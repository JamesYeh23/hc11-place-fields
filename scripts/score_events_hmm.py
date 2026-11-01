"""Step 10d: score sleep events with the receptive-field-free model.

    uv run python scripts/score_events_hmm.py Achilles_11012013
    uv run python scripts/score_events_hmm.py Achilles_11012013 --epochs PRE

Phase 1 scored each candidate event against place fields built from running.
This scores the same events against the **states** of the hidden Markov model
fitted in steps 10b and 10c, which never saw a position. Two orderings of those
states are used and both are always reported:

* **by position** -- states sorted by the position they map to, a map built
  from running behaviour. Comparable to Phase 1, and not receptive-field-free.
* **by the model** -- states sorted by the Fiedler vector of the symmetrised
  transition matrix, with no behaviour anywhere. This is Claim B, and it is the
  one that would mean a replay analysis can be done without place fields.

Four things make the numbers here readable, and three of them come from step
10c's results rather than from this script:

1. **Every score is read against the same event scored by a model fitted to
   shuffled training data** (D10.12), per session, never against zero. On
   Cicero_09012014, Cicero_09102014 and Gatsby_08022013 that control reached
   16-30 % overlap with the real fits on running, so a small margin there means
   nothing.
2. **Two sessions are marked vacuous in advance** (D10.14). Claim B on
   *running* -- where the structure is certainly present -- falls to 0.08 and
   0.06 on Gatsby_08022013 and Cicero_09012014, the two smallest ensembles. A
   null result on those sessions says nothing about sleep, and is reported as
   uninformative rather than as a negative, as in Phase 3.
3. **All five restarts, never the best.** Each event carries the median, min
   and max over the five fits, as everywhere else in this phase.
4. The state model is fitted on **every** running bout, not a training fold:
   sleep is the held-out data here, so there is nothing to hold out.

The scoring machinery is Phase 1's, unchanged -- `sequences.score_event`, the
same 20 ms bins, the same two within-event shuffles, the same maximisation over
circular offset applied identically to observation and shuffles. Only the
posterior it is handed is different.
"""

from __future__ import annotations

import argparse
import json
import logging
import time
from datetime import datetime, timezone

import numpy as np
import pandas as pd

from fit_hmm import roll_cells, running_sequences  # noqa: I001  (sibling script)
from hc11 import config, hmm, sequences
from place_fields import prepare

log = logging.getLogger("hmm-events")

EPOCHS = ("PRE", "POST")
ORDERINGS = ("position", "spectral")
DEFAULT_RESTARTS = 5
#: Shuffle kinds. Phase 1 reports both and they agree; the second is run on the
#: real model's first restart only, as a cross-check rather than a second grid.
PRIMARY_SHUFFLE = "column_cycle"
CROSS_CHECK_SHUFFLE = "time_permutation"


def fitted_models(session, track, running, pyr, n_states, bin_s, restarts):
    """``{(kind, seed): (model, order_by_position, order_by_model)}``.

    Fitted on every running bout. The shuffled models are the D10.12 control:
    the same fit on a circular shift of each cell along the whole running
    timeline, which is the one this project uses and not the within-bout roll
    that failed (METHODS, silent failure 8).
    """
    counts, positions = running_sequences(session, running, pyr, bin_s)
    out = {}
    for kind in ("real", "shuffled"):
        for seed in range(restarts):
            fit_on = (counts if kind == "real"
                      else roll_cells(counts, np.random.default_rng(5000 + seed)))
            model = hmm.fit(fit_on, n_states, bin_s, seed=seed)
            gamma = np.vstack([hmm.posterior(model, c) for c in counts])
            mapped = hmm.state_positions(gamma, np.concatenate(positions),
                                         n_states, circular=track.is_circular,
                                         extent=track.extent)
            out[(kind, seed)] = {
                "model": model,
                "position": np.argsort(np.nan_to_num(mapped, nan=np.inf)),
                "spectral": hmm.spectral_order(model.transition)}
    return out, counts


def read_retrying(path, attempts: int = 4, pause: float = 0.5) -> str:
    """Read a file, retrying a transient OSError rather than losing its contents.

    The repository lives in an iCloud-synced folder, where a file that was just
    deleted or written can briefly answer ``exists()`` and then fail the read
    with ``[Errno 35] Resource deadlock avoided``. Swallowing that and carrying
    on with an empty dictionary would silently drop every other session's
    provenance stamp, so this retries and then refuses.
    """
    for attempt in range(attempts):
        try:
            return path.read_text()
        except OSError as error:
            if attempt == attempts - 1:
                raise SystemExit(
                    f"could not read {path} after {attempts} attempts ({error}). "
                    "The folder is iCloud-synced; wait a moment and re-run rather "
                    "than deleting the file, which would lose the other sessions."
                ) from error
            time.sleep(pause)
    raise AssertionError("unreachable")


def event_counts(session, start_s, end_s, pyr, bin_s):
    """``(kept_counts, bin_index, n_bins_total)`` for one event, empty bins dropped."""
    counts, _ = sequences.bin_event(session, start_s, end_s, pyr, bin_s)
    if counts.shape[1] == 0:
        return np.zeros((0, len(pyr))), np.empty(0, dtype=int), 0
    per_bin = counts.sum(axis=0)
    keep = per_bin > 0
    return counts[:, keep].T, np.flatnonzero(keep), counts.shape[1]


def score(model_entry, ordering, counts, bin_index, grid, track, bin_s, shuffle,
          n_shuffles, rng):
    """One event, one model, one ordering: Phase 1's scorer on an HMM posterior."""
    model = model_entry["model"]
    gamma = hmm.memoryless_posterior(model, counts, bin_s=bin_s)
    event = sequences.EventPosteriors(
        posteriors={"state": gamma[:, model_entry[ordering]]},
        bin_index=bin_index, spikes_per_bin=counts.sum(axis=1),
        n_bins_total=len(bin_index), n_cells_active=int((counts.sum(axis=0) > 0).sum()))
    return sequences.score_event(event, grid, bin_s, circular=track.is_circular,
                                 shuffle=shuffle, n_shuffles=n_shuffles, rng=rng)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("session", nargs="?", default="Achilles_11012013")
    parser.add_argument("--epochs", nargs="*", default=list(EPOCHS))
    parser.add_argument("--restarts", type=int, default=DEFAULT_RESTARTS)
    parser.add_argument("--limit", type=int, default=None,
                        help="score only the first N events (timing runs)")
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(message)s")

    params = config.load_params()
    name = args.session
    fit_bin_s = params["decoding"]["time_bin_s"]
    bin_s = params["sequences"]["time_bin_s"]
    min_bins = params["sequences"]["min_bins"]
    n_shuffles = params["sequences"]["n_shuffles"]
    seed0 = params["sequences"]["shuffle_seed"]

    chosen = pd.read_csv(config.results_dir() / "hmm_state_choice.csv"
                         ).set_index("session")
    n_states = int(chosen.loc[name, "m_by_claim_a"])
    session, track, running, _ = prepare(name, params)
    pyr = np.asarray(session.units("pyr"))
    models, _ = fitted_models(session, track, running, pyr, n_states, fit_bin_s,
                              args.restarts)
    grid = np.arange(n_states) * track.extent / n_states
    log.info("%s: m = %d, %d pyramidal cells, %d models (%d real, %d shuffled)",
             name, n_states, len(pyr), len(models), args.restarts, args.restarts)

    candidates = pd.read_csv(config.results_dir() / "candidate_events.csv")
    rows = []
    for epoch in args.epochs:
        events = candidates[(candidates.session == name) & (candidates.epoch == epoch)]
        phase1 = pd.read_csv(
            config.results_dir() / "sequences" / f"{name}_{epoch}.csv")
        phase1_scored = set(np.round(phase1[phase1.scored].peak.to_numpy(), 6))
        if args.limit:
            events = events.head(args.limit)
        started, n_scored = time.time(), 0
        for position, event in enumerate(events.itertuples()):
            counts, bin_index, n_total = event_counts(session, event.start, event.end,
                                                      pyr, bin_s)
            if len(bin_index) < min_bins:
                continue
            n_scored += 1
            common = {"session": name, "epoch": epoch, "peak": event.peak,
                      "duration_ms": event.duration_ms, "n_bins_total": n_total,
                      "n_bins_used": len(bin_index),
                      "n_active_pyr": int((counts.sum(axis=0) > 0).sum()),
                      "scored_in_phase1": round(event.peak, 6) in phase1_scored}
            for ordering in ORDERINGS:
                for kind in ("real", "shuffled"):
                    values = []
                    for seed in range(args.restarts):
                        rng = np.random.default_rng(seed0 + 1000 * seed + position)
                        result = score(models[(kind, seed)], ordering, counts,
                                       bin_index, grid, track, bin_s,
                                       PRIMARY_SHUFFLE, n_shuffles, rng)
                        if result is not None:
                            values.append((result.rz, result.r))
                    if not values:
                        continue
                    rz = np.array([v[0] for v in values])
                    r = np.array([v[1] for v in values])
                    row = {**common, "ordering": ordering, "kind": kind,
                           "rz_median": float(np.median(rz)),
                           "rz_min": float(rz.min()), "rz_max": float(rz.max()),
                           "r_median": float(np.median(r)),
                           "n_restarts": len(values)}
                    if kind == "real":
                        rng = np.random.default_rng(seed0 + 99 + position)
                        cross = score(models[("real", 0)], ordering, counts,
                                      bin_index, grid, track, bin_s,
                                      CROSS_CHECK_SHUFFLE, n_shuffles, rng)
                        row["rz_time_permutation_seed0"] = (
                            float(cross.rz) if cross is not None else np.nan)
                    rows.append(row)
        log.info("  %s: %d of %d candidate events scorable (>= %d bins), %.0f s",
                 epoch, n_scored, len(events), min_bins, time.time() - started)

    path = config.results_dir() / "hmm_event_scores.csv"
    fresh = pd.DataFrame(rows)
    if path.exists():
        kept = pd.read_csv(path)
        kept = kept[~((kept.session == name) & (kept.epoch.isin(args.epochs)))]
        fresh = pd.concat([kept, fresh], ignore_index=True)
    fresh.sort_values(["session", "epoch", "peak", "ordering", "kind"]).round(5).to_csv(
        path, index=False)

    stamp_path = config.results_dir() / "hmm_event_scores.json"
    stamps = json.loads(read_retrying(stamp_path)) if stamp_path.exists() else {}
    stamps[name] = {"script": "scripts/score_events_hmm.py",
                    "git_commit": config.git_commit(),
                    "created": datetime.now(timezone.utc).isoformat(timespec="seconds"),
                    "n_states": n_states, "restarts": args.restarts,
                    "bin_s": bin_s, "fit_bin_s": fit_bin_s,
                    "n_shuffles": n_shuffles, "min_bins": min_bins,
                    "epochs": sorted(set(stamps.get(name, {}).get("epochs", []))
                                     | set(args.epochs))}
    stamp_path.write_text(json.dumps(stamps, indent=2, sort_keys=True) + "\n")
    log.info("wrote %s", path.relative_to(config.repo_root()))


if __name__ == "__main__":
    main()
