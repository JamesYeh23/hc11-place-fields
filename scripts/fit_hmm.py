"""Step 10b: fit the receptive-field-free decoder on running, and test Claim A.

    uv run python scripts/fit_hmm.py Achilles_11012013
    uv run python scripts/fit_hmm.py Achilles_11012013 --states 12 16 --restarts 5

Phase 4's method (Chen et al. 2016's DecodewoRF) is a hidden Markov model fitted
to binned population spiking with no behavioural measure anywhere in the fit.
One question decides whether it is worth building on: **do the states it
recovers correspond to positions on the track?** Chen et al. report a
correlation of 0.73-0.79 between the recovered state trajectory and the animal's
run trajectory. If we cannot get near that on the best-powered session we have,
the phase stops here -- the 10a go/no-go.

Two claims are kept apart from this milestone onward, because only the second is
about replay:

* **Claim A -- the states correspond to track positions.** Map each state to a
  position using *training* behaviour only, decode held-out running, compare.
  Behaviour is used for assessment, never for fitting. On a ring the statistic
  is the circular correlation, because a Pearson correlation there depends on
  where the ring is cut; the median decoding error is reported beside it against
  a permuted-map chance level.
* **Claim B -- events traverse states in orderly sequences.** No behaviour at
  all: order states by the Fiedler vector of the fitted transition matrix and
  take the weighted correlation across that ordering. The **gap** between the
  two orderings is reported wherever both exist.

Three things are reported that a single number would hide:

1. **All restarts, not the best.** EM has local optima here; on synthetic data
   three of five restarts recover the generating rates and two do not.
2. **Stickiness.** A hidden Markov model fitted to population spiking will
   learn a sticky transition matrix whether or not the states mean anything,
   and a sticky model makes any state path look smooth -- which is the failure
   mode 10a identified for this phase. The fitted dwell time is reported against
   two baselines: how long the animal actually spends crossing one position bin,
   and what the same fit produces from **shuffled** data.
3. **A fit on shuffled training data** (each cell's counts rolled independently
   within each running bout, preserving rate and burstiness, destroying
   cross-cell alignment). This is control 1 from the 10a scope, and it is run
   and read before the real numbers mean anything.
"""

from __future__ import annotations

import argparse
import json
import logging
from datetime import datetime, timezone

import numpy as np
import pandas as pd

from hc11 import config, hmm, sequences
from place_fields import prepare  # noqa: I001  (sibling script)

log = logging.getLogger("hmm")

DEFAULT_STATES = (8, 12, 16, 20, 24, 28)
DEFAULT_RESTARTS = 5
N_FOLDS = 3
N_MAP_PERMUTATIONS = 20
CLAIM_B_WINDOW_S = 2.0


def running_sequences(session, running, cluster_ids, bin_s):
    """One count matrix and position vector per running bout, never concatenated."""
    spikes = {int(c): session.spikes_for(int(c)) for c in cluster_ids}
    counts_out, positions_out = [], []
    for start, end in running.intervals:
        edges = np.arange(start, end + bin_s, bin_s)
        if len(edges) < 3:
            continue
        counts = np.empty((len(edges) - 1, len(cluster_ids)))
        for j, c in enumerate(cluster_ids):
            counts[:, j] = np.histogram(spikes[int(c)], bins=edges)[0]
        centres = 0.5 * (edges[:-1] + edges[1:])
        position = np.interp(centres, session.position_t, session.position_1d)
        if not np.isfinite(position).all():
            continue
        counts_out.append(counts)
        positions_out.append(position)
    return counts_out, positions_out


def roll_cells(sequences_in, rng):
    """Each cell circularly shifted along the whole running timeline, then re-split.

    The shift is across *all* bouts, not within one. Rolling inside a bout is
    not a control here and the first run of this script proved it: bouts run a
    median of six bins, the animal is in one place for all of them, so each
    cell's rate within a bout identifies the position whatever the order --
    and a model fitted to within-bout rolls decoded position *better* than the
    real fit. Shifting along the concatenated timeline is the project's
    existing idiom (the place-cell shuffle of step 4): it keeps every cell's
    rate and burst structure and destroys the relationship between a cell and
    where the animal was.
    """
    lengths = [len(c) for c in sequences_in]
    pooled = np.vstack(sequences_in)
    rolled = np.empty_like(pooled)
    for j in range(pooled.shape[1]):
        rolled[:, j] = np.roll(pooled[:, j], int(rng.integers(len(pooled))))
    out, start = [], 0
    for length in lengths:
        out.append(rolled[start:start + length])
        start += length
    return out


def folds(n_sequences, n_folds):
    """Contiguous blocks of bouts, so a fold is a stretch of the epoch."""
    edges = np.linspace(0, n_sequences, n_folds + 1).astype(int)
    for k in range(n_folds):
        test = np.arange(edges[k], edges[k + 1])
        train = np.setdiff1d(np.arange(n_sequences), test)
        if len(test) and len(train):
            yield train, test


def stickiness(model) -> dict:
    """Self-transition mass and the dwell time it implies.

    Expected dwell in a state is ``bin / (1 - P_ii)`` bins. A model whose states
    are positions should dwell about as long as the animal takes to cross one,
    and a model that has learned nothing but "the population rate is
    autocorrelated" will dwell much longer.
    """
    diagonal = np.diag(model.transition)
    dwell_bins = 1.0 / np.maximum(1.0 - diagonal, 1e-9)
    return {"median_self_transition": float(np.median(diagonal)),
            "median_dwell_bins": float(np.median(dwell_bins)),
            "median_dwell_s": float(np.median(dwell_bins) * model.bin_s),
            "max_self_transition": float(diagonal.max())}


def crossing_time_s(session, running, track, bin_s) -> float:
    """How long the animal actually spends crossing one position bin while running.

    The behavioural baseline for stickiness: if a state is a place on the track,
    this is roughly how long the model should stay in one.
    """
    width = track.extent / track.n_bins
    speeds = []
    for start, end in running.intervals:
        mask = (session.position_t >= start) & (session.position_t <= end)
        position = session.position_1d[mask]
        time = session.position_t[mask]
        finite = np.isfinite(position)
        position, time = position[finite], time[finite]
        if len(position) < 3:
            continue
        step = np.abs(np.diff(position))
        if track.is_circular:
            step = np.minimum(step, track.extent - step)
        speeds.append(step / np.maximum(np.diff(time), 1e-9))
    if not speeds:
        return float("nan")
    return float(width / np.median(np.concatenate(speeds)))


def claim_a(model, train_counts, train_positions, test_counts, test_positions,
            track, rng) -> dict:
    """Do the states correspond to positions? Behaviour only for assessment."""
    gamma = np.vstack([hmm.posterior(model, c) for c in train_counts])
    mapped = hmm.state_positions(gamma, np.concatenate(train_positions),
                                 model.n_states, circular=track.is_circular,
                                 extent=track.extent)
    states = np.concatenate([hmm.posterior(model, c).argmax(axis=1)
                             for c in test_counts])
    truth = np.concatenate(test_positions)
    decoded = mapped[states]
    keep = np.isfinite(decoded) & np.isfinite(truth)
    decoded, truth, states = decoded[keep], truth[keep], states[keep]
    if decoded.size < 10:
        return {"claim_a_r": np.nan, "median_error_cm": np.nan,
                "chance_error_cm": np.nan, "n_test_bins": int(decoded.size)}

    def error(values):
        if track.is_circular:
            return float(np.median(hmm.circular_distance(values, truth, track.extent)))
        return float(np.median(np.abs(values - truth)))

    chance = [error(mapped[rng.permutation(model.n_states)][states])
              for _ in range(N_MAP_PERMUTATIONS)]
    correlation = (hmm.circular_correlation(decoded, truth, track.extent)
                   if track.is_circular
                   else float(np.corrcoef(decoded, truth)[0, 1]))
    return {"claim_a_r": correlation, "median_error_cm": 100 * error(decoded),
            "chance_error_cm": 100 * float(np.median(chance)),
            "n_test_bins": int(decoded.size), "mapped": mapped}


def claim_b(model, windows, track, mapped) -> dict:
    """Weighted correlation of held-out running windows, under both orderings.

    The orderings differ only in how the posterior's columns are permuted, so
    the coordinate grid is identical and the two scores are comparable.
    """
    by_position = np.argsort(np.nan_to_num(mapped, nan=np.inf))
    by_transition = hmm.spectral_order(model.transition)
    grid = np.arange(model.n_states) * track.extent / model.n_states
    scores = {"position": [], "spectral": []}
    for counts, fine_bin_s in windows:
        # Memoryless: a transition matrix fitted at 250 ms is wrong at 20 ms.
        gamma = hmm.memoryless_posterior(model, counts, bin_s=fine_bin_s)
        times = np.arange(len(counts)) * fine_bin_s
        for label, order in (("position", by_position), ("spectral", by_transition)):
            r, _ = sequences.weighted_correlation(
                gamma[:, order], grid, times, track.is_circular)
            scores[label].append(abs(float(r)))
    out = {f"claim_b_{k}": float(np.median(v)) if v else np.nan
           for k, v in scores.items()}
    out["claim_b_gap"] = out["claim_b_position"] - out["claim_b_spectral"]
    out["n_claim_b_windows"] = len(windows)
    return out


def fine_windows(session, running, cluster_ids, test_indices, bin_s, window_s):
    """Held-out running windows at the scoring bin, for Claim B."""
    spikes = {int(c): session.spikes_for(int(c)) for c in cluster_ids}
    out = []
    for i in test_indices:
        if i >= len(running.intervals):
            continue
        start, end = running.intervals[i]
        for begin in np.arange(start, end - window_s, window_s):
            edges = np.arange(begin, begin + window_s + bin_s, bin_s)
            counts = np.empty((len(edges) - 1, len(cluster_ids)))
            for j, c in enumerate(cluster_ids):
                counts[:, j] = np.histogram(spikes[int(c)], bins=edges)[0]
            if counts.sum() > 0:
                out.append((counts, bin_s))
    return out


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("session", nargs="?", default="Achilles_11012013")
    parser.add_argument("--states", type=int, nargs="*", default=list(DEFAULT_STATES))
    parser.add_argument("--restarts", type=int, default=DEFAULT_RESTARTS)
    parser.add_argument("--no-shuffle-control", action="store_true")
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(message)s")

    params = config.load_params()
    name = args.session
    bin_s = params["decoding"]["time_bin_s"]
    fine_bin_s = params["sequences"]["time_bin_s"]
    session, track, running, laps = prepare(name, params)
    pyr = np.asarray(session.units("pyr"))
    counts, positions = running_sequences(session, running, pyr, bin_s)
    total_bins = sum(len(c) for c in counts)
    crossing = crossing_time_s(session, running, track, bin_s)
    log.info("%s: %d pyramidal cells, %d running bouts, %d bins of %g s, %d spikes",
             name, len(pyr), len(counts), total_bins, bin_s,
             int(sum(c.sum() for c in counts)))
    log.info("track %.2f m, %d position bins, circular=%s; bout length median %d bins, "
             "max %d", track.extent, track.n_bins, track.is_circular,
             int(np.median([len(c) for c in counts])), max(len(c) for c in counts))
    log.info("the animal crosses one position bin in %.2f s (%.1f bins of %g s) -- the "
             "behavioural baseline for how long a position-like state should last",
             crossing, crossing / bin_s, bin_s)

    rows = []
    for n_states in args.states:
        log.info("\n=== m = %d  (%d transition parameters against %d transitions)",
                 n_states, n_states * (n_states - 1), total_bins - len(counts))
        for fold_index, (train, test) in enumerate(folds(len(counts), N_FOLDS)):
            train_counts = [counts[i] for i in train]
            train_positions = [positions[i] for i in train]
            test_counts = [counts[i] for i in test]
            test_positions = [positions[i] for i in test]
            windows = fine_windows(session, running, pyr, test, fine_bin_s,
                                   CLAIM_B_WINDOW_S)
            for seed in range(args.restarts):
                for kind in ("real", "shuffled"):
                    if kind == "shuffled" and (args.no_shuffle_control
                                               or fold_index != 0):
                        continue
                    fit_on = (train_counts if kind == "real" else
                              roll_cells(train_counts,
                                         np.random.default_rng(5000 + seed)))
                    model = hmm.fit(fit_on, n_states, bin_s, seed=seed)
                    rng = np.random.default_rng(1000 + seed)
                    row = {"session": name, "kind": kind, "n_states": n_states,
                           "fold": fold_index, "seed": seed,
                           "train_log_likelihood": model.log_likelihood,
                           "held_out_log_likelihood": hmm.sequence_log_likelihood(
                               model, test_counts),
                           "n_iterations": model.n_iterations,
                           "converged": model.converged,
                           "n_train_bins": sum(len(c) for c in train_counts),
                           "crossing_time_s": crossing}
                    row.update(stickiness(model))
                    measured = claim_a(model, train_counts, train_positions,
                                       test_counts, test_positions, track, rng)
                    mapped = measured.pop("mapped", None)
                    row.update(measured)
                    if mapped is not None:
                        row.update(claim_b(model, windows, track, mapped))
                    rows.append(row)
            here = [r for r in rows if r["n_states"] == n_states
                    and r["fold"] == fold_index and r["kind"] == "real"]
            log.info("  fold %d  Claim A r over %d restarts: %s", fold_index, len(here),
                     "  ".join(f"{r['claim_a_r']:+.2f}" for r in here))
            log.info("          error %.0f-%.0f cm (chance %.0f) | Claim B position "
                     "%.2f-%.2f vs spectral %.2f-%.2f | dwell %.2f-%.2f s",
                     min(r["median_error_cm"] for r in here),
                     max(r["median_error_cm"] for r in here),
                     np.median([r["chance_error_cm"] for r in here]),
                     min(r["claim_b_position"] for r in here),
                     max(r["claim_b_position"] for r in here),
                     min(r["claim_b_spectral"] for r in here),
                     max(r["claim_b_spectral"] for r in here),
                     min(r["median_dwell_s"] for r in here),
                     max(r["median_dwell_s"] for r in here))
            shuffled = [r for r in rows if r["n_states"] == n_states
                        and r["fold"] == fold_index and r["kind"] == "shuffled"]
            if shuffled:
                log.info("          SHUFFLED control: Claim A r %s | dwell %.2f-%.2f s",
                         "  ".join(f"{r['claim_a_r']:+.2f}" for r in shuffled),
                         min(r["median_dwell_s"] for r in shuffled),
                         max(r["median_dwell_s"] for r in shuffled))

    # Merge on (session, m): the grid is run in pieces to fit the time budget,
    # so a later call must not delete the m values an earlier one computed.
    path = config.results_dir() / "hmm_running_fit.csv"
    fresh = pd.DataFrame(rows)
    if path.exists():
        kept = pd.read_csv(path)
        kept = kept[~((kept.session == name) & (kept.n_states.isin(args.states)))]
        fresh = pd.concat([kept, fresh], ignore_index=True)
    fresh.sort_values(["session", "n_states", "kind", "fold", "seed"]).round(5).to_csv(
        path, index=False)
    stamp = {"script": "scripts/fit_hmm.py", "session": name,
             "git_commit": config.git_commit(),
             "created": datetime.now(timezone.utc).isoformat(timespec="seconds"),
             "states": list(args.states), "restarts": args.restarts, "folds": N_FOLDS,
             "bin_s": bin_s, "fine_bin_s": fine_bin_s,
             "claim_b_window_s": CLAIM_B_WINDOW_S,
             "crossing_time_s": crossing, "chen_reference_r": [0.73, 0.79]}
    (config.results_dir() / "hmm_running_fit.json").write_text(
        json.dumps(stamp, indent=2) + "\n")
    log.info("\nwrote %s", (config.results_dir() / "hmm_running_fit.csv").relative_to(
        config.repo_root()))


if __name__ == "__main__":
    main()
