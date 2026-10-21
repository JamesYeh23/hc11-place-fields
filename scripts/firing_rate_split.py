"""Step 9b: the paper's Fig. 4D — do slow-firing cells gain replay and fast ones not?

    uv run python scripts/firing_rate_split.py
    uv run python scripts/firing_rate_split.py Achilles_10252013 --random-splits 20

Grosmark & Buzsáki report that place cells divided into equal subgroups by
off-line (sleep) firing rate behave differently: "Slow-firing, but not
fast-firing, neurons increased their contribution to neuronal sequences from
PRE to POST sleep." That claim needs no per-cell contribution score, invents
almost no parameters, and costs a minute to test, which is why it runs before
anything built on PCC.

The test. Each session's place cells are split at the median of their off-line
firing rate into equal halves (`hc11.subsets`). A decoding template is rebuilt
from each half alone, the same candidate events are re-decoded and re-scored
against it, and PRE is compared with POST **within each half**. The claim
predicts a POST-minus-PRE gain in the slow half and little or none in the fast
half.

Three things the measurement has to control, all of which would otherwise
produce a plausible number:

* **The event set must not move with the template.** Half an ensemble leaves
  more empty 20 ms bins, so fewer events clear the five-bin floor. The
  comparison therefore runs on the events scorable under *both* halves and
  under the full ensemble, so the two halves are scored on the same events.
* **A difference between halves is not evidence by itself.** Any split of an
  ensemble into halves produces some difference. `--random-splits` repeats the
  whole measurement on equal halves drawn at random, giving the spread of
  "difference between two arbitrary halves" against which the rate split is
  read.
* **The rate definition is ours, not the paper's.** `subsets.offline_rate_window`
  selects the whole sleep epochs or scored non-REM only; both splits are built
  and their agreement is reported, so a result that depends on that choice
  cannot be mistaken for one that does not.

Writes results/firing_rate_split.csv, results/firing_rate_split_cells.csv,
results/firing_rate_split_null.csv and results/firing_rate_split.png.
"""

from __future__ import annotations

import argparse
import json
import logging
from datetime import datetime, timezone

import numpy as np
import pandas as pd

from hc11 import config, plots, sequences, subsets
from hc11.mua import state_intervals_for_epoch
from hc11.results import write_session_table
from place_fields import prepare  # noqa: I001  (sibling script)
from replay_comparison import compare  # noqa: I001
from score_events import direction_tuning, place_cell_ids  # noqa: I001

log = logging.getLogger("subsets")

SESSIONS = list(config.LFP_SESSIONS)
SLEEP_EPOCHS = ("PRE", "POST")


def offline_rates(session, cluster_ids, params, window: str) -> np.ndarray:
    """Firing rate of every cell over the sleep epochs, or over scored non-REM."""
    if window == "epochs":
        intervals = np.array([session.epoch(e) for e in SLEEP_EPOCHS], dtype=float)
    elif window == "nrem":
        pieces = [state_intervals_for_epoch(session, e, params) for e in SLEEP_EPOCHS]
        intervals = np.vstack([p for p in pieces if p.size])
    else:
        raise ValueError(f"unknown offline_rate_window {window!r}")
    return np.array([subsets.rate_in_intervals(session, int(c), intervals)
                     for c in cluster_ids])


def score_with(session, track, events, cells, tuning, params) -> pd.DataFrame:
    """Score ``events`` against one template. One row per event, NaN if unscorable."""
    cfg = params["sequences"]
    bin_s = cfg["time_bin_s"]
    rng = np.random.default_rng(cfg["shuffle_seed"])
    rows = []
    for row in events.itertuples():
        event = sequences.event_posteriors(
            session, row.start, row.end, cells, tuning, bin_s,
            rate_floor_hz=cfg["rate_floor_hz"],
        )
        record = {"peak": row.peak, "scored": event.n_bins >= cfg["min_bins"],
                  "n_bins_used": event.n_bins}
        if record["scored"]:
            for kind in sequences.SHUFFLES:
                score = sequences.score_event(
                    event, track.bin_centers, bin_s, circular=track.is_circular,
                    shuffle=kind, n_shuffles=cfg["n_shuffles"], rng=rng,
                )
                record[f"rz_{kind}"] = score.rz
        rows.append(record)
    return pd.DataFrame(rows)


def half_scores(session, track, running, laps, events_by_epoch, cells, params):
    """``{epoch: DataFrame}`` for one subset of cells, or None if it has no template."""
    tuning = direction_tuning(session, track, running, laps, np.asarray(cells), params)
    if not tuning:
        return None
    return {epoch: score_with(session, track, events, np.asarray(cells), tuning, params)
            for epoch, events in events_by_epoch.items()}


def contrast(low, high, matched, kind: str) -> dict:
    """POST-minus-PRE in each half, on the same events, and the gap between them.

    An event whose score is undefined in either half (a degenerate posterior
    gives a zero-width null and so no z-score) is dropped from both halves, not
    from one: see :func:`hc11.subsets.finite_in_both`.
    """
    usable, dropped = {}, 0
    for epoch in SLEEP_EPOCHS:
        a = low[epoch].loc[matched[epoch], f"rz_{kind}"].to_numpy(dtype=float)
        b = high[epoch].loc[matched[epoch], f"rz_{kind}"].to_numpy(dtype=float)
        keep = subsets.finite_in_both(a, b)
        dropped += int((~keep).sum())
        usable[epoch] = (a[keep], b[keep])
    out = {"n_dropped_undefined": dropped}
    for label, index in (("low", 0), ("high", 1)):
        out[label] = compare(usable["PRE"][index], usable["POST"][index])
    out["gap"] = out["low"]["difference"] - out["high"]["difference"]
    return out


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("sessions", nargs="*", default=SESSIONS)
    parser.add_argument("--random-splits", type=int, default=None,
                        help="override subsets.n_random_splits (0 skips the control)")
    parser.add_argument("--no-figures", action="store_true")
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(message)s")

    params = config.load_params()
    cfg = params["subsets"]
    n_random = cfg["n_random_splits"] if args.random_splits is None else args.random_splits
    out_root = config.results_dir()
    candidates = pd.read_csv(out_root / "candidate_events.csv")
    fields_table = pd.read_csv(out_root / "place_fields_summary.csv")

    rows, cell_rows, null_rows = [], [], []
    for name in args.sessions:
        log.info("\n=== %s", name)
        session, track, running, laps = prepare(name, params)
        cells = place_cell_ids(fields_table, name)
        # Only events the full ensemble could score can enter the matched set, so
        # those are the only ones worth re-scoring: a four-fold saving that
        # changes no number, since `matched` intersects with them anyway.
        events_by_epoch, scored_full = {}, {}
        for epoch in SLEEP_EPOCHS:
            full = pd.read_csv(out_root / "sequences" / f"{name}_{epoch}.csv")
            keep = full[full.scored].peak.to_numpy()
            scored_full[epoch] = int(keep.size)
            here = candidates[(candidates.session == name)
                              & (candidates.epoch == epoch)]
            events_by_epoch[epoch] = here[here.peak.isin(keep)].reset_index(drop=True)
            if len(events_by_epoch[epoch]) != keep.size:
                raise SystemExit(
                    f"{name} {epoch}: {keep.size} scored events in the step-7 table but "
                    f"{len(events_by_epoch[epoch])} matched in candidate_events.csv")

        rates = offline_rates(session, cells, params, cfg["offline_rate_window"])
        split = subsets.equal_halves(rates, cells)
        other_window = "nrem" if cfg["offline_rate_window"] == "epochs" else "epochs"
        other = subsets.equal_halves(
            offline_rates(session, cells, params, other_window), cells)
        agreement = subsets.split_agreement(split, other)
        slow_rate = np.median(rates[np.isin(cells, split["low"])])
        fast_rate = np.median(rates[np.isin(cells, split["high"])])
        log.info("%d place cells -> %d slow / %d fast (%d dropped at the median)",
                 len(cells), len(split["low"]), len(split["high"]), len(split["dropped"]))
        log.info("  median off-line rate: slow %.2f Hz, fast %.2f Hz "
                 "(paper: 0.39 and 1.12 Hz); %.0f%% of cells land in the same half "
                 "under the %s definition", slow_rate, fast_rate, 100 * agreement,
                 other_window)
        for cluster_id, rate in zip(cells, rates, strict=True):
            cell_rows.append({
                "session": name, "cluster_id": int(cluster_id), "offline_rate_hz": rate,
                "half": ("slow" if cluster_id in split["low"]
                         else "fast" if cluster_id in split["high"] else "dropped"),
            })

        low = half_scores(session, track, running, laps, events_by_epoch,
                          split["low"], params)
        high = half_scores(session, track, running, laps, events_by_epoch,
                           split["high"], params)
        if low is None or high is None:
            log.warning("  no usable template for one half -- skipping %s", name)
            continue

        matched, counts = {}, {}
        for epoch in SLEEP_EPOCHS:
            both = low[epoch].scored.to_numpy() & high[epoch].scored.to_numpy()
            matched[epoch] = np.flatnonzero(both)
            counts[epoch] = {
                "full": scored_full[epoch],
                "slow": int(low[epoch].scored.sum()),
                "fast": int(high[epoch].scored.sum()),
                "matched": int(matched[epoch].size),
            }
        log.info("  scorable events  PRE %(full)d full / %(slow)d slow / %(fast)d fast "
                 "-> %(matched)d matched", counts["PRE"])
        log.info("  scorable events POST %(full)d full / %(slow)d slow / %(fast)d fast "
                 "-> %(matched)d matched", counts["POST"])

        for kind in sequences.SHUFFLES:
            result = contrast(low, high, matched, kind)
            if result["n_dropped_undefined"]:
                log.info("  %-18s %d event(s) dropped: no z-score in one half",
                         kind, result["n_dropped_undefined"])
            rows.append({
                "session": name, "shuffle": kind,
                "n_cells": len(cells), "n_per_half": len(split["low"]),
                "slow_rate_hz": slow_rate, "fast_rate_hz": fast_rate,
                "split_agreement": agreement,
                "n_matched_pre": counts["PRE"]["matched"],
                "n_matched_post": counts["POST"]["matched"],
                "n_dropped_undefined": result["n_dropped_undefined"],
                "n_scored_pre": result["low"]["n_pre"],
                "n_scored_post": result["low"]["n_post"],
                "slow_median_pre": result["low"]["median_pre"],
                "slow_median_post": result["low"]["median_post"],
                "slow_difference": result["low"]["difference"],
                "slow_rank_biserial": result["low"]["rank_biserial"],
                "slow_p": result["low"]["p"],
                "fast_median_pre": result["high"]["median_pre"],
                "fast_median_post": result["high"]["median_post"],
                "fast_difference": result["high"]["difference"],
                "fast_rank_biserial": result["high"]["rank_biserial"],
                "fast_p": result["high"]["p"],
                "gap_slow_minus_fast": result["gap"],
            })
            log.info("  %-18s slow %+.3f -> %+.3f (diff %+.3f, rb %+.3f, p %.3g)",
                     kind, result["low"]["median_pre"], result["low"]["median_post"],
                     result["low"]["difference"], result["low"]["rank_biserial"],
                     result["low"]["p"])
            log.info("  %-18s fast %+.3f -> %+.3f (diff %+.3f, rb %+.3f, p %.3g)",
                     "", result["high"]["median_pre"], result["high"]["median_post"],
                     result["high"]["difference"], result["high"]["rank_biserial"],
                     result["high"]["p"])
            log.info("  %-18s gap (slow - fast) %+.3f", "", result["gap"])

        if n_random:
            rng = np.random.default_rng(cfg["random_split_seed"])
            gaps = {kind: [] for kind in sequences.SHUFFLES}
            for draw in range(n_random):
                random_split = subsets.random_halves(cells, rng)
                a = half_scores(session, track, running, laps, events_by_epoch,
                                random_split["low"], params)
                b = half_scores(session, track, running, laps, events_by_epoch,
                                random_split["high"], params)
                if a is None or b is None:
                    continue
                random_matched = {
                    epoch: np.flatnonzero(a[epoch].scored.to_numpy()
                                          & b[epoch].scored.to_numpy())
                    for epoch in SLEEP_EPOCHS
                }
                for kind in sequences.SHUFFLES:
                    result = contrast(a, b, random_matched, kind)
                    gaps[kind].append(result["gap"])
                    null_rows.append({"session": name, "shuffle": kind, "draw": draw,
                                      "gap": result["gap"],
                                      "difference_a": result["low"]["difference"],
                                      "difference_b": result["high"]["difference"]})
            for row in rows:
                if row["session"] != name:
                    continue
                draws = np.abs(np.asarray(gaps[row["shuffle"]], dtype=float))
                draws = draws[np.isfinite(draws)]
                observed = abs(row["gap_slow_minus_fast"])
                row["n_random_splits"] = int(draws.size)
                row["random_gap_p95"] = (float(np.percentile(draws, 95))
                                         if draws.size else np.nan)
                row["p_vs_random_split"] = subsets.empirical_p(observed, draws)
                log.info("  %-18s |gap| %.3f against %d usable random splits: "
                         "95th percentile %.3f, p = %.3f", row["shuffle"], observed,
                         draws.size, row["random_gap_p95"], row["p_vs_random_split"])

    if not rows:
        return
    table = write_session_table(out_root / "firing_rate_split.csv", rows,
                                sessions=args.sessions, round_to=5)
    write_session_table(out_root / "firing_rate_split_cells.csv", cell_rows,
                        sessions=args.sessions, round_to=5)
    if null_rows:
        write_session_table(out_root / "firing_rate_split_null.csv", null_rows,
                            sessions=args.sessions, round_to=5)
    stamp = {
        "script": "scripts/firing_rate_split.py",
        "git_commit": config.git_commit(),
        "created": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "sessions": list(args.sessions), "omitted": config.LFP_OMITTED,
        "subsets": cfg, "sequences": params["sequences"],
    }
    (out_root / "firing_rate_split.json").write_text(json.dumps(stamp, indent=2) + "\n")
    if not args.no_figures:
        plots.plot_firing_rate_split(
            table.to_dict("records"),
            pd.read_csv(out_root / "firing_rate_split_null.csv").to_dict("records")
            if (out_root / "firing_rate_split_null.csv").exists() else [],
            out_root / "firing_rate_split.png", stamp)
    log.info("\nwrote %s", (out_root / "firing_rate_split.csv").relative_to(
        config.repo_root()))


if __name__ == "__main__":
    main()
