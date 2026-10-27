"""Step 10a-bis: what the MUA detector's boundary threshold costs, and what it buys.

    uv run python scripts/boundary_threshold_sweep.py
    uv run python scripts/boundary_threshold_sweep.py Achilles_10252013 --no-figures

Our candidate events are short. D7.5 records the consequence -- only 6-23 % of
them clear the five-bin floor at 20 ms -- and traces it to the detector's
boundary threshold rather than to anything in the scoring. Chen et al. (2016)
score sleep events of 11.5 bins where ours run 7.2-8.0, and a receptive-field-free
decoder has less to work with per event than ours does.

This sweep asks whether the shortfall is fixable. The **peak** threshold stays
at 3 SD throughout, because that is the one number Chen et al. actually state;
only the **boundary** threshold moves, and it is ours (D10.1). Lowering it
extends each event outward from its peak until the population rate falls back
to the boundary, so events get longer by construction. The question is whether
the added time carries signal or is drift into baseline, and the sweep reports
four things that answer it:

* **spikes per used bin** and the **empty-bin fraction** -- if duration grows
  while these worsen, the extra time is empty;
* **active cells per event** -- real extension should recruit cells, not just
  seconds;
* the **ripple-matched fraction** -- an event that has drifted into baseline
  stops coinciding with a sharp-wave ripple;
* the **multi-shank support** of the matched ripples, which is the independent
  check D6.19 established.

**Nothing here adopts a threshold.** The output is a trade-off curve; the
decision, if there is one, belongs in its own entry with its own evidence.
"""

from __future__ import annotations

import argparse
import copy
import json
import logging
from datetime import datetime, timezone

import numpy as np
import pandas as pd

from hc11 import config, plots
from hc11.events import nearest_within
from hc11.io import load_session
from hc11.mua import detect_mua_events
from hc11.preprocess import truncate_to_session
from hc11.results import write_session_table
from score_events import place_cell_ids  # noqa: I001  (sibling script)

log = logging.getLogger("boundary")

SESSIONS = list(config.LFP_SESSIONS)
EPOCHS = ("PRE", "POST")
#: Boundary thresholds to sweep, in standard deviations above the state's own
#: baseline mean. 1.0 is the pipeline's current value; the grid runs down to
#: where the boundary stops binding at all. [PROVISIONAL]
GRID = (1.0, 0.75, 0.5, 0.25, 0.0, -0.25, -0.5, -1.0)
#: Match window to a sharp-wave ripple, as in the step-6c cross-check (D6.19).
RIPPLE_WINDOW_S = 0.05
#: Scoring bin and floor, from `sequences` -- the quantities this sweep is about.
BIN_S = 0.02
MIN_BINS = 5


def per_event_structure(session, events, cluster_ids, bin_s=BIN_S) -> dict:
    """Bins, used bins and active cells for every event, vectorised over spikes."""
    ids = np.asarray(cluster_ids)
    order = np.argsort(session.spike_times, kind="stable")
    times = session.spike_times[order]
    owners = session.spike_ids[order]
    keep = np.isin(owners, ids)
    times, owners = times[keep], owners[keep]

    n_events = len(events.start)
    total_bins = np.maximum(np.round((events.end - events.start) / bin_s), 1).astype(int)
    used_bins = np.zeros(n_events, dtype=int)
    active = np.zeros(n_events, dtype=int)
    spikes = np.zeros(n_events, dtype=int)
    lo = np.searchsorted(times, events.start, "left")
    hi = np.searchsorted(times, events.end, "right")
    for i, (a, b) in enumerate(zip(lo, hi, strict=True)):
        if b <= a:
            continue
        spikes[i] = b - a
        bins = np.floor((times[a:b] - events.start[i]) / bin_s).astype(int)
        used_bins[i] = np.unique(bins).size
        active[i] = np.unique(owners[a:b]).size
    return {"total_bins": total_bins, "used_bins": used_bins,
            "active": active, "spikes": spikes}


def ripple_match(events, ripples: pd.DataFrame) -> tuple[np.ndarray, np.ndarray]:
    """For each event: is there a ripple peak within the window, and its support.

    Support is the *count* of other shanks showing the same transient, not the
    fraction: the denominator is constant within a session, which is all this
    sweep compares, and `n_comparison_shanks` is absent for the two sessions
    detected before step 8 added it.
    """
    if ripples.empty or len(events.start) == 0:
        return (np.zeros(len(events.start), dtype=bool),
                np.full(len(events.start), np.nan))
    peaks = np.sort(ripples.peak.to_numpy())
    order = np.argsort(ripples.peak.to_numpy())
    support = ripples.n_supporting_shanks.to_numpy()[order]
    nearest = nearest_within(events.peak, peaks, RIPPLE_WINDOW_S)
    matched = nearest >= 0
    matched_support = np.where(matched, support[np.clip(nearest, 0, None)], np.nan)
    return matched, matched_support


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("sessions", nargs="*", default=SESSIONS)
    parser.add_argument("--no-figures", action="store_true")
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(message)s")

    params = config.load_params()
    results = config.results_dir()
    fields_table = pd.read_csv(results / "place_fields_summary.csv")
    ripple_table = pd.read_csv(results / "ripple_events.csv")
    committed = pd.read_csv(results / "mua_events.csv")
    rows = []
    for name in args.sessions:
        log.info("\n=== %s", name)
        session = load_session(name)
        if params["preprocess"]["truncate_to_session"]:
            session, _ = truncate_to_session(session, verbose=False)
        place = place_cell_ids(fields_table, name)
        pyr = np.asarray(session.units("pyr"))
        log.info("%d place cells, %d pyramidal cells", len(place), len(pyr))
        for epoch in EPOCHS:
            ripples = ripple_table[(ripple_table.session == name)
                                   & (ripple_table.epoch == epoch)]
            # Absent for the two sessions detected before step 8 added it; the
            # sweep compares thresholds within a session, so it is context, not
            # a denominator anything here divides by.
            counts = ripples.get("n_comparison_shanks")
            comparison_shanks = (float(counts.median())
                                 if counts is not None and counts.notna().any()
                                 else float("nan"))
            log.info("  %s  %-7s %7s %9s %6s %7s %10s %11s %8s %7s %8s",
                     epoch, "low SD", "events", "dur (ms)", "used", "scor.",
                     "used|scor", "cells|scor", "spk/bin", "ripple", "shanks")
            for low_sd in GRID:
                tuned = copy.deepcopy(params)
                tuned["events"]["mua"]["threshold_low_sd"] = float(low_sd)
                events, intervals = detect_mua_events(session, epoch, tuned)
                if len(events.start) == 0:
                    # A row rather than a skip: below some boundary the trace
                    # never falls back, every candidate merges, and the merged
                    # events run past the duration cap and are dropped.
                    rows.append({"session": name, "epoch": epoch,
                                 "threshold_low_sd": float(low_sd), "n_events": 0})
                    log.info("  %s  %7.2f %9d   (no events: all merged past the "
                             "duration cap)", " " * len(epoch), low_sd, 0)
                    continue
                place_struct = per_event_structure(session, events, place)
                pyr_struct = per_event_structure(session, events, pyr)
                matched, support = ripple_match(events, ripples)
                scorable = place_struct["used_bins"] >= MIN_BINS
                duration_ms = 1000 * (events.end - events.start)
                state_min = events.state_seconds / 60
                row = {
                    "session": name, "epoch": epoch, "threshold_low_sd": float(low_sd),
                    "threshold_low_rate": float(events.params["threshold_low"]),
                    "baseline_mean": events.baseline_mean,
                    "baseline_sd": events.baseline_sd,
                    "n_events": len(events.start),
                    "events_per_min": len(events.start) / state_min if state_min else np.nan,
                    "median_duration_ms": float(np.median(duration_ms)),
                    "pct_at_duration_cap": float(100 * np.mean(
                        duration_ms >= 0.95 * params["events"]["mua"]["max_duration_ms"])),
                    "median_bins": float(np.median(place_struct["total_bins"])),
                    "median_used_bins_place": float(np.median(place_struct["used_bins"])),
                    "median_used_bins_pyr": float(np.median(pyr_struct["used_bins"])),
                    "pct_scorable": float(100 * np.mean(scorable)),
                    "n_scorable": int(scorable.sum()),
                    "median_active_place": float(np.median(place_struct["active"])),
                    "median_active_pyr": float(np.median(pyr_struct["active"])),
                    # Among scorable events only: the apples-to-apples figure
                    # against Chen et al.'s 11.5 bins, which is a median over
                    # the events they actually scored.
                    "median_used_bins_scorable": float(
                        np.median(pyr_struct["used_bins"][scorable]))
                    if scorable.any() else np.nan,
                    "median_active_pyr_scorable": float(
                        np.median(pyr_struct["active"][scorable]))
                    if scorable.any() else np.nan,
                    "median_spikes_per_used_bin": float(np.median(
                        pyr_struct["spikes"] / np.maximum(pyr_struct["used_bins"], 1))),
                    "median_empty_bin_frac": float(np.median(
                        1 - pyr_struct["used_bins"] / np.maximum(pyr_struct["total_bins"], 1))),
                    "pct_ripple_matched": float(100 * np.mean(matched)),
                    "n_comparison_shanks": comparison_shanks,
                    "median_support_matched": float(np.nanmedian(support))
                    if np.isfinite(support).any() else np.nan,
                    "median_peak_z": float(np.median(events.peak_z)),
                }
                if low_sd == params["events"]["mua"]["threshold_low_sd"]:
                    # The sweep at the pipeline's own threshold must reproduce
                    # the pipeline's own events. If it does not, the sweep is
                    # measuring something other than this detector.
                    here = committed[(committed.session == name)
                                     & (committed.epoch == epoch)]
                    if len(here):
                        gap = abs(float(np.median(here.duration_ms))
                                  - row["median_duration_ms"])
                        if len(here) != row["n_events"] or gap > 1.0:
                            raise SystemExit(
                                f"{name} {epoch}: at the pipeline's own threshold the "
                                f"sweep finds {row['n_events']} events of "
                                f"{row['median_duration_ms']:.1f} ms against the "
                                f"committed {len(here)} of "
                                f"{np.median(here.duration_ms):.1f} ms")
                        log.info("  %s  %7.2f   reproduces results/mua_events.csv "
                                 "exactly (%d events)", " " * len(epoch), low_sd,
                                 len(here))
                rows.append(row)
                log.info("  %s  %7.2f %9d %9.0f %6.1f %6.0f%% %10.1f %11.1f %8.2f "
                         "%6.0f%% %8.1f",
                         " " * len(epoch), low_sd, row["n_events"],
                         row["median_duration_ms"], row["median_used_bins_pyr"],
                         row["pct_scorable"], row["median_used_bins_scorable"],
                         row["median_active_pyr_scorable"],
                         row["median_spikes_per_used_bin"], row["pct_ripple_matched"],
                         row["median_support_matched"])

    if not rows:
        return
    table = write_session_table(results / "boundary_threshold_sweep.csv", rows,
                                sessions=args.sessions, round_to=4)
    stamp = {"script": "scripts/boundary_threshold_sweep.py",
             "git_commit": config.git_commit(),
             "created": datetime.now(timezone.utc).isoformat(timespec="seconds"),
             "grid": list(GRID), "bin_s": BIN_S, "min_bins": MIN_BINS,
             "note": "sweep only; no threshold adopted (D10.1)",
             "events": params["events"]["mua"]}
    (results / "boundary_threshold_sweep.json").write_text(
        json.dumps(stamp, indent=2) + "\n")
    if not args.no_figures:
        plots.plot_boundary_sweep(table.to_dict("records"),
                                  results / "boundary_threshold_sweep.png", stamp)
    log.info("\nwrote %s", (results / "boundary_threshold_sweep.csv").relative_to(
        config.repo_root()))


if __name__ == "__main__":
    main()
