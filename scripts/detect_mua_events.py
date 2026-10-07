"""Detect and characterise population-synchrony (MUA) events from spikes alone.

    uv run python scripts/detect_mua_events.py
    uv run python scripts/detect_mua_events.py Achilles_10252013 --no-figures

Writes results/mua_events.csv (one row per event), results/mua_summary.csv (one
row per session and epoch) and results/mua_threshold_sensitivity.csv, plus
per-session figures.
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
from hc11.behavior import running_mask_from_config
from hc11.events import participation
from hc11.fields import analyse_cell, select_running_samples
from hc11.io import load_session
from hc11.laps import DIRECTIONS, segment_laps_from_config
from hc11.mua import detect_mua_events
from hc11.preprocess import truncate_to_session
from hc11.track import Track

log = logging.getLogger("mua")

SESSIONS = ["Achilles_10252013", "Achilles_11012013"]
EPOCHS = ("PRE", "MAZE", "POST")
THRESHOLD_SWEEP = (2.5, 3.0, 4.0, 5.0)


def place_cell_ids(session, params) -> np.ndarray:
    """The decoding ensemble: place cells in either direction (step 5a)."""
    track = Track.from_session(session, params["place_fields"]["bin_size_cm"])
    running = running_mask_from_config(session, params)
    laps = segment_laps_from_config(session, track, params)
    rng = np.random.default_rng(params["place_fields"]["shuffle"]["random_seed"])
    found = set()
    for direction in DIRECTIONS:
        samples = select_running_samples(session, track, laps, running.mask, direction)
        if samples.n == 0:
            continue
        found |= {
            c.cluster_id
            for c in (analyse_cell(session, track, samples, cid, params, rng)
                      for cid in session.pyr_ids)
            if c.is_place_cell
        }
    return np.array(sorted(found))


def describe(events, part_pyr, part_pc, session, epoch, params, label="primary") -> dict:
    cfg = params["events"]["decodable"]
    dur_ms = events.duration * 1000 if len(events) else np.array([np.nan])
    n_active = part_pyr.n_active if len(events) else np.array([0])
    n_active_pc = part_pc.n_active if len(events) else np.array([0])
    return {
        "session": session.name, "epoch": epoch, "variant": label,
        "n_events": len(events),
        "state_min": events.state_seconds / 60,
        "rate_per_min": events.rate_per_min,
        "rate_per_s": len(events) / events.state_seconds if events.state_seconds else np.nan,
        "baseline_hz": events.baseline_mean, "baseline_sd_hz": events.baseline_sd,
        "dur_p10_ms": float(np.nanpercentile(dur_ms, 10)),
        "dur_median_ms": float(np.nanmedian(dur_ms)),
        "dur_p90_ms": float(np.nanpercentile(dur_ms, 90)),
        "n_pyr": part_pyr.n_cells, "n_place_cells": part_pc.n_cells,
        "active_pyr_median": float(np.median(n_active)),
        "active_pyr_p90": float(np.percentile(n_active, 90)),
        "active_frac_median": float(np.median(n_active) / max(part_pyr.n_cells, 1)),
        "spikes_median": float(np.median(part_pyr.n_spikes)) if len(events) else np.nan,
        "spikes_p90": float(np.percentile(part_pyr.n_spikes, 90)) if len(events) else np.nan,
        "active_pc_median": float(np.median(n_active_pc)),
        "pct_ge5_pyr": 100 * float(np.mean(n_active >= cfg["min_active_cells"])),
        "pct_ge10pct_pyr": 100 * float(
            np.mean(n_active >= cfg["min_active_fraction"] * part_pyr.n_cells)),
        "pct_ge5_pc": 100 * float(np.mean(n_active_pc >= cfg["min_active_cells"])),
        "pct_ge10pct_pc": 100 * float(
            np.mean(n_active_pc >= cfg["min_active_fraction"] * max(part_pc.n_cells, 1))),
        "n_ge5_pc": int(np.sum(n_active_pc >= cfg["min_active_cells"])),
    }


def main() -> None:
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("sessions", nargs="*")
    parser.add_argument("--no-figures", action="store_true")
    args = parser.parse_args()

    params = config.load_params()
    names = args.sessions or SESSIONS
    summary, per_event, sweep = [], [], []

    for name in names:
        session = load_session(name)
        if params["preprocess"]["truncate_to_session"]:
            session, _ = truncate_to_session(session, verbose=False)
        pyr = session.pyr_ids
        place = place_cell_ids(session, params)
        log.info("%s: %d pyramidal cells, %d place cells", name, len(pyr), len(place))

        figure_data = {}
        for epoch in EPOCHS:
            events, intervals = detect_mua_events(session, epoch, params)
            part_pyr = participation(session, events, pyr)
            part_pc = participation(session, events, place)
            summary.append(describe(events, part_pyr, part_pc, session, epoch, params))
            for i in range(len(events)):
                per_event.append({
                    "session": name, "epoch": epoch,
                    "start": events.start[i], "peak": events.peak[i], "end": events.end[i],
                    "duration_ms": 1000 * (events.end[i] - events.start[i]),
                    "peak_z": events.peak_z[i],
                    "n_active_pyr": part_pyr.n_active[i], "n_spikes": part_pyr.n_spikes[i],
                    "n_active_place_cells": part_pc.n_active[i],
                })
            figure_data[epoch] = (events, part_pyr, intervals)
            log.info("  %-5s %5d events, %5.2f/min over %6.1f min, median %3.0f ms, "
                     "median %2.0f/%d active", epoch, len(events), events.rate_per_min,
                     events.state_seconds / 60, np.median(events.duration * 1000) if len(events)
                     else np.nan, np.median(part_pyr.n_active) if len(events) else 0, len(pyr))

            # Threshold sensitivity, same pipeline, high threshold varied.
            for sd in THRESHOLD_SWEEP:
                p2 = copy.deepcopy(params)
                p2["events"]["mua"]["threshold_high_sd"] = sd
                ev2, _ = detect_mua_events(session, epoch, p2)
                pp2 = participation(session, ev2, pyr)
                sweep.append({
                    "session": name, "epoch": epoch, "threshold_high_sd": sd,
                    "n_events": len(ev2), "rate_per_s": (
                        len(ev2) / ev2.state_seconds if ev2.state_seconds else np.nan),
                    "dur_median_ms": float(np.median(ev2.duration * 1000)) if len(ev2) else np.nan,
                    "active_median": float(np.median(pp2.n_active)) if len(ev2) else np.nan,
                })

        if not args.no_figures:
            out = config.results_dir(name)
            meta = {
                "session": name, "git_commit": config.git_commit(),
                "created": datetime.now(timezone.utc).isoformat(timespec="seconds"),
                "params": params["events"], "n_pyr": len(pyr), "n_place_cells": len(place),
            }
            best = max(EPOCHS, key=lambda e: len(figure_data[e][0]))
            events, part, _ = figure_data[best]
            plots.plot_event_examples(session, events, part, pyr, best,
                                      out / "mua_examples.png", meta)
            plots.plot_event_rate_over_session(
                session, {e: figure_data[e][0] for e in EPOCHS},
                out / "mua_event_rate.png", meta)
            plots.plot_participation(
                {e: (figure_data[e][0], figure_data[e][1]) for e in EPOCHS}, len(pyr),
                name, out / "mua_participation.png", meta)

    results = config.results_dir()
    pd.DataFrame(summary).round(4).to_csv(results / "mua_summary.csv", index=False)
    pd.DataFrame(per_event).round(6).to_csv(results / "mua_events.csv", index=False)
    pd.DataFrame(sweep).round(4).to_csv(results / "mua_threshold_sensitivity.csv", index=False)
    (results / "mua_summary.json").write_text(json.dumps({
        "script": "scripts/detect_mua_events.py", "git_commit": config.git_commit(),
        "created": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "sessions": names, "params": params["events"],
    }, indent=2) + "\n")

    df = pd.DataFrame(summary)
    print()
    header = (f"{'session':<18} {'epoch':<5} {'state min':>9} {'events':>7} {'/min':>6} {'/s':>6} "
              f"{'dur ms':>7} {'act pyr':>8} {'>=5pc':>7} {'>=10%pc':>8}")
    print(header)
    print("-" * len(header))
    for _, r in df.iterrows():
        print(f"{r.session:<18} {r.epoch:<5} {r.state_min:9.1f} {r.n_events:7d} "
              f"{r.rate_per_min:6.1f} {r.rate_per_s:6.2f} {r.dur_median_ms:7.0f} "
              f"{r.active_pyr_median:8.0f} {r.pct_ge5_pc:6.0f}% {r.pct_ge10pct_pc:7.0f}%")
    print("-" * len(header))
    print("(act pyr = median active pyramidal cells per event; >=5pc / >=10%pc = share of "
          "events\n with at least 5 / 10% of place cells active -- the step 7 decodability "
          "numbers)")
    print(f"\nwrote {(results / 'mua_summary.csv').relative_to(config.repo_root())}")


if __name__ == "__main__":
    main()
