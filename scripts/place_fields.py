"""Compute one-dimensional place fields for one or all sessions.

    uv run python scripts/place_fields.py                     # all sessions
    uv run python scripts/place_fields.py Achilles_10252013    # one session
    uv run python scripts/place_fields.py --no-figures

Writes results/place_fields_summary.csv (one row per cell per direction),
results/lap_summary.csv, and per-session figures in results/<session>/.
"""

from __future__ import annotations

import argparse
import json
import logging
from datetime import datetime, timezone

import numpy as np
import pandas as pd

from hc11 import config, plots
from hc11.behavior import running_mask_from_config
from hc11.fields import analyse_cell, select_running_samples
from hc11.io import Session, list_sessions, load_session
from hc11.laps import DIRECTIONS, lap_summary, segment_laps_from_config
from hc11.preprocess import truncate_to_session
from hc11.track import Track

log = logging.getLogger("place_fields")

#: Below this many laps in a direction, treat that direction as under-sampled.
LAP_WARNING_THRESHOLD = 20


def prepare(session_name: str, params: dict):
    """Load a session and derive track, running mask and laps."""
    session = load_session(session_name)
    if params["preprocess"]["truncate_to_session"]:
        session, _ = truncate_to_session(session, verbose=False)
    track = Track.from_session(session, params["place_fields"]["bin_size_cm"])
    running = running_mask_from_config(session, params)
    laps = segment_laps_from_config(session, track, params)
    return session, track, running, laps


def analyse_session(
    session: Session, track: Track, running, laps, params: dict, rng: np.random.Generator
) -> list:
    cells = []
    for direction in DIRECTIONS:
        samples = select_running_samples(session, track, laps, running.mask, direction)
        if samples.n == 0:
            log.warning("%s %s: no running samples in this direction", session.name, direction)
            continue
        for cluster_id in session.pyr_ids:
            cells.append(analyse_cell(session, track, samples, cluster_id, params, rng))
    return cells


def _row(cell) -> dict:
    return {k: v for k, v in cell.__dict__.items() if k != "rate_map"}


def main() -> None:
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("sessions", nargs="*", help="session names (default: all)")
    parser.add_argument("--no-figures", action="store_true")
    args = parser.parse_args()

    params = config.load_params()
    names = args.sessions or list_sessions()
    seed = params["place_fields"]["shuffle"]["random_seed"]

    all_rows, lap_rows = [], []
    for name in names:
        session, track, running, laps = prepare(name, params)
        summary = lap_summary(laps)
        rng = np.random.default_rng(seed)
        cells = analyse_session(session, track, running, laps, params, rng)

        row = {
            "session": name, "animal": session.animal, "maze_type": session.maze_type,
            "maze_min": session.epoch_duration("MAZE") / 60,
            "running_min": running.seconds / 60,
            "track_extent_m": track.extent, "n_bins": track.n_bins,
        }
        for d in DIRECTIONS:
            row[f"n_laps_{d}"] = summary[d]["n_laps"]
            row[f"lap_s_{d}"] = summary[d]["total_s"]
            row[f"lap_median_dur_{d}"] = summary[d]["median_duration_s"]
            row[f"lap_median_cov_{d}"] = summary[d]["median_coverage"]
            sel = [c for c in cells if c.direction == d]
            row[f"n_place_cells_{d}"] = sum(c.is_place_cell for c in sel)
        unique = {c.cluster_id for c in cells if c.is_place_cell}
        pairs = sum(c.is_place_cell for c in cells)
        unique_nostab = {c.cluster_id for c in cells if c.is_place_cell_no_stability}
        row["n_pyr"] = len(session.pyr_ids)
        row["n_place_cells_unique"] = len(unique)
        row["n_place_cell_directions"] = int(pairs)
        row["n_place_cells_unique_no_stability"] = len(unique_nostab)
        row["under_sampled"] = any(
            0 < summary[d]["n_laps"] < LAP_WARNING_THRESHOLD for d in DIRECTIONS
        )
        lap_rows.append(row)
        all_rows.extend(_row(c) for c in cells)

        log.info(
            "%s: laps %d/%d (pos/neg), %d/%d pyramidal cells are place cells",
            name, summary["pos"]["n_laps"], summary["neg"]["n_laps"],
            len(unique), len(session.pyr_ids),
        )
        for d in DIRECTIONS:
            n = summary[d]["n_laps"]
            if 0 < n < LAP_WARNING_THRESHOLD:
                log.warning("%s %s: only %d laps -- under-sampled", name, d, n)

        if not args.no_figures:
            out = config.results_dir(name)
            meta = {
                "session": name, "git_commit": config.git_commit(),
                "created": datetime.now(timezone.utc).isoformat(timespec="seconds"),
                "params": {
                    k: params[k]
                    for k in ("behavior", "laps", "place_fields", "preprocess")
                },
                "track": str(track), "laps": summary,
            }
            plots.plot_laps(session, track, laps, running, out / "laps.png", meta)
            plots.plot_info_histogram(cells, name, out / "spatial_info_hist.png", meta)
            for d in DIRECTIONS:
                sel = [c for c in cells if c.direction == d]
                if not sel:
                    continue
                plots.plot_rate_map_heatmap(
                    sel, track, name, d, out / f"rate_maps_{d}.png", meta)
                plots.plot_field_grid(
                    sel, track, name, d, out / f"place_fields_{d}.png", meta)

    cells_df = pd.DataFrame(all_rows)
    laps_df = pd.DataFrame(lap_rows)
    results = config.results_dir()
    cells_df.round(6).to_csv(results / "place_fields_summary.csv", index=False)
    laps_df.round(4).to_csv(results / "lap_summary.csv", index=False)
    (results / "place_fields_summary.json").write_text(json.dumps({
        "script": "scripts/place_fields.py",
        "git_commit": config.git_commit(),
        "created": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "sessions": names,
        "params": params,
    }, indent=2) + "\n")

    print()
    header = (f"{'session':<18} {'maze':<17} {'laps +':>7} {'laps -':>7} {'pyr':>5} "
              f"{'PC':>5} {'PC%':>6} {'cell-dir':>9} {'run min':>8}")
    print(header)
    print("-" * len(header))
    for r in lap_rows:
        flag = " *" if r["under_sampled"] else ""
        print(f"{r['session']:<18} {r['maze_type']:<17} {r['n_laps_pos']:7d} {r['n_laps_neg']:7d} "
              f"{r['n_pyr']:5d} {r['n_place_cells_unique']:5d} "
              f"{100 * r['n_place_cells_unique'] / r['n_pyr']:5.0f}% "
              f"{r['n_place_cell_directions']:9d} {r['running_min']:8.1f}{flag}")
    print("-" * len(header))
    total_pyr = laps_df.n_pyr.sum()
    total_pc = laps_df.n_place_cells_unique.sum()
    print(f"{'total':<18} {'':<17} {laps_df.n_laps_pos.sum():7d} {laps_df.n_laps_neg.sum():7d} "
          f"{total_pyr:5d} {total_pc:5d} {100 * total_pc / total_pyr:5.0f}% "
          f"{laps_df.n_place_cell_directions.sum():9d}")
    print(f"\n* = fewer than {LAP_WARNING_THRESHOLD} laps in a direction")
    if total_pc and set(names) == set(config.SESSIONS):
        # Only meaningful for the complete set, and only against unique cells:
        # the paper's 491 counts neurons, not cell-direction pairs (D5.1).
        print(f"paper reports 491 place cells (unique neurons) across all 8 sessions; "
              f"ours is {total_pc} unique ({491 / total_pc:.2f}x). "
              f"See scripts/place_cell_sensitivity.py")


if __name__ == "__main__":
    main()
