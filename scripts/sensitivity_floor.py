"""The positive-control series: sequence score against ensemble size.

    uv run python scripts/sensitivity_floor.py

Each point is one session's median rZ for 2 s windows of real MAZE running —
the validation control written by scripts/validate_sequences.py, read back
here from results/sequence_validation.csv and plotted against the number of
place cells in that session's decoding ensemble (results/lap_summary.csv).

These windows are the one case where the ground truth is known: the animal ran
a trajectory, the decoder recovers it, so a score near zero is a failure of the
scoring, not an absence of sequence. Reading the five sessions in order of
ensemble size therefore says where the method's sensitivity floor lies, and the
replay effect sizes in results/replay_comparison.csv should be read against it
(D8.9).

Writes results/sensitivity_floor.csv and results/sensitivity_floor.png.
"""

from __future__ import annotations

from datetime import datetime, timezone

import pandas as pd

from hc11 import config, plots

WINDOW_LABEL = "MAZE running 2 s"


def main() -> None:
    params = config.load_params()
    results = config.results_dir()
    validation = pd.read_csv(results / "sequence_validation.csv")
    laps = pd.read_csv(results / "lap_summary.csv").set_index("session")

    windows = validation[(validation["control"] == "running")
                         & (validation["set"] == WINDOW_LABEL)]
    if windows.empty:
        raise SystemExit(
            f"no '{WINDOW_LABEL}' rows in results/sequence_validation.csv — "
            "run scripts/validate_sequences.py first")

    rows = [{
        "session": r.session,
        "n_place_cells": int(laps.loc[r.session].n_place_cells_unique),
        "geometry": ("circular" if "circ" in str(laps.loc[r.session].maze_type).lower()
                     else "linear"),
        "shuffle": r.shuffle,
        "n_windows": int(r.n),
        "median_rz": float(r.median_rz),
        "mean_rz": float(r.mean_rz),
        "pct_above_z": float(r.pct_above_z),
    } for r in windows.itertuples()]
    rows.sort(key=lambda r: (-r["n_place_cells"], r["shuffle"]))

    header = (f"{'session':<20}{'cells':>6}{'shuffle':>19}{'windows':>9}"
              f"{'median rZ':>11}{'% above z':>11}")
    print(header)
    print("-" * len(header))
    for r in rows:
        print(f"{r['session']:<20}{r['n_place_cells']:6d}{r['shuffle']:>19}"
              f"{r['n_windows']:9d}{r['median_rz']:+11.2f}{r['pct_above_z']:10.1f}%")
    print("-" * len(header))
    print("(2 s windows of real running: the decoder's own positive control)")

    out_csv = results / "sensitivity_floor.csv"
    pd.DataFrame(rows).round(4).to_csv(out_csv, index=False)

    figure_params = {
        "script": "scripts/sensitivity_floor.py",
        "git_commit": config.git_commit(),
        "created": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "window": WINDOW_LABEL,
        "sessions": list(config.LFP_SESSIONS),
        "omitted": config.LFP_OMITTED,
        "sequences": params["sequences"],
    }
    # plots._save writes the sidecar beside the figure, which is the same
    # results/sensitivity_floor.json the CSV would claim; one write, not two.
    out_png = plots.plot_sensitivity_floor(
        rows, results / "sensitivity_floor.png", figure_params)

    print(f"\nwrote {out_csv.relative_to(config.repo_root())}")
    print(f"wrote {out_png.relative_to(config.repo_root())}")


if __name__ == "__main__":
    main()
