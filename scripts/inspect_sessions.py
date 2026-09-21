"""Load every hc-11 session, print a one-line summary each, and save a CSV.

Usage::

    uv run python scripts/inspect_sessions.py            # all sessions
    uv run python scripts/inspect_sessions.py --warnings  # also list validation warnings

Writes ``results/session_overview.csv`` and a sidecar
``results/session_overview.json`` recording the git commit.
"""

from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone

import numpy as np
import pandas as pd

from hc11 import config
from hc11.io import Session, list_sessions, load_session
from hc11.validate import validate


def summarize(s: Session) -> dict:
    issues = validate(s)
    nan_2d = np.isnan(s.position_xy).any(axis=1)
    nan_1d = np.isnan(s.position_1d)
    return {
        "session": s.name,
        "animal": s.animal,
        "maze_type": s.maze_type,
        "pre_min": s.epoch_duration("PRE") / 60,
        "maze_min": s.epoch_duration("MAZE") / 60,
        "post_min": s.epoch_duration("POST") / 60,
        "n_pyr": len(s.pyr_ids),
        "n_int": len(s.int_ids),
        "n_shanks": len({s.shank_of_cluster[c] for c in s.cluster_ids}),
        "n_spikes": len(s.spike_times),
        "n_pos_samples": len(s.position_t),
        "pos_rate_hz": 1.0 / s.position_dt,
        "pct_nan_2d": 100 * nan_2d.mean(),
        "pct_nan_1d": 100 * nan_1d.mean(),
        "min_1d_valid": (~nan_1d).sum() * s.position_dt / 60,
        "n_errors": sum(i.severity == "error" for i in issues),
        "n_warnings": sum(i.severity == "warning" for i in issues),
    }, issues


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--warnings", action="store_true", help="print validation warnings")
    args = parser.parse_args()

    rows, all_issues = [], {}
    header = (f"{'session':<18} {'animal':<9} {'maze':<17} {'PRE':>6} {'MAZE':>6} {'POST':>6} "
              f"{'pyr':>4} {'int':>4} {'spikes':>11} {'NaN2D':>6} {'NaN1D':>6}")
    print(header)
    print("-" * len(header))
    for name in list_sessions():
        row, issues = summarize(load_session(name))
        rows.append(row)
        all_issues[name] = issues
        print(f"{row['session']:<18} {row['animal']:<9} {row['maze_type']:<17} "
              f"{row['pre_min']:6.1f} {row['maze_min']:6.1f} {row['post_min']:6.1f} "
              f"{row['n_pyr']:4d} {row['n_int']:4d} {row['n_spikes']:11,d} "
              f"{row['pct_nan_2d']:5.1f}% {row['pct_nan_1d']:5.1f}%")
    df = pd.DataFrame(rows)
    print("-" * len(header))
    print(f"{'total':<46} {'':>6} {'':>6} {'':>6} {df.n_pyr.sum():4d} {df.n_int.sum():4d} "
          f"{df.n_spikes.sum():11,d}")
    print("(epoch durations in minutes)")

    if args.warnings:
        print()
        for name, issues in all_issues.items():
            flagged = [i for i in issues if i.severity != "info"]
            if flagged:
                print(name)
                for issue in flagged:
                    print(f"  {issue}")

    out = config.results_dir() / "session_overview.csv"
    df.round(4).to_csv(out, index=False)
    sidecar = {
        "script": "scripts/inspect_sessions.py",
        "git_commit": config.git_commit(),
        "created": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "data_dir": str(config.data_dir()),
        "columns": {
            "pre_min/maze_min/post_min": "epoch durations, minutes",
            "pct_nan_2d": "% of position samples with TwoDLocation NaN",
            "pct_nan_1d": "% of position samples with OneDLocation NaN",
            "min_1d_valid": "minutes of MAZE with a defined OneDLocation",
        },
    }
    out.with_suffix(".json").write_text(json.dumps(sidecar, indent=2) + "\n")
    print(f"\nwrote {out.relative_to(config.repo_root())}")


if __name__ == "__main__":
    main()
