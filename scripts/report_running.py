"""Report how the running mask is built for every session.

Stage 1 is the authors' mask (OneDLocation defined); stage 2 applies our speed
threshold on top. Writes results/running_summary.csv.

    uv run python scripts/report_running.py
"""

from __future__ import annotations

import json
from datetime import datetime, timezone

import numpy as np
import pandas as pd

from hc11 import config
from hc11.behavior import running_mask_from_config
from hc11.io import list_sessions, load_session
from hc11.preprocess import truncate_to_session
from hc11.track import Track


def main() -> None:
    params = config.load_params()
    rows = []
    for name in list_sessions():
        session = load_session(name)
        if params["preprocess"]["truncate_to_session"]:
            session, _ = truncate_to_session(session, verbose=False)
        rm = running_mask_from_config(session, params)
        track = Track.from_session(session, params["place_fields"]["bin_size_cm"])
        speed_authors = rm.speed[rm.authors_mask]
        rows.append({
            "session": name,
            "maze": session.maze_type,
            "maze_min": session.epoch_duration("MAZE") / 60,
            "authors_min": rm.authors_seconds / 60,
            "final_min": rm.seconds / 60,
            "pct_of_maze": 100 * rm.seconds / session.epoch_duration("MAZE"),
            "pct_removed_by_speed": 100 * rm.fraction_removed_by_speed,
            "pct_below_5cm_s": 100 * float(np.mean(np.nan_to_num(speed_authors, nan=1e9) < 5)),
            "n_restored_bridging": rm.n_restored_by_bridging,
            "n_removed_short": rm.n_removed_short_epochs,
            "n_periods": len(rm.intervals),
            "median_speed_authors": rm.median_speed("authors"),
            "median_speed_final": rm.median_speed("final"),
            "track_extent_m": track.extent,
            "n_bins": track.n_bins,
        })

    df = pd.DataFrame(rows)
    header = (f"{'session':<18} {'maze':>6} {'authors':>8} {'final':>7} {'%maze':>6} "
              f"{'%cut by v':>10} {'%v<5':>6} {'periods':>8} {'med v':>6}")
    print(header)
    print("-" * len(header))
    for r in rows:
        print(f"{r['session']:<18} {r['maze_min']:6.1f} {r['authors_min']:8.1f} "
              f"{r['final_min']:7.1f} {r['pct_of_maze']:5.1f}% {r['pct_removed_by_speed']:9.1f}% "
              f"{r['pct_below_5cm_s']:5.1f}% {r['n_periods']:8d} {r['median_speed_final']:6.1f}")
    print("-" * len(header))
    print("(minutes; '%cut by v' = share of the authors' mask removed by the speed threshold)")

    out = config.results_dir() / "running_summary.csv"
    df.round(4).to_csv(out, index=False)
    out.with_suffix(".json").write_text(json.dumps({
        "script": "scripts/report_running.py",
        "git_commit": config.git_commit(),
        "created": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "params": {"behavior": params["behavior"], "preprocess": params["preprocess"]},
    }, indent=2) + "\n")
    print(f"\nwrote {out.relative_to(config.repo_root())}")


if __name__ == "__main__":
    main()
