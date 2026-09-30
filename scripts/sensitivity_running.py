"""Does the speed threshold change the place fields, or only the sample count?

Builds rate maps three ways for one session from each group identified in
decision D3.4 -- Achilles_10252013, where the authors' mask barely overlaps the
speed criterion (5.5 % removed), and Cicero_09012014, where it removes 51.9 % --
and compares what comes out:

    mask only          the authors' OneDLocation mask, no speed criterion
    mask + 15 cm/s     the configured analysis (decision D3.4)
    mask + 5 cm/s      a permissive threshold, excluding only near-immobility

Laps are segmented once per session from the authors' mask, so all three
conditions see identical laps and differ only in which samples are counted.

    uv run python scripts/sensitivity_running.py

Writes results/sensitivity_running.csv.
"""

from __future__ import annotations

import copy
import json
from datetime import datetime, timezone

import numpy as np
import pandas as pd

from hc11 import config
from hc11.behavior import running_mask_from_config
from hc11.fields import analyse_cell, select_running_samples
from hc11.io import load_session
from hc11.laps import DIRECTIONS, segment_laps_from_config
from hc11.preprocess import truncate_to_session
from hc11.track import Track

SESSIONS = ["Achilles_10252013", "Cicero_09012014"]
CONDITIONS = {"mask only": 0.0, "mask + 15 cm/s": 15.0, "mask + 5 cm/s": 5.0}


def main() -> None:
    params = config.load_params()
    rows = []
    for name in SESSIONS:
        session = load_session(name)
        if params["preprocess"]["truncate_to_session"]:
            session, _ = truncate_to_session(session, verbose=False)
        track = Track.from_session(session, params["place_fields"]["bin_size_cm"])
        laps = segment_laps_from_config(session, track, params)  # identical across conditions

        for label, threshold in CONDITIONS.items():
            p = copy.deepcopy(params)
            p["behavior"]["speed_threshold_cm_s"] = threshold
            running = running_mask_from_config(session, p)
            rng = np.random.default_rng(p["place_fields"]["shuffle"]["random_seed"])

            cells = []
            for direction in DIRECTIONS:
                samples = select_running_samples(session, track, laps, running.mask, direction)
                if samples.n == 0:
                    continue
                cells.extend(
                    analyse_cell(session, track, samples, c, p, rng) for c in session.pyr_ids
                )
            if not cells:
                continue
            pc = [c for c in cells if c.is_place_cell]
            info = np.array([c.info_bits_per_spike for c in pc], dtype=float)
            stab = np.array([c.stability for c in pc], dtype=float)
            all_stab = np.array([c.stability for c in cells], dtype=float)
            rows.append({
                "session": name,
                "condition": label,
                "speed_threshold_cm_s": threshold,
                "running_min": running.seconds / 60,
                "pct_of_authors_mask": 100 * running.fraction_kept,
                "n_place_cells": len({c.cluster_id for c in pc}),
                "n_place_cell_directions": len(pc),
                "median_info_place_cells": float(np.nanmedian(info)) if info.size else np.nan,
                "median_stability_place_cells": float(np.nanmedian(stab)) if stab.size else np.nan,
                "median_stability_all": float(np.nanmedian(all_stab)),
                "median_peak_rate_place_cells": float(
                    np.nanmedian([c.peak_rate for c in pc])) if pc else np.nan,
                "median_field_width_m": float(
                    np.nanmedian([c.largest_field_width_m for c in pc])) if pc else np.nan,
            })

    df = pd.DataFrame(rows)
    out = config.results_dir() / "sensitivity_running.csv"
    df.round(4).to_csv(out, index=False)
    out.with_suffix(".json").write_text(json.dumps({
        "script": "scripts/sensitivity_running.py",
        "git_commit": config.git_commit(),
        "created": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "conditions": CONDITIONS,
    }, indent=2) + "\n")

    header = (f"{'session':<18} {'condition':<15} {'run min':>8} {'% mask':>7} {'PC':>4} "
              f"{'med info':>9} {'med stab':>9} {'med peak':>9} {'width m':>8}")
    print(header)
    print("-" * len(header))
    for r in rows:
        print(f"{r['session']:<18} {r['condition']:<15} {r['running_min']:8.1f} "
              f"{r['pct_of_authors_mask']:6.0f}% {r['n_place_cells']:4d} "
              f"{r['median_info_place_cells']:9.2f} {r['median_stability_place_cells']:9.2f} "
              f"{r['median_peak_rate_place_cells']:9.1f} {r['median_field_width_m']:8.2f}")
    print(f"\nwrote {out.relative_to(config.repo_root())}")


if __name__ == "__main__":
    main()
