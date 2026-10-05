"""How many place cells, under which criteria and which counting convention.

Grosmark & Buzsaki report "n = 491 cells" across the eight sessions. That number
counts **unique neurons**, not cell-direction pairs: the same 491 is later split
by each cell's sleep firing rate (Fig. 4D) and 216 of them are singled out as
strongly contributing neurons (Fig. 2B, one point per neuron). So the like-for-
like comparison with our pipeline is the unique-cell column.

Variants (the primary one is whatever `config/params.yaml` holds):

    primary                 as configured
    cap relaxed             no maximum field width
    cap relaxed + 20 spikes also lowers the running-spike floor from 50 to 20
    peak rate only          the most permissive reading of the paper's own
                            criterion -- Chen et al. 2016's methods say "all
                            putative pyramidal neurons selected for analysis had
                            peak firing rate > 1 Hz" -- with no field, shuffle or
                            stability requirement at all

    uv run python scripts/place_cell_sensitivity.py

Writes results/place_cell_sensitivity.csv.
"""

from __future__ import annotations

import copy
import json
from datetime import datetime, timezone

import numpy as np
import pandas as pd

from hc11 import config
from place_fields import analyse_session, prepare  # noqa: I001  (sibling script)

PAPER_TOTAL = 491

VARIANTS: dict[str, dict] = {
    "primary": {},
    "cap relaxed": {"max_field_width_frac": 1.01},
    "cap relaxed + 20 spikes": {"max_field_width_frac": 1.01, "min_spikes": 20},
}


def _counts(cells, predicate) -> tuple[int, int]:
    """(unique cells, cell-direction pairs) satisfying ``predicate``."""
    selected = [c for c in cells if predicate(c)]
    return len({c.cluster_id for c in selected}), len(selected)


def main() -> None:
    base = config.load_params()
    rows = []

    for label, overrides in VARIANTS.items():
        params = copy.deepcopy(base)
        params["place_fields"]["criteria"].update(overrides)
        for name in config.SESSIONS:
            session, track, running, laps = prepare(name, params)
            rng = np.random.default_rng(params["place_fields"]["shuffle"]["random_seed"])
            cells = analyse_session(session, track, running, laps, params, rng)
            unique, pairs = _counts(cells, lambda c: c.is_place_cell)
            min_peak = params["place_fields"]["criteria"]["min_peak_rate_hz"]
            peak_u, peak_p = _counts(
                cells,
                lambda c, threshold=min_peak: np.isfinite(c.peak_rate)
                and c.peak_rate >= threshold,
            )
            rows.append({
                "variant": label, "session": name, "maze_type": session.maze_type,
                "n_pyr": len(session.pyr_ids),
                "unique_place_cells": unique, "cell_direction_pairs": pairs,
                "peak_rate_only_unique": peak_u, "peak_rate_only_pairs": peak_p,
            })

    df = pd.DataFrame(rows)
    out = config.results_dir() / "place_cell_sensitivity.csv"
    df.to_csv(out, index=False)
    out.with_suffix(".json").write_text(json.dumps({
        "script": "scripts/place_cell_sensitivity.py",
        "git_commit": config.git_commit(),
        "created": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "paper_total_unique_cells": PAPER_TOTAL,
        "variants": VARIANTS,
    }, indent=2) + "\n")

    n_pyr = int(df[df.variant == "primary"].n_pyr.sum())
    print(f"{'variant':<26}{'unique':>8}{'% pyr':>7}{'paper/ours':>12}{'cell-dir':>10}")
    print("-" * 63)
    for label in VARIANTS:
        g = df[df.variant == label]
        u, p = int(g.unique_place_cells.sum()), int(g.cell_direction_pairs.sum())
        print(f"{label:<26}{u:8d}{100 * u / n_pyr:6.0f}%{PAPER_TOTAL / u:12.2f}{p:10d}")
    g = df[df.variant == "primary"]
    u, p = int(g.peak_rate_only_unique.sum()), int(g.peak_rate_only_pairs.sum())
    print(f"{'peak rate only (>= 1 Hz)':<26}{u:8d}{100 * u / n_pyr:6.0f}%"
          f"{PAPER_TOTAL / u:12.2f}{p:10d}")
    print("-" * 63)
    print(f"{'paper (Grosmark & Buzsaki)':<26}{PAPER_TOTAL:8d}{100 * PAPER_TOTAL / n_pyr:6.0f}%")
    print(f"\npyramidal cells across the 8 sessions: {n_pyr}")

    print("\nper session (primary variant):")
    print(f"{'session':<20}{'pyr':>5}{'unique':>8}{'cell-dir':>10}{'peak-only':>11}")
    for _, r in df[df.variant == "primary"].iterrows():
        print(f"{r.session:<20}{r.n_pyr:5d}{r.unique_place_cells:8d}"
              f"{r.cell_direction_pairs:10d}{r.peak_rate_only_unique:11d}")
    print(f"\nwrote {out.relative_to(config.repo_root())}")


if __name__ == "__main__":
    main()
