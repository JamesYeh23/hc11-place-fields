"""DIAGNOSTIC: does the measurement base explain the place-cell gap?

Re-runs the primary place-cell criteria on a wider base -- our own linearisation
of the 2-D position over all MAZE time with valid tracking -- instead of scoring
only inside the authors' OneDLocation mask.

**This is not the replication base.** The authors' mask stays primary because it
is theirs; this only tests the hypothesis in decision D5.4. Results are written
to results/diagnostic_linearization.csv and are not used anywhere else.

    uv run python scripts/diagnostic_linearization.py
"""

from __future__ import annotations

import json
from datetime import datetime, timezone

import numpy as np
import pandas as pd

from hc11 import config
from hc11.behavior import running_mask_from_config
from hc11.fields import analyse_cell, select_running_samples
from hc11.io import load_session
from hc11.laps import DIRECTIONS, segment_laps_from_config
from hc11.linearize import coverage, with_own_linearization
from hc11.preprocess import truncate_to_session
from hc11.track import Track

PAPER_TOTAL = 491


def analyse(session, params):
    track = Track.from_session(session, params["place_fields"]["bin_size_cm"])
    running = running_mask_from_config(session, params)
    laps = segment_laps_from_config(session, track, params)
    rng = np.random.default_rng(params["place_fields"]["shuffle"]["random_seed"])
    cells = []
    for direction in DIRECTIONS:
        samples = select_running_samples(session, track, laps, running.mask, direction)
        if samples.n == 0:
            continue
        cells.extend(analyse_cell(session, track, samples, c, params, rng) for c in session.pyr_ids)
    unique = len({c.cluster_id for c in cells if c.is_place_cell})
    return {
        "unique": unique,
        "pairs": int(sum(c.is_place_cell for c in cells)),
        "running_min": running.seconds / 60,
        "n_laps": len(laps),
        "track_extent_m": track.extent,
    }


def main() -> None:
    params = config.load_params()
    rows = []
    for name in config.SESSIONS:
        raw = load_session(name)
        if params["preprocess"]["truncate_to_session"]:
            raw, _ = truncate_to_session(raw, verbose=False)
        cov = coverage(raw)
        authors = analyse(raw, params)
        own_session = with_own_linearization(raw)
        own = analyse(own_session, params)
        # How faithful our linearisation is to theirs, where both are defined.
        both = ~np.isnan(own_session.position_1d) & ~np.isnan(raw.position_1d)
        fidelity = (
            float(np.corrcoef(own_session.position_1d[both], raw.position_1d[both])[0, 1])
            if both.sum() > 10 else float("nan")
        )
        rows.append({
            "session": name, "maze_type": raw.maze_type, "n_pyr": len(raw.pyr_ids),
            "maze_kind": raw.maze_kind, "linearization_r": fidelity,
            "authors_cover_pct": 100 * cov["authors_frac"],
            "own_cover_pct": 100 * cov["own_frac"],
            "authors_run_min": authors["running_min"], "own_run_min": own["running_min"],
            "authors_laps": authors["n_laps"], "own_laps": own["n_laps"],
            "authors_unique": authors["unique"], "own_unique": own["unique"],
            "authors_pairs": authors["pairs"], "own_pairs": own["pairs"],
        })
        print(f"  {name} done", flush=True)

    df = pd.DataFrame(rows)
    out = config.results_dir() / "diagnostic_linearization.csv"
    df.round(4).to_csv(out, index=False)
    out.with_suffix(".json").write_text(json.dumps({
        "script": "scripts/diagnostic_linearization.py",
        "git_commit": config.git_commit(),
        "created": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "note": "diagnostic only; the authors' OneDLocation remains the replication base",
    }, indent=2) + "\n")

    header = (f"{'session':<18}{'pyr':>5}{'cover %':>17}{'running min':>16}{'laps':>12}"
              f"{'unique PC':>15}")
    print()
    print(header)
    print(f"{'':<18}{'':>5}{'theirs / ours':>17}{'theirs / ours':>16}"
          f"{'th / ours':>12}{'theirs / ours':>15}")
    print("-" * len(header))
    for _, r in df.iterrows():
        print(f"{r.session:<18}{r.n_pyr:5d}"
              f"{r.authors_cover_pct:8.0f}% /{r.own_cover_pct:5.0f}%"
              f"{r.authors_run_min:8.1f} /{r.own_run_min:6.1f}"
              f"{r.authors_laps:6d} /{r.own_laps:4d}"
              f"{r.authors_unique:8d} /{r.own_unique:5d}")
    print("-" * len(header))
    a, o, n = int(df.authors_unique.sum()), int(df.own_unique.sum()), int(df.n_pyr.sum())
    print(f"{'total':<18}{n:5d}{'':>17}{'':>16}{'':>12}{a:8d} /{o:5d}")
    print(f"\nall 8 sessions      : {a} -> {o} unique "
          f"({100 * a / n:.0f}% -> {100 * o / n:.0f}% of pyramidal)")

    # Our linearisation reproduces theirs exactly on linear tracks (r = 1.0000)
    # but only approximately on the ring, where the arc-length reconstruction is
    # cruder than their own -- so the circular rows measure our linearisation,
    # not the base. Report the linear sessions separately.
    lin = df[df.maze_kind == "linear"]
    la, lo, ln = (
        int(lin.authors_unique.sum()), int(lin.own_unique.sum()), int(lin.n_pyr.sum())
    )
    cir = df[df.maze_kind == "circular"]
    print(f"  linear sessions   : {la} -> {lo} unique "
          f"({100 * la / ln:.0f}% -> {100 * lo / ln:.0f}% of {ln} pyramidal), "
          f"linearisation r = {lin.linearization_r.min():.4f}-{lin.linearization_r.max():.4f}")
    print(f"  circular sessions : {int(cir.authors_unique.sum())} -> "
          f"{int(cir.own_unique.sum())} unique, linearisation r = "
          f"{cir.linearization_r.min():.3f}-{cir.linearization_r.max():.3f} "
          f"-- our arc-length reconstruction is worse than theirs, so these rows "
          f"measure our linearisation, not the base")
    print(f"\npaper               : {PAPER_TOTAL} unique "
          f"({100 * PAPER_TOTAL / n:.0f}% of pyramidal across all 8)")
    print(f"linear-session ceiling under this base: {100 * lo / ln:.0f}% vs the paper's 87%")
    print(f"\nwrote {out.relative_to(config.repo_root())}")


if __name__ == "__main__":
    main()
