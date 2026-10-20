"""Does the identity-shuffle residual track ensemble size rather than geometry?

    uv run python scripts/residual_vs_ensemble.py
    uv run python scripts/residual_vs_ensemble.py --event-set conjunction

D8.8 recorded that the negative control's residual is large in both circular
sessions *and* in one linear one, so track geometry does not explain it. The
obvious next suspect is ensemble size: Gatsby_08022013 is the smallest linear
ensemble as well as the largest linear residual, and if the residual follows
cell count then it is the same story the sensitivity-floor figure tells (D8.9)
and the two belong side by side.

This script answers that by ordering, not by a fit. **Five sessions cannot
establish a relationship**, so nothing here is regressed, no slope is quoted,
and the output says which readings are monotone with ensemble size and which
are not. Two quantities are reported, because they behave differently:

* the residual **as a fraction** of the real PRE/POST difference, which is how
  it is quoted everywhere else, and which is unstable wherever that difference
  is near zero (Gatsby_08022013's time permutation, a real gap of +0.014);
* the **absolute residual** in rZ units, which has no denominator and so says
  how noisy the control itself is.

Writes results/residual_vs_ensemble.csv and .png.
"""

from __future__ import annotations

import argparse
from datetime import datetime, timezone

import pandas as pd

from hc11 import config, plots
from hc11.results import leading_monotone_run

SHUFFLES = ("column_cycle", "time_permutation")


def gather(event_set: str) -> list[dict]:
    results = config.results_dir()
    comparison = pd.read_csv(results / "replay_comparison.csv")
    laps = pd.read_csv(results / "lap_summary.csv").set_index("session")
    per_session = comparison[~comparison.session.str.startswith("pooled")]
    here = per_session[per_session.event_set == event_set]
    if here.empty:
        raise SystemExit(f"no rows for event_set={event_set!r} in replay_comparison.csv")
    rows = [{
        "session": r.session,
        "n_place_cells": int(laps.loc[r.session].n_place_cells_unique),
        "geometry": r.geometry,
        "shuffle": r.shuffle,
        "difference": float(r.difference),
        "identity_residual": float(r.identity_residual),
        "residual_fraction": float(r.residual_fraction),
        "abs_residual": abs(float(r.identity_residual)),
        "fraction_interpretable": abs(float(r.difference)) >= plots.UNSTABLE_DIFFERENCE,
    } for r in here.itertuples()]
    rows.sort(key=lambda r: (-r["n_place_cells"], r["shuffle"]))
    return rows


def monotone_with_ensemble(rows: list[dict], shuffle: str, key: str) -> tuple[bool, int]:
    """Is ``key`` ordered by ensemble size, and over how many of the leading points?

    The run matters because one session can break an otherwise clean ordering,
    and that is a different statement from no ordering at all.
    """
    series = [r[key] for r in sorted((r for r in rows if r["shuffle"] == shuffle),
                                     key=lambda r: -r["n_place_cells"])]
    return leading_monotone_run(series)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--event-set", default="all", choices=("all", "conjunction"))
    args = parser.parse_args()

    rows = gather(args.event_set)
    header = (f"{'session':<20}{'cells':>6}{'geometry':>10}{'shuffle':>19}"
              f"{'real':>8}{'residual':>10}{'fraction':>10}")
    print(header)
    print("-" * len(header))
    for r in rows:
        flag = "" if r["fraction_interpretable"] else "  <- denominator too small"
        print(f"{r['session']:<20}{r['n_place_cells']:6d}{r['geometry']:>10}"
              f"{r['shuffle']:>19}{r['difference']:+8.3f}{r['identity_residual']:+10.3f}"
              f"{r['residual_fraction']:+10.3f}{flag}")
    print("-" * len(header))

    print("\nordered by ensemble size (largest first), does the residual rise as cells fall?")
    for shuffle in SHUFFLES:
        for key, label in (("residual_fraction", "fraction"), ("abs_residual", "|residual|")):
            whole, run = monotone_with_ensemble(rows, shuffle, key)
            verdict = ("monotone across all 5" if whole
                       else f"monotone over the {run} largest ensembles, then breaks")
            print(f"  {shuffle:<18}{label:<12}{verdict}")
    print("\nFive sessions cannot establish a relationship. No fit is made and none should be"
          "\nread into the figure; the ordering is reported, and that is all it is.")

    out_csv = config.results_dir() / f"residual_vs_ensemble_{args.event_set}.csv"
    pd.DataFrame(rows).round(4).to_csv(out_csv, index=False)
    out_png = plots.plot_residual_vs_ensemble(
        rows, config.results_dir() / f"residual_vs_ensemble_{args.event_set}.png",
        {"script": "scripts/residual_vs_ensemble.py",
         "git_commit": config.git_commit(),
         "created": datetime.now(timezone.utc).isoformat(timespec="seconds"),
         "event_set": args.event_set,
         "sessions": list(config.LFP_SESSIONS),
         "omitted": config.LFP_OMITTED,
         "note": "no curve is fitted; five sessions show an ordering at most"})
    print(f"\nwrote {out_csv.relative_to(config.repo_root())}")
    print(f"wrote {out_png.relative_to(config.repo_root())}")


if __name__ == "__main__":
    main()
