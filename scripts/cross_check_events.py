"""Step 6c: cross-check LFP ripples against population-synchrony (MUA) events.

    uv run python scripts/cross_check_events.py

Reads results/ripple_events.csv (step 6b) and results/mua_events.csv (the
earlier MUA run) and compares them as event lists: how many events each
detector finds, how many of each are matched by the other, and the distribution
of peak-to-peak offsets.

The two detectors see the same phenomenon through different instruments -- a
sharp-wave ripple is a transient in the 150-300 Hz field *and* a burst of
population firing -- so disagreement is informative in both directions. An
unmatched ripple is an LFP transient with no population burst behind it; an
unmatched MUA event is synchrony without a ripple. Neither is automatically an
error, and the point of the comparison is to find out which detector's
unmatched events look like the weak tail of the same distribution and which
look like a different population.

The multi-shank support measured in 6b is carried through, so the question the
brief asks -- whether the events few shanks see are the same ones the MUA
rejects -- is answered by a cross-tabulation rather than by inspection.

Writes results/event_cross_check.csv, results/cross_check_by_support.csv and
per-session figures.
"""

from __future__ import annotations

import argparse
import json
import logging
from datetime import datetime, timezone

import numpy as np
import pandas as pd

from hc11 import config, plots
from hc11.events import EventList, compare_events

log = logging.getLogger("crosscheck")

SESSIONS = ["Achilles_10252013", "Achilles_11012013"]
EPOCHS = ("PRE", "MAZE", "POST")


def event_list(table: pd.DataFrame, session: str, epoch: str, detector: str) -> EventList:
    """Rebuild an :class:`EventList` from a results CSV.

    Reconstructed rather than re-detected: the comparison has to be against the
    event sets that were actually reported, not against a fresh run that might
    differ because a parameter moved in between.
    """
    rows = table[(table.session == session) & (table.epoch == epoch)]
    return EventList(
        session=session, epoch=epoch, detector=detector,
        start=rows.start.to_numpy(), peak=rows.peak.to_numpy(), end=rows.end.to_numpy(),
        peak_z=rows.peak_z.to_numpy(), baseline_mean=float("nan"), baseline_sd=float("nan"),
        state_seconds=float("nan"),
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("sessions", nargs="*", default=SESSIONS)
    parser.add_argument("--no-figures", action="store_true")
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(message)s")

    params = config.load_params()
    window_s = params["events"]["ripples"]["match_window_ms"] / 1000.0
    out_root = config.results_dir()
    for required in ("ripple_events.csv", "mua_events.csv"):
        if not (out_root / required).exists():
            raise SystemExit(f"{out_root / required} is missing; run the earlier step first")
    ripples_table = pd.read_csv(out_root / "ripple_events.csv")
    mua_table = pd.read_csv(out_root / "mua_events.csv")

    stamp = {
        "generated": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "git_commit": config.git_commit(),
        "match_window_ms": params["events"]["ripples"]["match_window_ms"],
    }
    rows, support_rows, offset_store = [], [], {}

    for name in args.sessions or SESSIONS:
        log.info("\n=== %s", name)
        for epoch in EPOCHS:
            lfp_events = event_list(ripples_table, name, epoch, "lfp")
            mua_events = event_list(mua_table, name, epoch, "mua")
            if len(lfp_events) == 0 or len(mua_events) == 0:
                log.info("  %s: %d LFP, %d MUA -- nothing to compare",
                         epoch, len(lfp_events), len(mua_events))
                continue
            overlap = compare_events(lfp_events, mua_events, match_window_s=window_s)
            offsets_ms = 1000 * overlap.offsets_s
            offset_store[(name, epoch)] = offsets_ms
            rows.append({
                "session": name, "epoch": epoch,
                "n_lfp": overlap.n_a, "n_mua": overlap.n_b,
                "n_lfp_matched": overlap.n_a_matched, "n_mua_matched": overlap.n_b_matched,
                "fraction_lfp_matched": overlap.fraction_a_matched,
                "fraction_mua_matched": overlap.fraction_b_matched,
                "median_offset_ms": float(np.median(offsets_ms)) if offsets_ms.size else np.nan,
                "iqr_offset_ms": (
                    float(np.subtract(*np.percentile(offsets_ms, [75, 25])))
                    if offsets_ms.size else np.nan
                ),
                "match_window_ms": 1000 * window_s,
            })
            log.info(
                "  %-4s LFP %5d (%4.1f%% matched)  MUA %5d (%4.1f%% matched)  "
                "offset %+.1f ms (IQR %.1f)",
                epoch, overlap.n_a, 100 * overlap.fraction_a_matched,
                overlap.n_b, 100 * overlap.fraction_b_matched,
                np.median(offsets_ms) if offsets_ms.size else np.nan,
                np.subtract(*np.percentile(offsets_ms, [75, 25])) if offsets_ms.size else np.nan,
            )

            # --- does multi-shank support predict agreement with the MUA? ---
            rows_here = ripples_table[
                (ripples_table.session == name) & (ripples_table.epoch == epoch)
            ].reset_index(drop=True)
            matched = np.zeros(len(lfp_events), dtype=bool)
            for index in range(len(lfp_events)):
                hit = (
                    (mua_events.start < lfp_events.end[index])
                    & (mua_events.end > lfp_events.start[index])
                ) | (np.abs(mua_events.peak - lfp_events.peak[index]) <= window_s)
                matched[index] = hit.any()
            support = rows_here.n_supporting_shanks.to_numpy()
            for level in range(support.max() + 1):
                here = support == level
                if not here.any():
                    continue
                support_rows.append({
                    "session": name, "epoch": epoch, "n_supporting_shanks": level,
                    "n_events": int(here.sum()),
                    "fraction_of_events": float(here.mean()),
                    "fraction_matched_by_mua": float(matched[here].mean()),
                    "median_peak_z": float(np.median(rows_here.peak_z[here])),
                })
            low = support <= 2
            high = support >= 9
            log.info(
                "       support <=2: %d events, %.0f%% matched by MUA; "
                ">=9: %d events, %.0f%% matched",
                int(low.sum()), 100 * matched[low].mean() if low.any() else np.nan,
                int(high.sum()), 100 * matched[high].mean() if high.any() else np.nan,
            )

        if not args.no_figures:
            folder = config.results_dir(name)
            sidecar = {**stamp, "session": name}
            offsets = {e: offset_store[(name, e)] for e in EPOCHS if (name, e) in offset_store}
            if offsets:
                plots.plot_cross_check_offsets(
                    offsets, name, 1000 * window_s, folder / "cross_check_offsets.png", sidecar
                )
            here = [r for r in support_rows if r["session"] == name]
            if here:
                plots.plot_match_by_support(
                    here, name, folder / "cross_check_by_support.png", sidecar
                )

    pd.DataFrame(rows).round(4).to_csv(out_root / "event_cross_check.csv", index=False)
    pd.DataFrame(support_rows).round(4).to_csv(
        out_root / "cross_check_by_support.csv", index=False
    )
    (out_root / "event_cross_check.json").write_text(json.dumps(stamp, indent=2) + "\n")
    log.info("\nwrote %s", (out_root / "event_cross_check.csv").relative_to(config.repo_root()))


if __name__ == "__main__":
    main()
