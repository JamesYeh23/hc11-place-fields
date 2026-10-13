"""Step 6b: detect sharp-wave ripples on the LFP, on both candidate channels.

    uv run python scripts/detect_ripples.py
    uv run python scripts/detect_ripples.py Achilles_10252013 --no-figures

Reads results/channel_profiles.csv from step 6a-2, so run
scripts/select_ripple_channels.py first.

For each session this detects on the best channel under each of the two
rankings -- mean envelope level and p99/median excursion -- compares the two
event sets, and measures for every event on the primary channel how many other
CA1 shanks show the same ripple-band elevation. Shanks where the two rankings
disagree by more than a couple of channels are detected on as well, since those
are the cases where the ranking could actually matter.

Writes results/ripple_events.csv, results/ripple_summary.csv,
results/ripple_channel_comparison.csv and results/shank_support.csv, plus
per-session figures.
"""

from __future__ import annotations

import argparse
import json
import logging
from datetime import datetime, timezone

import numpy as np
import pandas as pd

from hc11 import config, plots, ripples
from hc11.events import compare_events
from hc11.io import load_session
from hc11.lfp import open_session_lfp
from hc11.mua import state_intervals_for_epoch
from hc11.preprocess import truncate_to_session

log = logging.getLogger("ripples")

SESSIONS = ["Achilles_10252013", "Achilles_11012013"]
EPOCHS = ("PRE", "MAZE", "POST")


def candidates(profiles: pd.DataFrame, session: str, stress_gap: int) -> dict:
    """The channels to detect on, from the step 6a-2 profile.

    ``level`` and ``excursion`` are the best channel overall under each
    ranking. ``stress`` holds the pair from any shank where the two rankings
    land more than ``stress_gap`` channels apart -- far enough that one of them
    is plausibly outside the pyramidal layer, which is the only case where the
    choice of ranking can change what is detected rather than shift it by a
    contact.
    """
    table = profiles[(profiles.session == session) & profiles.has_units].copy()
    if table.empty:
        raise SystemExit(f"{session}: no profiled channels on a shank with units")
    by_level = table.loc[table.groupby("shank").env_ripple.idxmax()]
    by_excursion = table.loc[table.groupby("shank").ripple_contrast.idxmax()]
    merged = by_level.merge(
        by_excursion[["shank", "channel"]], on="shank", suffixes=("_level", "_excursion")
    )
    stress = [
        (int(row.shank), int(row.channel_level), int(row.channel_excursion))
        for row in merged.itertuples()
        if abs(row.channel_level - row.channel_excursion) > stress_gap
    ]
    return {
        "level": int(by_level.loc[by_level.env_ripple.idxmax()].channel),
        "excursion": int(by_excursion.loc[by_excursion.ripple_contrast.idxmax()].channel),
        "per_shank_level": dict(zip(by_level.shank, by_level.channel, strict=True)),
        "per_shank_excursion": dict(zip(by_excursion.shank, by_excursion.channel, strict=True)),
        "stress": stress,
    }


def detect_on(lfp, session, channel: int, label: str, params: dict, intervals: dict):
    """Detect in every epoch on one channel; returns ``(envelope, {epoch: events})``."""
    cfg = params["events"]["ripples"]
    envelope = ripples.channel_envelope(
        lfp, channel, tuple(cfg["band"]),
        smoothing_sigma_ms=cfg["smoothing_sigma_ms"], label=label,
    )
    found = {}
    for epoch in EPOCHS:
        events = ripples.detect_ripples(envelope, intervals[epoch], session.name, epoch, params)
        found[epoch] = events
        log.info("  %-10s ch %3d  %s", label, channel, events)
    return envelope, found


def matched_indices(a, b, window_s: float) -> np.ndarray:
    """Indices of ``a``'s events that overlap, or peak close to, one of ``b``'s."""
    if len(a) == 0 or len(b) == 0:
        return np.empty(0, dtype=int)
    keep = [
        index for index in range(len(a))
        if (((b.start < a.end[index]) & (b.end > a.start[index]))
            | (np.abs(b.peak - a.peak[index]) <= window_s)).any()
    ]
    return np.array(keep, dtype=int)


def unmatched_indices(a, b, window_s: float) -> np.ndarray:
    """The complement of :func:`matched_indices` within ``a``."""
    matched = set(matched_indices(a, b, window_s).tolist())
    return np.array([i for i in range(len(a)) if i not in matched], dtype=int)


def subset_events(events, indices: np.ndarray):
    """A copy of ``events`` holding only ``indices``."""
    from dataclasses import replace

    indices = np.asarray(indices, dtype=int)
    return replace(
        events, start=events.start[indices], peak=events.peak[indices],
        end=events.end[indices], peak_z=events.peak_z[indices],
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("sessions", nargs="*", default=SESSIONS)
    parser.add_argument("--no-figures", action="store_true")
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(message)s")

    params = config.load_params()
    cfg = params["events"]["ripples"]
    out_root = config.results_dir()
    profile_path = out_root / "channel_profiles.csv"
    if not profile_path.exists():
        raise SystemExit(f"{profile_path} is missing; run scripts/select_ripple_channels.py first")
    profiles = pd.read_csv(profile_path)

    stamp = {
        "generated": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "git_commit": config.git_commit(),
        "ripples": cfg,
    }
    event_rows, summary_rows, comparison_rows, support_rows = [], [], [], []
    stress_rows = []

    for name in args.sessions or SESSIONS:
        log.info("\n=== %s", name)
        session = load_session(name)
        if params["preprocess"]["truncate_to_session"]:
            session, _ = truncate_to_session(session, verbose=False)
        lfp, meta = open_session_lfp(name)
        lfp.require_consistent_size(
            session.sess_duration, tolerance_s=params["lfp"]["size_check_tolerance_s"]
        )
        intervals = {e: state_intervals_for_epoch(session, e, params) for e in EPOCHS}
        picked = candidates(profiles, name, cfg["stress_gap_channels"])
        log.info(
            "candidates: level ch %d, excursion ch %d; stress pairs %s",
            picked["level"], picked["excursion"],
            picked["stress"] or "none",
        )

        runs = {}
        envelopes = {}
        for label, channel in (("level", picked["level"]), ("excursion", picked["excursion"])):
            envelopes[label], runs[label] = detect_on(
                lfp, session, channel, label, params, intervals
            )
        for shank, level_channel, excursion_channel in picked["stress"]:
            for label, channel in (
                (f"shank{shank}-level", level_channel),
                (f"shank{shank}-excursion", excursion_channel),
            ):
                envelopes[label], runs[label] = detect_on(
                    lfp, session, channel, label, params, intervals
                )

        # --- comparison of the two rankings, per epoch ---
        pairs = [("level", "excursion")] + [
            (f"shank{s}-level", f"shank{s}-excursion") for s, _, _ in picked["stress"]
        ]
        for first, second in pairs:
            for epoch in EPOCHS:
                a, b = runs[first][epoch], runs[second][epoch]
                overlap = compare_events(a, b, match_window_s=cfg["match_window_ms"] / 1000.0)
                comparison_rows.append({
                    "session": name, "epoch": epoch, "pair": f"{first} vs {second}",
                    "channel_a": a.params["channel"], "channel_b": b.params["channel"],
                    "n_a": overlap.n_a, "n_b": overlap.n_b,
                    "n_a_matched": overlap.n_a_matched, "n_b_matched": overlap.n_b_matched,
                    "fraction_a_matched": overlap.fraction_a_matched,
                    "fraction_b_matched": overlap.fraction_b_matched,
                    "median_offset_ms": (
                        1000 * float(np.median(overlap.offsets_s))
                        if overlap.offsets_s.size else float("nan")
                    ),
                    "rate_a_per_min": a.rate_per_min, "rate_b_per_min": b.rate_per_min,
                    "median_duration_a_ms": (
                        1000 * float(np.median(a.duration)) if len(a) else float("nan")
                    ),
                    "median_duration_b_ms": (
                        1000 * float(np.median(b.duration)) if len(b) else float("nan")
                    ),
                })
                log.info(
                    "  %s %s: %d vs %d events, %.0f%% of %s matched, median offset %s",
                    epoch, f"{first}/{second}", overlap.n_a, overlap.n_b,
                    100 * overlap.fraction_a_matched, first,
                    f"{1000 * np.median(overlap.offsets_s):+.1f} ms"
                    if overlap.offsets_s.size else "n/a",
                )

        # --- are the disagreeing events real? ---
        # The eye test is in the figure; this is the number behind it. An event
        # found on one channel and not the other is either a ripple the other
        # channel was too far from the layer to see, or noise. Multi-shank
        # support tells those apart without appealing to the trace.
        primary_per_shank = picked[f"per_shank_{cfg['primary_ranking']}"]
        for shank, _, _ in picked["stress"]:
            for first, second in ((f"shank{shank}-level", f"shank{shank}-excursion"),
                                  (f"shank{shank}-excursion", f"shank{shank}-level")):
                for epoch in EPOCHS:
                    a, b = runs[first][epoch], runs[second][epoch]
                    if len(a) == 0:
                        continue
                    window = cfg["match_window_ms"] / 1000.0
                    only = unmatched_indices(a, b, window)
                    if only.size == 0:
                        continue
                    channel = a.params["channel"]
                    others = [c for s2, c in primary_per_shank.items() if s2 != shank]
                    other_shanks = [s2 for s2 in primary_per_shank if s2 != shank]
                    subset = subset_events(a, only)
                    support = ripples.shank_support(
                        lfp, subset, others, other_shanks, intervals[epoch], params,
                    )
                    shared = subset_events(a, matched_indices(a, b, window))
                    shared_support = ripples.shank_support(
                        lfp, shared, others, other_shanks, intervals[epoch], params,
                    ) if len(shared) else None
                    stress_rows.append({
                        "session": name, "epoch": epoch, "shank": shank,
                        "channel": channel, "other_channel": b.params["channel"],
                        "ranking": "level" if first.endswith("level") else "excursion",
                        "n_events": len(a), "n_only_here": int(only.size),
                        "fraction_only_here": only.size / len(a),
                        "median_support_only": float(np.median(support.n_supporting)),
                        "median_support_shared": (
                            float(np.median(shared_support.n_supporting))
                            if shared_support is not None else float("nan")
                        ),
                        "fraction_only_with_majority_support": float(
                            (support.n_supporting >= support.n_channels / 2).mean()
                        ),
                        "median_peak_z_only": float(np.median(a.peak_z[only])),
                        "median_peak_z_all": float(np.median(a.peak_z)),
                    })
                    log.info(
                        "  %s shank %d ch %d only-events: %d of %d, median support %.0f/%d "
                        "(shared events %.0f/%d), median peak z %.2f vs %.2f overall",
                        epoch, shank, channel, only.size, len(a),
                        np.median(support.n_supporting), support.n_channels,
                        np.median(shared_support.n_supporting) if shared_support is not None
                        else float("nan"), support.n_channels,
                        np.median(a.peak_z[only]), np.median(a.peak_z),
                    )

        # --- multi-shank support for the primary channel ---
        primary_label = cfg["primary_ranking"]
        per_shank = picked[f"per_shank_{primary_label}"]
        primary_channel = picked[primary_label]
        comparison_channels = [c for c in per_shank.values() if c != primary_channel]
        comparison_shanks = [s for s, c in per_shank.items() if c != primary_channel]
        for epoch in EPOCHS:
            events = runs[primary_label][epoch]
            if len(events) == 0:
                continue
            support = ripples.shank_support(
                lfp, events, comparison_channels, comparison_shanks,
                intervals[epoch], params,
            )
            counts = support.distribution()
            if not np.isfinite(support.z).any():
                log.warning(
                    "  %s shank support could not be measured: no baseline for the "
                    "comparison channels in this state", epoch,
                )
            log.info(
                "  %s shank support (of %d other CA1 shanks): median %d, "
                "distribution %s",
                epoch, support.n_channels, int(np.median(support.n_supporting)),
                counts.tolist(),
            )
            support_rows.append({
                "session": name, "epoch": epoch, "channel": primary_channel,
                "n_events": len(events), "n_comparison_shanks": support.n_channels,
                "median_supporting": float(np.median(support.n_supporting)),
                "mean_supporting": float(support.n_supporting.mean()),
                "fraction_zero_support": float((support.n_supporting == 0).mean()),
                "fraction_majority_support": float(
                    (support.n_supporting >= support.n_channels / 2).mean()
                ),
                **{f"n_with_{k}": int(v) for k, v in enumerate(counts)},
            })
            for index in range(len(events)):
                event_rows.append({
                    "session": name, "epoch": epoch, "detector": events.detector,
                    "ranking": primary_label, "channel": primary_channel,
                    "start": events.start[index], "peak": events.peak[index],
                    "end": events.end[index],
                    "duration_ms": 1000 * events.duration[index],
                    "peak_z": events.peak_z[index],
                    "n_supporting_shanks": int(support.n_supporting[index]),
                    "median_support_z": (
                        float(np.nanmedian(support.z[index]))
                        if np.isfinite(support.z[index]).any() else float("nan")
                    ),
                })

        for label, found in runs.items():
            for epoch, events in found.items():
                summary_rows.append({
                    "session": name, "epoch": epoch, "ranking": label,
                    "channel": events.params["channel"], "n_events": len(events),
                    "state_min": events.state_seconds / 60,
                    "rate_per_min": events.rate_per_min,
                    "median_duration_ms": (
                        1000 * float(np.median(events.duration)) if len(events)
                        else float("nan")
                    ),
                    "median_peak_z": (
                        float(np.median(events.peak_z)) if len(events) else float("nan")
                    ),
                    "baseline_mean": events.baseline_mean,
                    "baseline_sd": events.baseline_sd,
                })

        if not args.no_figures:
            folder = config.results_dir(name)
            sidecar = {**stamp, "session": name, "candidates": {
                k: v for k, v in picked.items() if k in ("level", "excursion", "stress")
            }}
            plots.plot_ripple_rate_by_epoch(summary_rows, name, folder / "ripple_rates.png",
                                            sidecar)
            for first, second in pairs:
                a = runs[first]["POST"] if len(runs[first]["POST"]) else runs[first]["PRE"]
                b = runs[second]["POST"] if len(runs[second]["POST"]) else runs[second]["PRE"]
                examples = disagreement_examples(lfp, a, b, params, n=6)
                if examples:
                    tag = first.replace("-level", "")
                    plots.plot_disagreement_examples(
                        examples, meta.sampling_rate_hz, tuple(cfg["band"]), name,
                        a.params["channel"], b.params["channel"],
                        folder / f"ripple_disagreement_{tag}.png", sidecar,
                    )
            if support_rows:
                plots.plot_shank_support(
                    [r for r in support_rows if r["session"] == name], name,
                    folder / "ripple_shank_support.png", sidecar,
                )
        lfp.close()

    pd.DataFrame(event_rows).to_csv(out_root / "ripple_events.csv", index=False)
    pd.DataFrame(summary_rows).round(4).to_csv(out_root / "ripple_summary.csv", index=False)
    pd.DataFrame(comparison_rows).round(4).to_csv(
        out_root / "ripple_channel_comparison.csv", index=False
    )
    pd.DataFrame(support_rows).round(4).to_csv(out_root / "shank_support.csv", index=False)
    if stress_rows:
        pd.DataFrame(stress_rows).round(4).to_csv(
            out_root / "ripple_stress_pairs.csv", index=False
        )
    (out_root / "ripple_summary.json").write_text(json.dumps(stamp, indent=2) + "\n")
    log.info("\nwrote %s", (out_root / "ripple_summary.csv").relative_to(config.repo_root()))


def disagreement_examples(lfp, a, b, params: dict, n: int = 6):
    """Traces around events found on one channel and not the other.

    Both channels are shown for each example, so that "the other channel did
    not detect it" can be read off the trace rather than taken on trust.
    """
    cfg = params["events"]["ripples"]
    if len(a) == 0 or len(b) == 0:
        return []
    band = tuple(cfg["band"])
    pad = cfg["example_pad_ms"] / 1000.0
    window = cfg["match_window_ms"] / 1000.0
    channel_a, channel_b = a.params["channel"], b.params["channel"]

    def unmatched(first, second):
        keep = []
        for index in range(len(first)):
            start, end, peak = first.start[index], first.end[index], first.peak[index]
            hit = ((second.start < end) & (second.end > start)) | (
                np.abs(second.peak - peak) <= window
            )
            if not hit.any():
                keep.append(index)
        return np.array(keep, dtype=int)

    examples = []
    for first, second, owner in ((a, b, "a"), (b, a, "b")):
        missing = unmatched(first, second)
        if missing.size == 0:
            continue
        strongest = missing[np.argsort(first.peak_z[missing])[::-1][: n // 2]]
        for index in strongest:
            centre = float(first.peak[index])
            raw_a, filt_a = ripples.raw_and_filtered(lfp, channel_a, centre, pad, band)
            raw_b, filt_b = ripples.raw_and_filtered(lfp, channel_b, centre, pad, band)
            examples.append({
                "centre": centre, "found_on": owner, "epoch": first.epoch,
                "peak_z": float(first.peak_z[index]),
                "raw_a": raw_a, "filt_a": filt_a, "raw_b": raw_b, "filt_b": filt_b,
            })
    return examples


if __name__ == "__main__":
    main()
