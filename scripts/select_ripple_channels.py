"""Step 6a-2: pick the CA1 pyramidal-layer channel on each shank.

    uv run python scripts/select_ripple_channels.py
    uv run python scripts/select_ripple_channels.py Achilles_10252013 --no-figures

For each session this verifies the ``.eeg`` geometry against the XML and the
session clock, profiles every channel by band power over a fixed sample of
scored non-REM sleep, and selects the highest ripple-band channel on each
spike-sorted shank. Writes results/ripple_channels.csv, a per-session
results/<session>/channel_profile.csv, and two figures plus an example-trace
figure per session.

Nothing is detected here. The output is the channel each detector will run on
in step 6b, and the figures are there to be checked by eye first.
"""

from __future__ import annotations

import argparse
import json
import logging
from datetime import datetime, timezone

import numpy as np

from hc11 import config, plots
from hc11 import ripple_channels as rc
from hc11.io import load_session, shank_of
from hc11.lfp import open_session_lfp
from hc11.mua import state_intervals_for_epoch
from hc11.preprocess import truncate_to_session
from hc11.results import write_session_table

log = logging.getLogger("ripple_channels")

SESSIONS = list(config.LFP_SESSIONS)


def sleep_intervals(session, params) -> np.ndarray:
    """Scored non-REM in PRE and POST, which is where the layer is clearest.

    MAZE immobility is excluded from the *profiling* sample deliberately. It is
    a few minutes at most in these sessions (decision D6.3), and running theta
    bleeding into a short sample would tilt the band profile.
    """
    parts = [
        state_intervals_for_epoch(session, epoch, params) for epoch in ("PRE", "POST")
    ]
    parts = [p for p in parts if len(p)]
    return np.vstack(parts) if parts else np.empty((0, 2))


def shanks_with_units(session, meta) -> list[int]:
    """Spike groups that carry at least one sorted unit.

    hc-11's spike groups are not all shanks. Achilles' XML has 14 groups, and
    Channel_Orderings.pdf names two of them "Extra Channels" -- four channels
    each, outside the cell layer, with no sorted units. They are excluded here
    by the data rather than by their size: a group with no sorted pyramidal
    cell is not in *stratum pyramidale*, whatever the document calls it.
    """
    present = set(np.unique(shank_of(session.units())).tolist())
    return [s for s in sorted(meta.shanks) if s in present]


def profile_session(session_name: str, params: dict) -> dict:
    """Verify, profile and select, returning everything worth recording."""
    session = load_session(session_name)
    if params["preprocess"]["truncate_to_session"]:
        session, _ = truncate_to_session(session, verbose=False)
    lfp, meta = open_session_lfp(session_name)
    log.info("%s", meta)

    tolerance = params["lfp"]["size_check_tolerance_s"]
    check = lfp.require_consistent_size(session.sess_duration, tolerance_s=tolerance)
    log.info("size check: %s", check)

    cfg = params["lfp"]["channel_selection"]
    intervals = sleep_intervals(session, params)
    windows = rc.sample_windows(intervals, cfg["window_s"], cfg["n_windows"])
    if len(windows) == 0:
        raise SystemExit(f"{session_name}: no non-REM interval long enough to sample")
    log.info(
        "profiling on %d x %g s = %.1f min of non-REM (of %.1f min scored)",
        len(windows), cfg["window_s"], len(windows) * cfg["window_s"] / 60,
        (intervals[:, 1] - intervals[:, 0]).sum() / 60,
    )

    # Every channel is profiled, including the ones that cannot be selected:
    # seeing where the EMG and accelerometer sit on the profile is the check
    # that excluding them mattered.
    order = meta.probe_channels + meta.non_probe_channels
    block = rc.read_windows(lfp, windows, order)
    bands = {k: tuple(v) for k, v in cfg["bands"].items()}
    profile = rc.channel_profile(block, meta.sampling_rate_hz, channels=order, bands=bands)
    eligible = shanks_with_units(session, meta)
    skipped_groups = [s for s in sorted(meta.shanks) if s not in eligible]
    # What the excluded groups are depends on the probe. On the BUZ64SP
    # (Achilles, Cicero) they are the four-channel "Extra Channels" groups of
    # Channel_Orderings.pdf, which are not in the cell layer at all. On the
    # BUZ64 (Gatsby) every group is a real CA1 shank and the excluded ones
    # simply carry no sorted units. Either way the test is the same and it is
    # the data's, not the document's: no sorted pyramidal cell, no pyramidal
    # layer.
    log.info(
        "spike groups with sorted units: %s; without: %s (excluded -- no sorted "
        "units, so not in stratum pyramidale)",
        eligible, skipped_groups,
    )
    select = dict(skip_channels=cfg["skip_channels"], shanks=eligible)
    choices = rc.select_per_shank(profile, meta, statistic="mean", **select)
    sensitivity = rc.select_per_shank(profile, meta, statistic="contrast", **select)
    best = rc.best_channel(choices)
    for choice in choices:
        log.info("%s", choice)
    disagree = rc.agreement(choices, sensitivity)
    log.info(
        "contrast (p99/median) sensitivity: %s",
        "same channel on every shank" if not disagree
        else ", ".join(f"shank {s}: mean {a} vs contrast {b}" for s, (a, b) in disagree.items()),
    )
    log.info("best overall: shank %d, channel %d", best.shank, best.channel)

    return {
        "session": session_name,
        "lfp": lfp,
        "meta": meta,
        "check": check,
        "windows": windows,
        "profile": profile,
        "choices": choices,
        "sensitivity": sensitivity,
        "disagree": disagree,
        "eligible_shanks": eligible,
        "best": best,
        "n_nrem_s": float((intervals[:, 1] - intervals[:, 0]).sum()),
    }


def example_traces(result, params, n_examples: int, pad_ms: float):
    """``(time, raw, filtered)`` around the strongest ripple-band peaks.

    The peaks are located in the sampled non-REM block, then re-read from the
    recording at full context so each panel shows real neighbouring LFP rather
    than a window edge.
    """
    lfp, profile = result["lfp"], result["profile"]
    rate = result["meta"].sampling_rate_hz
    band = profile.bands["ripple"]
    channel = result["best"].channel

    block = rc.read_windows(lfp, result["windows"], [channel])[:, 0]
    envelope = rc.band_envelope(block.astype(float), rate, band)
    peaks = rc.top_envelope_peaks(envelope, rate, n_examples, refractory_s=1.0)

    # Map a position in the concatenated sample back onto the session clock.
    per_window = int(round(result["windows"][0][1] - result["windows"][0][0]) * rate)
    examples = []
    for index in peaks:
        window, within = divmod(int(index), per_window)
        centre = result["windows"][window][0] + within / rate
        pad = pad_ms / 1000.0
        _, raw = lfp.read(centre - pad, centre + pad, [channel])
        raw = raw[:, 0].astype(float)
        if len(raw) < 10:
            continue
        examples.append((centre, raw, rc.bandpass(raw, rate, band)))
    return examples, channel, band


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("sessions", nargs="*", default=SESSIONS)
    parser.add_argument("--no-figures", action="store_true")
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(message)s")

    params = config.load_params()
    cfg = params["lfp"]["channel_selection"]
    stamp = {
        "generated": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "git_commit": config.git_commit(),
        "channel_selection": cfg,
    }
    out_root = config.results_dir()
    rows, profile_rows = [], []

    names = args.sessions or SESSIONS
    for name in names:
        log.info("\n=== %s", name)
        result = profile_session(name, params)
        meta, profile = result["meta"], result["profile"]
        probe = set(meta.probe_channels)
        for shank, channels in sorted(meta.shanks.items()):
            for depth, channel in enumerate(channels):
                index = profile.index_of(channel)
                profile_rows.append({
                    "session": name, "shank": shank, "depth_in_shank": depth,
                    "channel": channel, "on_probe": True,
                    "sd": profile.sd[index],
                    **{f"env_{b}": profile.mean_envelope[b][index] for b in profile.bands},
                    "theta_delta": profile.theta_delta[index],
                    "ripple_z": profile.ripple_z()[index],
                    "ripple_contrast": profile.ripple_contrast[index],
                    "has_units": shank in result["eligible_shanks"],
                })
        for channel in meta.non_probe_channels:
            index = profile.index_of(channel)
            profile_rows.append({
                "session": name, "shank": np.nan, "depth_in_shank": np.nan,
                "channel": channel, "on_probe": False, "sd": profile.sd[index],
                **{f"env_{b}": profile.mean_envelope[b][index] for b in profile.bands},
                "theta_delta": profile.theta_delta[index],
                "ripple_z": profile.ripple_z()[index],
                "ripple_contrast": profile.ripple_contrast[index],
                "has_units": False,
            })
        for choice in result["choices"]:
            rows.append({
                "session": name, "shank": choice.shank, "channel": choice.channel,
                "depth_in_shank": meta.shanks[choice.shank].index(choice.channel),
                "env_ripple": choice.ripple, "ripple_z": choice.ripple_z,
                "theta_delta": choice.theta_delta, "sd": choice.sd,
                "runner_up": choice.runner_up, "margin": choice.margin,
                "ripple_contrast": profile.ripple_contrast[profile.index_of(choice.channel)],
                "contrast_choice": next(
                    (c.channel for c in result["sensitivity"] if c.shank == choice.shank),
                    None,
                ),
                "is_best_overall": choice.channel == result["best"].channel,
            })
        aux = [profile.ripple[profile.index_of(c)] for c in meta.non_probe_channels]
        log.info(
            "non-probe channels %s: ripple %s (max probe %.2f) -- %s",
            meta.non_probe_channels, [round(a, 1) for a in aux],
            profile.ripple[[profile.index_of(c) for c in probe]].max(),
            "excluded from selection",
        )

        if not args.no_figures:
            folder = config.results_dir(name)
            sidecar = {**stamp, "session": name, "size_check": str(result["check"]),
                       "n_windows": len(result["windows"]),
                       "n_nrem_s": result["n_nrem_s"]}
            plots.plot_channel_power_profile(profile, meta, result["choices"], name,
                                             folder / "lfp_channel_power.png", sidecar)
            plots.plot_band_profile(profile, meta, result["choices"], name,
                                    folder / "lfp_band_profile.png", sidecar)
            examples, channel, band = example_traces(
                result, params, cfg["n_example_events"], cfg["example_pad_ms"]
            )
            if examples:
                plots.plot_ripple_examples(
                    examples, meta.sampling_rate_hz, band, name, channel,
                    folder / "lfp_ripple_examples.png", sidecar,
                )
        result["lfp"].close()

    write_session_table(out_root / "ripple_channels.csv", rows, names)
    write_session_table(out_root / "channel_profiles.csv", profile_rows, names)
    (out_root / "ripple_channels.json").write_text(json.dumps(stamp, indent=2) + "\n")
    log.info("\nwrote %s", (out_root / "ripple_channels.csv").relative_to(config.repo_root()))


if __name__ == "__main__":
    main()
