"""Does the envelope smoothing kernel explain the MUA-to-ripple lead?

    uv run python scripts/offset_vs_smoothing.py

Step 6c reported population firing peaking about 5 ms before the ripple
envelope peak, with the caveat that the envelope is smoothed with an 8 ms
Gaussian before thresholding and that smoothing can move a peak. This settles
it by sweeping the kernel width.

A symmetric kernel cannot shift a symmetric peak. It *can* shift an asymmetric
one, towards the heavier side, by an amount that grows with the kernel. So the
prediction if the lead is an artefact is simple: the offset should scale with
kernel width and go to zero as the kernel does.

Two measurements, because they answer different questions:

**fixed** -- the event set detected at the configured kernel is held constant,
and only the *peak time* is re-located inside each event window at each kernel
width. This isolates the kernel's effect on peak timing from its effect on
which events exist, and is the direct test.

**redetected** -- detection is re-run at each width and events are re-matched
to the MUA list. This is the end-to-end number: it is what the reported offset
would have been had the pipeline used that kernel throughout.

Writes results/offset_vs_smoothing.csv. No figure: the output is 28 numbers in
two dimensions and the only question asked of it is whether one of them varies,
which a table answers more directly than a chart would. Once the question is
settled the offset itself belongs in the 6c figure, not here.
"""

from __future__ import annotations

import argparse
import json
import logging
from datetime import datetime, timezone

import numpy as np
import pandas as pd
from scipy.ndimage import gaussian_filter1d

from hc11 import config, ripples
from hc11.events import EventList
from hc11.io import load_session
from hc11.lfp import open_session_lfp
from hc11.mua import state_intervals_for_epoch
from hc11.preprocess import truncate_to_session

log = logging.getLogger("offset")

SESSIONS = ["Achilles_10252013", "Achilles_11012013"]
EPOCHS = ("PRE", "POST")  # MAZE has ~20 events; the offset there is not estimable
KERNELS_MS = (0.0, 1.0, 2.0, 4.0, 8.0, 16.0, 32.0)


def mua_event_list(table: pd.DataFrame, session: str, epoch: str) -> EventList:
    rows = table[(table.session == session) & (table.epoch == epoch)]
    return EventList(
        session=session, epoch=epoch, detector="mua",
        start=rows.start.to_numpy(), peak=rows.peak.to_numpy(), end=rows.end.to_numpy(),
        peak_z=rows.peak_z.to_numpy(), baseline_mean=float("nan"),
        baseline_sd=float("nan"), state_seconds=float("nan"),
    )


def nearest_mua_offsets(events: EventList, mua: EventList, window_s: float) -> np.ndarray:
    """``mua.peak - events.peak`` for each event with a partner, in seconds."""
    if len(events) == 0 or len(mua) == 0:
        return np.empty(0)
    offsets = []
    for index in range(len(events)):
        start, end, peak = events.start[index], events.end[index], events.peak[index]
        hits = ((mua.start < end) & (mua.end > start)) | (
            np.abs(mua.peak - peak) <= window_s
        )
        if hits.any():
            candidates = mua.peak[hits]
            offsets.append(candidates[np.argmin(np.abs(candidates - peak))] - peak)
    return np.array(offsets)


def repeak(envelope_values, rate_hz: float, events: EventList) -> EventList:
    """Same events, peak times re-located on a differently smoothed envelope."""
    from dataclasses import replace

    peaks = np.empty(len(events))
    for index, (start, end) in enumerate(zip(events.start, events.end, strict=True)):
        lo = max(int(round(start * rate_hz)), 0)
        hi = min(int(round(end * rate_hz)), len(envelope_values))
        if hi <= lo:
            peaks[index] = events.peak[index]
            continue
        peaks[index] = (lo + int(np.argmax(envelope_values[lo:hi]))) / rate_hz
    return replace(events, peak=peaks)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("sessions", nargs="*", default=SESSIONS)
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(message)s")

    params = config.load_params()
    cfg = params["events"]["ripples"]
    window_s = cfg["match_window_ms"] / 1000.0
    configured = cfg["smoothing_sigma_ms"]
    out_root = config.results_dir()
    mua_table = pd.read_csv(out_root / "mua_events.csv")
    rows = []

    for name in args.sessions or SESSIONS:
        log.info("\n=== %s", name)
        session = load_session(name)
        if params["preprocess"]["truncate_to_session"]:
            session, _ = truncate_to_session(session, verbose=False)
        lfp, meta = open_session_lfp(name)
        rate = meta.sampling_rate_hz
        channel = int(
            pd.read_csv(out_root / "ripple_summary.csv")
            .query("session == @name and ranking == 'level'").channel.iloc[0]
        )
        # One unsmoothed envelope per session; every kernel is a cheap filter of
        # it, rather than 7 passes over 12 GB.
        raw = ripples.channel_envelope(lfp, channel, tuple(cfg["band"]),
                                       smoothing_sigma_ms=0.0, label="level")
        log.info("channel %d, %d samples", channel, len(raw))

        for epoch in EPOCHS:
            intervals = state_intervals_for_epoch(session, epoch, params)
            mua = mua_event_list(mua_table, name, epoch)
            reference = None
            for kernel in KERNELS_MS:
                values = (
                    raw.values if kernel == 0 else
                    gaussian_filter1d(raw.values, sigma=kernel / 1000.0 * rate,
                                      mode="nearest").astype(np.float32)
                )
                smoothed = ripples.ChannelEnvelope(
                    channel=channel, rate_hz=rate, band=tuple(cfg["band"]),
                    values=values, smoothing_sigma_ms=kernel, label="level",
                )
                redetected = ripples.detect_ripples(
                    smoothed, intervals, name, epoch,
                    {"events": {"ripples": {**cfg, "smoothing_sigma_ms": kernel}}},
                )
                if kernel == configured:
                    reference = redetected
                offsets = nearest_mua_offsets(redetected, mua, window_s)
                rows.append({
                    "session": name, "epoch": epoch, "channel": channel,
                    "kernel_ms": kernel, "measurement": "redetected",
                    "n_events": len(redetected), "n_matched": offsets.size,
                    "median_offset_ms": 1000 * float(np.median(offsets)) if offsets.size
                    else np.nan,
                    "iqr_offset_ms": 1000 * float(
                        np.subtract(*np.percentile(offsets, [75, 25]))
                    ) if offsets.size else np.nan,
                })
            if reference is None:
                raise SystemExit(
                    f"configured kernel {configured} ms is not in the sweep; add it"
                )
            for kernel in KERNELS_MS:
                values = (
                    raw.values if kernel == 0 else
                    gaussian_filter1d(raw.values, sigma=kernel / 1000.0 * rate,
                                      mode="nearest").astype(np.float32)
                )
                moved = repeak(values, rate, reference)
                offsets = nearest_mua_offsets(moved, mua, window_s)
                rows.append({
                    "session": name, "epoch": epoch, "channel": channel,
                    "kernel_ms": kernel, "measurement": "fixed",
                    "n_events": len(moved), "n_matched": offsets.size,
                    "median_offset_ms": 1000 * float(np.median(offsets)) if offsets.size
                    else np.nan,
                    "iqr_offset_ms": 1000 * float(
                        np.subtract(*np.percentile(offsets, [75, 25]))
                    ) if offsets.size else np.nan,
                })
            for measurement in ("fixed", "redetected"):
                line = [
                    r for r in rows
                    if r["session"] == name and r["epoch"] == epoch
                    and r["measurement"] == measurement
                ]
                log.info(
                    "  %-4s %-10s %s", epoch, measurement,
                    "  ".join(f"{r['kernel_ms']:.0f}ms:{r['median_offset_ms']:+.1f}"
                              for r in line),
                )
        lfp.close()

    table = pd.DataFrame(rows)
    table.round(4).to_csv(out_root / "offset_vs_smoothing.csv", index=False)
    (out_root / "offset_vs_smoothing.json").write_text(json.dumps({
        "script": "scripts/offset_vs_smoothing.py", "git_commit": config.git_commit(),
        "created": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "kernels_ms": list(KERNELS_MS), "configured_kernel_ms": configured,
    }, indent=2) + "\n")

    unsmoothed = table.query("kernel_ms == 0 and measurement == 'fixed'").median_offset_ms
    log.info(
        "\nAt zero smoothing the offset is %s ms (median over session x epoch); "
        "the configured %g ms kernel is in the sweep for comparison.",
        f"{unsmoothed.median():+.1f}", configured,
    )
    log.info("wrote %s", (out_root / "offset_vs_smoothing.csv").relative_to(config.repo_root()))


if __name__ == "__main__":
    main()
