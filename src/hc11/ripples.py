"""Sharp-wave ripple detection from the local field potential (step 6b).

The detector is the one in :mod:`hc11.events`, fed a ripple-band amplitude
envelope instead of a population rate, so ripples and population-synchrony
events are defined by the same rule and can be compared directly in 6c.

Three things here are specific to the LFP and worth stating.

**The envelope is computed once over the whole recording, then sliced by
state.** Filtering each state interval separately would put a filter transient
at the edge of every interval, and non-REM intervals are short enough that
those edges would be a measurable fraction of the searched time. Thresholds are
still derived from, and detection still runs inside, one state at a time.

**Chunked, with overlap.** 43.6 million samples is 350 MB as float64 and the
Hilbert transform of that is 700 MB of complex128. :func:`channel_envelope`
works in chunks with an overlap that is discarded, which costs nothing in
accuracy -- the discarded margin is far longer than the filter's impulse
response -- and keeps peak memory at a few tens of MB.

**Multi-shank consistency is measured per event, not filtered on.** A real
sharp-wave ripple is a population event visible across neighbouring recording
sites; a movement artefact or a single bad channel is not. :func:`shank_support`
reports, for each event detected on the primary channel, how many *other* CA1
shanks also show an elevated ripple-band envelope within a few milliseconds.
Nothing is excluded on that basis: the count travels with the event so that it
can be compared against the MUA cross-check, which is the point of having it.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from scipy.ndimage import gaussian_filter1d

from hc11.events import EventList, detect_two_threshold
from hc11.ripple_channels import band_envelope, bandpass

#: Chunk and overlap for :func:`channel_envelope`, in samples. The overlap is
#: ~1.6 s at 1250 Hz, orders of magnitude longer than a 4th-order 150-300 Hz
#: filter's settling time, so discarding it removes every edge effect.
CHUNK_SAMPLES = 1 << 21
OVERLAP_SAMPLES = 1 << 11


@dataclass(frozen=True)
class ChannelEnvelope:
    """Ripple-band amplitude envelope of one channel, over a whole recording."""

    channel: int
    rate_hz: float
    band: tuple[float, float]
    values: np.ndarray  # (n_samples,) float32
    smoothing_sigma_ms: float = 0.0
    label: str = ""  # e.g. "level" or "excursion", for reporting

    def __len__(self) -> int:
        return len(self.values)

    @property
    def duration_s(self) -> float:
        return len(self.values) / self.rate_hz

    def slice(self, start_s: float, end_s: float) -> tuple[np.ndarray, np.ndarray]:
        """``(t, values)`` for one interval, clipped to the recording."""
        lo = max(int(round(start_s * self.rate_hz)), 0)
        hi = min(int(round(end_s * self.rate_hz)), len(self.values))
        if hi <= lo:
            return np.empty(0), np.empty(0)
        return (lo + np.arange(hi - lo)) / self.rate_hz, self.values[lo:hi]

    def baseline(self, intervals) -> tuple[float, float]:
        """Mean and s.d. of the envelope inside ``intervals``, pooled.

        Accumulated rather than concatenated: the intervals can run to hours and
        there is no reason to hold a copy of them.
        """
        total = total_sq = 0.0
        n = 0
        for start, end in np.atleast_2d(np.asarray(intervals, dtype=float)):
            _, values = self.slice(start, end)
            if values.size == 0:
                continue
            values = values.astype(np.float64)
            total += values.sum()
            total_sq += np.square(values).sum()
            n += values.size
        if n == 0:
            return float("nan"), float("nan")
        mean = total / n
        return mean, float(np.sqrt(max(total_sq / n - mean**2, 0.0)))


def channel_envelope(
    lfp,
    channel: int,
    band: tuple[float, float],
    smoothing_sigma_ms: float = 0.0,
    order: int = 4,
    chunk_samples: int = CHUNK_SAMPLES,
    overlap_samples: int = OVERLAP_SAMPLES,
    label: str = "",
) -> ChannelEnvelope:
    """Ripple-band envelope of one channel across the whole ``.eeg``."""
    rate = lfp.meta.sampling_rate_hz
    n = lfp.n_samples
    out = np.empty(n, dtype=np.float32)
    data = lfp.data
    position = 0
    while position < n:
        stop = min(position + chunk_samples, n)
        lo = max(position - overlap_samples, 0)
        hi = min(stop + overlap_samples, n)
        block = np.asarray(data[lo:hi, channel], dtype=np.float64)
        envelope = band_envelope(block, rate, band, order)
        out[position:stop] = envelope[position - lo : position - lo + (stop - position)]
        position = stop
    if smoothing_sigma_ms > 0:
        out = gaussian_filter1d(
            out, sigma=smoothing_sigma_ms / 1000.0 * rate, mode="nearest"
        ).astype(np.float32)
    return ChannelEnvelope(
        channel=int(channel), rate_hz=rate, band=tuple(band), values=out,
        smoothing_sigma_ms=smoothing_sigma_ms, label=label,
    )


def detect_ripples(
    envelope: ChannelEnvelope,
    intervals,
    session_name: str,
    epoch: str,
    params: dict,
) -> EventList:
    """Two-threshold detection on the envelope, inside ``intervals`` only.

    The baseline is pooled over all of the epoch's state intervals before any
    detection happens, then detection runs inside each interval separately so
    that no event can span a gap in the searched time -- the same two-pass
    structure as :func:`hc11.mua.detect_mua_events`, for the same reason.
    """
    cfg = params["events"]["ripples"]
    intervals = np.atleast_2d(np.asarray(intervals, dtype=float)) if len(intervals) else \
        np.empty((0, 2))
    state_seconds = float((intervals[:, 1] - intervals[:, 0]).sum()) if len(intervals) else 0.0
    mean, sd = envelope.baseline(intervals)
    detector = f"lfp:ch{envelope.channel}"
    if not np.isfinite(sd) or sd == 0:
        empty = np.empty(0)
        return EventList(
            session=session_name, epoch=epoch, detector=detector, start=empty,
            peak=empty.copy(), end=empty.copy(), peak_z=empty.copy(),
            baseline_mean=mean, baseline_sd=sd, state_seconds=state_seconds,
            params={**dict(cfg), "channel": envelope.channel},
        )

    high = mean + cfg["threshold_high_sd"] * sd
    low = mean + cfg["threshold_low_sd"] * sd
    starts, peaks, ends, peak_z = [], [], [], []
    for start, end in intervals:
        t, values = envelope.slice(start, end)
        if values.size == 0:
            continue
        s, p, e = detect_two_threshold(
            values, t, high=high, low=low,
            min_duration_s=cfg["min_duration_ms"] / 1000.0,
            max_duration_s=cfg["max_duration_ms"] / 1000.0,
            merge_gap_s=cfg["merge_gap_ms"] / 1000.0,
        )
        if len(s) == 0:
            continue
        starts.append(s)
        peaks.append(p)
        ends.append(e)
        peak_z.append((np.interp(p, t, values) - mean) / sd)

    def _cat(chunks):
        return np.concatenate(chunks) if chunks else np.empty(0)

    return EventList(
        session=session_name, epoch=epoch, detector=detector,
        start=_cat(starts), peak=_cat(peaks), end=_cat(ends), peak_z=_cat(peak_z),
        baseline_mean=mean, baseline_sd=sd, state_seconds=state_seconds,
        params={
            **dict(cfg), "channel": envelope.channel, "band": list(envelope.band),
            "threshold_high": high, "threshold_low": low,
            "channel_ranking": envelope.label,
        },
    )


# ---------------------------------------------------------------------------
# Multi-shank consistency
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class ShankSupport:
    """How many other CA1 shanks show the same event.

    ``z`` is ``(n_events, n_channels)``: the largest ripple-band envelope
    z-score on each comparison channel within ``window_s`` of the event's peak,
    against that channel's own baseline in the same brain state.
    """

    channels: np.ndarray  # (n_channels,) the comparison channels
    shanks: np.ndarray  # (n_channels,) their shank numbers
    z: np.ndarray  # (n_events, n_channels)
    window_s: float
    z_threshold: float
    baseline_mean: np.ndarray = field(default_factory=lambda: np.empty(0))
    baseline_sd: np.ndarray = field(default_factory=lambda: np.empty(0))

    @property
    def n_supporting(self) -> np.ndarray:
        """Per event: comparison channels above ``z_threshold``."""
        if self.z.size == 0:
            return np.zeros(self.z.shape[0], dtype=int)
        return (self.z >= self.z_threshold).sum(axis=1)

    @property
    def n_channels(self) -> int:
        return len(self.channels)

    def distribution(self) -> np.ndarray:
        """Counts of events by number of supporting shanks, length ``n+1``."""
        return np.bincount(self.n_supporting, minlength=self.n_channels + 1)


def windowed_peak_z(
    lfp,
    channels,
    times: np.ndarray,
    band: tuple[float, float],
    mean,
    sd,
    window_s: float,
    pad_s: float = 0.25,
    smoothing_sigma_ms: float = 0.0,
    order: int = 4,
) -> np.ndarray:
    """``(n_times, n_channels)`` largest envelope z within ``+/- window_s`` of each time.

    Reads one padded window per event covering *all* comparison channels at
    once and filters the whole block along time. Per event per channel would be
    the obvious loop and is an order of magnitude slower: the cost here is
    dominated by per-call overhead, not by arithmetic, and one read of
    ``(n, 11)`` is one call where eleven reads of ``(n, 1)`` are eleven.

    ``pad_s`` of context on each side is read and then discarded, so the
    filter's transient never lands inside the measured window.
    """
    times = np.asarray(times, dtype=float)
    channels = np.atleast_1d(np.asarray(channels, dtype=int))
    mean = np.atleast_1d(np.asarray(mean, dtype=float))
    sd = np.atleast_1d(np.asarray(sd, dtype=float))
    out = np.full((times.size, channels.size), np.nan)
    if times.size == 0:
        return out
    rate = lfp.meta.sampling_rate_hz
    margin = int(round(pad_s * rate))
    usable = np.isfinite(sd) & (sd > 0)
    for i, centre in enumerate(times):
        _, block = lfp.read(centre - pad_s - window_s, centre + pad_s + window_s, channels)
        if block.shape[0] < 4 * order + 1:
            continue
        envelope = band_envelope(block.astype(float), rate, band, order)
        if smoothing_sigma_ms > 0:
            envelope = gaussian_filter1d(
                envelope, sigma=smoothing_sigma_ms / 1000.0 * rate, axis=0, mode="nearest"
            )
        core = envelope[margin : envelope.shape[0] - margin] if margin > 0 else envelope
        if core.shape[0] == 0:
            continue
        peak = core.max(axis=0)
        out[i, usable] = (peak[usable] - mean[usable]) / sd[usable]
    return out


def baseline_windows(intervals, window_s: float, n_windows: int) -> np.ndarray:
    """Windows to estimate a comparison channel's baseline from.

    Normally a fixed, evenly spread sample of ``window_s`` windows. When no
    interval is that long -- MAZE immobility is fragments of under a second in
    these sessions -- the state itself is used instead, capped at the same
    total duration. The alternative is a silent failure: with no windows the
    baseline is NaN, every z-score is NaN, and "no other shank shows this
    event" is indistinguishable from "this was never measured".
    """
    from hc11 import ripple_channels as rc

    intervals = np.atleast_2d(np.asarray(intervals, dtype=float)) if len(intervals) else \
        np.empty((0, 2))
    if intervals.size == 0:
        return np.empty((0, 2))
    sampled = rc.sample_windows(intervals, window_s, n_windows)
    if len(sampled):
        return sampled
    budget = window_s * n_windows
    durations = intervals[:, 1] - intervals[:, 0]
    keep = np.cumsum(durations) <= budget
    if not keep.any():
        keep[0] = True
    return intervals[keep]


def shank_support(
    lfp,
    events: EventList,
    channels,
    shanks,
    intervals,
    params: dict,
    envelope_baselines=None,
) -> ShankSupport:
    """Per-event ripple-band elevation on other shanks' channels.

    ``envelope_baselines`` supplies ``(mean, sd)`` per comparison channel; when
    omitted they are estimated from ``intervals`` by sampling, because
    filtering every comparison channel over the whole recording costs minutes
    per session and buys nothing -- a baseline is a summary statistic.
    """
    from hc11 import ripple_channels as rc

    cfg = params["events"]["ripples"]
    band = tuple(cfg["band"])
    channels = np.asarray(channels, dtype=int)
    shanks = np.asarray(shanks, dtype=int)
    window_s = cfg["support_window_ms"] / 1000.0

    if envelope_baselines is None:
        windows = baseline_windows(
            intervals, cfg["baseline_window_s"], cfg["baseline_n_windows"]
        )
        block = rc.read_windows(lfp, windows, channels) if len(windows) else None
        means, sds = [], []
        for index in range(len(channels)):
            if block is None:
                means.append(float("nan"))
                sds.append(float("nan"))
                continue
            envelope = band_envelope(
                block[:, index].astype(float), lfp.meta.sampling_rate_hz, band
            )
            means.append(float(envelope.mean()))
            sds.append(float(envelope.std()))
        envelope_baselines = (np.array(means), np.array(sds))
    means, sds = (np.asarray(x, dtype=float) for x in envelope_baselines)

    z = windowed_peak_z(
        lfp, channels, events.peak, band, means, sds, window_s,
        smoothing_sigma_ms=cfg["smoothing_sigma_ms"],
    )
    return ShankSupport(
        channels=channels, shanks=shanks, z=z, window_s=window_s,
        z_threshold=cfg["support_z"], baseline_mean=means, baseline_sd=sds,
    )


def raw_and_filtered(lfp, channel: int, centre_s: float, pad_s: float, band) -> tuple:
    """``(raw, filtered)`` around one time, for an example figure."""
    _, block = lfp.read(centre_s - pad_s, centre_s + pad_s, [channel])
    raw = block[:, 0].astype(float)
    return raw, bandpass(raw, lfp.meta.sampling_rate_hz, tuple(band))
