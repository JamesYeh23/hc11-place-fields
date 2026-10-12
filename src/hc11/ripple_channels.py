"""Choosing the CA1 pyramidal-layer channel for ripple detection.

A 134-channel recording gives one ripple detector per shank at best, and the
shank's pyramidal layer has to be found before anything can be detected. The
standard criterion is power in the ripple band during non-REM sleep: ripples
are generated in *stratum pyramidale* and the band-passed envelope peaks there,
falling off above and below within a few hundred microns.

Two mistakes are easy enough here that the module is built to make them hard.

**Filtering the wrong axis.** The LFP array is ``(time, channel)``, and
``scipy.signal`` filters the *last* axis by default. Filtering across channels
instead of across time is not an error -- 134 channels is a long enough
"signal" for a 150-300 Hz band-pass to return numbers -- and the numbers it
returns look reasonable: a few tens of microvolts, slightly different per
channel. They are a spatial derivative of the probe geometry and carry no
information about where the pyramidal layer is, and the resulting profile is
nearly flat, which reads as "all channels are equivalent" rather than as a bug.
Every function here takes ``(time, channel)`` and passes ``axis=0``
explicitly, and :func:`channel_profile` loops one channel at a time over 1-D
arrays, where the ambiguity cannot arise.

**Letting a non-probe channel win.** The electromyogram has no pyramidal layer
but plenty of 150-300 Hz power. Selection runs over
:attr:`~hc11.lfp.LfpMetadata.probe_channels` only; see the ``hc11.lfp`` module
docstring for why that set is not the obvious one.

Band-pass is a zero-phase Butterworth in second-order-section form.
``sosfiltfilt`` rather than ``filtfilt`` because a narrow low-frequency band at
1250 Hz -- delta, 1-4 Hz -- is numerically unstable in transfer-function form
and returns NaN rather than an error.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.signal import butter, hilbert, sosfiltfilt

#: Bands used to describe a channel. ``ripple`` is the selection criterion;
#: the others are diagnostics that say whether the winner looks like a
#: pyramidal-layer channel rather than, say, a noisy one.
DEFAULT_BANDS: dict[str, tuple[float, float]] = {
    "delta": (1.0, 4.0),
    "theta": (5.0, 11.0),
    "ripple": (150.0, 300.0),
}


def bandpass(x: np.ndarray, rate_hz: float, band: tuple[float, float], order: int = 4):
    """Zero-phase Butterworth band-pass along **axis 0** (time).

    ``x`` may be 1-D or ``(time, channel)``. The axis is explicit and not the
    scipy default; see the module docstring.
    """
    low, high = band
    nyquist = 0.5 * rate_hz
    if not 0 < low < high < nyquist:
        raise ValueError(
            f"band {band} is not inside (0, {nyquist:g}) Hz for a {rate_hz:g} Hz signal"
        )
    sos = butter(order, [low / nyquist, high / nyquist], btype="bandpass", output="sos")
    return sosfiltfilt(sos, np.asarray(x, dtype=float), axis=0)


def band_envelope(
    x: np.ndarray, rate_hz: float, band: tuple[float, float], order: int = 4
) -> np.ndarray:
    """Analytic-signal amplitude envelope of ``x`` band-passed to ``band``."""
    return np.abs(hilbert(bandpass(x, rate_hz, band, order), axis=0))


@dataclass(frozen=True)
class ChannelProfile:
    """Per-channel summary of a stretch of LFP, used to pick a ripple channel.

    ``channels`` indexes the recording the same way the XML does. Every other
    array is parallel to it.
    """

    channels: np.ndarray  # (n_channels,) int, file numbering
    sd: np.ndarray  # (n_channels,) raw sample s.d., in file units
    mean_envelope: dict[str, np.ndarray]  # band -> (n_channels,) mean envelope
    median_envelope: dict[str, np.ndarray]  # band -> (n_channels,) median
    high_envelope: dict[str, np.ndarray]  # band -> (n_channels,) 99th percentile
    rate_hz: float
    n_samples: int
    bands: dict[str, tuple[float, float]]

    @property
    def ripple(self) -> np.ndarray:
        return self.mean_envelope["ripple"]

    @property
    def ripple_contrast(self) -> np.ndarray:
        """99th-percentile over median ripple envelope -- excursion, not level.

        The mean envelope answers "how much 150-300 Hz power is on this
        channel", which includes a baseline that is largely common to the
        probe. This answers "how far above its own baseline does this channel
        go", which is what a ripple actually is. The two agree when the layer
        is where the power is, and the sensitivity run in step 6a-2 reports
        when they do not.
        """
        return self.high_envelope["ripple"] / self.median_envelope["ripple"]

    @property
    def theta_delta(self) -> np.ndarray:
        """Theta/delta envelope ratio -- high in oriens/radiatum, low in the layer."""
        return self.mean_envelope["theta"] / self.mean_envelope["delta"]

    def ripple_z(self) -> np.ndarray:
        """Ripple power as a z-score across the channels in this profile."""
        r = self.ripple
        return (r - r.mean()) / r.std()

    def index_of(self, channel: int) -> int:
        hits = np.flatnonzero(self.channels == channel)
        if hits.size == 0:
            raise KeyError(f"channel {channel} is not in this profile")
        return int(hits[0])

    def restrict(self, channels) -> ChannelProfile:
        """A profile covering only ``channels`` (order preserved as given)."""
        keep = np.array([self.index_of(int(c)) for c in channels], dtype=int)
        return ChannelProfile(
            channels=self.channels[keep],
            sd=self.sd[keep],
            mean_envelope={k: v[keep] for k, v in self.mean_envelope.items()},
            median_envelope={k: v[keep] for k, v in self.median_envelope.items()},
            high_envelope={k: v[keep] for k, v in self.high_envelope.items()},
            rate_hz=self.rate_hz,
            n_samples=self.n_samples,
            bands=self.bands,
        )


def channel_profile(
    data: np.ndarray,
    rate_hz: float,
    channels=None,
    bands: dict[str, tuple[float, float]] | None = None,
    order: int = 4,
) -> ChannelProfile:
    """Mean band envelopes per channel for a ``(time, channel)`` block.

    One channel at a time, as a 1-D array: the per-channel loop is what keeps
    the time axis unambiguous, and it also keeps peak memory at one float64
    copy of a single channel rather than of the whole block.
    """
    data = np.asarray(data)
    if data.ndim != 2:
        raise ValueError(f"expected (time, channel), got shape {data.shape}")
    n_samples, n_channels = data.shape
    if n_samples < 3 * order + 1:
        raise ValueError(
            f"{n_samples} samples is too short to filter at order {order}"
        )
    bands = dict(bands if bands is not None else DEFAULT_BANDS)
    channels = (
        np.arange(n_channels) if channels is None else np.asarray(channels, dtype=int)
    )
    if channels.size != n_channels:
        raise ValueError(
            f"{channels.size} channel labels for {n_channels} columns of data"
        )

    sd = np.empty(n_channels)
    mean = {name: np.empty(n_channels) for name in bands}
    median = {name: np.empty(n_channels) for name in bands}
    high = {name: np.empty(n_channels) for name in bands}
    for col in range(n_channels):
        trace = np.asarray(data[:, col], dtype=float)
        sd[col] = trace.std()
        for name, band in bands.items():
            envelope = band_envelope(trace, rate_hz, band, order)
            mean[name][col] = envelope.mean()
            median[name][col] = np.median(envelope)
            high[name][col] = np.percentile(envelope, 99)
    return ChannelProfile(
        channels=channels,
        sd=sd,
        mean_envelope=mean,
        median_envelope=median,
        high_envelope=high,
        rate_hz=rate_hz,
        n_samples=n_samples,
        bands=bands,
    )


@dataclass(frozen=True)
class ShankChoice:
    """The channel chosen on one shank, with the numbers behind the choice."""

    shank: int
    channel: int
    ripple: float
    ripple_z: float
    theta_delta: float
    sd: float
    runner_up: int | None  # second-best channel on the shank
    margin: float  # winner / runner-up on the ranking statistic, 1.0 when alone
    statistic: str = "mean"  # which ranking produced this choice

    def __str__(self) -> str:
        runner = "-" if self.runner_up is None else str(self.runner_up)
        return (
            f"shank {self.shank:>2}: ch {self.channel:>3}  ripple {self.ripple:6.2f} "
            f"(z {self.ripple_z:+5.2f})  theta/delta {self.theta_delta:4.2f}  "
            f"sd {self.sd:6.1f}  runner-up {runner:>3} (x{self.margin:.2f})"
        )


#: Ranking statistics :func:`select_per_shank` accepts.
STATISTICS = {"mean": lambda p: p.ripple, "contrast": lambda p: p.ripple_contrast}


def select_per_shank(
    profile: ChannelProfile, meta, skip_channels=(), shanks=None, statistic: str = "mean"
) -> list[ShankChoice]:
    """Highest-ranking channel on each spike-sorted shank.

    Only ``meta.probe_channels`` are eligible, and only those present in
    ``profile``; a shank with no eligible channel is left out rather than
    filled with a guess. ``skip_channels`` drops individual channels on top of
    the XML's own skip flags. ``shanks`` restricts the shanks considered --
    hc-11's spike groups include non-shank "extra channel" groups that are not
    in the cell layer at all, and the caller is better placed than this
    function to say which groups those are.

    ``statistic`` is ``"mean"`` (mean ripple envelope, the stated criterion) or
    ``"contrast"`` (99th percentile over median, a sensitivity check).
    """
    if statistic not in STATISTICS:
        raise ValueError(f"unknown statistic {statistic!r}; expected one of {set(STATISTICS)}")
    rank = STATISTICS[statistic](profile)
    drop = set(int(c) for c in skip_channels) | set(meta.skipped_channels)
    available = set(int(c) for c in profile.channels)
    wanted = None if shanks is None else {int(s) for s in shanks}
    z = profile.ripple_z()
    choices: list[ShankChoice] = []
    for shank in sorted(meta.shanks):
        if wanted is not None and shank not in wanted:
            continue
        eligible = [c for c in meta.shanks[shank] if c not in drop and c in available]
        if not eligible:
            continue
        index = {c: profile.index_of(c) for c in eligible}
        ranked = sorted(eligible, key=lambda c: rank[index[c]], reverse=True)
        best = ranked[0]
        second = ranked[1] if len(ranked) > 1 else None
        margin = (
            1.0 if second is None else float(rank[index[best]] / rank[index[second]])
        )
        choices.append(
            ShankChoice(
                shank=shank,
                channel=int(best),
                ripple=float(profile.ripple[index[best]]),
                ripple_z=float(z[index[best]]),
                theta_delta=float(profile.theta_delta[index[best]]),
                sd=float(profile.sd[index[best]]),
                runner_up=None if second is None else int(second),
                margin=margin,
                statistic=statistic,
            )
        )
    return choices


def best_channel(choices: list[ShankChoice]) -> ShankChoice:
    """The strongest shank winner -- the single detection channel for a session."""
    if not choices:
        raise ValueError("no shank had an eligible channel")
    return max(choices, key=lambda c: c.ripple)


def agreement(first: list[ShankChoice], second: list[ShankChoice]) -> dict[int, tuple]:
    """Shanks where two selections differ, as ``{shank: (first, second)}``."""
    a = {c.shank: c.channel for c in first}
    b = {c.shank: c.channel for c in second}
    return {s: (a[s], b[s]) for s in sorted(set(a) & set(b)) if a[s] != b[s]}


# ---------------------------------------------------------------------------
# Sampling the recording
# ---------------------------------------------------------------------------


def sample_windows(intervals, window_s: float, n_windows: int) -> np.ndarray:
    """``n_windows`` windows of ``window_s``, spread evenly through ``intervals``.

    Ranking 134 channels does not need hours of LFP, and reading hours of 134
    channels is 10 GB. This picks a fixed sample instead: positions equally
    spaced along the *total* time inside ``intervals``, so long intervals
    contribute proportionally more windows and no window straddles a gap.

    Deterministic by construction -- no random seed to record, and two runs on
    the same session profile exactly the same samples. Intervals shorter than
    ``window_s`` are dropped, and fewer than ``n_windows`` come back if the
    state does not hold enough time.
    """
    intervals = np.atleast_2d(np.asarray(intervals, dtype=float))
    if intervals.size == 0:
        return np.empty((0, 2))
    usable = intervals[(intervals[:, 1] - intervals[:, 0]) >= window_s]
    if usable.size == 0:
        return np.empty((0, 2))
    # Each interval offers (duration - window_s) of possible start times.
    room = usable[:, 1] - usable[:, 0] - window_s
    total = room.sum()
    n = int(min(n_windows, np.floor(total / window_s) + len(usable)))
    if n <= 0:
        return np.empty((0, 2))
    edges = np.cumsum(room)
    # Midpoints of n equal slices of the available start-time range, so the
    # first and last windows sit inside the recording rather than on its edges.
    targets = (np.arange(n) + 0.5) * total / n if total > 0 else np.zeros(n)
    which = np.searchsorted(edges, targets, side="right").clip(0, len(usable) - 1)
    before = np.concatenate([[0.0], edges[:-1]])
    starts = usable[which, 0] + (targets - before[which])
    return np.column_stack([starts, starts + window_s])


def read_windows(lfp, windows, channels=None) -> np.ndarray:
    """Concatenate ``lfp.read`` over ``windows`` into one ``(time, channel)`` block.

    The windows are not contiguous in the recording, so the result has no
    meaningful time axis and must only be used for statistics that do not span
    window boundaries -- which a band-pass envelope does not, apart from a few
    samples of filter edge per window.
    """
    blocks = [lfp.read(start, stop, channels)[1] for start, stop in np.asarray(windows)]
    if not blocks:
        raise ValueError("no windows to read")
    return np.concatenate(blocks, axis=0)


def top_envelope_peaks(
    envelope: np.ndarray, rate_hz: float, n_peaks: int, refractory_s: float = 0.5
) -> np.ndarray:
    """Sample indices of the ``n_peaks`` largest envelope peaks, kept apart.

    A diagnostic for eyeballing whether the chosen channel's strongest moments
    look like ripples. It is deliberately *not* event detection: no threshold,
    no duration criterion, no state restriction -- just the largest peaks, so
    that what the figure shows is the best case rather than a thresholded
    sample of it.
    """
    envelope = np.asarray(envelope, dtype=float).ravel()
    gap = max(int(refractory_s * rate_hz), 1)
    order = np.argsort(envelope)[::-1]
    picked: list[int] = []
    for index in order:
        if len(picked) >= n_peaks:
            break
        if all(abs(index - other) >= gap for other in picked):
            picked.append(int(index))
    return np.array(sorted(picked), dtype=int)
