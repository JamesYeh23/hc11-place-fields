"""Tests for LFP ripple detection.

The detector itself is `hc11.events.detect_two_threshold`, already covered; what
is new here is everything around it — the chunked envelope, restricting
detection to a brain state, and the multi-shank support measurement. The
fixtures plant ripples at known times on known channels so the assertions are
about *which* events come back, not just how many.
"""

from __future__ import annotations

import numpy as np
import pytest

from hc11 import ripples
from hc11.lfp import LfpMetadata

RATE = 1250.0
BAND = (150.0, 300.0)

PARAMS = {
    "events": {
        "ripples": {
            "band": [150.0, 300.0],
            "smoothing_sigma_ms": 8.0,
            "threshold_high_sd": 3.0,
            "threshold_low_sd": 1.0,
            "min_duration_ms": 20,
            "max_duration_ms": 400,
            "merge_gap_ms": 15,
            "support_window_ms": 20.0,
            "support_z": 2.0,
            "baseline_window_s": 2.0,
            "baseline_n_windows": 4,
            "match_window_ms": 50.0,
            "example_pad_ms": 250.0,
            "primary_ranking": "level",
            "stress_gap_channels": 2,
        }
    }
}


class FakeLfp:
    """An in-memory stand-in for `LfpFile` with the same read interface."""

    def __init__(self, data: np.ndarray, rate: float = RATE):
        self.data = data
        self.n_samples = data.shape[0]
        self.meta = LfpMetadata(
            session="synthetic", n_channels_total=data.shape[1], sampling_rate_hz=rate,
            wideband_rate_hz=20000.0, n_bits=16, voltage_range=20.0, amplification=1000.0,
            shanks={1: list(range(data.shape[1]))},
        )

    def read(self, start_s, end_s, channels=None):
        lo = max(int(round(start_s * self.meta.sampling_rate_hz)), 0)
        hi = min(int(round(end_s * self.meta.sampling_rate_hz)), self.n_samples)
        picked = slice(None) if channels is None else np.atleast_1d(np.asarray(channels, int))
        if hi <= lo:
            width = self.data.shape[1] if channels is None else np.size(picked)
            return np.empty(0), np.empty((0, width))
        block = np.asarray(self.data[lo:hi, picked])
        if block.ndim == 1:
            block = block[:, None]
        return (lo + np.arange(hi - lo)) / self.meta.sampling_rate_hz, block


def plant(n_samples: int, n_channels: int, times_s, channels, amplitude=120.0,
          duration_s=0.06, freq=190.0, noise=12.0, seed=0) -> np.ndarray:
    """Noise with Gaussian-windowed 190 Hz bursts at given times and channels."""
    rng = np.random.default_rng(seed)
    data = rng.standard_normal((n_samples, n_channels)) * noise
    width = int(duration_s * RATE)
    t = np.arange(width) / RATE
    window = np.exp(-0.5 * ((t - duration_s / 2) / (duration_s / 5)) ** 2)
    burst = amplitude * window * np.sin(2 * np.pi * freq * t)
    for time_s in np.atleast_1d(times_s):
        start = int(time_s * RATE)
        for channel in np.atleast_1d(channels):
            data[start:start + width, channel] += burst
    return data


@pytest.fixture
def planted():
    """60 s, 4 channels; ripples on ch 0 at 10/20/30 s, one of them shared."""
    times = [10.0, 20.0, 30.0]
    data = plant(int(60 * RATE), 4, times, [0])
    data = data + plant(int(60 * RATE), 4, [20.0], [1, 2], seed=1) - \
        np.random.default_rng(1).standard_normal((int(60 * RATE), 4)) * 12.0
    return data, times


# --- envelope -------------------------------------------------------------


def test_envelope_peaks_at_the_planted_ripples(planted):
    data, times = planted
    envelope = ripples.channel_envelope(FakeLfp(data), 0, BAND)
    for time_s in times:
        at = int(time_s * RATE) + int(0.03 * RATE)
        assert envelope.values[at] > 5 * np.median(envelope.values)


def test_chunking_does_not_change_the_envelope(planted):
    """The chunk boundary must be invisible: that is what the overlap is for.

    The first and last second of the *recording* are excluded. The analytic
    signal is computed by FFT, which treats its input as periodic, so the two
    ends of a trace carry a wrap-around artefact whose size depends on the
    length transformed. That is a property of the signal's own ends, not of
    chunking: inside a 1 s margin the two agree to 3e-5 of the envelope's
    maximum, while at the very first sample they differ by 1.4e-2 of it.
    """
    data, _ = planted
    whole = ripples.channel_envelope(FakeLfp(data), 0, BAND, chunk_samples=1 << 20)
    chunked = ripples.channel_envelope(
        FakeLfp(data), 0, BAND, chunk_samples=20000, overlap_samples=2048
    )
    inside = slice(int(RATE), -int(RATE))
    assert np.allclose(whole.values[inside], chunked.values[inside],
                       atol=1e-4 * whole.values.max())


def test_internal_chunk_boundaries_are_clean(planted):
    """Pin the quantity that matters: agreement *at* the boundaries themselves."""
    data, _ = planted
    chunk = 20000
    whole = ripples.channel_envelope(FakeLfp(data), 0, BAND, chunk_samples=1 << 20)
    chunked = ripples.channel_envelope(
        FakeLfp(data), 0, BAND, chunk_samples=chunk, overlap_samples=2048
    )
    for boundary in range(chunk, len(whole.values), chunk):
        window = slice(boundary - 200, boundary + 200)
        assert np.allclose(whole.values[window], chunked.values[window],
                           atol=1e-4 * whole.values.max()), boundary


def test_a_chunk_boundary_on_a_ripple_is_still_clean(planted):
    """Pin the failure mode directly: boundary placed inside a planted event."""
    data, _ = planted
    boundary = int(10.0 * RATE) + 20
    chunked = ripples.channel_envelope(
        FakeLfp(data), 0, BAND, chunk_samples=boundary, overlap_samples=1024
    )
    whole = ripples.channel_envelope(FakeLfp(data), 0, BAND, chunk_samples=1 << 20)
    window = slice(boundary - 100, boundary + 100)
    assert np.allclose(whole.values[window], chunked.values[window],
                       atol=1e-3 * whole.values.max())


def test_envelope_slice_and_baseline_use_the_intervals_given(planted):
    data, _ = planted
    envelope = ripples.channel_envelope(FakeLfp(data), 0, BAND)
    t, values = envelope.slice(5.0, 6.0)
    assert values.size == int(RATE)
    assert t[0] == pytest.approx(5.0)
    quiet, _ = envelope.baseline([[0.0, 5.0]])
    loud, _ = envelope.baseline([[9.99, 10.07]])
    assert loud > 3 * quiet


def test_smoothing_widens_the_envelope(planted):
    data, _ = planted
    rough = ripples.channel_envelope(FakeLfp(data), 0, BAND, smoothing_sigma_ms=0.0)
    smooth = ripples.channel_envelope(FakeLfp(data), 0, BAND, smoothing_sigma_ms=8.0)
    assert smooth.values.max() < rough.values.max()
    assert np.std(np.diff(smooth.values)) < np.std(np.diff(rough.values))


# --- detection ------------------------------------------------------------


def test_detection_finds_the_planted_ripples_and_nothing_else(planted):
    data, times = planted
    envelope = ripples.channel_envelope(FakeLfp(data), 0, BAND, smoothing_sigma_ms=8.0)
    events = ripples.detect_ripples(envelope, [[0.0, 60.0]], "synthetic", "PRE", PARAMS)
    assert len(events) == len(times)
    assert np.allclose(sorted(events.peak), times, atol=0.05)


def test_detection_searches_only_the_state_intervals(planted):
    data, _ = planted
    envelope = ripples.channel_envelope(FakeLfp(data), 0, BAND, smoothing_sigma_ms=8.0)
    events = ripples.detect_ripples(
        envelope, [[0.0, 15.0], [25.0, 40.0]], "synthetic", "PRE", PARAMS
    )
    assert len(events) == 2
    assert not np.any(np.abs(events.peak - 20.0) < 0.5)
    assert events.state_seconds == pytest.approx(30.0)


def test_an_event_cannot_span_a_gap_between_intervals(planted):
    """A ripple straddling two intervals must not be reported as one long event."""
    data, _ = planted
    envelope = ripples.channel_envelope(FakeLfp(data), 0, BAND, smoothing_sigma_ms=8.0)
    events = ripples.detect_ripples(
        envelope, [[0.0, 10.03], [10.03, 60.0]], "synthetic", "PRE", PARAMS
    )
    assert (events.duration <= PARAMS["events"]["ripples"]["max_duration_ms"] / 1000).all()


def test_the_baseline_comes_from_the_state_not_the_recording(planted):
    """A noisier stretch outside the state must not raise the threshold."""
    data, _ = planted
    loud = data.copy()
    loud[int(40 * RATE):] *= 8.0
    envelope = ripples.channel_envelope(FakeLfp(loud), 0, BAND, smoothing_sigma_ms=8.0)
    quiet_only = ripples.detect_ripples(
        envelope, [[0.0, 35.0]], "synthetic", "PRE", PARAMS
    )
    everything = ripples.detect_ripples(envelope, [[0.0, 60.0]], "synthetic", "PRE", PARAMS)
    assert quiet_only.baseline_sd < everything.baseline_sd
    assert len(quiet_only) >= 2


def test_rate_is_per_minute_of_searched_state(planted):
    data, _ = planted
    envelope = ripples.channel_envelope(FakeLfp(data), 0, BAND, smoothing_sigma_ms=8.0)
    events = ripples.detect_ripples(envelope, [[0.0, 60.0]], "synthetic", "PRE", PARAMS)
    assert events.rate_per_min == pytest.approx(len(events) * 60 / 60.0)


def test_a_flat_channel_yields_no_events_rather_than_raising():
    flat = np.zeros((int(10 * RATE), 2))
    envelope = ripples.channel_envelope(FakeLfp(flat), 0, BAND)
    events = ripples.detect_ripples(envelope, [[0.0, 10.0]], "synthetic", "PRE", PARAMS)
    assert len(events) == 0


def test_events_record_the_channel_and_ranking_they_came_from(planted):
    data, _ = planted
    envelope = ripples.channel_envelope(FakeLfp(data), 0, BAND, smoothing_sigma_ms=8.0,
                                        label="level")
    events = ripples.detect_ripples(envelope, [[0.0, 60.0]], "synthetic", "PRE", PARAMS)
    assert events.detector == "lfp:ch0"
    assert events.params["channel"] == 0
    assert events.params["channel_ranking"] == "level"


# --- multi-shank support --------------------------------------------------


def test_support_counts_the_channels_that_share_an_event(planted):
    data, _ = planted
    lfp = FakeLfp(data)
    envelope = ripples.channel_envelope(lfp, 0, BAND, smoothing_sigma_ms=8.0)
    events = ripples.detect_ripples(envelope, [[0.0, 60.0]], "synthetic", "PRE", PARAMS)
    support = ripples.shank_support(
        lfp, events, [1, 2, 3], [2, 3, 4], [[0.0, 60.0]], PARAMS
    )
    shared = int(np.argmin(np.abs(events.peak - 20.0)))
    alone = int(np.argmin(np.abs(events.peak - 10.0)))
    assert support.n_supporting[shared] >= 2
    assert support.n_supporting[alone] == 0


def test_support_distribution_covers_every_count(planted):
    data, _ = planted
    lfp = FakeLfp(data)
    envelope = ripples.channel_envelope(lfp, 0, BAND, smoothing_sigma_ms=8.0)
    events = ripples.detect_ripples(envelope, [[0.0, 60.0]], "synthetic", "PRE", PARAMS)
    support = ripples.shank_support(lfp, events, [1, 2, 3], [2, 3, 4], [[0.0, 60.0]], PARAMS)
    counts = support.distribution()
    assert len(counts) == 4
    assert counts.sum() == len(events)


def test_windowed_peak_z_reads_all_channels_at_once(planted):
    """The multi-channel result must equal the per-channel one, channel by channel."""
    data, _ = planted
    lfp = FakeLfp(data)
    times = np.array([10.03, 20.03, 30.03])
    mean, sd = np.array([1.0, 1.0, 1.0]), np.array([1.0, 1.0, 1.0])
    together = ripples.windowed_peak_z(lfp, [1, 2, 3], times, BAND, mean, sd, 0.02)
    for index, channel in enumerate([1, 2, 3]):
        alone = ripples.windowed_peak_z(
            lfp, [channel], times, BAND, mean[:1], sd[:1], 0.02
        )
        assert np.allclose(together[:, index], alone[:, 0])


def test_windowed_peak_z_is_nan_where_the_baseline_is_degenerate(planted):
    data, _ = planted
    z = ripples.windowed_peak_z(
        FakeLfp(data), [1, 2], np.array([10.0]), BAND, [1.0, 1.0], [1.0, 0.0], 0.02
    )
    assert np.isfinite(z[0, 0])
    assert np.isnan(z[0, 1])


def test_support_of_no_events_is_empty(planted):
    data, _ = planted
    lfp = FakeLfp(data)
    envelope = ripples.channel_envelope(lfp, 0, BAND)
    events = ripples.detect_ripples(envelope, [], "synthetic", "PRE", PARAMS)
    support = ripples.shank_support(lfp, events, [1, 2], [2, 3], [[0.0, 60.0]], PARAMS)
    assert support.z.shape == (0, 2)
    assert support.n_supporting.size == 0


# --- baseline windows -----------------------------------------------------


def test_baseline_windows_sample_a_long_state():
    windows = ripples.baseline_windows([[0.0, 600.0]], window_s=2.0, n_windows=4)
    assert len(windows) == 4
    assert np.allclose(windows[:, 1] - windows[:, 0], 2.0)


def test_baseline_windows_fall_back_to_short_fragments():
    """MAZE immobility is fragments under a second; sampling returns nothing.

    Falling back matters because the failure is otherwise silent: no windows
    means a NaN baseline, which makes "no other shank saw this" and "this was
    never measured" the same number.
    """
    fragments = [[10.0, 10.4], [20.0, 20.3], [30.0, 30.5]]
    windows = ripples.baseline_windows(fragments, window_s=15.0, n_windows=40)
    assert len(windows) == 3
    assert np.allclose(windows, fragments)


def test_baseline_windows_respect_the_time_budget():
    fragments = [[t, t + 1.0] for t in range(0, 100)]
    windows = ripples.baseline_windows(fragments, window_s=15.0, n_windows=4)
    assert (windows[:, 1] - windows[:, 0]).sum() <= 60.0


def test_baseline_windows_of_nothing_is_empty():
    assert ripples.baseline_windows([], window_s=15.0, n_windows=4).shape == (0, 2)


def test_support_in_a_short_state_is_measured_not_zero(planted):
    """The MAZE case end to end: fragments too short to sample must still work."""
    data, _ = planted
    lfp = FakeLfp(data)
    envelope = ripples.channel_envelope(lfp, 0, BAND, smoothing_sigma_ms=8.0)
    fragments = [[19.5, 20.4]]
    events = ripples.detect_ripples(envelope, fragments, "synthetic", "MAZE", PARAMS)
    assert len(events) == 1
    support = ripples.shank_support(lfp, events, [1, 2, 3], [2, 3, 4], fragments, PARAMS)
    assert np.isfinite(support.z).any()
