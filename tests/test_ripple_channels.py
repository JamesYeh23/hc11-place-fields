"""Tests for ripple-band channel selection.

The point of most of these is the time axis. A band-pass applied across
channels instead of across time returns plausible numbers rather than an
error, so the tests are built so that an axis mistake changes the *answer*,
not just the arithmetic: synthetic data where one known channel carries a
ripple and the rest do not, and an assertion that that channel is the one
selected.
"""

from __future__ import annotations

import numpy as np
import pytest

from hc11 import ripple_channels as rc
from hc11.lfp import LfpMetadata

RATE = 1250.0


def pink_noise(n_samples: int, n_channels: int, rng) -> np.ndarray:
    """Broadband noise with a 1/f tilt, like real LFP and unlike white noise."""
    white = rng.standard_normal((n_samples, n_channels))
    spectrum = np.fft.rfft(white, axis=0)
    freqs = np.fft.rfftfreq(n_samples, 1 / RATE)
    scale = np.ones_like(freqs)
    scale[1:] = 1 / freqs[1:]
    return np.fft.irfft(spectrum * scale[:, None], n=n_samples, axis=0) * 200.0


def with_ripples(
    data: np.ndarray, channel: int, rate: float = RATE, n_events: int = 12,
    freq: float = 190.0, amplitude: float = 40.0, duration_s: float = 0.06,
    rng=None,
) -> np.ndarray:
    """Add Gaussian-windowed 190 Hz bursts to one channel, in place on a copy."""
    rng = rng if rng is not None else np.random.default_rng(0)
    out = data.copy()
    n = out.shape[0]
    width = int(duration_s * rate)
    for start in rng.choice(n - width, size=n_events, replace=False):
        t = np.arange(width) / rate
        window = np.exp(-0.5 * ((t - duration_s / 2) / (duration_s / 5)) ** 2)
        out[start:start + width, channel] += amplitude * window * np.sin(
            2 * np.pi * freq * t
        )
    return out


def fake_meta(n_shanks: int = 4, per_shank: int = 8, n_aux: int = 3) -> LfpMetadata:
    shanks = {
        s: list(range((s - 1) * per_shank, s * per_shank))
        for s in range(1, n_shanks + 1)
    }
    total = n_shanks * per_shank + n_aux
    return LfpMetadata(
        session="synthetic",
        n_channels_total=total,
        sampling_rate_hz=RATE,
        wideband_rate_hz=20000.0,
        n_bits=16,
        voltage_range=20.0,
        amplification=1000.0,
        shanks=shanks,
        anatomical_groups={**shanks, n_shanks + 1: list(range(n_shanks * per_shank, total))},
    )


@pytest.fixture
def noise():
    return pink_noise(12000, 35, np.random.default_rng(7))


# --- band-pass ------------------------------------------------------------


def test_bandpass_keeps_its_band_and_drops_the_rest():
    t = np.arange(5000) / RATE
    signal = np.sin(2 * np.pi * 200 * t) + 3 * np.sin(2 * np.pi * 7 * t)
    filtered = rc.bandpass(signal, RATE, (150, 300))
    # The 7 Hz component is three times larger in the raw signal and must not
    # survive; the 200 Hz one must come through close to unit amplitude.
    assert filtered.std() == pytest.approx(1 / np.sqrt(2), rel=0.1)


def test_bandpass_is_zero_phase():
    t = np.arange(5000) / RATE
    signal = np.sin(2 * np.pi * 200 * t)
    filtered = rc.bandpass(signal, RATE, (150, 300))
    mid = slice(1000, 4000)
    lag = np.argmax(np.correlate(filtered[mid], signal[mid], mode="same"))
    assert lag == len(signal[mid]) // 2


def test_bandpass_filters_time_not_channels():
    """A column-wise constant block has nothing to filter along time.

    Each channel here is a different constant, so there is strong variation
    across channels and none across time. Filtering time must return ~zero;
    filtering channels would return the spatial edges instead.
    """
    block = np.tile(np.arange(40.0), (3000, 1))
    assert np.abs(rc.bandpass(block, RATE, (150, 300))).max() < 1e-6


def test_bandpass_rejects_a_band_outside_nyquist():
    with pytest.raises(ValueError, match="not inside"):
        rc.bandpass(np.zeros(2000), RATE, (150, 900))


def test_delta_band_is_finite():
    """1-4 Hz at 1250 Hz is where a transfer-function filter returns NaN."""
    envelope = rc.band_envelope(pink_noise(8000, 1, np.random.default_rng(1))[:, 0],
                                RATE, (1.0, 4.0))
    assert np.isfinite(envelope).all()
    assert envelope.mean() > 0


# --- profiles -------------------------------------------------------------


def test_profile_finds_the_channel_carrying_ripples(noise):
    data = with_ripples(noise, channel=11)
    profile = rc.channel_profile(data, RATE)
    assert int(np.argmax(profile.ripple)) == 11


def test_profile_is_not_flat_across_channels(noise):
    """The symptom of filtering the channel axis is a near-constant profile."""
    profile = rc.channel_profile(with_ripples(noise, channel=11), RATE)
    spread = profile.ripple.std() / profile.ripple.mean()
    assert spread > 0.05


def test_profile_envelopes_are_all_finite(noise):
    profile = rc.channel_profile(noise, RATE)
    for band, values in profile.mean_envelope.items():
        assert np.isfinite(values).all(), band
        assert (values > 0).all(), band


def test_profile_channel_labels_need_not_start_at_zero(noise):
    labels = np.arange(100, 100 + noise.shape[1])
    profile = rc.channel_profile(with_ripples(noise, channel=4), RATE, channels=labels)
    assert profile.channels[int(np.argmax(profile.ripple))] == 104
    assert profile.index_of(104) == 4


def test_profile_rejects_mismatched_labels(noise):
    with pytest.raises(ValueError, match="channel labels"):
        rc.channel_profile(noise, RATE, channels=[0, 1, 2])


def test_profile_rejects_one_dimensional_input():
    with pytest.raises(ValueError, match=r"\(time, channel\)"):
        rc.channel_profile(np.zeros(4000), RATE)


def test_profile_rejects_a_block_too_short_to_filter():
    with pytest.raises(ValueError, match="too short"):
        rc.channel_profile(np.zeros((8, 4)), RATE)


def test_restrict_keeps_the_requested_channels(noise):
    profile = rc.channel_profile(noise, RATE)
    subset = profile.restrict([5, 2, 9])
    assert subset.channels.tolist() == [5, 2, 9]
    assert subset.ripple[0] == profile.ripple[5]


# --- selection ------------------------------------------------------------


def test_selection_picks_the_ripple_channel_on_each_shank(noise):
    meta = fake_meta()
    data = noise.copy()
    planted = {1: 3, 2: 9, 3: 20, 4: 25}
    for channel in planted.values():
        data = with_ripples(data, channel)
    profile = rc.channel_profile(data, RATE)
    choices = {c.shank: c.channel for c in rc.select_per_shank(profile, meta)}
    assert choices == planted


def test_selection_ignores_non_probe_channels(noise):
    """An EMG channel with the strongest ripple-band power must not be picked."""
    meta = fake_meta()
    aux = meta.non_probe_channels[0]
    data = with_ripples(with_ripples(noise, channel=3), aux, amplitude=400.0)
    profile = rc.channel_profile(data, RATE)
    chosen = [c.channel for c in rc.select_per_shank(profile, meta)]
    assert aux not in chosen
    assert rc.best_channel(rc.select_per_shank(profile, meta)).channel == 3


def test_selection_honours_the_xml_skip_flag(noise):
    meta = fake_meta()
    data = with_ripples(noise, channel=3)
    profile = rc.channel_profile(data, RATE)
    assert rc.select_per_shank(profile, meta)[0].channel == 3
    skipping = LfpMetadata(**{**meta.__dict__, "skipped_channels": (3,)})
    assert rc.select_per_shank(profile, skipping)[0].channel != 3


def test_selection_honours_an_explicit_skip_list(noise):
    meta = fake_meta()
    profile = rc.channel_profile(with_ripples(noise, channel=3), RATE)
    assert rc.select_per_shank(profile, meta, skip_channels=[3])[0].channel != 3


def test_a_shank_with_no_profiled_channel_is_left_out(noise):
    meta = fake_meta()
    profile = rc.channel_profile(noise, RATE).restrict(meta.shanks[1] + meta.shanks[2])
    assert [c.shank for c in rc.select_per_shank(profile, meta)] == [1, 2]


def test_margin_reports_how_clear_the_win_was(noise):
    meta = fake_meta()
    profile = rc.channel_profile(with_ripples(noise, channel=3, amplitude=200.0), RATE)
    winner = rc.select_per_shank(profile, meta)[0]
    assert winner.channel == 3
    assert winner.margin > 1.5
    assert winner.runner_up != 3


def test_best_channel_needs_a_choice():
    with pytest.raises(ValueError, match="no shank"):
        rc.best_channel([])


# --- sampling -------------------------------------------------------------


def test_sampled_windows_stay_inside_their_interval():
    windows = rc.sample_windows([[10, 100], [500, 560]], window_s=15, n_windows=6)
    for start, stop in windows:
        assert stop - start == pytest.approx(15)
        assert (10 <= start and stop <= 100) or (500 <= start and stop <= 560)


def test_sampling_is_deterministic():
    first = rc.sample_windows([[0, 1000]], window_s=10, n_windows=8)
    second = rc.sample_windows([[0, 1000]], window_s=10, n_windows=8)
    assert np.array_equal(first, second)


def test_longer_intervals_get_more_windows():
    windows = rc.sample_windows([[0, 900], [1000, 1100]], window_s=10, n_windows=10)
    in_long = ((windows[:, 0] >= 0) & (windows[:, 1] <= 900)).sum()
    assert in_long >= 7


def test_intervals_shorter_than_the_window_are_dropped():
    windows = rc.sample_windows([[0, 5], [100, 160]], window_s=15, n_windows=4)
    assert (windows[:, 0] >= 100).all()


def test_sampling_returns_nothing_when_no_interval_fits():
    assert rc.sample_windows([[0, 5]], window_s=15, n_windows=4).shape == (0, 2)
    assert rc.sample_windows([], window_s=15, n_windows=4).shape == (0, 2)


def test_sampling_is_capped_by_the_time_available():
    windows = rc.sample_windows([[0, 40]], window_s=15, n_windows=50)
    assert 0 < len(windows) <= 3


# --- peak picking ---------------------------------------------------------


def test_peaks_are_the_largest_and_respect_the_refractory_gap():
    envelope = np.zeros(10000)
    envelope[[1000, 1010, 4000, 7000]] = [10.0, 9.9, 8.0, 7.0]
    peaks = rc.top_envelope_peaks(envelope, RATE, n_peaks=3, refractory_s=0.5)
    assert peaks.tolist() == [1000, 4000, 7000]


def test_peaks_stop_when_there_are_not_enough():
    envelope = np.zeros(1000)
    envelope[500] = 1.0
    assert len(rc.top_envelope_peaks(envelope, RATE, n_peaks=5, refractory_s=0.5)) == 1


# --- ranking statistic ----------------------------------------------------


def test_contrast_ranks_excursion_not_level(noise):
    """A channel with constant band power loses to one with the same mean in bursts.

    Added 150-300 Hz noise on one channel and ripple bursts of matched total
    power on another: the mean envelope cannot tell them apart, and that is
    exactly the case the contrast statistic exists for.
    """
    meta = fake_meta()
    data = noise.copy()
    rng = np.random.default_rng(3)
    steady = rc.bandpass(rng.standard_normal((len(data), 1)), RATE, (150, 300))[:, 0]
    data[:, 1] += steady * 60
    data = with_ripples(data, channel=2, n_events=30, amplitude=90.0)
    profile = rc.channel_profile(data, RATE)
    # The steady channel is at least as loud on average ...
    assert profile.ripple[1] > 0.8 * profile.ripple[2]
    # ... but the bursty one goes much further above its own baseline.
    assert profile.ripple_contrast[2] > profile.ripple_contrast[1]
    by_contrast = rc.select_per_shank(profile, meta, statistic="contrast")
    assert by_contrast[0].channel == 2


def test_an_unknown_statistic_is_refused(noise):
    with pytest.raises(ValueError, match="unknown statistic"):
        rc.select_per_shank(rc.channel_profile(noise, RATE), fake_meta(), statistic="rms")


def test_choices_record_which_statistic_ranked_them(noise):
    profile = rc.channel_profile(noise, RATE)
    assert rc.select_per_shank(profile, fake_meta())[0].statistic == "mean"
    assert rc.select_per_shank(
        profile, fake_meta(), statistic="contrast"
    )[0].statistic == "contrast"


def test_shanks_can_be_restricted(noise):
    """The hc-11 'Extra Channels' groups are spike groups but not cell layer."""
    profile = rc.channel_profile(noise, RATE)
    chosen = rc.select_per_shank(profile, fake_meta(), shanks=[2, 4])
    assert [c.shank for c in chosen] == [2, 4]


def test_agreement_lists_only_the_shanks_that_differ():
    def choice(shank, channel):
        return rc.ShankChoice(shank, channel, 1.0, 0.0, 1.0, 1.0, None, 1.0)

    first = [choice(1, 10), choice(2, 20)]
    second = [choice(1, 10), choice(2, 21)]
    assert rc.agreement(first, second) == {2: (20, 21)}
