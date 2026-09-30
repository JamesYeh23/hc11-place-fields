"""Rate maps, Skaggs information, field detection, shuffle test, recovery."""

from __future__ import annotations

import numpy as np
import pytest

from conftest import SYNTH_NAME, requires_data, write_sessinfo
from hc11.fields import (
    RunningSamples,
    build_rate_map,
    find_fields,
    rate_map,
    select_running_samples,
    shuffle_test,
    skaggs_information,
    smooth_1d,
    spike_slots,
    split_half_stability,
)
from hc11.io import load_session
from hc11.laps import segment_laps
from hc11.track import Track

LINEAR = Track(kind="linear", lo=0.0, hi=1.6, bin_size=0.1, nominal_length=1.6)
CIRCLE = Track(kind="circular", lo=0.0, hi=2.9, bin_size=0.1)


class TestSkaggs:
    """Hand-computed cases. I = sum_i p_i (r_i/rbar) log2(r_i/rbar)."""

    def test_uniform_rate_carries_no_information(self):
        rate = np.full(8, 3.0)
        occ = np.full(8, 2.0)
        bits_spike, bits_sec = skaggs_information(rate, occ)
        assert bits_spike == pytest.approx(0.0, abs=1e-12)
        assert bits_sec == pytest.approx(0.0, abs=1e-12)

    def test_half_silent_half_active(self):
        """p = [0.5, 0.5], rates [0, 2]: rbar = 1, I = 0.5 * 2 * log2(2) = 1 bit/spike."""
        bits_spike, bits_sec = skaggs_information(np.array([0.0, 2.0]), np.array([1.0, 1.0]))
        assert bits_spike == pytest.approx(1.0)
        assert bits_sec == pytest.approx(1.0)  # rbar = 1 Hz, so bits/s == bits/spike

    def test_mild_modulation(self):
        """p = [0.5, 0.5], rates [1, 3]: rbar = 2.

        I = 0.5*0.5*log2(0.5) + 0.5*1.5*log2(1.5) = -0.25 + 0.4387 = 0.18872
        """
        bits_spike, bits_sec = skaggs_information(np.array([1.0, 3.0]), np.array([1.0, 1.0]))
        assert bits_spike == pytest.approx(0.188722, abs=1e-6)
        assert bits_sec == pytest.approx(2.0 * 0.188722, abs=1e-6)

    def test_uneven_occupancy_weights_the_bins(self):
        """p = [0.8, 0.2], rates [1, 5]: rbar = 1.8."""
        rate, occ = np.array([1.0, 5.0]), np.array([8.0, 2.0])
        p = occ / occ.sum()
        expected = float(np.sum(p * (rate / 1.8) * np.log2(rate / 1.8)))
        assert skaggs_information(rate, occ)[0] == pytest.approx(expected)

    def test_a_single_active_bin_of_n_gives_log2_n(self):
        """One of 8 equally visited bins active: I = log2(8) = 3 bits/spike."""
        rate = np.zeros(8)
        rate[3] = 4.0
        assert skaggs_information(rate, np.ones(8))[0] == pytest.approx(3.0)

    def test_nan_bins_are_excluded_and_p_renormalised(self):
        rate = np.array([0.0, 2.0, np.nan])
        occ = np.array([1.0, 1.0, 5.0])
        assert skaggs_information(rate, occ)[0] == pytest.approx(1.0)

    def test_silent_cell_is_nan(self):
        assert np.isnan(skaggs_information(np.zeros(5), np.ones(5))[0])

    def test_nothing_visited_is_nan(self):
        assert np.isnan(skaggs_information(np.full(5, np.nan), np.zeros(5))[0])


class TestSmoothing:
    def test_circular_smoothing_wraps(self):
        a = np.zeros(20)
        a[0] = 1.0
        out = smooth_1d(a, 1.0, circular=True)
        assert out[-1] == pytest.approx(out[1]), "mass must cross the wrap point"

    def test_linear_smoothing_does_not_wrap(self):
        a = np.zeros(20)
        a[0] = 1.0
        out = smooth_1d(a, 1.0, circular=False)
        assert out[-1] < 1e-6

    def test_sigma_zero_is_identity(self):
        a = np.arange(5.0)
        np.testing.assert_array_equal(smooth_1d(a, 0.0, circular=False), a)

    def test_ratio_is_unbiased_at_a_linear_edge(self):
        """Smoothing counts and occupancy separately keeps the edge rate correct."""
        occ = np.full(16, 2.0)
        counts = np.full(16, 6.0)  # a uniform 3 Hz cell
        r = rate_map(counts, occ, sigma_bins=1.0, circular=False, min_occupancy_s=0.1)
        np.testing.assert_allclose(r, 3.0, rtol=1e-6)

    def test_under_occupied_bins_are_nan(self):
        occ = np.full(8, 2.0)
        occ[3] = 0.01
        r = rate_map(np.ones(8), occ, sigma_bins=0.0, circular=False, min_occupancy_s=0.1)
        assert np.isnan(r[3]) and np.isfinite(r[0])


class TestFindFields:
    def test_single_field(self):
        rate = np.zeros(16)
        rate[6:9] = [2.0, 5.0, 2.0]
        fields = find_fields(rate, 0.1, threshold_frac=0.2, min_width_m=0.1,
                             max_width_m=1.2, circular=False)
        assert fields == [(6, 3)]

    def test_two_fields(self):
        rate = np.zeros(16)
        rate[2:4] = 5.0
        rate[10:12] = 4.0
        fields = find_fields(rate, 0.1, threshold_frac=0.2, min_width_m=0.1,
                             max_width_m=1.2, circular=False)
        assert len(fields) == 2

    def test_too_wide_is_rejected(self):
        rate = np.linspace(10.0, 9.0, 16)  # broad, nearly flat: 160 cm above threshold
        assert find_fields(rate, 0.1, threshold_frac=0.2, min_width_m=0.1,
                           max_width_m=1.2, circular=False) == []

    def test_field_wrapping_the_circular_boundary_is_one_field(self):
        rate = np.zeros(29)
        rate[-2:] = 4.0
        rate[:2] = 4.0
        wrapped = find_fields(rate, 0.1, threshold_frac=0.2, min_width_m=0.1,
                              max_width_m=1.2, circular=True)
        straight = find_fields(rate, 0.1, threshold_frac=0.2, min_width_m=0.1,
                               max_width_m=1.2, circular=False)
        assert len(wrapped) == 1 and wrapped[0][1] == 4
        assert len(straight) == 2, "without wrapping it would look like two fields"

    def test_silent_cell_has_no_fields(self):
        assert find_fields(np.zeros(16), 0.1, threshold_frac=0.2, min_width_m=0.1,
                           max_width_m=1.2, circular=False) == []


# ---------------------------------------------------------------------------
# End-to-end recovery from a synthetic Poisson place cell
# ---------------------------------------------------------------------------


def _place_cell_session(tmp_path, arrays, field_center_m, peak_hz=10.0, width_m=0.2,
                        n_laps=40, seed=0, baseline_hz=0.1, jitter=True):
    """A session with one Poisson place cell whose field we know exactly.

    Lap durations are jittered by default, as real running is. With
    ``jitter=False`` every lap takes exactly the same time, which makes the
    circular-shift null degenerate -- see
    ``test_stereotyped_running_weakens_the_shuffle``.
    """
    rng = np.random.default_rng(seed)
    dt = 0.0256
    lengths = (
        rng.integers(90, 150, n_laps) if jitter else np.full(n_laps, 120)
    )
    x = np.concatenate([
        np.linspace(0, 1.6, n) if i % 2 == 0 else np.linspace(1.6, 0, n)
        for i, n in enumerate(lengths)
    ])
    n = len(x)
    t = 100.0 + np.arange(n) * dt

    rate = baseline_hz + peak_hz * np.exp(-0.5 * ((x - field_center_m) / width_m) ** 2)
    counts = rng.poisson(rate * dt)
    spike_times = np.concatenate([
        t[i] + rng.uniform(0, dt, c) for i, c in enumerate(counts) if c
    ])
    spike_times = np.sort(spike_times)

    arrays = dict(arrays)
    arrays.update(
        SpikeTimes=spike_times[:, None],
        SpikeIDs=np.full((len(spike_times), 1), 102.0),
        PyrIDs=np.array([[102.0]]),
        IntIDs=np.array([[104.0]]),
        TimeStamps=t[None, :],
        TwoDLocation=np.column_stack([x, np.zeros_like(x)]),
        OneDLocation=x[:, None],
        MazeEpoch=np.array([[100.0, 100.0 + n * dt]]),
        PREEpoch=np.array([[0.0, 100.0]]),
        POSTEpoch=np.array([[100.0 + n * dt, 100.0 + n * dt + 50.0]]),
        sessDuration=np.array([[100.0 + n * dt + 50.0]]),
        Wake=np.array([[0.0, 100.0 + n * dt + 50.0]]),
        Drowsy=np.array([[0.0, 1.0]]), NREM=np.array([[1.0, 2.0]]),
        Intermediate=np.array([[2.0, 3.0]]), REM=np.array([[3.0, 4.0]]),
    )
    # IntIDs must exist in SpikeIDs for the loader's bookkeeping; give it one spike.
    arrays["SpikeTimes"] = np.vstack([arrays["SpikeTimes"], [[100.0 + n * dt - dt]]])
    arrays["SpikeIDs"] = np.vstack([arrays["SpikeIDs"], [[104.0]]])
    order = np.argsort(arrays["SpikeTimes"].ravel())
    arrays["SpikeTimes"] = arrays["SpikeTimes"][order]
    arrays["SpikeIDs"] = arrays["SpikeIDs"][order]
    return load_session(write_sessinfo(tmp_path / f"{SYNTH_NAME}_sessInfo.mat", arrays))


def _samples_for(session, track, direction="pos"):
    laps = segment_laps(
        session, track, position_smoothing_sigma_s=0.25, direction_threshold_cm_s=5.0,
        min_coverage_frac=0.5, min_lap_duration_s=0.5, max_lap_duration_s=60.0,
    )
    mask = np.ones(len(session.position_t), dtype=bool)
    return select_running_samples(session, track, laps, mask, direction)


class TestRecovery:
    @pytest.mark.parametrize("center", [0.35, 0.8, 1.25])
    def test_field_position_is_recovered(self, tmp_path, synthetic_arrays, center):
        s = _place_cell_session(tmp_path, synthetic_arrays, center)
        samples = _samples_for(s, LINEAR)
        rmap = build_rate_map(s, LINEAR, samples, 102, sigma_bins=1.0, min_occupancy_s=0.1)
        assert rmap.peak_position == pytest.approx(center, abs=0.12)
        assert rmap.peak_rate == pytest.approx(10.0, rel=0.4)

    def _info_and_p(self, session, n_shuffles=300, seed=1):
        samples = _samples_for(session, LINEAR)
        slots = spike_slots(session, 102, samples)
        rmap = build_rate_map(session, LINEAR, samples, 102, sigma_bins=1.0,
                              min_occupancy_s=0.1, slots=slots)
        info, _ = skaggs_information(rmap.rate, rmap.occupancy_s)
        null, _ = shuffle_test(samples, slots, sigma_bins=1.0, min_occupancy_s=0.1,
                               circular_track=False, n_shuffles=n_shuffles,
                               rng=np.random.default_rng(seed))
        p = (np.sum(null >= info) + 1) / (np.isfinite(null).sum() + 1)
        return info, float(np.nanmean(null)), float(p)

    def test_a_real_place_cell_is_significant(self, tmp_path, synthetic_arrays):
        s = _place_cell_session(tmp_path, synthetic_arrays, 0.8)
        info, null_mean, p = self._info_and_p(s)
        assert info > 0.5
        # The contrast against the null is what matters, not the absolute value:
        # smoothing and the Gaussian tails both bleed information out of a field.
        assert info > 3 * null_mean
        assert p < 0.05, f"a Poisson place cell should be detected (p={p})"

    def test_uniform_cells_are_called_significant_at_about_the_nominal_rate(
        self, tmp_path, synthetic_arrays
    ):
        """Calibration, not a single-case check.

        A spatially uniform cell has, by construction, a 5 % chance of beating the
        null at alpha = 0.05, so asserting that one particular uniform cell is
        non-significant would be a coin flip. Ten independent uniform cells should
        yield roughly zero or one false positive; three or more would mean the null
        is mis-calibrated and the place-cell counts are inflated.
        """
        flagged = 0
        for seed in range(10):
            s = _place_cell_session(
                tmp_path, synthetic_arrays, 0.8, peak_hz=0.0, baseline_hz=5.0, seed=seed
            )
            _, _, p = self._info_and_p(s, n_shuffles=200, seed=100 + seed)
            flagged += p < 0.05
        assert flagged <= 3, f"{flagged}/10 uniform cells called significant -- null too tight"

    def test_stereotyped_running_weakens_the_shuffle(self, tmp_path, synthetic_arrays):
        """A known limitation of the circular-shift null, pinned deliberately.

        If every lap takes exactly the same number of samples, shifting the spike
        train by k slots maps bin i onto bin i+k in *every* lap, so the shuffled
        rate map is a rotation of the true one and carries the same spatial
        information. The null then sits at the observed value and the test loses
        its power -- conservatively, so no cell is wrongly called a place cell.
        Real lap durations vary (Achilles_10252013: 2.0-6.9 s), which is what
        gives the shuffle its power on the actual data.
        """
        s = _place_cell_session(tmp_path, synthetic_arrays, 0.8, jitter=False)
        samples = _samples_for(s, LINEAR)
        slots = spike_slots(s, 102, samples)
        rmap = build_rate_map(s, LINEAR, samples, 102, sigma_bins=1.0,
                              min_occupancy_s=0.1, slots=slots)
        info, _ = skaggs_information(rmap.rate, rmap.occupancy_s)
        null, _ = shuffle_test(samples, slots, sigma_bins=1.0, min_occupancy_s=0.1,
                               circular_track=False, n_shuffles=200,
                               rng=np.random.default_rng(3))
        assert np.nanmean(null) > 0.5 * info, "documented degeneracy, not a passing test"

    def test_shuffle_preserves_the_spike_count(self, tmp_path, synthetic_arrays):
        s = _place_cell_session(tmp_path, synthetic_arrays, 0.8)
        samples = _samples_for(s, LINEAR)
        slots = spike_slots(s, 102, samples)
        per_slot = np.bincount(slots, minlength=samples.n)
        assert np.roll(per_slot, 137).sum() == per_slot.sum() == slots.size

    def test_split_half_stability_is_high_for_a_stable_cell(self, tmp_path, synthetic_arrays):
        s = _place_cell_session(tmp_path, synthetic_arrays, 0.8)
        samples = _samples_for(s, LINEAR)
        slots = spike_slots(s, 102, samples)
        r = split_half_stability(s, LINEAR, samples, 102, slots,
                                 sigma_bins=1.0, min_occupancy_s=0.1)
        assert r > 0.8

    def test_spikes_outside_the_running_mask_are_ignored(self, tmp_path, synthetic_arrays):
        s = _place_cell_session(tmp_path, synthetic_arrays, 0.8)
        samples = _samples_for(s, LINEAR)
        all_slots = spike_slots(s, 102, samples)
        # Restricting the mask to half the samples must not add spikes.
        half = np.zeros(len(s.position_t), dtype=bool)
        half[: len(half) // 2] = True
        laps = segment_laps(
            s, LINEAR, position_smoothing_sigma_s=0.25, direction_threshold_cm_s=5.0,
            min_coverage_frac=0.5, min_lap_duration_s=0.5, max_lap_duration_s=60.0,
        )
        fewer = select_running_samples(s, LINEAR, laps, half, "pos")
        assert spike_slots(s, 102, fewer).size <= all_slots.size


class TestSampleSelection:
    def test_slots_are_sequential_and_unique(self, tmp_path, synthetic_arrays):
        s = _place_cell_session(tmp_path, synthetic_arrays, 0.8)
        samples = _samples_for(s, LINEAR)
        assert isinstance(samples, RunningSamples)
        used = samples.slot_of_index[samples.slot_of_index >= 0]
        np.testing.assert_array_equal(np.sort(used), np.arange(samples.n))

    def test_bins_are_in_range(self, tmp_path, synthetic_arrays):
        s = _place_cell_session(tmp_path, synthetic_arrays, 0.8)
        samples = _samples_for(s, LINEAR)
        assert samples.bins.min() >= 0 and samples.bins.max() < LINEAR.n_bins

    def test_directions_are_disjoint(self, tmp_path, synthetic_arrays):
        s = _place_cell_session(tmp_path, synthetic_arrays, 0.8)
        pos = _samples_for(s, LINEAR, "pos")
        neg = _samples_for(s, LINEAR, "neg")
        assert not set(pos.indices) & set(neg.indices)


@requires_data
def test_real_session_produces_the_expected_staircase(real_sessions):
    """Place-cell peaks should tile the track, not pile up in one place."""
    from hc11 import config
    from hc11.behavior import running_mask_from_config
    from hc11.fields import analyse_cell
    from hc11.laps import segment_laps_from_config

    params = config.load_params()
    s = real_sessions("Achilles_10252013")
    track = Track.from_session(s, params["place_fields"]["bin_size_cm"])
    running = running_mask_from_config(s, params)
    laps = segment_laps_from_config(s, track, params)
    samples = select_running_samples(s, track, laps, running.mask, "pos")
    rng = np.random.default_rng(0)
    cells = [analyse_cell(s, track, samples, c, params, rng) for c in s.pyr_ids]
    peaks = np.array([c.peak_position_m for c in cells if c.is_place_cell])
    assert peaks.size > 20
    occupied = np.unique(track.digitize(peaks))
    assert occupied.size >= track.n_bins // 2, "peaks should spread across the track"
