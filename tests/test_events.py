"""The shared two-threshold detector, participation, and event comparison."""

from __future__ import annotations

import numpy as np
import pytest

from hc11.events import EventList, compare_events, detect_two_threshold, participation


def _trace(bumps, n=2000, dt=0.001, width=20):
    """A flat trace with Gaussian bumps at given (centre_index, height)."""
    t = np.arange(n) * dt
    signal = np.zeros(n)
    for centre, height in bumps:
        signal += height * np.exp(-0.5 * ((np.arange(n) - centre) / width) ** 2)
    return t, signal


KW = dict(high=3.0, low=1.0, min_duration_s=0.03, max_duration_s=0.4, merge_gap_s=0.03)


class TestDetectTwoThreshold:
    def test_finds_one_bump(self):
        t, signal = _trace([(1000, 5.0)])
        start, peak, end = detect_two_threshold(signal, t, **KW)
        assert len(start) == 1
        assert peak[0] == pytest.approx(1.0, abs=0.005)
        assert start[0] < peak[0] < end[0]

    def test_bump_below_the_high_threshold_is_ignored(self):
        t, signal = _trace([(1000, 2.0)])  # peaks at 2.0, below high=3
        assert len(detect_two_threshold(signal, t, **KW)[0]) == 0

    def test_boundaries_extend_to_the_low_threshold(self):
        t, signal = _trace([(1000, 5.0)])
        start, _, end = detect_two_threshold(signal, t, **KW)
        inside = (t >= start[0]) & (t < end[0])
        assert signal[inside].min() >= 1.0 - 1e-9, "inside the event, never below `low`"
        outside_before = signal[t < start[0]]
        assert outside_before.max() < 1.0 + 1e-9

    def test_two_separate_bumps_stay_separate(self):
        t, signal = _trace([(500, 5.0), (1500, 5.0)])
        assert len(detect_two_threshold(signal, t, **KW)[0]) == 2

    def test_close_bumps_are_merged(self):
        # Each bump clears `low` for ~43 ms, and the 18 ms trough between them
        # is shorter than merge_gap_s, so they become one event.
        t, signal = _trace([(1000, 5.0), (1060, 5.0)], width=12)
        merged = detect_two_threshold(signal, t, **KW)[0]
        separate = detect_two_threshold(signal, t, **{**KW, "merge_gap_s": 0.0})[0]
        assert len(merged) == 1
        assert len(separate) == 2, "without merging they are two events"

    def test_too_short_is_dropped(self):
        t, signal = _trace([(1000, 5.0)], width=2)  # a few ms wide
        assert len(detect_two_threshold(signal, t, **KW)[0]) == 0

    def test_too_long_is_dropped(self):
        t, signal = _trace([(1000, 5.0)], width=400)
        assert len(detect_two_threshold(signal, t, **KW)[0]) == 0

    def test_empty_and_flat_inputs(self):
        assert len(detect_two_threshold(np.empty(0), np.empty(0), **KW)[0]) == 0
        t, signal = _trace([])
        assert len(detect_two_threshold(signal, t, **KW)[0]) == 0

    def test_shape_mismatch_raises(self):
        with pytest.raises(ValueError, match="samples"):
            detect_two_threshold(np.zeros(5), np.zeros(4), **KW)

    def test_event_at_the_very_end_is_kept(self):
        t, signal = _trace([(1960, 5.0)])
        start, _, end = detect_two_threshold(signal, t, **KW)
        assert len(start) == 1 and end[0] <= t[-1] + 0.002


def _events(starts, ends, peaks=None, detector="mua"):
    starts, ends = np.asarray(starts, float), np.asarray(ends, float)
    peaks = 0.5 * (starts + ends) if peaks is None else np.asarray(peaks, float)
    return EventList(
        session="S", epoch="PRE", detector=detector, start=starts, peak=peaks, end=ends,
        peak_z=np.full(len(starts), 4.0), baseline_mean=1.0, baseline_sd=1.0,
        state_seconds=600.0,
    )


class TestEventList:
    def test_rate_and_duration(self):
        events = _events([0, 10], [0.1, 10.08])
        assert len(events) == 2
        np.testing.assert_allclose(events.duration, [0.1, 0.08])
        assert events.rate_per_min == pytest.approx(0.2)

    def test_intervals_shape(self):
        assert _events([1], [2]).intervals.shape == (1, 2)
        assert _events([], []).intervals.shape == (0, 2)


class TestParticipation:
    def test_counts_distinct_cells_and_spikes(self, synthetic_file):
        from hc11.io import load_session

        session = load_session(synthetic_file)
        window = (100.0, 160.0)
        events = _events([window[0]], [window[1]])
        part = participation(session, events, session.pyr_ids)
        times, ids = session.spikes_in(window)
        expected_cells = len(np.unique(ids[np.isin(ids, session.pyr_ids)]))
        assert part.n_active[0] == expected_cells
        assert part.n_spikes[0] == np.sum(np.isin(ids, session.pyr_ids))
        assert part.n_cells == len(session.pyr_ids)
        assert 0 <= part.active_fraction[0] <= 1

    def test_empty_event_window(self, synthetic_file):
        from hc11.io import load_session

        session = load_session(synthetic_file)
        part = participation(session, _events([299.98], [299.99]), session.pyr_ids)
        assert part.n_active[0] >= 0


class TestCompareEvents:
    def test_identical_lists_match_completely(self):
        a = _events([1.0, 5.0], [1.1, 5.1])
        overlap = compare_events(a, a)
        assert overlap.fraction_a_matched == 1.0 and overlap.fraction_b_matched == 1.0
        np.testing.assert_allclose(overlap.offsets_s, 0.0)

    def test_disjoint_lists_do_not_match(self):
        a = _events([1.0], [1.1])
        b = _events([50.0], [50.1], detector="lfp")
        overlap = compare_events(a, b)
        assert overlap.n_a_matched == 0 and overlap.n_b_matched == 0

    def test_overlapping_windows_match_even_with_offset_peaks(self):
        a = _events([1.00], [1.20], peaks=[1.02])
        b = _events([1.15], [1.30], peaks=[1.25], detector="lfp")
        overlap = compare_events(a, b, match_window_s=0.01)
        assert overlap.n_a_matched == 1
        assert overlap.offsets_s[0] == pytest.approx(0.23)

    def test_near_peaks_match_without_overlap(self):
        a = _events([1.00], [1.05], peaks=[1.02])
        b = _events([1.06], [1.10], peaks=[1.06], detector="lfp")
        assert compare_events(a, b, match_window_s=0.05).n_a_matched == 1

    def test_one_b_can_match_several_a(self):
        a = _events([1.0, 1.2], [1.1, 1.3])
        b = _events([0.95], [1.35], detector="lfp")
        overlap = compare_events(a, b)
        assert overlap.n_a_matched == 2 and overlap.n_b_matched == 1

    def test_empty_lists(self):
        overlap = compare_events(_events([], []), _events([1.0], [1.1]))
        assert overlap.n_a_matched == 0 and np.isnan(overlap.fraction_a_matched)
