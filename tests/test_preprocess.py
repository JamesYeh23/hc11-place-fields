"""Truncation, interval cleaning, and mask/interval conversion."""

from __future__ import annotations

import numpy as np
import pytest

from conftest import SYNTH_NAME, requires_data, write_sessinfo
from hc11.io import load_session
from hc11.preprocess import (
    clean_intervals,
    intervals_to_mask,
    mask_to_intervals,
    state_intervals,
    truncate_to_session,
)


class TestTruncate:
    def test_clean_session_is_returned_unchanged(self, synthetic_file):
        s = load_session(synthetic_file)
        out, report = truncate_to_session(s)
        assert out is s, "an unaffected session should not be copied"
        assert not report.changed and report.n_spikes_dropped == 0

    def test_drops_spikes_after_end(self, tmp_path, synthetic_arrays):
        arrays = dict(synthetic_arrays)
        arrays["SpikeTimes"] = np.append(arrays["SpikeTimes"], [[305.0], [310.0]], axis=0)
        arrays["SpikeIDs"] = np.append(arrays["SpikeIDs"], [[102.0], [205.0]], axis=0)
        s = load_session(write_sessinfo(tmp_path / f"{SYNTH_NAME}_sessInfo.mat", arrays))
        out, report = truncate_to_session(s)
        assert report.n_spikes_dropped == 2
        assert report.clusters_affected == 2
        assert report.seconds_beyond_end == pytest.approx(10.0)
        assert out.spike_times.max() <= out.sess_duration
        assert out.spike_times.size == s.spike_times.size - 2
        # The lookup index must be rebuilt, not inherited from the original.
        assert out.spikes_for(102).size == np.sum(out.spike_ids == 102)

    def test_clips_and_drops_state_intervals(self, tmp_path, synthetic_arrays):
        arrays = dict(synthetic_arrays)
        arrays["NREM"] = np.array([[30.0, 80.0], [290.0, 320.0], [400.0, 450.0]])
        s = load_session(write_sessinfo(tmp_path / f"{SYNTH_NAME}_sessInfo.mat", arrays))
        out, report = truncate_to_session(s)
        assert report.state_intervals_dropped["NREM"] == 1  # starts after the end
        assert report.state_intervals_clipped["NREM"] == 1  # 290-320 -> 290-300
        np.testing.assert_allclose(out.states["NREM"], [[30.0, 80.0], [290.0, 300.0]])
        assert report.state_seconds_removed == pytest.approx(20.0 + 50.0)

    def test_original_session_is_not_mutated(self, tmp_path, synthetic_arrays):
        arrays = dict(synthetic_arrays)
        arrays["SpikeTimes"] = np.append(arrays["SpikeTimes"], [[305.0]], axis=0)
        arrays["SpikeIDs"] = np.append(arrays["SpikeIDs"], [[102.0]], axis=0)
        s = load_session(write_sessinfo(tmp_path / f"{SYNTH_NAME}_sessInfo.mat", arrays))
        before = s.spike_times.size
        truncate_to_session(s)
        assert s.spike_times.size == before


class TestCleanIntervals:
    def test_drops_zero_and_negative_length(self):
        iv = np.array([[0.0, 1.0], [2.0, 2.0], [5.0, 4.0], [6.0, 7.0]])
        np.testing.assert_allclose(clean_intervals(iv), [[0.0, 1.0], [6.0, 7.0]])

    def test_empty_stays_0x2(self):
        assert clean_intervals(np.empty((0, 2))).shape == (0, 2)

    def test_warns_once_with_the_offending_rows(self, caplog):
        clean_intervals(np.array([[1.0, 1.0]]), label="REM")
        assert "REM" in caplog.text and "[[1.0, 1.0]]" in caplog.text

    def test_state_intervals_restricted_to_epoch(self, synthetic_file):
        s = load_session(synthetic_file)
        iv = state_intervals(s, "NREM", restrict_to=np.array([50.0, 200.0]))
        assert iv.min() >= 50.0 and iv.max() <= 200.0

    def test_unknown_state_raises(self, synthetic_file):
        with pytest.raises(KeyError):
            state_intervals(load_session(synthetic_file), "Dreaming")


class TestMasksAndIntervals:
    def test_round_trip(self):
        t = np.arange(0, 1.0, 0.1)
        mask = np.array([0, 1, 1, 0, 0, 1, 0, 0, 1, 1], dtype=bool)
        iv = mask_to_intervals(t, mask, dt=0.1)
        assert iv.shape == (3, 2)
        np.testing.assert_allclose(intervals_to_mask(t, iv), mask)

    def test_duration_counts_every_sample(self):
        t = np.arange(0, 1.0, 0.1)
        mask = np.array([0, 1, 1, 0, 0, 1, 0, 0, 0, 0], dtype=bool)
        iv = mask_to_intervals(t, mask, dt=0.1)
        assert np.sum(iv[:, 1] - iv[:, 0]) == pytest.approx(mask.sum() * 0.1)

    def test_run_touching_the_end(self):
        t = np.arange(5) * 0.1
        iv = mask_to_intervals(t, np.array([0, 0, 1, 1, 1], dtype=bool), dt=0.1)
        np.testing.assert_allclose(iv, [[0.2, 0.5]])

    def test_empty_mask(self):
        assert mask_to_intervals(np.arange(5) * 0.1, np.zeros(5, bool)).shape == (0, 2)

    def test_length_mismatch_raises(self):
        with pytest.raises(ValueError, match="samples"):
            mask_to_intervals(np.arange(5), np.zeros(4, bool))


@requires_data
def test_gatsby_truncation_matches_step2_numbers(real_sessions):
    """The one real session affected, with the counts reported in step 2."""
    s = real_sessions("Gatsby_08022013")
    out, report = truncate_to_session(s)
    assert report.n_spikes_dropped == 263_976
    assert report.fraction_spikes_dropped == pytest.approx(0.0467, abs=5e-4)
    assert report.clusters_affected == 80
    assert report.seconds_beyond_end == pytest.approx(1588.7, abs=0.1)
    assert out.spike_times.max() <= out.sess_duration
    for iv in out.states.values():
        assert iv.size == 0 or iv.max() <= out.sess_duration


@requires_data
def test_other_sessions_are_untouched_by_truncation(real_sessions):
    from hc11 import config

    for name in config.SESSIONS:
        if name == "Gatsby_08022013":
            continue
        _, report = truncate_to_session(real_sessions(name))
        assert not report.changed, f"{name} unexpectedly needed truncation"


@requires_data
def test_only_achilles_11012013_has_zero_length_intervals(real_sessions):
    from hc11 import config
    from hc11.io import STATE_NAMES

    for name in config.SESSIONS:
        s = real_sessions(name)
        dropped = sum(
            len(s.states[st]) - len(state_intervals(s, st, verbose=False)) for st in STATE_NAMES
        )
        assert dropped == (1 if name == "Achilles_11012013" else 0), name
