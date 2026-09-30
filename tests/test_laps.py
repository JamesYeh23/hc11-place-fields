"""Lap segmentation: direction, hysteresis, coverage, and the circular wrap."""

from __future__ import annotations

import numpy as np
import pytest

from conftest import SYNTH_NAME, requires_data, write_sessinfo
from hc11.io import load_session
from hc11.laps import lap_mask, lap_summary, segment_laps
from hc11.track import Track

KW = dict(
    position_smoothing_sigma_s=0.25,
    direction_threshold_cm_s=5.0,
    min_coverage_frac=0.5,
    min_lap_duration_s=0.5,
    max_lap_duration_s=60.0,
)


def _session_from_positions(tmp_path, arrays, x, dt=0.0256, maze_start=100.0):
    """Build a session whose 1-D position follows ``x`` (metres); NaN allowed."""
    n = len(x)
    t = maze_start + np.arange(n) * dt
    arrays = dict(arrays)
    arrays.update(
        TimeStamps=t[None, :],
        TwoDLocation=np.column_stack([x, np.zeros_like(x)]),
        OneDLocation=np.asarray(x)[:, None],
        MazeEpoch=np.array([[maze_start, maze_start + n * dt]]),
        PREEpoch=np.array([[0.0, maze_start]]),
        POSTEpoch=np.array([[maze_start + n * dt, maze_start + n * dt + 100.0]]),
        sessDuration=np.array([[maze_start + n * dt + 100.0]]),
        Wake=np.array([[0.0, maze_start + n * dt + 100.0]]),
        Drowsy=np.array([[0.0, 1.0]]), NREM=np.array([[1.0, 2.0]]),
        Intermediate=np.array([[2.0, 3.0]]), REM=np.array([[3.0, 4.0]]),
    )
    return load_session(write_sessinfo(tmp_path / f"{SYNTH_NAME}_sessInfo.mat", arrays))


def _triangle(n_laps, samples_per_lap, length=1.6):
    """Back-and-forth traversals of a linear track."""
    up = np.linspace(0, length, samples_per_lap)
    down = up[::-1]
    return np.concatenate([up if i % 2 == 0 else down for i in range(n_laps)])


class TestLinear:
    def test_counts_alternating_traversals(self, tmp_path, synthetic_arrays):
        x = _triangle(6, 120)
        s = _session_from_positions(tmp_path, synthetic_arrays, x)
        track = Track(kind="linear", lo=0.0, hi=1.6, bin_size=0.1, nominal_length=1.6)
        laps = segment_laps(s, track, **KW)
        summary = lap_summary(laps)
        assert summary["pos"]["n_laps"] == 3
        assert summary["neg"]["n_laps"] == 3
        assert all(lap.coverage > 0.9 for lap in laps)

    def test_does_not_assume_alternation(self, tmp_path, synthetic_arrays):
        """Two runs in the same direction in a row are two laps of that direction."""
        up = np.linspace(0, 1.6, 120)
        teleport = np.full(10, np.nan)  # tracking gap; no direction inferred across it
        x = np.concatenate([up, teleport, up, teleport, up[::-1]])
        s = _session_from_positions(tmp_path, synthetic_arrays, x)
        track = Track(kind="linear", lo=0.0, hi=1.6, bin_size=0.1)
        summary = lap_summary(segment_laps(s, track, **KW))
        assert summary["pos"]["n_laps"] == 2
        assert summary["neg"]["n_laps"] == 1

    def test_partial_traversal_is_rejected(self, tmp_path, synthetic_arrays):
        x = np.linspace(0.0, 0.4, 120)  # only 25 % of the track
        s = _session_from_positions(tmp_path, synthetic_arrays, x)
        track = Track(kind="linear", lo=0.0, hi=1.6, bin_size=0.1)
        assert segment_laps(s, track, **KW) == []

    def test_brief_slowdown_does_not_split_a_lap(self, tmp_path, synthetic_arrays):
        """The reason laps are segmented before the speed filter (decision D4.1)."""
        first = np.linspace(0, 0.7, 60)
        pause = np.full(20, 0.7)  # ~0.5 s stationary mid-traversal
        second = np.linspace(0.7, 1.6, 60)
        s = _session_from_positions(
            tmp_path, synthetic_arrays, np.concatenate([first, pause, second])
        )
        track = Track(kind="linear", lo=0.0, hi=1.6, bin_size=0.1)
        laps = segment_laps(s, track, **KW)
        assert len(laps) == 1, "a mid-traversal pause must not create two laps"
        assert laps[0].coverage > 0.9

    def test_wobble_does_not_flip_direction(self, tmp_path, synthetic_arrays):
        x = np.linspace(0, 1.6, 200)
        x[90:100] -= 0.01  # 1 cm backward jitter
        s = _session_from_positions(tmp_path, synthetic_arrays, x)
        track = Track(kind="linear", lo=0.0, hi=1.6, bin_size=0.1)
        laps = segment_laps(s, track, **KW)
        assert len(laps) == 1 and laps[0].direction == "pos"

    def test_too_long_is_rejected(self, tmp_path, synthetic_arrays):
        x = np.linspace(0, 1.6, 4000)  # ~102 s, over max_lap_duration_s
        s = _session_from_positions(tmp_path, synthetic_arrays, x)
        track = Track(kind="linear", lo=0.0, hi=1.6, bin_size=0.1)
        assert segment_laps(s, track, **KW) == []

    def test_lap_mask_selects_only_that_direction(self, tmp_path, synthetic_arrays):
        s = _session_from_positions(tmp_path, synthetic_arrays, _triangle(4, 120))
        track = Track(kind="linear", lo=0.0, hi=1.6, bin_size=0.1)
        laps = segment_laps(s, track, **KW)
        m_pos = lap_mask(s, laps, "pos")
        m_neg = lap_mask(s, laps, "neg")
        assert not np.any(m_pos & m_neg)
        assert np.all(lap_mask(s, laps) == (m_pos | m_neg))


class TestCircular:
    """One lap is one circuit, measured between crossings of the reward site at 0."""

    def test_one_full_circuit_from_the_reward_site_is_one_lap(self, tmp_path, synthetic_arrays):
        track = Track(kind="circular", lo=0.0, hi=2.9, bin_size=0.1)
        x = track.wrap(np.linspace(0.0, 2.9, 200, endpoint=False))
        s = _session_from_positions(tmp_path, synthetic_arrays, x)
        laps = segment_laps(s, track, **KW)
        assert len(laps) == 1
        assert laps[0].direction == "pos"
        assert laps[0].coverage == pytest.approx(1.0, abs=0.05)

    def test_continuous_running_splits_into_one_lap_per_circuit(self, tmp_path, synthetic_arrays):
        """Without the split at the reward site this would be a single 2-circuit 'lap'."""
        track = Track(kind="circular", lo=0.0, hi=2.9, bin_size=0.1)
        x = track.wrap(np.linspace(1.0, 1.0 + 5.8, 400))  # 2 circuits, starting mid-ring
        s = _session_from_positions(tmp_path, synthetic_arrays, x)
        laps = segment_laps(s, track, **KW)
        assert {lap.direction for lap in laps} == {"pos"}, "no spurious reversal at the wrap"
        assert all(lap.coverage <= 1.02 for lap in laps), "a lap must not span several circuits"
        assert max(lap.coverage for lap in laps) == pytest.approx(1.0, abs=0.05)
        # 1.0 -> 2.9 (0.66 of a circuit), 0 -> 2.9 (a full one), 0 -> 1.0 (0.34, too short)
        assert len(laps) == 2

    def test_wrap_is_not_read_as_a_direction_reversal(self, tmp_path, synthetic_arrays):
        track = Track(kind="circular", lo=0.0, hi=2.9, bin_size=0.1)
        x = track.wrap(np.linspace(2.5, 2.5 + 5.8, 400))
        s = _session_from_positions(tmp_path, synthetic_arrays, x)
        assert {lap.direction for lap in segment_laps(s, track, **KW)} == {"pos"}


@requires_data
def test_real_linear_session_laps(real_sessions):
    from hc11 import config

    s = real_sessions("Achilles_10252013")
    track = Track.from_session(s, 10.0)
    laps = segment_laps(s, track, **KW)
    summary = lap_summary(laps)
    assert summary["pos"]["n_laps"] > 20 and summary["neg"]["n_laps"] > 20
    assert summary["pos"]["median_coverage"] > 0.8
    for lap in laps:
        assert lap.t_end > lap.t_start
        assert config.SESSIONS  # keeps the import meaningful
