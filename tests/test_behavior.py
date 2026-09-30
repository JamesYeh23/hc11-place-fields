"""Speed estimation and the two-stage running mask."""

from __future__ import annotations

import numpy as np
import pytest

from conftest import SYNTH_NAME, requires_data, write_sessinfo
from hc11.behavior import running_mask, speed_cm_s
from hc11.io import load_session

PARAMS = dict(
    speed_threshold_cm_s=15.0,
    speed_smoothing_sigma_s=0.1,
    min_run_epoch_s=0.5,
    max_run_gap_s=0.2,
)


def _constant_speed_session(tmp_path, arrays, speed_cm_per_s, n=400):
    """A session whose animal moves at a known constant speed along x."""
    dt = 0.0256
    t = 100.0 + np.arange(n) * dt
    x = (speed_cm_per_s / 100.0) * np.arange(n) * dt
    arrays = dict(arrays)
    arrays["TimeStamps"] = t[None, :]
    arrays["TwoDLocation"] = np.column_stack([x, np.zeros_like(x)])
    arrays["OneDLocation"] = x[:, None]
    arrays["MazeEpoch"] = np.array([[100.0, 100.0 + n * dt]])
    arrays["PREEpoch"] = np.array([[0.0, 100.0]])
    arrays["POSTEpoch"] = np.array([[100.0 + n * dt, 300.0]])
    return load_session(write_sessinfo(tmp_path / f"{SYNTH_NAME}_sessInfo.mat", arrays))


class TestSpeed:
    @pytest.mark.parametrize("v", [5.0, 15.0, 60.0])
    def test_recovers_a_known_constant_speed(self, tmp_path, synthetic_arrays, v):
        s = _constant_speed_session(tmp_path, synthetic_arrays, v)
        out = speed_cm_s(s, smooth_sigma_s=0.1)
        np.testing.assert_allclose(np.median(out), v, rtol=1e-6)

    def test_nan_position_gives_nan_speed_and_does_not_spread(self, tmp_path, synthetic_arrays):
        s = _constant_speed_session(tmp_path, synthetic_arrays, 30.0)
        arrays = dict(synthetic_arrays)
        xy = s.position_xy.copy()
        xy[200] = np.nan  # a one-sample tracking dropout
        arrays.update(
            TimeStamps=s.position_t[None, :],
            TwoDLocation=xy,
            OneDLocation=s.position_1d[:, None],
            MazeEpoch=s.epoch("MAZE")[None, :],
            PREEpoch=s.epoch("PRE")[None, :],
            POSTEpoch=s.epoch("POST")[None, :],
        )
        s = load_session(write_sessinfo(tmp_path / f"{SYNTH_NAME}_sessInfo.mat", arrays))
        out = speed_cm_s(s, smooth_sigma_s=0.0)
        assert np.isnan(out[199:202]).all()  # the gradient touches its neighbours
        assert np.isfinite(out[:198]).all() and np.isfinite(out[203:]).all()

    def test_smoothing_does_not_shift_the_median(self, tmp_path, synthetic_arrays):
        s = _constant_speed_session(tmp_path, synthetic_arrays, 40.0)
        raw = np.nanmedian(speed_cm_s(s, 0.0))
        smoothed = np.nanmedian(speed_cm_s(s, 0.2))
        assert smoothed == pytest.approx(raw, rel=1e-6)


class TestRunningMask:
    def test_fast_session_is_kept_entirely(self, tmp_path, synthetic_arrays):
        s = _constant_speed_session(tmp_path, synthetic_arrays, 60.0)
        rm = running_mask(s, **PARAMS)
        assert rm.n_final == rm.n_authors
        assert rm.n_removed_by_speed == 0
        assert rm.fraction_kept == 1.0

    def test_slow_session_is_rejected_entirely(self, tmp_path, synthetic_arrays):
        s = _constant_speed_session(tmp_path, synthetic_arrays, 5.0)
        rm = running_mask(s, **PARAMS)
        assert rm.n_final == 0
        assert rm.fraction_removed_by_speed == 1.0
        assert rm.intervals.shape == (0, 2)

    def test_mask_never_exceeds_the_authors_mask(self, synthetic_file):
        s = load_session(synthetic_file)
        rm = running_mask(s, **PARAMS)
        assert not np.any(rm.mask & ~rm.authors_mask)

    def test_bridging_only_fills_short_enclosed_gaps(self, tmp_path, synthetic_arrays):
        dt, n = 0.0256, 400
        t = 100.0 + np.arange(n) * dt
        v = np.full(n, 0.60)  # m/s
        v[100:104] = 0.0  # 0.10 s dip -- shorter than max_run_gap_s
        v[200:250] = 0.0  # 1.28 s stop -- longer
        x = np.cumsum(v) * dt
        arrays = dict(synthetic_arrays)
        arrays["TimeStamps"] = t[None, :]
        arrays["TwoDLocation"] = np.column_stack([x, np.zeros_like(x)])
        arrays["OneDLocation"] = x[:, None]
        arrays["MazeEpoch"] = np.array([[100.0, 100.0 + n * dt]])
        arrays["PREEpoch"] = np.array([[0.0, 100.0]])
        arrays["POSTEpoch"] = np.array([[100.0 + n * dt, 300.0]])
        s = load_session(write_sessinfo(tmp_path / f"{SYNTH_NAME}_sessInfo.mat", arrays))
        rm = running_mask(s, **{**PARAMS, "speed_smoothing_sigma_s": 0.0})
        assert rm.n_restored_by_bridging > 0
        assert rm.mask[101], "a 0.1 s dip should be bridged"
        assert not rm.mask[225], "a 1.3 s stop should not be bridged"

    def test_short_epochs_are_dropped(self, tmp_path, synthetic_arrays):
        dt, n = 0.0256, 400
        t = 100.0 + np.arange(n) * dt
        v = np.zeros(n)
        v[50:55] = 0.60  # 0.13 s of running -- below min_run_epoch_s
        x = np.cumsum(v) * dt
        arrays = dict(synthetic_arrays)
        arrays["TimeStamps"] = t[None, :]
        arrays["TwoDLocation"] = np.column_stack([x, np.zeros_like(x)])
        arrays["OneDLocation"] = x[:, None]
        arrays["MazeEpoch"] = np.array([[100.0, 100.0 + n * dt]])
        arrays["PREEpoch"] = np.array([[0.0, 100.0]])
        arrays["POSTEpoch"] = np.array([[100.0 + n * dt, 300.0]])
        s = load_session(write_sessinfo(tmp_path / f"{SYNTH_NAME}_sessInfo.mat", arrays))
        rm = running_mask(s, **{**PARAMS, "speed_smoothing_sigma_s": 0.0})
        assert rm.n_removed_short_epochs > 0 and rm.n_final == 0

    def test_intervals_match_the_mask_duration(self, synthetic_file):
        s = load_session(synthetic_file)
        rm = running_mask(s, **PARAMS)
        assert np.sum(rm.intervals[:, 1] - rm.intervals[:, 0]) == pytest.approx(rm.seconds)


@requires_data
@pytest.mark.parametrize("name", ["Achilles_10252013", "Cicero_09012014"])
def test_running_mask_on_real_sessions(real_sessions, name):
    s = real_sessions(name)
    rm = running_mask(s, **PARAMS)
    assert 0 < rm.n_final <= rm.n_authors
    assert not np.any(rm.mask & np.isnan(s.position_1d))
    assert np.nanmedian(rm.speed[rm.mask]) >= 15.0
