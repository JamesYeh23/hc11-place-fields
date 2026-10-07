"""Our own linearisation of 2-D position (diagnostic path only)."""

from __future__ import annotations

import numpy as np
import pytest

from conftest import SYNTH_NAME, requires_data, write_sessinfo
from hc11.io import load_session
from hc11.linearize import coverage, linearize_2d, with_own_linearization


def _session_with_xy(tmp_path, arrays, xy, maze="1.6m Linear Maze"):
    n = len(xy)
    dt = 0.0256
    t = 100.0 + np.arange(n) * dt
    arrays = dict(arrays)
    arrays.update(
        TimeStamps=t[None, :], TwoDLocation=np.asarray(xy, dtype=float),
        OneDLocation=np.full((n, 1), np.nan), MazeType=maze,
        MazeEpoch=np.array([[100.0, 100.0 + n * dt]]),
        PREEpoch=np.array([[0.0, 100.0]]),
        POSTEpoch=np.array([[100.0 + n * dt, 100.0 + n * dt + 50.0]]),
        sessDuration=np.array([[100.0 + n * dt + 50.0]]),
    )
    return load_session(write_sessinfo(tmp_path / f"{SYNTH_NAME}_sessInfo.mat", arrays))


class TestLinear:
    def test_recovers_distance_along_a_diagonal_track(self, tmp_path, synthetic_arrays):
        """A track at 45 degrees: the coordinate must be arc length, not x or y."""
        d = np.linspace(0, 1.6, 300)
        xy = np.column_stack([d / np.sqrt(2) + 3.0, d / np.sqrt(2) - 7.0])
        s = _session_with_xy(tmp_path, synthetic_arrays, xy)
        out = linearize_2d(s)
        np.testing.assert_allclose(out, d, atol=1e-9)

    def test_origin_is_at_the_low_end(self, tmp_path, synthetic_arrays):
        d = np.linspace(0, 1.6, 200)
        s = _session_with_xy(tmp_path, synthetic_arrays, np.column_stack([d, np.zeros_like(d)]))
        assert np.nanmin(linearize_2d(s)) == pytest.approx(0.0)

    def test_nan_tracking_stays_nan(self, tmp_path, synthetic_arrays):
        d = np.linspace(0, 1.6, 200)
        xy = np.column_stack([d, np.zeros_like(d)])
        xy[50:60] = np.nan
        s = _session_with_xy(tmp_path, synthetic_arrays, xy)
        out = linearize_2d(s)
        assert np.isnan(out[50:60]).all()
        assert np.isfinite(out[:50]).all() and np.isfinite(out[60:]).all()

    def test_off_axis_wobble_is_projected_out(self, tmp_path, synthetic_arrays):
        rng = np.random.default_rng(0)
        d = np.linspace(0, 1.6, 400)
        xy = np.column_stack([d, rng.normal(0, 0.01, d.size)])
        s = _session_with_xy(tmp_path, synthetic_arrays, xy)
        np.testing.assert_allclose(linearize_2d(s), d, atol=0.02)


class TestCircular:
    def test_arc_length_increases_around_the_ring(self, tmp_path, synthetic_arrays):
        angle = np.linspace(0.01, 2 * np.pi - 0.01, 300)
        r = 0.5
        xy = np.column_stack([r * np.cos(angle) + 1.0, r * np.sin(angle) - 2.0])
        s = _session_with_xy(tmp_path, synthetic_arrays, xy, maze="Circular Maze")
        out = linearize_2d(s)
        assert np.all(np.diff(out) > 0)
        assert out.max() == pytest.approx(r * (2 * np.pi - 0.01), rel=0.02)

    def test_unknown_maze_raises(self, tmp_path, synthetic_arrays):
        d = np.linspace(0, 1.0, 100)
        with pytest.raises(ValueError):
            linearize_2d(
                _session_with_xy(tmp_path, synthetic_arrays,
                                 np.column_stack([d, d]), maze="Figure Eight Maze")
            )


def test_with_own_linearization_replaces_only_position_1d(tmp_path, synthetic_arrays):
    d = np.linspace(0, 1.6, 200)
    s = _session_with_xy(tmp_path, synthetic_arrays, np.column_stack([d, np.zeros_like(d)]))
    out = with_own_linearization(s)
    assert out is not s
    np.testing.assert_array_equal(out.position_xy, s.position_xy)
    np.testing.assert_array_equal(out.spike_times, s.spike_times)
    assert np.isnan(s.position_1d).all(), "the original must not be mutated"
    assert np.isfinite(out.position_1d).any()
    with pytest.raises(ValueError):
        out.position_1d[0] = 0.0  # still read-only


@requires_data
@pytest.mark.parametrize("name", ["Achilles_10252013", "Buddy_06272013", "Cicero_09012014",
                                  "Cicero_09172014", "Gatsby_08022013"])
def test_linear_sessions_reproduce_the_authors_linearization(real_sessions, name):
    """On linear tracks ours must be theirs up to an affine transform."""
    s = real_sessions(name)
    ours = linearize_2d(s)
    both = ~np.isnan(ours) & ~np.isnan(s.position_1d)
    assert both.sum() > 1000
    r = np.corrcoef(ours[both], s.position_1d[both])[0, 1]
    assert r > 0.9999, f"{name}: r = {r}"


@requires_data
def test_own_base_covers_more_than_the_authors_mask(real_sessions):
    for name in ["Achilles_10252013", "Buddy_06272013", "Gatsby_08022013"]:
        cov = coverage(real_sessions(name))
        assert cov["own_frac"] > cov["authors_frac"]
