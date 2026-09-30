"""Track geometry, with the circular wrap as the main target."""

from __future__ import annotations

import numpy as np
import pytest

from hc11.track import Track

LINEAR = Track(kind="linear", lo=0.0, hi=1.6, bin_size=0.1, nominal_length=1.6)
CIRCLE = Track(kind="circular", lo=0.0, hi=2.9, bin_size=0.1)


class TestBinning:
    def test_bin_count_and_edges(self):
        assert LINEAR.n_bins == 16
        assert LINEAR.bin_edges[0] == 0.0 and LINEAR.bin_edges[-1] == 1.6
        assert len(LINEAR.bin_centers) == 16
        np.testing.assert_allclose(LINEAR.actual_bin_size, 0.1)

    def test_bins_tile_a_non_multiple_extent(self):
        t = Track(kind="circular", lo=0.0, hi=2.8416, bin_size=0.1)
        assert t.n_bins == 28
        assert t.bin_edges[-1] == pytest.approx(2.8416)
        assert t.actual_bin_size == pytest.approx(2.8416 / 28)

    def test_digitize_within_range(self):
        idx = LINEAR.digitize(np.array([0.0, 0.05, 0.15, 1.59]))
        np.testing.assert_array_equal(idx, [0, 0, 1, 15])

    def test_digitize_nan_is_minus_one(self):
        np.testing.assert_array_equal(LINEAR.digitize(np.array([np.nan, 0.5])), [-1, 5])

    def test_digitize_endpoint_linear_vs_circular(self):
        assert LINEAR.digitize(np.array([1.6]))[0] == 15  # clipped into the last bin
        assert CIRCLE.digitize(np.array([2.9]))[0] == 0  # same place as 0.0


class TestWrap:
    """The reward site is at position 0, so 0 and `period` are the same place."""

    def test_distance_takes_the_short_way(self):
        assert CIRCLE.distance(0.05, 2.85) == pytest.approx(0.10)
        assert abs(0.05 - 2.85) == pytest.approx(2.80)  # what a naive subtraction gives

    def test_linear_track_does_not_wrap(self):
        assert LINEAR.distance(0.05, 1.55) == pytest.approx(1.50)

    def test_difference_keeps_its_sign(self):
        assert CIRCLE.difference(0.05, 2.85) == pytest.approx(0.10)  # forward past 0
        assert CIRCLE.difference(2.85, 0.05) == pytest.approx(-0.10)  # backward past 0

    def test_distance_is_symmetric_and_bounded(self):
        rng = np.random.default_rng(0)
        a, b = rng.uniform(0, 2.9, 500), rng.uniform(0, 2.9, 500)
        np.testing.assert_allclose(CIRCLE.distance(a, b), CIRCLE.distance(b, a))
        assert CIRCLE.distance(a, b).max() <= 2.9 / 2 + 1e-9

    def test_distance_to_self_is_zero(self):
        x = np.linspace(0, 2.9, 50)
        np.testing.assert_allclose(CIRCLE.distance(x, x), 0.0, atol=1e-12)

    def test_wrap_folds_into_range(self):
        np.testing.assert_allclose(CIRCLE.wrap(np.array([-0.1, 0.0, 2.9, 3.0])),
                                   [2.8, 0.0, 0.0, 0.1])

    def test_unwrap_makes_laps_continuous(self):
        # Two full laps crossing the reward site, sampled forwards.
        true = np.linspace(0.0, 5.8, 200)
        wrapped = CIRCLE.wrap(true)
        out = CIRCLE.unwrap(wrapped)
        np.testing.assert_allclose(out - out[0], true - true[0], atol=1e-9)
        assert np.all(np.diff(out) > 0), "unwrapped position must be monotonic here"

    def test_unwrap_does_not_bridge_a_tracking_gap(self):
        x = CIRCLE.wrap(np.linspace(0.0, 5.8, 200))
        x[80:120] = np.nan
        out = CIRCLE.unwrap(x)
        assert np.all(np.isnan(out[80:120]))
        # Each valid run is unwrapped on its own; no jump is inferred across the gap.
        assert np.all(np.diff(out[:80]) > 0)
        assert np.all(np.diff(out[120:]) > 0)

    def test_unwrap_is_identity_on_linear(self):
        x = np.array([0.1, 0.5, np.nan, 1.2])
        np.testing.assert_array_equal(LINEAR.unwrap(x), x)


def test_from_session_uses_observed_range(synthetic_file):
    from hc11.io import load_session

    s = load_session(synthetic_file)
    t = Track.from_session(s, bin_size_cm=10.0)
    assert t.kind == "linear"
    assert t.lo == pytest.approx(np.nanmin(s.position_1d))
    assert t.hi == pytest.approx(np.nanmax(s.position_1d))
    assert t.nominal_length == 1.6
