"""Tests for the per-cell contribution estimator (step 9c)."""

from __future__ import annotations

import numpy as np
import pytest

from hc11 import pcc


def test_shift_field_moves_one_cell_and_leaves_the_rest_alone():
    tuning = {"pos": np.array([[1.0, 2.0, 3.0, 4.0], [9.0, 9.0, 9.0, 9.0]])}
    out = pcc.shift_field(tuning, 0, 1)
    assert list(out["pos"][0]) == [4.0, 1.0, 2.0, 3.0]
    assert list(out["pos"][1]) == [9.0, 9.0, 9.0, 9.0]


def test_shift_field_uses_the_same_offset_in_every_direction():
    tuning = {"pos": np.array([[1.0, 2.0, 3.0]]), "neg": np.array([[4.0, 5.0, 6.0]])}
    out = pcc.shift_field(tuning, 0, 1)
    assert list(out["pos"][0]) == [3.0, 1.0, 2.0]
    assert list(out["neg"][0]) == [6.0, 4.0, 5.0]


def test_shift_field_preserves_the_cells_rate_and_shape():
    """A contribution measured against the shift must be about position, not rate."""
    field = np.array([[0.1, 5.0, 0.2, 0.0]])
    out = pcc.shift_field({"pos": field}, 0, 2)["pos"][0]
    assert sorted(out) == sorted(field[0])
    assert out.max() == field[0].max()


def test_shift_field_of_zero_is_the_identity_object():
    tuning = {"pos": np.ones((2, 3))}
    assert pcc.shift_field(tuning, 0, 0) is tuning


def test_standardise_needs_a_null_with_width():
    assert pcc.standardise(1.0, 0.5, 0.25) == pytest.approx(2.0)
    assert np.isnan(pcc.standardise(1.0, 0.5, 0.0))
    assert np.isnan(pcc.standardise(1.0, 0.5, float("nan")))


def test_contribution_is_the_gap_between_the_real_field_and_the_average_placement():
    result = pcc.contribution(1.0, np.array([0.0, 0.0, 0.0, 0.0]))
    assert result["pcc"] == pytest.approx(1.0)
    assert result["p"] == pytest.approx(1 / 5)


def test_contribution_reports_a_cell_that_is_no_better_than_its_shifts():
    shifts = np.array([0.0, 1.0, 2.0, 3.0, 4.0])
    result = pcc.contribution(2.0, shifts)          # exactly the mean placement
    assert result["pcc"] == pytest.approx(0.0)
    assert result["p"] == pytest.approx(1.0)        # every shift is at least as extreme


def test_contribution_z_scales_with_the_spread_of_placements():
    tight = pcc.contribution(1.0, np.array([0.0, 0.1, -0.1, 0.0]))
    loose = pcc.contribution(1.0, np.array([0.0, 1.0, -1.0, 0.0]))
    assert tight["z"] > loose["z"]


def test_contribution_survives_missing_shifts():
    result = pcc.contribution(1.0, np.array([0.0, np.nan, 0.0]))
    assert result["n_shifts"] == 2
    assert np.isfinite(result["pcc"])


def test_contribution_is_nan_without_enough_shifts():
    assert np.isnan(pcc.contribution(1.0, np.array([np.nan, 0.5]))["pcc"])
    assert np.isnan(pcc.contribution(float("nan"), np.arange(5.0))["pcc"])


def test_null_contributions_use_each_shift_as_a_stand_in_for_the_real_field():
    values = np.array([1.0, 2.0, 3.0])
    out = pcc.null_contributions(values)
    assert out[0] == pytest.approx(1.0 - 2.5)
    assert out[1] == pytest.approx(2.0 - 2.0)
    assert out[2] == pytest.approx(3.0 - 1.5)


def test_null_contributions_are_centred_on_zero_for_exchangeable_shifts():
    rng = np.random.default_rng(0)
    draws = rng.normal(0.4, 0.1, 500)
    assert abs(float(np.mean(pcc.null_contributions(draws)))) < 0.02


def test_null_contributions_match_the_real_statistic_in_form():
    """Leaving shift s out and scoring it is exactly what contribution() does."""
    values = np.array([0.3, 0.7, 0.1, 0.5])
    nulls = pcc.null_contributions(values)
    for s in range(values.size):
        others = np.delete(values, s)
        assert nulls[s] == pytest.approx(pcc.contribution(values[s], others)["pcc"])


def test_null_contributions_need_two_usable_shifts():
    assert np.all(np.isnan(pcc.null_contributions(np.array([1.0, np.nan]))))


def test_finite_everywhere_keeps_only_events_every_placement_scored():
    stack = np.array([[1.0, 2.0, 3.0], [1.0, np.nan, 3.0], [0.5, 2.5, 3.5]])
    assert list(pcc.finite_everywhere(stack)) == [True, False, True]


def test_finite_everywhere_handles_a_single_placement():
    assert list(pcc.finite_everywhere(np.array([1.0, np.nan]))) == [True, False]


def test_a_silent_cell_still_changes_the_posterior():
    """The flagged silent-failure case: a cell with no spikes is not zero-contribution.

    Its field never enters the spike term of the decoder, but it does enter the
    normalisation term ``-dt * sum_i f_i(x)``, so moving it still moves the
    posterior. If this ever stops being true, PCC for silent cells collapses to
    exactly zero and the all-events mean (the primary one) is quietly wrong.
    """
    # Two cells, four positions, two time bins. Cell 0 fires, cell 1 is silent.
    counts = np.array([[1, 0], [0, 0]])
    tuning = {"pos": np.array([[5.0, 1.0, 1.0, 1.0], [1.0, 1.0, 1.0, 9.0]])}
    positions = np.array([0.0, 1.0, 2.0, 3.0])
    times = np.array([0.0, 0.02])

    assert counts[1].sum() == 0, "cell 1 must be silent for this test to mean anything"
    real = pcc.best_abs_r(counts, tuning, positions, times, 0.02,
                          circular=False, rate_floor_hz=0.01)
    moved = pcc.best_abs_r(counts, pcc.shift_field(tuning, 1, 2), positions, times,
                           0.02, circular=False, rate_floor_hz=0.01)
    assert real != moved, "shifting a silent cell's field left the score unchanged"


def test_a_silent_cells_contribution_is_not_structurally_zero():
    """The same point one level up: the estimator must not return exactly 0.0."""
    counts = np.array([[1, 0, 1], [0, 0, 0]])
    tuning = {"pos": np.array([[6.0, 1.0, 1.0, 1.0], [1.0, 1.0, 1.0, 7.0]])}
    positions = np.array([0.0, 1.0, 2.0, 3.0])
    times = np.array([0.0, 0.02, 0.04])
    scores = [pcc.best_abs_r(counts, pcc.shift_field(tuning, 1, offset), positions,
                             times, 0.02, circular=False, rate_floor_hz=0.01)
              for offset in range(4)]
    assert len(set(scores)) > 1, "every placement of a silent cell's field scored alike"
    result = pcc.contribution(scores[0], np.array(scores[1:]))
    assert result["pcc"] != 0.0
