"""Tests for splitting an ensemble into halves, and for the random-split control."""

from __future__ import annotations

import numpy as np
import pytest

from hc11 import subsets


class FakeSession:
    """Just enough session to measure a firing rate."""

    def __init__(self, times_by_id):
        self._times = {int(k): np.sort(np.asarray(v, dtype=float))
                       for k, v in times_by_id.items()}

    def spikes_for(self, cluster_id):
        return self._times.get(int(cluster_id), np.empty(0))


def test_rate_counts_only_spikes_inside_the_intervals():
    session = FakeSession({1: [0.5, 1.5, 2.5, 10.5, 11.5]})
    rate = subsets.rate_in_intervals(session, 1, np.array([[0.0, 3.0], [10.0, 12.0]]))
    assert rate == pytest.approx(5 / 5.0)


def test_rate_is_nan_when_there_is_no_time_to_measure_in():
    session = FakeSession({1: [0.5]})
    assert np.isnan(subsets.rate_in_intervals(session, 1, np.empty((0, 2))))
    assert np.isnan(subsets.rate_in_intervals(session, 1, np.array([[1.0, 1.0]])))


def test_a_silent_cell_has_a_rate_of_zero_not_nan():
    """A place cell that never fires in sleep is slow, not missing."""
    session = FakeSession({1: [100.0]})
    assert subsets.rate_in_intervals(session, 1, np.array([[0.0, 10.0]])) == 0.0


def test_equal_halves_splits_evenly_and_keeps_order():
    ids = np.array([10, 20, 30, 40])
    split = subsets.equal_halves(np.array([0.4, 0.1, 0.9, 0.2]), ids)
    assert list(split["low"]) == [20, 40]
    assert list(split["high"]) == [10, 30]
    assert split["dropped"].size == 0


def test_equal_halves_drops_the_median_cell_when_the_count_is_odd():
    """The paper's subgroups are equal; the middle cell goes to neither side."""
    ids = np.array([1, 2, 3, 4, 5])
    split = subsets.equal_halves(np.array([1.0, 2.0, 3.0, 4.0, 5.0]), ids)
    assert list(split["low"]) == [1, 2]
    assert list(split["high"]) == [4, 5]
    assert list(split["dropped"]) == [3]
    assert len(split["low"]) == len(split["high"])


def test_equal_halves_breaks_ties_deterministically():
    ids = np.array([7, 3, 5, 1])
    first = subsets.equal_halves(np.ones(4), ids)
    second = subsets.equal_halves(np.ones(4), ids)
    assert list(first["low"]) == list(second["low"])


def test_equal_halves_rejects_mismatched_lengths():
    with pytest.raises(ValueError):
        subsets.equal_halves(np.ones(3), np.array([1, 2]))


def test_random_halves_are_equal_and_disjoint():
    ids = np.arange(11)
    split = subsets.random_halves(ids, np.random.default_rng(0))
    assert len(split["low"]) == len(split["high"]) == 5
    assert np.intersect1d(split["low"], split["high"]).size == 0
    assert split["dropped"].size == 1


def test_random_halves_actually_vary():
    ids = np.arange(20)
    rng = np.random.default_rng(0)
    draws = {tuple(subsets.random_halves(ids, rng)["low"]) for _ in range(20)}
    assert len(draws) > 1


def test_split_agreement_is_one_for_a_mirrored_split():
    a = {"low": np.array([1, 2]), "high": np.array([3, 4]), "dropped": np.array([])}
    b = {"low": np.array([3, 4]), "high": np.array([1, 2]), "dropped": np.array([])}
    assert subsets.split_agreement(a, b) == 1.0


def test_split_agreement_counts_only_cells_both_splits_kept():
    a = {"low": np.array([1, 2]), "high": np.array([3, 4]), "dropped": np.array([5])}
    b = {"low": np.array([1, 5]), "high": np.array([3, 4]), "dropped": np.array([2])}
    assert subsets.split_agreement(a, b) == pytest.approx(1.0)


def test_finite_in_both_drops_an_event_either_half_failed_on():
    a = np.array([1.0, np.nan, 3.0, 4.0])
    b = np.array([1.0, 2.0, np.nan, 4.0])
    assert list(subsets.finite_in_both(a, b)) == [True, False, False, True]


def test_finite_in_both_rejects_unpaired_arrays():
    with pytest.raises(ValueError):
        subsets.finite_in_both(np.ones(3), np.ones(2))


def test_empirical_p_counts_draws_at_least_as_extreme():
    assert subsets.empirical_p(0.5, [0.1, 0.2, 0.6, 0.9]) == pytest.approx(3 / 5)
    assert subsets.empirical_p(2.0, [0.1, 0.2, 0.6, 0.9]) == pytest.approx(1 / 5)


def test_empirical_p_is_two_sided_on_magnitude():
    assert subsets.empirical_p(-0.5, [0.6]) == subsets.empirical_p(0.5, [-0.6])


def test_empirical_p_is_nan_for_a_missing_observation():
    """A NaN must not come back as the smallest p the null can express."""
    assert np.isnan(subsets.empirical_p(float("nan"), [0.1, 0.2, 0.3]))


def test_empirical_p_ignores_non_finite_draws():
    assert subsets.empirical_p(0.5, [0.6, np.nan, np.inf]) == pytest.approx(2 / 2)


def test_empirical_p_is_nan_when_every_draw_is_missing():
    assert np.isnan(subsets.empirical_p(0.5, [np.nan, np.nan]))
