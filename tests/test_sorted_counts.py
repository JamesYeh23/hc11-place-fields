"""Tests for step 11a: the two baselines that probe the Phase 4 control's leak.

What matters here is not arithmetic but comparability. Both scripts exist to be
read against numbers produced elsewhere -- the shuffled-fit leak, the rate-only
decoder -- so the tests pin that they compute the *same* statistic rather than
one that merely looks like it.
"""

from __future__ import annotations

import numpy as np
import pytest

import population_rate_baseline as rate
import sorted_count_baseline as sorted_counts
import unfitted_model_baseline as unfitted
from hc11 import hmm


class Line:
    is_circular = False
    extent = 1.0


class Ring:
    is_circular = True
    extent = 1.0


def test_the_sorted_count_scorer_matches_the_rate_decoder_exactly():
    """"The same statistic as the rate-only decoder" is checked, not intended."""
    position = np.array([0.0, 0.25, 0.5, 0.75] * 10)
    level = np.array([0, 1, 2, 3] * 10)
    mine = sorted_counts.score(
        hmm.state_positions(np.eye(4)[level], position, 4, circular=False, extent=1.0),
        level, position, Line(), np.random.default_rng(0))
    theirs = rate.decode(level, position, level, position, 4, Line(),
                         np.random.default_rng(0))
    assert mine["r"] == pytest.approx(theirs[0])
    assert mine["median_error_cm"] == pytest.approx(theirs[1])
    assert mine["chance_error_cm"] == pytest.approx(theirs[2])


def test_sorting_destroys_identity_and_keeps_the_multiset():
    counts = [np.array([[3.0, 0.0, 1.0], [0.0, 2.0, 2.0]])]
    vectors = sorted_counts.sorted_vectors(counts, shape=False)
    assert list(vectors[0]) == [0.0, 1.0, 3.0]
    assert list(vectors[1]) == [0.0, 2.0, 2.0]


def test_the_shape_variant_removes_the_total_the_rate_test_already_covered():
    counts = [np.array([[2.0, 4.0], [1.0, 2.0]])]  # same shape, different total
    vectors = sorted_counts.sorted_vectors(counts, shape=True)
    assert np.allclose(vectors[0], vectors[1])


def test_the_shape_variant_survives_an_empty_bin():
    counts = [np.array([[0.0, 0.0], [1.0, 3.0]])]
    vectors = sorted_counts.sorted_vectors(counts, shape=True)
    assert np.isfinite(vectors).all()
    assert list(vectors[0]) == [0.0, 0.0]


def test_an_unfitted_model_is_marked_as_never_having_been_fitted():
    """Nothing downstream may mistake one of these for a fit."""
    counts = [np.random.default_rng(0).poisson(1.0, size=(40, 6)).astype(float)]
    model = unfitted.unfitted_model("exemplar", counts, 4, 0.25, 0.5,
                                    np.random.default_rng(0))
    assert model.n_iterations == 0
    assert model.seed == -1
    assert np.isnan(model.log_likelihood)


def test_exemplar_rates_are_real_training_bins_and_marginal_rates_are_not():
    """The gap between the two variants is cross-cell structure, so it has to exist."""
    pooled = np.array([[10.0, 0.0], [0.0, 10.0], [5.0, 5.0], [1.0, 1.0]])
    exemplar = unfitted.random_rates("exemplar", [pooled], 3, 1.0,
                                     np.random.default_rng(0))
    # Floored like a fitted model's rates, so a silent cell is 0.01 Hz and not
    # zero; the comparison is against the floored bins, not the raw ones.
    floored = np.maximum(pooled, hmm.RATE_FLOOR_HZ)
    assert all(any(np.allclose(row, bin_) for bin_ in floored) for row in exemplar)
    marginal = unfitted.random_rates("marginal", [pooled], 200, 1.0,
                                     np.random.default_rng(0))
    # Independent draws per cell must produce pairs that never co-occur.
    assert any(row[0] >= 10.0 and row[1] >= 10.0 for row in marginal)


def test_a_sticky_transition_reproduces_the_dwell_time_it_was_given():
    transition = unfitted.sticky_transition(5, dwell_s=1.0, bin_s=0.25)
    dwell = transition.diagonal()[0]
    assert 0.25 / (1.0 - dwell) == pytest.approx(1.0)
    assert transition.sum(axis=1) == pytest.approx(np.ones(5))


def test_a_uniform_transition_has_no_temporal_structure():
    transition = unfitted.uniform_transition(4)
    assert transition == pytest.approx(np.full((4, 4), 0.25))


def test_the_scorer_is_ring_aware():
    """On a ring, 0.98 and 0.02 are close; a linear scorer would call them far."""
    truth = np.array([0.98] * 20)
    mapped = np.array([0.02])
    label = np.zeros(20, dtype=int)
    ring = sorted_counts.score(mapped, label, truth, Ring(),
                               np.random.default_rng(0))
    line = sorted_counts.score(mapped, label, truth, Line(),
                               np.random.default_rng(0))
    assert ring["median_error_cm"] == pytest.approx(4.0)
    assert line["median_error_cm"] == pytest.approx(96.0)
