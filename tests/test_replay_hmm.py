"""Tests for step 10d: scoring sleep events with the state model.

The scoring machinery itself is Phase 1's and is tested in `test_sequences.py`.
What is new here is the adapter that hands it an HMM posterior instead of a
place-field posterior, and the comparison's rules about controls and vacuity.
"""

from __future__ import annotations

import numpy as np
import pytest

import replay_comparison_hmm as comparison
import score_events_hmm as scoring
from hc11 import hmm


class Ring:
    is_circular = True
    extent = 1.0


class Line:
    is_circular = False
    extent = 1.0


def ramp_model(n_states=6, n_cells=6, bin_s=0.02):
    """A model whose state k is the cell-k-fires state: a readable posterior."""
    rates = np.full((n_states, n_cells), 0.5)
    np.fill_diagonal(rates, 400.0)
    transition = np.full((n_states, n_states), 1.0 / n_states)
    return hmm.PoissonHMM(initial=np.full(n_states, 1.0 / n_states),
                          transition=transition, rates=rates, bin_s=bin_s,
                          log_likelihood=0.0, n_iterations=1, converged=True, seed=0)


def entry(model):
    order = np.arange(model.n_states)
    return {"model": model, "position": order, "spectral": order}


def test_a_clean_state_ramp_scores_far_above_its_shuffles():
    model = ramp_model()
    counts = np.eye(6)  # bin k has one spike from cell k: states sweep in order
    grid = np.arange(6) / 6.0
    result = scoring.score(entry(model), "position", counts, np.arange(6), grid,
                           Line(), 0.02, "column_cycle", 200,
                           np.random.default_rng(0))
    # A six-bin event cannot give a large z however clean it is: with so few
    # bins the column-cycle null itself reaches |r| ~ 0.37 on average. The
    # correlation is the strong assertion here and the z is the weak one.
    assert result.r > 0.9
    assert result.rz > 2.0


def test_a_scrambled_state_path_does_not():
    model = ramp_model()
    counts = np.eye(6)[[3, 0, 4, 1, 5, 2]]
    grid = np.arange(6) / 6.0
    result = scoring.score(entry(model), "position", counts, np.arange(6), grid,
                           Line(), 0.02, "column_cycle", 200,
                           np.random.default_rng(0))
    assert abs(result.r) < 0.6


def test_the_ordering_argument_actually_changes_the_posterior_columns():
    """A reversed ordering has to flip the sign of the correlation, not ignore it."""
    model = ramp_model()
    counts = np.eye(6)
    grid = np.arange(6) / 6.0
    straight = entry(model)
    reversed_order = {**straight, "spectral": np.arange(6)[::-1]}
    forward = scoring.score(straight, "position", counts, np.arange(6), grid, Line(),
                            0.02, "column_cycle", 50, np.random.default_rng(0))
    backward = scoring.score(reversed_order, "spectral", counts, np.arange(6), grid,
                             Line(), 0.02, "column_cycle", 50,
                             np.random.default_rng(0))
    assert forward.r == pytest.approx(-backward.r, abs=1e-9)


def test_the_scoring_posterior_is_memoryless():
    """A sticky transition matrix must not smooth the event's state path.

    The model is fitted at 250 ms and events are scored at 20 ms, so applying
    the fitted transitions here would let the state move twelve times too fast
    -- or, with a sticky matrix, not at all.
    """
    model = ramp_model()
    sticky = hmm.PoissonHMM(initial=model.initial,
                            transition=np.eye(model.n_states) * 0.999
                            + 0.001 / model.n_states,
                            rates=model.rates, bin_s=model.bin_s,
                            log_likelihood=0.0, n_iterations=1, converged=True,
                            seed=0)
    counts = np.eye(6)
    loose = hmm.memoryless_posterior(model, counts, bin_s=0.02)
    tight = hmm.memoryless_posterior(sticky, counts, bin_s=0.02)
    assert np.allclose(loose, tight)


def test_event_counts_drops_bins_with_no_spikes_and_reports_the_total():
    class FakeSession:
        spike_times = np.array([0.001, 0.0011, 0.081])
        spike_ids = np.array([1, 2, 1])

    counts, bin_index, n_total = scoring.event_counts(
        FakeSession(), 0.0, 0.1, np.array([1, 2]), 0.02)
    assert n_total == 5
    assert list(bin_index) == [0, 4]
    assert counts.shape == (2, 2)


def test_post_above_pre_is_a_positive_rank_biserial():
    result = comparison.compare(np.zeros(50), np.ones(50))
    assert result["rank_biserial"] == pytest.approx(1.0)
    assert result["difference"] == pytest.approx(1.0)


def test_pre_above_post_is_a_negative_rank_biserial():
    result = comparison.compare(np.ones(50), np.zeros(50))
    assert result["rank_biserial"] == pytest.approx(-1.0)


def test_subsampling_is_skipped_when_pre_is_not_the_larger_set():
    assert comparison.subsample(np.zeros(10), np.ones(20), 5,
                                np.random.default_rng(0)) == {}


def test_subsampling_cuts_pre_to_post_s_size():
    rng = np.random.default_rng(0)
    pre, post = rng.normal(size=400), rng.normal(size=40) + 1.0
    drawn = comparison.subsample(pre, post, 50, rng)
    assert drawn["subsample_difference_p5"] < drawn["subsample_median_difference"]
    assert drawn["subsample_median_rank_biserial"] > 0


def test_the_vacuous_sessions_are_the_two_smallest_ensembles():
    """Marked in advance from D10.14, not chosen after seeing the sleep result."""
    assert set(comparison.VACUOUS) == {"Gatsby_08022013", "Cicero_09012014"}
    assert all(value < 0.1 for value in comparison.VACUOUS.values())


def test_an_unreadable_stamp_refuses_rather_than_silently_starting_empty(tmp_path):
    """Swallowing the error would drop every other session's provenance."""
    missing = tmp_path / "nothing.json"
    with pytest.raises(SystemExit):
        scoring.read_retrying(missing, attempts=2, pause=0.0)
