"""Tests for event binning, the weighted correlation, the shuffles and rZ.

The weighted correlation is the piece everything else rests on, so it is pinned
against cases where the answer is known by construction: a posterior
concentrated on a perfect diagonal must give r = 1, the same posterior with its
time axis reversed must give r = -1, and a posterior with no position-time
relationship must give r near 0. The rZ tests then check that a real sequence
scores high and a structureless event scores near zero *through the shuffles*,
which is the only thing the Fig. 1C comparison actually depends on.
"""

from __future__ import annotations

from types import SimpleNamespace

import numpy as np
import pytest

from hc11 import sequences as sq

POSITIONS = np.arange(32) * 5.0  # 32 bins of 5 cm
BIN_S = 0.02


def diagonal_posterior(n_bins=8, n_position=32, spread=1.0, reverse=False, start=4,
                       step=3):
    """A posterior whose mass walks along the track at a constant rate."""
    centres = start + step * np.arange(n_bins)
    if reverse:
        centres = centres[::-1]
    grid = np.arange(n_position)[None, :]
    weights = np.exp(-0.5 * ((grid - centres[:, None]) / spread) ** 2)
    return weights / weights.sum(axis=1, keepdims=True)


def flat_posterior(n_bins=8, n_position=32):
    return np.full((n_bins, n_position), 1.0 / n_position)


def scattered_posterior(n_bins=8, n_position=32, seed=0):
    """Concentrated per bin, but at positions with no relation to bin order."""
    rng = np.random.default_rng(seed)
    centres = rng.permutation(n_position)[:n_bins]
    grid = np.arange(n_position)[None, :]
    weights = np.exp(-0.5 * ((grid - centres[:, None]) / 1.0) ** 2)
    return weights / weights.sum(axis=1, keepdims=True)


def times_for(posterior):
    return np.arange(posterior.shape[0]) * BIN_S


# --- binning --------------------------------------------------------------


def spike_session(spikes):
    spikes = sorted(spikes, key=lambda pair: pair[1])
    return SimpleNamespace(
        spike_ids=np.array([c for c, _ in spikes], dtype=np.int64),
        spike_times=np.array([t for _, t in spikes], dtype=float),
    )


def test_binning_counts_spikes_into_the_right_bins():
    session = spike_session([(1, 0.005), (1, 0.015), (2, 0.025), (3, 0.095)])
    counts, centres = sq.bin_event(session, 0.0, 0.1, np.array([1, 2, 3]), 0.02)
    assert counts.shape == (3, 5)
    assert counts[0].tolist() == [2, 0, 0, 0, 0]
    assert counts[1].tolist() == [0, 1, 0, 0, 0]
    assert counts[2].tolist() == [0, 0, 0, 0, 1]
    assert centres[0] == pytest.approx(0.01)


def test_binning_drops_a_final_partial_bin():
    """A short last bin would have a different expected count at the same rate."""
    session = spike_session([(1, 0.001)])
    counts, _ = sq.bin_event(session, 0.0, 0.055, np.array([1]), 0.02)
    assert counts.shape[1] == 2


def test_binning_ignores_cells_outside_the_ensemble():
    session = spike_session([(1, 0.005), (99, 0.006)])
    counts, _ = sq.bin_event(session, 0.0, 0.02, np.array([1]), 0.02)
    assert counts.sum() == 1


def test_an_event_shorter_than_one_bin_has_no_bins():
    session = spike_session([(1, 0.001)])
    counts, centres = sq.bin_event(session, 0.0, 0.01, np.array([1]), 0.02)
    assert counts.shape[1] == 0 and centres.size == 0


# --- weighted correlation -------------------------------------------------


def test_a_perfect_forward_diagonal_scores_one():
    posterior = diagonal_posterior(spread=0.4)
    r, offset = sq.weighted_correlation(posterior, POSITIONS, times_for(posterior))
    assert r == pytest.approx(1.0, abs=0.01)
    assert offset == 0


def test_a_perfect_reverse_diagonal_scores_minus_one():
    posterior = diagonal_posterior(spread=0.4, reverse=True)
    r, _ = sq.weighted_correlation(posterior, POSITIONS, times_for(posterior))
    assert r == pytest.approx(-1.0, abs=0.01)


def test_a_scattered_posterior_scores_near_zero():
    posterior = scattered_posterior(seed=3)
    r, _ = sq.weighted_correlation(posterior, POSITIONS, times_for(posterior))
    assert abs(r) < 0.6


def test_a_flat_posterior_has_no_correlation():
    posterior = flat_posterior()
    r, _ = sq.weighted_correlation(posterior, POSITIONS, times_for(posterior))
    assert r == pytest.approx(0.0, abs=1e-9)


def test_a_broad_posterior_scores_lower_than_a_sharp_one_on_the_same_path():
    sharp = diagonal_posterior(spread=0.4)
    broad = diagonal_posterior(spread=6.0)
    r_sharp, _ = sq.weighted_correlation(sharp, POSITIONS, times_for(sharp))
    r_broad, _ = sq.weighted_correlation(broad, POSITIONS, times_for(broad))
    assert abs(r_broad) < abs(r_sharp)


def test_the_correlation_is_vectorised_over_leading_axes():
    stack = np.stack([diagonal_posterior(spread=0.4),
                      diagonal_posterior(spread=0.4, reverse=True),
                      flat_posterior()])
    r, _ = sq.weighted_correlation(stack, POSITIONS, times_for(stack[0]))
    assert r.shape == (3,)
    assert r[0] == pytest.approx(1.0, abs=0.01)
    assert r[1] == pytest.approx(-1.0, abs=0.01)
    assert r[2] == pytest.approx(0.0, abs=1e-9)


# --- circular tracks ------------------------------------------------------


def test_a_wrapping_trajectory_scores_near_zero_on_a_linear_axis():
    """The case the circular handling exists for: a good sequence across the seam."""
    posterior = diagonal_posterior(n_bins=8, start=26, step=3)  # 26..47 -> wraps
    wrapped = np.roll(posterior, 0, axis=1)
    wrapped = np.stack([np.roll(row, 0) for row in wrapped])
    # rebuild with explicit wrap
    centres = (26 + 3 * np.arange(8)) % 32
    grid = np.arange(32)[None, :]
    delta = np.abs(grid - centres[:, None])
    delta = np.minimum(delta, 32 - delta)
    wrapped = np.exp(-0.5 * (delta / 0.4) ** 2)
    wrapped /= wrapped.sum(axis=1, keepdims=True)
    linear, _ = sq.weighted_correlation(wrapped, POSITIONS, times_for(wrapped),
                                        circular=False)
    circular, _ = sq.weighted_correlation(wrapped, POSITIONS, times_for(wrapped),
                                          circular=True)
    assert abs(linear) < 0.75
    assert abs(circular) > 0.98


def test_the_circular_offset_is_reported():
    centres = (26 + 3 * np.arange(8)) % 32
    grid = np.arange(32)[None, :]
    delta = np.abs(grid - centres[:, None])
    delta = np.minimum(delta, 32 - delta)
    posterior = np.exp(-0.5 * (delta / 0.4) ** 2)
    posterior /= posterior.sum(axis=1, keepdims=True)
    _, offset = sq.weighted_correlation(posterior, POSITIONS, times_for(posterior),
                                        circular=True)
    assert 0 < int(offset) < 32


def test_circular_scoring_leaves_a_non_wrapping_sequence_alone():
    posterior = diagonal_posterior(spread=0.4, start=4, step=2)
    linear, _ = sq.weighted_correlation(posterior, POSITIONS, times_for(posterior))
    circular, _ = sq.weighted_correlation(posterior, POSITIONS, times_for(posterior),
                                          circular=True)
    assert abs(circular) >= abs(linear) - 1e-9
    assert abs(circular) == pytest.approx(1.0, abs=0.01)


def test_offset_positions_has_one_column_when_linear():
    assert sq.offset_positions(POSITIONS, circular=False).shape == (32, 1)
    assert sq.offset_positions(POSITIONS, circular=True).shape == (32, 32)


# --- shuffles -------------------------------------------------------------


def test_column_cycle_preserves_each_bins_posterior_shape():
    rng = np.random.default_rng(0)
    posterior = diagonal_posterior()
    prepared = sq.shuffle_indices("column_cycle", 50, posterior.shape[0], 32, rng)
    shuffled = sq.apply_shuffle(posterior, prepared)
    assert shuffled.shape == (50, posterior.shape[0], 32)
    assert np.allclose(np.sort(shuffled, axis=2), np.sort(posterior, axis=1), atol=1e-12)
    assert np.allclose(shuffled.sum(axis=2), 1.0)


def test_time_permutation_preserves_the_set_of_bins():
    rng = np.random.default_rng(0)
    posterior = diagonal_posterior()
    prepared = sq.shuffle_indices("time_permutation", 50, posterior.shape[0], 32, rng)
    shuffled = sq.apply_shuffle(posterior, prepared)
    for realisation in shuffled:
        assert np.allclose(np.sort(realisation, axis=0), np.sort(posterior, axis=0))


def test_both_shuffles_destroy_a_perfect_sequence():
    rng = np.random.default_rng(1)
    posterior = diagonal_posterior(spread=0.4)
    times = times_for(posterior)
    for kind in sq.SHUFFLES:
        prepared = sq.shuffle_indices(kind, 400, posterior.shape[0], 32, rng)
        r, _ = sq.weighted_correlation(sq.apply_shuffle(posterior, prepared),
                                       POSITIONS, times)
        assert np.abs(r).mean() < 0.75


def test_an_unknown_shuffle_is_refused():
    with pytest.raises(ValueError, match="unknown shuffle"):
        sq.shuffle_indices("rotate", 10, 8, 32, np.random.default_rng(0))


# --- rZ -------------------------------------------------------------------


def one_event(posterior, n_bins_total=None):
    n_bins = posterior.shape[0]
    return sq.EventPosteriors(
        posteriors={"pos": posterior}, bin_index=np.arange(n_bins),
        spikes_per_bin=np.ones(n_bins), n_bins_total=n_bins_total or n_bins,
        n_cells_active=n_bins,
    )


@pytest.mark.parametrize("kind", sq.SHUFFLES)
def test_a_real_sequence_scores_high(kind):
    score = sq.score_event(
        one_event(diagonal_posterior(n_bins=10, spread=0.5)), POSITIONS, BIN_S,
        circular=False, shuffle=kind, n_shuffles=400, rng=np.random.default_rng(2),
    )
    assert score.rz > 3.0
    assert score.forward


@pytest.mark.parametrize("kind", sq.SHUFFLES)
def test_a_structureless_event_scores_near_zero(kind):
    scores = [
        sq.score_event(
            one_event(scattered_posterior(n_bins=10, seed=seed)), POSITIONS, BIN_S,
            circular=False, shuffle=kind, n_shuffles=400,
            rng=np.random.default_rng(100 + seed),
        ).rz
        for seed in range(25)
    ]
    assert abs(np.median(scores)) < 1.0


def test_a_reverse_sequence_scores_high_and_is_flagged_reverse():
    score = sq.score_event(
        one_event(diagonal_posterior(n_bins=10, spread=0.5, reverse=True)),
        POSITIONS, BIN_S, circular=False, shuffle="column_cycle", n_shuffles=400,
        rng=np.random.default_rng(3),
    )
    assert score.rz > 3.0
    assert not score.forward
    assert score.r < 0


def test_the_direction_with_the_larger_correlation_is_chosen():
    event = sq.EventPosteriors(
        posteriors={"pos": scattered_posterior(n_bins=10, seed=5),
                    "neg": diagonal_posterior(n_bins=10, spread=0.5)},
        bin_index=np.arange(10), spikes_per_bin=np.ones(10), n_bins_total=10,
        n_cells_active=10,
    )
    score = sq.score_event(event, POSITIONS, BIN_S, circular=False,
                           shuffle="column_cycle", n_shuffles=200,
                           rng=np.random.default_rng(4))
    assert score.direction == "neg"
    assert score.r_by_direction["neg"] > abs(score.r_by_direction["pos"])


def test_the_null_takes_the_same_maximum_over_directions():
    """Two directions must give a higher null than one, or rZ is inflated."""
    posterior = scattered_posterior(n_bins=10, seed=7)
    other = scattered_posterior(n_bins=10, seed=8)
    single = sq.score_event(one_event(posterior), POSITIONS, BIN_S, circular=False,
                            shuffle="column_cycle", n_shuffles=400,
                            rng=np.random.default_rng(9))
    both = sq.score_event(
        sq.EventPosteriors(posteriors={"pos": posterior, "neg": other},
                           bin_index=np.arange(10), spikes_per_bin=np.ones(10),
                           n_bins_total=10, n_cells_active=10),
        POSITIONS, BIN_S, circular=False, shuffle="column_cycle", n_shuffles=400,
        rng=np.random.default_rng(9),
    )
    assert both.null_mean > single.null_mean


def test_an_event_with_no_posterior_scores_nothing():
    empty = sq.EventPosteriors({}, np.empty(0, dtype=int), np.empty(0), 0, 0)
    assert sq.score_event(empty, POSITIONS, BIN_S, circular=False,
                          shuffle="column_cycle", n_shuffles=10,
                          rng=np.random.default_rng(0)) is None


def test_empty_bins_are_dropped_and_counted():
    """A silent bin says something about the place fields, not about the animal."""
    session = spike_session([(1, 0.001), (1, 0.085)])
    tuning = {"pos": np.full((1, 32), 5.0)}
    event = sq.event_posteriors(session, 0.0, 0.1, np.array([1]), tuning, 0.02)
    assert event.n_bins == 2
    assert event.n_bins_total == 5
    assert event.bin_index.tolist() == [0, 4]
    assert event.empty_fraction == pytest.approx(0.6)


def test_bin_index_keeps_real_timing_through_a_gap():
    """Kept bins carry their original index, so a gap stays a gap in the score."""
    session = spike_session([(1, 0.001), (1, 0.085)])
    tuning = {"pos": np.full((1, 32), 5.0)}
    event = sq.event_posteriors(session, 0.0, 0.1, np.array([1]), tuning, 0.02)
    assert (event.bin_index * 0.02).tolist() == [0.0, pytest.approx(0.08)]
