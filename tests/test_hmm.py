"""The Poisson hidden Markov model: recovery, arithmetic, and the two orderings."""

from __future__ import annotations

import numpy as np
import pytest

from hc11 import hmm


def simulate(n_bins, rates, transition, bin_s, rng, initial=None):
    """Draw a state path and Poisson counts from a known model."""
    n_states = rates.shape[0]
    initial = np.full(n_states, 1 / n_states) if initial is None else initial
    states = np.empty(n_bins, dtype=int)
    states[0] = rng.choice(n_states, p=initial)
    for t in range(1, n_bins):
        states[t] = rng.choice(n_states, p=transition[states[t - 1]])
    counts = rng.poisson(rates[states] * bin_s)
    return states, counts.astype(float)


def ring_model(n_states=4, n_cells=8, peak=20.0, stay=0.6):
    """Each state fires two cells hard and the rest barely; states form a ring."""
    rates = np.full((n_states, n_cells), 0.2)
    for state in range(n_states):
        rates[state, (2 * state) % n_cells] = peak
        rates[state, (2 * state + 1) % n_cells] = peak
    transition = np.full((n_states, n_states), (1 - stay) / (n_states - 1))
    np.fill_diagonal(transition, stay)
    return rates, transition


def test_forward_backward_matches_brute_force_on_a_tiny_model():
    """Two states, three bins: enumerate all eight paths and compare."""
    rng = np.random.default_rng(0)
    rates = np.array([[5.0, 0.5], [0.5, 5.0]])
    transition = np.array([[0.7, 0.3], [0.4, 0.6]])
    initial = np.array([0.6, 0.4])
    counts = rng.poisson(rates[[0, 1, 0]] * 0.25).astype(float)
    log_b = hmm.log_emissions(counts, rates, 0.25)
    gamma, _, log_likelihood = hmm.forward_backward(log_b, initial, transition)

    total, occupancy = 0.0, np.zeros((3, 2))
    for a in range(2):
        for b in range(2):
            for c in range(2):
                path = (a, b, c)
                p = initial[a] * transition[a, b] * transition[b, c]
                p *= np.exp(log_b[0, a] + log_b[1, b] + log_b[2, c])
                total += p
                for t, state in enumerate(path):
                    occupancy[t, state] += p
    assert log_likelihood == pytest.approx(np.log(total), rel=1e-9)
    assert gamma == pytest.approx(occupancy / total, rel=1e-9)


def test_em_likelihood_never_decreases():
    rng = np.random.default_rng(1)
    rates, transition = ring_model()
    _, counts = simulate(400, rates, transition, 0.25, rng)
    history = []
    for n_iterations in range(1, 9):
        model = hmm.fit([counts], 4, 0.25, seed=3, n_iterations=n_iterations, tol=0.0)
        history.append(model.log_likelihood)
    assert all(b >= a - 1e-6 for a, b in zip(history, history[1:], strict=False))


def recovery_error(model, rates):
    """Worst relative rate error, matching each true state to its nearest fit.

    States are permutation-ambiguous, so the comparison has to be by content.
    """
    error = []
    for true_rate in rates:
        distance = np.abs(model.rates - true_rate[None, :]).sum(axis=1)
        error.append(distance.min() / true_rate.sum())
    return max(error)


def test_em_recovers_rates_it_generated_given_restarts():
    """The best of five restarts recovers the model; a single restart may not."""
    rng = np.random.default_rng(2)
    rates, transition = ring_model()
    _, counts = simulate(3000, rates, transition, 0.25, rng)
    fits = [hmm.fit([counts], 4, 0.25, seed=seed, n_iterations=200) for seed in range(5)]
    best = max(fits, key=lambda m: m.log_likelihood)
    assert recovery_error(best, rates) < 0.1
    assert best.n_states == 4 and best.n_cells == 8


def test_restarts_can_disagree_which_is_why_all_five_are_reported():
    """EM on this model has local optima: reporting only the best hides that."""
    rng = np.random.default_rng(2)
    rates, transition = ring_model()
    _, counts = simulate(3000, rates, transition, 0.25, rng)
    fits = [hmm.fit([counts], 4, 0.25, seed=seed, n_iterations=200) for seed in range(5)]
    likelihoods = [m.log_likelihood for m in fits]
    assert max(likelihoods) - min(likelihoods) > 1.0, "no spread to report"
    errors = [recovery_error(m, rates) for m in fits]
    assert max(errors) > 0.1 > min(errors), "some restarts should miss and some should not"


def test_sequences_are_never_concatenated():
    """A transition must never be counted across the gap between two sequences."""
    rng = np.random.default_rng(4)
    rates, transition = ring_model(n_states=2, n_cells=4)
    pieces = [simulate(50, rates, transition, 0.25, rng)[1] for _ in range(6)]
    joined = hmm.fit([np.vstack(pieces)], 2, 0.25, seed=0, n_iterations=30)
    apart = hmm.fit(pieces, 2, 0.25, seed=0, n_iterations=30)
    assert not np.allclose(joined.transition, apart.transition)


def test_single_bin_sequences_are_dropped_not_crashed_on():
    rng = np.random.default_rng(5)
    rates, transition = ring_model(n_states=2, n_cells=4)
    good = simulate(60, rates, transition, 0.25, rng)[1]
    model = hmm.fit([good, good[:1]], 2, 0.25, seed=0, n_iterations=10)
    assert model.n_states == 2


def test_fit_refuses_when_nothing_is_long_enough():
    with pytest.raises(ValueError):
        hmm.fit([np.zeros((1, 3))], 2, 0.25)


def test_a_large_ensemble_stays_finite():
    rng = np.random.default_rng(6)
    rates, transition = ring_model(n_states=6, n_cells=90, peak=30.0)
    _, counts = simulate(300, rates, transition, 0.25, rng)
    gamma, _, log_likelihood = hmm.forward_backward(
        hmm.log_emissions(counts, rates, 0.25), np.full(6, 1 / 6), transition)
    assert np.isfinite(log_likelihood)
    assert np.allclose(gamma.sum(axis=1), 1.0)


def test_the_per_bin_offset_is_added_back_exactly():
    """Emissions are scaled per bin before exponentiating; the constant must return.

    Shifting every bin's log-emissions by the same amount must leave the state
    posteriors untouched and move the log likelihood by exactly n_bins x shift.
    If the offset were dropped, a 90-cell ensemble would silently report a
    likelihood hundreds of nats wrong and every model comparison with it.
    """
    rng = np.random.default_rng(11)
    rates, transition = ring_model(n_states=3, n_cells=6)
    _, counts = simulate(120, rates, transition, 0.25, rng)
    log_b = hmm.log_emissions(counts, rates, 0.25)
    initial = np.full(3, 1 / 3)
    gamma, _, base = hmm.forward_backward(log_b, initial, transition)
    shifted_gamma, _, shifted = hmm.forward_backward(log_b - 50.0, initial, transition)
    assert shifted == pytest.approx(base - 50.0 * len(counts), rel=1e-9)
    assert shifted_gamma == pytest.approx(gamma, rel=1e-9)


def test_held_out_likelihood_beats_a_structureless_model():
    """Four states must explain held-out data better than one does."""
    rng = np.random.default_rng(7)
    rates, transition = ring_model(n_states=4, n_cells=8)
    train = simulate(2000, rates, transition, 0.25, rng)[1]
    test = simulate(800, rates, transition, 0.25, rng)[1]
    fits = [hmm.fit([train], 4, 0.25, seed=seed, n_iterations=100) for seed in range(3)]
    structured = max(fits, key=lambda m: m.log_likelihood)
    flat = hmm.fit([train], 1, 0.25, seed=0, n_iterations=100)
    assert (hmm.sequence_log_likelihood(structured, [test])
            > hmm.sequence_log_likelihood(flat, [test]))


def test_spectral_order_follows_a_chain_and_uses_no_behaviour():
    transition = np.zeros((5, 5))
    for i in range(5):
        for j in range(5):
            if abs(i - j) == 1:
                transition[i, j] = 0.4
        transition[i, i] = 1 - transition[i].sum()
    order = hmm.spectral_order(transition)
    forward = list(order) == sorted(range(5))
    backward = list(order) == sorted(range(5))[::-1]
    assert forward or backward, f"chain came back as {list(order)}"


def test_state_positions_average_a_ring_the_short_way():
    gamma = np.array([[1.0, 0.0], [1.0, 0.0], [0.0, 1.0]])
    positions = np.array([0.1, 2.8, 1.5])
    linear = hmm.state_positions(gamma, positions, 2, circular=False)
    ring = hmm.state_positions(gamma, positions, 2, circular=True, extent=2.9)
    assert linear[0] == pytest.approx(1.45)
    assert min(ring[0], 2.9 - ring[0]) < 0.2, "the ring mean took the long way round"


def test_circular_distance_and_correlation():
    assert hmm.circular_distance(0.1, 2.8, 2.9) == pytest.approx(0.2)
    # A distribution with a defined mean direction: a uniform ring has none,
    # and the correlation is numerically meaningless there.
    truth = np.linspace(0.2, 1.6, 200)
    assert hmm.circular_correlation(truth, truth, 2.9) == pytest.approx(1.0, abs=1e-6)
    assert hmm.circular_correlation(truth, (truth + 0.4) % 2.9, 2.9) == pytest.approx(
        1.0, abs=1e-6), "a constant offset on a ring is still a perfect correlation"
    assert hmm.circular_correlation(truth, (2.9 - truth) % 2.9, 2.9) == pytest.approx(
        -1.0, abs=1e-6), "a reversal should come back as -1"
    rng = np.random.default_rng(8)
    assert abs(hmm.circular_correlation(truth, rng.permutation(truth), 2.9)) < 0.3
