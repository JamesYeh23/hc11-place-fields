"""A Poisson hidden Markov model, fitted to binned population spiking.

Phase 4 asks whether a decoder that needs no place fields can find structure
where place fields ran out of resolution. Chen, Grosmark, Penagos & Wilson
(2016) call that decoder DecodewoRF: an *m*-state hidden Markov model fitted
directly to spike counts, with a Poisson observation model and a stationary
transition matrix, and no behavioural measure anywhere in the fit.

**Where this diverges from them, and why** (D10.6). They use a hierarchical
Dirichlet process prior and Markov chain Monte Carlo, which infers *m* as part
of the posterior. This is a fixed-*m* model fitted by expectation-maximisation,
with *m* chosen by held-out likelihood. The reason is not principle but
dependency risk: the published HDP-HMM implementations are unmaintained, and a
load-bearing library this project cannot test is worse than a simpler model it
can. A likelihood curve over *m* is also more legible than a number that falls
out of a sampler.

Three things in here exist because of how the project has been bitten before:

* sequences are fitted as **separate sequences**, never concatenated, so a
  transition can never span a gap in the searched time (the rule D6.2 set for
  detection, applied to fitting);
* the forward-backward recursion is scaled per bin and the per-bin maximum is
  factored out of the emissions before exponentiating, because 90-odd cells
  make raw emission likelihoods underflow to zero long before the arithmetic
  looks wrong;
* every fit reports the **whole set of restarts**, not the best one. A model
  whose five restarts disagree is a finding, not a fit to be quietly cleaned up.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

#: Floor on a fitted rate, in Hz. A state in which a cell never fired would
#: otherwise give log(0) and veto that state for any bin the cell spikes in --
#: the same reasoning as `sequences.rate_floor_hz`, and the same value.
RATE_FLOOR_HZ = 0.01


@dataclass(frozen=True)
class PoissonHMM:
    """A fitted model: initial distribution, transition matrix, rates per state.

    ``rates`` is ``(n_states, n_cells)`` in Hz. ``transition[i, j]`` is
    p(state j at t+1 | state i at t). ``log_likelihood`` is the training
    likelihood of the sequences it was fitted on, and ``converged`` says
    whether EM stopped because it had converged rather than because it hit the
    iteration cap -- which is worth keeping, since a model that ran out of
    iterations is not the model the likelihood curve claims.
    """

    initial: np.ndarray
    transition: np.ndarray
    rates: np.ndarray
    bin_s: float
    log_likelihood: float
    n_iterations: int
    converged: bool
    seed: int

    @property
    def n_states(self) -> int:
        return self.rates.shape[0]

    @property
    def n_cells(self) -> int:
        return self.rates.shape[1]


def log_emissions(counts: np.ndarray, rates: np.ndarray, bin_s: float) -> np.ndarray:
    """``(n_bins, n_states)`` log p(counts | state), up to a constant in state.

    The ``log(n!)`` term of the Poisson likelihood is dropped: it does not
    depend on the state, so it cancels in every posterior and shifts every
    model's likelihood by the same amount. Likelihoods from this module are
    therefore comparable to each other and not to a textbook figure.
    """
    counts = np.asarray(counts, dtype=float)
    expected = np.maximum(rates, RATE_FLOOR_HZ) * bin_s  # (n_states, n_cells)
    return counts @ np.log(expected).T - expected.sum(axis=1)[None, :]


def forward_backward(log_b: np.ndarray, initial: np.ndarray, transition: np.ndarray):
    """Scaled forward-backward. Returns ``(gamma, xi_sum, log_likelihood)``.

    ``gamma`` is ``(n_bins, n_states)`` posterior state occupancy, ``xi_sum``
    the ``(n_states, n_states)`` expected transition counts summed over the
    sequence.
    """
    n_bins, n_states = log_b.shape
    offset = log_b.max(axis=1)
    b = np.exp(log_b - offset[:, None])

    alpha = np.empty((n_bins, n_states))
    scale = np.empty(n_bins)
    alpha[0] = initial * b[0]
    scale[0] = alpha[0].sum()
    if scale[0] <= 0:
        raise FloatingPointError("the first bin has zero probability under every state")
    alpha[0] /= scale[0]
    for t in range(1, n_bins):
        alpha[t] = (alpha[t - 1] @ transition) * b[t]
        scale[t] = alpha[t].sum()
        if scale[t] <= 0:
            raise FloatingPointError(f"bin {t} has zero probability under every state")
        alpha[t] /= scale[t]

    beta = np.ones((n_bins, n_states))
    xi_sum = np.zeros((n_states, n_states))
    for t in range(n_bins - 2, -1, -1):
        weighted = b[t + 1] * beta[t + 1]
        beta[t] = (transition @ weighted) / scale[t + 1]
        xi_sum += (alpha[t][:, None] * transition * weighted[None, :]) / scale[t + 1]

    gamma = alpha * beta
    gamma /= gamma.sum(axis=1, keepdims=True)
    return gamma, xi_sum, float(np.sum(np.log(scale)) + offset.sum())


def _initialise(sequences, n_states, bin_s, rng):
    """k-means on the pooled population vectors, as a starting point for EM."""
    pooled = np.vstack(sequences)
    centres = pooled[rng.choice(len(pooled), size=n_states, replace=False)]
    for _ in range(25):
        distance = ((pooled[:, None, :] - centres[None, :, :]) ** 2).sum(axis=2)
        labels = distance.argmin(axis=1)
        moved = False
        for state in range(n_states):
            here = pooled[labels == state]
            if len(here):
                centre = here.mean(axis=0)
                moved |= not np.allclose(centre, centres[state])
                centres[state] = centre
        if not moved:
            break
    rates = np.maximum(centres / bin_s, RATE_FLOOR_HZ)
    transition = np.full((n_states, n_states), 1.0 / n_states)
    initial = np.full(n_states, 1.0 / n_states)
    return initial, transition, rates


def fit(sequences, n_states: int, bin_s: float, *, seed: int = 0,
        n_iterations: int = 200, tol: float = 1e-6) -> PoissonHMM:
    """Fit by expectation-maximisation over a list of independent sequences.

    ``sequences`` is a list of ``(n_bins, n_cells)`` count arrays. They are
    never concatenated: each contributes its own forward-backward pass, so no
    transition is ever counted across the gap between two of them.
    """
    sequences = [np.asarray(s, dtype=float) for s in sequences if len(s) > 1]
    if not sequences:
        raise ValueError("no sequence has more than one bin")
    rng = np.random.default_rng(seed)
    initial, transition, rates = _initialise(sequences, n_states, bin_s, rng)

    previous = -np.inf
    converged = False
    iteration = 0
    while iteration < n_iterations:
        iteration += 1
        start = np.zeros(n_states)
        trans = np.zeros((n_states, n_states))
        weight = np.zeros(n_states)
        weighted_counts = np.zeros_like(rates)
        total = 0.0
        for counts in sequences:
            gamma, xi_sum, log_likelihood = forward_backward(
                log_emissions(counts, rates, bin_s), initial, transition)
            start += gamma[0]
            trans += xi_sum
            weight += gamma.sum(axis=0)
            weighted_counts += gamma.T @ counts
            total += log_likelihood
        initial = start / start.sum()
        transition = trans / np.maximum(trans.sum(axis=1, keepdims=True), 1e-300)
        empty = trans.sum(axis=1) <= 0
        transition[empty] = 1.0 / n_states          # a state nothing enters
        rates = np.maximum(weighted_counts / np.maximum(weight, 1e-12)[:, None] / bin_s,
                           RATE_FLOOR_HZ)
        if np.isfinite(previous) and abs(total - previous) <= tol * abs(previous):
            previous = total
            converged = True
            break
        previous = total

    return PoissonHMM(initial=initial, transition=transition, rates=rates, bin_s=bin_s,
                      log_likelihood=float(previous), n_iterations=iteration,
                      converged=converged, seed=seed)


def posterior(model: PoissonHMM, counts: np.ndarray) -> np.ndarray:
    """``(n_bins, n_states)`` posterior state occupancy for one sequence."""
    gamma, _, _ = forward_backward(
        log_emissions(counts, model.rates, model.bin_s), model.initial, model.transition)
    return gamma


def memoryless_posterior(model: PoissonHMM, counts: np.ndarray,
                         bin_s: float | None = None) -> np.ndarray:
    """Per-bin posterior over states from the emissions alone, uniform prior.

    The model is fitted at the running bin (250 ms) but events are scored at
    20 ms, and a transition matrix estimated for 250 ms steps is simply wrong
    at 20 ms -- it would let the state move twelve times too fast. Chen et al.
    score with "no temporal prior", which is this: the state field supplies the
    likelihood and nothing supplies a dynamics. It is also exactly what the
    receptive-field decoder in `hc11.decode` does, so the two methods are
    scored the same way and differ only in where their fields came from.
    """
    log_b = log_emissions(counts, model.rates, model.bin_s if bin_s is None else bin_s)
    log_b -= log_b.max(axis=1, keepdims=True)
    weights = np.exp(log_b)
    return weights / weights.sum(axis=1, keepdims=True)


def sequence_log_likelihood(model: PoissonHMM, sequences) -> float:
    """Total log likelihood of held-out sequences under a fitted model."""
    total = 0.0
    for counts in sequences:
        counts = np.asarray(counts, dtype=float)
        if len(counts) < 1:
            continue
        _, _, value = forward_backward(
            log_emissions(counts, model.rates, model.bin_s),
            model.initial, model.transition)
        total += value
    return float(total)


def spectral_order(transition: np.ndarray) -> np.ndarray:
    """Order states along the transition structure, using no behaviour at all.

    The Fiedler vector of the graph Laplacian of the symmetrised transition
    matrix: states that pass probability to each other end up adjacent. This is
    the behaviour-free ordering Claim B needs, and the point of having it is
    that it can disagree with the behaviour-derived one.
    """
    weights = 0.5 * (transition + transition.T)
    np.fill_diagonal(weights, 0.0)
    degree = weights.sum(axis=1)
    laplacian = np.diag(degree) - weights
    values, vectors = np.linalg.eigh(laplacian)
    fiedler = vectors[:, np.argsort(values)[1]] if len(values) > 1 else vectors[:, 0]
    return np.argsort(fiedler)


def state_positions(gamma: np.ndarray, positions: np.ndarray, n_states: int,
                    *, circular: bool = False, extent: float = 1.0) -> np.ndarray:
    """Each state's position, as its posterior-weighted mean over training bins.

    This is the state-to-position map. **It uses behaviour**, and only for
    assessment: nothing in the fit has seen a position. On a ring the mean is
    circular, because averaging 0.1 m and 2.8 m on a 2.9 m track gives 1.45 m,
    which is the far side.
    """
    weight = gamma.sum(axis=0)
    if not circular:
        total = gamma.T @ positions
        return np.where(weight > 0, total / np.maximum(weight, 1e-12), np.nan)
    angle = 2 * np.pi * positions / extent
    cos = gamma.T @ np.cos(angle)
    sin = gamma.T @ np.sin(angle)
    mean = np.arctan2(sin, cos) % (2 * np.pi)
    return np.where(weight > 0, mean * extent / (2 * np.pi), np.nan)


def circular_distance(a: np.ndarray, b: np.ndarray, extent: float) -> np.ndarray:
    """Shortest distance between two positions on a ring of length ``extent``."""
    gap = np.abs(np.asarray(a) - np.asarray(b)) % extent
    return np.minimum(gap, extent - gap)


def circular_correlation(a: np.ndarray, b: np.ndarray, extent: float) -> float:
    """Jammalamadaka-SenGupta circular correlation, the ring's analogue of r.

    Chen et al. quote a Pearson correlation between the state trajectory and
    the run trajectory (0.73-0.79). Their tracks are circular, and Pearson on a
    circular variable depends on where the origin is cut, so this is reported
    in its place and labelled as what it is.
    """
    a = 2 * np.pi * np.asarray(a, dtype=float) / extent
    b = 2 * np.pi * np.asarray(b, dtype=float) / extent
    keep = np.isfinite(a) & np.isfinite(b)
    a, b = a[keep], b[keep]
    if a.size < 3:
        return float("nan")
    a_mean = np.arctan2(np.sin(a).mean(), np.cos(a).mean())
    b_mean = np.arctan2(np.sin(b).mean(), np.cos(b).mean())
    numerator = np.sum(np.sin(a - a_mean) * np.sin(b - b_mean))
    denominator = np.sqrt(np.sum(np.sin(a - a_mean) ** 2)
                          * np.sum(np.sin(b - b_mean) ** 2))
    return float(numerator / denominator) if denominator > 0 else float("nan")
