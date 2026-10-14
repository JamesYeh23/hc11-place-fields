"""Sequence scores for candidate events (step 7).

Each candidate event is binned at 20 ms, every bin is decoded to a posterior
over position with the step-4 place fields, and the resulting
``(time, position)`` matrix is scored by how well decoded position tracks bin
index. The score is normalised against within-event shuffles, which is what
makes it comparable between events of different length and different cell
counts.

**The weighted correlation.** Treat the posterior ``P(x, t)`` as a weight on
each ``(position, time)`` pair and take the ordinary Pearson correlation under
those weights::

    E[x]  = sum_{x,t} P(x,t) x      / sum P
    E[t]  = sum_{x,t} P(x,t) t      / sum P
    E[xt] = sum_{x,t} P(x,t) x t    / sum P
    cov   = E[xt] - E[x] E[t]
    r     = cov / sqrt( (E[x^2] - E[x]^2) (E[t^2] - E[t]^2) )

No position is read off the posterior: a bin whose posterior is spread over
half the track contributes its spread, not a point estimate. That matters here
because a 20 ms bin carries about one spike (step 5), so most single-bin
posteriors are broad and a maximum-a-posteriori trajectory would be mostly
noise. The signal is in the trajectory across bins.

**Empty bins are dropped.** A bin with no spikes still has a posterior -- the
memoryless Poisson decoder returns ``-dt * sum_i f_i(x)``, which favours
positions where the ensemble is quiet. That is a statement about the place
fields, not about the animal, so empty bins are excluded and the event is
required to retain a minimum number of them.

**Circular tracks.** Achilles_11012013 is a circular maze, where position has
no origin: a trajectory crossing the start of the linearisation would score
near zero despite being a perfect sequence. The correlation is therefore
evaluated at every circular offset of the position axis and the best taken --
"is there an origin that makes this trajectory linear in time", which is the
question that has meaning on a ring. The same maximisation is applied inside
the shuffles, so the normalised score is not inflated by it. On a linear track
there is one offset and the code path is identical.

**Field direction.** On a linear track the place fields are directional, so an
event is scored against each direction's templates and the larger ``|r|``
kept. That inflates the raw correlation and the inflation is not small. It is
harmless for the normalised score *only* because every shuffle goes through
the same maximisation: :func:`score_event` applies one shuffle realisation to
both directions' posteriors and takes the same maximum, so the null carries
the same bias as the observation. Raw ``r`` is reported too and should not be
compared against a nominal significance level.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from hc11.decode import decode_log_posterior

#: Shuffle kinds :func:`score_event` accepts.
SHUFFLES = ("column_cycle", "time_permutation")


# ---------------------------------------------------------------------------
# Binning and decoding one event
# ---------------------------------------------------------------------------


def bin_event(
    session, start_s: float, end_s: float, cluster_ids: np.ndarray, bin_s: float
) -> tuple[np.ndarray, np.ndarray]:
    """``(counts, centres)`` for one event: ``(n_cells, n_bins)`` and bin centres.

    Bin edges start at ``start_s``; a final partial bin is dropped rather than
    counted at a shorter width, which would give it a different expected count.
    """
    cluster_ids = np.asarray(cluster_ids)
    n_bins = int(np.floor((end_s - start_s) / bin_s))
    if n_bins <= 0 or cluster_ids.size == 0:
        return np.zeros((cluster_ids.size, 0)), np.empty(0)
    edges = start_s + np.arange(n_bins + 1) * bin_s
    lo, hi = np.searchsorted(session.spike_times, [edges[0], edges[-1]])
    times, ids = session.spike_times[lo:hi], session.spike_ids[lo:hi]
    counts = np.zeros((cluster_ids.size, n_bins))
    index = {int(c): i for i, c in enumerate(cluster_ids)}
    which = np.clip(np.searchsorted(edges, times, side="right") - 1, 0, n_bins - 1)
    for cluster, bin_index in zip(ids, which, strict=True):
        row = index.get(int(cluster))
        if row is not None:
            counts[row, bin_index] += 1
    return counts, edges[:-1] + bin_s / 2


@dataclass(frozen=True)
class EventPosteriors:
    """One event decoded against each available field direction."""

    posteriors: dict[str, np.ndarray]  # direction -> (n_kept_bins, n_position)
    bin_index: np.ndarray  # (n_kept_bins,) index of each kept bin within the event
    spikes_per_bin: np.ndarray  # (n_kept_bins,) spikes from the decoding ensemble
    n_bins_total: int  # before empty bins were dropped
    n_cells_active: int

    @property
    def n_bins(self) -> int:
        return len(self.bin_index)

    @property
    def empty_fraction(self) -> float:
        return 1.0 - self.n_bins / self.n_bins_total if self.n_bins_total else float("nan")


def event_posteriors(
    session, start_s: float, end_s: float, cluster_ids: np.ndarray,
    tuning_by_direction: dict[str, np.ndarray], bin_s: float,
    rate_floor_hz: float = 0.01,
) -> EventPosteriors:
    """Decode one event against each direction's templates, dropping empty bins."""
    counts, _ = bin_event(session, start_s, end_s, cluster_ids, bin_s)
    if counts.shape[1] == 0:
        return EventPosteriors({}, np.empty(0, dtype=int), np.empty(0), 0, 0)
    per_bin = counts.sum(axis=0)
    keep = per_bin > 0
    posteriors = {}
    for direction, tuning in tuning_by_direction.items():
        log_posterior = decode_log_posterior(
            counts[:, keep], tuning, bin_s, rate_floor_hz=rate_floor_hz
        )
        posteriors[direction] = np.exp(log_posterior).T  # (n_bins, n_position)
    return EventPosteriors(
        posteriors=posteriors,
        bin_index=np.flatnonzero(keep),
        spikes_per_bin=per_bin[keep],
        n_bins_total=counts.shape[1],
        n_cells_active=int((counts.sum(axis=1) > 0).sum()),
    )


# ---------------------------------------------------------------------------
# The weighted correlation
# ---------------------------------------------------------------------------


def offset_positions(positions: np.ndarray, circular: bool) -> np.ndarray:
    """``(n_position, n_offsets)`` position values under each circular offset.

    Column ``k`` is the position axis with its origin moved to bin ``k``. One
    column when ``circular`` is false, so callers have a single code path.
    """
    positions = np.asarray(positions, dtype=float)
    size = positions.size
    if not circular:
        return positions[:, None]
    step = float(np.median(np.diff(positions))) if size > 1 else 1.0
    index = np.arange(size)
    return ((index[:, None] - index[None, :]) % size) * step


def weighted_correlation(
    weights: np.ndarray, positions: np.ndarray, times: np.ndarray, circular: bool = False
) -> tuple[np.ndarray, np.ndarray]:
    """Weighted Pearson correlation of position against time, over the posterior.

    ``weights`` is ``(..., n_time, n_position)``; leading axes are independent
    replicates (shuffles). Returns ``(r, offset)`` with the leading shape of
    ``weights``: on a circular track ``r`` is the correlation at the offset that
    maximises ``|r|`` and ``offset`` is that offset's index, otherwise the
    single value and zero. The formula is in the module docstring.
    """
    weights = np.asarray(weights, dtype=float)
    times = np.asarray(times, dtype=float)
    grid = offset_positions(positions, circular)  # (n_position, n_offsets)

    total = weights.sum(axis=(-2, -1))
    by_position = weights.sum(axis=-2)  # (..., n_position)
    time_weighted = (weights * times[:, None]).sum(axis=-2)  # (..., n_position)
    sum_t = (weights * times[:, None]).sum(axis=(-2, -1))
    sum_tt = (weights * (times**2)[:, None]).sum(axis=(-2, -1))

    sum_x = by_position @ grid  # (..., n_offsets)
    sum_xx = by_position @ (grid**2)
    sum_xt = time_weighted @ grid

    with np.errstate(invalid="ignore", divide="ignore"):
        mean_x = sum_x / total[..., None]
        mean_t = sum_t / total
        var_x = sum_xx / total[..., None] - mean_x**2
        var_t = sum_tt / total - mean_t**2
        covariance = sum_xt / total[..., None] - mean_x * mean_t[..., None]
        r = covariance / np.sqrt(var_x * var_t[..., None])
    r = np.where(np.isfinite(r), r, 0.0)
    best = np.argmax(np.abs(r), axis=-1)
    return np.take_along_axis(r, best[..., None], axis=-1)[..., 0], best


# ---------------------------------------------------------------------------
# Shuffles
# ---------------------------------------------------------------------------


def shuffle_indices(kind: str, n_shuffles: int, n_bins: int, n_position: int, rng):
    """Index arrays describing ``n_shuffles`` realisations of one shuffle kind.

    Returned as indices rather than shuffled data so that the *same* realisation
    can be applied to each field direction's posterior. Applying independent
    shuffles per direction would give the maximisation over directions more
    freedom in the null than the observation had, and the normalised score would
    come out too low.
    """
    if kind == "column_cycle":
        shifts = rng.integers(0, n_position, size=(n_shuffles, n_bins))
        columns = (np.arange(n_position)[None, None, :] - shifts[:, :, None]) % n_position
        return ("column_cycle", columns)
    if kind == "time_permutation":
        order = np.argsort(rng.random((n_shuffles, n_bins)), axis=1)
        return ("time_permutation", order)
    raise ValueError(f"unknown shuffle {kind!r}; expected one of {SHUFFLES}")


def apply_shuffle(posterior: np.ndarray, prepared) -> np.ndarray:
    """``(n_shuffles, n_time, n_position)`` from one ``(n_time, n_position)`` posterior."""
    kind, index = prepared
    if kind == "column_cycle":
        return np.take_along_axis(
            np.broadcast_to(posterior, index.shape), index, axis=2
        )
    return posterior[index]


# ---------------------------------------------------------------------------
# Scoring
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class SequenceScore:
    """One event's score under one shuffle kind."""

    r: float  # signed weighted correlation, at the chosen direction and offset
    rz: float  # |r| standardised against the shuffle distribution of |r|
    direction: str  # field direction that gave the larger |r|
    offset: int  # circular offset chosen (0 on a linear track)
    null_mean: float
    null_sd: float
    n_shuffles: int
    shuffle: str
    r_by_direction: dict[str, float] = field(default_factory=dict)

    @property
    def forward(self) -> bool:
        """Decoded position increasing with time."""
        return self.r > 0


def score_event(
    event: EventPosteriors,
    positions: np.ndarray,
    bin_s: float,
    *,
    circular: bool,
    shuffle: str,
    n_shuffles: int,
    rng,
) -> SequenceScore | None:
    """Score one decoded event, normalised against ``n_shuffles`` shuffles.

    Returns ``None`` when the event has no posterior to score. The maximisation
    over field direction and, on a circular track, over position offset is
    applied identically to the observation and to every shuffle.
    """
    if not event.posteriors or event.n_bins == 0:
        return None
    times = event.bin_index * bin_s
    observed = {}
    for direction, posterior in event.posteriors.items():
        r, offset = weighted_correlation(posterior, positions, times, circular)
        observed[direction] = (float(r), int(offset))
    best_direction = max(observed, key=lambda d: abs(observed[d][0]))
    best_r, best_offset = observed[best_direction]

    prepared = shuffle_indices(
        shuffle, n_shuffles, event.n_bins, len(positions), rng
    )
    null = None
    for posterior in event.posteriors.values():
        r, _ = weighted_correlation(
            apply_shuffle(posterior, prepared), positions, times, circular
        )
        null = np.abs(r) if null is None else np.maximum(null, np.abs(r))

    mean, sd = float(null.mean()), float(null.std())
    rz = (abs(best_r) - mean) / sd if sd > 0 else float("nan")
    return SequenceScore(
        r=best_r, rz=rz, direction=best_direction, offset=best_offset,
        null_mean=mean, null_sd=sd, n_shuffles=n_shuffles, shuffle=shuffle,
        r_by_direction={d: value for d, (value, _) in observed.items()},
    )
