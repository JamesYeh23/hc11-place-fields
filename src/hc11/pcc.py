"""Per-cell contribution (PCC): how much one cell's place field carries a sequence.

Grosmark & Buzsáki define PCC as "the neuron's mean contribution across all PRE
and POST putative events as determined by a cell-specific shuffling technique".
The technique itself is in the Science supplementary materials, which this
project does not have. **The one procedural constraint we do have comes from
the caption of their Fig. 2D** — "templates in which an increasing percentage of
neuron's place fields were shuffled" — so the perturbation is applied to the
decoding template, not to the spike train. Everything else here is our choice,
recorded in D9.3, and a replication resting on a figure caption should say so.

The estimator, in one line: shift one cell's place field around the track,
re-score every event with that altered template, and ask how much the score
falls. Concretely, for cell *i* and epoch *E*:

    PCC(i, E) = mean_e rZ(e) − mean_s mean_e rZ(e | field i shifted by s)

The shifts are **every** non-zero circular offset of the position axis, so the
estimator is deterministic — no seed, no sampling noise — and the same shifts
provide the null: if the cell's real field were just another placement, PCC
would sit inside the spread of what the other placements give.

That spread is also what makes the **null label** possible at no extra cost. Take
one shifted field as though it were the cell's real one, leave the rest as its
null, and the identical statistic comes out of data in which the cell's field
carries no information about position at all. Comparing the distribution of real
ΔPCC with the distribution of those null ΔPCC values says directly whether the
rigid/plastic split can be told apart from what the estimator produces from
noise. That comparison is the headline of step 9c, not a check on it.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from hc11.decode import decode_log_posterior
from hc11.sequences import weighted_correlation


@dataclass(frozen=True)
class PreparedEvent:
    """One event's binned spikes and the step-7 numbers needed to re-score it.

    ``counts`` has empty time bins already dropped, matching what
    :func:`hc11.sequences.event_posteriors` scores. ``null_mean`` and
    ``null_sd`` are the within-event shuffle statistics computed in step 7 and
    reused here rather than recomputed per perturbation (D9.4): that reuse is
    the single approximation in the estimator, and it is validated against fresh
    nulls before being relied on.
    """

    peak: float
    counts: np.ndarray
    times: np.ndarray
    present: np.ndarray
    observed_abs_r: float
    null_mean: dict
    null_sd: dict

    @property
    def n_cells(self) -> int:
        return self.counts.shape[0]


def shift_field(tuning: dict, row: int, offset: int) -> dict:
    """One cell's field moved around the track, by the same offset in every direction.

    A circular shift keeps the cell's peak rate, field width and shape and
    destroys only *where* the field is, so a contribution measured against it is
    a contribution of position rather than of excitability. The same offset is
    applied in both running directions because a cell is one cell: shifting its
    two directional fields independently would also destroy the relationship
    between them, which is a second change on top of the one being tested.
    """
    if offset == 0:
        return tuning
    out = {}
    for direction, values in tuning.items():
        moved = values.copy()
        moved[row] = np.roll(values[row], offset)
        out[direction] = moved
    return out


def best_abs_r(
    counts: np.ndarray, tuning: dict, positions: np.ndarray, times: np.ndarray,
    bin_s: float, *, circular: bool, rate_floor_hz: float,
) -> float:
    """Largest |weighted correlation| over running directions, as step 7 scores it.

    The direction (and, on a ring, the position offset) is maximised over in the
    same way for the real template and for every perturbed one, so a perturbation
    cannot win or lose by changing which direction is chosen.
    """
    best = 0.0
    for values in tuning.values():
        posterior = np.exp(
            decode_log_posterior(counts, values, bin_s, rate_floor_hz=rate_floor_hz)
        ).T
        r, _ = weighted_correlation(posterior, positions, times, circular)
        best = max(best, abs(float(r)))
    return best


def standardise(abs_r: float, null_mean: float, null_sd: float) -> float:
    """rZ from an already-computed null. NaN when the null has no width."""
    if not np.isfinite(null_sd) or null_sd <= 0:
        return float("nan")
    return (abs_r - null_mean) / null_sd


def finite_everywhere(stack: np.ndarray) -> np.ndarray:
    """Events usable at *every* placement of the field, as a mask over events.

    ``stack`` is ``(n_placements, n_events)``. A placement that leaves one event
    unscorable — a degenerate posterior gives a zero-width null and no z-score —
    must not be compared against placements that scored it, or the contribution
    becomes a difference between two event sets rather than between two field
    positions. One mask, applied to all placements including the real one.
    """
    stack = np.atleast_2d(np.asarray(stack, dtype=float))
    return np.all(np.isfinite(stack), axis=0)


def contribution(observed_mean: float, shift_means: np.ndarray) -> dict:
    """PCC and its within-cell significance, from the mean score at every shift.

    ``shift_means[s]`` is the epoch's mean rZ with this cell's field placed at
    offset *s*. PCC is how far the real placement sits above the average
    placement; the spread across placements is the null, so the test is exact
    and its finest p is ``1 / (n_shifts + 1)`` — on a 16-bin track that is 0.06,
    which is why the standardised ``z`` rather than ``p`` decides which cells
    count as strongly contributing (D9.3).
    """
    shift_means = np.asarray(shift_means, dtype=float)
    usable = shift_means[np.isfinite(shift_means)]
    if not np.isfinite(observed_mean) or usable.size < 2:
        return {"pcc": float("nan"), "z": float("nan"), "p": float("nan"),
                "n_shifts": int(usable.size)}
    centre = float(usable.mean())
    spread = float(usable.std(ddof=1))
    pcc = observed_mean - centre
    deviations = np.abs(usable - centre)
    return {
        "pcc": float(pcc),
        "z": float(pcc / spread) if spread > 0 else float("nan"),
        "p": float((1 + int(np.sum(deviations >= abs(pcc)))) / (usable.size + 1)),
        "n_shifts": int(usable.size),
    }


def null_contributions(shift_means: np.ndarray) -> np.ndarray:
    """The same statistic with each shifted field standing in for the real one.

    For shift *s*: ``shift_means[s] − mean(shift_means[s' != s])``. These are the
    PCC values the estimator produces when the field it is measuring carries no
    information about position — the null labels, one per shift, at no extra
    computational cost because every term is already in hand.
    """
    values = np.asarray(shift_means, dtype=float)
    finite = np.isfinite(values)
    out = np.full(values.shape, np.nan)
    if finite.sum() < 2:
        return out
    total = values[finite].sum()
    n = int(finite.sum())
    out[finite] = values[finite] - (total - values[finite]) / (n - 1)
    return out
