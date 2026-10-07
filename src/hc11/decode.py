r"""Memoryless Bayesian decoding of position from place-cell spike counts.

This is the object the paper's sleep analysis rests on: a spatial decoder built
from place-cell firing-rate vectors during track running, later applied to
candidate ripple events. Here it is built and validated on running data, where
the true position is known.

The model
---------
For a time bin of width :math:`\Delta` in which cell *i* fired :math:`n_i` spikes,
with tuning curve :math:`f_i(x)` in Hz,

.. math::

    P(x \mid \mathbf{n}) \;\propto\; P(x) \prod_i
        \frac{(f_i(x)\Delta)^{n_i}}{n_i!}\, e^{-f_i(x)\Delta}

and with a uniform prior, dropping the :math:`n_i!` terms (constant in *x*),

.. math::

    \log P(x \mid \mathbf{n}) \;=\; \sum_i n_i \log f_i(x)
        \;-\; \Delta \sum_i f_i(x) \;+\; \text{const}.

Everything is computed in log space and normalised with a log-sum-exp, so a bin
with many spikes cannot underflow the likelihood to zero.

Assumptions, and where they fail
--------------------------------
*Poisson firing.* Real place cells fire in complex-spike bursts, so the spike
count variance exceeds the mean. The decoder is therefore overconfident: the
posterior is sharper than the evidence warrants, and a single bursting cell can
pull the estimate toward its own field. This hurts the posterior's calibration
more than the peak position, which is why decoding *error* stays usable.

*Conditional independence across cells.* Place cells are not independent given
position; they are organised by theta phase and by assembly membership, and
pairs with overlapping fields co-fire more than independence predicts. Shared
variability is thus counted repeatedly as independent evidence, which again
sharpens the posterior more than it shifts the peak.

*Memorylessness.* No transition model links successive bins, so the decoder has
no notion that the animal cannot teleport. This is deliberate: it is what makes
the same decoder applicable to ripple events, where the "trajectory" is virtual
and compressed roughly 10-20x.

*Bin size.* At 250 ms on a running animal the assumptions are mild. At 20 ms
they are not: most cells contribute zero spikes, so each bin is decoded from a
handful of spikes, and within-bin theta sequences mean the animal's spikes
encode positions ahead of and behind its true location rather than at it. The
20 ms results here are reported to quantify that cost, not as a competing
estimate.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from scipy.special import logsumexp

from hc11.fields import RunningSamples, occupancy_map, rate_map, spike_map, spike_slots
from hc11.io import Session
from hc11.laps import Lap
from hc11.track import Track

# ---------------------------------------------------------------------------
# Core decoder
# ---------------------------------------------------------------------------


def decode_log_posterior(
    counts: np.ndarray,
    tuning: np.ndarray,
    dt: float,
    *,
    rate_floor_hz: float = 0.01,
    log_prior: np.ndarray | None = None,
) -> np.ndarray:
    """Log posterior over position for each time bin.

    Parameters
    ----------
    counts
        ``(n_cells, n_time)`` spike counts.
    tuning
        ``(n_cells, n_bins)`` firing rates in Hz. NaNs (unvisited bins) are
        treated as ``rate_floor_hz``.
    dt
        Time-bin width in seconds.
    rate_floor_hz
        Floor applied to the tuning curves, so that one spike from a cell with a
        zero-rate estimate at some position cannot veto that position outright.
    log_prior
        ``(n_bins,)``; uniform when omitted.

    Returns
    -------
    ``(n_bins, n_time)`` log posterior, normalised so each column sums to 1 in
    probability.
    """
    counts = np.asarray(counts, dtype=np.float64)
    tuning = np.asarray(tuning, dtype=np.float64)
    if counts.ndim != 2 or tuning.ndim != 2:
        raise ValueError("counts must be (n_cells, n_time) and tuning (n_cells, n_bins)")
    if counts.shape[0] != tuning.shape[0]:
        raise ValueError(
            f"{counts.shape[0]} cells in counts but {tuning.shape[0]} in tuning"
        )

    rates = np.where(np.isfinite(tuning), tuning, rate_floor_hz)
    rates = np.maximum(rates, rate_floor_hz)

    # (n_bins, n_time): sum_i n_i log f_i(x)  -  dt sum_i f_i(x)
    log_like = np.log(rates).T @ counts - dt * rates.sum(axis=0)[:, None]
    if log_prior is not None:
        log_like = log_like + np.asarray(log_prior, dtype=np.float64)[:, None]
    return log_like - logsumexp(log_like, axis=0, keepdims=True)


def posterior_mode(log_posterior: np.ndarray, bin_centers: np.ndarray) -> np.ndarray:
    """Maximum-a-posteriori position per time bin."""
    return np.asarray(bin_centers)[np.argmax(log_posterior, axis=0)]


# ---------------------------------------------------------------------------
# Building tuning curves and binned spike counts
# ---------------------------------------------------------------------------


def tuning_from_laps(
    session: Session,
    track: Track,
    samples: RunningSamples,
    cluster_ids: np.ndarray,
    lap_subset: np.ndarray,
    *,
    sigma_bins: float,
    min_occupancy_s: float,
    slots_by_cell: dict[int, np.ndarray] | None = None,
) -> np.ndarray:
    """``(n_cells, n_bins)`` tuning curves from a subset of this direction's laps.

    ``lap_subset`` is a boolean mask over ``samples`` selecting the training
    samples. Bins never visited in the training laps come back NaN and are
    floored by the decoder.
    """
    occupancy = occupancy_map(samples, lap_subset)
    out = np.empty((len(cluster_ids), track.n_bins), dtype=np.float64)
    for row, cluster_id in enumerate(cluster_ids):
        slots = (
            slots_by_cell[int(cluster_id)]
            if slots_by_cell is not None
            else spike_slots(session, int(cluster_id), samples)
        )
        counts = spike_map(slots, samples, lap_subset)
        out[row] = rate_map(
            counts, occupancy,
            sigma_bins=sigma_bins, circular=track.is_circular,
            min_occupancy_s=min_occupancy_s,
        )
    return out


def bin_lap(
    session: Session,
    samples: RunningSamples,
    lap: Lap,
    cluster_ids: np.ndarray,
    time_bin_s: float,
    *,
    slots_by_cell: dict[int, np.ndarray],
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Time-bin one lap: spike counts, true position, and bin centre times.

    Only samples that passed the running mask contribute a true position, so a
    decoded bin is always compared against a position the animal held while
    running. Bins with no such sample are dropped.
    """
    edges = np.arange(lap.t_start, lap.t_end + time_bin_s, time_bin_s)
    if len(edges) < 2:
        empty = np.empty((len(cluster_ids), 0))
        return empty, np.empty(0), np.empty(0)
    n_time = len(edges) - 1

    counts = np.zeros((len(cluster_ids), n_time), dtype=np.float64)
    for row, cluster_id in enumerate(cluster_ids):
        times = session.spikes_for(int(cluster_id))
        lo, hi = np.searchsorted(times, [edges[0], edges[-1]])
        if hi > lo:
            counts[row] = np.histogram(times[lo:hi], bins=edges)[0]

    # True position: mean of the running samples falling in each bin.
    t_run = session.position_t[samples.indices]
    x_run = session.position_1d[samples.indices]
    which = np.digitize(t_run, edges) - 1
    inside = (which >= 0) & (which < n_time)
    sums = np.bincount(which[inside], weights=x_run[inside], minlength=n_time)
    n_in = np.bincount(which[inside], minlength=n_time)
    with np.errstate(invalid="ignore"):
        true_pos = np.where(n_in > 0, sums / np.maximum(n_in, 1), np.nan)

    keep = n_in > 0
    centres = 0.5 * (edges[:-1] + edges[1:])
    return counts[:, keep], true_pos[keep], centres[keep]


# ---------------------------------------------------------------------------
# Leave-one-lap-out cross-validation
# ---------------------------------------------------------------------------


@dataclass
class DecodeResult:
    """Decoded positions and errors for one session, direction and bin size."""

    session: str
    direction: str
    time_bin_s: float
    cluster_ids: np.ndarray
    bin_centers: np.ndarray
    true_pos: np.ndarray  # (n_time,) m
    decoded_pos: np.ndarray  # (n_time,) m, NaN where undecodable
    error_cm: np.ndarray  # (n_time,) wrap-aware absolute error
    lap_of_bin: np.ndarray  # (n_time,) which held-out lap
    n_spikes_per_bin: np.ndarray
    n_laps: int
    circular: bool
    posteriors: dict = field(default_factory=dict, repr=False)  # lap index -> (n_bins, n_time)

    @property
    def n_cells(self) -> int:
        return len(self.cluster_ids)

    @property
    def decodable(self) -> np.ndarray:
        return np.isfinite(self.error_cm)

    @property
    def median_error_cm(self) -> float:
        e = self.error_cm[self.decodable]
        return float(np.median(e)) if e.size else float("nan")

    @property
    def mean_error_cm(self) -> float:
        e = self.error_cm[self.decodable]
        return float(np.mean(e)) if e.size else float("nan")

    @property
    def fraction_undecodable(self) -> float:
        n = self.error_cm.size
        return float(np.mean(~self.decodable)) if n else float("nan")

    def error_percentiles(self, qs=(10, 25, 50, 75, 90)) -> dict[int, float]:
        e = self.error_cm[self.decodable]
        return {q: (float(np.percentile(e, q)) if e.size else float("nan")) for q in qs}

    def error_by_position(self, n_groups: int = 8) -> tuple[np.ndarray, np.ndarray]:
        """Median error as a function of true position (are the track ends worse?)."""
        ok = self.decodable & np.isfinite(self.true_pos)
        if not ok.any():
            return np.empty(0), np.empty(0)
        lo, hi = self.bin_centers[0], self.bin_centers[-1]
        edges = np.linspace(lo, hi, n_groups + 1)
        idx = np.clip(np.digitize(self.true_pos[ok], edges) - 1, 0, n_groups - 1)
        centres = 0.5 * (edges[:-1] + edges[1:])
        med = np.array([
            np.median(self.error_cm[ok][idx == g]) if np.any(idx == g) else np.nan
            for g in range(n_groups)
        ])
        return centres, med

    def confusion_matrix(self, track: Track, normalize: bool = True) -> np.ndarray:
        """``(n_bins, n_bins)`` counts of true bin (rows) vs decoded bin (columns)."""
        ok = self.decodable & np.isfinite(self.true_pos)
        n = len(self.bin_centers)
        out = np.zeros((n, n), dtype=np.float64)
        if not ok.any():
            return out
        true_bin = track.digitize(self.true_pos[ok])
        dec_bin = track.digitize(self.decoded_pos[ok])
        np.add.at(out, (true_bin, dec_bin), 1.0)
        if normalize:
            with np.errstate(invalid="ignore"):
                out = out / np.maximum(out.sum(axis=1, keepdims=True), 1)
        return out


@dataclass
class LapData:
    """One lap's binned spike counts and true positions."""

    index: int
    direction: str
    counts: np.ndarray  # (n_cells, n_time)
    true_pos: np.ndarray  # (n_time,)
    centres: np.ndarray  # (n_time,) bin centre times


@dataclass
class DecoderInputs:
    """Everything the decoder needs, computed once per session/direction/bin size.

    ``tuning_by_lap[i]`` is the ``(n_cells, n_bins)`` template for decoding lap
    ``i``, built from every lap *except* ``i``. Precomputing these makes the
    chance baseline and the ensemble-size sweep cheap: both reuse the same
    templates, one by permuting their rows and the other by slicing them.
    """

    session: str
    direction: str
    time_bin_s: float
    cluster_ids: np.ndarray
    bin_centers: np.ndarray
    laps: list[LapData]
    tuning_by_lap: dict[int, np.ndarray]
    full_tuning: np.ndarray  # built from all laps; leaky, for diagnostics only
    circular: bool

    @property
    def n_cells(self) -> int:
        return len(self.cluster_ids)


def prepare_decoder(
    session: Session,
    track: Track,
    samples: RunningSamples,
    laps: list[Lap],
    cluster_ids: np.ndarray,
    params: dict,
    *,
    time_bin_s: float,
) -> DecoderInputs:
    """Bin every lap and build its leave-one-out template."""
    pf = params["place_fields"]
    cluster_ids = np.asarray(cluster_ids)
    slots_by_cell = {int(c): spike_slots(session, int(c), samples) for c in cluster_ids}

    wanted = (
        laps if samples.direction == "both"
        else [lap for lap in laps if lap.direction == samples.direction]
    )
    lap_data, tuning_by_lap = [], {}
    for lap in wanted:
        held_out = samples.lap_of_sample == lap.index
        if held_out.sum() < 2 or (~held_out).sum() < 2:
            continue
        counts, true_pos, centres = bin_lap(
            session, samples, lap, cluster_ids, time_bin_s, slots_by_cell=slots_by_cell
        )
        if counts.shape[1] == 0:
            continue
        lap_data.append(LapData(lap.index, lap.direction, counts, true_pos, centres))
        tuning_by_lap[lap.index] = tuning_from_laps(
            session, track, samples, cluster_ids, ~held_out,
            sigma_bins=pf["smoothing_sigma_bins"], min_occupancy_s=pf["min_occupancy_s"],
            slots_by_cell=slots_by_cell,
        )

    full_tuning = tuning_from_laps(
        session, track, samples, cluster_ids, np.ones(samples.n, dtype=bool),
        sigma_bins=pf["smoothing_sigma_bins"], min_occupancy_s=pf["min_occupancy_s"],
        slots_by_cell=slots_by_cell,
    )
    return DecoderInputs(
        session=session.name, direction=samples.direction, time_bin_s=time_bin_s,
        cluster_ids=cluster_ids, bin_centers=track.bin_centers, laps=lap_data,
        tuning_by_lap=tuning_by_lap, full_tuning=full_tuning, circular=track.is_circular,
    )


def decode(
    inputs: DecoderInputs,
    track: Track,
    params: dict,
    *,
    cell_subset: np.ndarray | None = None,
    tuning_override: np.ndarray | None = None,
    permutation: np.ndarray | None = None,
    keep_posteriors: bool = False,
) -> DecodeResult:
    """Decode every lap against its leave-one-out template.

    ``cell_subset`` restricts the ensemble (ensemble-size sweep);
    ``permutation`` reassigns tuning curves between cells (chance baseline);
    ``tuning_override`` replaces the cross-validated templates entirely and is
    for diagnostics and tests only -- passing ``inputs.full_tuning`` there is
    deliberately leaky.
    """
    dec = params["decoding"]
    rows = np.arange(inputs.n_cells) if cell_subset is None else np.asarray(cell_subset)
    cluster_ids = inputs.cluster_ids[rows]

    true_all, dec_all, lap_all, spk_all, posteriors = [], [], [], [], {}
    for lap in inputs.laps:
        if tuning_override is not None:
            tuning = tuning_override[rows]
        else:
            tuning = inputs.tuning_by_lap[lap.index]
            tuning = tuning[permutation][rows] if permutation is not None else tuning[rows]
        counts = lap.counts[rows]
        log_post = decode_log_posterior(
            counts, tuning, inputs.time_bin_s, rate_floor_hz=dec["rate_floor_hz"]
        )
        decoded = posterior_mode(log_post, inputs.bin_centers)
        total = counts.sum(axis=0)
        decoded[total < dec["min_spikes_per_bin"]] = np.nan
        true_all.append(lap.true_pos)
        dec_all.append(decoded)
        lap_all.append(np.full(len(lap.true_pos), lap.index))
        spk_all.append(total)
        if keep_posteriors:
            posteriors[lap.index] = log_post

    if not true_all:
        empty = np.empty(0)
        return DecodeResult(
            session=inputs.session, direction=inputs.direction, time_bin_s=inputs.time_bin_s,
            cluster_ids=cluster_ids, bin_centers=inputs.bin_centers, true_pos=empty,
            decoded_pos=empty, error_cm=empty, lap_of_bin=empty, n_spikes_per_bin=empty,
            n_laps=0, circular=inputs.circular,
        )

    true_pos = np.concatenate(true_all)
    decoded = np.concatenate(dec_all)
    return DecodeResult(
        session=inputs.session, direction=inputs.direction, time_bin_s=inputs.time_bin_s,
        cluster_ids=cluster_ids, bin_centers=inputs.bin_centers, true_pos=true_pos,
        decoded_pos=decoded,
        error_cm=np.abs(track.difference(decoded, true_pos)) * 100.0,
        lap_of_bin=np.concatenate(lap_all), n_spikes_per_bin=np.concatenate(spk_all),
        n_laps=len(inputs.laps), circular=inputs.circular, posteriors=posteriors,
    )


# ---------------------------------------------------------------------------
# Baselines and ensemble size
# ---------------------------------------------------------------------------


def chance_baseline(
    inputs: DecoderInputs, track: Track, params: dict, rng: np.random.Generator
) -> np.ndarray:
    """Median error when tuning curves are randomly reassigned between cells.

    Keeps both the spike counts and the set of fields intact and destroys only
    which field belongs to which cell, so it measures what the track geometry
    and firing statistics alone deliver. The permutation is applied to the same
    cross-validated templates, so this baseline is cross-validated too.
    """
    n_shuffles = params["decoding"]["n_chance_shuffles"]
    out = np.full(n_shuffles, np.nan)
    for i in range(n_shuffles):
        perm = rng.permutation(inputs.n_cells)
        out[i] = decode(inputs, track, params, permutation=perm).median_error_cm
    return out


def error_vs_ensemble_size(
    inputs: DecoderInputs, track: Track, params: dict, rng: np.random.Generator
) -> list[dict]:
    """Median decoding error for random subsets of the ensemble."""
    dec = params["decoding"]
    sizes = [s for s in dec["ensemble_sizes"] if s <= inputs.n_cells]
    if inputs.n_cells not in sizes:
        sizes.append(inputs.n_cells)
    rows = []
    for size in sizes:
        draws = 1 if size == inputs.n_cells else dec["n_ensemble_draws"]
        errors = []
        for _ in range(draws):
            subset = (
                np.arange(inputs.n_cells) if size == inputs.n_cells
                else rng.choice(inputs.n_cells, size=size, replace=False)
            )
            errors.append(decode(inputs, track, params, cell_subset=subset).median_error_cm)
        rows.append({
            "n_cells": size,
            "median_error_cm": float(np.nanmedian(errors)),
            "p25_cm": float(np.nanpercentile(errors, 25)),
            "p75_cm": float(np.nanpercentile(errors, 75)),
            "n_draws": draws,
        })
    return rows
