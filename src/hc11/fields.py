"""One-dimensional place fields: rate maps, per-cell metrics, shuffle test.

Pipeline for one session and one running direction:

1. Select samples: inside a lap of that direction (:mod:`hc11.laps`) **and** in
   the running mask (:mod:`hc11.behavior`). Laps are segmented before the speed
   filter; the speed filter enters here (decision D4.1).
2. Accumulate occupancy and spike counts per spatial bin.
3. Smooth occupancy and spike counts **separately**, then divide. Smoothing the
   ratio instead would let a single spike in a briefly visited bin dominate its
   neighbourhood. Circular tracks smooth with wraparound.
4. Score each cell: Skaggs spatial information, peak and mean rate, field count
   and width, split-half stability, laps active.
5. Test spatial information against a circular-shift null (:func:`shuffle_test`).

Spikes are assigned to the nearest position sample rather than to an interpolated
position. At 39.06 Hz and the observed running speeds this places a spike within
~1.5 cm of its interpolated position, an order of magnitude finer than the 10 cm
bins, and it keeps every spike on the same sample grid the shuffle operates on.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from scipy.ndimage import gaussian_filter1d

from hc11.io import Session
from hc11.laps import Lap, lap_mask
from hc11.track import Track

# ---------------------------------------------------------------------------
# Sample selection
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class RunningSamples:
    """The position samples that feed one direction's rate maps.

    ``slot`` numbers the selected samples 0..n-1 in time order: this is the
    "compressed running time" the circular-shift shuffle operates on.
    """

    session: str
    direction: str
    indices: np.ndarray  # (n,) indices into session.position_t
    bins: np.ndarray  # (n,) spatial bin of each selected sample
    lap_of_sample: np.ndarray  # (n,) lap index of each selected sample
    slot_of_index: np.ndarray  # (M,) slot per position sample, -1 if unselected
    dt: float
    n_bins: int

    @property
    def n(self) -> int:
        return len(self.indices)

    @property
    def seconds(self) -> float:
        return self.n * self.dt


def select_running_samples(
    session: Session,
    track: Track,
    laps: list[Lap],
    speed_mask: np.ndarray,
    direction: str,
) -> RunningSamples:
    """Samples inside a lap of ``direction`` that also pass the speed filter."""
    in_lap = lap_mask(session, laps, direction)
    bins_all = track.digitize(session.position_1d)
    keep = in_lap & np.asarray(speed_mask, dtype=bool) & (bins_all >= 0)
    indices = np.flatnonzero(keep)

    lap_of_sample = np.full(len(session.position_t), -1, dtype=np.int64)
    for lap in laps:
        if lap.direction == direction:
            lap_of_sample[lap.start_sample : lap.stop_sample] = lap.index

    slot_of_index = np.full(len(session.position_t), -1, dtype=np.int64)
    slot_of_index[indices] = np.arange(len(indices))

    return RunningSamples(
        session=session.name,
        direction=direction,
        indices=indices,
        bins=bins_all[indices],
        lap_of_sample=lap_of_sample[indices],
        slot_of_index=slot_of_index,
        dt=session.position_dt,
        n_bins=track.n_bins,
    )


def spike_slots(session: Session, cluster_id: int, samples: RunningSamples) -> np.ndarray:
    """Slot index of each of one cluster's spikes that lands on a selected sample.

    Each spike is attached to the nearest position sample; spikes falling on
    unselected samples (not running, or not in a lap of this direction) are
    dropped.
    """
    t_pos = session.position_t
    times = session.spikes_for(cluster_id)
    if times.size == 0 or t_pos.size == 0:
        return np.empty(0, dtype=np.int64)
    right = np.searchsorted(t_pos, times)
    left = np.clip(right - 1, 0, len(t_pos) - 1)
    right = np.clip(right, 0, len(t_pos) - 1)
    nearest = np.where(
        np.abs(times - t_pos[left]) <= np.abs(t_pos[right] - times), left, right
    )
    slots = samples.slot_of_index[nearest]
    return slots[slots >= 0]


# ---------------------------------------------------------------------------
# Rate maps
# ---------------------------------------------------------------------------


def smooth_1d(
    a: np.ndarray, sigma_bins: float, circular: bool, truncate: float = 4.0
) -> np.ndarray:
    """Gaussian smoothing along the track; wraps on a circular track.

    On a linear track the array is padded with zeros (``mode="constant"``), which
    is the right choice here because occupancy and spike counts are smoothed with
    the same kernel and then divided: both are attenuated identically at the
    edges, so the ratio stays unbiased.
    """
    if sigma_bins <= 0:
        return np.asarray(a, dtype=np.float64).copy()
    mode = "wrap" if circular else "constant"
    return gaussian_filter1d(
        np.asarray(a, dtype=np.float64), sigma_bins, mode=mode, cval=0.0, truncate=truncate
    )


@dataclass(frozen=True)
class RateMap:
    """One cell's directional rate map, with the pieces it was built from."""

    cluster_id: int
    direction: str
    bin_centers: np.ndarray  # (n_bins,) m
    occupancy_s: np.ndarray  # (n_bins,) raw seconds per bin
    spike_counts: np.ndarray  # (n_bins,) raw spike counts
    rate: np.ndarray  # (n_bins,) Hz, smoothed; NaN where under-occupied
    circular: bool
    bin_size_m: float
    n_spikes: int

    @property
    def valid(self) -> np.ndarray:
        return ~np.isnan(self.rate)

    @property
    def mean_rate(self) -> float:
        """Occupancy-weighted mean rate over visited bins (Hz)."""
        v = self.valid
        if not v.any() or self.occupancy_s[v].sum() == 0:
            return float("nan")
        return float(np.sum(self.rate[v] * self.occupancy_s[v]) / self.occupancy_s[v].sum())

    @property
    def peak_rate(self) -> float:
        v = self.valid
        return float(np.nanmax(self.rate[v])) if v.any() else float("nan")

    @property
    def peak_position(self) -> float:
        v = self.valid
        if not v.any():
            return float("nan")
        return float(self.bin_centers[v][np.argmax(self.rate[v])])


def occupancy_map(samples: RunningSamples, lap_subset: np.ndarray | None = None) -> np.ndarray:
    """Seconds spent in each spatial bin."""
    bins = samples.bins if lap_subset is None else samples.bins[lap_subset]
    return np.bincount(bins, minlength=samples.n_bins).astype(np.float64) * samples.dt


def spike_map(
    slots: np.ndarray, samples: RunningSamples, lap_subset: np.ndarray | None = None
) -> np.ndarray:
    """Spike count per spatial bin, from slot indices."""
    if slots.size == 0:
        return np.zeros(samples.n_bins, dtype=np.float64)
    keep = slots if lap_subset is None else slots[lap_subset[slots]]
    if keep.size == 0:
        return np.zeros(samples.n_bins, dtype=np.float64)
    return np.bincount(samples.bins[keep], minlength=samples.n_bins).astype(np.float64)


def rate_map(
    counts: np.ndarray,
    occupancy: np.ndarray,
    *,
    sigma_bins: float,
    circular: bool,
    min_occupancy_s: float,
    truncate: float = 4.0,
) -> np.ndarray:
    """Smoothed firing rate per bin (Hz); NaN where raw occupancy is too low.

    Occupancy and counts are smoothed separately and then divided -- never the
    ratio.
    """
    num = smooth_1d(counts, sigma_bins, circular, truncate)
    den = smooth_1d(occupancy, sigma_bins, circular, truncate)
    with np.errstate(invalid="ignore", divide="ignore"):
        rate = np.where(den > 0, num / den, np.nan)
    rate[occupancy < min_occupancy_s] = np.nan
    return rate


def build_rate_map(
    session: Session,
    track: Track,
    samples: RunningSamples,
    cluster_id: int,
    *,
    sigma_bins: float,
    min_occupancy_s: float,
    slots: np.ndarray | None = None,
) -> RateMap:
    slots = spike_slots(session, cluster_id, samples) if slots is None else slots
    occ = occupancy_map(samples)
    counts = spike_map(slots, samples)
    return RateMap(
        cluster_id=int(cluster_id),
        direction=samples.direction,
        bin_centers=track.bin_centers,
        occupancy_s=occ,
        spike_counts=counts,
        rate=rate_map(
            counts, occ,
            sigma_bins=sigma_bins, circular=track.is_circular, min_occupancy_s=min_occupancy_s,
        ),
        circular=track.is_circular,
        bin_size_m=track.actual_bin_size,
        n_spikes=int(slots.size),
    )


# ---------------------------------------------------------------------------
# Metrics
# ---------------------------------------------------------------------------


def skaggs_information(rate: np.ndarray, occupancy_s: np.ndarray) -> tuple[float, float]:
    r"""Skaggs spatial information, in bits/spike and bits/second.

    With :math:`p_i` the fraction of time spent in bin :math:`i`, :math:`\lambda_i`
    the firing rate there, and :math:`\bar\lambda = \sum_i p_i \lambda_i` the
    overall mean rate:

    .. math::

        I_{\text{spike}} = \sum_i p_i \frac{\lambda_i}{\bar\lambda}
                           \log_2 \frac{\lambda_i}{\bar\lambda}
        \qquad\text{(bits/spike)}

        I_{\text{sec}} = \sum_i p_i \lambda_i \log_2 \frac{\lambda_i}{\bar\lambda}
                       = \bar\lambda \, I_{\text{spike}}
        \qquad\text{(bits/second)}

    Bins with zero rate contribute zero (the limit of :math:`x \log x` as
    :math:`x \to 0`). Bins that are NaN in ``rate`` (under-occupied) are excluded
    and :math:`p_i` is renormalised over the rest. Returns ``(nan, nan)`` if the
    cell never fires or nothing was visited.

    Reference: Skaggs, McNaughton, Gothard & Markus (1993), *NIPS* 5, 1030-1037.
    """
    rate = np.asarray(rate, dtype=np.float64)
    occ = np.asarray(occupancy_s, dtype=np.float64)
    valid = np.isfinite(rate) & (occ > 0)
    if not valid.any():
        return float("nan"), float("nan")
    p = occ[valid] / occ[valid].sum()
    lam = rate[valid]
    mean_rate = float(np.sum(p * lam))
    if mean_rate <= 0:
        return float("nan"), float("nan")
    ratio = lam / mean_rate
    with np.errstate(divide="ignore", invalid="ignore"):
        terms = np.where(lam > 0, p * ratio * np.log2(ratio), 0.0)
    bits_per_spike = float(np.sum(terms))
    return bits_per_spike, bits_per_spike * mean_rate


def _contiguous_runs(mask: np.ndarray, circular: bool) -> list[tuple[int, int]]:
    """Contiguous ``True`` runs as ``(start, length)``; wraps if ``circular``."""
    if not mask.any():
        return []
    if mask.all():
        return [(0, len(mask))]
    padded = np.concatenate([[False], mask, [False]])
    edges = np.flatnonzero(np.diff(padded.astype(np.int8)))
    runs = [(int(s), int(e - s)) for s, e in zip(edges[0::2], edges[1::2], strict=True)]
    if circular and len(runs) > 1 and mask[0] and mask[-1]:
        (_, l0), (sl, ll) = runs[0], runs[-1]
        runs = runs[1:-1] + [(sl, ll + l0)]  # merge the run across the wrap point
    return runs


def find_fields(
    rate: np.ndarray,
    bin_size_m: float,
    *,
    threshold_frac: float,
    min_width_m: float,
    max_width_m: float,
    circular: bool,
) -> list[tuple[int, int]]:
    """Place fields as contiguous runs above ``threshold_frac`` of the peak rate."""
    valid = np.isfinite(rate)
    if not valid.any():
        return []
    peak = np.nanmax(rate)
    if not np.isfinite(peak) or peak <= 0:
        return []
    above = np.where(valid, rate >= threshold_frac * peak, False)
    fields = []
    for start, length in _contiguous_runs(above, circular):
        width = length * bin_size_m
        if min_width_m <= width <= max_width_m:
            fields.append((start, length))
    return fields


def split_half_stability(
    session: Session,
    track: Track,
    samples: RunningSamples,
    cluster_id: int,
    slots: np.ndarray,
    *,
    sigma_bins: float,
    min_occupancy_s: float,
) -> float:
    """Pearson correlation between rate maps from odd- and even-numbered laps.

    ``nan`` when either half has no valid bins or a constant map (e.g. a cell that
    fires in only one half).
    """
    lap_ids = np.unique(samples.lap_of_sample)
    lap_ids = lap_ids[lap_ids >= 0]
    if lap_ids.size < 2:
        return float("nan")
    order = {lap: i for i, lap in enumerate(lap_ids)}
    parity = np.array([order[lap] % 2 for lap in samples.lap_of_sample])

    maps = []
    for half in (0, 1):
        subset = parity == half
        occ = occupancy_map(samples, subset)
        counts = spike_map(slots, samples, subset)
        maps.append(rate_map(
            counts, occ, sigma_bins=sigma_bins, circular=track.is_circular,
            min_occupancy_s=min_occupancy_s,
        ))
    both = np.isfinite(maps[0]) & np.isfinite(maps[1])
    if both.sum() < 3:
        return float("nan")
    a, b = maps[0][both], maps[1][both]
    if np.std(a) == 0 or np.std(b) == 0:
        return float("nan")
    return float(np.corrcoef(a, b)[0, 1])


def laps_active(samples: RunningSamples, slots: np.ndarray) -> int:
    """Number of laps in which the cell fired at least one spike."""
    if slots.size == 0:
        return 0
    return int(len(np.unique(samples.lap_of_sample[slots])))


# ---------------------------------------------------------------------------
# Shuffle test
# ---------------------------------------------------------------------------


def shuffle_test(
    samples: RunningSamples,
    slots: np.ndarray,
    *,
    sigma_bins: float,
    min_occupancy_s: float,
    circular_track: bool,
    n_shuffles: int,
    rng: np.random.Generator,
    min_shift_slots: int = 1,
) -> tuple[np.ndarray, float]:
    """Null distribution of Skaggs information under circular shifts.

    The spike train is shifted circularly **in compressed running time**: the
    selected running samples are numbered 0..n-1 (``slot``), each spike sits in a
    slot, and a shuffle rolls the per-slot spike counts by a random offset. This
    preserves each cell's spike count, its firing statistics, and the animal's
    trajectory, while destroying the pairing between the two -- and, because the
    shift happens inside the running samples only, every shuffled spike still
    lands at a position the animal actually occupied while running.

    Returns the null information values (bits/spike) and their mean occupancy,
    which the caller compares against the observed value.
    """
    occ = occupancy_map(samples)
    n_slots = samples.n
    null = np.full(n_shuffles, np.nan)
    if slots.size == 0 or n_slots < 2:
        return null, float("nan")

    per_slot = np.bincount(slots, minlength=n_slots).astype(np.float64)
    max_shift = n_slots - min_shift_slots
    shifts = rng.integers(min_shift_slots, max(max_shift, min_shift_slots + 1), size=n_shuffles)
    for i, shift in enumerate(shifts):
        rolled = np.roll(per_slot, shift)
        counts = np.bincount(samples.bins, weights=rolled, minlength=samples.n_bins)
        r = rate_map(
            counts, occ, sigma_bins=sigma_bins, circular=circular_track,
            min_occupancy_s=min_occupancy_s,
        )
        null[i], _ = skaggs_information(r, occ)
    return null, float(np.nanmean(null))


# ---------------------------------------------------------------------------
# Per-cell analysis
# ---------------------------------------------------------------------------


@dataclass
class CellResult:
    """Everything computed for one cell in one direction."""

    session: str
    animal: str
    maze_type: str
    direction: str
    cluster_id: int
    shank: int
    n_spikes: int
    mean_rate: float
    peak_rate: float
    peak_position_m: float
    info_bits_per_spike: float
    info_bits_per_second: float
    info_null_mean: float
    info_z: float
    info_p: float
    n_fields: int
    field_width_m: float  # total width of all fields
    largest_field_width_m: float
    in_field_mean_rate: float
    out_field_mean_rate: float
    stability: float
    n_laps: int
    n_laps_active: int
    is_place_cell: bool
    is_place_cell_no_stability: bool
    rate_map: RateMap = field(repr=False, default=None)


def analyse_cell(
    session: Session,
    track: Track,
    samples: RunningSamples,
    cluster_id: int,
    params: dict,
    rng: np.random.Generator,
) -> CellResult:
    """Rate map, metrics, shuffle test and classification for one cell."""
    pf = params["place_fields"]
    crit = pf["criteria"]
    sigma = pf["smoothing_sigma_bins"]
    min_occ = pf["min_occupancy_s"]
    n_shuffles = pf["shuffle"]["n_shuffles"]
    alpha = pf["shuffle"]["alpha"]

    slots = spike_slots(session, cluster_id, samples)
    rmap = build_rate_map(
        session, track, samples, cluster_id,
        sigma_bins=sigma, min_occupancy_s=min_occ, slots=slots,
    )
    info_spike, info_sec = skaggs_information(rmap.rate, rmap.occupancy_s)

    null, null_mean = shuffle_test(
        samples, slots,
        sigma_bins=sigma, min_occupancy_s=min_occ, circular_track=track.is_circular,
        n_shuffles=n_shuffles, rng=rng,
    )
    finite_null = null[np.isfinite(null)]
    if finite_null.size and np.isfinite(info_spike):
        # +1 in numerator and denominator: a Monte-Carlo p-value can never be 0.
        p = (np.sum(finite_null >= info_spike) + 1) / (finite_null.size + 1)
        sd = float(np.std(finite_null))
        z = (info_spike - float(np.mean(finite_null))) / sd if sd > 0 else float("nan")
    else:
        p, z = float("nan"), float("nan")

    fields = find_fields(
        rmap.rate, rmap.bin_size_m,
        threshold_frac=crit["field_threshold_frac_of_peak"],
        min_width_m=crit["min_field_width_cm"] / 100.0,
        max_width_m=crit["max_field_width_cm"] / 100.0,
        circular=track.is_circular,
    )
    in_field = np.zeros(track.n_bins, dtype=bool)
    for start, length in fields:
        in_field[(np.arange(start, start + length)) % track.n_bins] = True
    valid = rmap.valid
    in_rate = float(np.nanmean(rmap.rate[in_field & valid])) if (in_field & valid).any() else np.nan
    out_rate = (
        float(np.nanmean(rmap.rate[~in_field & valid])) if (~in_field & valid).any() else np.nan
    )

    stability = split_half_stability(
        session, track, samples, cluster_id, slots, sigma_bins=sigma, min_occupancy_s=min_occ,
    )

    passes_core = bool(
        rmap.n_spikes >= crit["min_spikes"]
        and np.isfinite(rmap.peak_rate)
        and rmap.peak_rate >= crit["min_peak_rate_hz"]
        and np.isfinite(p)
        and p < alpha
        and len(fields) >= 1
    )
    passes_stability = bool(np.isfinite(stability) and stability >= crit["min_stability_r"])

    n_laps = int(len(np.unique(samples.lap_of_sample[samples.lap_of_sample >= 0])))
    return CellResult(
        session=session.name,
        animal=session.animal,
        maze_type=session.maze_type,
        direction=samples.direction,
        cluster_id=int(cluster_id),
        shank=int(session.shank_of_cluster[int(cluster_id)]),
        n_spikes=rmap.n_spikes,
        mean_rate=rmap.mean_rate,
        peak_rate=rmap.peak_rate,
        peak_position_m=rmap.peak_position,
        info_bits_per_spike=info_spike,
        info_bits_per_second=info_sec,
        info_null_mean=null_mean,
        info_z=z,
        info_p=p,
        n_fields=len(fields),
        field_width_m=float(sum(length for _, length in fields) * rmap.bin_size_m),
        largest_field_width_m=(
            float(max(length for _, length in fields) * rmap.bin_size_m) if fields else float("nan")
        ),
        in_field_mean_rate=in_rate,
        out_field_mean_rate=out_rate,
        stability=stability,
        n_laps=n_laps,
        n_laps_active=laps_active(samples, slots),
        is_place_cell=passes_core and passes_stability,
        is_place_cell_no_stability=passes_core,
        rate_map=rmap,
    )
