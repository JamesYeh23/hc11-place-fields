"""Segment MAZE running into directional traversals ("laps").

Segmentation runs on the **authors' mask blocks** -- the contiguous runs of samples
where ``OneDLocation`` is defined -- and *not* on the speed-filtered mask
(decision D4.1). A traversal that dips below the speed threshold for a moment
would otherwise be cut into several short laps, none of which covers enough of
the track to qualify. The speed filter is applied afterwards, when occupancy and
spikes are accumulated (see :mod:`hc11.fields`).

Direction comes from the sign of a smoothed derivative of the linearised position
(unwrapped first on a circular track), with hysteresis: the direction only flips
when the smoothed speed exceeds a threshold in the opposite sense, so a stationary
wobble does not produce a direction change. No alternation is assumed -- two
outbound runs in a row are perfectly possible and are segmented as such.

On a circular track a run is split again at every crossing of the reward site
(position 0, the wrap point), so that one lap is one circuit. Without this an
animal running continuously round the ring would yield a single "lap" spanning
many circuits.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.ndimage import gaussian_filter1d

from hc11.io import Session
from hc11.track import Track

#: Direction labels. "pos" = increasing linearised position, "neg" = decreasing.
DIRECTIONS: tuple[str, str] = ("pos", "neg")


@dataclass(frozen=True)
class Lap:
    """One directional traversal of the track."""

    index: int
    direction: str  # "pos" or "neg"
    start_sample: int  # inclusive index into session.position_t
    stop_sample: int  # exclusive
    t_start: float
    t_end: float
    x_start: float
    x_end: float
    coverage: float  # |net displacement| / track extent
    duration: float

    @property
    def n_samples(self) -> int:
        return self.stop_sample - self.start_sample

    @property
    def interval(self) -> np.ndarray:
        return np.array([self.t_start, self.t_end])


def _mask_blocks(mask: np.ndarray) -> list[tuple[int, int]]:
    """Contiguous ``True`` runs of ``mask`` as ``(start, stop)`` index pairs."""
    padded = np.concatenate([[False], mask, [False]])
    edges = np.flatnonzero(np.diff(padded.astype(np.int8)))
    return list(zip(edges[0::2], edges[1::2], strict=True))


def _direction_with_hysteresis(velocity: np.ndarray, threshold: float) -> np.ndarray:
    """Direction (+1 / -1 / 0) from a velocity trace, holding through slow patches.

    A sample is ``+1`` once velocity rises above ``threshold`` and stays ``+1``
    until velocity falls below ``-threshold``. Leading samples before the first
    clear movement are ``0``.
    """
    out = np.zeros(len(velocity), dtype=np.int8)
    state = 0
    for i, v in enumerate(velocity):
        if v > threshold:
            state = 1
        elif v < -threshold:
            state = -1
        out[i] = state
    return out


def segment_laps(
    session: Session,
    track: Track,
    *,
    position_smoothing_sigma_s: float,
    direction_threshold_cm_s: float,
    min_coverage_frac: float,
    min_lap_duration_s: float,
    max_lap_duration_s: float,
    base_mask: np.ndarray | None = None,
) -> list[Lap]:
    """Split the MAZE epoch into directional traversals.

    ``base_mask`` defaults to "``OneDLocation`` is defined". Pass a different mask
    only to test alternatives; passing the speed-filtered mask will fragment laps.
    """
    x_raw = session.position_1d
    mask = ~np.isnan(x_raw) if base_mask is None else np.asarray(base_mask, dtype=bool)
    t = session.position_t
    dt = session.position_dt
    sigma = max(position_smoothing_sigma_s / dt, 1e-6)
    x_unwrapped = track.unwrap(x_raw)

    laps: list[Lap] = []
    for start, stop in _mask_blocks(mask):
        if stop - start < 3:
            continue
        x = x_unwrapped[start:stop]
        if np.any(np.isnan(x)):
            continue
        smoothed = gaussian_filter1d(x, sigma, mode="nearest")
        velocity = np.gradient(smoothed, dt) * 100.0  # cm/s, signed along the track
        direction = _direction_with_hysteresis(velocity, direction_threshold_cm_s)

        # Split the block wherever the direction changes, and -- on a circular
        # track -- additionally at every crossing of the reward site at `lo`,
        # so that one lap is one circuit rather than one uninterrupted run.
        cuts = list(np.flatnonzero(np.diff(direction)) + 1)
        if track.is_circular:
            circuit = np.floor((x - track.lo) / track.extent)
            cuts.extend((np.flatnonzero(np.diff(circuit)) + 1).tolist())
        cuts = np.unique(np.asarray(cuts, dtype=np.int64)) if cuts else np.empty(0, dtype=np.int64)
        for seg_start, seg_stop in zip(
            np.concatenate([[0], cuts]), np.concatenate([cuts, [len(direction)]]), strict=True
        ):
            sign = direction[seg_start]
            if sign == 0:
                continue
            i0, i1 = start + seg_start, start + seg_stop
            duration = (i1 - i0) * dt
            if not (min_lap_duration_s <= duration <= max_lap_duration_s):
                continue
            net = abs(x_unwrapped[i1 - 1] - x_unwrapped[i0])
            coverage = net / track.extent
            if coverage < min_coverage_frac:
                continue
            laps.append(
                Lap(
                    index=len(laps),
                    direction="pos" if sign > 0 else "neg",
                    start_sample=int(i0),
                    stop_sample=int(i1),
                    t_start=float(t[i0]),
                    t_end=float(t[i1 - 1] + dt),
                    x_start=float(x_raw[i0]),
                    x_end=float(x_raw[i1 - 1]),
                    coverage=float(coverage),
                    duration=float(duration),
                )
            )
    return laps


def segment_laps_from_config(
    session: Session, track: Track, params: dict, **overrides
) -> list[Lap]:
    """:func:`segment_laps` with arguments from ``config/params.yaml``."""
    cfg = params["laps"]
    kwargs = {
        "position_smoothing_sigma_s": cfg["position_smoothing_sigma_s"],
        "direction_threshold_cm_s": cfg["direction_threshold_cm_s"],
        "min_coverage_frac": cfg["min_coverage_frac"],
        "min_lap_duration_s": cfg["min_lap_duration_s"],
        "max_lap_duration_s": cfg["max_lap_duration_s"],
    }
    kwargs.update(overrides)
    return segment_laps(session, track, **kwargs)


def lap_mask(session: Session, laps: list[Lap], direction: str | None = None) -> np.ndarray:
    """Boolean sample mask covering the given laps (optionally one direction)."""
    mask = np.zeros(len(session.position_t), dtype=bool)
    for lap in laps:
        if direction is None or lap.direction == direction:
            mask[lap.start_sample : lap.stop_sample] = True
    return mask


def lap_summary(laps: list[Lap]) -> dict:
    """Counts and durations per direction, for reporting."""
    out: dict = {}
    for d in DIRECTIONS:
        sel = [lap for lap in laps if lap.direction == d]
        durations = np.array([lap.duration for lap in sel]) if sel else np.array([])
        out[d] = {
            "n_laps": len(sel),
            "total_s": float(durations.sum()) if sel else 0.0,
            "median_duration_s": float(np.median(durations)) if sel else float("nan"),
            "median_coverage": (
                float(np.median([lap.coverage for lap in sel])) if sel else float("nan")
            ),
        }
    out["n_laps_total"] = len(laps)
    return out
