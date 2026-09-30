"""Speed and running periods during the MAZE epoch.

The running mask is built in two stages (decision D3.4):

1. **The authors' mask is primary.** ``OneDLocation`` is NaN outside track
   running, so a defined linearised position *is* the released definition of "on
   the track, running". Nothing else can be used for a 1-D rate map anyway, since
   a sample without a linearised position has no bin to go in.
2. **Our speed threshold is applied on top**, as a check rather than a filter. If
   the authors' mask already encodes a speed criterion, thresholding at
   15 cm/s should remove very little; a large removal would mean we have
   misunderstood what the mask is. :func:`running_mask` reports exactly how much
   each stage removes.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from scipy.ndimage import gaussian_filter1d

from hc11.io import Session
from hc11.preprocess import mask_to_intervals


def _nan_aware_gaussian(x: np.ndarray, sigma_samples: float) -> np.ndarray:
    """Gaussian smoothing that ignores NaNs instead of spreading them.

    Normalised convolution: smooth the zero-filled signal and the validity mask
    with the same kernel, then divide. Samples that are NaN in the input stay NaN.
    """
    valid = np.isfinite(x)
    if sigma_samples <= 0 or not valid.any():
        return x.copy()
    filled = np.where(valid, x, 0.0)
    num = gaussian_filter1d(filled, sigma_samples, mode="nearest")
    den = gaussian_filter1d(valid.astype(np.float64), sigma_samples, mode="nearest")
    with np.errstate(invalid="ignore", divide="ignore"):
        out = np.where(den > 0, num / den, np.nan)
    out[~valid] = np.nan
    return out


def speed_cm_s(session: Session, smooth_sigma_s: float = 0.0) -> np.ndarray:
    """Instantaneous running speed from 2-D tracking, in cm/s.

    Speed is differentiated from the raw 2-D position and then smoothed, rather
    than the other way round, so that a tracking gap never produces a spurious
    jump: differences spanning a NaN are NaN, and the smoothing step ignores them.
    ``NaN`` is returned wherever the speed is undefined.
    """
    t = session.position_t
    xy = session.position_xy_cm
    with np.errstate(invalid="ignore"):
        dxy = np.gradient(xy, t, axis=0)
    v = np.hypot(dxy[:, 0], dxy[:, 1])
    if smooth_sigma_s > 0:
        v = _nan_aware_gaussian(v, smooth_sigma_s / session.position_dt)
    return v


@dataclass(frozen=True)
class RunningMask:
    """A per-sample running mask over the MAZE epoch, with provenance."""

    session: str
    mask: np.ndarray  # (M,) bool, samples used for place fields
    speed: np.ndarray  # (M,) cm/s, NaN where undefined
    authors_mask: np.ndarray  # (M,) bool, OneDLocation is defined
    intervals: np.ndarray  # (n, 2) contiguous running periods, s
    dt: float
    speed_threshold_cm_s: float
    n_removed_by_speed: int
    n_restored_by_bridging: int
    n_removed_short_epochs: int
    params: dict = field(default_factory=dict)

    @property
    def n_authors(self) -> int:
        return int(self.authors_mask.sum())

    @property
    def n_final(self) -> int:
        return int(self.mask.sum())

    @property
    def seconds(self) -> float:
        return self.n_final * self.dt

    @property
    def authors_seconds(self) -> float:
        return self.n_authors * self.dt

    @property
    def fraction_removed_by_speed(self) -> float:
        """Share of the authors' mask that our speed threshold rejects."""
        return self.n_removed_by_speed / self.n_authors if self.n_authors else float("nan")

    @property
    def fraction_kept(self) -> float:
        return self.n_final / self.n_authors if self.n_authors else float("nan")

    def median_speed(self, which: str = "final") -> float:
        m = self.mask if which == "final" else self.authors_mask
        return float(np.nanmedian(self.speed[m])) if m.any() else float("nan")

    def __str__(self) -> str:
        return (
            f"{self.session}: authors' mask {self.authors_seconds / 60:.1f} min "
            f"({self.n_authors:,} samples); speed > {self.speed_threshold_cm_s:g} cm/s removes "
            f"{self.n_removed_by_speed:,} ({self.fraction_removed_by_speed:.1%}); "
            f"bridging restores {self.n_restored_by_bridging:,}; short epochs remove "
            f"{self.n_removed_short_epochs:,}; final {self.seconds / 60:.1f} min "
            f"({self.fraction_kept:.1%} kept) in {len(self.intervals)} periods; "
            f"median speed {self.median_speed():.1f} cm/s"
        )


def _bridge_short_gaps(
    mask: np.ndarray, allowed: np.ndarray, max_gap: int
) -> tuple[np.ndarray, int]:
    """Fill ``False`` runs of at most ``max_gap`` samples that sit between ``True`` runs.

    Only samples where ``allowed`` holds may be filled, so bridging never invents a
    running sample at a time with no linearised position.
    """
    if max_gap <= 0 or not mask.any():
        return mask.copy(), 0
    out = mask.copy()
    padded = np.concatenate([[True], mask, [True]])
    edges = np.flatnonzero(np.diff(padded.astype(np.int8)))
    restored = 0
    for start, stop in zip(edges[0::2], edges[1::2], strict=True):  # False runs
        if start == 0 or stop == len(mask):  # not enclosed by running on both sides
            continue
        if stop - start <= max_gap and allowed[start:stop].all():
            out[start:stop] = True
            restored += stop - start
    return out, restored


def _drop_short_runs(mask: np.ndarray, min_len: int) -> tuple[np.ndarray, int]:
    """Remove ``True`` runs shorter than ``min_len`` samples."""
    if min_len <= 1 or not mask.any():
        return mask.copy(), 0
    out = mask.copy()
    padded = np.concatenate([[False], mask, [False]])
    edges = np.flatnonzero(np.diff(padded.astype(np.int8)))
    removed = 0
    for start, stop in zip(edges[0::2], edges[1::2], strict=True):  # True runs
        if stop - start < min_len:
            out[start:stop] = False
            removed += stop - start
    return out, removed


def running_mask(
    session: Session,
    *,
    speed_threshold_cm_s: float,
    speed_smoothing_sigma_s: float,
    min_run_epoch_s: float,
    max_run_gap_s: float,
    require_linearized: bool = True,
) -> RunningMask:
    """Samples to use for 1-D place fields: on the track, running.

    With ``require_linearized=True`` (the default, and the only option that can
    feed a 1-D rate map) the mask starts from the samples where ``OneDLocation``
    is defined. Set it to ``False`` only for the sensitivity analysis that asks
    what a speed criterion alone would select.
    """
    dt = session.position_dt
    speed = speed_cm_s(session, speed_smoothing_sigma_s)
    authors = ~np.isnan(session.position_1d)
    base = authors if require_linearized else np.ones_like(authors)

    with np.errstate(invalid="ignore"):
        fast = np.nan_to_num(speed, nan=-1.0) >= speed_threshold_cm_s
    after_speed = base & fast
    n_removed_by_speed = int(base.sum() - after_speed.sum())

    bridged, n_restored = _bridge_short_gaps(
        after_speed, allowed=base, max_gap=int(round(max_run_gap_s / dt))
    )
    final, n_short = _drop_short_runs(bridged, min_len=int(round(min_run_epoch_s / dt)))

    return RunningMask(
        session=session.name,
        mask=final,
        speed=speed,
        authors_mask=authors,
        intervals=mask_to_intervals(session.position_t, final, dt=dt),
        dt=dt,
        speed_threshold_cm_s=speed_threshold_cm_s,
        n_removed_by_speed=n_removed_by_speed,
        n_restored_by_bridging=n_restored,
        n_removed_short_epochs=n_short,
        params={
            "speed_threshold_cm_s": speed_threshold_cm_s,
            "speed_smoothing_sigma_s": speed_smoothing_sigma_s,
            "min_run_epoch_s": min_run_epoch_s,
            "max_run_gap_s": max_run_gap_s,
            "require_linearized": require_linearized,
        },
    )


def running_mask_from_config(session: Session, params: dict, **overrides) -> RunningMask:
    """:func:`running_mask` with arguments taken from ``config/params.yaml``."""
    behavior = params["behavior"]
    kwargs = {
        "speed_threshold_cm_s": behavior["speed_threshold_cm_s"],
        "speed_smoothing_sigma_s": behavior["speed_smoothing_sigma_s"],
        "min_run_epoch_s": behavior["min_run_epoch_s"],
        "max_run_gap_s": behavior["max_run_gap_s"],
        "require_linearized": behavior.get("require_linearized", True),
    }
    kwargs.update(overrides)
    return running_mask(session, **kwargs)
