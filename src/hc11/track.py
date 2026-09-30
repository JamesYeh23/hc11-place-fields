"""Track geometry: spatial binning, and the circular wrap made explicit.

A :class:`Track` is built from the linearised positions actually observed in a
session, **not** from an assumed track length (decision D3.3). On the circular
sessions the observed 1-D range (2.84-2.90 m) is shorter than the 3.14 m
circumference implied by the documented 1 m diameter, for reasons we have not
established, so an assumed circumference would put every bin edge in the wrong
place and would make the wrap discontinuity land at the wrong position.

Every position comparison goes through :meth:`Track.difference` or
:meth:`Track.distance`, which take the short way around on a circular track. The
reward site is at position 0 on the circular maze, so 0 and ``period`` are the
same place: a plain ``a - b`` would report the two sides of the reward site as
being a full lap apart.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Literal

import numpy as np

from hc11.io import Session

TrackKind = Literal["linear", "circular"]


@dataclass(frozen=True)
class Track:
    """Spatial extent and binning of one session's linearised track.

    ``lo``/``hi`` are the observed range of ``OneDLocation`` in metres. For a
    circular track the period is ``hi - lo``: positions are dense, so this
    underestimates the true circumference by at most one inter-sample step
    (< 1.5 cm at the observed running speeds), which is well inside a 10 cm bin.
    """

    kind: TrackKind
    lo: float
    hi: float
    bin_size: float
    nominal_length: float | None = None  # from the maze name; linear only
    session: str | None = None

    # --- construction -------------------------------------------------------

    @classmethod
    def from_session(cls, session: Session, bin_size_cm: float) -> Track:
        x = session.position_1d
        if np.all(np.isnan(x)):
            raise ValueError(f"{session.name}: OneDLocation is entirely NaN")
        return cls(
            kind=session.maze_kind,
            lo=float(np.nanmin(x)),
            hi=float(np.nanmax(x)),
            bin_size=bin_size_cm / 100.0,
            nominal_length=session.track_length_m,
            session=session.name,
        )

    # --- extent -------------------------------------------------------------

    @property
    def extent(self) -> float:
        """Observed span of the track, in metres."""
        return self.hi - self.lo

    @property
    def period(self) -> float | None:
        """Wrap period in metres for a circular track; ``None`` if linear."""
        return self.extent if self.kind == "circular" else None

    @property
    def is_circular(self) -> bool:
        return self.kind == "circular"

    # --- binning ------------------------------------------------------------

    @property
    def n_bins(self) -> int:
        return max(1, int(round(self.extent / self.bin_size)))

    @property
    def bin_edges(self) -> np.ndarray:
        """``n_bins + 1`` edges tiling ``[lo, hi]`` exactly."""
        return np.linspace(self.lo, self.hi, self.n_bins + 1)

    @property
    def bin_centers(self) -> np.ndarray:
        edges = self.bin_edges
        return 0.5 * (edges[:-1] + edges[1:])

    @property
    def actual_bin_size(self) -> float:
        """Bin width actually used (``extent / n_bins``), in metres."""
        return self.extent / self.n_bins

    def digitize(self, x: np.ndarray) -> np.ndarray:
        """Bin index for each position; ``-1`` where the position is NaN.

        Positions exactly at ``hi`` fall in the last bin (linear) or in bin 0
        (circular, where ``hi`` is the same place as ``lo``).
        """
        x = np.asarray(x, dtype=np.float64)
        out = np.full(x.shape, -1, dtype=np.int64)
        ok = ~np.isnan(x)
        if not ok.any():
            return out
        idx = np.floor((x[ok] - self.lo) / self.actual_bin_size).astype(np.int64)
        if self.is_circular:
            idx %= self.n_bins
        else:
            idx = np.clip(idx, 0, self.n_bins - 1)
        out[ok] = idx
        return out

    # --- wrap-aware arithmetic ---------------------------------------------

    def difference(self, a: np.ndarray | float, b: np.ndarray | float) -> np.ndarray | float:
        """Signed displacement ``a - b``, taking the short way around if circular.

        The result lies in ``[-period/2, period/2)`` for a circular track, so the
        sign still means "forward" or "backward" for steps smaller than half a lap.
        """
        d = np.subtract(a, b, dtype=np.float64)
        if not self.is_circular:
            return d
        period = self.extent
        return (d + period / 2.0) % period - period / 2.0

    def distance(self, a: np.ndarray | float, b: np.ndarray | float) -> np.ndarray | float:
        """Absolute distance, taking the short way around if circular."""
        return np.abs(self.difference(a, b))

    def wrap(self, x: np.ndarray | float) -> np.ndarray | float:
        """Fold positions into ``[lo, hi)`` (circular) or leave them alone (linear)."""
        if not self.is_circular:
            return np.asarray(x, dtype=np.float64)
        return self.lo + np.mod(np.subtract(x, self.lo, dtype=np.float64), self.extent)

    def unwrap(self, x: np.ndarray) -> np.ndarray:
        """Undo wrapping along a trajectory, giving a continuous cumulative position.

        Needed for direction and lap detection on a circular track, where a lap
        crossing the reward site would otherwise look like a huge jump backwards.
        NaNs are preserved; each run of valid samples is unwrapped on its own, so
        no jump is inferred across a tracking gap.
        """
        x = np.asarray(x, dtype=np.float64)
        if not self.is_circular:
            return x.copy()
        out = np.full(x.shape, np.nan)
        ok = ~np.isnan(x)
        if not ok.any():
            return out
        # Unwrap within each contiguous run of valid samples.
        edges = np.flatnonzero(np.diff(np.concatenate([[False], ok, [False]]).astype(np.int8)))
        for start, stop in zip(edges[0::2], edges[1::2], strict=True):
            seg = x[start:stop]
            steps = self.difference(seg[1:], seg[:-1])
            out[start:stop] = seg[0] + np.concatenate([[0.0], np.cumsum(steps)])
        return out

    def __str__(self) -> str:
        nominal = f", nominal {self.nominal_length} m" if self.nominal_length else ""
        return (
            f"Track({self.session}, {self.kind}, observed [{self.lo:.3f}, {self.hi:.3f}] m"
            f"{nominal}, {self.n_bins} bins of {self.actual_bin_size * 100:.2f} cm)"
        )
