"""Analysis-layer cleanup of the raw :class:`~hc11.io.Session`.

:func:`hc11.io.load_session` returns exactly what the ``.mat`` file contains. The
fixes for the anomalies found in step 2 live here instead, so every one of them is
an explicit, reversible call in analysis code rather than a hidden change at load
time (decision D2.3):

* :func:`truncate_to_session` — drop spikes and clip state intervals that fall
  after ``sess_duration`` (Gatsby_08022013; decision D3.1).
* :func:`clean_intervals` / :func:`state_intervals` — drop zero-length state
  intervals (Achilles_11012013; decision D3.2).
"""

from __future__ import annotations

import dataclasses
import logging
from dataclasses import dataclass

import numpy as np

from hc11.io import STATE_NAMES, Session

log = logging.getLogger(__name__)


@dataclass(frozen=True)
class TruncationReport:
    """What :func:`truncate_to_session` removed."""

    session: str
    sess_duration: float
    last_spike_time: float
    n_spikes_dropped: int
    n_spikes_total: int
    seconds_beyond_end: float
    clusters_affected: int
    state_intervals_dropped: dict[str, int]
    state_intervals_clipped: dict[str, int]
    state_seconds_removed: float

    @property
    def changed(self) -> bool:
        return bool(
            self.n_spikes_dropped
            or sum(self.state_intervals_dropped.values())
            or sum(self.state_intervals_clipped.values())
        )

    @property
    def fraction_spikes_dropped(self) -> float:
        return self.n_spikes_dropped / self.n_spikes_total if self.n_spikes_total else 0.0

    def __str__(self) -> str:
        if not self.changed:
            return f"{self.session}: nothing beyond sess_duration={self.sess_duration:.3f} s"
        dropped = {k: v for k, v in self.state_intervals_dropped.items() if v}
        clipped = {k: v for k, v in self.state_intervals_clipped.items() if v}
        return (
            f"{self.session}: dropped {self.n_spikes_dropped:,} spikes "
            f"({self.fraction_spikes_dropped:.2%} of {self.n_spikes_total:,}) from "
            f"{self.clusters_affected} clusters, spanning up to "
            f"{self.seconds_beyond_end:.1f} s past sess_duration="
            f"{self.sess_duration:.3f} s; state intervals dropped={dropped or '{}'}, "
            f"clipped={clipped or '{}'} ({self.state_seconds_removed:.1f} s of state time)"
        )


def truncate_to_session(
    session: Session, *, verbose: bool = True
) -> tuple[Session, TruncationReport]:
    """Restrict ``session`` to ``[0, sess_duration]``.

    Spikes after the end are dropped; state intervals are clipped to the end, and
    dropped entirely if they start after it. Position is untouched (it lies inside
    the MAZE epoch by construction).

    Only Gatsby_08022013 is affected in the released data: ``sess_duration``,
    ``POSTEpoch`` and ``Sessions_Recordings_Summary.pdf`` all agree that the session
    ends at 30 413.628 s, while spikes and some state intervals continue for another
    ~1 589 s. Returns the new session and a :class:`TruncationReport`; the report is
    returned even when nothing changed, so callers can log it unconditionally.
    """
    end = session.sess_duration
    keep = session.spike_times <= end
    n_dropped = int((~keep).sum())
    clusters = int(len(np.unique(session.spike_ids[~keep]))) if n_dropped else 0

    states, dropped, clipped = {}, {}, {}
    removed_s = 0.0
    for name in STATE_NAMES:
        iv = session.states[name]
        inside = iv[:, 0] < end
        dropped[name] = int((~inside).sum())
        removed_s += float(np.sum(iv[~inside, 1] - iv[~inside, 0]))
        kept = iv[inside].copy()
        over = kept[:, 1] > end
        clipped[name] = int(over.sum())
        removed_s += float(np.sum(kept[over, 1] - end))
        kept[over, 1] = end
        states[name] = kept

    report = TruncationReport(
        session=session.name,
        sess_duration=end,
        last_spike_time=(
            float(session.spike_times[-1]) if session.spike_times.size else float("nan")
        ),
        n_spikes_dropped=n_dropped,
        n_spikes_total=int(session.spike_times.size),
        seconds_beyond_end=float(session.spike_times.max() - end) if n_dropped else 0.0,
        clusters_affected=clusters,
        state_intervals_dropped=dropped,
        state_intervals_clipped=clipped,
        state_seconds_removed=removed_s,
    )
    if verbose and report.changed:
        log.warning("truncate_to_session: %s", report)

    if not report.changed:
        return session, report

    from types import MappingProxyType

    truncated = dataclasses.replace(
        session,
        spike_times=session.spike_times[keep],
        spike_ids=session.spike_ids[keep],
        states=MappingProxyType({k: _readonly(v) for k, v in states.items()}),
    )
    return truncated, report


def _readonly(a: np.ndarray) -> np.ndarray:
    a = np.ascontiguousarray(a)
    a.setflags(write=False)
    return a


def clean_intervals(
    intervals: np.ndarray, *, label: str = "intervals", verbose: bool = True
) -> np.ndarray:
    """Drop zero-length (and any negative-length) intervals.

    Achilles_11012013 has one zero-length REM interval, ``[13281, 13281]``. It
    contributes nothing to a duration sum, but code that assumes ``start < end`` --
    interval intersection, sample masking, midpoint interpolation -- can produce an
    empty or degenerate result from it, so it is removed up front (decision D3.2).
    """
    iv = np.asarray(intervals, dtype=np.float64)
    if iv.size == 0:
        return np.empty((0, 2), dtype=np.float64)
    good = iv[:, 1] > iv[:, 0]
    n_bad = int((~good).sum())
    if n_bad and verbose:
        log.warning(
            "clean_intervals: dropped %d zero- or negative-length %s: %s",
            n_bad, label, iv[~good].tolist(),
        )
    return iv[good]


def state_intervals(
    session: Session, state: str, *, restrict_to: np.ndarray | None = None, verbose: bool = True
) -> np.ndarray:
    """Cleaned ``(n, 2)`` intervals for one behavioural state.

    ``restrict_to`` optionally clips the result to an epoch such as
    ``session.epoch("MAZE")``; intervals outside it are dropped.
    """
    if state not in session.states:
        raise KeyError(f"Unknown state {state!r}; expected one of {STATE_NAMES}")
    iv = clean_intervals(session.states[state], label=f"{session.name} {state}", verbose=verbose)
    if restrict_to is None or iv.size == 0:
        return iv
    lo, hi = float(restrict_to[0]), float(restrict_to[1])
    iv = iv[(iv[:, 1] > lo) & (iv[:, 0] < hi)].copy()
    np.clip(iv, lo, hi, out=iv)
    return clean_intervals(iv, label=f"{session.name} {state} (restricted)", verbose=False)


def intervals_to_mask(t: np.ndarray, intervals: np.ndarray) -> np.ndarray:
    """Boolean mask of samples ``t`` falling inside any ``[start, end)`` interval."""
    mask = np.zeros(len(t), dtype=bool)
    for start, end in np.asarray(intervals).reshape(-1, 2):
        mask |= (t >= start) & (t < end)
    return mask


def mask_to_intervals(t: np.ndarray, mask: np.ndarray, *, dt: float | None = None) -> np.ndarray:
    """Contiguous ``True`` runs of ``mask`` as ``(n, 2)`` intervals on the ``t`` grid.

    Each run spans from its first sample to one sample interval past its last, so
    that interval durations sum to ``mask.sum() * dt`` rather than losing one sample
    per run.
    """
    t = np.asarray(t, dtype=np.float64)
    mask = np.asarray(mask, dtype=bool)
    if mask.size != t.size:
        raise ValueError(f"mask has {mask.size} samples but t has {t.size}")
    if not mask.any():
        return np.empty((0, 2), dtype=np.float64)
    step = float(np.median(np.diff(t))) if dt is None else float(dt)
    padded = np.concatenate([[False], mask, [False]])
    edges = np.flatnonzero(np.diff(padded.astype(np.int8)))
    starts, ends = edges[0::2], edges[1::2] - 1
    return np.column_stack([t[starts], t[ends] + step])
