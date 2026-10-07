"""Shared machinery for candidate-event detection.

Both detectors — population synchrony from spikes (:mod:`hc11.mua`) and
sharp-wave ripples from the local field potential (LFP, ``hc11.ripples``, step
6b) — reduce their signal to one trace per brain state and then hand it to the
same :func:`detect_two_threshold` here. That is deliberate: the two event lists
have to be compared directly in the 6c cross-check, and a difference in how
boundaries are defined would show up as a difference in the events themselves.

The detection rule, from Chen et al. (2016) and standard sharp-wave ripple
(SPW-R) practice:

1. a candidate must **peak** above a high threshold (default 3 SD above the
   state's own baseline);
2. its **boundaries** extend outwards to where the signal falls back below a low
   threshold (default 1 SD), so that duration reflects the whole event rather
   than the part that happened to clear the peak threshold;
3. candidates closer together than ``merge_gap`` are merged, then
4. the survivors are filtered on duration.

Thresholds are computed **within a brain state**, never session-wide: baseline
population rate and ripple power both differ between non-REM sleep and quiet
waking, so a single session-wide threshold would hand one epoch systematically
more events than another and make the PRE/MAZE/POST comparison meaningless.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from hc11.io import Session


@dataclass(frozen=True)
class EventList:
    """Candidate events from one detector, for one session and epoch.

    Both detectors return this shape, so the 6c cross-check is a comparison of
    two ``EventList`` objects and nothing else.
    """

    session: str
    epoch: str  # "PRE" / "MAZE" / "POST"
    detector: str  # "mua" / "lfp"
    start: np.ndarray  # (n,) s
    peak: np.ndarray  # (n,) s, time of the maximum
    end: np.ndarray  # (n,) s
    peak_z: np.ndarray  # (n,) peak height in baseline SD
    baseline_mean: float
    baseline_sd: float
    state_seconds: float  # time searched, i.e. the state's duration in this epoch
    params: dict = field(default_factory=dict)

    def __len__(self) -> int:
        return len(self.start)

    @property
    def duration(self) -> np.ndarray:
        return self.end - self.start

    @property
    def rate_per_min(self) -> float:
        return 60.0 * len(self) / self.state_seconds if self.state_seconds > 0 else float("nan")

    @property
    def intervals(self) -> np.ndarray:
        return np.column_stack([self.start, self.end]) if len(self) else np.empty((0, 2))

    def __str__(self) -> str:
        return (
            f"{self.session} {self.epoch} [{self.detector}]: {len(self)} events, "
            f"{self.rate_per_min:.1f}/min over {self.state_seconds / 60:.1f} min, "
            f"median duration {np.median(self.duration) * 1000:.0f} ms"
            if len(self) else
            f"{self.session} {self.epoch} [{self.detector}]: no events"
        )


def _runs(mask: np.ndarray) -> list[tuple[int, int]]:
    """Contiguous ``True`` runs of ``mask`` as ``(start, stop)`` index pairs."""
    if not mask.any():
        return []
    padded = np.concatenate([[False], mask, [False]])
    edges = np.flatnonzero(np.diff(padded.astype(np.int8)))
    return list(zip(edges[0::2], edges[1::2], strict=True))


def detect_two_threshold(
    signal: np.ndarray,
    t: np.ndarray,
    *,
    high: float,
    low: float,
    min_duration_s: float,
    max_duration_s: float,
    merge_gap_s: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Two-threshold event detection on one contiguous trace.

    Returns ``(start, peak, end)`` times. ``signal`` and ``t`` must be a single
    contiguous stretch: callers run this once per state interval so that an
    event can never span a gap in the searched time.
    """
    signal = np.asarray(signal, dtype=np.float64)
    t = np.asarray(t, dtype=np.float64)
    if signal.shape != t.shape:
        raise ValueError(f"signal has {signal.shape} samples but t has {t.shape}")
    if signal.size == 0:
        empty = np.empty(0)
        return empty, empty.copy(), empty.copy()
    dt = float(np.median(np.diff(t))) if t.size > 1 else 0.0

    candidates = [(s, e) for s, e in _runs(signal >= low) if signal[s:e].max() >= high]
    if not candidates:
        empty = np.empty(0)
        return empty, empty.copy(), empty.copy()

    merged = [candidates[0]]
    for start, stop in candidates[1:]:
        prev_start, prev_stop = merged[-1]
        if t[start] - t[prev_stop - 1] < merge_gap_s:
            merged[-1] = (prev_start, stop)
        else:
            merged.append((start, stop))

    starts, peaks, ends = [], [], []
    for start, stop in merged:
        begin, finish = t[start], t[stop - 1] + dt
        if not (min_duration_s <= finish - begin <= max_duration_s):
            continue
        starts.append(begin)
        peaks.append(t[start + int(np.argmax(signal[start:stop]))])
        ends.append(finish)
    return np.array(starts), np.array(peaks), np.array(ends)


# ---------------------------------------------------------------------------
# Per-event participation
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class Participation:
    """Which cells fired in each event, and how much."""

    n_active: np.ndarray  # (n_events,) distinct cells with >= 1 spike
    n_spikes: np.ndarray  # (n_events,) total spikes from the cell set
    n_cells: int  # size of the cell set the fractions are against

    @property
    def active_fraction(self) -> np.ndarray:
        return self.n_active / self.n_cells if self.n_cells else np.zeros_like(self.n_active)


def participation(
    session: Session, events: EventList, cluster_ids: np.ndarray
) -> Participation:
    """Count distinct active cells and total spikes in each event window."""
    cluster_ids = np.asarray(cluster_ids)
    wanted = np.isin(session.spike_ids, cluster_ids)
    times = session.spike_times[wanted]
    ids = session.spike_ids[wanted]
    n_active = np.zeros(len(events), dtype=np.int64)
    n_spikes = np.zeros(len(events), dtype=np.int64)
    for i, (start, end) in enumerate(zip(events.start, events.end, strict=True)):
        lo, hi = np.searchsorted(times, [start, end])
        n_spikes[i] = hi - lo
        n_active[i] = len(np.unique(ids[lo:hi])) if hi > lo else 0
    return Participation(n_active=n_active, n_spikes=n_spikes, n_cells=len(cluster_ids))


# ---------------------------------------------------------------------------
# Comparing two event lists (the 6c cross-check)
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class EventOverlap:
    """How two event lists correspond."""

    n_a: int
    n_b: int
    n_a_matched: int
    n_b_matched: int
    offsets_s: np.ndarray  # peak(b) - peak(a) for matched pairs
    match_window_s: float

    @property
    def fraction_a_matched(self) -> float:
        return self.n_a_matched / self.n_a if self.n_a else float("nan")

    @property
    def fraction_b_matched(self) -> float:
        return self.n_b_matched / self.n_b if self.n_b else float("nan")


def compare_events(a: EventList, b: EventList, *, match_window_s: float = 0.05) -> EventOverlap:
    """Match events between two detectors by overlap, then by peak proximity.

    An event of ``a`` counts as matched if any event of ``b`` overlaps it in
    time, or if their peaks fall within ``match_window_s``. Offsets are reported
    for the nearest-peak partner of each matched ``a`` event.
    """
    if len(a) == 0 or len(b) == 0:
        return EventOverlap(len(a), len(b), 0, 0, np.empty(0), match_window_s)

    a_matched = np.zeros(len(a), dtype=bool)
    b_matched = np.zeros(len(b), dtype=bool)
    offsets = []
    for i, (start, end, peak) in enumerate(zip(a.start, a.end, a.peak, strict=True)):
        overlaps = (b.start < end) & (b.end > start)
        near = np.abs(b.peak - peak) <= match_window_s
        hits = overlaps | near
        if hits.any():
            a_matched[i] = True
            b_matched |= hits
            offsets.append(b.peak[hits][np.argmin(np.abs(b.peak[hits] - peak))] - peak)
    return EventOverlap(
        n_a=len(a), n_b=len(b),
        n_a_matched=int(a_matched.sum()), n_b_matched=int(b_matched.sum()),
        offsets_s=np.array(offsets), match_window_s=match_window_s,
    )
