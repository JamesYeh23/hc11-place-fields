"""Population-synchrony (multi-unit activity, MUA) event detection from spikes.

This detector uses spike times only — no local field potential (LFP). It exists
as an independent check on ripple detection: the ripple band alone also picks up
electrical artifacts, which have no reason to coincide with a burst of
population firing. Where the two agree, the event is real.

Method: bin every pyramidal spike at 1 ms, smooth with a Gaussian, and threshold
the resulting population rate with the shared two-threshold rule in
:mod:`hc11.events`.

Two things about how the baseline is computed matter more than the thresholds:

* **Per state, never session-wide.** Baseline population rate in non-REM sleep
  differs from quiet waking on the maze, so one session-wide threshold would
  hand one epoch systematically more events than another — and the whole point
  of Fig. 1C is to compare event content across PRE, MAZE and POST.
* **Per state interval, never across gaps.** The rate trace is built and
  searched inside each state interval separately, so no event can span a gap in
  the searched time, and smoothing never bleeds activity across a boundary.
"""

from __future__ import annotations

import numpy as np
from scipy.ndimage import gaussian_filter1d

from hc11.behavior import bridge_short_gaps, speed_cm_s
from hc11.events import EventList, detect_two_threshold
from hc11.io import Session
from hc11.preprocess import mask_to_intervals, state_intervals


def state_intervals_for_epoch(
    session: Session, epoch: str, params: dict, *, verbose: bool = False
) -> np.ndarray:
    """The periods of ``epoch`` in which candidate events may be detected.

    PRE and POST use scored non-REM sleep. MAZE uses immobility — speed below
    ``maze_immobility_speed_cm_s`` — which is Chen et al.'s quiet-wakefulness
    criterion, and is where sharp-wave ripples occur on the track. Samples with
    no speed estimate (tracking lost) are excluded rather than assumed still.
    """
    cfg = params["events"]["states"]
    if epoch.upper() not in ("PRE", "MAZE", "POST"):
        raise ValueError(f"Unknown epoch {epoch!r}; expected PRE, MAZE or POST")
    window = session.epoch(epoch)
    if epoch.upper() in ("PRE", "POST"):
        return state_intervals(session, cfg["sleep_state"], restrict_to=window, verbose=verbose)

    speed = speed_cm_s(session, cfg["speed_smoothing_sigma_s"])
    measured = np.isfinite(speed)
    still = measured & (speed < cfg["maze_immobility_speed_cm_s"])
    # Position differentiated from 39 Hz LED tracking is noisy enough that a
    # stationary animal's apparent speed crosses 2 cm/s constantly, which chops
    # one quiet period into dozens of sub-threshold fragments. Bridge brief
    # excursions before applying the minimum duration -- the same treatment
    # running periods get in hc11.behavior. Samples with no speed estimate
    # (tracking lost) are never bridged into immobility: we cannot tell whether
    # the animal was still.
    still, _ = bridge_short_gaps(
        still, allowed=measured,
        max_gap=int(round(cfg["max_immobility_gap_s"] / session.position_dt)),
    )
    intervals = mask_to_intervals(session.position_t, still, dt=session.position_dt)
    if intervals.size == 0:
        return intervals
    intervals = intervals[
        (intervals[:, 1] > window[0]) & (intervals[:, 0] < window[1])
    ].copy()
    np.clip(intervals, window[0], window[1], out=intervals)
    keep = intervals[:, 1] - intervals[:, 0] >= cfg["min_immobility_s"]
    return intervals[keep]


def population_rate(
    session: Session,
    interval: np.ndarray,
    cluster_ids: np.ndarray,
    *,
    bin_s: float,
    sigma_s: float,
) -> tuple[np.ndarray, np.ndarray]:
    """Smoothed population firing rate (spikes/s, summed over cells) on one interval.

    Returns ``(t, rate)`` on a regular ``bin_s`` grid covering ``interval``.
    """
    start, end = float(interval[0]), float(interval[1])
    n_bins = max(int(round((end - start) / bin_s)), 1)
    edges = start + np.arange(n_bins + 1) * bin_s
    times, ids = session.spikes_in((start, end))
    if times.size:
        times = times[np.isin(ids, cluster_ids)]
    counts = np.histogram(times, bins=edges)[0].astype(np.float64) if times.size else (
        np.zeros(n_bins)
    )
    rate = gaussian_filter1d(counts, sigma_s / bin_s, mode="nearest") / bin_s
    return edges[:-1] + bin_s / 2, rate


def detect_mua_events(
    session: Session,
    epoch: str,
    params: dict,
    *,
    cluster_ids: np.ndarray | None = None,
) -> tuple[EventList, np.ndarray]:
    """Detect population-synchrony events in one epoch.

    Returns the events and the state intervals they were searched in. The
    baseline mean and SD are pooled across all of the epoch's state intervals,
    in one pass, before any detection happens.
    """
    cfg = params["events"]["mua"]
    ids = session.units(cfg["unit_kind"]) if cluster_ids is None else np.asarray(cluster_ids)
    intervals = state_intervals_for_epoch(session, epoch, params)
    bin_s = cfg["bin_ms"] / 1000.0
    sigma_s = cfg["smoothing_sigma_ms"] / 1000.0

    # Pass 1: pooled baseline across the epoch's state intervals.
    total = 0.0
    total_sq = 0.0
    n = 0
    traces = []
    for interval in intervals:
        t, rate = population_rate(session, interval, ids, bin_s=bin_s, sigma_s=sigma_s)
        traces.append((t, rate))
        total += rate.sum()
        total_sq += np.square(rate).sum()
        n += rate.size
    state_seconds = float(np.sum(intervals[:, 1] - intervals[:, 0])) if len(intervals) else 0.0
    if n == 0:
        empty = np.empty(0)
        return EventList(
            session=session.name, epoch=epoch, detector="mua", start=empty,
            peak=empty.copy(), end=empty.copy(), peak_z=empty.copy(),
            baseline_mean=float("nan"), baseline_sd=float("nan"),
            state_seconds=state_seconds, params=dict(cfg),
        ), intervals

    mean = total / n
    sd = float(np.sqrt(max(total_sq / n - mean**2, 0.0)))
    high = mean + cfg["threshold_high_sd"] * sd
    low = mean + cfg["threshold_low_sd"] * sd

    # Pass 2: detect inside each interval separately.
    starts, peaks, ends, peak_z = [], [], [], []
    for t, rate in traces:
        s, p, e = detect_two_threshold(
            rate, t, high=high, low=low,
            min_duration_s=cfg["min_duration_ms"] / 1000.0,
            max_duration_s=cfg["max_duration_ms"] / 1000.0,
            merge_gap_s=cfg["merge_gap_ms"] / 1000.0,
        )
        if len(s) == 0:
            continue
        starts.append(s)
        peaks.append(p)
        ends.append(e)
        at_peak = np.interp(p, t, rate)
        peak_z.append((at_peak - mean) / sd if sd > 0 else np.full(len(p), np.nan))

    def _cat(chunks):
        return np.concatenate(chunks) if chunks else np.empty(0)

    events = EventList(
        session=session.name, epoch=epoch, detector="mua",
        start=_cat(starts), peak=_cat(peaks), end=_cat(ends), peak_z=_cat(peak_z),
        baseline_mean=mean, baseline_sd=sd, state_seconds=state_seconds,
        params={**dict(cfg), "n_cells": len(ids), "threshold_high": high, "threshold_low": low},
    )
    return events, intervals
