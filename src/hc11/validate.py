"""Consistency checks on a loaded :class:`~hc11.io.Session`.

:func:`validate` never raises on bad data: it returns a list of :class:`Issue`
records, each with a severity:

``error``
    The session violates an invariant that downstream code relies on (e.g.
    unsorted spike times, overlapping PRE/MAZE epochs, a unit labelled both
    pyramidal and interneuron). Analyses should not proceed.
``warning``
    A real anomaly in the data that analyses must account for, but which does
    not make the session unusable (e.g. spikes recorded after the documented
    end of the session, zero-length state intervals).
``info``
    A descriptive measurement reported for the record (sampling interval,
    fraction of missing position, linearised-position range).
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Literal

import numpy as np

from hc11.io import EPOCH_NAMES, STATE_NAMES, Session

Severity = Literal["error", "warning", "info"]

#: Expected position sampling interval (s). The description quotes ~39 Hz.
NOMINAL_POSITION_DT_S = 0.0256
#: Relative tolerance on the median sampling interval before warning.
POSITION_DT_RTOL = 0.02
#: Floating-point slack for boundary comparisons (s). Epoch boundaries are
#: stored to 0.1 ms or coarser; spike times to 50 us (20 kHz).
TIME_EPS_S = 1e-6
#: Allowed overshoot of linearised position beyond the nominal track length.
TRACK_LENGTH_TOL_M = 0.01


@dataclass(frozen=True)
class Issue:
    severity: Severity
    code: str
    message: str

    def __str__(self) -> str:
        return f"[{self.severity.upper():7s}] {self.code}: {self.message}"


def _interval_overlap_s(a: np.ndarray, b: np.ndarray) -> tuple[int, float]:
    """Number of overlapping (row of a, row of b) pairs and total overlap (s)."""
    if len(a) == 0 or len(b) == 0:
        return 0, 0.0
    lo = np.maximum(a[:, None, 0], b[None, :, 0])
    hi = np.minimum(a[:, None, 1], b[None, :, 1])
    ov = np.clip(hi - lo, 0.0, None)
    return int(np.count_nonzero(ov > 0)), float(ov.sum())


def _check_spikes(s: Session) -> list[Issue]:
    out: list[Issue] = []
    t = s.spike_times
    if len(t) != len(s.spike_ids):
        out.append(Issue("error", "spike_length_mismatch",
                         f"{len(t)} spike times vs {len(s.spike_ids)} spike IDs"))
    if not np.all(np.isfinite(t)):
        out.append(Issue("error", "spike_times_nonfinite",
                         f"{np.count_nonzero(~np.isfinite(t))} non-finite spike times"))
    if np.any(np.diff(t) < 0):
        out.append(Issue("error", "spike_times_unsorted",
                         f"{np.count_nonzero(np.diff(t) < 0)} decreasing steps in spike_times"))
    before = t < -TIME_EPS_S
    after = t > s.sess_duration + TIME_EPS_S
    if before.any():
        out.append(Issue("warning", "spikes_before_session_start",
                         f"{before.sum():,} spikes before t=0 (earliest {t.min():.4f} s)"))
    if after.any():
        out.append(Issue(
            "warning", "spikes_after_session_end",
            f"{after.sum():,} spikes ({after.mean():.2%}) after sess_duration="
            f"{s.sess_duration:.3f} s, spanning {t[after].min():.3f}-{t.max():.3f} s "
            f"({t.max() - s.sess_duration:.0f} s beyond the end); "
            f"{len(np.unique(s.spike_ids[after]))} clusters affected",
        ))
    return out


def _check_units(s: Session) -> list[Issue]:
    out: list[Issue] = []
    present = s.cluster_ids
    both = np.intersect1d(s.pyr_ids, s.int_ids)
    if both.size:
        out.append(Issue("error", "pyr_int_overlap",
                         f"IDs labelled both pyramidal and interneuron: {both.tolist()}"))
    for label, ids in (("pyr_ids", s.pyr_ids), ("int_ids", s.int_ids)):
        if len(np.unique(ids)) != len(ids):
            out.append(Issue("error", f"{label}_duplicates", f"{label} contains duplicate IDs"))
        silent = np.setdiff1d(ids, present)
        if silent.size:
            out.append(Issue("warning", f"{label}_without_spikes",
                             f"{label} not present in spike_ids: {silent.tolist()}"))
    unclassified = np.setdiff1d(present, np.union1d(s.pyr_ids, s.int_ids))
    if unclassified.size:
        n_spk = np.isin(s.spike_ids, unclassified).sum()
        out.append(Issue("warning", "unclassified_clusters",
                         f"{unclassified.size} clusters are neither pyr nor int "
                         f"({n_spk:,} spikes): {unclassified.tolist()}"))
    within = present % 100
    if np.any(within < 2):
        out.append(Issue("warning", "noise_or_mua_cluster_present",
                         "cluster numbers 0 (noise) or 1 (multi-unit) found: "
                         f"{present[within < 2].tolist()}"))
    return out


def _check_epochs(s: Session) -> list[Issue]:
    out: list[Issue] = []
    ep = {k: s.epoch(k) for k in EPOCH_NAMES}
    for k, (a, b) in ep.items():
        if not a < b:
            out.append(Issue("error", "epoch_nonpositive", f"{k} epoch has start >= end: {a}, {b}"))
        if a < -TIME_EPS_S or b > s.sess_duration + TIME_EPS_S:
            out.append(Issue("error", "epoch_outside_session",
                             f"{k} [{a}, {b}] outside [0, {s.sess_duration}]"))
    for first, second in zip(EPOCH_NAMES, EPOCH_NAMES[1:], strict=False):
        gap = ep[second][0] - ep[first][1]
        if gap < -TIME_EPS_S:
            out.append(Issue("error", "epochs_overlap",
                             f"{first} ends {-gap:.3f} s after {second} starts"))
        elif gap > TIME_EPS_S:
            out.append(Issue("info", "epoch_gap", f"{gap:.3f} s gap between {first} and {second}"))
    if abs(ep["POST"][1] - s.sess_duration) > TIME_EPS_S:
        out.append(Issue("warning", "post_end_not_session_end",
                         f"POST ends at {ep['POST'][1]:.3f} s but sess_duration is "
                         f"{s.sess_duration:.3f} s"))
    return out


def _check_position(s: Session) -> list[Issue]:
    out: list[Issue] = []
    t, xy, x1 = s.position_t, s.position_xy, s.position_1d
    maze = s.epoch("MAZE")
    if len(t) < 2:
        return [Issue("error", "position_too_short", f"only {len(t)} position samples")]
    dt = np.diff(t)
    if np.any(dt <= 0):
        out.append(Issue("error", "position_t_not_monotonic",
                         f"{np.count_nonzero(dt <= 0)} non-increasing timestamp steps"))
    if t[0] < maze[0] - TIME_EPS_S or t[-1] > maze[1] + TIME_EPS_S:
        out.append(Issue("error", "position_outside_maze",
                         f"position spans [{t[0]:.3f}, {t[-1]:.3f}] s but MAZE is "
                         f"[{maze[0]:.3f}, {maze[1]:.3f}] s"))
    med = float(np.median(dt))
    out.append(Issue("info", "position_dt",
                     f"median {med * 1e3:.4f} ms ({1 / med:.3f} Hz), "
                     f"range {dt.min() * 1e3:.3f}-{dt.max() * 1e3:.3f} ms"))
    if abs(med - NOMINAL_POSITION_DT_S) > POSITION_DT_RTOL * NOMINAL_POSITION_DT_S:
        out.append(Issue("warning", "position_dt_unexpected",
                         f"median dt {med * 1e3:.3f} ms, "
                         f"expected ~{NOMINAL_POSITION_DT_S * 1e3} ms"))
    gaps = dt > 1.5 * med
    if gaps.any():
        out.append(Issue("warning", "position_time_gaps",
                         f"{gaps.sum()} gaps longer than 1.5x the median interval "
                         f"(longest {dt.max():.3f} s)"))

    nan2 = np.isnan(xy)
    partial = nan2.any(axis=1) & ~nan2.all(axis=1)
    if partial.any():
        out.append(Issue("warning", "xy_partial_nan",
                         f"{partial.sum()} samples with only x or y NaN"))
    nan2 = nan2.any(axis=1)
    nan1 = np.isnan(x1)
    out.append(Issue("info", "position_nan",
                     f"2-D NaN {nan2.mean():.1%}; 1-D NaN {nan1.mean():.1%} "
                     f"(1-D NaN while 2-D valid: {(nan1 & ~nan2).mean():.1%}); "
                     f"1-D valid for {(~nan1).sum() * med / 60:.1f} of "
                     f"{(maze[1] - maze[0]) / 60:.1f} min"))
    if np.any(~nan1 & nan2):
        out.append(Issue("warning", "1d_valid_where_2d_nan",
                         f"{np.count_nonzero(~nan1 & nan2)} samples have 1-D but no 2-D position"))

    if np.all(nan1):
        out.append(Issue("error", "position_1d_all_nan", "OneDLocation is entirely NaN"))
        return out
    lo, hi = float(np.nanmin(x1)), float(np.nanmax(x1))
    if s.maze_kind == "linear":
        length = s.track_length_m
        out.append(Issue("info", "position_1d_range",
                         f"linear track: 1-D range [{lo:.4f}, {hi:.4f}] m, "
                         f"nominal length {length} m"))
        if lo < -TRACK_LENGTH_TOL_M or (length is not None and hi > length + TRACK_LENGTH_TOL_M):
            out.append(Issue("warning", "position_1d_out_of_track",
                             f"1-D range [{lo:.4f}, {hi:.4f}] m exceeds [0, {length}] m"))
    else:
        out.append(Issue("info", "position_1d_range",
                         f"circular track: 1-D range [{lo:.4f}, {hi:.4f}] m; position wraps at 0 "
                         f"(reward site), implied circumference >= {hi:.3f} m"))
        if lo < -TRACK_LENGTH_TOL_M:
            out.append(Issue("warning", "position_1d_negative", f"1-D minimum {lo:.4f} m < 0"))
    return out


def _check_states(s: Session) -> list[Issue]:
    out: list[Issue] = []
    for name in STATE_NAMES:
        iv = s.states[name]
        if len(iv) == 0:
            out.append(Issue("info", "state_empty", f"{name} has no intervals"))
            continue
        zero = np.where(iv[:, 1] == iv[:, 0])[0]
        neg = np.where(iv[:, 1] < iv[:, 0])[0]
        if zero.size:
            out.append(Issue("warning", "state_zero_length",
                             f"{name} rows {zero.tolist()} have start == end: "
                             f"{iv[zero].tolist()}"))
        if neg.size:
            out.append(Issue("error", "state_negative_length",
                             f"{name} rows {neg.tolist()} have start > end: {iv[neg].tolist()}"))
        if np.any(np.diff(iv[:, 0]) < 0):
            out.append(Issue("warning", "state_unsorted", f"{name} intervals are not sorted"))
        n_self, _ = _interval_overlap_s(iv, iv)
        # Each positive-length row overlaps itself exactly once; anything more is
        # a genuine overlap between two different rows.
        if n_self > np.count_nonzero(iv[:, 1] > iv[:, 0]):
            out.append(Issue("warning", "state_self_overlap",
                             f"{name} intervals overlap each other"))
        outside = (iv[:, 0] < -TIME_EPS_S) | (iv[:, 1] > s.sess_duration + TIME_EPS_S)
        if outside.any():
            out.append(Issue("warning", "state_outside_session",
                             f"{outside.sum()} {name} intervals extend beyond "
                             f"[0, {s.sess_duration:.3f}] s (latest end {iv[:, 1].max():.1f} s)"))
    for i, a in enumerate(STATE_NAMES):
        for b in STATE_NAMES[i + 1:]:
            n, total = _interval_overlap_s(s.states[a], s.states[b])
            if n:
                out.append(Issue("warning", "states_overlap",
                                 f"{a} and {b} overlap in {n} interval pairs, {total:.1f} s total"))
    return out


def validate(session: Session) -> list[Issue]:
    """Run every check on ``session`` and return the issues found."""
    issues: list[Issue] = []
    for check in (_check_spikes, _check_units, _check_epochs, _check_position, _check_states):
        issues.extend(check(session))
    return issues


def errors(issues: list[Issue]) -> list[Issue]:
    return [i for i in issues if i.severity == "error"]


def warnings(issues: list[Issue]) -> list[Issue]:
    return [i for i in issues if i.severity == "warning"]
