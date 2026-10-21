"""Splitting a place-cell ensemble into halves, and the control for doing so.

Step 9b asks the paper's Fig. 4D question — do slow-firing cells increase their
contribution to sequences from PRE to POST while fast-firing cells do not — by
rebuilding the decoding template from one half of the ensemble at a time. That
needs three things kept separate and individually testable: a firing rate
measured over the off-line (sleep) epochs, an equal split of the cells by that
rate, and equal splits made at random, which is the control. A difference
between two halves of an ensemble is not evidence until the spread of
differences between *arbitrary* halves of the same ensemble is known.
"""

from __future__ import annotations

import numpy as np


def rate_in_intervals(session, cluster_id: int, intervals: np.ndarray) -> float:
    """Firing rate of one cluster inside ``intervals`` (Hz), or NaN if they are empty.

    ``intervals`` is ``(n, 2)`` of [start, end] seconds. Overlapping intervals
    would double-count both spikes and time; the caller is expected to pass the
    disjoint intervals the rest of the pipeline uses.
    """
    intervals = np.atleast_2d(np.asarray(intervals, dtype=float))
    if intervals.size == 0:
        return float("nan")
    duration = float(np.sum(intervals[:, 1] - intervals[:, 0]))
    if duration <= 0:
        return float("nan")
    times = session.spikes_for(int(cluster_id))
    lo = np.searchsorted(times, intervals[:, 0], side="left")
    hi = np.searchsorted(times, intervals[:, 1], side="right")
    return float(np.sum(hi - lo) / duration)


def equal_halves(values: np.ndarray, cluster_ids: np.ndarray) -> dict:
    """Split ``cluster_ids`` into equal low and high halves by ``values``.

    The paper's subgroups are equal, so with an odd number of cells the median
    cell is **dropped** rather than assigned to a side: assigning it would make
    one group larger and would do so by a rule nobody could check. Ties are
    broken by cluster ID, which is arbitrary but deterministic.

    Returns ``{"low": ids, "high": ids, "dropped": ids}``.
    """
    values = np.asarray(values, dtype=float)
    cluster_ids = np.asarray(cluster_ids)
    if values.shape != cluster_ids.shape:
        raise ValueError(
            f"{values.size} values for {cluster_ids.size} cells")
    order = np.lexsort((cluster_ids, values))
    ranked = cluster_ids[order]
    half = len(ranked) // 2
    dropped = ranked[half:half + 1] if len(ranked) % 2 else ranked[:0]
    return {
        "low": np.sort(ranked[:half]),
        "high": np.sort(ranked[len(ranked) - half:]),
        "dropped": np.sort(dropped),
    }


def random_halves(cluster_ids: np.ndarray, rng) -> dict:
    """Two equal halves drawn at random, in the same shape as :func:`equal_halves`.

    The control for any claim of the form "these two halves of the ensemble
    behave differently": it answers how different two halves are when the split
    carries no information at all.
    """
    cluster_ids = np.asarray(cluster_ids)
    shuffled = rng.permutation(cluster_ids)
    half = len(shuffled) // 2
    dropped = shuffled[half:half + 1] if len(shuffled) % 2 else shuffled[:0]
    return {
        "low": np.sort(shuffled[:half]),
        "high": np.sort(shuffled[len(shuffled) - half:]),
        "dropped": np.sort(dropped),
    }


def split_agreement(a: dict, b: dict) -> float:
    """Fraction of cells that land in the same half under two splits.

    Symmetric in the labels: a split that is the exact mirror of another agrees
    completely, because which half is called "low" is a naming choice.
    """
    shared = np.intersect1d(
        np.concatenate([a["low"], a["high"]]),
        np.concatenate([b["low"], b["high"]]),
    )
    if shared.size == 0:
        return float("nan")
    in_a_low = np.isin(shared, a["low"])
    in_b_low = np.isin(shared, b["low"])
    same = float(np.mean(in_a_low == in_b_low))
    return max(same, 1.0 - same)


def finite_in_both(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    """Mask of entries finite in *both* arrays.

    Two halves of an ensemble are compared on the same events, so an event one
    half cannot score — a degenerate posterior gives a shuffle distribution of
    zero width and an undefined z-score — has to leave the comparison for both
    halves, not just for the half that failed on it. Dropping it from one side
    only would compare two different event sets and call the difference a
    result.
    """
    a, b = np.asarray(a, dtype=float), np.asarray(b, dtype=float)
    if a.shape != b.shape:
        raise ValueError(f"{a.size} and {b.size} values cannot be paired")
    return np.isfinite(a) & np.isfinite(b)


def empirical_p(observed: float, draws) -> float:
    """Two-sided p for ``observed`` against a null of ``draws``, NaN-safe.

    ``(1 + #{|draw| >= |observed|}) / (n + 1)``: the plus-one keeps p away from
    zero, which an empirical null can never justify.

    Non-finite draws are discarded, and a non-finite ``observed`` returns NaN
    rather than a p-value. Without that guard a NaN observation compares false
    against every draw and comes back as the smallest p the null can express —
    a missing number reported as the strongest possible result.
    """
    draws = np.asarray(list(draws), dtype=float)
    draws = np.abs(draws[np.isfinite(draws)])
    if not np.isfinite(observed) or draws.size == 0:
        return float("nan")
    return float((1 + int(np.sum(draws >= abs(observed)))) / (draws.size + 1))
