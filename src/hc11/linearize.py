"""Our own linearisation of the 2-D position, for a diagnostic only.

The replication base is and remains the authors' ``OneDLocation``: it is their
linearisation, and reproducing their analysis means using it. This module exists
to answer one question raised in decision D5.4 — whether the place-cell gap comes
from the *measurement base* rather than from our criteria, since ``OneDLocation``
is defined for only 7-35 % of each MAZE epoch while 2-D tracking covers 46-98 %.

The coordinate produced here has its own origin and scale and is **not** aligned
to ``OneDLocation``. That is deliberate: place-cell counts do not depend on where
the origin sits, and calibrating against the authors' values would reintroduce
the very restriction the diagnostic is trying to lift.

* **Linear tracks**: project the 2-D position onto its first principal axis. This
  is exactly what the authors' linearisation does — on every linear session
  ``OneDLocation`` correlates with that projection at |r| = 1.0000
  (``tests/test_io.py::test_linear_1d_is_projection_of_2d``).
* **Circular tracks**: arc length about the fitted centre, wrapped into
  ``[0, 2*pi*R)`` with ``R`` the median radius.
"""

from __future__ import annotations

import dataclasses

import numpy as np

from hc11.io import Session


def linearize_2d(session: Session) -> np.ndarray:
    """A 1-D coordinate (metres) wherever 2-D tracking is valid, NaN elsewhere.

    Uses the maze type to choose the projection; raises for an unknown type.
    """
    xy = session.position_xy
    valid = ~np.isnan(xy).any(axis=1)
    out = np.full(len(xy), np.nan)
    if valid.sum() < 3:
        return out

    if session.maze_kind == "linear":
        centred = xy[valid] - xy[valid].mean(axis=0)
        axis = np.linalg.svd(centred, full_matrices=False)[2][0]
        projection = centred @ axis
        # Orient so the coordinate increases along the same direction as the
        # authors' where both are defined, then put the origin at the low end.
        both = valid & ~np.isnan(session.position_1d)
        if both.sum() > 10:
            theirs = session.position_1d[both]
            ours = (xy[both] - xy[valid].mean(axis=0)) @ axis
            if np.corrcoef(ours, theirs)[0, 1] < 0:
                projection = -projection
        out[valid] = projection - projection.min()
        return out

    if session.maze_kind == "circular":
        centre = np.median(xy[valid], axis=0)
        centred = xy[valid] - centre
        radius = float(np.median(np.hypot(centred[:, 0], centred[:, 1])))
        angle = np.arctan2(centred[:, 1], centred[:, 0])
        out[valid] = np.mod(angle, 2 * np.pi) * radius
        return out

    raise ValueError(f"{session.name}: unsupported maze kind {session.maze_kind!r}")


def with_own_linearization(session: Session) -> Session:
    """A copy of ``session`` whose ``position_1d`` is :func:`linearize_2d`.

    Everything downstream (``Track``, running mask, laps, rate maps) then works
    unchanged on the wider base, because they all read ``position_1d``.
    """
    x = linearize_2d(session)
    x.setflags(write=False)
    return dataclasses.replace(session, position_1d=x)


def coverage(session: Session) -> dict[str, float]:
    """Fraction of MAZE samples covered by each linearisation, for reporting."""
    n = len(session.position_t)
    theirs = int(np.sum(~np.isnan(session.position_1d)))
    ours = int(np.sum(~np.isnan(linearize_2d(session))))
    return {
        "n_samples": n,
        "authors_frac": theirs / n,
        "own_frac": ours / n,
        "authors_min": theirs * session.position_dt / 60,
        "own_min": ours * session.position_dt / 60,
    }
