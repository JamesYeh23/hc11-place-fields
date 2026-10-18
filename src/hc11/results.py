"""Writing per-session result tables without losing the sessions you did not run.

Every analysis script writes one CSV covering the sessions it was asked for.
Naively that means running one session rewrites the table and deletes the other
seven, which is both a data-loss footgun and the reason a long run could not be
split across several invocations.

:func:`write_session_table` merges instead: rows for the sessions just computed
replace their previous versions, rows for every other session are carried
through untouched. Running the pipeline session by session therefore
accumulates, and re-running one session refreshes only that session.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd

#: Column identifying the session in a results table.
SESSION_COLUMN = "session"


def write_session_table(
    path: str | Path,
    rows,
    sessions=None,
    *,
    round_to: int | None = None,
    sort: bool = True,
) -> pd.DataFrame:
    """Write ``rows`` to ``path``, keeping rows for sessions not in this run.

    ``sessions`` is the set of sessions this run computed; their old rows are
    dropped before the new ones are added. When omitted it is taken from the
    new rows themselves, which is right whenever a run that produced no rows
    for a session should leave that session's old rows alone -- and wrong if a
    run is meant to *clear* a session, which no caller wants.

    The merged table is returned so a caller can report on the whole of it
    rather than only on what it just computed.
    """
    path = Path(path)
    new = pd.DataFrame(rows)
    if new.empty and not path.exists():
        new.to_csv(path, index=False)
        return new

    if sessions is None:
        sessions = set(new[SESSION_COLUMN].unique()) if SESSION_COLUMN in new else set()
    else:
        sessions = set(sessions)

    frames = []
    if path.exists():
        old = pd.read_csv(path)
        if SESSION_COLUMN in old.columns:
            # With no sessions named, nothing is being replaced, so everything
            # already there is kept -- an empty run must not empty the table.
            kept = old[~old[SESSION_COLUMN].isin(sessions)] if sessions else old
            if not kept.empty:
                frames.append(kept)
        # A table with no session column cannot be merged by session, so the
        # only available meaning is to replace it; silently merging is worse.
    if not new.empty:
        frames.append(new)
    merged = pd.concat(frames, ignore_index=True) if frames else new

    if sort and SESSION_COLUMN in merged.columns:
        order = [c for c in (SESSION_COLUMN, "epoch", "variant", "ranking", "shuffle")
                 if c in merged.columns]
        merged = merged.sort_values(order, kind="stable").reset_index(drop=True)
    if round_to is not None:
        merged = merged.round(round_to)
    path.parent.mkdir(parents=True, exist_ok=True)
    merged.to_csv(path, index=False)
    return merged
