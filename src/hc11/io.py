"""Loader for CRCNS hc-11 ``<session>_sessInfo.mat`` files.

The files are MATLAB v7.3, i.e. HDF5 with a 512-byte MATLAB userblock. They are
read with :mod:`h5py` in two layers:

1. :func:`read_matlab` — a *generic* reader that turns any MATLAB v7.3 object
   into plain Python: structs → ``dict``, numeric arrays → ``ndarray`` in
   **MATLAB orientation**, char arrays → ``str``, cell arrays / object references
   → ``list``, and MATLAB empties → zero-size arrays of the declared shape.
2. :func:`load_session` — the *hc-11-specific* assembly that normalises every
   field (vectors to 1-D, interval matrices to ``(n, 2)``, IDs to ``int64``) and
   returns an immutable :class:`Session`.

Orientation rule: HDF5 stores MATLAB arrays with their axes reversed, so the
MATLAB array is always ``h5py_array.T``. This is applied unconditionally rather
than by guessing from the shape, which would silently mis-read a ``2 × 2``
interval matrix. The actual on-disk layouts are documented in
``docs/data_schema.md``.

The raw file is opened read-only and never modified.
"""

from __future__ import annotations

import re
from collections.abc import Mapping
from dataclasses import dataclass, field
from pathlib import Path
from types import MappingProxyType
from typing import Any, Literal

import h5py
import numpy as np

from hc11 import config

#: Names of the three recording epochs, in chronological order.
EPOCH_NAMES: tuple[str, ...] = ("PRE", "MAZE", "POST")
#: Mapping from our epoch names to the field names used in the .mat file.
_EPOCH_FIELDS: dict[str, str] = {"PRE": "PREEpoch", "MAZE": "MazeEpoch", "POST": "POSTEpoch"}
#: The five scored behavioural states, as named in the .mat file.
STATE_NAMES: tuple[str, ...] = ("Wake", "Drowsy", "NREM", "Intermediate", "REM")

_SESSINFO_SUFFIX = "_sessInfo.mat"
UnitKind = Literal["pyr", "int", "all"]


# ---------------------------------------------------------------------------
# Generic MATLAB v7.3 helpers
# ---------------------------------------------------------------------------


def _attr_str(obj: h5py.HLObject, name: str) -> str | None:
    """Read a string-valued HDF5 attribute (MATLAB stores these as bytes)."""
    if name not in obj.attrs:
        return None
    value = obj.attrs[name]
    if isinstance(value, bytes):
        return value.decode("ascii")
    if isinstance(value, np.ndarray) and value.dtype.kind == "S":
        return b"".join(value.ravel()).decode("ascii")
    return str(value)


def matlab_class(obj: h5py.HLObject) -> str | None:
    """The ``MATLAB_class`` attribute of an HDF5 object (e.g. ``'double'``)."""
    return _attr_str(obj, "MATLAB_class")


def is_matlab_empty(dataset: h5py.Dataset) -> bool:
    """True if ``dataset`` is MATLAB's placeholder for an empty array.

    MATLAB v7.3 cannot write a zero-size HDF5 dataset, so an empty array is
    stored as a tiny dataset holding the array's *dimensions*, flagged with a
    ``MATLAB_empty`` attribute. Its values must never be used as data.
    """
    if "MATLAB_empty" not in dataset.attrs:
        return False
    flag = np.asarray(dataset.attrs["MATLAB_empty"])
    return bool(flag.ravel()[0]) if flag.size else True


def _empty_from_placeholder(dataset: h5py.Dataset) -> np.ndarray:
    """Zero-size array with the MATLAB shape recorded in an empty placeholder."""
    dims = np.asarray(dataset[()]).ravel()
    shape = tuple(int(d) for d in dims) if dims.size >= 2 else (0, 0)
    if 0 not in shape:  # defensive: a placeholder must describe an empty array
        shape = (0, 0)
    return np.empty(shape, dtype=np.float64)


def decode_matlab_string(codes: np.ndarray) -> str | list[str]:
    """Decode a MATLAB char array (stored as uint16 UTF-16 code units).

    ``codes`` must already be in MATLAB orientation. A single row (``1 × n``)
    decodes to one ``str``; an ``m × n`` char matrix decodes to ``m`` strings
    with MATLAB's right-padding stripped.
    """
    arr = np.asarray(codes)
    if arr.size == 0:
        return ""
    if arr.dtype.kind not in "ui":
        raise TypeError(f"MATLAB char data must be integer code units, got {arr.dtype}")
    if arr.ndim == 1:
        arr = arr[np.newaxis, :]
    rows = [
        np.asarray(row, dtype=np.uint16).tobytes().decode("utf-16-le").rstrip("\x00")
        for row in arr
    ]
    return rows[0] if len(rows) == 1 else [r.rstrip() for r in rows]


def read_matlab(obj: h5py.HLObject, file: h5py.File | None = None) -> Any:
    """Recursively convert a MATLAB v7.3 HDF5 object to plain Python.

    Numeric arrays are returned in MATLAB orientation (HDF5 axes reversed).
    Object-reference datasets (cell arrays, and some struct arrays) are
    dereferenced against ``file``.
    """
    file = file if file is not None else obj.file
    if isinstance(obj, h5py.Group):
        return {key: read_matlab(obj[key], file) for key in obj}

    if not isinstance(obj, h5py.Dataset):
        raise TypeError(f"Unsupported HDF5 object at {obj.name!r}: {type(obj)!r}")

    if is_matlab_empty(obj):
        return "" if matlab_class(obj) == "char" else _empty_from_placeholder(obj)

    if obj.dtype == h5py.ref_dtype:
        refs = np.asarray(obj[()]).T  # MATLAB orientation of the cell array
        return [read_matlab(file[ref], file) for ref in refs.ravel(order="F")]

    data = np.asarray(obj[()]).T
    if matlab_class(obj) == "char":
        return decode_matlab_string(data)
    if matlab_class(obj) == "logical":
        return data.astype(bool)
    return data


# ---------------------------------------------------------------------------
# Field normalisation
# ---------------------------------------------------------------------------


def as_vector(arr: np.ndarray, name: str = "array") -> np.ndarray:
    """Flatten a MATLAB row/column vector (``1 × n``, ``n × 1``) to 1-D.

    Raises ``ValueError`` for anything that is genuinely two-dimensional, so a
    matrix is never silently flattened.
    """
    a = np.asarray(arr)
    if a.ndim > 2 or (a.ndim == 2 and min(a.shape) > 1):
        raise ValueError(f"{name}: expected a vector, got shape {a.shape}")
    return a.reshape(-1)


def as_intervals(arr: np.ndarray, name: str = "intervals") -> np.ndarray:
    """Normalise a MATLAB ``[start, end]`` matrix to shape ``(n_intervals, 2)``.

    ``arr`` must be in MATLAB orientation (as returned by :func:`read_matlab`):
    ``n × 2``, or a ``1 × 2`` / ``2``-element vector for a single interval.
    Empty input returns ``np.empty((0, 2))``.
    """
    a = np.asarray(arr, dtype=np.float64)
    if a.size == 0:
        return np.empty((0, 2), dtype=np.float64)
    if a.ndim == 1:
        if a.size != 2:
            raise ValueError(f"{name}: a 1-D interval needs 2 elements, got {a.size}")
        return a.reshape(1, 2)
    if a.ndim == 2 and a.shape[1] == 2:
        return a
    raise ValueError(
        f"{name}: expected an n x 2 [start, end] matrix in MATLAB orientation, got {a.shape}"
    )


def as_int_ids(arr: np.ndarray, name: str = "ids") -> np.ndarray:
    """Cast MATLAB double-precision IDs to ``int64`` after checking they are whole."""
    a = as_vector(arr, name)
    if a.size == 0:
        return np.empty(0, dtype=np.int64)
    if not np.all(np.isfinite(a)):
        raise ValueError(f"{name}: contains non-finite values")
    rounded = np.round(a)
    if not np.array_equal(a, rounded):
        bad = a[a != rounded][:5]
        raise ValueError(f"{name}: non-integer IDs, e.g. {bad.tolist()}")
    return rounded.astype(np.int64)


def shank_of(cluster_id: int | np.ndarray) -> int | np.ndarray:
    """Shank number of a cluster: the hundreds digit of its ID (``floor(ID / 100)``)."""
    return np.asarray(cluster_id) // 100 if np.ndim(cluster_id) else int(cluster_id) // 100


def parse_animal(session_name: str) -> str:
    """``'Achilles_10252013'`` → ``'Achilles'``."""
    match = re.fullmatch(r"([A-Za-z]+)_(\d{8})", session_name)
    if match is None:
        raise ValueError(f"Session name {session_name!r} is not of the form Animal_MMDDYYYY")
    return match.group(1)


def _readonly(a: np.ndarray) -> np.ndarray:
    a = np.ascontiguousarray(a)
    a.setflags(write=False)
    return a


# ---------------------------------------------------------------------------
# Session
# ---------------------------------------------------------------------------


@dataclass(frozen=True, eq=False, repr=False)
class Session:
    """One hc-11 recording session, as stored in its ``sessInfo.mat``.

    All times are in seconds on the session clock (PRE starts at 0). All
    positions are in **metres**, as stored; use :attr:`position_xy_cm` /
    :attr:`position_1d_cm` for centimetres. Arrays are read-only.

    Values are exactly what the file contains: nothing is clipped, dropped or
    interpolated here. Use :func:`hc11.validate.validate` to find anomalies.
    """

    name: str
    animal: str
    path: Path
    spike_times: np.ndarray  # (n_spikes,) float64, s
    spike_ids: np.ndarray  # (n_spikes,) int64, cluster ID per spike
    pyr_ids: np.ndarray  # (n_pyr,) int64, sorted
    int_ids: np.ndarray  # (n_int,) int64, sorted
    shank_of_cluster: Mapping[int, int]  # cluster ID -> shank number
    position_t: np.ndarray  # (n_samples,) float64, s
    position_xy: np.ndarray  # (n_samples, 2) float64, m; NaN = tracking lost
    position_1d: np.ndarray  # (n_samples,) float64, m; NaN outside track running
    maze_type: str
    epochs: Mapping[str, np.ndarray]  # "PRE"/"MAZE"/"POST" -> (2,) [start, end], s
    states: Mapping[str, np.ndarray]  # state name -> (n, 2) [start, end], s
    sess_duration: float  # s
    _spike_order: np.ndarray = field(init=False, repr=False)  # spikes sorted by (id, time)
    _id_bounds: Mapping[int, tuple[int, int]] = field(init=False, repr=False)

    def __post_init__(self) -> None:
        order = np.lexsort((self.spike_times, self.spike_ids))
        sorted_ids = self.spike_ids[order]
        ids, first = np.unique(sorted_ids, return_index=True)
        last = np.append(first[1:], len(sorted_ids))
        bounds = {int(i): (int(a), int(b)) for i, a, b in zip(ids, first, last, strict=True)}
        object.__setattr__(self, "_spike_order", _readonly(order))
        object.__setattr__(self, "_id_bounds", MappingProxyType(bounds))

    # --- derived descriptors -------------------------------------------------

    @property
    def maze_kind(self) -> Literal["linear", "circular"]:
        """``'linear'`` or ``'circular'``, parsed from :attr:`maze_type`."""
        lowered = self.maze_type.lower()
        if "linear" in lowered:
            return "linear"
        if "circular" in lowered:
            return "circular"
        raise ValueError(f"Unrecognised maze type {self.maze_type!r}")

    @property
    def track_length_m(self) -> float | None:
        """Nominal linear-track length from the maze name (``'1.6m Linear Maze'`` → 1.6).

        ``None`` for the circular maze, whose name carries no length.
        """
        match = re.search(r"(\d+(?:\.\d+)?)\s*m\b", self.maze_type)
        return float(match.group(1)) if match and self.maze_kind == "linear" else None

    @property
    def position_xy_cm(self) -> np.ndarray:
        """Two-dimensional position in centimetres (a new array)."""
        return self.position_xy * 100.0

    @property
    def position_1d_cm(self) -> np.ndarray:
        """Linearised position in centimetres (a new array)."""
        return self.position_1d * 100.0

    @property
    def position_dt(self) -> float:
        """Median position sampling interval in seconds."""
        return float(np.median(np.diff(self.position_t)))

    @property
    def cluster_ids(self) -> np.ndarray:
        """Every cluster ID that has at least one spike, sorted."""
        return np.fromiter(self._id_bounds.keys(), dtype=np.int64)

    def epoch(self, name: str) -> np.ndarray:
        """``[start, end]`` of ``'PRE'``, ``'MAZE'`` or ``'POST'``."""
        try:
            return self.epochs[name.upper()]
        except KeyError:
            raise KeyError(f"Unknown epoch {name!r}; expected one of {EPOCH_NAMES}") from None

    def epoch_duration(self, name: str) -> float:
        start, end = self.epoch(name)
        return float(end - start)

    # --- convenience queries ------------------------------------------------

    def units(self, kind: UnitKind = "all") -> np.ndarray:
        """Cluster IDs of putative pyramidal cells, interneurons, or both."""
        if kind == "pyr":
            return self.pyr_ids
        if kind == "int":
            return self.int_ids
        if kind == "all":
            return np.union1d(self.pyr_ids, self.int_ids)
        raise ValueError(f"kind must be 'pyr', 'int' or 'all', got {kind!r}")

    def spikes_for(self, cluster_id: int) -> np.ndarray:
        """Sorted spike times of one cluster (empty if it never fired)."""
        start, stop = self._id_bounds.get(int(cluster_id), (0, 0))
        return self.spike_times[self._spike_order[start:stop]]

    def spikes_in(
        self, interval: np.ndarray | tuple[float, float], cluster_id: int | None = None
    ) -> tuple[np.ndarray, np.ndarray]:
        """Spikes with ``start <= t < end``, as ``(times, ids)``.

        If ``cluster_id`` is given, only that cluster's spikes are returned.
        """
        start, end = float(interval[0]), float(interval[1])
        if cluster_id is not None:
            times = self.spikes_for(cluster_id)
            lo, hi = np.searchsorted(times, [start, end], side="left")
            return times[lo:hi], np.full(hi - lo, int(cluster_id), dtype=np.int64)
        lo, hi = np.searchsorted(self.spike_times, [start, end], side="left")
        return self.spike_times[lo:hi], self.spike_ids[lo:hi]

    def __repr__(self) -> str:
        return (
            f"Session({self.name!r}, {self.maze_type!r}, "
            f"{len(self.pyr_ids)} pyr / {len(self.int_ids)} int, "
            f"{len(self.spike_times):,} spikes, {self.sess_duration / 3600:.2f} h)"
        )


# ---------------------------------------------------------------------------
# Public entry points
# ---------------------------------------------------------------------------


def list_sessions(data_root: str | Path | None = None) -> list[str]:
    """Names of all sessions with a ``*_sessInfo.mat`` under ``data_root``.

    Searches recursively, so it works both for the flat
    ``NoveltySessInfoMatFiles/`` layout and for per-session folders produced by
    unpacking the full ``X.tar.gz`` archives. Defaults to
    :func:`hc11.config.data_dir`.
    """
    root = Path(data_root) if data_root is not None else config.data_dir()
    names = {p.name[: -len(_SESSINFO_SUFFIX)] for p in root.rglob(f"*{_SESSINFO_SUFFIX}")}
    known = [s for s in config.SESSIONS if s in names]
    return known + sorted(names - set(config.SESSIONS))


def _find_session_file(session: str | Path) -> Path:
    candidate = Path(session)
    if candidate.suffix == ".mat":
        if candidate.is_file():
            return candidate
        raise FileNotFoundError(f"No such file: {candidate}")
    return config.session_path(str(session))


def load_session(path_or_name: str | Path) -> Session:
    """Load one ``sessInfo.mat`` into a :class:`Session`.

    ``path_or_name`` may be a path to the ``.mat`` file or a bare session name
    such as ``'Achilles_10252013'`` (resolved through :mod:`hc11.config`).
    """
    path = _find_session_file(path_or_name)
    name = path.name.removesuffix(_SESSINFO_SUFFIX) if path.suffix == ".mat" else path.stem

    with h5py.File(path, "r") as f:
        if "sessInfo" not in f:
            raise KeyError(f"{path.name}: no top-level 'sessInfo' struct (found {list(f)})")
        raw = read_matlab(f["sessInfo"])

    def get(struct: str, key: str) -> Any:
        try:
            return raw[struct][key]
        except KeyError:
            raise KeyError(f"{path.name}: missing field sessInfo.{struct}.{key}") from None

    spike_times = as_vector(get("Spikes", "SpikeTimes"), "SpikeTimes").astype(np.float64)
    spike_ids = as_int_ids(get("Spikes", "SpikeIDs"), "SpikeIDs")
    if spike_times.shape != spike_ids.shape:
        raise ValueError(
            f"{path.name}: {spike_times.size} spike times but {spike_ids.size} spike IDs"
        )
    pyr_ids = np.sort(as_int_ids(get("Spikes", "PyrIDs"), "PyrIDs"))
    int_ids = np.sort(as_int_ids(get("Spikes", "IntIDs"), "IntIDs"))
    all_ids = np.union1d(np.unique(spike_ids), np.union1d(pyr_ids, int_ids))
    shank_map = MappingProxyType({int(c): int(c) // 100 for c in all_ids})

    position_t = as_vector(get("Position", "TimeStamps"), "TimeStamps").astype(np.float64)
    xy = np.asarray(get("Position", "TwoDLocation"), dtype=np.float64)
    if xy.ndim != 2 or xy.shape[1] != 2:
        raise ValueError(f"{path.name}: TwoDLocation should be n x 2, got {xy.shape}")
    position_1d = as_vector(get("Position", "OneDLocation"), "OneDLocation").astype(np.float64)
    if not (len(position_t) == len(xy) == len(position_1d)):
        raise ValueError(
            f"{path.name}: position lengths disagree: TimeStamps={len(position_t)}, "
            f"TwoDLocation={len(xy)}, OneDLocation={len(position_1d)}"
        )
    maze_type = get("Position", "MazeType")
    if not isinstance(maze_type, str):
        raise TypeError(f"{path.name}: MazeType decoded to {type(maze_type).__name__}, not str")

    epochs = {}
    for ours, theirs in _EPOCH_FIELDS.items():
        iv = as_intervals(get("Epochs", theirs), theirs)
        if iv.shape != (1, 2):
            raise ValueError(f"{path.name}: {theirs} should hold one interval, got {iv.shape}")
        epochs[ours] = _readonly(iv[0])
    states = {s: _readonly(as_intervals(get("Epochs", s), s)) for s in STATE_NAMES}
    sess_duration = float(as_vector(get("Epochs", "sessDuration"), "sessDuration")[0])

    return Session(
        name=name,
        animal=parse_animal(name),
        path=path,
        spike_times=_readonly(spike_times),
        spike_ids=_readonly(spike_ids),
        pyr_ids=_readonly(pyr_ids),
        int_ids=_readonly(int_ids),
        shank_of_cluster=shank_map,
        position_t=_readonly(position_t),
        position_xy=_readonly(xy),
        position_1d=_readonly(position_1d),
        maze_type=maze_type.strip(),
        epochs=MappingProxyType(epochs),
        states=MappingProxyType(states),
        sess_duration=sess_duration,
    )
