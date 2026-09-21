"""Shared fixtures: synthetic MATLAB v7.3 files and cached real sessions."""

from __future__ import annotations

from pathlib import Path

import h5py
import numpy as np
import pytest

from hc11 import config

# ---------------------------------------------------------------------------
# Writing MATLAB-v7.3-shaped HDF5 by hand
# ---------------------------------------------------------------------------


def write_matlab(group: h5py.Group, name: str, matlab_array, cls: str = "double") -> h5py.Dataset:
    """Write ``matlab_array`` the way MATLAB v7.3 would: axes reversed on disk.

    ``matlab_array`` is given in MATLAB orientation (e.g. an n x 2 interval
    matrix). 1-D input is treated as a MATLAB row vector (1 x n).
    """
    a = np.asarray(matlab_array)
    if a.ndim == 1:
        a = a[np.newaxis, :]
    ds = group.create_dataset(name, data=a.T)
    ds.attrs["MATLAB_class"] = np.bytes_(cls)
    return ds


def write_matlab_string(group: h5py.Group, name: str, text: str) -> h5py.Dataset:
    codes = np.frombuffer(text.encode("utf-16-le"), dtype=np.uint16)
    ds = group.create_dataset(name, data=codes.reshape(-1, 1))  # 1 x n in MATLAB
    ds.attrs["MATLAB_class"] = np.bytes_("char")
    ds.attrs["MATLAB_int_decode"] = np.int32(2)
    return ds


def write_matlab_empty(group: h5py.Group, name: str, dims=(0, 0), cls: str = "double"):
    """MATLAB's placeholder for an empty array: the dims, flagged MATLAB_empty."""
    ds = group.create_dataset(name, data=np.asarray(dims, dtype=np.uint64))
    ds.attrs["MATLAB_class"] = np.bytes_(cls)
    ds.attrs["MATLAB_empty"] = np.uint8(1)
    return ds


def make_struct(parent: h5py.Group, name: str) -> h5py.Group:
    g = parent.create_group(name)
    g.attrs["MATLAB_class"] = np.bytes_("struct")
    return g


# ---------------------------------------------------------------------------
# A small, internally consistent synthetic session
# ---------------------------------------------------------------------------

SYNTH_NAME = "Testrat_01012020"


def synthetic_session_arrays(seed: int = 0) -> dict:
    """Arrays for a clean synthetic session, in MATLAB orientation."""
    rng = np.random.default_rng(seed)
    pre, maze, post = (0.0, 100.0), (100.0, 160.0), (160.0, 300.0)
    pyr = np.array([102, 103, 205])
    inter = np.array([104, 210])
    ids = np.concatenate([pyr, inter])
    times = np.sort(rng.uniform(0.01, 299.9, 2000))
    spike_ids = rng.choice(ids, size=times.size)
    spike_ids[: ids.size] = ids  # every unit fires at least once
    dt = 0.0256
    t = np.arange(maze[0], maze[1] - dt / 2, dt)
    x1 = 0.8 + 0.8 * np.sin(2 * np.pi * (t - t[0]) / 10.0)
    xy = np.column_stack([x1, np.full_like(x1, 0.3)])
    xy[50:60] = np.nan
    x1 = x1.copy()
    x1[40:70] = np.nan
    return {
        "SpikeTimes": times[:, None],  # n x 1, like MATLAB
        "SpikeIDs": spike_ids.astype(float)[:, None],
        "PyrIDs": pyr.astype(float)[:, None],
        "IntIDs": inter.astype(float)[:, None],
        "TwoDLocation": xy,  # n x 2
        "OneDLocation": x1[:, None],  # n x 1
        "TimeStamps": t[None, :],  # 1 x n (as documented)
        "MazeType": "1.6m Linear Maze",
        "PREEpoch": np.array([pre]),
        "MazeEpoch": np.array([maze]),
        "POSTEpoch": np.array([post]),
        "sessDuration": np.array([[post[1]]]),
        "Wake": np.array([[0.0, 20.0], [100.0, 160.0]]),
        "Drowsy": np.array([[20.0, 30.0]]),
        "NREM": np.array([[30.0, 80.0], [170.0, 250.0]]),
        "Intermediate": np.array([[80.0, 82.0]]),
        "REM": np.array([[82.0, 95.0], [250.0, 290.0]]),
    }


def write_sessinfo(path: Path, arrays: dict, empty_fields: tuple[str, ...] = ()) -> Path:
    """Write a synthetic ``sessInfo.mat`` with the real file's layout."""
    layout = {
        "Spikes": ("SpikeTimes", "SpikeIDs", "PyrIDs", "IntIDs"),
        "Position": ("TwoDLocation", "OneDLocation", "TimeStamps", "MazeType"),
        "Epochs": ("PREEpoch", "MazeEpoch", "POSTEpoch", "sessDuration",
                   "Wake", "Drowsy", "NREM", "Intermediate", "REM"),
    }
    with h5py.File(path, "w", userblock_size=512) as f:
        root = make_struct(f, "sessInfo")
        for struct, fields in layout.items():
            g = make_struct(root, struct)
            for name in fields:
                if name in empty_fields:
                    write_matlab_empty(g, name, dims=(0, 2))
                elif name == "MazeType":
                    write_matlab_string(g, name, arrays[name])
                else:
                    write_matlab(g, name, arrays[name])
    with open(path, "r+b") as fh:  # MATLAB's human-readable header
        fh.write(b"MATLAB 7.3 MAT-file, Platform: TEST, HDF5 schema 1.00 .".ljust(128))
    return path


@pytest.fixture
def synthetic_arrays() -> dict:
    return synthetic_session_arrays()


@pytest.fixture
def synthetic_file(tmp_path, synthetic_arrays) -> Path:
    return write_sessinfo(tmp_path / f"{SYNTH_NAME}_sessInfo.mat", synthetic_arrays)


# ---------------------------------------------------------------------------
# Real data (skipped when absent)
# ---------------------------------------------------------------------------


def _real_data_available() -> bool:
    try:
        config.data_dir()
    except FileNotFoundError:
        return False
    return True


requires_data = pytest.mark.skipif(
    not _real_data_available(), reason="hc-11 sessInfo files not found (see README)"
)


@pytest.fixture(scope="session")
def real_sessions():
    """Lazily load and cache real sessions for the whole test run."""
    from hc11.io import load_session

    cache: dict = {}

    def get(name: str):
        if name not in cache:
            cache[name] = load_session(name)
        return cache[name]

    return get
