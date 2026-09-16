"""Checks that the raw hc-11 data is present and complete.

These are skipped automatically when the data has not been downloaded, so the
suite still passes on a fresh clone.
"""

from __future__ import annotations

import pytest

from hc11 import config

pytestmark = pytest.mark.requires_data

#: HDF5 superblock signature. In a MATLAB v7.3 file it sits after a 512-byte
#: userblock that holds a human-readable description string.
_HDF5_MAGIC = b"\x89HDF\r\n\x1a\n"
_USERBLOCK_BYTES = 512


@pytest.fixture(scope="module")
def data_directory():
    try:
        return config.data_dir()
    except FileNotFoundError as exc:
        pytest.skip(str(exc).splitlines()[0])


def test_all_eight_sessions_present(data_directory):
    missing = [s for s in config.SESSIONS if s not in config.available_sessions()]
    assert not missing, f"missing session files: {missing}"


def test_no_unexpected_session_files(data_directory):
    unexpected = [s for s in config.available_sessions() if s not in config.SESSIONS]
    assert not unexpected, f"unexpected session files: {unexpected}"


@pytest.mark.parametrize("session", config.SESSIONS)
def test_session_file_is_matlab_v73_hdf5(data_directory, session):
    """v7.3 .mat files are HDF5; a v7 file here would break the h5py loader."""
    path = config.session_path(session)
    assert path.stat().st_size > 1_000_000, "session file is suspiciously small"
    with open(path, "rb") as fh:
        userblock = fh.read(_USERBLOCK_BYTES)
        magic = fh.read(len(_HDF5_MAGIC))
    description = userblock[:128].decode("latin-1")
    assert "MATLAB 7.3" in description, (
        f"{path.name} does not declare itself a MATLAB v7.3 file "
        f"(header says: {description.strip()!r}); scipy.io.loadmat may be needed instead."
    )
    assert magic == _HDF5_MAGIC, (
        f"{path.name} has no HDF5 superblock at offset {_USERBLOCK_BYTES}; "
        "the loader's h5py assumption would not hold."
    )
