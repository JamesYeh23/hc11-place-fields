"""Tests for resumable extraction of a .eeg from a CRCNS tarball.

The real archives are 9 GB and extraction is measured in minutes, so these
build a small tarball with the same shape -- one member, a path of
``<session>/<session>.eeg`` -- and check the parts that would corrupt 12 GB of
data if they were wrong: the declared size, resuming at the right offset, and
refusing to resume onto something that is not a prefix.
"""

from __future__ import annotations

import io
import tarfile

import pytest

from hc11 import eeg_archive

SESSION = "Rat_01012000"
MEMBER = f"{SESSION}/{SESSION}.eeg"


def make_archive(path, payload: bytes, name: str = MEMBER, extra: bool = True):
    with tarfile.open(path, mode="w:gz") as tar:
        if extra:  # the real archives carry the .xml alongside the recording
            xml = b"<parameters/>"
            info = tarfile.TarInfo(f"{SESSION}/{SESSION}.xml")
            info.size = len(xml)
            tar.addfile(info, io.BytesIO(xml))
        info = tarfile.TarInfo(name)
        info.size = len(payload)
        tar.addfile(info, io.BytesIO(payload))
    return path


@pytest.fixture
def payload():
    return bytes(range(256)) * 400  # 102,400 bytes


@pytest.fixture
def archive(tmp_path, payload):
    return make_archive(tmp_path / "session_eeg.tar.gz", payload)


def test_member_info_reads_the_declared_size(archive, payload):
    info = eeg_archive.member_info(archive)
    assert info.name == MEMBER
    assert info.filename == f"{SESSION}.eeg"
    assert info.size == len(payload)


def test_member_info_raises_when_there_is_no_eeg(tmp_path, payload):
    archive = make_archive(tmp_path / "wrong.tar.gz", payload, name=f"{SESSION}/notes.txt",
                           extra=False)
    with pytest.raises(eeg_archive.ArchiveError, match="no .eeg member"):
        eeg_archive.member_info(archive)


def test_extract_writes_the_whole_member(archive, tmp_path, payload):
    dest = tmp_path / "out" / f"{SESSION}.eeg"
    progress = eeg_archive.extract(archive, dest)
    assert progress.complete
    assert progress.bytes_written == len(payload)
    assert dest.read_bytes() == payload


def test_extract_resumes_from_a_partial_file(archive, tmp_path, payload):
    dest = tmp_path / f"{SESSION}.eeg"
    dest.write_bytes(payload[:30000])
    progress = eeg_archive.extract(archive, dest)
    assert progress.bytes_at_start == 30000
    assert progress.bytes_written == len(payload) - 30000
    assert dest.read_bytes() == payload


def test_a_deadline_stops_early_and_the_next_call_finishes(archive, tmp_path, payload):
    dest = tmp_path / f"{SESSION}.eeg"
    first = eeg_archive.extract(archive, dest, deadline_s=0.0, chunk_bytes=1024)
    assert not first.complete
    assert first.bytes_written == 1024  # the deadline is checked after a chunk
    second = eeg_archive.extract(archive, dest, chunk_bytes=1024)
    assert second.complete
    assert dest.read_bytes() == payload


def test_extracting_a_complete_file_is_a_no_op(archive, tmp_path, payload):
    dest = tmp_path / f"{SESSION}.eeg"
    dest.write_bytes(payload)
    progress = eeg_archive.extract(archive, dest)
    assert progress.complete
    assert progress.bytes_written == 0


def test_a_longer_destination_is_refused_rather_than_truncated(archive, tmp_path, payload):
    dest = tmp_path / f"{SESSION}.eeg"
    dest.write_bytes(payload + b"tail")
    with pytest.raises(eeg_archive.ArchiveError, match="exceeds"):
        eeg_archive.extract(archive, dest)
    assert dest.stat().st_size == len(payload) + 4  # untouched


def test_progress_reports_fraction_and_state(archive, tmp_path, payload):
    dest = tmp_path / f"{SESSION}.eeg"
    dest.write_bytes(payload[: len(payload) // 2])
    progress = eeg_archive.extract(archive, dest, deadline_s=0.0, chunk_bytes=4096)
    assert progress.fraction == pytest.approx(0.5)
    assert "incomplete" in str(progress)


def test_a_resume_that_spends_its_budget_skipping_writes_nothing(archive, tmp_path,
                                                                 payload):
    """The deadline is checked while skipping too, and a timed-out skip is safe.

    It would be worse to write after a short skip than to make no progress: the
    bytes would land at the wrong offset. Nothing is written and the next call
    starts over from the same place on disk.
    """
    dest = tmp_path / f"{SESSION}.eeg"
    dest.write_bytes(payload[:40000])
    progress = eeg_archive.extract(archive, dest, deadline_s=0.0, chunk_bytes=1024)
    assert progress.bytes_written == 0
    assert dest.stat().st_size == 40000
    assert eeg_archive.extract(archive, dest).complete
    assert dest.read_bytes() == payload
