"""Resumable extraction of a ``.eeg`` recording from a CRCNS ``.tar.gz``.

The hc-11 LFP downloads are per-session gzipped tarballs holding a single
member, ``<session>/<session>.eeg``, of 11-12 GB. Two things about them shape
this module:

**The uncompressed size is in the tar header, at the front of the stream.**
:func:`member_info` reads it without decompressing the body, so the file-size
consistency check (:func:`hc11.lfp.size_check`) can run against the archive
*before* spending minutes on extraction. If the declared size is not a whole
number of frames at the channel count in the XML, nothing downstream is worth
doing, and we find that out in under a second.

**Extraction has to be interruptible and resumable.** Writing 12 GB across a
network-backed mount takes minutes, and the surrounding automation cannot hold
a single call open that long. :func:`extract` therefore appends to the
destination and, given a ``deadline_s``, stops cleanly at a chunk boundary and
reports how far it got; the next call skips the bytes already on disk and
carries on. Resuming is safe because a single-member tar is a pure byte stream:
the first *n* bytes of the member are always the same *n* bytes.

The destination is never truncated and never overwritten in place. A
destination already *larger* than the member is treated as an error rather than
silently rewound, because that means the file on disk is not a prefix of this
archive and continuing would splice two different recordings together.
"""

from __future__ import annotations

import tarfile
import time
from dataclasses import dataclass
from pathlib import Path

#: Read/write unit. 8 MiB keeps the write syscall count low on a slow mount
#: without making the deadline check coarse (one chunk is well under a second).
CHUNK_BYTES = 8 << 20


class ArchiveError(ValueError):
    """The archive is not the single-member ``.eeg`` tarball we expect."""


@dataclass(frozen=True)
class ArchiveMember:
    """The one file inside a session's LFP archive."""

    archive: Path
    name: str  # path as stored in the tar, e.g. "Achilles_10252013/...eeg"
    size: int  # uncompressed bytes, from the tar header

    @property
    def filename(self) -> str:
        return self.name.rsplit("/", 1)[-1]

    def __str__(self) -> str:
        return (
            f"{self.archive.name}: {self.name}, {self.size:,} bytes "
            f"({self.size / 1e9:.2f} GB uncompressed)"
        )


def member_info(archive: str | Path, suffix: str = ".eeg") -> ArchiveMember:
    """The ``suffix`` member's name and uncompressed size, read from the header.

    Streams the archive (``r|gz``) and stops at the first matching member, so
    only the leading few kilobytes are decompressed.
    """
    archive = Path(archive)
    with tarfile.open(archive, mode="r|gz") as tar:
        for member in tar:
            if member.isfile() and member.name.endswith(suffix):
                return ArchiveMember(archive=archive, name=member.name, size=member.size)
    raise ArchiveError(f"{archive.name}: no {suffix} member found")


@dataclass(frozen=True)
class ExtractProgress:
    """How far one :func:`extract` call got."""

    member: ArchiveMember
    dest: Path
    bytes_at_start: int
    bytes_written: int
    elapsed_s: float

    @property
    def bytes_on_disk(self) -> int:
        return self.bytes_at_start + self.bytes_written

    @property
    def complete(self) -> bool:
        return self.bytes_on_disk >= self.member.size

    @property
    def fraction(self) -> float:
        return self.bytes_on_disk / self.member.size if self.member.size else 1.0

    @property
    def rate_mb_s(self) -> float:
        return self.bytes_written / 1e6 / self.elapsed_s if self.elapsed_s > 0 else 0.0

    def __str__(self) -> str:
        state = "complete" if self.complete else "incomplete"
        return (
            f"{self.dest.name}: {self.bytes_on_disk:,} / {self.member.size:,} bytes "
            f"({100 * self.fraction:.1f}%, {state}); wrote "
            f"{self.bytes_written / 1e9:.2f} GB in {self.elapsed_s:.0f} s "
            f"({self.rate_mb_s:.0f} MB/s)"
        )


def extract(
    archive: str | Path,
    dest: str | Path,
    deadline_s: float | None = None,
    suffix: str = ".eeg",
    chunk_bytes: int = CHUNK_BYTES,
    member: ArchiveMember | None = None,
) -> ExtractProgress:
    """Append the archive's ``suffix`` member to ``dest``, resuming if partial.

    Bytes already in ``dest`` are skipped in the decompressed stream rather than
    re-written. With ``deadline_s`` set, the call returns at the first chunk
    boundary past that many seconds; :attr:`ExtractProgress.complete` says
    whether to call again. Skipping counts against the deadline too, so a resume
    that spends its whole budget seeking makes no progress but also loses
    nothing -- it never writes at the wrong offset.
    """
    archive, dest = Path(archive), Path(dest)
    info = member if member is not None else member_info(archive, suffix)
    dest.parent.mkdir(parents=True, exist_ok=True)
    already = dest.stat().st_size if dest.exists() else 0
    started = time.monotonic()

    if already > info.size:
        raise ArchiveError(
            f"{dest.name}: {already:,} bytes on disk exceeds the {info.size:,} bytes "
            f"in {archive.name}. It is not a prefix of this archive, so resuming "
            "would splice two recordings together. Move or delete it first."
        )
    if already == info.size:
        return ExtractProgress(info, dest, already, 0, time.monotonic() - started)

    with tarfile.open(archive, mode="r|gz") as tar:
        stream = None
        for candidate in tar:
            if candidate.isfile() and candidate.name == info.name:
                stream = tar.extractfile(candidate)
                break
        if stream is None:
            raise ArchiveError(f"{archive.name}: member {info.name} not found on re-read")

        remaining_skip = already
        while remaining_skip > 0:
            block = stream.read(min(chunk_bytes, remaining_skip))
            if not block:
                raise ArchiveError(
                    f"{archive.name}: stream ended {remaining_skip:,} bytes into the "
                    f"{already:,} bytes already on disk; the archive and {dest.name} "
                    "do not agree."
                )
            remaining_skip -= len(block)
            if deadline_s is not None and time.monotonic() - started > deadline_s:
                return ExtractProgress(info, dest, already, 0, time.monotonic() - started)

        written = 0
        with open(dest, "ab") as out:
            while True:
                block = stream.read(chunk_bytes)
                if not block:
                    break
                out.write(block)
                written += len(block)
                if deadline_s is not None and time.monotonic() - started > deadline_s:
                    break

    return ExtractProgress(info, dest, already, written, time.monotonic() - started)
