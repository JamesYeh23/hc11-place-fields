"""Neuroscope ``.xml`` metadata and memory-mapped ``.eeg`` access.

The local field potential (LFP) file is a flat binary of ``int16`` samples,
channel-interleaved, one frame per sample time::

    [ch0 t0][ch1 t0]...[chN-1 t0][ch0 t1][ch1 t1]...

Two channel counts matter here and they are **not** the same number, so they are
kept as separate, differently named quantities throughout:

``n_channels_total``
    Every channel in the file, including the electromyogram (EMG) and
    accelerometer channels the hc-11 description mentions alongside the
    hippocampal LFP. **This is the stride.** Every byte offset, the file-size
    check and the memory map all use it. Getting it wrong by even one does not
    merely misreport the duration — it mis-strides every frame, so each
    "channel" becomes a slowly drifting mixture of the real ones. That produces
    a plausible-looking trace with plausible-looking spectra and no error
    anywhere, which is why :meth:`LfpFile.check_size` is strict rather than
    advisory.

``probe_channels``
    Only the silicon-probe channels, grouped by shank, taken from the XML's
    ``anatomicalDescription`` / ``spikeDetection`` groups. **Channel selection
    for ripple detection operates on these alone.** An EMG channel has no
    pyramidal layer and no ripples, but it does have power in the 150-300 Hz
    band, so letting one into the selection pool would quietly win.

Neuroscope numbers channels from 0, and so do we.
"""

from __future__ import annotations

import os
import xml.etree.ElementTree as ET
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

#: Environment variable pointing at the directory holding the LFP recordings.
LFP_DIR_ENV_VAR = "HC11_LFP_DIR"

#: Bytes per sample in a Neuroscope ``.eeg`` file.
BYTES_PER_SAMPLE = 2
EEG_DTYPE = np.int16


class LfpFormatError(ValueError):
    """The XML or the .eeg file is not shaped the way the format requires."""


# ---------------------------------------------------------------------------
# Path resolution
# ---------------------------------------------------------------------------


def lfp_dir(required: bool = True) -> Path:
    """Directory holding the per-session LFP recordings.

    Resolution order: ``$HC11_LFP_DIR`` → ``<repo>/data/raw/lfp`` →
    ``<repo>/../hc-11_LFP``. The files are tens of gigabytes, so they usually
    live outside the repository and outside any synced folder; the environment
    variable is the intended way to point at them.
    """
    from hc11 import config

    root = config.repo_root()
    candidates = []
    env = os.environ.get(LFP_DIR_ENV_VAR)
    if env:
        candidates.append(Path(env).expanduser())
    candidates += [root / "data" / "raw" / "lfp", root.parent / "hc-11_LFP"]
    for candidate in candidates:
        if candidate.is_dir() and (
            any(candidate.glob("*.eeg")) or any(candidate.glob("*/*.eeg"))
        ):
            return candidate
    if not required:
        return candidates[0]
    listed = "\n  ".join(str(c) for c in candidates)
    raise FileNotFoundError(
        f"No .eeg files found. Looked in:\n  {listed}\n"
        f"Set {LFP_DIR_ENV_VAR} to the directory holding them."
    )


def session_lfp_paths(session: str, root: Path | None = None) -> tuple[Path, Path]:
    """``(eeg_path, xml_path)`` for one session, flat or in a per-session folder."""
    base = Path(root) if root is not None else lfp_dir()
    for folder in (base / session, base):
        eeg, xml = folder / f"{session}.eeg", folder / f"{session}.xml"
        if eeg.is_file() and xml.is_file():
            return eeg, xml
    raise FileNotFoundError(
        f"{session}: need both {session}.eeg and {session}.xml under {base} "
        f"(directly or in a {session}/ subfolder)"
    )


# ---------------------------------------------------------------------------
# XML metadata
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class LfpMetadata:
    """What the Neuroscope ``.xml`` says about a recording."""

    session: str
    n_channels_total: int  # the stride: every channel in the file
    sampling_rate_hz: float  # of the .eeg (lfpSamplingRate)
    wideband_rate_hz: float | None
    n_bits: int
    voltage_range: float | None
    amplification: float | None
    shanks: dict[int, list[int]] = field(default_factory=dict)  # shank -> channels
    skipped_channels: tuple[int, ...] = ()  # Neuroscope "skip" flag
    path: Path | None = None

    @property
    def probe_channels(self) -> list[int]:
        """Silicon-probe channels only, in shank order. Never the stride."""
        return [c for shank in sorted(self.shanks) for c in self.shanks[shank]]

    @property
    def n_probe_channels(self) -> int:
        return len(self.probe_channels)

    @property
    def non_probe_channels(self) -> list[int]:
        """Channels present in the file but on no shank: EMG, accelerometer, sync."""
        probe = set(self.probe_channels)
        return [c for c in range(self.n_channels_total) if c not in probe]

    @property
    def n_shanks(self) -> int:
        return len(self.shanks)

    def good_channels(self, shank: int) -> list[int]:
        """A shank's channels with the skipped ones removed."""
        skipped = set(self.skipped_channels)
        return [c for c in self.shanks[shank] if c not in skipped]

    def __str__(self) -> str:
        return (
            f"{self.session}: {self.n_channels_total} channels total "
            f"({self.n_probe_channels} probe on {self.n_shanks} shanks, "
            f"{len(self.non_probe_channels)} non-probe), "
            f"{self.sampling_rate_hz:g} Hz, {self.n_bits}-bit"
            + (f", {len(self.skipped_channels)} skipped" if self.skipped_channels else "")
        )


def _first_float(root: ET.Element, *paths: str) -> float | None:
    for path in paths:
        node = root.find(path)
        if node is not None and node.text and node.text.strip():
            try:
                return float(node.text)
            except ValueError:
                continue
    return None


def parse_xml(path: str | Path, session: str | None = None) -> LfpMetadata:
    """Read a Neuroscope ``.xml`` parameter file.

    Raises :class:`LfpFormatError` if the channel count or sampling rate is
    missing: without them the ``.eeg`` cannot be indexed at all, and a guess
    would silently mis-stride every sample.
    """
    path = Path(path)
    root = ET.parse(path).getroot()
    name = session if session is not None else path.stem

    n_channels = _first_float(root, "acquisitionSystem/nChannels")
    if n_channels is None:
        raise LfpFormatError(f"{path.name}: no acquisitionSystem/nChannels")
    rate = _first_float(root, "fieldPotentials/lfpSamplingRate")
    if rate is None:
        raise LfpFormatError(
            f"{path.name}: no fieldPotentials/lfpSamplingRate; refusing to assume one"
        )

    shanks: dict[int, list[int]] = {}
    for container in ("anatomicalDescription/channelGroups", "spikeDetection/channelGroups"):
        node = root.find(container)
        if node is None:
            continue
        for index, group in enumerate(node.findall("group"), start=1):
            channels = [
                int(c.text) for c in group.iter("channel") if c.text and c.text.strip()
            ]
            if channels:
                shanks.setdefault(index, channels)
        if shanks:
            break  # anatomical grouping wins when both are present

    skipped = tuple(
        int(c.text)
        for c in root.iter("channel")
        if c.get("skip") == "1" and c.text and c.text.strip()
    )

    meta = LfpMetadata(
        session=name,
        n_channels_total=int(n_channels),
        sampling_rate_hz=rate,
        wideband_rate_hz=_first_float(root, "acquisitionSystem/samplingRate"),
        n_bits=int(_first_float(root, "acquisitionSystem/nBits") or 16),
        voltage_range=_first_float(root, "acquisitionSystem/voltageRange"),
        amplification=_first_float(root, "acquisitionSystem/amplification"),
        shanks=shanks,
        skipped_channels=skipped,
        path=path,
    )
    out_of_range = [
        c for c in meta.probe_channels if not 0 <= c < meta.n_channels_total
    ]
    if out_of_range:
        raise LfpFormatError(
            f"{path.name}: shank channels outside 0..{meta.n_channels_total - 1}: "
            f"{out_of_range[:8]}"
        )
    return meta


# ---------------------------------------------------------------------------
# The .eeg file
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class SizeCheck:
    """Comparison of the file's size against the metadata and the session clock."""

    path: Path
    n_bytes: int
    n_channels_total: int
    bytes_per_frame: int
    n_samples: int
    remainder_bytes: int
    duration_s: float
    expected_duration_s: float | None

    @property
    def exact(self) -> bool:
        """True when the file is a whole number of frames at this channel count."""
        return self.remainder_bytes == 0

    @property
    def duration_difference_s(self) -> float | None:
        if self.expected_duration_s is None:
            return None
        return self.duration_s - self.expected_duration_s

    def __str__(self) -> str:
        text = (
            f"{self.path.name}: {self.n_bytes:,} bytes / "
            f"({self.n_channels_total} ch x {BYTES_PER_SAMPLE} B) = "
            f"{self.n_samples:,} samples = {self.duration_s:.3f} s"
        )
        if not self.exact:
            text += f"  [{self.remainder_bytes} trailing bytes -- NOT a whole frame]"
        if self.expected_duration_s is not None:
            text += (
                f"; sessInfo says {self.expected_duration_s:.3f} s "
                f"(difference {self.duration_difference_s:+.3f} s)"
            )
        return text


class LfpFile:
    """Memory-mapped read access to one ``.eeg`` recording.

    The file is never read as a whole: :meth:`read` maps it and returns a copy
    of the requested window only, shaped ``(time, channel)``.
    """

    def __init__(self, path: str | Path, meta: LfpMetadata):
        self.path = Path(path)
        self.meta = meta
        if meta.n_channels_total <= 0:
            raise LfpFormatError(f"{self.path.name}: non-positive channel count")
        self._memmap: np.memmap | None = None

    # --- geometry -----------------------------------------------------------

    @property
    def bytes_per_frame(self) -> int:
        """One sample across **all** channels — the stride."""
        return self.meta.n_channels_total * BYTES_PER_SAMPLE

    @property
    def n_bytes(self) -> int:
        return self.path.stat().st_size

    @property
    def n_samples(self) -> int:
        return self.n_bytes // self.bytes_per_frame

    @property
    def duration_s(self) -> float:
        return self.n_samples / self.meta.sampling_rate_hz

    def check_size(self, expected_duration_s: float | None = None) -> SizeCheck:
        """Compare the file size against the **full** channel count from the XML.

        A file whose size is not a whole number of frames is the signature of a
        wrong channel count, and :attr:`SizeCheck.exact` reports it. Callers
        that cannot tolerate a mis-strided read should use
        :meth:`require_consistent_size`.
        """
        n_bytes = self.n_bytes
        return SizeCheck(
            path=self.path, n_bytes=n_bytes,
            n_channels_total=self.meta.n_channels_total,
            bytes_per_frame=self.bytes_per_frame,
            n_samples=n_bytes // self.bytes_per_frame,
            remainder_bytes=n_bytes % self.bytes_per_frame,
            duration_s=(n_bytes // self.bytes_per_frame) / self.meta.sampling_rate_hz,
            expected_duration_s=expected_duration_s,
        )

    def require_consistent_size(
        self, expected_duration_s: float | None = None, tolerance_s: float = 1.0
    ) -> SizeCheck:
        """:meth:`check_size`, but raise unless the geometry holds up.

        Raises if the file is not a whole number of frames, or if its duration
        disagrees with ``expected_duration_s`` by more than ``tolerance_s``.
        """
        check = self.check_size(expected_duration_s)
        if not check.exact:
            raise LfpFormatError(
                f"{self.path.name}: {check.n_bytes:,} bytes is not a whole number of "
                f"{check.bytes_per_frame}-byte frames ({check.remainder_bytes} left over). "
                f"The channel count of {check.n_channels_total} is probably wrong; a wrong "
                "count mis-strides every sample and yields plausible-looking garbage."
            )
        if (
            expected_duration_s is not None
            and abs(check.duration_difference_s) > tolerance_s
        ):
            raise LfpFormatError(
                f"{self.path.name}: file is {check.duration_s:.3f} s at "
                f"{check.n_channels_total} channels but the session is "
                f"{expected_duration_s:.3f} s "
                f"({check.duration_difference_s:+.3f} s, tolerance {tolerance_s} s)."
            )
        return check

    # --- reading ------------------------------------------------------------

    @property
    def data(self) -> np.memmap:
        """The whole file as a lazy ``(time, channel)`` memory map."""
        if self._memmap is None:
            self._memmap = np.memmap(
                self.path, dtype=EEG_DTYPE, mode="r",
                shape=(self.n_samples, self.meta.n_channels_total),
            )
        return self._memmap

    def sample_index(self, t: float) -> int:
        return int(round(t * self.meta.sampling_rate_hz))

    def read(
        self,
        start_s: float,
        end_s: float,
        channels: int | list[int] | np.ndarray | None = None,
    ) -> tuple[np.ndarray, np.ndarray]:
        """Return ``(t, data)`` for one interval; ``data`` is ``(time, channel)``.

        Only the requested window is materialised. ``channels`` indexes the file
        — i.e. the same numbering the XML uses — and defaults to all of them.
        """
        if end_s < start_s:
            raise ValueError(f"end_s ({end_s}) is before start_s ({start_s})")
        lo = max(self.sample_index(start_s), 0)
        hi = min(self.sample_index(end_s), self.n_samples)
        if channels is None:
            picked = slice(None)
            n_out = self.meta.n_channels_total
        else:
            picked = np.atleast_1d(np.asarray(channels, dtype=np.int64))
            bad = picked[(picked < 0) | (picked >= self.meta.n_channels_total)]
            if bad.size:
                raise IndexError(
                    f"channel(s) {bad.tolist()} outside "
                    f"0..{self.meta.n_channels_total - 1}"
                )
            n_out = picked.size
        if hi <= lo:
            return np.empty(0), np.empty((0, n_out), dtype=EEG_DTYPE)
        block = np.array(self.data[lo:hi, picked])
        if block.ndim == 1:
            block = block[:, np.newaxis]
        t = (lo + np.arange(hi - lo)) / self.meta.sampling_rate_hz
        return t, block

    def close(self) -> None:
        self._memmap = None

    def __enter__(self) -> LfpFile:
        return self

    def __exit__(self, *exc) -> None:
        self.close()

    def __repr__(self) -> str:
        return (
            f"LfpFile({self.path.name}, {self.n_samples:,} samples x "
            f"{self.meta.n_channels_total} ch, {self.duration_s / 3600:.2f} h)"
        )


def open_session_lfp(session: str, root: Path | None = None) -> tuple[LfpFile, LfpMetadata]:
    """Open one session's ``.eeg`` with the metadata from its ``.xml``."""
    eeg_path, xml_path = session_lfp_paths(session, root)
    meta = parse_xml(xml_path, session=session)
    return LfpFile(eeg_path, meta), meta
