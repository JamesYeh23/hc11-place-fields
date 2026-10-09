"""Neuroscope XML parsing and memory-mapped .eeg access.

Fixtures mimic the real thing: 8 shanks of 16 probe channels plus EMG and
accelerometer channels, 136 channels in total, so that the two channel counts
(the full stride vs the probe-only selection pool) are genuinely different
numbers and a confusion between them cannot pass unnoticed.
"""

from __future__ import annotations

import numpy as np
import pytest

from hc11.lfp import (
    BYTES_PER_SAMPLE,
    LFP_DIR_ENV_VAR,
    LfpFile,
    LfpFormatError,
    lfp_dir,
    open_session_lfp,
    parse_xml,
    session_lfp_paths,
)

N_SHANKS = 8
PER_SHANK = 16
N_PROBE = N_SHANKS * PER_SHANK  # 128
N_AUX = 8  # 1 EMG + 3 accelerometer + 4 spare, present in the file, on no shank
N_TOTAL = N_PROBE + N_AUX  # 136 -- the stride
RATE = 1250.0


def write_xml(path, n_channels=N_TOTAL, rate=RATE, n_shanks=N_SHANKS,
              per_shank=PER_SHANK, skip=(), omit_rate=False, omit_channels=False):
    groups = []
    channel = 0
    for _ in range(n_shanks):
        lines = []
        for _ in range(per_shank):
            flag = ' skip="1"' if channel in skip else ' skip="0"'
            lines.append(f"      <channel{flag}>{channel}</channel>")
            channel += 1
        groups.append("    <group>\n      <channels>\n"
                      + "\n".join(lines)
                      + "\n      </channels>\n    </group>")
    acquisition = "".join([
        "  <acquisitionSystem>\n",
        "" if omit_channels else f"    <nChannels>{n_channels}</nChannels>\n",
        "    <samplingRate>20000</samplingRate>\n",
        "    <nBits>16</nBits>\n",
        "    <voltageRange>20</voltageRange>\n",
        "    <amplification>1000</amplification>\n",
        "  </acquisitionSystem>\n",
    ])
    potentials = "" if omit_rate else (
        f"  <fieldPotentials>\n    <lfpSamplingRate>{rate}</lfpSamplingRate>\n"
        "  </fieldPotentials>\n"
    )
    path.write_text(
        "<?xml version='1.0'?>\n<parameters>\n"
        + acquisition + potentials
        + "  <anatomicalDescription>\n    <channelGroups>\n"
        + "\n".join(groups)
        + "\n    </channelGroups>\n  </anatomicalDescription>\n</parameters>\n"
    )
    return path


def write_eeg(path, n_samples, n_channels=N_TOTAL, extra_bytes=0, seed=0):
    """A file whose value encodes (sample, channel), so mis-striding is visible."""
    rng = np.random.default_rng(seed)
    data = np.empty((n_samples, n_channels), dtype=np.int16)
    for c in range(n_channels):
        data[:, c] = (c * 100 + rng.integers(-5, 5, n_samples)).astype(np.int16)
    with open(path, "wb") as fh:
        fh.write(data.tobytes())
        if extra_bytes:
            fh.write(b"\x00" * extra_bytes)
    return path


@pytest.fixture
def recording(tmp_path):
    """A small but realistically shaped session: 60 s at 1250 Hz, 136 channels."""
    folder = tmp_path / "lfp"
    folder.mkdir()
    n_samples = int(60 * RATE)
    write_xml(folder / "Testrat_01012020.xml")
    write_eeg(folder / "Testrat_01012020.eeg", n_samples)
    return folder, n_samples


# ---------------------------------------------------------------------------
# XML
# ---------------------------------------------------------------------------


class TestParseXml:
    def test_reads_the_basics(self, recording):
        folder, _ = recording
        meta = parse_xml(folder / "Testrat_01012020.xml")
        assert meta.n_channels_total == N_TOTAL
        assert meta.sampling_rate_hz == RATE
        assert meta.wideband_rate_hz == 20000
        assert meta.n_bits == 16
        assert meta.voltage_range == 20 and meta.amplification == 1000

    def test_the_two_channel_counts_differ(self, recording):
        """The whole point of keeping them separate."""
        folder, _ = recording
        meta = parse_xml(folder / "Testrat_01012020.xml")
        assert meta.n_probe_channels == N_PROBE == 128
        assert meta.n_channels_total == N_TOTAL == 136
        assert meta.n_probe_channels != meta.n_channels_total
        assert len(meta.non_probe_channels) == N_AUX
        assert set(meta.non_probe_channels) == set(range(N_PROBE, N_TOTAL))

    def test_shank_grouping(self, recording):
        folder, _ = recording
        meta = parse_xml(folder / "Testrat_01012020.xml")
        assert meta.n_shanks == N_SHANKS
        assert sorted(meta.shanks) == list(range(1, N_SHANKS + 1))
        assert meta.shanks[1] == list(range(PER_SHANK))
        assert meta.shanks[2] == list(range(PER_SHANK, 2 * PER_SHANK))
        assert meta.probe_channels == list(range(N_PROBE))

    def test_skipped_channels(self, tmp_path):
        path = write_xml(tmp_path / "s.xml", skip=(3, 17))
        meta = parse_xml(path)
        assert meta.skipped_channels == (3, 17)
        assert 3 not in meta.good_channels(1)
        assert meta.good_channels(1) == [c for c in range(PER_SHANK) if c != 3]
        assert len(meta.good_channels(2)) == PER_SHANK - 1

    def test_missing_channel_count_raises(self, tmp_path):
        path = write_xml(tmp_path / "s.xml", omit_channels=True)
        with pytest.raises(LfpFormatError, match="nChannels"):
            parse_xml(path)

    def test_missing_sampling_rate_raises(self, tmp_path):
        """Never assume 1250 Hz: a wrong rate silently rescales every timestamp."""
        path = write_xml(tmp_path / "s.xml", omit_rate=True)
        with pytest.raises(LfpFormatError, match="lfpSamplingRate"):
            parse_xml(path)

    def test_shank_channel_outside_the_file_raises(self, tmp_path):
        path = write_xml(tmp_path / "s.xml", n_channels=64)  # 128 probe channels declared
        with pytest.raises(LfpFormatError, match="outside"):
            parse_xml(path)

    def test_session_name_defaults_to_the_stem(self, recording):
        folder, _ = recording
        assert parse_xml(folder / "Testrat_01012020.xml").session == "Testrat_01012020"
        assert parse_xml(folder / "Testrat_01012020.xml", session="X").session == "X"


# ---------------------------------------------------------------------------
# File size -- the stride check
# ---------------------------------------------------------------------------


class TestSizeCheck:
    def test_size_arithmetic_uses_the_full_channel_count(self, recording):
        folder, n_samples = recording
        lfp, meta = open_session_lfp("Testrat_01012020", folder)
        check = lfp.check_size(expected_duration_s=n_samples / RATE)
        assert check.n_channels_total == N_TOTAL, "must be the full count, not the probe count"
        assert check.bytes_per_frame == N_TOTAL * BYTES_PER_SAMPLE
        assert check.n_samples == n_samples
        assert check.exact and check.remainder_bytes == 0
        assert check.duration_difference_s == pytest.approx(0.0)
        # The probe count would have given a different, wrong answer.
        assert n_samples != lfp.n_bytes // (meta.n_probe_channels * BYTES_PER_SAMPLE)

    def test_a_wrong_channel_count_is_detected_not_accepted(self, recording):
        """The headline assertion: off-by-N must fail loudly.

        At 135 or 137 channels the byte count is no longer a whole number of
        frames. Without this check the file would still 'read' -- every frame
        shifted by one channel, so each trace becomes a slow mixture of the
        real ones, with plausible spectra and no error anywhere.
        """
        folder, n_samples = recording
        _, meta = open_session_lfp("Testrat_01012020", folder)
        for wrong in (N_TOTAL - 1, N_TOTAL + 1, meta.n_probe_channels):
            import dataclasses

            bad_meta = dataclasses.replace(meta, n_channels_total=wrong)
            bad = LfpFile(folder / "Testrat_01012020.eeg", bad_meta)
            check = bad.check_size()
            assert not check.exact, f"{wrong} channels should not divide the file evenly"
            with pytest.raises(LfpFormatError, match="whole number"):
                bad.require_consistent_size()

    def test_trailing_bytes_are_detected(self, tmp_path):
        folder = tmp_path / "lfp"
        folder.mkdir()
        write_xml(folder / "S_01012020.xml")
        write_eeg(folder / "S_01012020.eeg", 1000, extra_bytes=7)
        lfp, _ = open_session_lfp("S_01012020", folder)
        check = lfp.check_size()
        assert check.remainder_bytes == 7 and not check.exact
        with pytest.raises(LfpFormatError, match="whole number"):
            lfp.require_consistent_size()

    def test_duration_mismatch_raises(self, recording):
        folder, n_samples = recording
        lfp, _ = open_session_lfp("Testrat_01012020", folder)
        lfp.require_consistent_size(n_samples / RATE, tolerance_s=1.0)  # fine
        with pytest.raises(LfpFormatError, match="session is"):
            lfp.require_consistent_size(n_samples / RATE + 30.0, tolerance_s=1.0)

    def test_check_size_without_an_expectation(self, recording):
        folder, _ = recording
        lfp, _ = open_session_lfp("Testrat_01012020", folder)
        check = lfp.check_size()
        assert check.expected_duration_s is None
        assert check.duration_difference_s is None
        assert "bytes" in str(check)


# ---------------------------------------------------------------------------
# Reading
# ---------------------------------------------------------------------------


class TestRead:
    def test_returns_time_by_channel(self, recording):
        folder, _ = recording
        lfp, _ = open_session_lfp("Testrat_01012020", folder)
        t, data = lfp.read(1.0, 2.0)
        assert data.shape == (int(RATE), N_TOTAL)
        assert t.shape == (int(RATE),)
        assert t[0] == pytest.approx(1.0) and t[-1] < 2.0

    def test_channel_values_are_not_mixed_between_channels(self, recording):
        """Each fixture channel sits near c*100; a mis-stride would blur that."""
        folder, _ = recording
        lfp, _ = open_session_lfp("Testrat_01012020", folder)
        _, data = lfp.read(5.0, 6.0)
        for c in (0, 1, 63, N_PROBE - 1, N_PROBE, N_TOTAL - 1):
            assert abs(float(np.median(data[:, c])) - c * 100) <= 5

    def test_channel_subset(self, recording):
        folder, _ = recording
        lfp, _ = open_session_lfp("Testrat_01012020", folder)
        _, subset = lfp.read(1.0, 1.5, channels=[0, 130])
        assert subset.shape == (int(0.5 * RATE), 2)
        assert abs(float(np.median(subset[:, 1])) - 130 * 100) <= 5

    def test_single_channel_still_two_dimensional(self, recording):
        folder, _ = recording
        lfp, _ = open_session_lfp("Testrat_01012020", folder)
        _, data = lfp.read(1.0, 1.1, channels=7)
        assert data.ndim == 2 and data.shape[1] == 1

    def test_out_of_range_channel_raises(self, recording):
        folder, _ = recording
        lfp, _ = open_session_lfp("Testrat_01012020", folder)
        with pytest.raises(IndexError, match="outside"):
            lfp.read(0.0, 1.0, channels=[0, N_TOTAL])

    def test_request_beyond_the_file_is_clipped(self, recording):
        folder, n_samples = recording
        lfp, _ = open_session_lfp("Testrat_01012020", folder)
        t, data = lfp.read(59.5, 120.0)
        assert data.shape[0] == n_samples - lfp.sample_index(59.5)
        assert t[-1] < 60.0

    def test_empty_and_reversed_windows(self, recording):
        folder, _ = recording
        lfp, _ = open_session_lfp("Testrat_01012020", folder)
        _, data = lfp.read(5.0, 5.0)
        assert data.shape[0] == 0 and data.shape[1] == N_TOTAL
        with pytest.raises(ValueError, match="before"):
            lfp.read(5.0, 4.0)

    def test_context_manager_and_repr(self, recording):
        folder, _ = recording
        lfp, _ = open_session_lfp("Testrat_01012020", folder)
        with lfp as handle:
            handle.read(0.0, 0.1)
        assert "LfpFile" in repr(lfp)


class TestLargeFileStaysOutOfMemory:
    def test_reads_a_slice_of_a_file_larger_than_we_would_load(self, tmp_path):
        """A 544 MB fixture: reading a 1 s window must touch ~0.3 MB, not 544.

        Proportionally this is the real thing -- Achilles_10252013 at 136
        channels is about 11 GB -- and it fails outright if the reader ever
        materialises the file.
        """
        folder = tmp_path / "lfp"
        folder.mkdir()
        write_xml(folder / "Big_01012020.xml")
        n_samples = int(1250 * 1600)  # 1600 s
        path = folder / "Big_01012020.eeg"
        with open(path, "wb") as fh:  # write in chunks; never hold it all
            chunk = np.tile(
                (np.arange(N_TOTAL, dtype=np.int16) * 100), (12500, 1)
            )
            for _ in range(n_samples // 12500):
                fh.write(chunk.tobytes())
        size_mb = path.stat().st_size / 1e6
        assert size_mb > 500, f"fixture should be large, got {size_mb:.0f} MB"

        lfp, _ = open_session_lfp("Big_01012020", folder)
        lfp.require_consistent_size(n_samples / RATE)
        t, data = lfp.read(800.0, 801.0, channels=[5])
        assert data.shape == (1250, 1)
        assert float(np.median(data)) == pytest.approx(500)
        assert data.nbytes < 10_000, "only the window may be materialised"
        assert lfp.duration_s == pytest.approx(1600.0)


# ---------------------------------------------------------------------------
# Path resolution
# ---------------------------------------------------------------------------


class TestPaths:
    def test_env_var_wins(self, recording, monkeypatch):
        folder, _ = recording
        monkeypatch.setenv(LFP_DIR_ENV_VAR, str(folder))
        assert lfp_dir() == folder

    def test_per_session_subfolder_layout(self, tmp_path):
        root = tmp_path / "lfp"
        (root / "S_01012020").mkdir(parents=True)
        write_xml(root / "S_01012020" / "S_01012020.xml")
        write_eeg(root / "S_01012020" / "S_01012020.eeg", 100)
        eeg, xml = session_lfp_paths("S_01012020", root)
        assert eeg.parent.name == "S_01012020" and xml.suffix == ".xml"

    def test_flat_layout(self, recording):
        folder, _ = recording
        eeg, xml = session_lfp_paths("Testrat_01012020", folder)
        assert eeg.parent == folder and xml.parent == folder

    def test_missing_xml_names_both_files(self, tmp_path):
        folder = tmp_path / "lfp"
        folder.mkdir()
        write_eeg(folder / "S_01012020.eeg", 100)
        with pytest.raises(FileNotFoundError, match=r"S_01012020\.xml"):
            session_lfp_paths("S_01012020", folder)

    def test_absent_directory_lists_what_was_tried(self, tmp_path, monkeypatch):
        monkeypatch.setenv(LFP_DIR_ENV_VAR, str(tmp_path / "nope"))
        with pytest.raises(FileNotFoundError, match=LFP_DIR_ENV_VAR):
            lfp_dir()
        assert lfp_dir(required=False) == tmp_path / "nope"
