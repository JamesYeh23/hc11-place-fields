"""Tests for hc11.io: MATLAB v7.3 helpers, the Session object, and real data."""

from __future__ import annotations

import dataclasses

import h5py
import numpy as np
import pytest

from conftest import (
    SYNTH_NAME,
    make_struct,
    requires_data,
    write_matlab,
    write_matlab_empty,
    write_matlab_string,
    write_sessinfo,
)
from hc11 import config
from hc11.io import (
    STATE_NAMES,
    as_int_ids,
    as_intervals,
    as_vector,
    decode_matlab_string,
    is_matlab_empty,
    list_sessions,
    load_session,
    parse_animal,
    read_matlab,
    shank_of,
)
from hc11.validate import errors, validate

# ===========================================================================
# Unit tests: generic MATLAB v7.3 reading
# ===========================================================================


@pytest.fixture
def h5(tmp_path):
    with h5py.File(tmp_path / "t.mat", "w") as f:
        yield f


class TestOrientation:
    def test_n_by_2_matrix_round_trips(self, h5):
        m = np.array([[1.0, 2.0], [3.0, 4.0], [5.0, 6.0]])  # 3 intervals
        ds = write_matlab(h5, "m", m)
        assert ds.shape == (2, 3), "fixture should store it transposed, like MATLAB"
        np.testing.assert_array_equal(read_matlab(ds), m)

    def test_square_2x2_is_not_ambiguous(self, h5):
        # A shape-based heuristic ("transpose if it has 2 rows") gets this wrong.
        m = np.array([[10.0, 20.0], [30.0, 40.0]])  # rows are [start, end]
        iv = as_intervals(read_matlab(write_matlab(h5, "m", m)))
        np.testing.assert_array_equal(iv, m)

    @pytest.mark.parametrize("shape", [(1, 5), (5, 1)])
    def test_row_and_column_vectors_become_1d(self, h5, shape):
        v = np.arange(5.0).reshape(shape)
        out = as_vector(read_matlab(write_matlab(h5, "v", v)))
        assert out.shape == (5,)
        np.testing.assert_array_equal(out, np.arange(5.0))

    def test_as_vector_rejects_matrix(self):
        with pytest.raises(ValueError, match="expected a vector"):
            as_vector(np.zeros((3, 2)))

    def test_single_interval_from_1x2(self, h5):
        iv = as_intervals(read_matlab(write_matlab(h5, "e", np.array([[0.0, 18079.5]]))))
        assert iv.shape == (1, 2)
        np.testing.assert_array_equal(iv[0], [0.0, 18079.5])

    def test_as_intervals_rejects_wrong_width(self):
        with pytest.raises(ValueError, match="n x 2"):
            as_intervals(np.zeros((4, 3)))


class TestStrings:
    def test_decode_row(self, h5):
        ds = write_matlab_string(h5, "s", "Circular Maze")
        assert ds.dtype == np.uint16 and ds.shape == (13, 1)
        assert read_matlab(ds) == "Circular Maze"

    def test_decode_non_ascii(self, h5):
        assert read_matlab(write_matlab_string(h5, "s", "Buzsáki")) == "Buzsáki"

    def test_decode_char_matrix(self):
        rows = np.array([[ord(c) for c in "ab "], [ord(c) for c in "cde"]], dtype=np.uint16)
        assert decode_matlab_string(rows) == ["ab", "cde"]

    def test_decode_rejects_floats(self):
        with pytest.raises(TypeError):
            decode_matlab_string(np.array([1.5, 2.5]))


class TestEmpty:
    def test_placeholder_is_detected(self, h5):
        ds = write_matlab_empty(h5, "e", dims=(0, 2))
        assert is_matlab_empty(ds)
        assert not is_matlab_empty(write_matlab(h5, "x", np.array([[0.0, 2.0]])))

    def test_placeholder_values_are_never_returned(self, h5):
        # The placeholder stores [0, 2] -- which looks exactly like an interval.
        out = read_matlab(write_matlab_empty(h5, "e", dims=(0, 2)))
        assert out.size == 0
        iv = as_intervals(out)
        assert iv.shape == (0, 2)

    def test_empty_char_is_empty_string(self, h5):
        assert read_matlab(write_matlab_empty(h5, "e", cls="char")) == ""


class TestReferences:
    def test_cell_array_of_refs_is_dereferenced(self, h5):
        a = write_matlab(h5, "a", np.array([[1.0, 2.0]]))
        b = write_matlab_string(h5, "b", "hello")
        refs = h5.create_dataset("cell", data=np.array([[a.ref], [b.ref]]), dtype=h5py.ref_dtype)
        refs.attrs["MATLAB_class"] = np.bytes_("cell")
        out = read_matlab(refs)
        assert isinstance(out, list) and len(out) == 2
        np.testing.assert_array_equal(out[0], [[1.0, 2.0]])
        assert out[1] == "hello"

    def test_struct_group_becomes_dict(self, h5):
        g = make_struct(h5, "s")
        write_matlab(g, "x", np.array([[1.0]]))
        assert set(read_matlab(g)) == {"x"}


class TestIds:
    def test_whole_floats_cast_to_int64(self):
        out = as_int_ids(np.array([[102.0, 1319.0]]))
        assert out.dtype == np.int64
        np.testing.assert_array_equal(out, [102, 1319])

    @pytest.mark.parametrize("bad", [[102.5], [np.nan], [np.inf]])
    def test_non_integer_ids_raise(self, bad):
        with pytest.raises(ValueError):
            as_int_ids(np.array(bad))

    def test_shank_is_hundreds_digit(self):
        assert shank_of(1319) == 13
        np.testing.assert_array_equal(shank_of(np.array([102, 205, 1003])), [1, 2, 10])

    def test_parse_animal(self):
        assert parse_animal("Achilles_10252013") == "Achilles"
        with pytest.raises(ValueError):
            parse_animal("not-a-session")


# ===========================================================================
# Unit tests: load_session on a synthetic file
# ===========================================================================


class TestSyntheticSession:
    def test_loads_with_correct_shapes(self, synthetic_file, synthetic_arrays):
        s = load_session(synthetic_file)
        n_spk = synthetic_arrays["SpikeTimes"].size
        n_pos = synthetic_arrays["TimeStamps"].size
        assert s.name == SYNTH_NAME and s.animal == "Testrat"
        assert s.spike_times.shape == (n_spk,) and s.spike_ids.shape == (n_spk,)
        assert s.spike_ids.dtype == np.int64
        assert s.position_t.shape == (n_pos,)
        assert s.position_xy.shape == (n_pos, 2)
        assert s.position_1d.shape == (n_pos,)
        assert s.maze_type == "1.6m Linear Maze"
        assert s.maze_kind == "linear" and s.track_length_m == 1.6
        assert set(s.epochs) == {"PRE", "MAZE", "POST"}
        assert all(v.shape == (2,) for v in s.epochs.values())
        assert set(s.states) == set(STATE_NAMES)
        assert all(v.ndim == 2 and v.shape[1] == 2 for v in s.states.values())
        np.testing.assert_array_equal(s.states["NREM"], synthetic_arrays["NREM"])
        assert s.sess_duration == 300.0

    def test_values_preserved_exactly(self, synthetic_file, synthetic_arrays):
        s = load_session(synthetic_file)
        np.testing.assert_array_equal(s.spike_times, synthetic_arrays["SpikeTimes"].ravel())
        np.testing.assert_array_equal(s.position_xy, synthetic_arrays["TwoDLocation"])
        np.testing.assert_array_equal(s.position_t, synthetic_arrays["TimeStamps"].ravel())

    def test_cm_helpers(self, synthetic_file):
        s = load_session(synthetic_file)
        np.testing.assert_allclose(s.position_xy_cm, s.position_xy * 100)
        np.testing.assert_allclose(s.position_1d_cm, s.position_1d * 100)

    def test_is_immutable(self, synthetic_file):
        s = load_session(synthetic_file)
        with pytest.raises(dataclasses.FrozenInstanceError):
            s.name = "x"
        with pytest.raises(ValueError):
            s.spike_times[0] = -1.0
        with pytest.raises(TypeError):
            s.epochs["PRE"] = np.zeros(2)

    def test_units_and_shanks(self, synthetic_file):
        s = load_session(synthetic_file)
        np.testing.assert_array_equal(s.units("pyr"), [102, 103, 205])
        np.testing.assert_array_equal(s.units("int"), [104, 210])
        np.testing.assert_array_equal(s.units("all"), [102, 103, 104, 205, 210])
        assert s.shank_of_cluster[205] == 2
        with pytest.raises(ValueError):
            s.units("granule")

    def test_spikes_for_matches_boolean_mask(self, synthetic_file):
        s = load_session(synthetic_file)
        for cid in s.units("all"):
            expected = s.spike_times[s.spike_ids == cid]
            np.testing.assert_array_equal(s.spikes_for(cid), expected)
        assert s.spikes_for(999).size == 0

    def test_spikes_in_interval_is_half_open(self, synthetic_file):
        s = load_session(synthetic_file)
        start, end = s.epoch("MAZE")
        t, ids = s.spikes_in((start, end))
        mask = (s.spike_times >= start) & (s.spike_times < end)
        np.testing.assert_array_equal(t, s.spike_times[mask])
        np.testing.assert_array_equal(ids, s.spike_ids[mask])
        t1, ids1 = s.spikes_in((start, end), cluster_id=102)
        np.testing.assert_array_equal(t1, s.spike_times[mask & (s.spike_ids == 102)])
        assert np.all(ids1 == 102)

    def test_empty_state_field_loads_as_0x2(self, tmp_path, synthetic_arrays):
        path = write_sessinfo(
            tmp_path / f"{SYNTH_NAME}_sessInfo.mat", synthetic_arrays, empty_fields=("Drowsy",)
        )
        s = load_session(path)
        assert s.states["Drowsy"].shape == (0, 2)

    def test_mismatched_spike_lengths_raise(self, tmp_path, synthetic_arrays):
        arrays = dict(synthetic_arrays, SpikeIDs=synthetic_arrays["SpikeIDs"][:-1])
        path = write_sessinfo(tmp_path / f"{SYNTH_NAME}_sessInfo.mat", arrays)
        with pytest.raises(ValueError, match="spike times but"):
            load_session(path)

    def test_list_sessions_finds_nested_files(self, tmp_path, synthetic_arrays):
        # Unpacking the full per-session archives gives one folder per session.
        (tmp_path / SYNTH_NAME).mkdir()
        write_sessinfo(tmp_path / SYNTH_NAME / f"{SYNTH_NAME}_sessInfo.mat", synthetic_arrays)
        assert list_sessions(tmp_path) == [SYNTH_NAME]

    def test_synthetic_session_validates_clean(self, synthetic_file):
        issues = validate(load_session(synthetic_file))
        assert [i for i in issues if i.severity != "info"] == []


# ===========================================================================
# Integration tests: every real session
# ===========================================================================

#: Transcribed from Sessions_Recordings_Summary.pdf (durations rounded to 1 s).
RECORDING_SUMMARY = {
    #  session             maze                  PRE     MAZE   POST   pyr  int
    "Buddy_06272013":    ("1.6m Linear Maze", 10718, 2328, 8009, 48, 20),
    "Gatsby_08022013":   ("1.6m Linear Maze", 17079, 2335, 10999, 66, 14),
    "Gatsby_08282013":   ("Circular Maze",    16761, 2693, 15894, 41, 10),
    "Achilles_10252013": ("1.6m Linear Maze", 18080, 2068, 14714, 120, 17),
    "Achilles_11012013": ("Circular Maze",    19330, 2674, 14625, 92, 12),
    "Cicero_09012014":   ("1.6m Linear Maze", 15532, 5561, 13664, 55, 18),
    "Cicero_09102014":   ("Circular Maze",    17246, 3030, 16223, 81, 24),
    "Cicero_09172014":   ("2m Linear Maze",   16001, 3078, 14482, 59, 13),
}

#: Anomalies found during step 2 inspection (docs/data_schema.md). Pinned so
#: that any *new* warning -- or a known one disappearing -- fails loudly.
KNOWN_WARNINGS = {
    "Achilles_10252013": set(),
    "Achilles_11012013": {"state_zero_length"},
    "Buddy_06272013": set(),
    "Cicero_09012014": set(),
    "Cicero_09102014": set(),
    "Cicero_09172014": set(),
    "Gatsby_08022013": {"spikes_after_session_end", "state_outside_session"},
    "Gatsby_08282013": set(),
}


@requires_data
def test_all_sessions_listed():
    assert list_sessions() == list(config.SESSIONS)


@requires_data
@pytest.mark.parametrize("name", config.SESSIONS)
class TestRealSessions:
    def test_no_validation_errors(self, real_sessions, name):
        issues = validate(real_sessions(name))
        assert errors(issues) == [], "\n".join(map(str, errors(issues)))

    def test_known_warnings_only(self, real_sessions, name):
        codes = {i.code for i in validate(real_sessions(name)) if i.severity == "warning"}
        assert codes == KNOWN_WARNINGS[name]

    def test_matches_recording_summary(self, real_sessions, name):
        s = real_sessions(name)
        maze, pre, mz, post, n_pyr, n_int = RECORDING_SUMMARY[name]
        assert s.maze_type == maze
        assert abs(s.epoch_duration("PRE") - pre) <= 1
        assert abs(s.epoch_duration("MAZE") - mz) <= 1
        assert abs(s.epoch_duration("POST") - post) <= 1
        assert len(s.pyr_ids) == n_pyr
        assert len(s.int_ids) == n_int

    def test_orientation_is_logical(self, real_sessions, name):
        s = real_sessions(name)
        assert s.position_xy.ndim == 2 and s.position_xy.shape[1] == 2
        assert s.position_xy.shape[0] > 50_000
        assert s.position_t.ndim == 1 and s.position_1d.ndim == 1
        assert s.spike_times.ndim == 1 and s.spike_times.size > 1_000_000
        for iv in s.states.values():
            assert iv.shape[1] == 2 and iv.shape[0] >= 1

    def test_spikes_for_is_consistent(self, real_sessions, name):
        s = real_sessions(name)
        total = sum(s.spikes_for(c).size for c in s.cluster_ids)
        assert total == s.spike_times.size


@requires_data
def test_documented_example_matches_file(real_sessions):
    """The data description's worked example is Gatsby_08282013."""
    s = real_sessions("Gatsby_08282013")
    assert s.spike_times.size == 5_775_679
    assert len(s.pyr_ids) == 41 and len(s.int_ids) == 10
    assert s.position_t.size == 105_180
    assert s.maze_type == "Circular Maze"


# ===========================================================================
# Independent cross-checks of the loader
# ===========================================================================


@requires_data
@pytest.mark.parametrize("name", config.SESSIONS)
def test_agrees_with_pymatreader(real_sessions, name):
    """A second, independently written v7.3 reader must give identical values."""
    pymatreader = pytest.importorskip("pymatreader")
    m = pymatreader.read_mat(str(config.session_path(name)))["sessInfo"]
    s = real_sessions(name)
    np.testing.assert_array_equal(m["Spikes"]["SpikeTimes"].ravel(), s.spike_times)
    np.testing.assert_array_equal(m["Spikes"]["SpikeIDs"].ravel(), s.spike_ids)
    np.testing.assert_array_equal(np.sort(m["Spikes"]["PyrIDs"].ravel()), s.pyr_ids)
    np.testing.assert_array_equal(np.sort(m["Spikes"]["IntIDs"].ravel()), s.int_ids)
    np.testing.assert_array_equal(m["Position"]["TimeStamps"].ravel(), s.position_t)
    np.testing.assert_array_equal(m["Position"]["TwoDLocation"], s.position_xy)
    np.testing.assert_array_equal(m["Position"]["OneDLocation"].ravel(), s.position_1d)
    assert m["Position"]["MazeType"] == s.maze_type
    for ours, theirs in (("PRE", "PREEpoch"), ("MAZE", "MazeEpoch"), ("POST", "POSTEpoch")):
        np.testing.assert_array_equal(m["Epochs"][theirs].ravel(), s.epochs[ours])
    for state in STATE_NAMES:
        np.testing.assert_array_equal(np.atleast_2d(m["Epochs"][state]), s.states[state])
    assert float(m["Epochs"]["sessDuration"]) == s.sess_duration


@requires_data
@pytest.mark.parametrize("name", config.SESSIONS)
def test_linear_1d_is_projection_of_2d(real_sessions, name):
    """On a linear track, OneDLocation must be a linear function of 2-D position.

    This is a physical check that x and y are paired per sample and that the
    position arrays were not transposed or reordered: a scrambled layout would
    destroy the correlation.
    """
    s = real_sessions(name)
    if s.maze_kind != "linear":
        pytest.skip("circular maze")
    ok = ~np.isnan(s.position_1d) & ~np.isnan(s.position_xy).any(axis=1)
    centred = s.position_xy[ok] - s.position_xy[ok].mean(axis=0)
    pc1 = np.linalg.svd(centred, full_matrices=False)[2][0]
    r = np.corrcoef(centred @ pc1, s.position_1d[ok])[0, 1]
    assert abs(r) > 0.999
