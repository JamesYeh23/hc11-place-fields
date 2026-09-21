"""Each validate() check fires on a deliberately broken synthetic session."""

from __future__ import annotations

import numpy as np
import pytest

from conftest import SYNTH_NAME, write_sessinfo
from hc11.io import load_session
from hc11.validate import validate


def _issues(tmp_path, arrays):
    path = write_sessinfo(tmp_path / f"{SYNTH_NAME}_sessInfo.mat", arrays)
    return {(i.severity, i.code) for i in validate(load_session(path))}


def _perturb(arrays, **changes):
    out = {k: (v.copy() if isinstance(v, np.ndarray) else v) for k, v in arrays.items()}
    out.update(changes)
    return out


def test_clean_session_has_only_info(tmp_path, synthetic_arrays):
    assert {sev for sev, _ in _issues(tmp_path, synthetic_arrays)} == {"info"}


@pytest.mark.parametrize(
    ("name", "change", "expected"),
    [
        (
            "spikes past end",
            lambda a: {"SpikeTimes": np.append(a["SpikeTimes"][:-1], [[350.0]], axis=0)},
            ("warning", "spikes_after_session_end"),
        ),
        (
            "unsorted spikes",
            lambda a: {"SpikeTimes": a["SpikeTimes"][::-1].copy()},
            ("error", "spike_times_unsorted"),
        ),
        (
            "pyr/int overlap",
            lambda a: {"IntIDs": np.array([[104.0], [102.0]])},
            ("error", "pyr_int_overlap"),
        ),
        (
            "unclassified cluster",
            lambda a: {"IntIDs": np.array([[104.0]])},  # 210 fires but is unlabelled
            ("warning", "unclassified_clusters"),
        ),
        (
            "epochs overlap",
            lambda a: {"PREEpoch": np.array([[0.0, 110.0]])},
            ("error", "epochs_overlap"),
        ),
        (
            "position outside maze",
            lambda a: {"MazeEpoch": np.array([[105.0, 160.0]]),
                       "PREEpoch": np.array([[0.0, 105.0]])},
            ("error", "position_outside_maze"),
        ),
        (
            "zero-length state",
            lambda a: {"REM": np.array([[82.0, 95.0], [250.0, 250.0]])},
            ("warning", "state_zero_length"),
        ),
        (
            "negative-length state",
            lambda a: {"REM": np.array([[95.0, 82.0]])},
            ("error", "state_negative_length"),
        ),
        (
            "cross-state overlap",
            lambda a: {"Drowsy": np.array([[15.0, 30.0]])},  # Wake is [0, 20]
            ("warning", "states_overlap"),
        ),
        (
            "state past end",
            lambda a: {"NREM": np.array([[30.0, 80.0], [170.0, 320.0]])},
            ("warning", "state_outside_session"),
        ),
        (
            "1-D beyond track",
            lambda a: {"OneDLocation": np.where(np.isnan(a["OneDLocation"]),
                                                np.nan, a["OneDLocation"] + 0.5)},
            ("warning", "position_1d_out_of_track"),
        ),
    ],
)
def test_check_fires(tmp_path, synthetic_arrays, name, change, expected):
    arrays = _perturb(synthetic_arrays, **change(synthetic_arrays))
    assert expected in _issues(tmp_path, arrays), name


def test_zero_length_row_does_not_mask_self_overlap(tmp_path, synthetic_arrays):
    rem = np.array([[82.0, 95.0], [90.0, 100.0], [250.0, 250.0]])
    issues = _issues(tmp_path, _perturb(synthetic_arrays, REM=rem))
    assert ("warning", "state_self_overlap") in issues
