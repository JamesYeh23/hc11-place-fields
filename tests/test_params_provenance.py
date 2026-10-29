"""Every externally attributed parameter, frozen against its citation.

Twice in this project the config credited a parameter to a paper that does not
state it, while `docs/decisions.md` recorded it correctly (the multi-unit
boundary threshold, D10.1; the ripple envelope threshold, D10.3). Once it ran
the other way, and the decision log was wrong too (`min_peak_rate_hz`, which is
Chen et al.'s and which D8.1's table called ours).

So "the decision log is authoritative" is not a check, it is a hope. This is the
check: the set of `[CHEN2016]`, `[GB2016]` and `[DATA]` tags in
`config/params.yaml` is frozen here, and **a parameter cannot acquire, lose or
change an attribution without someone editing this file**. Editing this file is
the moment to open the paper and read the sentence, which is the only thing that
would have caught any of the three.

It deliberately does not check the *values*: those move for good reasons and
have their own tests. Only who is said to have chosen them.
"""

from __future__ import annotations

import re
from pathlib import Path

#: Parameter path -> the source it is attributed to, as of D10.3.
#: Changing a line here means you have read the paper again. That is the point.
ATTRIBUTED = {
    "behavior.require_linearized": "DATA",
    "behavior.speed_threshold_cm_s": "CHEN2016",
    "dataset.position_rate_hz": "DATA",
    "dataset.spike_time_units": "DATA",
    "decoding.prior": "CHEN2016",
    "decoding.time_bin_s": "CHEN2016",
    "decoding.time_bins_s": "CHEN2016",
    "events.mua.max_duration_ms": "CHEN2016",
    "events.mua.threshold_high_sd": "CHEN2016",
    "events.mua.unit_kind": "GB2016",
    "events.ripples.band": "CHEN2016",
    "events.ripples.max_duration_ms": "CHEN2016",
    "events.states.maze_immobility_speed_cm_s": "CHEN2016",
    "events.states.sleep_state": "DATA",
    "lfp.channel_count_override.Gatsby_08022013": "DATA",
    "lfp.channel_selection.bands.ripple": "CHEN2016",
    "place_fields.bin_size_cm": "CHEN2016",
    "place_fields.criteria.min_peak_rate_hz": "CHEN2016",
    "place_fields.expected_total_place_cells": "GB2016",
    "place_fields.separate_directions": "GB2016",
    "sequences.min_bins": "CHEN2016",
    "sequences.significance_z": "CHEN2016",
    "sequences.time_bin_s": "CHEN2016",
    "units.include_types": "GB2016",
}

EXTERNAL = ("CHEN2016", "GB2016", "DATA")


def tags_in_params() -> dict[str, str]:
    """``{dotted.parameter.path: TAG}`` for every externally attributed line."""
    text = (Path(__file__).resolve().parents[1] / "config" / "params.yaml").read_text()
    path: list[tuple[int, str]] = []
    found: dict[str, str] = {}
    for raw in text.splitlines():
        if not raw.strip() or raw.lstrip().startswith("#"):
            continue
        match = re.match(r"(\s*)([A-Za-z_][A-Za-z0-9_]*)\s*:(.*)$", raw)
        if not match:
            continue
        indent = len(match.group(1))
        while path and path[-1][0] >= indent:
            path.pop()
        path.append((indent, match.group(2)))
        tag = re.search(r"#\s*\[([A-Z0-9]+)\]", match.group(3))
        if tag and tag.group(1) in EXTERNAL:
            found[".".join(key for _, key in path)] = tag.group(1)
    return found


def test_no_parameter_gains_or_loses_an_attribution_unnoticed():
    found = tags_in_params()
    added = set(found) - set(ATTRIBUTED)
    removed = set(ATTRIBUTED) - set(found)
    assert not added, (
        f"newly attributed to a paper: {sorted(added)} -- open the paper, find the "
        "sentence, then add it to ATTRIBUTED")
    assert not removed, (
        f"no longer attributed: {sorted(removed)} -- if that is deliberate, remove it "
        "from ATTRIBUTED and say why in decisions.md")


def test_no_attribution_changes_which_paper_it_credits():
    found = tags_in_params()
    moved = {k: (ATTRIBUTED[k], found[k]) for k in set(found) & set(ATTRIBUTED)
             if ATTRIBUTED[k] != found[k]}
    assert not moved, f"attribution changed: {moved}"


def test_the_parser_finds_a_tag_at_every_nesting_depth():
    """Guard the guard: a path bug would silently shrink what is frozen."""
    found = tags_in_params()
    depths = {len(key.split(".")) for key in found}
    assert depths >= {2, 3, 4}, f"only found tags at depths {sorted(depths)}"
    assert found["lfp.channel_selection.bands.ripple"] == "CHEN2016"
    assert found["events.mua.threshold_high_sd"] == "CHEN2016"
    assert found["events.ripples.band"] == "CHEN2016"


def test_the_two_thresholds_that_were_wrong_are_not_attributed():
    """The boundary and ripple thresholds are ours; they must carry no citation."""
    found = tags_in_params()
    assert "events.mua.threshold_low_sd" not in found
    assert "events.ripples.threshold_high_sd" not in found
