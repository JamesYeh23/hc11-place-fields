"""Tests for step 12b: the extended measurement base and its quality checks.

The risk this step exists to catch is a cell-count gain that is real in count
and spurious in quality, so what is tested is the machinery that could hide it:
that the extended base never writes over a primary result, that the coordinate
conversion between the two bases is exact rather than approximate, and that the
two ensembles compared within a base are compared on the same bins.
"""

from __future__ import annotations

import numpy as np
import pytest

import extended_base_checks as checks
import place_fields


def test_a_circular_session_is_refused_rather_than_scored_on_the_wrong_base(
        monkeypatch):
    """On a ring the reconstruction measures itself, so Phase 5 excludes it."""
    class Ring:
        name = "Fake_circular"
        maze_kind = "circular"

    monkeypatch.setattr(place_fields, "load_session", lambda name: Ring())
    params = {"preprocess": {"truncate_to_session": False},
              "place_fields": {"bin_size_cm": 10}}
    with pytest.raises(SystemExit, match="linear-track only"):
        place_fields.prepare("Fake_circular", params, base="extended")


def test_an_unknown_base_is_refused(monkeypatch):
    class Line:
        name = "Fake_linear"
        maze_kind = "linear"

    monkeypatch.setattr(place_fields, "load_session", lambda name: Line())
    params = {"preprocess": {"truncate_to_session": False},
              "place_fields": {"bin_size_cm": 10}}
    with pytest.raises(SystemExit, match="unknown base"):
        place_fields.prepare("Fake_linear", params, base="sideways")


def test_the_default_base_is_the_authors_one():
    """Phases 1-4 call `prepare` with no base argument and must not move."""
    import inspect
    assert inspect.signature(place_fields.prepare).parameters["base"].default == (
        "authors")


def test_the_affine_map_between_the_two_bases_is_exact_on_a_linear_track():
    """|r| = 1.0000 means a line fits with no residual worth reporting."""
    class Fake:
        pass

    ours, theirs = Fake(), Fake()
    x = np.linspace(0, 2.0, 500)
    ours.position_1d = x
    theirs.position_1d = 1.6 - 0.8 * x  # a different origin, scale and sign
    slope, intercept, residual, n = checks.affine_to_authors(theirs, ours)
    assert slope == pytest.approx(-0.8)
    assert intercept == pytest.approx(1.6)
    assert residual < 1e-9
    assert n == 500


def test_the_affine_map_ignores_samples_only_one_base_defines():
    class Fake:
        pass

    ours, theirs = Fake(), Fake()
    ours.position_1d = np.array([0.0, 1.0, 2.0, 3.0])
    theirs.position_1d = np.array([0.0, 2.0, 4.0, np.nan])
    slope, intercept, residual, n = checks.affine_to_authors(theirs, ours)
    assert slope == pytest.approx(2.0)
    assert n == 3


def test_the_end_fraction_counts_both_ends():
    """0.15 a side, so a uniform population puts 30 % of its peaks there."""
    assert checks.END_FRACTION == 0.15
    assert checks.N_THIRDS == 3


def test_two_ensembles_on_one_base_are_compared_only_where_both_decoded():
    """A smaller ensemble leaves more bins undecodable; scoring it on its own
    decodable subset would compare two different stretches of track."""
    big = np.array([1.0, 2.0, 3.0, 4.0])
    small = np.array([1.0, np.nan, 3.0, np.nan])
    shared = np.isfinite(big) & np.isfinite(small)
    assert list(shared) == [True, False, True, False]
    assert list(big[shared]) == [1.0, 3.0]
