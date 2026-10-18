"""Tests for merging per-session result tables."""

from __future__ import annotations

import pandas as pd
import pytest

from hc11.results import write_session_table


def rows_for(session, value):
    return [{"session": session, "epoch": "PRE", "value": value}]


def test_a_first_write_just_writes(tmp_path):
    path = tmp_path / "t.csv"
    write_session_table(path, rows_for("A", 1))
    assert pd.read_csv(path).session.tolist() == ["A"]


def test_a_second_session_is_added_not_substituted(tmp_path):
    """The footgun this exists for: running one session must not delete the rest."""
    path = tmp_path / "t.csv"
    write_session_table(path, rows_for("A", 1))
    write_session_table(path, rows_for("B", 2))
    table = pd.read_csv(path)
    assert sorted(table.session) == ["A", "B"]
    assert table.set_index("session").value.to_dict() == {"A": 1, "B": 2}


def test_rerunning_a_session_replaces_only_its_rows(tmp_path):
    path = tmp_path / "t.csv"
    write_session_table(path, rows_for("A", 1) + rows_for("B", 2))
    write_session_table(path, rows_for("A", 99))
    table = pd.read_csv(path).set_index("session")
    assert table.value.to_dict() == {"A": 99, "B": 2}


def test_an_explicit_session_set_clears_those_rows_even_with_none_produced(tmp_path):
    path = tmp_path / "t.csv"
    write_session_table(path, rows_for("A", 1) + rows_for("B", 2))
    write_session_table(path, [], sessions=["A"])
    assert pd.read_csv(path).session.tolist() == ["B"]


def test_without_an_explicit_set_an_empty_run_changes_nothing(tmp_path):
    path = tmp_path / "t.csv"
    write_session_table(path, rows_for("A", 1))
    write_session_table(path, [])
    assert pd.read_csv(path).session.tolist() == ["A"]


def test_rows_come_back_sorted_by_session(tmp_path):
    path = tmp_path / "t.csv"
    write_session_table(path, rows_for("B", 1))
    merged = write_session_table(path, rows_for("A", 2))
    assert merged.session.tolist() == ["A", "B"]


def test_a_table_with_no_session_column_is_replaced_wholesale(tmp_path):
    path = tmp_path / "t.csv"
    pd.DataFrame([{"k": 1}]).to_csv(path, index=False)
    write_session_table(path, [{"k": 2}])
    assert pd.read_csv(path).k.tolist() == [2]


def test_rounding_applies_to_the_merged_table(tmp_path):
    path = tmp_path / "t.csv"
    write_session_table(path, [{"session": "A", "value": 1.23456}], round_to=2)
    assert pd.read_csv(path).value.tolist() == [pytest.approx(1.23)]
