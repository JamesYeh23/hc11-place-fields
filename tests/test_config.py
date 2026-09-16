"""Tests for path resolution and parameter loading."""

from __future__ import annotations

import pytest

from hc11 import config


def test_repo_root_contains_pyproject():
    assert (config.repo_root() / "pyproject.toml").is_file()


def test_config_file_exists():
    assert config.config_path().is_file()


def test_sessions_inventory():
    assert len(config.SESSIONS) == 8
    assert len(set(config.SESSIONS)) == 8
    assert len(set(config.SESSION_ANIMAL.values())) == 4


def test_params_load_and_have_required_sections():
    params = config.load_params()
    for section in (
        "dataset",
        "behavior",
        "units",
        "place_fields",
        "decoding",
        "laps",
        "output",
    ):
        assert section in params, f"missing config section: {section}"


@pytest.mark.parametrize(
    ("key", "expected"),
    [
        ("behavior.speed_threshold_cm_s", 15.0),
        ("place_fields.bin_size_cm", 10.0),
        ("decoding.prior", "uniform"),
        ("place_fields.expected_total_place_cells", 491),
    ],
)
def test_param_lookup(key, expected):
    assert config.param(key) == expected


def test_param_missing_key_raises_with_full_path():
    with pytest.raises(KeyError, match="place_fields.not_a_parameter"):
        config.param("place_fields.not_a_parameter")
    with pytest.raises(KeyError, match="no_such_section.x"):
        config.param("no_such_section.x")


def test_param_accepts_explicit_params_dict():
    params = {"a": {"b": 3}}
    assert config.param("a.b", params) == 3


def test_data_dir_not_required_returns_a_path():
    # Must not raise even if the data is absent.
    assert config.data_dir(required=False).name


def test_git_commit_is_a_string():
    assert isinstance(config.git_commit(), str)


def test_results_dir_is_under_repo_root():
    path = config.results_dir("unit_test_session", create=False)
    assert path.parent.name == "results"
    assert config.repo_root() in path.parents
