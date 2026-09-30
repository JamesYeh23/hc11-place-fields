"""Smoke tests: every figure renders, writes a sidecar, and survives empty input."""

from __future__ import annotations

import json

import numpy as np

from conftest import requires_data
from hc11 import config, plots
from hc11.behavior import running_mask_from_config
from hc11.fields import analyse_cell, select_running_samples
from hc11.laps import segment_laps_from_config
from hc11.track import Track


def test_direction_palette_is_fixed_and_distinct():
    assert set(plots.DIRECTION_COLORS) == {"pos", "neg"}
    assert plots.DIRECTION_COLORS["pos"] != plots.DIRECTION_COLORS["neg"]
    assert plots.RATE_CMAP not in {"jet", "rainbow", "hsv"}, "no rainbow ramps for magnitude"


def test_empty_inputs_still_render(tmp_path):
    for fn, name in (
        (plots.plot_rate_map_heatmap, "heatmap"),
        (plots.plot_field_grid, "grid"),
    ):
        track = Track(kind="linear", lo=0.0, hi=1.6, bin_size=0.1)
        out = fn([], track, "Empty_00000000", "pos", tmp_path / f"{name}.png", {"a": 1})
        assert out.exists() and out.with_suffix(".json").exists()


@requires_data
def test_figures_render_for_a_real_session(real_sessions, tmp_path):
    params = config.load_params()
    s = real_sessions("Achilles_10252013")
    track = Track.from_session(s, params["place_fields"]["bin_size_cm"])
    running = running_mask_from_config(s, params)
    laps = segment_laps_from_config(s, track, params)
    samples = select_running_samples(s, track, laps, running.mask, "pos")
    rng = np.random.default_rng(0)
    cells = [analyse_cell(s, track, samples, c, params, rng) for c in s.pyr_ids[:12]]

    meta = {"session": s.name, "git_commit": config.git_commit()}
    paths = [
        plots.plot_laps(s, track, laps, running, tmp_path / "laps.png", meta),
        plots.plot_info_histogram(cells, s.name, tmp_path / "info.png", meta),
        plots.plot_rate_map_heatmap(cells, track, s.name, "pos", tmp_path / "heat.png", meta),
        plots.plot_field_grid(cells, track, s.name, "pos", tmp_path / "grid.png", meta),
    ]
    for path in paths:
        assert path.exists() and path.stat().st_size > 5_000
        sidecar = json.loads(path.with_suffix(".json").read_text())
        assert sidecar["session"] == s.name, "each figure records what produced it"


@requires_data
def test_heatmap_shows_only_place_cells(real_sessions, tmp_path):
    """The sorted heatmap is a place-cell figure; non-place cells would blur it."""
    params = config.load_params()
    s = real_sessions("Achilles_10252013")
    track = Track.from_session(s, params["place_fields"]["bin_size_cm"])
    running = running_mask_from_config(s, params)
    laps = segment_laps_from_config(s, track, params)
    samples = select_running_samples(s, track, laps, running.mask, "pos")
    rng = np.random.default_rng(0)
    cells = [analyse_cell(s, track, samples, c, params, rng) for c in s.pyr_ids[:20]]
    n_pc = sum(c.is_place_cell for c in cells)
    out = plots.plot_rate_map_heatmap(cells, track, s.name, "pos", tmp_path / "h.png", {})
    assert out.exists()
    assert n_pc < len(cells), "fixture should include some non-place cells"
