"""Bayesian position decoding during running, cross-validated lap by lap.

    uv run python scripts/decode_position.py                    # all sessions
    uv run python scripts/decode_position.py Achilles_10252013
    uv run python scripts/decode_position.py --no-figures

For every session, every bin size in `decoding.time_bins_s`, and every
condition (each running direction against its own fields, plus a
direction-agnostic decoder whose fields are merged across directions), this
decodes each lap against fields built from the other laps and reports the
error, a chance baseline, and the dependence on ensemble size.

The ensemble is the session's unique place cells, identical in all three
conditions, so the directional and merged rows differ only in the template.

Writes results/decoding_summary.csv, results/decoding_ensemble.csv and
per-session figures in results/<session>/.
"""

from __future__ import annotations

import argparse
import json
import logging
from datetime import datetime, timezone

import numpy as np
import pandas as pd

from hc11 import config, plots
from hc11.decode import chance_baseline, decode, error_vs_ensemble_size, prepare_decoder
from hc11.fields import analyse_cell, select_running_samples
from hc11.io import list_sessions
from hc11.laps import DIRECTIONS
from place_fields import prepare  # noqa: I001  (sibling script)

log = logging.getLogger("decode")


def place_cells(session, track, running, laps, params) -> tuple[np.ndarray, dict]:
    """Unique place-cell IDs for the session, plus the per-direction counts."""
    rng = np.random.default_rng(params["place_fields"]["shuffle"]["random_seed"])
    unique, per_direction = set(), {}
    for direction in DIRECTIONS:
        samples = select_running_samples(session, track, laps, running.mask, direction)
        if samples.n == 0:
            per_direction[direction] = 0
            continue
        cells = [analyse_cell(session, track, samples, c, params, rng) for c in session.pyr_ids]
        hits = {c.cluster_id for c in cells if c.is_place_cell}
        per_direction[direction] = len(hits)
        unique |= hits
    return np.array(sorted(unique)), per_direction


def main() -> None:
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("sessions", nargs="*")
    parser.add_argument("--no-figures", action="store_true")
    args = parser.parse_args()

    params = config.load_params()
    dec = params["decoding"]
    names = args.sessions or list_sessions()
    rows, ensemble_rows = [], []

    for name in names:
        session, track, running, laps = prepare(name, params)
        cells, per_direction = place_cells(session, track, running, laps, params)
        if len(cells) < dec["min_units_for_decoding"]:
            log.warning("%s: only %d place cells, below min_units_for_decoding=%d -- skipped",
                        name, len(cells), dec["min_units_for_decoding"])
            continue

        conditions = [(d, d) for d in DIRECTIONS] + [("merged", None)]
        figures: dict = {"results": [], "chance": {}, "curves": {}}
        for bin_s in dec["time_bins_s"]:
            for label, direction in conditions:
                samples = select_running_samples(session, track, laps, running.mask, direction)
                if samples.n == 0:
                    continue
                inputs = prepare_decoder(
                    session, track, samples, laps, cells, params, time_bin_s=bin_s
                )
                if not inputs.laps:
                    continue
                result = decode(inputs, track, params, keep_posteriors=True)
                chance = chance_baseline(inputs, track, params, np.random.default_rng(7))
                curve = error_vs_ensemble_size(inputs, track, params, np.random.default_rng(8))
                leaky = decode(inputs, track, params, tuning_override=inputs.full_tuning)

                pcts = result.error_percentiles()
                positions, by_position = result.error_by_position()
                rows.append({
                    "session": name, "maze_type": session.maze_type, "condition": label,
                    "time_bin_s": bin_s, "n_cells": result.n_cells,
                    "n_place_cells_this_direction": per_direction.get(label, len(cells)),
                    "n_laps": result.n_laps, "n_time_bins": int(result.error_cm.size),
                    "median_error_cm": result.median_error_cm,
                    "mean_error_cm": result.mean_error_cm,
                    **{f"p{q}_error_cm": v for q, v in pcts.items()},
                    "fraction_undecodable": result.fraction_undecodable,
                    "median_spikes_per_bin": float(np.median(result.n_spikes_per_bin)),
                    "chance_median_error_cm": float(np.median(chance)),
                    "chance_min_cm": float(np.min(chance)), "chance_max_cm": float(np.max(chance)),
                    "leaky_median_error_cm": leaky.median_error_cm,
                    "error_end_zones_cm": float(np.nanmean(by_position[[0, -1]]))
                    if by_position.size else np.nan,
                    "error_middle_cm": float(np.nanmean(by_position[1:-1]))
                    if by_position.size > 2 else np.nan,
                    "track_extent_m": track.extent,
                })
                for entry in curve:
                    ensemble_rows.append({
                        "session": name, "condition": label, "time_bin_s": bin_s, **entry
                    })
                if label == "pos" or (label == "merged" and not per_direction.get("pos")):
                    figures["results"].append(result)
                    figures["chance"][bin_s] = float(np.median(chance))
                    figures["curves"][bin_s] = curve
                    if bin_s == dec["time_bin_s"]:
                        figures["primary"] = (result, inputs)

        log.info("%s: %d place cells, %s", name, len(cells),
                 ", ".join(f"{r['condition']}@{r['time_bin_s'] * 1000:.0f}ms "
                           f"{r['median_error_cm']:.1f}cm"
                           for r in rows if r["session"] == name))

        if not args.no_figures and figures["results"]:
            out = config.results_dir(name)
            meta = {
                "session": name, "git_commit": config.git_commit(),
                "created": datetime.now(timezone.utc).isoformat(timespec="seconds"),
                "params": {"decoding": dec, "place_fields": params["place_fields"]},
                "n_place_cells": len(cells),
            }
            plots.plot_error_distribution(
                figures["results"], figures["chance"], name,
                out / "decoding_error.png", meta)
            plots.plot_error_vs_ensemble(
                figures["curves"], name, out / "decoding_ensemble.png", meta)
            if "primary" in figures:
                result, inputs = figures["primary"]
                plots.plot_decoding_examples(
                    result, inputs, track, out / "decoding_examples.png", meta)
                plots.plot_confusion(result, track, out / "decoding_confusion.png", meta)

    df = pd.DataFrame(rows)
    results = config.results_dir()
    df.round(4).to_csv(results / "decoding_summary.csv", index=False)
    pd.DataFrame(ensemble_rows).round(4).to_csv(results / "decoding_ensemble.csv", index=False)
    (results / "decoding_summary.json").write_text(json.dumps({
        "script": "scripts/decode_position.py", "git_commit": config.git_commit(),
        "created": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "sessions": names, "params": params,
    }, indent=2) + "\n")

    print()
    header = (f"{'session':<18} {'cond':<7} {'bin':>6} {'cells':>6} {'laps':>5} "
              f"{'median':>8} {'mean':>7} {'p90':>7} {'chance':>8} {'undec':>7}")
    print(header)
    print("-" * len(header))
    for _, r in df.iterrows():
        print(f"{r.session:<18} {r.condition:<7} {r.time_bin_s * 1000:5.0f}ms {r.n_cells:6d} "
              f"{r.n_laps:5d} {r.median_error_cm:7.1f} {r.mean_error_cm:6.1f} "
              f"{r.p90_error_cm:6.1f} {r.chance_median_error_cm:7.1f} "
              f"{100 * r.fraction_undecodable:6.1f}%")
    print("-" * len(header))
    print("(errors in cm; 'chance' = tuning curves permuted between cells)")
    print(f"\nwrote {(results / 'decoding_summary.csv').relative_to(config.repo_root())}")


if __name__ == "__main__":
    main()
