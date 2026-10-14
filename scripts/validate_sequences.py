"""Two controls on the sequence score, before any of it is believed.

    uv run python scripts/validate_sequences.py

**Negative control — identity-shuffled place fields.** Run
``scripts/score_events.py --identity-shuffle`` first; this script reads the
result. Permuting which cell owns which field leaves every other number in the
pipeline alone: same events, same bins, same spikes, same shuffles. If a
PRE/POST difference survives that, it is not replay — it is in the event sets or
in the scoring, and the Fig. 1C comparison means nothing.

**Positive control — MAZE running.** Decode actual running the same way and
score it. The animal's position is known and the decoder was shown in step 5b to
recover it to within 5 cm, so these must score high; if they do not, the scoring
is broken whatever the sleep events do.

Running is decoded with **leave-one-lap-out templates**, as in step 5b: a window
inside lap *i* is decoded against fields built from every lap but *i*. Scoring
running against fields that include it would be the leakage step 5b went to
some trouble to avoid.

Two window lengths, and the reason is not pedantry. At roughly 30 cm/s a window
matched to a candidate event -- about 200 ms -- covers some 6 cm, which is
under one 10 cm position bin. Such a window *cannot* show a trajectory however
good the decoder is, so a low score there would say nothing. The 2 s window
covers most of the track and is the control that can actually fail.

Writes results/sequence_validation.csv and results/sequence_running_scores.csv.
"""

from __future__ import annotations

import argparse
import json
import logging
from datetime import datetime, timezone

import numpy as np
import pandas as pd

from hc11 import config, sequences
from hc11.decode import tuning_from_laps
from hc11.fields import select_running_samples
from hc11.laps import DIRECTIONS
from place_fields import prepare  # noqa: I001  (sibling script)
from score_events import place_cell_ids  # noqa: I001  (sibling script)

log = logging.getLogger("validate")

SESSIONS = ["Achilles_10252013", "Achilles_11012013"]
EPOCHS = ("PRE", "POST")


def running_windows(session, track, running, laps, params, window_s):
    """Yield ``(lap_index, direction, start, end)`` for contiguous running windows.

    Windows are cut inside one lap, never across a lap boundary: a window
    spanning the turnaround would contain two trajectories and score as neither.
    """
    minimum = params["sequences"]["min_laps_per_direction"]
    for direction in DIRECTIONS:
        samples = select_running_samples(session, track, laps, running.mask, direction)
        if samples.n == 0:
            continue
        lap_ids = np.unique(samples.lap_of_sample)
        if len(lap_ids) < minimum:
            continue
        for lap in lap_ids:
            here = samples.lap_of_sample == lap
            times = session.position_t[samples.indices[here]]
            if times.size < 2:
                continue
            start, stop = float(times[0]), float(times[-1])
            edge = start
            while edge + window_s <= stop:
                yield int(lap), direction, edge, edge + window_s
                edge += window_s


def leave_one_out_tuning(session, track, running, laps, cells, params):
    """``{direction: {lap: (n_cells, n_position)}}`` plus the per-direction samples."""
    pf = params["place_fields"]
    minimum = params["sequences"]["min_laps_per_direction"]
    out, sample_by_direction = {}, {}
    for direction in DIRECTIONS:
        samples = select_running_samples(session, track, laps, running.mask, direction)
        if samples.n == 0:
            continue
        lap_ids = np.unique(samples.lap_of_sample)
        if len(lap_ids) < minimum:
            continue
        sample_by_direction[direction] = samples
        out[direction] = {}
        for lap in lap_ids:
            held_out = samples.lap_of_sample == lap
            if held_out.sum() < 2 or (~held_out).sum() < 2:
                continue
            out[direction][int(lap)] = tuning_from_laps(
                session, track, samples, cells, ~held_out,
                sigma_bins=pf["smoothing_sigma_bins"],
                min_occupancy_s=pf["min_occupancy_s"],
            )
    return out, sample_by_direction


def median_rz(rows: pd.DataFrame, label: str, epoch: str) -> float:
    """Median rZ for one set and epoch, or NaN when that row is absent."""
    hit = rows[(rows.set == label) & (rows.epoch == epoch)]
    return float(hit.median_rz.iloc[0]) if len(hit) else float("nan")


def summarise(scores: pd.DataFrame, label: str, kind: str) -> dict:
    column = f"rz_{kind}"
    values = scores[column].dropna()
    return {
        "set": label, "shuffle": kind, "n": int(values.size),
        "median_rz": float(values.median()) if values.size else np.nan,
        "mean_rz": float(values.mean()) if values.size else np.nan,
        "pct_above_z": 100 * float((values > 1.65).mean()) if values.size else np.nan,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("sessions", nargs="*", default=SESSIONS)
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(message)s")

    params = config.load_params()
    cfg = params["sequences"]
    bin_s = cfg["time_bin_s"]
    out_root = config.results_dir()
    sequence_dir = out_root / "sequences"
    fields_table = pd.read_csv(out_root / "place_fields_summary.csv")
    rows, running_rows = [], []

    for name in args.sessions or SESSIONS:
        log.info("\n=== %s", name)

        # --- negative control: identity-shuffled place fields ----------------
        for epoch in EPOCHS:
            real = sequence_dir / f"{name}_{epoch}.csv"
            shuffled = sequence_dir / f"{name}_{epoch}_identity_shuffle.csv"
            if not (real.exists() and shuffled.exists()):
                log.warning("  %s: missing %s -- run score_events.py --identity-shuffle",
                            epoch, shuffled.name)
                continue
            for path, label in ((real, "real fields"), (shuffled, "identity-shuffled")):
                table = pd.read_csv(path)
                table = table[table.scored]
                for kind in sequences.SHUFFLES:
                    rows.append({
                        "session": name, "epoch": epoch, "control": "identity",
                        **summarise(table, label, kind),
                    })
        for kind in sequences.SHUFFLES:
            here = [r for r in rows if r["session"] == name and r["shuffle"] == kind
                    and r["control"] == "identity"]
            for label in ("real fields", "identity-shuffled"):
                values = {r["epoch"]: r["median_rz"] for r in here if r["set"] == label}
                if values:
                    log.info(
                        "  %-18s %-16s median rZ  PRE %+0.3f  POST %+0.3f  "
                        "(POST-PRE %+0.3f)",
                        kind, label, values.get("PRE", np.nan),
                        values.get("POST", np.nan),
                        values.get("POST", np.nan) - values.get("PRE", np.nan),
                    )

        # --- positive control: MAZE running ----------------------------------
        session, track, running, laps = prepare(name, params)
        cells = place_cell_ids(fields_table, name)
        tuning, _ = leave_one_out_tuning(session, track, running, laps, cells, params)
        for window_s in cfg["running_window_s"]:
            rng = np.random.default_rng(cfg["shuffle_seed"])
            scored = []
            for lap, direction, start, end in running_windows(
                session, track, running, laps, params, window_s
            ):
                templates = {
                    other: by_lap[lap]
                    for other, by_lap in tuning.items() if lap in by_lap
                }
                if not templates:
                    continue
                event = sequences.event_posteriors(
                    session, start, end, cells, templates, bin_s,
                    rate_floor_hz=cfg["rate_floor_hz"],
                )
                if event.n_bins < cfg["min_bins"]:
                    continue
                record = {
                    "session": name, "window_s": window_s, "lap": lap,
                    "direction": direction, "start": start,
                    "n_bins_used": event.n_bins,
                    "spikes_per_bin": float(np.mean(event.spikes_per_bin)),
                }
                for kind in sequences.SHUFFLES:
                    score = sequences.score_event(
                        event, track.bin_centers, bin_s, circular=track.is_circular,
                        shuffle=kind, n_shuffles=cfg["n_shuffles"], rng=rng,
                    )
                    record[f"r_{kind}"] = score.r
                    record[f"rz_{kind}"] = score.rz
                scored.append(record)
            if not scored:
                log.warning("  running %.1f s: no windows survived the bin floor", window_s)
                continue
            table = pd.DataFrame(scored)
            running_rows.extend(scored)
            for kind in sequences.SHUFFLES:
                summary = summarise(table, f"MAZE running {window_s:g} s", kind)
                rows.append({"session": name, "epoch": "MAZE_RUN",
                             "control": "running", **summary})
                log.info(
                    "  running %4.1f s  %-18s n %4d  median rZ %+0.2f  %.0f%% above 1.65",
                    window_s, kind, summary["n"], summary["median_rz"],
                    summary["pct_above_z"],
                )

    table = pd.DataFrame(rows)
    table.round(4).to_csv(out_root / "sequence_validation.csv", index=False)
    pd.DataFrame(running_rows).round(6).to_csv(
        out_root / "sequence_running_scores.csv", index=False)

    # The verdict the controls exist to give, stated rather than left to the reader.
    log.info("\n--- verdict ---")
    for name in args.sessions or SESSIONS:
        for kind in sequences.SHUFFLES:
            rows_here = table[(table.session == name) & (table.shuffle == kind)
                              & (table.control == "identity")]
            real_gap = (median_rz(rows_here, "real fields", "POST")
                        - median_rz(rows_here, "real fields", "PRE"))
            null_gap = (median_rz(rows_here, "identity-shuffled", "POST")
                        - median_rz(rows_here, "identity-shuffled", "PRE"))
            log.info(
                "%s %-18s POST-PRE: real %+0.3f, identity-shuffled %+0.3f -- %s",
                name, kind, real_gap, null_gap,
                "control passes" if abs(null_gap) < 0.25 * abs(real_gap) or
                np.sign(null_gap) != np.sign(real_gap)
                else "CONTROL FAILS: the gap survives identity shuffling",
            )

    (out_root / "sequence_validation.json").write_text(json.dumps({
        "script": "scripts/validate_sequences.py", "git_commit": config.git_commit(),
        "created": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "sequences": cfg,
    }, indent=2) + "\n")
    log.info("wrote %s",
             (out_root / "sequence_validation.csv").relative_to(config.repo_root()))


if __name__ == "__main__":
    main()
