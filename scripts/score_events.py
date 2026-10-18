"""Steps 7a and 7b: decode candidate events and score their sequences.

    uv run python scripts/score_events.py                          # everything
    uv run python scripts/score_events.py Achilles_10252013 --epochs PRE
    uv run python scripts/score_events.py --identity-shuffle       # the control

Reads results/candidate_events.csv (step 6d) and the place-cell list from
results/place_fields_summary.csv. Each event is binned at
`sequences.time_bin_s`, decoded against the step-4 place fields of each running
direction, and scored by the weighted correlation between decoded position and
bin index, normalised against within-event shuffles of both kinds.

Output goes to results/sequences/<session>_<epoch>[_<tag>].csv, one file per
session and epoch, so a long run can be done in pieces and resumed. 7c reads
whatever is there.

``--identity-shuffle`` permutes which cell owns which place field before
decoding, leaving every other number in the pipeline untouched. It is the
control that says whether a PRE/POST difference lives in the replay or in the
event sets: with the fields detached from the cells, nothing should survive.
"""

from __future__ import annotations

import argparse
import json
import logging
import time
from datetime import datetime, timezone

import numpy as np
import pandas as pd

from hc11 import config, sequences
from hc11.decode import tuning_from_laps
from hc11.fields import select_running_samples
from hc11.laps import DIRECTIONS
from place_fields import prepare  # noqa: I001  (sibling script)

log = logging.getLogger("sequences")

SESSIONS = list(config.LFP_SESSIONS)
EPOCHS = ("PRE", "MAZE", "POST")


def place_cell_ids(table: pd.DataFrame, session: str) -> np.ndarray:
    """Place cells in either direction -- the same ensemble step 5b decoded with."""
    rows = table[(table.session == session) & table.is_place_cell]
    return np.array(sorted(rows.cluster_id.unique()))


def direction_tuning(session, track, running, laps, cluster_ids, params) -> dict:
    """``{direction: (n_cells, n_position)}`` from all MAZE running laps.

    Built from every lap, not leave-one-out: the held-out data here is sleep,
    which contributed nothing to the fields, so there is no lap to hold out.

    A direction with fewer than ``sequences.min_laps_per_direction`` laps is
    left out. Achilles_11012013 is a circular maze run one way, and its "pos"
    direction holds exactly one lap of 226 samples -- a template built from that
    is noise, and because the score maximises over directions, a noise template
    costs more than it looks: the shuffles get the same extra maximisation, so
    the null rises and every event's normalised score falls.
    """
    pf = params["place_fields"]
    minimum = params["sequences"]["min_laps_per_direction"]
    tuning = {}
    for direction in DIRECTIONS:
        samples = select_running_samples(session, track, laps, running.mask, direction)
        if samples.n == 0:
            continue
        n_laps = len({int(lap) for lap in samples.lap_of_sample})
        if n_laps < minimum:
            log.info(
                "  %s: %d lap(s) only, below the %d-lap floor -- not used as a template",
                direction, n_laps, minimum,
            )
            continue
        tuning[direction] = tuning_from_laps(
            session, track, samples, cluster_ids, np.ones(samples.n, dtype=bool),
            sigma_bins=pf["smoothing_sigma_bins"], min_occupancy_s=pf["min_occupancy_s"],
        )
    return tuning


def permute_identities(tuning: dict, rng) -> dict:
    """Give every cell someone else's place field, the same way in each direction.

    One permutation across all directions, because a cell is one cell: permuting
    independently per direction would also destroy the relationship between a
    cell's two directional fields, which is a second change on top of the one
    being tested.
    """
    any_tuning = next(iter(tuning.values()))
    order = rng.permutation(any_tuning.shape[0])
    return {direction: values[order] for direction, values in tuning.items()}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("sessions", nargs="*", default=SESSIONS)
    parser.add_argument("--epochs", nargs="*", default=list(EPOCHS))
    parser.add_argument("--identity-shuffle", action="store_true")
    parser.add_argument("--limit", type=int, default=0,
                        help="score only the first N events of each epoch (smoke test)")
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(message)s")

    params = config.load_params()
    cfg = params["sequences"]
    bin_s = cfg["time_bin_s"]
    tag = "identity_shuffle" if args.identity_shuffle else ""
    out_root = config.results_dir() / "sequences"
    out_root.mkdir(parents=True, exist_ok=True)
    candidates = pd.read_csv(config.results_dir() / "candidate_events.csv")
    fields_table = pd.read_csv(config.results_dir() / "place_fields_summary.csv")

    for name in args.sessions or SESSIONS:
        log.info("\n=== %s%s", name, "  [identity-shuffled fields]"
                 if args.identity_shuffle else "")
        session, track, running, laps = prepare(name, params)
        cells = place_cell_ids(fields_table, name)
        tuning = direction_tuning(session, track, running, laps, cells, params)
        if args.identity_shuffle:
            tuning = permute_identities(
                tuning, np.random.default_rng(cfg["identity_shuffle_seed"])
            )
        log.info(
            "%d place cells, directions %s, %d position bins of %g cm, circular=%s",
            len(cells), sorted(tuning), track.n_bins, 100 * track.bin_size,
            track.is_circular,
        )

        for epoch in args.epochs:
            rows_in = candidates[
                (candidates.session == name) & (candidates.epoch == epoch)
            ].reset_index(drop=True)
            if args.limit:
                rows_in = rows_in.iloc[: args.limit]
            if rows_in.empty:
                continue
            started = time.monotonic()
            rng = np.random.default_rng(cfg["shuffle_seed"])
            rows_out = []
            for row in rows_in.itertuples():
                event = sequences.event_posteriors(
                    session, row.start, row.end, cells, tuning, bin_s,
                    rate_floor_hz=cfg["rate_floor_hz"],
                )
                record = {
                    "session": name, "epoch": epoch, "peak": row.peak,
                    "start": row.start, "end": row.end,
                    "duration_ms": row.duration_ms,
                    "n_active_place": row.n_active_place,
                    "ripple_matched": row.ripple_matched,
                    "n_supporting_shanks": row.n_supporting_shanks,
                    "n_bins_total": event.n_bins_total,
                    "n_bins_used": event.n_bins,
                    "empty_bin_fraction": event.empty_fraction,
                    "spikes_per_bin": (
                        float(np.mean(event.spikes_per_bin)) if event.n_bins else np.nan
                    ),
                    "scored": False,
                }
                if event.n_bins >= cfg["min_bins"]:
                    record["scored"] = True
                    for kind in sequences.SHUFFLES:
                        score = sequences.score_event(
                            event, track.bin_centers, bin_s,
                            circular=track.is_circular, shuffle=kind,
                            n_shuffles=cfg["n_shuffles"], rng=rng,
                        )
                        record[f"r_{kind}"] = score.r
                        record[f"rz_{kind}"] = score.rz
                        record[f"null_mean_{kind}"] = score.null_mean
                        record[f"null_sd_{kind}"] = score.null_sd
                        record[f"direction_{kind}"] = score.direction
                        record[f"offset_{kind}"] = score.offset
                rows_out.append(record)

            table = pd.DataFrame(rows_out)
            stem = f"{name}_{epoch}" + (f"_{tag}" if tag else "")
            table.round(6).to_csv(out_root / f"{stem}.csv", index=False)
            scored = table[table.scored]
            log.info(
                "  %-4s %5d events, %5d scored (%.0f%%): bins/event %.1f (used %.1f), "
                "spikes/bin %.2f, empty bins %.0f%%  [%.0f s]",
                epoch, len(table), len(scored),
                100 * len(scored) / max(len(table), 1),
                table.n_bins_total.mean(), table.n_bins_used.mean(),
                scored.spikes_per_bin.mean() if len(scored) else np.nan,
                100 * table.empty_bin_fraction.mean(),
                time.monotonic() - started,
            )
            if len(scored):
                for kind in sequences.SHUFFLES:
                    log.info(
                        "       %-18s median rZ %+.3f, %.0f%% above 1.65, "
                        "%.0f%% forward",
                        kind, scored[f"rz_{kind}"].median(),
                        100 * (scored[f"rz_{kind}"] > 1.65).mean(),
                        100 * (scored[f"r_{kind}"] > 0).mean(),
                    )

    (out_root / f"run{'_' + tag if tag else ''}.json").write_text(json.dumps({
        "script": "scripts/score_events.py", "git_commit": config.git_commit(),
        "created": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "sequences": cfg, "identity_shuffle": args.identity_shuffle,
    }, indent=2) + "\n")
    log.info("\nwrote %s", out_root.relative_to(config.repo_root()))


if __name__ == "__main__":
    main()
