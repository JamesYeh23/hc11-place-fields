"""Step 6d: characterise the candidate events, and split the unmatched MUA events.

    uv run python scripts/characterise_events.py
    uv run python scripts/characterise_events.py Achilles_10252013 --no-figures

Three jobs.

**Characterisation.** Rate, duration and participation for each candidate-event
set, per epoch, with the number of events that carry enough cells to decode
reported separately for the full MUA set and for the strict conjunction (an MUA
event with a ripple partner). That is the choice step 7 inherits, so the cost of
the strict variant is given in usable events rather than in percentages alone.

**The unmatched MUA events.** Roughly 45 % of MUA events have no ripple partner.
Three explanations make different predictions and all three are tested:

1. *low-synchrony bursts* -- systematically fewer active cells than matched
   events;
2. *events on shanks the primary channel does not cover* -- normal
   participation, but active cells concentrated on particular shanks;
3. *detection noise* -- a higher failure rate against the spike-jitter null.

Prediction 2 needs care, because it is confounded with prediction 1: an event
with six active cells cannot occupy twelve shanks, so fewer cells forces higher
concentration whatever the physiology. Concentration is therefore compared
against a per-event null -- the same number of cells drawn at random from the
recorded population, respecting how many cells each shank actually has -- and it
is the *excess* over that null that is reported.

Reads results/mua_events.csv, results/ripple_events.csv,
results/mua_jitter_per_event.csv and results/place_fields_summary.csv. Writes
results/candidate_events.csv -- the set step 7 consumes, one row per MUA event
carrying its ripple match and multi-shank support -- plus
results/event_characterisation.csv, results/event_decodability.csv,
results/unmatched_mua_hypotheses.csv and per-session figures.
"""

from __future__ import annotations

import argparse
import json
import logging
from datetime import datetime, timezone

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu

from hc11 import config, plots
from hc11.events import EventList, participation, participation_by_shank
from hc11.io import load_session, shank_of
from hc11.preprocess import truncate_to_session
from hc11.results import write_session_table

log = logging.getLogger("characterise")

SESSIONS = list(config.LFP_SESSIONS)
EPOCHS = ("PRE", "MAZE", "POST")
#: Monte Carlo draws for the per-event concentration null.
N_NULL_DRAWS = 400


def event_list(table: pd.DataFrame, session: str, epoch: str, detector: str) -> EventList:
    rows = table[(table.session == session) & (table.epoch == epoch)]
    return EventList(
        session=session, epoch=epoch, detector=detector,
        start=rows.start.to_numpy(), peak=rows.peak.to_numpy(), end=rows.end.to_numpy(),
        peak_z=rows.peak_z.to_numpy(), baseline_mean=float("nan"),
        baseline_sd=float("nan"),
        state_seconds=float(rows.state_seconds.iloc[0]) if "state_seconds" in rows
        and len(rows) else float("nan"),
    )


def matched_mask(a: EventList, b: EventList, window_s: float) -> np.ndarray:
    """Per event of ``a``: does any event of ``b`` overlap it or peak within the window."""
    out = np.zeros(len(a), dtype=bool)
    if len(a) == 0 or len(b) == 0:
        return out
    for index in range(len(a)):
        out[index] = (
            ((b.start < a.end[index]) & (b.end > a.start[index]))
            | (np.abs(b.peak - a.peak[index]) <= window_s)
        ).any()
    return out


def subset(events: EventList, mask) -> EventList:
    from dataclasses import replace

    mask = np.asarray(mask)
    return replace(
        events, start=events.start[mask], peak=events.peak[mask], end=events.end[mask],
        peak_z=events.peak_z[mask],
    )


def place_cell_ids(table: pd.DataFrame, session: str) -> np.ndarray:
    """Cluster IDs that are a place cell in either direction (step 5a's ensemble)."""
    rows = table[(table.session == session) & table.is_place_cell]
    return np.array(sorted(rows.cluster_id.unique()))


def compare(label: str, values, mask, higher_means: str) -> dict:
    """Matched against unmatched on one quantity, with an effect size.

    Reports the rank-biserial correlation rather than a bare p-value: with
    thousands of events almost anything is significant, and the question is
    whether the difference is large enough to be an explanation.
    """
    values = np.asarray(values, dtype=float)
    yes, no = values[mask], values[~mask]
    finite_yes, finite_no = yes[np.isfinite(yes)], no[np.isfinite(no)]
    if finite_yes.size < 2 or finite_no.size < 2:
        return {"quantity": label, "n_matched": int(finite_yes.size),
                "n_unmatched": int(finite_no.size), "median_matched": np.nan,
                "median_unmatched": np.nan, "mean_matched": np.nan,
                "mean_unmatched": np.nan, "difference": np.nan,
                "rank_biserial": np.nan, "p": np.nan, "higher_means": higher_means}
    statistic, p = mannwhitneyu(finite_unmatched := finite_no, finite_yes,
                                alternative="two-sided")
    rank_biserial = 2 * statistic / (finite_unmatched.size * finite_yes.size) - 1
    return {
        "quantity": label,
        "n_matched": int(finite_yes.size), "n_unmatched": int(finite_no.size),
        "median_matched": float(np.median(finite_yes)),
        "median_unmatched": float(np.median(finite_no)),
        # The mean matters for the binary quantities, where a median of 1.0 on
        # both sides hides the whole difference.
        "mean_matched": float(finite_yes.mean()),
        "mean_unmatched": float(finite_no.mean()),
        "difference": float(np.median(finite_no) - np.median(finite_yes)),
        # Positive: unmatched events score higher than matched ones.
        "rank_biserial": float(rank_biserial), "p": float(p),
        "higher_means": higher_means,
    }


def concentration_null(n_active: np.ndarray, shank_sizes: np.ndarray, rng) -> np.ndarray:
    """Expected busiest-shank share if the same cells were drawn at random.

    Drawing respects how many cells each shank has, so a shank with twice the
    cells is twice as likely to supply one. Cached per distinct active-cell
    count: there are at most a few dozen of those and thousands of events.
    """
    population = np.repeat(np.arange(shank_sizes.size), shank_sizes)
    cache: dict[int, float] = {}
    out = np.empty(n_active.size)
    for index, k in enumerate(n_active.astype(int)):
        if k <= 0:
            out[index] = np.nan
            continue
        if k not in cache:
            draws = min(k, population.size)
            shares = [
                np.bincount(rng.choice(population, size=draws, replace=False),
                            minlength=shank_sizes.size).max() / draws
                for _ in range(N_NULL_DRAWS)
            ]
            cache[k] = float(np.mean(shares))
        out[index] = cache[k]
    return out


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("sessions", nargs="*", default=SESSIONS)
    parser.add_argument("--no-figures", action="store_true")
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(message)s")

    params = config.load_params()
    window_s = params["events"]["ripples"]["match_window_ms"] / 1000.0
    decodable = params["events"]["decodable"]
    out_root = config.results_dir()
    mua_table = pd.read_csv(out_root / "mua_events.csv")
    ripple_table = pd.read_csv(out_root / "ripple_events.csv")
    fields_table = pd.read_csv(out_root / "place_fields_summary.csv")
    jitter_table = pd.read_csv(out_root / "mua_jitter_per_event.csv")
    summary_table = pd.read_csv(out_root / "mua_summary.csv")
    channels_table = pd.read_csv(out_root / "ripple_channels.csv")

    stamp = {
        "generated": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "git_commit": config.git_commit(),
        "match_window_ms": params["events"]["ripples"]["match_window_ms"],
        "decodable": decodable,
    }
    char_rows, decode_rows, hypothesis_rows, shank_profile_rows = [], [], [], []
    candidate_rows = []

    names = args.sessions or SESSIONS
    for name in names:
        log.info("\n=== %s", name)
        session = load_session(name)
        if params["preprocess"]["truncate_to_session"]:
            session, _ = truncate_to_session(session, verbose=False)
        pyr = session.pyr_ids
        place = place_cell_ids(fields_table, name)
        # The shank the detection channel sits on, from step 6a-2's own table --
        # not inferred from the channel number, which would be wrong for the
        # four-channel "Extra Channels" groups.
        best = channels_table.query("session == @name and is_best_overall")
        primary_shank = int(best.shank.iloc[0]) if len(best) else -1
        log.info("%d pyramidal cells, %d place cells", len(pyr), len(place))

        for epoch in EPOCHS:
            mua = event_list(mua_table, name, epoch, "mua")
            lfp = event_list(ripple_table, name, epoch, "lfp")
            if len(mua) == 0:
                continue
            state_min = float(
                summary_table.query(
                    "session == @name and epoch == @epoch and variant == 'primary'"
                ).state_min.iloc[0]
            )
            matched = matched_mask(mua, lfp, window_s)
            variants = {
                "mua": mua,
                "conjunction": subset(mua, matched),
                "mua_without_ripple": subset(mua, ~matched),
            }

            for label, events in variants.items():
                if len(events) == 0:
                    continue
                part_pyr = participation(session, events, pyr)
                part_pc = participation(session, events, place)
                enough_pc = (part_pc.n_active >= decodable["min_active_cells"]) & (
                    part_pc.active_fraction >= decodable["min_active_fraction"]
                )
                char_rows.append({
                    "session": name, "epoch": epoch, "variant": label,
                    "n_events": len(events), "state_min": state_min,
                    "rate_per_min": len(events) / state_min if state_min else np.nan,
                    "median_duration_ms": 1000 * float(np.median(events.duration)),
                    "median_active_pyr": float(np.median(part_pyr.n_active)),
                    "median_active_place": float(np.median(part_pc.n_active)),
                    "median_peak_z": float(np.median(events.peak_z)),
                })
                decode_rows.append({
                    "session": name, "epoch": epoch, "variant": label,
                    "n_events": len(events),
                    "n_pyr": len(pyr), "n_place_cells": len(place),
                    "min_active_cells": decodable["min_active_cells"],
                    "min_active_fraction": decodable["min_active_fraction"],
                    "n_decodable": int(enough_pc.sum()),
                    "pct_decodable": 100 * float(enough_pc.mean()),
                    "decodable_per_min": (
                        enough_pc.sum() / state_min if state_min else np.nan
                    ),
                })

            # The set step 7 consumes: every MUA event, with the ripple match and
            # the multi-shank support alongside it, so the conjunction is a
            # filter on a column rather than a second pipeline.
            part_pyr_all = participation(session, mua, pyr)
            part_pc_all = participation(session, mua, place)
            support_lookup = dict(zip(
                ripple_table.query("session == @name and epoch == @epoch").peak.round(4),
                ripple_table.query(
                    "session == @name and epoch == @epoch"
                ).n_supporting_shanks,
                strict=True,
            ))
            for index in range(len(mua)):
                partner = np.nan
                support = np.nan
                if matched[index] and len(lfp):
                    near = np.abs(lfp.peak - mua.peak[index])
                    closest = int(np.argmin(near))
                    partner = float(lfp.peak[closest] - mua.peak[index])
                    support = support_lookup.get(round(float(lfp.peak[closest]), 4), np.nan)
                candidate_rows.append({
                    "session": name, "epoch": epoch,
                    "start": mua.start[index], "peak": mua.peak[index],
                    "end": mua.end[index],
                    "duration_ms": 1000 * mua.duration[index],
                    "peak_z": mua.peak_z[index],
                    "n_active_pyr": int(part_pyr_all.n_active[index]),
                    "n_active_place": int(part_pc_all.n_active[index]),
                    "decodable": bool(
                        part_pc_all.n_active[index] >= decodable["min_active_cells"]
                        and part_pc_all.active_fraction[index]
                        >= decodable["min_active_fraction"]
                    ),
                    "ripple_matched": bool(matched[index]),
                    "ripple_peak_offset_s": partner,
                    "n_supporting_shanks": support,
                })

            full = next(r for r in decode_rows
                        if r["session"] == name and r["epoch"] == epoch
                        and r["variant"] == "mua")
            strict = next((r for r in decode_rows
                           if r["session"] == name and r["epoch"] == epoch
                           and r["variant"] == "conjunction"), None)
            log.info(
                "  %-4s decodable: MUA %d/%d (%.0f%%)  conjunction %s  -- strict keeps %s",
                epoch, full["n_decodable"], full["n_events"], full["pct_decodable"],
                f"{strict['n_decodable']}/{strict['n_events']} ({strict['pct_decodable']:.0f}%)"
                if strict else "n/a",
                f"{100 * strict['n_decodable'] / max(full['n_decodable'], 1):.0f}% "
                f"of usable events" if strict else "n/a",
            )

            # --- the three hypotheses, PRE and POST only (MAZE has ~20 events) ---
            if epoch == "MAZE":
                continue
            part_pyr = participation(session, mua, pyr)
            part_pc = participation(session, mua, place)
            by_shank = participation_by_shank(session, mua, pyr)
            rng = np.random.default_rng(params["output"]["random_seed"])
            shank_sizes = np.array([
                int((shank_of(pyr) == s).sum()) for s in by_shank.shanks
            ])
            expected = concentration_null(part_pyr.n_active, shank_sizes, rng)
            primary_share = by_shank.fraction_on(primary_shank)
            expected_primary = (
                shank_sizes[np.flatnonzero(by_shank.shanks == primary_shank)[0]]
                / shank_sizes.sum()
                if (by_shank.shanks == primary_shank).any() else np.nan
            )

            jitter = jitter_table[
                (jitter_table.session == name) & (jitter_table.epoch == epoch)
            ]
            passes = np.full(len(mua), np.nan)
            if len(jitter):
                lookup = dict(zip(jitter.peak.round(4), jitter.exceeds_null_p95, strict=True))
                passes = np.array([
                    float(lookup.get(round(p, 4), np.nan)) for p in mua.peak
                ])

            tests = [
                compare("active_pyramidal_cells", part_pyr.n_active, matched,
                        "unmatched events are larger"),
                compare("active_place_cells", part_pc.n_active, matched,
                        "unmatched events are larger"),
                compare("duration_ms", 1000 * mua.duration, matched,
                        "unmatched events are longer"),
                compare("peak_z", mua.peak_z, matched, "unmatched events are stronger"),
                compare("shanks_with_an_active_cell", by_shank.n_shanks_active, matched,
                        "unmatched events are spread wider"),
                compare("busiest_shank_share", by_shank.concentration, matched,
                        "unmatched events are more concentrated"),
                compare("busiest_shank_share_excess", by_shank.concentration - expected,
                        matched, "unmatched events are more concentrated than chance"),
                compare("share_on_primary_shank", primary_share, matched,
                        "unmatched events favour the primary channel's shank"),
                compare("jitter_null_passed", passes, matched,
                        "unmatched events beat the jitter null more often"),
            ]
            for row in tests:
                hypothesis_rows.append({"session": name, "epoch": epoch, **row})
                log.info(
                    "    %-32s median %6.3f / %-6.3f  mean %6.3f / %-6.3f  "
                    "rb %+.3f  p %.2g",
                    row["quantity"], row["median_matched"], row["median_unmatched"],
                    row["mean_matched"], row["mean_unmatched"],
                    row["rank_biserial"], row["p"],
                )
            log.info("    (expected share on the primary shank by cell count: %.3f)",
                     expected_primary)

            for index, shank in enumerate(by_shank.shanks):
                total = by_shank.counts.sum(axis=1)
                with np.errstate(invalid="ignore"):
                    share = np.where(total > 0, by_shank.counts[:, index] / total, np.nan)
                shank_profile_rows.append({
                    "session": name, "epoch": epoch, "shank": int(shank),
                    "n_pyr_on_shank": int(shank_sizes[index]),
                    "is_primary_shank": int(shank) == primary_shank,
                    "share_matched": float(np.nanmean(share[matched])),
                    "share_unmatched": float(np.nanmean(share[~matched])),
                    "expected_share": float(shank_sizes[index] / shank_sizes.sum()),
                })

        if not args.no_figures:
            folder = config.results_dir(name)
            sidecar = {**stamp, "session": name}
            rows_here = [r for r in decode_rows if r["session"] == name]
            if rows_here:
                plots.plot_decodable_by_variant(
                    rows_here, name, folder / "event_decodability.png", sidecar
                )
            profile_here = [r for r in shank_profile_rows if r["session"] == name]
            if profile_here:
                plots.plot_shank_profile(
                    profile_here, name, folder / "unmatched_shank_profile.png", sidecar
                )

    write_session_table(out_root / "candidate_events.csv", candidate_rows, names, round_to=6)
    write_session_table(out_root / "event_characterisation.csv", char_rows, names, round_to=4)
    write_session_table(out_root / "event_decodability.csv", decode_rows, names, round_to=4)
    write_session_table(out_root / "unmatched_mua_hypotheses.csv", hypothesis_rows,
                        names, round_to=5)
    write_session_table(out_root / "unmatched_shank_profile.csv", shank_profile_rows,
                        names, round_to=5)
    (out_root / "event_characterisation.json").write_text(json.dumps(stamp, indent=2) + "\n")
    log.info("\nwrote %s and %s",
             (out_root / "candidate_events.csv").relative_to(config.repo_root()),
             (out_root / "event_characterisation.csv").relative_to(config.repo_root()))


if __name__ == "__main__":
    main()
