"""Step 9c: per-cell contribution (PCC) on one session, with its null labels first.

    uv run python scripts/per_cell_contribution.py Achilles_10252013 --chunk 1/4
    uv run python scripts/per_cell_contribution.py Achilles_10252013 --summarise
    uv run python scripts/per_cell_contribution.py Achilles_10252013 --validate 40

The estimator and the reasoning behind every choice are in `hc11.pcc` and D9.3.
In short: shift one cell's place field to every other position on the track,
re-score every event of an epoch with that template, and take how far the real
placement sits above the average placement as that cell's contribution. Running
it for both sleep epochs gives PCC(PRE) and PCC(POST), and their difference is
the ΔPCC the paper splits at the median into rigid and plastic cells.

**The null labels are computed first and reported first.** Each shifted field
can stand in for the real one, which yields a ΔPCC the estimator produced from a
field carrying no positional information. If the distribution of those null ΔPCC
values is indistinguishable from the real one, the rigid/plastic split is
reproducible from noise and nothing downstream of it means anything -- and that
statement depends on no later choice, which is why it leads.

Three stages, because the computation is 77 cells x 15 shifts x 2471 events and
a single call has a time limit:

* ``--chunk k/n`` computes the shift means for a slice of the cells and writes
  results/pcc/<session>_shift_means_<k>of<n>.csv. Chunks are independent.
* ``--summarise`` reads every chunk, checks the coverage is complete, and
  derives PCC, the null labels, the two ΔPCC distributions, the rigid/plastic
  split and the figures.
* ``--validate N`` re-scores N random (cell, shift, event) triples with fresh
  500-shuffle nulls and compares them against the reused nulls, which is the one
  approximation in the estimator (D9.4).
"""

from __future__ import annotations

import argparse
import json
import logging
from datetime import datetime, timezone

import numpy as np
import pandas as pd
from scipy.stats import ks_2samp, mannwhitneyu

from hc11 import config, pcc, sequences, subsets
from place_fields import prepare  # noqa: I001  (sibling script)
from score_events import direction_tuning, place_cell_ids  # noqa: I001

log = logging.getLogger("pcc")

EPOCHS = ("PRE", "POST")
#: Largest disagreement tolerated between re-scoring the real template here and
#: the step-7 table. The step-7 CSV is written with six decimals, so a few times
#: 1e-6 is the floor; anything larger means the two code paths have diverged.
STEP7_TOLERANCE = 1e-3


def load_events(session, name, epoch, cells, params):
    """Prepared events for one epoch: binned counts plus step 7's scores and nulls."""
    cfg = params["sequences"]
    table = pd.read_csv(config.results_dir() / "sequences" / f"{name}_{epoch}.csv")
    table = table[table.scored].reset_index(drop=True)
    events = []
    for row in table.itertuples():
        counts, _ = sequences.bin_event(session, row.start, row.end, cells,
                                        cfg["time_bin_s"])
        keep = counts.sum(axis=0) > 0
        counts = counts[:, keep]
        events.append(pcc.PreparedEvent(
            peak=float(row.peak),
            counts=counts,
            times=np.flatnonzero(keep) * cfg["time_bin_s"],
            present=counts.sum(axis=1) > 0,
            observed_abs_r=abs(float(row.r_column_cycle)),
            null_mean={k: float(getattr(row, f"null_mean_{k}"))
                       for k in sequences.SHUFFLES},
            null_sd={k: float(getattr(row, f"null_sd_{k}"))
                     for k in sequences.SHUFFLES},
        ))
    observed = {k: table[f"rz_{k}"].to_numpy(dtype=float) for k in sequences.SHUFFLES}
    return events, observed


def placement_means_for_cell(events, observed, tuning, row, offsets, track, params):
    """Mean rZ over the epoch's events for every placement of one cell's field.

    Placement 0 is the cell's real field, and it is computed **through the same
    path as every shifted one**. The alternative -- taking the real mean from
    step 7's table and the shifted means from here -- would put any difference
    between the two code paths straight into PCC, for every cell, in the same
    direction. Step 7's numbers are instead used as an independent check:
    ``shift0_vs_step7_<kind>`` is the largest disagreement over the epoch's
    events, and it should be at the rounding floor of the step-7 CSV.

    One event mask serves every placement (:func:`hc11.pcc.finite_everywhere`),
    so a contribution is a difference between field positions and never between
    event sets.
    """
    cfg = params["sequences"]
    placements = [0, *offsets]
    by_kind = {kind: [] for kind in sequences.SHUFFLES}
    for offset in placements:
        moved = pcc.shift_field(tuning, row, offset)
        abs_r = np.array([
            pcc.best_abs_r(event.counts, moved, track.bin_centers, event.times,
                           cfg["time_bin_s"], circular=track.is_circular,
                           rate_floor_hz=cfg["rate_floor_hz"])
            for event in events
        ])
        for kind in sequences.SHUFFLES:
            by_kind[kind].append(np.array([
                pcc.standardise(r, e.null_mean[kind], e.null_sd[kind])
                for r, e in zip(abs_r, events, strict=True)
            ]))

    present = np.array([bool(event.present[row]) for event in events])
    records = [{"shift": int(offset)} for offset in placements]
    for kind in sequences.SHUFFLES:
        stack = np.vstack(by_kind[kind])
        keep = pcc.finite_everywhere(np.vstack([stack, observed[kind][None, :]]))
        with_spikes = keep & present
        for index, record in enumerate(records):
            values = stack[index]
            record[f"mean_rz_{kind}"] = (
                float(values[keep].mean()) if keep.any() else np.nan)
            record[f"mean_rz_present_{kind}"] = (
                float(values[with_spikes].mean()) if with_spikes.any() else np.nan)
            record[f"n_events_{kind}"] = int(keep.sum())
            record[f"n_present_{kind}"] = int(with_spikes.sum())
        deviation = (float(np.max(np.abs(stack[0][keep] - observed[kind][keep])))
                     if keep.any() else np.nan)
        for record in records:
            record[f"shift0_vs_step7_{kind}"] = deviation
    return records


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("session")
    parser.add_argument("--chunk", default=None, help="k/n: compute cells k of n chunks")
    parser.add_argument("--summarise", action="store_true")
    parser.add_argument("--validate", type=int, default=0)
    parser.add_argument("--no-figures", action="store_true")
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(message)s")

    params = config.load_params()
    name = args.session
    out_root = config.results_dir() / "pcc"
    out_root.mkdir(parents=True, exist_ok=True)

    if args.chunk or args.validate:
        session, track, running, laps = prepare(name, params)
        fields_table = pd.read_csv(config.results_dir() / "place_fields_summary.csv")
        cells = place_cell_ids(fields_table, name)
        tuning = direction_tuning(session, track, running, laps, cells, params)
        offsets = list(range(1, track.n_bins))
        log.info("%s: %d place cells, %d position bins, %d shifts per cell, "
                 "directions %s", name, len(cells), track.n_bins, len(offsets),
                 sorted(tuning))
        prepared = {epoch: load_events(session, name, epoch, cells, params)
                    for epoch in EPOCHS}
        for epoch in EPOCHS:
            log.info("  %s: %d scored events", epoch, len(prepared[epoch][0]))

    if args.chunk:
        k, n = (int(x) for x in args.chunk.split("/"))
        rows_out, worst = [], 0.0
        for row, cluster_id in enumerate(cells):
            if row % n != (k - 1) % n:
                continue
            for epoch in EPOCHS:
                events, observed = prepared[epoch]
                records = placement_means_for_cell(events, observed, tuning, row,
                                                   offsets, track, params)
                for record in records:
                    rows_out.append({"session": name, "cluster_id": int(cluster_id),
                                     "epoch": epoch, **record})
                deviations = [records[0][f"shift0_vs_step7_{kind}"]
                              for kind in sequences.SHUFFLES]
                worst = max([worst, *(d for d in deviations if np.isfinite(d))])
            log.info("  cell %d (%d/%d of this chunk) done", cluster_id,
                     len([r for r in rows_out if r["epoch"] == "PRE"]) // (len(offsets) + 1),
                     len(cells[(k - 1) % n::n]))
        if worst > STEP7_TOLERANCE:
            raise SystemExit(
                f"{name}: re-scoring the real template disagrees with step 7 by "
                f"{worst:.2e} rZ, above the {STEP7_TOLERANCE:.0e} tolerance -- the two "
                "code paths have diverged, so no contribution computed from them means "
                "anything")
        log.info("  real template reproduces step 7 to %.2e rZ "
                 "(the step-7 CSV's own rounding floor is 1e-6)", worst)
        path = out_root / f"{name}_shift_means_{k}of{n}.csv"
        pd.DataFrame(rows_out).round(6).to_csv(path, index=False)
        log.info("wrote %s", path.relative_to(config.repo_root()))

    if args.validate:
        rng = np.random.default_rng(params["sequences"]["shuffle_seed"])
        cfg = params["sequences"]
        rows_out = []
        for _ in range(args.validate):
            epoch = EPOCHS[rng.integers(len(EPOCHS))]
            events, _ = prepared[epoch]
            event = events[rng.integers(len(events))]
            row = int(rng.integers(len(cells)))
            offset = int(offsets[rng.integers(len(offsets))])
            moved = pcc.shift_field(tuning, row, offset)
            posteriors = {d: np.exp(
                sequences.decode_log_posterior(
                    event.counts, values, cfg["time_bin_s"],
                    rate_floor_hz=cfg["rate_floor_hz"])).T
                for d, values in moved.items()}
            perturbed = sequences.EventPosteriors(
                posteriors=posteriors,
                bin_index=np.rint(event.times / cfg["time_bin_s"]).astype(int),
                spikes_per_bin=event.counts.sum(axis=0),
                n_bins_total=event.counts.shape[1], n_cells_active=int(event.present.sum()),
            )
            record = {"session": name, "epoch": epoch, "peak": event.peak,
                      "cluster_id": int(cells[row]), "shift": offset}
            for kind in sequences.SHUFFLES:
                fresh = sequences.score_event(
                    perturbed, track.bin_centers, cfg["time_bin_s"],
                    circular=track.is_circular, shuffle=kind,
                    n_shuffles=cfg["n_shuffles"], rng=rng)
                reused = pcc.standardise(abs(fresh.r), event.null_mean[kind],
                                         event.null_sd[kind])
                record[f"fresh_rz_{kind}"] = fresh.rz
                record[f"reused_rz_{kind}"] = reused
                record[f"fresh_null_mean_{kind}"] = fresh.null_mean
                record[f"reused_null_mean_{kind}"] = event.null_mean[kind]
                record[f"fresh_null_sd_{kind}"] = fresh.null_sd
                record[f"reused_null_sd_{kind}"] = event.null_sd[kind]
            rows_out.append(record)
        table = pd.DataFrame(rows_out)
        path = config.results_dir() / f"pcc_null_reuse_{name}.csv"
        table.round(6).to_csv(path, index=False)
        log.info("\nnull-reuse validation on %d (cell, shift, event) triples", len(table))
        for kind in sequences.SHUFFLES:
            fresh, reused = table[f"fresh_rz_{kind}"], table[f"reused_rz_{kind}"]
            ok = np.isfinite(fresh) & np.isfinite(reused)
            difference = (reused - fresh)[ok]
            log.info("  %-18s r = %.4f, bias %+.4f rZ (median %+.4f), "
                     "95%% of |error| below %.4f, max %.4f",
                     kind, float(np.corrcoef(fresh[ok], reused[ok])[0, 1]),
                     float(difference.mean()), float(difference.median()),
                     float(np.percentile(np.abs(difference), 95)),
                     float(np.abs(difference).max()))
            sd_ratio = (table[f"reused_null_sd_{kind}"] / table[f"fresh_null_sd_{kind}"])[ok]
            log.info("  %-18s null sd ratio reused/fresh: median %.4f, "
                     "5-95%% %.4f-%.4f", "", float(sd_ratio.median()),
                     float(np.percentile(sd_ratio, 5)), float(np.percentile(sd_ratio, 95)))
            # A constant bias cancels in delta-PCC; one that differs by epoch does not.
            for epoch in EPOCHS:
                here = ok & (table.epoch == epoch)
                if here.sum() >= 2:
                    log.info("  %-18s %s bias %+.4f rZ over %d triples", "", epoch,
                             float((reused - fresh)[here].mean()), int(here.sum()))
        log.info("wrote %s", path.relative_to(config.repo_root()))

    if args.summarise:
        summarise(name, params, out_root, no_figures=args.no_figures)


def summarise(name, params, out_root, *, no_figures=False) -> None:
    """Derive PCC, the null labels and the two ΔPCC distributions from the chunks."""
    from hc11 import plots

    chunks = sorted(out_root.glob(f"{name}_shift_means_*.csv"))
    if not chunks:
        raise SystemExit(f"no shift-mean chunks for {name}; run --chunk first")
    table = pd.concat([pd.read_csv(p) for p in chunks], ignore_index=True)
    cells = sorted(table.cluster_id.unique())
    shifts = sorted(table.loc[table["shift"] != 0, "shift"].unique())
    expected = len(cells) * len(EPOCHS) * (len(shifts) + 1)
    if len(table) != expected:
        raise SystemExit(
            f"{name}: {len(table)} rows but {expected} expected "
            f"({len(cells)} cells x {len(EPOCHS)} epochs x {len(shifts) + 1} placements)"
            " -- a chunk is missing or was computed with different settings")
    log.info("%s: %d cells x %d shifts, complete", name, len(cells), len(shifts))

    rows, null_rows = [], []
    for kind in sequences.SHUFFLES:
        for cluster_id in cells:
            per_epoch, per_epoch_null = {}, {}
            for epoch in EPOCHS:
                here = table[(table.cluster_id == cluster_id) & (table.epoch == epoch)]
                real = here.loc[here["shift"] == 0, f"mean_rz_{kind}"].iloc[0]
                moved = here[here["shift"] != 0].sort_values("shift")
                values = moved[f"mean_rz_{kind}"].to_numpy(dtype=float)
                per_epoch[epoch] = pcc.contribution(real, values)
                at_zero = here["shift"] == 0
                per_epoch[epoch]["n_events"] = int(
                    here.loc[at_zero, f"n_events_{kind}"].iloc[0])
                per_epoch[epoch]["n_present"] = int(
                    here.loc[at_zero, f"n_present_{kind}"].iloc[0])
                per_epoch_null[epoch] = pcc.null_contributions(values)
            rows.append({
                "session": name, "shuffle": kind, "cluster_id": int(cluster_id),
                "pcc_pre": per_epoch["PRE"]["pcc"], "z_pre": per_epoch["PRE"]["z"],
                "p_pre": per_epoch["PRE"]["p"],
                "pcc_post": per_epoch["POST"]["pcc"], "z_post": per_epoch["POST"]["z"],
                "p_post": per_epoch["POST"]["p"],
                "delta_pcc": per_epoch["POST"]["pcc"] - per_epoch["PRE"]["pcc"],
                "n_events_pre": per_epoch["PRE"]["n_events"],
                "n_events_post": per_epoch["POST"]["n_events"],
                "n_present_pre": per_epoch["PRE"]["n_present"],
                "n_present_post": per_epoch["POST"]["n_present"],
            })
            for index, shift in enumerate(shifts):
                null_rows.append({
                    "session": name, "shuffle": kind, "cluster_id": int(cluster_id),
                    "shift": int(shift),
                    "pcc_pre": per_epoch_null["PRE"][index],
                    "pcc_post": per_epoch_null["POST"][index],
                    "delta_pcc": (per_epoch_null["POST"][index]
                                  - per_epoch_null["PRE"][index]),
                })
    real = pd.DataFrame(rows)
    nulls = pd.DataFrame(null_rows)
    real.round(6).to_csv(config.results_dir() / f"pcc_{name}.csv", index=False)
    nulls.round(6).to_csv(config.results_dir() / f"pcc_null_{name}.csv", index=False)

    stamp = {"script": "scripts/per_cell_contribution.py", "session": name,
             "git_commit": config.git_commit(),
             "created": datetime.now(timezone.utc).isoformat(timespec="seconds"),
             "pcc": params["pcc"], "sequences": params["sequences"]}
    (config.results_dir() / f"pcc_{name}.json").write_text(
        json.dumps(stamp, indent=2) + "\n")
    properties = ("session_rate_hz", "info_bits_per_spike", "n_fields", "peak_rate_hz")
    props = cell_properties(name, params)
    report(name, real, nulls, params)
    log.info("\n=== %s: rigid/plastic comparisons, null labels beside real ones", name)
    control = label_control(name, real, nulls, props, params, properties)
    if control:
        pd.DataFrame(control).round(6).to_csv(
            config.results_dir() / f"pcc_label_control_{name}.csv", index=False)
    if not no_figures:
        plots.plot_pcc_null_first(
            real.to_dict("records"), nulls.to_dict("records"),
            config.results_dir() / f"pcc_{name}.png", stamp)


def cell_properties(name, params) -> pd.DataFrame:
    """Per-cell quantities the paper's Fig. 3A and Fig. 4A compare, by cluster ID.

    ``session_rate_hz`` is the whole-session rate the paper uses for its
    rigid-versus-plastic firing-rate claim, not the on-maze running rate, which
    is a different number about a different thing.
    """
    session, track, running, laps = prepare(name, params)
    fields = pd.read_csv(config.results_dir() / "place_fields_summary.csv")
    here = fields[(fields.session == name) & fields.is_place_cell]
    whole = np.array([[0.0, session.sess_duration]])
    rows = []
    for cluster_id, group in here.groupby("cluster_id"):
        rows.append({
            "cluster_id": int(cluster_id),
            "session_rate_hz": subsets.rate_in_intervals(session, int(cluster_id), whole),
            "offline_rate_hz": subsets.rate_in_intervals(
                session, int(cluster_id),
                np.array([session.epoch(e) for e in EPOCHS], dtype=float)),
            "info_bits_per_spike": float(group.info_bits_per_spike.max()),
            "n_fields": float(group.n_fields.max()),
            "peak_rate_hz": float(group.peak_rate.max()),
        })
    return pd.DataFrame(rows).set_index("cluster_id")


def split_comparison(delta, ids, props, properties) -> dict:
    """Median split of ``delta`` into rigid and plastic, and the group differences."""
    delta = np.asarray(delta, dtype=float)
    median = float(np.median(delta))
    rigid, plastic = ids[delta <= median], ids[delta > median]
    out = {"median": median, "n_rigid": len(rigid), "n_plastic": len(plastic),
           "sum_plastic": float(delta[delta > median].sum()),
           "sum_rigid": float(delta[delta <= median].sum())}
    for name in properties:
        a = props[name].reindex(rigid).to_numpy(dtype=float)
        b = props[name].reindex(plastic).to_numpy(dtype=float)
        ok_a, ok_b = a[np.isfinite(a)], b[np.isfinite(b)]
        if ok_a.size < 2 or ok_b.size < 2:
            out[name] = (np.nan, np.nan, np.nan, np.nan)
            continue
        statistic, p = mannwhitneyu(ok_b, ok_a, alternative="two-sided")
        out[name] = (float(np.median(ok_a)), float(np.median(ok_b)),
                     float(2 * statistic / (ok_a.size * ok_b.size) - 1), float(p))
    return out


def label_control(name, real, nulls, props, params, properties) -> list[dict]:
    """Every rigid/plastic comparison, null labels beside real ones, never alone.

    The null labels come from the same estimator run on fields that carry no
    positional information, so they say what a median split of this many cells
    produces from noise. Where the null's effects are as large as the real
    ones, the real ones are not evidence -- and that statement is made before
    the real numbers are read, not after.
    """
    threshold = params["pcc"]["strong_z"]
    out = []
    for kind in sequences.SHUFFLES:
        here = real[real.shuffle == kind]
        strong = here[(here.z_pre >= threshold) | (here.z_post >= threshold)]
        ids = strong.cluster_id.to_numpy()
        if len(ids) < 6:
            log.info("  %s: only %d strongly contributing cells -- no split attempted",
                     kind, len(ids))
            continue
        observed = split_comparison(
            strong.delta_pcc.to_numpy(), ids,
            props, properties)
        null_here = nulls[(nulls.shuffle == kind) & nulls.cluster_id.isin(ids)]
        draws = []
        for _, group in null_here.groupby("shift"):
            ordered = group.set_index("cluster_id").reindex(ids)
            draws.append(split_comparison(
                ordered.delta_pcc.to_numpy(), ids, props, properties))
        log.info("\n  --- %s: %d strongly contributing cells, %d rigid / %d plastic",
                 kind, len(ids), observed["n_rigid"], observed["n_plastic"])
        log.info("      %-22s %-28s %s", "", "REAL labels", "NULL labels (15 shifts)")
        log.info("      %-22s plastic %+.4f rigid %+.4f    median plastic %+.4f "
                 "(range %+.4f..%+.4f)", "sum of delta-PCC",
                 observed["sum_plastic"], observed["sum_rigid"],
                 float(np.median([d["sum_plastic"] for d in draws])),
                 float(np.min([d["sum_plastic"] for d in draws])),
                 float(np.max([d["sum_plastic"] for d in draws])))
        for prop in properties:
            rigid_median, plastic_median, rb, p = observed[prop]
            null_rb = np.array([d[prop][2] for d in draws], dtype=float)
            null_p = np.array([d[prop][3] for d in draws], dtype=float)
            log.info("      %-22s %6.3f vs %6.3f  rb %+.3f p %.3f    |rb| up to %.3f, "
                     "%d of %d shifts reach p<0.05", prop, rigid_median, plastic_median,
                     rb, p, float(np.nanmax(np.abs(null_rb))),
                     int(np.sum(null_p < 0.05)), null_p.size)
            out.append({
                "session": name, "shuffle": kind, "property": prop,
                "n_rigid": observed["n_rigid"], "n_plastic": observed["n_plastic"],
                "real_rigid_median": rigid_median, "real_plastic_median": plastic_median,
                "real_rank_biserial": rb, "real_p": p,
                "null_rank_biserial_median": float(np.nanmedian(null_rb)),
                "null_rank_biserial_p95": float(np.nanpercentile(np.abs(null_rb), 95)),
                "null_rank_biserial_max": float(np.nanmax(np.abs(null_rb))),
                "null_shifts_significant": int(np.sum(null_p < 0.05)),
                "n_null_draws": int(null_p.size),
                "real_exceeds_null": bool(abs(rb) > np.nanpercentile(np.abs(null_rb), 95)),
            })
    return out


def report(name, real, nulls, params) -> None:
    """Print the null-label comparison first, the real numbers second."""

    threshold = params["pcc"]["strong_z"]
    for kind in sequences.SHUFFLES:
        r = real[real.shuffle == kind]
        n = nulls[nulls.shuffle == kind]
        real_delta = r.delta_pcc.to_numpy(dtype=float)
        null_delta = n.delta_pcc.to_numpy(dtype=float)
        real_delta = real_delta[np.isfinite(real_delta)]
        null_delta = null_delta[np.isfinite(null_delta)]
        log.info("\n=== %s, %s", name, kind)
        log.info("  NULL LABELS FIRST -- delta-PCC produced by shifted fields")
        log.info("    null  n=%5d  median %+.4f  IQR %+.4f..%+.4f  sd %.4f",
                 null_delta.size, np.median(null_delta),
                 np.percentile(null_delta, 25), np.percentile(null_delta, 75),
                 null_delta.std(ddof=1))
        log.info("    real  n=%5d  median %+.4f  IQR %+.4f..%+.4f  sd %.4f",
                 real_delta.size, np.median(real_delta),
                 np.percentile(real_delta, 25), np.percentile(real_delta, 75),
                 real_delta.std(ddof=1))
        stat, p = ks_2samp(real_delta, null_delta)
        u, pu = mannwhitneyu(real_delta, null_delta, alternative="two-sided")
        log.info("    two-sample KS D = %.4f, p = %.3g | Mann-Whitney p = %.3g, "
                 "rank-biserial %+.3f", stat, p, pu,
                 2 * u / (real_delta.size * null_delta.size) - 1)
        log.info("    spread ratio real/null sd = %.3f",
                 real_delta.std(ddof=1) / null_delta.std(ddof=1))

        # "Strongly contributing" is signed: a cell whose real field scores
        # *better* than its own alternative placements. Taking |z| would also
        # admit cells whose real placement scores worse than its shifts, which
        # is a real effect but not the one the paper's phrase describes. The
        # two-sided count is logged beside it so the choice is visible (D9.3).
        strong = r[(r.z_pre >= threshold) | (r.z_post >= threshold)]
        either_sign = r[(r.z_pre.abs() >= threshold) | (r.z_post.abs() >= threshold)]
        log.info("  then the real labels")
        log.info("    %d of %d cells reach z >= %.2f in either epoch "
                 "(%d if the sign is ignored); paper: 216 of 491 = 44%%",
                 len(strong), len(r), threshold, len(either_sign))
        if len(strong) >= 4:
            median = strong.delta_pcc.median()
            rigid = strong[strong.delta_pcc <= median]
            plastic = strong[strong.delta_pcc > median]
            log.info("    rigid %d cells (delta-PCC <= %+.4f), plastic %d",
                     len(rigid), median, len(plastic))
        contributing_pre = r[r.z_pre >= threshold]
        both = contributing_pre[contributing_pre.z_post >= threshold]
        if len(contributing_pre):
            log.info("    of %d cells contributing to PRE, %d (%.0f%%) also contribute "
                     "to POST (paper: 73%%)", len(contributing_pre), len(both),
                     100 * len(both) / len(contributing_pre))


if __name__ == "__main__":
    main()
