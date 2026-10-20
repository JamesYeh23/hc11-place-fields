"""Step 7c: the Fig. 1C comparison — rZ distributions for PRE, MAZE and POST.

    uv run python scripts/replay_comparison.py
    uv run python scripts/replay_comparison.py --no-figures

Reads results/sequences/<session>_<epoch>.csv and reports, for each session,
each shuffle kind and each candidate-event set:

* the cumulative distribution of rZ with medians marked;
* the identity-shuffle residual beside every effect size — the same PRE/POST
  difference measured with place fields detached from the cells that own them.
  It is large in three of the five sessions, both circular ones and one linear
  one, so no effect size should be read without it (D8.8). It is measured once
  per session and shuffle, on the full event set, so the conjunction rows carry
  the same residual over a different denominator;
* PRE against POST by Mann-Whitney U with a rank-biserial effect size;
* the same test on PRE subsampled to POST's event count, repeated
  ``sequences.n_subsamples`` times, so the unequal-n objection is answered with
  a distribution rather than an argument;
* significant events per minute of searched state, never counts — the amount of
  non-REM to search differs by up to a factor of two between epochs, and not
  always in the same direction (Cicero_09102014 has more in POST than in PRE).

Sessions are reported separately and **grouped by track geometry**. Linear and
circular tracks put rZ on different scales — a ring maximises the correlation
over 28-29 position offsets where a line maximises over one, so the circular
null is broader and rZ compressed (D8.4). Both circular sessions also leave a
large identity-shuffle residual, but so does one linear session, and that
residual's size orders by ensemble count better than by geometry (D8.8, D8.10).
Magnitudes are therefore pooled only within a geometry. What *is* pooled across
all five is the **direction** of the effect, because a sign has no scale: the
headline is how many sessions show POST > PRE, under each shuffle and each
event set.

Writes results/replay_comparison.csv, results/replay_subsamples.csv and
results/replay_distributions.csv, plus per-session figures.
"""

from __future__ import annotations

import argparse
import json
import logging
from datetime import datetime, timezone

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu

from hc11 import config, plots, results, sequences

log = logging.getLogger("replay")

SESSIONS = list(config.LFP_SESSIONS)
EPOCHS = ("PRE", "MAZE", "POST")
EVENT_SETS = {
    "all": lambda table: table,
    "conjunction": lambda table: table[table.ripple_matched],
}

#: Track geometry per session, which decides how rZ may be compared (D8.4).
GEOMETRY = {
    "Achilles_10252013": "linear",
    "Achilles_11012013": "circular",
    "Cicero_09012014": "linear",
    "Cicero_09102014": "circular",
    "Gatsby_08022013": "linear",
}

#: Below this many scored events in an epoch, a PRE/POST test is underpowered
#: and is reported as such rather than as a negative result.
MIN_SCORED_FOR_TEST = 100


def load_scores(session: str, epochs=EPOCHS) -> pd.DataFrame:
    """Scored events for one session, all epochs, as one table."""
    folder = config.results_dir() / "sequences"
    frames = []
    for epoch in epochs:
        path = folder / f"{session}_{epoch}.csv"
        if not path.exists():
            continue
        table = pd.read_csv(path)
        frames.append(table[table.scored])
    if not frames:
        raise SystemExit(f"{session}: no scored events; run scripts/score_events.py first")
    return pd.concat(frames, ignore_index=True)


def identity_residual(validation: pd.DataFrame, session: str, shuffle: str) -> float:
    """POST minus PRE median rZ with cell identities permuted, or NaN.

    The negative control's own effect size. It belongs beside every real effect
    size rather than in a separate table: a difference of +0.20 means something
    different when the control leaves +0.12 behind than when it leaves -0.03.
    """
    rows = validation[(validation.session == session) & (validation.shuffle == shuffle)
                      & (validation.control == "identity")]
    if rows.empty:
        return float("nan")

    def median(label, epoch):
        hit = rows[(rows.set == label) & (rows.epoch == epoch)]
        return float(hit.median_rz.iloc[0]) if len(hit) else float("nan")

    return median("identity-shuffled", "POST") - median("identity-shuffled", "PRE")


def compare(pre: np.ndarray, post: np.ndarray) -> dict:
    """POST against PRE: medians, Mann-Whitney U, rank-biserial effect size.

    Rank-biserial rather than a bare p-value: with thousands of events a
    trivial shift is significant, and the replication question is whether the
    shift is the size the paper reports, not whether it is non-zero.
    """
    if pre.size < 2 or post.size < 2:
        return {"n_pre": int(pre.size), "n_post": int(post.size),
                "median_pre": np.nan, "median_post": np.nan, "difference": np.nan,
                "rank_biserial": np.nan, "p": np.nan}
    statistic, p = mannwhitneyu(post, pre, alternative="two-sided")
    return {
        "n_pre": int(pre.size), "n_post": int(post.size),
        "median_pre": float(np.median(pre)), "median_post": float(np.median(post)),
        "difference": float(np.median(post) - np.median(pre)),
        # Positive means POST events score higher than PRE events.
        "rank_biserial": float(2 * statistic / (pre.size * post.size) - 1),
        "p": float(p),
    }


def subsample(pre: np.ndarray, post: np.ndarray, n_draws: int, rng) -> dict:
    """``compare`` repeated with PRE cut at random to POST's size."""
    if pre.size <= post.size or post.size < 2:
        return {}
    differences, effects, ps = [], [], []
    for _ in range(n_draws):
        drawn = rng.choice(pre, size=post.size, replace=False)
        result = compare(drawn, post)
        differences.append(result["difference"])
        effects.append(result["rank_biserial"])
        ps.append(result["p"])
    return {
        "n_draws": n_draws, "n_pre_drawn": int(post.size),
        "median_difference": float(np.median(differences)),
        "difference_p5": float(np.percentile(differences, 5)),
        "difference_p95": float(np.percentile(differences, 95)),
        "median_rank_biserial": float(np.median(effects)),
        "pct_p_below_05": 100 * float(np.mean(np.array(ps) < 0.05)),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("sessions", nargs="*", default=SESSIONS)
    parser.add_argument("--no-figures", action="store_true")
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(message)s")

    params = config.load_params()
    cfg = params["sequences"]
    threshold = cfg["significance_z"]
    out_root = config.results_dir()
    state = pd.read_csv(out_root / "mua_summary.csv").query("variant == 'primary'")
    validation_path = out_root / "sequence_validation.csv"
    validation = (pd.read_csv(validation_path) if validation_path.exists()
                  else pd.DataFrame(columns=["session", "shuffle", "control", "set",
                                             "epoch", "median_rz"]))
    stamp = {
        "generated": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "git_commit": config.git_commit(), "sequences": cfg,
    }
    rows, subsample_rows, distribution_rows = [], [], []
    names = args.sessions or SESSIONS
    pooled = {}

    for name in names:
        log.info("\n=== %s", name)
        scores = load_scores(name)
        pooled[name] = scores
        for set_name, select in EVENT_SETS.items():
            table = select(scores)
            for kind in sequences.SHUFFLES:
                column = f"rz_{kind}"
                # Drop unscorable events once, for both columns together: an
                # event whose shuffles were degenerate has a NaN rZ but a
                # perfectly good r, and taking them from different masks
                # silently misaligns the two.
                scored = {
                    epoch: table[(table.epoch == epoch) & table[column].notna()]
                    for epoch in EPOCHS
                }
                by_epoch = {e: f[column].to_numpy() for e, f in scored.items()}
                for epoch, values in by_epoch.items():
                    minutes = state.query(
                        "session == @name and epoch == @epoch"
                    ).state_min
                    minutes = float(minutes.iloc[0]) if len(minutes) else np.nan
                    forward = scored[epoch][f"r_{kind}"].to_numpy() > 0
                    distribution_rows.append({
                        "session": name, "event_set": set_name, "shuffle": kind,
                        "epoch": epoch, "n_scored": int(values.size),
                        "state_min": minutes,
                        "median_rz": float(np.median(values)) if values.size else np.nan,
                        "mean_rz": float(values.mean()) if values.size else np.nan,
                        "p75_rz": float(np.percentile(values, 75)) if values.size else np.nan,
                        "n_significant": int((values > threshold).sum()),
                        "pct_significant": (
                            100 * float((values > threshold).mean()) if values.size
                            else np.nan
                        ),
                        "significant_per_min": (
                            (values > threshold).sum() / minutes
                            if minutes and np.isfinite(minutes) else np.nan
                        ),
                        "pct_forward": 100 * float(forward.mean()) if forward.size else np.nan,
                        "median_rz_forward": (
                            float(np.median(values[forward])) if forward.any() else np.nan
                        ),
                        "median_rz_reverse": (
                            float(np.median(values[~forward])) if (~forward).any()
                            else np.nan
                        ),
                    })
                result = compare(by_epoch["PRE"], by_epoch["POST"])
                residual = identity_residual(validation, name, kind)
                difference = result["difference"]
                rows.append({
                    "session": name, "geometry": GEOMETRY.get(name, "unknown"),
                    "event_set": set_name, "shuffle": kind, **result,
                    "identity_residual": residual,
                    # How much of the difference the negative control reproduces.
                    # Reported as a fraction and left to the reader: there is no
                    # principled cutoff and inventing one turns a graded result
                    # into a binary one.
                    "residual_fraction": (
                        residual / difference
                        if np.isfinite(residual) and difference not in (0, np.nan)
                        and np.isfinite(difference) and difference != 0
                        else np.nan
                    ),
                    "post_gt_pre": bool(difference > 0),
                    "underpowered": bool(
                        min(result["n_pre"], result["n_post"]) < MIN_SCORED_FOR_TEST
                    ),
                })
                rng = np.random.default_rng(cfg["subsample_seed"])
                drawn = subsample(by_epoch["PRE"], by_epoch["POST"],
                                  cfg["n_subsamples"], rng)
                if drawn:
                    subsample_rows.append({"session": name, "event_set": set_name,
                                           "shuffle": kind, **drawn})
                log.info(
                    "  %-11s %-18s PRE %+0.3f (n=%d)  POST %+0.3f (n=%d)  "
                    "diff %+0.3f  rb %+0.3f  p %.2g   control leaves %+0.3f (%s)%s",
                    set_name, kind, result["median_pre"], result["n_pre"],
                    result["median_post"], result["n_post"], result["difference"],
                    result["rank_biserial"], result["p"], residual,
                    f"{100 * residual / difference:.0f}% of it"
                    if np.isfinite(residual) and difference else "n/a",
                    f"   [matched-n {drawn['median_difference']:+0.3f}]" if drawn else "",
                )

        if not args.no_figures:
            folder = config.results_dir(name)
            sidecar = {**stamp, "session": name}
            for set_name, select in EVENT_SETS.items():
                table = select(scores)
                suffix = "" if set_name == "all" else f"_{set_name}"
                plots.plot_sequence_cdf(
                    table, name, set_name, threshold,
                    folder / f"replay_cdf{suffix}.png", sidecar,
                )

    # --- pooled within a geometry; the direction count across all of them -----
    if len(names) > 1:
        for geometry in ("linear", "circular"):
            members = [n for n in names if GEOMETRY.get(n) == geometry]
            if len(members) < 2:
                continue
            log.info("\n=== pooled, %s tracks only (%d sessions: %s)",
                     geometry, len(members), ", ".join(members))
            combined = pd.concat([pooled[n] for n in members], ignore_index=True)
            for set_name, select in EVENT_SETS.items():
                table = select(combined)
                for kind in sequences.SHUFFLES:
                    column = f"rz_{kind}"
                    pre = table[table.epoch == "PRE"][column].dropna().to_numpy()
                    post = table[table.epoch == "POST"][column].dropna().to_numpy()
                    result = compare(pre, post)
                    rows.append({"session": f"pooled_{geometry}", "geometry": geometry,
                                 "event_set": set_name, "shuffle": kind, **result,
                                 "post_gt_pre": bool(result["difference"] > 0),
                                 "underpowered": False})
                    rng = np.random.default_rng(cfg["subsample_seed"])
                    drawn = subsample(pre, post, cfg["n_subsamples"], rng)
                    if drawn:
                        subsample_rows.append({"session": f"pooled_{geometry}",
                                               "event_set": set_name,
                                               "shuffle": kind, **drawn})
                    log.info(
                        "  %-11s %-18s PRE %+0.3f  POST %+0.3f  diff %+0.3f  "
                        "rb %+0.3f  p %.2g",
                        set_name, kind, result["median_pre"], result["median_post"],
                        result["difference"], result["rank_biserial"], result["p"],
                    )

        # The headline: a sign has no scale, so it pools across geometries.
        log.info("\n=== how many sessions show POST > PRE (n=%d)", len(names))
        per_session = [r for r in rows if r["session"] in names]
        for set_name in EVENT_SETS:
            for kind in sequences.SHUFFLES:
                here = [r for r in per_session
                        if r["event_set"] == set_name and r["shuffle"] == kind]
                positive = [r for r in here if r["post_gt_pre"]]
                significant = [r for r in positive if r["p"] < 0.05]
                weak = [r for r in here if r["underpowered"]]
                log.info(
                    "  %-11s %-18s %d/%d positive, %d of those at p<0.05%s",
                    set_name, kind, len(positive), len(here), len(significant),
                    f"; {len(weak)} underpowered" if weak else "",
                )

    results.round_for_csv(pd.DataFrame(rows)).to_csv(
        out_root / "replay_comparison.csv", index=False)
    results.round_for_csv(pd.DataFrame(subsample_rows)).to_csv(
        out_root / "replay_subsamples.csv", index=False)
    results.round_for_csv(pd.DataFrame(distribution_rows)).to_csv(
        out_root / "replay_distributions.csv", index=False)
    (out_root / "replay_comparison.json").write_text(json.dumps(stamp, indent=2) + "\n")
    log.info("\nwrote %s",
             (out_root / "replay_comparison.csv").relative_to(config.repo_root()))


if __name__ == "__main__":
    main()
