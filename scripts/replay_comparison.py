"""Step 7c: the Fig. 1C comparison — rZ distributions for PRE, MAZE and POST.

    uv run python scripts/replay_comparison.py
    uv run python scripts/replay_comparison.py --no-figures

Reads results/sequences/<session>_<epoch>.csv and reports, for each session,
each shuffle kind and each candidate-event set:

* the cumulative distribution of rZ with medians marked;
* PRE against POST by Mann-Whitney U with a rank-biserial effect size;
* the same test on PRE subsampled to POST's event count, repeated
  ``sequences.n_subsamples`` times, so the unequal-n objection is answered with
  a distribution rather than an argument;
* significant events per minute of searched state, never counts — PRE holds
  roughly twice the non-REM of POST in both sessions.

Sessions are reported separately. They are different mazes with different
decodable fractions, and two of them is an internal replication rather than a
sample to pool; the pooled result is reported second and marked as secondary.

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

from hc11 import config, plots, sequences

log = logging.getLogger("replay")

SESSIONS = ["Achilles_10252013", "Achilles_11012013"]
EPOCHS = ("PRE", "MAZE", "POST")
EVENT_SETS = {
    "all": lambda table: table,
    "conjunction": lambda table: table[table.ripple_matched],
}


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
                by_epoch = {
                    epoch: table[table.epoch == epoch][column].dropna().to_numpy()
                    for epoch in EPOCHS
                }
                for epoch, values in by_epoch.items():
                    minutes = state.query(
                        "session == @name and epoch == @epoch"
                    ).state_min
                    minutes = float(minutes.iloc[0]) if len(minutes) else np.nan
                    forward = table[(table.epoch == epoch)][f"r_{kind}"].to_numpy() > 0
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
                rows.append({"session": name, "event_set": set_name, "shuffle": kind,
                             **result})
                rng = np.random.default_rng(cfg["subsample_seed"])
                drawn = subsample(by_epoch["PRE"], by_epoch["POST"],
                                  cfg["n_subsamples"], rng)
                if drawn:
                    subsample_rows.append({"session": name, "event_set": set_name,
                                           "shuffle": kind, **drawn})
                log.info(
                    "  %-11s %-18s PRE %+0.3f (n=%d)  POST %+0.3f (n=%d)  "
                    "diff %+0.3f  rb %+0.3f  p %.2g%s",
                    set_name, kind, result["median_pre"], result["n_pre"],
                    result["median_post"], result["n_post"], result["difference"],
                    result["rank_biserial"], result["p"],
                    f"   [matched-n diff {drawn['median_difference']:+0.3f}, "
                    f"{drawn['pct_p_below_05']:.0f}% of draws p<0.05]" if drawn else "",
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

    # --- pooled, secondary ---------------------------------------------------
    if len(names) > 1:
        log.info("\n=== pooled (secondary)")
        combined = pd.concat(pooled.values(), ignore_index=True)
        for set_name, select in EVENT_SETS.items():
            table = select(combined)
            for kind in sequences.SHUFFLES:
                column = f"rz_{kind}"
                pre = table[table.epoch == "PRE"][column].dropna().to_numpy()
                post = table[table.epoch == "POST"][column].dropna().to_numpy()
                result = compare(pre, post)
                rows.append({"session": "pooled", "event_set": set_name,
                             "shuffle": kind, **result})
                rng = np.random.default_rng(cfg["subsample_seed"])
                drawn = subsample(pre, post, cfg["n_subsamples"], rng)
                if drawn:
                    subsample_rows.append({"session": "pooled", "event_set": set_name,
                                           "shuffle": kind, **drawn})
                log.info(
                    "  %-11s %-18s PRE %+0.3f  POST %+0.3f  diff %+0.3f  rb %+0.3f  p %.2g",
                    set_name, kind, result["median_pre"], result["median_post"],
                    result["difference"], result["rank_biserial"], result["p"],
                )

    pd.DataFrame(rows).round(5).to_csv(out_root / "replay_comparison.csv", index=False)
    pd.DataFrame(subsample_rows).round(5).to_csv(
        out_root / "replay_subsamples.csv", index=False)
    pd.DataFrame(distribution_rows).round(5).to_csv(
        out_root / "replay_distributions.csv", index=False)
    (out_root / "replay_comparison.json").write_text(json.dumps(stamp, indent=2) + "\n")
    log.info("\nwrote %s",
             (out_root / "replay_comparison.csv").relative_to(config.repo_root()))


if __name__ == "__main__":
    main()
