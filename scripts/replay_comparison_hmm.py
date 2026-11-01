"""Step 10d: PRE against POST, scored by the receptive-field-free model.

    uv run python scripts/replay_comparison_hmm.py

Reads `results/hmm_event_scores.csv` (step 10d's scoring) and asks the Phase 1
question of it: do POST sleep events score higher than PRE? The statistics are
Phase 1's -- medians, Mann-Whitney U, rank-biserial, and PRE subsampled to
POST's event count -- so the two phases' answers are comparable.

Three things are different, and all three come from this phase's own controls:

* **Every comparison is reported as a pair**, the real fit beside the model
  fitted to shuffled training data (D10.12). The shuffled model is scored on the
  *same events* with the *same shuffles*, so the difference between the two
  PRE/POST effects is the quantity that means something. On Cicero_09012014,
  Cicero_09102014 and Gatsby_08022013 the shuffled fit overlapped the real one
  on running by 16-30 %, so a small margin on those sessions means nothing.
* **Two sessions are marked vacuous in advance** (D10.14), not after the fact.
  Claim B on *running* -- where the structure is certainly there -- is 0.08 on
  Gatsby_08022013 and 0.06 on Cicero_09012014, the two smallest ensembles. A
  null on those sessions carries no information about sleep. They are reported
  as uninformative, as Phase 3 reported its own unpowered cells, rather than
  counted as negatives.
* **Both state orderings are reported everywhere**, and their gap with them.
  Only the model's own ordering is receptive-field-free; the position ordering
  is the comparison point and uses behaviour.
"""

from __future__ import annotations

import json
import logging
from datetime import datetime, timezone

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu

from hc11 import config, plots, results

log = logging.getLogger("hmm-replay")

#: Sessions whose running Claim B is too small for a sleep null to mean anything
#: (D10.14), with the running value that makes them vacuous.
VACUOUS = {"Gatsby_08022013": 0.08, "Cicero_09012014": 0.06}


def compare(pre: np.ndarray, post: np.ndarray) -> dict:
    """POST against PRE. Phase 1's `replay_comparison.compare`, unchanged."""
    if pre.size < 2 or post.size < 2:
        return {"n_pre": int(pre.size), "n_post": int(post.size),
                "median_pre": np.nan, "median_post": np.nan, "difference": np.nan,
                "rank_biserial": np.nan, "p": np.nan}
    statistic, p = mannwhitneyu(post, pre, alternative="two-sided")
    return {"n_pre": int(pre.size), "n_post": int(post.size),
            "median_pre": float(np.median(pre)), "median_post": float(np.median(post)),
            "difference": float(np.median(post) - np.median(pre)),
            "rank_biserial": float(2 * statistic / (pre.size * post.size) - 1),
            "p": float(p)}


def subsample(pre: np.ndarray, post: np.ndarray, n_draws: int, rng) -> dict:
    """`compare` repeated with PRE cut at random to POST's size."""
    if pre.size <= post.size or post.size < 2:
        return {}
    differences, effects = [], []
    for _ in range(n_draws):
        result = compare(rng.choice(pre, size=post.size, replace=False), post)
        differences.append(result["difference"])
        effects.append(result["rank_biserial"])
    return {"subsample_median_difference": float(np.median(differences)),
            "subsample_difference_p5": float(np.percentile(differences, 5)),
            "subsample_difference_p95": float(np.percentile(differences, 95)),
            "subsample_median_rank_biserial": float(np.median(effects))}


def main() -> None:
    logging.basicConfig(level=logging.INFO, format="%(message)s")
    params = config.load_params()
    cfg = params["sequences"]
    scores = pd.read_csv(config.results_dir() / "hmm_event_scores.csv")
    rows = []
    for (name, ordering), part in scores.groupby(["session", "ordering"]):
        row = {"session": name, "ordering": ordering,
               "vacuous_by_running_claim_b": name in VACUOUS,
               "running_claim_b": VACUOUS.get(name, np.nan)}
        for kind in ("real", "shuffled"):
            here = part[part.kind == kind]
            pre = here[here.epoch == "PRE"].rz_median.dropna().to_numpy()
            post = here[here.epoch == "POST"].rz_median.dropna().to_numpy()
            result = compare(pre, post)
            rng = np.random.default_rng(cfg["subsample_seed"])
            result.update(subsample(pre, post, cfg["n_subsamples"], rng))
            result["pct_significant_pre"] = 100 * float(
                np.mean(pre >= cfg["significance_z"])) if pre.size else np.nan
            result["pct_significant_post"] = 100 * float(
                np.mean(post >= cfg["significance_z"])) if post.size else np.nan
            row.update({f"{kind}_{key}": value for key, value in result.items()})
        row["margin_difference"] = row["real_difference"] - row["shuffled_difference"]
        row["margin_rank_biserial"] = (row["real_rank_biserial"]
                                       - row["shuffled_rank_biserial"])
        rows.append(row)
    table = pd.DataFrame(rows).sort_values(["session", "ordering"])
    # The Claim A / Claim B gap, carried from 10b onward: how much the
    # behaviour-free ordering gives up against the behaviour-derived one.
    gaps = table.pivot(index="session", columns="ordering", values="real_median_post")
    table = table.merge(
        (gaps["position"] - gaps["spectral"]).rename("claim_ab_gap_post"),
        left_on="session", right_index=True)
    # Phase 1's answer on the same sessions and the same event sets, so the two
    # methods are compared in the artifact rather than only in the write-up.
    # Phase 1 decoded place fields from a place-cell ensemble; this phase decodes
    # states from every pyramidal cell, so the absolute rz levels are not
    # comparable and the PRE-to-POST differences are.
    phase1 = pd.read_csv(config.results_dir() / "replay_comparison.csv")
    phase1 = phase1[(phase1.event_set == "all")
                    & (phase1.shuffle == "column_cycle")]
    table = table.merge(
        phase1[["session", "difference", "rank_biserial", "p"]].rename(columns={
            "difference": "phase1_difference",
            "rank_biserial": "phase1_rank_biserial", "p": "phase1_p"}),
        on="session", how="left")
    results.round_for_csv(table).to_csv(
        config.results_dir() / "hmm_replay_comparison.csv", index=False)

    for _, row in table.iterrows():
        mark = (f"  [VACUOUS: running Claim B {row.running_claim_b:.2f}]"
                if row.vacuous_by_running_claim_b else "")
        log.info("%-20s %-9s real  PRE %.2f (n=%d) -> POST %.2f (n=%d)  "
                 "diff %+.3f  rb %+.3f  p %.2g%s",
                 row.session, row.ordering, row.real_median_pre, row.real_n_pre,
                 row.real_median_post, row.real_n_post, row.real_difference,
                 row.real_rank_biserial, row.real_p, mark)
        log.info("%-20s %-9s shuf  PRE %.2f -> POST %.2f  diff %+.3f  rb %+.3f  "
                 "p %.2g   | MARGIN diff %+.3f, rb %+.3f", "", "",
                 row.shuffled_median_pre, row.shuffled_median_post,
                 row.shuffled_difference, row.shuffled_rank_biserial,
                 row.shuffled_p, row.margin_difference, row.margin_rank_biserial)
    plots.plot_hmm_replay(table,
                          config.results_dir() / "hmm_replay_comparison.png")
    (config.results_dir() / "hmm_replay_comparison.json").write_text(json.dumps(
        {"script": "scripts/replay_comparison_hmm.py",
         "git_commit": config.git_commit(),
         "created": datetime.now(timezone.utc).isoformat(timespec="seconds"),
         "n_subsamples": cfg["n_subsamples"], "significance_z": cfg["significance_z"],
         "vacuous_sessions": sorted(VACUOUS)}, indent=2) + "\n")
    log.info("wrote results/hmm_replay_comparison.csv")


if __name__ == "__main__":
    main()
