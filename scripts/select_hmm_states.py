"""Step 10c: choose the number of states, by Claim A rather than by likelihood.

    uv run python scripts/select_hmm_states.py

Reads ``results/hmm_running_fit.csv`` -- the cross-validated sweep over m, each
session capped at its step 10a limit -- and writes the selection table and
figure. Nothing is fitted here.

**Why Claim A and not held-out likelihood.** The two criteria measure different
things. Held-out likelihood measures how well the model predicts spikes, and it
rises close to monotonically as more states carve the observation space more
finely: there is almost always a little more predictable structure to buy. Claim
A measures whether the states correspond to *positions*, which saturates at
roughly one state per position bin and then degrades, because a single position
starts to fragment across several states and the state-to-position map stops
being a map. The step 10a cap arithmetic already landed at about one state per
position bin, and step 10b's Claim A peak at m = 12 against 29 position bins is
consistent with that. This phase asks whether a receptive-field-free method
recovers positional structure, so the phase's own criterion selects the model.

The likelihood curve is reported beside it at every m, and where the two
criteria disagree the table says so. Where they agree, that is evidence the
choice did not matter.

**Medians, never the best.** Every selection statistic is a median across the
five restarts and the three folds. EM has local optima here and the restart
spread within one (m, fold) cell reaches 0.5; taking the best restart would be
choosing the model by its luckiest initialisation. Each fold's own preferred m
is reported too, and a session whose folds disagree is called unstable rather
than given a number.
"""

from __future__ import annotations

import json
import logging
from datetime import datetime, timezone

import numpy as np
import pandas as pd

from fit_hmm import SESSION_CAPS  # noqa: I001  (sibling script)
from hc11 import config, plots, results

log = logging.getLogger("select")

#: How far apart two folds' preferred m may sit, in grid steps, and still agree.
STABLE_WITHIN_STEPS = 1


def summarise(frame: pd.DataFrame, grid: list[int]) -> pd.DataFrame:
    """One row per m: the pooled medians, the per-fold medians, the control."""
    real = frame[frame.kind == "real"]
    null = frame[frame.kind == "shuffled"]
    rows = []
    for m in grid:
        here = real[real.n_states == m]
        there = null[null.n_states == m]
        row = {"n_states": m,
               "claim_a_median": here.claim_a_r.median(),
               "claim_a_min": here.claim_a_r.min(),
               "claim_a_max": here.claim_a_r.max(),
               "claim_a_iqr": here.claim_a_r.quantile(0.75)
               - here.claim_a_r.quantile(0.25),
               "shuffled_median": there.claim_a_r.median(),
               "shuffled_p95": there.claim_a_r.quantile(0.95),
               "held_out_ll_median": here.held_out_log_likelihood.median(),
               "median_error_cm": here.median_error_cm.median(),
               "chance_error_cm": here.chance_error_cm.median(),
               "shuffled_error_cm": there.median_error_cm.median(),
               "dwell_s": here.median_dwell_s.median(),
               "shuffled_dwell_s": there.median_dwell_s.median(),
               "claim_b_position": here.claim_b_position.median(),
               "claim_b_spectral": here.claim_b_spectral.median(),
               "converged": here.converged.mean(),
               "n_fits": len(here)}
        row["margin_over_shuffled"] = row["claim_a_median"] - row["shuffled_median"]
        row["claim_b_gap"] = row["claim_b_position"] - row["claim_b_spectral"]
        for fold, part in here.groupby("fold"):
            row[f"claim_a_fold{int(fold)}"] = part.claim_a_r.median()
        rows.append(row)
    return pd.DataFrame(rows)


def choose(table: pd.DataFrame, frame: pd.DataFrame) -> dict:
    """Claim A picks m; the likelihood's pick is recorded beside it, not used."""
    grid = list(table.n_states)
    by_claim_a = int(table.loc[table.claim_a_median.idxmax(), "n_states"])
    by_likelihood = int(table.loc[table.held_out_ll_median.idxmax(), "n_states"])
    fold_columns = [c for c in table.columns if c.startswith("claim_a_fold")]
    fold_choice = {c: int(table.loc[table[c].idxmax(), "n_states"])
                   for c in sorted(fold_columns)}
    indices = [grid.index(v) for v in fold_choice.values()]
    unstable = (max(indices) - min(indices)) > STABLE_WITHIN_STEPS
    chosen = table[table.n_states == by_claim_a].iloc[0]
    monotone, run = results.leading_monotone_run(list(table.held_out_ll_median))
    prominence = float(table.claim_a_median.max() - table.claim_a_median.min())
    # Does the chosen m sit within each fold's own restart noise of that fold's
    # best m? A looser, and more honest, reading of stability than the argmax:
    # the argmax of a flat curve moves on nothing.
    within_noise = all(
        float(table.loc[table.n_states == by_claim_a, column].iloc[0])
        >= float(table[column].max()) - float(chosen.claim_a_iqr)
        for column in fold_columns)
    real = frame[(frame.kind == "real") & (frame.n_states == by_claim_a)].claim_a_r
    null = frame[(frame.kind == "shuffled")
                 & (frame.n_states == by_claim_a)].claim_a_r
    exceeded = float(np.mean(null.to_numpy()[:, None] >= real.to_numpy()[None, :]))
    return {"m_by_claim_a": by_claim_a, "m_by_likelihood": by_likelihood,
            "criteria_agree": by_claim_a == by_likelihood,
            # An agreement is only evidence if the likelihood actually turned
            # over inside the grid. Where it is still rising at the cap, the two
            # criteria coincide because the grid ran out, not because they agree.
            "likelihood_rising_at_cap": by_likelihood == max(grid),
            "likelihood_monotone": monotone, "likelihood_monotone_run": run,
            "fold_choices": "/".join(str(fold_choice[c]) for c in sorted(fold_choice)),
            "folds_unstable": bool(unstable),
            "choice_within_fold_noise": bool(within_noise),
            "claim_a_at_choice": float(chosen.claim_a_median),
            "claim_a_curve_range": prominence,
            "restart_iqr_at_choice": float(chosen.claim_a_iqr),
            "curve_flatter_than_noise": prominence < float(chosen.claim_a_iqr),
            "shuffled_at_choice": float(chosen.shuffled_median),
            "shuffled_p95_at_choice": float(chosen.shuffled_p95),
            "margin_at_choice": float(chosen.margin_over_shuffled),
            "shuffled_exceeds_real_frac": exceeded,
            "median_error_cm": float(chosen.median_error_cm),
            "shuffled_error_cm": float(chosen.shuffled_error_cm),
            "chance_error_cm": float(chosen.chance_error_cm),
            "dwell_s": float(chosen.dwell_s),
            "shuffled_dwell_s": float(chosen.shuffled_dwell_s),
            "claim_b_position": float(chosen.claim_b_position),
            "claim_b_spectral": float(chosen.claim_b_spectral),
            "claim_b_gap": float(chosen.claim_b_gap)}


def main() -> None:
    logging.basicConfig(level=logging.INFO, format="%(message)s")
    frame = pd.read_csv(config.results_dir() / "hmm_running_fit.csv")
    phase = pd.read_csv(config.results_dir() / "hmm_bout_phase_baseline.csv")
    tables, choices = {}, []
    for name, part in frame.groupby("session"):
        grid = sorted(part[part.kind == "real"].n_states.unique())
        table = summarise(part, grid)
        table.insert(0, "session", name)
        tables[name] = table
        decision = {"session": name, **choose(table, part)}
        decision["bout_phase_r"] = float(phase[phase.session == name].phase_r.median())
        decision["crossing_time_s"] = float(part.crossing_time_s.median())
        decision["largest_m_swept"] = int(max(grid))
        decision["state_cap_10a"] = SESSION_CAPS[name]
        choices.append(decision)
    per_m = pd.concat(tables.values(), ignore_index=True)
    results.round_for_csv(per_m).to_csv(
        config.results_dir() / "hmm_state_selection.csv", index=False)
    decisions = pd.DataFrame(choices)
    results.round_for_csv(decisions).to_csv(
        config.results_dir() / "hmm_state_choice.csv", index=False)

    for _, row in decisions.iterrows():
        agreement = ("agree" if row.criteria_agree else "DIVERGE")
        if row.criteria_agree and row.likelihood_rising_at_cap:
            agreement = "agree, but only because the grid ran out"
        log.info("%-20s m=%2d by Claim A (likelihood would choose %2d -- %s)",
                 row.session, row.m_by_claim_a, row.m_by_likelihood, agreement)
        log.info("%-20s folds choose %s (%s); curve spans %.2f against a restart "
                 "IQR of %.2f%s", "", row.fold_choices,
                 "unstable" if row.folds_unstable else "stable",
                 row.claim_a_curve_range, row.restart_iqr_at_choice,
                 " -- flatter than the noise, so m barely matters"
                 if row.curve_flatter_than_noise else "")
        log.info("%-20s Claim A %.2f vs shuffled fit %.2f (p95 %.2f, margin %+.2f); "
                 "%.0f%% of shuffled fits beat a real one", "",
                 row.claim_a_at_choice, row.shuffled_at_choice,
                 row.shuffled_p95_at_choice, row.margin_at_choice,
                 100 * row.shuffled_exceeds_real_frac)

    plots.plot_hmm_state_selection(per_m, decisions,
                                   config.results_dir() / "hmm_state_selection.png")
    (config.results_dir() / "hmm_state_selection.json").write_text(json.dumps(
        {"script": "scripts/select_hmm_states.py", "git_commit": config.git_commit(),
         "created": datetime.now(timezone.utc).isoformat(timespec="seconds"),
         "criterion": "median Claim A across restarts and folds",
         "reported_alongside": "median held-out log likelihood",
         "stable_within_steps": STABLE_WITHIN_STEPS}, indent=2) + "\n")
    log.info("wrote results/hmm_state_selection.csv and results/hmm_state_choice.csv")


if __name__ == "__main__":
    main()
