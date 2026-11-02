"""Step 10e: the two methods on one axis, against the number of cells available.

    uv run python scripts/method_comparison.py

Step 10d settled the question this step was built to ask, so this does not
re-derive it: **the receptive-field-free method does not clear the Phase 2
sensitivity floor.** The two sessions where it works (77 and 73 place cells) are
the two where place fields already worked; of the three below the floor, two are
uninformative on their own running Claim B (D10.21) and the third is a negative.
The method sees 1.3-3.4x more cells per session than Phase 1 did, and the extra
cells did not convert into resolution on any session that needed them.

What this step adds is the picture that makes that one statement checkable: every
sleep effect, every method's own control, and both methods' running sensitivity,
plotted against the ensemble each method actually had. The relationship between
available cells and whether either method resolves anything is the result.

Three readings the figure is built to support, each of which is also a column in
`results/method_comparison.csv`:

* A session moves **right** between the two methods -- more cells -- and does not
  move **up**. The connecting line is the point.
* Each effect is drawn against its own control: Phase 1's identity shuffle
  (fields detached from cells) and Phase 4's shuffled fit (the fit destroyed,
  each cell's emission rates kept). They destroy different things, so neither is
  a check on the other and the larger is the one to believe (D10.19).
* Both methods' running sensitivity -- where the structure is certainly present
  -- orders the five sessions identically, and by ensemble size.
"""

from __future__ import annotations

import json
import logging
from datetime import datetime, timezone

import numpy as np
import pandas as pd

from hc11 import config, plots, results
from replay_comparison_hmm import VACUOUS, compare  # noqa: I001  (sibling script)

log = logging.getLogger("methods")

SHUFFLE = "column_cycle"


def phase1_effects(session: str) -> dict:
    """Phase 1's POST-over-PRE effect and its identity-shuffle control.

    Recomputed here from the per-event scores rather than read from
    `replay_comparison.csv`, because that table carries the control as a
    difference in medians and this step needs the control's rank-biserial on
    the same footing as the real effect.
    """
    folder = config.results_dir() / "sequences"
    out = {}
    for label, suffix in (("phase1", ""), ("phase1_control", "_identity_shuffle")):
        by_epoch = {}
        for epoch in ("PRE", "POST"):
            path = folder / f"{session}_{epoch}{suffix}.csv"
            if not path.exists():
                by_epoch[epoch] = np.empty(0)
                continue
            table = pd.read_csv(path)
            table = table[table.scored]
            by_epoch[epoch] = table[f"rz_{SHUFFLE}"].dropna().to_numpy()
        result = compare(by_epoch["PRE"], by_epoch["POST"])
        out[f"{label}_rank_biserial"] = result["rank_biserial"]
        out[f"{label}_difference"] = result["difference"]
        out[f"{label}_p"] = result["p"]
        out[f"{label}_n_pre"] = result["n_pre"]
        out[f"{label}_n_post"] = result["n_post"]
    return out


def main() -> None:
    logging.basicConfig(level=logging.INFO, format="%(message)s")
    fields = pd.read_csv(config.results_dir() / "place_fields_summary.csv")
    floor = pd.read_csv(config.results_dir() / "sensitivity_floor.csv")
    floor = floor[floor.shuffle == SHUFFLE].set_index("session")
    hmm_sleep = pd.read_csv(config.results_dir() / "hmm_replay_comparison.csv")
    hmm_running = pd.read_csv(config.results_dir() / "hmm_state_selection.csv")
    choice = pd.read_csv(config.results_dir() / "hmm_state_choice.csv").set_index(
        "session")

    rows = []
    for name in sorted(hmm_sleep.session.unique()):
        spectral = hmm_sleep[(hmm_sleep.session == name)
                             & (hmm_sleep.ordering == "spectral")].iloc[0]
        position = hmm_sleep[(hmm_sleep.session == name)
                             & (hmm_sleep.ordering == "position")].iloc[0]
        chosen = choice.loc[name]
        running = hmm_running[(hmm_running.session == name)
                              & (hmm_running.n_states == chosen.m_by_claim_a)].iloc[0]
        place_cells = int(fields[(fields.session == name)
                                 & fields.is_place_cell].cluster_id.nunique())
        row = {"session": name,
               "n_place_cells": place_cells,
               "phase4_rank_biserial": spectral.real_rank_biserial,
               "phase4_control_rank_biserial": spectral.shuffled_rank_biserial,
               "phase4_p": spectral.real_p,
               "phase4_position_rank_biserial": position.real_rank_biserial,
               "phase4_position_control": position.shuffled_rank_biserial,
               "running_rz_phase1": float(floor.loc[name, "median_rz"])
               if name in floor.index else np.nan,
               # The behaviour-free Claim B is the measure matched to the
               # behaviour-free sleep effect; the position one is what D10.14
               # and D10.21 used to mark the uninformative pair. Both are kept
               # because they do not order the sessions the same way.
               "running_claim_b_phase4": running.claim_b_spectral,
               "running_claim_b_position": running.claim_b_position,
               "claim_ab_gap_post": spectral.claim_ab_gap_post,
               "informative": name not in VACUOUS}
        row.update(phase1_effects(name))
        rows.append(row)
    table = pd.DataFrame(rows)

    # Two different counts, both worth having. ``n_pyramidal`` is the ensemble
    # Phase 4 decoded with -- every pyramidal cell in the session, which is what
    # the x-axis of the figure means. ``max_active_per_event`` is how many of
    # them ever fired inside one event, which is far smaller and is the reason
    # a larger ensemble does not buy as much as it looks like it should.
    overview = pd.read_csv(config.results_dir() / "session_overview.csv")
    table["n_pyramidal"] = table.session.map(
        overview.set_index("session").n_pyr.to_dict())
    scores = pd.read_csv(config.results_dir() / "hmm_event_scores.csv")
    table["max_active_per_event"] = table.session.map(
        scores.groupby("session").n_active_pyr.max().to_dict())
    table["median_active_per_event"] = table.session.map(
        scores.groupby("session").n_active_pyr.median().to_dict())
    table["cells_ratio"] = table.n_pyramidal / table.n_place_cells
    table["phase1_margin"] = (table.phase1_rank_biserial
                              - table.phase1_control_rank_biserial)
    table["phase4_margin"] = (table.phase4_rank_biserial
                              - table.phase4_control_rank_biserial)
    table["phase4_margin_ratio"] = (table.phase4_rank_biserial
                                    / table.phase4_control_rank_biserial)
    table = table.sort_values("n_place_cells", ascending=False)
    results.round_for_csv(table).to_csv(
        config.results_dir() / "method_comparison.csv", index=False)

    for _, row in table.iterrows():
        log.info("%-20s %3d place cells -> %3d pyramidal (%.1fx)%s",
                 row.session, row.n_place_cells, row.n_pyramidal, row.cells_ratio,
                 "" if row.informative else "   [uninformative]")
        log.info("%-20s   Phase 1 rb %+.3f (control %+.3f, margin %+.3f) | "
                 "Phase 4 rb %+.3f (control %+.3f, margin %+.3f)", "",
                 row.phase1_rank_biserial, row.phase1_control_rank_biserial,
                 row.phase1_margin, row.phase4_rank_biserial,
                 row.phase4_control_rank_biserial, row.phase4_margin)
        log.info("%-20s   running: Phase 2 window rz %5.2f | Phase 4 Claim B "
                 "%.2f behaviour-free, %.2f by position", "", row.running_rz_phase1,
                 row.running_claim_b_phase4, row.running_claim_b_position)

    plots.plot_method_comparison(
        table, config.results_dir() / "method_comparison.png")
    (config.results_dir() / "method_comparison.json").write_text(json.dumps(
        {"script": "scripts/method_comparison.py", "git_commit": config.git_commit(),
         "created": datetime.now(timezone.utc).isoformat(timespec="seconds"),
         "shuffle": SHUFFLE, "uninformative_sessions": sorted(VACUOUS),
         "phase4_effect": "behaviour-free (spectral) state ordering"},
        indent=2) + "\n")
    log.info("wrote results/method_comparison.csv")


if __name__ == "__main__":
    main()
