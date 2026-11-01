"""Step 10c follow-up: does population firing rate alone predict position?

    uv run python scripts/population_rate_baseline.py

The fit-on-shuffled-data control of step 10c is not a null on four of five
sessions (D10.12): a model fitted to spike trains with the cross-cell alignment
destroyed still decodes held-out position at r = 0.34-0.51 on the linear-track
sessions. The bout-timing explanation was tested and rejected
(``bout_phase_baseline.py``). This script tests the remaining cheap candidate.

**The hypothesis.** The shuffle destroys *which* cell fires. It does not stop a
fitted model from having states that differ in *total* population rate, and the
posterior is then evaluated on real spikes with a map built from real training
behaviour. So if population rate varies systematically with position -- near
reward sites, at turnarounds, wherever the animal slows -- the whole decode is
available from rate alone, with no cell identity involved at any point. That
would also predict the observed geometry split, since reward structure differs
between a linear track with two ends and a ring with one reward site.

**The test, with no model in it.** Two measurements per session:

1. The plain association between binned population rate and position: Pearson
   on a linear track, the circular-linear coefficient on a ring, where a
   Pearson correlation would depend on where the ring was cut.
2. A rate-only decoder scored exactly like Claim A: quantile-bin the rate into
   as many levels as that session's chosen *m*, map each level to the mean
   training position, decode held-out bins, and take the same statistic and the
   same permuted-map chance level. This is the like-for-like comparison, because
   Claim A's number is a decode and a raw correlation is not.

If the rate-only decoder reaches roughly the size of the leak, that is the
mechanism. If it does not, the last cheap explanation is eliminated and
"undetermined" becomes a stronger statement than it was.
"""

from __future__ import annotations

import json
import logging
from datetime import datetime, timezone

import numpy as np
import pandas as pd

from fit_hmm import N_FOLDS, N_MAP_PERMUTATIONS, folds, running_sequences  # noqa: I001
from hc11 import config, hmm
from place_fields import prepare

log = logging.getLogger("rate")


def rate_levels(train_rate, test_rate, n_levels):
    """Quantile bins of population rate, cut on the training bins only."""
    edges = np.quantile(train_rate, np.linspace(0, 1, n_levels + 1)[1:-1])
    return (np.searchsorted(edges, train_rate, side="right"),
            np.searchsorted(edges, test_rate, side="right"))


def decode(train_level, train_position, test_level, test_position, n_levels,
           track, rng):
    """Score a level-to-position map the way Claim A scores a state-to-position map."""
    mapped = hmm.state_positions(np.eye(n_levels)[train_level], train_position,
                                 n_levels, circular=track.is_circular,
                                 extent=track.extent)
    decoded = mapped[test_level]
    keep = np.isfinite(decoded) & np.isfinite(test_position)
    decoded, truth, level = decoded[keep], test_position[keep], test_level[keep]

    def error(values):
        if track.is_circular:
            return float(np.median(hmm.circular_distance(values, truth, track.extent)))
        return float(np.median(np.abs(values - truth)))

    correlation = (hmm.circular_correlation(decoded, truth, track.extent)
                   if track.is_circular
                   else float(np.corrcoef(decoded, truth)[0, 1]))
    chance = [error(mapped[rng.permutation(n_levels)][level])
              for _ in range(N_MAP_PERMUTATIONS)]
    return correlation, 100 * error(decoded), 100 * float(np.median(chance))


def main() -> None:
    logging.basicConfig(level=logging.INFO, format="%(message)s")
    params = config.load_params()
    bin_s = params["decoding"]["time_bin_s"]
    chosen = pd.read_csv(config.results_dir() / "hmm_state_choice.csv"
                         ).set_index("session")
    rows = []
    for name in config.LFP_SESSIONS:
        session, track, running, _ = prepare(name, params)
        pyr = np.asarray(session.units("pyr"))
        counts, positions = running_sequences(session, running, pyr, bin_s)
        n_levels = int(chosen.loc[name, "m_by_claim_a"])
        pooled_rate = np.concatenate([c.sum(axis=1) for c in counts]) / bin_s
        pooled_position = np.concatenate(positions)
        association = (hmm.circular_linear_correlation(pooled_position, pooled_rate,
                                                       track.extent)
                       if track.is_circular
                       else abs(float(np.corrcoef(pooled_rate, pooled_position)[0, 1])))
        for fold_index, (train, test) in enumerate(folds(len(counts), N_FOLDS)):
            train_rate = np.concatenate([counts[i].sum(axis=1) for i in train]) / bin_s
            test_rate = np.concatenate([counts[i].sum(axis=1) for i in test]) / bin_s
            train_level, test_level = rate_levels(train_rate, test_rate, n_levels)
            correlation, error, chance = decode(
                train_level, np.concatenate([positions[i] for i in train]),
                test_level, np.concatenate([positions[i] for i in test]),
                n_levels, track, np.random.default_rng(7000 + fold_index))
            rows.append({"session": name, "fold": fold_index, "n_levels": n_levels,
                         "rate_position_association": association,
                         "rate_only_r": correlation, "median_error_cm": error,
                         "chance_error_cm": chance,
                         "circular": bool(track.is_circular)})
    frame = pd.DataFrame(rows)
    frame.round(5).to_csv(config.results_dir() / "hmm_rate_baseline.csv", index=False)

    choice = chosen.reset_index()
    for name, part in frame.groupby("session"):
        row = choice[choice.session == name].iloc[0]
        log.info("%-20s rate~position %.2f | rate-only decode r %s (median %+.2f), "
                 "error %.0f cm vs %.0f chance", name,
                 part.rate_position_association.iloc[0],
                 "  ".join(f"{v:+.2f}" for v in part.rate_only_r),
                 float(part.rate_only_r.median()),
                 float(part.median_error_cm.median()),
                 float(part.chance_error_cm.median()))
        log.info("%-20s   for comparison: shuffled-fit Claim A %+.2f, real %+.2f",
                 "", row.shuffled_at_choice, row.claim_a_at_choice)
    (config.results_dir() / "hmm_rate_baseline.json").write_text(json.dumps(
        {"script": "scripts/population_rate_baseline.py",
         "git_commit": config.git_commit(),
         "created": datetime.now(timezone.utc).isoformat(timespec="seconds"),
         "folds": N_FOLDS, "bin_s": bin_s,
         "levels": "the session's chosen m"}, indent=2) + "\n")
    log.info("wrote results/hmm_rate_baseline.csv")


if __name__ == "__main__":
    main()
