"""Step 11a: how much of Claim A survives when the fit is removed entirely.

    uv run python scripts/unfitted_model_baseline.py

The fit-on-shuffled-data control is **not** a per-bin identity permutation and
not a global relabelling: `fit_hmm.roll_cells` gives each cell its own circular
time shift along the concatenated running timeline (D10.26). So the only thing
the fit ever receives from the real data under that control is **each cell's own
marginal count distribution**; every cross-cell and cell-to-position relation is
destroyed. Measured on the five sessions, the per-bin count multiset survives in
under 0.2 % of bins and the per-bin total is uncorrelated with the real one.

That narrows the candidate mechanism for the leak to one that is not about the
shuffle at all. Claim A computes its state-to-position map from **real**
training behaviour and evaluates the posterior on **real** spikes; only the rate
vectors come from the shuffled fit. So an arbitrary m-way partition of real
population vectors may be position-informative however it was obtained.

This script removes the fit entirely and changes nothing else -- same folds,
same `claim_a`, same statistic, same permuted-map chance level -- using rate
matrices that no optimiser ever saw:

* **exemplar** -- *m* real training population vectors, picked at random, used
  as the rate matrix. Keeps the cross-cell structure of real bins; keeps no
  information about which bins they were.
* **marginal** -- each cell's *m* rates drawn independently from that cell's own
  training counts. This is exactly what the shuffled fit was given: real
  per-cell marginals, no co-occurrence.

The two fail differently and the gap between them is the point: whatever
*exemplar* has that *marginal* does not is cross-cell structure. The result and
what it says about Claim A are D10.27 and D10.28.

A third variant is a robustness check rather than a result. Removing the fit
also removes the fitted stickiness, and a negative could otherwise be blamed on
that; **exemplar_sticky** repeats *exemplar* with a transition matrix whose
self-transition reproduces the fitted model's median dwell time on that session
(D10.13), leaving the rates unfitted.
"""

from __future__ import annotations

import json
import logging
from datetime import datetime, timezone

import numpy as np
import pandas as pd

from fit_hmm import (  # noqa: I001  (sibling script)
    N_FOLDS,
    claim_a,
    folds,
    running_sequences,
)
from hc11 import config, hmm, results
from place_fields import prepare

log = logging.getLogger("unfitted")

DEFAULT_DRAWS = 5
VARIANTS = ("exemplar", "marginal", "exemplar_sticky")


def uniform_transition(n_states: int) -> np.ndarray:
    """No temporal structure at all: the fit is removed, stickiness with it."""
    return np.full((n_states, n_states), 1.0 / n_states)


def sticky_transition(n_states: int, dwell_s: float, bin_s: float) -> np.ndarray:
    """Self-transition reproducing a target dwell time, the rest spread evenly.

    Expected dwell is ``bin / (1 - P_ii)``, so ``P_ii = 1 - bin / dwell``.
    """
    stay = float(np.clip(1.0 - bin_s / max(dwell_s, bin_s), 0.0, 1.0 - 1e-6))
    transition = np.full((n_states, n_states), (1.0 - stay) / max(n_states - 1, 1))
    np.fill_diagonal(transition, stay)
    return transition


def random_rates(variant: str, train_counts, n_states: int, bin_s: float, rng):
    """A rate matrix no optimiser produced, in Hz, floored like a fitted one."""
    pooled = np.vstack(train_counts)
    if variant.startswith("exemplar"):
        chosen = rng.choice(len(pooled), size=n_states,
                            replace=len(pooled) < n_states)
        rates = pooled[chosen] / bin_s
    elif variant == "marginal":
        # Each cell independently: its own counts, resampled. Real marginals,
        # no co-occurrence -- what the shuffled fit was handed.
        rates = np.empty((n_states, pooled.shape[1]))
        for cell in range(pooled.shape[1]):
            rates[:, cell] = rng.choice(pooled[:, cell], size=n_states) / bin_s
    else:
        raise ValueError(variant)
    return np.maximum(rates, hmm.RATE_FLOOR_HZ)


def unfitted_model(variant, train_counts, n_states, bin_s, dwell_s, rng):
    transition = (sticky_transition(n_states, dwell_s, bin_s)
                  if variant.endswith("_sticky") else uniform_transition(n_states))
    return hmm.PoissonHMM(
        initial=np.full(n_states, 1.0 / n_states), transition=transition,
        rates=random_rates(variant, train_counts, n_states, bin_s, rng),
        bin_s=bin_s, log_likelihood=float("nan"), n_iterations=0,
        converged=True, seed=-1)


def main() -> None:
    logging.basicConfig(level=logging.INFO, format="%(message)s")
    params = config.load_params()
    bin_s = params["decoding"]["time_bin_s"]
    choice = pd.read_csv(config.results_dir() / "hmm_state_choice.csv").set_index(
        "session")
    rows = []
    for name in config.LFP_SESSIONS:
        session, track, running, _ = prepare(name, params)
        pyr = np.asarray(session.units("pyr"))
        counts, positions = running_sequences(session, running, pyr, bin_s)
        n_states = int(choice.loc[name, "m_by_claim_a"])
        dwell_s = float(choice.loc[name, "dwell_s"])
        for fold_index, (train, test) in enumerate(folds(len(counts), N_FOLDS)):
            train_counts = [counts[i] for i in train]
            train_positions = [positions[i] for i in train]
            test_counts = [counts[i] for i in test]
            test_positions = [positions[i] for i in test]
            for variant in VARIANTS:
                for draw in range(DEFAULT_DRAWS):
                    model = unfitted_model(
                        variant, train_counts, n_states, bin_s, dwell_s,
                        np.random.default_rng(90000 + 100 * draw + fold_index))
                    measured = claim_a(model, train_counts, train_positions,
                                       test_counts, test_positions, track,
                                       np.random.default_rng(1000 + draw))
                    measured.pop("mapped", None)
                    rows.append({"session": name, "variant": variant,
                                 "n_states": n_states, "fold": fold_index,
                                 "draw": draw, **measured})
    frame = pd.DataFrame(rows)
    results.round_for_csv(frame).to_csv(
        config.results_dir() / "hmm_unfitted_baseline.csv", index=False)

    for name, part in frame.groupby("session"):
        row = choice.loc[name]
        log.info("%-20s m=%d | real fit %+.2f (%.0f cm) | shuffled fit %+.2f (%.0f cm) "
                 "| chance %.0f cm", name, int(row.m_by_claim_a),
                 row.claim_a_at_choice, row.median_error_cm,
                 row.shuffled_at_choice, row.shuffled_error_cm, row.chance_error_cm)
        for variant in VARIANTS:
            here = part[part.variant == variant]
            log.info("%-20s   %-16s %+.2f  (%+.2f to %+.2f)   %.0f cm", "", variant,
                     here.claim_a_r.median(), here.claim_a_r.min(),
                     here.claim_a_r.max(), here.median_error_cm.median())
    (config.results_dir() / "hmm_unfitted_baseline.json").write_text(json.dumps(
        {"script": "scripts/unfitted_model_baseline.py",
         "git_commit": config.git_commit(),
         "created": datetime.now(timezone.utc).isoformat(timespec="seconds"),
         "variants": list(VARIANTS), "draws": DEFAULT_DRAWS, "folds": N_FOLDS,
         "n_states": "each session's chosen m (D10.10)"}, indent=2) + "\n")
    log.info("wrote results/hmm_unfitted_baseline.csv")


if __name__ == "__main__":
    main()
