"""Step 11a: how much position a bin's count multiset carries, with identity gone.

    uv run python scripts/sorted_count_baseline.py

**This is a bound on a control this project does not have, not an explanation of
the one it does.** The fit-on-shuffled-data control is a per-cell circular time
shift (D10.26) and it destroys the per-bin count multiset -- the multiset
survives in under 0.2 % of bins. So nothing measured here can account for that
control's leak; the matched test for it is `unfitted_model_baseline.py`, which explains it
(D10.27). The question this one answers, and its answer, are D10.29.

What this answers is whether a **per-bin identity permutation** would be worth
building as a second, stricter control. Such a shuffle would destroy which cell
fired while leaving each bin's multiset of counts intact: how many cells fired,
how the counts are spread across them, how peaked or flat the population vector
is. The population-rate test collapsed that multiset to its sum, which is a weak
summary -- three cells firing four spikes each and twelve cells firing one each
share a total and look nothing alike.

So: decode held-out position from the **sorted count vector** alone, identity
destroyed by construction, through the same folds and the same statistic as the
rate-only decoder (`population_rate_baseline.decode`), so the two are directly
comparable. Three variants:

* **sorted_raw** -- the sorted count vector, partitioned into as many groups as
  that session's chosen *m* by k-means on the training bins. The matched
  comparison: an *m*-way partition, as the model has.
* **sorted_shape** -- the same, with each vector divided by its own total. This
  removes the population-rate information already tested in D10.16 and leaves
  the sparsity profile on its own.
* **nearest_neighbour** -- each held-out bin takes the position of the training
  bin whose sorted vector is closest. Not matched to anything the model does: it
  uses the training set at full resolution instead of *m* groups. It is **not**
  a strict upper bound, and the results show why -- averaging the positions
  inside a cluster suppresses variance, so the *m*-way partition beats the
  single nearest neighbour on two sessions. The two are reported together and
  the larger is the one to read.

Reading it: if the sorted counts carry little position, a per-bin identity
permutation would be a clean and stricter control than ours and is worth
building. If they carry a lot, that control would leak too and is not worth it.
"""

from __future__ import annotations

import json
import logging
from datetime import datetime, timezone

import numpy as np
import pandas as pd
from scipy.cluster.vq import kmeans2

from fit_hmm import N_FOLDS, N_MAP_PERMUTATIONS, folds, running_sequences  # noqa: I001
from hc11 import config, hmm, results
from place_fields import prepare

log = logging.getLogger("sorted")

DEFAULT_DRAWS = 5
VARIANTS = ("sorted_raw", "sorted_shape", "nearest_neighbour")


def score(mapped, test_label, test_position, track, rng) -> dict:
    """A label-to-position map, scored exactly as Claim A and the rate test score one.

    `tests/test_sorted_counts.py` pins this against
    `population_rate_baseline.decode` so "the same statistic" is a checked
    claim rather than an intention.
    """
    decoded = mapped[test_label]
    keep = np.isfinite(decoded) & np.isfinite(test_position)
    decoded, truth, label = decoded[keep], test_position[keep], test_label[keep]
    if decoded.size < 10:
        return {"r": np.nan, "median_error_cm": np.nan, "chance_error_cm": np.nan}

    def error(values):
        if track.is_circular:
            return float(np.median(hmm.circular_distance(values, truth, track.extent)))
        return float(np.median(np.abs(values - truth)))

    chance = [error(mapped[rng.permutation(len(mapped))][label])
              for _ in range(N_MAP_PERMUTATIONS)]
    r = (hmm.circular_correlation(decoded, truth, track.extent) if track.is_circular
         else float(np.corrcoef(decoded, truth)[0, 1]))
    return {"r": r, "median_error_cm": 100 * error(decoded),
            "chance_error_cm": 100 * float(np.median(chance))}


def sorted_vectors(counts, shape: bool):
    """Each bin's counts in sorted order; identity is gone by construction."""
    vectors = np.sort(np.vstack(counts), axis=1).astype(float)
    if shape:
        total = vectors.sum(axis=1, keepdims=True)
        vectors = np.divide(vectors, total, out=np.zeros_like(vectors),
                            where=total > 0)
    return vectors


def partition(train, test, n_groups, seed):
    """k-means on the training vectors; held-out bins join the nearest centroid."""
    centroids, train_label = kmeans2(train, n_groups, minit="++", seed=seed,
                                     missing="warn")
    test_label = np.argmin(
        ((test[:, None, :] - centroids[None, :, :]) ** 2).sum(axis=2), axis=1)
    return train_label, test_label, len(centroids)


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
        for fold_index, (train, test) in enumerate(folds(len(counts), N_FOLDS)):
            train_position = np.concatenate([positions[i] for i in train])
            test_position = np.concatenate([positions[i] for i in test])
            for variant in VARIANTS:
                shape = variant == "sorted_shape"
                train_vectors = sorted_vectors([counts[i] for i in train], shape)
                test_vectors = sorted_vectors([counts[i] for i in test], shape)
                draws = 1 if variant == "nearest_neighbour" else DEFAULT_DRAWS
                for draw in range(draws):
                    if variant == "nearest_neighbour":
                        # Every training bin is its own label, mapping to its
                        # own position: the tightest bound the data allows.
                        test_label = np.argmin(
                            ((test_vectors[:, None, :]
                              - train_vectors[None, :, :]) ** 2).sum(axis=2), axis=1)
                        mapped = train_position
                    else:
                        train_label, test_label, n_groups = partition(
                            train_vectors, test_vectors, n_states, 11000 + draw)
                        mapped = hmm.state_positions(
                            np.eye(n_groups)[train_label], train_position, n_groups,
                            circular=track.is_circular, extent=track.extent)
                    measured = score(mapped, test_label, test_position, track,
                                     np.random.default_rng(1000 + draw))
                    rows.append({"session": name, "variant": variant,
                                 "n_groups": n_states if not shape else n_states,
                                 "fold": fold_index, "draw": draw, **measured})
    frame = pd.DataFrame(rows)
    results.round_for_csv(frame).to_csv(
        config.results_dir() / "hmm_sorted_count_baseline.csv", index=False)

    for name, part in frame.groupby("session"):
        row = choice.loc[name]
        log.info("%-20s real fit %+.2f (%.0f cm) | shuffled fit %+.2f (%.0f cm) | "
                 "chance %.0f cm", name, row.claim_a_at_choice, row.median_error_cm,
                 row.shuffled_at_choice, row.shuffled_error_cm, row.chance_error_cm)
        for variant in VARIANTS:
            here = part[part.variant == variant]
            log.info("%-20s   %-18s %+.2f  (%+.2f to %+.2f)   %.0f cm vs %.0f chance",
                     "", variant, here.r.median(), here.r.min(), here.r.max(),
                     here.median_error_cm.median(), here.chance_error_cm.median())
    (config.results_dir() / "hmm_sorted_count_baseline.json").write_text(json.dumps(
        {"script": "scripts/sorted_count_baseline.py",
         "git_commit": config.git_commit(),
         "created": datetime.now(timezone.utc).isoformat(timespec="seconds"),
         "variants": list(VARIANTS), "draws": DEFAULT_DRAWS, "folds": N_FOLDS,
         "answers": "whether a per-bin identity permutation would be a clean control",
         "does_not_explain": "the per-cell time-shift control actually in use"},
        indent=2) + "\n")
    log.info("wrote results/hmm_sorted_count_baseline.csv")


if __name__ == "__main__":
    main()
