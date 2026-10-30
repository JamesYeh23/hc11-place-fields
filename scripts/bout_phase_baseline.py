"""Step 10c: how much of a position decode is available from bout timing alone.

    uv run python scripts/bout_phase_baseline.py

The fit-on-shuffled-data control of step 10b scored near zero on
Achilles_11012013 (Claim A r 0.05-0.16 against the real fit's 0.41). On the
other four sessions it does not: on Achilles_10252013 and Cicero_09012014 a
model fitted to shuffled spike trains decodes held-out position at r = 0.4-0.6.
Before any real 10c number is read, that has to be explained or the control is
useless.

The hypothesis this script tests is that the leak is **bout timing**, not
spiking. Running comes in bouts, and a bout is a directed run: the animal is at
one end when it starts and the other when it ends, so the bin's *index within
its bout* predicts position whatever the spikes say. A shuffled fit that has
learned nothing but "the population is autocorrelated on the bout timescale"
can ride that.

So: decode held-out position from the within-bout index alone -- mean training
position at each index, with indices past the cap pooled -- using the same
folds and the same statistic as Claim A, and no spikes anywhere. If this
reproduces the shuffled control, the control's leak is bout timing and Claim A
on those sessions is measuring the track's geometry rather than the model.
"""

from __future__ import annotations

import json
import logging
from datetime import datetime, timezone

import numpy as np
import pandas as pd

from fit_hmm import N_FOLDS, folds, running_sequences  # noqa: I001  (sibling script)
from hc11 import config, hmm
from place_fields import prepare

log = logging.getLogger("phase")

MAX_INDEX = 12  # indices at or past this are pooled; bout medians run 5-9 bins


def phase_map(positions, n_index, circular, extent):
    """Mean training position at each within-bout index, on a ring if circular."""
    index = np.concatenate([np.minimum(np.arange(len(p)), n_index - 1)
                            for p in positions])
    return hmm.state_positions(np.eye(n_index)[index], np.concatenate(positions),
                               n_index, circular=circular, extent=extent)


def main() -> None:
    logging.basicConfig(level=logging.INFO, format="%(message)s")
    params = config.load_params()
    bin_s = params["decoding"]["time_bin_s"]
    rows = []
    for name in config.LFP_SESSIONS:
        session, track, running, _ = prepare(name, params)
        pyr = np.asarray(session.units("pyr"))
        counts, positions = running_sequences(session, running, pyr, bin_s)
        for fold_index, (train, test) in enumerate(folds(len(counts), N_FOLDS)):
            mapped = phase_map([positions[i] for i in train], MAX_INDEX,
                               track.is_circular, track.extent)
            index = np.concatenate([np.minimum(np.arange(len(positions[i])),
                                               MAX_INDEX - 1) for i in test])
            truth = np.concatenate([positions[i] for i in test])
            decoded = mapped[index]
            keep = np.isfinite(decoded) & np.isfinite(truth)
            decoded, truth = decoded[keep], truth[keep]
            if track.is_circular:
                r = hmm.circular_correlation(decoded, truth, track.extent)
                error = float(np.median(hmm.circular_distance(decoded, truth,
                                                              track.extent)))
            else:
                r = float(np.corrcoef(decoded, truth)[0, 1])
                error = float(np.median(np.abs(decoded - truth)))
            rows.append({"session": name, "fold": fold_index, "phase_r": r,
                         "median_error_cm": 100 * error,
                         "n_test_bins": int(decoded.size),
                         "median_bout_bins": int(np.median([len(c) for c in counts]))})
    frame = pd.DataFrame(rows)
    path = config.results_dir() / "hmm_bout_phase_baseline.csv"
    frame.round(5).to_csv(path, index=False)
    for name, part in frame.groupby("session"):
        log.info("%-20s phase-only r per fold: %s  (median %+.2f), error %.0f cm",
                 name, "  ".join(f"{v:+.2f}" for v in part.phase_r),
                 float(part.phase_r.median()), float(part.median_error_cm.median()))
    (config.results_dir() / "hmm_bout_phase_baseline.json").write_text(json.dumps(
        {"script": "scripts/bout_phase_baseline.py", "git_commit": config.git_commit(),
         "created": datetime.now(timezone.utc).isoformat(timespec="seconds"),
         "max_index": MAX_INDEX, "folds": N_FOLDS, "bin_s": bin_s}, indent=2) + "\n")
    log.info("wrote %s", path.relative_to(config.repo_root()))


if __name__ == "__main__":
    main()
