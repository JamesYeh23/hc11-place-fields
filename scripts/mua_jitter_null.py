"""Is the population synchrony real? A spike-jitter null.

Each cell's spikes are displaced independently by U(-100, +100) ms, which keeps
every cell's rate and the slow population envelope but destroys co-firing on the
timescale of an event. Two comparisons are reported per epoch:

**Re-detected** — the whole detector is re-run on jittered spikes. If the
observed events were only the slow envelope crossing a threshold, the null would
produce events of the same size at the same rate.

**Same windows** — participation is recounted inside the *observed* event
windows using jittered spikes. This is the per-event version: it asks, for each
real event, how many cells would have been active there by chance.

    uv run python scripts/mua_jitter_null.py

Writes results/mua_jitter_null.csv and results/mua_jitter_per_event.csv.
"""

from __future__ import annotations

import json
from datetime import datetime, timezone

import numpy as np
import pandas as pd

from hc11 import config
from hc11.events import participation
from hc11.io import load_session
from hc11.mua import detect_mua_events, jitter_spikes
from hc11.preprocess import truncate_to_session

SESSIONS = ["Achilles_10252013", "Achilles_11012013"]
EPOCHS = ("PRE", "POST")  # MAZE has almost no searchable time (D6.3)


def main() -> None:
    params = config.load_params()
    cfg = params["events"]["mua"]
    window_s = cfg["jitter_window_ms"] / 1000.0
    n_repeats = cfg["n_jitter_repeats"]
    rows, per_event_rows = [], []

    for name in SESSIONS:
        session = load_session(name)
        if params["preprocess"]["truncate_to_session"]:
            session, _ = truncate_to_session(session, verbose=False)
        pyr = session.pyr_ids

        for epoch in EPOCHS:
            observed, _ = detect_mua_events(session, epoch, params)
            obs_part = participation(session, observed, pyr)

            redetected_counts, redetected_active, same_window_active = [], [], []
            for repeat in range(n_repeats):
                rng = np.random.default_rng(1000 + repeat)
                shuffled = jitter_spikes(session, rng, window_s=window_s, cluster_ids=pyr)
                null_events, _ = detect_mua_events(shuffled, epoch, params)
                null_part = participation(shuffled, null_events, pyr)
                redetected_counts.append(len(null_events))
                redetected_active.append(null_part.n_active)
                same_window_active.append(participation(shuffled, observed, pyr).n_active)
                print(f"  {name} {epoch} repeat {repeat + 1}/{n_repeats}: "
                      f"{len(null_events)} null events", flush=True)

            null_active = np.concatenate(redetected_active) if redetected_active else np.array([0])
            same_window = np.vstack(same_window_active)  # (n_repeats, n_observed)
            null_per_event_mean = same_window.mean(axis=0)
            null_per_event_p95 = np.percentile(same_window, 95, axis=0)
            exceeds = obs_part.n_active > null_per_event_p95

            rows.append({
                "session": name, "epoch": epoch, "n_pyr": len(pyr),
                "jitter_window_ms": cfg["jitter_window_ms"], "n_repeats": n_repeats,
                "observed_events": len(observed),
                "observed_rate_per_s": len(observed) / observed.state_seconds,
                "null_events_mean": float(np.mean(redetected_counts)),
                "null_rate_per_s": float(np.mean(redetected_counts)) / observed.state_seconds,
                "event_ratio": len(observed) / max(float(np.mean(redetected_counts)), 1e-9),
                "observed_active_median": float(np.median(obs_part.n_active)),
                "null_active_median": float(np.median(null_active)),
                "observed_active_p90": float(np.percentile(obs_part.n_active, 90)),
                "null_active_p90": float(np.percentile(null_active, 90)),
                "same_window_null_median": float(np.median(null_per_event_mean)),
                "pct_events_above_null_p95": 100 * float(np.mean(exceeds)),
                "median_excess_cells": float(np.median(obs_part.n_active - null_per_event_mean)),
            })
            for i in range(len(observed)):
                per_event_rows.append({
                    "session": name, "epoch": epoch, "peak": observed.peak[i],
                    "observed_active": obs_part.n_active[i],
                    "null_active_mean": null_per_event_mean[i],
                    "null_active_p95": null_per_event_p95[i],
                    "exceeds_null_p95": bool(exceeds[i]),
                })

    df = pd.DataFrame(rows)
    results = config.results_dir()
    df.round(4).to_csv(results / "mua_jitter_null.csv", index=False)
    pd.DataFrame(per_event_rows).round(4).to_csv(
        results / "mua_jitter_per_event.csv", index=False)
    (results / "mua_jitter_null.json").write_text(json.dumps({
        "script": "scripts/mua_jitter_null.py", "git_commit": config.git_commit(),
        "created": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "jitter_window_ms": cfg["jitter_window_ms"], "n_repeats": n_repeats,
    }, indent=2, default=str) + "\n")

    print()
    header = (f"{'session':<18}{'epoch':<6}{'obs ev':>8}{'null ev':>9}{'ratio':>7}"
              f"{'obs act':>9}{'null act':>9}{'same-win':>10}{'>null p95':>11}")
    print(header)
    print("-" * len(header))
    for _, r in df.iterrows():
        print(f"{r.session:<18}{r.epoch:<6}{r.observed_events:8d}{r.null_events_mean:9.0f}"
              f"{r.event_ratio:7.1f}{r.observed_active_median:9.0f}{r.null_active_median:9.0f}"
              f"{r.same_window_null_median:10.1f}{r.pct_events_above_null_p95:10.0f}%")
    print("-" * len(header))
    print("(act = median active pyramidal cells per event; same-win = null participation\n"
          " recounted inside the observed event windows)")


if __name__ == "__main__":
    main()
