# Parameter and design decisions

The full methods of Grosmark & Buzsáki (2016) are in an online supplement not available to
this project. Every parameter whose value is not stated in a source available to us is
marked **provisional** here and in `config/params.yaml`.

Provenance tags used in the config file:

| Tag | Meaning |
|-----|---------|
| `[CHEN2016]` | Stated in Chen, Grosmark, Penagos & Wilson 2016, *Sci Rep* 6:32193 |
| `[GB2016]` | Stated in the main text of Grosmark & Buzsáki 2016, *Science* 351:1440 |
| `[DATA]` | To be confirmed from the `.mat` files themselves |
| `[PROVISIONAL]` | Chosen by us; rationale below |

---

## Step 1 — repository and environment

### D1.1 — Repository sits beside the data, not around it

The repo root is `Place Field/hc11-place-fields/`, a sibling of `hc-11_Data/` and
`Article/`. The 290 MB of `.mat` files and the copyrighted paper PDFs stay outside the
repository entirely, so no `.gitignore` mistake can ever commit them. The documented
public layout is still `<repo>/data/raw/NoveltySessInfoMatFiles/`; `hc11.config.data_dir`
resolves `$HC11_DATA_DIR` first, then that public path, then the local sibling layout, so
the repo works unmodified for someone who clones it and follows the README.

### D1.2 — `pyproject.toml` + uv rather than `environment.yml`

Phase 1 needs no non-Python binaries (`h5py` wheels ship their own HDF5), so conda buys
nothing. uv gives a committed `uv.lock` for exact reproducibility, and a `src/` layout
package that a reviewer can `pip install -e .` without uv if they prefer.

### D1.3 — One config file, no defaults in function signatures

Analysis functions take their parameters explicitly; the values come from
`config/params.yaml`. A parameter that appears as a default in a function signature is a
parameter that will silently disagree with the config file six months from now.
`hc11.config.param` raises `KeyError` on an unknown dotted key rather than returning a
fallback, so a renamed parameter fails loudly.

### D1.4 — Figures carry their own provenance

`output.write_sidecar_json: true`. Each figure gets `<name>.json` recording the parameters
used and the git commit (with a `-dirty` suffix when the working tree is not clean). This
is what makes "which version of the code made this figure" answerable after the fact.

---

## Step 2 — loader

### D2.1 — Transpose unconditionally, never by shape

The loader reverses HDF5 axes for *every* numeric dataset (MATLAB array = `h5py.T`) and
then normalises each field to its documented logical shape, raising if it does not fit.
Guessing orientation from shape ("if it has 2 rows, transpose") silently mis-reads a
`2 × 2` interval matrix, and the real files are not even internally consistent
(`TimeStamps` is a row vector, `OneDLocation` a column). A test pins the `2 × 2` case.

### D2.2 — Two layers: generic MATLAB reader, then hc-11 assembly

`read_matlab` knows nothing about hc-11 (structs, chars, empties, references, logicals);
`load_session` knows nothing about HDF5. Each can be tested on its own, and the generic
layer handles `MATLAB_empty` placeholders and object references even though no current
file uses them.

### D2.3 — The loader never modifies data

No clipping, no dropping, no interpolation, no unit conversion of stored arrays. Metres
stay metres (`position_xy_cm` / `position_1d_cm` are derived views). Anomalies are
reported by `validate()`, not "fixed" in the loader, so that each fix is an explicit,
reviewable analysis decision.

### D2.4 — `validate()` returns issues with a severity instead of raising

`error` = an invariant downstream code relies on is broken; `warning` = a real anomaly to
account for; `info` = a measurement for the record. The integration tests assert zero
errors **and** pin the exact set of warning codes per session, so a new anomaly (or a
disappearing one, e.g. after a data re-download) fails the suite.

Spikes outside `[0, sess_duration]` are a *warning*, not an error: the ask was to
report them, and every Phase 1 analysis restricts spikes to an explicit epoch anyway,
which already excludes them.

### D2.5 — No cache

Loading a session takes 0.2–0.7 s straight from the gzip-compressed HDF5. An `.npz`
cache would save ~0.3 s per session at the cost of an invalidation problem. Not worth it;
revisit if Phase 2 shuffles need thousands of reloads (they won't — load once, shuffle
in memory).

### D2.6 — Cross-check with an independent reader

`pymatreader` (a dev dependency only) is used in `tests/test_io.py` to re-read every
session and compare every field for exact equality with `hc11.io`. A physical check adds
a second, reader-independent line of evidence: on linear tracks `OneDLocation` must be a
linear function of `TwoDLocation` (|r| > 0.999 against the first principal axis), which
would fail if the position arrays were transposed or mis-paired.

---

## Step 3 — anomaly resolution and behaviour

James's decisions, 2026-09-30. All three are implemented in an analysis layer
(`preprocess.py`, `track.py`, `behavior.py`); `load_session` still returns raw values.

### D3.1 — Truncate Gatsby_08022013 at `sessDuration`

`truncate_to_session(session)` drops spikes after `sess_duration` and clips state
intervals to it, returning a new `Session` plus a `TruncationReport` that is logged.
Only Gatsby_08022013 is affected: 263 976 spikes (4.67 %) from 80 clusters, up to
1 588.7 s past the end, plus 4 state intervals dropped and 1 clipped (1 578.4 s).

Rationale: `sessDuration`, `POSTEpoch` and `Sessions_Recordings_Summary.pdf` agree that
the session ends at 30 413.628 s, so the extra data lies outside every documented epoch
and its provenance is unknown. Keeping it in would inflate whole-session firing rates by
~5 % for that session only. Implemented as a separate function, not inside the loader, so
that reversing the decision means not calling it.

### D3.2 — Drop zero-length state intervals

`clean_intervals` / `state_intervals` drop rows with `end <= start`, logging what was
dropped. Only Achilles_11012013 REM row 12 (`[13281, 13281]`) is affected. A zero-length
interval contributes nothing to a duration sum but can yield a degenerate result in
interval intersection or sample masking.

### D3.3 — Circular tracks: observed range, explicit wrap **(open question)**

On the three circular sessions, `OneDLocation` spans 2.84–2.90 m, while the documented
1 m-diameter platform implies a 3.14 m circumference, and the median 2-D radius
(0.45–0.53 m) implies 2.82–3.34 m. **The cause is unknown.** Possibilities not
distinguished by the released data: the linearisation may cover only the portion of the
ring the animal actually ran, it may be scaled to a nominal path radius rather than the
measured one, or the reward zone may be excluded.

Consequences for the code:

- `Track.from_session` takes the extent from the **observed** `OneDLocation` range per
  session, never from an assumed circumference. An assumed 3.14 m would misplace every
  bin edge and put the wrap discontinuity at a position the animal never occupied.
- The wrap is explicit: `Track.difference` / `Track.distance` take the short way around,
  so two positions either side of the reward site at 0 are correctly ~0.1 m apart rather
  than a full lap. `Track.unwrap` gives continuous position for direction and lap
  detection. Unit-tested in `tests/test_track.py::TestWrap`.
- The circular period is taken as `hi - lo`, which underestimates the true circumference
  by at most one inter-sample step (< 1.5 cm at observed speeds) — well inside a 10 cm bin.

This stays an open question; it does not affect Phase 1 linear-track place fields.

### D3.4 — Running = the authors' mask ∩ our speed threshold (a substantive choice)

Stage 1 is `OneDLocation` being defined; stage 2 is speed > 15 cm/s (2-D tracking,
differentiated then smoothed with σ = 0.1 s). Short dips (≤ 0.2 s) are bridged, and runs
under 0.5 s dropped.

**This is a methodological choice, not a confirmation of the released mask.**
`OneDLocation` marks *being on the linearised track*; how strictly it was curated varies
from session to session, so it does not reliably carry a running criterion, and the speed
threshold supplies one the data does not otherwise provide. Treating stage 2 as a
formality would mean accepting a definition of "running" that silently differs between
sessions — the opposite of a controlled comparison across the eight of them.

The evidence: the speed threshold removes

| Session | Removed by speed | Samples < 5 cm/s | Mask block duration (median / max) |
|---|---|---|---|
| Buddy_06272013 | 3.3 % | 0.6 % | 1.6 s / 5.7 s |
| Achilles_10252013 | 5.5 % | 0.3 % | 1.9 s / 6.9 s |
| Gatsby_08022013 | 11.0 % | 1.1 % | 2.8 s / 15.4 s |
| Gatsby_08282013 | 30.8 % | 8.4 % | 1.4 s / 30.6 s |
| Achilles_11012013 | 36.2 % | 9.0 % | 0.7 s / 30.2 s |
| Cicero_09102014 | 39.5 % | 10.6 % | 0.4 s / 23.2 s |
| Cicero_09172014 | 50.4 % | 12.1 % | 2.2 s / 29.7 s |
| Cicero_09012014 | 51.9 % | 16.5 % | 0.8 s / 58.7 s |

In the first three sessions the mask is made of short blocks the length of a single
traversal, and 71–87 % of the slow samples sit in the end zones — i.e. the acceleration
and deceleration at the ends of a run, which any speed threshold trims. In the other five
the mask contains blocks of up to 23–59 s with 8–17 % of samples below 5 cm/s spread
across the middle of the track, i.e. genuine pauses.

In the first three sessions the mask is made of short blocks the length of a single
traversal; in the other five it contains blocks of up to 23–59 s with genuine mid-track
pauses. The stage-1/stage-2 split is kept precisely so this stays visible per session
(`results/running_summary.csv`).

What it does **not** change is the place fields themselves — see D4.5.

### D3.5 — Refractory-period violations cannot be a quality metric

The minimum inter-spike interval within a cluster is **exactly 0.85 ms in all 8
sessions** — an implausible coincidence unless a refractory censoring step was applied
before release. Standard cluster-quality measures based on ISI violations (refractory
contamination rate, fraction of ISIs below 2 ms) are therefore uninformative here: they
measure the upstream cleaning, not the isolation quality of the cluster. If unit quality
needs assessing in Phase 2, it has to come from waveform or amplitude data in the
per-session archives, which this phase does not download.

---

## Step 4 — place fields

### D4.1 — Laps are segmented before the speed filter, not after

`segment_laps` runs on the authors' mask blocks; the speed filter enters later, in
`select_running_samples`, when occupancy and spikes are accumulated. Segmenting on the
speed-filtered mask would cut a traversal with a brief mid-track slowdown into several
fragments, each covering too little of the track to qualify as a lap, so a session's lap
count would depend on the speed threshold. Pinned by
`tests/test_laps.py::test_brief_slowdown_does_not_split_a_lap`.

### D4.2 — Direction by hysteresis, with no assumption of alternation

Direction is the sign of a smoothed (σ = 0.25 s) derivative of the linearised position,
with a ±5 cm/s hysteresis band: the direction flips only when the animal clearly moves
the other way. Runs are then split wherever the direction changes, so two traversals in
the same direction in a row are two laps of that direction — nothing assumes the animal
alternates. A lap must cover ≥ 50 % of the track extent end to end and last 0.5–60 s.

### D4.3 — On a circular track, one lap is one circuit

A constant-direction run is split again at every crossing of the reward site at position
0. Without this, an animal running continuously round the ring yields a single "lap"
spanning many circuits, with a coverage of 2, 3, … Found by a test written for the wrap,
not by inspection of the data.

Consequence for the circular sessions: their laps are strongly unidirectional (1/75,
0/19, 1/79 for pos/neg), which matches the data description — the animals were
"gently encouraged to run unidirectionally" on the circular platform. Their median lap
coverage is 0.66–0.99: the linearised position is undefined for much of the ring in two
of them, so many segments are partial circuits.

### D4.4 — Circular-shift shuffle in compressed running time

The null shifts each cell's spike train circularly through the *selected running samples*
rather than through session time. This preserves the spike count and the trajectory,
destroys the pairing between them, and keeps every shuffled spike at a position the animal
actually occupied while running. Spikes are attached to the nearest position sample
(≈ 1.5 cm at observed speeds, versus 10 cm bins), which also keeps them on the grid the
shuffle operates on.

**Known limitation, pinned as a test.** The shuffle's power depends on lap-to-lap
variability. If every lap took exactly the same time, a shift of *k* slots would map bin
*i* onto bin *i + k* in every lap, so the shuffled map would be a rotation of the true one
and would carry the same spatial information — the null would sit at the observed value
and no cell would be significant. Real lap durations vary (Achilles_10252013: 2.0–6.9 s),
which is what gives the test its power. The failure mode is conservative, never
anti-conservative. See `tests/test_fields.py::test_stereotyped_running_weakens_the_shuffle`
and the calibration test alongside it, which checks that spatially uniform cells are
flagged at roughly the nominal 5 % rate rather than more often.

### D4.5 — The speed threshold barely changes the place fields

`scripts/sensitivity_running.py` rebuilds the maps three ways for one session from each
group in D3.4 (laps held fixed, so only the samples differ):

| Session | Condition | Running | % of mask | Place cells | Median info | Median stability | Median peak |
|---|---|---|---|---|---|---|---|
| Achilles_10252013 | mask only | 4.3 min | 96 % | 76 | 0.62 | 0.97 | 5.5 Hz |
| Achilles_10252013 | mask + 15 cm/s | 4.0 min | 90 % | 76 | 0.62 | 0.97 | 5.5 Hz |
| Achilles_10252013 | mask + 5 cm/s | 4.2 min | 96 % | 76 | 0.62 | 0.97 | 5.5 Hz |
| Cicero_09012014 | mask only | 19.6 min | 97 % | 16 | 0.54 | 0.97 | 5.6 Hz |
| Cicero_09012014 | mask + 15 cm/s | 8.4 min | 41 % | 16 | 0.54 | 0.97 | 6.9 Hz |
| Cicero_09012014 | mask + 5 cm/s | 16.0 min | 79 % | 16 | 0.53 | 0.97 | 5.8 Hz |

**Mask-only fields are not noticeably worse.** Even in Cicero_09012014, where the
threshold discards 59 % of the samples, the place-cell count, median spatial information
and median stability are unchanged to two decimal places. Cell membership is identical in
Achilles_10252013 and 76 % overlapping in Cicero_09012014 (3 cell-directions swap each
way). The one thing that does move is the median peak rate (5.6 → 6.9 Hz in
Cicero_09012014), which is what you would expect from removing slow, low-rate occupancy
from the denominator.

So the threshold is not doing necessary work for *identifying* place cells. It is kept as
primary anyway, for the reason in D3.4 — it makes "running" mean the same thing in all
eight sessions — and because rate estimates do depend on it. But the fields themselves are
robust to it, which is a reassuring negative result rather than a justification.

### D4.6 — Place-cell criteria, and why our count is below the paper's

Criteria (all in `config/params.yaml`, all provisional): ≥ 50 spikes during running in that
direction, peak rate ≥ 1 Hz, spatial information significant against the shuffle null at
α = 0.05, at least one detected field, and split-half (odd vs even laps) r ≥ 0.3. A cell
counts once per session if it qualifies in either direction. Counts are reported both with
and without the stability criterion, because the paper states no threshold for it; in
practice it changes nothing (cells that pass the other criteria are highly stable).

**Result: 296 place cells across the 8 sessions, against the paper's 491 — a factor of
1.66, outside the factor-1.5 tolerance James set.** Per the working rules, the criteria
have *not* been adjusted to close the gap. The diagnosis:

| Criterion | Cell-directions failing it alone |
|---|---|
| no detected field | 101 |
| fewer than 50 spikes | 32 |
| peak rate < 1 Hz | 6 |
| unstable (r < 0.3) | 3 |
| information not significant | 2 |

The dominant constraint is field *detection*, not spatial tuning: 105 cell-directions have
significant spatial information, a peak above 1 Hz and enough spikes, yet no field, because
the region above 20 % of their peak is wider than the provisional 120 cm cap — on a 160 cm
track that cap excludes any cell tuned more broadly than 75 % of the track. Relaxing only
that cap gives 335 place cells (1.47x); additionally lowering the spike floor to 20 gives
356 (1.38x). Both are inside the tolerance, and both are single-parameter changes to
numbers we invented rather than took from a source — which is exactly why they need
James's decision rather than ours.

Provisional choices already recorded in `config/params.yaml`, to be revisited once real
maps exist:

| Parameter | Value | Reasoning |
|-----------|-------|-----------|
| `behavior.speed_smoothing_sigma_s` | 0.1 s | Speed differentiated from ~39 Hz tracking is noisy enough to chatter across a 15 cm/s threshold; ~4 samples of Gaussian smoothing suppresses that without blurring real run onsets. |
| `behavior.min_run_epoch_s` | 0.5 s | Shorter "runs" contribute too few spikes to inform a tuning curve and mostly reflect tracking jitter. |
| `place_fields.smoothing_sigma_bins` | 1.0 | Chen et al. describe smoothing over ~5 bins. Read as ~5 bins of kernel support, i.e. σ ≈ 1 bin with truncation at 4σ↓. An alternative reading is σ = 5 bins (50 cm), which would be unusually heavy for a linear track. Flagged for sensitivity check. |
| `place_fields.min_occupancy_s` | 0.1 s | Bins visited for less than ~4 tracking samples give a rate estimate dominated by a single spike; treated as unvisited (NaN) rather than 0 Hz, which would otherwise bias both spatial information and the decoder. |
| `place_fields.criteria.*` | see config | Standard CA1 place-cell criteria. The paper's own count (n = 491 place cells across 8 sessions) is the target used to check whether these are too strict or too loose. |
| `units.min_spikes_maze` | 100 | Below this, a tuning curve is not estimable in any useful sense. |

## Step 5 — decoding

| Parameter | Value | Reasoning |
|-----------|-------|-----------|
| `decoding.time_bin_s` | 0.25 s | Conventional bin for decoding *behaviour* (as opposed to the 20 ms bins used for replay in Phase 2). At 15 cm/s the animal moves ≈ 3.75 cm per bin, well under the 10 cm spatial bin, so binning does not itself dominate decoding error. |
| `decoding.rate_floor_hz` | 0.01 Hz | Without a floor, a single spike from a cell with a 0 Hz estimate at some position assigns that position zero posterior probability outright. |
| `decoding.cv_scheme` | leave-one-lap-out | Laps are the natural independent unit on a linear track; random bin-wise splits leak information through temporal autocorrelation of position. |
| `laps.end_zone_frac` | 0.15 | Enough of the track end to capture the turnaround reliably despite tracking dropout at the reward sites. |
