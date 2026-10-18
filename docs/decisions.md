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

## Step 5b — Bayesian decoding

### D5.6 — Memoryless Poisson decoder, uniform prior, log space

`log P(x | n) = sum_i n_i log f_i(x) - dt sum_i f_i(x) + const`, normalised with a
log-sum-exp. The `-dt sum_i f_i(x)` term matters and is easy to drop by accident: it is
what makes *too few* spikes evidence *against* a high-rate position (pinned by
`test_too_few_spikes_is_evidence_against_a_high_rate_position`). Tuning curves are floored
at `rate_floor_hz = 0.01` so one spike cannot veto a position outright, and NaN bins
(never visited in the training laps) take the same floor.

The assumptions — Poisson firing, conditional independence across cells, no transition
model — are documented in the module docstring together with where they fail. In short:
bursting and assembly structure both make the posterior sharper than the evidence
warrants, which damages calibration more than the peak position, so decoding *error*
stays usable while posterior probabilities should not be read as confidences.

### D5.7 — Leave-one-lap-out, with the templates precomputed

Each lap is decoded against fields built from every lap except itself — no spikes and no
occupancy from the held-out lap reach its own template. Three tests enforce this: the
template for lap *i* is compared against an independently computed one; spikes fabricated
into lap *i* change every other template and the all-laps template but not lap *i*'s own;
and the honest error must be no better than a deliberately leaky template built from all
laps. On the real data that penalty is **+0.89 cm on average** — small, because one lap of
40 barely moves a template, but consistently in the right direction.

Templates are precomputed once per session/direction/bin size, so the chance baseline
(which permutes their rows) and the ensemble sweep (which slices them) are nearly free.

### D5.8 — Chance baseline by permuting tuning curves between cells

Rather than shuffling spikes, the baseline reassigns which field belongs to which cell,
keeping the spike counts and the set of fields intact and destroying only the
correspondence. It is applied to the same cross-validated templates, so the baseline is
cross-validated too. It lands at **36–71 cm**, against a uniform-guess expectation of
53 cm on the 1.6 m track — i.e. permuted fields are no better than guessing, as they
should be.

### D5.9 — Results on running data

Primary condition, 250 ms bins, each direction against its own fields:

| Session | Cells | Median error | Mean | p90 | Chance |
|---|---|---|---|---|---|
| Achilles_10252013 | 77 | **4.7–5.1 cm** | 7.2–8.0 | 13.6–16.2 | 40–43 |
| Achilles_11012013 | 73 | 4.9 | 10.7 | 16.0 | 67 |
| Buddy_06272013 | 16 | 8.8 | 14.5–16.2 | 31–34 | 36–45 |
| Cicero_09012014 | 16 | 10.6 | 17.7–21.4 | 44–68 | 38–53 |
| Cicero_09102014 | 36 | 4.9 | 13.3 | 16.6 | 71 |
| Cicero_09172014 | 25 | 9.0–9.3 | 20.5–21.2 | 54–58 | 49–63 |
| Gatsby_08022013 | 30 | 6.6–7.9 | 11.9–16.7 | 26–43 | 38–43 |
| Gatsby_08282013 | 31 | 5.4 | 16.3 | 45.5 | 68 |

**The median error of 4.7–10.6 cm is better than the "low tens of centimetres" we
expected, and that is not leakage.** Three things explain it. The cross-validation checks
above all pass. The spatial bin is 10 cm and the decoder returns a bin centre, so even
perfect bin identification yields a median error near 2.5 cm — our 5 cm means the right
bin or an adjacent one. And a 1.6 m linear track with 77 cells is a far easier problem
than the 2-D open fields where Chen et al. report 8.5–12.5 cm with 49 neurons.

Error is **lower at the track ends** (6.1 cm) than in the middle (8.3 cm), the opposite of
the usual expectation. The ends are where the animal slows and turns, so more time bins
accumulate there per unit distance, and the boundary itself restricts where the posterior
can place the animal.

### D5.10 — The directional split is worth 6–34 %

Same cells, same laps, only the template differs:

| Session | Directional | Merged | Cost of merging |
|---|---|---|---|
| Achilles_10252013 | 4.9 cm | 5.9 cm | +19 % |
| Buddy_06272013 | 8.8 | 9.3 | +6 % |
| Cicero_09012014 | 10.6 | 12.6 | +19 % |
| Cicero_09172014 | 9.2 | 10.1 | +10 % |
| Gatsby_08022013 | 7.2 | 9.7 | +34 % |

On the three circular sessions the two are identical to within 0.3 cm, as they must be:
those sessions are unidirectional, so there is nothing to merge.

### D5.11 — 20 ms bins cost a factor of two, and more than that in the tail

| | 250 ms | 20 ms |
|---|---|---|
| Median error | 4.7–10.6 cm | 10.0–24.0 cm |
| Undecodable bins | 0–3 % | 5–60 % |
| Median spikes per bin | 10–40 | **0–4** |

At 20 ms the median spike count per bin is between 0 and 4 across sessions, so a large
share of bins carry no information at all (60 % undecodable in Cicero_09012014, which has
16 cells). The p90 error reaches 84–126 cm, at or above the chance level, meaning the
worst decile is pure guesswork. This is expected rather than a defect: 20 ms is the bin
Grosmark & Buzsáki use for *ripple* events, where ~10–20× temporal compression packs a
whole trajectory into a few hundred milliseconds and each bin is read as part of a
sequence rather than as an independent position estimate. The number to carry into
Phase 2 is that a 20 ms bin on this dataset is worth roughly one spike, so sequence
scoring must lean on the trajectory across bins, not on any single posterior.

### D5.12 — Ensemble size

Achilles_10252013, 250 ms: 12.1 cm with 5 cells, 9.7 with 10, 7.1 with 20, 5.9 with 30,
5.5 with 40, 5.2 with all 55. The curve flattens past ~30 cells, which is why the
16-cell sessions (Buddy, Cicero_09012014) sit at 9–11 cm while the 73–77-cell sessions
reach 5 cm. This is the dependence Chen et al. emphasise, and it means our place-cell
shortfall (D5.4) costs decoding accuracy directly: sessions where we admit fewer cells
decode measurably worse.

| Parameter | Value | Reasoning |
|-----------|-------|-----------|
| `decoding.time_bin_s` | 0.25 s | Conventional bin for decoding *behaviour* (as opposed to the 20 ms bins used for replay in Phase 2). At 15 cm/s the animal moves ≈ 3.75 cm per bin, well under the 10 cm spatial bin, so binning does not itself dominate decoding error. |
| `decoding.rate_floor_hz` | 0.01 Hz | Without a floor, a single spike from a cell with a 0 Hz estimate at some position assigns that position zero posterior probability outright. |
| `decoding.cv_scheme` | leave-one-lap-out | Laps are the natural independent unit on a linear track; random bin-wise splits leak information through temporal autocorrelation of position. |
| `laps.end_zone_frac` | 0.15 | Enough of the track end to capture the turnaround reliably despite tracking dropout at the reward sites. |

---

## Step 5a — place-cell count: convention, criteria, and the residual gap

### D5.1 — The paper's 491 counts unique neurons, so ours must too

Grosmark & Buzsáki: *"a spatial Bayesian decoder, constructed from the firing-rate vectors
of place cells (n = 491 cells) during track running (eight novel exploration sessions)"*.
Three things fix the convention as **unique neurons**, not cell-direction pairs:

1. The same 491 is later divided *"into equal subgroups by either 'off-line' sleep-firing
   rate or ripple-rate gain"* — properties a neuron has once, not once per direction.
2. *"n = 216 neurons"* is described as a subset of them (Fig. 2B, one point per neuron).
3. Per-cell contribution (PCC) is defined per neuron.

Our pipeline reports both. The ratio between them is only 1.21 overall, not the ~2 a
linear track might suggest, because the three circular sessions are unidirectional (ratio
exactly 1.00) and on the linear tracks many cells qualify in one direction only
(1.31–1.57). **The convention does not close the gap**: 304 unique vs 370 cell-direction
pairs, against the paper's 491 unique.

### D5.2 — Field-width cap as a fraction of track extent

`criteria.max_field_width_frac: 0.75` replaces the fixed `max_field_width_cm: 120`, which
was inconsistent across mazes — 75 % of a 1.6 m track, 60 % of the 2 m track, but only
41 % of a 2.84–2.90 m ring, so the same cell shape was judged differently depending on
which maze it was recorded on.

Why 0.75 and not something tighter:

- A field is defined as the contiguous run above 20 % of the peak rate. A textbook CA1
  field with σ = 20 cm already spans 2σ·√(2 ln 5) ≈ 72 cm at that threshold — 45 % of a
  1.6 m track before the 10 cm smoothing kernel widens it further. A cap at 0.5 would
  reject ordinary fields.
- A spatially uniform cell spans 100 % by construction.
- 0.75 sits between the widest plausible genuine field and total non-localisation, and it
  reproduces the previous behaviour on the 1.6 m tracks, so this change is about
  consistency across track lengths rather than about loosening the criterion.

Effect: +8 unique place cells (296 → 304), all on the 2 m and circular sessions.
Provisional, like every other criterion we invented.

### D5.3 — The 50-spike floor stays primary

`min_spikes: 50` during running, per direction. The 20-spike variant is reported in
`results/place_cell_sensitivity.csv` as a sensitivity row only and is **not** adopted.

### D5.4 — The gap is real and cannot be closed by relaxing our criteria

`scripts/place_cell_sensitivity.py`, all counts across the 8 sessions (562 pyramidal
cells):

| Variant | Unique cells | % of pyramidal | 491 / ours | Cell-direction pairs |
|---|---|---|---|---|
| **primary** | **304** | 54 % | **1.62** | 370 |
| cap relaxed (no max width) | 346 | 62 % | 1.42 | 458 |
| cap relaxed + 20-spike floor | 367 | 65 % | 1.34 | 489 |
| peak rate ≥ 1 Hz only, nothing else | 392 | 70 % | 1.25 | 618 |
| Grosmark & Buzsáki | 491 | 87 % | — | — |

The last row of our own results is the important one. Chen et al. (2016), whose methods
section describes the same recordings, state that *"all putative pyramidal neurons selected
for analysis had peak firing rate > 1 Hz"*. Applying **only** that — no field requirement,
no shuffle test, no stability, no spike floor — still yields 392, not 491. Measuring on
the most inclusive base available (every sample with a defined `OneDLocation`, directions
merged) gives 353. **There is no setting of our criteria that reaches 87 % of pyramidal
cells**, so the difference is not a threshold we have mis-set.

Where it plausibly comes from, none of which we can test with the released files:

- **The released `OneDLocation` covers only part of the behaviour.** It is defined for
  7–35 % of each MAZE epoch (docs/data_schema.md), so cells are scored on a fraction of
  the running the authors had. 186 of 562 pyramidal cells fire fewer than 50 spikes in
  their best direction within our running periods, and 10 fire none at all. Place-cell
  fraction and running time rise together across the eight sessions (r = 0.63), which is
  *consistent with* this explanation but is not evidence for it: with n = 8 that
  correlation is not statistically meaningful (p ≈ 0.09), and the sessions differ in
  several other ways at the same time. Tested directly in D5.5.
- **Field construction may differ.** Two-dimensional fields, or fields built over all MAZE
  time rather than running periods, would admit cells that our 1-D running-only maps
  cannot evaluate.
- **The supplement is unavailable**, and with it the actual inclusion criteria.

Per the working rules the criteria were **not** tuned toward 491. The gap is recorded as a
finding: our place-cell population is a conservative subset of theirs, and step 5b's
decoder is built on that subset.

### D5.5 — Diagnostic: Fig. 1A's 77 cells, and the measurement base tested directly

Two checks on D5.4's reasoning. Both are **diagnostics**; the authors' `OneDLocation`
remains the replication base in every primary result, because it is their linearisation.

**Fig. 1A is not a match.** The legend reads *"Simultaneous recording of 77 place cells
(rightward runs) used to generate a sequence template"* — 77 in a **single** direction. The
figure is almost certainly Achilles_10252013: only three sessions have ≥ 77 pyramidal
cells (120, 92, 81) and the other two are circular, where "rightward" has no meaning. Our
counts for that session:

| | Place cells |
|---|---|
| increasing position | 55 |
| decreasing position | 52 |
| **unique (either direction)** | **77** |
| both directions | 30 |

Neither direction approaches 77; our *unique* count equalling their *one-direction* count
is a coincidence and should not be read as a match. What the figure does give us is a
per-session, per-direction anchor that is much tighter than the 491 aggregate: on their
best session they admitted 77 of 120 pyramidal cells (64 %) in one direction, where we
admit 55 (46 %). They are roughly 1.4× more inclusive per direction — the same factor as
the overall gap, which at least says the discrepancy is uniform rather than concentrated
in the aggregate.

**The measurement base explains part of the gap, not all of it.**
`scripts/diagnostic_linearization.py` rebuilds everything on our own linearisation of the
2-D position over all MAZE time with valid tracking (`src/hc11/linearize.py`), with the
primary criteria unchanged:

| Session | Maze | Coverage theirs → ours | Running min | Unique place cells |
|---|---|---|---|---|
| Achilles_10252013 | 1.6 m linear | 13 % → 79 % | 4.0 → 11.4 | 77 → 83 |
| Buddy_06272013 | 1.6 m linear | 7 % → 97 % | 2.6 → 17.4 | 16 → **32** |
| Cicero_09012014 | 1.6 m linear | 22 % → 97 % | 8.4 → 17.3 | 16 → **29** |
| Cicero_09172014 | 2 m linear | 32 % → 68 % | 7.3 → 8.6 | 25 → 28 |
| Gatsby_08022013 | 1.6 m linear | 19 % → 98 % | 6.3 → 26.7 | 30 → 35 |
| **linear total** | | | | **164 → 207** (47 % → 59 % of 348) |
| Achilles_11012013 | circular | 58 % → 81 % | 15.4 → 19.9 | 73 → 65 |
| Cicero_09102014 | circular | 29 % → 46 % | 7.6 → 9.1 | 36 → 32 |
| Gatsby_08282013 | circular | 36 % → 95 % | 10.9 → 19.8 | 31 → 14 |

Read the linear rows only. On those five sessions our linearisation reproduces the
authors' exactly (r = 1.0000 against `OneDLocation` wherever both are defined), so the
comparison isolates the base. On the ring our arc-length reconstruction is merely close
(r = 0.972–0.997) and its extent disagrees with theirs, so the circular rows measure the
quality of our linearisation rather than the effect of the base — which is why the
all-session total (304 → 318) understates the effect.

Conclusion: widening the base from 7–32 % to 68–98 % of each linear MAZE epoch raises the
yield from 47 % to 59 % of pyramidal cells. Buddy_06272013 doubles (16 → 32), exactly the
session where the released mask was thinnest at 2.6 min of running. So the restricted
`OneDLocation` mask is a real contributor to the gap — but 59 % is still well short of the
paper's 87 %, so it is not the whole story, and the remainder stays unexplained.

---

## Step 6c (first half) — population-synchrony detection from spikes

Built before the LFP arrived, so the two detectors can be compared as soon as it does.

### D6.1 — One detection rule, shared by both detectors

`events.detect_two_threshold` is used by the MUA detector now and by the ripple detector
later: peak above a high threshold, boundaries extended out to a low threshold, merge,
then filter on duration. Both return the same `EventList`, so the 6c cross-check is a
comparison of two event lists and nothing else. A difference in how boundaries are
defined would otherwise show up as a difference in the events.

### D6.2 — Baselines per state and per interval

Thresholds come from the baseline of the epoch's own brain state, pooled across that
state's intervals, never session-wide: baseline population rate differs between non-REM
sleep and quiet waking, and a session-wide threshold would hand one epoch systematically
more events than another — which is exactly the comparison Fig. 1C rests on. Detection
then runs inside each interval separately, so no event spans a gap in the searched time
and smoothing never bleeds across a state boundary. Pinned by
`test_baseline_is_computed_within_the_state`.

### D6.3 — MAZE immobility cannot be defined from the released position data **(blocked)**

Chen et al.'s quiet-wakefulness criterion is speed < 2 cm/s. On these two sessions that
yields **0.4 min of searchable time out of 34–45 min**, which is not a threshold problem:

| | Achilles_10252013 | Achilles_11012013 |
|---|---|---|
| MAZE epoch | 34.5 min | 44.6 min |
| samples below 2 cm/s | 22 % | — |
| but: contiguous runs below 2 cm/s | 0.6 min total, median 0.18 s | 0.6 min, median 0.18 s |
| tracking lost | 7.6 min (235 runs, max 26 s) | 8.7 min (480 runs, max 338 s) |
| searchable immobility found | 0.4 min | 0.4 min |

The sub-2 cm/s samples are overwhelmingly **turnaround zero-crossings** — the speed passes
through zero as the animal reverses — not pauses. Bridging brief excursions (which is
needed anyway, and is applied) barely helps. Meanwhile the genuinely quiet periods are plausibly
the long tracking-lost runs: the animal sits at a reward site and the head-mounted LED is
occluded or out of frame.

**Those runs are an open hypothesis, not a ruled-out one.** They cannot be used *yet*,
because a lost LED is indistinguishable from an animal that left the tracked area. But
the elevated population rate during them (259–389 Hz, against a non-REM baseline of
60–77 Hz) is **ambiguous evidence, not disqualifying**: awake sharp-wave ripples at a
reward site, with the LED occluded, would produce exactly that signature – and those are
precisely the MAZE events the paper cares about. Reading the high rate as "the animal was
active, so these are not quiet periods" would beg the question.

Either channel in the `.eeg` settles it. A high theta/delta ratio during those runs means
locomotion or attentive waking; a low ratio with ripple-band bursts means awake
sharp-wave ripples at the reward site. The EMG channel answers the same question
independently of the hippocampal LFP.

The scored states do not rescue it either: Wake covers 96–100 % of both MAZE epochs, with
1.4 min of Drowsy in one session and none in the other.

**So MAZE events are blocked until the LFP arrives.** The theta/delta ratio and the EMG
channel both live in the `.eeg` file, and either would define immobility without relying
on the LED. PRE and POST are unaffected — non-REM comes from the scored states.

### D6.4 — The mean + 3 SD threshold is self-referential

The baseline SD is computed over a trace that contains the events, so more events means a
larger SD means a higher threshold. On the real data the SD (109 Hz) exceeds the mean
(77 Hz) for this reason, and the rate trace is strongly skewed (median 37 Hz, p99 537 Hz,
max 1323 Hz; only 32 % of time is above the mean). Two consequences worth stating:

- A flat Poisson background still produces threshold crossings — about 0.35–0.68/s in a
  20-cell synthetic session. Event *count* alone is therefore not evidence of synchrony.
  What separates a real burst from a noise crossing is **participation**: a crossing
  involves a few cells firing once each, a burst involves many. Pinned by
  `test_injected_bursts_stand_out_from_background_noise`.
- Robust statistics would make this worse, not better: median + 3 MAD-equivalent SD is
  200 Hz against mean + 3 SD's 411 Hz, so it would roughly double the event count.

The threshold stays at 3 SD for now — it is the value Chen et al. give for the *ripple
envelope*, and we have no published value for a MUA detector — with the sweep in
`results/mua_threshold_sensitivity.csv` recorded so the choice can be revisited once the
LFP cross-check says which crossings are real.

### D6.5 – The events carry real fine-timescale synchrony (spike-jitter null)

The flat-Poisson comparison in D6.4 is about threshold arithmetic and nothing else: it
destroys the slow population envelope along with the synchrony, so it cannot distinguish
"real co-firing" from "a bursty envelope crossing a threshold". The null that can is a
**spike jitter**: displace each cell's spikes independently by U(-100, +100) ms, which
preserves every cell's rate and the envelope on timescales above ~200 ms while destroying
co-firing on the 57–68 ms scale of an event. 5 repeats, `mua.jitter_spikes`.

Three comparisons, all pointing the same way:

| | Achilles_10252013 PRE | POST | Achilles_11012013 PRE | POST |
|---|---|---|---|---|
| Observed events | 7046 | 3955 | 4709 | 3328 |
| Jitter-null events | 3619 | 1830 | 2781 | 1853 |
| **Observed / null** | **1.9x** | **2.2x** | **1.7x** | **1.8x** |
| Observed active cells, per event | 18 | 18 | 11 | 12 |
| Null active cells, *same windows* | 10.8 | 10.0 | 6.4 | 6.8 |
| **Events above the null's own 95th percentile** | **89 %** | **90 %** | **84 %** | **86 %** |

The per-event "same windows" row is the decisive one and the only apples-to-apples
comparison: participation is recounted inside the *observed* event windows using jittered
spikes, so duration and timing are held fixed and only the co-firing changes. Observed
participation is ~1.7x the null, and **84–90 % of individual events exceed the 95th
percentile of their own null** – these are not a slow envelope crossing a threshold.

**A trap worth recording.** Comparing the median participation of *re-detected* null
events against observed gives the opposite-looking answer – 22 vs 18 cells, i.e. the null
appears *more* synchronous. It is an artefact of detection, not of the data: jittering
smooths the rate trace, so null events are nearly twice as long (median 124 ms vs 68 ms)
and simply accumulate more distinct cells over their longer span. Normalised for
duration, the observed events are denser (2.50 vs 1.85 active cells per 10 ms). Any
comparison of participation between two event sets has to control for duration.

### D6.6 – POST has less non-REM sleep than PRE, in both sessions

| Session | Epoch | Duration | Wake | Drowsy | **NREM** | Intermediate | REM |
|---|---|---|---|---|---|---|---|
| Achilles_10252013 | PRE | 301.3 min | 7 % | 25 % | **56 %** | 2 % | 10 % |
| Achilles_10252013 | POST | 245.2 min | 17 % | 36 % | **37 %** | 2 % | 6 % |
| Achilles_11012013 | PRE | 322.2 min | 9 % | 45 % | **32 %** | 1 % | 13 % |
| Achilles_11012013 | POST | 243.8 min | 9 % | 59 % | **30 %** | 0 % | 2 % |

In absolute terms non-REM falls from 168.5 to 89.7 min and from 104.0 to 73.0 min. Part
of that is simply that the POST epoch is shorter (245 vs 301 and 244 vs 322 min), but the
*proportion* also falls in Achilles_10252013 (56 % to 37 %), with Drowsy and Wake taking
up the difference. Achilles_11012013 holds its non-REM share (32 % to 30 %) but its REM
collapses, from 27 % of scored sleep to 5 %.

This is the opposite of a naive post-learning expectation, and the most likely reason is
mundane: the animal has already slept for roughly five hours in PRE, so sleep pressure is
largely discharged by the time POST begins. We cannot test that here.

**Consequence for Fig. 1C.** Any epoch comparison has to be normalised by searched time,
never by raw event count – PRE would otherwise win on sleep duration alone. Event *rates*
are already reported per minute of state. The sequence-score distributions themselves are
unaffected (they compare per-event scores, not counts), but POST rests on roughly half as
many events as PRE in Achilles_10252013, so its distribution is the noisier of the two and
any test across epochs should account for the unequal n.


---

## Step 6a — LFP file access

### D6.7 — Two channel counts, never one variable

`LfpMetadata` exposes `n_channels_total` and `probe_channels` as separate, differently
named quantities, and nothing reuses one for the other.

- **`n_channels_total` is the stride.** Every byte offset, the memory map and the
  file-size check use it. The hc-11 `.eeg` holds "hippocampal LFP, EMG and accelerometer
  data", so it is strictly larger than the probe count — 136 vs 128 in the test fixture.
- **`probe_channels` is the selection pool.** Ripple-channel selection operates on these
  alone. An electromyogram channel has no pyramidal layer and no ripples, but it does
  have power in the 150–300 Hz band, so it would quietly win a power-based contest.

### D6.8 — A wrong channel count fails loudly

`LfpFile.check_size` divides the file size by `n_channels_total × 2` bytes and reports
the remainder; `require_consistent_size` raises unless the file is a whole number of
frames and its duration matches `sessDuration` within a tolerance.

This is strict rather than advisory because the failure mode is silent. An off-by-one
channel count does not merely misreport the duration: it shifts every frame by one
channel, so each "channel" becomes a slowly drifting mixture of the real ones. The result
has plausible amplitudes, plausible spectra, and no error anywhere — it would survive
ripple detection and produce events. `test_a_wrong_channel_count_is_detected_not_accepted`
pins this for 135, 137 and 128 channels against a 136-channel file.

A missing `lfpSamplingRate` likewise raises rather than defaulting to 1250 Hz: a wrong
rate rescales every timestamp, silently decoupling LFP events from spike times.

### D6.9 — Memory-mapped reads only

`LfpFile.read(start_s, end_s, channels)` maps the file and copies out only the requested
window, shaped `(time, channel)`. The full session is never materialised — at 136
channels Achilles_10252013 is about 11 GB. A 544 MB fixture in the test suite checks that
a one-second read touches ~0.3 MB; the test fails outright if the reader ever loads the
file.

### D6.10 — The `.eeg` files are extracted from the archives, not streamed

CRCNS ships the LFP as one gzipped tarball per session holding a single member,
`<session>/<session>.eeg`, 11.7 and 12.3 GB uncompressed. Two routes were available:
decompress on the fly and keep only the channels in use, or extract once to disk and
memory-map.

We extract. The streaming route was tried first and works — it reads a whole archive in
about 75 s — but it gives up random access, which every later step needs: detection runs
over the full recording on a handful of channels, and the characterisation figures need
to jump to an arbitrary event time. Extraction costs 24 GB of disk once and leaves
`LfpFile`'s memory map as the only access path in the codebase, which is the simpler
thing to reason about and the thing the tests cover.

`src/hc11/eeg_archive.py` does the extraction, and is resumable: it appends to the
destination and can stop at a chunk boundary on a deadline, with the next call skipping
the bytes already on disk. That is not incidental. A 12 GB write to this machine's data
volume takes about 170 s, which does not reliably fit in one run, and the first two
attempts at a plain `tar -x` were killed partway and left 4.6 GB and 3.5 GB fragments
with no indication that anything was wrong. Resuming is safe because a single-member tar
is a pure byte stream; a destination *larger* than the member raises instead of being
truncated, since that means it is not a prefix of this archive.

`member_info` reads the uncompressed size from the tar header — the front of the stream —
so `size_check` can run against the declared size before committing minutes to an
extraction.

### D6.11 — Both `.eeg` files match the XML and the session clock exactly

| Session | bytes | ÷ (134 × 2) | duration | `sessDuration` | difference |
|---|---|---|---|---|---|
| Achilles_10252013 | 11,678,469,572 | 43,576,379 frames, remainder 0 | 34,861.103 s | 34,861.103 s | +0.000 s |
| Achilles_11012013 | 12,270,786,020 | 45,786,515 frames, remainder 0 | 36,629.212 s | 36,629.212 s | +0.000 s |

`require_consistent_size` passes on both. The remainder of zero is the stronger of the
two results: 134 is the only channel count near the right value that divides either file
exactly, so the stride is confirmed by the data and not only by the XML.

### D6.12 — Channel grouping verified against `Channel_Orderings.pdf`

The XML's `spikeDetection` grouping for both Achilles sessions matches the dataset's own
channel-order document line for line, once its 1-based numbering is shifted to the
0-based numbering Neuroscope and this codebase use:

- 12 CA1 shanks of 10 channels — six on a left-hemisphere BUZ64SP probe, six on a right —
  which are spike groups 1–6 and 8–13.
- Two groups of four, spike groups 7 and 14, which the document calls **"Extra
  Channels"**. They are not shanks in the cell layer.
- Six channels outside the probe: neck-muscle EMG (0-based 128–130) and head-mounted
  accelerometer (131–133). These are exactly the channels absent from `spikeDetection`,
  which is how D6.7 identifies them without a label.

Two independent checks agree with the document. The EMG channels carry 9.8–21.8 of
150–300 Hz envelope against 2.4–5.0 for the accelerometer — muscle is broadband, an
accelerometer is not — and the "Extra Channels" groups carry no sorted units in
`sessInfo` and the lowest ripple-band power on the probe.

### D6.13 — Selection: highest mean 150–300 Hz envelope during non-REM, on a shank with units

The ripple channel for each shank is the one with the highest mean ripple-band envelope
over a sample of scored non-REM sleep. 150–300 Hz is Chen et al.'s band. The rest is
ours and provisional:

- **Sample, not the whole recording.** 40 windows of 15 s = 10 min, spread evenly through
  the non-REM of PRE and POST, deterministic (positions are equally spaced along total
  non-REM time, so there is no seed to record). Ranking 134 channels does not need the
  3–4 h of non-REM these sessions contain.
- **Non-REM only, PRE and POST.** MAZE immobility is a few minutes at most here (D6.3)
  and running theta bleeding into a short sample would tilt the profile.
- **Shanks with sorted units only.** The "Extra Channels" groups of D6.12 are excluded,
  by the data rather than by their size: a group with no sorted pyramidal cell is not in
  *stratum pyramidale*, whatever the document calls it. Were they included, spike group 7
  of Achilles_10252013 would have been "selected" at channel 62 with less ripple power
  than the EMG.

Chosen channels (0-based), with the two sessions side by side:

| Spike group | Anatomy | 10252013 | 11012013 |
|---|---|---|---|
| 1 | Left CA1 shank 1 | 2 | 3 |
| 2 | Left CA1 shank 2 | 10 | 10 |
| 3 | Left CA1 shank 3 | 20 | 20 |
| 4 | Left CA1 shank 4 | 31 | 31 |
| 5 | Left CA1 shank 5 | 41 | 41 |
| 6 | Left CA1 shank 6 | 56 | 57 |
| 8 | Right CA1 shank 1 | 69 | 71 |
| 9 | Right CA1 shank 2 | 77 | 79 |
| 10 | Right CA1 shank 3 | 84 | 86 |
| 11 | Right CA1 shank 4 | 95 | 95 |
| 12 | Right CA1 shank 5 | 105 | 107 |
| 13 | Right CA1 shank 6 | 116 | 118 |

The primary detection channel is the strongest of these: **channel 2** in
Achilles_10252013 and **channel 3** in Achilles_11012013, both on left CA1 shank 1.

### D6.14 — Why the selection is believed, given that the margins are small

The winner beats the next channel on its own shank by 1–14 % of mean envelope, which on
its own would not be worth much. Three things make the choice credible anyway, and none
of them is the margin.

1. **The depth profile is smooth and peaked.** On left shank 1 of Achilles_10252013 the
   envelope runs 40.3, 44.1, **45.8**, 43.5, 41.7, 37.7, 37.8, 32.7, 32.5, 29.8 down the
   shank. Every shank has that shape. A profile driven by per-channel noise would not be
   monotonic away from a single maximum, and a small margin between two adjacent channels
   20 µm apart is what a real layer looks like at this electrode spacing.
2. **Theta runs anti-phase to the ripple band.** Within each shank the theta and delta
   envelopes peak several channels *deeper* than the ripple envelope
   (`lfp_band_profile.png`). That is the CA1 depth signature — ripples maximal in
   *stratum pyramidale*, theta maximal towards *lacunosum-moleculare* — and it is not
   something a noisy channel can produce.
3. **The two sessions agree.** Seven of twelve shanks select the identical channel a week
   apart, and the rest differ by one or two. The right-hemisphere probe is the one that
   moved, consistently by two channels (40 µm).

**Sensitivity.** Ranking instead by *contrast* — 99th percentile over median envelope,
which measures excursion above a channel's own baseline rather than level — moves the
choice by one channel on 8 of 12 shanks in Achilles_10252013 and 5 of 12 in
Achilles_11012013. One disagreement is larger: spike group 3 of Achilles_10252013, mean
picks channel 20 and contrast picks 27, seven channels (140 µm) apart. Both rankings are
recorded in `results/ripple_channels.csv`; the mean is primary because it is the stated
criterion. Whether the choice matters at all is testable in step 6b and should be
tested there rather than argued here.

### D6.15 — The filter axis, and a bug worth recording

The first attempt at this analysis reported a mean ripple envelope of 44–46 on *every one*
of the 134 channels, with the theta/delta ratio NaN. It looked like corrupt data, and
several hours went into checking the recording: raw values varied smoothly, per-channel
s.d. ranged 646–944, neighbouring channels correlated 0.56–0.86, the Welch spectrum fell
47.6 dB (1–4 Hz) to 7.1 dB (400–600 Hz), and little-endian was confirmed against
big-endian by a 25-fold difference in s.d. The data was fine. The analysis was indexing
`data[channel]` — a row, one time point across 134 channels — instead of
`data[:, channel]`, so it was band-passing the probe's spatial profile rather than time.

The reason this is a decision and not just a fixed bug is that **it produced no error and
plausible numbers**. 134 samples is enough for a 150–300 Hz filter to return something,
and what it returns is a few tens of units, slightly different per channel — a flat
profile reading as "all channels are equivalent" rather than as a fault. The 1–4 Hz band
is the only part that complained, and it complained with NaN rather than an exception,
because a narrow low-frequency Butterworth in transfer-function form is numerically
unstable at this length.

Three things in `src/hc11/ripple_channels.py` now make this hard to repeat: every filter
call passes `axis=0` explicitly rather than taking scipy's default of `-1`;
`channel_profile` loops one channel at a time over 1-D arrays, where the axis cannot be
ambiguous; and the band-pass is `sosfiltfilt` on second-order sections rather than
`filtfilt` on a transfer function, so the delta band is numerically stable. Two tests
pin the behaviour by *answer* rather than by arithmetic —
`test_bandpass_filters_time_not_channels` and
`test_profile_finds_the_channel_carrying_ripples`, the latter planting a 190 Hz burst on
one known channel — and a third, `test_profile_is_not_flat_across_channels`, pins the
symptom directly.

---

## Step 6b — Ripple detection

### D6.16 — One envelope per channel for the whole recording, then sliced by state

`channel_envelope` band-passes and Hilbert-transforms a channel across the
entire `.eeg` in chunks, and detection slices that envelope by state interval.
The alternative — filter each state interval separately — puts a filter
transient at the edge of every interval, and non-REM intervals here are often
under a minute, so those edges would be a measurable fraction of the searched
time.

Chunking is exact where it matters. Internal chunk boundaries agree with a
whole-signal transform to 3 × 10⁻⁵ of the envelope maximum; the only deviation
is in the first and last second of the *recording*, where the FFT's periodicity
assumption produces a wrap-around artefact whose size depends on the length
transformed. That is a property of the signal's ends, not of chunking, and it
is pinned in two tests — one for the interior, one placing a chunk boundary
inside a planted ripple.

Thresholds are still per state, as in the MUA detector and for the same reason:
baseline ripple power differs between non-REM sleep and quiet waking, so one
session-wide threshold would hand one epoch systematically more events.

### D6.17 — Primary channel: level ranking, against the prior expectation

The brief's instruction was to prefer the excursion ranking (99th percentile
over median) if the two rankings disagreed materially on event quality, on the
reasoning that mean envelope level partly measures channel noise while
excursion measures transients against the channel's own background. That
reasoning is sound a priori. The data does not support acting on it here, so
the config sets `primary_ranking: level`.

**The two rankings barely differ where it was expected to matter.** The best
channels under each ranking are adjacent contacts (ch 2 vs 1, ch 3 vs 4), and
their event sets match 94–98 % with a median peak offset of 0.0 ms. Whichever
is chosen, the same events come out.

**Where they differ by 140 µm, the level ranking is the better one.** Spike
group 3 of Achilles_10252013 is the one shank where the rankings land seven
contacts apart: ch 20 (level) against ch 27 (excursion).

| | ch 20 (level) | ch 27 (excursion) |
|---|---|---|
| PRE events | 4579 | 4144 |
| POST events | 3166 | 2969 |
| events found only here (PRE) | 780 | 271 |
| median other-shank support of those | 7 / 11 | 6 / 11 |
| median other-shank support of shared events | 11 / 11 | 11 / 11 |
| median peak z of only-here events | 3.41 | 3.41 |
| median peak z, all events | 5.69 | 6.17 |

Read together: both channels' unique events are weak ones, sitting just above
the 3 SD threshold, and **both sets are real** — a median of 6–7 other CA1
shanks show the same transient, against 11/11 for the events the two channels
share. So neither ranking is manufacturing noise. The level-ranked channel
simply finds more of the same thing, and its extra events are if anything
better corroborated.

`ripple_disagreement_shank3.png` shows why, and it is the opposite of the
prior worry. Around events only ch 20 detects, ch 27's band-passed trace is
oscillating more or less continuously rather than in discrete bursts: it has a
high, near-constant 150–300 Hz background, which is what a contact away from
*stratum pyramidale* picks up. A discrete ripple does not stand out above that
background, and dividing by a median does not rescue a channel whose background
is itself spiky.

This is recorded as a decision rather than a result because the prior was
reasonable and the evidence went the other way. Both rankings are detected on
and both are in `results/ripple_channel_comparison.csv`; `primary_ranking` is
one config line, so the sensitivity is a re-run rather than a rewrite.

### D6.18 — Multi-shank support, measured per event and not filtered on

For every event on the primary channel, `shank_support` reports the largest
ripple-band envelope z-score within ±20 ms on each of the other 11 CA1 shanks'
chosen channels, against that channel's own baseline in the same state. A
shank counts as supporting at z ≥ 2. Nothing is excluded on that basis — the
count travels with the event in `ripple_events.csv`.

The distributions are strongly bimodal towards full agreement. In
Achilles_10252013 POST, 2366 of 3284 events are seen on all 11 other shanks and
only 13 on none; PRE and the other session are the same shape. Median support
is 11/11 in five of the six session × epoch combinations and 10/11 in the
sixth.

Two implementation notes, both of which were bugs first:

- The comparison channels' baselines come from a sample of the same state, not
  a full-recording pass. Eleven extra passes cost about four minutes per
  session and a baseline is a summary statistic.
- When no interval is long enough to sample — MAZE immobility is fragments of
  under a second — the state itself is used instead, capped at the same total
  duration. Without that fallback the baseline is NaN, every z is NaN, and
  **"no other shank saw this event" is indistinguishable from "this was never
  measured"**. The first run reported exactly that: all 18 MAZE events with
  zero support. With the fallback they have a median of 10/11.

### D6.19 — Cross-check: the LFP and MUA detectors agree, and agree about which events are weak

| Session | Epoch | ripples | matched by MUA | MUA events | matched by a ripple | offset (MUA − ripple) |
|---|---|---|---|---|---|---|
| 10252013 | PRE | 4761 | 80.1 % | 7046 | 51.1 % | −4.3 ms (IQR 17.2) |
| 10252013 | MAZE | 18 | 72.2 % | 20 | 60.0 % | −5.9 ms (IQR 30.7) |
| 10252013 | POST | 3284 | 76.9 % | 3955 | 58.4 % | −5.1 ms (IQR 14.6) |
| 11012013 | PRE | 3324 | 69.5 % | 4709 | 47.0 % | −4.7 ms (IQR 16.0) |
| 11012013 | MAZE | 25 | 56.0 % | 20 | 65.0 % | −5.5 ms (IQR 25.0) |
| 11012013 | POST | 2869 | 72.0 % | 3328 | 57.9 % | −5.9 ms (IQR 15.1) |

The offset is consistent across every epoch and both sessions: population
firing peaks about 5 ms **before** the ripple envelope peak, with an
interquartile range of 15–17 ms. Part of that is physiology and part may be
ours — the envelope is smoothed with an 8 ms Gaussian before thresholding,
which can shift its peak. The sign and the tightness are the result; the exact
5 ms should not be quoted without re-running at a smaller smoothing kernel.

The more useful result is the cross-tabulation the brief asked for. Agreement
with the MUA detector rises monotonically with multi-shank support:

| other shanks showing the event | matched by MUA, 10252013 PRE | POST |
|---|---|---|
| 0 | 37 % | 38 % |
| ≤ 2 | 53 % | 52 % |
| ≥ 9 | 85 % | 80 % |
| 11 | 87 % | 82 % |

The two checks are independent — one is the field at other recording sites, the
other is spiking at the same site — and they pick out the same events as weak.
That is the strongest evidence so far that the low-support events are a real
category and not an artefact of either method.

### D6.20 — Recommendation: MUA primary, ripples as corroboration

For the Fig. 1C candidate-event set, the MUA detector should be primary and
ripple matching should ride along as a per-event property. Three reasons, in
order of weight:

1. **Coverage.** LFP exists locally for 2 of the 8 sessions. Each one is a
   9.3–9.7 GB download and about 12 GB extracted. An LFP-primary pipeline
   cannot run on the full dataset without another ~76 GB; an MUA-primary one
   runs today on all eight, with the LFP available as a check on two.
2. **The sequence analysis is a spike analysis.** Scores are computed from the
   decoded posterior over spikes. A candidate-event definition drawn from the
   same signal keeps the detection and the measurement on one instrument.
3. **Ripples are close to a subset.** 70–80 % of ripples have an MUA partner,
   against 47–58 % the other way, so taking MUA as primary discards less.

Chen et al. use the conjunction — ripple-band power *combined with* multi-unit
activity above mean + 3 SD — and that remains available as the strict variant,
since every ripple event now carries its MUA match and every MUA event can be
asked for a ripple. The three candidate sets (MUA, LFP, conjunction) differ by
enough to matter and the Fig. 1C distributions should be reported for the
primary with at least the conjunction as a sensitivity row.

What this does **not** settle: ~45 % of MUA events have no ripple. They may be
lower-synchrony population bursts, events on a shank the primary channel does
not cover, or detection noise. That is worth a look in 6d, where the
per-event participation numbers are already to hand.

### D6.21 — MAZE events are few, but they are real ripples

MAZE yields 18 and 25 events over 0.4 min of scored immobility in the two
sessions — the immobility problem of D6.3, unchanged. What is new is that these
events are not junk: median multi-shank support is 10/11 and 11/11, and 56–72 %
have an MUA partner, both comparable to sleep. So the shortage is a shortage of
*detected immobility*, not of awake sharp-wave ripples, and it supports the
open hypothesis in D6.3 rather than the alternative that the tracking-lost
periods are artefact. Resolving it still needs EMG or a theta/delta ratio, both
of which are now reachable: the EMG channels are identified (0-based 128–130)
and the band machinery is in place.

---

## Step 6d — Characterising the candidate events

### D6.22 — The MUA-to-ripple lead survives the smoothing sweep; it is not an artefact

Step 6c reported population firing peaking about 5 ms before the ripple
envelope peak, with the caveat that the envelope is smoothed with an 8 ms
Gaussian and that smoothing can move a peak. `scripts/offset_vs_smoothing.py`
sweeps the kernel from 0 to 32 ms. Median offset in ms (MUA peak minus ripple
peak), events held fixed and peaks re-located:

| | 0 ms | 1 | 2 | 4 | 8 | 16 | 32 |
|---|---|---|---|---|---|---|---|
| 10252013 PRE | −4.5 | −4.3 | −4.3 | −3.9 | −4.3 | −4.7 | −5.1 |
| 10252013 POST | −4.9 | −4.7 | −4.5 | −4.9 | −5.1 | −5.3 | −5.5 |
| 11012013 PRE | −4.9 | −4.9 | −5.0 | −4.7 | −4.7 | −5.3 | −5.9 |
| 11012013 POST | −5.5 | −5.5 | −5.5 | −5.5 | −5.9 | −6.3 | −6.3 |

Re-detecting at each width rather than holding the event set fixed gives the
same picture (both are in `results/offset_vs_smoothing.csv`).

**The offset is already there at zero smoothing.** Across the range that
matters — 0 to 8 ms — it moves by at most 0.6 ms, which is less than one
sample at 1250 Hz. Only at 32 ms, four times the configured kernel, does it
grow by a consistent ~0.7 ms, in the direction a kernel would push an
asymmetric peak. So the kernel contributes a fraction of a millisecond to a
5 ms effect.

The lead is therefore reported as a result: population firing peaks **4.5–6 ms
before** the ripple-band envelope peak, in both sessions and both sleep epochs.
Two limits on how hard to lean on it. The measurement granularity is 0.8 ms
(one LFP sample), so the quantity is good to about ±1 ms, not better. And a
distribution with an interquartile range of 15–17 ms has a median that is well
determined but describes a broad relationship, not a fixed latency. With the
question settled, the offset now appears in `cross_check_offsets.png`.

### D6.23 — Why ~45 % of MUA events have no ripple: they are smaller, not elsewhere

Three explanations were tested against the ~45 % of MUA events with no ripple
partner. The table is Achilles_10252013 PRE; POST and the second session give
the same pattern and are in `results/unmatched_mua_hypotheses.csv`. "rb" is the
rank-biserial correlation, positive meaning the unmatched events score higher.

| quantity | matched | unmatched | rb |
|---|---|---|---|
| active pyramidal cells (mean) | 25.3 | 15.3 | **−0.61** |
| active place cells (median) | 15 | 10 | −0.58 |
| duration (median, ms) | 82 | 57 | −0.48 |
| peak z (median) | 4.62 | 3.78 | −0.39 |
| shanks with an active cell (median) | 9 | 8 | −0.51 |
| busiest-shank share (median) | 0.200 | 0.241 | +0.38 |
| busiest-shank share, **excess over the cell-count null** | −0.007 | −0.007 | **+0.03** |
| share on the detection channel's shank (mean) | 0.124 | 0.138 | +0.09 |
| passes the spike-jitter null (mean) | 0.918 | 0.867 | **−0.05** |

**1. Low-synchrony bursts — supported, and this is the whole story.** Unmatched
events have about 60 % of the active cells of matched ones (15.3 against 25.3;
10.0 against 16.3 in the other session), with a rank-biserial of −0.58 to
−0.63. They are also shorter and weaker. This is a large effect by any
standard and it is consistent across both sessions and both sleep epochs.

**2. Events on shanks the detection channel does not cover — rejected.** The
raw busiest-shank share *is* higher for unmatched events, which looks like
support until the confound is removed: an event with nine active cells cannot
occupy twelve shanks, so fewer cells forces higher concentration whatever the
physiology. Compared against a per-event null — the same number of cells drawn
at random, weighted by how many cells each shank actually has — the excess is
−0.007 for both groups, rank-biserial +0.03. The prediction was that unmatched
events would favour shanks *away* from the detection channel; the observed
share on that shank is 0.138 against 0.124, i.e. marginally *towards* it, and
both sit at the 0.117 expected from cell counts alone. In Achilles_11012013
PRE the direction even reverses (rb −0.018, p = 0.28).
`unmatched_shank_profile.png` shows the two series tracking each other and the
cell-count expectation across all twelve shanks.

**3. Detection noise — a small contribution, not the explanation.** Unmatched
events fail the spike-jitter null more often, but only by 5–9 percentage
points: 86.7 % of them still beat their own null against 91.8 % of matched
events (80.0 % against 88.7 % in the second session). So the great majority of
unmatched events are genuine co-firing, just less of it.

**Conclusion.** The unmatched MUA events are the weak tail of the same
distribution, not a separate population and not a blind spot in the LFP
coverage. That agrees with the 6c result that multi-shank support predicts
MUA agreement: both detectors are measuring the same underlying strength and
disagreeing about where to put the threshold on it.

Two consequences. The LFP is not adding coverage of a part of the recording the
spikes miss, which removes the main argument for an LFP-primary or
LFP-supplemented event set. And the conjunction is best understood as a
*stricter threshold* on event strength rather than as a different kind of
filter — which is what makes its cost in D6.24 the relevant number.

### D6.24 — What the strict variant costs in usable events

"Usable" means at least 5 active place cells and at least 10 % of the place-cell
ensemble, the `events.decodable` criterion from step 5.

| Session | Epoch | MUA events | usable | conjunction events | usable | strict keeps |
|---|---|---|---|---|---|---|
| 10252013 | PRE | 7046 | 5827 (83 %) | 3601 | 3433 (95 %) | **59 %** |
| 10252013 | POST | 3955 | 3368 (85 %) | 2309 | 2201 (95 %) | **65 %** |
| 10252013 | MAZE | 20 | 17 (85 %) | 12 | 12 (100 %) | 71 % |
| 11012013 | PRE | 4709 | 3198 (68 %) | 2213 | 1924 (87 %) | **60 %** |
| 11012013 | POST | 3328 | 2465 (74 %) | 1927 | 1717 (89 %) | **70 %** |
| 11012013 | MAZE | 20 | 16 (80 %) | 13 | 11 (85 %) | 69 % |

The strict variant costs **30–41 % of the usable events** and buys a rise in the
decodable fraction from 68–85 % to 87–95 %. Given D6.23 — the events it drops
are weak, not different in kind — that is a straightforward sensitivity knob
rather than a correctness question, and it is why the primary stays the full
MUA set with the conjunction reported alongside.

Both counts stay comfortably above what Fig. 1C needs; even the smallest cell,
Achilles_11012013 POST under the conjunction, is 1717 events. The real cost is
in the tails of the score distribution, where the strict variant will have
roughly two thirds of the events to resolve them with.

### D6.25 — `results/candidate_events.csv` is the handover to step 7

One row per MUA event, carrying `ripple_matched`, `ripple_peak_offset_s`,
`n_supporting_shanks`, `n_active_pyr`, `n_active_place` and `decodable`. The
conjunction is therefore `ripple_matched == True` applied to a column, not a
second pipeline, and the sensitivity variants step 7 needs — full set, strict
conjunction, decodable-only, support thresholds — are all filters on this one
file. 19,078 rows across both sessions and all three epochs.

`n_supporting_shanks` is NaN for an event with no ripple partner, which is the
honest value: support is measured on the ripple, and an event without one was
never measured. It must not be read as zero.

---

## Step 7 — Sequence scores and the Fig. 1C comparison

### D7.1 — The score is a weighted correlation over the whole posterior

No position is read off the posterior. The decoded distribution is used as a
weight on each `(position, time)` pair and the ordinary Pearson correlation is
taken under those weights (formula in the `hc11.sequences` docstring). A bin
whose posterior is spread over half the track contributes its spread rather
than a point estimate, which matters because a maximum-a-posteriori trajectory
through broad posteriors is mostly noise.

Measured during events, a 20 ms bin carries **5–6 spikes** from the place-cell
ensemble, not the ~1 spike step 5 found during running. That is the difference
between a 15 cm/s animal and a population burst, and it means the posteriors
are sharper than expected. Fewer than 1 % of bins are empty.

### D7.2 — Maximising over field direction, and paying for it in the null

On a linear track the fields are directional, so an event is scored against
each direction's templates and the larger `|r|` kept. This inflates the raw
correlation and the inflation is not small.

It is harmless for `rZ` **only** because the same maximisation runs inside the
null: one shuffle realisation is applied to both directions' posteriors and the
same maximum taken. Shuffling the directions independently would give the null
more freedom than the observation had and would push every score up.
`test_the_null_takes_the_same_maximum_over_directions` pins this. Raw `r` is
reported but must not be compared against a nominal significance level.

### D7.3 — A direction needs five laps before its fields are used

Achilles_11012013 is a circular maze run one way. Its "pos" direction holds
**exactly one lap of 226 samples**, against 75 laps and 11,730 samples for
"neg". A template built from one lap is noise.

A noise template is worse here than no template, and the reason is D7.2: the
score maximises over directions, so the shuffles get the same extra freedom,
the null rises, and *every* event's normalised score falls. Removing it moved
that session's median `rZ` from −0.060 to +0.110 in PRE and from +0.148 to
+0.307 in POST — a change of the same size as the PRE/POST effect being
measured, from a template nobody would have defended if asked.

Floor set at five laps. Achilles_10252013 has 40 and 42.

### D7.4 — Circular position: maximise over the origin

On the circular maze, position has no origin, and a trajectory crossing the
seam of the linearisation scores near zero despite being a perfect sequence
(`test_a_wrapping_trajectory_scores_near_zero_on_a_linear_axis`: 0.75 linear
against 0.98 circular on the same posterior). The correlation is therefore
evaluated at every circular offset of the position axis and the best taken —
"is there an origin that makes this trajectory linear in time", which is the
question that has meaning on a ring. The maximisation runs inside the shuffles
too, for the same reason as D7.2. Linear tracks take the same code path with a
single offset.

### D7.5 — The five-bin floor is the largest attrition step in the pipeline

Bins with no spikes are dropped before scoring: the memoryless Poisson decoder
returns `-dt * sum_i f_i(x)` for an empty bin, which describes where the
ensemble is quiet, not where the animal is. Chen et al. require more than five
bins, so the floor is five *kept* bins.

At 20 ms that needs a 100 ms event, and our events are shorter than that:

| | median duration | ≥3 bins | ≥4 | **≥5** | ≥6 |
|---|---|---|---|---|---|
| 10252013 PRE | 68 ms | 62 % | 37 % | **24 %** | 17 % |
| 10252013 POST | 62 ms | 54 % | 32 % | **21 %** | 16 % |
| 11012013 PRE | 57 ms | 46 % | 25 % | **15 %** | 10 % |
| 11012013 POST | 58 ms | 47 % | 28 % | **19 %** | 13 % |

So 15–24 % of candidate events are scored, and this is a larger loss than every
other filter in the pipeline combined. It is stated rather than tuned away:
lowering the floor to three bins would mean estimating a correlation from three
time points, which is not an estimate. The real cause is upstream — the MUA
detector's 1 SD boundaries give a median event of 55–68 ms, where the replay
literature usually works with events twice that — and it is a candidate for
revisiting if step 7's effects turn out to be marginal.

The surviving events are the longer ones, which also carry more cells. The
PRE/POST comparison is unaffected, since the same filter applies to both, but
the absolute rate of significant events is a rate within a selected subset.

### D7.6 — The two shuffles answer different questions, and the difference is measurable

**Column-cycle** rolls each time bin's posterior independently. It destroys the
position-time relationship and keeps each bin's posterior shape — but it also
spreads the posterior mass uniformly across position, which raises the variance
of position in the null. Any event whose decoded position is merely *localised*
across bins therefore scores above this null, sequence or no sequence.

**Time-bin permutation** reorders the bins. It preserves the set of posteriors
exactly, so the position variance is unchanged and only the covariance moves.
It isolates sequence from localisation.

The identity-shuffle control (D7.7) shows this is not a theoretical worry. With
place fields detached from cells, median `rZ` is still **+0.59** under
column-cycle in Achilles_10252013, because the same cells firing across
neighbouring bins still produce a localised posterior. Under time-permutation
it is +0.25. Neither is zero.

**Consequence: the absolute level of `rZ` is not interpretable, and no claim
should rest on "x % of events exceed z = 1.65".** The PRE/POST *difference* is
what survives the control, and that is what the Fig. 1C comparison reports.

### D7.7 — Both controls pass

**Negative control: identity-shuffled place fields.** Permuting which cell owns
which field, leaving the events, bins, spikes and shuffles identical. Median
`rZ`, POST minus PRE:

| Session | shuffle | real fields | identity-shuffled |
|---|---|---|---|
| 10252013 | column cycle | **+0.302** | −0.005 |
| 10252013 | time permutation | **+0.220** | −0.033 |
| 11012013 | column cycle | **+0.197** | +0.049 |
| 11012013 | time permutation | **+0.433** | +0.071 |

The level survives identity shuffling (D7.6); the PRE/POST gap does not. Three
of four shuffled gaps are under a sixth of the real one, and one is zero.

**Positive control: MAZE running.** Real trajectories, decoded the same way with
leave-one-lap-out templates as in step 5b:

| Session | window | shuffle | n | median rZ | above 1.65 |
|---|---|---|---|---|---|
| 10252013 | 2.0 s | column cycle | 62 | **+13.23** | 95 % |
| 10252013 | 2.0 s | time permutation | 62 | **+10.68** | 94 % |
| 11012013 | 2.0 s | column cycle | 145 | **+7.50** | 82 % |
| 11012013 | 2.0 s | time permutation | 145 | **+7.91** | 85 % |
| 10252013 | 0.2 s | column cycle | 984 | −0.10 | 11 % |
| 11012013 | 0.2 s | column cycle | 1606 | −0.67 | 2 % |

The scoring works: a known trajectory scores an order of magnitude above
anything in sleep.

The 0.2 s rows are not a failure and are the reason both lengths are reported.
At roughly 30 cm/s, 200 ms of running covers about 6 cm — **under one 10 cm
position bin**. Such a window cannot contain a trajectory however good the
decoder is, so its score is uninformative by construction. That sleep events of
that duration score at all is the point: replay is time-compressed, so 60 ms of
replay traverses what takes the animal a second or more. It also sets the
sensitivity ceiling: the method can only see a sequence when the decoded
position actually moves across bins.

### D7.8 — The result: POST exceeds PRE, in every cell of the design

Median `rZ`, PRE against POST, with the rank-biserial effect size and the
Mann-Whitney p. "matched-n" repeats the test 200 times with PRE cut at random
to POST's event count.

| Session | Event set | Shuffle | PRE | POST | diff | rb | p | matched-n diff |
|---|---|---|---|---|---|---|---|---|
| 10252013 | all | column cycle | +0.548 | +0.850 | **+0.302** | +0.114 | 3.4e−06 | +0.302 (100 % of draws p<0.05) |
| 10252013 | all | time permutation | +0.188 | +0.407 | **+0.220** | +0.131 | 1.2e−07 | +0.217 (100 %) |
| 10252013 | conjunction | column cycle | +0.529 | +0.873 | **+0.344** | +0.136 | 5.4e−07 | +0.336 (100 %) |
| 10252013 | conjunction | time permutation | +0.188 | +0.438 | **+0.250** | +0.147 | 6.0e−08 | +0.248 (100 %) |
| 11012013 | all | column cycle | +0.110 | +0.307 | **+0.197** | +0.071 | 0.026 | +0.197 (77 %) |
| 11012013 | all | time permutation | +0.291 | +0.723 | **+0.433** | +0.178 | 3.1e−08 | +0.433 (100 %) |
| 11012013 | conjunction | column cycle | +0.129 | +0.350 | **+0.221** | +0.091 | 0.0091 | +0.221 (100 %) |
| 11012013 | conjunction | time permutation | +0.251 | +0.776 | **+0.525** | +0.189 | 6.4e−08 | +0.516 (100 %) |

Eight of eight positive. The subsampled tests return the same difference to
three decimals, which is what one expects when the difference is in the
distributions rather than in the sample sizes.

Per minute of searched state, rather than per event:

| Session | shuffle | PRE | POST |
|---|---|---|---|
| 10252013 | column cycle | 2.01 /min | 2.69 /min |
| 10252013 | time permutation | 0.83 /min | 1.58 /min |
| 11012013 | column cycle | 0.45 /min | 0.81 /min |
| 11012013 | time permutation | 0.70 /min | 1.72 /min |

POST wins on a rate basis by a wider relative margin than on medians, which is
the stronger form of the claim: POST has roughly half the non-REM of PRE and
still produces more significant events per minute.

**How much to make of it.** The effect is consistent — two sessions, two
shuffles, two event sets, all positive, with the subsampling and the identity
control both closing the obvious objections — but it is *small*. A rank-biserial
of +0.07 to +0.19 means that drawing one PRE and one POST event at random, POST
scores higher about 54–59 % of the time. The direction replicates; the
magnitude is modest, and D7.5's 15–24 % scoring rate and D7.6's uninterpretable
absolute level are the two things that would have to be fixed before the
magnitude could be compared against the paper's.

Forward and reverse traversals split about evenly (49–52 % forward in every
epoch of both sessions) and score similarly, with no consistent asymmetry
between them.

The conjunction set gives a *larger* effect than the full set in all four
comparisons, which fits D6.23: it is a stricter threshold on event strength, and
stronger events carry more sequence.

---

## Step 8a — audit before scaling to eight sessions

### D8.1 — Which parameters are session-general, and which are an Achilles-only assumption

Everything in `config/params.yaml` has been *applied* to whatever sessions its stage ran
on. The question here is narrower: for each parameter, has its value been **checked**
against more than the two Achilles sessions?

| Parameter | Value | Provenance | Checked on | Verdict |
|---|---|---|---|---|
| `place_fields.bin_size_cm` | 10 cm | [CHEN2016] | 8 | session-general |
| `criteria.max_field_width_frac` | 0.75 | [PROVISIONAL] | 8 (`place_cell_sensitivity.csv`) | session-general **by construction** — it replaced a fixed 120 cm precisely because that was 75 % of one maze and 41 % of another (D5.2) |
| `criteria.min_spikes` / `min_peak_rate_hz` / `alpha` / `min_stability_r` | 50 / 1 Hz / 0.05 / 0.3 | [PROVISIONAL] | 8 | applied to 8, but the *yield* ranges 29–79 % of pyramidal cells, so the criteria are uniform and their consequences are not |
| `behavior.speed_threshold_cm_s` | 15 | [CHEN2016] | 8, with a sensitivity run on Achilles_10252013 **and Cicero_09012014** | session-general; the one sensitivity check that already includes a non-Achilles session |
| `sequences.min_laps_per_direction` | 5 | [PROVISIONAL] | 2 | **not Achilles-specific — it is a circular/linear discriminator** (below) |
| `events.mua.threshold_high_sd` | 3.0 | [CHEN2016], but stated there for the *ripple envelope*, not for a population rate | 2 | **untested on 6** |
| `events.mua.threshold_low_sd`, `min/max_duration_ms`, `merge_gap_ms` | 1.0, 30/400 ms, 30 ms | [PROVISIONAL] | 2 | **untested on 6** |
| `events.ripples.*` | band [CHEN2016], rest [PROVISIONAL] | 2 | **untested on 6** |
| `lfp.channel_selection.*` | 40 × 15 s NREM, bands | [PROVISIONAL] | 2 | **untested on 6, and the probe geometry differs** (D8.2) |
| `sequences.min_bins` | 5 | [CHEN2016] | 2 | **the riskiest parameter in the config** (D8.3) |
| `events.decodable` | ≥ 5 cells and ≥ 10 % of the ensemble | [CHEN2016] | 2 | binds far harder on a 16-cell ensemble than on a 77-cell one (D8.3) |

**The five-lap floor is not an Achilles artefact.** Laps per direction across all eight
sessions: every linear session has 24–44 in both directions; every circular session has 0
or 1 in the direction the animal did not run (Achilles_11012013: 1 and 75;
Cicero_09102014: 0 and 19; Gatsby_08282013: 1 and 79). The floor therefore separates
circular from linear with a margin of 23 laps on one side and 4 on the other. It was
*discovered* on the one circular session we happened to run, but it generalises, and it
will fire on exactly the two circular sessions still to come.

**Nothing hard-codes a channel count or a session name in `src/hc11/`** (checked): the
only occurrences are docstrings describing Achilles as an example. Eleven scripts default
to a two-session list, which 8b must widen. That is a default, not a method.

### D8.2 — There are two probe geometries in this dataset, not one

From `Channel_Orderings.pdf`, read for all four animals rather than only the one we had
processed:

| Animal | Probe | CA1 shanks | Channels/shank | Extra-channel groups | Misc. channels | Total | Skipped (1-based) |
|---|---|---|---|---|---|---|---|
| Achilles | BUZ64SP 6-shank × 2 | 12 | 10 | 2 × 4 | EMG 129–131, accel 132–134 | 134 | none |
| Cicero | BUZ64SP 6-shank × 2 | 12 | 10 | 2 × 4 | EMG 129–131, accel 132–134 | 134 | none |
| Buddy | BUZ64 8-shank × 2 | **16** | **8** | **none** | **accel 129–131, EMG 132** | **132** | **24, 27, 58** |
| Gatsby | BUZ64 8-shank × 2 | **16** | **8** | **none** | EMG 129–131, accel 132–134 | 134 | **35, 45, 47** (session 2 only) |

Four consequences, none of which requires a methodological change but all of which
change what the step-6 output will look like:

1. **The channel count is not 134 for every session.** Buddy is 132. The code takes it
   from the XML and `require_consistent_size` would catch any mismatch, so this is a
   prediction to verify rather than a risk.
2. **The "extra channel" exclusion has nothing to exclude on Buddy and Gatsby.** Our
   filter keeps spike groups that carry sorted units (D6.13), which is data-driven and
   adapts; on Buddy and Gatsby all 16 groups are genuine CA1 shanks and the filter will
   simply drop the ones without units.
3. **The multi-shank support denominator varies by session.** Shanks carrying units:
   Achilles 12 and 12, Cicero 11/12/9, Gatsby 8 and 6, Buddy 10. Support was reported as
   "of 11 other CA1 shanks" on Achilles; on Gatsby_08282013 it will be "of 5". The
   *distribution* is therefore not comparable across sessions without normalising to a
   fraction, and 8b should report the fraction alongside the count.
4. **The XML skip flag has never run on real data.** Achilles has no skipped channels, so
   `parse_xml`'s skip handling and `good_channels()` have only ever been exercised by
   unit tests. Buddy and Gatsby_08282013 will exercise them for the first time. If the
   skipped channels do not appear in the XML with `skip="1"` the code will silently
   include three bad channels in the selection pool — worth checking explicitly in 8b
   against the three numbers above rather than trusting it.

### D8.3 — The two validated sessions are the two best-case sessions in the dataset

This is the finding that matters most before committing to 64 GB of download. Achilles
has the two largest place-cell ensembles of the eight. From `results/decoding_summary.csv`,
at the **20 ms bin the sequence scoring uses**, measured during running:

| Session | Place cells | Median spikes per 20 ms bin | Empty bins | Median error (20 ms) |
|---|---|---|---|---|
| **Achilles_10252013** | **77** | **3–4** | **5–9 %** | 10.0–10.4 cm |
| **Achilles_11012013** | **73** | **2** | **13 %** | 10.6 cm |
| Cicero_09102014 | 36 | 1 | 38 % | 12.2 cm |
| Gatsby_08282013 | 31 | 1 | 42 % | 17.8 cm |
| Gatsby_08022013 | 30 | 1 | 27–36 % | 14.2–16.0 cm |
| Cicero_09172014 | 25 | 1 | 42–43 % | 19.2–21.2 cm |
| Buddy_06272013 | 16 | 0–1 | **47–51 %** | 18.0–18.6 cm |
| Cicero_09012014 | 16 | 0 | **58–60 %** | 19.2–19.9 cm |

**What this predicts for step 7 on the six unseen sessions.** Spike density inside events
runs 1.6–2.4× the running density (Achilles: 5.55 and 4.88 in-event against 3–4 and 2
running), and under 1 % of in-event bins were empty on Achilles. Scaling that ratio, the
six unseen sessions should see roughly 0–2.5 spikes per in-event bin and an empty-bin
fraction in the tens of percent rather than under one.

The five-bin floor counts **non-empty** bins. On Achilles it was effectively "≥ 5 raw
bins", i.e. a 100 ms event, and it already discarded 76–85 % of candidates. On a session
where 30 % of in-event bins are empty, an event needs about 7 raw bins — 140 ms — to
yield 5 scorable ones, and the ≥ 7-bin fraction on Achilles was 10–17 %. On Buddy and
Cicero_09012014, with 16 place cells each, the scorable fraction could be low single
digits.

Separately, `events.decodable` requires ≥ 5 active place cells. On a 16-cell ensemble
that is 31 % of the ensemble active in a single event; on Achilles's 77 it is 6 %. The
same criterion is a far stronger filter on the small sessions.

**This is a prediction, not a result, and it is not a reason to change anything.** It is a
reason to expect that two to four of the six new sessions may return too few scorable
events for a PRE/POST comparison, and to plan to report that as the outcome rather than
treat it as a failure. If a session returns fewer than roughly 100 scorable events per
epoch, its individual PRE/POST test will be underpowered and should be reported with its
n rather than as a negative result.

### D8.4 — The circular sessions belong in a separate column, and the reason is the score, not the decoder

Of the six sessions still to come, **two are circular** (Cicero_09102014,
Gatsby_08282013), not three; the third circular session, Achilles_11012013, is one of the
two already analysed. The final split is five linear and three circular.

**The decoder is fine on circular tracks.** At 250 ms the circular sessions decode to
4.9, 4.9 and 5.4 cm median error against chance baselines of 67–72 cm — as good as the
linear sessions. The approximate linearisation of D5.5 (r = 0.972–0.997 against the
authors' `OneDLocation`, where our linear reconstruction is r = 1.0000) does **not** enter
any primary result: every primary analysis uses the authors' own `OneDLocation`, and our
reconstruction was built only for the measurement-base diagnostic. The unresolved extent
mismatch (D3.3: `OneDLocation` spans 2.84–2.90 m where a 1 m-diameter platform implies
3.14 m) shifts the absolute scale of position but not the decoder's accuracy in bins,
because `Track` takes the extent from the observed range per session.

**The sequence score is not on the same scale between geometries, and that is a
methodological consequence of D7.4.** On a ring the weighted correlation is maximised over
all circular offsets of the position axis; on a linear track there is one offset. The null
gets the same freedom, so `rZ` stays calibrated — but a null maximised over 28–29 offsets
is systematically broader than one maximised over a single offset, which compresses `rZ`
for the same quality of sequence. The two sessions we have already show the signature:

| | column cycle PRE → POST | time permutation PRE → POST |
|---|---|---|
| Achilles_10252013 (linear) | +0.548 → +0.850 | +0.188 → +0.407 |
| Achilles_11012013 (circular) | +0.110 → +0.307 | +0.291 → +0.723 |

The circular session's column-cycle medians are roughly a fifth of the linear session's
while its time-permutation medians are nearly double. That is a scale difference between
geometries interacting with the two shuffles, not a difference in biology, and it means
**absolute `rZ` must not be pooled across geometries**.

Recommendation for 8c, to be applied regardless of how the numbers come out:

- Report linear and circular sessions in separate columns throughout.
- Pool the **direction** of the PRE/POST effect — how many of eight sessions show
  POST > PRE — because a sign is scale-free and that count is the headline the brief
  asks for.
- Pool **magnitudes** only within a geometry, and say so wherever a pooled median `rZ`
  appears. A pooled median across five linear and three circular sessions is a mixture of
  two scales and should not be quoted as one number.

### D8.5 — Runtime and disk for the eight-session run

Expected `.eeg` size is `sessDuration × 1250 Hz × n_channels × 2 bytes`. The model
reproduces both known files to the byte (11.68 and 12.27 GB), and `sessDuration` equals
the POST epoch end in all eight sessions, so there is no truncation ambiguity in the
prediction — including for Gatsby_08022013, whose anomaly is spikes *past* `sessDuration`
(D3.1), not a short recording.

| Session | Channels | Predicted `.eeg` | Archive (≈ 0.79 ×) |
|---|---|---|---|
| Buddy_06272013 | 132 | 6.95 GB | ~5.5 GB |
| Cicero_09012014 | 134 | 11.64 GB | ~9.2 GB |
| Cicero_09102014 | 134 | 12.23 GB | ~9.7 GB |
| Cicero_09172014 | 134 | 11.24 GB | ~8.9 GB |
| Gatsby_08022013 | 134 | 10.19 GB | ~8.1 GB |
| Gatsby_08282013 | 134 | 11.84 GB | ~9.4 GB |
| **six-session total** | | **64.1 GB** | **~51 GB** |

Disk: 197 GB free on the data volume, 45 GB already in use by the two Achilles sessions.
Peak requirement is ~115 GB if the archives are kept alongside the extracted files, which
fits with ~80 GB to spare. Deleting each archive after extraction would hold the peak at
~64 GB.

Compute, scaled from measured Achilles timings:

| Stage | Per session | Six sessions |
|---|---|---|
| extraction (write-bound, ~70 MB/s) | 100–175 s | ~13 min |
| channel selection | ~40 s | ~4 min |
| MUA detection (recomputes place cells) | ~35 s | ~4 min |
| spike-jitter null (5 repeats) | ~60 s | ~6 min |
| ripple detection + multi-shank support | ~80 s | ~8 min |
| cross-check, characterisation | ~35 s | ~4 min |
| sequence scoring (+ identity-shuffle control) | ~20 s | ~2 min |
| running positive control | ~60 s | ~6 min |
| **total** | **~7–8 min** | **~45–50 min** |

Plus re-running the two Achilles sessions through any stage whose output table is written
whole rather than per session, and the pooled 8c analysis, which is seconds. Call it an
hour of compute once the data is on disk, spread over the 180 s call limit.

### D8.6 — Gatsby_08022013: the XML says 141 channels, the recording is 134 **(blocked)**

The file-size check refuses this session, and the refusal is correct. What it refuses is
worth stating precisely, because the obvious reading of the evidence is the wrong one.

**What the XML says.** `acquisitionSystem/nChannels` is **141**. The anatomical grouping
holds 18 groups: 16 CA1 shanks of 8 channels (0–127) and two groups of three (134–136,
137–139). `spikeDetection` holds 16 groups covering 127 channels — shank 8 has seven,
omitting channel 59. Channels 128–133 and 140 appear in no grouping at all.

**What the data description says.** `Channel_Orderings.pdf` gives Gatsby 128 probe
channels, EMG on 1-based 129–131 and accelerometer on 132–134: **134 total**, with the
auxiliary block immediately after the probe. The XML's auxiliary groups sit at 0-based
134–139, six channels further on.

**What the file says.** 10,720,804,000 bytes. At 134 channels that is exactly 40,003,000
frames with no remainder. At 141 it is 38,017,035 frames with **130 bytes left over** —
under half a frame.

The durations point opposite ways, which is what makes this worth a decision rather than
a fix:

| | frames | duration | against `sessDuration` (30,413.5 s) |
|---|---|---|---|
| 134 channels | 40,003,000 exactly | 32,002.4 s | **+1,588.9 s (26.5 min too long)** |
| 141 channels | 38,017,035 + 130 B | 30,413.6 s | **+0.1 s** |

So 141 matches the session clock and leaves a partial frame; 134 divides exactly and runs
26 minutes long. Taken on arithmetic alone, 141 looks right.

**It is not.** A direct stride test on 60 MB of the recording, 150 s in, settles it:

| | per-channel s.d. (ch 0–7) | s.d. range over ch 0–127 | neighbour correlation (ch 0–63) | ch 0 vs ch 60 |
|---|---|---|---|---|
| **stride 134** | 422, 448, 454, 464, 1362, 507, 474, 580 | 85–2543 | median **+0.925**, min +0.169 | +0.343 |
| stride 141 | 624, 623, 624, 624, 625, 624, 625, 623 | 620–626 | median +0.513, **min +0.509** | +0.187 |

At 134 the channels differ in amplitude the way real electrodes do, one channel is
visibly hot (s.d. 1362 among neighbours near 450), neighbours on a shank correlate at
0.93 and distant channels much less. At 141 every channel has the same s.d. to within
1 %, and every neighbour pair the same correlation to within 0.004. That uniformity is
the signature D6.7 describes: at the wrong stride each "channel" is the same drifting
mixture of the real ones, so they all end up with identical statistics.

**The 26 extra minutes are not an anomaly, they are the known one.** D3.1 records that
this session alone carries 263,976 spikes (4.67 % of its total) after `sessDuration`.
4.67 % of 8.89 h is 25 minutes. The recording is 8.89 h, the authors defined the session
as the first 8.4482 h of it, and both the spike file and the LFP run past that mark. We
already truncate spikes at `sessDuration`; all three epochs end at or before it, so the
extra time is never read.

**Conclusion: the `.eeg` is 134 channels and 8.89 h, and the XML's `nChannels` does not
describe it.** We cannot tell from the released files why — the XML may describe a
wideband acquisition that was reduced before the LFP was written.

**This is blocked rather than fixed, because every way forward is a change.** Using the
session needs two of them:

1. **Override the channel count for this session**, against its own XML. Everything in
   `hc11.lfp` takes the stride from the XML on purpose (D6.7), and a per-session override
   is exactly the kind of exception that makes a pipeline untrustworthy — unless it is
   declared in config, justified by the table above, and applied to nothing else.
2. **Accept a recording longer than `sessDuration`.** `require_consistent_size` currently
   treats any disagreement beyond 1 s as fatal. A file *shorter* than the session is
   dangerous; a file *longer* is not, because every analysis window ends at
   `sessDuration`. Separating the two is defensible on its own merits, but it is still a
   change made while looking at a session that fails.

The alternative is to drop Gatsby_08022013 and report four sessions. That costs the only
BUZ64-probe session in the analysis and the only test of the
data-past-`sessDuration` path, both of which were reasons it was chosen.

**One amendment resolves itself here.** The brief asked for an assertion on the parsed
skip set for Gatsby_08022013, on the understanding that it was the only downloaded
session exercising the XML skip flag. It has **no skip flags at all** — and neither does
the PDF show any for this session; the red channels are Buddy's (24, 27, 58) and
Gatsby_08282013's (35, 45, 47), both now omitted. Gatsby_08022013 excludes its bad
channel by a different mechanism: channel 59 is simply **absent from shank 8's
spikeDetection group**, which gives it seven channels instead of eight. That is the
mechanism channel selection actually uses — it iterates `meta.shanks[shank]` — so the
channel is correctly excluded whatever the skip flag says. The consequence to record is
that **no session in this analysis exercises the `skip="1"` path**, and it remains
covered only by unit tests.
