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
