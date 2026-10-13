# Worklog

One entry per completed step: what was done, what came out of it, what is still open.

## 2026-09-16 — Step 1: repository and environment setup

**Done**

- Created the repository at `Place Field/hc11-place-fields/`, a sibling of the existing
  `hc-11_Data/` and `Article/` folders, so that no data or paper PDF is ever inside the
  git tree (decision D1.1).
- Package layout: `src/hc11/`, `config/`, `scripts/`, `notebooks/`, `tests/`, `docs/`,
  `results/`.
- `pyproject.toml` (hatchling build, `src/` layout) with runtime deps numpy, scipy,
  pandas, h5py, matplotlib, pyyaml, tqdm and a `dev` extra (pytest, pytest-cov, ruff,
  jupyterlab, ipykernel). `uv.lock` committed, 113 packages resolved.
- `.gitignore` excludes `data/`, `*.mat`, `*.pdf`, `results/**`, and the usual Python and
  macOS noise.
- `config/params.yaml`: every analysis parameter for steps 3–5 in one file, each tagged
  `[CHEN2016]`, `[GB2016]`, `[DATA]` or `[PROVISIONAL]`.
- `src/hc11/config.py`: repo-root detection, data-directory resolution
  (`$HC11_DATA_DIR` → `<repo>/data/raw/NoveltySessInfoMatFiles` →
  `<repo>/../hc-11_Data/NoveltySessInfoMatFiles`), session inventory, YAML parameter
  loading with a dotted-key accessor that raises on unknown keys, and `git_commit()` with
  a `-dirty` suffix for figure provenance.
- `docs/decisions.md` seeded with the step-1 decisions and the rationale for every
  provisional parameter already in the config.
- `docs/data_schema.md` stubbed, to be filled from the actual HDF5 tree in step 2.
- Tests: 23 passing (`tests/test_config.py`, `tests/test_data_present.py`). Ruff clean.

**Results**

- All 8 session files are present and resolve correctly:
  Achilles ×2, Buddy ×1, Cicero ×3, Gatsby ×2; 287 MB total.
- Confirmed from the file headers that every session file is a genuine MATLAB v7.3 file:
  128-byte description reads `MATLAB 7.3 MAT-file, Platform: PCWIN64, Created on:
  Fri Nov 25 17:35:00 2016 HDF5 schema 1.00`, with the HDF5 superblock signature at byte
  offset 512 (after the userblock). So `h5py` is the right reader and `scipy.io.loadmat`
  would indeed fail, as expected.

**Open questions**

- `dataset.position_rate_hz` (39.06) and the unit-type label spellings (`pyr` / `int`) are
  placeholders taken from the data description; step 2 must confirm all of them from the
  file and the loader must warn rather than trust them.
- `place_fields.smoothing_sigma_bins`: Chen et al.'s "smoothing over ~5 bins" is
  ambiguous. Currently read as σ ≈ 1 bin (≈5 bins of kernel support). Worth a sensitivity
  check once real maps exist — the alternative reading, σ = 5 bins = 50 cm, would be very
  heavy for a linear track.
- The `results/` sidecar-JSON writer is specified in the config but not yet implemented;
  it lands with the first figure-producing step.

**Time**: ~45 min.

## 2026-09-21 — Step 2: sessInfo loader, validation, schema documentation

**Done**

- Walked the full HDF5 tree of all 8 files (groups, datasets, shapes, dtypes,
  `MATLAB_class` / `MATLAB_empty` attributes). The tree is identical across sessions.
  Documented every field in `docs/data_schema.md`.
- `src/hc11/io.py`: generic MATLAB v7.3 reader (`read_matlab`: structs, char → str,
  `MATLAB_empty` placeholders → zero-size arrays, object references dereferenced,
  unconditional axis reversal) plus hc-11 assembly into an immutable `Session`
  dataclass with read-only arrays. Helpers: `as_vector`, `as_intervals`, `as_int_ids`,
  `decode_matlab_string`, `is_matlab_empty`, `shank_of`, `parse_animal`. Methods:
  `units()`, `spikes_for()`, `spikes_in()`, `epoch()`, `epoch_duration()`,
  `maze_kind`, `track_length_m`, `position_xy_cm`, `position_1d_cm`.
  `list_sessions()` searches recursively, so per-session folders work too.
- `src/hc11/validate.py`: `validate(session) -> list[Issue]` with error / warning / info
  severities covering spikes, unit labels, epochs, position, states, 1-D range.
- Tests: 124 passing, 3 skipped (circular sessions in a linear-only check).
  Synthetic-HDF5 unit tests for every helper (including the 2 × 2 orientation trap,
  empty placeholders, references, non-ASCII strings); one test per validation check;
  integration tests over all 8 real sessions (zero errors, exact known-warning set,
  agreement with `Sessions_Recordings_Summary.pdf`, bit-for-bit agreement with
  `pymatreader`, and a physical 1-D-vs-2-D projection check).
- `scripts/inspect_sessions.py` → `results/session_overview.csv` (+ sidecar JSON).
- README "Loading the data" section; decisions D2.1–D2.6; config `dataset` section
  updated with confirmed values (label-string placeholders removed — the file uses ID
  vectors, not labels).

**Results**

- 8 sessions, 4 animals, 562 putative pyramidal cells and 128 interneurons,
  64.1 M spikes. 5 linear-track sessions (4 × 1.6 m, 1 × 2 m) and 3 circular.
- Every cluster in every session is labelled either pyramidal or interneuron — zero
  unclassified clusters, no overlap, no labelled unit without spikes.
- Position timestamps are an exact, gap-free 25.6 ms grid (39.0625 Hz) spanning the
  MAZE epoch precisely.
- `OneDLocation` is defined for only 2.8–26 min per session (42–93 % NaN) and appears to
  be pre-restricted to track running (median speed 14–55 cm/s where defined vs
  4–22 cm/s where not). On linear tracks it is an exact linear projection of the 2-D
  position (|r| = 1.0000).
- Validation: **0 errors** in all sessions. Warnings only in two sessions (below).

**Discrepancies with the data description / recording summary**

1. **Gatsby_08022013:** 263 976 spikes (4.67 %, all 80 clusters) and 5 state intervals
   extend up to 1 589 s past `sessDuration` (30 413.6 s). `sessDuration`, `POSTEpoch`
   and the recording summary (POST = 10 999 s) all agree with each other; the extra data
   lies outside every documented epoch.
2. **Achilles_11012013:** REM row 12 is zero-length (`[13281, 13281]`).
3. The description says v7.3 files open with `scipy.io` — they don't.
4. `OneDLocation`'s NaN mask is a running mask, not only a tracking-loss mask (the
   description only documents NaN for `TwoDLocation`).
5. Circular maze: 1-D range (2.84–2.90 m) is shorter than a 1 m-diameter circle's
   circumference (3.14 m) and than what the 2-D radius implies (2.82–3.34 m).
6. File names differ from those cited in the description
   (`Sessions_Recordings_Summary.pdf`, `Channel_Orderings.pdf`).

**Open questions (need James's decision)**

- **Q1 — Gatsby_08022013 tail.** Clip to `sessDuration` (i.e. trust the documented
  epochs, which agree with the recording summary) or treat POST as extending to the last
  spike? Irrelevant for Phase 1 MAZE analyses; matters for whole-session firing rates in
  step 3 and for POST events in Phase 2. Recommendation: clip — both the file's epoch
  fields and the summary PDF put the end at 30 413.6 s.
- **Q2 — Running definition for step 4.** `OneDLocation` is already restricted to
  track running by the authors. Options: (a) use the authors' mask *and* our
  speed > 15 cm/s threshold; (b) authors' mask only; (c) ignore it, compute running from
  2-D speed alone and linearise ourselves. Recommendation: (a) as primary, with (c) as a
  sensitivity check. Note Buddy_06272013 has only 2.8 min of defined 1-D position.
- **Q3 — Zero-length REM row.** Drop zero-length intervals when building state masks
  (harmless either way for durations)? Recommendation: drop, with a warning.
- Fig. 1A of the paper shows 77 simultaneously recorded place cells in one session; only
  Achilles_10252013 (120 pyr), Achilles_11012013 (92) and Cicero_09102014 (81) have
  enough pyramidal units — a useful cross-check for step 4 place-cell criteria.

**Time**: ~2 h.

## 2026-09-30 — Step 3 (part A): anomaly decisions implemented

**Done**

- `src/hc11/preprocess.py`: `truncate_to_session` (+ `TruncationReport`),
  `clean_intervals`, `state_intervals`, `intervals_to_mask`, `mask_to_intervals`.
- `src/hc11/track.py`: `Track` built from the observed 1-D range per session, with
  wrap-aware `difference` / `distance` / `wrap` / `unwrap` and spatial binning.
- `src/hc11/behavior.py`: NaN-aware speed from 2-D tracking, and a two-stage
  `running_mask` (authors' mask, then speed threshold) reporting what each stage does.
- `scripts/report_running.py` → `results/running_summary.csv`.
- Config: new `preprocess` and `track` sections; `behavior.require_linearized`.
- Decisions D3.1–D3.5; 45 new tests (177 total, 3 skipped).

**Results**

- Truncation affects only Gatsby_08022013: 263 976 spikes (4.67 %), 4 state intervals
  dropped, 1 clipped. All other sessions are untouched (asserted in tests).
- Zero-length intervals affect only Achilles_11012013 (1 REM row), also asserted.
- Running time after both stages: 2.6–15.4 min per session, 6.7–34.6 % of the MAZE epoch.

**The speed threshold is not a formality (D3.4)**

Removal ranges from 3.3 % (Buddy) to 51.9 % (Cicero_09012014). In the three
low-removal sessions the authors' mask is built of single-traversal blocks and the slow
samples are at the track ends; in the other five it contains blocks up to 23–59 s with
8–17 % of samples below 5 cm/s in mid-track. `OneDLocation` marks the linearised track,
not running, and how strictly it was restricted varies by session.

**Open questions**

- The circular 1-D range vs. circumference mismatch (D3.3) remains unexplained.
- Cicero_09012014 loses over half its linearised time to the speed threshold, leaving
  8.4 min of running from a 92.7 min MAZE epoch. Worth watching in step 4.

**Time**: ~1 h.

## 2026-09-30 — Step 4 (part B): directional place fields

**Done**

- `src/hc11/laps.py`: direction from a smoothed derivative with hysteresis, traversals
  segmented on the authors' mask blocks (D4.1), split at the reward site on circular
  tracks so one lap is one circuit (D4.3).
- `src/hc11/fields.py`: sample selection, occupancy and spike maps smoothed separately
  (wraparound on circular tracks), Skaggs information in bits/spike and bits/second,
  field detection, split-half stability, laps active, and a circular-shift shuffle in
  compressed running time (D4.4).
- `src/hc11/plots.py`: lap plot, sorted rate-map heatmap, field grid, information
  histogram. Two fixed colourblind-safe hues for direction (checked: worst-case protan
  ΔE 21.9), a perceptually uniform ramp for rate.
- `scripts/place_fields.py`, `scripts/sensitivity_running.py`.
- 41 new tests (218 total, 3 skipped): hand-computed Skaggs cases, synthetic Poisson
  place-cell recovery, shuffle calibration, lap segmentation including the circular wrap.

**Laps per direction** (the adequacy number)

| Session | Maze | pos | neg | Pyr | Place cells |
|---|---|---|---|---|---|
| Achilles_10252013 | 1.6 m linear | 40 | 42 | 120 | 76 (63 %) |
| Achilles_11012013 | circular | 1 | 75 | 92 | 72 (78 %) |
| Buddy_06272013 | 1.6 m linear | 27 | 24 | 48 | 16 (33 %) |
| Cicero_09012014 | 1.6 m linear | 44 | 44 | 55 | 16 (29 %) |
| Cicero_09102014 | circular | 0 | 19 | 81 | 35 (43 %) |
| Cicero_09172014 | 2 m linear | 28 | 27 | 59 | 23 (39 %) |
| Gatsby_08022013 | 1.6 m linear | 41 | 42 | 66 | 29 (44 %) |
| Gatsby_08282013 | circular | 1 | 79 | 41 | 29 (71 %) |
| **Total** | | **182** | **352** | **562** | **296 (53 %)** |

Under-sampled directions: the three circular sessions, which are unidirectional by
design; Cicero_09102014 has only 19 laps in its one direction. Every linear session has
24–44 laps per direction. Buddy_06272013, flagged as thin on running time (2.6 min),
turns out to have 27/24 laps — adequate by lap count, but its cells fire so few spikes
during running (median 12 per cell-direction) that the 50-spike floor excludes most.

**Validation**: Achilles_10252013's sorted heatmap reproduces the staircase of Fig. 1A,
and individual fields are single-peaked with peaks of 1.7–26.6 Hz.

**Open questions**

- 296 place cells vs the paper's 491 (1.66x), outside the agreed tolerance. Criteria were
  not tuned. The binding constraint is the provisional 120 cm field-width cap: 105
  cell-directions have significant spatial information and a peak above 1 Hz but no
  detected field. Relaxing that cap alone gives 335 (1.47x); also lowering the spike floor
  to 20 gives 356 (1.38x). Needs James's decision — see D4.6.
- The speed threshold barely changes the fields (D4.5): with laps held fixed, place-cell
  counts and median information and stability are unchanged even where it discards 59 %
  of samples. Kept as primary for cross-session comparability, not because the fields
  need it.

**Time**: ~3 h.

## 2026-10-05 — Step 5a: place-cell count resolved

**Done**

- Established that the paper's 491 counts unique neurons (D5.1); the pipeline now reports
  unique cells and cell-direction pairs side by side everywhere.
- Replaced the fixed 120 cm field-width cap with `max_field_width_frac: 0.75` of the
  usable track extent (D5.2), so the 1.6 m, 2 m and circular tracks are treated alike.
- `scripts/place_cell_sensitivity.py` → `results/place_cell_sensitivity.csv`.

**Results**

| Variant | Unique | % pyr | 491/ours | Cell-dir |
|---|---|---|---|---|
| primary | 304 | 54 % | 1.62 | 370 |
| cap relaxed | 346 | 62 % | 1.42 | 458 |
| cap relaxed + 20 spikes | 367 | 65 % | 1.34 | 489 |
| peak ≥ 1 Hz only | 392 | 70 % | 1.25 | 618 |
| paper | 491 | 87 % | — | — |

The fraction cap added 8 unique cells (296 → 304), all on the 2 m and circular sessions.

**The gap is not a threshold we mis-set.** Applying only Chen et al.'s stated criterion
(peak > 1 Hz) with nothing else still gives 392; the most inclusive base available (all
samples with a defined OneDLocation, directions merged) gives 353. 186 of 562 pyramidal
cells fire < 50 spikes in their best direction during the released running periods.
Recorded in D5.4 and left as a finding, not tuned away.

**Time**: ~1 h.

## 2026-10-06 — Step 5a diagnostics (D5.5)

**Done**

- `src/hc11/linearize.py` + `scripts/diagnostic_linearization.py`: our own linearisation
  of the 2-D position (principal axis for linear tracks, arc length for the ring) and a
  rerun of the primary criteria on that wider base. Diagnostic only — the authors'
  `OneDLocation` stays the replication base everywhere else.
- Softened the r = 0.63 claim in D5.4 (n = 8, p ≈ 0.09: consistent with, not evidence for).
- 12 new tests (230 total, 3 skipped), including that our linearisation reproduces the
  authors' on all five linear sessions at r > 0.9999.

**Results**

- Fig. 1A shows 77 place cells on *rightward runs* — one direction. Achilles_10252013
  gives 55 and 52 per direction, 77 unique. Neither direction is near 77, so our unique
  count matching their one-direction count is a coincidence. It does give a tighter
  anchor than the 491 aggregate: 77/120 = 64 % per direction for them, 55/120 = 46 % for
  us, a factor of 1.4 — the same as the overall gap.
- Wider base, linear sessions only (where our linearisation is exact): coverage 7–32 % →
  68–98 % of the MAZE epoch, place cells 164 → 207 (47 % → 59 % of 348 pyramidal cells).
  Buddy_06272013 doubles, 16 → 32, which is the session whose released mask was thinnest.
- Circular sessions go the other way (140 → 111), but our arc-length reconstruction there
  is only r = 0.972–0.997 against theirs and disagrees on extent, so those rows measure
  our linearisation rather than the base. The all-session total (304 → 318) understates
  the effect for that reason.

**Conclusion**: the restricted `OneDLocation` mask is a real contributor to the gap, worth
roughly 12 percentage points on the linear sessions, but 59 % is still well short of the
paper's 87 %. The remainder is unexplained.

**Time**: ~1 h.

## 2026-10-06 — Step 5b: Bayesian position decoding

**Done**

- `src/hc11/decode.py`: memoryless Poisson decoder with a uniform prior, computed in log
  space; leave-one-lap-out cross-validation with precomputed per-lap templates; chance
  baseline by permuting tuning curves between cells; ensemble-size sweep.
- `src/hc11/plots.py`: posterior heatmaps with the true trajectory, error distributions,
  error vs ensemble size, confusion matrices.
- `scripts/decode_position.py` → `results/decoding_summary.csv`,
  `results/decoding_ensemble.csv`, per-session figures.
- `select_running_samples(direction=None)` for the direction-agnostic decoder.
- 20 new tests (250 total, 3 skipped), including three that enforce no leakage.

**Results (250 ms, directional)**: median error 4.7–10.6 cm against a chance level of
36–71 cm. Achilles_10252013 reaches 4.7–5.1 cm with 77 cells; the 16-cell sessions sit at
8.8–10.6 cm. Error is lower at the track ends (6.1 cm) than in the middle (8.3 cm).

**Better than the expected "low tens of cm", and verified not to be leakage**: the honest
error is 0.89 cm worse than a deliberately leaky all-laps template, a 10 cm spatial bin
puts a floor near 2.5 cm on any median error, and a 1-D track with 77 cells is easier than
the 2-D open fields behind Chen et al.'s 8.5–12.5 cm.

**Directional split** is worth 6–34 % on linear sessions, and nothing on the circular ones,
which are unidirectional.

**20 ms bins**: median 10–24 cm, 5–60 % of bins undecodable, median spikes per bin 0–4, p90
at or above chance. Expected, and the useful number for Phase 2: a 20 ms bin here is worth
about one spike, so sequence scoring must rely on the trajectory across bins.

**Ensemble size**: 12.1 cm at 5 cells → 5.2 cm at 55, flattening past ~30. Our place-cell
shortfall (D5.4) therefore costs decoding accuracy directly.

**Time**: ~2 h.

## 2026-10-07 — Step 6c (first half): MUA event detection

**Done**

- `src/hc11/events.py`: shared `EventList`, two-threshold detector, per-event
  `participation`, and `compare_events` for the forthcoming LFP cross-check.
- `src/hc11/mua.py`: per-state-interval population rate (1 ms bins, 10 ms Gaussian),
  pooled per-state baseline, detection via the shared rule.
- `src/hc11/plots.py`: event rasters, event rate over the session, participation.
- `scripts/detect_mua_events.py` → `mua_summary.csv`, `mua_events.csv`,
  `mua_threshold_sensitivity.csv`.
- 31 new tests (281 total, 3 skipped).

**Results (PRE/POST, non-REM)**

| Session | Epoch | State min | Events | /s | Median dur | Median active |
|---|---|---|---|---|---|---|
| Achilles_10252013 | PRE | 168.5 | 7046 | 0.70 | 68 ms | 18/120 |
| Achilles_10252013 | POST | 89.7 | 3955 | 0.74 | 62 ms | 18/120 |
| Achilles_11012013 | PRE | 104.0 | 4709 | 0.76 | 57 ms | 11/92 |
| Achilles_11012013 | POST | 73.0 | 3328 | 0.76 | 58 ms | 12/92 |

- Durations match expectation: median 57–68 ms, about half in the 50–100 ms band.
- Participation matches the paper's "10–15 % of the pyramidal population per ripple":
  median 15 % and 13 %. The distribution is cleanly lognormal — raw skew +1.5/+1.8,
  skew of the log +0.17 — i.e. many small events and few large ones, as Chen et al. note.
- **Decodability for step 7**: 97 %/94 % of events have ≥ 5 place cells active, 70 %/51 %
  have ≥ 10, and the median event has 12/10 place cells active. There is plenty to score.
- Rate is 0.70–0.76/s, above the expected 0.2–0.5/s. The threshold sweep puts 4 SD at
  0.39–0.41/s, inside the range. Not changed — see D6.4 and the report to James.

**Blocked**: MAZE events (D6.3). Immobility cannot be defined from the released position
data; it needs the theta/delta ratio or EMG from the `.eeg` file.

**Time**: ~2 h.

## 2026-10-08 — Step 6c addenda: jitter null, state asymmetry, MAZE hypothesis

**Done**

- `mua.jitter_spikes` + `scripts/mua_jitter_null.py` → `mua_jitter_null.csv`,
  `mua_jitter_per_event.csv`. Per-cell independent jitter, U(±100 ms), 5 repeats.
- Revised D6.3: the MAZE tracking-lost runs are an open hypothesis for the LFP to
  resolve, not ruled out. Elevated population rate there is ambiguous — awake
  sharp-wave ripples at an occluded reward site look exactly like that.
- D6.6: per-epoch scored-state breakdown for both sessions.

**Results**

- The events are real synchrony. Observed events outnumber jitter-null events by
  1.7–2.2x, and in the per-event comparison (participation recounted in the *observed*
  windows using jittered spikes) observed participation is ~1.7x the null, with
  **84–90 % of individual events above the 95th percentile of their own null**.
- Recorded a trap: comparing median participation of *re-detected* null events against
  observed suggests the null is more synchronous (22 vs 18 cells). That is a duration
  artefact — jittered events are nearly twice as long (124 vs 68 ms). Per unit duration
  the observed events are denser (2.50 vs 1.85 cells per 10 ms).
- POST has less non-REM than PRE in both sessions: 168.5 → 89.7 min and 104.0 → 73.0 min.
  The proportion falls in one session (56 % → 37 %) and holds in the other, where REM
  collapses instead (27 % → 5 % of scored sleep). Fig. 1C must therefore normalise by
  searched time, and POST rests on about half the events PRE does in Achilles_10252013.

**Time**: ~1 h.

## 2026-10-08 — Step 6a infrastructure: XML parser and memory-mapped .eeg reader

**Done** (written against synthetic fixtures; no real LFP yet)

- `src/hc11/lfp.py`: `parse_xml` (channel count, LFP sampling rate, wideband rate, bit
  depth, voltage range, amplification, shank grouping, skipped channels),
  `LfpMetadata` with `n_channels_total` and `probe_channels` kept separate,
  `LfpFile` with memory-mapped `read()`, `check_size()` and `require_consistent_size()`,
  and `$HC11_LFP_DIR` path resolution supporting flat and per-session layouts.
- Config `lfp` section; README section on where to put the `.eeg` files.
- Decisions D6.7–D6.9.
- 27 new tests (315 total, 3 skipped). Fixtures are 8 shanks × 16 probe channels plus 8
  EMG/accelerometer channels = 136 total, so the two channel counts are different
  numbers and a confusion between them cannot pass. Channel values encode their own
  index, so a mis-strided read is visible. One fixture is 544 MB, proportionally
  equivalent to the real 11 GB file, and asserts a 1 s read materialises < 10 kB.

**Not done**: channel selection by ripple-band power (6a-2) needs real LFP.

**Time**: ~1 h.


## 2026-10-11 — Step 6a-2: LFP in place, file geometry verified, ripple channels selected

**Done**

- Extracted both `.eeg` files from the CRCNS archives (11.7 and 12.3 GB) with a new
  resumable extractor, `src/hc11/eeg_archive.py`. Two earlier plain `tar -x` attempts had
  been killed partway and left 4.6 and 3.5 GB fragments with no error; the extractor
  appends, stops on a deadline at a chunk boundary, and resumes by skipping what is
  already on disk. `member_info` reads the uncompressed size from the tar header without
  decompressing the body.
- **`require_consistent_size` passes on both files, exactly.** 11,678,469,572 and
  12,270,786,020 bytes are both whole numbers of 268-byte frames with remainder zero, and
  both durations match `sessDuration` to +0.000 s (34,861.103 s and 36,629.212 s). 134 is
  the only plausible channel count that divides either file exactly, so the stride is
  confirmed by the data rather than only by the XML.
- Verified the channel grouping against `Channel_Orderings.pdf`: 12 CA1 shanks of 10
  channels (spike groups 1–6, 8–13), two four-channel "Extra Channels" groups (7 and 14),
  EMG on 0-based 128–130 and accelerometer on 131–133. Every line matches. The EMG
  channels carry 9.8–21.8 of ripple-band envelope against 2.4–5.0 for the accelerometer,
  which is an independent confirmation of which is which — and a reminder that EMG
  channel 130 (20.0) would have outranked two real spike groups in a naive selection.
- Selected one channel per shank by mean 150–300 Hz envelope over 10 min of sampled
  non-REM, restricted to shanks carrying sorted units. Primary detection channel: **ch 2**
  (Achilles_10252013) and **ch 3** (Achilles_11012013), both left CA1 shank 1. Full table
  in D6.13; `results/ripple_channels.csv` and `results/channel_profiles.csv`.
- Figures per session: `lfp_channel_power.png` (profile across channels and collapsed
  onto depth), `lfp_band_profile.png` (delta/theta/ripple together), and
  `lfp_ripple_examples.png` (raw and band-passed LFP around the nine strongest
  ripple-band peaks).
- New modules `src/hc11/eeg_archive.py` and `src/hc11/ripple_channels.py`, new script
  `scripts/select_ripple_channels.py`, `lfp.size_check`/`require_size_consistency`
  factored out so the archive path shares the file path's validation, config section
  `lfp.channel_selection`, decisions D6.10–D6.15. 43 new tests (358 total, 3 skipped).

**Results worth looking at**

- The ripple-band profile is a clean sawtooth: each shank peaks two to three channels
  from its top and decays monotonically with depth. Theta and delta peak several channels
  *deeper* than the ripple band, which is the CA1 depth signature and is the strongest
  evidence that the selection is tracking anatomy.
- Seven of twelve shanks select the identical channel across two sessions a week apart;
  the rest differ by one or two, all on the right-hemisphere probe.
- The example traces are textbook: 40–80 ms spindle-shaped bursts at roughly 180–200 Hz
  on a slow deflection, often in pairs 100–200 ms apart.

**Open questions**

- Within-shank margins are 1–14 %, so the winner and its neighbour are near-equivalent.
  Ranking by excursion (p99/median) instead of level moves the choice by one channel on
  8 of 12 shanks in one session, and by seven channels on spike group 3. Whether that
  matters is a step-6b question: detect on both and compare the event sets.
- Whether to detect on one channel or combine shanks is still open. Everything is
  recorded per shank so either is available.

**Cost of the bug** (D6.15): the first run reported a uniform 44–46 envelope on all 134
channels with NaN theta/delta, from indexing `data[channel]` (a time point) instead of
`data[:, channel]`. It produced no error and plausible numbers, and several hours went
into verifying the recording before the fault was located in the analysis. The fix is
structural, not a patch: explicit `axis=0`, a per-channel loop over 1-D arrays, and
`sosfiltfilt` instead of `filtfilt`.

**Time**: ~3 h, of which about half was the extraction and the bug above.

## 2026-10-13 — Steps 6b and 6c: ripple detection and the MUA cross-check

**Done**

- `src/hc11/ripples.py`: whole-recording chunked ripple-band envelope, two-threshold
  detection per state through the shared `hc11.events` machinery, per-event multi-shank
  support, and the windowed z-score measurement behind it.
  `scripts/detect_ripples.py` (6b) and `scripts/cross_check_events.py` (6c).
  Config section `events.ripples`. Decisions D6.16–D6.21. 23 new tests (381 total,
  3 skipped).
- Detected on both candidate channels in both sessions, plus both channels of spike
  group 3 of Achilles_10252013, where the two rankings land 140 µm apart.

**Results**

- Event rates, primary channel: PRE 28.3 and 32.0 per minute of non-REM, POST 36.6 and
  39.3, MAZE 44.5 and 63.6 per minute of (very little) immobility. Median duration
  54–64 ms. Rates are per minute of searched state throughout, never counts.
- **The two channel rankings are interchangeable at the primary**: the best channels are
  adjacent contacts and their event sets match 94–98 % at 0.0 ms median offset.
- **Where they differ by seven contacts, the level ranking wins** (D6.17). Both channels'
  unique events are weak but real — median 6–7 of 11 other shanks show them, against
  11/11 for shared events — and the level-ranked channel finds more of them. The
  disagreement figure shows the excursion-ranked channel oscillating near-continuously
  in the ripple band rather than in discrete bursts. Config set to `primary_ranking:
  level`, against the brief's suggested default, with the evidence in D6.17.
- **Multi-shank support is strongly bimodal towards agreement**: median 11/11 in five of
  six session × epoch combinations. In Achilles_10252013 POST, 2366 of 3284 events are on
  all 11 other shanks and 13 on none.
- **Cross-check**: 70–80 % of ripples have an MUA partner, 47–58 % the other way. The
  population-rate peak leads the ripple envelope peak by about 5 ms in every epoch of
  both sessions (IQR 15–17 ms), part of which may be the 8 ms envelope smoothing.
- **The two independent checks agree about which events are weak.** Agreement with the
  MUA detector rises monotonically with multi-shank support: 37–38 % at zero support,
  52–53 % at ≤ 2, 80–85 % at ≥ 9.
- Recommendation (D6.20): **MUA primary, ripple matching as a per-event property**, with
  Chen et al.'s conjunction available as the strict variant. LFP exists for 2 of 8
  sessions; the sequence analysis is a spike analysis; and ripples are close to a subset
  of the MUA set.

**Throughput and extraction — what it actually costs**

Correcting an earlier assumption. The 70 MB/s figure recorded in D6.10 is the **write**
throughput of the connected folder. Reads are about 35× faster:

| operation | measured |
|---|---|
| sequential read from the connected folder | 2.0–2.6 GB/s |
| write to the connected folder | ~70 MB/s |
| local VM scratch disk | 5.6 GB free — too small for one `.eeg` |
| one channel across a whole 11.7 GB `.eeg` (memory-mapped, strided) | ~5 s |
| band-pass + Hilbert envelope of one whole channel | ~15 s |
| **a full single-channel pass, end to end** | **~20 s** |
| extracting one `.eeg` from its archive (write-bound) | ~170 s |

So a second pass over the full LFP is cheap *if it is a read*. Twelve channels is about
four minutes per session; all 134 would be about 45 minutes. What is expensive is writing
— which is why extraction is resumable (D6.10) and why nothing else writes a large file
to that folder.

Two environment constraints shape the code and are not visible in it:

- **Background jobs do not survive.** Work started with `&`, `nohup` or `setsid` is killed
  at the end of the call that started it, silently and part-way. Two 12 GB extractions
  died that way leaving 4.6 and 3.5 GB fragments with no error. Every long operation must
  therefore finish inside a single ≤180 s call, which is why `eeg_archive.extract` takes a
  deadline and resumes.
- **The local scratch disk cannot hold one recording**, so nothing can be staged whole.

What this means for later steps: the streaming approach that was used for channel
selection (decompress the archive and keep only sampled windows, ~75 s per session) is no
longer needed now that the files are extracted, and the memory map is both simpler and
faster. If the remaining six sessions' LFP is wanted, budget ~9.5 GB download and ~170 s
extraction each, about 76 GB in total, and the analysis time is minutes, not hours.

**Open**

- ~45 % of MUA events have no ripple. Worth characterising in 6d — they may be
  lower-synchrony bursts, events on shanks the primary channel does not cover, or noise.
- The 5 ms MUA-to-ripple lead should be re-measured with a smaller smoothing kernel
  before it is quoted as physiology.
- D6.3 (MAZE immobility) is unchanged but better supported: the few MAZE events that are
  detected look like real awake ripples (median support 10–11 of 11), so the shortage is
  of detected immobility rather than of events. EMG channels 128–130 are now identified
  and the band machinery exists, so a theta/delta or EMG state score is reachable.

**Time**: ~2 h.

## 2026-10-13 (later) — Step 6d: characterisation, and the three loose ends

**Done**

- `scripts/offset_vs_smoothing.py`: sweeps the envelope smoothing kernel from 0 to 32 ms
  and re-measures the MUA-to-ripple offset, both with the event set held fixed and
  re-detecting at each width.
- `scripts/characterise_events.py`: per-epoch characterisation of each candidate-event
  set, the decodability comparison between the full MUA set and the strict conjunction,
  and the three-way test on the unmatched MUA events.
- `hc11.events.participation_by_shank` and `ShankParticipation`, with the concentration
  and per-shank share measures the shank-coverage hypothesis needs. 9 new tests
  (390 total, 3 skipped). Decisions D6.22–D6.24.

**1. The offset is real.** Median MUA-minus-ripple offset barely moves across kernel
widths — −4.5 ms at zero smoothing against −4.3 ms at the configured 8 ms in
Achilles_10252013 PRE, and the same flatness in all four session × epoch combinations.
Only at 32 ms, four times the configured kernel, does it grow by a consistent ~0.7 ms.
So the kernel contributes a fraction of a millisecond to a 5 ms effect. Reported as a
result: population firing peaks 4.5–6 ms before the ripple envelope peak. Good to about
±1 ms — one LFP sample is 0.8 ms — and the IQR of 15–17 ms means it describes a broad
relationship, not a latency. With the question settled the offset now appears in
`cross_check_offsets.png`.

**2. The unmatched MUA events are smaller, not elsewhere.**

- *Low-synchrony bursts*: **supported, and it is the whole story.** Unmatched events
  carry about 60 % of the active cells of matched ones (mean 15.3 against 25.3), with a
  rank-biserial of −0.58 to −0.63 across both sessions. They are also shorter (57 against
  82 ms) and weaker (peak z 3.78 against 4.62).
- *Shank coverage*: **rejected.** Raw concentration is higher for unmatched events, but
  that is the cell-count confound — fewer cells forces higher concentration. Against a
  per-event null drawing the same number of cells at random, weighted by each shank's cell
  count, the excess is −0.007 for both groups (rb +0.03). The prediction was that the
  active cells would sit *away* from the detection channel's shank; the observed share on
  that shank is 0.138 against 0.124 and the expected share from cell counts is 0.117, so
  both groups are at chance.
- *Detection noise*: **a small contribution.** Unmatched events fail the spike-jitter null
  5–9 points more often (86.7 % pass against 91.8 %), so most of them are still genuine
  co-firing.

This agrees with the 6c finding that multi-shank support predicts MUA agreement. Both
detectors measure the same underlying strength and put the threshold in different places.
Consequence for the plan: the LFP is not covering a part of the recording the spikes
miss, so there is no coverage argument for an LFP-primary set, and the conjunction is a
stricter threshold rather than a different filter.

**3. What the strict variant costs.**

| Session | Epoch | MUA usable | conjunction usable | strict keeps |
|---|---|---|---|---|
| 10252013 | PRE | 5827 / 7046 (83 %) | 3433 / 3601 (95 %) | 59 % |
| 10252013 | POST | 3368 / 3955 (85 %) | 2201 / 2309 (95 %) | 65 % |
| 11012013 | PRE | 3198 / 4709 (68 %) | 1924 / 2213 (87 %) | 60 % |
| 11012013 | POST | 2465 / 3328 (74 %) | 1717 / 1927 (89 %) | 70 % |

30–41 % of usable events, for a rise in the decodable fraction from 68–85 % to 87–95 %.
Even the smallest cell is 1717 events, so both variants are workable for Fig. 1C; the
cost falls on the tails of the score distribution.

**Also worth noting**

- Achilles_11012013 has a markedly lower decodable fraction than Achilles_10252013 (68–74 %
  against 83–85 %) on the full MUA set. It has fewer place cells, and its events have
  fewer active cells (median 14–15 against 21–22). Not a problem, but it means the two
  sessions will not contribute equally to any pooled distribution.
- MAZE still yields 20 events per session. Unchanged from D6.3 and still waiting on a
  state score from EMG or theta/delta.

**Open for step 7**

- Primary candidate-event set: full MUA, with `ripple_matched` and `n_supporting_shanks`
  carried per event so the conjunction is a filter rather than a re-run.
- Whether to pool the two sessions or report them separately, given the participation
  difference above.

**Time**: ~1.5 h.
