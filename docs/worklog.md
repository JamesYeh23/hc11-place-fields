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

## 2026-10-14 — Step 7: sequence scores and the Fig. 1C comparison

**Done**

- `src/hc11/sequences.py`: event binning, decoding against each direction's fields, the
  weighted correlation over the whole posterior (vectorised over a leading shuffle axis),
  both shuffle kinds, and `rZ`. `scripts/score_events.py` (7a, 7b, plus the
  `--identity-shuffle` control), `scripts/validate_sequences.py` (both controls),
  `scripts/replay_comparison.py` (7c). Config section `sequences`. Decisions D7.1–D7.8.
  28 new tests (418 total, 3 skipped).

**The headline: POST exceeds PRE in all eight cells of the design** — two sessions × two
shuffle kinds × two event sets. Median `rZ` differences +0.197 to +0.525, rank-biserial
+0.071 to +0.189, p from 0.026 to 6e−08. Subsampling PRE to POST's event count 200 times
returns the same difference to three decimals. Per minute of searched state the margin is
wider still: POST has about half the non-REM of PRE and produces more significant events
per minute (2.01 → 2.69 and 0.83 → 1.58 in one session, 0.45 → 0.81 and 0.70 → 1.72 in
the other).

**Both controls pass.** Identity-shuffled place fields collapse the PRE/POST gap to
−0.005, −0.033, +0.049 and +0.071 against real gaps of +0.302, +0.220, +0.197 and +0.433.
MAZE running at 2 s windows scores a median `rZ` of +7.5 to +13.2 with 82–95 % above 1.65,
so the scoring works on a trajectory that is known.

**Three things found along the way that mattered more than expected**

1. *A one-lap template was costing the whole effect.* Achilles_11012013's "pos" direction
   has one lap of 226 samples against 75 laps for "neg". Because the score maximises over
   directions, a noise template gives the *shuffles* the same extra freedom, the null
   rises, and every event's score falls. Removing it moved that session's medians from
   −0.060 to +0.110 (PRE) and +0.148 to +0.307 (POST) — a shift the size of the effect
   being measured. Now a five-lap floor in config (D7.3).
2. *The absolute level of `rZ` is not interpretable.* With place fields detached from
   cells, median `rZ` is still +0.59 under column-cycle. The column-cycle shuffle spreads
   posterior mass uniformly over position, so it rewards an event whose decoded position
   is merely localised, sequence or not. Time-permutation preserves each bin's posterior
   and is the cleaner test, but it is not zero either (+0.25). Only the PRE/POST
   difference survives the control, and no claim here rests on "x % of events exceed
   1.65" (D7.6).
3. *A 200 ms window of real running scores zero.* At ~30 cm/s it covers under one 10 cm
   position bin, so it cannot contain a trajectory whatever the decoder does. That is why
   the running control is reported at two lengths, and it sets the method's sensitivity
   ceiling (D7.7).

**The uncomfortable number**: the five-bin floor at 20 ms needs a 100 ms event, and our
events have a median duration of 55–68 ms, so **only 15–24 % of candidate events are
scored**. That is a larger loss than every other filter combined, and it traces back to
the MUA detector's 1 SD boundaries rather than to anything in step 7. Lowering the floor
to three bins would mean fitting a correlation to three points, which is not an estimate.
Recorded in D7.5 rather than tuned.

**Measured, worth keeping**

- During events a 20 ms bin carries 5–6 spikes from the place-cell ensemble, not the ~1
  spike step 5 found during running, and under 1 % of bins are empty. The posteriors are
  sharper than the brief anticipated.
- Forward and reverse traversals split 49–52 % in every epoch of both sessions, with no
  consistent score asymmetry.
- The conjunction event set gives a *larger* PRE/POST effect than the full set in all four
  comparisons, consistent with D6.23's finding that it is a stricter threshold on event
  strength.

**Open**

- The effect is consistent but small: rank-biserial +0.07 to +0.19 is a 54–59 % chance
  that a random POST event outscores a random PRE event. Before the magnitude can be
  compared against the paper's, D7.5 (scoring rate) and D7.6 (uninterpretable level)
  would both have to be addressed.
- Longer candidate events would address D7.5 directly. That means revisiting the MUA
  boundary definition — a step-6 change, not a step-7 one.

**Time**: ~2.5 h.

## 2026-10-16 — Step 8a: audit before scaling to eight sessions

**Done.** No analysis run and no parameter changed. Everything below comes from results
already in the repo, plus a read of `Channel_Orderings.pdf` for the three animals we had
not processed. Decisions D8.1–D8.5.

**The parameter audit (D8.1).** Place-field criteria, the 10 cm bin and the 15 cm/s speed
threshold have all been applied to eight sessions, and the field-width fraction was
introduced *because* a fixed centimetre value was maze-dependent — those are
session-general. Everything downstream of place fields has only ever been checked on two
sessions: the MUA thresholds and duration limits, every ripple parameter, the channel
selection sample, the five-bin scoring floor and the decodability criterion.

The five-lap floor turns out not to be an Achilles artefact: lap counts per direction are
24–44 on every linear session and 0 or 1 on the unused direction of every circular one, so
it separates the two geometries with a wide margin either side. It will fire on exactly
the two circular sessions still to come.

**Two probe geometries, not one (D8.2).** Achilles and Cicero use a BUZ64SP with 12 CA1
shanks of 10 channels plus two four-channel "extra channel" groups. Buddy and Gatsby use a
BUZ64 with **16 shanks of 8 and no extra groups**; Buddy has **132** channels, not 134,
and its accelerometer and EMG are in the opposite order with a single EMG channel. Three
consequences for step 6: the multi-shank support denominator varies from 5 to 11 across
sessions, so the count must be reported as a fraction too; and **the XML skip flag has
never run on real data** — Achilles skips nothing, while Buddy skips 1-based 24, 27, 58
and Gatsby_08282013 skips 35, 45, 47. If those do not come through as `skip="1"` the code
will silently keep three bad channels in the selection pool, so 8b should check against
those numbers rather than trust it.

**The two validated sessions are the two best-case sessions (D8.3).** Achilles has the
largest two ensembles in the dataset, 77 and 73 place cells; the rest have 16 to 36. At
the 20 ms bin the scoring uses, measured during running, Achilles carries 2–4 spikes per
bin with 5–13 % of bins empty, while Buddy and Cicero_09012014 carry a median of **zero**
with 47–60 % empty. Since the five-bin floor counts non-empty bins and already discarded
76–85 % of candidate events on Achilles, the scorable fraction on the weakest sessions
could be low single digits. Prediction, not result, and not a reason to change anything —
but a reason to expect two to four of the six new sessions to return too few scorable
events for a meaningful PRE/POST test, and to plan on reporting that with its n rather
than as a negative finding.

**Circular sessions need their own column (D8.4).** Only **two** of the six remaining are
circular (Cicero_09102014, Gatsby_08282013), not three — the third, Achilles_11012013, is
already analysed; the final split is five linear, three circular. The decoder is not the
problem: the circular sessions decode to 4.9–5.4 cm at 250 ms, as good as the linear ones,
and the approximate linearisation never enters a primary result because every primary
analysis uses the authors' own `OneDLocation`. The problem is the score. On a ring the
weighted correlation is maximised over all 28–29 circular offsets and on a line over one,
so the circular null is broader and `rZ` is compressed for the same quality of sequence.
The two sessions already in hand show it: the circular session's column-cycle medians are
about a fifth of the linear session's while its time-permutation medians are nearly
double. Absolute `rZ` must not be pooled across geometries; the *count* of sessions with
POST > PRE can be, because a sign has no scale.

**Cost (D8.5).** Predicted `.eeg` sizes from `sessDuration × 1250 × n_channels × 2`
reproduce both known files to the byte, so the six remaining come to **64.1 GB extracted,
~51 GB of archives**. 197 GB free; peak ~115 GB if archives are kept. Compute is about
7–8 minutes per session, ~45–50 minutes for six, an hour all in.

**One code change 8b will need**, flagged rather than made: eleven scripts default to a
two-session list. Widening a default is not a methodological change, but it is a change,
so it belongs in 8b's commit rather than this one.

**Time**: ~45 min.

## 2026-10-17 — Step 8b (part 1): four sessions through the pipeline, and a blocked fifth

**Scope.** Five sessions, not eight: Buddy_06272013, Cicero_09172014 and Gatsby_08282013
have no LFP downloaded, for storage reasons rather than methodological ones. Recorded in
`config.LFP_SESSIONS` and `config.LFP_OMITTED` so any result table can say what is
missing. Four of the five are now through; Gatsby_08022013 is blocked on a decision
(below).

**Infrastructure, no methodology.**

- `config.LFP_SESSIONS` replaces the two-session default hard-coded in ten scripts.
- `src/hc11/results.py`: `write_session_table` merges per-session result tables instead
  of rewriting them. Every script wrote its whole CSV, so running one session would have
  silently deleted the other seven — and it made a five-session run impossible to split
  across the 180 s call limit. 8 new tests (426 total, 3 skipped).
- `detect_ripples.py` gains `--no-stress` / `--stress-only`. The stress-pair sensitivity
  re-detects on both channels of every shank where the two rankings disagree;
  Cicero_09012014 has five such shanks and the analysis alone takes 2m45s. Splitting it
  changes what runs when, not what is computed.
- Multi-shank support is now reported as a fraction as well as a count, since the
  denominator ranges from 10 to 12 across these sessions and would otherwise make "7 of
  10" and "7 of 12" look like the same number.

**Channel selection replicates.** Both Cicero sessions pick **channel 68 on shank 8** —
the same channel, nine days apart, as Achilles picked the same channels a week apart.
Margins are slightly wider than Achilles (1.03–1.16 against 1.01–1.09). The EMG warning
from D6.12 is much sharper here: in Cicero_09012014 the EMG channel scores **35.1** of
ripple-band envelope, higher than six of the eleven shank winners (29.07–33.88). Had the
non-probe channels not been excluded, the EMG would have been selected over half the
probe.

**Cross-check replicates.** 67–75 % of ripples have an MUA partner against 46–57 % the
other way; population firing leads the ripple envelope by 2.9–7.8 ms; and agreement rises
monotonically with multi-shank support (31–36 % at ≤ 2 shanks, 78–83 % at ≥ 9) in both
sessions, exactly as on Achilles.

**Two session-level surprises.**

- **Cicero_09102014 has more non-REM in POST than PRE** (138.9 against 92.8 min) — the
  opposite of both Achilles sessions, where POST had roughly half of PRE. The PRE/POST
  normalisation by searched time matters in both directions, and this session is the
  check that it is not quietly encoding "POST is shorter".
- **Cicero_09012014 has 7.1 minutes of searchable MAZE immobility**, against 0.4 minutes
  on each Achilles session, and 279 MAZE events against 20. Its MAZE epoch is 92.7 min,
  more than twice the others. D6.3's immobility problem is therefore not universal: it is
  severe where the MAZE epoch is short, and this session gives the first MAZE column with
  enough events to be worth looking at.

**The scoring floor behaves as 8a predicted, and the limit session shows the limit.**

| Session | Place cells | Spikes/bin in events | Empty bins | Scored PRE / POST |
|---|---|---|---|---|
| Achilles_10252013 | 77 | 5.55 | 1 % | 1646 (23 %) / 825 (21 %) |
| Achilles_11012013 | 73 | 4.88 | 0 % | 695 (15 %) / 605 (18 %) |
| Cicero_09102014 | 36 | 3.28 | 3 % | 562 (15 %) / 836 (14 %) |
| Cicero_09012014 | 16 | 1.84 | 28 % | 310 (6 %) / 435 (9 %) |

The 8a prediction was that in-event density would be 1.6–2.4× the running density and the
scorable fraction would fall; both held. Cicero_09012014 scores 6–9 % of its events, but
310 and 435 events are still enough to test.

**Every comparison so far is POST > PRE in sign — sixteen of sixteen** across four
sessions, two shuffles and two event sets. Significance is not uniform: Cicero_09012014's
full-set tests are p = 0.091 and 0.15 while its conjunction tests are p = 0.005 and
0.0072, and Cicero_09102014's conjunction column-cycle test is p = 0.61 with a
rank-biserial of +0.019.

**The controls are no longer unanimous.**

| Session | shuffle | real POST−PRE | identity-shuffled | ratio |
|---|---|---|---|---|
| Achilles_10252013 | column cycle | +0.302 | −0.005 | — |
| Achilles_10252013 | time permutation | +0.220 | −0.033 | — |
| Cicero_09012014 | column cycle | +0.202 | **−0.369** | — |
| Cicero_09012014 | time permutation | +0.185 | −0.114 | — |
| Achilles_11012013 | column cycle | +0.197 | +0.049 | 0.25 |
| Achilles_11012013 | time permutation | +0.433 | +0.071 | 0.16 |
| **Cicero_09102014** | column cycle | +0.209 | **+0.087** | **0.42** |
| **Cicero_09102014** | time permutation | +0.230 | **+0.060** | **0.26** |

The automated verdict calls Cicero_09102014 a failure, on a threshold (residual under a
quarter of the real gap) that we invented in step 7 and that has no external standing.
What is not arbitrary is the pattern: **the two linear sessions leave zero or negative
residual, and the two circular sessions leave +0.05 to +0.09.** That is consistent with
D8.4 — on a ring the score maximises over 28–29 position offsets, which gives the
observed data more freedom to find apparent structure that permuting cell identities does
not remove. It is a reason to weight the circular sessions' PRE/POST differences less, and
a reason the geometries must not be pooled.

**The positive control degrades with ensemble size**, which is the honest reading of the
limit session. Median rZ for 2 s running windows: Achilles +13.2 and +7.5,
Cicero_09102014 +5.6, **Cicero_09012014 +1.8** with only 52 % of windows above 1.65. With
16 place cells the scoring can still recognise a known trajectory, but barely.

**Blocked: Gatsby_08022013.** Its XML says 141 channels; the `.eeg` is 134. See D8.6. No
work done on it pending a decision.

**Time**: ~2 h.

## 2026-10-19 — Step 8b (part 2) and 8c: the fifth session, and the five-session comparison

**Gatsby_08022013 unblocked (D8.7).** Both changes D8.6 identified were authorised and
made. The channel count for this one session comes from the data (134) rather than its XML
(141), declared in `config/params.yaml` as `lfp.channel_count_override` and recorded on
every `LfpMetadata` as `channel_count_source`; the stride evidence that pins it to a single
value is in D8.6. The file-size check is now asymmetric — a recording *shorter* than
`sessDuration` is still fatal, a *longer* one is accepted — because every window this
project reads lies inside an epoch and every epoch ends at or before `sessDuration`.

That second change was made while looking at a failing session, which is worth saying
plainly rather than burying. What bounds it: Gatsby is the only long file (+1,588.8 s, the
known D3.1 anomaly); the other four match the session clock to 0.0 s, so the check's input
and verdict are unchanged for them; and the change touches no threshold, no scoring
parameter and no epoch boundary. 11 new tests in `test_lfp.py`.

**A coverage gap, stated rather than discovered later.** Gatsby_08022013 turns out to have
no `skip="1"` flags at all — it excludes channel 59 by leaving it out of shank 8's
`spikeDetection` group, which is the mechanism channel selection actually uses. The two
sessions whose XML does carry skip flags are both omitted for storage. **No session in this
analysis exercises the skip-flag path**; it is covered by unit tests and nothing else. Now
in METHODS.md.

**Gatsby in numbers.** 30 place cells, 8 of 16 spike groups carrying units (support
denominator 7, the lowest in the set), 671 PRE and 408 POST scorable events.

**8c: five sessions, 20 comparisons, all positive.** Every session shows POST > PRE under
both shuffles and both event sets — 20/20 — with 3, 3, 4 and 4 of the five significant at
p < 0.05 in the four cells. Normalising by searched non-REM minutes gives the same 20/20.
Pooled magnitudes within geometry: linear +0.137 to +0.298 (rb +0.083 to +0.120, p ≤ 4.2 ×
10⁻⁶), circular +0.177 to +0.333 (rb +0.050 to +0.126, p 0.045 to 4.2 × 10⁻⁷).

**The clean geometry story from part 1 did not survive the fifth session (D8.8).** With
four sessions the identity-shuffle residual looked like a circular-track artefact: linear
sessions left nothing, circular left a sixth to a quarter. Gatsby_08022013 is linear and
leaves **+0.121 against a real gap of +0.195 — 62 %, the largest in the set**. Both
circular sessions still leave +16 % to +42 % in all four of their readings, so the circular
caution stands; but geometry is not sufficient as an explanation, and whether the residual
and the circular rZ compression of D8.4 share a cause is now recorded as an open question
rather than guessed at. The residual is printed beside every effect size and quoted in the
README; the invented 25 % pass/fail threshold is gone from `validate_sequences.py`.

Two ratios in that table are uninterpretable and are flagged as such: Gatsby's −180 % and
Cicero_09012014's −183 % are percentages of a near-zero or oppositely-signed denominator.

**Sensitivity floor as a figure in its own right.** `scripts/sensitivity_floor.py` and
`plots.plot_sensitivity_floor` write `results/sensitivity_floor.png` and `.csv`: median rZ
for 2 s windows of real running against place-cell count — +13.2 at 77 cells, +7.5 at 73,
+5.6 at 36, +3.3 at 30, +1.8 at 16, where half of that last session's windows fall below
z = 1.65. Circular sessions are drawn hollow because rZ is not comparable across geometry
(D8.4), and ensemble size and geometry are confounded in these five sessions, so the figure
is a floor rather than a calibration curve. 2 new tests.

**Other findings worth keeping.** Cicero_09102014 has *more* non-REM in POST than PRE
(138.9 vs 92.8 min), reversing both Achilles sessions — the "POST is half of PRE"
observation from step 7 is an Achilles fact, not a dataset fact. Cicero_09012014 has 7.1
min of searchable MAZE immobility and 279 MAZE events, so D6.3's MAZE immobility problem is
not universal either; MAZE remains underpowered in the other four (1–5 scorable events).

**Write-up.** README and METHODS.md now carry the five-session result in place of the
two-session one, name the three omitted sessions and the storage reason wherever the scope
appears, and state the skip-flag coverage gap.

**A reporting bug found while checking the write-up.** `results/replay_comparison.csv`
rounded every column to five decimals, so a p-value of 3.4 × 10⁻⁶ was written as `0.0` —
73 of the 90 rows of `unmatched_mua_hypotheses.csv` were affected the same way. The printed
report was right and the saved table was not, which is the worse of the two failure modes.
`results.round_for_csv` now rounds p-value columns to three significant figures and
everything else to the fixed places, and `write_session_table` uses it, so every table
written through it is fixed at once. 4 new tests; the affected tables were regenerated.

**Tests**: 442 passing, 3 skipped. **Time**: ~2.5 h.

## 2026-10-20 — Two additions before pushing: the residual against ensemble size, and the silent failures as a section

**Does the residual track cells rather than geometry? On one shuffle, yes — suggestively
(D8.10).** `scripts/residual_vs_ensemble.py` and `plots.plot_residual_vs_ensemble` plot the
identity-shuffle residual against place-cell count, hollow markers for circular tracks, two
panels (signed fraction and absolute residual), **no curve fitted**. The signed fraction
rises across the four largest ensembles — 77, 73, 36, 30 — which alternate linear,
circular, circular, linear, so cell count orders them and geometry does not; the
16-cell session then reverses the sign. The *absolute* residual is monotone across all five
(0.005, 0.049, 0.087, 0.121, 0.369) on the column-cycle shuffle, in both event sets. Read
as noise rather than bias, that is the sensitivity-floor story again: the negative control
weakens as the ensemble shrinks, and at 16 cells it is noisier than the effect it controls.
Time permutation shows neither pattern, which fits — column cycling is the shuffle that
rewards positional localisation.

Five points, with geometry and ensemble size confounded, establish nothing. The script
reports how far each series is monotone and refuses to go further; `leading_monotone_run`
in `hc11.results` is the whole statistic, and it is four lines long. The geometry
explanation is not retired: both circular sessions still sit above both large-ensemble
linear ones.

Also recorded: the residual is computed once per session and shuffle on the **full** event
set, so the conjunction rows carry the same residual over a different denominator. That is
another reason the absolute residual is the series to read.

**The silent failures are now a section of their own in METHODS.md**, not a list inside the
controls section. Lead paragraph on what they share — every one returned a plausible number
instead of an error — then a table of what each produced and what caught it, then a
paragraph each. The point the table makes is that in all five cases the detector was a
cross-check against something the result had to agree with (a second measurement, a
physical expectation, a control, a document), never the code failing. Several produced
results that *supported* the hypothesis under test, which is the dangerous shape.

Stale text in `replay_comparison.py`'s docstring fixed on the way: it still said the
residual was "a sixth to a half" on circular tracks and that PRE holds twice the non-REM of
POST "in both sessions".

**Tests**: 448 passing, 3 skipped (6 new: 4 for `leading_monotone_run`, 2 for the figure).
**Time**: ~1 h.

## 2026-10-21 — Step 9b: the firing-rate split, and the control that ate it

**Why this ran first.** Phase 3 is the paper's Fig. 2–4 material: a per-cell contribution
(PCC) score, a rigid/plastic split at the median of its PRE-to-POST change, and the claims
that follow. The scoping report (9a) found that Fig. 4D — slow-firing cells increase their
sequence contribution from PRE to POST, fast-firing cells do not — needs no PCC at all, is
fully specified in the main text, and costs a minute of compute. It is the paper's own
independent confirmation of the PCC result, so it is the cheapest decisive thing to run
before building anything.

**New code.** `hc11.subsets` (off-line firing rate over intervals, equal halves with the
median cell dropped, random halves, split agreement, NaN-safe pairing and empirical p) with
18 tests, and `scripts/firing_rate_split.py`, which rebuilds the decoding template from one
half of the ensemble at a time and re-scores the same events. Four provisional parameters
in `config/params.yaml` under `subsets`, reasoning in D9.1.

**The halves land where the paper's do**: median off-line rate 0.34–0.50 Hz slow against
0.81–1.39 Hz fast, versus the paper's 0.39 and 1.12 Hz. 67–88 % of cells keep their half
under the non-REM-only rate definition.

**The result (D9.2).** The slow half gains more than the fast half in 7 of 10
session-by-shuffle cells, including both shuffles of both Achilles sessions. But the fast
half also gains in 8 of 10 — so the "but not fast-firing" half of the claim does not hold —
and, decisively, **in 9 of 10 cells the slow-minus-fast gap sits inside the spread of gaps
produced by splitting the same ensemble at random**. The one exception is p = 0.095.
Rebuilding a template from half an ensemble is a large perturbation that produces gaps of
this size routinely; the firing-rate split does not produce a larger one.

Reported as a non-replication of Fig. 4D as stated. Nothing was tuned and no variant was
searched for one that would separate.

**A sixth silent failure, caught by reading the line above the answer.** Cicero_09102014's
time-permutation gap came back as `p = 0.048` — the smallest value a 20-draw null can
express, and a result that would have been the strongest in the table. It was NaN: one
event in the slow half had a degenerate posterior, so its null had zero width, its z-score
was undefined, and the median inherited it. `NaN >= draw` is false for every draw, so the
count of more-extreme draws was zero. A missing measurement had become the most significant
number on the page. Fixed by `subsets.empirical_p` (NaN in, NaN out; non-finite draws
discarded) and `subsets.finite_in_both` (an event either half cannot score leaves the
comparison for both halves, so the two are never compared on different event sets). Both
unit-tested on exactly these cases, and added to the METHODS catalogue, which is now six.

While adding it I corrected an overstatement of my own from the step-8c write-up: that
section said three of the five failures produced results supporting the hypothesis under
test. Counting them honestly it was two of five, and it is three of six now that the NaN
comparison joins them. The names are in METHODS.md so the claim can be checked.

**Stopping here.** The brief says proceed to 9c if 9b holds. It does not hold, so no PCC
code was written. What 9b does hand the PCC phase is a calibrated expectation: if dividing
38 cells by a threefold firing-rate difference produces a gap indistinguishable from a
random division of the same cells, the null-label control for the rigid/plastic split is
very likely to produce effects the size of the real ones — which is why the brief puts it
first, and why it should be run before the real labels are ever looked at.

**Tests**: 466 passing, 3 skipped. **Time**: ~2 h.

## 2026-10-24 — Step 9c: per-cell contribution on Achilles_10252013, null labels first

**Provenance, because it matters for a repository whose history is part of the point.**
The estimator (`src/hc11/pcc.py`), its script and its tests were **drafted in a separate
session on 23 October** and were sitting uncommitted in the working tree at the start of
this one. They were adopted deliberately rather than silently: every line was reviewed
before anything was run, five changes were made (below), and the commit says so.

**What the review changed.**

1. *The real template is now scored through the same code path as every shifted one.* The
   draft took the unperturbed mean from the step-7 table and the shifted means from its own
   loop. Any difference between those paths would have landed in every cell's PCC in the
   same direction, undetectably. Placement 0 now goes through the same function, and step 7
   serves as an independent check instead: agreement is 1.64e-05 rZ, the step-7 CSV's own
   rounding floor, and the run aborts above 1e-3.
2. *One event mask across all placements* (`pcc.finite_everywhere`), so a contribution is
   never a difference between event sets. The same bug shape 9b found a week ago.
3. *The silent-cell case is tested, not assumed.* Two tests assert that shifting a silent
   cell's field changes the score, and that its contribution is not structurally zero. It
   was already handled correctly — the field enters the decoder's normalisation term — but
   "appears to handle it" is the standard that produced six silent failures.
4. *"Strongly contributing" made signed* rather than two-sided, with the two-sided count
   logged beside it (30 against 38).
5. `np.rint` instead of `int()` on a float bin index in the validation path.

**The null-reuse approximation passes.** Validation raised from 600 to **40,000 triples**
(68 s): correlation 0.996, overall bias +0.0004 ± 0.0006 rZ, and — the number that matters,
since a bias common to both epochs cancels in ΔPCC — epoch-differential bias bounded at
**|bias| < 0.0023 rZ**. What it cannot bound is a per-cell bias below ~0.01 rZ; that is
recorded in D9.4 along with why it does not threaten this particular conclusion.

**The result (D9.5): the rigid/plastic split is reproducible from noise.**

- Real ΔPCC against null-label ΔPCC: Kolmogorov–Smirnov p = 0.39 and 0.50, Mann-Whitney
  p = 0.090 and 0.104. The real distribution is a third wider and shifted +0.001. **Not
  distinguishable.** Per cell, 10 and 12 of 77 fall outside their own null at p < 0.1,
  against 7.7 by chance.
- Every rigid/plastic property difference is smaller than what the null labels produce at
  their largest, and on peak rate the null labels reach p < 0.05 in **6 of 15 draws** where
  the real labels reach p = 0.21.
- The paper's claims: rigid cells firing faster — absent, wrong sign. Plastic cells more
  spatially specific — right direction, inside the null envelope. Plastic cells with fewer
  fields — absent. PRE-to-POST contributor overlap **38 % and 29 % against their 73 %**.

**What the estimator does get right**, so the failure is not overread: the 77 per-cell
contributions sum to +0.325 and +0.296 rZ, against the population POST−PRE differences of
+0.302 and +0.220 measured in step 8c. It is additive and it recovers the effect it
decomposes. It is the per-cell split it cannot support.

**Calibration carried in from 9b**, as the brief asked: a threefold rate split across 38
cells could not be separated from a random one (gaps +0.148, +0.208 against random-split
95th percentiles 0.302, 0.233), so a population claim here needs ~0.25–0.30 rZ to be
evidence; and 15-against-15 median splits of noise produce |rank-biserial| up to 0.43–0.56,
so a property difference below ~0.45 is not evidence either. Every real difference measured
in 9c is below both floors.

Cost, for the record: 77 cells x 16 placements x 2471 events ran in 20 chunks of ~17 s.
The whole estimator is about six minutes of compute per session.

**Tests**: 484 passing, 3 skipped (22 in `test_pcc.py`). **Time**: ~2 h.

## 2026-10-24 (later) — Step 9d: per-cell contribution on all five sessions

The 9c estimator run unchanged on the remaining four sessions, null labels first
throughout. About 10 minutes of compute for the shift means (chunked per session) and
another 5 for 40,000 null-reuse validation triples each. Four additions the brief asked
for, all in D9.6.

**Additivity, established per session rather than assumed.** Per-cell contributions sum to
the step-8c population difference within a factor of two everywhere, and within 15 % on
Achilles_10252013 (1.07) and Cicero_09012014 (1.14). The estimator is not broken; the
loosest ratios are the circular sessions (1.7–2.0), consistent with D8.4's rZ compression.

**The two circular sessions *do* separate from their nulls — in the wrong way.** This was
the surprise, and it is the reason the write-up needs care. Achilles_11012013 and
Cicero_09102014 reject at p = 0.006 and 0.002, but the separation is entirely a shift in
**location**: their spread ratios are 1.07 and 0.99 (Levene p = 0.64, 0.74 after centring).
A constant added to every cell is the population POST>PRE effect showing up per cell, and
it is precisely what a within-session median split throws away, since the split ranks cells
against each other. The quantity the split needs is excess *spread*, and only
Gatsby_08022013 has it (ratio 1.58–1.60, p ≤ 0.0016).

**The overlap claim is a function of ensemble size.** 29–38 % at 77 cells, undefined at 73
and 36 (one qualifying cell each), **77–80 % at 30 and 16 cells** — the paper reports 73 %.
Our numbers bracket theirs, and the match lands in the two smallest ensembles, where the
statistic rests on 13 and 5 cells. A cell in a 16-cell ensemble carries a sixteenth of the
posterior and clears a threshold easily in both epochs; in a 77-cell ensemble it carries a
seventy-seventh and does not. This is the cleanest cross-session finding in Phase 3,
because it needs no control at all.

**Three of five sessions cannot test a property difference.** Null-label |rank-biserial|
ceilings: 0.56 at 15v15, 0.55 at 9v8, but 0.76–0.81 at 5v5 and 4v4, 1.00 at 4v4, and
Cicero_09102014 produces only two strongly contributing cells so there is no split to make.
Marked vacuous rather than listed as further negatives. Of 32 comparisons, **one** exceeds
its own null 95th percentile — Gatsby time-permutation place-field count, +0.639 against a
ceiling of 0.639, so it ties rather than beats, in the cell whose population difference is
+0.014 at p = 0.52, and in the direction opposite to the paper's.

**The per-cell bias bound carried forward.** Achilles_11012013's real ΔPCC spread is only
1.1–1.2 times what an undetected per-cell null-reuse bias could account for; the others run
1.4 to 10.5. No conclusion rests on it — a bias cannot manufacture a negative — but it is
recorded so that session's numbers are never read as if the bound had been closed.

`results/pcc_cross_session.png` puts the two findings side by side: the overlap against
ensemble size with the paper's 73 % marked, and every rigid/plastic effect against its own
null ceiling, 31 of 32 below the diagonal.

**Phase 3 stops here for building.** 9e and 9f would measure an estimator whose resolution
is now bounded below the effects they would look for. Next is the write-up decision.

**Tests**: 484 passing, 3 skipped. **Time**: ~1.5 h.

## 2026-10-25 — Phase 3 write-up, and the public documents restructured around three phases

No new analysis. The documents now carry all three phases rather than Phase 2 with a
Phase 3 appendix.

**Framing first, and fixed in the decision log (D9.7).** Phase 3 is a measurement claim,
not a refutation: 77 cells in the best session against the paper's 491, five sessions of
eight, and a procedure reconstructed from four sentences and a figure caption. The claim is
that this estimator on these ensembles cannot resolve the per-cell structure, with two
measured floors saying by how much. D9.7 also fixes three requirements that keep the
framing honest in the documents themselves — the caveat goes early rather than in a trailing
limitations list, vacuous sessions are not counted as negatives, and additivity is reported
beside every negative so "below resolution" is never mistaken for "broken".

**README restructured around the arc.** A reader should get it in the first fifty lines:
Phase 1 replicated the direction of the published effect (20 of 20), Phase 2 found the
method's sensitivity floor (+13.2 down to +1.8 across ensembles), Phase 3 bounded what the
per-cell estimator can resolve. Phase 3's section leads with the overlap result, because it
needs no control and is the cleanest cross-session comparison in the project: 29 % to 80 %,
rising as the ensemble shrinks, with the two matches to the published 73 % being the two
smallest linear ensembles, on 13 and 5 cells. Then the null labels with the location/spread
decomposition, additivity as the reason the negative is credible, the two floors, the three
vacuous sessions, and Achilles_11012013's unresolved per-cell bias bound.

**METHODS.md gained a stage 6** in the same form as the others — the seven choices, the
parameters, the validation (placement-0 path check at 1.6e-05 rZ, null reuse at 40,000
triples, the two floors), the outcome.

**The silent-failure catalogue stays at six, and gained a near-miss with a different
detection mode.** The placement-0 path error was caught by *reading code*, not by reading
output: it would have shifted every cell's contribution in the same direction, which is the
one shape no downstream control can see. The other six were all caught because something
independent disagreed. The honest count of hypothesis-supporting failures is kept where it
belongs — two of the first five, three of six now — with a note saying why it is spelled
out, since that is the kind of number that drifts upward in retelling.

Tagged **v0.2.0**. **Tests**: 484 passing, 3 skipped. **Time**: ~1 h.

## 2026-10-26 — Step 10a-bis: is the event-length shortfall fixable? Mostly, and not the part that matters

Phase 4 asks whether a receptive-field-free decoder works where place fields ran out of
resolution. Before building one on events of 7–8 bins when Chen et al. scored 11.5, this
step asks whether the shortfall is the data or a parameter.

**It is a parameter, and the config lied about whose it was.** `threshold_low_sd` was
tagged `[CHEN2016]`. Chen et al. state one MUA threshold — peak above mean + 3 SD — and
nothing about boundaries. decisions.md and METHODS.md both had it as provisional all
along; only `config/params.yaml` was wrong, and a parameter believed to be the paper's is
one nobody sweeps. It survived the step-8a audit because that audit asked which parameters
had been *checked on more than two sessions*, not which provenance tags were true. Fixed in
place with the correction noted (D10.1).

**The sweep** (`scripts/boundary_threshold_sweep.py`) moves the boundary from 1.0 to −1.0
SD, leaves the 3 SD peak alone, and reports duration, bins, scorable fraction, active cells,
spike density and the two independent quality checks (ripple match, multi-shank support).
Two guards: at the pipeline's own threshold it must reproduce `results/mua_events.csv`
exactly — it does, 10 of 10 session-epochs on count and median duration, and it raises
rather than writing a table if that ever stops being true — and the ripple match goes
through `events.nearest_within`, now unit-tested, rather than a hand-rolled search.

**What came back (D10.2).** At −0.5 SD the median scorable event runs 9.0–10.5 bins against
6.0–7.0 today and Chen's 11.5, and 1.9–3.3× more events clear the five-bin floor. The extra
duration is not drift into baseline: the ripple-matched fraction is flat to three points
higher and multi-shank support does not move at all. **But active cells per event is flat or
slightly lower in all five sessions** — the extension buys time, not participants — and
spikes per used bin falls 29–34 %.

That is the answer this step existed for. Phase 4 would start from better footing on bins
and *unchanged* footing on cells, and Chen et al.'s detection power scales with the
active-cell count, where we sit at 15–27 against the ~40 their own Fig. 4B needed. So if
10d comes back null, the attribution is the ensemble, not the event length — which is what
we wanted to know before building rather than after.

Two things the sweep found that nobody asked for. Below −0.5 SD **every session returns
zero events**: the rate never falls back, everything merges, and the merged events exceed
the 400 ms cap and are dropped. And in Cicero_09012014 the PRE-minus-POST gap in the
scorable fraction **widens from −4.8 to −14.6 points** as the boundary drops, where the
other four sessions hold steady — a parameter change that helps one epoch more than the
other, in the session with the smallest ensemble. Nothing is adopted; that asymmetry would
have to be understood before anything could be.

**Tests**: 489 passing, 3 skipped (5 new for `nearest_within`). **Time**: ~1.5 h.

## 2026-10-27 — Auditing every citation tag, and a seventh failure mode

D10.1 found one parameter credited to a paper that does not state it. This checked the
other twenty-four tags in `config/params.yaml` against the sentences they cite.

**Five were wrong, in three directions** (D10.3). A duration limit credited to Grosmark &
Buzsáki when the "50 to 400 milliseconds" sentence is Chen et al.'s introduction. The
ripple envelope's 3 SD threshold credited to Chen et al. when their 3 SD sits in a
parenthesis attached to *multi-unit activity* — and decisions.md had said exactly that
since D8.1, so the config contradicted the decision log and the config was the file being
read. A sleep-state label credited to a paper that never describes the scoring, when the
labels come from the released `sessInfo`. A `require_linearized` flag tagged "[GB2016]
implicit", where implicit is not stated. And the opposite error: `min_peak_rate_hz: 1.0`
is **Chen et al.'s** stated criterion and we had been calling it our own choice.

Three more overstate the agreement without being wrong, and were reworded rather than
retagged. The sharpest is place-field smoothing: Chen et al. write "a Gaussian template
(5×1 …) **with a half SD**", which reads as σ = 0.5 bins, and we use 1.0 — twice as wide.
The old comment quoted the "5×1" and omitted the half-SD clause, which is the part of the
sentence that fixes the width and disagrees with us. The value is unchanged; changing it
would move every downstream number and is a decision with its own evidence, not a
correction.

No parameter value moved. Seventeen tags were accurate, including every one the headline
rests on.

**METHODS.md now lists seven silent failures rather than six**, and the lead says plainly
that three different things caught them: six by a cross-check that disagreed, one by
reading code (the placement-0 path, which would have moved every number together), and
this one by reading a tag against its paper. Nothing about a mislabelled tag is detectable
from output — the code is right, the controls pass, the number is reasonable. What it
changes is whether anyone ever questions the number, and the boundary threshold went
un-swept until Phase 4 needed longer events.

**D10.4** records the Cicero_09012014 asymmetry as a condition on future work: the boundary
stays at 1.0 SD, and if a later step proposes loosening it, the PRE/POST scorable-fraction
gap widening from −4.8 to −14.6 points in that one session has to be explained first, on
all five sessions, before it is applied.

**Tests**: 489 passing, 3 skipped. **Time**: ~45 min.

## 2026-10-29 — Step 10b: the receptive-field-free model recovers position, with caveats

**Three corrections first.** `min_peak_rate_hz: 1.0` is Chen et al.'s stated criterion and
METHODS.md's parameter table and D8.1's audit table both called it ours; both corrected in
place (README was already right). The Gaussian smoothing divergence is now a named note in
METHODS.md's place-field stage rather than only a config comment — their "(5×1 …) with a
half SD" reads as σ = 0.5 bins and we use 1.0, twice as wide, unchanged because changing it
re-runs the project. And D10.5 records that `decisions.md` is authoritative for provenance
with `tests/test_params_provenance.py` as the enforcement: all twenty-four external tags
frozen, so no parameter can gain, lose or change an attribution without someone editing a
file whose only purpose is to make them open the paper. The decision log was itself wrong
once, which is why the rule needed a test rather than a designation.

**`scripts/fit_hmm.py` had not landed** before the connection dropped — checked before
rewriting rather than clobbering.

**A control that was not controlling anything.** The first shuffle control rolled each
cell's counts within each bout and scored **+0.47 to +0.52 against the real fit's +0.41** —
a null beating the data. Bouts run a median of six bins and the animal is in one place for
all of them, so a cell's rate within a bout identifies the position whatever the order is;
rolling inside a bout destroys nothing. Replaced with a circular shift along the whole
running timeline, the same idiom as the step-4 place-cell shuffle, which drops the null to
0.05–0.16. Believed as first written, it would have reported that the model learns nothing.

**Claim A passes the go/no-go.** Fitted on 92 pyramidal cells, 403 bouts, 3019 bins,
91,681 spikes: held-out decoding error **10.7–15.7 cm against ~71 cm chance**, median
correlation 0.71 at m = 12, best restarts reaching **0.80–0.87** against Chen et al.'s
0.73–0.79, shuffled-fit null at 0.05–0.16.

**Three caveats that travel with it.** Restarts within one (m, fold) cell differ by up to
**0.52** — all converged, so this is EM finding different optima, and a best-fit report
would have shown 0.85 where the median is 0.61. One fold of three is markedly worse (fold
medians 0.41–0.71, 0.75–0.82, 0.36–0.53), so which third is held out matters more than m
does. And the model dwells **1.5–2.3× longer per state than the animal takes to cross one
position bin** (0.48–0.77 s against 0.33 s), with the shuffled fit at 0.36–0.50 s — so much
of the stickiness is available with no population structure, which is the failure mode 10a
named. Every sleep score in 10d gets read against the shuffled-fit model on the same
events, not against zero.

**Claim A and Claim B agree on running** (D10.8): the behaviour-free spectral ordering
scores within 0.012–0.033 of the behaviour-derived one, on a quantity whose own value is
0.35. That is the licence for scoring sleep events without any behavioural measure.

**A tension that 10c has to confront**: held-out likelihood rises monotonically with m
across the whole grid while Claim A peaks at m = 12. The criterion that will select m does
not agree with the control that says whether the states are places.

**Tests**: 506 passing, 3 skipped (13 for the HMM, 4 for provenance). **Time**: ~2.5 h.
