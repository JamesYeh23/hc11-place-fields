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

