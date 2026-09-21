# hc-11 `sessInfo.mat` schema

Derived from walking the HDF5 tree of all 8 session files with `h5py` (step 2,
2026-09-21), **not** copied from the CRCNS data description. Where the description and
the files disagree, the files win and the disagreement is listed at the bottom.

## File format

- MATLAB v7.3 = HDF5 with a 512-byte userblock. Header string in every file:
  `MATLAB 7.3 MAT-file, Platform: PCWIN64, Created on: Fri Nov 25 17:35:00 2016 HDF5 schema 1.00`.
- Datasets are gzip-compressed, chunked along the long axis.
- **All 8 files have an identical tree**: same groups, same fields, same dtypes, same
  `MATLAB_class` attributes. Only array lengths differ.
- **No** field is a MATLAB empty (`MATLAB_empty` attribute absent everywhere) and **no**
  field uses HDF5 object references. The loader still handles both (tested on synthetic
  files) in case a future file revision introduces them.
- One top-level group, `/sessInfo` (`MATLAB_class='struct'`), with three struct fields:
  `Spikes`, `Epochs`, `Position`.

### Orientation

HDF5 stores MATLAB arrays with axes reversed; the MATLAB array is always `h5py_array.T`.
The loader applies this unconditionally rather than guessing from shape (a guess would
mis-read a 2 × 2 interval matrix). Note the fields are **not** stored in a consistent
orientation relative to each other — e.g. `TimeStamps` is a MATLAB row vector while
`OneDLocation` is a column vector — so field-by-field normalisation is required.

In the tables below, *on-disk* is the `h5py` shape, *MATLAB* is the logical shape, and
*loaded* is what `hc11.io.Session` exposes. `N` = number of spikes, `M` = number of
position samples, `k` = number of intervals.

## `sessInfo.Spikes`

| Field | On-disk | MATLAB | dtype | Loaded as | Unit | Meaning / quirks |
|---|---|---|---|---|---|---|
| `SpikeTimes` | (1, N) | N × 1 | float64 | `spike_times` (N,) float64 | s | Session clock; PRE starts at 0. Globally sorted. Resolution 50 µs (20 kHz). Many exact ties across clusters (36k–455k per session, expected at 20 kHz with 50–137 units); **no** duplicate spikes within a cluster. Minimum within-cluster ISI is exactly **0.85 ms in every session**, suggesting a refractory de-duplication step was applied upstream. |
| `SpikeIDs` | (1, N) | N × 1 | float64 | `spike_ids` (N,) int64 | — | Cluster ID per spike. All whole numbers; cast after checking. `shank = ID // 100`, `cluster_in_shank = ID % 100`. No cluster 0 (noise) or 1 (multi-unit) present. |
| `PyrIDs` | (1, n_pyr) | n_pyr × 1 | float64 | `pyr_ids` sorted int64 | — | Putative pyramidal cells. |
| `IntIDs` | (1, n_int) | n_int × 1 | float64 | `int_ids` sorted int64 | — | Putative interneurons. |

Verified in all 8 sessions: `PyrIDs ∩ IntIDs = ∅`; `PyrIDs ∪ IntIDs` = exactly the set of IDs
in `SpikeIDs` (**zero unclassified clusters**); every labelled unit has spikes.
Derived: `shank_of_cluster` (mapping ID → shank). Shanks in use range 1–15; not every
shank carries a unit in every session.

## `sessInfo.Position`

| Field | On-disk | MATLAB | dtype | Loaded as | Unit | Meaning / quirks |
|---|---|---|---|---|---|---|
| `TimeStamps` | (M, 1) | 1 × M | float64 | `position_t` (M,) | s | Exactly **25.6 ms** spacing (39.0625 Hz) with **no gaps** — min = max = 25.600 ms in every session, i.e. tracking was resampled onto a uniform grid. First sample = MAZE start exactly; last = MAZE end − 25.6 ms. |
| `TwoDLocation` | (2, M) | M × 2 | float64 | `position_xy` (M, 2); `position_xy_cm` | m | Camera coordinates, columns (x, y). NaN in both columns together when tracking fails (never only one). Arbitrary origin; the track runs diagonally in some sessions. |
| `OneDLocation` | (1, M) | M × 1 | float64 | `position_1d` (M,); `position_1d_cm` | m | Linearised position. **Defined only during track running** — NaN for 42–93 % of samples, most of them with valid 2-D tracking (see below). Linear: range exactly [0, L] with L = 1.6 or 2.0 m, and an exact linear projection of `TwoDLocation` (\|r\| = 1.0000 against the first principal axis). Circular: [0, 2.84–2.90] m, 0 = reward site. |
| `MazeType` | (n_char, 1) | 1 × n_char | uint16 | `maze_type` str | — | `MATLAB_class='char'`, `MATLAB_int_decode=2`: UTF-16 code units. Values: `'1.6m Linear Maze'` (4), `'2m Linear Maze'` (1), `'Circular Maze'` (3). |

Derived on `Session`: `maze_kind` (`'linear'`/`'circular'`), `track_length_m` (parsed from
the name; `None` for circular), `position_dt`.

### `OneDLocation` is a running mask, not just a tracking-loss mask

| Session | Maze | 2-D NaN | 1-D NaN | 1-D NaN with valid 2-D | 1-D defined | Median speed, 1-D defined / undefined |
|---|---|---|---|---|---|---|
| Achilles_10252013 | 1.6 m linear | 21.4 % | 87.2 % | 65.8 % | 4.4 of 34.5 min | 53.0 / 11.8 cm/s |
| Achilles_11012013 | circular | 18.5 % | 41.5 % | 23.0 % | 26.1 of 44.6 min | 21.2 / 11.2 cm/s |
| Buddy_06272013 | 1.6 m linear | 2.8 % | 92.8 % | 89.9 % | 2.8 of 38.8 min | 54.8 / 11.4 cm/s |
| Cicero_09012014 | 1.6 m linear | 2.6 % | 78.2 % | 75.6 % | 20.2 of 92.7 min | 14.4 / 8.2 cm/s |
| Cicero_09102014 | circular | 53.7 % | 71.3 % | 17.6 % | 14.5 of 50.5 min | 20.7 / 6.7 cm/s |
| Cicero_09172014 | 2 m linear | 32.4 % | 67.9 % | 35.4 % | 16.5 of 51.3 min | 14.5 / 4.5 cm/s |
| Gatsby_08022013 | 1.6 m linear | 1.7 % | 81.3 % | 79.6 % | 7.3 of 38.9 min | 42.3 / 21.7 cm/s |
| Gatsby_08282013 | circular | 4.5 % | 63.7 % | 59.1 % | 16.3 of 44.9 min | 27.6 / 9.9 cm/s |

`OneDLocation` is never defined where `TwoDLocation` is NaN. Speed is from finite
differences of unsmoothed 2-D position, so the absolute values are noisy, but the contrast
is consistent: the authors appear to have already restricted `OneDLocation` to
track-running segments. This matters for step 4 — see open questions in the worklog.

## `sessInfo.Epochs`

| Field | On-disk | MATLAB | dtype | Loaded as | Unit | Meaning / quirks |
|---|---|---|---|---|---|---|
| `PREEpoch` | (2, 1) | 1 × 2 | float64 | `epochs['PRE']` (2,) | s | Always starts at 0. |
| `MazeEpoch` | (2, 1) | 1 × 2 | float64 | `epochs['MAZE']` (2,) | s | Starts exactly where PRE ends. |
| `POSTEpoch` | (2, 1) | 1 × 2 | float64 | `epochs['POST']` (2,) | s | Starts exactly where MAZE ends; ends at `sessDuration`. |
| `sessDuration` | (1, 1) | 1 × 1 | float64 | `sess_duration` float | s | Total concatenated duration. |
| `Wake` | (2, k) | k × 2 | float64 | `states['Wake']` (k, 2) | s | Active waking. Covers 96–100 % of the MAZE epoch. |
| `Drowsy` | (2, k) | k × 2 | float64 | `states['Drowsy']` (k, 2) | s | Drowsy wake / light sleep. |
| `NREM` | (2, k) | k × 2 | float64 | `states['NREM']` (k, 2) | s | Non-rapid-eye-movement (non-REM) sleep. |
| `Intermediate` | (2, k) | k × 2 | float64 | `states['Intermediate']` (k, 2) | s | Short NREM→REM transition state with high spindle power. |
| `REM` | (2, k) | k × 2 | float64 | `states['REM']` (k, 2) | s | Rapid-eye-movement sleep. |

State matrices: boundaries are **whole seconds**; rows sorted; no state overlaps another
state or itself in any session; together they cover 95.8–99.4 % of each session (short
unscored gaps). k per state ranges 6–72.

Epoch boundaries and cell counts match `Sessions_Recordings_Summary.pdf` for all 8
sessions to within 1 s (pinned in `tests/test_io.py::RECORDING_SUMMARY`).

## Anomalies (pinned in `tests/test_io.py::KNOWN_WARNINGS`)

1. **Gatsby_08022013 — data beyond the end of the session.** `sessDuration` =
   POST end = 30 413.628 s, and the recording summary agrees (POST = 10 999 s). But
   **263 976 spikes (4.67 %)** from all 80 clusters fall between 30 413.63 and 32 002.37 s,
   i.e. up to **1 589 s (26.5 min) after the declared end**, and 5 state intervals
   (Drowsy ×2, NREM, Intermediate, REM) also extend to 31 473–31 996 s. The loader keeps
   these values untouched; `validate()` flags them.
2. **Achilles_11012013 — zero-length REM interval**: row 12 is `[13281, 13281]`.
   Harmless for duration sums; could matter for code that assumes `start < end`.

## Discrepancies with the CRCNS data description (v0.8, 7 Dec 2016)

| # | Description says | File shows |
|---|---|---|
| D1 | The v7.3 files "can be opened using the Python function 'scipy.io'". | False: `scipy.io.loadmat` does not read v7.3/HDF5. Read with `h5py` (or `pymatreader`). |
| D2 | `OneDLocation` is a "linearization of the track-running". NaN is described only for `TwoDLocation` (tracking loss). | Consistent in wording, but the practical consequence is unstated: `OneDLocation` is NaN for 42–93 % of MAZE samples, most of which have valid 2-D tracking. Only 2.8–26 min per session carry a 1-D position. |
| D3 | All times given in seconds within `sessDuration` (implied). | Gatsby_08022013 has spikes and state intervals up to 1 589 s beyond `sessDuration` (anomaly 1). |
| D4 | Circular maze is "1 m diameter". | 1-D range tops out at 2.84–2.90 m, below the π m ≈ 3.14 m a 1 m circle implies; the median 2-D radius (0.45–0.53 m) implies 2.82–3.34 m. Not resolved; circular sessions are out of scope for Phase 1 place fields. |
| D5 | Tracking at "~39 Hz". | Exactly 39.0625 Hz (25.6 ms), gap-free — resampled. Consistent, just more precise. |
| D6 | Summary PDF named `HC-11_Recording_Summary.pdf`, channel map `HC-11_Channel_Orders.pdf`. | Downloaded files are named `Sessions_Recordings_Summary.pdf` and `Channel_Orderings.pdf`. |

Shapes in the description's MATLAB listing (`TwoDLocation [105180x2]`,
`TimeStamps [1x105180]`, …) **do** match the files exactly; its worked example is
Gatsby_08282013 (5 775 679 spikes, 41 pyr, 10 int), and that is pinned as a test.
