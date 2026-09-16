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

_Not started._

---

## Step 3 — quality control

_Not started._

---

## Step 4 — place fields

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
