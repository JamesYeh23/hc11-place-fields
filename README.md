# hc11-place-fields

Replication of hippocampal place-field and Bayesian position-decoding analyses on the
public CRCNS **hc-11** dataset.

## Goal

Reproduce the core spatial-coding results underlying:

- **Grosmark, A. D. & Buzsáki, G. (2016)**, *Diversity in neural firing dynamics supports
  both rigid and learned hippocampal sequences*, **Science** 351:1440–1443.
- **Chen, Z., Grosmark, A. D., Penagos, H. & Wilson, M. A. (2016)**, *Uncovering
  representations of sleep-associated hippocampal ensemble spike activity*,
  **Scientific Reports** 6:32193 — used here as the source of concrete decoding parameters.

The dataset is dorsal CA1 extracellular recordings from 4 rats across 8 novel-maze
sessions, each structured as PRE sleep → MAZE exploration → POST sleep.

**Phase 1** (this repo's current scope):

1. Repository and environment setup
2. A loader for the `*_sessInfo.mat` session files, with the HDF5 schema documented
3. A quality-control (QC) pass across all 8 sessions
4. Directional place fields on the linear-track sessions
5. Memoryless Bayesian position decoding during running, cross-validated

**Phase 2** (not started): ripple-candidate event detection, 20 ms-bin decoding of
PRE/MAZE/POST events, weighted-correlation sequence scores with shuffles, per-cell
contribution scores, and rigid vs. plastic cell classification.

## Status

| Step | Description | State |
|------|-------------|-------|
| 1 | Repository setup | done |
| 2 | Session loader + schema documentation | not started |
| 3 | QC across all 8 sessions | not started |
| 4 | Place fields (linear track) | not started |
| 5 | Bayesian decoding | not started |

## Getting the data

The raw data is **not** in this repository and never will be — it is ~290 MB and is
distributed by CRCNS under its own terms. Download it yourself:

1. Register for a free account at <https://crcns.org/> and request access to the
   `hc-11` dataset (<https://crcns.org/data-sets/hc/hc-11>).
2. Download, at minimum:
   - `NoveltySessInfoMatFiles.tar.gz` — one `<session>_sessInfo.mat` per session
   - `Sessions_Recordings_Summary.pdf` — per-session recording metadata
   - `Channel_Orderings.pdf` — probe/shank channel maps
   The local field potential (LFP) archives (`<animal>_eeg.tar.gz`) are **not** needed
   for Phase 1.
3. Extract so that the session files land in `data/raw/NoveltySessInfoMatFiles/`:

   ```bash
   mkdir -p data/raw
   tar -xzf NoveltySessInfoMatFiles.tar.gz -C data/raw/
   ```

Expected contents (8 sessions, 4 animals):

```
data/raw/NoveltySessInfoMatFiles/
├── Achilles_10252013_sessInfo.mat
├── Achilles_11012013_sessInfo.mat
├── Buddy_06272013_sessInfo.mat
├── Cicero_09012014_sessInfo.mat
├── Cicero_09102014_sessInfo.mat
├── Cicero_09172014_sessInfo.mat
├── Gatsby_08022013_sessInfo.mat
└── Gatsby_08282013_sessInfo.mat
```

### Pointing the code somewhere else

If your data lives elsewhere, set the `HC11_DATA_DIR` environment variable to the
directory holding the `*_sessInfo.mat` files:

```bash
export HC11_DATA_DIR="/path/to/NoveltySessInfoMatFiles"
```

Resolution order is: `HC11_DATA_DIR` → `<repo>/data/raw/NoveltySessInfoMatFiles` →
`<repo>/../hc-11_Data/NoveltySessInfoMatFiles`.

Verify with:

```bash
python -m hc11.config
```

## Installation

This project uses [uv](https://docs.astral.sh/uv/). Python ≥ 3.10.

```bash
uv sync --extra dev        # creates .venv and installs hc11 in editable mode
uv run pytest              # run the test suite
```

Without uv:

```bash
python -m venv .venv && source .venv/bin/activate
pip install -e ".[dev]"
pytest
```

> **Note:** this repository lives inside an iCloud Drive folder. A virtual environment is
> thousands of small files and should not be synced. `.venv/` is gitignored, but if iCloud
> churn becomes a problem, put the environment outside the synced tree:
> `export UV_PROJECT_ENVIRONMENT=~/.venvs/hc11`.

## Layout

```
src/hc11/        importable package (loader, place fields, decoder)
config/          analysis parameters — every tunable number lives here
scripts/         command-line entry points that produce results/
notebooks/       exploratory work only; nothing load-bearing
tests/           pytest suite
docs/            data schema, parameter decisions, worklog
results/         generated figures and tables (gitignored)
data/            raw data (gitignored, see above)
```

## Reproducibility notes

- Every analysis parameter lives in `config/params.yaml`. Nothing is hard-coded in
  analysis functions.
- The full methods of the Science paper are in an online supplement not available to this
  project. Parameters chosen without a published value are labelled **provisional** in
  both the config file and `docs/decisions.md`, with the reasoning.
- Figures are written to `results/<session>/` with a sidecar JSON recording the exact
  parameters and the git commit used to generate them.

## License

Code: MIT (see `LICENSE`). The hc-11 dataset is not covered by this license and remains
subject to the CRCNS data-sharing terms.
