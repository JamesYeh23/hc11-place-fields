# hc11-place-fields

An independent replication of the hippocampal replay analysis in **Grosmark & Buzsáki
(2016)**, built from scratch in Python against the public CRCNS **hc-11** dataset — from
raw MATLAB session files through place fields, a Bayesian position decoder and
candidate-event detection, to the sequence-score distributions of the paper's Figure 1C.

**The direction of the published effect replicates: sequence scores are higher in
post-experience sleep than in pre-experience sleep, in both sessions analysed, under
both shuffle procedures and both candidate-event definitions — eight comparisons out of
eight.** The magnitude is not established: the effect is small, the paper reports no
comparable effect size, and two measurement limitations described below would have to be
resolved before the sizes could be compared at all.

![Cumulative distributions of sequence score by epoch](docs/figures/replay_cdf.png)

*Sequence score (rZ) for candidate events in Achilles_10252013, by epoch. POST sleep
(green) sits to the right of PRE sleep (blue) across the whole distribution, not only in
the tail. The two panels are two different shuffle procedures; both are reported because
they answer different questions (see below). MAZE has only three scorable events and is
shown for completeness.*

---

## Abbreviations

| | |
|---|---|
| **CA1** | *cornu ammonis* field 1, the hippocampal subregion recorded here |
| **LFP** | local field potential — the low-frequency extracellular voltage |
| **MUA** | multi-unit activity — here, the summed firing rate of the recorded population |
| **NREM** | non-rapid-eye-movement sleep |
| **SPW-R** | sharp-wave ripple, a 150–300 Hz LFP transient accompanying population bursts |
| **rZ** | the sequence score: a weighted correlation standardised against within-event shuffles |
| **EMG** | electromyogram — neck-muscle activity, recorded alongside the probe |
| **PRE / MAZE / POST** | the three epochs of each session: sleep, track running, sleep |

---

## What was found

Candidate events were detected in each epoch, decoded at 20 ms resolution against the
place fields, and scored by how well decoded position tracks time within the event. The
score is standardised against shuffles of the same event, giving **rZ**.

Median rZ, PRE against POST. "rb" is the rank-biserial correlation — an effect size
bounded by ±1. "matched-n" repeats the test 200 times with PRE cut at random to POST's
event count.

| Session | Event set | Shuffle | PRE | POST | difference | rb | p | matched-n |
|---|---|---|---|---|---|---|---|---|
| Achilles_10252013 | all | column cycle | +0.548 | +0.850 | **+0.302** | +0.114 | 3.4 × 10⁻⁶ | +0.302 (100 % of draws p < 0.05) |
| Achilles_10252013 | all | time permutation | +0.188 | +0.407 | **+0.220** | +0.131 | 1.2 × 10⁻⁷ | +0.217 (100 %) |
| Achilles_10252013 | conjunction | column cycle | +0.529 | +0.873 | **+0.344** | +0.136 | 5.4 × 10⁻⁷ | +0.336 (100 %) |
| Achilles_10252013 | conjunction | time permutation | +0.188 | +0.438 | **+0.250** | +0.147 | 6.0 × 10⁻⁸ | +0.248 (100 %) |
| Achilles_11012013 | all | column cycle | +0.110 | +0.307 | **+0.197** | +0.071 | 0.026 | +0.197 (77 %) |
| Achilles_11012013 | all | time permutation | +0.291 | +0.723 | **+0.433** | +0.178 | 3.1 × 10⁻⁸ | +0.433 (100 %) |
| Achilles_11012013 | conjunction | column cycle | +0.129 | +0.350 | **+0.221** | +0.091 | 0.0091 | +0.221 (100 %) |
| Achilles_11012013 | conjunction | time permutation | +0.251 | +0.776 | **+0.525** | +0.189 | 6.4 × 10⁻⁸ | +0.516 (100 %) |

A rank-biserial correlation of +0.07 to +0.19 means that **if you draw one PRE event and
one POST event at random, the POST event scores higher 54–59 % of the time.** That is a
real shift and a consistent one, but it is not a large one.

Normalising by searched time rather than by event — POST holds roughly half the NREM
sleep of PRE in both sessions — widens the margin:

| Session | Shuffle | PRE | POST |
|---|---|---|---|
| Achilles_10252013 | column cycle | 2.01 significant events/min | 2.69 /min |
| Achilles_10252013 | time permutation | 0.83 /min | 1.58 /min |
| Achilles_11012013 | column cycle | 0.45 /min | 0.81 /min |
| Achilles_11012013 | time permutation | 0.70 /min | 1.72 /min |

**Why the magnitude cannot be compared to the paper.** Grosmark & Buzsáki report the
Figure 1C distributions graphically and give no effect size, no median rZ, and no
event-count denominator, so there is no published number to compare ours against. Two
further obstacles are ours rather than theirs: only 15–24 % of our candidate events are
long enough to score, and the absolute level of rZ carries a bias that makes
"fraction of events above threshold" meaningless here. Both are described below.

Forward and reverse traversals split evenly (49–52 % forward in every epoch of both
sessions) with no consistent difference in score between them.

---

## How we know it is not an artifact

A sequence score can look impressive for reasons that have nothing to do with replay, so
the result rests on two controls run through the identical pipeline.

### Negative control: place fields detached from the cells that own them

Permuting which cell owns which place field leaves everything else untouched — the same
events, the same bins, the same spikes, the same shuffles — and destroys the only thing
that could carry replay. If the PRE/POST difference survives that, it is a property of
the event sets or of the scoring, not of the hippocampus.

| Session | Shuffle | POST − PRE, real fields | POST − PRE, identities permuted |
|---|---|---|---|
| Achilles_10252013 | column cycle | **+0.302** | −0.005 |
| Achilles_10252013 | time permutation | **+0.220** | −0.033 |
| Achilles_11012013 | column cycle | **+0.197** | +0.049 |
| Achilles_11012013 | time permutation | **+0.433** | +0.071 |

The gap collapses. Three of the four permuted gaps are under a sixth of the real one and
one is zero.

The control also returned something we did not expect and that changes how the result may
be read: **with identities permuted, the median rZ is still +0.59**, not zero. The
column-cycle shuffle spreads posterior mass uniformly across position, so it rewards any
event whose decoded position is merely *localised*, sequence or not — and permuting
identities does not stop a burst of co-firing cells from producing a localised posterior.
So the absolute level of rZ is not interpretable, and nothing in this repository rests on
a count of events exceeding a significance threshold. Only the PRE/POST difference
survives the control, and only that is reported as a result.

### Positive control: score real running the same way

If the scoring cannot recognise a trajectory the animal actually ran, it cannot be
trusted on sleep. Running periods were decoded with leave-one-lap-out templates — fields
built from every lap except the one being scored — and scored identically:

| Session | Window | Shuffle | n | median rZ | above z = 1.65 |
|---|---|---|---|---|---|
| Achilles_10252013 | 2.0 s | column cycle | 62 | **+13.23** | 95 % |
| Achilles_10252013 | 2.0 s | time permutation | 62 | **+10.68** | 94 % |
| Achilles_11012013 | 2.0 s | column cycle | 145 | **+7.50** | 82 % |
| Achilles_11012013 | 2.0 s | time permutation | 145 | **+7.91** | 85 % |
| Achilles_10252013 | 0.2 s | column cycle | 984 | −0.10 | 11 % |
| Achilles_11012013 | 0.2 s | column cycle | 1606 | −0.67 | 2 % |

A known trajectory scores an order of magnitude above anything seen in sleep.

The 0.2 s rows are not a failure, and they are the reason both window lengths are
reported. At roughly 30 cm/s, 200 ms of running covers about 6 cm — less than one 10 cm
position bin. Such a window cannot contain a trajectory however good the decoder is, so a
low score there measures nothing. That sleep events of that duration score at all is the
point: replay is time-compressed, so 60 ms of replay traverses what takes the animal a
second or more. It also fixes the method's ceiling — a sequence is only visible when the
decoded position actually moves between bins.

---

## The pipeline

Each stage has a module in `src/hc11/`, a command-line script in `scripts/`, and its
parameters in `config/params.yaml`. [`docs/METHODS.md`](docs/METHODS.md) is the full
narrative version; this is the two-minute walkthrough.

**1. Load.** `hc11.io` reads the per-session MATLAB v7.3 files into an immutable
`Session` holding spikes, cluster identities, position, epochs and scored sleep states.
The loader never modifies what it reads; anomalies are reported by `hc11.validate` as
severity-tagged issues and resolved in a separate analysis layer (`hc11.preprocess`), so
a reader can always tell data from decision. Every field is documented in
[`docs/data_schema.md`](docs/data_schema.md).

**2. Behaviour and place fields.** `hc11.behavior` and `hc11.laps` define running periods
and segment them into laps; `hc11.fields` builds directional rate maps, smoothing
occupancy and spike counts separately before dividing, and classifies place cells by
spatial information against a circular-shift shuffle, field detection and split-half
stability. **304 place cells** across the eight sessions.

![Sorted rate maps](docs/figures/rate_maps_pos.png)

*Place-field rate maps for one session and one running direction, each cell normalised to
its own peak and the cells sorted by where that peak falls. The staircase is the signature
of a population tiling the track.*

**3. Decoder.** `hc11.decode` implements the memoryless Poisson decoder,
`log P(x | n) = Σᵢ nᵢ log fᵢ(x) − Δ Σᵢ fᵢ(x)`, normalised in log space. Validated by
decoding each lap against fields built from every *other* lap: **median error 4.7–10.6 cm**
across the eight sessions against a chance baseline of 36–71 cm.

![Decoded against true position](docs/figures/decoding_examples.png)

*Decoded posterior against the animal's true position for example laps. The decoder is
cross-validated — the fields used here were built without the lap being decoded.*

**4. Candidate events.** Two independent detectors, sharing one two-threshold rule in
`hc11.events` so their outputs are directly comparable. `hc11.mua` detects population
synchrony from spikes alone; `hc11.ripples` detects 150–300 Hz sharp-wave ripples on a
CA1 pyramidal-layer channel chosen per shank by `hc11.ripple_channels`. The MUA set is
primary — the LFP exists locally for only two of the eight sessions — with the ripple
match carried as a per-event property, so the stricter conjunction of the two is a filter
on a column rather than a second pipeline.

![Agreement with the MUA detector against multi-shank support](docs/figures/cross_check_by_support.png)

*Two independent checks on the same events. The x axis counts how many other CA1 shanks
show the same ripple-band transient; the y axis is the fraction of those events that the
spike-based detector also found. The two agree about which events are weak — 37 % at zero
support rising to 87 % at eleven — which is the main evidence that the low-support events
are a real category rather than an artifact of either method.*

**5. Sequence scores.** `hc11.sequences` bins each event at 20 ms, decodes every bin, and
takes the weighted correlation between decoded position and bin index over the whole
posterior, standardised against 500 within-event shuffles of two kinds. No position is
read off the posterior: a bin whose posterior is spread over half the track contributes
its spread, not a point estimate.

---

## Where this diverges from the paper, and what turned up on the way

This section is the reason to trust the rest. None of it was tuned away.

### We find 304 place cells; the paper reports 491

A factor of 1.62, outside the tolerance set at the start, and it is **not** a threshold we
have mis-set. Relaxing our criteria one at a time:

| Variant | Unique place cells | % of 562 pyramidal cells |
|---|---|---|
| **our primary criteria** | **304** | 54 % |
| field-width cap removed | 346 | 62 % |
| cap removed + spike floor lowered to 20 | 367 | 65 % |
| peak rate ≥ 1 Hz and nothing else | 392 | 70 % |
| Grosmark & Buzsáki | 491 | 87 % |

The fourth row is the decisive one. Chen et al. (2016), describing the same recordings,
state that all pyramidal neurons selected for analysis had a peak firing rate above 1 Hz.
Applying *only* that — no field requirement, no shuffle test, no stability criterion, no
spike floor — still yields 392, not 491. There is no setting of our criteria that reaches
87 % of pyramidal cells.

Part of the gap is traceable. The released `OneDLocation` variable marks *being on the
linearised track*, not running, and how strictly it was curated varies by session: it is
defined for only 7–35 % of each MAZE epoch. Rebuilding everything on our own linearisation
of the raw two-dimensional position over all tracked MAZE time widens the base to 68–98 %
and raises the yield from 47 % to 59 % of pyramidal cells on the linear tracks — one
session doubles, the one where the released mask was thinnest at 2.6 minutes of running.
That is a real contributor, but 59 % is still well short of 87 %, and the remainder is
unexplained. The authors' mask remains the base for every primary result, because it is
theirs.

### Absolute rZ is not interpretable

Covered under the negative control above: identity-permuted fields still give a median rZ
of +0.59 under the column-cycle shuffle, because that shuffle rewards positional
localisation as well as sequence. Threshold-crossing fractions are therefore meaningless
in this repository and are never quoted as a result. Only the PRE/POST difference, which
the control collapses to zero, carries information.

### A single-lap directional template cost one session its entire effect

Achilles_11012013 is a circular maze, run in one direction. Its other direction contains
**exactly one lap of 226 position samples**, against 75 laps for the direction the animal
actually ran. A place-field template built from one lap is noise — and because the score
takes the better of the two directions, a noise template hands the *shuffles* the same
extra freedom, raising the null and lowering every event's score. Removing it moved that
session's median rZ from −0.060 to +0.110 in PRE and from +0.148 to +0.307 in POST: a
shift the size of the effect being measured, from a template nobody would have defended
if asked. Directions now need five laps before their fields are used.

### Only 15–24 % of candidate events are long enough to score

Scoring needs at least five 20 ms bins with spikes in them, which needs a 100 ms event.
Our events have a median duration of 55–68 ms:

| | median duration | ≥ 3 bins | ≥ 4 | **≥ 5** |
|---|---|---|---|---|
| Achilles_10252013 PRE | 68 ms | 62 % | 37 % | **24 %** |
| Achilles_10252013 POST | 62 ms | 54 % | 32 % | **21 %** |
| Achilles_11012013 PRE | 57 ms | 46 % | 25 % | **15 %** |
| Achilles_11012013 POST | 58 ms | 47 % | 28 % | **19 %** |

This is a larger loss than every other filter in the pipeline combined, and it traces to
the event detector's boundary threshold rather than to anything in the scoring. It is
recorded rather than tuned: a three-bin floor would mean fitting a correlation to three
points, which is not an estimate. The same filter applies to PRE and POST, so the
comparison is unaffected; what it limits is any absolute rate.

### POST has less NREM sleep than PRE, in both sessions

168.5 → 89.7 minutes in one session, 104.0 → 73.0 in the other. Part of that is simply a
shorter POST epoch, but the *proportion* also falls in one session (56 % → 37 %) and the
other's REM sleep collapses instead (27 % → 5 % of scored sleep). The likely reason is
mundane — the animal has already slept for some five hours before POST begins — and it
cannot be tested here. Every rate in this repository is therefore normalised by searched
time, and the PRE/POST test is reported both on all events and on PRE subsampled to
POST's event count.

### Quiet wakefulness on the maze cannot be defined from the released position data

Chen et al.'s criterion is speed below 2 cm/s, which yields 0.4 minutes of searchable time
out of 34–45 minutes. This is not a threshold problem: the sub-threshold samples are
turnaround zero-crossings, not pauses, with a median run of 0.18 s. The genuinely quiet
periods are plausibly the 7.6–8.7 minutes of lost tracking, when the animal sits at a
reward site with its head-mounted light out of view — an open hypothesis, not a ruled-out
one, and one the EMG or theta/delta channels of the LFP would settle.

---

## Scope and limitations

- **Two of eight sessions** carry the full analysis. Steps 1–5 (loading through decoding)
  run on all eight; event detection and sequence scoring need the LFP, which is roughly
  9.5 GB compressed per session, and two were downloaded. Both are from the same animal.
- **The Science paper's methods are in an online supplement that was not available.**
  Parameters came from Chen et al. (2016) where that paper states one, and are otherwise
  our own choice, labelled **provisional** in `config/params.yaml` with the reasoning in
  `docs/decisions.md`.
- The two sessions are reported separately as the primary result, pooled only as a
  secondary one. They are different mazes with different place-cell yields, and two
  sessions are an internal replication rather than a sample.
- Nothing here addresses the paper's rigid-versus-plastic cell classification or its
  per-cell contribution scores. The replication stops at Figure 1C.

---

## Reproducing it

### Install

[uv](https://docs.astral.sh/uv/), Python ≥ 3.10:

```bash
uv sync --extra dev
uv run pytest
```

Or with pip: `python -m venv .venv && source .venv/bin/activate && pip install -e ".[dev]" && pytest`.

> This repository lives in a synced folder in its original location, and a virtual
> environment is thousands of small files. `.venv/` is gitignored; if sync churn becomes a
> problem, put the environment elsewhere with
> `export UV_PROJECT_ENVIRONMENT=~/.venvs/hc11`.

### Get the data

The raw data is not in this repository and never will be. Register at
<https://crcns.org/> and request the [hc-11](https://crcns.org/data-sets/hc/hc-11)
dataset, then download:

| File | Size | Needed for |
|---|---|---|
| `NoveltySessInfoMatFiles.tar.gz` | ~290 MB | everything |
| `Channel_Orderings.pdf`, `Sessions_Recordings_Summary.pdf` | small | verifying the probe layout |
| `<session>_eeg.tar.gz` | ~9.5 GB each | ripple detection only (steps 6–7) |

Extract the session files to `data/raw/NoveltySessInfoMatFiles/`, or point at them:

```bash
export HC11_DATA_DIR=/path/to/NoveltySessInfoMatFiles   # → data/raw/... → ../hc-11_Data/...
export HC11_LFP_DIR=/path/to/lfp                        # → data/raw/lfp → ../hc-11_LFP
```

The LFP directory needs both `<session>.eeg` and `<session>.xml` per session, flat or in
per-session subfolders; without the XML the binary cannot be indexed at all. CRCNS ships
each recording as a gzipped tarball, and `hc11.eeg_archive` unpacks it: it reads the
uncompressed size from the tar header so the file-size check can run before the
extraction starts, and it is resumable, because a 12 GB write does not always finish in
one go. Check the setup with `python -m hc11.config`.

### Run it

```bash
uv run python scripts/inspect_sessions.py --warnings   # per-session overview and QC
uv run python scripts/place_fields.py                  # place fields, all 8 sessions   ~40 s
uv run python scripts/decode_position.py               # cross-validated decoding       ~1 min
uv run python scripts/detect_mua_events.py             # population-synchrony events
uv run python scripts/mua_jitter_null.py               # is the synchrony real?
uv run python scripts/select_ripple_channels.py        # LFP: one channel per shank     ~1 min
uv run python scripts/detect_ripples.py                # sharp-wave ripples             ~2 min
uv run python scripts/cross_check_events.py            # ripples against MUA
uv run python scripts/characterise_events.py           # event properties, decodability
uv run python scripts/score_events.py                  # sequence scores                ~15 s
uv run python scripts/score_events.py --identity-shuffle   # the negative control
uv run python scripts/validate_sequences.py            # both controls
uv run python scripts/replay_comparison.py             # the Fig. 1C comparison
```

Everything lands in `results/` (gitignored): one CSV per analysis and per-session figures,
each figure with a sidecar JSON recording the parameters and the git commit that made it.
The whole pipeline is a few minutes on a laptop once the LFP is extracted; extraction
itself is the slow part, at roughly three minutes per session.

### Tests

```
418 passed, 3 skipped
```

The three skips are a parametrised test that does not apply to circular mazes. **Sixteen
test functions, plus all of `tests/test_data_present.py`, are marked `requires_data` and
skip silently when the hc-11 session files are absent** — so a passing run on a fresh
clone with no data is *not* full coverage. Everything else runs against synthetic
fixtures: HDF5 files built to the same schema, synthetic place fields with known peaks,
posteriors with known correlations, and a 544 MB LFP fixture that checks a one-second read
touches under 10 kB.

---

## Layout

```
src/hc11/        the package: loader, behaviour, fields, decoder, events, ripples, sequences
config/          params.yaml — every tunable number, with its provenance
scripts/         command-line entry points; each writes to results/
tests/           pytest suite, synthetic fixtures plus data-dependent checks
docs/            METHODS.md, decisions.md (the decision log), worklog.md, data_schema.md
results/         generated figures and tables (gitignored)
data/            raw data (gitignored)
```

`docs/decisions.md` records every methodological choice with its reasoning, numbered
D1.1 through D7.8 and cross-referenced from `docs/METHODS.md`. `docs/worklog.md` is the
running diary, including what went wrong and when.

---

## Citations

> Grosmark, A. D. & Buzsáki, G. (2016). Diversity in neural firing dynamics supports both
> rigid and learned hippocampal sequences. *Science* **351**(6280), 1440–1443.
> doi:[10.1126/science.aad1935](https://doi.org/10.1126/science.aad1935)

> Chen, Z., Grosmark, A. D., Penagos, H. & Wilson, M. A. (2016). Uncovering representations
> of sleep-associated hippocampal ensemble spike activity. *Scientific Reports* **6**,
> 32193. doi:[10.1038/srep32193](https://doi.org/10.1038/srep32193)

> Grosmark, A. D., Long, J. & Buzsáki, G. (2016). Recordings from hippocampal area CA1,
> PRE, during and POST novel spatial learning. CRCNS.org.
> doi:[10.6080/K0862DC5](https://doi.org/10.6080/K0862DC5)

## License

Code: MIT, see `LICENSE`. The hc-11 dataset is not covered by that license and remains
subject to the CRCNS data-sharing terms.
