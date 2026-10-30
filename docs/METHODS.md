# Methods

How the pipeline works and why, in the order it runs. This is the narrative version:
[`decisions.md`](decisions.md) is the decision log, numbered D1.1 to D8.9 and referenced
throughout here; [`worklog.md`](worklog.md) is the running diary;
[`data_schema.md`](data_schema.md) documents the raw file format field by field. Results
and effect sizes are in the [README](../README.md).

Every parameter named below lives in `config/params.yaml` and carries a provenance tag
there:

- **[CHEN2016]** — stated in Chen, Grosmark, Penagos & Wilson (2016), *Sci. Rep.* 6:32193,
  which describes the same recordings and is the only published source of concrete
  numbers available to this project.
- **[GB2016]** — stated in the main text of Grosmark & Buzsáki (2016), *Science* 351:1440.
- **[DATA]** — measured from the files themselves.
- **[PROVISIONAL]** — our choice. The Science paper's full methods are in an online
  supplement that was not available, so wherever it would have supplied a number, we
  chose one and recorded the reasoning.

Abbreviations are spelled out at first use in the README; the frequent ones here are LFP
(local field potential), MUA (multi-unit activity), NREM (non-rapid-eye-movement sleep),
SPW-R (sharp-wave ripple) and rZ (the standardised sequence score).

---

## 1. Loading — `hc11.io`, `hc11.validate`, `hc11.preprocess`

### What it does

Each session is one MATLAB v7.3 file, which is HDF5 with a 512-byte userblock in front of
the magic number. `read_matlab` reads it generically and `Session` assembles the hc-11
structure on top: spike times and cluster identities, putative pyramidal and interneuron
ID lists, two-dimensional and linearised position, the PRE/MAZE/POST epoch boundaries,
and the scored sleep states.

Two rules make the rest trustworthy.

**The loader never modifies data (D2.3).** `Session` is immutable, holds exactly what the
file contains, and clips, interpolates and truncates nothing. Anything that needs fixing
is fixed in `hc11.preprocess`, by analysis code, visibly.

**Arrays are transposed unconditionally, never by shape (D2.1).** MATLAB and h5py disagree
about axis order, and a shape-dependent transpose silently does the wrong thing on any
array that happens to be square or one-dimensional. Every array is `h5py_array.T`, always.

`validate()` returns a list of severity-tagged issues rather than raising (D2.4), so a
caller sees every problem at once instead of the first one.

### Parameters

| | Value | Source |
|---|---|---|
| position sample rate | 39.0625 Hz | [DATA] — an exact, gap-free 25.6 ms grid in all 8 sessions |
| spike time units | seconds on the session clock | [DATA] |
| truncate to `sessDuration` | yes | [PROVISIONAL], D3.1 |
| drop zero-length state intervals | yes | [PROVISIONAL], D3.2 |

### Validation

- Cross-checked against an independent MATLAB reader on every field (D2.6).
- Two sessions carry anomalies, both reported as warnings and resolved in the analysis
  layer: Gatsby_08022013 has 263,976 spikes (4.67 %) after `sessDuration`, truncated per
  D3.1; Achilles_11012013 has one zero-length REM interval, dropped per D3.2.
- A finding that closed off a planned quality metric: **the minimum inter-spike interval
  within a cluster is exactly 0.85 ms in all eight sessions** (D3.5). That is not a
  coincidence — a refractory censoring step was applied before release — so
  ISI-violation measures of cluster quality measure the upstream cleaning, not the
  isolation. Unit quality would have to come from waveform data this project does not
  download.

---

## 2. Behaviour and place fields — `hc11.behavior`, `hc11.laps`, `hc11.fields`, `hc11.track`

### What it does

**Running periods** are the intersection of two criteria (D3.4). Stage 1 is the authors'
own `OneDLocation` variable being defined. Stage 2 is speed above 15 cm/s, computed by
differentiating the two-dimensional position and then smoothing — in that order, with a
NaN-aware kernel, so lost tracking does not leak into neighbouring samples.

The two stages are kept separate on purpose. `OneDLocation` marks *being on the linearised
track*, not running, and how strictly it was curated varies between sessions: the speed
filter removes 3.3 % of its samples in one session and 51.9 % in another. In the first
three sessions the mask is short blocks the length of a single traversal; in the other
five it contains blocks of up to 23–59 s with genuine mid-track pauses. Treating stage 2
as a formality would mean accepting a definition of running that silently differs by
session.

**Laps** are segmented on the authors' mask-level continuous runs, *before* the speed
filter is applied, never after (D4.1) — a lap is a behavioural unit, and chopping it at
every deceleration would invent laps. Direction comes from hysteresis on unwrapped
position, with no assumption that directions alternate (D4.2). On a circular track one lap
is one circuit, split at the reward site (D4.3).

**Rate maps** smooth the occupancy map and the spike-count map separately and then divide,
which is the order that keeps a briefly visited bin from inventing a high rate. Bins
visited for under `min_occupancy_s` become NaN rather than 0 Hz, since a 0 Hz estimate
would bias both the spatial information and the decoder.

**Place cells** must satisfy all of: at least 50 spikes while running in that direction, a
peak rate of at least 1 Hz, spatial information significant against a shuffle null at
α = 0.05, at least one detected field, and split-half stability (odd against even laps)
r ≥ 0.3. A cell counts once per session if it qualifies in either direction (D5.1).

The shuffle is a **circular shift in compressed running time** (D4.4): the selected
running samples are numbered 0…n−1 and the spike train is rotated within that numbering,
which preserves each cell's rate and the animal's trajectory while destroying the
spike-to-position relationship. Rotating in wall-clock time instead would mostly move
spikes into periods the cell was not being measured in.

### Parameters

| | Value | Source |
|---|---|---|
| speed threshold | 15 cm/s | [CHEN2016] |
| speed smoothing σ | 0.1 s | [PROVISIONAL] — ~4 samples; suppresses chatter across the threshold |
| bridge short dips / minimum run | 0.2 s / 0.5 s | [PROVISIONAL] |
| spatial bin | 10 cm | [CHEN2016] |
| rate-map smoothing σ | 1 bin | [PROVISIONAL] — **twice Chen et al.'s stated width**; see the divergence note below |
| minimum occupancy | 0.1 s | [PROVISIONAL] |
| place cell: peak rate | 1 Hz | **[CHEN2016]** — "all putative pyramidal neurons selected for analysis had peak firing rate > 1 Hz"; this project called it provisional until D10.3 |
| place cell: min spikes / α / stability | 50 / 0.05 / 0.3 | [PROVISIONAL], D4.6, D5.3 |
| maximum field width | 0.75 of track extent | [PROVISIONAL], D5.2 |

The field-width cap is a fraction rather than a fixed length because a fixed 120 cm was
75 % of the 1.6 m track but only 41 % of the circular one, so the same cell shape was
judged differently by which maze it was recorded on (D5.2).

**A known divergence from Chen et al., stated here rather than buried in a config
comment.** Their methods read: place fields "smoothed using a Gaussian template (5×1 for
one-dimensional …) **with a half SD**". The literal reading of that is σ = 0.5 spatial
bins, i.e. 5 cm at our 10 cm binning. **We use σ = 1 bin — twice as wide.** The value was
chosen by reading the "5×1" as the kernel's full support and inferring σ ≈ 1 from it,
which is a defensible reading of half the sentence and not of the other half. It is left
unchanged because every rate map, place-cell count, decoding error and sequence score in
Phases 1 to 3 was computed with it, so changing it is a re-run of the project rather than
a correction — but a reader comparing our place-field widths or our 304-against-491 cell
count with theirs should know that our maps are smoothed twice as hard as the paper
specifies (D10.3).

### Validation

- **The speed threshold barely changes the place fields** (D4.5), which is the reassuring
  direction: the stage-1/stage-2 split matters for how much data is used, not for what the
  fields look like.
- Lap segmentation on circular tracks was wrong in an early version — laps spanned
  multiple circuits, giving coverage of 2.0 and 3.0 — and was caught by a unit test
  written for the wrap, not by inspection. Fixed by splitting constant-direction runs at
  every crossing of the reward site (D4.3).
- The shuffle degenerates when laps are perfectly periodic: a circular shift becomes a
  rotation of the rate map and preserves the spatial information. Pinned as a test
  (`test_stereotyped_running_weakens_the_shuffle`) rather than hidden.

### Outcome

**304 place cells** across the eight sessions, against the paper's 491. The gap is a
finding, not a tuning target; see *Known limitations* below and D5.4.

---

## 3. Position decoder — `hc11.decode`

### What it does

A memoryless Poisson decoder with a uniform prior, in log space (D5.6):

```
log P(x | n) = Σᵢ nᵢ log fᵢ(x) − Δ Σᵢ fᵢ(x) + const
```

normalised with a log-sum-exp. The `−Δ Σᵢ fᵢ(x)` term is what makes *too few* spikes
evidence *against* a high-rate position, and it is easy to drop by accident; it has its
own test. Tuning curves are floored at 0.01 Hz so that one spike from a cell whose
estimated rate is zero somewhere cannot veto that position outright.

Cross-validation is **leave-one-lap-out** (D5.7), with the templates precomputed once per
session so the chance baseline and the ensemble-size sweep reuse them. The chance baseline
permutes tuning curves between cells, which preserves the shape of the population's
tuning while destroying the correspondence to position (D5.8).

### Parameters

| | Value | Source |
|---|---|---|
| time bins | 250 ms primary, 20 ms secondary | [CHEN2016] |
| rate floor | 0.01 Hz | [PROVISIONAL] |
| prior | uniform | [PROVISIONAL] |

### Validation

Median error on held-out laps, 250 ms bins, each direction against its own fields:

| Session | Cells | Median error | Chance |
|---|---|---|---|
| Achilles_10252013 | 77 | 4.7–5.1 cm | 40–43 cm |
| Achilles_11012013 | 73 | 4.9 cm | 67 cm |
| Buddy_06272013 | 16 | 8.8 cm | 36–45 cm |
| Cicero_09012014 | 16 | 10.6 cm | 38–53 cm |
| Cicero_09102014 | 36 | 4.9 cm | 71 cm |
| Cicero_09172014 | 25 | 9.0–9.3 cm | 49–63 cm |
| Gatsby_08022013 | 30 | 6.6–7.9 cm | 38–43 cm |
| Gatsby_08282013 | 31 | 5.4 cm | 68 cm |

**This is better than the "low tens of centimetres" expected, and it was checked for
leakage before being believed** (D5.9). Three explanations, all of which hold: the
cross-validation checks pass (honest 5.18 cm against a deliberately leaky 4.56 cm;
occupancy and spikes partition exactly between held-in and held-out laps; three separate
leakage tests); the spatial bin is 10 cm and the decoder returns a bin centre, so even
perfect bin identification gives a median error near 2.5 cm; and a 1.6 m linear track with
77 cells is a far easier problem than the two-dimensional open fields where Chen et al.
report 8.5–12.5 cm with 49 neurons.

Two further results worth keeping: the directional split is worth 6–34 % of accuracy
against a direction-merged template (D5.10), and 20 ms bins cost roughly a factor of two
in median error and more than that in the tail (D5.11) — which is why 20 ms is used only
where it must be, inside candidate events.

---

## 4. Candidate events — `hc11.events`, `hc11.mua`, `hc11.lfp`, `hc11.ripple_channels`, `hc11.ripples`

### What it does

Two detectors, deliberately sharing one rule.

**The shared rule (D6.1).** A candidate must peak above a high threshold (3 standard
deviations above the state's own baseline); its boundaries extend outward to where the
signal falls back below a low threshold (1 SD), so duration describes the whole event
rather than the part that cleared the peak; candidates closer than a merge gap are joined;
survivors are filtered on duration. Thresholds are computed **within a brain state, never
session-wide** (D6.2) — baseline population rate and ripple power both differ between NREM
sleep and quiet waking, so one session-wide threshold would hand one epoch systematically
more events and make the PRE/POST comparison meaningless.

**MUA detector** (`hc11.mua`): the smoothed population firing rate of the pyramidal
ensemble, binned at 1 ms and smoothed with a 10 ms Gaussian.

**Ripple detector** (`hc11.ripples`): the 150–300 Hz amplitude envelope of one CA1
pyramidal-layer channel. The envelope is computed once across the whole recording in
overlapping chunks and then sliced by state (D6.16), because filtering each state interval
separately would put a filter transient at the edge of every one and NREM intervals are
often under a minute. Internal chunk boundaries agree with a whole-signal transform to
3 × 10⁻⁵ of the envelope maximum.

**Choosing the ripple channel** (`hc11.ripple_channels`, D6.13): the highest mean 150–300 Hz
envelope over a fixed, deterministic 10-minute sample of scored NREM, on each shank that
carries sorted units. Two things make this harder than it sounds and both are handled
explicitly:

- **Two channel counts, never one variable** (D6.7). `n_channels_total` (134) is the file
  stride and includes the EMG and accelerometer channels; `probe_channels` (128) is the
  selection pool. An off-by-one in the stride does not merely misreport the duration, it
  shifts every frame by one channel so each "channel" becomes a drifting mixture of the
  real ones — plausible amplitudes, plausible spectra, no error anywhere. `check_size` is
  therefore strict rather than advisory (D6.8).
- **Nothing in the hc-11 XML labels the EMG or accelerometer channels.** The only thing
  separating them from the probe is that `anatomicalDescription` contains every channel
  while `spikeDetection` contains the spike-sorted shanks alone. Taking the anatomical
  grouping as the probe definition — the obvious reading — puts the EMG into the
  ripple-channel contest, where it would have outranked two of the fourteen spike groups
  on ripple-band power.

**The handover** (D6.25): `results/candidate_events.csv`, one row per MUA event carrying
`ripple_matched`, `n_supporting_shanks`, `n_active_place` and `decodable`. The strict
conjunction is a filter on a column, not a second pipeline. `n_supporting_shanks` is NaN,
not zero, for an event with no ripple partner — support is measured on the ripple, and an
event without one was never measured.

### Parameters

| | Value | Source |
|---|---|---|
| ripple band | 150–300 Hz | [CHEN2016] |
| peak threshold | mean + 3 SD | [CHEN2016] |
| boundary threshold | mean + 1 SD | [PROVISIONAL] |
| MUA bin / smoothing | 1 ms / 10 ms | [PROVISIONAL] |
| ripple envelope smoothing | 8 ms | [PROVISIONAL] |
| duration limits | 30–400 ms (MUA), 20–400 ms (ripple) | 400 ms [CHEN2016], floors [PROVISIONAL] |
| minimum cell activation | ≥ 5 cells and ≥ 10 % of the ensemble | [CHEN2016] |
| channel profiling sample | 40 windows × 15 s of NREM | [PROVISIONAL] |
| multi-shank support | max envelope z within ±20 ms, supporting at z ≥ 2 | [PROVISIONAL] |

### Validation

- **File geometry.** All five `.eeg` files are a whole number of 268-byte frames with
  remainder zero, and four of the five match `sessDuration` to +0.0 s (D6.11, D8.7). 134 is
  the only plausible channel count that divides any of them exactly, so the stride is
  confirmed by the data and not only by the XML.
- **One session's XML contradicts its own recording, and the data wins** (D8.6, D8.7).
  Gatsby_08022013 declares 141 channels; the file is an exact number of 134-channel frames
  and 130 bytes short of a 141-channel one. A direct stride test settles it: at 134 the
  channels have the statistics of a real probe (per-channel s.d. 85–2543, one visibly hot
  channel, neighbouring contacts correlating at +0.93), while at 141 every channel has the
  same s.d. to within 1 % and every neighbour pair the same correlation to within 0.004 —
  the mis-stride signature described below. The override lives in `config/params.yaml`
  under `lfp.channel_count_override`, applies to that session alone, and every
  `LfpMetadata` records whether its stride came from the XML or from an override. That
  session's recording also runs 26.5 minutes past `sessDuration`, which is the known
  spikes-past-session anomaly of D3.1 seen in the LFP; a recording *longer* than the
  session is accepted because every analysis window ends at `sessDuration`, while a
  *shorter* one remains fatal, since it can be read past its end.
- **The `skip="1"` channel flag is never exercised by this analysis** (D8.7). The parser
  handles it and unit tests cover it, but the two sessions whose XML uses it
  (Buddy_06272013, Gatsby_08282013) are among the three omitted for storage.
  Gatsby_08022013 excludes its bad channel by a different route — channel 59 is absent from
  shank 8's `spikeDetection` group, which is the grouping channel selection iterates — so
  it is correctly dropped without the flag being read. This is a real coverage gap and is
  stated here rather than left to be discovered.
- **Channel grouping against the dataset's own `Channel_Orderings.pdf`** (D6.12, D8.2):
  for the Achilles and Cicero sessions, 12 CA1 shanks of 10 channels, two four-channel
  "extra channel" groups, EMG on 128–130 and accelerometer on 131–133, every line matching.
  Gatsby_08022013 uses a different probe (BUZ64, 16 shanks of 8, of which 8 carry sorted
  units), which is why group membership rather than a fixed geometry drives the selection
  and why multi-shank support is reported as a fraction: the denominator is 7 there and
  10–12 elsewhere. Two independent confirmations: the EMG
  channels carry 9.8–21.8 of ripple-band envelope against 2.4–5.0 for the accelerometer
  (muscle is broadband, an accelerometer is not), and the extra-channel groups carry no
  sorted units and the lowest ripple-band power on the probe.
- **The channel selection tracks anatomy, not noise** (D6.14), on three grounds that are
  not the margin: each shank's depth profile is smooth and single-peaked; theta and delta
  peak several channels *deeper* than the ripple band, which is the CA1 depth signature;
  and seven of twelve shanks select the identical channel in two Achilles sessions a week
  apart. The Cicero pair repeats it nine days apart, both sessions selecting channel 68
  on shank 8 (D8.9).
- The per-shank margins are only 1–14 %, and ranking by excursion (99th percentile over
  median) instead of level moves the choice by one contact on most shanks and by seven on
  one. Both rankings were detected on. Where they differ most, the level ranking is
  marginally better, and the config says so against the prior expectation (D6.17).

### Outcome

Across the five sessions: MUA event rates of 41–46 per minute of NREM, median duration
55–69 ms. **61–80 % of ripples have an MUA partner and 45–58 % the other way**, with
population firing peaking 3–6 ms before the ripple envelope in four sessions (D6.19,
D6.22). Gatsby_08022013 is the exception at −0.3 to +0.5 ms — no measurable lead — which
is unexplained and may follow from its different probe and its much smaller set of
unit-bearing shanks.

---

## 5. Sequence scores — `hc11.sequences`

### What it does

Each candidate event is binned at 20 ms, every bin decoded against the place fields, and
the resulting `(time, position)` posterior scored by a **weighted correlation**: the
posterior is a weight on each `(position, time)` pair and the ordinary Pearson correlation
is taken under those weights (D7.1).

```
E[x] = Σ P(x,t)·x / Σ P        E[t] = Σ P(x,t)·t / Σ P
r    = (E[xt] − E[x]E[t]) / sqrt( (E[x²] − E[x]²)(E[t²] − E[t]²) )
```

No position is read off the posterior. A bin whose posterior is spread over half the track
contributes its spread rather than a point estimate — which matters, because a
maximum-a-posteriori trajectory through broad posteriors is mostly noise.

The score is standardised against 500 within-event shuffles, giving **rZ**. Two shuffle
kinds are run and both reported, because they answer different questions (D7.6):

- **column cycle** rolls each time bin's posterior independently. It destroys the
  position–time relationship and keeps each bin's posterior shape — but it also spreads
  posterior mass uniformly across position, raising the variance of position in the null,
  so any event whose decoded position is merely *localised* scores above it.
- **time-bin permutation** reorders the bins. It preserves the set of posteriors exactly,
  so position variance is unchanged and only the covariance moves. It isolates sequence
  from localisation.

Three maximisations are involved, and **all three are applied identically to the
observation and to every shuffle**, which is the only reason rZ is interpretable:

1. **Field direction** (D7.2). An event is scored against each running direction's
   templates and the larger `|r|` kept. One shuffle realisation is applied to both
   directions and the same maximum taken; shuffling them independently would give the null
   more freedom than the observation had.
2. **Circular position** (D7.4). On the circular maze, position has no origin and a
   trajectory crossing the seam of the linearisation scores near zero despite being a
   perfect sequence. The correlation is evaluated at every circular offset and the best
   kept — "is there an origin that makes this linear in time", the question that has
   meaning on a ring. Linear tracks take the same code path with one offset.
3. The same applies to taking `|r|` rather than signed `r`, so that forward and reverse
   replay do not cancel. Traversal direction is recorded separately as the sign.

**Empty bins are dropped** before scoring (D7.5). A bin with no spikes still has a
posterior — the memoryless decoder returns `−Δ Σᵢ fᵢ(x)`, which favours positions where
the ensemble is quiet — and that describes the place fields, not the animal.

### Parameters

| | Value | Source |
|---|---|---|
| time bin | 20 ms | [CHEN2016] |
| minimum bins kept | 5 | [CHEN2016] |
| shuffles per event per kind | 500 | [PROVISIONAL] |
| minimum laps for a directional template | 5 | [PROVISIONAL], D7.3 |
| significance threshold quoted | z = 1.65 | [CHEN2016] |
| PRE subsamples | 200 | [PROVISIONAL] |

### Validation

Both controls are in the README with their numbers. The short version is no longer a
single sentence, and that is the result rather than a reporting problem:

- **Identity-permuted place fields** collapse the PRE/POST gap to approximately zero in two
  of the five sessions, leave +16 % to +42 % of it in both circular sessions, and leave
  **+62 %** of it in one linear session, Gatsby_08022013 (D8.8). Geometry is a strong
  predictor of a large residual but not a sufficient explanation, and whether the residual
  shares a cause with the circular rZ compression of D8.4 is an open question. Every
  effect size is quoted with its residual attached.
- **Real running** scores a median rZ of +1.8 to +13.2 across the five sessions, ordered by
  ensemble size (D8.9). The low end is the session with 16 place cells, where only half the
  windows clear z = 1.65. `results/sensitivity_floor.png` is this series, and it is the
  context every effect size needs: the scoring's power varies nearly tenfold across
  sessions the pipeline otherwise treats identically.
- **The two controls are probably one finding.** Plotted against place-cell count rather
  than against session, the *absolute* identity-shuffle residual is monotone with ensemble
  size across all five sessions on the column-cycle shuffle — 0.005, 0.049, 0.087, 0.121,
  0.369 — while the four largest ensembles order the same way on the signed fraction with
  geometry alternating between them (D8.10). Read as noise rather than as bias, that says
  the negative control gets weaker as the ensemble shrinks, exactly as the positive control
  does. `results/residual_vs_ensemble_all.png` plots it with no curve fitted: five
  sessions, with geometry and ensemble size moving together, can show an ordering and
  nothing more.

### Measured along the way

- **During events a 20 ms bin carries 5–6 spikes** from the place-cell ensemble, not the
  ~1 spike found during running, and under 1 % of bins are empty. The posteriors are
  sharper than anticipated — a 15 cm/s animal and a population burst are different
  regimes.
- Forward and reverse traversals split 45–53 % in every sleep epoch of every session, with
  no consistent difference in score.
- The conjunction event set gives a *larger* PRE/POST effect than the full set in **9 of
  the 10** session-by-shuffle comparisons, consistent with it being a stricter threshold on
  event strength rather than a different kind of filter (D6.23). The single exception is
  Cicero_09102014's column cycle, +0.149 against +0.209 — the same cell where that
  session's rank-biserial is +0.019 at p = 0.61.

---

## 6. Per-cell contribution — `hc11.pcc`, `hc11.subsets`

### What it does

The paper's Figures 2 to 4 rest on a per-cell contribution score (PCC): each neuron's mean
contribution to PRE and POST sequences, cells with significant scores split at the median
of their PRE-to-POST change (ΔPCC) into *rigid* and *plastic*. **The main text gives four
sentences and defers the procedure to supplementary materials this project does not have.**
One figure caption — Fig. 2D's "templates in which an increasing percentage of neuron's
place fields were shuffled" — fixes the only procedural detail available: the perturbation
is applied to the decoding template, not to the spike train. Everything else below is our
choice, and a replication resting on a figure caption should say so plainly.

The estimator, for cell *i* and epoch *E*:

    PCC(i, E) = mean_e rZ(e) − mean_s mean_e rZ(e | field i placed at offset s)

Seven choices, all provisional, all in D9.3:

1. **The field is moved by circular shift**, the same offset in both running directions. It
   keeps the cell's peak rate, width and shape and destroys only *where* the field is, so a
   contribution is a contribution of position rather than of excitability.
2. **Every non-zero offset is used**, not a random sample — 15 placements on a 16-bin
   linear track, 27–28 on a circular one. This makes the estimator deterministic and gives
   the null labels for free. It also makes the within-cell permutation test coarse: with 15
   alternatives the finest two-sided p is 1/16 = 0.0625, so **no cell on a linear track can
   reach p < 0.05** and significance is decided on a standardised score (`pcc.strong_z`,
   1.65) with the exact p reported beside it.
3. **The perturbed score reuses the event's own shuffle null** rather than drawing a fresh
   one, which is the estimator's single approximation; it is validated below.
4. **All events of the epoch enter the mean**, not only those the cell fired in. A silent
   cell is not zero-contribution: its field never enters the spike term but does enter the
   decoder's normalisation, so moving it still moves the posterior. The present-only mean
   is kept as a diagnostic.
5. **Significance is signed**: a cell counts as strongly contributing when its real
   placement scores *better* than its own alternatives, not merely differently.
6. **One event mask across all placements** (`pcc.finite_everywhere`), so a contribution is
   never a difference between event sets.
7. **The null labels** take each shifted placement in turn as though it were the real
   field, with the remaining placements as its null — the same statistic, same events, same
   spike counts, same epoch asymmetry, everything except the cell owning its field.

`hc11.subsets` is the simpler sibling used for the paper's Fig. 4D: it splits an ensemble
into equal halves by off-line firing rate, and — the part that matters — into equal halves
at random, which is the control for splitting at all.

### Parameters

| | Value | Source |
|---|---|---|
| placements per cell | every non-zero circular offset | [PROVISIONAL] |
| strongly contributing | signed z ≥ 1.65 in PRE or POST | [PROVISIONAL], matches `sequences.significance_z` |
| rigid/plastic boundary | median ΔPCC within session | [GB2016] |
| null-reuse validation sample | 40,000 (cell, placement, event) triples | [PROVISIONAL] |
| vacuous-test threshold | null-label \|rank-biserial\| ≥ 0.6 | [PROVISIONAL] |
| random ensemble splits (Fig. 4D) | 20 | [PROVISIONAL] |

### Validation

- **The real template is scored through the same code path as every shifted one.** The
  first draft took the unperturbed mean from the step-5 results table and the perturbed
  means from its own loop; any difference between those paths would have landed in every
  cell's PCC in the same direction. Step 5's numbers now serve as an independent check
  instead: agreement is **1.6–1.8 × 10⁻⁵ rZ**, the rounding floor of the step-5 CSV, and a
  run aborts above 10⁻³.
- **The null-reuse approximation**, on 40,000 re-scored triples per session: correlation
  **0.996**, overall bias **+0.0004 ± 0.0006 rZ**, and — the quantity that does not cancel
  in a PRE-to-POST difference — epoch-differential bias bounded at **|bias| < 0.0023 rZ**.
  A *per-cell* bias below ~0.010 rZ is beyond what that sample can see, which matters for
  one session and is recorded in D9.4 and in the limitations below.
- **The silent-cell case is tested, not assumed**: two tests assert that shifting a silent
  cell's field changes the score and that its contribution is not structurally zero.
- **Two calibration floors**, both measured rather than asserted. A threefold firing-rate
  split of 38 cells against 38 produced a population gap indistinguishable from splitting
  the same cells at random, so a subset claim needs **~0.25–0.30 rZ** to be evidence
  (D9.2). A median split of *noise* at 15 cells against 15 produces rank-biserial effects
  up to **0.43–0.56**, and at 4 against 4 up to **1.00** (D9.5, D9.6).

### Outcome

The estimator is additive: per-cell contributions sum to the population PRE/POST difference
within a factor of two on every session and within 15 % on two of them. The split built on
it is not supported: of 32 rigid-versus-plastic comparisons **one** exceeds its own null's
95th percentile, and it ties rather than beats, in the session and shuffle with no
population effect, in the direction opposite to the paper's. Three of five sessions cannot
test the split at all. The PRE-to-POST contributor overlap — the paper's 73 % — runs from
29 % to 80 % across our sessions and **rises as the ensemble shrinks**, which is what a
threshold on a per-cell score does when each cell carries a larger share of the posterior.

---

## Controls and validation, in one place

| Control | Question | Result | Where |
|---|---|---|---|
| Leave-one-lap-out decoding | Does the decoder work without seeing the lap? | Median error 4.7–10.6 cm against 36–71 cm chance | D5.9 |
| Leakage checks on the decoder | Is the low error leakage? | Honest 5.18 cm against leaky 4.56 cm; occupancy and spikes partition exactly; three tests | D5.9 |
| Tuning-curve permutation | What is chance for this ensemble? | 36–71 cm | D5.8 |
| Flat-Poisson threshold check | Does the 3 SD rule produce events on noise? | Yes, 0.35–0.68/s — so event *count* alone proves nothing | D6.4 |
| **Spike-jitter null** | Is the population synchrony real? | 84–90 % of events exceed the 95th percentile of their own null; participation 1.7× the null | D6.5 |
| **Multi-shank support vs MUA agreement** | Do two independent checks pick out the same weak events? | Agreement rises monotonically, 37 % at zero support to 87 % at eleven | D6.19 |
| Cell-count-matched shank null | Are unmatched MUA events on shanks the detector misses? | No — concentration excess over the null is identical for matched and unmatched | D6.23 |
| Smoothing-kernel sweep | Is the 5 ms MUA-to-ripple lead an artifact of smoothing? | No — flat from 0 to 8 ms, less than one sample | D6.22 |
| **Identity-shuffled place fields** | Is the PRE/POST difference replay? | Collapses in 2 of 5 sessions; leaves +16–42 % in both circular sessions and +62 % in one linear one | D7.7, D8.8 |
| **MAZE running positive control** | Can the scoring recognise a known trajectory? | Median rZ +1.8 to +13.2 (46–95 % above 1.65), ordered by ensemble size | D7.7, D8.9 |
| **Per-session reporting with no pooling across geometry** | Is a pooled rZ a mixture of two scales? | Yes — a ring maximises over 28–29 offsets, a line over one | D8.4 |
| Residual against ensemble size | Does the control's residual follow cell count rather than geometry? | The absolute residual is monotone with ensemble size across all five sessions on the column-cycle shuffle; suggestive only, no fit | D8.10 |
| **Random ensemble splits** | Is a difference between two halves of an ensemble more than splitting it at all produces? | No — 9 of 10 observed gaps sit inside the random-split spread | D9.2 |
| **Null labels for the rigid/plastic split** | What does a median split of ΔPCC produce from fields carrying no position? | Effects as large as the real ones, and larger in three sessions | D9.5, D9.6 |
| **Additivity of the per-cell contribution** | Does the estimator decompose the population effect it claims to? | Yes — the per-cell sum matches within a factor of two on every session | D9.6 |
| Null-reuse validation | Is reusing each event's shuffle null safe for a perturbed template? | Bias below 0.0023 rZ on the quantity that matters; a per-cell bias below ~0.010 rZ is unresolved | D9.4 |

The five in bold are the ones the headline result depends on.

---

## Eight silent failures

Every one of these produced plausible output instead of an error. None raised an
exception and none made a number look obviously wrong. **Two different things caught the
eight, and the differences matter more than the count.** Seven were caught by a
**cross-check against something the result had to agree with** — a second measurement, a
physical expectation, a control, a document — which is the argument for the controls in
the previous section: they are not ceremony, they are this project's main detector. The
seventh entry was caught by **reading a citation against its paper**, because nothing about
it was wrong except who got the credit for choosing it — and a parameter that looks like
the authors' is a parameter nobody sweeps. A third detection mode, **reading the code**,
appears only among the near-misses below, where the error would have moved every number in
the same direction and so no cross-check could have seen it.

| # | What broke | What it produced instead of an error | What caught it |
|---|---|---|---|
| 1 | Band-passing across the probe's spatial profile, not across time (D6.15) | A flat channel-power profile: a uniform 44–46 on every one of 134 channels, reading as "all channels are equivalent" | Physical implausibility — real electrodes at different depths cannot all carry the same ripple power |
| 2 | NaN multi-shank support read as a clean zero (D6.18) | "No other shank saw this event", on all 18 MAZE events | The result being *too* clean — exactly zero everywhere, where real data is never that tidy |
| 3 | A duration confound in the spike-jitter null (D6.5) | The null looking *more* synchronous than the data, 22 active cells against 18 | A second comparison that controlled for duration, run because the first answer was the wrong sign |
| 4 | A cell-count confound in the shank-concentration test (D6.23) | Apparent confirmation: unmatched events did have more concentrated activity | A per-event null drawing the same number of cells, built before accepting a result that agreed with the hypothesis |
| 5 | p-values below 1 × 10⁻⁵ rounding to zero in saved tables (step 8c) | `p = 0.0` in two results CSVs while the printed report showed 3.4 × 10⁻⁶ | A verification script comparing the README's numbers against the CSVs, which flagged six "mismatches" that were the file's fault |
| 6 | An undefined statistic compared against a null (step 9b) | `p = 0.048`, the smallest value a 20-draw null can express, for a gap that was NaN | Reading the printed line beside the medians it came from, which were NaN too |
| 7 | A parameter credited to a paper that does not state it (step 10a-bis) | Nothing. Correct code, reasonable value, every test and control passing — and a number nobody thought to sweep | Reading the citation tag against the paper; five tags of twenty-five were wrong (D10.3) |
| 8 | A shuffle that preserved the thing it was meant to destroy (step 10b) | A null that beat the data: Claim A +0.47 and +0.52 against the real fit's +0.41, reading as "the model has learned nothing" | The sign — a null outscoring the data is the null's fault until shown otherwise |

**1. Band-passing across the probe instead of across time.** The first channel-power
analysis indexed `data[channel]` — a single time point across 134 channels — instead of
`data[:, channel]`. A 150–300 Hz filter runs happily on a 134-sample "signal" and returns
a few tens of units, slightly different per channel. Only the 1–4 Hz band complained, and
it complained with NaN rather than an error, because a narrow low-frequency Butterworth in
transfer-function form is unstable at that length. Several hours went into verifying the
recording before the fault was located in the analysis. *Prevented now by* explicit
`axis=0` on every filter call, a per-channel loop over one-dimensional arrays where the
axis cannot be ambiguous, `sosfiltfilt` instead of `filtfilt`, and three tests that pin
the behaviour by *answer* rather than by arithmetic — including one that plants a 190 Hz
burst on a known channel and asserts that channel is the one selected.

**2. "No shank saw this event", which actually meant "never measured".** Multi-shank
support needs a baseline per comparison channel, sampled from the same brain state. MAZE
immobility comes in fragments under a second, too short to sample from, so the baseline
was NaN, every z-score was NaN, and the count of supporting shanks came out as a clean
zero. *Prevented now by* a fallback that uses the state itself when no window fits, and by
a warning when nothing finite was measured. With the fallback, MAZE support is a median of
10 of 11.

**3. A duration confound in the jitter null.** Comparing the median participation of
*re-detected* null events against observed suggested the null was more synchronous.
It is an artifact of detection: jittering smooths the rate trace, so null events are
nearly twice as long (124 ms against 68 ms) and accumulate more distinct cells over their
longer span. Normalised for duration the observed events are denser (2.50 against 1.85
active cells per 10 ms). *Prevented now by* the rule that any comparison of participation
between two event sets must control for duration, which is why the decisive comparison
recounts participation inside the *observed* windows.

**4. A cell-count confound in the shank-concentration test.** Testing whether unmatched
MUA events have their active cells concentrated on particular shanks, the raw
concentration was higher for unmatched events — which looks like confirmation. But an
event with nine active cells cannot occupy twelve shanks, so fewer cells forces higher
concentration whatever the physiology, and unmatched events have fewer cells. *Prevented
now by* comparing against a per-event null that draws the same number of cells at random,
weighted by how many cells each shank actually has. Against that null the excess is
identical for both groups and the hypothesis is rejected.

**5. A results table that rounded its own p-values to zero.** Every results CSV was
written with `frame.round(5)`, so a p-value of 3.4 × 10⁻⁶ was saved as `0.0` — in
`replay_comparison.csv` and in 73 of the 90 rows of `unmatched_mua_hypotheses.csv`. The
printed report showed the real value, so nothing looked wrong while the analysis ran; only
the saved artifact was damaged, which is the harder version to notice because it surfaces
when someone reads the table months later. *Prevented now by* `results.round_for_csv`,
which rounds p-value columns to significant figures rather than decimal places, used by
`write_session_table` so every table written through it is covered.

**6. A missing number reported as the strongest possible result.** The random-split control
of step 9b compares an observed gap against 20 draws from a null. In one session the gap
was undefined — one half of the ensemble produced a degenerate posterior on one event, so
its z-score was NaN and the median inherited it — and `NaN >= draw` is false for every
draw, so the count of more-extreme draws was zero and the empirical p came back as 1/21 =
0.048. A missing measurement had turned into the most significant result the null could
produce. *Prevented now by* `subsets.empirical_p`, which returns NaN for a non-finite
observation and discards non-finite draws, and by `subsets.finite_in_both`, which drops an
event from *both* halves when either half cannot score it — so the two halves are always
compared on one event set. Both are unit-tested on exactly these cases.

**A seventh entry belongs here, and it is a different kind of failure from the six
above.** `config/params.yaml` attributed the multi-unit detector's boundary threshold to
Chen et al., who state only the 3 SD peak and say nothing about where a candidate's
boundaries extend to. The number was ours. An audit of all twenty-five citation tags
(D10.3) found four more: a duration limit credited to the wrong paper, a ripple threshold
carried over from a sentence about multi-unit activity, a sleep-state label credited to a
paper that does not describe the scoring, and — the opposite error — a 1 Hz peak-rate
criterion that **is** Chen et al.'s and that we had been calling our own choice.

**8. A shuffle that preserved the thing it was meant to destroy.** The fit-on-shuffled-data
control of step 10b (D10.6) rolled each cell's spike counts circularly *within* each running
bout, meaning to keep every cell's rate and burstiness while destroying the alignment
between cells. It destroyed nothing. Running bouts are a median of five to nine 250 ms bins
and the animal is in essentially one place for all of them, so a cell's rate within a bout
identifies the position whatever order the bins arrive in. The control scored Claim A at
+0.47 and +0.52 against the real fit's +0.41 — a null outscoring the data, which taken at
face value says the model has learned nothing and the phase should stop. *Prevented now by*
rolling along the whole concatenated running timeline instead, which is this project's
existing place-cell shuffle idiom from step 4; the null then falls to 0.05–0.16 on that
session. Caught by the sign, not by a test: a null is not supposed to beat the data, and
when one does the null is the first suspect.

**This is the first entry in the list whose failure pointed toward a negative, and that is
worth saying.** Every other one that had a direction — the cell-count confound, the rounded
p-values, the NaN comparison — made a hypothesis look *better* supported than it was, and
the habit of mind those three train is to be most suspicious of a result you like. This one
runs the other way: believed as written, it would have discarded a real effect because its
control was broken in the data's favour. A control that kills a result deserves exactly the
same scrutiny as a result that confirms one, and nothing in the first seven entries would
have prompted it.

**No test, control or cross-check could have found any of them.** The code was correct,
every number was reasonable, nothing disagreed with anything. The output of a pipeline is
identical whether a parameter is the authors' or ours; what differs is whether anyone ever
questions it. A parameter that looks settled does not get swept, and the boundary threshold
went un-swept for that reason until Phase 4 needed longer events and went looking. The
step-8a audit could not have caught it either, because it asked which parameters had been
*checked on more than two sessions*, not which attributions were true — a coverage gap in
the audit rather than a fault in the output.

The detection mode is therefore unique in this list: not a second measurement, not a
control, not reading the code, but **reading the tag against the paper it cites**. That is
now done for every tag, and the result is D10.3.

**Two near-misses belong in the same section, because they were caught differently.**


The first was caught by **reading code rather than output**, which is the detection mode
none of the eight above had. In the first draft of the per-cell contribution estimator the
real template was scored through one code path and the shifted templates through another:
the unperturbed mean came from the step-5 results table, the perturbed means from the
estimator's own loop. Nothing it produced would have looked wrong — the difference between
two correct-looking paths is a small constant, and it would have landed in *every* cell's
contribution in the *same* direction, which is precisely the shape no downstream control
can see. It was found in review, before the estimator was ever run on more than four cells,
and the fix was to compute the real placement through the same function as the other
fifteen and keep the step-5 table as an independent check (agreement: 1.6 × 10⁻⁵ rZ). The
lesson is narrow and worth stating: the six failures above were caught because something
independent disagreed, and this one could not have been, because the error would have moved
every number together.

The second: **a one-lap directional template** on the circular session (D7.3) depressed that session's sequence scores by as
much as the effect being measured, because the score maximises over directions and a noise
template hands the shuffles the same freedom. It was caught by checking lap counts per
direction before trusting the output — again a cross-check, not an error.

What they have in common is worth stating plainly, because it generalises past this
project: **a numerical pipeline fails by returning a number** — or, in the seventh case,
by returning the right number for a reason nobody checked. Three of the eight — the
cell-count confound, the rounded p-values and the NaN comparison — produced results that
*supported* the hypothesis under test or made it look stronger, which is the most dangerous
shape a silent failure can take, and the only reason any of them was found is that
something independent had to agree and did not. (Of the first five, two had that shape; the
NaN comparison makes three. The count is kept honest here because it is the kind of number
that drifts upward in retelling.) **Exactly one of the eight, the newest, pointed the other
way** — toward abandoning something real — which is why the catalogue is kept by detection
mode rather than by whether the error flattered the hypothesis.

---

## Known limitations and open questions

### The place-cell shortfall: 304 against 491

**What would settle it:** the Science supplementary methods, which would state the actual
inclusion criteria; or the authors' own linearisation over all MAZE time rather than the
released `OneDLocation`.

Applying only Chen et al.'s stated criterion — peak rate above 1 Hz, nothing else — still
gives 392, so this is not a threshold we have mis-set (D5.4). Widening the measurement base
to our own linearisation of the raw position raises the linear-track yield from 47 % to
59 % of pyramidal cells, so the restricted released mask is a real contributor, but 59 %
is still short of 87 % (D5.5). The remainder is unexplained. Our population is a
conservative subset of theirs, and every downstream result is built on that subset.

### Circular-track extent does not match the documented platform

**What would settle it:** the authors' linearisation procedure, or the raw tracking
calibration.

`OneDLocation` spans 2.84–2.90 m on the circular sessions, while a documented 1 m-diameter
platform implies a 3.14 m circumference and the median two-dimensional radius implies
2.82–3.34 m (D3.3). The cause is unknown: the linearisation may cover only the arc the
animal ran, may be scaled to a nominal rather than measured radius, or may exclude the
reward zone. The code takes the extent from the observed range per session rather than
from an assumed circumference, because an assumed 3.14 m would misplace every bin edge and
put the wrap discontinuity where the animal never was.

### Roughly 45 % of MUA events have no ripple partner

**Settled, and recorded here because the answer constrains what the LFP adds.** Three
explanations were tested (D6.23). They are *smaller*, not elsewhere: unmatched events carry
about 60 % of the active cells of matched ones (rank-biserial −0.58 to −0.63), are shorter
and weaker. The shank-coverage explanation is rejected against the cell-count-matched null
above. Detection noise contributes a little — they fail the jitter null 5–9 points more
often — but 80–87 % of them still beat it. They are the weak tail of one distribution.

The consequence: the LFP is not covering events the spikes miss, so there is no coverage
argument for an LFP-primary event set, and the conjunction is best understood as a
stricter threshold.

### How much sleep there is to search varies by session, in both directions

**What would settle it:** more sessions, or sleep-pressure measurements the dataset does
not contain.

Searched non-REM minutes, PRE → POST: 168.5 → 89.7, 104.0 → 73.0, 112.9 → 106.5, 100.7 →
63.8 — and **92.8 → 138.9 in Cicero_09102014** (D6.6). The "POST holds about half the sleep
of PRE" statement from the two-session stage was an Achilles fact; one session has half
again as much non-REM in POST. Where POST is shorter, the likely reason — the animal has
already slept some five hours before POST begins — cannot be tested here. Every rate is
normalised by searched time, and the PRE/POST test is reported both on all events and on
PRE subsampled to POST's count wherever PRE is the larger set (which it is not in either
Cicero session).

### Only 6–23 % of candidate events are long enough to score

**What would settle it:** a longer candidate event, which means revisiting the MUA
detector's boundary threshold — a step-4 change, not a step-5 one.

The five-bin floor at 20 ms needs a 100 ms event and ours run 55–69 ms (D7.5). The
scorable fraction is 14–23 % in four sessions and **6–9 % in Cicero_09012014**, whose
16-cell ensemble rarely puts spikes from five cells into five consecutive bins; that
session still contributes 310 PRE and 435 POST scorable events, and no comparison in this
repository rests on fewer than 100 per epoch. This is a
larger loss than every other filter combined. It is recorded rather than tuned: a
three-bin floor would mean fitting a correlation to three points. The same filter applies
to PRE and POST so the comparison is unaffected, but any absolute rate is a rate within a
selected subset, and the surviving events are systematically the longer ones with more
cells.

### Quiet wakefulness on the maze is undefined

**What would settle it:** the theta/delta ratio or the EMG channel of the LFP, both of
which are now reachable — the EMG channels are identified as 128–130 and the band
machinery exists.

Chen et al.'s speed-below-2 cm/s criterion yields 0.4 minutes of searchable time out of
34–45 (D6.3). The sub-threshold samples are turnaround zero-crossings with a median run of
0.18 s, not pauses. The genuinely quiet periods are plausibly the 7.6–8.7 minutes of lost
tracking, when the animal sits at a reward site with its head-mounted light occluded. The
elevated population rate during those periods is *ambiguous*, not disqualifying: awake
sharp-wave ripples at a reward site would produce exactly that signature. The few MAZE
events that are detected look like real ripples — median support 10–11 of 11 other shanks
— which supports the hypothesis rather than the alternative (D6.21).

### The absolute level of rZ carries a bias

**What would settle it:** a shuffle that preserves both each bin's posterior and the
population's positional localisation. Time-bin permutation is closer than column cycling
but is not zero either.

Identity-permuted fields still give a median rZ of +0.59 under column cycling (D7.6).
Threshold-crossing fractions are therefore not reported as results anywhere in this
repository, and only the PRE/POST difference — which the control collapses — is treated as
information.

### Five sessions of eight, three animals of four

Event detection and sequence scoring need the LFP, which is roughly 9.5 GB compressed per
session. Five were downloaded: Achilles_10252013, Achilles_11012013, Cicero_09012014,
Cicero_09102014 and Gatsby_08022013. **Buddy_06272013, Cicero_09172014 and Gatsby_08282013
are omitted because their LFP was not downloaded, for local storage reasons** — not because
of anything in the data and not because they failed any check (D8.9). Steps 1–3 above run
on all eight sessions, and the place-cell total of 304 is all eight.

Each session is reported separately as the primary result. Pooling is secondary and happens
only within a track geometry, because rZ is not on the same scale on a ring as on a line
(D8.4). The omissions have one methodological consequence beyond sample size: both sessions
whose XML carries a `skip="1"` channel flag are among them, so that path is untested
against real data (D8.7).

### The sensitivity floor, and what it does to effect sizes

**What would settle it:** sessions that vary ensemble size independently of track geometry
and event count — which these five do not.

The positive control runs from +13.2 at 77 place cells to +1.8 at 16 (D8.9). Whatever the
replay effect sizes mean, they are not measurements on a common instrument, and the
smallest-ensemble session is also the one with the lowest scorable fraction. Ensemble size,
geometry and event count all move together here, so `results/sensitivity_floor.png` is a
floor — a statement of where the method stops working — and not a calibration curve. No rZ
in this repository is corrected by it.
