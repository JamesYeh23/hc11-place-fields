"""Step 12b: is the extended base's cell gain real in quality, or only in count?

    uv run python scripts/extended_base_checks.py

Widening the measurement base from the authors' `OneDLocation` to our own
linearisation raises the place-cell yield on the three linear LFP sessions:
77 -> 83, 16 -> 29, 30 -> 35 (D11.1). This asks whether those extra cells are
worth having, **before** anything downstream is allowed to use them.

The reason to ask is in the scope measurement. The authors' mask is fast
traversals of the whole track; what the extension adds is slow time at the ends
-- median speed 21-25 cm/s against 28-61, and 60-98 % of the extra running
samples in the outer 15 % of the track, which on a 1.6 m linear maze is where
the reward sites are. A cell that fires during slow approach and turnaround at a
reward site produces a map peak at a track end that can pass a spatial-
information test without being useful to a decoder.

**No cross-check currently in this project would catch that.** Phase 2's
sensitivity floor is indexed by cell *count*, and a bad cell counts the same as
a good one, so a gain that is real in count and spurious in quality would make
every downstream number look better while the instrument got worse. That is why
check 2 exists.

Three checks, and the second decides:

1. **Where the gained cells' fields peak**, against the cells that were already
   there. The retained-cell baseline is 46 / 55 / 46 % of peaks in the outer
   15 % of the track against 30 % expected under uniform coverage, so the
   comparison is gained-against-retained, not against an absolute threshold.
2. **Decoding error by true position, both bases.** A real gain lowers error
   across the track; an end-field gain lowers it in the outer thirds and leaves
   the middle flat. The middle is where replay sequences have to resolve.
3. **The gained cells' in-field speed.** A field whose spikes arrive almost
   entirely below about 25 cm/s is a reward-site response rather than a place
   field.

The two bases are reported side by side throughout and the authors' mask stays
primary. Nothing here changes a Phase 1-4 number.
"""

from __future__ import annotations

import json
import logging
from datetime import datetime, timezone

import numpy as np
import pandas as pd

from hc11 import config, plots, results
from hc11.behavior import speed_cm_s
from hc11.decode import decode, prepare_decoder
from hc11.fields import select_running_samples
from place_fields import prepare  # noqa: I001  (sibling script)

log = logging.getLogger("extended")

LINEAR_LFP = ("Achilles_10252013", "Cicero_09012014", "Gatsby_08022013")
#: Fraction of the track at each end counted as "the ends". 0.15 a side, so a
#: uniform population would put 30 % of its peaks there. [PROVISIONAL] -- chosen
#: in the 12a scope before any of these numbers existed, and not revisited.
END_FRACTION = 0.15
#: Thirds of the track for check 2: a real gain improves all three.
N_THIRDS = 3
#: A spike counts as in-field if it falls within half a field width of the peak.
#: [PROVISIONAL] -- the summary table carries a width and a peak and no explicit
#: field boundary, so this is the simplest reading of the two numbers it has.
FIELD_HALF_WIDTHS = 0.5


def affine_to_authors(session_authors, session_extended):
    """``(a, b)`` with ``theirs ~= a * ours + b``, fitted where both are defined.

    On a linear track the two linearisations are the same projection up to
    origin and sign (|r| = 1.0000), so this is exact rather than a model, and it
    lets the extended base's errors be read on the authors' 1.6 m scale instead
    of on its own longer one.
    """
    theirs, ours = session_authors.position_1d, session_extended.position_1d
    both = np.isfinite(theirs) & np.isfinite(ours)
    a, b = np.polyfit(ours[both], theirs[both], 1)
    residual = float(np.max(np.abs(theirs[both] - (a * ours[both] + b))))
    return float(a), float(b), residual, int(both.sum())


def decoded_pairs(session, track, running, laps, cluster_ids, params, time_bin_s):
    """``(true, decoded)`` in metres for the direction-merged decoder.

    ``direction=None`` is the merged condition, as in `decode_position.py`:
    one template built across both running directions, which is the condition
    Phase 2's ensemble curves use.
    """
    samples = select_running_samples(session, track, laps, running.mask, None)
    inputs = prepare_decoder(session, track, samples, laps, cluster_ids, params,
                             time_bin_s=time_bin_s)
    result = decode(inputs, track, params)
    # Unfiltered, and in bin order: two ensembles on the same base produce the
    # same bins in the same order, so the caller can intersect their decodable
    # sets and compare them bin for bin. Filtering here would hide that a
    # smaller ensemble leaves more bins undecodable and silently compare two
    # different stretches of track.
    return result.true_pos, result.decoded_pos


def in_field_speed(session, running, cells, track):
    """Median speed at the in-field spikes of each cell, in cm/s."""
    speed = speed_cm_s(session, 0.0)
    position = session.position_1d
    time = session.position_t
    out = {}
    for cluster_id, peak, width in cells:
        if not np.isfinite(peak) or not np.isfinite(width) or width <= 0:
            continue
        spikes = session.spikes_for(int(cluster_id))
        if spikes.size == 0:
            continue
        index = np.clip(np.searchsorted(time, spikes), 0, len(time) - 1)
        inside = running.mask[index] & (
            np.abs(position[index] - peak) <= FIELD_HALF_WIDTHS * width)
        if inside.sum() < 10:
            continue
        out[int(cluster_id)] = float(np.nanmedian(speed[index[inside]]))
    return out


def main() -> None:
    logging.basicConfig(level=logging.INFO, format="%(message)s")
    params = config.load_params()
    time_bin_s = params["decoding"]["time_bin_s"]
    authors_fields = pd.read_csv(config.results_dir() / "place_fields_summary.csv")
    extended_fields = pd.read_csv(
        config.results_dir() / "place_fields_summary_extended.csv")

    peak_rows, error_rows, speed_rows, summary = [], [], [], []
    for name in LINEAR_LFP:
        a_session, a_track, a_running, a_laps = prepare(name, params)
        e_session, e_track, e_running, e_laps = prepare(name, params, base="extended")
        slope, intercept, residual, n_both = affine_to_authors(a_session, e_session)
        log.info("%s: extended -> authors' coordinate is %.4f x + %.4f "
                 "(max residual %.4f m over %d samples)",
                 name, slope, intercept, residual, n_both)

        authors_cells = set(authors_fields[(authors_fields.session == name)
                                           & authors_fields.is_place_cell].cluster_id)
        extended_rows = extended_fields[(extended_fields.session == name)
                                        & extended_fields.is_place_cell]
        extended_cells = set(extended_rows.cluster_id)
        gained = extended_cells - authors_cells

        # --- check 1: where the fields peak, on the authors' scale ------------
        # Three groups, not two. "retained, authors' base" is the same cells
        # scored the old way: without it, a shift in the retained cells' own
        # peaks would be invisible, and it is not -- the extension rebuilds
        # their fields over end-inclusive time too.
        low = float(np.nanmin(a_session.position_1d))

        authors_rows = authors_fields[(authors_fields.session == name)
                                      & authors_fields.is_place_cell]
        wanted_peaks = [
            (row, "gained" if row.cluster_id in gained else "retained",
             slope * row.peak_position_m + intercept)
            for row in extended_rows.itertuples()
        ] + [
            (row, "retained, authors' base", row.peak_position_m)
            for row in authors_rows.itertuples() if row.cluster_id in extended_cells
        ]
        for row, group, peak in wanted_peaks:
            relative = (peak - low) / a_track.extent
            peak_rows.append({
                "session": name, "cluster_id": int(row.cluster_id),
                "direction": row.direction, "group": group,
                "peak_on_authors_scale_m": peak, "relative_position": relative,
                "in_outer_ends": bool(relative < END_FRACTION
                                      or relative > 1 - END_FRACTION)})

        # --- check 2: decoding error by position, both bases ------------------
        # A third row isolates the cells from the base: the extended base's
        # templates with only the cells that were already there. If it matches
        # "extended", the gained cells contribute nothing of their own.
        raw = {}
        for base, session, track, running, laps, cells in (
                ("authors", a_session, a_track, a_running, a_laps,
                 np.array(sorted(authors_cells))),
                ("extended", e_session, e_track, e_running, e_laps,
                 np.array(sorted(extended_cells))),
                ("extended, retained only", e_session, e_track, e_running, e_laps,
                 np.array(sorted(extended_cells & authors_cells)))):
            raw[base] = (cells, *decoded_pairs(session, track, running, laps, cells,
                                               params, time_bin_s))
        # The two extended rows share a base, so they share their bins: compare
        # them only where *both* ensembles decoded, or the smaller one is scored
        # on an easier subset of the track.
        shared = np.isfinite(raw["extended"][2]) & np.isfinite(
            raw["extended, retained only"][2])
        for base, (cells, true_pos, decoded) in raw.items():
            usable = (np.isfinite(decoded) if base == "authors" else shared)
            true_pos, decoded = true_pos[usable], decoded[usable]
            if base.startswith("extended"):
                true_pos = slope * true_pos + intercept
                decoded = slope * decoded + intercept
            relative = (true_pos - low) / a_track.extent
            # Only the stretch of track the authors' base also covers. The
            # extended coordinate runs about 2.0 m where theirs runs 1.6, so
            # clipping the overhang into the outer thirds would charge the
            # extended base for decoding track the comparison does not include
            # -- which on the first run turned Achilles_10252013's outer-low
            # error into 25 cm and would have been read as the gain failing.
            inside = (relative >= 0) & (relative <= 1)
            dropped = int((~inside).sum())
            relative, error = relative[inside], np.abs(
                decoded - true_pos)[inside] * 100
            third = np.minimum((relative * N_THIRDS).astype(int), N_THIRDS - 1)
            for index in range(N_THIRDS):
                here = third == index
                error_rows.append({
                    "session": name, "base": base, "n_cells": len(cells),
                    "third": ["outer low", "middle", "outer high"][index],
                    "n_bins": int(here.sum()), "n_bins_off_authors_track": dropped,
                    "median_error_cm": float(np.median(error[here]))
                    if here.any() else np.nan})

        # --- check 3: in-field speed of the gained cells ----------------------
        wanted = [(r.cluster_id, r.peak_position_m, r.largest_field_width_m)
                  for r in extended_rows.itertuples()]
        speeds = in_field_speed(e_session, e_running, wanted, e_track)
        for cluster_id, value in speeds.items():
            speed_rows.append({
                "session": name, "cluster_id": cluster_id,
                "group": "gained" if cluster_id in gained else "retained",
                "median_in_field_speed_cm_s": value})

        summary.append({"session": name, "n_authors": len(authors_cells),
                        "n_extended": len(extended_cells), "n_gained": len(gained),
                        "n_lost": len(authors_cells - extended_cells),
                        "affine_slope": slope, "affine_max_residual_m": residual})

    peaks = pd.DataFrame(peak_rows)
    errors = pd.DataFrame(error_rows)
    speeds = pd.DataFrame(speed_rows)
    out = config.results_dir()
    results.round_for_csv(peaks).to_csv(out / "extended_base_peaks.csv", index=False)
    results.round_for_csv(errors).to_csv(out / "extended_base_errors.csv", index=False)
    results.round_for_csv(speeds).to_csv(out / "extended_base_speeds.csv", index=False)
    results.round_for_csv(pd.DataFrame(summary)).to_csv(
        out / "extended_base_summary.csv", index=False)

    log.info("\nCHECK 1 -- where the fields peak (%% in the outer %.0f %% of the track; "
             "uniform would be %.0f %%)", 100 * END_FRACTION, 200 * END_FRACTION)
    for name in LINEAR_LFP:
        part = peaks[peaks.session == name]
        line = f"  {name:20s}"
        for group in ("retained, authors' base", "retained", "gained"):
            here = part[part.group == group]
            label = "retained (authors)" if "authors" in group else group
            line += (f" | {label} {100 * here.in_outer_ends.mean():4.0f}% "
                     f"(n={len(here):3d})" if len(here) else f" | {label}  --")
        log.info(line)

    log.info("\nCHECK 2 -- decoding error by third of the track, cm "
             "(authors' coordinate, %g s bins)", time_bin_s)
    for name in LINEAR_LFP:
        part = errors[errors.session == name]
        for base in ("authors", "extended", "extended, retained only"):
            here = part[part.base == base].set_index("third")
            cells = int(here.n_cells.iloc[0])
            log.info("  %-20s %-23s %2d cells | %s | %d off-track bins",
                     name, base, cells,
                     "  ".join(f"{t}: {here.median_error_cm[t]:5.1f}"
                               for t in ("outer low", "middle", "outer high")),
                     int(here.n_bins_off_authors_track.iloc[0]))

    log.info("\nCHECK 3 -- median in-field speed, cm/s")
    for name in LINEAR_LFP:
        part = speeds[speeds.session == name]
        line = f"  {name:20s}"
        for group in ("retained", "gained"):
            here = part[part.group == group]
            line += (f" {group} {here.median_in_field_speed_cm_s.median():5.1f} "
                     f"(n={len(here):3d})" if len(here) else f" {group}   --")
        log.info(line)

    plots.plot_extended_base(peaks, errors, speeds,
                             out / "extended_base_checks.png")
    (out / "extended_base_checks.json").write_text(json.dumps(
        {"script": "scripts/extended_base_checks.py", "git_commit": config.git_commit(),
         "created": datetime.now(timezone.utc).isoformat(timespec="seconds"),
         "sessions": list(LINEAR_LFP), "end_fraction": END_FRACTION,
         "n_thirds": N_THIRDS, "field_half_widths": FIELD_HALF_WIDTHS,
         "time_bin_s": time_bin_s,
         "retained_baseline_pct": "46 / 55 / 46 (D11.1)"}, indent=2) + "\n")
    log.info("\nwrote results/extended_base_*.csv")


if __name__ == "__main__":
    main()
