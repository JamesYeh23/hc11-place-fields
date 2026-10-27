"""Figures for the place-field analysis.

Colour choices, once, here:

* **Directions are categorical**, so they get two fixed hues from the Okabe-Ito
  colourblind-safe set -- blue ``#0072B2`` for increasing position, vermillion
  ``#D55E00`` for decreasing -- assigned in a fixed order and never cycled. The
  pair was checked for colour-vision-deficiency separation (worst-case protan
  Delta-E 21.9, normal-vision 31.2, both well above the 8/15 floors). Direction is
  also carried by the legend and by panel titles, never by colour alone.
* **Firing rate is magnitude**, so it gets a sequential, perceptually uniform
  ramp (``magma``, light to dark) rather than a rainbow map like ``jet``, whose
  uneven lightness invents boundaries that are not in the data.
* Axes and grids are recessive; text stays in neutral ink rather than taking the
  series colour.
"""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402

#: Fixed categorical hues for running direction (Okabe-Ito).
DIRECTION_COLORS = {"pos": "#0072B2", "neg": "#D55E00"}
DIRECTION_LABELS = {"pos": "increasing position", "neg": "decreasing position"}
#: Sequential ramp for firing rate (perceptually uniform, CVD-safe).
RATE_CMAP = "magma"

INK = "#1a1a1a"
MUTED = "#6b6b6b"
GRID = "#e6e6e6"


def _style(ax) -> None:
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(MUTED)
        ax.spines[side].set_linewidth(0.8)
    ax.tick_params(colors=MUTED, labelsize=8, width=0.8)
    for label in ax.get_xticklabels() + ax.get_yticklabels():
        label.set_color(INK)
    ax.grid(True, color=GRID, linewidth=0.6, alpha=0.9)
    ax.set_axisbelow(True)


def _save(fig, path: Path, params: dict, dpi: int = 150) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=dpi, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    path.with_suffix(".json").write_text(json.dumps(params, indent=2, default=str) + "\n")
    return path


def plot_laps(session, track, laps, running, path: Path, params: dict) -> Path:
    """Linearised position against time, with laps coloured by direction."""
    fig, ax = plt.subplots(figsize=(13, 4))
    t0 = session.epoch("MAZE")[0]
    t = session.position_t - t0
    x = session.position_1d

    ax.plot(t, x, color="#c9c9c9", linewidth=0.6, label="on track (authors' mask)", zorder=1)
    seen = set()
    for lap in laps:
        sl = slice(lap.start_sample, lap.stop_sample)
        ax.plot(
            t[sl], x[sl],
            color=DIRECTION_COLORS[lap.direction], linewidth=1.6, zorder=3,
            label=(None if lap.direction in seen else DIRECTION_LABELS[lap.direction]),
        )
        seen.add(lap.direction)
    slow = running.authors_mask & ~running.mask
    ax.plot(t[slow], x[slow], ".", color=MUTED, markersize=1.6, alpha=0.5, zorder=2,
            label=f"below {running.speed_threshold_cm_s:g} cm/s")

    n_pos = sum(lap.direction == "pos" for lap in laps)
    n_neg = sum(lap.direction == "neg" for lap in laps)
    ax.set_xlabel("time from MAZE onset (s)", color=INK, fontsize=9)
    ax.set_ylabel("linearised position (m)", color=INK, fontsize=9)
    ax.set_title(
        f"{session.name} — {session.maze_type} — {n_pos} laps increasing, {n_neg} decreasing",
        color=INK, fontsize=11, loc="left", pad=28,
    )
    ax.legend(loc="lower center", bbox_to_anchor=(0.5, 1.02), fontsize=8, frameon=False,
              ncols=4, labelcolor=INK)
    _style(ax)
    return _save(fig, path, params)


def plot_rate_map_heatmap(cells, track, session_name, direction, path: Path, params: dict) -> Path:
    """Place cells as rows, sorted by peak position, each normalised to its own peak."""
    cells = [c for c in cells if c.is_place_cell and np.isfinite(c.peak_position_m)]
    fig, ax = plt.subplots(figsize=(6.5, max(3.0, 0.11 * max(len(cells), 8))))
    if not cells:
        ax.text(0.5, 0.5, "no place cells", ha="center", va="center", color=MUTED,
                transform=ax.transAxes)
        _style(ax)
        return _save(fig, path, params)

    order = np.argsort([c.peak_position_m for c in cells])
    stack = np.vstack([cells[i].rate_map.rate for i in order])
    with np.errstate(invalid="ignore"):
        stack = stack / np.nanmax(stack, axis=1, keepdims=True)

    im = ax.imshow(
        stack, aspect="auto", cmap=RATE_CMAP, vmin=0, vmax=1, interpolation="nearest",
        extent=(track.lo, track.hi, len(cells) - 0.5, -0.5),
    )
    cbar = fig.colorbar(im, ax=ax, pad=0.02, fraction=0.04)
    cbar.set_label("firing rate / cell peak", color=INK, fontsize=8)
    cbar.ax.tick_params(colors=MUTED, labelsize=7)
    ax.set_xlabel("linearised position (m)", color=INK, fontsize=9)
    ax.set_ylabel("place cell (sorted by peak position)", color=INK, fontsize=9)
    ax.set_title(
        f"{session_name} — {DIRECTION_LABELS[direction]} — {len(cells)} place cells",
        color=INK, fontsize=11, loc="left",
    )
    ax.grid(False)
    ax.tick_params(colors=MUTED, labelsize=8)
    return _save(fig, path, params)


def plot_field_grid(cells, track, session_name, direction, path: Path, params: dict,
                    n_show: int = 20) -> Path:
    """Individual rate maps for the most spatially informative cells."""
    cells = [c for c in cells if c.is_place_cell]
    cells = sorted(cells, key=lambda c: -c.info_bits_per_spike)[:n_show]
    if not cells:
        fig, ax = plt.subplots(figsize=(4, 2))
        ax.text(0.5, 0.5, "no place cells", ha="center", va="center", color=MUTED)
        ax.axis("off")
        return _save(fig, path, params)

    ncols = 5
    nrows = int(np.ceil(len(cells) / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(2.5 * ncols, 1.7 * nrows), squeeze=False)
    color = DIRECTION_COLORS[direction]
    for ax, cell in zip(axes.ravel(), cells, strict=False):
        rate = cell.rate_map.rate
        ax.plot(track.bin_centers, rate, color=color, linewidth=1.8)
        ax.fill_between(track.bin_centers, 0, np.nan_to_num(rate), color=color, alpha=0.18)
        ax.set_title(
            f"unit {cell.cluster_id} · {cell.info_bits_per_spike:.2f} bits/spk\n"
            f"{cell.peak_rate:.1f} Hz · r={cell.stability:.2f}",
            fontsize=7, color=INK,
        )
        ax.set_xlim(track.lo, track.hi)
        ax.set_ylim(0, None)
        _style(ax)
        ax.tick_params(labelsize=6)
    for ax in axes.ravel()[len(cells):]:
        ax.axis("off")
    fig.suptitle(
        f"{session_name} — {DIRECTION_LABELS[direction]} — top {len(cells)} by spatial information",
        fontsize=11, color=INK, x=0.01, ha="left",
    )
    fig.supxlabel("linearised position (m)", fontsize=9, color=INK)
    fig.supylabel("firing rate (Hz)", fontsize=9, color=INK)
    fig.tight_layout()
    return _save(fig, path, params)


def plot_info_histogram(cells, session_name, path: Path, params: dict) -> Path:
    """Spatial information, place cells against the rest."""
    fig, ax = plt.subplots(figsize=(6, 3.4))
    info = np.array([c.info_bits_per_spike for c in cells], dtype=float)
    is_pc = np.array([c.is_place_cell for c in cells], dtype=bool)
    finite = np.isfinite(info)
    if finite.any():
        bins = np.histogram_bin_edges(info[finite], bins=30)
        ax.hist(info[finite & ~is_pc], bins=bins, color="#bdbdbd",
                label=f"not a place cell (n={int((finite & ~is_pc).sum())})")
        ax.hist(info[finite & is_pc], bins=bins, color=DIRECTION_COLORS["pos"], alpha=0.85,
                label=f"place cell (n={int((finite & is_pc).sum())})")
        median = float(np.nanmedian(info[finite & is_pc])) if (finite & is_pc).any() else np.nan
        if np.isfinite(median):
            ax.axvline(median, color=INK, linewidth=1.2, linestyle="--")
            ax.text(median, ax.get_ylim()[1] * 0.95, f"  median {median:.2f}",
                    color=INK, fontsize=8, va="top")
    ax.set_xlabel("Skaggs spatial information (bits/spike)", color=INK, fontsize=9)
    ax.set_ylabel("cell-directions", color=INK, fontsize=9)
    ax.set_title(f"{session_name} — spatial information", color=INK, fontsize=11, loc="left")
    ax.legend(fontsize=8, frameon=False, labelcolor=INK)
    _style(ax)
    return _save(fig, path, params)


# ---------------------------------------------------------------------------
# Decoding figures
# ---------------------------------------------------------------------------

#: Bin sizes are ordered, not categorical, so they get steps of one hue.
BIN_SIZE_COLORS = {0.25: "#0072B2", 0.02: "#56B4E9"}
CHANCE_COLOR = "#6b6b6b"


def plot_decoding_examples(result, inputs, track, path: Path, params: dict,
                           n_laps: int = 4) -> Path:
    """Posterior as a heatmap with the true trajectory drawn over it."""
    laps = [lap for lap in inputs.laps if lap.index in result.posteriors][:n_laps]
    if not laps:
        fig, ax = plt.subplots(figsize=(4, 2))
        ax.text(0.5, 0.5, "no decoded laps", ha="center", va="center", color=MUTED)
        ax.axis("off")
        return _save(fig, path, params)

    fig, axes = plt.subplots(1, len(laps), figsize=(3.4 * len(laps), 3.2), squeeze=False)
    for ax, lap in zip(axes[0], laps, strict=False):
        posterior = np.exp(result.posteriors[lap.index])
        t_rel = lap.centres - lap.centres[0]
        im = ax.imshow(
            posterior, aspect="auto", origin="lower", cmap=RATE_CMAP,
            extent=(t_rel[0], t_rel[-1] + result.time_bin_s, track.lo, track.hi),
            vmin=0, vmax=max(float(posterior.max()), 1e-9),
        )
        ax.plot(t_rel + result.time_bin_s / 2, lap.true_pos, color="#ffffff",
                linewidth=2.4, alpha=0.95)
        ax.plot(t_rel + result.time_bin_s / 2, lap.true_pos, color=INK,
                linewidth=1.2, label="true position")
        ax.set_title(f"lap {lap.index} ({lap.direction})", fontsize=8, color=INK)
        ax.set_xlabel("time in lap (s)", fontsize=8, color=INK)
        ax.tick_params(colors=MUTED, labelsize=7)
        ax.grid(False)
    axes[0][0].set_ylabel("position (m)", fontsize=9, color=INK)
    axes[0][0].legend(fontsize=7, frameon=False, loc="upper left", labelcolor=INK)
    cbar = fig.colorbar(im, ax=axes[0].tolist(), pad=0.015, fraction=0.03)
    cbar.set_label("posterior probability", color=INK, fontsize=8)
    cbar.ax.tick_params(colors=MUTED, labelsize=7)
    fig.suptitle(
        f"{result.session} — {DIRECTION_LABELS.get(result.direction, result.direction)} — "
        f"{result.n_cells} cells, {result.time_bin_s * 1000:.0f} ms bins, "
        f"median error {result.median_error_cm:.1f} cm",
        fontsize=10, color=INK, x=0.01, ha="left",
    )
    return _save(fig, path, params)


def plot_error_distribution(results, chance_by_bin, session_name, path: Path,
                            params: dict) -> Path:
    """Cumulative error distribution, one line per bin size, with the chance level."""
    fig, ax = plt.subplots(figsize=(6, 3.6))
    for result in results:
        errors = np.sort(result.error_cm[result.decodable])
        if errors.size == 0:
            continue
        color = BIN_SIZE_COLORS.get(result.time_bin_s, INK)
        ax.plot(errors, np.linspace(0, 1, errors.size), color=color, linewidth=2,
                label=f"{result.time_bin_s * 1000:.0f} ms bins "
                      f"(median {result.median_error_cm:.1f} cm)")
        ax.plot([result.median_error_cm], [0.5], "o", color=color, markersize=5)
    for chance in chance_by_bin.values():
        if np.isfinite(chance):
            ax.axvline(chance, color=CHANCE_COLOR, linestyle="--", linewidth=1.2)
            ax.text(chance, 0.04, f" chance {chance:.0f} cm", color=CHANCE_COLOR, fontsize=8)
            break
    ax.set_xlabel("absolute decoding error (cm)", color=INK, fontsize=9)
    ax.set_ylabel("cumulative fraction of time bins", color=INK, fontsize=9)
    ax.set_title(f"{session_name} — decoding error", color=INK, fontsize=11, loc="left")
    ax.set_ylim(0, 1)
    ax.legend(fontsize=8, frameon=False, loc="lower right", labelcolor=INK)
    _style(ax)
    return _save(fig, path, params)


def plot_error_vs_ensemble(curves, session_name, path: Path, params: dict) -> Path:
    """Median error against the number of cells, one line per bin size."""
    fig, ax = plt.subplots(figsize=(6, 3.6))
    for bin_s, rows in sorted(curves.items(), reverse=True):
        if not rows:
            continue
        n = [r["n_cells"] for r in rows]
        med = [r["median_error_cm"] for r in rows]
        color = BIN_SIZE_COLORS.get(bin_s, INK)
        ax.fill_between(n, [r["p25_cm"] for r in rows], [r["p75_cm"] for r in rows],
                        color=color, alpha=0.15)
        ax.plot(n, med, "-o", color=color, linewidth=2, markersize=5,
                label=f"{bin_s * 1000:.0f} ms bins")
        ax.annotate(f"{med[-1]:.1f} cm", (n[-1], med[-1]), textcoords="offset points",
                    xytext=(6, 0), fontsize=8, color=INK, va="center")
    ax.set_xlabel("number of place cells in the ensemble", color=INK, fontsize=9)
    ax.set_ylabel("median decoding error (cm)", color=INK, fontsize=9)
    ax.set_title(f"{session_name} — error vs ensemble size", color=INK, fontsize=11, loc="left")
    ax.set_ylim(0, None)
    ax.legend(fontsize=8, frameon=False, labelcolor=INK)
    _style(ax)
    return _save(fig, path, params)


def plot_confusion(result, track, path: Path, params: dict) -> Path:
    """True position (rows) against decoded position (columns), row-normalised."""
    matrix = result.confusion_matrix(track)
    fig, ax = plt.subplots(figsize=(4.6, 4.0))
    im = ax.imshow(matrix, origin="lower", cmap=RATE_CMAP, vmin=0, vmax=1,
                   extent=(track.lo, track.hi, track.lo, track.hi), interpolation="nearest")
    ax.plot([track.lo, track.hi], [track.lo, track.hi], color="#ffffff", linewidth=0.8,
            alpha=0.6)
    cbar = fig.colorbar(im, ax=ax, pad=0.02, fraction=0.045)
    cbar.set_label("fraction of bins", color=INK, fontsize=8)
    cbar.ax.tick_params(colors=MUTED, labelsize=7)
    ax.set_xlabel("decoded position (m)", color=INK, fontsize=9)
    ax.set_ylabel("true position (m)", color=INK, fontsize=9)
    ax.set_title(
        f"{result.session} — {result.time_bin_s * 1000:.0f} ms bins, "
        f"median {result.median_error_cm:.1f} cm",
        color=INK, fontsize=10, loc="left",
    )
    ax.grid(False)
    ax.tick_params(colors=MUTED, labelsize=8)
    return _save(fig, path, params)


# ---------------------------------------------------------------------------
# Candidate-event figures (step 6)
# ---------------------------------------------------------------------------

EPOCH_COLORS = {"PRE": "#0072B2", "MAZE": "#009E73", "POST": "#D55E00"}


def plot_event_examples(session, events, part, cluster_ids, epoch, path: Path,
                        params: dict, n_show: int = 6) -> Path:
    """Spike rasters around events spanning the participation range."""
    if len(events) == 0:
        fig, ax = plt.subplots(figsize=(4, 2))
        ax.text(0.5, 0.5, "no events", ha="center", va="center", color=MUTED)
        ax.axis("off")
        return _save(fig, path, params)

    order = np.argsort(part.n_active)
    picks = order[np.linspace(0, len(order) - 1, min(n_show, len(order))).astype(int)]
    cluster_ids = np.asarray(cluster_ids)
    row_of = {int(c): i for i, c in enumerate(cluster_ids)}

    fig, axes = plt.subplots(1, len(picks), figsize=(2.6 * len(picks), 3.6), squeeze=False)
    for ax, idx in zip(axes[0], picks, strict=True):
        centre = events.peak[idx]
        lo, hi = centre - 0.25, centre + 0.25
        times, ids = session.spikes_in((lo, hi))
        keep = np.isin(ids, cluster_ids)
        rows = np.array([row_of[int(i)] for i in ids[keep]]) if keep.any() else np.empty(0)
        ax.plot((times[keep] - centre) * 1000, rows, "|", color=INK, markersize=3,
                markeredgewidth=0.8)
        ax.axvspan((events.start[idx] - centre) * 1000, (events.end[idx] - centre) * 1000,
                   color=EPOCH_COLORS.get(epoch, "#0072B2"), alpha=0.16, zorder=0)
        ax.set_title(f"{part.n_active[idx]}/{part.n_cells} cells\n"
                     f"{1000 * (events.end[idx] - events.start[idx]):.0f} ms, "
                     f"z={events.peak_z[idx]:.1f}", fontsize=7, color=INK)
        ax.set_xlim(-250, 250)
        ax.set_ylim(-1, len(cluster_ids))
        ax.set_xlabel("ms from peak", fontsize=8, color=INK)
        _style(ax)
        ax.tick_params(labelsize=6)
    axes[0][0].set_ylabel("pyramidal cell", fontsize=9, color=INK)
    fig.suptitle(f"{session.name} — {epoch} — events across the participation range",
                 fontsize=10, color=INK, x=0.01, ha="left")
    fig.tight_layout()
    return _save(fig, path, params)


def plot_event_rate_over_session(session, events_by_epoch, path: Path, params: dict,
                                 window_min: float = 5.0) -> Path:
    """Event rate through the session, one colour per epoch."""
    fig, ax = plt.subplots(figsize=(11, 3.2))
    for epoch, events in events_by_epoch.items():
        if len(events) == 0:
            continue
        start, end = session.epoch(epoch)
        edges = np.arange(start, end + window_min * 60, window_min * 60)
        counts = np.histogram(events.peak, bins=edges)[0] / window_min
        centres = 0.5 * (edges[:-1] + edges[1:]) / 60
        ax.step(centres, counts, where="mid", color=EPOCH_COLORS.get(epoch, INK),
                linewidth=1.8, label=f"{epoch} ({len(events)} events)")
    for epoch in ("PRE", "MAZE", "POST"):
        ax.axvline(session.epoch(epoch)[0] / 60, color=GRID, linewidth=1, zorder=0)
    ax.set_xlabel("time in session (min)", color=INK, fontsize=9)
    ax.set_ylabel("events per minute", color=INK, fontsize=9)
    ax.set_title(f"{session.name} — candidate-event rate", color=INK, fontsize=11, loc="left")
    ax.legend(fontsize=8, frameon=False, ncols=3, labelcolor=INK)
    _style(ax)
    return _save(fig, path, params)


def plot_participation(by_epoch, n_cells, session_name, path: Path, params: dict) -> Path:
    """Active-cell counts per event, per epoch, with the decodability cut marked."""
    fig, axes = plt.subplots(1, 2, figsize=(10, 3.4))
    bins = np.arange(0, n_cells + 2) - 0.5
    for epoch, (events, part) in by_epoch.items():
        if len(events) == 0:
            continue
        color = EPOCH_COLORS.get(epoch, INK)
        axes[0].hist(part.n_active, bins=bins, histtype="step", linewidth=1.8, color=color,
                     density=True, label=f"{epoch} (n={len(events)})")
        counts = np.sort(part.n_active)
        axes[1].plot(counts, 1 - np.linspace(0, 1, counts.size), color=color, linewidth=1.8,
                     label=epoch)
    axes[0].set_xlabel(f"active pyramidal cells per event (of {n_cells})", color=INK, fontsize=9)
    axes[0].set_ylabel("fraction of events", color=INK, fontsize=9)
    axes[0].set_title(f"{session_name} — participation", color=INK, fontsize=11, loc="left")
    axes[0].legend(fontsize=8, frameon=False, labelcolor=INK)
    axes[1].axvline(5, color=MUTED, linestyle="--", linewidth=1.2)
    axes[1].text(5, 0.95, "  5 cells", color=MUTED, fontsize=8)
    axes[1].set_xlabel("active cells per event", color=INK, fontsize=9)
    axes[1].set_ylabel("fraction of events at least this large", color=INK, fontsize=9)
    axes[1].set_xscale("symlog", linthresh=10)
    axes[1].set_title("survival (log x) — is the tail lognormal-ish?", color=INK,
                      fontsize=10, loc="left")
    axes[1].legend(fontsize=8, frameon=False, labelcolor=INK)
    for ax in axes:
        _style(ax)
    fig.tight_layout()
    return _save(fig, path, params)


# ---------------------------------------------------------------------------
# LFP channel selection (step 6a)
# ---------------------------------------------------------------------------


def plot_channel_power_profile(profile, meta, choices, session_name, path: Path,
                               params: dict) -> Path:
    """Ripple-band power per channel, with the per-shank winners marked.

    Shanks are not given 14 colours: they are an ordered spatial sequence, so
    they are separated along the x axis instead and labelled on the top edge.
    The lower panel collapses the shanks onto depth within the shank, which is
    where the pyramidal layer shows as a peak if the selection is picking up
    anatomy rather than noise.
    """
    chosen = {c.channel for c in choices}
    probe = set(meta.probe_channels)
    order = list(profile.channels)
    x = np.arange(len(order))
    ripple = profile.ripple

    fig, (top, bottom) = plt.subplots(
        2, 1, figsize=(11, 7), gridspec_kw={"height_ratios": [2, 1.4]}
    )

    top.plot(x, ripple, color=MUTED, linewidth=0.9, zorder=2)
    non_probe = [i for i, c in enumerate(order) if c not in probe]
    if non_probe:
        top.scatter([x[i] for i in non_probe], ripple[non_probe], s=26, marker="x",
                    color="#CC79A7", linewidth=1.2, zorder=4,
                    label="not on a sorted shank (ineligible)")
    picked = [i for i, c in enumerate(order) if c in chosen]
    if picked:
        top.scatter([x[i] for i in picked], ripple[picked], s=46, color="#0072B2",
                    zorder=5, label="selected for this shank")
    for shank in sorted(meta.shanks):
        channels = [c for c in meta.shanks[shank] if c in set(order)]
        if not channels:
            continue
        left = order.index(channels[0]) - 0.5
        right = order.index(channels[-1]) + 0.5
        top.axvline(right, color=GRID, linewidth=0.8, zorder=1)
        top.text((left + right) / 2, 1.02, str(shank), transform=top.get_xaxis_transform(),
                 ha="center", va="bottom", fontsize=7, color=MUTED)
    top.text(-0.045, 1.015, "shank", transform=top.transAxes, fontsize=7, color=MUTED,
             ha="right")
    top.set_xlabel("channel (file numbering, in shank order)")
    top.set_ylabel("mean 150-300 Hz envelope\n(a.d. units)")
    top.set_xlim(-1, len(order))
    top.set_xticks(x[::8])
    top.set_xticklabels([str(order[i]) for i in range(0, len(order), 8)])
    top.legend(frameon=False, fontsize=8, loc="upper right")
    top.set_title(
        f"{session_name}: ripple-band power by channel during non-REM "
        f"({profile.n_samples / profile.rate_hz / 60:.0f} min sampled)",
        color=INK, fontsize=10, pad=26,
    )
    _style(top)

    selected_shanks = {c.shank for c in choices}
    for shank in sorted(meta.shanks):
        channels = [c for c in meta.shanks[shank] if c in set(order)]
        if len(channels) < 2:
            continue
        values = [ripple[profile.index_of(c)] for c in channels]
        # Shanks left out of the selection are drawn, so that what was excluded
        # is visible, but dotted and unmarked so they cannot be misread as
        # candidates.
        included = shank in selected_shanks
        bottom.plot(range(len(channels)), values, color=MUTED, linewidth=0.8,
                    alpha=0.55 if included else 0.4,
                    linestyle="-" if included else ":", zorder=2)
        if not included:
            continue
        picked_channel = next(c.channel for c in choices if c.shank == shank)
        depth = channels.index(picked_channel)
        bottom.scatter([depth], [values[depth]], s=22, color="#0072B2", zorder=4)
    bottom.set_xlabel("position within shank (0 = first channel listed)")
    bottom.set_ylabel("mean 150-300 Hz\nenvelope")
    bottom.set_title(
        "collapsed onto depth within shank; one line per shank "
        "(dotted = no sorted units, not selected from)",
        color=INK, fontsize=9,
    )
    _style(bottom)

    fig.tight_layout()
    return _save(fig, path, params)


def plot_band_profile(profile, meta, choices, session_name, path: Path,
                      params: dict) -> Path:
    """Each band's mean envelope per channel, normalised to its own maximum.

    Ripple power alone cannot distinguish the pyramidal layer from a channel
    that is simply noisy; a noisy channel is elevated in every band at once,
    and this panel makes that visible.
    """
    order = list(profile.channels)
    x = np.arange(len(order))
    chosen = {c.channel for c in choices}
    fig, ax = plt.subplots(figsize=(11, 4.2))
    styles = {"delta": (":", MUTED), "theta": ("--", "#009E73"), "ripple": ("-", "#0072B2")}
    for name in ("delta", "theta", "ripple"):
        if name not in profile.mean_envelope:
            continue
        values = profile.mean_envelope[name]
        style, colour = styles[name]
        low, high = profile.bands[name]
        ax.plot(x, values / values.max(), style, color=colour, linewidth=1.1,
                label=f"{name} ({low:g}-{high:g} Hz)")
    picked = [i for i, c in enumerate(order) if c in chosen]
    ax.scatter([x[i] for i in picked], np.full(len(picked), 1.04), marker="v", s=22,
               color="#0072B2", clip_on=False, label="selected")
    ax.set_xlabel("channel (file numbering, in shank order)")
    ax.set_ylabel("mean envelope\n(fraction of band maximum)")
    ax.set_xlim(-1, len(order))
    ax.set_ylim(0, 1.1)
    ax.legend(frameon=False, fontsize=8, ncol=4, loc="lower center",
              bbox_to_anchor=(0.5, -0.42))
    ax.set_title(f"{session_name}: band profiles across channels", color=INK, fontsize=10)
    _style(ax)
    fig.tight_layout()
    return _save(fig, path, params)


def plot_ripple_examples(examples, rate_hz, band, session_name, channel, path: Path,
                         params: dict) -> Path:
    """Raw and band-passed LFP around the strongest ripple-band peaks.

    ``examples`` is a list of ``(centre_time_s, raw, filtered)`` triples, each
    trace the same length and centred on its peak. Raw sits above the filtered
    trace in each panel at its own scale, since a 40 uV ripple on a 1 mV
    sharp wave is invisible if the two share an axis.
    """
    n = len(examples)
    columns = min(3, n)
    rows = int(np.ceil(n / columns))
    fig, axes = plt.subplots(rows, columns, figsize=(4.0 * columns, 2.4 * rows),
                             squeeze=False)
    for index, ax in enumerate(axes.ravel()):
        if index >= n:
            ax.axis("off")
            continue
        centre, raw, filtered = examples[index]
        t = (np.arange(len(raw)) - len(raw) // 2) / rate_hz * 1000.0
        raw = np.asarray(raw, dtype=float)
        raw = raw - raw.mean()
        span = max(np.abs(raw).max(), 1.0)
        offset = -1.6 * span
        ax.plot(t, raw, color=INK, linewidth=0.7)
        ax.plot(t, offset + filtered * (span / max(np.abs(filtered).max(), 1e-9)) * 0.8,
                color="#0072B2", linewidth=0.7)
        ax.axvline(0, color=MUTED, linewidth=0.6, linestyle=":")
        ax.set_title(f"t = {centre:.2f} s", color=INK, fontsize=8)
        ax.set_yticks([])
        ax.set_xlabel("ms from peak", fontsize=8)
        if index % columns == 0:
            ax.text(-0.02, 0.78, "raw", transform=ax.transAxes, ha="right", fontsize=7,
                    color=INK)
            ax.text(-0.02, 0.18, f"{band[0]:g}-{band[1]:g} Hz", transform=ax.transAxes,
                    ha="right", fontsize=7, color="#0072B2")
        _style(ax)
        ax.grid(False)
    fig.suptitle(f"{session_name}: strongest ripple-band peaks on channel {channel}",
                 color=INK, fontsize=10)
    fig.tight_layout()
    return _save(fig, path, params)


# ---------------------------------------------------------------------------
# Ripple detection (step 6b)
# ---------------------------------------------------------------------------

#: Epochs are categorical and ordered; a sequential ramp would imply the gap
#: between PRE and MAZE is the same kind of thing as the gap between MAZE and
#: POST, which it is not. Three fixed Okabe-Ito hues instead.
EPOCH_COLORS = {"PRE": "#56B4E9", "MAZE": "#E69F00", "POST": "#009E73"}
EPOCH_ORDER = ("PRE", "MAZE", "POST")


def plot_ripple_rate_by_epoch(summary_rows, session_name, path: Path, params: dict) -> Path:
    """Event rate per minute of searched state, by epoch and channel ranking.

    Rate per minute of *state*, never raw counts: PRE holds roughly twice the
    non-REM of POST in these sessions, so counts would compare sleep duration
    rather than ripple incidence.
    """
    rows = [r for r in summary_rows if r["session"] == session_name]
    rankings = sorted({r["ranking"] for r in rows})
    fig, ax = plt.subplots(figsize=(1.6 + 1.5 * len(rankings), 3.6))
    width = 0.8 / max(len(EPOCH_ORDER), 1)
    x = np.arange(len(rankings))
    for index, epoch in enumerate(EPOCH_ORDER):
        values = [
            next((r["rate_per_min"] for r in rows
                  if r["ranking"] == ranking and r["epoch"] == epoch), np.nan)
            for ranking in rankings
        ]
        ax.bar(x + index * width - 0.4 + width / 2, values, width * 0.9,
               color=EPOCH_COLORS[epoch], label=epoch)
    ax.set_xticks(x)
    ax.set_xticklabels(
        [f"{r}\nch {next(q['channel'] for q in rows if q['ranking'] == r)}" for r in rankings],
        fontsize=8,
    )
    ax.set_ylabel("events per minute of searched state")
    ax.set_title(f"{session_name}: ripple rate by epoch", color=INK, fontsize=10)
    ax.legend(frameon=False, fontsize=8, ncol=3)
    _style(ax)
    fig.tight_layout()
    return _save(fig, path, params)


def plot_disagreement_examples(examples, rate_hz, band, session_name, channel_a, channel_b,
                               path: Path, params: dict) -> Path:
    """Both candidate channels around events that only one of them detected.

    Each panel is the same stretch of time on both channels, so whether the
    channel that missed the event missed a ripple or correctly ignored noise is
    visible rather than asserted.
    """
    n = len(examples)
    columns = min(3, n)
    rows = int(np.ceil(n / columns))
    fig, axes = plt.subplots(rows, columns, figsize=(4.2 * columns, 2.9 * rows), squeeze=False)
    for index, ax in enumerate(axes.ravel()):
        if index >= n:
            ax.axis("off")
            continue
        example = examples[index]
        traces = [
            (example["raw_a"], example["filt_a"], f"ch {channel_a}", "#0072B2"),
            (example["raw_b"], example["filt_b"], f"ch {channel_b}", "#D55E00"),
        ]
        spans = [max(np.abs(raw - raw.mean()).max(), 1.0) for raw, _, _, _ in traces]
        step = 2.4 * max(spans)
        for row, (raw, filtered, label, colour) in enumerate(traces):
            t = (np.arange(len(raw)) - len(raw) // 2) / rate_hz * 1000.0
            offset = -row * step
            ax.plot(t, offset + (raw - raw.mean()), color=colour, linewidth=0.7)
            scale = 0.45 * step / max(np.abs(filtered).max(), 1e-9)
            ax.plot(t, offset - 0.55 * step + filtered * scale, color=colour,
                    linewidth=0.7, alpha=0.75)
            ax.text(-0.02, 1 - (row + 0.25) / len(traces), label, transform=ax.transAxes,
                    ha="right", fontsize=7, color=colour)
        ax.axvline(0, color=MUTED, linewidth=0.6, linestyle=":")
        owner = channel_a if example["found_on"] == "a" else channel_b
        ax.set_title(
            f"{example['epoch']} t={example['centre']:.2f} s - only ch {owner} "
            f"(z {example['peak_z']:.1f})",
            color=INK, fontsize=8,
        )
        ax.set_yticks([])
        ax.set_xlabel("ms from peak", fontsize=8)
        _style(ax)
        ax.grid(False)
    fig.suptitle(
        f"{session_name}: events found on one candidate channel and not the other "
        f"(raw above, {band[0]:g}-{band[1]:g} Hz below)",
        color=INK, fontsize=10,
    )
    fig.tight_layout()
    return _save(fig, path, params)


def plot_shank_support(support_rows, session_name, path: Path, params: dict) -> Path:
    """How many other CA1 shanks show each event, as a distribution per epoch."""
    fig, ax = plt.subplots(figsize=(6.4, 3.8))
    for row in support_rows:
        counts = np.array([
            row[key] for key in sorted(
                (k for k in row if k.startswith("n_with_")),
                key=lambda k: int(k.rsplit("_", 1)[1]),
            )
        ], dtype=float)
        if counts.sum() == 0:
            continue
        ax.plot(np.arange(len(counts)), counts / counts.sum(), marker="o", markersize=4,
                linewidth=1.2, color=EPOCH_COLORS.get(row["epoch"], MUTED),
                label=f"{row['epoch']} (n={row['n_events']})")
    ax.set_xlabel("other CA1 shanks with an elevated ripple-band envelope")
    ax.set_ylabel("fraction of events")
    ax.set_title(
        f"{session_name}: multi-shank consistency of detected events",
        color=INK, fontsize=10,
    )
    ax.legend(frameon=False, fontsize=8)
    _style(ax)
    fig.tight_layout()
    return _save(fig, path, params)


def plot_cross_check_offsets(offsets_by_epoch, session_name, window_ms, path: Path,
                             params: dict) -> Path:
    """Peak-to-peak offset between matched LFP and MUA events.

    Signed, LFP peak minus MUA peak, so a systematic lead or lag is visible
    rather than folded away by taking the absolute value — the two detectors
    measure different things and are not expected to peak at the same instant.
    """
    fig, ax = plt.subplots(figsize=(6.4, 3.6))
    edges = np.linspace(-window_ms, window_ms, 41)
    for epoch, offsets in offsets_by_epoch.items():
        offsets = np.asarray(offsets, dtype=float)
        if offsets.size == 0:
            continue
        ax.hist(offsets, bins=edges, histtype="step", linewidth=1.3, density=True,
                color=EPOCH_COLORS.get(epoch, MUTED),
                label=f"{epoch} (n={offsets.size}, median {np.median(offsets):+.1f} ms)")
    ax.axvline(0, color=MUTED, linewidth=0.8, linestyle=":")
    ax.set_xlabel("MUA peak minus ripple peak (ms)")
    ax.set_ylabel("density")
    ax.set_title(f"{session_name}: timing of matched events", color=INK, fontsize=10)
    ax.legend(frameon=False, fontsize=8)
    _style(ax)
    fig.tight_layout()
    return _save(fig, path, params)


def plot_match_by_support(support_rows, session_name, path: Path, params: dict) -> Path:
    """Agreement with the MUA detector as a function of multi-shank support.

    The question this answers: are the ripples that only one or two shanks see
    the same ones the population-rate detector does not corroborate? A rising
    curve says the two independent checks agree about which events are weak.
    """
    fig, ax = plt.subplots(figsize=(6.4, 3.8))
    for epoch in EPOCH_ORDER:
        rows = sorted(
            (r for r in support_rows if r["epoch"] == epoch),
            key=lambda r: r["n_supporting_shanks"],
        )
        if not rows:
            continue
        x = [r["n_supporting_shanks"] for r in rows]
        y = [100 * r["fraction_matched_by_mua"] for r in rows]
        sizes = np.array([r["n_events"] for r in rows], dtype=float)
        ax.plot(x, y, linewidth=1.2, color=EPOCH_COLORS.get(epoch, MUTED), label=epoch,
                zorder=2)
        ax.scatter(x, y, s=8 + 60 * sizes / sizes.max(),
                   color=EPOCH_COLORS.get(epoch, MUTED), zorder=3, alpha=0.85)
    ax.set_xlabel("other CA1 shanks showing the event")
    ax.set_ylabel("% of ripples matched by an MUA event")
    ax.set_ylim(0, 102)
    ax.set_title(
        f"{session_name}: two independent checks, compared\n"
        "(marker area is the number of events at that support level)",
        color=INK, fontsize=10,
    )
    ax.legend(frameon=False, fontsize=8, loc="lower right")
    _style(ax)
    fig.tight_layout()
    return _save(fig, path, params)


def plot_decodable_by_variant(decode_rows, session_name, path: Path, params: dict) -> Path:
    """Events with enough place cells to decode, full MUA set against the conjunction.

    Counts, not percentages: the question this answers is what step 7 has left
    to work with, and a percentage of a smaller set hides the cost.
    """
    variants = ["mua", "conjunction"]
    labels = {"mua": "MUA (full)", "conjunction": "conjunction\n(MUA + ripple)"}
    epochs = [e for e in EPOCH_ORDER
              if any(r["epoch"] == e and r["variant"] in variants for r in decode_rows)]
    fig, ax = plt.subplots(figsize=(1.8 + 1.4 * len(epochs), 3.8))
    width = 0.36
    x = np.arange(len(epochs))
    for index, variant in enumerate(variants):
        total = [
            next((r["n_events"] for r in decode_rows
                  if r["epoch"] == e and r["variant"] == variant), 0) for e in epochs
        ]
        usable = [
            next((r["n_decodable"] for r in decode_rows
                  if r["epoch"] == e and r["variant"] == variant), 0) for e in epochs
        ]
        offset = (index - 0.5) * width
        colour = "#0072B2" if variant == "mua" else "#009E73"
        ax.bar(x + offset, total, width * 0.92, color=colour, alpha=0.28,
               label=f"{labels[variant]}: all events" if index == 0 else None)
        ax.bar(x + offset, usable, width * 0.92, color=colour,
               label=f"{labels[variant]}: enough cells to decode")
        for position, (whole, part) in zip(x + offset, zip(total, usable, strict=True),
                                           strict=True):
            if whole:
                ax.text(position, part, f"{100 * part / whole:.0f}%", ha="center",
                        va="bottom", fontsize=7, color=INK)
    ax.set_xticks(x)
    ax.set_xticklabels(epochs)
    ax.set_ylabel("events")
    ax.set_title(
        f"{session_name}: how many candidate events are usable\n"
        "(pale = all events, solid = enough active place cells)",
        color=INK, fontsize=10,
    )
    ax.legend(frameon=False, fontsize=7.5, loc="upper right")
    _style(ax)
    fig.tight_layout()
    return _save(fig, path, params)


def plot_shank_profile(profile_rows, session_name, path: Path, params: dict) -> Path:
    """Where the active cells sat, for MUA events with and without a ripple partner.

    Each point is one shank's mean share of the active cells in an event. The
    grey line is the share expected from how many cells each shank has, which
    is the only thing the two series would track if the unmatched events were
    just smaller versions of the matched ones.
    """
    present = {r["epoch"] for r in profile_rows}
    epochs = [e for e in EPOCH_ORDER if e in present]
    fig, axes = plt.subplots(1, len(epochs), figsize=(4.6 * len(epochs), 3.6),
                             squeeze=False, sharey=True)
    for ax, epoch in zip(axes[0], epochs, strict=True):
        rows = sorted((r for r in profile_rows if r["epoch"] == epoch),
                      key=lambda r: r["shank"])
        x = np.arange(len(rows))
        ax.plot(x, [r["expected_share"] for r in rows], color=MUTED, linewidth=1.0,
                linestyle="--", label="expected from cell counts")
        ax.plot(x, [r["share_matched"] for r in rows], marker="o", markersize=4,
                linewidth=1.2, color="#0072B2", label="MUA with a ripple")
        ax.plot(x, [r["share_unmatched"] for r in rows], marker="s", markersize=4,
                linewidth=1.2, color="#D55E00", label="MUA without a ripple")
        for index, row in enumerate(rows):
            if row["is_primary_shank"]:
                ax.axvline(index, color=GRID, linewidth=6, zorder=0)
                ax.text(index, ax.get_ylim()[1], "detection\nshank", ha="center",
                        va="top", fontsize=7, color=MUTED)
        ax.set_xticks(x)
        ax.set_xticklabels([r["shank"] for r in rows], fontsize=7)
        ax.set_xlabel("shank")
        ax.set_title(epoch, color=INK, fontsize=9)
        _style(ax)
    axes[0][0].set_ylabel("mean share of an event's active cells")
    axes[0][-1].legend(frameon=False, fontsize=7.5)
    fig.suptitle(f"{session_name}: where the active cells sat", color=INK, fontsize=10)
    fig.tight_layout()
    return _save(fig, path, params)


def plot_sequence_cdf(scores, session_name, set_label, threshold, path: Path,
                      params: dict) -> Path:
    """Cumulative distribution of rZ by epoch, one panel per shuffle kind.

    A cumulative distribution rather than a histogram: the comparison is a
    question about whether one distribution sits to the right of another, and
    a CDF answers it by eye without a bin width to argue about. Medians are
    marked because the test reported alongside is on medians and ranks.
    """
    kinds = [k for k in ("column_cycle", "time_permutation")
             if f"rz_{k}" in scores.columns]
    fig, axes = plt.subplots(1, len(kinds), figsize=(5.4 * len(kinds), 4.0),
                             squeeze=False, sharey=True)
    for ax, kind in zip(axes[0], kinds, strict=True):
        for epoch in EPOCH_ORDER:
            values = scores[scores.epoch == epoch][f"rz_{kind}"].dropna().to_numpy()
            if values.size < 2:
                continue
            ordered = np.sort(values)
            fraction = np.arange(1, ordered.size + 1) / ordered.size
            colour = EPOCH_COLORS[epoch]
            ax.plot(ordered, fraction, linewidth=1.4, color=colour,
                    label=f"{epoch} (n={ordered.size}, median {np.median(values):+.2f})")
            median = float(np.median(values))
            ax.plot([median], [0.5], marker="o", markersize=5, color=colour, zorder=4)
        ax.axvline(threshold, color=MUTED, linewidth=0.8, linestyle=":")
        ax.text(threshold, 0.02, f" z={threshold:g}", fontsize=7, color=MUTED,
                ha="left", va="bottom")
        ax.set_xlabel("sequence score rZ")
        ax.set_xlim(-3, 8)
        ax.set_title(kind.replace("_", " "), color=INK, fontsize=9)
        ax.legend(frameon=False, fontsize=7.5, loc="lower right")
        _style(ax)
    axes[0][0].set_ylabel("cumulative fraction of events")
    fig.suptitle(f"{session_name}: sequence scores, {set_label} events",
                 color=INK, fontsize=10)
    fig.tight_layout()
    return _save(fig, path, params)



def plot_sensitivity_floor(rows, path: Path, params: dict) -> Path:
    """Positive-control score against ensemble size — where the method stops working.

    Each point is one session's median rZ for 2 s windows of *real running*, a
    trajectory the animal demonstrably ran and the decoder demonstrably
    recovers. It is the cleanest statement in the project of how much the
    method depends on how many cells it has: the question is not whether these
    windows contain a sequence, it is whether the scoring can still see one.

    Plotted against place cells rather than against session identity, because
    the ensemble is the variable that moves. Circular sessions are drawn hollow:
    a ring maximises the correlation over 28-29 position offsets where a line
    maximises over one, so the circular null is broader and rZ compressed
    (D8.4). Ensemble size and geometry are therefore not fully separable in
    these five sessions, and the trend should be read as a floor, not a
    calibration curve.
    """
    fig, ax = plt.subplots(figsize=(6.8, 4.4))
    for kind, colour, marker in (("column_cycle", "#0072B2", "o"),
                                 ("time_permutation", "#D55E00", "s")):
        here = sorted((r for r in rows if r["shuffle"] == kind),
                      key=lambda r: r["n_place_cells"])
        if not here:
            continue
        x = [r["n_place_cells"] for r in here]
        y = [r["median_rz"] for r in here]
        ax.plot(x, y, linewidth=1.4, color=colour, zorder=2,
                label=kind.replace("_", " "))
        for r in here:
            circular = str(r.get("geometry", "")).lower().startswith("circ")
            ax.plot(r["n_place_cells"], r["median_rz"], marker=marker,
                    markersize=6.5, color=colour, zorder=3,
                    markerfacecolor="white" if circular else colour,
                    markeredgecolor=colour, markeredgewidth=1.4)
    labelled = sorted((r for r in rows if r["shuffle"] == "column_cycle"),
                      key=lambda r: r["n_place_cells"])
    for i, r in enumerate(labelled):
        crowded = any(abs(o["n_place_cells"] - r["n_place_cells"]) < 10
                      for j, o in enumerate(labelled) if j < i)
        ax.annotate(r["session"].replace("_", "\n"),
                    (r["n_place_cells"], r["median_rz"]), textcoords="offset points",
                    xytext=(0, -26 if crowded else 12), fontsize=6.5, color=MUTED,
                    ha="center", va="bottom" if not crowded else "bottom")
    ax.axhline(1.65, color=MUTED, linewidth=0.8, linestyle=":", zorder=1)
    ax.margins(x=0.09, y=0.12)
    ax.text(ax.get_xlim()[1], 1.65, " z = 1.65 ", fontsize=7, color=MUTED,
            ha="right", va="bottom")
    ax.set_xlabel("place cells in the decoding ensemble")
    ax.set_ylabel("median rZ, 2 s windows of real running")
    ax.set_title(
        "Where the sequence score stops working\n"
        "positive control against ensemble size, five sessions",
        color=INK, fontsize=10,
    )
    handles = [Line2D([], [], color="#0072B2", marker="o", markersize=6.5,
                      linewidth=1.4, label="column cycle"),
               Line2D([], [], color="#D55E00", marker="s", markersize=6.5,
                      linewidth=1.4, label="time permutation"),
               Line2D([], [], color=MUTED, marker="o", markersize=6.5, linewidth=0,
                      markerfacecolor="white", markeredgecolor=MUTED,
                      label="circular track (rZ compressed, D8.4)")]
    ax.legend(handles=handles, frameon=False, fontsize=7.5, loc="upper left")
    _style(ax)
    fig.tight_layout()
    return _save(fig, path, params)


#: Below this, a PRE/POST difference is too small for a residual *fraction* of it
#: to mean anything: the denominator is noise, so the ratio is noise over noise.
UNSTABLE_DIFFERENCE = 0.05


def plot_residual_vs_ensemble(rows, path: Path, params: dict) -> Path:
    """Identity-shuffle residual against ensemble size — does it track cells, not geometry?

    Two panels, because the quantity asked for and the quantity that behaves are
    not the same one. Left: the residual as a fraction of the real PRE/POST
    difference, which is how it is quoted everywhere else and which inherits its
    denominator's instability — points whose real difference is under
    :data:`UNSTABLE_DIFFERENCE` are ringed, not hidden. Right: the absolute
    residual in rZ units, which has no denominator at all.

    Circular sessions are drawn hollow. **No line is fitted.** Five sessions
    cannot establish a relationship, and the figure exists to show the ordering,
    not to measure a slope.
    """
    fig, axes = plt.subplots(1, 2, figsize=(10.4, 4.3))
    panels = (
        (axes[0], "residual_fraction", "residual ÷ real PRE−POST difference", True),
        (axes[1], "abs_residual", "|identity-shuffle residual| (rZ)", False),
    )
    for ax, key, ylabel, mark_unstable in panels:
        for kind, colour, marker in (("column_cycle", "#0072B2", "o"),
                                     ("time_permutation", "#D55E00", "s")):
            here = sorted((r for r in rows if r["shuffle"] == kind),
                          key=lambda r: r["n_place_cells"])
            if not here:
                continue
            # Markers only: a connecting line between five sessions would read as
            # a trend, and no trend is claimed.
            for r in here:
                circular = str(r.get("geometry", "")).lower().startswith("circ")
                ax.plot(r["n_place_cells"], r[key], marker=marker, markersize=6.5,
                        color=colour, zorder=3, linestyle="none",
                        markerfacecolor="white" if circular else colour,
                        markeredgecolor=colour, markeredgewidth=1.4)
                if mark_unstable and abs(r["difference"]) < UNSTABLE_DIFFERENCE:
                    ax.plot(r["n_place_cells"], r[key], marker="o", markersize=13,
                            markerfacecolor="none", markeredgecolor=MUTED,
                            markeredgewidth=0.9, linestyle="none", zorder=4)
        ax.axhline(0.0, color=MUTED, linewidth=0.8, linestyle=":", zorder=1)
        ax.set_xlabel("place cells in the decoding ensemble")
        ax.set_ylabel(ylabel)
        ax.margins(x=0.12, y=0.14)
        _style(ax)
    for r in rows:
        if r["shuffle"] != "column_cycle":
            continue
        for ax, key in ((axes[0], "residual_fraction"), (axes[1], "abs_residual")):
            ax.annotate(r["session"].replace("_", "\n"),
                        (r["n_place_cells"], r[key]), textcoords="offset points",
                        xytext=(0, 11), fontsize=6.5, color=MUTED, ha="center")
    handles = [Line2D([], [], color="#0072B2", marker="o", markersize=6.5, linewidth=0,
                      label="column cycle"),
               Line2D([], [], color="#D55E00", marker="s", markersize=6.5, linewidth=0,
                      label="time permutation"),
               Line2D([], [], color=MUTED, marker="o", markersize=6.5, linewidth=0,
                      markerfacecolor="white", markeredgecolor=MUTED,
                      label="circular track"),
               Line2D([], [], color=MUTED, marker="o", markersize=11, linewidth=0,
                      markerfacecolor="none", markeredgecolor=MUTED,
                      label="real difference < 0.05: fraction not interpretable")]
    fig.suptitle("Does the negative control's residual follow ensemble size?\n"
                 "five sessions, no line fitted", color=INK, fontsize=10)
    fig.tight_layout(rect=(0, 0.09, 1, 1))
    fig.legend(handles=handles, frameon=False, fontsize=7.5, ncol=2,
               loc="lower center", bbox_to_anchor=(0.5, 0.0))
    return _save(fig, path, params)


def plot_firing_rate_split(rows, null_rows, path: Path, params: dict) -> Path:
    """The paper's Fig. 4D question, with the control that makes it readable.

    Left column: the POST-minus-PRE change in median rZ for the slow-firing and
    fast-firing halves of each session's ensemble, scored on the same events.
    Right column: the gap between those two halves against the spread of gaps
    produced by splitting the same ensemble at random, drawn as a grey band.
    A rate split that lands inside its own random band has not shown anything,
    however large the gap looks on the left.
    """
    kinds = [k for k in ("column_cycle", "time_permutation")
             if any(r["shuffle"] == k for r in rows)]
    fig, axes = plt.subplots(len(kinds), 2, figsize=(11.0, 2.9 * len(kinds) + 1.1),
                             squeeze=False)
    sessions = sorted({r["session"] for r in rows},
                      key=lambda s: -max(r["n_cells"] for r in rows if r["session"] == s))
    y = np.arange(len(sessions))
    for row_index, kind in enumerate(kinds):
        here = {r["session"]: r for r in rows if r["shuffle"] == kind}
        left, right = axes[row_index]

        for label, key, colour, filled in (("slow half", "slow_difference", "#0072B2", True),
                                           ("fast half", "fast_difference", "#D55E00", False)):
            left.plot([here[s][key] for s in sessions], y, marker="o", markersize=7,
                      linestyle="none", color=colour, label=label,
                      markerfacecolor=colour if filled else "white",
                      markeredgecolor=colour, markeredgewidth=1.5)
        left.axvline(0.0, color=MUTED, linewidth=0.8, linestyle=":")
        left.set_xlabel("POST − PRE, median rZ (same events)")
        left.set_yticks(y)
        left.set_yticklabels([s.replace("_", "\n") for s in sessions], fontsize=7)
        left.set_title(f"{kind.replace('_', ' ')}: each half", fontsize=9, color=INK)
        if row_index == 0:
            left.legend(frameon=False, fontsize=8, loc="best")

        for i, session in enumerate(sessions):
            draws = [abs(n["gap"]) for n in null_rows
                     if n["session"] == session and n["shuffle"] == kind]
            if draws:
                right.barh(i, 2 * np.percentile(draws, 95), left=-np.percentile(draws, 95),
                           height=0.45, color=MUTED, alpha=0.22, zorder=1)
            right.plot(here[session]["gap_slow_minus_fast"], i, marker="D", markersize=7,
                       color="#0072B2", zorder=3, linestyle="none")
        right.axvline(0.0, color=MUTED, linewidth=0.8, linestyle=":")
        right.set_xlabel("gap: slow minus fast")
        right.set_yticks(y)
        right.set_yticklabels([])
        right.set_title(f"{kind.replace('_', ' ')}: against random splits",
                        fontsize=9, color=INK)
        if row_index == 0 and null_rows:
            right.text(0.98, 0.04, "grey: 95th percentile of |gap| over random halves",
                       transform=right.transAxes, fontsize=7, color=MUTED, ha="right")
        for ax in (left, right):
            ax.margins(y=0.18)
            _style(ax)
    fig.suptitle("Do slow-firing cells gain replay from PRE to POST, and fast-firing not?\n"
                 "the paper's Fig. 4D, with a random-split control",
                 color=INK, fontsize=10)
    fig.tight_layout()
    return _save(fig, path, params)


def plot_pcc_null_first(real_rows, null_rows, path: Path, params: dict) -> Path:
    """The ΔPCC the estimator produces from noise, against the one it produces from data.

    Null first, in grey, and drawn behind: the question this figure answers is
    not how large the real spread is but whether it can be told apart from the
    spread of ΔPCC values obtained by treating a shifted field as if it were the
    cell's own. If the two distributions sit on top of each other, every
    rigid/plastic comparison downstream is a comparison of noise, and no later
    choice can repair that.
    """
    kinds = [k for k in ("column_cycle", "time_permutation")
             if any(r["shuffle"] == k for r in real_rows)]
    fig, axes = plt.subplots(1, len(kinds), figsize=(5.4 * len(kinds), 4.2), squeeze=False)
    for ax, kind in zip(axes[0], kinds, strict=True):
        real = np.array([r["delta_pcc"] for r in real_rows if r["shuffle"] == kind],
                        dtype=float)
        null = np.array([r["delta_pcc"] for r in null_rows if r["shuffle"] == kind],
                        dtype=float)
        real, null = real[np.isfinite(real)], null[np.isfinite(null)]
        span = np.percentile(np.concatenate([real, null]), [1, 99])
        edges = np.linspace(span[0], span[1], 41)
        ax.hist(null, bins=edges, density=True, color=MUTED, alpha=0.35,
                label=f"null labels (shifted fields, n={null.size})")
        ax.hist(real, bins=edges, density=True, histtype="step", linewidth=1.6,
                color="#0072B2", label=f"real fields (n={real.size})")
        for values, colour, style in ((null, MUTED, "--"), (real, "#0072B2", "-")):
            ax.axvline(float(np.median(values)), color=colour, linewidth=1.0,
                       linestyle=style)
        ax.set_xlabel("ΔPCC  (POST − PRE per-cell contribution)")
        ax.set_ylabel("density")
        ax.set_title(kind.replace("_", " "), fontsize=9, color=INK)
        ax.legend(frameon=False, fontsize=7.5, loc="upper left")
        _style(ax)
    fig.suptitle("What the estimator produces from noise, and from data\n"
                 "null labels first: if these overlap, the rigid/plastic split is noise",
                 color=INK, fontsize=10)
    fig.tight_layout()
    return _save(fig, path, params)


def plot_pcc_cross_session(summary_rows, control_rows, path: Path, params: dict) -> Path:
    """What the per-cell contribution survives across five sessions.

    Left: the PRE-to-POST contributor overlap the paper reports as 73 %, against
    the size of the ensemble it was measured in. It is not a constant of the
    hippocampus in our data; it rises as the ensemble shrinks, which is what a
    threshold on a noisy per-cell statistic does.

    Right: every rigid/plastic property difference against the largest effect
    the *null labels* reach in the same session. Points below the diagonal are
    differences a median split of noise can already produce.
    """
    fig, axes = plt.subplots(1, 2, figsize=(10.6, 4.4))

    linear = {"Achilles_10252013", "Gatsby_08022013", "Cicero_09012014"}
    ax = axes[0]
    ax.axhline(73, color=MUTED, linewidth=0.9, linestyle="--")
    ax.text(0.99, 0.955, " paper: 73 % ", transform=ax.transAxes, fontsize=7.5,
            color=MUTED, ha="right", va="top")
    for row in summary_rows:
        if not np.isfinite(row.get("overlap_pct", np.nan)):
            continue
        circular = row["session"] not in linear
        marker = "o" if row["shuffle"] == "column_cycle" else "s"
        colour = "#0072B2" if row["shuffle"] == "column_cycle" else "#D55E00"
        ax.plot(row["n_cells"], row["overlap_pct"], marker=marker, markersize=7,
                linestyle="none", color=colour,
                markerfacecolor="white" if circular else colour,
                markeredgecolor=colour, markeredgewidth=1.5)
        if row["shuffle"] == "column_cycle":
            ax.annotate(f"{row['session'].split('_')[0]}\n{row['n_contributing_pre']} cells",
                        (row["n_cells"], row["overlap_pct"]), textcoords="offset points",
                        xytext=(0, 11), fontsize=6.5, color=MUTED, ha="center")
    ax.set_xlabel("place cells in the ensemble")
    ax.set_ylabel("PRE contributors that also contribute to POST (%)")
    ax.set_title("The overlap claim, against ensemble size", fontsize=9, color=INK)
    ax.margins(x=0.16, y=0.22)
    _style(ax)

    ax = axes[1]
    limit = 1.05
    ax.plot([0, limit], [0, limit], color=MUTED, linewidth=0.9, linestyle=":")
    for row in control_rows:
        ceiling = row.get("null_rank_biserial_p95", np.nan)
        if not np.isfinite(ceiling):
            continue
        beats = abs(row["real_rank_biserial"]) > ceiling
        ax.plot(ceiling, abs(row["real_rank_biserial"]), marker="D" if beats else "o",
                markersize=7 if beats else 5.5, linestyle="none",
                color="#D55E00" if beats else "#0072B2",
                markerfacecolor="#D55E00" if beats else "white",
                markeredgecolor="#D55E00" if beats else "#0072B2", markeredgewidth=1.3)
    ax.text(0.03, 0.95, "above the line: larger than noise labels reach\n"
            "below: a median split of noise does this too",
            transform=ax.transAxes, fontsize=7.5, color=MUTED, va="top")
    ax.set_xlabel("null labels' 95th percentile |rank-biserial|")
    ax.set_ylabel("real labels' |rank-biserial|")
    ax.set_title("Every rigid/plastic difference against its own null",
                 fontsize=9, color=INK)
    ax.set_xlim(0, limit)
    ax.set_ylim(0, limit)
    _style(ax)

    fig.suptitle("Per-cell contribution across five sessions: what survives its null",
                 color=INK, fontsize=10)
    fig.tight_layout()
    return _save(fig, path, params)


def plot_boundary_sweep(rows, path: Path, params: dict) -> Path:
    """What lowering the detector's boundary threshold buys, and what it costs.

    Four panels, two of them the purchase and two the price. Bins per scorable
    event and the scorable fraction both rise steeply; active cells per event
    does not move, and spike density per bin falls. A boundary that lifts the
    first pair while leaving the second flat is buying **time, not
    participants** — which matters, because the detection power of both
    decoders scales with how many cells are active in an event.

    PRE solid, POST dashed: the two must move together, or a change here would
    shift the PRE/POST comparison for a reason that is not biology.
    """
    fig, axes = plt.subplots(2, 2, figsize=(10.6, 7.4))
    sessions = sorted({r["session"] for r in rows})
    colours = ["#0072B2", "#D55E00", "#009E73", "#CC79A7", "#56B4E9"]
    panels = (
        (axes[0][0], "median_used_bins_scorable", "bins with spikes, scorable events"),
        (axes[0][1], "pct_scorable", "events clearing the five-bin floor (%)"),
        (axes[1][0], "median_active_pyr_scorable", "active pyramidal cells per event"),
        (axes[1][1], "pct_ripple_matched", "events matching a ripple (%)"),
    )
    for ax, key, ylabel in panels:
        for session, colour in zip(sessions, colours, strict=False):
            for epoch, style in (("PRE", "-"), ("POST", "--")):
                here = sorted((r for r in rows if r["session"] == session
                               and r["epoch"] == epoch and r.get("n_events")),
                              key=lambda r: -r["threshold_low_sd"])
                here = [r for r in here if np.isfinite(r.get(key, np.nan))]
                if not here:
                    continue
                ax.plot([r["threshold_low_sd"] for r in here], [r[key] for r in here],
                        style, color=colour, linewidth=1.3, marker="o", markersize=3.5,
                        label=session.replace("_", " ") if epoch == "PRE" else None)
        ax.axvline(1.0, color=MUTED, linewidth=0.9, linestyle=":")
        ax.set_xlabel("boundary threshold (SD above baseline)")
        ax.set_ylabel(ylabel)
        ax.invert_xaxis()
        _style(ax)
    axes[0][0].axhline(11.5, color=MUTED, linewidth=0.9, linestyle="-.")
    axes[0][0].text(0.02, 0.95, "-·- Chen et al., 11.5 bins", transform=axes[0][0].transAxes,
                    fontsize=7, color=MUTED, va="top")
    axes[1][0].text(0.02, 0.06, "flat: the extra duration adds time, not cells",
                    transform=axes[1][0].transAxes, fontsize=7, color=MUTED)
    axes[1][1].text(0.02, 0.06, "flat: the extra duration is still ripple-locked",
                    transform=axes[1][1].transAxes, fontsize=7, color=MUTED)
    axes[0][1].text(0.98, 0.06, "dotted: the current 1 SD\nsolid PRE, dashed POST",
                    transform=axes[0][1].transAxes, fontsize=7, color=MUTED, ha="right")
    axes[0][0].legend(frameon=False, fontsize=7, loc="lower right")
    fig.suptitle("Lowering the boundary threshold: what it buys and what it costs\n"
                 "peak threshold fixed at 3 SD throughout; no threshold adopted (D10.1)",
                 color=INK, fontsize=10)
    fig.tight_layout()
    return _save(fig, path, params)
