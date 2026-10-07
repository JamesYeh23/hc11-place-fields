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
