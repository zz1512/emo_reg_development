"""Fan plots and pairwise difference matrices for the two anatomical groupings."""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import nibabel as nib
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.cm import ScalarMappable
from nilearn import plotting

from analysis.fig03.gam_supplement.prepare_inputs import SYSTEM_ORDER, YEO_ORDER

PROJECT_ROOT = Path(__file__).resolve().parents[3]

SURFACE_TEMPLATE = PROJECT_ROOT / "resources" / "templates" / "fslr_32k"
SURFACES = {
    "lh": SURFACE_TEMPLATE / "Q1-Q6_R440.L.very_inflated.32k_fs_LR.surf.gii",
    "rh": SURFACE_TEMPLATE / "Q1-Q6_R440.R.very_inflated.32k_fs_LR.surf.gii",
}
DIVERGING = LinearSegmentedColormap.from_list(
    "signed_partial_r2",
    ["#3b0a74", "#8d6cab", "#f4f1ea", "#e6c36a", "#f6e27a"],
)
INCREASING_COLOR = "#3e9d57"
DECREASING_COLOR = "#e8872d"
INCREASING_CMAP = LinearSegmentedColormap.from_list("increasing_parcels", [INCREASING_COLOR, INCREASING_COLOR])
DECREASING_CMAP = LinearSegmentedColormap.from_list("decreasing_parcels", [DECREASING_COLOR, DECREASING_COLOR])
GROUP_ORDERS = {"yeo": YEO_ORDER, "system": SYSTEM_ORDER}
VIEWS = (("lh", "lateral"), ("rh", "lateral"), ("lh", "medial"), ("rh", "medial"))
SHORT_GROUP_LABELS = {
    "reappraisal_only": "reappraisal\nonly",
    "common_appraisal": "common\nappraisal",
    "nonmodifiable_emotion": "nonmodifiable\nemotion",
    "modifiable_emotion": "modifiable\nemotion",
}


def _dark_axes(ax: plt.Axes) -> None:
    ax.set_facecolor("black")
    ax.tick_params(colors="white", labelsize=8)
    for spine in ax.spines.values():
        spine.set_color("white")
    ax.xaxis.label.set_color("white")
    ax.yaxis.label.set_color("white")
    ax.title.set_color("white")


def _hemisphere_labels(dlabel_path: Path, hemi: str) -> np.ndarray:
    """Parcel labels on one fsLR hemisphere, in the surface vertex order."""
    image = nib.load(str(dlabel_path))
    data = np.asarray(image.get_fdata()).reshape(-1)
    structure = "CIFTI_STRUCTURE_CORTEX_LEFT" if hemi == "lh" else "CIFTI_STRUCTURE_CORTEX_RIGHT"
    brain_axis = image.header.get_axis(1)
    for name, row_slice, _sub_axis in brain_axis.iter_structures():
        if name == structure:
            return np.rint(data[row_slice]).astype(np.int32)
    raise ValueError(f"{dlabel_path} has no {structure} vertices")


def _vertex_map(stat_by_parcel: dict[int, float], dlabel_path: Path, hemi: str) -> np.ndarray:
    labels = _hemisphere_labels(dlabel_path, hemi)
    values = np.full(labels.shape, np.nan, dtype=float)
    for parcel, value in stat_by_parcel.items():
        if np.isfinite(value):
            values[labels == int(parcel)] = value
    return values


def _style_surface_axis(ax: plt.Axes) -> None:
    ax.set_facecolor("black")
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)


def _paint_surface(ax: plt.Axes, hemi: str, view: str, vertex_values: np.ndarray,
                   cmap, vmin: float, vmax: float) -> None:
    hemi_name = "left" if hemi == "lh" else "right"
    mesh = str(SURFACES[hemi])
    if np.isfinite(vertex_values).any():
        plotting.plot_surf_stat_map(
            mesh, vertex_values, hemi=hemi_name, view=view, cmap=cmap,
            vmin=vmin, vmax=vmax, colorbar=False, axes=ax, bg_map=None,
        )
    else:
        plotting.plot_surf(mesh, hemi=hemi_name, view=view, axes=ax, bg_map=None, colorbar=False)
    _style_surface_axis(ax)


def _draw_brains(fig: plt.Figure, spec, stat_by_parcel: dict[int, float], dlabel_path: Path,
                 vmax: float, colorbar_label: str = "partial R$^2$") -> None:
    block = spec.subgridspec(2, 1, height_ratios=[1.0, 0.07], hspace=0.08)
    inner = block[0].subgridspec(2, 2, wspace=0.02, hspace=0.02)
    for index, (hemi, view) in enumerate(VIEWS):
        ax = fig.add_subplot(inner[index // 2, index % 2], projection="3d")
        _paint_surface(
            ax, hemi, view, _vertex_map(stat_by_parcel, dlabel_path, hemi),
            DIVERGING, -vmax, vmax,
        )
    color_ax = fig.add_subplot(block[1])
    colorbar = fig.colorbar(
        ScalarMappable(norm=Normalize(vmin=-vmax, vmax=vmax), cmap=DIVERGING),
        cax=color_ax, orientation="horizontal",
    )
    colorbar.ax.xaxis.set_tick_params(color="white", labelcolor="white", labelsize=7)
    colorbar.set_label(colorbar_label, color="white", fontsize=8)
    colorbar.outline.set_edgecolor("white")


def _overlay_brains(fig: plt.Figure, host: plt.Axes, stat_by_parcel: dict[int, float],
                    dlabel_path: Path, cmap, corner: str) -> None:
    """Small lateral and medial views sitting on one trajectory panel."""
    fig.canvas.draw()
    pos = host.get_position()
    width = pos.width * 0.40
    height = pos.height * 0.46
    pad_x = pos.width * 0.035
    if corner == "lower right":
        x = pos.x1 - width - pad_x
        y = pos.y0 + pos.height * 0.10
    else:
        x = pos.x1 - width - pad_x
        y = pos.y1 - height - pos.height * 0.08
    cell_w = width * 0.48
    cell_h = height * 0.48
    gap_x = width * 0.04
    gap_y = height * 0.04
    for index, (hemi, view) in enumerate(VIEWS):
        column = index % 2
        row = index // 2
        ax = fig.add_axes(
            [
                x + column * (cell_w + gap_x),
                y + (1 - row) * (cell_h + gap_y),
                cell_w, cell_h,
            ],
            projection="3d",
        )
        _paint_surface(
            ax, hemi, view, _vertex_map(stat_by_parcel, dlabel_path, hemi),
            cmap, 0.0, 1.0,
        )


def _draw_curves(ax: plt.Axes, curves: pd.DataFrame, values: pd.Series, norm: Normalize,
                 *, title: str, color: str | None = None) -> None:
    for parcel, rows in curves.groupby("parcel", sort=False):
        magnitude = values.get(int(parcel), np.nan)
        line_color = color if color is not None else DIVERGING(norm(0.0 if not np.isfinite(magnitude) else magnitude))
        ax.plot(rows["age"], rows["centered_effect"], color=line_color, lw=1.0, alpha=0.9)
    ax.axhline(0, color="white", lw=0.6, alpha=0.7)
    ax.set_title(title, fontsize=10)
    ax.set_xlabel("Age")
    ax.set_ylabel("Centered age effect")
    _dark_axes(ax)


def plot_parcel_suite(dimension: str, parcel_info: pd.DataFrame, stats: pd.DataFrame,
                      curves: pd.DataFrame, dlabel_path: Path, output_path: Path) -> None:
    """One reference-style fan of the parcels that belong to this grouping."""
    info = parcel_info.copy()
    info["parcel"] = info["parcel"].astype(int)
    if dimension == "system":
        # A blank system cell is read back from CSV as NaN, and str(NaN) is "nan".
        label = info["system"].fillna("").astype(str).str.strip()
        info = info[~label.str.lower().isin({"", "nan", "none"})]
    keep = set(info["parcel"].astype(int))
    stats = stats[stats["parcel"].astype(int).isin(keep)].copy()
    curves = curves[curves["parcel"].astype(int).isin(keep)].copy()
    signed = stats.set_index(stats["parcel"].astype(int))["signed_partial_r2"]
    finite = signed.to_numpy(dtype=float)
    finite = finite[np.isfinite(finite)]
    vmax = float(np.max(np.abs(finite))) if finite.size else 1.0
    if vmax <= 0:
        vmax = 1.0
    norm = Normalize(vmin=-vmax, vmax=vmax)
    shapes = stats.set_index(stats["parcel"].astype(int))["shape"]
    increasing = shapes[shapes == "increasing"].index
    decreasing = shapes[shapes == "decreasing"].index

    selected = {int(parcel): float(value) for parcel, value in signed.items()}
    fig = plt.figure(figsize=(13.2, 8.4), facecolor="black")
    grid = fig.add_gridspec(2, 2, width_ratios=[1.0, 1.25], wspace=0.24, hspace=0.30,
                            left=0.07, right=0.86, top=0.90, bottom=0.07)
    _draw_brains(fig, grid[0, 0], selected, dlabel_path, vmax)
    fan = fig.add_subplot(grid[0, 1])
    _draw_curves(fan, curves, signed, norm, title=f"{dimension} parcels (n={stats['parcel'].nunique()})")
    colorbar = fig.colorbar(ScalarMappable(norm=norm, cmap=DIVERGING), ax=fan, fraction=0.046, pad=0.02)
    colorbar.ax.yaxis.set_tick_params(color="white", labelcolor="white")
    colorbar.outline.set_edgecolor("white")
    colorbar.set_label("partial R$^2$", color="white")
    up = fig.add_subplot(grid[1, 0])
    down = fig.add_subplot(grid[1, 1])
    _draw_curves(up, curves[curves["parcel"].isin(increasing)], signed, norm,
                 title=f"Monotone increasing trajectories (n={len(increasing)})",
                 color=INCREASING_COLOR)
    _draw_curves(down, curves[curves["parcel"].isin(decreasing)], signed, norm,
                 title=f"Monotone decreasing trajectories (n={len(decreasing)})",
                 color=DECREASING_COLOR)
    _overlay_brains(fig, up, {int(parcel): 1.0 for parcel in increasing}, dlabel_path,
                    INCREASING_CMAP, "lower right")
    _overlay_brains(fig, down, {int(parcel): 1.0 for parcel in decreasing}, dlabel_path,
                    DECREASING_CMAP, "upper right")
    fig.suptitle(
        "Centered f(age), anchored at the mean age. Color is the descriptive partial "
        "R$^2$ signed by the net age change.",
        color="white", fontsize=11, y=0.98,
    )
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=140, facecolor=fig.get_facecolor())
    plt.close(fig)


def plot_difference_matrix(dimension: str, stats: pd.DataFrame, curves: pd.DataFrame,
                           output_path: Path) -> None:
    """Grid of d(age) = f_row(age) - f_column(age)."""
    order = list(GROUP_ORDERS[dimension])
    subset = stats[stats["dimension"] == dimension].copy()
    curve_subset = curves[curves["dimension"] == dimension].copy()
    signed = {
        (row.group_a, row.group_b): float(row.signed_partial_r2)
        for row in subset.itertuples(index=False)
    }
    finite = np.array([value for value in signed.values() if np.isfinite(value)], dtype=float)
    vmax = float(np.max(np.abs(finite))) if finite.size else 1.0
    if vmax <= 0:
        vmax = 1.0
    norm = Normalize(vmin=-vmax, vmax=vmax)
    effects = []
    for rows in curve_subset.groupby(["group_a", "group_b"], sort=False):
        effects.append(rows[1]["centered_effect"].to_numpy(dtype=float))
    stacked = np.vstack(effects) if effects else np.zeros((1, 2))
    y_limit = float(np.nanmax(np.abs(stacked))) if np.isfinite(stacked).any() else 1.0
    if y_limit <= 0:
        y_limit = 1.0

    size = 2.05 * len(order) + 1.4
    fig, axes = plt.subplots(len(order), len(order), figsize=(size, size),
                             sharex=True, sharey=True, squeeze=False, facecolor="black")
    fig.subplots_adjust(left=0.12, right=0.90, top=0.90, bottom=0.08, wspace=0.15, hspace=0.15)
    lookup = {(a, b): rows for (a, b), rows in curve_subset.groupby(["group_a", "group_b"], sort=False)}
    for row_index, group_a in enumerate(order):
        for column_index, group_b in enumerate(order):
            ax = axes[row_index, column_index]
            ax.set_facecolor("black")
            if group_a == group_b:
                ax.text(0.5, 0.5, SHORT_GROUP_LABELS.get(group_a, group_a), color="white",
                        ha="center", va="center", fontsize=7, transform=ax.transAxes)
                ax.set_xticks([])
                ax.set_yticks([])
            else:
                rows = lookup[(group_a, group_b)]
                color = DIVERGING(norm(signed[(group_a, group_b)]))
                ax.plot(rows["age"], rows["centered_effect"], color=color, lw=1.4)
                ax.axhline(0, color="white", lw=0.4, alpha=0.6)
                ax.set_ylim(-y_limit, y_limit)
                _dark_axes(ax)
                ax.tick_params(labelsize=6)
            if row_index == len(order) - 1:
                ax.set_xlabel(SHORT_GROUP_LABELS.get(group_b, group_b), color="white", fontsize=7)
            if column_index == 0 and group_a != group_b:
                ax.set_ylabel(SHORT_GROUP_LABELS.get(group_a, group_a), color="white", fontsize=7)
            for spine in ax.spines.values():
                spine.set_color("#666666")
    colorbar = fig.colorbar(ScalarMappable(norm=norm, cmap=DIVERGING), ax=axes.ravel().tolist(),
                            fraction=0.02, pad=0.02)
    colorbar.ax.yaxis.set_tick_params(color="white", labelcolor="white")
    colorbar.set_label("signed partial R$^2$ of d(age)", color="white")
    fig.suptitle(
        f"{dimension}: d(age) = f_row(age) - f_column(age), centered at the mean age",
        color="white", fontsize=12,
    )
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=160, facecolor=fig.get_facecolor())
    plt.close(fig)


def write_square_matrices(stats: pd.DataFrame, output_dir: Path) -> list[Path]:
    """Write one labeled square file per comparison statistic.

    Rows are group A and columns are group B. The cell is f_A(age) - f_B(age).
    """
    columns = (
        "signed_partial_r2", "p_linear", "q_linear", "p_smooth", "q_smooth", "shape",
    )
    output_dir.mkdir(parents=True, exist_ok=True)
    written = []
    for dimension, order in GROUP_ORDERS.items():
        subset = stats[stats["dimension"] == dimension]
        for column in columns:
            grid = pd.DataFrame(np.nan, index=list(order), columns=list(order), dtype=object)
            for row in subset.itertuples(index=False):
                grid.loc[row.group_a, row.group_b] = getattr(row, column)
            path = output_dir / f"{dimension}_{column}_matrix.csv"
            grid.to_csv(path, index_label="group")
            written.append(path)
    return written


def write_linearity_table(result_dir: Path) -> Path:
    """One table of the h(age) test for parcels, Yeo networks, and the four systems."""
    parcels = pd.read_csv(result_dir / "parcel_gam_statistics.csv")
    groups = pd.read_csv(result_dir / "group_gam_statistics.csv")
    parcel_rows = pd.DataFrame({
        "level": "parcel",
        "name": parcels["parcel"].astype(int).astype(str),
        "gamma_h": parcels["gamma_h"],
        "p_smooth": parcels["p_smooth"],
        "q_smooth": parcels["q_smooth"],
        "p_linear": parcels["p_linear"],
        "q_linear": parcels["q_linear"],
        "shape": parcels["shape"],
        "signed_partial_r2": parcels["signed_partial_r2"],
    })
    group_rows = pd.DataFrame({
        "level": groups["dimension"],
        "name": groups["group"],
        "gamma_h": groups["gamma_h"],
        "p_smooth": groups["p_smooth"],
        "q_smooth": groups["q_smooth"],
        "p_linear": groups["p_linear"],
        "q_linear": groups["q_linear"],
        "shape": groups["shape"],
        "signed_partial_r2": groups["signed_partial_r2"],
    })
    table = pd.concat([parcel_rows, group_rows], ignore_index=True)
    table["departs_from_linear"] = table["q_smooth"] < 0.05
    path = result_dir / "linearity_deviation.csv"
    table.to_csv(path, index=False)
    return path


def plot_linearity(result_dir: Path, figure_dir: Path, dlabel_path: Path) -> list[Path]:
    """Brain map of parcel nonlinearity, and curves for networks and systems."""
    figure_dir.mkdir(parents=True, exist_ok=True)
    parcels = pd.read_csv(result_dir / "parcel_gam_statistics.csv")
    significant = parcels[parcels["q_smooth"] < 0.05]
    stat_by_parcel = {
        int(row.parcel): (1.0 if row.gamma_h > 0 else -1.0)
        for row in significant.itertuples(index=False)
        if np.isfinite(row.gamma_h) and row.gamma_h != 0
    }
    parcel_path = figure_dir / "linearity_deviation_parcels.png"
    fig = plt.figure(figsize=(8.2, 6.4), facecolor="black")
    _draw_brains(fig, fig.add_gridspec(1, 1)[0], stat_by_parcel, dlabel_path, 1.0,
                 "sign of gamma for h(age)")
    fig.suptitle(
        f"Parcels departing from linearity (q_smooth < 0.05, n={len(significant)}). "
        "Gold is a positive gamma, purple is negative.",
        color="white", fontsize=11,
    )
    fig.savefig(parcel_path, dpi=140, facecolor=fig.get_facecolor())
    plt.close(fig)

    groups = pd.read_csv(result_dir / "group_gam_statistics.csv")
    curves = pd.read_csv(result_dir / "group_curves.csv")
    group_path = figure_dir / "linearity_deviation_groups.png"
    fig, axes = plt.subplots(1, 2, figsize=(12.5, 5.2), facecolor="black")
    for ax, dimension, title in (
        (axes[0], "yeo", "Yeo 7 networks"),
        (axes[1], "system", "Four systems"),
    ):
        subset = groups[groups["dimension"] == dimension]
        q_by_group = dict(zip(subset["group"], subset["q_smooth"]))
        gamma_by_group = dict(zip(subset["group"], subset["gamma_h"]))
        for group, rows in curves[curves["dimension"] == dimension].groupby("group", sort=False):
            q_value = float(q_by_group[group])
            gamma = float(gamma_by_group[group])
            if np.isfinite(q_value) and q_value < 0.05:
                color = "#f6e27a" if gamma > 0 else "#8d6cab"
            else:
                color = "#8d8d8d"
            ax.plot(rows["age"], rows["centered_effect"], color=color, lw=1.8,
                    label=f"{group}  q={q_value:.3g}")
        ax.axhline(0, color="white", lw=0.6, alpha=0.7)
        ax.set_title(title, fontsize=11)
        ax.set_xlabel("Age")
        ax.set_ylabel("Centered f(age)")
        _dark_axes(ax)
        legend = ax.legend(frameon=False, fontsize=7, loc="best")
        for text in legend.get_texts():
            text.set_color("white")
    fig.suptitle(
        "Gold or purple: h(age) departs from linearity at q < 0.05. Gray: the smooth term does not.",
        color="white", fontsize=11,
    )
    fig.tight_layout()
    fig.savefig(group_path, dpi=140, facecolor=fig.get_facecolor())
    plt.close(fig)
    return [parcel_path, group_path]


def write_parcel_pair_matrices(stats: pd.DataFrame, output_dir: Path) -> list[Path]:
    """Square files for the 200 parcel-by-parcel trajectory comparisons."""
    subset = stats[stats["dimension"] == "parcel"]
    if subset.empty:
        raise ValueError("difference_statistics.csv has no parcel pairs")
    parcels = sorted({int(value) for value in subset["group_a"]} | {int(value) for value in subset["group_b"]})
    labels = [str(parcel) for parcel in parcels]
    written = []
    for column in ("signed_partial_r2", "p_linear", "q_linear", "p_smooth", "q_smooth", "net_change", "shape"):
        grid = pd.DataFrame(np.nan, index=labels, columns=labels, dtype=object)
        for row in subset.itertuples(index=False):
            grid.loc[str(int(row.group_a)), str(int(row.group_b))] = getattr(row, column)
        path = output_dir / f"parcel_{column}_matrix.csv"
        grid.to_csv(path, index_label="parcel")
        written.append(path)
    return written


def plot_parcel_pair_heatmaps(stats: pd.DataFrame, figure_dir: Path) -> list[Path]:
    """Heatmaps of d(age) for every ordered pair of Schaefer parcels."""
    subset = stats[stats["dimension"] == "parcel"].copy()
    parcels = sorted({int(value) for value in subset["group_a"]} | {int(value) for value in subset["group_b"]})
    index = {parcel: position for position, parcel in enumerate(parcels)}
    figure_dir.mkdir(parents=True, exist_ok=True)
    written = []
    panels = (
        ("net_change", "Net change of d(age), row minus column", True),
        ("q_smooth", "q_smooth for a nonlinear difference in d(age)", False),
    )
    for column, title, diverging in panels:
        grid = np.full((len(parcels), len(parcels)), np.nan)
        for row in subset.itertuples(index=False):
            grid[index[int(row.group_a)], index[int(row.group_b)]] = float(getattr(row, column))
        fig, ax = plt.subplots(figsize=(8.4, 7.2), facecolor="black")
        ax.set_facecolor("black")
        if diverging:
            finite = grid[np.isfinite(grid)]
            vmax = float(np.max(np.abs(finite))) if finite.size else 1.0
            if vmax <= 0:
                vmax = 1.0
            image = ax.imshow(grid, cmap=DIVERGING, vmin=-vmax, vmax=vmax, interpolation="nearest")
        else:
            image = ax.imshow(np.log10(np.clip(grid, 1e-300, 1)), cmap="magma_r",
                              vmin=-5, vmax=0, interpolation="nearest")
        ax.set_xlabel("Parcel B")
        ax.set_ylabel("Parcel A")
        ax.set_title(title, color="white", fontsize=11)
        ticks = np.linspace(0, len(parcels) - 1, 5, dtype=int)
        ax.set_xticks(ticks)
        ax.set_yticks(ticks)
        ax.set_xticklabels([str(parcels[tick]) for tick in ticks], color="white", fontsize=8)
        ax.set_yticklabels([str(parcels[tick]) for tick in ticks], color="white", fontsize=8)
        for spine in ax.spines.values():
            spine.set_color("white")
        colorbar = fig.colorbar(image, ax=ax, fraction=0.046, pad=0.04)
        colorbar.ax.yaxis.set_tick_params(color="white", labelcolor="white")
        colorbar.set_label("net change" if diverging else "log10 q_smooth", color="white")
        path = figure_dir / f"parcel_pair_{column}.png"
        fig.savefig(path, dpi=140, facecolor=fig.get_facecolor())
        plt.close(fig)
        written.append(path)
    return written


def plot_all(input_dir: Path, result_dir: Path, figure_dir: Path, dlabel_path: Path) -> list[Path]:
    parcel_info = pd.read_csv(input_dir / "parcel_info.csv")
    stats = pd.read_csv(result_dir / "parcel_gam_statistics.csv")
    curves = pd.read_csv(result_dir / "parcel_curves.csv")
    difference_stats = pd.read_csv(result_dir / "difference_statistics.csv")
    difference_curves = pd.read_csv(result_dir / "difference_curves.csv")
    written = []
    for dimension in ("yeo", "system"):
        fan_path = figure_dir / f"{dimension}_parcel_trajectories.png"
        matrix_path = figure_dir / f"{dimension}_difference_matrix.png"
        plot_parcel_suite(dimension, parcel_info, stats, curves, dlabel_path, fan_path)
        plot_difference_matrix(dimension, difference_stats, difference_curves, matrix_path)
        written.extend([fan_path, matrix_path])
    return written
