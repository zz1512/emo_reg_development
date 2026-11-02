"""Curves of the linear model and the smooth model on one age axis."""

from __future__ import annotations

import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.cm import ScalarMappable
from matplotlib.colors import LinearSegmentedColormap, Normalize

from analysis.fig03.gam_supplement.prepare_inputs import SYSTEM_ORDER, YEO_ORDER

DIVERGING = LinearSegmentedColormap.from_list(
    "m3_slope",
    ["#3b0a74", "#8d6cab", "#f4f1ea", "#e6c36a", "#f6e27a"],
)
SHORT_GROUP_LABELS = {
    "reappraisal_only": "reappraisal\nonly",
    "common_appraisal": "common\nappraisal",
    "nonmodifiable_emotion": "nonmodifiable\nemotion",
    "modifiable_emotion": "modifiable\nemotion",
}
PAIR_ORDERS = {"yeo": YEO_ORDER, "system": SYSTEM_ORDER}

DIMENSIONS = ("yeo", "system", "amygdala")


def _save_figure(figure: plt.Figure, path: Path, **kwargs) -> None:
    """Write SVG with text left as text, so labels stay editable."""
    if path.suffix.lower() == ".svg":
        plt.rcParams["svg.fonttype"] = "none"
    path.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(path, **kwargs)


def plot_curve_panels(statistics: pd.DataFrame, curves: pd.DataFrame, figure_dir: Path) -> list[Path]:
    """One panel per network, system, and amygdala ROI. Parcels stay in the table."""
    figure_dir.mkdir(parents=True, exist_ok=True)
    written = []
    for dimension in DIMENSIONS:
        stats = statistics[statistics["dimension"] == dimension].copy()
        subset = curves[curves["dimension"] == dimension].copy()
        if stats.empty or subset.empty:
            raise ValueError(f"No {dimension} rows to plot")
        names = stats["name"].astype(str).tolist()
        columns = 3 if dimension == "amygdala" else 4
        rows = math.ceil(len(names) / columns)
        figure, axes = plt.subplots(rows, columns, figsize=(columns * 3.4, rows * 2.6), squeeze=False)
        flat_axes = axes.ravel()
        for index, name in enumerate(names):
            axis = flat_axes[index]
            piece = subset[subset["name"].astype(str) == name]
            for model, color, style in (("m2", "#4d4d4d", "--"), ("m3", "#b35806", "-")):
                line = piece[piece["model"] == model].sort_values("age")
                axis.plot(line["age"], line["centered_effect"], color=color, linestyle=style,
                          linewidth=1.6, label=model)
            row = stats[stats["name"].astype(str) == name].iloc[0]
            axis.set_title(f"{name}\nq(m3 vs m2)={row['q_m3_vs_m2']:.3g}", fontsize=9)
            axis.axhline(0, color="0.7", linewidth=0.6)
            axis.set_xlabel("Age")
            axis.set_ylabel("Centered age effect")
            if index == 0:
                axis.legend(frameon=False, fontsize=8)
        for axis in flat_axes[len(names):]:
            axis.axis("off")
        figure.tight_layout()
        path = figure_dir / f"{dimension}_m2_m3_curves.svg"
        _save_figure(figure, path)
        plt.close(figure)
        written.append(path)
    return written


def _dark_axes(axis: plt.Axes) -> None:
    axis.set_facecolor("black")
    axis.tick_params(colors="white", labelsize=6)
    for spine in axis.spines.values():
        spine.set_color("#666666")
    axis.xaxis.label.set_color("white")
    axis.yaxis.label.set_color("white")


def write_slope_matrices(stats: pd.DataFrame, output_dir: Path) -> list[Path]:
    """Square files for the m3 slope of group A minus group B. The diagonal is empty."""
    output_dir.mkdir(parents=True, exist_ok=True)
    written = []
    for dimension, order in PAIR_ORDERS.items():
        subset = stats[stats["dimension"] == dimension]
        for column in ("slope", "p_slope", "q_slope"):
            grid = pd.DataFrame(np.nan, index=list(order), columns=list(order), dtype=float)
            for row in subset.itertuples(index=False):
                grid.loc[row.group_a, row.group_b] = getattr(row, column)
            path = output_dir / f"{dimension}_m3_{column}_matrix.csv"
            grid.to_csv(path, index_label="group")
            written.append(path)
    return written


def plot_m3_slope_matrices(stats: pd.DataFrame, curves: pd.DataFrame, figure_dir: Path) -> list[Path]:
    """One grid per dimension. Each panel is the m3 age curve of row minus column."""
    figure_dir.mkdir(parents=True, exist_ok=True)
    written = []
    for dimension, order in PAIR_ORDERS.items():
        subset = stats[stats["dimension"] == dimension].copy()
        curve_subset = curves[curves["dimension"] == dimension].copy()
        if subset.empty or curve_subset.empty:
            raise ValueError(f"No {dimension} pairwise m3 slopes to plot")
        slopes = {
            (row.group_a, row.group_b): float(row.slope)
            for row in subset.itertuples(index=False)
        }
        finite = np.array(list(slopes.values()), dtype=float)
        finite = finite[np.isfinite(finite)]
        vmax = float(np.max(np.abs(finite))) if finite.size else 1.0
        if vmax <= 0:
            vmax = 1.0
        norm = Normalize(vmin=-vmax, vmax=vmax)
        stacked = curve_subset["centered_effect"].to_numpy(dtype=float)
        y_limit = float(np.nanmax(np.abs(stacked))) if np.isfinite(stacked).any() else 1.0
        if y_limit <= 0:
            y_limit = 1.0
        lookup = {
            (group_a, group_b): rows
            for (group_a, group_b), rows in curve_subset.groupby(["group_a", "group_b"], sort=False)
        }
        size = 2.05 * len(order) + 1.4
        figure, axes = plt.subplots(
            len(order), len(order), figsize=(size, size),
            sharex=True, sharey=True, squeeze=False, facecolor="black",
        )
        figure.subplots_adjust(left=0.12, right=0.90, top=0.90, bottom=0.08, wspace=0.15, hspace=0.15)
        order = list(order)
        for row_index, group_a in enumerate(order):
            for column_index, group_b in enumerate(order):
                axis = axes[row_index, column_index]
                axis.set_facecolor("black")
                if group_a == group_b:
                    axis.text(
                        0.5, 0.5, SHORT_GROUP_LABELS.get(group_a, group_a),
                        color="white", ha="center", va="center", fontsize=7, transform=axis.transAxes,
                    )
                    axis.set_xticks([])
                    axis.set_yticks([])
                else:
                    rows = lookup[(group_a, group_b)]
                    axis.plot(
                        rows["age"], rows["centered_effect"],
                        color=DIVERGING(norm(slopes[(group_a, group_b)])), lw=1.4,
                    )
                    axis.axhline(0, color="white", lw=0.4, alpha=0.6)
                    axis.set_ylim(-y_limit, y_limit)
                    _dark_axes(axis)
                if row_index == len(order) - 1:
                    axis.set_xlabel(SHORT_GROUP_LABELS.get(group_b, group_b), color="white", fontsize=7)
                if column_index == 0 and group_a != group_b:
                    axis.set_ylabel(SHORT_GROUP_LABELS.get(group_a, group_a), color="white", fontsize=7)
                for spine in axis.spines.values():
                    spine.set_color("#666666")
        colorbar = figure.colorbar(
            ScalarMappable(norm=norm, cmap=DIVERGING), ax=axes.ravel().tolist(),
            fraction=0.02, pad=0.02,
        )
        colorbar.ax.yaxis.set_tick_params(color="white", labelcolor="white")
        colorbar.set_label("m3 slope of d(age), per year", color="white")
        figure.suptitle(
            f"{dimension}: m3 d(age) = f_row(age) - f_column(age), centered at the mean age",
            color="white", fontsize=12,
        )
        path = figure_dir / f"{dimension}_m3_slope_matrix.svg"
        _save_figure(figure, path, facecolor=figure.get_facecolor())
        plt.close(figure)
        written.append(path)
    return written


def significant_m3_vs_m1_parcels(statistics: pd.DataFrame, parcel_info: pd.DataFrame,
                                 alpha: float = 0.05) -> pd.DataFrame:
    """Parcels whose m3 versus m1 comparison survives BH within parcels."""
    stats = statistics[statistics["dimension"].astype(str) == "parcel"].copy()
    if "q_m3_vs_m1" not in stats.columns:
        raise ValueError("nested_model_comparison.csv has no q_m3_vs_m1 column")
    stats["parcel"] = stats["name"].astype(int)
    selected = stats[stats["q_m3_vs_m1"].astype(float) < alpha].copy()
    info = parcel_info.loc[:, ["parcel", "name", "hemi", "network"]].copy()
    info["parcel"] = info["parcel"].astype(int)
    info = info.rename(columns={"name": "parcel_name"})
    selected = selected.merge(info, on="parcel", how="left")
    return selected.sort_values(["q_m3_vs_m1", "parcel"]).reset_index(drop=True)


def plot_m3_m1_parcel_trajectories(statistics: pd.DataFrame, curves: pd.DataFrame,
                                   parcel_info: pd.DataFrame, dlabel_path: Path,
                                   figure_path: Path, alpha: float = 0.05) -> pd.DataFrame:
    """Draw m3 curves for parcels with q(m3 vs m1) < alpha, one color per parcel."""
    from analysis.fig03.gam_pygam.plot_trajectories import plot_significant_trajectories

    selected = significant_m3_vs_m1_parcels(statistics, parcel_info, alpha)
    if selected.empty:
        return selected
    parcel_curves = curves[
        (curves["dimension"].astype(str) == "parcel") & (curves["model"].astype(str) == "m3")
    ].copy()
    if parcel_curves.empty:
        raise ValueError("nested_curves.csv has no parcel m3 curves")
    font_name = None
    from matplotlib import font_manager
    for name in ("Microsoft YaHei", "Noto Sans CJK SC", "Source Han Sans SC", "WenQuanYi Micro Hei"):
        if any(found.name == name for found in font_manager.fontManager.ttflist):
            font_name = name
            break
    if font_name:
        curve_title = f"q(m3 vs m1) < {alpha:g} 的脑区（n={len(selected)}）"
        figure_title = "m3 的中心化年龄曲线。协变量取均值，零点是平均年龄处的预测值。"
    else:
        curve_title = f"Parcels with q(m3 vs m1) < {alpha:g} (n={len(selected)})"
        figure_title = "Centered m3 age curve. Covariates are at their means. Zero is the prediction at the mean age."
    plot_significant_trajectories(
        selected, parcel_curves, dlabel_path, figure_path, alpha,
        curve_title=curve_title, figure_title=figure_title,
    )
    return selected


def plot_age_parent_edu(design: pd.DataFrame, figure_path: Path) -> pd.DataFrame:
    """Scatter of age against parent_edu, with Pearson and Spearman summaries."""
    from analysis.fig03.gam_pygam.age_parent_edu import age_parent_edu_association

    summary = age_parent_edu_association(design)
    font_name = None
    from matplotlib import font_manager
    for name in ("Microsoft YaHei", "Noto Sans CJK SC", "Source Han Sans SC", "WenQuanYi Micro Hei"):
        if any(found.name == name for found in font_manager.fontManager.ttflist):
            font_name = name
            break
    if font_name:
        plt.rcParams["font.sans-serif"] = [font_name]
        plt.rcParams["axes.unicode_minus"] = False
        xlabel, ylabel = "父母受教育水平", "年龄"
    else:
        xlabel, ylabel = "Parental education", "Age"
    age = design["age"].to_numpy(dtype=float)
    education = design["parent_edu"].to_numpy(dtype=float)
    row = summary.iloc[0]
    slope, intercept = np.polyfit(education, age, 1)
    x_line = np.linspace(float(np.min(education)), float(np.max(education)), 50)
    figure, axis = plt.subplots(figsize=(6.2, 4.6))
    axis.scatter(education, age, s=12, c="#3d6f99", alpha=0.45, linewidths=0)
    axis.plot(x_line, intercept + slope * x_line, color="#c45c26", lw=1.6)
    axis.set_xlabel(xlabel)
    axis.set_ylabel(ylabel)
    axis.set_title(
        f"n={int(row['n'])}   Pearson r={row['pearson_r']:.3f} (p={row['pearson_p']:.3g})"
        f"\nSpearman r={row['spearman_r']:.3f} (p={row['spearman_p']:.3g})",
        fontsize=11,
    )
    figure.tight_layout()
    _save_figure(figure, figure_path)
    plt.close(figure)
    return summary
