"""Age trajectories of Schaefer parcels with q_age < 0.05.

Reads the historical pygam tables written by run_pygam.py. Significance is
the parcel BH q already stored in pygam_age_comparison.csv. The curve comes
from pygam_curves.csv: covariates are held at their means, then the 1% to 99%
age-quantile curve is mean-centered. A positive net change means that centered
Reappraise minus Look effect is higher at the older end of this grid.

Run from the repository root:

    python -m analysis.fig03.gam_pygam.plot_trajectories --output "outputs/fig03/gam_pygam_volume"

The 200-parcel fsLR dlabel is the default. Parcel ids that are absent from it,
including Kong 2022 400-parcel ids, are drawn on
Schaefer2018_400Parcels_Kong2022_17Networks_order.dlabel.nii. That surface file
uses the same parcel ids and names as
Schaefer2018_400Parcels_Kong2022_17Networks_order_info.txt.
"""

from __future__ import annotations

import argparse
import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import nibabel as nib
import numpy as np
import pandas as pd
from matplotlib.colors import ListedColormap
from matplotlib.patches import Patch

from analysis.fig03.gam_supplement.plot_suites import (
    VIEWS,
    _hemisphere_labels,
    _paint_surface,
    _vertex_map,
)

PROJECT_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_OUTPUT = PROJECT_ROOT / "outputs" / "fig03" / "gam_pygam_volume"
DEFAULT_DLABEL = PROJECT_ROOT / "resources/atlases/schaefer2018/Schaefer2018_200Parcels_7Networks_order.dlabel.nii"
DEFAULT_DLABEL_400 = (
    PROJECT_ROOT / "resources/atlases/schaefer2018"
    / "Schaefer2018_400Parcels_Kong2022_17Networks_order.dlabel.nii"
)
DEFAULT_PARCEL_INFO_400 = (
    PROJECT_ROOT / "resources/atlases/schaefer2018"
    / "Schaefer2018_400Parcels_Kong2022_17Networks_order_info.txt"
)
LOCATION_COLORS = ("#E1A100", "#5C4B8A", "#1F8A70")
PARCEL_CURVE_COLORS = (
    "#4C9BE8", "#E07A3D", "#3CB371", "#E15B64", "#C57BDB", "#E6C15A",
    "#2EC4B6", "#F2A2C0", "#8FD14F", "#7AA2FF", "#FF8C6B", "#6ED0F6",
)


def resolve_project_path(path: Path) -> Path:
    if path.is_absolute():
        return path
    from_cwd = Path.cwd() / path
    if from_cwd.exists():
        return from_cwd
    from_root = PROJECT_ROOT / path
    if from_root.exists():
        return from_root
    return from_cwd


def _find_comparison(statistics_dir: Path) -> Path:
    preferred = statistics_dir / "pygam_age_comparison.csv"
    if preferred.is_file():
        return preferred
    matches = [
        path for path in statistics_dir.glob("*.csv")
        if "age" in path.stem.lower() and "compar" in path.stem.lower()
    ]
    if len(matches) == 1:
        return matches[0]
    raise FileNotFoundError(
        f"Could not find pygam_age_comparison.csv in {statistics_dir}"
    )


def _chinese_font() -> str | None:
    from matplotlib import font_manager
    for name in ("Microsoft YaHei", "Noto Sans CJK SC", "Source Han Sans SC", "WenQuanYi Micro Hei"):
        if any(found.name == name for found in font_manager.fontManager.ttflist):
            return name
    return None


def _net_change(curve: pd.DataFrame) -> float:
    ordered = curve.sort_values("age")
    effect = ordered["centered_effect"].to_numpy(dtype=float)
    finite = np.isfinite(effect)
    if finite.sum() < 2:
        return np.nan
    return float(effect[finite][-1] - effect[finite][0])


def significant_parcels(comparison: pd.DataFrame, curves: pd.DataFrame,
                        parcel_info: pd.DataFrame, alpha: float) -> pd.DataFrame:
    """Parcels whose BH q is below alpha, with the sign of the plotted curve."""
    parcels = comparison[comparison["dimension"].astype(str) == "parcel"].copy()
    parcels["parcel"] = parcels["name"].astype(int)
    kept = parcels[parcels["q_age"].astype(float) < alpha].copy()
    changes = []
    for parcel, rows in curves[curves["dimension"].astype(str) == "parcel"].groupby("name", sort=False):
        changes.append({"parcel": int(parcel), "net_change": _net_change(rows)})
    change_table = pd.DataFrame(changes)
    if change_table.empty:
        kept["net_change"] = np.nan
    else:
        kept = kept.merge(change_table, on="parcel", how="left")
    info = parcel_info.copy()
    info["parcel"] = info["parcel"].astype(int)
    info = info.rename(columns={"name": "parcel_name"})
    keep_columns = [column for column in ("parcel", "parcel_name", "hemi", "network") if column in info.columns]
    kept = kept.drop(columns=["name"], errors="ignore")
    kept = kept.merge(info.loc[:, keep_columns], on="parcel", how="left")
    kept["direction"] = np.where(kept["net_change"] >= 0, "higher_at_older_age", "higher_at_younger_age")
    kept.loc[~np.isfinite(kept["net_change"]), "direction"] = "unknown"
    return kept.sort_values(["direction", "q_age", "parcel"]).reset_index(drop=True)


def _style_axis(ax: plt.Axes) -> None:
    ax.set_facecolor("black")
    ax.tick_params(colors="white", labelsize=8)
    for spine in ax.spines.values():
        spine.set_color("white")
    ax.xaxis.label.set_color("white")
    ax.yaxis.label.set_color("white")
    ax.title.set_color("white")


def _parcel_colors(parcel_ids: list[int]) -> dict[int, str]:
    return {
        int(parcel): PARCEL_CURVE_COLORS[index % len(PARCEL_CURVE_COLORS)]
        for index, parcel in enumerate(parcel_ids)
    }


def _discrete_paint(colors: dict[int, str]) -> tuple[dict[int, float], ListedColormap, float, float]:
    parcels = list(colors)
    stat = {parcel: float(index + 1) for index, parcel in enumerate(parcels)}
    swatches = [colors[parcel] for parcel in parcels]
    if len(swatches) <= 1:
        swatch = swatches[0] if swatches else "#bbbbbb"
        return stat, ListedColormap([swatch, swatch]), 1.0, 2.0
    return stat, ListedColormap(swatches), 1.0, float(len(swatches))


def _draw_curves(ax: plt.Axes, curves: pd.DataFrame, selected: pd.DataFrame,
                 title: str, colors: dict[int, str]) -> None:
    parcel_curves = curves[curves["dimension"].astype(str) == "parcel"].copy()
    if not parcel_curves.empty:
        parcel_curves["parcel"] = parcel_curves["name"].astype(int)
    for parcel in selected["parcel"].astype(int):
        rows = parcel_curves[parcel_curves["parcel"] == int(parcel)] if not parcel_curves.empty else parcel_curves
        if rows.empty:
            continue
        ordered = rows.sort_values("age")
        ax.plot(
            ordered["age"], ordered["centered_effect"],
            color=colors[int(parcel)], lw=1.6, alpha=0.95,
        )
    ax.axhline(0, color="white", lw=0.6, alpha=0.7)
    ax.set_title(title, fontsize=11)
    ax.set_xlabel("Age")
    ax.set_ylabel("Centered age effect")
    _style_axis(ax)


def _surface_covers_parcels(dlabel_path: Path, parcels: list[int]) -> bool:
    """True when every requested parcel id is a cortical label on this dlabel."""
    if not parcels or not dlabel_path.is_file():
        return False
    present: set[int] = set()
    for hemi in ("lh", "rh"):
        labels = _hemisphere_labels(dlabel_path, hemi)
        present.update(int(value) for value in np.unique(labels) if int(value) > 0)
    return set(int(parcel) for parcel in parcels).issubset(present)


def _label_names(dlabel_path: Path) -> dict[int, str]:
    image = nib.load(str(dlabel_path))
    table = image.header.get_axis(0).label[0]
    names = {}
    for key, value in table.items():
        parcel = int(key)
        if parcel <= 0:
            continue
        names[parcel] = value[0] if isinstance(value, (tuple, list)) else str(value)
    return names


def _names_match(dlabel_path: Path, expected: dict[int, str]) -> bool:
    """Parcel ids 1–200 exist on both atlases, so the name has to agree."""
    if not expected or not dlabel_path.is_file():
        return False
    table = _label_names(dlabel_path)
    return all(table.get(int(parcel)) == name for parcel, name in expected.items())


def surface_dlabel_for(parcels: list[int], requested: Path,
                       names: dict[int, str] | None = None) -> Path:
    """Use the requested dlabel, or the Kong 2022 400-parcel surface when needed.

    The same integer is a different parcel on the 200-parcel dlabel. When names
    are available, the surface file is accepted only if those names match.
    """
    candidates = [requested]
    if DEFAULT_DLABEL_400.resolve() != requested.resolve():
        candidates.append(DEFAULT_DLABEL_400)
    if names:
        for candidate in candidates:
            if _names_match(candidate, names):
                return candidate
        return requested
    for candidate in candidates:
        if _surface_covers_parcels(candidate, parcels):
            return candidate
    return requested


def _draw_brains(fig: plt.Figure, spec, colors: dict[int, str], dlabel_path: Path) -> None:
    stat, cmap, vmin, vmax = _discrete_paint(colors)
    inner = spec.subgridspec(2, 2, wspace=0.02, hspace=0.02)
    for index, (hemi, view) in enumerate(VIEWS):
        ax = fig.add_subplot(inner[index // 2, index % 2], projection="3d")
        _paint_surface(ax, hemi, view, _vertex_map(stat, dlabel_path, hemi), cmap, vmin, vmax)


def plot_significant_trajectories(selected: pd.DataFrame, curves: pd.DataFrame,
                                  dlabel_path: Path, output_path: Path, alpha: float,
                                  curve_title: str | None = None,
                                  figure_title: str | None = None) -> None:
    font = _chinese_font()
    if font:
        plt.rcParams["font.sans-serif"] = [font]
        plt.rcParams["axes.unicode_minus"] = False
    selected = selected.sort_values("parcel").reset_index(drop=True)
    if font:
        if curve_title is None:
            curve_title = f"q < {alpha:g} 的脑区（n={len(selected)}）"
        if figure_title is None:
            figure_title = "Reappraise − Look 的中心化年龄曲线。零点是这条曲线自身的平均值。"
        axis_labels = ("年龄", "中心化年龄效应")
    else:
        if curve_title is None:
            curve_title = f"Parcels with q < {alpha:g} (n={len(selected)})"
        if figure_title is None:
            figure_title = "Centered Reappraise − Look curve. Zero is the mean of that curve."
        axis_labels = ("Age", "Centered age effect")
    legend_rows = max(len(selected), 1)
    legend_columns = 1 if legend_rows <= 8 else 2 if legend_rows <= 24 else 3
    shown_rows = math.ceil(legend_rows / legend_columns)
    fig = plt.figure(figsize=(12.4, 5.6 + max(0, shown_rows - 3) * 0.28), facecolor="black")
    parcel_ids = selected["parcel"].astype(int).tolist()
    colors = _parcel_colors(parcel_ids)
    names = None
    if "parcel_name" in selected.columns and parcel_ids:
        named = {
            int(row.parcel): str(row.parcel_name)
            for row in selected.itertuples(index=False)
            if isinstance(row.parcel_name, str) and row.parcel_name.strip()
        }
        if len(named) == len(parcel_ids):
            names = named
    dlabel_path = surface_dlabel_for(parcel_ids, dlabel_path, names)
    names_agree = names is None or _names_match(dlabel_path, names)
    draw_brain = names_agree and _surface_covers_parcels(dlabel_path, parcel_ids)
    if draw_brain:
        print(f"Surface labels: {dlabel_path}")
    bottom = min(0.08 + 0.034 * shown_rows, 0.62)
    if draw_brain:
        grid = fig.add_gridspec(1, 2, width_ratios=[1.0, 1.25], wspace=0.20,
                                left=0.06, right=0.98, top=0.84, bottom=bottom)
        _draw_brains(fig, grid[0, 0], colors, dlabel_path)
        curve_ax = fig.add_subplot(grid[0, 1])
    else:
        if parcel_ids:
            print("Surface labels do not cover these parcels. The figure shows trajectories only.")
        grid = fig.add_gridspec(1, 1, left=0.08, right=0.98, top=0.84, bottom=bottom)
        curve_ax = fig.add_subplot(grid[0, 0])
    _draw_curves(curve_ax, curves, selected, curve_title, colors)
    curve_ax.set_xlabel(axis_labels[0])
    curve_ax.set_ylabel(axis_labels[1])
    handles = []
    name_by_parcel = {}
    if "parcel_name" in selected.columns:
        name_by_parcel = {
            int(row.parcel): str(row.parcel_name)
            for row in selected.itertuples(index=False)
            if isinstance(row.parcel_name, str) and row.parcel_name.strip()
        }
    for parcel in parcel_ids:
        label = f"{parcel}  {name_by_parcel[parcel]}" if parcel in name_by_parcel else str(parcel)
        handles.append(Patch(facecolor=colors[parcel], edgecolor="none", label=label))
    if handles:
        fig.legend(
            handles=handles, loc="lower center", frameon=False, fontsize=8,
            labelcolor="white", ncol=legend_columns, bbox_to_anchor=(0.5, 0.01),
        )
    fig.suptitle(figure_title, color="white", fontsize=12)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    if output_path.suffix.lower() == ".svg":
        plt.rcParams["svg.fonttype"] = "none"
    fig.savefig(output_path, dpi=160, facecolor=fig.get_facecolor())
    plt.close(fig)


def plot_group_age_curves(statistics: pd.DataFrame, curves: pd.DataFrame, figure_dir: Path) -> list[Path]:
    """One SVG panel per Yeo network, system, and amygdala ROI that was fit."""
    figure_dir.mkdir(parents=True, exist_ok=True)
    written = []
    for dimension in ("yeo", "system", "amygdala"):
        stats = statistics[statistics["dimension"].astype(str) == dimension].copy()
        subset = curves[curves["dimension"].astype(str) == dimension].copy()
        if stats.empty or subset.empty:
            continue
        names = stats["name"].astype(str).tolist()
        columns = 3 if dimension == "amygdala" else 4
        rows = math.ceil(len(names) / columns)
        figure, axes = plt.subplots(rows, columns, figsize=(columns * 3.4, rows * 2.6), squeeze=False)
        flat_axes = axes.ravel()
        for index, name in enumerate(names):
            axis = flat_axes[index]
            piece = subset[subset["name"].astype(str) == name].sort_values("age")
            axis.plot(piece["age"], piece["centered_effect"], color="#b35806", linewidth=1.6)
            row = stats[stats["name"].astype(str) == name].iloc[0]
            q_value = float(row["q_age"]) if "q_age" in stats.columns else np.nan
            axis.set_title(f"{name}\nq_age={q_value:.3g}", fontsize=9)
            axis.axhline(0, color="0.7", linewidth=0.6)
            axis.set_xlabel("Age")
            axis.set_ylabel("Centered age effect")
        for axis in flat_axes[len(names):]:
            axis.axis("off")
        figure.tight_layout()
        path = figure_dir / f"{dimension}_age_curves.svg"
        plt.rcParams["svg.fonttype"] = "none"
        figure.savefig(path)
        plt.close(figure)
        written.append(path)
    return written


def plot_labeled_parcels(parcel_ids: list[int], dlabel_path: Path, output_path: Path,
                         info_path: Path = DEFAULT_PARCEL_INFO_400) -> None:
    """Paint Kong 2022 parcels on the fsLR inflated cortex.

    Colors are distinct so neighboring parcels stay separable. The info table
    supplies the legend names; vertex placement comes only from the dlabel.
    """
    from analysis.fig03.gam_supplement.prepare_inputs import read_schaefer_order_info

    ids = [int(parcel) for parcel in parcel_ids]
    if not ids:
        raise ValueError("No parcels to draw")
    info = read_schaefer_order_info(info_path).set_index("parcel")
    missing = [parcel for parcel in ids if parcel not in info.index]
    if missing:
        raise ValueError(f"Parcel ids are absent from {info_path}: {missing}")
    expected = {parcel: str(info.loc[parcel, "name"]) for parcel in ids}
    dlabel_path = surface_dlabel_for(ids, dlabel_path, expected)
    if not _names_match(dlabel_path, expected):
        raise ValueError(f"Surface label names do not match {info_path.name}: {dlabel_path}")
    font = _chinese_font()
    if font:
        plt.rcParams["font.sans-serif"] = [font]
        plt.rcParams["axes.unicode_minus"] = False
    colors = [LOCATION_COLORS[index % len(LOCATION_COLORS)] for index in range(len(ids))]
    stat = {parcel: float(index + 1) for index, parcel in enumerate(ids)}
    cmap = ListedColormap(colors)
    fig = plt.figure(figsize=(8.6, 5.4), facecolor="black")
    grid = fig.add_gridspec(2, 2, wspace=0.02, hspace=0.02, left=0.02, right=0.98, top=0.90, bottom=0.14)
    for index, (hemi, view) in enumerate(VIEWS):
        ax = fig.add_subplot(grid[index // 2, index % 2], projection="3d")
        _paint_surface(ax, hemi, view, _vertex_map(stat, dlabel_path, hemi), cmap, 1.0, float(len(ids)))
    handles = []
    for parcel, color in zip(ids, colors):
        row = info.loc[parcel]
        handles.append(Patch(
            facecolor=color, edgecolor="none",
            label=f"{parcel}  {row['name']}",
        ))
    fig.legend(
        handles=handles, loc="lower center", ncol=1, frameon=False,
        fontsize=8, labelcolor="white", bbox_to_anchor=(0.5, 0.01),
    )
    fig.suptitle("Schaefer 400，Kong 2022", color="white", fontsize=13)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=180, facecolor=fig.get_facecolor())
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description="Plot significant pygam parcel age trajectories")
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--alpha", type=float, default=0.05)
    parser.add_argument("--dlabel", type=Path, default=DEFAULT_DLABEL)
    args = parser.parse_args()
    output_dir = resolve_project_path(args.output)
    statistics_dir = output_dir / "statistics"
    comparison = pd.read_csv(_find_comparison(statistics_dir))
    curves = pd.read_csv(statistics_dir / "pygam_curves.csv")
    parcel_info = pd.read_csv(output_dir / "inputs" / "parcel_info.csv")
    selected = significant_parcels(comparison, curves, parcel_info, args.alpha)
    list_path = statistics_dir / "significant_parcels.csv"
    columns = [column for column in (
        "parcel", "parcel_name", "hemi", "network", "q_age", "partial_r2", "net_change", "direction",
    ) if column in selected.columns]
    selected.loc[:, columns].to_csv(list_path, index=False)
    figure_path = output_dir / "figures" / "significant_parcel_trajectories.svg"
    plot_significant_trajectories(
        selected, curves, resolve_project_path(args.dlabel), figure_path, args.alpha)
    print(f"Significant parcels: {len(selected)}")
    print(f"Figure: {figure_path}")
    print(f"List: {list_path}")


if __name__ == "__main__":
    main()
