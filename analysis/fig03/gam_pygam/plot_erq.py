"""Conditional predictions and ERQ slopes from one continuous age-by-ERQ model.

The caller supplies units passing either pre-specified sequential FDR test. The three ERQ
reference curves are predictions from the same fit, not separately fitted
groups. Their common centering reference preserves the ERQ main effect.
"""

from __future__ import annotations

from pathlib import Path
import re

import matplotlib

matplotlib.use("Agg")
matplotlib.rcParams["svg.fonttype"] = "none"
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from analysis.fig03.gam_pygam.plot_region_age import AGE_TICKS, BOX_ASPECT, zoom_y_limits


LEVELS = (
    ("low", "Low (P16)", "#00A8B5"),
    ("median", "Median (P50)", "#E5A000"),
    ("high", "High (P84)", "#E76F51"),
)
MANIFEST_COLUMNS = (
    "dimension", "name", "label", "q_full_vs_baseline", "q_joint_erq", "q_interaction",
    "q_additive_vs_baseline", "q_full_vs_additive", "age_figure_file", "slope_figure_file",
)


def _with_keys(table: pd.DataFrame, required: set[str], label: str) -> pd.DataFrame:
    missing = sorted(required - set(table.columns))
    if missing:
        raise ValueError(f"{label} is missing columns: {', '.join(missing)}")
    result = table.copy()
    if result[["dimension", "name"]].isna().any().any():
        raise ValueError(f"{label} has missing dimension/name keys")
    result["dimension"] = result["dimension"].astype(str)
    result["name"] = result["name"].astype(str)
    parcel = result["dimension"] == "parcel"
    if parcel.any():
        ids = pd.to_numeric(result.loc[parcel, "name"], errors="raise")
        if not np.isfinite(ids).all() or not (ids == np.floor(ids)).all():
            raise ValueError(f"{label} has invalid parcel identifiers")
        result.loc[parcel, "name"] = ids.astype(int).astype(str)
    return result


def _interval_arrays(table: pd.DataFrame, value: str, label: str):
    age, effect, low, high = (table[column].to_numpy(dtype=float)
                              for column in ("age", value, "ci_low", "ci_high"))
    if not all(np.isfinite(array).all() for array in (age, effect, low, high)):
        raise ValueError(f"{label} contains non-finite predictions or CI values")
    if np.any(low > high) or np.any(effect < low) or np.any(effect > high):
        raise ValueError(f"{label} has an invalid confidence interval")
    if table["age"].duplicated().any():
        raise ValueError(f"{label} has duplicate ages")
    return age, effect, low, high


def _axes(axis, age: np.ndarray, effect: np.ndarray,
          low: np.ndarray, high: np.ndarray, ylabel: str) -> None:
    axis.set_xticks(AGE_TICKS, labels=[str(tick) for tick in AGE_TICKS])
    xmin = min(float(age.min()), AGE_TICKS[0])
    xmax = max(float(age.max()), AGE_TICKS[-1])
    margin = 0.05 * (xmax - xmin)
    axis.set_xlim(xmin - margin, xmax + margin)
    axis.set_ylim(*zoom_y_limits(effect, low, high))
    axis.set_box_aspect(BOX_ASPECT)
    axis.set_xlabel("Age")
    axis.set_ylabel(ylabel)
    axis.spines["top"].set_visible(False)
    axis.spines["right"].set_visible(False)


def _filename(dimension: str, name: str) -> str:
    slug = f"{int(name):04d}" if dimension == "parcel" else re.sub(
        r"[^A-Za-z0-9_.-]+", "_", name).strip("._")
    prefix = re.sub(r"[^A-Za-z0-9_.-]+", "_", dimension).strip("._")
    if not slug or not prefix:
        raise ValueError(f"Cannot construct an ERQ figure filename for {dimension!r}, {name!r}")
    return f"{prefix}_{slug}.svg"


def write_erq_figures(curves: pd.DataFrame, slopes: pd.DataFrame,
                      statistics: pd.DataFrame, figure_dir: Path) -> pd.DataFrame:
    """Write two SVGs per sequentially selected unit and return file paths.

    CI bounds are drawn as supplied; neither selection nor additional FDR
    correction is performed here. ``q_interaction`` is the value already
    corrected over the full analysis family before selection.
    An empty selection returns an empty manifest without creating figures.
    """
    if statistics.empty:
        return pd.DataFrame(columns=MANIFEST_COLUMNS)
    statistics = _with_keys(statistics, {
        "dimension", "name", "label", "q_full_vs_baseline", "q_joint_erq", "q_interaction",
        "q_additive_vs_baseline", "q_full_vs_additive",
    }, "ERQ statistics")
    curves = _with_keys(curves, {
        "dimension", "name", "age", "erq_level", "erq_value", "effect", "ci_low", "ci_high",
    }, "ERQ conditional curves")
    slopes = _with_keys(slopes, {
        "dimension", "name", "age", "erq_slope", "ci_low", "ci_high",
    }, "ERQ slopes")
    if statistics.duplicated(["dimension", "name"]).any():
        raise ValueError("ERQ statistics has duplicate dimension/name rows")
    q_values = statistics[["q_additive_vs_baseline", "q_full_vs_additive"]].to_numpy(dtype=float)
    valid = np.isfinite(q_values) & (q_values >= 0) & (q_values < .05)
    if np.any(~valid.any(axis=1)):
        raise ValueError("ERQ figures require at least one sequential q < 0.05 per unit")
    figure_dir = Path(figure_dir)
    filenames: set[str] = set()
    manifest: list[dict] = []
    for row in statistics.itertuples(index=False):
        key = f"{row.dimension} {row.name}"
        conditional = curves[(curves["dimension"] == row.dimension) & (curves["name"] == row.name)]
        slope = slopes[(slopes["dimension"] == row.dimension) & (slopes["name"] == row.name)].sort_values("age")
        if set(conditional["erq_level"]) != {item[0] for item in LEVELS} or slope.empty:
            raise ValueError(f"Incomplete ERQ prediction grids for {key}")
        level_arrays = []
        for level, description, color in LEVELS:
            values = conditional[conditional["erq_level"] == level].sort_values("age")
            erq_values = values["erq_value"].to_numpy(dtype=float)
            if not np.isfinite(erq_values).all() or np.unique(erq_values).size != 1:
                raise ValueError(f"{key} {level} must use one finite ERQ reference value")
            arrays = _interval_arrays(values, "effect", f"{key} {level}")
            level_arrays.append((description, color, erq_values[0], arrays))
        slope_arrays = _interval_arrays(slope, "erq_slope", key)
        filename = _filename(row.dimension, row.name)
        if filename in filenames:
            raise ValueError(f"ERQ units produce duplicate figure filename: {filename}")
        filenames.add(filename)
        age_path = Path("age_by_erq") / filename
        slope_path = Path("erq_slope_by_age") / filename
        for path in (age_path, slope_path):
            (figure_dir / path).parent.mkdir(parents=True, exist_ok=True)

        figure, axis = plt.subplots(figsize=(5.0, 6.8))
        try:
            for description, color, erq_value, arrays in level_arrays:
                age, effect, low, high = arrays
                axis.fill_between(age, low, high, color=color, alpha=0.18, linewidth=0, zorder=1)
                axis.plot(age, effect, color=color, lw=1.8, zorder=2,
                          label=f"{description}: ERQ = {erq_value:.3g}")
            envelope = [np.concatenate([item[3][index] for item in level_arrays]) for index in range(4)]
            _axes(axis, *envelope, ylabel="Centered model prediction")
            axis.set_title(f"{row.label}\nq(ERQ)={row.q_additive_vs_baseline:.3g}; "
                           f"q(interaction)={row.q_full_vs_additive:.3g}\n95% CI (with mean)", fontsize=9)
            figure.legend(*axis.get_legend_handles_labels(), loc="lower center", frameon=False,
                          title="ERQ reference values (one continuous model)",
                          fontsize=8, title_fontsize=8, bbox_to_anchor=(0.5, 0.015))
            figure.tight_layout(rect=(0, 0.18, 1, 1))
            figure.savefig(figure_dir / age_path, format="svg")
        finally:
            plt.close(figure)

        figure, axis = plt.subplots(figsize=(5.0, 6.2))
        try:
            age, effect, low, high = slope_arrays
            axis.fill_between(age, low, high, color="#FFD45C", alpha=0.4, linewidth=0, zorder=1)
            axis.plot(age, effect, color="#E5A000", lw=1.8, zorder=2)
            axis.axhline(0, color="#777777", lw=0.7, linestyle="--", zorder=0)
            _axes(axis, *slope_arrays, ylabel="ERQ slope (activation per score point)")
            q_interaction = f"{row.q_interaction:.3g}" if np.isfinite(row.q_interaction) else "NA"
            axis.set_title(f"{row.label}\nq(ERQ)={row.q_additive_vs_baseline:.3g}\nq(interaction)={q_interaction}"
                           "\n95% pointwise CI", fontsize=9)
            figure.text(0.5, 0.015, "Smooth age x ERQ interaction; q corrected over full family",
                        ha="center", fontsize=7)
            figure.tight_layout(rect=(0, 0.035, 1, 1))
            figure.savefig(figure_dir / slope_path, format="svg")
        finally:
            plt.close(figure)
        manifest.append({
            "dimension": row.dimension, "name": row.name, "label": row.label,
            "q_full_vs_baseline": row.q_full_vs_baseline,
            "q_joint_erq": row.q_joint_erq, "q_interaction": row.q_interaction,
            "q_additive_vs_baseline": row.q_additive_vs_baseline, "q_full_vs_additive": row.q_full_vs_additive,
            "age_figure_file": age_path.as_posix(), "slope_figure_file": slope_path.as_posix(),
        })
    return pd.DataFrame(manifest, columns=MANIFEST_COLUMNS)
