"""Overview and zoom SVGs for each parcel, network, system, or Juelich ROI.

The points are m3 residuals plus its centered age effect. The line and its 95%
pointwise CI describe that same effect. Significant parcels and parcels 156/157 in
the 432 atlas use -1.5..1.5; all other units use -1..1. The plotting box has
a 3:4 width-to-height ratio. Points outside the unit's range remain in the
saved partial-residual table and are counted in the trajectory list. A separate
zoom shows only the same curve and CI, with limits based on their actual range.
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

CURVE_COLOR = "#E5A000"
BAND_COLOR = "#FFD45C"
SCATTER_COLOR = "#666666"
Y_LIMITS = (-1.0, 1.0)
Y_LIMITS_432_SELECTED = (-1.5, 1.5)
BOX_ASPECT = 4.0 / 3.0
DIMENSIONS = ("parcel", "yeo", "system", "amygdala")
TRAJECTORY_COLUMNS = (
    "dimension", "name", "parcel", "parcel_name", "q_age",
    "q_m3_vs_m1", "q_m3_vs_m2", "partial_r2", "partial_r2_m3_vs_m1",
    "partial_r2_m3_vs_m2", "reason",
)


def curve_style_metadata() -> dict:
    return {
        "y_limits": {
            "default": list(Y_LIMITS),
            "parcel432_significant_or_156_157": list(Y_LIMITS_432_SELECTED),
        },
        "significant": "q_m3_vs_m1 < 0.05 or q_m3_vs_m2 < 0.05",
        "axes_width_to_height": "3:4",
        "color": "bright gold",
        "views": {
            "overview": "age_curves/<unit>.svg, m3 partial residuals and centered age curve",
            "zoom": "age_curves_zoom/<unit>_zoom.svg, same curve and CI without scatter",
        },
        "zoom_y_limits": "effect and CI envelope plus 15% of its span at each end",
    }


def _number(row: pd.Series, name: str, fallback: str | None = None) -> float:
    value = row.get(name, row.get(fallback, np.nan) if fallback else np.nan)
    return float(value) if pd.notna(value) else np.nan


def _unit_name(dimension: str, name) -> str:
    return str(int(name)) if dimension == "parcel" else str(name)


def trajectory_rows(statistics: pd.DataFrame, parcel_info: pd.DataFrame,
                    marginal: list[int], amygdala: list[int]) -> pd.DataFrame:
    """Selected parcels and every fitted network, system, and Juelich ROI."""
    info = parcel_info.copy()
    info["parcel"] = info["parcel"].astype(int)
    if info["parcel"].duplicated().any():
        raise ValueError("parcel_info.csv has duplicate parcel ids")
    names = info.set_index("parcel")["name"].to_dict() if "name" in info else {}
    marginal_ids = {int(parcel) for parcel in marginal}
    amygdala_ids = {int(parcel) for parcel in amygdala}
    rows = []
    for dimension in DIMENSIONS:
        units = statistics[statistics["dimension"].astype(str) == dimension].copy()
        if dimension == "parcel":
            units = units.assign(_parcel=units["name"].astype(int)).sort_values("_parcel")
        for _, row in units.iterrows():
            name = _unit_name(dimension, row["name"])
            parcel = int(name) if dimension == "parcel" else None
            q_any = _number(row, "q_m3_vs_m1", "q_age")
            q_beyond = _number(row, "q_m3_vs_m2")
            reasons = []
            if np.isfinite(q_any) and q_any < 0.05:
                reasons.append("significant_m3_vs_m1")
            if np.isfinite(q_beyond) and q_beyond < 0.05:
                reasons.append("significant_m3_vs_m2")
            if parcel in marginal_ids:
                reasons.append("marginal")
            if parcel in amygdala_ids:
                reasons.append("tian_amygdala")
            if dimension != "parcel":
                reasons.append(f"all_{dimension}")
            if not reasons:
                continue
            parcel_name = names.get(parcel, "") if parcel is not None else name
            rows.append({
                "dimension": dimension,
                "name": name,
                "parcel": parcel,
                "parcel_name": parcel_name if isinstance(parcel_name, str) else "",
                "q_age": q_any,
                "q_m3_vs_m1": q_any,
                "q_m3_vs_m2": q_beyond,
                "partial_r2": _number(row, "partial_r2_m3_vs_m1", "partial_r2"),
                "partial_r2_m3_vs_m1": _number(row, "partial_r2_m3_vs_m1", "partial_r2"),
                "partial_r2_m3_vs_m2": _number(row, "partial_r2_m3_vs_m2"),
                "reason": "+".join(reasons),
            })
    result = pd.DataFrame(rows, columns=TRAJECTORY_COLUMNS)
    if result.duplicated(["dimension", "name"]).any():
        raise ValueError("Age statistics have duplicate dimension/name rows")
    return result


def plot_one_region(age: np.ndarray, residual: np.ndarray, curve_age: np.ndarray,
                    effect: np.ndarray, ci_low: np.ndarray, ci_high: np.ndarray,
                    title: str, output_path: Path,
                    y_limits: tuple[float, float] = Y_LIMITS,
                    show_scatter: bool = True) -> None:
    """One view of the centered effect, with a 3:4 plotting box."""
    output_path.parent.mkdir(parents=True, exist_ok=True)
    figure, axis = plt.subplots(figsize=(4.8, 5.8))
    observed = np.isfinite(age) & np.isfinite(residual)
    if show_scatter:
        axis.scatter(
            age[observed], residual[observed],
            s=8, c=SCATTER_COLOR, alpha=0.4, linewidths=0, zorder=1, clip_on=True,
        )
    order = np.argsort(curve_age)
    band = np.isfinite(ci_low) & np.isfinite(ci_high)
    if band.any():
        axis.fill_between(
            curve_age[order], ci_low[order], ci_high[order],
            color=BAND_COLOR, alpha=0.4, linewidth=0, zorder=2,
        )
    finite_effect = np.isfinite(effect)
    if finite_effect.any():
        axis.plot(curve_age[order], effect[order], color=CURVE_COLOR, lw=1.8, zorder=3)
    axis.set_ylim(*y_limits)
    ages = np.concatenate([age[np.isfinite(age)], curve_age[np.isfinite(curve_age)]])
    if ages.size:
        xmin, xmax = float(ages.min()), float(ages.max())
        padding = 0.05 * (xmax - xmin) if xmax > xmin else 0.5
        axis.set_xlim(xmin - padding, xmax + padding)
    axis.set_box_aspect(BOX_ASPECT)
    axis.set_title(title, fontsize=9)
    axis.set_xlabel("Age")
    axis.set_ylabel("Partial residual (m3)" if show_scatter else "Centered age effect")
    axis.spines["top"].set_visible(False)
    axis.spines["right"].set_visible(False)
    figure.tight_layout()
    figure.savefig(output_path, format="svg")
    plt.close(figure)


def zoom_y_limits(effect: np.ndarray, ci_low: np.ndarray,
                  ci_high: np.ndarray) -> tuple[float, float]:
    """Fit the complete curve/CI envelope, including a small plotting margin."""
    values = np.concatenate([effect, ci_low, ci_high])
    values = values[np.isfinite(values)]
    if not values.size:
        return Y_LIMITS
    low, high = float(values.min()), float(values.max())
    span = max(high - low, 0.05 * max(abs(low), abs(high)), 1e-6)
    return low - 0.15 * span, high + 0.15 * span


def write_region_curves(statistics: pd.DataFrame, curves: pd.DataFrame,
                        residuals: pd.DataFrame, parcel_info: pd.DataFrame,
                        marginal: list[int], amygdala: list[int],
                        figure_dir: Path) -> pd.DataFrame:
    """Write both SVG views and return the selection plus limits/counts.

    New residual tables are long: sub_id, age, dimension, name, residual. Old
    wide parcel tables remain readable for redraws of earlier parcel fits.
    """
    selected = trajectory_rows(statistics, parcel_info, marginal, amygdala)
    required = {"effect", "ci_low", "ci_high", "dimension", "name", "age"}
    missing = sorted(required - set(curves.columns))
    if missing:
        raise ValueError(f"pygam_curves.csv is missing {', '.join(missing)}. Refit this run.")
    unit_curves = curves.copy()
    unit_curves["name"] = [
        _unit_name(str(dimension), name)
        for dimension, name in zip(unit_curves["dimension"], unit_curves["name"])
    ]
    long_residuals = {"dimension", "name", "residual"} <= set(residuals.columns)
    residual_table = residuals.copy()
    if long_residuals:
        residual_table["name"] = [
            _unit_name(str(dimension), name)
            for dimension, name in zip(residual_table["dimension"], residual_table["name"])
        ]
    else:
        selected = selected[selected["dimension"] == "parcel"].reset_index(drop=True)
        print("Old wide residual table: redrawing parcels only. Refit to save network, system, and Juelich residuals.")
    selected["n_observed"] = 0
    selected["n_outside_y_limits"] = 0
    selected["y_min"] = Y_LIMITS[0]
    selected["y_max"] = Y_LIMITS[1]
    selected["figure_file"] = ""
    selected["zoom_figure_file"] = ""
    selected["zoom_y_min"] = np.nan
    selected["zoom_y_max"] = np.nan
    is_432 = set(parcel_info["parcel"].astype(int)) == set(range(1, 433))
    output_dir = figure_dir / "age_curves"
    written = set()
    for row in selected.itertuples():
        curve = unit_curves[
            (unit_curves["dimension"].astype(str) == row.dimension)
            & (unit_curves["name"] == row.name)
        ].sort_values("age")
        if curve.empty:
            raise ValueError(f"No saved age curve for {row.dimension} {row.name}. Refit this run.")
        if long_residuals:
            observations = residual_table[
                (residual_table["dimension"].astype(str) == row.dimension)
                & (residual_table["name"] == row.name)
            ]
            if observations.empty:
                raise ValueError(f"No saved residuals for {row.dimension} {row.name}. Refit this run.")
            age = observations["age"].to_numpy(dtype=float)
            residual = observations["residual"].to_numpy(dtype=float)
        elif row.dimension == "parcel" and row.name in residual_table:
            age = residual_table["age"].to_numpy(dtype=float)
            residual = residual_table[row.name].to_numpy(dtype=float)
        else:
            raise ValueError(f"No saved residuals for {row.dimension} {row.name}. Refit this run.")
        significant = any(np.isfinite(q) and q < 0.05
                          for q in (row.q_m3_vs_m1, row.q_m3_vs_m2))
        wide_range = is_432 and row.dimension == "parcel" and (
            significant or int(row.name) in (156, 157))
        y_limits = Y_LIMITS_432_SELECTED if wide_range else Y_LIMITS
        observed = np.isfinite(age) & np.isfinite(residual)
        outside = observed & ((residual < y_limits[0]) | (residual > y_limits[1]))
        selected.loc[row.Index, "n_observed"] = int(observed.sum())
        selected.loc[row.Index, "n_outside_y_limits"] = int(outside.sum())
        selected.loc[row.Index, "y_min"] = y_limits[0]
        selected.loc[row.Index, "y_max"] = y_limits[1]
        if row.dimension == "parcel":
            filename = f"parcel_{int(row.name):04d}.svg"
            label = f"{int(row.name)}  {row.parcel_name}".strip()
        else:
            slug = re.sub(r"[^A-Za-z0-9_.-]+", "_", row.name).strip("._") or "unit"
            prefix = "amygdala_Juelich" if row.dimension == "amygdala" else row.dimension
            filename = f"{prefix}_{slug}.svg"
            label = f"Juelich {row.name}" if row.dimension == "amygdala" else f"{row.dimension} {row.name}"
        if filename in written:
            raise ValueError(f"Age-curve names produce the same filename: {filename}")
        written.add(filename)
        selected.loc[row.Index, "figure_file"] = f"age_curves/{filename}"
        q_any = f"{row.q_m3_vs_m1:.3g}" if np.isfinite(row.q_m3_vs_m1) else "na"
        q_beyond = f"{row.q_m3_vs_m2:.3g}" if np.isfinite(row.q_m3_vs_m2) else "na"
        curve_age = curve["age"].to_numpy(dtype=float)
        effect = curve["effect"].to_numpy(dtype=float)
        ci_low = curve["ci_low"].to_numpy(dtype=float)
        ci_high = curve["ci_high"].to_numpy(dtype=float)
        title = f"{label}\nq(m3-m1)={q_any}   q(m3-m2)={q_beyond}   95% CI"
        plot_one_region(
            age, residual, curve_age, effect, ci_low, ci_high,
            title=title,
            output_path=output_dir / filename,
            y_limits=y_limits,
        )
        zoom_name = f"{Path(filename).stem}_zoom.svg"
        zoom_limits = zoom_y_limits(effect, ci_low, ci_high)
        selected.loc[row.Index, "zoom_figure_file"] = f"age_curves_zoom/{zoom_name}"
        selected.loc[row.Index, "zoom_y_min"] = zoom_limits[0]
        selected.loc[row.Index, "zoom_y_max"] = zoom_limits[1]
        plot_one_region(
            age, residual, curve_age, effect, ci_low, ci_high,
            title=title + "\nCurve detail",
            output_path=figure_dir / "age_curves_zoom" / zoom_name,
            y_limits=zoom_limits, show_scatter=False,
        )
    print(f"Age curves: {len(written)} in {output_dir}")
    print(f"Curve details: {len(written)} in {figure_dir / 'age_curves_zoom'}")
    print(f"Residual points outside each unit's y limits: {int(selected['n_outside_y_limits'].sum())}")
    return selected
