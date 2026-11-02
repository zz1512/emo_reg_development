"""One age-curve SVG per parcel.

The points are residuals after the intercept and the four covariates. The line
is the full-model age effect at the mean of those covariates, with its 95%
pointwise confidence interval. The top and right spines are omitted.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
matplotlib.rcParams["svg.fonttype"] = "none"
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

CURVE_COLOR = "#4C78A8"


def trajectory_rows(statistics: pd.DataFrame, parcel_info: pd.DataFrame,
                    marginal: list[int], amygdala: list[int]) -> pd.DataFrame:
    """Parcels that each get their own figure: significant, marginal, and amygdala."""
    parcels = statistics[statistics["dimension"].astype(str) == "parcel"].copy()
    parcels["parcel"] = parcels["name"].astype(int)
    info = parcel_info.copy()
    info["parcel"] = info["parcel"].astype(int)
    name_column = "name" if "name" in info.columns else None
    if name_column is not None:
        info = info.rename(columns={"name": "parcel_name"})
    keep = [column for column in ("parcel", "parcel_name", "hemi", "structure") if column in info.columns]
    parcels = parcels.merge(info.loc[:, keep], on="parcel", how="left")
    marginal_ids = {int(parcel) for parcel in marginal}
    amygdala_ids = {int(parcel) for parcel in amygdala}
    rows = []
    for row in parcels.itertuples():
        reasons = []
        q_value = float(row.q_age) if np.isfinite(row.q_age) else np.nan
        if np.isfinite(q_value) and q_value < 0.05:
            reasons.append("significant")
        if int(row.parcel) in marginal_ids:
            reasons.append("marginal")
        if int(row.parcel) in amygdala_ids:
            reasons.append("amygdala")
        if not reasons:
            continue
        parcel_name = ""
        if hasattr(row, "parcel_name") and isinstance(row.parcel_name, str):
            parcel_name = row.parcel_name
        rows.append({
            "parcel": int(row.parcel),
            "parcel_name": parcel_name,
            "q_age": q_value,
            "partial_r2": float(row.partial_r2) if np.isfinite(row.partial_r2) else np.nan,
            "reason": "+".join(reasons),
        })
    if not rows:
        return pd.DataFrame(columns=["parcel", "parcel_name", "q_age", "partial_r2", "reason"])
    return pd.DataFrame(rows).sort_values("parcel").reset_index(drop=True)


def plot_one_region(age: np.ndarray, residual: np.ndarray, curve_age: np.ndarray,
                    effect: np.ndarray, ci_low: np.ndarray, ci_high: np.ndarray,
                    title: str, output_path: Path) -> None:
    """Scatter, age curve, and 95% interval. No zero line. No top or right spine."""
    output_path.parent.mkdir(parents=True, exist_ok=True)
    figure, axis = plt.subplots(figsize=(4.4, 3.5))
    observed = np.isfinite(age) & np.isfinite(residual)
    axis.scatter(
        age[observed], residual[observed],
        s=10, c=CURVE_COLOR, alpha=0.35, linewidths=0, zorder=2,
    )
    order = np.argsort(curve_age)
    band = np.isfinite(ci_low) & np.isfinite(ci_high)
    if band.any():
        axis.fill_between(
            curve_age[order], ci_low[order], ci_high[order],
            color=CURVE_COLOR, alpha=0.2, linewidth=0, zorder=1,
        )
    finite_effect = np.isfinite(effect)
    if finite_effect.any():
        axis.plot(curve_age[order], effect[order], color=CURVE_COLOR, lw=1.8, zorder=3)
    axis.set_title(title, fontsize=10)
    axis.set_xlabel("Age")
    axis.set_ylabel("Residual after covariates")
    axis.spines["top"].set_visible(False)
    axis.spines["right"].set_visible(False)
    figure.tight_layout()
    figure.savefig(output_path, format="svg")
    plt.close(figure)


def write_region_curves(statistics: pd.DataFrame, curves: pd.DataFrame,
                        residuals: pd.DataFrame, parcel_info: pd.DataFrame,
                        marginal: list[int], amygdala: list[int],
                        figure_dir: Path) -> pd.DataFrame:
    """Write one SVG for every selected parcel and return the parcel list."""
    selected = trajectory_rows(statistics, parcel_info, marginal, amygdala)
    if "effect" not in curves.columns:
        raise ValueError("pygam_curves.csv has no age-effect column. Refit this Tian run.")
    parcel_curves = curves[curves["dimension"].astype(str) == "parcel"].copy()
    parcel_curves["parcel"] = parcel_curves["name"].astype(int)
    output_dir = figure_dir / "age_curves"
    for row in selected.itertuples():
        curve = parcel_curves[parcel_curves["parcel"] == int(row.parcel)].sort_values("age")
        column = str(int(row.parcel))
        if column not in residuals.columns or curve.empty:
            print(f"No residual curve for parcel {column}.")
            continue
        q_text = f"{row.q_age:.3g}" if np.isfinite(row.q_age) else "na"
        name = f"  {row.parcel_name}" if row.parcel_name else ""
        plot_one_region(
            residuals["age"].to_numpy(dtype=float),
            residuals[column].to_numpy(dtype=float),
            curve["age"].to_numpy(dtype=float),
            curve["effect"].to_numpy(dtype=float),
            curve["ci_low"].to_numpy(dtype=float),
            curve["ci_high"].to_numpy(dtype=float),
            title=f"{int(row.parcel)}{name}   q={q_text}   95% CI",
            output_path=output_dir / f"parcel_{int(row.parcel):04d}.svg",
        )
    print(f"Region curves: {len(selected)} in {output_dir}")
    return selected
