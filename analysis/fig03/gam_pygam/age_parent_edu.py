"""Correlation of age and parental education in the 400-parcel design.

The design uses the same subject rules as the Kong 2022 parcel-only run:
valid ER list, minus second-level removals, then parental education, then
both effectSize volumes, then mean FD for every remaining subject. Parcel
means are not read. Age and parent_edu in that table do not depend on the
atlas, so this is the design the 400-parcel model uses.

If inputs/design.csv is already present, that file is used. Pass --rebuild
to write it again from the source tables.

    python -m analysis.fig03.gam_pygam.age_parent_edu --output "outputs/fig03/gam_pygam_volume_400"

To build the design on the server when that file is not there yet:

    python -m analysis.fig03.gam_pygam.age_parent_edu --input "/public/home/dingrui/fmri_analysis/data/beh/participants_in_tfmri_demographics.csv" --subjects "/public/home/dingrui/fmri_analysis/data/beh/sub_list_valid_4_ER.txt" --excluded "/public/home/dingrui/fmri_analysis/data/fmri/subs_removed_from_2ndlevel.txt" --parent-edu "/public/home/dingrui/fmri_analysis/data/beh/env_family_student_cbcl.tsv" --meanfd "/public/home/dingrui/fmri_analysis/data/fmri/qc_res_subs_runs_fmri.csv" --fig3-volume-dir "/public/home/dingrui/fmri_analysis/data/fmri/res_1stlevel" --output "outputs/fig03/gam_pygam_volume_400"
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

from analysis.fig03.gam_pygam.run_pygam import resolve_project_path
from analysis.fig03.gam_supplement.prepare_inputs import (
    ORIGINAL_DEMOGRAPHICS,
    ORIGINAL_EXCLUDED,
    ORIGINAL_MEANFD,
    ORIGINAL_PARENT_EDU,
    ORIGINAL_SUBJECTS,
    ORIGINAL_VOLUME_DIR,
    PARENT_EDU_ID_COLUMN,
    PARENT_EDU_VALUE_COLUMN,
    attach_mean_fd,
    attach_parent_education,
    build_design_from_sources,
    effect_volume_paths,
)

PROJECT_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_OUTPUT = PROJECT_ROOT / "outputs" / "fig03" / "gam_pygam_volume_400"


def _two_sided(result) -> tuple[float, float]:
    if hasattr(result, "statistic"):
        return float(result.statistic), float(result.pvalue)
    return float(result[0]), float(result[1])


def age_parent_edu_association(design: pd.DataFrame) -> pd.DataFrame:
    """Pearson and Spearman association of age with parent_edu."""
    if "parent_edu" not in design.columns:
        raise ValueError(
            "The design table has no parent_edu column. "
            f"Columns: {list(design.columns)}"
        )
    age = pd.to_numeric(design["age"], errors="coerce").to_numpy(dtype=float)
    education = pd.to_numeric(design["parent_edu"], errors="coerce").to_numpy(dtype=float)
    finite = np.isfinite(age) & np.isfinite(education)
    if int(finite.sum()) < 3:
        raise ValueError("Fewer than 3 subjects have both age and parent_edu")
    age = age[finite]
    education = education[finite]
    pearson_r, pearson_p = _two_sided(stats.pearsonr(age, education))
    spearman_r, spearman_p = _two_sided(stats.spearmanr(age, education))
    return pd.DataFrame([{
        "n": int(age.size),
        "pearson_r": pearson_r,
        "pearson_p": pearson_p,
        "spearman_r": spearman_r,
        "spearman_p": spearman_p,
        "age_mean": float(np.mean(age)),
        "age_sd": float(np.std(age, ddof=1)),
        "parent_edu_mean": float(np.mean(education)),
        "parent_edu_sd": float(np.std(education, ddof=1)),
    }])


def prepare_design(demographics_path: Path, subjects_path: Path, excluded_path: Path,
                   parent_edu_path: Path, meanfd_path: Path, volume_dir: Path,
                   input_dir: Path) -> pd.DataFrame:
    """Write the 400-parcel subject table without reading voxel values."""
    metadata = build_design_from_sources(demographics_path, subjects_path, excluded_path)
    metadata, dropped_parent_edu = attach_parent_education(metadata, parent_edu_path)
    input_dir.mkdir(parents=True, exist_ok=True)
    dropped_parent_edu.to_csv(input_dir / "dropped_parent_edu.csv", index=False)
    (input_dir / "parent_edu_source.txt").write_text(
        "\n".join([
            f"file: {parent_edu_path}",
            f"id_column: {PARENT_EDU_ID_COLUMN}",
            f"value_column: {PARENT_EDU_VALUE_COLUMN}",
            "design_column: parent_edu",
            f"n_dropped: {len(dropped_parent_edu)}",
        ]) + "\n",
        encoding="utf-8",
    )
    metadata, _rpsl_paths, _lkng_paths, missing_volume = effect_volume_paths(volume_dir, metadata)
    metadata = attach_mean_fd(metadata, meanfd_path)
    (input_dir / "meanfd_source.txt").write_text(
        "\n".join([
            f"file: {meanfd_path}",
            "layout: run-level mean_fd averaged over EMO and SOC within sub_id",
            "design_column: meanFD",
            f"n_subjects: {len(metadata)}",
        ]) + "\n",
        encoding="utf-8",
    )
    (input_dir / "dropped_no_volume.txt").write_text(
        "\n".join(missing_volume) + ("\n" if missing_volume else ""), encoding="utf-8")
    design = metadata.loc[:, ["sub_id", "age", "gender", "site_id", "parent_edu", "meanFD"]].copy()
    design.to_csv(input_dir / "design.csv", index=False)
    print(f"Design: {input_dir / 'design.csv'} (n={len(design)})")
    return design


def _chinese_font() -> str | None:
    from matplotlib import font_manager
    for name in ("Microsoft YaHei", "Noto Sans CJK SC", "Source Han Sans SC", "WenQuanYi Micro Hei"):
        if any(found.name == name for found in font_manager.fontManager.ttflist):
            return name
    return None


def plot_association(design: pd.DataFrame, summary: pd.DataFrame, output_path: Path) -> None:
    font = _chinese_font()
    if font:
        plt.rcParams["font.sans-serif"] = [font]
        plt.rcParams["axes.unicode_minus"] = False
        xlabel, ylabel = "父母受教育水平", "年龄"
    else:
        xlabel, ylabel = "Parental education", "Age"
    age = design["age"].to_numpy(dtype=float)
    education = design["parent_edu"].to_numpy(dtype=float)
    row = summary.iloc[0]
    slope, intercept = np.polyfit(education, age, 1)
    x_line = np.linspace(float(np.min(education)), float(np.max(education)), 50)
    fig, ax = plt.subplots(figsize=(6.2, 4.6))
    ax.scatter(education, age, s=12, c="#3d6f99", alpha=0.45, linewidths=0)
    ax.plot(x_line, intercept + slope * x_line, color="#c45c26", lw=1.6)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(
        f"n={int(row['n'])}   Pearson r={row['pearson_r']:.3f} (p={row['pearson_p']:.3g})"
        f"\nSpearman r={row['spearman_r']:.3f} (p={row['spearman_p']:.3g})",
        fontsize=11,
    )
    fig.tight_layout()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    if output_path.suffix.lower() == ".svg":
        plt.rcParams["svg.fonttype"] = "none"
    fig.savefig(output_path, dpi=160)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description="Age and parent_edu correlation for the 400-parcel design")
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--rebuild", action="store_true",
                        help="Rebuild design.csv from the source tables even if it already exists.")
    parser.add_argument("--input", type=Path, default=ORIGINAL_DEMOGRAPHICS)
    parser.add_argument("--subjects", type=Path, default=ORIGINAL_SUBJECTS)
    parser.add_argument("--excluded", type=Path, default=ORIGINAL_EXCLUDED)
    parser.add_argument("--parent-edu", type=Path, default=ORIGINAL_PARENT_EDU)
    parser.add_argument("--meanfd", type=Path, default=ORIGINAL_MEANFD)
    parser.add_argument("--fig3-volume-dir", type=Path, default=ORIGINAL_VOLUME_DIR)
    args = parser.parse_args()
    output_dir = resolve_project_path(args.output)
    input_dir = output_dir / "inputs"
    design_path = input_dir / "design.csv"
    if design_path.is_file() and not args.rebuild:
        design = pd.read_csv(design_path)
        print(f"Using saved design: {design_path} (n={len(design)})")
    else:
        design = prepare_design(
            resolve_project_path(args.input),
            resolve_project_path(args.subjects),
            resolve_project_path(args.excluded),
            resolve_project_path(args.parent_edu),
            resolve_project_path(args.meanfd),
            resolve_project_path(args.fig3_volume_dir),
            input_dir,
        )
    summary = age_parent_edu_association(design)
    statistics_dir = output_dir / "statistics"
    statistics_dir.mkdir(parents=True, exist_ok=True)
    summary_path = statistics_dir / "age_parent_edu_correlation.csv"
    summary.to_csv(summary_path, index=False)
    figure_path = output_dir / "figures" / "age_parent_edu.svg"
    plot_association(design, summary, figure_path)
    row = summary.iloc[0]
    print(
        f"age 与 parent_edu：n={int(row['n'])}，"
        f"Pearson r={row['pearson_r']:.4f}，p={row['pearson_p']:.4g}，"
        f"Spearman r={row['spearman_r']:.4f}，p={row['spearman_p']:.4g}"
    )
    print(f"Table: {summary_path}")
    print(f"Figure: {figure_path}")


if __name__ == "__main__":
    main()
