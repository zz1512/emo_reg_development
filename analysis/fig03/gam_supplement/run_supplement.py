"""Two separate steps for the k=3 mgcv age supplement.

Run from the repository root. Step 1 reports the linear age term and the
smooth term. Step 2 builds the comparison matrices and the trajectory figures.
It reads the tables written by step 1.

Real data, using the effectSize volumes read by the system-identification script:

    python -m analysis.fig03.gam_supplement.run_supplement --step terms --input "/public/home/dingrui/fmri_analysis/data/beh/participants_in_tfmri_demographics.csv" --subjects "/public/home/dingrui/fmri_analysis/data/beh/sub_list_valid_4_ER.txt" --excluded "/public/home/dingrui/fmri_analysis/data/fmri/subs_removed_from_2ndlevel.txt" --parent-edu "/public/home/dingrui/fmri_analysis/data/beh/env_family_student_cbcl.tsv" --meanfd "/public/home/dingrui/fmri_analysis/data/fmri/qc_res_subs_runs_fmri.csv" --fig3-volume-dir "/public/home/dingrui/fmri_analysis/data/fmri/res_1stlevel" --fig3-volume-atlas "resources/atlases/schaefer2018/Schaefer2018_200Parcels_7Networks_order_FSLMNI152_2mm.nii.gz" --schaefer-atlas "resources/atlases/schaefer2018/Schaefer2018_200Parcels_7Networks_order.dlabel.nii" --system-mask-dir "/public/home/dingrui/fmri_analysis/res/res_bayes" --output "outputs/fig03/gam_supplement_volume"

    python -m analysis.fig03.gam_supplement.run_supplement --step matrix --schaefer-atlas "resources/atlases/schaefer2018/Schaefer2018_200Parcels_7Networks_order.dlabel.nii" --output "outputs/fig03/gam_supplement_volume"

A finished design table can replace the demographics and subject lists. Pass it
with --input; it must contain sub_id, age, gender, site_id, parent_edu, and meanFD.
Parental education otherwise comes from env_family_student_cbcl.tsv, column
父母受教育水平_父母受教育水平, aligned by 被试编号. meanFD comes from
qc_res_subs_runs_fmri.csv. meanFD is the mean of EMO and SOC mean_fd within sub_id.
"""

from __future__ import annotations

import argparse
import json
import shutil
import subprocess
from pathlib import Path

import pandas as pd

from analysis.fig03.gam_supplement.plot_suites import (
    plot_all,
    plot_linearity,
    plot_parcel_pair_heatmaps,
    write_linearity_table,
    write_parcel_pair_matrices,
    write_square_matrices,
)
from analysis.fig03.gam_supplement.prepare_inputs import (
    ORIGINAL_DEMOGRAPHICS,
    ORIGINAL_EXCLUDED,
    ORIGINAL_MEANFD,
    ORIGINAL_PARENT_EDU,
    ORIGINAL_SUBJECTS,
    ORIGINAL_SURFACE_DIR,
    ORIGINAL_SYSTEM_MASKS,
    ORIGINAL_VOLUME_DIR,
    PARENT_EDU_ID_COLUMN,
    PARENT_EDU_VALUE_COLUMN,
    attach_mean_fd,
    attach_parent_education,
    build_design_from_sources,
    effect_volume_paths,
    make_toy_volumes,
    prepare_from_effect_volumes,
    prepare_from_volumes,
    read_design_tsv,
    schaefer_parcel_info,
    write_toy_system_masks,
)

PROJECT_ROOT = Path(__file__).resolve().parents[3]
R_SCRIPT = Path(__file__).resolve().parent / "mgcv_age_decomposition.R"
DEFAULT_ATLAS = PROJECT_ROOT / "resources/atlases/schaefer2018/Schaefer2018_200Parcels_7Networks_order_FSLMNI152_2mm.nii.gz"
DEFAULT_DLABEL = PROJECT_ROOT / "resources/atlases/schaefer2018/Schaefer2018_200Parcels_7Networks_order.dlabel.nii"
DEFAULT_RSCRIPT = Path(r"C:\Program Files\R\R-4.3.2\bin\Rscript.exe")
DEFAULT_OUTPUT = PROJECT_ROOT / "outputs" / "fig03" / "gam_supplement"


def find_rscript(explicit: Path | None) -> Path:
    if explicit is not None:
        if not explicit.is_file():
            raise FileNotFoundError(f"Rscript not found: {explicit}")
        return explicit
    discovered = shutil.which("Rscript")
    if discovered:
        return Path(discovered)
    if DEFAULT_RSCRIPT.is_file():
        return DEFAULT_RSCRIPT
    raise FileNotFoundError("Rscript was not found on PATH or in R 4.3.2")


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


def fit_models(rscript: Path, input_dir: Path, result_dir: Path, mode: str) -> None:
    result_dir.mkdir(parents=True, exist_ok=True)
    completed = subprocess.run(
        [str(rscript), str(R_SCRIPT), str(input_dir), str(result_dir), mode],
        check=False,
    )
    if completed.returncode != 0:
        raise RuntimeError(f"mgcv {mode} step failed with exit code {completed.returncode}")


def prepare_real_inputs(input_path: Path | None, subjects_path: Path, excluded_path: Path,
                        parent_edu_path: Path, meanfd_path: Path, volume_dir: Path,
                        volume_atlas: Path, dlabel_path: Path, system_mask_dir: Path,
                        input_dir: Path) -> dict[str, int]:
    source = ORIGINAL_DEMOGRAPHICS if input_path is None else input_path
    if subjects_path.is_file() and excluded_path.is_file():
        metadata = build_design_from_sources(source, subjects_path, excluded_path)
        metadata, dropped_parent_edu = attach_parent_education(metadata, parent_edu_path)
    else:
        metadata = read_design_tsv(source)
        dropped_parent_edu = pd.DataFrame(columns=["sub_id", "reason"])
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
    metadata, rpsl_paths, lkng_paths, missing_volume = effect_volume_paths(volume_dir, metadata)
    if subjects_path.is_file() and excluded_path.is_file():
        metadata = attach_mean_fd(metadata, meanfd_path)
    elif "meanFD" not in metadata.columns:
        raise ValueError("A finished design table must contain meanFD")
    (input_dir / "meanfd_source.txt").write_text(
        "\n".join([
            f"file: {meanfd_path}",
            "layout: run-level mean_fd averaged over EMO and SOC within sub_id",
            "design_column: meanFD",
            f"n_subjects: {len(metadata)}",
        ]) + "\n",
        encoding="utf-8",
    )
    input_dir.mkdir(parents=True, exist_ok=True)
    (input_dir / "dropped_no_volume.txt").write_text(
        "\n".join(missing_volume) + ("\n" if missing_volume else ""), encoding="utf-8")
    counts = prepare_from_effect_volumes(
        metadata, rpsl_paths, lkng_paths, volume_atlas, dlabel_path, system_mask_dir, input_dir)
    counts["n_dropped_no_volume"] = len(missing_volume)
    counts["n_dropped_parent_edu"] = int(len(dropped_parent_edu))
    return counts


def prepare_toy_inputs(volume_atlas: Path, dlabel_path: Path, input_dir: Path) -> dict[str, int]:
    import tempfile

    import nibabel as nib

    with tempfile.TemporaryDirectory(prefix="gam-supplement-") as temp_dir:
        metadata, paths = make_toy_volumes(Path(temp_dir), volume_atlas)
        mask_dir = Path(temp_dir) / "system_masks"
        write_toy_system_masks(nib.load(str(volume_atlas)), schaefer_parcel_info(dlabel_path), mask_dir)
        return prepare_from_volumes(
            metadata, paths["rpsl"], paths["lkng"], volume_atlas, dlabel_path, mask_dir, input_dir)


def write_step_metadata(output_dir: Path, step: str, payload: dict) -> None:
    path = output_dir / "run_metadata.json"
    existing = {}
    if path.is_file():
        existing = json.loads(path.read_text(encoding="utf-8"))
    steps = existing.get("steps", {})
    steps[step] = payload
    existing.update({
        "model": "f(age) = beta*age + gamma*h(age) + gender + site_id + parent_edu + meanFD",
        "h_age": "direction in the cubic-regression basis s(age, k=3) orthogonal to intercept and age",
        "penalty": "cubic-regression penalty of that direction, REML",
        "inference": "Wald p-values of beta and gamma; BH FDR within parcels and within each difference family",
        "anova": False,
        "covariates": ["gender", "site_id", "parent_edu", "meanFD"],
        "steps": steps,
    })
    path.write_text(json.dumps(existing, indent=2), encoding="utf-8")


def run_terms(*, toy: bool, input_path: Path | None, subjects_path: Path, excluded_path: Path,
              parent_edu_path: Path, meanfd_path: Path, volume_dir: Path, volume_atlas: Path,
              dlabel_path: Path, system_mask_dir: Path, output_dir: Path, rscript: Path) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    input_dir = output_dir / "inputs"
    result_dir = output_dir / "statistics"
    if toy:
        counts = prepare_toy_inputs(volume_atlas, dlabel_path, input_dir)
    else:
        counts = prepare_real_inputs(
            input_path, subjects_path, excluded_path, parent_edu_path, meanfd_path,
            volume_dir, volume_atlas, dlabel_path, system_mask_dir, input_dir)
    fit_models(rscript, input_dir, result_dir, "terms")
    linearity_table = write_linearity_table(result_dir)
    linearity_figures = plot_linearity(result_dir, output_dir / "figures", dlabel_path)
    write_step_metadata(output_dir, "terms", {
        "counts": counts,
        "statistics": [
            str(result_dir / "parcel_gam_statistics.csv"),
            str(result_dir / "group_gam_statistics.csv"),
            str(linearity_table),
        ],
        "figures": [str(path) for path in linearity_figures],
    })
    print(f"Linearity table: {linearity_table}")
    print(f"Linearity figures: {output_dir / 'figures'}")


def run_matrix(*, output_dir: Path, dlabel_path: Path, rscript: Path) -> None:
    input_dir = output_dir / "inputs"
    result_dir = output_dir / "statistics"
    figure_dir = output_dir / "figures"
    required = (
        input_dir / "design.csv",
        input_dir / "group_response.csv",
        input_dir / "parcel_info.csv",
        result_dir / "parcel_gam_statistics.csv",
        result_dir / "parcel_curves.csv",
    )
    missing = [str(path) for path in required if not path.is_file()]
    if missing:
        raise FileNotFoundError(
            "The matrix step needs the tables from --step terms. Missing: " + ", ".join(missing)
        )
    fit_models(rscript, input_dir, result_dir, "matrix")
    difference_stats = pd.read_csv(result_dir / "difference_statistics.csv")
    matrix_files = write_square_matrices(difference_stats, result_dir)
    matrix_files.extend(write_parcel_pair_matrices(difference_stats, result_dir))
    figures = plot_all(input_dir, result_dir, figure_dir, dlabel_path)
    figures.extend(plot_parcel_pair_heatmaps(difference_stats, figure_dir))
    write_step_metadata(output_dir, "matrix", {
        "matrices": [str(path) for path in matrix_files],
        "figures": [str(path) for path in figures],
    })
    print(f"Comparison matrices: {result_dir}")
    print(f"Trajectory figures: {figure_dir}")


def main() -> None:
    parser = argparse.ArgumentParser(description="mgcv k=3 age decomposition, in two steps")
    parser.add_argument("--step", choices=("terms", "matrix"), required=True,
                        help="terms reports beta and h(age); matrix builds comparisons and figures")
    parser.add_argument("--toy", action="store_true")
    parser.add_argument("--input", type=Path,
                        help="Fig. 3 design table, or the original demographics CSV")
    parser.add_argument("--subjects", type=Path, default=ORIGINAL_SUBJECTS,
                        help="Included subject list from the original Figure 3 script")
    parser.add_argument("--excluded", type=Path, default=ORIGINAL_EXCLUDED,
                        help="Subjects removed from the second-level analysis")
    parser.add_argument("--parent-edu", type=Path, default=ORIGINAL_PARENT_EDU,
                        help="Family table env_family_student_cbcl.tsv; 被试编号 aligned to sub_id")
    parser.add_argument("--meanfd", type=Path, default=ORIGINAL_MEANFD,
                        help="Run-level QC csv; mean_fd is averaged over EMO and SOC for each sub_id")
    parser.add_argument("--fig3-volume-dir", type=Path, default=ORIGINAL_VOLUME_DIR,
                        help="res_1stlevel folders with effectSize-reg/pas_neg_Observe_across_runs.nii.gz")
    parser.add_argument("--fig3-surface-dir", type=Path, default=ORIGINAL_SURFACE_DIR,
                        help="Accepted for the old command. The real step reads effectSize volumes")
    parser.add_argument("--fig3-volume-atlas", type=Path, default=DEFAULT_ATLAS,
                        help="Schaefer volume atlas, nearest-resampled once onto the effect-size grid")
    parser.add_argument("--schaefer-atlas", type=Path, default=DEFAULT_DLABEL,
                        help="Schaefer dlabel used for parcel names, Yeo networks, and surface figures")
    parser.add_argument("--system-mask-dir", type=Path, default=ORIGINAL_SYSTEM_MASKS,
                        help="Directory of the four mask_<system>.nii.gz files")
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--rscript", type=Path)
    args = parser.parse_args()
    if args.step == "matrix" and args.toy:
        parser.error("--toy prepares data inside --step terms. Run terms first, then matrix on the same --output")
    rscript = find_rscript(args.rscript)
    output_dir = resolve_project_path(args.output)
    if args.step == "terms":
        run_terms(
            toy=args.toy,
            input_path=None if args.input is None else resolve_project_path(args.input),
            subjects_path=resolve_project_path(args.subjects),
            excluded_path=resolve_project_path(args.excluded),
            parent_edu_path=resolve_project_path(args.parent_edu),
            meanfd_path=resolve_project_path(args.meanfd),
            volume_dir=resolve_project_path(args.fig3_volume_dir),
            volume_atlas=resolve_project_path(args.fig3_volume_atlas),
            dlabel_path=resolve_project_path(args.schaefer_atlas),
            system_mask_dir=resolve_project_path(args.system_mask_dir),
            output_dir=output_dir,
            rscript=rscript,
        )
    else:
        run_matrix(
            output_dir=output_dir,
            dlabel_path=resolve_project_path(args.schaefer_atlas),
            rscript=rscript,
        )


if __name__ == "__main__":
    main()
