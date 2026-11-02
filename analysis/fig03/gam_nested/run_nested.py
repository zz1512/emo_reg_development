"""Fit m1, m2, and m3 with mgcv on volume Reappraise minus Look.

m1 has the four covariates and no age term. m2 adds linear age. m3 replaces
that line with s(age, k=3, bs="cr"). Nested F tests compare m2 with m1, m3
with m1, and m3 with m2. BH correction stays inside parcels, Yeo networks,
systems, and the six amygdala ROIs, and each comparison is corrected on its own.

Run from the repository root:

    python -m analysis.fig03.gam_nested.run_nested --input "/public/home/dingrui/fmri_analysis/data/beh/participants_in_tfmri_demographics.csv" --subjects "/public/home/dingrui/fmri_analysis/data/beh/sub_list_valid_4_ER.txt" --excluded "/public/home/dingrui/fmri_analysis/data/fmri/subs_removed_from_2ndlevel.txt" --parent-edu "/public/home/dingrui/fmri_analysis/data/beh/env_family_student_cbcl.tsv" --meanfd "/public/home/dingrui/fmri_analysis/data/fmri/qc_res_subs_runs_fmri.csv" --fig3-volume-dir "/public/home/dingrui/fmri_analysis/data/fmri/res_1stlevel" --fig3-volume-atlas "resources/atlases/schaefer2018/Schaefer2018_200Parcels_7Networks_order_FSLMNI152_2mm.nii.gz" --schaefer-atlas "resources/atlases/schaefer2018/Schaefer2018_200Parcels_7Networks_order.dlabel.nii" --system-mask-dir "/public/home/dingrui/fmri_analysis/res/res_bayes" --output "outputs/fig03/gam_nested_volume"

Schaefer 200 x Yeo 7 with no meanFD magnitude filter. Defaults stay REML,
cr, and k=3. The result is outputs/fig03/gam_nested_volume_200. This does not
write gam_nested_volume or the _FD05, _FD02, and _shared_inputs directories.

    python -m analysis.fig03.gam_nested.run_nested --input "/public/home/dingrui/fmri_analysis/data/beh/participants_in_tfmri_demographics.csv" --subjects "/public/home/dingrui/fmri_analysis/data/beh/sub_list_valid_4_ER.txt" --excluded "/public/home/dingrui/fmri_analysis/data/fmri/subs_removed_from_2ndlevel.txt" --parent-edu "/public/home/dingrui/fmri_analysis/data/beh/env_family_student_cbcl.tsv" --meanfd "/public/home/dingrui/fmri_analysis/data/fmri/qc_res_subs_runs_fmri.csv" --fig3-volume-dir "/public/home/dingrui/fmri_analysis/data/fmri/res_1stlevel" --fig3-volume-atlas "resources/atlases/schaefer2018/Schaefer2018_200Parcels_7Networks_order_FSLMNI152_2mm.nii.gz" --schaefer-atlas "resources/atlases/schaefer2018/Schaefer2018_200Parcels_7Networks_order.dlabel.nii" --system-mask-dir "/public/home/dingrui/fmri_analysis/res/res_bayes" --output "outputs/fig03/gam_nested_volume_200"

Schaefer 200 x Yeo 7 with meanFD ceilings 0.5 and 0.2. Defaults stay REML,
cr, and k=3. Volumes are read once. This does not write gam_nested_volume.

    python -m analysis.fig03.gam_nested.run_nested --meanfd-max 0.5 0.2 --input "/public/home/dingrui/fmri_analysis/data/beh/participants_in_tfmri_demographics.csv" --subjects "/public/home/dingrui/fmri_analysis/data/beh/sub_list_valid_4_ER.txt" --excluded "/public/home/dingrui/fmri_analysis/data/fmri/subs_removed_from_2ndlevel.txt" --parent-edu "/public/home/dingrui/fmri_analysis/data/beh/env_family_student_cbcl.tsv" --meanfd "/public/home/dingrui/fmri_analysis/data/fmri/qc_res_subs_runs_fmri.csv" --fig3-volume-dir "/public/home/dingrui/fmri_analysis/data/fmri/res_1stlevel" --fig3-volume-atlas "resources/atlases/schaefer2018/Schaefer2018_200Parcels_7Networks_order_FSLMNI152_2mm.nii.gz" --schaefer-atlas "resources/atlases/schaefer2018/Schaefer2018_200Parcels_7Networks_order.dlabel.nii" --system-mask-dir "/public/home/dingrui/fmri_analysis/res/res_bayes" --output "outputs/fig03/gam_nested_volume_200"

Schaefer 400 x Yeo 7 with no meanFD magnitude filter. Defaults stay REML, cr,
and k=3. The result is outputs/fig03/gam_nested_volume_400. This does not write
the _FD05, _FD02, _shared_inputs, gcv, gcv_ps5, reml, or 17net directories.

    python -m analysis.fig03.gam_nested.run_nested --input "/public/home/dingrui/fmri_analysis/data/beh/participants_in_tfmri_demographics.csv" --subjects "/public/home/dingrui/fmri_analysis/data/beh/sub_list_valid_4_ER.txt" --excluded "/public/home/dingrui/fmri_analysis/data/fmri/subs_removed_from_2ndlevel.txt" --parent-edu "/public/home/dingrui/fmri_analysis/data/beh/env_family_student_cbcl.tsv" --meanfd "/public/home/dingrui/fmri_analysis/data/fmri/qc_res_subs_runs_fmri.csv" --fig3-volume-dir "/public/home/dingrui/fmri_analysis/data/fmri/res_1stlevel" --fig3-volume-atlas "resources/atlases/schaefer2018/Schaefer2018_400Parcels_7Networks_order_FSLMNI152_2mm.nii.gz" --schaefer-atlas "resources/atlases/schaefer2018/Schaefer2018_400Parcels_7Networks_order.dlabel.nii" --system-mask-dir "/public/home/dingrui/fmri_analysis/res/res_bayes" --output "outputs/fig03/gam_nested_volume_400"

The 400-parcel, four-covariate run keeps parent_edu and writes one directory
per meanFD threshold. Subjects with meanFD above the threshold are dropped.
Volumes are read once into outputs/fig03/gam_nested_volume_400_shared_inputs.

    python -m analysis.fig03.gam_nested.run_nested --meanfd-max 0.5 0.2 --input "/public/home/dingrui/fmri_analysis/data/beh/participants_in_tfmri_demographics.csv" --subjects "/public/home/dingrui/fmri_analysis/data/beh/sub_list_valid_4_ER.txt" --excluded "/public/home/dingrui/fmri_analysis/data/fmri/subs_removed_from_2ndlevel.txt" --parent-edu "/public/home/dingrui/fmri_analysis/data/beh/env_family_student_cbcl.tsv" --meanfd "/public/home/dingrui/fmri_analysis/data/fmri/qc_res_subs_runs_fmri.csv" --fig3-volume-dir "/public/home/dingrui/fmri_analysis/data/fmri/res_1stlevel" --fig3-volume-atlas "resources/atlases/schaefer2018/Schaefer2018_400Parcels_7Networks_order_FSLMNI152_2mm.nii.gz" --schaefer-atlas "resources/atlases/schaefer2018/Schaefer2018_400Parcels_7Networks_order.dlabel.nii" --system-mask-dir "/public/home/dingrui/fmri_analysis/res/res_bayes" --output "outputs/fig03/gam_nested_volume_400"

Kong 2022 17-network 400 parcels are parcel means only. Pass --parcels-only.
That run does not fit the 17 networks, the four systems, or the amygdala.

    python -m analysis.fig03.gam_nested.run_nested --parcels-only --meanfd-max 0.5 0.2 --input "/public/home/dingrui/fmri_analysis/data/beh/participants_in_tfmri_demographics.csv" --subjects "/public/home/dingrui/fmri_analysis/data/beh/sub_list_valid_4_ER.txt" --excluded "/public/home/dingrui/fmri_analysis/data/fmri/subs_removed_from_2ndlevel.txt" --parent-edu "/public/home/dingrui/fmri_analysis/data/beh/env_family_student_cbcl.tsv" --meanfd "/public/home/dingrui/fmri_analysis/data/fmri/qc_res_subs_runs_fmri.csv" --fig3-volume-dir "/public/home/dingrui/fmri_analysis/data/fmri/res_1stlevel" --fig3-volume-atlas "resources/atlases/schaefer2018/Schaefer2018_400Parcels_Kong2022_17Networks_order_FSLMNI152_2mm.nii.gz" --parcel-info "resources/atlases/schaefer2018/Schaefer2018_400Parcels_Kong2022_17Networks_order_info.txt" --schaefer-atlas "resources/atlases/schaefer2018/Schaefer2018_400Parcels_Kong2022_17Networks_order.dlabel.nii" --output "outputs/fig03/gam_nested_volume_400_17net"
"""

from __future__ import annotations

import argparse
import json
import shutil
import subprocess
from pathlib import Path

import pandas as pd

from analysis.fig03.gam_nested.amygdala import (
    AMYGDALA_ORDER,
    MIN_AMYGDALA_VOXELS,
    amygdala_masks,
    load_juelich_atlas,
)
from analysis.fig03.gam_nested.plot_nested import (
    plot_age_parent_edu,
    plot_curve_panels,
    plot_m3_m1_parcel_trajectories,
    plot_m3_slope_matrices,
    write_slope_matrices,
)
from analysis.fig03.gam_supplement.prepare_inputs import (
    ORIGINAL_DEMOGRAPHICS,
    ORIGINAL_EXCLUDED,
    ORIGINAL_MEANFD,
    ORIGINAL_PARENT_EDU,
    ORIGINAL_SUBJECTS,
    ORIGINAL_SYSTEM_MASKS,
    ORIGINAL_VOLUME_DIR,
    PARENT_EDU_ID_COLUMN,
    PARENT_EDU_VALUE_COLUMN,
    analysis_design,
    apply_meanfd_exclusion,
    attach_mean_fd,
    attach_parent_education,
    build_design_from_sources,
    effect_volume_paths,
    ensure_subject_roster,
    n_analysis_subjects,
    write_subject_roster,
    prepare_from_effect_volumes,
    read_design_tsv,
    read_schaefer_order_info,
)

PROJECT_ROOT = Path(__file__).resolve().parents[3]
R_SCRIPT = Path(__file__).resolve().parent / "mgcv_nested_models.R"
SLOPE_SCRIPT = Path(__file__).resolve().parent / "mgcv_m3_pair_slopes.R"
DEFAULT_ATLAS = PROJECT_ROOT / "resources/atlases/schaefer2018/Schaefer2018_200Parcels_7Networks_order_FSLMNI152_2mm.nii.gz"
DEFAULT_DLABEL = PROJECT_ROOT / "resources/atlases/schaefer2018/Schaefer2018_200Parcels_7Networks_order.dlabel.nii"
DEFAULT_RSCRIPT = Path(r"C:\Program Files\R\R-4.3.2\bin\Rscript.exe")
DEFAULT_OUTPUT = PROJECT_ROOT / "outputs" / "fig03" / "gam_nested_volume"


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


def parcel_info_for(volume_atlas: Path, requested: Path | None) -> Path:
    """Kong 2022 names live in order_info.txt. The 7-network dlabel parser rejects them."""
    if requested is not None:
        if not requested.is_file():
            raise FileNotFoundError(f"Parcel info not found: {requested}")
        return requested
    sibling = Path(str(volume_atlas).replace("_FSLMNI152_2mm.nii.gz", "_info.txt"))
    if sibling.is_file():
        return sibling
    raise FileNotFoundError(
        "Parcel-only analysis needs --parcel-info pointing at the Schaefer order_info.txt"
    )


def prepare_real_inputs(input_path: Path | None, subjects_path: Path, excluded_path: Path,
                        parent_edu_path: Path, meanfd_path: Path, volume_dir: Path,
                        volume_atlas: Path, dlabel_path: Path, system_mask_dir: Path,
                        juelich_atlas: Path | None, juelich_labels: Path | None,
                        input_dir: Path, parcels_only: bool = False,
                        parcel_info_path: Path | None = None) -> dict[str, int]:
    source = ORIGINAL_DEMOGRAPHICS if input_path is None else input_path
    if subjects_path.is_file() and excluded_path.is_file():
        metadata = build_design_from_sources(source, subjects_path, excluded_path)
        write_subject_roster(input_dir, metadata, parent_edu_path, volume_dir, meanfd_path)
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
    (input_dir / "dropped_no_volume.txt").write_text(
        "\n".join(missing_volume) + ("\n" if missing_volume else ""), encoding="utf-8")
    parcel_table = None
    extra_masks = None
    if parcels_only:
        info_path = parcel_info_for(volume_atlas, parcel_info_path)
        parcel_table = read_schaefer_order_info(info_path)
        print(f"Parcel-only atlas names: {len(parcel_table)} rows from {info_path}")
    else:
        atlas_image, labels, atlas_source = load_juelich_atlas(juelich_atlas, juelich_labels)
        extra_masks = amygdala_masks(atlas_image, labels)
        if tuple(extra_masks) != AMYGDALA_ORDER:
            raise RuntimeError(f"Amygdala ROI order is {tuple(extra_masks)}, expected {AMYGDALA_ORDER}")
        (input_dir / "amygdala_source.txt").write_text(
            "\n".join([
                f"atlas: {atlas_source}",
                "atlas_name: Juelich-maxprob-thr25-2mm",
                "nuclei: LB, CM, SF from Amunts et al. 2005",
                "hemisphere: MNI x < 0 is left and x > 0 is right when the label is bilateral",
                "contrast: effectSize-reg_neg_Observe_across_runs minus effectSize-pas_neg_Observe_across_runs",
                f"minimum_voxels: {MIN_AMYGDALA_VOXELS}",
            ]) + "\n",
            encoding="utf-8",
        )
    counts = prepare_from_effect_volumes(
        metadata, rpsl_paths, lkng_paths, volume_atlas, dlabel_path, system_mask_dir, input_dir,
        extra_masks=extra_masks, min_extra_voxels=MIN_AMYGDALA_VOXELS,
        parcel_info_table=parcel_table, include_groups=not parcels_only)
    counts["n_dropped_no_volume"] = len(missing_volume)
    counts["n_dropped_parent_edu"] = int(len(dropped_parent_edu))
    counts["parcels_only"] = bool(parcels_only)
    return counts


def mgcv_method_name(requested: str) -> str:
    """Map the CLI name onto the string mgcv::gam accepts."""
    key = requested.strip().lower()
    if key in {"gcv", "gcv.cp"}:
        return "GCV.Cp"
    if key == "reml":
        return "REML"
    raise ValueError(f"mgcv method must be gcv or reml, got {requested}")


def mgcv_basis_k(requested: int) -> int:
    """k is the smooth basis size. mgcv needs at least 3 for cr and ps."""
    if int(requested) < 3:
        raise ValueError(f"basis k must be an integer >= 3, got {requested}")
    return int(requested)


def mgcv_basis_name(requested: str) -> str:
    """Map the CLI name onto mgcv's bs argument. ps is a P-spline."""
    key = requested.strip().lower()
    if key in {"cr", "ps"}:
        return key
    raise ValueError(f"mgcv basis must be cr or ps, got {requested}")


def fit_models(rscript: Path, input_dir: Path, result_dir: Path, method: str, basis_k: int,
               basis: str) -> None:
    result_dir.mkdir(parents=True, exist_ok=True)
    completed = subprocess.run(
        [str(rscript), str(R_SCRIPT), str(input_dir), str(result_dir), method, str(basis_k), basis],
        check=False,
    )
    if completed.returncode != 0:
        raise RuntimeError(f"mgcv nested models failed with exit code {completed.returncode}")


def fit_pair_slopes(rscript: Path, input_dir: Path, result_dir: Path, method: str, basis_k: int,
                    basis: str) -> None:
    if not (input_dir / "design.csv").is_file() or not (input_dir / "group_response.csv").is_file():
        raise FileNotFoundError(
            f"{input_dir} does not contain the saved design and group means. "
            "Run the nested models once before the slope comparison."
        )
    result_dir.mkdir(parents=True, exist_ok=True)
    completed = subprocess.run(
        [str(rscript), str(SLOPE_SCRIPT), str(input_dir), str(result_dir), method, str(basis_k), basis],
        check=False,
    )
    if completed.returncode != 0:
        raise RuntimeError(f"m3 pairwise slopes failed with exit code {completed.returncode}")


def write_metadata(output_dir: Path, counts: dict, method: str, basis_k: int, basis: str) -> None:
    payload = {
        "models": {
            "m1": "y ~ gender + site_id + parent_edu + meanFD",
            "m2": "y ~ age + gender + site_id + parent_edu + meanFD",
            "m3": f"y ~ s(age, k={basis_k}, bs={basis}) + gender + site_id + parent_edu + meanFD",
        },
        "method": method,
        "basis_k": basis_k,
        "basis": basis,
        "comparison": ["m2 vs m1", "m3 vs m1", "m3 vs m2"],
        "bh": "within parcel only" if counts.get("parcels_only") else "within parcel, yeo, system, and amygdala, separately for each comparison",
        "parcels_only": bool(counts.get("parcels_only")),
        "contrast": "Reappraise minus Look negative",
        "meanfd_rule": "drop subjects with meanFD greater than meanfd_max" if "meanfd_max" in counts else "no meanFD magnitude filter",
        "counts": counts,
    }
    (output_dir / "run_metadata.json").write_text(json.dumps(payload, indent=2), encoding="utf-8")


def meanfd_directory_suffix(threshold: float) -> str:
    """0.5 becomes _FD05 and 0.2 becomes _FD02."""
    scaled = float(threshold) * 10.0
    if threshold <= 0 or abs(scaled - round(scaled)) > 1e-6:
        raise ValueError(
            f"meanFD threshold must be a positive tenth of a millimeter, got {threshold}"
        )
    return f"_FD{int(round(scaled)):02d}"


def subset_saved_inputs(source_dir: Path, dest_dir: Path, max_meanfd: float) -> dict[str, int]:
    """Copy saved volume means and keep subjects with meanFD at or below the threshold."""
    if source_dir.resolve() == dest_dir.resolve():
        raise ValueError("The meanFD subset cannot overwrite the shared inputs")
    if dest_dir.exists():
        shutil.rmtree(dest_dir)
    shutil.copytree(source_dir, dest_dir)
    design = pd.read_csv(dest_dir / "design.csv")
    design["sub_id"] = design["sub_id"].astype(str)
    updated, dropped = apply_meanfd_exclusion(design, max_meanfd)
    kept = analysis_design(updated)
    if len(kept) < 12:
        raise ValueError(f"meanFD <= {max_meanfd:g} leaves {len(kept)} subjects")
    updated.to_csv(dest_dir / "design.csv", index=False)
    dropped.to_csv(dest_dir / "dropped_meanfd.csv", index=False)
    ids = kept.loc[:, ["sub_id"]]
    for name in ("parcel_response.csv", "amygdala_response.csv"):
        if name == "amygdala_response.csv" and not (dest_dir / name).is_file():
            continue
        table = pd.read_csv(dest_dir / name)
        table["sub_id"] = table["sub_id"].astype(str)
        ordered = ids.merge(table, on="sub_id", how="left", validate="one_to_one")
        if ordered.drop(columns=["sub_id"]).isna().any().any():
            raise ValueError(f"{name} lost subjects while applying the meanFD filter")
        ordered.to_csv(dest_dir / name, index=False)
    groups = pd.read_csv(dest_dir / "group_response.csv")
    groups["sub_id"] = groups["sub_id"].astype(str)
    groups = groups[groups["sub_id"].isin(set(ids["sub_id"].astype(str)))].copy()
    groups.to_csv(dest_dir / "group_response.csv", index=False)
    note = dest_dir / "meanfd_source.txt"
    previous = note.read_text(encoding="utf-8").rstrip() if note.is_file() else ""
    note.write_text(
        previous + "\n" + "\n".join([
            f"max_meanFD: {max_meanfd:g}",
            "rule: drop subjects with meanFD greater than max_meanFD",
            f"n_dropped_meanfd: {len(dropped)}",
            f"n_kept: {len(kept)}",
        ]) + "\n",
        encoding="utf-8",
    )
    return {
        "n_subjects": int(len(kept)),
        "n_dropped_meanfd": int(len(dropped)),
        "meanfd_max": float(max_meanfd),
    }


def write_followup_figures(output_dir: Path, dlabel_path: Path) -> None:
    """m3-versus-m1 parcel curves and the age versus parent_edu scatter, both SVG."""
    input_dir = output_dir / "inputs"
    result_dir = output_dir / "statistics"
    figure_dir = output_dir / "figures"
    statistics = pd.read_csv(result_dir / "nested_model_comparison.csv")
    curves = pd.read_csv(result_dir / "nested_curves.csv")
    design = analysis_design(pd.read_csv(input_dir / "design.csv"))
    parcel_info = pd.read_csv(input_dir / "parcel_info.csv")
    association = plot_age_parent_edu(design, figure_dir / "age_parent_edu.svg")
    association.to_csv(result_dir / "age_parent_edu_correlation.csv", index=False)
    selected = plot_m3_m1_parcel_trajectories(
        statistics, curves, parcel_info, dlabel_path,
        figure_dir / "significant_m3_vs_m1_parcels.svg",
    )
    columns = [column for column in (
        "parcel", "parcel_name", "hemi", "network", "q_m3_vs_m1", "p_m3_vs_m1", "F_m3_vs_m1",
    ) if column in selected.columns]
    selected.loc[:, columns].to_csv(result_dir / "significant_m3_vs_m1_parcels.csv", index=False)
    row = association.iloc[0]
    print(
        f"age 与 parent_edu：n={int(row['n'])}，"
        f"Pearson r={row['pearson_r']:.4f}，p={row['pearson_p']:.4g}，"
        f"Spearman r={row['spearman_r']:.4f}，p={row['spearman_p']:.4g}"
    )
    print(f"Parcels with q(m3 vs m1) < 0.05: {len(selected)}")
    print(f"Association figure: {figure_dir / 'age_parent_edu.svg'}")
    if selected.empty:
        print("No parcel m3-versus-m1 figure was written.")
    else:
        print(f"Parcel trajectories: {figure_dir / 'significant_m3_vs_m1_parcels.svg'}")


def main() -> None:
    parser = argparse.ArgumentParser(description="Nested m1, m2, and m3 mgcv age models")
    parser.add_argument("--input", type=Path, default=None)
    parser.add_argument("--subjects", type=Path, default=ORIGINAL_SUBJECTS)
    parser.add_argument("--excluded", type=Path, default=ORIGINAL_EXCLUDED)
    parser.add_argument("--parent-edu", type=Path, default=ORIGINAL_PARENT_EDU)
    parser.add_argument("--meanfd", type=Path, default=ORIGINAL_MEANFD)
    parser.add_argument("--fig3-volume-dir", type=Path, default=ORIGINAL_VOLUME_DIR)
    parser.add_argument("--fig3-volume-atlas", type=Path, default=DEFAULT_ATLAS)
    parser.add_argument("--schaefer-atlas", type=Path, default=DEFAULT_DLABEL)
    parser.add_argument("--system-mask-dir", type=Path, default=ORIGINAL_SYSTEM_MASKS)
    parser.add_argument("--juelich-atlas", type=Path, default=None,
                        help="Juelich-maxprob-thr25-2mm.nii.gz. Defaults to FSLDIR, then nilearn.")
    parser.add_argument("--juelich-labels", type=Path, default=None,
                        help="FSL Juelich.xml with the integer labels for that NIfTI.")
    parser.add_argument("--parcels-only", action="store_true",
                        help="Fit parcel means only. Required for the Kong 2022 17-network 400-parcel atlas.")
    parser.add_argument("--parcel-info", type=Path, default=None,
                        help="Schaefer order_info.txt. Required with --parcels-only when it is not beside the volume atlas.")
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--meanfd-max", type=float, nargs="+", default=None,
                        help="Drop subjects with meanFD above each value. 0.5 0.2 writes _FD05 and _FD02.")
    parser.add_argument("--rscript", type=Path, default=None)
    parser.add_argument("--method", choices=("gcv", "reml"), default="reml",
                        help="mgcv smoothing-parameter method. gcv is GCV.Cp. Default is REML.")
    parser.add_argument("--basis-k", type=int, default=3,
                        help="Basis dimension k of s(age). Default 3. Must be >= 3.")
    parser.add_argument("--basis", choices=("cr", "ps"), default="cr",
                        help="mgcv smooth basis. cr is a cubic regression spline. ps is a P-spline. Default is cr.")
    parser.add_argument("--step", choices=("all", "slopes"), default="all",
                        help="slopes compares Yeo and system m3 slopes from an existing output directory.")
    args = parser.parse_args()
    smooth_method = mgcv_method_name(args.method)
    basis_k = mgcv_basis_k(args.basis_k)
    smooth_basis = mgcv_basis_name(args.basis)
    output_dir = resolve_project_path(args.output)
    if args.meanfd_max:
        output_dirs = [
            Path(str(output_dir) + meanfd_directory_suffix(threshold))
            for threshold in args.meanfd_max
        ]
    else:
        output_dirs = [output_dir]
    rscript = find_rscript(args.rscript)
    dlabel_path = resolve_project_path(args.schaefer_atlas)
    volume_atlas = resolve_project_path(args.fig3_volume_atlas)
    parcel_info_path = (
        parcel_info_for(volume_atlas, None if args.parcel_info is None else resolve_project_path(args.parcel_info))
        if args.parcels_only else None
    )
    if args.step == "slopes":
        if args.parcels_only:
            raise RuntimeError(
                "--parcels-only fits parcel means only. The slope step needs Yeo networks and the four systems."
            )
        for one_output in output_dirs:
            group_path = one_output / "inputs" / "group_response.csv"
            if not group_path.is_file() or pd.read_csv(group_path).empty:
                raise RuntimeError(
                    f"{one_output} has no network or system means. "
                    "The slope step does not run on a parcel-only fit."
                )
            input_dir = one_output / "inputs"
            result_dir = one_output / "statistics"
            fit_pair_slopes(rscript, input_dir, result_dir, smooth_method, basis_k, smooth_basis)
            pair_stats = pd.read_csv(result_dir / "m3_pair_slopes.csv")
            pair_curves = pd.read_csv(result_dir / "m3_pair_curves.csv")
            matrices = write_slope_matrices(pair_stats, result_dir)
            figures = plot_m3_slope_matrices(pair_stats, pair_curves, one_output / "figures")
            for dimension in ("yeo", "system"):
                subset = pair_stats[pair_stats["dimension"] == dimension]
                n_significant = int((subset["q_slope"] < 0.05).sum())
                print(f"{one_output.name} {dimension}: {len(subset)} ordered pairs, {n_significant} with q_slope < 0.05")
            print(f"Slope table: {result_dir / 'm3_pair_slopes.csv'}")
            print(f"Matrices: {', '.join(path.name for path in matrices)}")
            print(f"Figures: {', '.join(path.name for path in figures)}")
        return
    if args.meanfd_max:
        shared_inputs = Path(str(output_dir) + "_shared_inputs") / "inputs"
        stamp = shared_inputs.parent / "atlas.txt"
        stamp_lines = [str(volume_atlas.resolve()), str(dlabel_path.resolve())]
        if args.parcels_only:
            stamp_lines.extend(["parcels_only", str(parcel_info_path.resolve())])
        stamp_text = "\n".join(stamp_lines) + "\n"
        saved_means = shared_inputs / "parcel_response.csv"
        if saved_means.is_file() and stamp.is_file() and stamp.read_text(encoding="utf-8") == stamp_text:
            print(f"Using saved volume means: {shared_inputs}")
            ensure_subject_roster(
                shared_inputs,
                ORIGINAL_DEMOGRAPHICS if args.input is None else resolve_project_path(args.input),
                resolve_project_path(args.subjects),
                resolve_project_path(args.excluded),
                resolve_project_path(args.parent_edu),
                resolve_project_path(args.meanfd),
                resolve_project_path(args.fig3_volume_dir),
            )
            shared_counts = {
                "n_subjects": n_analysis_subjects(pd.read_csv(shared_inputs / "design.csv")),
                "input": "saved volume means",
                "parcels_only": bool(args.parcels_only),
            }
        else:
            shared_counts = prepare_real_inputs(
                None if args.input is None else resolve_project_path(args.input),
                resolve_project_path(args.subjects),
                resolve_project_path(args.excluded),
                resolve_project_path(args.parent_edu),
                resolve_project_path(args.meanfd),
                resolve_project_path(args.fig3_volume_dir),
                volume_atlas,
                dlabel_path,
                resolve_project_path(args.system_mask_dir),
                None if args.juelich_atlas is None else resolve_project_path(args.juelich_atlas),
                None if args.juelich_labels is None else resolve_project_path(args.juelich_labels),
                shared_inputs,
                parcels_only=args.parcels_only,
                parcel_info_path=parcel_info_path,
            )
            stamp.parent.mkdir(parents=True, exist_ok=True)
            stamp.write_text(stamp_text, encoding="utf-8")
        prepared = []
        for threshold, one_output in zip(args.meanfd_max, output_dirs):
            filtered = subset_saved_inputs(shared_inputs, one_output / "inputs", threshold)
            counts = {**shared_counts, **filtered}
            print(
                f"meanFD <= {threshold:g}: {counts['n_subjects']} subjects, "
                f"dropped {counts['n_dropped_meanfd']} -> {one_output}"
            )
            prepared.append((one_output, counts))
    else:
        counts = prepare_real_inputs(
            None if args.input is None else resolve_project_path(args.input),
            resolve_project_path(args.subjects),
            resolve_project_path(args.excluded),
            resolve_project_path(args.parent_edu),
            resolve_project_path(args.meanfd),
            resolve_project_path(args.fig3_volume_dir),
            volume_atlas,
            dlabel_path,
            resolve_project_path(args.system_mask_dir),
            None if args.juelich_atlas is None else resolve_project_path(args.juelich_atlas),
            None if args.juelich_labels is None else resolve_project_path(args.juelich_labels),
            output_dir / "inputs",
            parcels_only=args.parcels_only,
            parcel_info_path=parcel_info_path,
        )
        prepared = [(output_dir, counts)]
    for one_output, counts in prepared:
        input_dir = one_output / "inputs"
        result_dir = one_output / "statistics"
        fit_models(rscript, input_dir, result_dir, smooth_method, basis_k, smooth_basis)
        statistics = pd.read_csv(result_dir / "nested_model_comparison.csv")
        curves = pd.read_csv(result_dir / "nested_curves.csv")
        if (statistics["dimension"].astype(str) == "yeo").any():
            figures = plot_curve_panels(statistics, curves, one_output / "figures")
        else:
            figures = []
        write_metadata(one_output, counts, smooth_method, basis_k, smooth_basis)
        write_followup_figures(one_output, dlabel_path)
        for dimension in ("parcel", "yeo", "system", "amygdala"):
            subset = statistics[statistics["dimension"].astype(str) == dimension]
            if subset.empty:
                continue
            n_nonlinear = int((subset["q_m3_vs_m2"] < 0.05).sum())
            n_any_age = int((subset["q_m3_vs_m1"] < 0.05).sum())
            print(
                f"{one_output.name} {dimension}: {len(subset)} models, "
                f"{n_any_age} with q(m3 vs m1) < 0.05, "
                f"{n_nonlinear} with q(m3 vs m2) < 0.05"
            )
        print(f"Comparison table: {result_dir / 'nested_model_comparison.csv'}")
        if figures:
            print(f"Figures: {', '.join(path.name for path in figures)}")


if __name__ == "__main__":
    main()
