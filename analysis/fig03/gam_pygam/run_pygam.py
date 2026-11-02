"""Historical pygam age test on volume Reappraise minus Look.

The full model is intercept + s(age) + gender + site_id + parent_edu + meanFD.
s(age) is a penalized B-spline with n_splines=5. The smoothing parameter is
chosen by GCV on 12 values from 1e-5 to 1e5. The reduced model is intercept
plus the same four covariates and is nested in the full model. p_age is the
nested F test of that pair. It asks whether the age smooth reduces residual
error, and it does not separate a straight line from a bend.

Subject volumes stay on their own grid. The Schaefer atlas is nearest-neighbor
labeled onto that grid once. This command does not run the vertex GLM in
original_code/brain_act_uni/brain_act_er_age.py.

The 200-parcel run also fits Yeo 7 networks, the four systems, and six amygdala ROIs. Pass
--parcels-only to fit only parcel means from the Kong 2022 400-parcel atlas.
That path reads order_info.txt for names and does not fit networks or systems.

Run from the repository root, in the environment that can import pygam:

    python -m analysis.fig03.gam_pygam.run_pygam --input "/public/home/dingrui/fmri_analysis/data/beh/participants_in_tfmri_demographics.csv" --subjects "/public/home/dingrui/fmri_analysis/data/beh/sub_list_valid_4_ER.txt" --excluded "/public/home/dingrui/fmri_analysis/data/fmri/subs_removed_from_2ndlevel.txt" --parent-edu "/public/home/dingrui/fmri_analysis/data/beh/env_family_student_cbcl.tsv" --meanfd "/public/home/dingrui/fmri_analysis/data/fmri/qc_res_subs_runs_fmri.csv" --fig3-volume-dir "/public/home/dingrui/fmri_analysis/data/fmri/res_1stlevel" --fig3-volume-atlas "resources/atlases/schaefer2018/Schaefer2018_200Parcels_7Networks_order_FSLMNI152_2mm.nii.gz" --schaefer-atlas "resources/atlases/schaefer2018/Schaefer2018_200Parcels_7Networks_order.dlabel.nii" --system-mask-dir "/public/home/dingrui/fmri_analysis/res/res_bayes" --output "outputs/fig03/gam_pygam_volume"

Schaefer 200 x Yeo 7 with no meanFD magnitude filter. The age smooth stays
5 B-splines and the 12-point GCV grid. The result is
outputs/fig03/gam_pygam_volume_200. This does not write gam_pygam_volume or
the _FD05, _FD02, and _shared_inputs directories.

    python -m analysis.fig03.gam_pygam.run_pygam --input "/public/home/dingrui/fmri_analysis/data/beh/participants_in_tfmri_demographics.csv" --subjects "/public/home/dingrui/fmri_analysis/data/beh/sub_list_valid_4_ER.txt" --excluded "/public/home/dingrui/fmri_analysis/data/fmri/subs_removed_from_2ndlevel.txt" --parent-edu "/public/home/dingrui/fmri_analysis/data/beh/env_family_student_cbcl.tsv" --meanfd "/public/home/dingrui/fmri_analysis/data/fmri/qc_res_subs_runs_fmri.csv" --fig3-volume-dir "/public/home/dingrui/fmri_analysis/data/fmri/res_1stlevel" --fig3-volume-atlas "resources/atlases/schaefer2018/Schaefer2018_200Parcels_7Networks_order_FSLMNI152_2mm.nii.gz" --schaefer-atlas "resources/atlases/schaefer2018/Schaefer2018_200Parcels_7Networks_order.dlabel.nii" --system-mask-dir "/public/home/dingrui/fmri_analysis/res/res_bayes" --output "outputs/fig03/gam_pygam_volume_200"

Schaefer 200 x Yeo 7 with meanFD ceilings 0.5 and 0.2. The age smooth stays
5 B-splines and the 12-point GCV grid. This does not write gam_pygam_volume.

    python -m analysis.fig03.gam_pygam.run_pygam --meanfd-max 0.5 0.2 --input "/public/home/dingrui/fmri_analysis/data/beh/participants_in_tfmri_demographics.csv" --subjects "/public/home/dingrui/fmri_analysis/data/beh/sub_list_valid_4_ER.txt" --excluded "/public/home/dingrui/fmri_analysis/data/fmri/subs_removed_from_2ndlevel.txt" --parent-edu "/public/home/dingrui/fmri_analysis/data/beh/env_family_student_cbcl.tsv" --meanfd "/public/home/dingrui/fmri_analysis/data/fmri/qc_res_subs_runs_fmri.csv" --fig3-volume-dir "/public/home/dingrui/fmri_analysis/data/fmri/res_1stlevel" --fig3-volume-atlas "resources/atlases/schaefer2018/Schaefer2018_200Parcels_7Networks_order_FSLMNI152_2mm.nii.gz" --schaefer-atlas "resources/atlases/schaefer2018/Schaefer2018_200Parcels_7Networks_order.dlabel.nii" --system-mask-dir "/public/home/dingrui/fmri_analysis/res/res_bayes" --output "outputs/fig03/gam_pygam_volume_200"

Schaefer 400 x Yeo 7 with no meanFD magnitude filter. The age smooth stays
5 B-splines and the 12-point GCV grid. The result is
outputs/fig03/gam_pygam_volume_400_7net. This does not write gam_pygam_volume,
gam_pygam_volume_400, or the _FD05, _FD02, and _shared_inputs directories.

    python -m analysis.fig03.gam_pygam.run_pygam --input "/public/home/dingrui/fmri_analysis/data/beh/participants_in_tfmri_demographics.csv" --subjects "/public/home/dingrui/fmri_analysis/data/beh/sub_list_valid_4_ER.txt" --excluded "/public/home/dingrui/fmri_analysis/data/fmri/subs_removed_from_2ndlevel.txt" --parent-edu "/public/home/dingrui/fmri_analysis/data/beh/env_family_student_cbcl.tsv" --meanfd "/public/home/dingrui/fmri_analysis/data/fmri/qc_res_subs_runs_fmri.csv" --fig3-volume-dir "/public/home/dingrui/fmri_analysis/data/fmri/res_1stlevel" --fig3-volume-atlas "resources/atlases/schaefer2018/Schaefer2018_400Parcels_7Networks_order_FSLMNI152_2mm.nii.gz" --schaefer-atlas "resources/atlases/schaefer2018/Schaefer2018_400Parcels_7Networks_order.dlabel.nii" --system-mask-dir "/public/home/dingrui/fmri_analysis/res/res_bayes" --output "outputs/fig03/gam_pygam_volume_400_7net"

Schaefer 400 x Yeo 7, four systems, and six amygdala ROIs. Volumes are read
once. meanFD above 0.5 and above 0.2 are dropped into separate directories.
Figures are SVG with editable text.

    python -m analysis.fig03.gam_pygam.run_pygam --meanfd-max 0.5 0.2 --input "/public/home/dingrui/fmri_analysis/data/beh/participants_in_tfmri_demographics.csv" --subjects "/public/home/dingrui/fmri_analysis/data/beh/sub_list_valid_4_ER.txt" --excluded "/public/home/dingrui/fmri_analysis/data/fmri/subs_removed_from_2ndlevel.txt" --parent-edu "/public/home/dingrui/fmri_analysis/data/beh/env_family_student_cbcl.tsv" --meanfd "/public/home/dingrui/fmri_analysis/data/fmri/qc_res_subs_runs_fmri.csv" --fig3-volume-dir "/public/home/dingrui/fmri_analysis/data/fmri/res_1stlevel" --fig3-volume-atlas "resources/atlases/schaefer2018/Schaefer2018_400Parcels_7Networks_order_FSLMNI152_2mm.nii.gz" --schaefer-atlas "resources/atlases/schaefer2018/Schaefer2018_400Parcels_7Networks_order.dlabel.nii" --system-mask-dir "/public/home/dingrui/fmri_analysis/res/res_bayes" --output "outputs/fig03/gam_pygam_volume_400_7net"

Schaefer 400 plus Tian S2, and Schaefer 200 plus Tian S2. Both keep Yeo 7 and
the four systems. Tian labels 1..32 are appended, lamy and mamy supply the
four amygdala parcels inside the parcel BH family, and only meanFD above 0.5
is dropped. See the README for the commands and the map and curve outputs.

400 parcels, Kong 2022 parcel means only:

    python -m analysis.fig03.gam_pygam.run_pygam --parcels-only --input "/public/home/dingrui/fmri_analysis/data/beh/participants_in_tfmri_demographics.csv" --subjects "/public/home/dingrui/fmri_analysis/data/beh/sub_list_valid_4_ER.txt" --excluded "/public/home/dingrui/fmri_analysis/data/fmri/subs_removed_from_2ndlevel.txt" --parent-edu "/public/home/dingrui/fmri_analysis/data/beh/env_family_student_cbcl.tsv" --meanfd "/public/home/dingrui/fmri_analysis/data/fmri/qc_res_subs_runs_fmri.csv" --fig3-volume-dir "/public/home/dingrui/fmri_analysis/data/fmri/res_1stlevel" --fig3-volume-atlas "resources/atlases/schaefer2018/Schaefer2018_400Parcels_Kong2022_17Networks_order_FSLMNI152_2mm.nii.gz" --parcel-info "resources/atlases/schaefer2018/Schaefer2018_400Parcels_Kong2022_17Networks_order_info.txt" --output "outputs/fig03/gam_pygam_volume_400"
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

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
    SYSTEM_ORDER,
    YEO_ORDER,
    analysis_design,
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
    schaefer_parcel_info,
)

PROJECT_ROOT = Path(__file__).resolve().parents[3]
MYLIB_ROOT = PROJECT_ROOT / "src" / "prj_analysis" / "src"
DEFAULT_ATLAS = PROJECT_ROOT / "resources/atlases/schaefer2018/Schaefer2018_200Parcels_7Networks_order_FSLMNI152_2mm.nii.gz"
DEFAULT_DLABEL = PROJECT_ROOT / "resources/atlases/schaefer2018/Schaefer2018_200Parcels_7Networks_order.dlabel.nii"
DEFAULT_ATLAS_400 = PROJECT_ROOT / "resources/atlases/schaefer2018/Schaefer2018_400Parcels_Kong2022_17Networks_order_FSLMNI152_2mm.nii.gz"
DEFAULT_PARCEL_INFO_400 = PROJECT_ROOT / "resources/atlases/schaefer2018/Schaefer2018_400Parcels_Kong2022_17Networks_order_info.txt"
DEFAULT_OUTPUT = PROJECT_ROOT / "outputs" / "fig03" / "gam_pygam_volume"
COVARIATES = ("gender", "site_id", "parent_edu", "meanFD")
K_SPLINES = 5
LAM_LOG_RANGE = (-5.0, 5.0, 12)
GRID_SIZE = 101


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


def load_historical_stats():
    """Import the pygam helpers from this repository, ahead of any installed copy."""
    stats_path = MYLIB_ROOT / "mylib" / "stat" / "stats.py"
    if not stats_path.is_file():
        raise FileNotFoundError(
            "Historical GAM helpers were not found at "
            f"{stats_path}. Copy src/prj_analysis/src/mylib/stat/stats.py there."
        )
    if str(MYLIB_ROOT) not in sys.path:
        sys.path.insert(0, str(MYLIB_ROOT))
    try:
        from mylib.stat import stats as stat_utils
    except ModuleNotFoundError as exc:
        raise ModuleNotFoundError(
            "The historical GAM needs pygam. Activate the er-manuscript environment, then rerun this command."
        ) from exc
    return stat_utils


def _parcel_info_path(volume_atlas: Path, requested: Path | None) -> Path:
    if requested is not None:
        return requested
    sibling = Path(str(volume_atlas).replace("_FSLMNI152_2mm.nii.gz", "_info.txt"))
    if sibling.is_file():
        return sibling
    if DEFAULT_PARCEL_INFO_400.is_file():
        return DEFAULT_PARCEL_INFO_400
    raise FileNotFoundError(
        "Parcel-only analysis needs --parcel-info pointing at the Schaefer order_info.txt"
    )


def _require_tian_columns(input_dir: Path, n_cortical: int) -> None:
    columns = set(pd.read_csv(input_dir / "parcel_response.csv", nrows=0).columns)
    missing = [str(n_cortical + label) for label in range(1, 33) if str(n_cortical + label) not in columns]
    if missing:
        raise ValueError(
            "Tian labels did not survive resampling onto the effect-size grid: "
            + ", ".join(missing)
        )


def prepare_volume_inputs(input_path: Path | None, subjects_path: Path, excluded_path: Path,
                          parent_edu_path: Path, meanfd_path: Path, volume_dir: Path,
                          volume_atlas: Path, dlabel_path: Path, system_mask_dir: Path,
                          input_dir: Path, parcels_only: bool = False,
                          parcel_info_path: Path | None = None,
                          juelich_atlas: Path | None = None,
                          juelich_labels: Path | None = None,
                          tian_atlas: Path | None = None) -> dict[str, int]:
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
    extra_parcel_info = None
    atlas_for_volumes = volume_atlas
    min_extra_voxels = 5
    tian_note = None
    if parcels_only and tian_atlas is not None:
        raise ValueError("--tian-atlas cannot be combined with --parcels-only")
    if parcels_only:
        info_path = _parcel_info_path(volume_atlas, parcel_info_path)
        parcel_table = read_schaefer_order_info(info_path)
        print(f"Parcel-only atlas names: {len(parcel_table)} rows from {info_path}")
    elif tian_atlas is not None:
        from analysis.fig03.gam_pygam.tian_combined import build_combined_atlas, marginal_parcels_for

        input_dir.mkdir(parents=True, exist_ok=True)
        cortical_info = schaefer_parcel_info(dlabel_path)
        n_cortical = int(cortical_info["parcel"].nunique())
        if n_cortical not in (200, 400):
            raise ValueError(f"Tian combined atlas expects 200 or 400 Schaefer parcels, got {n_cortical}")
        cortical_ids = set(cortical_info["parcel"].astype(int))
        marginal = marginal_parcels_for(n_cortical, cortical_ids)
        combined_path = input_dir / "combined_atlas.nii.gz"
        extra_parcel_info, amygdala_table = build_combined_atlas(
            volume_atlas, tian_atlas, n_cortical, combined_path)
        atlas_for_volumes = combined_path
        tian_note = {
            "n_cortical": n_cortical,
            "n_parcels": n_cortical + 32,
            "marginal_parcels": marginal,
            "amygdala_parcels": amygdala_table["parcel"].astype(int).tolist(),
            "combined_atlas": combined_path.name,
            "tian_atlas": str(tian_atlas),
            "partial_r2": "max(0, 1 - RSS_full / RSS_reduced)",
        }
        amygdala_table.to_csv(input_dir / "amygdala_parcels.csv", index=False)
        (input_dir / "tian_source.txt").write_text(
            "\n".join([
                f"tian_atlas: {tian_atlas}",
                f"schaefer_atlas: {volume_atlas}",
                f"n_cortical: {n_cortical}",
                "tian_labels: 1..32 shifted by n_cortical",
                "amygdala: lamy and mamy, left and right, inside the parcel family",
                "marginal_parcels: " + (", ".join(str(parcel) for parcel in marginal) or "none"),
            ]) + "\n",
            encoding="utf-8",
        )
    else:
        from analysis.fig03.gam_nested.amygdala import (
            AMYGDALA_ORDER,
            MIN_AMYGDALA_VOXELS,
            amygdala_masks,
            load_juelich_atlas,
        )
        atlas_image, labels, atlas_source = load_juelich_atlas(juelich_atlas, juelich_labels)
        extra_masks = amygdala_masks(atlas_image, labels)
        if tuple(extra_masks) != AMYGDALA_ORDER:
            raise RuntimeError(f"Amygdala ROI order is {tuple(extra_masks)}, expected {AMYGDALA_ORDER}")
        min_extra_voxels = MIN_AMYGDALA_VOXELS
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
        metadata, rpsl_paths, lkng_paths, atlas_for_volumes, dlabel_path, system_mask_dir, input_dir,
        extra_masks=extra_masks, min_extra_voxels=min_extra_voxels,
        parcel_info_table=parcel_table, include_groups=not parcels_only,
        extra_parcel_info=extra_parcel_info)
    counts["n_dropped_no_volume"] = len(missing_volume)
    counts["n_dropped_parent_edu"] = int(len(dropped_parent_edu))
    counts["parcels_only"] = bool(parcels_only)
    if tian_note is not None:
        _require_tian_columns(input_dir, int(tian_note["n_cortical"]))
        (input_dir / "tian_combined.json").write_text(json.dumps(tian_note, indent=2), encoding="utf-8")
        counts["tian_combined"] = True
        counts["n_cortical"] = int(tian_note["n_cortical"])
        counts["n_amygdala"] = len(tian_note["amygdala_parcels"])
        counts["marginal_parcels"] = list(tian_note["marginal_parcels"])
    return counts


def _age_grid(age: np.ndarray) -> np.ndarray:
    finite = age[np.isfinite(age)]
    quantiles = np.linspace(0.01, 0.99, GRID_SIZE)
    return np.unique(np.quantile(finite, quantiles))


def age_effect_on_residual_scale(model, y: np.ndarray, covariates: np.ndarray,
                                 finite: np.ndarray, grid: np.ndarray,
                                 covariate_mean: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Age curve on the scale of residuals after the intercept and covariates.

    The residual subtracts the ordinary-least-squares fit of the intercept and
    the four covariates. The curve is the full-model prediction at the mean of
    those covariates minus that reduced prediction. The band is the full model's
    95% pointwise confidence interval shifted by the same reduced prediction.
    """
    y = np.asarray(y, dtype=float)
    covariates = np.asarray(covariates, dtype=float)
    if covariates.ndim == 1:
        covariates = covariates.reshape(-1, 1)
    finite = np.asarray(finite, dtype=bool)
    grid = np.asarray(grid, dtype=float)
    covariate_mean = np.asarray(covariate_mean, dtype=float).reshape(-1)
    reduced_design = np.column_stack([np.ones(int(finite.sum())), covariates[finite]])
    beta = np.linalg.lstsq(reduced_design, y[finite], rcond=None)[0]
    residual = np.full(y.shape[0], np.nan, dtype=float)
    residual[finite] = y[finite] - reduced_design @ beta
    design = np.column_stack([
        grid,
        np.repeat(covariate_mean.reshape(1, -1), grid.size, axis=0),
    ])
    fitted = np.asarray(model.predict(design), dtype=float).reshape(-1)
    reduced_at_mean = float(beta[0] + np.dot(covariate_mean, beta[1:]))
    effect = fitted - reduced_at_mean
    intervals = np.asarray(model.confidence_intervals(design, width=0.95), dtype=float)
    if intervals.shape == (2, grid.size) and grid.size != 2:
        intervals = intervals.T
    if intervals.shape != (grid.size, 2):
        raise ValueError(f"confidence_intervals returned shape {intervals.shape}, expected {(grid.size, 2)}")
    return effect, intervals[:, 0] - reduced_at_mean, intervals[:, 1] - reduced_at_mean, residual


def _derivative(curve: np.ndarray, grid: np.ndarray) -> np.ndarray:
    slope = np.full(grid.size, np.nan, dtype=np.float64)
    if grid.size < 3:
        return slope
    slope[1:-1] = (curve[2:] - curve[:-2]) / (grid[2:] - grid[:-2])
    slope[0] = (curve[1] - curve[0]) / (grid[1] - grid[0])
    slope[-1] = (curve[-1] - curve[-2]) / (grid[-1] - grid[-2])
    return slope


def _fit_one(stat_utils, y: np.ndarray, age: np.ndarray, covariates: np.ndarray,
             grid: np.ndarray, lam_grid: np.ndarray,
             store_partial: bool = False) -> tuple[dict, np.ndarray, dict | None]:
    finite = np.isfinite(y) & np.isfinite(age) & np.all(np.isfinite(covariates), axis=1)
    n_obs = int(finite.sum())
    empty_curve = np.full(grid.size, np.nan)

    def _empty_partial() -> dict | None:
        if not store_partial:
            return None
        return {
            "effect": empty_curve.copy(),
            "ci_low": empty_curve.copy(),
            "ci_high": empty_curve.copy(),
            "residual": np.full(y.shape[0], np.nan, dtype=float),
        }

    if n_obs < max(12, K_SPLINES + 2):
        return {"n_obs": n_obs}, empty_curve, _empty_partial()
    rss_full, edf_full, rss_reduced, edf_reduced, model = stat_utils.fit_gam_full_and_reduced(
        y[finite],
        age[finite, None],
        covariates[finite],
        k_splines=K_SPLINES,
        lam_grid=lam_grid,
        max_iter=1000,
        verbose=False,
    )
    f_value, p_value, df1, df2, partial_r2 = stat_utils.nested_anova_gaussian(
        rss_full=rss_full,
        rss_reduced=rss_reduced,
        n_obs=n_obs,
        edf_full=edf_full,
        edf_reduced=edf_reduced,
    )
    covariate_mean = covariates[finite].mean(axis=0)
    design = np.column_stack([grid, np.repeat(covariate_mean.reshape(1, -1), grid.size, axis=0)])
    fitted = np.asarray(model.predict(design), dtype=np.float64).reshape(-1)
    centered = fitted - np.nanmean(fitted)
    slope = _derivative(centered, grid)
    finite_slope = np.isfinite(slope)
    age_of_max_change = float(grid[finite_slope][np.nanargmax(np.abs(slope[finite_slope]))]) if finite_slope.any() else np.nan
    partial = None
    if store_partial:
        effect, ci_low, ci_high, residual = age_effect_on_residual_scale(
            model, y, covariates, finite, grid, covariate_mean)
        partial = {
            "effect": effect,
            "ci_low": ci_low,
            "ci_high": ci_high,
            "residual": residual,
        }
    return {
        "n_obs": n_obs,
        "rss_full": rss_full,
        "rss_reduced": rss_reduced,
        "edf_full": edf_full,
        "edf_reduced": edf_reduced,
        "F_age": f_value,
        "df1": df1,
        "df2": df2,
        "p_age": p_value,
        "partial_r2": partial_r2,
        "age_of_max_change": age_of_max_change,
    }, centered, partial


def _correct(stat_utils, rows: list[dict]) -> pd.DataFrame:
    table = pd.DataFrame(rows)
    if "p_age" not in table.columns:
        table["p_age"] = np.nan
    table["q_age"] = np.nan
    table["reject_age_fdr"] = False
    finite = np.isfinite(table["p_age"].to_numpy(dtype=float))
    if finite.any():
        corrected = stat_utils.fdr_bh(table.loc[finite, "p_age"].to_numpy(dtype=float), alpha=0.05)
        table.loc[finite, "q_age"] = corrected["qvals"]
        table.loc[finite, "reject_age_fdr"] = corrected["reject"]
    return table


def fit_saved_means(stat_utils, input_dir: Path,
                    store_partial: bool = False) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame | None]:
    design = analysis_design(pd.read_csv(input_dir / "design.csv"))
    design["sub_id"] = design["sub_id"].astype(str)
    age = design["age"].to_numpy(dtype=float)
    missing = [name for name in COVARIATES if name not in design.columns]
    if missing:
        raise ValueError(f"design.csv is missing covariates: {', '.join(missing)}")
    covariates = design.loc[:, list(COVARIATES)].to_numpy(dtype=float)
    grid = _age_grid(age)
    lam_grid = 10.0 ** np.linspace(*LAM_LOG_RANGE)
    parcels = pd.read_csv(input_dir / "parcel_response.csv")
    parcels["sub_id"] = parcels["sub_id"].astype(str)
    parcels = design[["sub_id"]].merge(parcels, on="sub_id", how="left", validate="one_to_one")
    groups = pd.read_csv(input_dir / "group_response.csv")
    groups["sub_id"] = groups["sub_id"].astype(str)
    rows: list[dict] = []
    curves: list[pd.DataFrame] = []
    residual_columns: dict[str, np.ndarray] = {}

    def add_unit(dimension: str, name: str, values: np.ndarray, index: int, total: int) -> None:
        print(f"pygam {dimension} {name} ({index}/{total})", flush=True)
        statistics, centered, partial = _fit_one(
            stat_utils, values, age, covariates, grid, lam_grid,
            store_partial=store_partial and dimension == "parcel")
        statistics.update({"dimension": dimension, "name": str(name)})
        rows.append(statistics)
        curve = {
            "dimension": dimension,
            "name": str(name),
            "age": grid,
            "centered_effect": centered,
        }
        if partial is not None:
            curve["effect"] = partial["effect"]
            curve["ci_low"] = partial["ci_low"]
            curve["ci_high"] = partial["ci_high"]
            if dimension == "parcel":
                residual_columns[str(name)] = partial["residual"]
        curves.append(pd.DataFrame(curve))

    parcel_names = [column for column in parcels.columns if column != "sub_id"]
    parcel_names = sorted(parcel_names, key=lambda name: int(name))
    fitted_dimensions = ["parcel"]
    for index, name in enumerate(parcel_names, start=1):
        add_unit("parcel", name, parcels[name].to_numpy(dtype=float), index, len(parcel_names))
    for dimension, order in (("yeo", YEO_ORDER), ("system", SYSTEM_ORDER)):
        if "dimension" not in groups.columns or groups.empty:
            continue
        subset = groups[groups["dimension"] == dimension]
        if subset.empty:
            continue
        fitted_dimensions.append(dimension)
        for index, name in enumerate(order, start=1):
            matched = design[["sub_id"]].merge(
                subset.loc[subset["group"] == name, ["sub_id", "y"]],
                on="sub_id", how="left", validate="one_to_one")
            add_unit(dimension, name, matched["y"].to_numpy(dtype=float), index, len(order))
    amygdala_path = input_dir / "amygdala_response.csv"
    if amygdala_path.is_file():
        amygdala = pd.read_csv(amygdala_path)
        amygdala["sub_id"] = amygdala["sub_id"].astype(str)
        order_path = input_dir / "amygdala_order.csv"
        if order_path.is_file():
            amygdala_names = pd.read_csv(order_path)["roi"].astype(str).tolist()
        else:
            amygdala_names = [column for column in amygdala.columns if column != "sub_id"]
        missing_rois = [name for name in amygdala_names if name not in amygdala.columns]
        if missing_rois:
            raise ValueError(f"amygdala_response.csv is missing: {', '.join(missing_rois)}")
        fitted_dimensions.append("amygdala")
        for index, name in enumerate(amygdala_names, start=1):
            matched = design[["sub_id"]].merge(
                amygdala.loc[:, ["sub_id", name]],
                on="sub_id", how="left", validate="one_to_one")
            add_unit("amygdala", name, matched[name].to_numpy(dtype=float), index, len(amygdala_names))
    statistics = pd.concat([_correct(stat_utils, [row for row in rows if row["dimension"] == dimension])
                            for dimension in fitted_dimensions], ignore_index=True)
    column_order = [
        "dimension", "name", "n_obs", "F_age", "df1", "df2", "p_age", "q_age",
        "reject_age_fdr", "partial_r2", "rss_full", "rss_reduced",
        "edf_full", "edf_reduced", "age_of_max_change",
    ]
    residuals = None
    if store_partial:
        residuals = pd.DataFrame({"sub_id": design["sub_id"].to_numpy(), "age": age})
        for name in sorted(residual_columns, key=lambda value: int(value)):
            residuals[name] = residual_columns[name]
    return statistics.reindex(columns=column_order), pd.concat(curves, ignore_index=True), residuals


def write_metadata(output_dir: Path, counts: dict) -> None:
    covariate_terms = " + ".join(COVARIATES)
    payload = {
        "full_model": f"intercept + s(age, n_splines=5) + {covariate_terms}",
        "reduced_model": f"intercept + {covariate_terms}",
        "smoothing": "GCV grid from 1e-5 to 1e5, 12 values",
        "comparison": "nested F of the full model against the reduced model; p_age",
        "curve": (
            "residual after the intercept and four covariates; "
            "the age curve is the full model at the mean covariates minus that reduced prediction, "
            "with the full model's 95% pointwise confidence interval"
            if counts.get("tian_combined") else
            "covariates at their means; the 1% to 99% age-quantile curve is then mean-centered"
        ),
        "bh": (
            "within the combined parcels, within yeo, and within system"
            if counts.get("tian_combined") else
            "within parcel" if counts.get("parcels_only") else
            "within parcel, within yeo, within system, and within amygdala"
        ),
        "contrast": "effectSize-reg_neg_Observe_across_runs minus effectSize-pas_neg_Observe_across_runs",
        "meanfd_rule": "drop subjects with meanFD greater than meanfd_max" if "meanfd_max" in counts else "no meanFD magnitude filter",
        "counts": counts,
    }
    (output_dir / "run_metadata.json").write_text(json.dumps(payload, indent=2), encoding="utf-8")


def saved_counts(input_dir: Path) -> dict:
    """Describe an inputs directory that was written by an earlier prepare step."""
    design = pd.read_csv(input_dir / "design.csv")
    groups = pd.read_csv(input_dir / "group_response.csv")
    has_groups = (not groups.empty) and {"yeo", "system"}.intersection(set(groups["dimension"].astype(str)))
    counts = {
        "n_subjects": n_analysis_subjects(design),
        "input": "saved volume means",
        "parcels_only": not bool(has_groups),
        "n_amygdala": 6 if (input_dir / "amygdala_response.csv").is_file() else 0,
    }
    note_path = input_dir / "tian_combined.json"
    if note_path.is_file():
        note = json.loads(note_path.read_text(encoding="utf-8"))
        counts["tian_combined"] = True
        counts["n_cortical"] = int(note["n_cortical"])
        counts["n_amygdala"] = len(note.get("amygdala_parcels", []))
        counts["marginal_parcels"] = list(note.get("marginal_parcels", []))
        counts["parcels_only"] = False
    return counts


def _write_tian_followup(output_dir: Path, statistics: pd.DataFrame, curves: pd.DataFrame,
                         dlabel_path: Path) -> None:
    """Per-parcel residual curves and the partial-R² / mask brain files."""
    from analysis.fig03.gam_pygam.brain_maps import write_brain_maps
    from analysis.fig03.gam_pygam.plot_region_age import write_region_curves

    input_dir = output_dir / "inputs"
    result_dir = output_dir / "statistics"
    figure_dir = output_dir / "figures"
    note = json.loads((input_dir / "tian_combined.json").read_text(encoding="utf-8"))
    residuals = pd.read_csv(result_dir / "covariate_residuals.csv")
    parcel_info = pd.read_csv(input_dir / "parcel_info.csv")
    selected = write_region_curves(
        statistics, curves, residuals, parcel_info,
        list(note.get("marginal_parcels", [])),
        list(note.get("amygdala_parcels", [])),
        figure_dir,
    )
    selected.to_csv(result_dir / "trajectory_parcels.csv", index=False)
    write_brain_maps(
        input_dir / str(note["combined_atlas"]),
        statistics,
        list(note.get("marginal_parcels", [])),
        output_dir / "maps",
        dlabel_path,
    )
    print(f"Trajectory parcels: {result_dir / 'trajectory_parcels.csv'}")


def rewrite_saved_maps(output_dir: Path, dlabel_path: Path) -> None:
    """Rebuild the MNI volumes and fsLR copies from a finished Tian fit."""
    from analysis.fig03.gam_pygam.brain_maps import write_brain_maps

    statistics_path = output_dir / "statistics" / "pygam_age_comparison.csv"
    note_path = output_dir / "inputs" / "tian_combined.json"
    if not statistics_path.is_file() or not note_path.is_file():
        raise FileNotFoundError(
            f"{output_dir} has no saved Tian fit. "
            "statistics/pygam_age_comparison.csv and inputs/tian_combined.json are required."
        )
    note = json.loads(note_path.read_text(encoding="utf-8"))
    write_brain_maps(
        output_dir / "inputs" / str(note["combined_atlas"]),
        pd.read_csv(statistics_path),
        list(note.get("marginal_parcels", [])),
        output_dir / "maps",
        dlabel_path,
    )


def write_followup_figures(output_dir: Path, statistics: pd.DataFrame, curves: pd.DataFrame,
                           dlabel_path: Path) -> None:
    """Parcel trajectories, group age curves, and the age versus parent_edu scatter, all SVG."""
    from analysis.fig03.gam_nested.plot_nested import plot_age_parent_edu
    from analysis.fig03.gam_pygam.plot_trajectories import (
        plot_group_age_curves,
        plot_significant_trajectories,
        significant_parcels,
    )

    input_dir = output_dir / "inputs"
    result_dir = output_dir / "statistics"
    figure_dir = output_dir / "figures"
    design = analysis_design(pd.read_csv(input_dir / "design.csv"))
    association = plot_age_parent_edu(design, figure_dir / "age_parent_edu.svg")
    association.to_csv(result_dir / "age_parent_edu_correlation.csv", index=False)
    if (input_dir / "tian_combined.json").is_file():
        _write_tian_followup(output_dir, statistics, curves, dlabel_path)
        print(f"Association figure: {figure_dir / 'age_parent_edu.svg'}")
        return
    parcel_info = pd.read_csv(input_dir / "parcel_info.csv")
    selected = significant_parcels(statistics, curves, parcel_info, 0.05)
    list_path = result_dir / "significant_parcels.csv"
    selected.to_csv(list_path, index=False)
    if selected.empty:
        print("No parcels with q_age < 0.05. The trajectory figure was not written.")
    else:
        font_name = None
        from matplotlib import font_manager
        for name in ("Microsoft YaHei", "Noto Sans CJK SC", "Source Han Sans SC", "WenQuanYi Micro Hei"):
            if any(found.name == name for found in font_manager.fontManager.ttflist):
                font_name = name
                break
        if font_name:
            curve_title = f"q_age < 0.05 的脑区（n={len(selected)}）"
            figure_title = "全模型的中心化年龄曲线。协变量取均值，零点是这条曲线自身的平均值。"
        else:
            curve_title = f"Parcels with q_age < 0.05 (n={len(selected)})"
            figure_title = "Centered full-model age curve. Covariates are at their means. Zero is the mean of that curve."
        figure_path = figure_dir / "significant_parcel_trajectories.svg"
        plot_significant_trajectories(
            selected, curves, dlabel_path, figure_path, 0.05,
            curve_title=curve_title, figure_title=figure_title,
        )
        print(f"Trajectory figure: {figure_path}")
    for path in plot_group_age_curves(statistics, curves, figure_dir):
        print(f"Group figure: {path}")
    print(f"Association figure: {figure_dir / 'age_parent_edu.svg'}")
    print(f"Significant parcels: {list_path}")


def run_fit(stat_utils, output_dir: Path, counts: dict, dlabel_path: Path) -> None:
    input_dir = output_dir / "inputs"
    result_dir = output_dir / "statistics"
    store_partial = (input_dir / "tian_combined.json").is_file()
    statistics, curves, residuals = fit_saved_means(stat_utils, input_dir, store_partial=store_partial)
    result_dir.mkdir(parents=True, exist_ok=True)
    statistics_path = result_dir / "pygam_age_comparison.csv"
    curves_path = result_dir / "pygam_curves.csv"
    statistics.to_csv(statistics_path, index=False)
    curves.to_csv(curves_path, index=False)
    if residuals is not None:
        residuals.to_csv(result_dir / "covariate_residuals.csv", index=False)
    write_metadata(output_dir, counts)
    write_followup_figures(output_dir, statistics, curves, dlabel_path)
    for dimension in ("parcel", "yeo", "system", "amygdala"):
        subset = statistics[statistics["dimension"].astype(str) == dimension]
        if subset.empty:
            continue
        n_significant = int((subset["q_age"] < 0.05).sum())
        print(f"{output_dir.name} {dimension}: {len(subset)} units, {n_significant} with q_age < 0.05")
    print(f"Subjects: {counts['n_subjects']}")
    print(f"Statistics: {statistics_path}")
    print(f"Curves: {curves_path}")


def main() -> None:
    parser = argparse.ArgumentParser(description="Historical pygam age model on volume means")
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
                        help="Fit parcel means only. Use this for the 400-parcel Kong 2022 atlas.")
    parser.add_argument("--parcel-info", type=Path, default=None,
                        help="Schaefer order_info.txt. Required for names when --parcels-only is set.")
    parser.add_argument("--tian-atlas", type=Path, default=None,
                        help="Tian Subcortex S2. Labels 1..32 are appended after the Schaefer parcels.")
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--meanfd-max", type=float, nargs="+", default=None,
                        help="Drop subjects with meanFD above each value. 0.5 0.2 writes _FD05 and _FD02.")
    parser.add_argument("--step", choices=("all", "fit", "maps"), default="all",
                        help="fit refits saved means. maps rewrites Tian brain files and does not refit.")
    args = parser.parse_args()
    from analysis.fig03.gam_nested.run_nested import meanfd_directory_suffix, subset_saved_inputs

    output_dir = resolve_project_path(args.output)
    if args.meanfd_max:
        output_dirs = [Path(str(output_dir) + meanfd_directory_suffix(threshold)) for threshold in args.meanfd_max]
    else:
        output_dirs = [output_dir]
    dlabel_path = resolve_project_path(args.schaefer_atlas)
    tian_atlas = None if args.tian_atlas is None else resolve_project_path(args.tian_atlas)
    if args.step == "maps":
        for one_output in output_dirs:
            rewrite_saved_maps(one_output, dlabel_path)
        return
    stat_utils = load_historical_stats()
    if args.step == "fit":
        prepared = []
        for one_output in output_dirs:
            input_dir = one_output / "inputs"
            required = ("design.csv", "parcel_response.csv", "group_response.csv")
            missing = [name for name in required if not (input_dir / name).is_file()]
            if missing:
                raise FileNotFoundError(
                    f"{input_dir} is missing {', '.join(missing)}. Run without --step fit first."
                )
            prepared.append((one_output, saved_counts(input_dir)))
    elif args.meanfd_max:
        shared_inputs = Path(str(output_dir) + "_shared_inputs") / "inputs"
        volume_atlas = resolve_project_path(args.fig3_volume_atlas)
        stamp = shared_inputs.parent / "atlas.txt"
        stamp_text = f"{volume_atlas.resolve()}\n{dlabel_path.resolve()}\n"
        if tian_atlas is not None:
            stamp_text += f"tian:{tian_atlas.resolve()}\n"
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
            shared_counts = saved_counts(shared_inputs)
        else:
            shared_counts = prepare_volume_inputs(
                None if args.input is None else resolve_project_path(args.input),
                resolve_project_path(args.subjects),
                resolve_project_path(args.excluded),
                resolve_project_path(args.parent_edu),
                resolve_project_path(args.meanfd),
                resolve_project_path(args.fig3_volume_dir),
                volume_atlas,
                dlabel_path,
                resolve_project_path(args.system_mask_dir),
                shared_inputs,
                parcels_only=args.parcels_only,
                parcel_info_path=None if args.parcel_info is None else resolve_project_path(args.parcel_info),
                juelich_atlas=None if args.juelich_atlas is None else resolve_project_path(args.juelich_atlas),
                juelich_labels=None if args.juelich_labels is None else resolve_project_path(args.juelich_labels),
                tian_atlas=tian_atlas,
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
        counts = prepare_volume_inputs(
            None if args.input is None else resolve_project_path(args.input),
            resolve_project_path(args.subjects),
            resolve_project_path(args.excluded),
            resolve_project_path(args.parent_edu),
            resolve_project_path(args.meanfd),
            resolve_project_path(args.fig3_volume_dir),
            resolve_project_path(args.fig3_volume_atlas),
            dlabel_path,
            resolve_project_path(args.system_mask_dir),
            output_dir / "inputs",
            parcels_only=args.parcels_only,
            parcel_info_path=None if args.parcel_info is None else resolve_project_path(args.parcel_info),
            juelich_atlas=None if args.juelich_atlas is None else resolve_project_path(args.juelich_atlas),
            juelich_labels=None if args.juelich_labels is None else resolve_project_path(args.juelich_labels),
            tian_atlas=tian_atlas,
        )
        prepared = [(output_dir, counts)]
    for one_output, counts in prepared:
        run_fit(stat_utils, one_output, counts, dlabel_path)


if __name__ == "__main__":
    main()
