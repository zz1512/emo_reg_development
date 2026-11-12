"""Smooth age x habitual-reappraisal models on saved Tian432 response means.

Refit baseline, additive-ERQ and full ERQ-by-age models on identical rows.
Inference uses an unpenalized reference in the same fixed spline space;
predictions and their covariance come from the GCV-selected pyGAM model.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import f as f_distribution, t

from analysis.fig03.gam_nested.amygdala import AMYGDALA_ORDER
from analysis.fig03.gam_pygam.run_pygam import (
    COVARIATES, K_SPLINES, LAM_LOG_RANGE, SYSTEM_ORDER, YEO_ORDER, _age_grid, analysis_design,
    resolve_project_path, validate_saved_inputs,
)

ANALYSIS_KIND = "erq_age_smooth_pygam_v3"
PREVIOUS_ANALYSIS_KIND = "erq_age_smooth_pygam_v2"
LAM_GRID = 10.0 ** np.linspace(*LAM_LOG_RANGE)
TESTS = ("joint_erq", "joint_age_erq", "age", "erq", "interaction", "additive_vs_baseline")
# Retain the historical joint-test name as a compatibility alias only.
PRIMARY_TEST = "full_vs_baseline"
SEQUENTIAL_TESTS = ("additive_vs_baseline", "full_vs_additive")
COMPARISON_FIELDS = ("F", "df1", "df2", "p", "partial_r2_reference", "rss_reference_reduced")
PAIRED_MODEL_COLUMNS = (
    "n_obs_baseline", "n_obs_full", "rss_baseline", "rss_full", "edf_baseline", "edf_full",
    "r2_baseline", "r2_full", "gcv_baseline", "gcv_full", "lambda_age_baseline", "lambda_age_full",
    "lambda_interaction_full", "delta_rss_full_vs_baseline", "delta_r2_full_vs_baseline",
    "rss_reference_baseline", "rss_reference_full", "rank_baseline_reference", "rank_full_reference",
    "n_obs_additive", "rss_additive", "edf_additive", "r2_additive", "gcv_additive",
    "lambda_age_additive", "rss_reference_additive", "rank_additive_reference",
    "delta_rss_additive_vs_baseline", "delta_r2_additive_vs_baseline",
    "delta_rss_full_vs_additive", "delta_r2_full_vs_additive",
)
DIMENSIONS = ("parcel", "yeo", "system", "amygdala")
ALPHA = 0.05
CURVE_COLUMNS = ("dimension", "name", "age", "erq_level", "erq_value", "effect", "ci_low", "ci_high")
SLOPE_COLUMNS = ("dimension", "name", "age", "erq_slope", "ci_low", "ci_high")


def _subject_keys(table: pd.DataFrame, column: str, source: str) -> pd.Series:
    def normalize(value) -> str:
        if pd.isna(value):
            return ""
        value = str(value).strip().casefold()
        if value in {"", "na", "nan", "none", "null", "<na>"}:
            return ""
        if value.startswith("sub-"):
            value = value[4:]
        if value.endswith(".0") and value[:-2].isdigit():
            value = value[:-2]
        return value

    keys = table[column].map(normalize)
    if (keys == "").any():
        raise ValueError(f"{source} has missing subject IDs")
    if keys.duplicated().any():
        duplicates = keys[keys.duplicated(keep=False)].unique()[:10].tolist()
        raise ValueError(f"{source} has duplicate subject IDs after normalization: {duplicates}")
    return keys


def _rowwise_real_ids(path: Path, n_rows: int) -> pd.Series:
    """Restore IDs only under the user's verified row-correspondence contract."""
    lines = [line.strip() for line in path.read_text(encoding="utf-8-sig").splitlines() if line.strip()]
    records = lines[1:]  # The source table has one header row.
    if len(records) != n_rows:
        raise ValueError(f"ERQ ID source row count mismatch: {path} has {len(records)} data rows, "
                         f"but the ERQ score table has {n_rows}. Rowwise ID restoration requires identical row counts.")
    ids = pd.Series([line.lstrip('"')[:6] for line in records], dtype=str)
    invalid = ~ids.str.fullmatch(r"[A-Za-z0-9]{6}")
    if invalid.any():
        examples = [{"data_row": int(index) + 1, "prefix": ids.iloc[index]}
                    for index in np.flatnonzero(invalid.to_numpy())[:10]]
        raise ValueError(f"ERQ ID source requires six alphanumeric ID characters at each data-row start: {examples}")
    _subject_keys(pd.DataFrame({"real_id": ids}), "real_id", "ERQ restored real IDs")
    return ids


def prepare_sample(input_dir: Path, erq_path: Path, score_column: str = "ERQ_CR",
                   id_column: str | None = None, id_source_path: Path | None = None
                   ) -> tuple[pd.DataFrame, pd.DataFrame, dict, pd.DataFrame]:
    """Merge ERQ after the original inclusion filter; retain an exclusion audit."""
    roster = pd.read_csv(input_dir / "design.csv", dtype={"sub_id": str, "site_id": str})
    required = {"sub_id", "age", *COVARIATES}
    if missing := required - set(roster.columns):
        raise ValueError(f"design.csv is missing: {sorted(missing)}")
    roster["_subject_key"] = _subject_keys(roster, "sub_id", "design.csv")
    profiles = pd.read_csv(erq_path, dtype=str, encoding="utf-8-sig", sep=None, engine="python", keep_default_na=False)
    if id_column is None:
        id_column = next((column for column in ("sub_id", "sub_ids", "subject") if column in profiles), None)
    if id_column not in profiles:
        raise ValueError(
            f"ERQ subject ID column {id_column!r} not found; auto-detection accepts sub_id, sub_ids or subject. "
            f"Available columns: {profiles.columns.tolist()!r}. Set --erq-id-column to the actual ID column.")
    if score_column not in profiles:
        hint = ("Use --erq-column ERQ_CR." if "ERQ_CR" in profiles else
                "Set --erq-column to the actual habitual-reappraisal score column.")
        raise ValueError(
            f"ERQ score column {score_column!r} not found. Available columns: {profiles.columns.tolist()!r}. "
            f"{hint}")
    original_keys = _subject_keys(profiles, id_column, "ERQ table")
    profiles["erq_source_row"] = np.arange(1, len(profiles) + 1)
    profiles["erq_real_sub_id"] = (_rowwise_real_ids(id_source_path, len(profiles))
                                   if id_source_path is not None else profiles[id_column])
    profiles["_subject_key"] = (_subject_keys(profiles, "erq_real_sub_id", "ERQ restored real IDs")
                                if id_source_path is not None else original_keys)
    profiles["erq_score_raw"] = profiles[score_column]
    profiles["erq"] = pd.to_numeric(profiles[score_column], errors="coerce")
    empty_score = profiles[score_column].str.strip().str.casefold().isin(
        ["", "na", "nan", "none", "null", "<na>", "n/a"])
    profiles["erq_score_status"] = np.select(
        [np.isfinite(profiles["erq"]), empty_score, profiles["erq"].isna()],
        ["valid", "missing", "not_numeric"], default="nonfinite")
    profiles["_erq_source_id"] = profiles[id_column]
    profiles["_erq_matched"] = True
    roster = roster.drop(columns=["erq"], errors="ignore").merge(
        profiles[["_subject_key", "erq", "_erq_matched", "_erq_source_id", "erq_source_row",
                  "erq_real_sub_id", "erq_score_raw", "erq_score_status"]],
        on="_subject_key", how="left", validate="one_to_one", sort=False)
    eligible_ids = set(analysis_design(roster)["sub_id"].astype(str))
    eligible = roster["sub_id"].astype(str).isin(eligible_ids)
    matched = roster["_erq_matched"].eq(True)
    roster["erq_subject_key"] = roster["_subject_key"]
    roster["erq_subject_matched"] = matched
    roster["original_in_analysis"] = eligible.astype(int)
    roster["original_excluded_by"] = roster.get("excluded_by", pd.Series("", index=roster.index)).fillna("")
    roster["filter_failed_previously_excluded"] = ~eligible
    roster["erq_excluded_by"] = np.where(eligible, "", "previously_excluded")
    stages = {"original_in_analysis": len(eligible_ids)}
    filter_steps = [{"filter": "previously_excluded", "n_before": len(roster),
                     "n_newly_excluded": int((~eligible).sum()), "n_remaining": int(eligible.sum())}]

    def exclude(mask: pd.Series, reason: str) -> None:
        before = int((roster["erq_excluded_by"] == "").sum())
        failed = eligible & mask.fillna(False)
        roster[f"filter_failed_{reason}"] = failed
        removed = (roster["erq_excluded_by"] == "") & failed
        roster.loc[removed, "erq_excluded_by"] = reason
        remaining = int((roster["erq_excluded_by"] == "").sum())
        stages[f"after_{reason}"] = remaining
        filter_steps.append({"filter": reason, "n_before": before,
                             "n_newly_excluded": int(removed.sum()), "n_remaining": remaining})

    exclude(roster["_erq_matched"].isna(), "ERQ_subject_not_found")
    for column in ("erq", "age", "gender", "parent_edu", "meanFD"):
        roster[column] = pd.to_numeric(roster[column], errors="coerce")
        missing_value = ~np.isfinite(roster[column])
        exclude(matched & missing_value if column == "erq" else missing_value, f"{column}_missing_or_nonfinite")
    exclude(~roster["gender"].isin([0, 1]), "gender_not_0_or_1")
    roster["site_id"] = roster["site_id"].astype("string").str.strip()
    exclude(roster["site_id"].isna() | (roster["site_id"] == ""), "site_id_missing")
    roster["in_erq_analysis"] = (roster["erq_excluded_by"] == "").astype(int)
    sample = roster.loc[roster["in_erq_analysis"] == 1,
                        ["sub_id", "age", *COVARIATES, "erq"]].reset_index(drop=True)
    erq_rows = profiles[["erq_source_row", "_erq_source_id", "erq_real_sub_id", "_subject_key",
                         "erq_score_raw", "erq_score_status", "erq"]].rename(columns={
        "_erq_source_id": "erq_source_sub_id", "_subject_key": "erq_subject_key"})
    erq_rows.insert(4, "erq_score_column", score_column)
    imaging = roster[["_subject_key", "sub_id", "original_in_analysis", "original_excluded_by",
                      "in_erq_analysis", "erq_excluded_by"]].rename(columns={
        "_subject_key": "erq_subject_key", "sub_id": "imaging_sub_id"})
    erq_rows = erq_rows.merge(imaging, on="erq_subject_key", how="left", validate="one_to_one", sort=False)
    erq_rows["in_imaging_roster"] = erq_rows["imaging_sub_id"].notna()
    erq_rows["erq_excluded_by"] = erq_rows["erq_excluded_by"].fillna("not_in_imaging_roster")
    erq_rows["original_excluded_by"] = erq_rows["original_excluded_by"].fillna("")
    for column in ("original_in_analysis", "in_erq_analysis"):
        erq_rows[column] = erq_rows[column].fillna(0).astype(int)
    metadata = {
        "erq_file": str(erq_path.resolve()), "erq_column": score_column, "erq_id_column": id_column,
        "n_roster": len(roster), "n_original_in_analysis": len(eligible_ids),
        "n_erq_rows": len(profiles), "n_ids_matched_all_roster": int(matched.sum()),
        "n_original_with_erq_id": int((eligible & matched).sum()),
        "n_original_with_finite_erq": int((eligible & matched & np.isfinite(roster["erq"])).sum()),
        "n_complete_predictors": len(sample),
        "stage_counts": stages, "filter_steps": filter_steps,
        "original_exclusion_counts": roster.loc[~eligible, "original_excluded_by"].value_counts().to_dict(),
        "exclusions": roster.loc[roster["in_erq_analysis"] == 0, "erq_excluded_by"].value_counts().to_dict(),
        "id_examples": {"imaging": roster[["sub_id", "_subject_key"]].head(5).rename(
            columns={"_subject_key": "normalized_id"}).to_dict(orient="records"),
            "erq": profiles[[id_column, "erq_real_sub_id", "_subject_key"]].head(5).rename(
                columns={id_column: "source_id", "_subject_key": "normalized_id"}).to_dict(orient="records")},
        "id_restoration": {"enabled": id_source_path is not None,
                           "source_file": str(id_source_path.resolve()) if id_source_path is not None else None,
                           "method": "user-verified same data-row order, first six characters after one header; equal row counts and unique real IDs required"},
        "erq_aggregation": {"unit": "one score-table row per subject", "operation": "numeric conversion and one-to-one ID join; no score averaging",
                            "n_source_rows": len(profiles), "n_source_rows_outside_imaging": int((~erq_rows["in_imaging_roster"]).sum()),
                            "n_source_rows_in_imaging": int(erq_rows["in_imaging_roster"].sum()),
                            "n_source_rows_originally_excluded": int((erq_rows["in_imaging_roster"] & erq_rows["original_in_analysis"].eq(0)).sum()),
                            "n_source_rows_originally_eligible": int(erq_rows["original_in_analysis"].sum()),
                            "n_source_rows_retained": int(erq_rows["in_erq_analysis"].sum()),
                            "score_status_counts": profiles["erq_score_status"].value_counts().to_dict(),
                            "duplicate_policy": "reject ambiguous duplicate IDs"},
        "id_matching": "strings, whitespace/case normalized, optional sub- prefix and numeric .0 removed; leading zeros preserved; duplicates rejected",
        "age_and_covariate_source": "saved imaging design.csv; profile age/sex never overwrite the imaging design",
    }
    return sample, roster.drop(columns=["_subject_key", "_erq_matched"]).rename(
        columns={"_erq_source_id": "erq_source_sub_id"}), metadata, erq_rows


def build_covariates(sample: pd.DataFrame) -> tuple[np.ndarray, list[str], dict]:
    """Four covariate types; site is a categorical fixed effect."""
    sites = sorted(sample["site_id"].astype(str).unique())
    site = pd.Categorical(sample["site_id"].astype(str), categories=sites)
    encoded = pd.get_dummies(site, prefix="site", drop_first=True, dtype=float)
    numeric = sample[["gender", "parent_edu", "meanFD"]].astype(float).reset_index(drop=True)
    table = pd.concat([numeric, encoded.reset_index(drop=True)], axis=1)
    return table.to_numpy(dtype=float), table.columns.tolist(), {
        "site_levels": sites, "site_reference": sites[0], "columns": table.columns.tolist(),
        "site_encoding": "categorical fixed effects, treatment dummies; arbitrary numeric site distances are not used",
        "centering": "age, ERQ and all nuisance design columns centered within each unit's common complete sample; original units retained",
    }


def _aligned(table: pd.DataFrame, sample: pd.DataFrame, columns: list[str], source: str) -> pd.DataFrame:
    table = table.copy()
    table["_subject_key"] = _subject_keys(table, "sub_id", source)
    requested = sample[["sub_id"]].copy()
    requested["_subject_key"] = _subject_keys(requested, "sub_id", "analysis sample")
    return requested.merge(table[["_subject_key", *columns]], on="_subject_key", how="left",
                           validate="one_to_one", sort=False)


def load_units(input_dir: Path, sample: pd.DataFrame) -> list[tuple[str, str, str, np.ndarray]]:
    units = []
    parcels = pd.read_csv(input_dir / "parcel_response.csv", dtype={"sub_id": str})
    names = [str(number) for number in range(1, 433)]
    matched = _aligned(parcels, sample, names, "parcel_response.csv")
    info = pd.read_csv(input_dir / "parcel_info.csv")
    if info["parcel"].duplicated().any():
        raise ValueError("parcel_info.csv has duplicate parcel IDs")
    labels = info.set_index("parcel")["name"].astype(str).to_dict()
    for name in names:
        units.append(("parcel", name, f"{name} {labels.get(int(name), '')}".strip(),
                      pd.to_numeric(matched[name], errors="coerce").to_numpy(dtype=float)))
    groups = pd.read_csv(input_dir / "group_response.csv", dtype={"sub_id": str})
    for dimension, order in (("yeo", YEO_ORDER), ("system", SYSTEM_ORDER)):
        for name in order:
            subset = groups.loc[(groups["dimension"] == dimension) & (groups["group"] == name)]
            matched = _aligned(subset, sample, ["y"], f"{dimension} {name}")
            units.append((dimension, str(name), f"{dimension} {name}",
                          pd.to_numeric(matched["y"], errors="coerce").to_numpy(dtype=float)))
    amygdala = pd.read_csv(input_dir / "amygdala_response.csv", dtype={"sub_id": str})
    matched = _aligned(amygdala, sample, list(AMYGDALA_ORDER), "amygdala_response.csv")
    for name in AMYGDALA_ORDER:
        units.append(("amygdala", name, f"Juelich {name}",
                      pd.to_numeric(matched[name], errors="coerce").to_numpy(dtype=float)))
    return units


def fit_unit(y: np.ndarray, age: np.ndarray, erq: np.ndarray,
             covariates: np.ndarray, lam_grid: np.ndarray | None = None) -> tuple[dict, dict | None]:
    """Fit three nested GAMs, with fixed-spline-space sequential F tests.

    For inference, center the spline basis over valid observed ages and omit
    its final column. Its four independent columns define the age function;
    their products with centered ERQ define the interaction. This removes
    the constant-by-ERQ alias with the ERQ main effect. It is a reference
    fit, not a test based on the GCV model's approximate EDF or bare term p's.
    """
    from pygam import LinearGAM, l, s

    y, age, erq = (np.asarray(values, dtype=float) for values in (y, age, erq))
    covariates = np.asarray(covariates, dtype=float)
    if covariates.ndim == 1:
        covariates = covariates.reshape(-1, 1)
    finite = np.isfinite(y) & np.isfinite(age) & np.isfinite(erq) & np.all(np.isfinite(covariates), axis=1)
    n = int(finite.sum())
    statistics = {"status": "insufficient_observations", "n_obs": n,
                  "rss_full": np.nan, "rss_reference_full": np.nan, "r2_full": np.nan,
                  "edf_full": np.nan, "lambda_age": np.nan, "lambda_interaction": np.nan}
    statistics.update({column: np.nan for column in PAIRED_MODEL_COLUMNS})
    statistics.update(n_obs_baseline=n, n_obs_additive=n, n_obs_full=n)
    for test in (*TESTS, "full_vs_additive"):
        statistics.update({f"{field}_{test}": np.nan for field in COMPARISON_FIELDS})
    statistics.update({f"{field}_{PRIMARY_TEST}": np.nan for field in COMPARISON_FIELDS})
    for coefficient in ("erq_reference", "erq_additive_reference"):
        for field in ("beta", "se", "ci_low", "ci_high", "t", "df", "p"):
            statistics[f"{field}_{coefficient}"] = np.nan
    if n < 12:
        return statistics, None
    response, observed_age, observed_erq = y[finite], age[finite], erq[finite]
    if np.ptp(observed_age) == 0:
        statistics["status"] = "rank_deficient_targets"
        return statistics, None
    age_center, erq_center = float(observed_age.mean()), float(observed_erq.mean())
    ac, ec = observed_age - age_center, observed_erq - erq_center
    nuisance = covariates[finite] - covariates[finite].mean(axis=0)
    features = np.column_stack([ac, ec, nuisance])
    age_term = s(0, n_splines=K_SPLINES)
    age_term.compile(features)
    raw_basis = age_term.build_columns(features).toarray()
    centered_basis = raw_basis[:, :-1] - raw_basis[:, :-1].mean(axis=0)
    interaction = ec[:, None] * centered_basis
    full = np.column_stack([np.ones(n), centered_basis, ec, interaction, nuisance])
    reduced_cov = np.column_stack([np.ones(n), nuisance])
    reduced_age = np.column_stack([reduced_cov, centered_basis])
    beta, _, rank, _ = np.linalg.lstsq(full, response, rcond=None)
    cov_rank = int(np.linalg.matrix_rank(reduced_cov))
    statistics.update(age_center=age_center, erq_center=erq_center, rank_full=int(rank), df_residual=n - int(rank),
                      age_min=float(observed_age.min()), age_max=float(observed_age.max()),
                      erq_min=float(observed_erq.min()), erq_max=float(observed_erq.max()),
                      erq_sd=float(observed_erq.std(ddof=1)), n_response_excluded=int(len(y) - n))
    if rank - cov_rank != 2 * (K_SPLINES - 1) + 1:
        statistics["status"] = "rank_deficient_targets"
        return statistics, None
    if n <= rank:
        return statistics, None
    residual = response - full @ beta
    rss = float(residual @ residual)
    variation = float(np.sum((response - response.mean()) ** 2))
    if np.ptp(response) == 0 or rss <= np.finfo(float).eps * max(variation, 1e-300):
        statistics["status"] = "zero_residual_variance"
        return statistics, None
    terms = age_term + l(1, penalties=None) + s(0, n_splines=K_SPLINES, by=1)
    for index in range(2, features.shape[1]):
        terms += l(index, penalties=None)
    penalties = np.asarray(LAM_GRID if lam_grid is None else lam_grid, dtype=float).reshape(-1)
    if not penalties.size or np.any(~np.isfinite(penalties) | (penalties <= 0)):
        raise ValueError("Smoothing parameters must be finite positive values")
    # Refit baseline m3 on exactly the full model's complete observations.
    # Independent term instances prevent the two grid searches sharing state.
    baseline_features = np.column_stack([ac, nuisance])
    baseline_terms = s(0, n_splines=K_SPLINES, edge_knots=age_term.edge_knots_.copy())
    for index in range(1, baseline_features.shape[1]):
        baseline_terms += l(index, penalties=None)
    baseline_search = np.array([[age_lam, *([0] * nuisance.shape[1])] for age_lam in penalties])
    baseline_model = LinearGAM(baseline_terms, fit_intercept=True, max_iter=1000).gridsearch(
        baseline_features, response, lam=baseline_search, objective="GCV", progress=False)
    if not baseline_model._is_fitted:
        raise RuntimeError("No baseline m3 model could be fitted by GCV")
    additive_terms = s(0, n_splines=K_SPLINES, edge_knots=age_term.edge_knots_.copy()) + l(1, penalties=None)
    for index in range(2, features.shape[1]):
        additive_terms += l(index, penalties=None)
    additive_search = np.array([[age_lam, 0, *([0] * nuisance.shape[1])] for age_lam in penalties])
    additive_model = LinearGAM(additive_terms, fit_intercept=True, max_iter=1000).gridsearch(
        features, response, lam=additive_search, objective="GCV", progress=False)
    if not additive_model._is_fitted:
        raise RuntimeError("No additive ERQ model could be fitted by GCV")
    search = np.array([[age_lam, 0, interaction_lam, *([0] * nuisance.shape[1])]
                       for age_lam in penalties for interaction_lam in penalties])
    model = LinearGAM(terms, fit_intercept=True, max_iter=1000).gridsearch(
        features, response, lam=search, objective="GCV", progress=False)
    if not model._is_fitted:
        raise RuntimeError("No full age-by-ERQ model could be fitted by GCV")
    prediction = np.asarray(model.predict(features), dtype=float)
    # Verify that the inference reference and the pyGAM model span the same
    # fitted-function space, despite their different identifiability coding.
    pygam_basis = model._modelmat(features).toarray()
    basis_beta, _, basis_rank, _ = np.linalg.lstsq(pygam_basis, response, rcond=None)
    if basis_rank != rank or not np.allclose(pygam_basis @ basis_beta, full @ beta, rtol=1e-7, atol=1e-8):
        raise RuntimeError("The inference reference does not span the fitted pyGAM spline space")
    baseline_basis = baseline_model._modelmat(baseline_features).toarray()
    baseline_basis_beta, _, baseline_basis_rank, _ = np.linalg.lstsq(baseline_basis, response, rcond=None)
    baseline_beta, _, baseline_rank, _ = np.linalg.lstsq(reduced_age, response, rcond=None)
    if baseline_basis_rank != baseline_rank or not np.allclose(
            baseline_basis @ baseline_basis_beta, reduced_age @ baseline_beta, rtol=1e-7, atol=1e-8):
        raise RuntimeError("The baseline reference does not span the fitted baseline pyGAM spline space")
    df_residual = n - int(rank)
    inverse = np.linalg.pinv(full, rcond=np.finfo(float).eps * max(full.shape))
    covariance = (rss / df_residual) * (inverse @ inverse.T)
    critical = float(t.ppf(0.975, df_residual))
    rss_pygam = float(np.sum((response - prediction) ** 2))
    baseline_prediction = np.asarray(baseline_model.predict(baseline_features), dtype=float)
    rss_baseline = float(np.sum((response - baseline_prediction) ** 2))
    baseline_residual = response - reduced_age @ baseline_beta
    statistics.update(status="ok", rss_full=rss_pygam, rss_reference_full=rss,
                      r2_full=1 - rss_pygam / variation, edf_full=float(model.statistics_["edof"]),
                      lambda_age=float(model.lam[0][0]), lambda_interaction=float(model.lam[2][0]))
    statistics.update(
        rss_baseline=rss_baseline, r2_baseline=1 - rss_baseline / variation,
        edf_baseline=float(baseline_model.statistics_["edof"]),
        gcv_baseline=float(baseline_model.statistics_["GCV"]), gcv_full=float(model.statistics_["GCV"]),
        lambda_age_baseline=float(baseline_model.lam[0][0]),
        lambda_age_full=statistics["lambda_age"], lambda_interaction_full=statistics["lambda_interaction"],
        delta_rss_full_vs_baseline=rss_baseline - rss_pygam,
        delta_r2_full_vs_baseline=(rss_baseline - rss_pygam) / variation,
        rss_reference_baseline=float(baseline_residual @ baseline_residual),
        rank_baseline_reference=int(baseline_rank), rank_full_reference=int(rank))
    # ERQ's reference coefficient is the age-averaged ERQ slope, because the
    # interaction basis was centered over the valid subjects' actual ages.
    coefficient = float(beta[K_SPLINES])
    se = float(np.sqrt(max(covariance[K_SPLINES, K_SPLINES], 0.0)))
    statistics.update(beta_erq_reference=coefficient, se_erq_reference=se,
                      ci_low_erq_reference=coefficient - critical * se,
                      ci_high_erq_reference=coefficient + critical * se,
                      t_erq_reference=coefficient / se, df_erq_reference=df_residual,
                      p_erq_reference=float(2 * t.sf(abs(coefficient / se), df_residual)))
    reduced_additive = np.column_stack([reduced_age, ec])
    additive_beta, _, additive_rank, _ = np.linalg.lstsq(reduced_additive, response, rcond=None)
    additive_residual = response - reduced_additive @ additive_beta
    rss_additive_reference = float(additive_residual @ additive_residual)
    additive_basis = additive_model._modelmat(features).toarray()
    additive_basis_beta, _, additive_basis_rank, _ = np.linalg.lstsq(additive_basis, response, rcond=None)
    if additive_basis_rank != additive_rank or not np.allclose(
            additive_basis @ additive_basis_beta, reduced_additive @ additive_beta, rtol=1e-7, atol=1e-8):
        raise RuntimeError("The additive reference does not span the fitted additive pyGAM space")
    additive_df = n - int(additive_rank)
    additive_inverse = np.linalg.pinv(reduced_additive)
    additive_covariance = rss_additive_reference / additive_df * (additive_inverse @ additive_inverse.T)
    additive_se = float(np.sqrt(additive_covariance[-1, -1]))
    additive_coefficient = float(additive_beta[-1])
    additive_critical = float(t.ppf(.975, additive_df))
    additive_rss = float(np.sum((response - additive_model.predict(features)) ** 2))
    statistics.update(
        rss_additive=additive_rss, edf_additive=float(additive_model.statistics_["edof"]),
        r2_additive=1 - additive_rss / variation, gcv_additive=float(additive_model.statistics_["GCV"]),
        lambda_age_additive=float(additive_model.lam[0][0]),
        rss_reference_additive=rss_additive_reference, rank_additive_reference=int(additive_rank),
        delta_rss_additive_vs_baseline=rss_baseline - additive_rss,
        delta_r2_additive_vs_baseline=(rss_baseline - additive_rss) / variation,
        delta_rss_full_vs_additive=additive_rss - rss_pygam,
        delta_r2_full_vs_additive=(additive_rss - rss_pygam) / variation,
        beta_erq_additive_reference=additive_coefficient, se_erq_additive_reference=additive_se,
        ci_low_erq_additive_reference=additive_coefficient - additive_critical * additive_se,
        ci_high_erq_additive_reference=additive_coefficient + additive_critical * additive_se,
        t_erq_additive_reference=additive_coefficient / additive_se, df_erq_additive_reference=additive_df,
        p_erq_additive_reference=float(2 * t.sf(abs(additive_coefficient / additive_se), additive_df)))
    reduced_no_age = np.column_stack([reduced_cov, ec, interaction])
    reduced_no_erq = np.column_stack([reduced_age, interaction])
    for test, reduced in (("joint_erq", reduced_age), ("joint_age_erq", reduced_cov),
                          ("age", reduced_no_age), ("erq", reduced_no_erq),
                          ("interaction", reduced_additive)):
        reduced_beta, _, reduced_rank, _ = np.linalg.lstsq(reduced, response, rcond=None)
        reduced_residual = response - reduced @ reduced_beta
        reduced_rss = float(reduced_residual @ reduced_residual)
        df1 = int(rank) - int(reduced_rank)
        gain = max(0.0, reduced_rss - rss)
        statistic = (gain / df1) / (rss / df_residual)
        statistics.update({f"F_{test}": statistic, f"df1_{test}": df1, f"df2_{test}": df_residual,
                           f"p_{test}": float(f_distribution.sf(statistic, df1, df_residual)),
                           f"partial_r2_reference_{test}": gain / reduced_rss,
                           f"rss_reference_reduced_{test}": reduced_rss})
    for field in COMPARISON_FIELDS:
        statistics[f"{field}_{PRIMARY_TEST}"] = statistics[f"{field}_joint_erq"]
        statistics[f"{field}_full_vs_additive"] = statistics[f"{field}_interaction"]
    baseline_reference_rss = statistics["rss_reference_baseline"]
    gain_additive = max(0.0, baseline_reference_rss - rss_additive_reference)
    added_df = int(additive_rank) - int(baseline_rank)
    additive_f = (gain_additive / added_df) / (rss_additive_reference / additive_df)
    statistics.update(F_additive_vs_baseline=additive_f, df1_additive_vs_baseline=added_df,
                      df2_additive_vs_baseline=additive_df,
                      p_additive_vs_baseline=float(f_distribution.sf(additive_f, added_df, additive_df)),
                      partial_r2_reference_additive_vs_baseline=gain_additive / baseline_reference_rss,
                      rss_reference_reduced_additive_vs_baseline=baseline_reference_rss)
    prediction_df = n - float(model.statistics_["edof"])
    if prediction_df <= 0:
        raise RuntimeError("Full pyGAM prediction CI requires positive residual degrees of freedom")
    return statistics, {
        "model": model, "critical": float(t.ppf(0.975, prediction_df)),
        "age_center": age_center, "erq_center": erq_center,
        "age": observed_age, "erq": observed_erq, "n_covariates": nuisance.shape[1],
        "reference_beta": beta, "reference_covariance": covariance, "reference_df": df_residual,
        "age_term": age_term, "age_basis_mean": raw_basis[:, :-1].mean(axis=0),
    }


def correct_families(table: pd.DataFrame) -> pd.DataFrame:
    """BH on each complete dimension/test family, before follow-up selection."""
    table = table.copy()
    for test in TESTS:
        table[f"q_{test}"] = np.nan
        table[f"reject_{test}"] = False
        for dimension in DIMENSIONS:
            mask = (table["dimension"] == dimension) & np.isfinite(table[f"p_{test}"])
            values = table.loc[mask, f"p_{test}"].to_numpy(dtype=float)
            if not values.size:
                continue
            order = np.argsort(values)
            adjusted = np.minimum.accumulate((values[order] * len(values) / np.arange(1, len(values) + 1))[::-1])[::-1]
            q = np.empty_like(adjusted)
            q[order] = np.clip(adjusted, 0, 1)
            table.loc[mask, f"q_{test}"] = q
            table.loc[mask, f"reject_{test}"] = q < ALPHA
    for target, source in ((PRIMARY_TEST, "joint_erq"), ("full_vs_additive", "interaction")):
        for field in (*COMPARISON_FIELDS, "q", "reject"):
            if f"{field}_{source}" in table:
                table[f"{field}_{target}"] = table[f"{field}_{source}"]
    return table


def selected_units(statistics: pd.DataFrame) -> pd.DataFrame:
    """Either pre-specified sequential test may select a descriptive plot."""
    return statistics.loc[(statistics["q_additive_vs_baseline"] < ALPHA)
                          | (statistics["q_full_vs_additive"] < ALPHA)].copy()


def reference_slopes(dimension: str, name: str, details: dict) -> pd.DataFrame:
    """Pointwise reference ERQ slopes at in-range ages; exploratory, no q test."""
    ages = np.array([6., 8., 10., 12., 14., 16., 18.])
    ages = ages[(ages >= details["age"].min()) & (ages <= details["age"].max())]
    features = np.zeros((len(ages), 2 + details["n_covariates"]))
    features[:, 0] = ages - details["age_center"]
    basis = details["age_term"].build_columns(features).toarray()[:, :-1] - details["age_basis_mean"]
    contrast = np.zeros((len(ages), len(details["reference_beta"])))
    contrast[:, K_SPLINES] = 1
    contrast[:, K_SPLINES + 1:2 * K_SPLINES] = basis
    effect = contrast @ details["reference_beta"]
    se = np.sqrt(np.maximum(np.einsum("ij,jk,ik->i", contrast, details["reference_covariance"], contrast), 0))
    critical = float(t.ppf(.975, details["reference_df"]))
    return pd.DataFrame({"dimension": dimension, "name": name, "age": ages, "beta_erq": effect,
                         "se": se, "ci_low": effect - critical * se, "ci_high": effect + critical * se,
                         "df": details["reference_df"], "inference": "OLS reference pointwise 95% CI; exploratory"})


def primary_comparison_table(statistics: pd.DataFrame) -> pd.DataFrame:
    """A compact, explicit paired-model view derived from the full table."""
    columns = ["dimension", "name", "label", "status", "n_obs", *PAIRED_MODEL_COLUMNS,
               *(f"{field}_{PRIMARY_TEST}" for field in (*COMPARISON_FIELDS, "q", "reject"))]
    result = statistics.loc[:, columns].copy()
    result.insert(3, "comparison", PRIMARY_TEST)
    return result


def conditional_predictions(dimension: str, name: str, details: dict) -> tuple[pd.DataFrame, pd.DataFrame]:
    grid = _age_grid(details["age"])
    model = details["model"]
    beta, covariance = np.asarray(model.coef_), np.asarray(model.statistics_["cov"])
    ac = grid - details["age_center"]
    at_mean_erq = np.column_stack([ac, np.zeros((len(grid), 1 + details["n_covariates"]))])
    observed_at_mean = np.column_stack([
        details["age"] - details["age_center"],
        np.zeros((len(details["age"]), 1 + details["n_covariates"]))])
    reference = model._modelmat(observed_at_mean).toarray().mean(axis=0)
    curves = []
    for level, score in zip(("low", "median", "high"), np.quantile(details["erq"], [.16, .5, .84])):
        ec = float(score - details["erq_center"])
        features = at_mean_erq.copy()
        features[:, 1] = ec
        design = model._modelmat(features).toarray()
        effect = (design - reference) @ beta
        half_width = details["critical"] * np.sqrt(np.maximum(np.einsum("ij,jk,ik->i", design, covariance, design), 0))
        curves.append(pd.DataFrame({"dimension": dimension, "name": name, "age": grid,
                                    "erq_level": level, "erq_value": float(score), "effect": effect,
                                    "ci_low": effect - half_width, "ci_high": effect + half_width}))
    plus_one_erq = at_mean_erq.copy()
    plus_one_erq[:, 1] = 1.0
    contrast = model._modelmat(plus_one_erq).toarray() - model._modelmat(at_mean_erq).toarray()
    slope = contrast @ beta
    half_width = details["critical"] * np.sqrt(np.maximum(np.einsum("ij,jk,ik->i", contrast, covariance, contrast), 0))
    slopes = pd.DataFrame({"dimension": dimension, "name": name, "age": grid, "erq_slope": slope,
                           "ci_low": slope - half_width, "ci_high": slope + half_width})
    return pd.concat(curves, ignore_index=True), slopes


def _write_figures(output_dir: Path, statistics: pd.DataFrame, curves: pd.DataFrame, slopes: pd.DataFrame) -> None:
    from analysis.fig03.gam_pygam.plot_erq import write_erq_figures

    manifest_path = output_dir / "statistics" / "figure_manifest.csv"
    previous = pd.read_csv(manifest_path) if manifest_path.is_file() else pd.DataFrame()
    selected = selected_units(statistics)
    manifest = write_erq_figures(curves, slopes, selected, output_dir / "figures")
    figure_columns = ("age_figure_file", "slope_figure_file")
    new_paths = {str(value) for column in figure_columns for value in manifest.get(column, [])}
    root = (output_dir / "figures").resolve()
    for column in figure_columns:
        for value in previous.get(column, []):
            if str(value) in new_paths:
                continue
            path = (root / str(value)).resolve()
            if root not in path.parents:
                raise ValueError(f"Old figure manifest points outside its figures directory: {value}")
            if path.is_file():
                path.unlink()
    manifest.to_csv(manifest_path, index=False)


def _write_sample_audits(directory: Path, roster: pd.DataFrame, erq_rows: pd.DataFrame, metadata: dict) -> None:
    """Keep source-row provenance and all failed filters beside sequential counts."""
    directory.mkdir(parents=True, exist_ok=True)
    roster.to_csv(directory / "erq_sample_roster.csv", index=False)
    erq_rows.to_csv(directory / "erq_rows.csv", index=False)
    erq_rows[["erq_source_row", "erq_source_sub_id", "erq_real_sub_id", "imaging_sub_id", "in_imaging_roster"]].to_csv(
        directory / "erq_id_crosswalk.csv", index=False)
    summary = pd.DataFrame(metadata["filter_steps"])
    summary["n_failed_criterion"] = [int(roster[f"filter_failed_{name}"].sum()) for name in summary["filter"]]
    summary.to_csv(directory / "erq_filter_summary.csv", index=False)
    columns = ["sub_id", "original_excluded_by", "erq_source_row", "erq_source_sub_id", "erq_real_sub_id",
               "erq_score_raw", "erq_score_status", "erq", "erq_excluded_by"]
    failures = []
    for name in summary["filter"]:
        part = roster.loc[roster[f"filter_failed_{name}"], columns].copy()
        part.insert(0, "filter", name)
        part.insert(1, "is_first_exclusion", part["erq_excluded_by"].eq(name))
        failures.append(part)
    pd.concat(failures, ignore_index=True).to_csv(directory / "erq_exclusions_by_filter.csv", index=False)
    original = roster.loc[roster["original_in_analysis"] == 0, columns].copy()
    original.insert(0, "filter", original["original_excluded_by"].replace("", "original_reason_not_recorded"))
    original.to_csv(directory / "original_exclusions_by_filter.csv", index=False)
    original.groupby("filter").size().rename("n_excluded").reset_index().to_csv(
        directory / "original_exclusion_summary.csv", index=False)
    (directory / "erq_sample_summary.json").write_text(json.dumps(metadata, indent=2, ensure_ascii=False), encoding="utf-8")


def prepare_checked_sample(input_dir: Path, erq_path: Path, output_dir: Path,
                           erq_column: str = "ERQ_CR", erq_id_column: str | None = None,
                           erq_id_source: Path | None = None
                           ) -> tuple[pd.DataFrame, pd.DataFrame, dict, pd.DataFrame]:
    """Validate and audit sample matching before any model fitting."""
    input_dir, output_dir = input_dir.resolve(), output_dir.resolve()
    if output_dir in (input_dir, input_dir.parent) or input_dir in output_dir.parents:
        raise ValueError("ERQ output must be separate from the source age-analysis directory and its inputs")
    metadata_path = output_dir / "run_metadata.json"
    if metadata_path.is_file() and json.loads(metadata_path.read_text(encoding="utf-8")).get("analysis_kind") not in (
            ANALYSIS_KIND, PREVIOUS_ANALYSIS_KIND, "erq_age_smooth_pygam_v1"):
        raise ValueError("Output contains another analysis; choose a separate ERQ output directory")
    note = validate_saved_inputs(input_dir)
    if note is None or int(note["n_cortical"]) != 400:
        raise ValueError("This ERQ analysis requires complete Tian432 + Yeo7 + system4 + Juelich6 saved inputs")
    sample, roster, sample_metadata, erq_rows = prepare_sample(input_dir, erq_path, erq_column, erq_id_column, erq_id_source)
    diagnostic_dir = output_dir / "diagnostics"
    _write_sample_audits(diagnostic_dir, roster, erq_rows, sample_metadata)
    summary_path = diagnostic_dir / "erq_sample_summary.json"
    print(f"ERQ columns: ID={sample_metadata['erq_id_column']!r}; score={sample_metadata['erq_column']!r}", flush=True)
    if erq_id_source is not None:
        print(f"ERQ real IDs restored from first six characters of verified corresponding rows: {erq_id_source}", flush=True)
    print(f"Sample counts: roster={sample_metadata['n_roster']}; original included={sample_metadata['n_original_in_analysis']}; "
          f"ERQ rows={sample_metadata['n_erq_rows']}; original ID matched={sample_metadata['n_original_with_erq_id']}; "
          f"finite ERQ matched={sample_metadata['n_original_with_finite_erq']}; complete={len(sample)}", flush=True)
    print(f"First exclusion reasons: {sample_metadata['exclusions']}; diagnostics: {summary_path}", flush=True)
    if sample_metadata["n_original_in_analysis"] == 0:
        raise ValueError(f"No subjects marked in_analysis=1 in the source design. Diagnostics: {summary_path}")
    if sample_metadata["n_original_with_erq_id"] == 0:
        raise ValueError(
            "ERQ subject ID overlap is 0 among originally included subjects. Imaging and ERQ tables may use different "
            "ID systems; use an ERQ file with matching original IDs or a verified mapping. "
            "For a verified row-correspondence table, specify --erq-id-source; do not infer row correspondence. "
            f"ID examples: {sample_metadata['id_examples']!r}. Diagnostics: {summary_path}")
    if len(sample) < 12:
        raise ValueError(f"Only {len(sample)} complete subjects after ERQ matching; at least 12 are required. "
                         f"Exclusions: {sample_metadata['exclusions']!r}. Diagnostics: {summary_path}")
    if (sample["meanFD"] > .5).any():
        raise ValueError("This ERQ run is restricted to the existing FD <= 0.5 cache; no extra silent filtering is applied")
    return sample, roster, sample_metadata, erq_rows


def run_analysis(input_dir: Path, erq_path: Path, output_dir: Path,
                 erq_column: str = "ERQ_CR", erq_id_column: str | None = None,
                 erq_id_source: Path | None = None, skip_docx: bool = False) -> pd.DataFrame:
    from analysis.fig03.gam_pygam.result_tables import check_docx_dependency, write_erq_tables

    input_dir, output_dir = input_dir.resolve(), output_dir.resolve()
    metadata_path = output_dir / "run_metadata.json"
    sample, roster, sample_metadata, erq_rows = prepare_checked_sample(
        input_dir, erq_path, output_dir, erq_column, erq_id_column, erq_id_source)
    if not skip_docx:
        check_docx_dependency()
    covariates, _, encoding = build_covariates(sample)
    units = load_units(input_dir, sample)
    print(f"Complete predictor sample: {len(sample)}; units: {len(units)}; "
          f"GCV candidates per valid unit: baseline={len(LAM_GRID)}, additive={len(LAM_GRID)}, full={len(LAM_GRID) ** 2}", flush=True)
    rows, fitted, response_exclusions = [], {}, []
    for index, (dimension, name, label, response) in enumerate(units, start=1):
        response_exclusions.extend({"dimension": dimension, "name": name, "label": label, "sub_id": subject,
                                    "reason": "contrast_missing_or_nonfinite", "models": "baseline_additive_full"}
                                   for subject in sample.loc[~np.isfinite(response), "sub_id"])
        row, details = fit_unit(response, sample["age"].to_numpy(), sample["erq"].to_numpy(), covariates)
        row.update(dimension=dimension, name=name, label=label)
        rows.append(row)
        if details is not None:
            fitted[(dimension, name)] = details
        if index % 50 == 0 or index == len(units):
            print(f"ERQ baseline/additive/full triples: {index}/{len(units)}", flush=True)
    statistics = correct_families(pd.DataFrame(rows))
    leading_columns = ["dimension", "name", "label", "status", "n_obs"]
    statistics = statistics[leading_columns + [column for column in statistics if column not in leading_columns]]
    joint_selected = statistics.loc[statistics[f"q_{PRIMARY_TEST}"] < ALPHA].copy()
    selected = selected_units(statistics)
    curve_frames, slope_frames = [], []
    for row in selected.itertuples():
        curves, slopes = conditional_predictions(row.dimension, row.name, fitted[(row.dimension, row.name)])
        curve_frames.append(curves)
        slope_frames.append(slopes)
    curves = pd.concat(curve_frames, ignore_index=True) if curve_frames else pd.DataFrame(columns=CURVE_COLUMNS)
    slopes = pd.concat(slope_frames, ignore_index=True) if slope_frames else pd.DataFrame(columns=SLOPE_COLUMNS)
    result_dir, audit_dir = output_dir / "statistics", output_dir / "inputs"
    result_dir.mkdir(parents=True, exist_ok=True)
    audit_dir.mkdir(parents=True, exist_ok=True)
    sample.to_csv(audit_dir / "analysis_design.csv", index=False)
    _write_sample_audits(audit_dir, roster, erq_rows, sample_metadata)
    pd.DataFrame(response_exclusions, columns=["dimension", "name", "label", "sub_id", "reason", "models"]).to_csv(
        result_dir / "unit_response_exclusions.csv", index=False)
    statistics.to_csv(result_dir / "erq_age_statistics.csv", index=False)
    comparison = primary_comparison_table(statistics)
    comparison.to_csv(result_dir / "erq_full_vs_baseline.csv", index=False)
    comparison.loc[comparison[f"q_{PRIMARY_TEST}"] < ALPHA].to_csv(
        result_dir / "full_vs_baseline_significant.csv", index=False)
    joint_selected.to_csv(result_dir / "joint_erq_significant.csv", index=False)
    selected.to_csv(result_dir / "sequential_significant.csv", index=False)
    for test in SEQUENTIAL_TESTS:
        statistics.loc[statistics[f"q_{test}"] < ALPHA].to_csv(result_dir / f"{test}_significant.csv", index=False)
    reference_frames = [reference_slopes(dimension, name, details) for (dimension, name), details in fitted.items()]
    references = (pd.concat(reference_frames, ignore_index=True) if reference_frames else pd.DataFrame(
        columns=["dimension", "name", "age", "beta_erq", "se", "ci_low", "ci_high", "df", "inference"]))
    references.to_csv(result_dir / "erq_reference_slopes_at_ages.csv", index=False)
    write_erq_tables(statistics, result_dir, include_docx=not skip_docx)
    curves.to_csv(result_dir / "conditional_age_curves.csv", index=False)
    slopes.to_csv(result_dir / "erq_slopes_by_age.csv", index=False)
    _write_figures(output_dir, statistics, curves, slopes)
    metadata = {
        "analysis_kind": ANALYSIS_KIND, "source_inputs": str(input_dir),
        "baseline_model": "intercept + s(age, n_splines=5) + gender + parent_edu + meanFD + categorical(site_id)",
        "additive_model": "baseline + ERQ_c",
        "full_model": "intercept + s(age, n_splines=5) + ERQ_c + ERQ_c*g(age, n_splines=5) + gender + parent_edu + meanFD + categorical(site_id)",
        "estimation": "pyGAM: s(0,n_splines=5) + l(1,penalties=None) + s(0,n_splines=5,by=1) + unpenalized nuisance terms",
        "paired_models": "baseline, additive and full refitted on identical complete rows, age spline space, centered nuisance columns and categorical site coding; historical m3 statistics are not reused",
        "primary_comparisons": list(SEQUENTIAL_TESTS),
        "primary_comparison": PRIMARY_TEST, "primary_comparison_file": "statistics/erq_full_vs_baseline.csv",
        "compatibility_alias": "joint_erq and full_vs_baseline are identical tests, including F/p/q/df/partial_r2_reference",
        "smoothing": {"method": "GCV", "baseline_combinations": len(LAM_GRID),
                      "additive_combinations": len(LAM_GRID),
                      "full_combinations": len(LAM_GRID) ** 2, "age_lambda_grid": LAM_GRID.tolist(),
                      "interaction_lambda_grid": LAM_GRID.tolist(), "combinations": len(LAM_GRID) ** 2,
                      "note": "independent penalties for age and interaction; ERQ/nuisance/intercept unpenalized"},
        "inference": "nested F tests of an unpenalized OLS reference in the full fixed spline space; matrix-rank df, not GCV EDF tests or pyGAM term p-values; independent homoskedastic Gaussian errors",
        "tests": {"full_vs_baseline": "full spline space vs baseline s(age) + covariates; ERQ main + varying coefficient (normally 5 df); joint_erq is a compatibility alias",
                  "additive_vs_baseline": "additive ERQ model vs baseline; ERQ slope without interaction (1 df), using additive residual df",
                  "full_vs_additive": "full vs additive ERQ model; interaction block (normally 4 df); interaction is a compatibility alias",
                  "joint_age_erq": "full spline space vs covariates; age + ERQ + varying coefficient (normally 9 df)",
                  "age": "age function at mean ERQ, interaction retained (normally 4 df)",
                  "erq": "age-averaged ERQ slope in the centered-basis reference (1 df)",
                  "interaction": "full spline space vs additive s(age) + ERQ + covariates (normally 4 df)"},
        "reference_coding": "center age B-spline columns over observed ages, omit the final column, multiply remaining 4 columns by centered ERQ; removes constant-by-ERQ alias",
        "effect_sizes": "rss/r2/edf/gcv_baseline and _full use independently GCV-selected pyGAM fits; delta_rss = RSS_baseline - RSS_full and delta_r2 = R2_full - R2_baseline, signed and unclipped; rss_reference_* and partial_r2_reference_* use only OLS inference references, not the penalized fits",
        "bh": "each test separately within parcel432, Yeo7, system4 and Juelich6; all valid units before selection; not across-family FDR control",
        "followup_selection": "either q_additive_vs_baseline < 0.05 or q_full_vs_additive < 0.05; no gate between the two pre-specified tests; model-implied plots, no subgroup fitting",
        "tables": "all 449 units retained, including failures; primary ERQ coefficient from additive OLS reference; interaction represented by a block F test, not a scalar beta; signed penalized fit gains kept separate",
        "bayes_factors": "not computed for ERQ",
        "prediction_ci": "95% pointwise full-mean CI, conditional on GCV lambdas, shifted by the same average prediction over valid observed ages at mean ERQ/covariates for all ERQ values",
        "slope_ci": "95% pointwise fitted-model contrast CI for d(activation)/d(ERQ) = beta_ERQ + g(age), conditional on GCV lambdas",
        "sample": sample_metadata, "covariate_encoding": encoding,
        "unit_response_exclusion_records": len(response_exclusions),
        "n_units": statistics.groupby("dimension").size().to_dict(),
        "statuses": statistics["status"].value_counts().to_dict(),
        "n_joint_erq_significant": joint_selected.groupby("dimension").size().to_dict(),
        "n_full_vs_baseline_significant": joint_selected.groupby("dimension").size().reindex(DIMENSIONS, fill_value=0).to_dict(),
        "n_sequential_significant": {test: statistics.loc[statistics[f"q_{test}"] < ALPHA].groupby(
            "dimension").size().reindex(DIMENSIONS, fill_value=0).to_dict() for test in SEQUENTIAL_TESTS},
    }
    metadata_path.write_text(json.dumps(metadata, indent=2, ensure_ascii=False), encoding="utf-8")
    print(f"Complete predictor sample: {len(sample)}; valid models: {(statistics['status'] == 'ok').sum()}/{len(statistics)}")
    print(f"Sequential FDR significant (union): {len(selected)}; manuscript tables: {result_dir}")
    return statistics


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--inputs", type=Path, help="Existing Tian432 age-analysis inputs directory")
    parser.add_argument("--erq", type=Path, default=Path("/public/home/dingrui/fmri_analysis/data/beh/er_profiles_new.csv"))
    parser.add_argument("--erq-column", default="ERQ_CR", help="Habitual-reappraisal score column (default: ERQ_CR)")
    parser.add_argument("--erq-id-column", default=None, help="Subject ID column; auto-detect sub_id, sub_ids or subject")
    parser.add_argument("--erq-id-source", type=Path, default=None,
                        help="User-verified row-corresponding table: first six characters of each data row restore real IDs")
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--step", choices=("fit", "figures", "check", "tables"), default="fit",
                        help="fit three models, redraw figures, export saved tables, or audit matching")
    parser.add_argument("--skip-docx", action="store_true", help="Write CSV only; use --step tables later for Word export")
    args = parser.parse_args()
    output_dir = resolve_project_path(args.output)
    if args.step in ("figures", "tables"):
        metadata = json.loads((output_dir / "run_metadata.json").read_text(encoding="utf-8"))
        if metadata.get("analysis_kind") != ANALYSIS_KIND:
            raise ValueError("The output lacks the three-model schema; rerun --step fit to fit all three models")
        result_dir = output_dir / "statistics"
        if args.step == "tables":
            from analysis.fig03.gam_pygam.result_tables import write_erq_tables

            write_erq_tables(pd.read_csv(result_dir / "erq_age_statistics.csv", dtype={"name": str}),
                             result_dir, include_docx=not args.skip_docx)
            return
        tables = [pd.read_csv(result_dir / name, dtype={"name": str}) for name in (
            "erq_age_statistics.csv", "conditional_age_curves.csv", "erq_slopes_by_age.csv")]
        _write_figures(output_dir, *tables)
        return
    if args.inputs is None:
        parser.error("--inputs is required for --step fit or check")
    if args.step == "check":
        prepare_checked_sample(resolve_project_path(args.inputs), resolve_project_path(args.erq), output_dir,
                               args.erq_column, args.erq_id_column,
                               resolve_project_path(args.erq_id_source) if args.erq_id_source is not None else None)
        return
    run_analysis(resolve_project_path(args.inputs), resolve_project_path(args.erq), output_dir,
                 args.erq_column, args.erq_id_column,
                 resolve_project_path(args.erq_id_source) if args.erq_id_source is not None else None,
                 skip_docx=args.skip_docx)


if __name__ == "__main__":
    main()
