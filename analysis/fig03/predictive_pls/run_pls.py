"""Nested site holdout for the locked four-model PLS comparison.

The outer loop holds out one site. Component count is chosen by a second
site holdout inside the outer training sites. Each subject receives one
outer-test prediction per model.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
import sklearn

from analysis.fig03.predictive_pls.models import (
    MODELS,
    NEURAL_MODELS,
    PredictionFailure,
    comparison_rows,
    fit_fold,
    point_performance,
    predict_fold,
    regression_metrics,
    select_k,
    stored_text,
)
from analysis.fig03.predictive_pls.precheck import PARCELS, project_path, run_precheck

CONFIG_PATH = Path(__file__).with_name("config_main.json")
MEMBERSHIP_RULE = "midthickness_area_fraction>=0.50"
MEMBERSHIP_SOURCE = "Donahue2018_liberal_bilateral"
ACCEPTED_MEMBERSHIP = {
    "Donahue2018_liberal_bilateral": ("midthickness_area_fraction>=0.50", 0.50),
    "Buhle2014_reappraisal_PFC": (
        "schaefer2mm_voxel_fraction>=0.25_inside_buhle_z2.3_and_harvard_oxford_pfc", 0.25),
}


def file_sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_config() -> dict:
    return json.loads(CONFIG_PATH.read_text(encoding="utf-8"))


def load_membership(path: Path) -> pd.DataFrame:
    """Accept only the anatomical area rule, including every Tian parcel as non-PFC."""
    table = pd.read_csv(path)
    required = {"parcel", "name", "hemisphere", "coverage_fraction", "is_pfc",
                "boundary", "review_band", "definition_source", "rule"}
    if missing := required - set(table):
        raise PredictionFailure(f"PFC membership is missing columns: {sorted(missing)}")
    table["parcel"] = pd.to_numeric(table["parcel"], errors="raise").astype(int)
    if table["parcel"].duplicated().any() or set(table["parcel"]) != set(range(1, 433)):
        raise PredictionFailure("PFC membership must contain each parcel 1..432 once")
    table = table.sort_values("parcel").reset_index(drop=True)
    sources = set(table["definition_source"])
    rules = set(table["rule"])
    if len(sources) != 1 or len(rules) != 1 or next(iter(sources)) not in ACCEPTED_MEMBERSHIP:
        raise PredictionFailure(
            "PFC membership must use Donahue 2018 liberal cortex or the Buhle 2014 reappraisal PFC map. "
            "Name-based candidate lists are not accepted.")
    expected_rule, threshold = ACCEPTED_MEMBERSHIP[next(iter(sources))]
    if rules != {expected_rule}:
        raise PredictionFailure(f"PFC membership rule must be {expected_rule}")
    coverage = table["coverage_fraction"].to_numpy(dtype=float)
    is_pfc = table["is_pfc"].astype(str).str.lower().isin(("1", "true", "t")).to_numpy()
    if not np.isfinite(coverage).all() or np.any((coverage < 0) | (coverage > 1)):
        raise PredictionFailure("PFC coverage fractions must be finite and inside [0, 1]")
    cortical = table["parcel"].le(400).to_numpy()
    expected = cortical & (coverage >= threshold)
    if not np.array_equal(is_pfc, expected) or is_pfc[~cortical].any() or np.any(coverage[~cortical] != 0):
        raise PredictionFailure(
            f"Tian parcels 401..432 must be non-PFC with zero coverage; cortical labels must follow coverage >= {threshold}")
    if not is_pfc.any() or is_pfc.all():
        raise PredictionFailure("PFC and non-PFC sets must both be nonempty")
    table["is_pfc"] = is_pfc
    return table


def _raw_from_frame(frame: pd.DataFrame) -> dict:
    return {column: frame[column].to_numpy(dtype=float) for column in ("age", "gender", "parent_edu", "meanFD")}


def _slice_raw(raw: dict, index: np.ndarray) -> dict:
    return {key: values[index] for key, values in raw.items()}


def nested_loso(raw: dict, features: dict, y: np.ndarray, sites: np.ndarray, task: str,
                sub_id: np.ndarray, allow_repeated_rows: bool = False) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Return outer predictions, inner scores, and fold records aligned to sub_id."""
    y = np.asarray(y, dtype=float)
    sites = np.asarray(sites).astype(str)
    sub_id = np.asarray(sub_id).astype(str)
    if len({*raw["age"].shape, y.shape[0], sites.shape[0], sub_id.shape[0]}) != 1:
        raise PredictionFailure("Target, sites, and covariates have different lengths")
    identity = pd.DataFrame({"sub_id": sub_id, "site_id": sites}).drop_duplicates()
    if identity["sub_id"].duplicated().any():
        raise PredictionFailure("A subject maps to more than one site")
    if not allow_repeated_rows and pd.Series(sub_id).duplicated().any():
        raise PredictionFailure("Duplicate subject IDs must be resolved before LOSO; this cache has no family-group column")
    unique_sites = pd.unique(sites)
    if len(unique_sites) < 3:
        raise PredictionFailure(f"{task}: fewer than 3 sites, nested LOSO is not defined")
    prediction = {model: np.full(y.size, np.nan) for model in MODELS}
    selected = {model: np.full(y.size, -1, dtype=int) for model in MODELS}
    inner_rows, fold_rows = [], []
    for held in unique_sites:
        test = np.flatnonzero(sites == held)
        train = np.flatnonzero(sites != held)
        if np.intersect1d(train, test).size:
            raise PredictionFailure(f"{task}: outer train and test rows overlap at site {held}")
        train_sites = pd.unique(sites[train])
        if len(train_sites) < 2:
            raise PredictionFailure(f"{task}: outer site {held} leaves fewer than 2 training sites")
        prepared = {}
        limits = []
        limit_notes = []
        for inner in train_sites:
            inner_train = train[sites[train] != inner]
            raw_train = _slice_raw(raw, inner_train)
            for model in NEURAL_MODELS:
                matrix, parcels = features[model]
                state = fit_fold(raw_train["age"], raw_train["gender"], raw_train["parent_edu"], raw_train["meanFD"],
                                 y[inner_train], task, matrix[inner_train], parcels)
                prepared[(inner, model)] = state
                limits.append(state.max_k)
                if state.max_k < 10:
                    limit_notes.append(f"{inner}:{model}:{state.limit_name}={state.max_k}")
        common_max = int(min(limits)) if limits else 0
        chosen = {}
        for model in NEURAL_MODELS:
            scores = []
            for component in range(common_max + 1):
                sse = 0.0
                n_obs = 0
                for inner in train_sites:
                    inner_test = train[sites[train] == inner]
                    raw_test = _slice_raw(raw, inner_test)
                    yhat = predict_fold(prepared[(inner, model)], raw_test["age"], raw_test["gender"],
                                        raw_test["parent_edu"], raw_test["meanFD"],
                                        features[model][0][inner_test], component)
                    if yhat.shape != y[inner_test].shape:
                        raise PredictionFailure("Inner prediction length does not match the held-out site")
                    sse += float(np.sum((y[inner_test] - yhat) ** 2))
                    n_obs += int(inner_test.size)
                mse = sse / n_obs
                scores.append(mse)
                inner_rows.append({"task": task, "outer_site": held, "model": model, "k": component,
                                   "pooled_mse": mse, "n": n_obs})
            chosen[model] = select_k(scores)
        raw_train = _slice_raw(raw, train)
        raw_test = _slice_raw(raw, test)
        covariate_state = fit_fold(raw_train["age"], raw_train["gender"], raw_train["parent_edu"],
                                   raw_train["meanFD"], y[train], task)
        prediction["M0_covariates"][test] = predict_fold(
            covariate_state, raw_test["age"], raw_test["gender"], raw_test["parent_edu"], raw_test["meanFD"], None, 0)
        selected["M0_covariates"][test] = 0
        fold_rows.append(_fold_record(task, held, "M0_covariates", covariate_state, 0, train, train_sites, common_max, limit_notes))
        for model in NEURAL_MODELS:
            matrix, parcels = features[model]
            state = fit_fold(raw_train["age"], raw_train["gender"], raw_train["parent_edu"], raw_train["meanFD"],
                             y[train], task, matrix[train], parcels, components_through=chosen[model])
            if chosen[model] > state.max_k:
                raise PredictionFailure(
                    f"{task}: outer site {held} {model} cannot fit the inner-selected k={chosen[model]} (max {state.max_k})")
            prediction[model][test] = predict_fold(
                state, raw_test["age"], raw_test["gender"], raw_test["parent_edu"], raw_test["meanFD"],
                matrix[test], chosen[model])
            selected[model][test] = chosen[model]
            fold_rows.append(_fold_record(task, held, model, state, chosen[model], train, train_sites, common_max, limit_notes))
    if any(not np.isfinite(values).all() for values in prediction.values()):
        raise PredictionFailure(f"{task}: outer predictions are incomplete")
    rows = pd.DataFrame({"task": task, "sub_id": sub_id, "site_id": sites, "y_true": y})
    for model in MODELS:
        rows[f"pred_{model}"] = prediction[model]
        rows[f"error_{model}"] = y - prediction[model]
        rows[f"k_{model}"] = selected[model]
    rows["outer_site"] = sites
    return rows, pd.DataFrame(inner_rows), pd.DataFrame(fold_rows)


def _fold_record(task, held, model, state, k, train, train_sites, common_max, limit_notes) -> dict:
    return {"task": task, "outer_site": held, "model": model, "k": int(k), "n_train": int(train.size),
            "n_training_sites": int(len(train_sites)), "rank_C": state.rank_c, "rank_X": state.rank_x,
            "n_features_input": int(state.parcel_ids.size), "n_features_used": int(state.keep.sum()),
            "dropped_parcels": ";".join(map(str, state.dropped_parcels)),
            "common_k_max": int(common_max), "limit_name": state.limit_name,
            "truncation": ";".join(limit_notes),
            "selected_k_converged": True if int(k) == 0 else bool(state.converged.get(int(k), False))}


def performance_tables(predictions: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    enriched = predictions.copy()
    return point_performance(enriched), comparison_rows(enriched), site_performance(enriched)


def site_performance(predictions: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for (task, site), frame in predictions.groupby(["task", "site_id"], sort=True):
        y = frame["y_true"].to_numpy(dtype=float)
        constant = y.size < 2 or np.std(y) <= 1e-12
        for model in MODELS:
            metrics = regression_metrics(y, frame[f"pred_{model}"])
            if constant:
                metrics["r2"] = np.nan
                metrics["pearson_r"] = np.nan
            rows.append({"task": task, "model": model, "site_id": site, "n": metrics["n"],
                         "y_min": float(np.min(y)), "y_max": float(np.max(y)), "y_constant": bool(constant),
                         "r2": metrics["r2"], "mae": metrics["mae"], "rmse": metrics["rmse"],
                         "pearson_r": metrics["pearson_r"]})
    table = pd.DataFrame(rows)
    deltas = []
    for (task, site), frame in table.groupby(["task", "site_id"], sort=True):
        by_model = frame.set_index("model")
        deltas.append({"task": task, "site_id": site,
                       "delta_r2_M2_minus_M1": by_model.at["M2_distributed432", "r2"] - by_model.at["M1_PFC", "r2"],
                       "delta_mae_M2_minus_M1": by_model.at["M2_distributed432", "mae"] - by_model.at["M1_PFC", "mae"],
                       "delta_rmse_M2_minus_M1": by_model.at["M2_distributed432", "rmse"] - by_model.at["M1_PFC", "rmse"]})
    return table.merge(pd.DataFrame(deltas), on=["task", "site_id"], how="left")


def _aligned_matrix(response: pd.DataFrame, sample: pd.DataFrame) -> np.ndarray:
    response = response.copy()
    response["subject_key"] = response["sub_id"].map(lambda value: str(value).strip().casefold().removeprefix("sub-"))
    sample = sample.copy()
    sample["subject_key"] = sample["sub_id"].map(lambda value: str(value).strip().casefold().removeprefix("sub-"))
    merged = sample[["subject_key"]].merge(response, on="subject_key", how="left", validate="one_to_one")
    if merged[PARCELS].isna().any().any():
        raise PredictionFailure("Analysis subjects are missing parcel means")
    return merged[PARCELS].to_numpy(dtype=float)


def run_analysis(inputs: Path, behavior: Path, membership: Path, output: Path,
                 id_column="sub_ids", score_column="reappraisal_success",
                 id_source: Path | None = None) -> dict:
    inputs, behavior, membership, output = (path.resolve() for path in (inputs, behavior, membership, output))
    output.mkdir(parents=True, exist_ok=True)
    precheck = run_precheck(inputs, behavior, output / "precheck", id_column, score_column, id_source)
    table = load_membership(membership)
    mask = table["is_pfc"].to_numpy(dtype=bool)
    response = pd.read_csv(inputs / "parcel_response.csv", dtype={"sub_id": str})
    config = load_config()
    predictions, inner_tables, fold_tables = [], [], []
    frozen = output / "frozen"
    frozen.mkdir(exist_ok=True)
    for task in ("success", "age"):
        sample = pd.read_csv(output / "precheck" / f"{task}_sample.csv")
        if sample["site_id"].nunique() < 3:
            raise PredictionFailure(f"{task}: precheck retained fewer than 3 sites")
        raw = _raw_from_frame(sample)
        y_column = "success" if task == "success" else "age"
        y = sample[y_column].to_numpy(dtype=float)
        matrix = _aligned_matrix(response, sample)
        parcels = np.arange(1, 433)
        features = {
            "M1_PFC": (matrix[:, mask], parcels[mask]),
            "M2_distributed432": (matrix, parcels),
            "M3_nonPFC": (matrix[:, ~mask], parcels[~mask]),
        }
        oof, inner, folds = nested_loso(raw, features, y, sample["site_id"].to_numpy(), task, sample["sub_id"].to_numpy())
        predictions.append(oof)
        inner_tables.append(inner)
        fold_tables.append(folds)
        np.savez_compressed(
            frozen / f"{task}_arrays.npz", age=raw["age"], gender=raw["gender"], parent_edu=raw["parent_edu"],
            meanFD=raw["meanFD"], y=y, sites=stored_text(sample["site_id"]),
            sub_id=stored_text(sample["sub_id"]), X=matrix, pfc_mask=mask, parcel_ids=parcels)
    prediction_table = pd.concat(predictions, ignore_index=True)
    if prediction_table.duplicated(["task", "sub_id"]).any():
        raise PredictionFailure("A subject received more than one outer prediction")
    prediction_table.to_csv(output / "oof_predictions.csv", index=False)
    pd.concat(inner_tables, ignore_index=True).to_csv(output / "inner_cv_scores.csv", index=False)
    pd.concat(fold_tables, ignore_index=True).to_csv(output / "fold_models.csv", index=False)
    point, comparisons, sites = performance_tables(prediction_table)
    feature_counts = {"M0_covariates": 0, "M1_PFC": int(mask.sum()),
                      "M2_distributed432": 432, "M3_nonPFC": int((~mask).sum())}
    point["n_features"] = point["model"].map(feature_counts)
    point.to_csv(output / "model_performance.csv", index=False)
    comparisons.to_csv(output / "model_comparisons.csv", index=False)
    sites.to_csv(output / "site_performance.csv", index=False)
    table.to_csv(frozen / "pfc_parcel_membership.csv", index=False)
    manifest = {
        "analysis": config["analysis"],
        "sklearn_requested": config["sklearn_requested"],
        "sklearn_runtime": sklearn.__version__,
        "primary_comparison": config["primary_comparison"],
        "n_pfc": int(mask.sum()),
        "n_non_pfc": int((~mask).sum()),
        "precheck": precheck["tasks"],
        "family_structure": "no family or repeated-measure column in the input cache; duplicate subject IDs are rejected",
        "source_sha256": {
            str(inputs / "design.csv"): file_sha256(inputs / "design.csv"),
            str(inputs / "parcel_response.csv"): file_sha256(inputs / "parcel_response.csv"),
            str(membership): file_sha256(membership),
            str(behavior): file_sha256(behavior),
        },
        "score_direction": "original reappraisal_success score; direction and formula still require verification before prose claims",
    }
    if id_source is not None:
        manifest["source_sha256"][str(id_source)] = file_sha256(id_source)
    (frozen / "config_snapshot.json").write_text(json.dumps(manifest, indent=2, ensure_ascii=False), encoding="utf-8")
    (frozen / "config_main.json").write_text(CONFIG_PATH.read_text(encoding="utf-8"), encoding="utf-8")
    print(f"OOF predictions: {output / 'oof_predictions.csv'}", flush=True)
    return manifest


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--inputs", type=Path, required=True)
    parser.add_argument("--behavior", type=Path, required=True)
    parser.add_argument("--behavior-id-column", default="sub_ids")
    parser.add_argument("--score-column", default="reappraisal_success")
    parser.add_argument("--behavior-id-source", type=Path)
    parser.add_argument("--membership", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    run_analysis(project_path(args.inputs), project_path(args.behavior), project_path(args.membership),
                 project_path(args.output), args.behavior_id_column, args.score_column,
                 project_path(args.behavior_id_source) if args.behavior_id_source else None)


if __name__ == "__main__":
    main()
