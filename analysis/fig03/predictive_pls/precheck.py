"""Read-only sample/site audit for planned Tian432 cross-site prediction.

Requires only NumPy and pandas. No image loading, preprocessing fit, PLS,
feature selection, ERQ completeness filter, or finalized anatomical PFC mask.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import re

import numpy as np
import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[3]
PARCELS = [str(index) for index in range(1, 433)]
COLUMNS = ["sub_id", "age", "gender", "site_id", "parent_edu", "meanFD"]
SCORE_MISSING = {"", "na", "nan", "none", "null", "<na>", "n/a"}


def project_path(path: Path) -> Path:
    return (path if path.is_absolute() else PROJECT_ROOT / path).resolve()


def subject_keys(values: pd.Series, source: str) -> pd.Series:
    def normalize(value):
        value = str(value).strip().casefold()
        if value in SCORE_MISSING:
            return ""
        if value.startswith("sub-"):
            value = value[4:]
        if value.endswith(".0") and value[:-2].isdigit():
            value = value[:-2]
        return value

    keys = values.map(normalize)
    if keys.eq("").any() or keys.duplicated().any():
        examples = values.loc[keys.eq("") | keys.duplicated(keep=False)].head(10).tolist()
        raise ValueError(f"{source}: missing or duplicate normalized subject IDs: {examples}")
    return keys


def read_behavior(path: Path, id_column: str, score_column: str, id_source: Path | None):
    table = pd.read_csv(path, dtype=str, encoding="utf-8-sig", sep=None,
                        engine="python", keep_default_na=False)
    if missing := {id_column, score_column} - set(table):
        raise ValueError(f"Behavior columns missing: {sorted(missing)}; available: {table.columns.tolist()}")
    subject_keys(table[id_column], "Behavior source IDs")
    real_ids = table[id_column].copy()
    if id_source is not None:
        lines = [line.strip() for line in id_source.read_text(encoding="utf-8-sig").splitlines() if line.strip()]
        lines = lines[1:]  # One header, then the user-verified corresponding data rows.
        if len(lines) != len(table):
            raise ValueError("Rowwise ID source and behavior table have different data-row counts")
        real_ids = pd.Series([line.lstrip('"')[:6] for line in lines], index=table.index)
        if not real_ids.str.fullmatch(r"[A-Za-z0-9]{6}").all():
            raise ValueError("Rowwise source must begin every data row with a six-character alphanumeric ID")
    raw = table[score_column]
    score = pd.to_numeric(raw, errors="coerce")
    status = np.select([np.isfinite(score), raw.str.strip().str.casefold().isin(SCORE_MISSING), score.isna()],
                       ["valid", "missing", "not_numeric"], default="nonfinite")
    rows = pd.DataFrame({"source_row": np.arange(1, len(table) + 1),
                         "source_sub_id": table[id_column], "real_sub_id": real_ids,
                         "subject_key": subject_keys(real_ids, "Behavior real IDs"),
                         "score_column": score_column, "score_raw": raw,
                         "success": score.where(np.isfinite(score)), "score_status": status})
    return rows


def pfc_candidates(info: pd.DataFrame) -> pd.DataFrame:
    """Label-based candidate sets, deliberately not anatomical ground truth."""
    info = info[["parcel", "name"]].copy()
    info["parcel"] = pd.to_numeric(info["parcel"], errors="raise")
    if (info["parcel"].duplicated().any() or info["name"].isna().any()
            or set(info["parcel"]) != set(range(1, 433))):
        raise ValueError("parcel_info.csv must uniquely name every parcel 1..432")
    info["parcel"] = info["parcel"].astype(int)
    info = info.sort_values("parcel").reset_index(drop=True)
    cortical = info["parcel"].le(400)
    info["hemisphere"] = info["name"].str.extract(r"_(LH|RH)_", expand=False).fillna("")
    info["candidate_pfc_named"] = cortical & info["name"].str.contains("PFC", case=False, regex=False)
    info["candidate_ofc_named"] = cortical & info["name"].str.contains("OFC", case=False, regex=False)
    info["candidate_pfc_plus_ofc"] = info["candidate_pfc_named"] | info["candidate_ofc_named"]
    info["boundary_review"] = cortical & info["name"].str.contains(r"FrOper|Ins|FEF|Cing|_Med_", case=False, regex=True)
    info["definition_status"] = "label_candidates_only_not_final_anatomical_PFC"
    return info


def audit_filters(roster: pd.DataFrame, rules: list, task: str):
    reason = pd.Series("", index=roster.index, dtype=object)
    summaries, exclusions = [], []
    for name, failed in rules:
        failed = pd.Series(failed, index=roster.index).fillna(True).astype(bool)
        first = reason.eq("") & failed
        summaries.append({"task": task, "filter": name, "n_before": int(reason.eq("").sum()),
                          "n_newly_excluded": int(first.sum()), "n_remaining": int((reason.eq("") & ~failed).sum()),
                          "n_failed_in_original_sample": int((failed & roster["original_in_analysis"]).sum())})
        for index in roster.index[failed & (roster["original_in_analysis"] | (name == "previously_excluded"))]:
            exclusions.append({"task": task, "filter": name, "sub_id": roster.at[index, "sub_id"],
                               "site_id": roster.at[index, "site_id"], "is_first_exclusion": bool(first.at[index])})
        reason.loc[first] = name
    return reason, summaries, exclusions


def site_tables(sample: pd.DataFrame, task: str):
    summaries, splits, bins = [], [], []
    for site, test in sample.groupby("site_id", sort=True):
        train = sample.loc[sample["site_id"].ne(site)]
        train_counts = train.groupby("site_id").size()
        inner_min = len(train) - int(train_counts.max()) if len(train_counts) >= 2 else 0
        summaries.append({"task": task, "site_id": site, "n": len(test),
                          "age_min": test["age"].min(), "age_max": test["age"].max(),
                          "age_mean": test["age"].mean(), "age_sd": test["age"].std(),
                          "success_n_valid": int(test["success"].notna().sum()),
                          "success_mean": test["success"].mean(), "success_sd": test["success"].std(),
                          "gender_1_n": int(test["gender"].eq(1).sum()), "gender_0_n": int(test["gender"].eq(0).sum()),
                          "meanFD_mean": test["meanFD"].mean()})
        splits.append({"task": task, "held_out_site": site, "n_test": len(test), "n_train": len(train),
                       "n_training_sites": len(train_counts), "inner_loso_possible": len(train_counts) >= 2,
                       "min_inner_training_n": inner_min,
                       "test_age_outside_training_range_n": int((~test["age"].between(train["age"].min(), train["age"].max())).sum()),
                       "test_target_constant": test["age" if task == "age" else "success"].nunique() < 2})
        lower = (np.floor(test["age"] / 2) * 2).astype(int)
        for start, count in lower.value_counts().sort_index().items():
            bins.append({"task": task, "site_id": site, "age_bin_left_closed": f"[{start},{start + 2})", "n": int(count)})
    return summaries, splits, bins


def run_precheck(inputs: Path, behavior: Path, output: Path, id_column="sub_ids",
                 score_column="reappraisal_success", id_source: Path | None = None):
    inputs, behavior, output = inputs.resolve(), behavior.resolve(), output.resolve()
    sources = [inputs / name for name in ("design.csv", "parcel_response.csv", "parcel_info.csv", "tian_combined.json")]
    sources += [behavior] + ([id_source.resolve()] if id_source is not None else [])
    if output in (inputs, inputs.parent) or inputs in output.parents or any(path.parent == output for path in sources):
        raise ValueError("Use a separate precheck output directory, outside source inputs and source-file directories")
    if output.exists() and any(output.iterdir()) and not (output / "precheck_summary.json").is_file():
        raise ValueError("Output is nonempty and is not a previous prediction precheck directory")
    note = json.loads((inputs / "tian_combined.json").read_text(encoding="utf-8"))
    if int(note["n_cortical"]) != 400:
        raise ValueError("This precheck is for the planned Tian432 prediction analysis")
    design = pd.read_csv(inputs / "design.csv", dtype=str, keep_default_na=False)
    if missing := set(COLUMNS) - set(design):
        raise ValueError(f"design.csv missing columns: {sorted(missing)}")
    design["subject_key"] = subject_keys(design["sub_id"], "Imaging design")
    design["original_in_analysis"] = (pd.to_numeric(design["in_analysis"], errors="coerce").eq(1)
                                       if "in_analysis" in design else True)
    design["original_excluded_by"] = design.get("excluded_by", pd.Series("", index=design.index))
    for column in ("age", "gender", "parent_edu", "meanFD"):
        design[column] = pd.to_numeric(design[column], errors="coerce")
    design["site_id"] = design["site_id"].str.strip()
    brain = pd.read_csv(inputs / "parcel_response.csv", dtype={"sub_id": str}, keep_default_na=False)
    if set(brain) != {"sub_id", *PARCELS}:
        raise ValueError("parcel_response.csv must contain sub_id and exactly the 432 parcel columns")
    brain_keys = subject_keys(brain["sub_id"], "Brain responses")
    values = brain[PARCELS].apply(pd.to_numeric, errors="coerce")
    brain_audit = pd.DataFrame({"subject_key": brain_keys, "brain_matched": True,
                                "brain_nonfinite_count": (~np.isfinite(values)).sum(axis=1)})
    roster = design.merge(brain_audit,
                          on="subject_key", how="left", validate="one_to_one", sort=False)
    rows = read_behavior(behavior, id_column, score_column, id_source)
    roster = roster.merge(rows[["subject_key", "source_row", "source_sub_id", "score_raw", "score_status", "success"]],
                          on="subject_key", how="left", validate="one_to_one", sort=False)
    base_rules = [("previously_excluded", ~roster["original_in_analysis"])]
    for column in ("age", "gender", "parent_edu", "meanFD"):
        base_rules.append((f"{column}_missing_or_nonfinite", ~np.isfinite(roster[column])))
    base_rules += [("gender_not_0_or_1", ~roster["gender"].isin([0, 1])),
                   ("meanFD_above_0_5", roster["meanFD"].gt(.5)),
                   ("site_id_missing", roster["site_id"].str.casefold().isin(SCORE_MISSING)),
                   ("brain_subject_not_found", roster["brain_matched"].isna()),
                   ("brain_has_nonfinite_parcels", roster["brain_nonfinite_count"].gt(0))]
    task_rules = {"age": base_rules, "success": base_rules + [
        ("behavior_subject_not_found", roster["source_row"].isna()),
        ("success_missing_or_invalid", roster["source_row"].notna() & roster["success"].isna())]}
    output.mkdir(parents=True, exist_ok=True)
    filter_rows, exclusion_rows, sites, splits, bins = [], [], [], [], []
    task_counts = {}
    for task, rules in task_rules.items():
        reason, summary, excluded = audit_filters(roster, rules, task)
        roster[f"excluded_by_{task}"] = reason
        roster[f"in_{task}_analysis"] = reason.eq("")
        sample = roster.loc[reason.eq(""), [*COLUMNS, "success"]]
        sample.to_csv(output / f"{task}_sample.csv", index=False)
        a, b, c = site_tables(sample, task)
        sites.extend(a)
        splits.extend(b)
        bins.extend(c)
        filter_rows.extend(summary)
        exclusion_rows.extend(excluded)
        task_counts[task] = {"n": len(sample), "n_sites": int(sample["site_id"].nunique())}
    roster.to_csv(output / "sample_roster.csv", index=False)
    pd.DataFrame(filter_rows).to_csv(output / "filter_summary.csv", index=False)
    pd.DataFrame(exclusion_rows, columns=["task", "filter", "sub_id", "site_id", "is_first_exclusion"]).to_csv(
        output / "exclusions_by_filter.csv", index=False)
    pd.DataFrame(sites, columns=["task", "site_id", "n", "age_min", "age_max", "age_mean", "age_sd",
                               "success_n_valid", "success_mean", "success_sd", "gender_1_n", "gender_0_n", "meanFD_mean"]).to_csv(
        output / "site_summary.csv", index=False)
    pd.DataFrame(splits, columns=["task", "held_out_site", "n_test", "n_train", "n_training_sites",
                                "inner_loso_possible", "min_inner_training_n", "test_age_outside_training_range_n",
                                "test_target_constant"]).to_csv(output / "loso_feasibility.csv", index=False)
    pd.DataFrame(bins, columns=["task", "site_id", "age_bin_left_closed", "n"]).to_csv(output / "site_age_bins.csv", index=False)
    rows.merge(roster[["subject_key", "sub_id", "original_in_analysis", "in_age_analysis", "in_success_analysis"]],
               on="subject_key", how="left", validate="one_to_one").to_csv(output / "behavior_rows.csv", index=False)
    candidates = pfc_candidates(pd.read_csv(inputs / "parcel_info.csv"))
    candidates.to_csv(output / "parcel_feature_candidates.csv", index=False)
    for flag, filename in (("candidate_pfc_named", "pfc_named_candidates.csv"),
                           ("candidate_pfc_plus_ofc", "pfc_plus_ofc_candidates.csv"),
                           ("boundary_review", "pfc_boundary_review.csv")):
        candidates.loc[candidates[flag]].to_csv(output / filename, index=False)
    quality = pd.DataFrame({"parcel": range(1, 433),
                            "n_nonfinite_in_response_table": (~np.isfinite(values)).sum().to_numpy(),
                            "n_distinct_finite_in_response_table": values.where(np.isfinite(values)).nunique().to_numpy()})
    quality.to_csv(output / "parcel_quality.csv", index=False)
    warnings = []
    for task, counts in task_counts.items():
        if counts["n_sites"] < 3:
            warnings.append(f"{task}: fewer than 3 retained sites; nested outer/inner LOSO is not possible")
        elif counts["n_sites"] < 6:
            warnings.append(f"{task}: only {counts['n_sites']} retained sites; cross-site uncertainty needs caution (a diagnostic flag, not an exclusion rule)")
    if any(row["test_target_constant"] or row["n_test"] < 2 for row in splits):
        warnings.append("Some test sites have a constant target or fewer than 2 participants; site R2 is undefined")
    if any(row["test_age_outside_training_range_n"] for row in splits):
        warnings.append("Some held-out ages are outside the training age range; inspect loso_feasibility.csv")
    if quality["n_distinct_finite_in_response_table"].lt(2).any():
        warnings.append("Some parcel responses are constant/nonfinite; inspect parcel_quality.csv. No feature has been selected or dropped")
    warnings.append("PFC files are label-based candidates; no final anatomical PFC/non-PFC partition has been approved")
    metadata = {"analysis_kind": "prediction_sample_precheck_v1", "tasks": task_counts,
                "score_column": score_column, "rowwise_id_restoration": id_source is not None,
                "rowwise_id_assumption": "user-verified identical data-row order; source supplies IDs only",
                "max_meanFD": .5, "brain_features": 432, "cerebellum_included": False,
                "sample_policy": "same all-432-finite sample across the four models within each task; age does not require success or ERQ",
                "pfc_named_candidates": int(candidates["candidate_pfc_named"].sum()),
                "pfc_plus_ofc_candidates": int(candidates["candidate_pfc_plus_ofc"].sum()),
                "warnings": warnings,
                "source_sha256": {str(path): hashlib.sha256(path.read_bytes()).hexdigest() for path in sources}}
    (output / "precheck_summary.json").write_text(json.dumps(metadata, indent=2, ensure_ascii=False), encoding="utf-8")
    report = ["# PLS 预测分析：样本预检", "", "本报告未拟合模型、未标准化或插补、未选取显著脑区，也未确定最终 PFC 边界。", "",
              f"行为指标：`{score_column}`；原始评分方向和计算公式仍须核实。", ""]
    for task, counts in task_counts.items():
        report.append(f"- {task}：{counts['n']} 人，{counts['n_sites']} 个站点。")
    report += ["", f"PFC 名称候选 {metadata['pfc_named_candidates']} 区；PFC+OFC 名称候选 {metadata['pfc_plus_ofc_candidates']} 区。",
               "这些是候选集，不能将其补集直接宣称为解剖学非前额叶。", "", "## 需检查", ""]
    report += [f"- {warning}" for warning in warnings]
    report += ["", "优先阅读 site_summary.csv、loso_feasibility.csv、filter_summary.csv；逐被试原因见 exclusions_by_filter.csv，行为整合见 behavior_rows.csv。",
               "筛选人数按 n_newly_excluded 相加；同一被试可能不满足多项条件，不能把全部条件失败人数相加。",
               "全脑均值缓存及源行为表只读；输出中的质量统计仅用于诊断，不可拿来代替训练折内处理。"]
    (output / "REPORT.md").write_text("\n".join(report) + "\n", encoding="utf-8")
    print("\n".join(report[:10]), flush=True)
    print(f"Full precheck report: {output / 'REPORT.md'}", flush=True)
    return metadata


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--inputs", type=Path, required=True)
    parser.add_argument("--behavior", type=Path, required=True)
    parser.add_argument("--behavior-id-column", default="sub_ids")
    parser.add_argument("--score-column", default="reappraisal_success")
    parser.add_argument("--behavior-id-source", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    run_precheck(project_path(args.inputs), project_path(args.behavior), project_path(args.output),
                 args.behavior_id_column, args.score_column,
                 project_path(args.behavior_id_source) if args.behavior_id_source else None)


if __name__ == "__main__":
    main()
