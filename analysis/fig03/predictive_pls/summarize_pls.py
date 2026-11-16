"""Read-only PLS reviewer tables from saved outer-site predictions and intervals."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from analysis.fig03 import reporting as io

# Reporting must also run without sklearn or the fitting modules installed.
MODELS = ("M0_covariates", "M1_PFC", "M2_distributed432", "M3_nonPFC")
PRIMARY = "M2_distributed432-M1_PFC"
COMPARISONS = (PRIMARY, "M1_PFC-M0_covariates", "M2_distributed432-M0_covariates",
               "M3_nonPFC-M0_covariates", "M2_distributed432-M3_nonPFC")
TASKS = ("success", "age")
MODEL_NAMES = dict(zip(MODELS, ("Covariates", "PFC", "Distributed", "Non-PFC")))


def _json(path, sources):
    if not path.is_file():
        return {}
    sources.append(path)
    return json.loads(path.read_text(encoding="utf-8"))


def load_results(source, mock=False):
    paths = [source / name for name in ("model_performance.csv", "model_comparisons.csv", "oof_predictions.csv")]
    performance, comparisons, predictions = (io.read_csv(path) for path in paths)
    sources = paths.copy()
    snapshot = _json(source / "frozen/config_snapshot.json", sources)
    bootstrap = _json(source / "frozen/bootstrap_status.json", sources)
    if snapshot.get("simulated", False) and not mock:
        raise ValueError("Simulated input requires --mock")
    io.required(performance, ["task", "model", "n", "n_sites", "n_features", "median_k", "r2", "mae", "rmse", "pearson_r",
                              "prediction_error_improvement_vs_m0", "ci_low_r2", "ci_high_r2", "n_bootstrap_valid"], paths[0])
    io.required(comparisons, ["task", "comparison", "role", "delta_r2", "delta_mae", "delta_rmse", "ci_low_delta_r2",
                              "ci_high_delta_r2", "n_bootstrap_valid", "inference_valid", "supported"], paths[1])
    io.required(predictions, ["task", "sub_id", "site_id", "outer_site", "y_true", *[f"pred_{model}" for model in MODELS]], paths[2])
    for frame, field, expected in ((performance, "model", MODELS), (comparisons, "comparison", COMPARISONS)):
        if frame.duplicated(["task", field]).any() or set(frame.task) != set(TASKS):
            raise ValueError("Duplicate or incomplete PLS tasks/results")
        for task in TASKS:
            if set(frame.loc[frame.task.eq(task), field]) != set(expected):
                raise ValueError(f"Missing PLS {field} rows for {task}")
    if set(predictions.task) != set(TASKS) or predictions.duplicated(["task", "sub_id"]).any():
        raise ValueError("Each subject must have exactly one outer-test prediction per task")
    if predictions[["sub_id", "site_id", "outer_site"]].isna().any().any() or not predictions.site_id.equals(predictions.outer_site):
        raise ValueError("Invalid held-out site identities in saved predictions")
    if predictions.groupby("sub_id").site_id.nunique().gt(1).any():
        raise ValueError("A subject has inconsistent site membership across tasks")
    if predictions.sub_id.str.contains("SIMULATED", case=False).any() and not mock:
        raise ValueError("Simulated subject IDs require --mock")
    membership_path = source / "frozen/pfc_parcel_membership.csv"
    if not membership_path.is_file():
        raise FileNotFoundError("Need frozen/pfc_parcel_membership.csv to identify the actual fitted PFC definition")
    membership = io.read_csv(membership_path)
    io.required(membership, ["parcel", "is_pfc", "definition_source", "rule"], membership_path)
    if membership.parcel.duplicated().any() or set(membership.parcel) != set(range(1, 433)):
        raise ValueError("PFC membership must contain 432 unique parcels")
    pfc = membership.is_pfc.map(io.boolean)
    if pfc.isna().any() or membership.loc[pfc.eq(True), "parcel"].gt(400).any():
        raise ValueError("Invalid PFC membership or subcortical PFC assignments")
    n_pfc = int(pfc.eq(True).sum())
    if not 0 < n_pfc <= 400 or membership.definition_source.nunique() != 1 or membership.rule.nunique() != 1:
        raise ValueError("Missing/ambiguous PFC definition")
    definition, rule = str(membership.definition_source.iloc[0]), str(membership.rule.iloc[0])
    sources.append(membership_path)
    if snapshot and (snapshot.get("n_pfc") != n_pfc or snapshot.get("n_non_pfc") != 432 - n_pfc):
        raise ValueError("Frozen PFC definition disagrees with config snapshot")
    if snapshot.get("primary_comparison", PRIMARY) != PRIMARY:
        raise ValueError("Unexpected prespecified PLS comparison")
    features = dict(zip(MODELS, (0, n_pfc, 432, 432 - n_pfc)))
    index = performance.set_index(["task", "model"])
    # Recalculate simple metrics only to validate saved tables, never refit models.
    for task, group in predictions.groupby("task"):
        y = group.y_true.to_numpy(dtype=float)
        if not np.isfinite(y).all() or group.site_id.nunique() < 3:
            raise ValueError("Invalid target or fewer than three held-out sites")
        total = np.sum((y - y.mean()) ** 2)
        sse_m0 = np.sum((y - group[f"pred_{MODELS[0]}"]) ** 2)
        for model in MODELS:
            predicted = group[f"pred_{model}"].to_numpy(dtype=float)
            if not np.isfinite(predicted).all():
                raise ValueError("Nonfinite outer-test predictions")
            error = y - predicted
            sse = np.sum(error ** 2)
            r = np.nan if np.std(y) <= 1e-12 or np.std(predicted) <= 1e-12 else np.corrcoef(y, predicted)[0, 1]
            expected = [len(y), group.site_id.nunique(), features[model],
                        1 - sse / total if total > 1e-12 else np.nan, np.mean(np.abs(error)), np.sqrt(sse / len(y)), r,
                        1 - sse / sse_m0 if sse_m0 > 1e-12 else np.nan]
            io.close(index.loc[(task, model), ["n", "n_sites", "n_features", "r2", "mae", "rmse", "pearson_r",
                                              "prediction_error_improvement_vs_m0"]], expected, f"OOF vs saved {task}/{model}")
    comparisons = comparisons.copy()
    for field in ("inference_valid", "supported"):
        comparisons[field] = comparisons[field].map(io.boolean)
        if comparisons[field].isna().any():
            raise ValueError(f"Missing {field} flags")
    for row in comparisons.itertuples():
        left, right = row.comparison.split("-")
        io.close([row.delta_r2, row.delta_mae, row.delta_rmse],
                 [index.at[(row.task, left), field] - index.at[(row.task, right), field] for field in ("r2", "mae", "rmse")],
                 f"comparison difference {row.task}/{row.comparison}")
        expected_role = "primary" if row.task == "success" and row.comparison == PRIMARY else "support"
        if row.role != expected_role:
            raise ValueError("PLS primary/support role disagrees with the locked analysis")
    metrics_path = source / "bootstrap_metrics.csv"
    metrics = io.read_csv(metrics_path) if metrics_path.is_file() else None
    if metrics is not None:
        sources.append(metrics_path)
        io.required(metrics, ["task", "replicate", "valid"], metrics_path)
        if metrics.duplicated(["task", "replicate"]).any() or set(metrics.task) != set(TASKS):
            raise ValueError("Duplicate/incomplete bootstrap task records")
        metrics["valid"] = metrics.valid.map(io.boolean)
        if metrics.valid.isna().any():
            raise ValueError("Missing bootstrap validity flags")
    requested = bootstrap.get("requested_replicates", 0)
    if not isinstance(requested, int) or requested < 0:
        raise ValueError("Invalid requested bootstrap count")
    report_states, requested_counts, completed_counts, fractions = [], [], [], []
    for row in comparisons.itertuples():
        valid_count = int(row.n_bootstrap_valid)
        if valid_count < 0 or valid_count != row.n_bootstrap_valid:
            raise ValueError("Invalid valid-bootstrap count")
        if not bootstrap:
            if valid_count or np.isfinite([row.ci_low_delta_r2, row.ci_high_delta_r2]).any():
                raise ValueError("Intervals exist without frozen/bootstrap_status.json; verify the source run")
            if row.inference_valid or row.supported:
                raise ValueError("Inference flags cannot be true without completed bootstrap results")
            state, completed, fraction = "bootstrap_not_run", 0, np.nan
        else:
            if requested <= 0 or valid_count > requested or bootstrap.get("valid_by_task", {}).get(row.task) != valid_count:
                raise ValueError("Bootstrap counts disagree with status metadata")
            fraction = valid_count / requested
            quality = requested >= 2000 and valid_count >= 1000 and fraction >= .95
            if bool(row.inference_valid) != quality:
                raise ValueError("inference_valid disagrees with the prespecified bootstrap quality rules")
            lower, upper = row.ci_low_delta_r2, row.ci_high_delta_r2
            if valid_count and (not np.isfinite([lower, upper]).all() or lower > upper):
                raise ValueError("Missing or reversed bootstrap interval")
            if not valid_count and not pd.isna([lower, upper]).all():
                raise ValueError("Bootstrap interval exists without valid repetitions")
            should_support = quality and row.task == "success" and row.comparison == PRIMARY and np.isfinite(lower) and lower > 0
            if bool(row.supported) != bool(should_support):
                raise ValueError("supported disagrees with the prespecified primary endpoint rule")
            completed = np.nan
            if metrics is not None:
                group = metrics.loc[metrics.task.eq(row.task)]
                if not group.replicate.between(0, requested - 1).all() or not group.replicate.eq(group.replicate.astype(int)).all():
                    raise ValueError("Invalid bootstrap replicate IDs")
                completed = len(group)
                valid = group.loc[group.valid.eq(True)]
                if len(valid) != valid_count:
                    raise ValueError("Saved bootstrap counts disagree with replicate records")
                column = f"delta_r2_{row.comparison}"
                if valid_count:
                    io.required(valid, [column], metrics_path)
                values = valid[column].to_numpy(dtype=float) if valid_count else np.array([])
                if not np.isfinite(values).all():
                    raise ValueError("A valid bootstrap replicate has a nonfinite comparison metric")
                expected = np.percentile(values, [2.5, 97.5]) if len(values) else [np.nan, np.nan]
                io.close([lower, upper], expected, "bootstrap comparison interval")
            state = ("bootstrap_completion_unverified" if pd.isna(completed) else
                     "bootstrap_incomplete" if completed != requested else
                     "bootstrap_quality_insufficient" if not quality else "inference_available")
        report_states.append(state)
        requested_counts.append(requested)
        completed_counts.append(completed)
        fractions.append(fraction)
    comparisons["report_status"] = report_states
    comparisons["report_n_bootstrap_requested"] = requested_counts
    comparisons["report_n_bootstrap_completed"] = completed_counts
    comparisons["report_bootstrap_fraction_valid"] = fractions
    comparisons["confirmation_applicable"] = comparisons.task.eq("success") & comparisons.comparison.eq(PRIMARY)
    comparisons["supported_for_report"] = comparisons.supported & comparisons.report_status.eq("inference_available")
    for column in ("n", "n_sites"):
        comparisons[column] = [index.at[(task, "M1_PFC"), column] for task in comparisons.task]
    comparisons["pfc_definition"], comparisons["pfc_rule"] = definition, rule
    performance = performance.assign(pfc_definition=definition, pfc_rule=rule)
    for task in TASKS:
        rows = comparisons.loc[comparisons.task.eq(task)]
        if rows.n_bootstrap_valid.nunique() != 1:
            raise ValueError("Comparisons from one task must share bootstrap repetitions")
        matched = performance.task.eq(task)
        if not performance.loc[matched, "n_bootstrap_valid"].eq(rows.n_bootstrap_valid.iloc[0]).all():
            raise ValueError("Model performance and comparison tables have different bootstrap counts")
        if metrics is not None and bootstrap:
            valid = metrics.loc[metrics.task.eq(task) & metrics.valid.eq(True)]
            for model in MODELS:
                column = f"r2_{model}"
                if len(valid):
                    io.required(valid, [column], metrics_path)
                values = valid[column].to_numpy(dtype=float) if len(valid) else np.array([])
                if not np.isfinite(values).all():
                    raise ValueError("A valid bootstrap replicate has nonfinite model R2")
                interval = np.percentile(values, [2.5, 97.5]) if len(values) else [np.nan, np.nan]
                io.close(index.loc[(task, model), ["ci_low_r2", "ci_high_r2"]], interval, "model R2 interval")
    return performance, comparisons, sources, definition, rule, snapshot


def run_summary(source, output=None, *, mock=False, docx=True):
    source = Path(source).resolve()
    performance, comparisons, sources, definition, rule, snapshot = load_results(source, mock)
    main = comparisons.loc[comparisons.comparison.eq(PRIMARY)].copy()
    model_index = performance.set_index(["task", "model"])
    for column in ("n", "n_sites"):
        main[column] = [model_index.at[(task, "M1_PFC"), column] for task in main.task]
    for model in MODELS:
        main[f"r2_{model}"] = [model_index.at[(task, model), "r2"] for task in main.task]
    main["pfc_definition"], main["pfc_rule"] = definition, rule
    states = {"bootstrap_not_run": "Bootstrap not run", "bootstrap_completion_unverified": "Completion unverified",
              "bootstrap_incomplete": "Bootstrap incomplete", "bootstrap_quality_insufficient": "Insufficient valid bootstrap",
              "inference_available": "Available"}
    headers = ["Task", "N", "Sites", "Delta R2", "95% CI", "Valid B", "Status"]
    main_rows = [[r.task, str(int(r.n)), str(int(r.n_sites)), io.number(r.delta_r2),
                  io.interval(r.ci_low_delta_r2, r.ci_high_delta_r2), str(int(r.n_bootstrap_valid)), states[r.report_status]]
                 for r in main.itertuples()]
    pheaders = ["Task", "Model", "N", "Features", "R2 [95% CI]", "MAE", "RMSE", "r"]
    prows = [[r.task, MODEL_NAMES[r.model], str(int(r.n)), str(int(r.n_features)),
              f"{io.number(r.r2)} [{io.interval(r.ci_low_r2, r.ci_high_r2)}]", io.number(r.mae), io.number(r.rmse), io.number(r.pearson_r)]
             for r in performance.itertuples()]
    paragraphs = []
    for row in main.itertuples():
        values = (f"{row.task}：N={int(row.n)}，{int(row.n_sites)}站点；全脑相对PFC的"
                  f"ΔR²={io.number(row.delta_r2)}，95%CI={io.interval(row.ci_low_delta_r2, row.ci_high_delta_r2)}。")
        if row.report_status != "inference_available":
            values += f"推断尚不可用于确认性结论（{states[row.report_status]}）；False不能解释为无增量。"
        elif row.task == "age":
            values += "这是次要分析，区间作描述；supported=False不表示年龄预测无差异。"
        elif row.supported_for_report:
            values += "达到预设主比较的正向增量判据；还需结合M2−M0判断是否超过协变量。"
        else:
            values += "未达到预设的正向增量证据判据，不能据此声称两模型等效。"
        baseline = comparisons.loc[comparisons.task.eq(row.task) & comparisons.comparison.eq("M2_distributed432-M0_covariates")].iloc[0]
        values += f" M2−M0的ΔR²={io.number(baseline.delta_r2)}，95%CI={io.interval(baseline.ci_low_delta_r2, baseline.ci_high_delta_r2)}。"
        paragraphs.append(values)
    notes = ("Outer and inner cross-validation use site holdout. Delta R2 is distributed minus PFC. "
             "Intervals use paired site bootstrap with full refitting. Confirmation applies only to success M2 minus M1, "
             "with 2000 requested replicates, at least 1000 valid and at least 95 percent valid, and a positive lower bound. "
             "This report also checks completion against saved replicate records. Age and other comparisons are descriptive. "
             "Negative R2 is retained; negative delta MAE or RMSE indicates lower error. No p, q or BF is inferred from these intervals.")
    text = ("# PLS跨站点预测汇报表\n\n" + f"实际拟合PFC定义：{definition}；规则：{rule}。\n\n" +
            io.markdown_table(headers, main_rows) + "\n\n## 四模型表现\n\n" + io.markdown_table(pheaders, prows) +
            "\n\n## 数值总结\n\n" + "\n\n".join(paragraphs) + "\n\n" + notes +
            "\n\n两个任务可有不同N，同一任务的四模型使用相同被试和划分。"
            "主表2行、模型表现8行、全部预设比较10行。不同PFC定义分别运行，不能拼接成同一分析。"
            "\n\n重评分数保持原方向；在评分方向与公式尚未核实时，不自动写成能力更强或更弱。"
            "年龄预测改善不等于解释重评能力的年龄差异，也不直接证明网络重组。")
    frames = {"pls_main_summary.csv": main, "pls_model_performance.csv": performance,
              "pls_all_comparisons.csv": comparisons}
    sections = ("Cross site PLS prediction results", f"Fitted PFC definition {definition} with rule {rule}",
                [("Distributed versus PFC prediction", headers, main_rows, [.65, .45, .45, .75, 1.5, .65, 2.2], notes),
                 ("Performance of the four models", pheaders, prows, [.65, 1.2, .45, .55, 1.75, .75, .75, .65],
                  "Performance is pooled across outer held-out sites. Success uses age covariates; age prediction does not include age as a covariate.")])
    return io.export(source, output or source / "reviewer_summary", "pls_reviewer_summary_v1", frames, text,
                     sections, sources, mock=mock, docx=docx,
                     extra={"pfc_definition": definition, "pfc_rule": rule, "source_snapshot_checked": bool(snapshot),
                            "score_direction": snapshot.get("score_direction", "not verified"),
                            "bootstrap_report_status": dict(zip(main.task, main.report_status))})


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--analysis-dir", type=Path, required=True)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--skip-docx", action="store_true")
    parser.add_argument("--mock", action="store_true")
    args = parser.parse_args()
    run_summary(io.project_path(args.analysis_dir), io.project_path(args.output) if args.output else None,
                mock=args.mock, docx=not args.skip_docx)


if __name__ == "__main__":
    main()
