"""Extract age BF reviewer tables from completed Tian432/232 Bayesian fits."""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd

from analysis.fig03 import reporting as io
from analysis.fig03.gam_pygam.result_tables import age_long_table

COMPARISONS = ("m3_vs_m1", "m3_vs_m2", "m2_vs_m1")
KEY = ["dimension", "name"]
LABELS = {"parcel": "Whole-brain parcels", "yeo": "Yeo networks", "system": "Systems", "amygdala": "Amygdala subdivisions"}
BF_NAMES = dict(zip(COMPARISONS, ("BF31", "BF32", "BF21")))


def load_results(source, mock=False):
    path = source / "age_bayes_statistics.csv"
    wide = io.read_csv(path)
    sources = [path]
    fields = [*KEY, "label", "status", "n_obs", "integration_converged", "quadrature_relative_error", "quadrature_check_log_difference"]
    fields += [f"{field}_{comparison}" for comparison in COMPARISONS for field in
               ("bf10", "bf01", "log_bf10", "evidence", "sensitivity_log_bf10_min", "sensitivity_log_bf10_max", "sensitivity_direction_stable")]
    io.required(wide, fields, path)
    if wide[KEY].isna().any().any() or wide.duplicated(KEY).any():
        raise ValueError("Missing or duplicate Bayesian unit identities")
    metadata_path = source / "run_metadata.json"
    metadata = json.loads(metadata_path.read_text(encoding="utf-8")) if metadata_path.is_file() else {}
    if metadata_path.is_file():
        sources.append(metadata_path)
        if metadata.get("analysis_kind") != "age_bayesian_spline_v1":
            raise ValueError("Expected age_bayesian_spline_v1 results")
        if metadata.get("max_meanFD") != .5 or metadata.get("prior", {}).get("main_scales") != [.5, .5]:
            raise ValueError("Expected the FD<=0.5 analysis and main prior scales 0.5/0.5")
    if (metadata.get("simulated", False) or wide.label.astype(str).str.contains("SIMULATED", case=False).any()) and not mock:
        raise ValueError("Simulated input requires --mock")
    n_parcels = int(wide.dimension.eq("parcel").sum())
    families = {"parcel": n_parcels}
    if n_parcels == 432:
        families.update(yeo=7, system=4, amygdala=6)
    elif n_parcels != 232:
        raise ValueError("Expected 432+7+4+6 or 232 units")
    if set(wide.dimension) != set(families) or wide.groupby("dimension").size().to_dict() != families:
        raise ValueError("Missing/extra Bayesian analysis units")
    if set(wide.loc[wide.dimension.eq("parcel"), "name"]) != set(map(str, range(1, n_parcels + 1))):
        raise ValueError("Incorrect parcel IDs")
    if metadata and metadata.get("n_cortical") != n_parcels - 32:
        raise ValueError("Atlas count differs from source metadata")
    for row in wide.to_dict("records"):
        if row["status"] != "ok":
            continue
        if not np.isfinite(row["n_obs"]) or row["n_obs"] <= 0 or row["n_obs"] != int(row["n_obs"]):
            raise ValueError("Invalid n_obs in an ok Bayesian fit")
        if io.boolean(row["integration_converged"]) is not True:
            raise ValueError("An ok Bayesian fit has unconverged integration")
        for field in ("quadrature_relative_error", "quadrature_check_log_difference"):
            if not np.isfinite(row[field]) or not 0 <= row[field] <= 1e-6:
                raise ValueError(f"Bayesian numerical integration failed: {field}")
        logs = [row[f"log_bf10_{comparison}"] for comparison in COMPARISONS]
        if not np.isfinite(logs).all():
            raise ValueError("Nonfinite log BF in an ok fit; do not substitute zero")
        io.close(logs[0] - logs[1], logs[2], "log BF31 - log BF32 = log BF21")
        for comparison, log_bf in zip(COMPARISONS, logs):
            expected_label = ("strong_support_full" if log_bf >= np.log(10) else
                              "support_full" if log_bf >= np.log(3) else
                              "strong_support_reduced" if log_bf <= -np.log(10) else
                              "support_reduced" if log_bf <= -np.log(3) else "inconclusive")
            if row[f"evidence_{comparison}"] != expected_label:
                raise ValueError("Evidence label disagrees with the saved log BF")
            for field, sign in (("bf10", 1), ("bf01", -1)):
                bf = row[f"{field}_{comparison}"]
                if pd.isna(bf) or bf < 0 or (bf == 0 and sign * log_bf > -700) or (np.isinf(bf) and sign * log_bf < 700):
                    raise ValueError("Inconsistent linear/log Bayes factors")
                if np.finfo(float).tiny <= bf < math.inf:
                    io.close(np.log(bf), sign * log_bf, f"{field}/{comparison}")
            low, high = (row[f"sensitivity_log_bf10_{bound}_{comparison}"] for bound in ("min", "max"))
            if not np.isfinite([low, high]).all() or low > log_bf + 1e-8 or high < log_bf - 1e-8:
                raise ValueError("Invalid prior sensitivity bounds")
            if io.boolean(row[f"sensitivity_direction_stable_{comparison}"]) != (low > 0 or high < 0):
                raise ValueError("Prior sensitivity direction flag disagrees with its bounds")
    long = age_long_table(wide)
    saved = source / "age_bayes_comparisons.csv"
    if saved.is_file():
        other = io.read_csv(saved)
        keys = [*KEY, "comparison"]
        io.required(other, [*keys, "status", "n_obs", "bf10", "bf01", "log_bf10"], saved)
        if other.duplicated(keys).any():
            raise ValueError("Duplicate saved BF comparison rows")
        a, b = long.set_index(keys).sort_index(), other.set_index(keys).sort_index()
        if not a.index.equals(b.index) or not a.status.equals(b.status):
            raise ValueError("Stale or incomplete BF comparison export")
        for field in ("n_obs", "bf10", "bf01", "log_bf10", "sensitivity_log_bf10_min", "sensitivity_log_bf10_max"):
            io.required(other, [field], saved)
            io.close(a[field], b[field], f"wide/long {field}")
        sources.append(saved)
    prior_path = source / "age_bayes_prior_sensitivity.csv"
    if prior_path.is_file():
        prior = io.read_csv(prior_path)
        io.required(prior, [*KEY, "linear_prior_scale", "nonlinear_prior_scale", "integration_converged",
                            "quadrature_relative_error", "quadrature_check_log_difference"], prior_path)
        if prior.duplicated([*KEY, "linear_prior_scale", "nonlinear_prior_scale"]).any():
            raise ValueError("Duplicate prior sensitivity combinations")
        expected_keys = set(map(tuple, wide.loc[wide.status.eq("ok"), KEY].to_numpy()))
        if set(map(tuple, prior[KEY].to_numpy())) != expected_keys:
            raise ValueError("Prior sensitivity unit set differs from valid fits")
        for key, group in prior.groupby(KEY):
            if set(zip(group.linear_prior_scale, group.nonlinear_prior_scale)) != {(a, b) for a in (.25, .5, 1.) for b in (.25, .5, 1.)}:
                raise ValueError(f"Incomplete nine-prior sensitivity grid for {key}")
            if not group.integration_converged.map(io.boolean).eq(True).all():
                raise ValueError("Unconverged prior sensitivity integration")
            for field in ("quadrature_relative_error", "quadrature_check_log_difference"):
                if not group[field].between(0, 1e-6).all():
                    raise ValueError("Prior sensitivity numerical tolerance failed")
            original = wide.loc[wide.dimension.eq(key[0]) & wide.name.eq(key[1])].iloc[0]
            for comparison in COMPARISONS:
                values = group[f"log_bf10_{comparison}"]
                if not np.isfinite(values).all():
                    raise ValueError("Nonfinite prior log BF")
                io.close([values.min(), values.max()], [original[f"sensitivity_log_bf10_min_{comparison}"],
                                                       original[f"sensitivity_log_bf10_max_{comparison}"]], "prior summary")
                main_prior = group.loc[group.linear_prior_scale.eq(.5) & group.nonlinear_prior_scale.eq(.5)].iloc[0]
                io.close(main_prior[f"log_bf10_{comparison}"], original[f"log_bf10_{comparison}"], "main prior evidence")
        sources.append(prior_path)
    long["valid_test"] = long.status.eq("ok")
    long["invalid_reason"] = np.where(long.valid_test, "", "source_status=" + long.status.astype(str))
    return long, families, sources, bool(metadata), prior_path.is_file()


def aggregate(long, families):
    rows = []
    for dimension, expected in families.items():
        for comparison in COMPARISONS:
            group = long.loc[long.dimension.eq(dimension) & long.comparison.eq(comparison)]
            valid = group.loc[group.valid_test]
            logs = valid.log_bf10
            row = dict(dimension=dimension, analysis_scope=LABELS[dimension], expected_units=expected,
                       comparison=comparison, bf_name=BF_NAMES[comparison], n_valid=len(valid), n_invalid=len(group)-len(valid),
                       n_support_first=int(logs.ge(np.log(3)).sum()), n_support_second=int(logs.le(-np.log(3)).sum()),
                       n_inconclusive=int((logs.gt(-np.log(3)) & logs.lt(np.log(3))).sum()),
                       n_strong_first=int(logs.ge(np.log(10)).sum()), n_strong_second=int(logs.le(-np.log(10)).sum()),
                       n_prior_direction_stable=int(valid.sensitivity_direction_stable.map(io.boolean).eq(True).sum()))
            for field in ("n_obs", "bf10", "log_bf10"):
                row.update({f"{field}_min": valid[field].min(), f"{field}_max": valid[field].max()})
            rows.append(row)
    return pd.DataFrame(rows)


def run_summary(source, output=None, *, mock=False, docx=True):
    source = Path(source).resolve()
    long, families, sources, metadata_checked, priors_checked = load_results(source, mock)
    summary = aggregate(long, families)
    main = summary.loc[summary.comparison.ne(COMPARISONS[2])]
    joint = summary.loc[summary.comparison.eq(COMPARISONS[2])]
    headers = ["Analysis scope", "N", "Comparison", "Valid", "BF >= 3", "1/3 < BF < 3", "BF <= 1/3"]
    def display(frame):
        return [[f"{r.analysis_scope} ({r.expected_units})", io.interval(r.n_obs_min, r.n_obs_max, True),
                 r.bf_name, str(r.n_valid), str(r.n_support_first), str(r.n_inconclusive), str(r.n_support_second)]
                for r in frame.itertuples()]
    notes = ("BF31 compares M3 with M1; BF32 compares M3 with M2; BF21 compares M2 with M1. "
             "Counts describe evidence at the main prior scales 0.5/0.5, without BH adjustment. "
             "BF32 favoring M2 supports linearity relative to M3, not absence of age associations. "
             "Invalid fits are excluded from valid counts. Sensitivity bounds are not confidence intervals.")
    statements = []
    for comparison in COMPARISONS:
        group = summary.loc[summary.comparison.eq(comparison)]
        statements.append(f"{BF_NAMES[comparison]}：{int(group.n_valid.sum())}项有效比较中，"
                          f"{int(group.n_support_first.sum())}项BF≥3、{int(group.n_inconclusive.sum())}项证据不足、"
                          f"{int(group.n_support_second.sum())}项BF≤1/3；无效{int(group.n_invalid.sum())}项。")
    text = ("# 年龄贝叶斯因子汇报表\n\n" + io.markdown_table(headers, display(main)) +
            "\n\n## 补充线性年龄比较\n\n" + io.markdown_table(headers, display(joint)) +
            "\n\n## 数值总结\n\n" + "\n\n".join(statements) +
            "\n\nBF31：年龄函数相对无年龄；BF32：非线性相对线性；BF21：线性年龄相对无年龄。"
            "支持第一个或第二个模型的数量不是FDR显著脑区数；不平均各脑区BF，不从BF大小推断关联方向。"
            "完整明细保留各单元BF、log BF、先验敏感性和可匹配的原频率学统计。\n\n" + notes +
            f"\n\n来源metadata核验：{metadata_checked}；九套先验明细核验：{priors_checked}。未提供的来源检查不能视为通过。")
    frames = {"bf_main_summary.csv": main, "bf_linear_summary.csv": joint, "bf_all_summary.csv": summary,
              "bf_primary_details.csv": long.loc[long.comparison.ne(COMPARISONS[2])],
              "bf_linear_details.csv": long.loc[long.comparison.eq(COMPARISONS[2])],
              "bf_invalid_tests.csv": long.loc[~long.valid_test]}
    sections = ("Bayesian evidence for age associations", "Completed age models with shared covariates and prespecified priors",
                [("Main age comparisons", headers, display(main), [1.8, .6, .75, .5, .8, 1.35, .9], notes),
                 ("Supplementary linear age comparison", headers, display(joint), [1.8, .6, .75, .5, .8, 1.35, .9],
                  "BF21 describes the linear age model relative to the model without age.")])
    return io.export(source, output or source / "reviewer_summary", "age_bf_reviewer_summary_v1", frames, text,
                     sections, sources, mock=mock, docx=docx,
                     extra={"source_metadata_checked": metadata_checked, "prior_grid_checked": priors_checked,
                            "n_invalid_units": int((~long.valid_test).sum() // 3), "atlas_parcels": families["parcel"]})


COMPARISON_LABELS = {
    "m3_vs_m1": "Age spline vs. no age",
    "m3_vs_m2": "Age spline vs. linear age",
    "m2_vs_m1": "Linear age vs. no age",
}
SELECTED_FIELDS = dict(zip(COMPARISONS, ("age_spline_vs_no_age", "age_spline_vs_linear_age", "linear_age_vs_no_age")))


def run_selected_summary(source, parcels, output=None, *, mock=False, docx=True):
    """Report individual highlighted Tian432 parcels without altering all-unit inference."""
    source = Path(source).resolve()
    if not parcels or len(parcels) != len(set(parcels)):
        raise ValueError("Choose at least one parcel, without duplicate IDs")
    long, families, sources, metadata_checked, priors_checked = load_results(source, mock)
    if families["parcel"] != 432:
        raise ValueError("Selected reporting uses Schaefer400/Tian432 IDs; do not reuse these IDs in the 232 atlas")
    if any(not 1 <= parcel <= 432 for parcel in parcels):
        raise ValueError("Selected parcel IDs must be inside 1..432")
    selected = long.loc[long.dimension.eq("parcel") & long.name.isin(list(map(str, parcels)))].copy()
    selected["model_comparison"] = selected.comparison.map(COMPARISON_LABELS)
    selected["_order"] = selected.name.map({str(parcel): i for i, parcel in enumerate(parcels)})
    selected = selected.sort_values(["_order", "comparison"], kind="stable").drop(columns="_order")
    records = []
    for parcel in parcels:
        group = selected.loc[selected.name.eq(str(parcel))].set_index("comparison")
        first = group.loc[COMPARISONS[0]]
        record = dict(parcel=parcel, label=first.label, n_obs=first.n_obs, status=first.status,
                      valid_unit=bool(group.valid_test.all()))
        for comparison, suffix in SELECTED_FIELDS.items():
            for field in ("bf10", "log_bf10", "sensitivity_log_bf10_min", "sensitivity_log_bf10_max"):
                record[f"{field}_{suffix}"] = group.at[comparison, field]
        records.append(record)
    compact = pd.DataFrame(records)
    headers = ["Parcel", "Region", "N", *COMPARISON_LABELS.values()]
    def bf_display(value, log_value):
        return f"log BF = {io.number(log_value)}" if value == 0 or np.isinf(value) else io.number(value)
    display = []
    for row in compact.to_dict("records"):
        display.append([str(row["parcel"]), str(row["label"]).removeprefix("7Networks_"),
                        "NA" if pd.isna(row["n_obs"]) else str(int(row["n_obs"])),
                        *(bf_display(row[f"bf10_{suffix}"], row[f"log_bf10_{suffix}"]) if row["valid_unit"] else "NA"
                          for suffix in SELECTED_FIELDS.values())])
    valid = compact.loc[compact.valid_unit]
    sample = io.interval(valid.n_obs.min(), valid.n_obs.max(), True)
    english = [f"We examined Bayesian evidence for age associations in {len(parcels)} parcels highlighted by the preceding regional analyses "
               f"(parcel IDs: {', '.join(map(str, parcels))}; N = {sample} across valid parcels)."]
    comparisons_en = {
        "m3_vs_m1": ("the age spline model", "the model without age terms"),
        "m3_vs_m2": ("the age spline model", "the linear age model"),
        "m2_vs_m1": ("the linear age model", "the model without age terms"),
    }
    for comparison, suffix in SELECTED_FIELDS.items():
        first, second = comparisons_en[comparison]
        groups = {"first": [], "second": [], "inconclusive": []}
        for row in valid.to_dict("records"):
            log_bf = row[f"log_bf10_{suffix}"]
            category = "first" if log_bf >= np.log(3) else "second" if log_bf <= -np.log(3) else "inconclusive"
            value = row[f"bf10_{suffix}"]
            estimate = bf_display(value, log_bf)
            groups[category].append(f"{row['parcel']} ({estimate if value == 0 or np.isinf(value) else 'BF = ' + estimate})")
        parts = []
        if groups["first"]:
            parts.append(f"evidence favored {first} in parcels {', '.join(groups['first'])}")
        if groups["second"]:
            parts.append(f"evidence favored {second} in parcels {', '.join(groups['second'])}")
        if groups["inconclusive"]:
            parts.append(f"evidence was inconclusive in parcels {', '.join(groups['inconclusive'])}")
        if parts:
            english.append(f"When comparing {first} with {second}, " + "; ".join(parts) + ".")
    invalid = compact.loc[~compact.valid_unit, "parcel"].tolist()
    if invalid:
        english.append(f"Bayesian comparisons were not estimable for parcels {', '.join(map(str, invalid))}.")
    english.append("These results are follow-up characterizations of selected parcels, rather than an independent whole-brain confirmation.")
    notes = ("Entries are Bayes factors for the first model relative to the second. All models include the same covariates. "
             "The spline age model allows both linear and nonlinear associations. BF >= 3 favors the first model; BF <= 1/3 favors the second; "
             "intermediate values are inconclusive. Main prior scales are 0.5/0.5. Selected reporting does not change priors or recompute any BF. "
             "These parcels were highlighted in earlier analyses; selection is not an independent confirmation. Invalid fits are displayed as NA. "
             "The sensitivity ranges in CSV summarize priors, not confidence intervals. See log BF if a linear-scale value overflows.")
    text = ("# 关注脑区的逐区贝叶斯结果\n\n" + io.markdown_table(headers, display) +
            "\n\n## 英文结果段落\n\n" + " ".join(english) + "\n\n## 读取说明\n\n" + notes +
            "\n\n本表为既有432版结果的摘选，不重新计算BF，不在6区内重新校正。"
            "不把156、157自动称为显著或边缘显著，不把其余4区预设为BF支持年龄关联。"
            "具体支持方向只由各区真实BF决定；BF支持线性相对样条，不等于没有年龄关联。")
    frames = {"selected_region_summary.csv": compact, "selected_region_comparisons.csv": selected,
              "selected_invalid_tests.csv": selected.loc[~selected.valid_test]}
    sections = ("Bayesian evidence in highlighted brain parcels", "Individual regional estimates from the completed 432 parcel analysis",
                [("Regional model comparisons", headers, display, [.45, 1.8, .45, 1.5, 1.5, 1.15], notes)])
    return io.export(source, output or source / "reviewer_selected_parcels", "age_bf_selected_summary_v1", frames, text,
                     sections, sources, mock=mock, docx=docx,
                     extra={"selected_parcels": list(parcels), "source_metadata_checked": metadata_checked,
                            "prior_grid_checked": priors_checked, "n_invalid_units": len(invalid), "atlas_parcels": 432})


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--analysis-dir", type=Path, required=True, help="Completed bayes/ directory")
    parser.add_argument("--output", type=Path)
    parser.add_argument("--skip-docx", action="store_true")
    parser.add_argument("--mock", action="store_true")
    parser.add_argument("--parcels", type=int, nargs="+", help="Report individual Schaefer400/Tian432 parcels in one table")
    args = parser.parse_args()
    source, output = io.project_path(args.analysis_dir), io.project_path(args.output) if args.output else None
    if args.parcels is not None:
        run_selected_summary(source, args.parcels, output, mock=args.mock, docx=not args.skip_docx)
    else:
        run_summary(source, output, mock=args.mock, docx=not args.skip_docx)


if __name__ == "__main__":
    main()
