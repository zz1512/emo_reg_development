"""Extract reviewer R1-6 tables from saved ERQ v3 results, without refitting.

Only pandas is required for CSV/Markdown; Word additionally needs python-docx.
The saved p/q values are copied, never recalculated or corrected again.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path

import pandas as pd

from analysis.fig03.gam_pygam import result_tables as tables

ROOT = Path(__file__).resolve().parents[3]
KIND = "erq_reviewer_summary_v1"
ALPHA = .05
FAMILIES = {"parcel": ("Whole-brain parcels", 432), "yeo": ("Yeo networks", 7),
            "system": ("Systems", 4), "amygdala": ("Amygdala subdivisions", 6)}
PRIMARY = ("additive_vs_baseline", "full_vs_additive")
COMPARISONS = (*PRIMARY, "full_vs_baseline")
NUMERIC = ("n_obs", "beta", "se", "ci_low", "ci_high", "t", "F", "df1", "df2", "p", "q", "partial_r2_reference")
KEY = ["dimension", "name", "comparison"]


def _boolean(value):
    if pd.isna(value) or str(value).strip() == "":
        return None
    text = str(value).strip().casefold()
    if text in ("true", "1", "1.0"):
        return True
    if text in ("false", "0", "0.0"):
        return False
    raise ValueError(f"Unrecognized reject value {value!r}; expected true/false or 1/0")


def _read(path):
    return pd.read_csv(path, dtype={"name": str}, float_precision="round_trip", encoding="utf-8-sig")


def load_results(analysis_dir: Path, mock=False):
    directory = analysis_dir / "statistics"
    metadata_path = analysis_dir / "run_metadata.json"
    metadata = json.loads(metadata_path.read_text(encoding="utf-8")) if metadata_path.is_file() else {}
    if metadata and metadata.get("analysis_kind") != "erq_age_smooth_pygam_v3":
        raise ValueError("ERQ three-model v3 results are required; old full-vs-baseline results cannot supply the two separate effects")
    sources = [metadata_path] if metadata else []
    wide_path, long_path = directory / "erq_age_statistics.csv", directory / "erq_sequential_comparisons.csv"
    wide = None
    if wide_path.is_file():
        wide_source = _read(wide_path)
        required = {"dimension", "name", "label", "status", "n_obs", "beta_erq_additive_reference"}
        required.update(f"{field}_{comparison}" for comparison in COMPARISONS
                        for field in ("F", "df1", "df2", "p", "q", "reject", "partial_r2_reference"))
        required.update(f"{field}_erq_additive_reference" for field in ("se", "ci_low", "ci_high", "t"))
        if missing := required - set(wide_source):
            raise ValueError(f"Incomplete/old ERQ wide table; missing columns: {sorted(missing)}")
        wide = tables.erq_long_table(wide_source)
        sources.append(wide_path)
    if long_path.is_file():
        frame = _read(long_path)
        sources.append(long_path)
    elif wide is not None:
        frame = wide.copy()
    else:
        raise FileNotFoundError(f"No v3 ERQ result tables under {directory}")
    required = {*KEY, "label", "status", "reject", *NUMERIC}
    if missing := required - set(frame):
        raise ValueError(f"Incomplete ERQ long table; missing columns: {sorted(missing)}")
    # Catch partially updated/stale exports instead of silently combining runs.
    if wide is not None and long_path.is_file():
        for source in (frame, wide):
            if source.duplicated(KEY).any():
                raise ValueError("Duplicate dimension/name/comparison rows in ERQ results")
        left, right = frame.set_index(KEY).sort_index(), wide.set_index(KEY).sort_index()
        if not left.index.equals(right.index):
            raise ValueError("Wide and long tables contain different analysis units/comparisons")
        for field in (*NUMERIC, "status", "reject", "label"):
            for a, b in zip(left[field], right[field]):
                if pd.isna(a) and pd.isna(b):
                    continue
                equal = (math.isclose(float(a), float(b), rel_tol=1e-9, abs_tol=1e-12)
                         if field in NUMERIC else _boolean(a) == _boolean(b) if field == "reject" else str(a) == str(b))
                if not equal:
                    raise ValueError(f"Wide/long result mismatch in {field}; regenerate the v3 table export first")
    definition_path = directory / "erq_table_definitions.json"
    definitions = json.loads(definition_path.read_text(encoding="utf-8")) if definition_path.is_file() else {}
    if definition_path.is_file():
        sources.append(definition_path)
    simulated = (metadata.get("simulated", False) or definitions.get("simulated", False)
                 or frame["label"].astype(str).str.contains("SIMULATED", case=False).any())
    if simulated and not mock:
        raise ValueError("Source is explicitly simulated; pass --mock to produce a labelled demonstration only")
    return frame, sources, metadata


def validate_results(frame: pd.DataFrame):
    frame = frame.copy()
    if frame.duplicated(KEY).any():
        raise ValueError("Duplicate dimension/name/comparison rows")
    if frame[[*KEY, "label", "status"]].isna().any().any():
        raise ValueError("Missing unit identifiers, labels or statuses")
    if set(frame["dimension"]) != set(FAMILIES) or set(frame["comparison"]) != set(COMPARISONS):
        raise ValueError("All four analysis families and all three v3 comparisons are required")
    for dimension, (_, expected) in FAMILIES.items():
        unit_sets = [set(frame.loc[(frame["dimension"] == dimension) & (frame["comparison"] == comparison), "name"])
                     for comparison in COMPARISONS]
        if any(len(units) != expected or units != unit_sets[0] for units in unit_sets):
            raise ValueError(f"{dimension}: expected the same {expected} units in each comparison, including failed units")
    if set(frame.loc[frame["dimension"] == "parcel", "name"]) != {str(i) for i in range(1, 433)}:
        raise ValueError("Parcel IDs must be exactly 1..432")
    frame["reject"] = frame["reject"].map(_boolean)
    for column in NUMERIC:
        original = frame[column]
        converted = pd.to_numeric(original, errors="coerce")
        bad = original.notna() & original.astype(str).str.strip().ne("") & converted.isna()
        if bad.any():
            raise ValueError(f"Non-numeric values in {column}: {original[bad].head(3).tolist()}")
        frame[column] = converted
    counts = frame.groupby(["dimension", "name"])["n_obs"].nunique()
    if counts.gt(1).any():
        raise ValueError("The three comparisons have different n_obs for the same unit")
    reasons = []
    for row in frame.to_dict("records"):
        failures = []
        if row["status"] != "ok":
            failures.append(f"status={row['status']}")
        for column in ("n_obs", "F", "df1", "df2", "p", "q", "partial_r2_reference"):
            if not math.isfinite(row[column]):
                failures.append(f"{column}_missing_or_nonfinite")
        for column in ("p", "q", "partial_r2_reference"):
            if math.isfinite(row[column]) and not 0 <= row[column] <= 1:
                raise ValueError(f"{column} outside [0,1] for {row['dimension']}/{row['name']}")
        for column in ("n_obs", "df1", "df2"):
            if math.isfinite(row[column]) and (row[column] <= 0 or row[column] != int(row[column])):
                failures.append(f"{column}_not_positive_integer")
        if math.isfinite(row["F"]) and row["F"] < 0:
            raise ValueError("F cannot be negative")
        if pd.isna(row["reject"]):
            failures.append("reject_missing")
        elif math.isfinite(row["q"]) and row["reject"] != (row["q"] < ALPHA):
            raise ValueError(f"reject/q inconsistency for {row['dimension']}/{row['name']}/{row['comparison']}")
        if row["comparison"] == PRIMARY[0]:
            for column in ("beta", "se", "ci_low", "ci_high", "t"):
                if not math.isfinite(row[column]):
                    failures.append(f"{column}_missing_or_nonfinite")
            if row["se"] <= 0 or row["ci_low"] > row["ci_high"]:
                failures.append("invalid_coefficient_interval")
        reasons.append(";".join(failures))
    frame["invalid_reason"] = reasons
    frame["valid_test"] = frame["invalid_reason"].eq("")
    frame["significant_from_q"] = frame["valid_test"] & frame["q"].lt(ALPHA)
    order = {family: i for i, family in enumerate(FAMILIES)}
    frame["_family_order"] = frame["dimension"].map(order)
    frame["_comparison_order"] = frame["comparison"].map({test: i for i, test in enumerate(COMPARISONS)})
    frame["_unit_order"] = pd.to_numeric(frame["name"], errors="coerce")
    return frame.sort_values(["_family_order", "_unit_order", "name", "_comparison_order"], kind="stable").drop(
        columns=["_family_order", "_unit_order", "_comparison_order"]).reset_index(drop=True)


def aggregate(frame):
    rows = []
    for dimension, (label, expected) in FAMILIES.items():
        for comparison in COMPARISONS:
            subset = frame.loc[(frame["dimension"] == dimension) & (frame["comparison"] == comparison)]
            valid = subset.loc[subset["valid_test"]]
            row = {"dimension": dimension, "analysis_scope": label, "expected_units": expected,
                   "comparison": comparison, "n_units": len(subset), "n_valid": len(valid),
                   "n_invalid": len(subset) - len(valid), "n_significant": int(valid["significant_from_q"].sum()),
                   "n_unadjusted_p_lt_0_05": int(valid["p"].lt(ALPHA).sum())}
            for field in ("n_obs", "p", "q", "F", "df1", "df2", "partial_r2_reference"):
                row[f"{field}_min"], row[f"{field}_max"] = valid[field].min(), valid[field].max()
            rows.append(row)
    summary = pd.DataFrame(rows)
    compact = []
    for dimension, (label, expected) in FAMILIES.items():
        row = {"dimension": dimension, "analysis_scope": label, "expected_units": expected}
        valid = frame.loc[(frame["dimension"] == dimension) & frame["comparison"].isin(PRIMARY) & frame["valid_test"]]
        row.update(n_obs_min=valid["n_obs"].min(), n_obs_max=valid["n_obs"].max())
        for comparison, prefix in zip(PRIMARY, ("erq_main", "interaction")):
            source = summary.loc[(summary["dimension"] == dimension) & (summary["comparison"] == comparison)].iloc[0]
            for field in ("n_units", "n_valid", "n_invalid", "n_significant", "n_unadjusted_p_lt_0_05",
                          "n_obs_min", "n_obs_max", "p_min", "p_max", "q_min", "q_max",
                          "partial_r2_reference_min", "partial_r2_reference_max"):
                row[f"{prefix}_{field}"] = source[field]
        compact.append(row)
    return pd.DataFrame(compact), summary


def _number(value):
    if pd.isna(value):
        return "NA"
    return "<0.001" if 0 <= value < .001 else f"{value:.3f}"


def _range(low, high, integer=False):
    if pd.isna(low) or pd.isna(high):
        return "NA"
    render = (lambda value: str(int(value))) if integer else _number
    return render(low) if low == high else f"{render(low)} to {render(high)}"


def display_rows(compact):
    result = []
    for row in compact.to_dict("records"):
        result.append([f"{row['analysis_scope']} ({row['expected_units']})",
                       _range(row["n_obs_min"], row["n_obs_max"], True),
                       _range(row["erq_main_q_min"], row["erq_main_q_max"]),
                       f"{row['erq_main_n_significant']}/{row['erq_main_n_valid']}",
                       _range(row["interaction_q_min"], row["interaction_q_max"]),
                       f"{row['interaction_n_significant']}/{row['interaction_n_valid']}"])
    return result


def result_text(summary):
    parts, chinese = [], []
    for comparison, description, zh in ((PRIMARY[0], "constant ERQ associations", "恒定ERQ关联"),
                                         (PRIMARY[1], "age by ERQ interactions", "年龄×ERQ交互"),
                                         (COMPARISONS[2], "joint ERQ and interaction contributions", "ERQ与交互的联合贡献")):
        group = summary.loc[summary["comparison"] == comparison]
        valid, significant, invalid = (int(group[field].sum()) for field in ("n_valid", "n_significant", "n_invalid"))
        minimum = group["q_min"].min()
        if not valid:
            parts.append(f"No valid tests were available for {description}.")
            chinese.append(f"{zh}没有有效检验，不能判定是否显著。")
        elif not significant:
            parts.append(f"For {description}, none of the {valid} valid tests survived BH correction (minimum q = {_number(minimum)}).")
            chinese.append(f"{zh}的{valid}项有效检验中，0项通过BH校正（最小q={_number(minimum)}）。")
        else:
            parts.append(f"For {description}, {significant} of {valid} valid tests survived BH correction (minimum q = {_number(minimum)}).")
            chinese.append(f"{zh}的{valid}项有效检验中，{significant}项通过BH校正（最小q={_number(minimum)}）。")
        if invalid:
            parts.append(f"A further {invalid} tests were not estimable or failed validation and were not counted as nonsignificant.")
            chinese.append(f"另有{invalid}项不可计算或未通过结果校验，不计作不显著。")
    common = "BH correction was performed separately for each comparison within the parcel432, Yeo7, system4 and Juelich6 families."
    caution = "Nonsignificant results do not establish absence of associations, developmental independence, or mediation of age-related neural differences."
    return " ".join([common, *parts, caution]), "\n".join(chinese)


NOTES = ("E1 includes smooth age and gender, categorical site, parental education and mean FD. "
         "E2 adds centered ERQ; E3 adds the smooth ERQ by age interaction. The main comparison is E2 versus E1; "
         "the interaction comparison is E3 versus E2. Each comparison was BH corrected separately within each family. "
         "Sig/valid gives q < .05 tests divided by valid tests; invalid results are not nonsignificant. "
         "N is the common sample size or its range across valid units. Ranges summarize units, not confidence intervals. "
         "Display values are rounded; counts use unrounded q. Coefficients and inference use the fixed-spline OLS reference, "
         "not the independently GCV-selected penalized fits.")


def write_word(compact, summary, output, mock):
    document = tables._document("ERQ associations with neural activation",
                                "Four analysis families with all prespecified units retained", mock)
    document.add_paragraph("Main comparisons").runs[0].bold = True
    headers = ["Analysis scope", "N", "ERQ main\nq range", "Sig/valid", "Interaction\nq range", "Sig/valid"]
    table = tables._three_line_table(document, headers, [1.55, .55, 1.1, .7, 1.1, .7])
    for row in display_rows(compact):
        tables._append_row(table, row)
    document.add_paragraph("Notes  " + NOTES)
    document.add_paragraph("Supplementary joint comparison").runs[0].bold = True
    table = tables._three_line_table(document, ["Analysis scope", "N", "Joint q range", "Sig/valid"], [2.5, .7, 1.5, 1.1])
    for row in summary.loc[summary["comparison"] == COMPARISONS[2]].to_dict("records"):
        tables._append_row(table, [f"{row['analysis_scope']} ({row['expected_units']})",
                                  _range(row["n_obs_min"], row["n_obs_max"], True),
                                  _range(row["q_min"], row["q_max"]), f"{row['n_significant']}/{row['n_valid']}"])
    document.add_paragraph("The supplementary E3 versus E1 comparison tests ERQ and its age interaction jointly; it is not the ERQ main effect alone.")
    invalid = int(summary["n_invalid"].sum())
    if invalid:
        document.add_paragraph(f"Nonestimable or invalid test rows: {invalid}. See erq_invalid_tests.csv before making an all-units claim.")
    document.save(output / "erq_response_table.docx")


def run_summary(analysis_dir: Path, output: Path | None = None, *, include_docx=True,
                expect_all_nonsignificant=False, mock=False):
    analysis_dir = analysis_dir.resolve()
    output = (output or analysis_dir / "reviewer_R1_6").resolve()
    statistics_dir = analysis_dir / "statistics"
    if output == analysis_dir or output in analysis_dir.parents or output == statistics_dir or statistics_dir in output.parents:
        raise ValueError("Use a separate reviewer summary directory, outside the source statistics directory")
    manifest = output / "summary_metadata.json"
    previous = json.loads(manifest.read_text(encoding="utf-8")) if manifest.is_file() else {}
    if output.exists() and any(output.iterdir()) and previous.get("analysis_kind") != KIND:
        raise ValueError("Output is nonempty and is not an owned ERQ reviewer-summary directory")
    frame, sources, metadata = load_results(analysis_dir, mock)
    frame = validate_results(frame)
    # Even an invalid row with finite significant q contradicts an all-false assertion.
    if expect_all_nonsignificant and frame["q"].lt(ALPHA).any():
        examples = frame.loc[frame["q"].lt(ALPHA), KEY + ["q"]].head().to_dict("records")
        raise ValueError(f"--expect-all-nonsignificant contradicted by q < .05: {examples}")
    compact, summary = aggregate(frame)
    if include_docx:
        tables.check_docx_dependency()
    english, chinese = result_text(summary)
    output.mkdir(parents=True, exist_ok=True)
    exports = {"erq_response_summary.csv": compact,
               "erq_comparison_summary.csv": summary,
               "erq_joint_summary.csv": summary.loc[summary["comparison"] == COMPARISONS[2]],
               "erq_primary_details.csv": frame.loc[frame["comparison"].isin(PRIMARY)],
               "erq_joint_details.csv": frame.loc[frame["comparison"] == COMPARISONS[2]],
               "erq_invalid_tests.csv": frame.loc[~frame["valid_test"]]}
    for filename, data in exports.items():
        data.to_csv(output / filename, index=False, float_format="%.17g", encoding="utf-8-sig")
    if include_docx:
        write_word(compact, summary, output, mock)
    elif "erq_response_table.docx" in previous.get("outputs", []) and (output / "erq_response_table.docx").is_file():
        (output / "erq_response_table.docx").unlink()  # Only this exporter-owned stale display table.
    lines = ["# ERQ审稿回复结果汇总", ""]
    if mock:
        lines += ["**SIMULATED DATA ONLY：以下均为模拟结果，不能用于论文填写。**", ""]
    lines += ["## 回复信主表", "", "建议4行数据，不含表头；每范围一行，同时列两项主要比较。", "",
              "| Analysis scope | N | ERQ main q range | Sig/valid | Interaction q range | Sig/valid |",
              "|---|---|---|---|---|---|"]
    lines += ["| " + " | ".join(map(str, row)) + " |" for row in display_rows(compact)]
    lines += ["", "联合比较另见erq_joint_summary.csv（4行），不能称为单独ERQ主效应。",
              "主效应和交互明细898行；联合明细449行。失败单元保留，不能填成0效应。", "",
              "## 数值总结", "", chinese, "", "## 英文结果段落", "", english, "",
              "## 填写与解释", "", "不显著仍报告有效单元数、N、q范围和显著数；不将全部449单元平均成一个beta、F或p。",
              "详细beta、SE、95%CI、t/F、自由度、p、q及partial R²保存在两份明细表，不因不显著而删除。",
              "主效应beta来自E2；平滑交互没有单个beta。q区间是跨单元的最小到最大值，不是置信区间。",
              "若有无效检验，仅能说有效单元中未发现校正后显著结果；不对无效单元作阴性结论。",
              "只读已有结果，不重新拟合、不重新BH；未校正p<0.05不等于校正后显著。", "", NOTES, "",
              "完整来源、哈希、行数和校验状态见summary_metadata.json。"]
    (output / "erq_result_summary.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    if not metadata:
        lines_note = "No run_metadata.json found; v3 table schema and all 449 units were validated"
    else:
        lines_note = "Source run_metadata.json identified ERQ v3"
    result = {"analysis_kind": KIND, "simulated": mock, "alpha": ALPHA, "analysis_dir": str(analysis_dir),
              "source_metadata_check": lines_note,
              "source_sha256": {str(path): hashlib.sha256(path.read_bytes()).hexdigest() for path in sources},
              "all_tests_valid": bool(frame["valid_test"].all()), "n_invalid_tests": int((~frame["valid_test"]).sum()),
              "all_available_q_nonsignificant": not bool(frame["q"].lt(ALPHA).any()),
              "all_tests_valid_and_nonsignificant": bool(frame["valid_test"].all()) and not bool(frame["q"].lt(ALPHA).any()),
              "row_counts": {filename: len(data) for filename, data in exports.items()},
              "source_changes": "none; no model fitting or BH recalculation",
              "outputs": [*exports, "erq_result_summary.md", "summary_metadata.json"] + (["erq_response_table.docx"] if include_docx else [])}
    manifest.write_text(json.dumps(result, indent=2, ensure_ascii=False), encoding="utf-8")
    print(chinese, flush=True)
    print(f"Tables and reading notes: {output}", flush=True)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--analysis-dir", type=Path, required=True, help="ERQ v3 output root, containing statistics/")
    parser.add_argument("--output", type=Path, help="Default: ANALYSIS_DIR/reviewer_R1_6")
    parser.add_argument("--skip-docx", action="store_true")
    parser.add_argument("--expect-all-nonsignificant", action="store_true", help="Assert no saved q < .05 across all three comparisons")
    parser.add_argument("--mock", action="store_true", help="Label all outputs as simulated demonstrations")
    args = parser.parse_args()
    resolve = lambda path: (path if path.is_absolute() else ROOT / path).resolve()
    run_summary(resolve(args.analysis_dir), resolve(args.output) if args.output else None,
                include_docx=not args.skip_docx, expect_all_nonsignificant=args.expect_all_nonsignificant, mock=args.mock)


if __name__ == "__main__":
    main()
