"""Full-precision CSV and editable Word three-line tables, without refitting.

Standalone export deliberately needs only pandas and python-docx, allowing
analysis and document rendering to use separate Python environments.
"""

from __future__ import annotations

import argparse
import importlib.util
import json
import math
from pathlib import Path

import pandas as pd

FAMILIES = {"parcel": "Parcel", "yeo": "Yeo network", "system": "System", "amygdala": "Juelich amygdala"}
EFFECTS = (("additive_vs_baseline", "ERQ main effect", "E2 vs E1"),
           ("full_vs_additive", "Age by ERQ", "E3 vs E2"),
           ("full_vs_baseline", "Joint ERQ", "E3 vs E1"))


def check_docx_dependency():
    if importlib.util.find_spec("docx") is None:
        raise ImportError("Word export requires python-docx: python -m pip install python-docx; "
                          "or use --skip-docx and export saved CSV with result_tables later")


def erq_long_table(statistics: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for row in statistics.to_dict("records"):
        for comparison, effect, models in EFFECTS:
            main = comparison == "additive_vs_baseline"
            record = {key: row.get(key) for key in ("dimension", "name", "label", "status", "n_obs")}
            record.update(comparison=comparison, effect=effect, models=models,
                          estimate_source="fixed-spline OLS reference", beta_units="contrast per ERQ score point" if main else "",
                          statistic_type="t" if main else "F")
            for field in ("beta", "se", "ci_low", "ci_high", "t"):
                record[field] = row.get(f"{field}_erq_additive_reference", float("nan")) if main else float("nan")
            for field in ("F", "df1", "df2", "p", "q", "reject", "partial_r2_reference", "rss_reference_reduced"):
                record[field] = row.get(f"{field}_{comparison}", float("nan"))
            record["statistic"] = record["t"] if main else record["F"]
            record["delta_r2_pygam"] = row.get(f"delta_r2_{comparison}", float("nan"))
            record["delta_rss_pygam"] = row.get(f"delta_rss_{comparison}", float("nan"))
            record["rss_reference_target"] = row.get("rss_reference_additive" if main else "rss_reference_full", float("nan"))
            rows.append(record)
    return pd.DataFrame(rows)


def age_long_table(statistics: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for row in statistics.to_dict("records"):
        for comparison in ("m3_vs_m1", "m3_vs_m2", "m2_vs_m1"):
            record = {key: row.get(key) for key in ("dimension", "name", "label", "status", "n_obs")}
            record["comparison"] = comparison
            for field in ("bf10", "bf01", "log_bf10", "evidence", "sensitivity_log_bf10_min",
                          "sensitivity_log_bf10_max", "sensitivity_direction_stable", "sensitivity_evidence_labels"):
                record[field] = row.get(f"{field}_{comparison}", float("nan"))
            for field in ("F", "df1", "df2", "p", "q", "partial_r2"):
                record[f"frequentist_{field}"] = row.get(f"{field}_{comparison}", float("nan"))
            rows.append(record)
    return pd.DataFrame(rows)


def _number(value, probability=False):
    try:
        value = float(value)
    except (ValueError, TypeError):
        return "—"
    if math.isnan(value):
        return "—"
    if math.isinf(value):
        return "overflow; see log BF"
    if probability:
        return "<.001" if value < .001 else f"{value:.3f}"
    if value != 0 and (abs(value) >= 10000 or abs(value) < .001):
        return f"{value:.2e}"
    return f"{value:.3f}"


def _document(title, note, mock):
    check_docx_dependency()
    from docx import Document
    from docx.shared import Inches, Pt, RGBColor
    from docx.oxml import OxmlElement
    from docx.oxml.ns import qn

    document = Document()
    section = document.sections[0]
    section.page_width, section.page_height = Inches(8.5), Inches(11)
    section.top_margin = section.bottom_margin = Inches(.55)
    section.left_margin = section.right_margin = Inches(.55)
    normal = document.styles["Normal"]
    normal.font.name, normal.font.size = "Times New Roman", Pt(10)
    normal.paragraph_format.space_after = Pt(4)
    document.styles["Title"].font.name = "Times New Roman"
    document.styles["Title"].font.size = Pt(17)
    document.styles["Title"].font.color.rgb = RGBColor(0, 0, 0)
    # Bundled/default Word templates may decorate Title with an accent rule
    # and theme fonts. Explicitly remove that residue for journal tables.
    for style_name in ("Normal", "Title", "Title Char"):
        if style_name not in document.styles:
            continue
        style = document.styles[style_name]
        ppr = style.element.find(qn("w:pPr"))
        if ppr is not None:
            for border in list(ppr.findall(qn("w:pBdr"))):
                ppr.remove(border)
        rpr = style.element.get_or_add_rPr()
        fonts = rpr.find(qn("w:rFonts"))
        if fonts is not None:
            for attribute in list(fonts.attrib):
                if "theme" in attribute.lower():
                    del fonts.attrib[attribute]
            for attribute in ("ascii", "hAnsi", "eastAsia", "cs"):
                fonts.set(qn(f"w:{attribute}"), "Times New Roman")
    document.add_paragraph(title, style="Title")
    if mock:
        document.add_paragraph("SIMULATED DATA  Method and layout demonstration only").runs[0].bold = True
    document.add_paragraph(note)
    footer = section.footer.paragraphs[0]
    footer.alignment = 2
    footer.add_run("Page ")
    field = OxmlElement("w:fldSimple")
    field.set(qn("w:instr"), "PAGE")
    footer._p.append(field)
    return document


def _three_line_table(document, headers, widths):
    from docx.shared import Inches, Pt
    from docx.oxml import OxmlElement
    from docx.oxml.ns import qn

    table = document.add_table(rows=1, cols=len(headers))
    table.autofit = False
    for column, width in zip(table.columns, widths):
        column.width = Inches(width)
    borders = OxmlElement("w:tblBorders")
    for edge in ("top", "bottom", "left", "right", "insideH", "insideV"):
        element = OxmlElement(f"w:{edge}")
        element.set(qn("w:val"), "single" if edge in ("top", "bottom") else "nil")
        element.set(qn("w:sz"), "8")
        element.set(qn("w:color"), "000000")
        borders.append(element)
    table._tbl.tblPr.append(borders)
    margins = OxmlElement("w:tblCellMar")
    for edge, width in (("top", 30), ("bottom", 30), ("left", 70), ("right", 70)):
        margin = OxmlElement(f"w:{edge}")
        margin.set(qn("w:w"), str(width))
        margin.set(qn("w:type"), "dxa")
        margins.append(margin)
    table._tbl.tblPr.append(margins)
    repeat = OxmlElement("w:tblHeader")
    table.rows[0]._tr.get_or_add_trPr().append(repeat)
    for cell, text, width in zip(table.rows[0].cells, headers, widths):
        cell.width = Inches(width)
        cell.text = text
        p = cell.paragraphs[0]
        p.alignment = 1
        p.paragraph_format.keep_with_next = True
        p.paragraph_format.space_after = Pt(4)
        for run in p.runs:
            run.bold = True
        border = OxmlElement("w:tcBorders")
        bottom = OxmlElement("w:bottom")
        bottom.set(qn("w:val"), "single")
        bottom.set(qn("w:sz"), "4")
        border.append(bottom)
        cell._tc.get_or_add_tcPr().append(border)
    return table


def _append_row(table, texts, keep=False):
    from docx.shared import Pt
    from docx.oxml import OxmlElement
    from docx.enum.table import WD_CELL_VERTICAL_ALIGNMENT

    row = table.add_row()
    row._tr.get_or_add_trPr().append(OxmlElement("w:cantSplit"))
    for index, (cell, text) in enumerate(zip(row.cells, texts)):
        cell.text = str(text)
        cell.vertical_alignment = WD_CELL_VERTICAL_ALIGNMENT.CENTER
        p = cell.paragraphs[0]
        p.alignment = 0 if index == 0 else 1
        p.paragraph_format.space_before = Pt(3)
        p.paragraph_format.space_after = Pt(3)
        p.paragraph_format.keep_with_next = keep
        for run in p.runs:
            run.font.size = Pt(9)
    return row


def write_erq_tables(statistics: pd.DataFrame, output_dir: Path, include_docx=True, mock=False):
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    long = erq_long_table(statistics)
    long.to_csv(output_dir / "erq_sequential_comparisons.csv", index=False, float_format="%.17g")
    counts = {dimension: {test: int(((statistics["dimension"] == dimension) & statistics[f"p_{test}"].notna()).sum())
                          for test, _, _ in EFFECTS} for dimension in FAMILIES}
    definition = {
        "units": len(statistics), "table_rows": len(long), "bh_valid_tests": counts,
        "inference": "Fixed five-spline-space OLS reference; independent homoskedastic Gaussian errors",
        "main_effect": "E2 ERQ coefficient, contrast units per one raw ERQ score point; t test equals E2-vs-E1 F with one numerator df",
        "interaction": "E3-vs-E2 block F test; no scalar interaction beta exists for g(age)",
        "joint": "E3-vs-E1 supplementary block F test",
        "multiplicity": "Each comparison separately within parcel432, Yeo7, system4, Juelich6; aliases not counted twice; not joint FDR across comparisons/families",
        "all_units": "Retain nonsignificant and failed fits; inspect status; invalid results are missing, not zeros",
        "full_model_erq_coefficient": "beta_erq_reference and SE/CI/t/p_erq/q_erq in erq_age_statistics.csv describe age-averaged ERQ slope in E3, distinct from the E2 coefficient",
        "fit_metrics": "RSS/R2/EDF/GCV/lambdas and signed delta metrics in erq_age_statistics.csv describe the three penalized pyGAM fits; inference RSS and partial R2 use OLS reference fits",
        "rounding": "Word is rounded for display; CSV has 17 significant digits; p/q shown as <.001 are not exact zero",
        "bayes_factors": "not calculated for ERQ", "simulated": mock,
    }
    (output_dir / "erq_table_definitions.json").write_text(json.dumps(definition, indent=2), encoding="utf-8")
    if not include_docx:
        return long
    document = _document("ERQ associations and age interactions",
                         f"All {len(statistics)} analysis units are retained. E1 includes smooth age and covariates; "
                         "E2 adds ERQ; E3 adds the smooth age by ERQ interaction. "
                         "Each comparison is BH corrected separately within each analysis family. "
                         "Rows show OLS reference inference in the same fixed spline space as the GAM.", mock)
    table = _three_line_table(document,
        ["Region", "N", "Effect", "Beta\n(SE)", "95% CI", "t or F\n(df)", "p", "q", "Partial\nR²"],
        [1.65, .35, .85, .85, .95, 1.05, .55, .55, .6])
    for record in long.to_dict("records"):
        main = record["comparison"] == "additive_vs_baseline"
        beta = f"{_number(record['beta'])}\n({_number(record['se'])})" if main else "—"
        ci = f"[{_number(record['ci_low'])},\n{_number(record['ci_high'])}]" if main else "—"
        df2 = _number(record["df2"]).split(".")[0]
        df1 = _number(record["df1"]).split(".")[0]
        stat = (f"t={_number(record['statistic'])}\n({df2})" if main
                else f"F={_number(record['statistic'])}\n({df1}, {df2})")
        label = f"{FAMILIES.get(record['dimension'], record['dimension'])} {record['label']}"
        if record["status"] != "ok":
            label += f"\n{record['status']}"
        _append_row(table, [label, record["n_obs"], record["effect"], beta, ci, stat,
                           _number(record["p"], True), _number(record["q"], True), _number(record["partial_r2_reference"])],
                    keep=record["comparison"] != "full_vs_baseline")
    document.add_paragraph("Notes  Beta is the E2 constant ERQ slope. Smooth interaction and joint tests have no single beta. "
                           "CI denotes a 95% confidence interval. A dash denotes not applicable or an invalid fit; consult status in CSV. "
                           "Pointwise E3 slopes and the age-averaged E3 ERQ coefficient are separate from the E2 coefficient. "
                           "The full-precision CSV and definitions file retain the numerical fields and valid BH family sizes.")
    document.save(output_dir / "erq_manuscript_table.docx")
    return long


def write_age_tables(statistics: pd.DataFrame, output_dir: Path, include_docx=True, mock=False):
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    long = age_long_table(statistics)
    long.to_csv(output_dir / "age_bayes_comparisons.csv", index=False, float_format="%.17g")
    if not include_docx:
        return long
    document = _document("Bayesian evidence for age associations",
                         f"All {len(statistics)} units are retained. Main linear and nonlinear prior scales are 0.5 residual SD. "
                         "The sensitivity interval spans nine combinations of scales 0.25, 0.5 and 1. "
                         "BF10 favors the first model; log BF is the natural logarithm. "
                         "Raw Bayes factors are not multiplicity adjusted.", mock)
    table = _three_line_table(document,
        ["Region", "N", "Models", "BF10", "log BF10", "Sensitivity\nlog BF10", "Evidence"],
        [1.7, .35, .7, .8, .75, 1.25, 1.15])
    for row in long.to_dict("records"):
        label = f"{FAMILIES.get(row['dimension'], row['dimension'])} {row['label']}"
        if row["status"] != "ok":
            label += f"\n{row['status']}"
        _append_row(table, [label, row["n_obs"], row["comparison"].replace("_vs_", " vs "),
                           _number(row["bf10"]), _number(row["log_bf10"]),
                           f"[{_number(row['sensitivity_log_bf10_min'])},\n{_number(row['sensitivity_log_bf10_max'])}]",
                           str(row["evidence"]).replace("strong_", "strong\n").replace("_", " ")
                           if row["status"] == "ok" else "—"],
                    keep=row["comparison"] != "m2_vs_m1")
    document.add_paragraph("Notes  Evidence labels at BF 3 and 10 are descriptive conventions. "
                           "M3 vs M2 addresses nonlinearity, not the absence of an age association. "
                           "A log BF sensitivity interval crossing zero indicates prior-dependent direction. "
                           "Source frequentist statistics, when available, are retained separately in the full CSV.")
    document.save(output_dir / "age_bayes_manuscript_table.docx")
    return long


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--kind", choices=("erq", "age"), required=True)
    parser.add_argument("--statistics", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--mock", action="store_true", help="Mark synthetic examples explicitly")
    args = parser.parse_args()
    statistics = pd.read_csv(args.statistics, dtype={"name": str})
    writer = write_erq_tables if args.kind == "erq" else write_age_tables
    writer(statistics, args.output, mock=args.mock)


if __name__ == "__main__":
    main()
