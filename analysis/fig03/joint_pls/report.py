"""Full precision summaries, SVG figures and Word three-line tables."""
from __future__ import annotations

import json
from pathlib import Path
import numpy as np
import pandas as pd

from .io import MODELS, SCHEMES, own_output, write_csv
from .comparison_summary import summarize
from .visualize import probability

LABELS = {"M0_covariates": "Covariates", "M1_PFC": "PFC", "M2_distributed432": "Distributed", "M3_nonPFC": "Non-PFC"}
TITLES = {"joint": "Joint age and success prediction", "shared": "Common reference behavior score prediction"}


def number(value):
    return "NA" if pd.isna(value) else f"{0. if abs(float(value)) < .0005 else float(value):.3f}"


def markdown_table(headers, rows):
    return "\n".join(["| " + " | ".join(headers) + " |", "|" + "---|" * len(headers),
                      *("| " + " | ".join(map(str, r)) + " |" for r in rows)])


def report(root: Path, docx=True, plots=True, preview=False, model_ci=False):
    root = root.resolve()
    output = root / "report"
    own_output(output, "joint_pls_report_v1")
    manifest = json.loads((root / "analysis_manifest.json").read_text(encoding="utf-8"))
    sample = json.loads((root / "precheck" / "sample_metadata.json").read_text(encoding="utf-8"))
    definition = sample.get("pfc_definition", {})
    labels = LABELS.copy()
    meta_pfc = definition.get("source") == "Buhle2014_reappraisal_PFC"
    if meta_pfc:
        labels.update(M1_PFC="Meta-PFC", M3_nonPFC="Remaining regions")
    definition_note = (f"PFC definition: {definition['source']}; {definition['n_pfc']} PFC regions, "
                       f"{definition['n_remaining']} remaining regions." if definition else
                       "PFC definition was not recorded in this older precheck; inspect the original membership source.")
    performance = pd.read_csv(root / "model_performance.csv")
    comparisons, permutations = summarize(root, output)
    primary = comparisons.loc[comparisons.comparison.eq("distributed_minus_pfc")].copy()
    write_csv(performance, output / "model_performance_full_precision.csv")
    texts = []
    if manifest.get("simulated"):
        texts.append("**SIMULATED DATA ONLY — do not use these values in a manuscript.**")
    texts.append(definition_note)
    if definition:
        texts.append(f"Membership rule: `{definition['rule']}`.")
    if meta_pfc:
        texts.append("Remaining regions are the complement of the selected meta-PFC set; "
                     "they may include anatomical prefrontal regions outside the meta-analysis mask.")
    texts.append(f"Common N = {sample['n']}; sites = {sample['sites']}. Both outcomes are modeled in one PLS2 fit. "
                 "Gender and meanFD adjustment and X/Y scaling are fitted only on training rows. "
                 "Parent education is neither a covariate nor an eligibility criterion. "
                 "Site labels group inner/outer LOSO and restrict permutation; no test-site regression is fitted.")
    texts.append("All five signed pairwise comparisons are displayed in this order: M1-M0, M2-M0, M3-M0, "
                 "M2-M1 and M3-M1. M3-M2 / M2-M3 is omitted because it is not an analysis target. "
                 "Joint and shared routes have separate figures. "
                 "The PFC comparisons are M2 minus M1 (all regions versus PFC) and M3 minus M1 "
                 "(remaining regions versus PFC). The first comparison evaluates "
                 "adding remaining regions to PFC; the second contrasts two separate feature sets and does "
                 "not itself establish incremental information after combining them. "
                 "The new contrast is recovered within each original paired repeat as (M3-M0)-(M1-M0); "
                 "neither model fits nor P values from different contrasts are reused as if equivalent.")
    sections = []
    for scheme in SCHEMES:
        p = performance.loc[performance.scheme.eq(scheme)].set_index("model")
        rows = [[labels[m], int(p.at[m, "n"]), int(p.at[m, "n_features"]),
                 number(p.at[m, "r2"]), number(p.at[m, "mae"]), number(p.at[m, "rmse"])] for m in MODELS]
        headers = ["Model", "N", "Regions", "R2", "MAE", "RMSE"]
        comp = comparisons.loc[comparisons.scheme.eq(scheme)]
        comp_rows = [[f"{labels[r.full_model]} vs {labels[r.reduced_model]}", number(r.delta_r2),
                      probability(r.p_permutation), probability(r.q_bh),
                      str(r.inference)] for r in comp.itertuples()]
        comp_headers = ["Comparison", "Delta R2", "P", "BH q", "Inference status"]
        note = ("R2 is the equally weighted mean of the two R2 values for cross-fitted, covariate-adjusted, "
                "training-standardized outcomes. MAE/RMSE pool both standardized outcomes. This is joint prediction, not pure shared variance."
                if scheme == "joint" else
                "The target is the first Y direction of a full432 reference PLS fitted within each training fold. "
                "Every model uses the same held-out target in that fold. Target construction is repeated inside inner CV. "
                "The score is training-unit-variance and age-oriented; it is brain-informed and not an independent behavioral factor.")
        texts += [f"\n## {TITLES[scheme]}", markdown_table(headers, rows), note,
                  markdown_table(comp_headers, comp_rows)]
        row = primary.loc[primary.scheme.eq(scheme)].iloc[0]
        if pd.isna(row.q_primary):
            sentence = (f"The distributed-minus-PFC difference in R2 was {number(row.delta_r2)}. "
                        "An inferential conclusion is pending completion of the prespecified conditional permutation test.")
        elif row.delta_r2 > 0 and row.q_primary < .05:
            sentence = (f"The distributed model improved prediction relative to {labels['M1_PFC']} (delta R2 = {number(row.delta_r2)}, "
                        f"permutation P = {probability(row.p_permutation)}, q = {probability(row.q_primary)}), under the specified reduced-model null.")
        else:
            sentence = (f"The distributed-minus-PFC difference was {number(row.delta_r2)} (permutation P = "
                        f"{probability(row.p_permutation)}, q = {probability(row.q_primary)}); the conditional permutation test "
                        "did not support a positive incremental contribution after the prespecified correction.")
        secondary = comp.loc[comp.comparison.eq("nonpfc_minus_pfc")].iloc[0]
        sentence += (f" The remaining-regions-minus-PFC difference was {number(secondary.delta_r2)} "
                     f"(conditional permutation P = {probability(secondary.p_permutation)}, "
                     f"secondary BH q = {probability(secondary.q_secondary)}).")
        texts.append(sentence)
        sections.append((scheme, headers, rows, comp_headers, comp_rows, note, sentence))
    warning = ("Permutation inference uses an approximate Freedman-Lane reduced linear model with site-restricted paired "
               "bivariate residual permutations. The PFC null includes all PFC features, gender, meanFD and categorical site. "
               "It preserves modeled PFC associations and tests additional conditional information; it is not a generic test "
               "of equality of prediction algorithms. Two primary M2-minus-M1 route comparisons are jointly BH-corrected. "
               "The two newly requested M3-minus-M1 route comparisons form a separate secondary BH family. "
               "The three model-minus-M0 contrasts across both routes form an additional exploratory BH family of six tests; "
               "these require the covariates-null permutations. This is not a combined ten-test correction. "
               "The q_bh column reports the appropriate within-family adjusted P. "
               "M3-M2 / M2-M3 is omitted from the reported comparisons. "
               "Bootstrap is an independent optional branch and is never nested in permutation. "
               "Negative R2 is retained. A shared score may be dominated by a single behavior; inspect reference_axes.csv.")
    texts.append(warning)
    (output / "result_summary.md").write_text("\n\n".join(texts), encoding="utf-8")
    if docx:
        from analysis.fig03.gam_pygam.result_tables import _document, _three_line_table, _append_row
        document = _document("Joint and common score PLS prediction results",
                             f"Common sample N = {sample['n']}; sites = {sample['sites']}. "
                             "Nested leave-one-site-out validation. No parent education covariate or education completeness filter. "
                             + definition_note,
                             manifest.get("simulated", False))
        for index, (scheme, headers, rows, ch, cr, note, sentence) in enumerate(sections):
            if index:
                document.add_page_break()
            document.add_paragraph(TITLES[scheme]).runs[0].bold = True
            table = _three_line_table(document, headers, [1.45, .55, .7, .85, .85, .85])
            for row in rows:
                _append_row(table, row)
            document.add_paragraph(note)
            document.add_paragraph("Direct model comparisons").runs[0].bold = True
            # Long implementation status remains in CSV; use a compact Word table.
            table = _three_line_table(document, ch[:4], [3.35, .85, .7, .7])
            for row in cr:
                _append_row(table, row[:4])
            document.add_paragraph(sentence)
            document.add_paragraph("P is available only for a fully completed appropriate reduced-model null. "
                                   "BH q corrects M2 vs M1 across the two routes as the original primary family, "
                                   "and M3 vs M1 across routes as a separate secondary family. "
                                   "The three baseline comparisons across both routes form a separate six-test exploratory BH family. "
                                   "These are three separate correction families, not one combined ten-test correction. "
                                   "The reference score may primarily represent success; inspect its behavior loadings.")
        document.save(output / "response_tables.docx")
    if plots:
        from .visualize import figures
        model_intervals = None
        if model_ci:
            from .bootstrap_endpoints import read_intervals
            model_intervals = read_intervals(root)
            if model_intervals is None or not comparisons.bootstrap_complete.all():
                raise ValueError("--model-ci requires complete recovered model intervals; run bootstrap_endpoints first")
        figures(performance, comparisons, permutations, output, labels,
                simulated=manifest.get("simulated", False), preview=preview,
                model_intervals=model_intervals)
    return output


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Standalone CSV/Word export; no model fitting or sklearn required")
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--skip-docx", action="store_true")
    parser.add_argument("--no-plots", action="store_true")
    parser.add_argument("--preview", action="store_true", help="Also render PNG previews under report/preview; formal figures remain SVG")
    parser.add_argument("--model-ci", action="store_true", help="Show optional single-model R2 confidence intervals; requires bootstrap_endpoints output")
    args = parser.parse_args()
    print(report(args.output, not args.skip_docx, not args.no_plots, args.preview, args.model_ci))
