"""Full precision summaries, SVG figures and Word three-line tables."""
from __future__ import annotations

import json
from pathlib import Path
import numpy as np
import pandas as pd

from .io import MODELS, SCHEMES, own_output, write_csv

LABELS = {"M0_covariates": "Covariates", "M1_PFC": "PFC", "M2_distributed432": "Distributed", "M3_nonPFC": "Non-PFC"}
TITLES = {"joint": "Joint age and success prediction", "shared": "Common reference behavior score prediction"}


def number(value):
    return "NA" if pd.isna(value) else f"{float(value):.3f}"


def markdown_table(headers, rows):
    return "\n".join(["| " + " | ".join(headers) + " |", "|" + "---|" * len(headers),
                      *("| " + " | ".join(map(str, r)) + " |" for r in rows)])


def report(root: Path, docx=True, plots=True):
    root = root.resolve()
    output = root / "report"
    own_output(output, "joint_pls_report_v1")
    manifest = json.loads((root / "analysis_manifest.json").read_text(encoding="utf-8"))
    sample = json.loads((root / "precheck" / "sample_metadata.json").read_text(encoding="utf-8"))
    performance = pd.read_csv(root / "model_performance.csv")
    comparisons = pd.read_csv(root / "model_comparisons.csv")
    axes = pd.read_csv(root / "reference_axes.csv")
    primary = comparisons.loc[comparisons.comparison.eq("distributed_minus_pfc")].copy()
    write_csv(performance, output / "model_performance_full_precision.csv")
    write_csv(comparisons, output / "model_comparisons_full_precision.csv")
    write_csv(primary, output / "primary_comparisons.csv")
    texts = []
    if manifest.get("simulated"):
        texts.append("**SIMULATED DATA ONLY — do not use these values in a manuscript.**")
    texts.append(f"Common N = {sample['n']}; sites = {sample['sites']}. Both outcomes are modeled in one PLS2 fit. "
                 "Gender and meanFD adjustment and X/Y scaling are fitted only on training rows. "
                 "Parent education is neither a covariate nor an eligibility criterion. "
                 "Site labels group inner/outer LOSO and restrict permutation; no test-site regression is fitted.")
    sections = []
    for scheme in SCHEMES:
        p = performance.loc[performance.scheme.eq(scheme)].set_index("model")
        rows = [[LABELS[m], int(p.at[m, "n"]), int(p.at[m, "n_features"]),
                 number(p.at[m, "r2"]), number(p.at[m, "mae"]), number(p.at[m, "rmse"])] for m in MODELS]
        headers = ["Model", "N", "Regions", "R2", "MAE", "RMSE"]
        comp = comparisons.loc[comparisons.scheme.eq(scheme)]
        comp_rows = [[f"{LABELS[r.full_model]} vs {LABELS[r.reduced_model]}", number(r.delta_r2),
                      number(r.p_permutation), number(r.q_primary), str(r.inference)] for r in comp.itertuples()]
        comp_headers = ["Comparison", "Delta R2", "P", "Primary q", "Inference status"]
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
            sentence = (f"The distributed model improved prediction relative to PFC (delta R2 = {number(row.delta_r2)}, "
                        f"permutation P = {number(row.p_permutation)}, q = {number(row.q_primary)}), under the specified reduced-model null.")
        else:
            sentence = (f"The distributed-minus-PFC difference was {number(row.delta_r2)} (permutation P = "
                        f"{number(row.p_permutation)}, q = {number(row.q_primary)}); the conditional permutation test "
                        "did not support a positive incremental contribution after the prespecified correction.")
        texts.append(sentence)
        sections.append((scheme, headers, rows, comp_headers, comp_rows, note, sentence))
    warning = ("Permutation inference uses an approximate Freedman-Lane reduced linear model with site-restricted paired "
               "bivariate residual permutations. The PFC null includes all PFC features, gender, meanFD and categorical site. "
               "It preserves modeled PFC associations and tests additional conditional information; it is not a generic test "
               "of equality of prediction algorithms. Two primary route comparisons are jointly BH-corrected. "
               "Other comparisons are descriptive unless their corresponding null is separately run. "
               "Bootstrap is an independent optional branch and is never nested in permutation. "
               "Negative R2 is retained. A shared score may be dominated by a single behavior; inspect reference_axes.csv.")
    texts.append(warning)
    (output / "result_summary.md").write_text("\n\n".join(texts), encoding="utf-8")
    if docx:
        from analysis.fig03.gam_pygam.result_tables import _document, _three_line_table, _append_row
        document = _document("Joint and common score PLS prediction results",
                             f"Common sample N = {sample['n']}; sites = {sample['sites']}. "
                             "Nested leave-one-site-out validation. No parent education covariate or education completeness filter.",
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
                                   "Primary q corrects the two prespecified distributed versus PFC comparisons. "
                                   "Other P values, when requested, are exploratory and unadjusted.")
        document.save(output / "response_tables.docx")
    if plots:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        for scheme in SCHEMES:
            frame = performance.loc[performance.scheme.eq(scheme)].set_index("model").loc[list(MODELS)]
            fig, ax = plt.subplots(figsize=(6.6, 4.5))
            ax.bar([LABELS[m] for m in MODELS], frame.r2, color=["#999999", "#4e9770", "#8b79ae", "#c49b69"])
            ax.axhline(0, color="black", lw=.7)
            ax.set(ylabel="Cross-validated R²", title=TITLES[scheme])
            ax.spines[["top", "right"]].set_visible(False)
            if manifest.get("simulated"):
                fig.text(.5, .01, "SIMULATED DATA ONLY", ha="center", fontsize=9)
            fig.tight_layout(rect=(0, .035, 1, 1))
            fig.savefig(output / f"{scheme}_model_performance.svg")
            plt.close(fig)
    return output


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Standalone CSV/Word export; no model fitting or sklearn required")
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--skip-docx", action="store_true")
    parser.add_argument("--no-plots", action="store_true")
    args = parser.parse_args()
    print(report(args.output, not args.skip_docx, not args.no_plots))
