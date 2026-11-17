"""Publication SVGs for held-out performance and paired PFC comparisons."""
from __future__ import annotations

from pathlib import Path
import numpy as np
import pandas as pd

from .io import SCHEMES, MODELS, write_csv

COLORS = {"M0_covariates": "#0072FF", "M1_PFC": "#00A651",
          "M2_distributed432": "#8B2CFF", "M3_nonPFC": "#FF7A00"}
ROUTES = {"joint": "Joint age + success prediction", "shared": "Reference score prediction"}
FOCUS = ("distributed_minus_pfc", "nonpfc_minus_pfc")


def probability(value):
    return "pending" if pd.isna(value) else f"{float(value):.3g}"


def _style():
    import matplotlib
    matplotlib.use("Agg")
    return {"font.family": "Arial", "font.size": 11, "axes.titlesize": 12,
            "axes.labelsize": 11, "svg.fonttype": "none", "axes.spines.top": False,
            "axes.spines.right": False, "axes.linewidth": .7,
            "xtick.labelsize": 10, "ytick.labelsize": 10, "savefig.facecolor": "white"}


def _save(fig, path, preview):
    fig.savefig(path, bbox_inches="tight")
    if preview:
        folder = path.parent / "preview"
        folder.mkdir(exist_ok=True)
        fig.savefig(folder / f"{path.stem}.png", dpi=180, bbox_inches="tight")


def dumbbell_figures(performance, comparisons, output, simulated=False, preview=False, source_note=None,
                     pfc_label="PFC"):
    """Connect observed model R²; show paired-difference CIs on their own axis.

    Repeat archives contain differences, not absolute model R². Those intervals
    must not be attached to either model endpoint or reconstructed as model CIs.
    """
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    from matplotlib.ticker import MaxNLocator, FormatStrFormatter

    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    pfc_count = int(performance.loc[performance.model.eq("M1_PFC"), "n_features"].iloc[0])
    remaining_count = int(performance.loc[performance.model.eq("M3_nonPFC"), "n_features"].iloc[0])
    n, sites = int(performance["n"].iloc[0]), int(performance.sites.iloc[0])
    model_names = {"M2_distributed432": "All regions", "M3_nonPFC": "Remaining regions"}
    route_names = {"joint": "Joint age + success", "shared": "Reference score"}
    selections = [("pfc_remaining_dumbbell", ("nonpfc_minus_pfc",)),
                  ("pfc_all_comparisons_dumbbell", FOCUS)]
    for filename, names in selections:
        order = [(s, name) for s in SCHEMES for name in names]
        index = comparisons.set_index(["scheme", "comparison"])
        selected = index.loc[order].copy()
        selected = selected.reset_index()
        models = performance.set_index(["scheme", "model"])
        table = selected.copy()
        table["pfc_r2"] = [models.at[(s, "M1_PFC"), "r2"] for s, _ in order]
        table["other_r2"] = [models.at[(r.scheme, r.full_model), "r2"] for r in selected.itertuples()]
        table["difference_interval"] = "paired within-site bootstrap 95% percentile interval for delta_r2 only"
        write_csv(table, output / f"{filename}_data.csv")
        with plt.rc_context(_style()):
            fig, axes = plt.subplots(1, 3, figsize=(10.8, 3.9 if len(names) == 1 else 5.2),
                                     gridspec_kw={"width_ratios": [1.2, 1, .57], "wspace": .42})
            fig.subplots_adjust(left=.205, right=.97, bottom=.29 if len(names) == 1 else .24, top=.72 if len(names) == 1 else .77)
            model_ax, delta_ax, stats_ax = axes
            positions = np.arange(len(order))[::-1]
            limits = (table[["pfc_r2", "other_r2"]].to_numpy(float),
                      table[["delta_r2", "ci_low", "ci_high"]].to_numpy(float))
            for ax, values in zip(axes[:2], limits):
                values = values[np.isfinite(values)]
                low, high = min(0., values.min()), max(0., values.max())
                span = max(high-low, .04)
                ax.set_xlim(low-.14*span, high+.17*span)
                ax.axvline(0, lw=.8, ls=(0, (3, 3)), color="#ACB1B5", zorder=0)
                ax.spines["left"].set_visible(False)
                ax.tick_params(axis="y", length=0)
                ax.xaxis.set_major_locator(MaxNLocator(4))
                ax.xaxis.set_major_formatter(FormatStrFormatter("%.2f"))
            for ax in axes:
                ax.set_ylim(-.5, len(order)-.35)
            labels = []
            for pos, row in zip(positions, table.itertuples()):
                color = COLORS[row.full_model]
                model_ax.plot([row.pfc_r2, row.other_r2], [pos, pos], color="#B8BEC3", lw=2.5, zorder=1)
                model_ax.scatter(row.pfc_r2, pos, s=82, color=COLORS["M1_PFC"], edgecolors="white", linewidths=.8, zorder=3)
                model_ax.scatter(row.other_r2, pos, s=82, color=color, edgecolors="white", linewidths=.8, zorder=3)
                for value in (row.pfc_r2, row.other_r2):
                    rounded = 0. if abs(value) < .0005 else value
                    model_ax.annotate(f"{rounded:.3f}", (value, pos), xytext=(0, 12), textcoords="offset points",
                                      ha="center", va="bottom", fontsize=10, color="#41474C")
                interval_ready = pd.notna(row.ci_low) and pd.notna(row.ci_high)
                if interval_ready:
                    delta_ax.plot([row.ci_low, row.ci_high], [pos, pos], lw=2.1, color=color, zorder=1)
                    delta_ax.scatter([row.ci_low, row.ci_high], [pos, pos], marker="|", color=color, s=56, zorder=2)
                delta_ax.scatter(row.delta_r2, pos, s=65, color=color, edgecolors="white", linewidths=.6, zorder=3)
                interval_text = f"[{row.ci_low:.3f}, {row.ci_high:.3f}]" if interval_ready else "95% CI pending"
                delta_ax.text(.5, pos-.24, f"{row.delta_r2:.3f}  {interval_text}", transform=delta_ax.get_yaxis_transform(),
                              ha="center", va="center", fontsize=9, color="#51575C")
                q = row.q_primary if row.comparison == FOCUS[0] else row.q_secondary
                stats_ax.text(.03, pos+.09, f"P = {probability(row.p_permutation)}", fontsize=10, va="center")
                stats_ax.text(.03, pos-.13, f"q = {probability(q)}", fontsize=10, va="center", color="#51575C")
                label = route_names[row.scheme]
                if len(names) > 1:
                    label += "\n" + model_names[row.full_model]
                labels.append(label)
            model_ax.set(yticks=positions, yticklabels=labels, xlabel="Cross-validated R²")
            model_ax.tick_params(axis="y", pad=13)
            delta_ax.set(yticks=[], xlabel="ΔR² relative to PFC")
            stats_ax.set(xticks=[], yticks=[])
            stats_ax.set_axis_off()
            model_ax.set_title("Model performance", loc="left", pad=23, fontweight="bold", fontsize=11)
            delta_ax.set_title("Paired difference · 95% CI", loc="left", pad=23, fontweight="bold", fontsize=11)
            stats_ax.set_title("Permutation test", loc="left", pad=23, fontweight="bold", fontsize=11)
            if len(names) > 1:
                for ax in axes:
                    ax.axhline(1.5, color="#E2E5E7", lw=.7, zorder=0)
            legend_models = ["M1_PFC", "M3_nonPFC"] if len(names) == 1 else ["M1_PFC", "M2_distributed432", "M3_nonPFC"]
            legend_names = {"M1_PFC": f"{pfc_label} ({pfc_count})", "M2_distributed432": "All regions (432)",
                            "M3_nonPFC": f"Remaining regions ({remaining_count})"}
            handles = [Line2D([], [], marker="o", ls="none", color=COLORS[m], markersize=7, label=legend_names[m])
                       for m in legend_models]
            heading = "Remaining regions versus PFC" if len(names) == 1 else "Brain model comparisons relative to PFC"
            fig.suptitle(heading, x=.205, ha="left", y=.985, fontsize=14, fontweight="bold")
            fig.legend(handles=handles, loc="upper left", bbox_to_anchor=(.198, .92),
                       ncol=len(handles), frameon=False, fontsize=10, handletextpad=.4, columnspacing=1.6)
            fig.text(.97, .956, f"N = {n:,} · {sites} sites", ha="right", fontsize=10, color="#51575C")
            foot_y = .095 if len(names) == 1 else .08
            fig.text(.205, foot_y+.065, "Intervals describe ΔR², not individual model R². P and q are from conditional permutation tests.",
                     fontsize=8.4, color="#62686D")
            source = ("SIMULATED DATA — layout demonstration only." if simulated else
                      "Primary and secondary comparisons use separate BH correction families.")
            if source_note:
                source = "Screenshot point estimates; bootstrap intervals and new M3−M1 P/q require the server archives."
            fig.text(.205, foot_y+.02, source, fontsize=8.4, color="#62686D")
            _save(fig, output / f"{filename}.svg", preview)
            if filename == "pfc_remaining_dumbbell":
                _save(fig, output / "model_comparison_overview.svg", preview)
            plt.close(fig)


def legacy_figures(performance, comparisons, permutations, output, labels, simulated=False,
            preview=False, source_note=None):
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    from matplotlib.ticker import MaxNLocator, FormatStrFormatter
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    focus = comparisons.loc[comparisons.comparison.isin(FOCUS)]
    plot_data = focus.merge(performance[["scheme", "model", "r2"]].rename(columns={"model": "full_model", "r2": "full_r2"}),
                            on=["scheme", "full_model"], validate="many_to_one")
    plot_data = plot_data.merge(performance[["scheme", "model", "r2"]].rename(columns={"model": "reduced_model", "r2": "pfc_r2"}),
                                on=["scheme", "reduced_model"], validate="many_to_one")
    write_csv(plot_data, output / "figure_comparison_data.csv")
    dumbbell_figures(performance, comparisons, output, simulated=simulated, preview=preview,
                     source_note=source_note, pfc_label=labels["M1_PFC"])
    footer = ("SIMULATED DATA ONLY" if simulated else
              source_note or "Held-out sites; paired conditional permutation tests. Primary and secondary BH families are separate.")
    with plt.rc_context(_style()):
        fig, axes = plt.subplots(2, 2, figsize=(12.8, 7.6), gridspec_kw={"width_ratios": [1, 1.55]})
        fig.subplots_adjust(left=.15, right=.97, top=.90, bottom=.17, wspace=.70, hspace=.69)
        observed_values = performance.r2.to_numpy()
        low, high = min(observed_values.min(), 0), max(observed_values.max(), 0)
        span = max(high-low, .04)
        for row, scheme in enumerate(SCHEMES):
            ax, effect = axes[row]
            p = performance.loc[performance.scheme.eq(scheme)].set_index("model").loc[list(MODELS)]
            y = np.arange(4)[::-1]
            ax.axvline(0, color="#A4A4A4", lw=.8, ls="--", zorder=0)
            for pos, model in zip(y, MODELS):
                val = p.at[model, "r2"]
                ax.scatter(val, pos, s=62, color=COLORS[model], zorder=3)
                rounded = 0. if abs(val) < .0005 else val
                ax.annotate(f"{rounded:.3f}", (val, pos), xytext=(7, 0), textcoords="offset points", va="center", fontsize=10)
            ax.set(yticks=y, yticklabels=[f"{labels[m]} ({int(p.at[m, 'n_features'])})" for m in MODELS],
                   ylim=(-.6, 3.6), xlim=(low-.1*span, high+.38*span), xlabel="Cross-validated R²")
            ax.xaxis.set_major_locator(MaxNLocator(4))
            ax.xaxis.set_major_formatter(FormatStrFormatter("%.2f"))
            ax.spines["left"].set_visible(False)
            ax.tick_params(axis="y", length=0)
            ax.set_title(f"{'A' if row == 0 else 'C'}   {ROUTES[scheme]}", loc="left", pad=16, fontweight="bold")
            selected = focus.loc[focus.scheme.eq(scheme)].set_index("comparison").loc[list(FOCUS)]
            extent = selected[["delta_r2", "ci_low", "ci_high"]].to_numpy(float)
            finite = extent[np.isfinite(extent)]
            elow, ehigh = min(finite.min(), 0), max(finite.max(), 0)
            espan = max(ehigh-elow, .03)
            effect.axvline(0, color="#A4A4A4", lw=.8, ls="--", zorder=0)
            for pos, name in zip([1, 0], FOCUS):
                r = selected.loc[name]
                color = COLORS[r.full_model]
                if pd.notna(r.ci_low) and pd.notna(r.ci_high):
                    effect.plot([r.ci_low, r.ci_high], [pos, pos], lw=2, color=color, zorder=2)
                    effect.scatter([r.ci_low, r.ci_high], [pos, pos], marker="|", s=50, color=color)
                effect.scatter(r.delta_r2, pos, s=70, color=color, zorder=3)
                qcol = "q_primary" if name == FOCUS[0] else "q_secondary"
                stats = (f"ΔR² = {r.delta_r2:.3f}\nP = {probability(r.p_permutation)}; "
                         f"q = {probability(r[qcol])}")
                effect.text(.985, pos+.13, stats, transform=effect.get_yaxis_transform(),
                            ha="right", va="center", fontsize=10, linespacing=1.6,
                            bbox=dict(facecolor="white", edgecolor="none", alpha=.92, pad=1.5))
            effect.set(yticks=[1, 0], yticklabels=[f"{labels['M2_distributed432']} − PFC", f"{labels['M3_nonPFC']} − PFC"],
                       ylim=(-.5, 1.5), xlim=(elow-.08*espan, ehigh+.88*espan),
                       xlabel="ΔR² relative to PFC")
            effect.xaxis.set_major_locator(MaxNLocator(4))
            effect.spines["left"].set_visible(False)
            effect.tick_params(axis="y", length=0)
            effect.set_title(f"{'B' if row == 0 else 'D'}   Paired model comparisons", loc="left", pad=16, fontweight="bold")
        n, sites = int(performance["n"].iloc[0]), int(performance.sites.iloc[0])
        fig.suptitle(f"PLS model comparisons   |   N = {n:,}; {sites} held-out sites", x=.15, ha="left", fontsize=15, fontweight="bold")
        has_intervals = focus.ci_low.notna().all()
        note = ("Lines: 95% within-site bootstrap intervals; P/q come from permutation, not interval overlap." if has_intervals else
                "Bootstrap intervals are not shown until every requested repeat is valid.")
        fig.text(.15, .085, note, fontsize=9, color="#555555")
        fig.text(.15, .052, "Joint R² averages two outcomes. The reference score can be dominated by one behavior; inspect its loadings.", fontsize=9, color="#555555")
        fig.text(.15, .019, footer, fontsize=9, color="#555555")
        _save(fig, output / "model_comparison_dotplot.svg", preview)
        plt.close(fig)

        # Individual original per-route filenames remain available as dot plots.
        for scheme in SCHEMES:
            fig, ax = plt.subplots(figsize=(6.7, 3.5))
            p = performance.loc[performance.scheme.eq(scheme)].set_index("model").loc[list(MODELS)]
            ax.axvline(0, color="#AAAAAA", ls="--", lw=.8)
            ax.scatter(p.r2, np.arange(4)[::-1], c=[COLORS[m] for m in MODELS], s=65)
            ax.set(yticks=np.arange(4)[::-1], yticklabels=[labels[m] for m in MODELS],
                   xlabel="Cross-validated R²", title=ROUTES[scheme], ylim=(-.5, 3.5))
            fig.tight_layout()
            _save(fig, output / f"{scheme}_model_performance.svg", False)
            plt.close(fig)

        fig, axes = plt.subplots(2, 2, figsize=(10.8, 7.4))
        fig.subplots_adjust(hspace=.7, wspace=.35, bottom=.16, top=.88)
        for row, scheme in enumerate(SCHEMES):
            for col, name in enumerate(FOCUS):
                ax = axes[row, col]
                r = focus.loc[focus.scheme.eq(scheme) & focus.comparison.eq(name)].iloc[0]
                values = permutations.loc[permutations.scheme.eq(scheme) & permutations.comparison.eq(name) &
                                          permutations["null"].eq("pfc"), "delta_r2"].to_numpy(float)
                if len(values):
                    bins = min(35, max(5, int(np.sqrt(len(values)))))
                    ax.hist(values, bins=bins, color="#D5DAE0", edgecolor="white", linewidth=.6)
                    ax.axvline(r.delta_r2, color=COLORS[r.full_model], lw=2, label=f"Observed ΔR² = {r.delta_r2:.3f}")
                    ax.axvline(0, color="#888888", ls=":", lw=.8)
                    ax.legend(frameon=False, fontsize=9, loc="upper right")
                else:
                    ax.text(.5, .5, "No raw permutation records available", ha="center", va="center", transform=ax.transAxes)
                ax.set(xlabel="ΔR² under the PFC-reduced null", ylabel="Permutation count")
                ax.set_title(f"{ROUTES[scheme]}\n{labels[r.full_model]} − PFC", loc="left", fontsize=11)
                ax.text(.02, .97, f"Valid: {len(values)}/{int(r.permutations_requested)}\nP = {probability(r.p_permutation)}",
                        transform=ax.transAxes, va="top", fontsize=9)
        fig.suptitle("Conditional permutation tests of PFC comparisons", fontsize=15, fontweight="bold")
        fig.text(.08, .07, "The same paired repeats are used for both contrasts. A nonzero null mean is retained.", fontsize=9, color="#555555")
        fig.text(.08, .035, footer, fontsize=9, color="#555555")
        _save(fig, output / "pfc_comparison_permutation_nulls.svg", preview)
        plt.close(fig)


def figures(performance, comparisons, permutations, output, labels, simulated=False,
            preview=False, source_note=None, model_intervals=None):
    """Two-panel five-comparison figures in a compact side-by-side layout."""
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    from matplotlib.ticker import MaxNLocator, FormatStrFormatter
    from .comparison_summary import COMPARISONS

    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    index = performance.set_index(["scheme", "model"])
    interval_index = None
    if model_intervals is not None:
        if (len(model_intervals) != 8 or model_intervals.duplicated(["scheme", "model"]).any() or
                not model_intervals.complete.eq(True).all()):
            raise ValueError("Expected eight unique completed absolute-model bootstrap intervals")
        interval_index = model_intervals.set_index(["scheme", "model"])
    names = [row[0] for row in COMPARISONS]
    all_plot_rows = []
    for scheme in SCHEMES:
        table = comparisons.loc[comparisons.scheme.eq(scheme)].set_index("comparison").loc[names].reset_index()
        for role, key in (("compared", "full_model"), ("reference", "reduced_model")):
            table[f"{role}_r2"] = [index.at[(scheme, model), "r2"] for model in table[key]]
            for edge in ("low", "high"):
                table[f"{role}_ci_{edge}"] = ([interval_index.at[(scheme, model), f"ci_{edge}"] for model in table[key]]
                                               if interval_index is not None else np.nan)
        table["endpoint_interval"] = "absolute model R2, paired bootstrap percentile 95% CI" if interval_index is not None else "not_shown_optional"
        table["reference_color"] = table.reduced_model.map(COLORS)
        table["compared_color"] = table.full_model.map(COLORS)
        all_plot_rows.append(table)
        write_csv(table, output / f"{scheme}_all_model_comparisons_data.csv")
        with plt.rc_context(_style()):
            n_rows = len(table)
            fig, axes = plt.subplots(1, 2, figsize=(5.8, 4.15),
                                     gridspec_kw={"width_ratios": [1.08, 1], "wspace": .08})
            fig.subplots_adjust(left=.18, right=.99, top=.80, bottom=.18)
            model_ax, delta_ax = axes
            positions = np.arange(n_rows)[::-1]
            r2_values = table[["reference_r2", "compared_r2", "reference_ci_low", "reference_ci_high",
                               "compared_ci_low", "compared_ci_high"]].to_numpy(float)
            delta_values = table[["delta_r2", "ci_low", "ci_high"]].to_numpy(float)
            xlabels = ("Model R² (95% CI)" if interval_index is not None else "Model R²",
                       "ΔR² (95% CI)")
            for ax, values, xlab in zip(axes, (r2_values, delta_values), xlabels):
                values = values[np.isfinite(values)]
                low, high = min(0., values.min()), max(0., values.max())
                span = max(high-low, .025)
                ax.set_xlim(low-.10*span, high+.12*span)
                ax.spines["left"].set_visible(False)
                ax.tick_params(axis="y", length=0, pad=3, labelsize=14)
                ax.tick_params(axis="x", labelsize=14, pad=3, width=.7, length=3)
                ax.xaxis.set_major_locator(MaxNLocator(3))
                ax.xaxis.set_major_formatter(FormatStrFormatter("%.2f"))
                ax.set_xlabel(xlab, fontsize=15, labelpad=5)
                ax.set_ylim(-.42, n_rows - 1 + .42)
            row_labels = []
            for pos, row in zip(positions, table.itertuples()):
                full = row.full_model.split("_")[0]
                reduced = row.reduced_model.split("_")[0]
                row_labels.append(f"{full} − {reduced}")
                a, b = row.reference_r2, row.compared_r2
                model_ax.plot([a, b], [pos, pos], color="#AAB1B7", lw=1.5, zorder=1)
                for value, lo, hi, color in (
                        (a, row.reference_ci_low, row.reference_ci_high, COLORS[row.reduced_model]),
                        (b, row.compared_ci_low, row.compared_ci_high, COLORS[row.full_model])):
                    if pd.notna(lo) and pd.notna(hi):
                        model_ax.plot([lo, hi], [pos, pos], color=color, lw=1.4, zorder=2)
                        model_ax.scatter([lo, hi], [pos, pos], color=color, marker="|", s=28, zorder=2)
                    model_ax.scatter(value, pos, s=28, color=color, edgecolors="white", linewidths=.5, zorder=3)
                if pd.notna(row.ci_low) and pd.notna(row.ci_high):
                    delta_ax.plot([row.ci_low, row.ci_high], [pos, pos], color="#515C66", lw=1.4, zorder=2)
                    delta_ax.scatter([row.ci_low, row.ci_high], [pos, pos], color="#515C66", marker="|", s=26, zorder=2)
                delta_ax.scatter(row.delta_r2, pos, s=22, color="#515C66", zorder=3)
            model_ax.set(yticks=positions, yticklabels=row_labels)
            delta_ax.set(yticks=[])
            legend_labels = {"M0_covariates": "M0: Covariates", "M1_PFC": f"M1: {labels['M1_PFC']}",
                             "M2_distributed432": "M2: Whole brain", "M3_nonPFC": "M3: Remaining regions"}
            handles = [Line2D([], [], color=COLORS[m], marker="o", markersize=7, ls="none", label=legend_labels[m])
                       for m in MODELS]
            fig.legend(handles=handles, loc="lower center", bbox_to_anchor=(.585, .835),
                       ncol=2, frameon=False, fontsize=14, columnspacing=1.6,
                       handletextpad=.4, borderaxespad=0, labelspacing=.35)
            if simulated:
                fig.text(.18, .01, "SIMULATED DATA", fontsize=8, color="#616A72")
            _save(fig, output / f"{scheme}_all_model_comparisons.svg", preview)
            plt.close(fig)
    write_csv(pd.concat(all_plot_rows, ignore_index=True), output / "figure_comparison_data.csv")
