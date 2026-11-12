"""Tables and SVG figures for one finished PLS output directory."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
matplotlib.rcParams["svg.fonttype"] = "none"
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from analysis.fig03.predictive_pls.models import MODELS, PRIMARY_COMPARISON
from analysis.fig03.predictive_pls.precheck import project_path

GUIDE = """# PLS 结果读取说明

本目录的数字来自同一版外层留站点预测。CSV 保留全精度，Word 三线表使用同一套格式化数字。

## 阅读顺序

1. `frozen/config_snapshot.json`：输入哈希、PFC 区数、scikit-learn 实际版本、分数方向尚未改写。
2. `precheck/`：两个任务的人数、站点和逐项排除原因。
3. `frozen/pfc_parcel_membership.csv`：432 区的面积比例和 PFC 标记。Tian 401–432 为非 PFC。
4. `oof_predictions.csv`：每名被试一次外层测试预测。
5. `model_performance.csv`、`model_comparisons.csv`：R²、MAE、RMSE、Pearson r，以及 bootstrap 区间（若已合并）。
6. `site_performance.csv`：各站点结果。目标恒定或人数少于 2 时，该站点 R² 为 NA，MAE 和 RMSE 仍保留。
7. `inner_cv_scores.csv`、`fold_models.csv`：内层候选 MSE 和选定的成分数。
8. `bootstrap_metrics.csv`：每次站点重抽样的有效性、原因和指标。

## 指标方向

- R² 可以是负数，表示外层预测比用测试集总均值更差。不要把负数改成 0。
- 主比较是重评成功的 `M2_distributed432-M1_PFC`。`delta_r2` 大于 0 表示全部分区的外层 R² 高于 PFC 分区。
- `delta_mae` 和 `delta_rmse` 用同一顺序相减。负数表示前一个模型的绝对误差更小。
- `prediction_error_improvement_vs_m0` 等于 `1 - SSE_M / SSE_M0`，是预测误差改善，不是因果解释比例。
- Pearson r 只描述预测和实测的线性关系，不用它的普通相关 p 值代替模型比较。

## 区间

区间是站点整组重抽样、每次完整重拟合后，有效重复的 2.5% 和 97.5% 分位数。
`inference_valid` 为真需要请求满 2000 次、有效不少于 1000 次且有效比例不少于 95%。
只有重评成功的 M2−M1，在上述条件成立且区间下界大于 0 时，`supported` 才为真。
年龄任务和其他比较保留区间，但不作为确认性发现。有效比例不足时仍可阅读描述性区间，此时 `supported` 为假。
无效重复保留失败原因，不改记成 0 效应，也不为了凑满次数而重新抽样。

## 不能这样写

- 未发现增量，不能写成 PFC 与全脑等效。
- 年龄预测变好，不能写成脑活动解释了重评能力的发育，也不能写成网络重组。
- M2 优于一个很差的 M1、但没有优于协变量模型 M0 时，不能写成神经指标提供了可靠的行为增量。
- 跨界分区的均值没有拆成 PFC 和非 PFC 体素。结论是超出预定义 PFC 分区集合的增量预测信息。
- 重评成功分数在方向核实前只称原分数，不写成越高能力越强。
"""


def _format(value) -> str:
    if pd.isna(value):
        return "NA"
    if isinstance(value, (bool, np.bool_)):
        return "true" if value else "false"
    if isinstance(value, (int, np.integer)):
        return str(int(value))
    if isinstance(value, (float, np.floating)):
        return format(float(value), ".12g")
    return str(value)


def _three_line(document, frame: pd.DataFrame):
    from docx.oxml import OxmlElement
    from docx.oxml.ns import qn
    from docx.shared import Pt

    headers = [str(column) for column in frame.columns]
    table = document.add_table(rows=1, cols=len(headers))
    borders = OxmlElement("w:tblBorders")
    for edge in ("top", "bottom", "left", "right", "insideH", "insideV"):
        element = OxmlElement(f"w:{edge}")
        element.set(qn("w:val"), "single" if edge in ("top", "bottom") else "nil")
        element.set(qn("w:sz"), "8")
        element.set(qn("w:color"), "000000")
        borders.append(element)
    table._tbl.tblPr.append(borders)
    for cell, text in zip(table.rows[0].cells, headers):
        cell.text = text
        for paragraph in cell.paragraphs:
            for run in paragraph.runs:
                run.bold = True
                run.font.size = Pt(8)
    for record in frame.itertuples(index=False):
        cells = table.add_row().cells
        for cell, value in zip(cells, record):
            cell.text = _format(value)
            for paragraph in cell.paragraphs:
                for run in paragraph.runs:
                    run.font.size = Pt(8)
    header_bottom = OxmlElement("w:tcBorders")
    line = OxmlElement("w:bottom")
    line.set(qn("w:val"), "single")
    line.set(qn("w:sz"), "4")
    header_bottom.append(line)
    for cell in table.rows[0].cells:
        cell._tc.get_or_add_tcPr().append(header_bottom)
    return table


def write_docx(analysis: Path, performance: pd.DataFrame, comparisons: pd.DataFrame) -> None:
    from docx import Document
    from docx.shared import Pt

    document = Document()
    style = document.styles["Normal"]
    style.font.name = "Times New Roman"
    style.font.size = Pt(10)
    document.add_heading("跨站点 PLS：宽范围解剖 PFC", level=1)
    document.add_paragraph(
        "重评成功是确认性终点，年龄是预先指定的次要分析。主比较为全部 432 分区相对 PFC 分区的外层 R² 增量。"
        "区间来自站点整组重抽样。负 R² 保留。Word 单元格与 CSV 使用同一格式。")
    document.add_heading("模型表现", level=2)
    _three_line(document, performance)
    document.add_heading("模型比较", level=2)
    _three_line(document, comparisons)
    document.add_paragraph(
        "supported 只可能出现在重评成功的 M2−M1，并且需要满额有效 bootstrap。"
        "MAE 与 RMSE 的增量为前一模型减后一模型，负数表示前一模型误差更小。")
    document.save(analysis / "tables_main.docx")


def _figures(analysis: Path, performance: pd.DataFrame, comparisons: pd.DataFrame,
             predictions: pd.DataFrame, sites: pd.DataFrame) -> None:
    figure_dir = analysis / "figures"
    figure_dir.mkdir(exist_ok=True)
    figure, axis = plt.subplots(figsize=(7.2, 4.2))
    offsets = np.arange(len(MODELS))
    width = 0.35
    for shift, task in ((-width / 2, "success"), (width / 2, "age")):
        frame = performance.loc[performance["task"].eq(task)].set_index("model").loc[list(MODELS)]
        yerr = None
        if frame["ci_low_r2"].notna().any():
            yerr = np.vstack([frame["r2"] - frame["ci_low_r2"], frame["ci_high_r2"] - frame["r2"]])
        axis.bar(offsets + shift, frame["r2"], width=width, label=task, yerr=yerr, capsize=3)
    axis.axhline(0, color="#444444", lw=0.6)
    axis.set_xticks(offsets, MODELS, rotation=15, ha="right")
    axis.set_ylabel("Outer-test R2")
    axis.legend()
    axis.spines["top"].set_visible(False)
    axis.spines["right"].set_visible(False)
    figure.tight_layout()
    figure.savefig(figure_dir / "model_r2.svg", format="svg")
    plt.close(figure)
    for task in ("success", "age"):
        frame = sites.loc[sites["task"].eq(task) & sites["model"].eq("M2_distributed432")]
        figure, axis = plt.subplots(figsize=(6.5, 3.8))
        axis.bar(frame["site_id"].astype(str), frame["delta_r2_M2_minus_M1"], color="#4C78A8")
        axis.axhline(0, color="#444444", lw=0.6)
        axis.set_ylabel("Site R2, M2 minus M1")
        axis.set_xlabel("Site")
        axis.spines["top"].set_visible(False)
        axis.spines["right"].set_visible(False)
        figure.tight_layout()
        figure.savefig(figure_dir / f"site_delta_{task}.svg", format="svg")
        plt.close(figure)
        observed = predictions.loc[predictions["task"].eq(task)]
        figure, axis = plt.subplots(figsize=(4.6, 4.2))
        for site, group in observed.groupby("site_id"):
            axis.scatter(group["y_true"], group["pred_M2_distributed432"], s=12, alpha=0.75, label=str(site))
        axis.set_xlabel("Observed")
        axis.set_ylabel("Outer-test M2 prediction")
        axis.legend(frameon=False, fontsize=7)
        axis.spines["top"].set_visible(False)
        axis.spines["right"].set_visible(False)
        figure.tight_layout()
        figure.savefig(figure_dir / f"oof_scatter_{task}.svg", format="svg")
        plt.close(figure)
    primary = comparisons.loc[comparisons["comparison"].eq(PRIMARY_COMPARISON), "role"]
    if primary.empty:
        raise RuntimeError("Model comparison table has no M2-M1 rows")


def write_report(analysis: Path, skip_docx: bool = False) -> None:
    analysis = analysis.resolve()
    performance = pd.read_csv(analysis / "model_performance.csv")
    comparisons = pd.read_csv(analysis / "model_comparisons.csv")
    predictions = pd.read_csv(analysis / "oof_predictions.csv")
    sites = pd.read_csv(analysis / "site_performance.csv")
    (analysis / "README_RESULTS.md").write_text(GUIDE, encoding="utf-8")
    _figures(analysis, performance, comparisons, predictions, sites)
    if not skip_docx:
        write_docx(analysis, performance, comparisons)
    print(f"Result guide: {analysis / 'README_RESULTS.md'}", flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--analysis", type=Path, required=True)
    parser.add_argument("--skip-docx", action="store_true")
    args = parser.parse_args()
    write_report(project_path(args.analysis), args.skip_docx)


if __name__ == "__main__":
    main()
