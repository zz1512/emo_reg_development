**SIMULATED DATA ONLY：仅为模拟，不可填写论文。**

# PLS跨站点预测汇报表

实际拟合PFC定义：SIMULATED anatomical PFC；规则：SIMULATED rule。

| Task | N | Sites | Delta R2 | 95% CI | Valid B | Status |
|---|---|---|---|---|---|---|
| success | 24 | 3 | 0.750 | 0.669 to 0.835 | 1950 | Available |
| age | 24 | 3 | 0.188 | 0.103 to 0.272 | 1950 | Available |

## 四模型表现

| Task | Model | N | Features | R2 [95% CI] | MAE | RMSE | r |
|---|---|---|---|---|---|---|---|
| success | Covariates | 24 | 0 | -0.691 [-0.747 to -0.633] | 1.039 | 1.217 | 0.219 |
| success | PFC | 24 | 100 | -0.584 [-0.646 to -0.526] | 0.948 | 1.177 | 0.230 |
| success | Distributed | 24 | 432 | 0.166 [0.107 to 0.225] | 0.711 | 0.854 | 0.546 |
| success | Non-PFC | 24 | 332 | -0.204 [-0.262 to -0.145] | 0.794 | 1.027 | 0.418 |
| age | Covariates | 24 | 0 | -0.943 [-1.003 to -0.885] | 1.244 | 1.542 | -0.283 |
| age | PFC | 24 | 100 | -0.086 [-0.147 to -0.026] | 0.892 | 1.153 | 0.264 |
| age | Distributed | 24 | 432 | 0.102 [0.046 to 0.158] | 0.849 | 1.048 | 0.480 |
| age | Non-PFC | 24 | 332 | 0.275 [0.218 to 0.332] | 0.826 | 0.942 | 0.600 |

## 数值总结

success：N=24，3站点；全脑相对PFC的ΔR²=0.750，95%CI=0.669 to 0.835。达到预设主比较的正向增量判据；还需结合M2−M0判断是否超过协变量。 M2−M0的ΔR²=0.858，95%CI=0.774 to 0.939。

age：N=24，3站点；全脑相对PFC的ΔR²=0.188，95%CI=0.103 to 0.272。这是次要分析，区间作描述；supported=False不表示年龄预测无差异。 M2−M0的ΔR²=1.045，95%CI=0.963 to 1.124。

Outer and inner cross-validation use site holdout. Delta R2 is distributed minus PFC. Intervals use paired site bootstrap with full refitting. Confirmation applies only to success M2 minus M1, with 2000 requested replicates, at least 1000 valid and at least 95 percent valid, and a positive lower bound. This report also checks completion against saved replicate records. Age and other comparisons are descriptive. Negative R2 is retained; negative delta MAE or RMSE indicates lower error. No p, q or BF is inferred from these intervals.

两个任务可有不同N，同一任务的四模型使用相同被试和划分。主表2行、模型表现8行、全部预设比较10行。不同PFC定义分别运行，不能拼接成同一分析。

重评分数保持原方向；在评分方向与公式尚未核实时，不自动写成能力更强或更弱。年龄预测改善不等于解释重评能力的年龄差异，也不直接证明网络重组。
