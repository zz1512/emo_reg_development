# 两套联合 PLS：结果读取说明

执行与方法见 [README](README.md)。正式根目录默认为 `outputs/fig03/joint_pls_432_FD05_both/`；`mock` 路径下全部数值均为模拟，不能用于论文。

## 先核实样本

`precheck/sample_metadata.json` 的 `n`、`sites` 是重新预检后的共同样本量和站点数，两个结局、四个模型及两套方案完全一致。`parent_education_used` 必须为 false。不要把旧分析的 N 直接抄到新表。

| 文件 | 读取用途 |
|---|---|
| `precheck/upstream_cohort_audit.csv` | 从原始来源重建时：人口学表与有效名单、二级排除名单的纳入关系 |
| `precheck/sample_roster.csv` | 每位有效名单被试的旧状态、年龄、成功、FD、站点及全部排除原因 |
| `precheck/filter_summary.csv` | 每项过滤命中总数、按顺序新排除数和剩余数；只有新排除数可累加 |
| `precheck/exclusions_by_filter.csv` | 哪一项具体排除了谁；同一人可命中多项 |
| `precheck/behavior_alignment.csv` | 行号、旧 ID、真实 ID、原始成功分数及数值解析状态；行为表未匹配的人也保留 |
| `precheck/common_sample.csv` | 最终共同分析被试名单，含年龄＋成功，不含脑区大矩阵 |
| `precheck/site_summary.csv` | 各站点 N、年龄范围、成功分布 |
| `precheck/requires_brain_recovery.csv` | 其余条件通过但旧缓存缺少脑区值的人 |
| `precheck/brain_recovery_log.csv` | 若发生恢复：原始 NIfTI、图谱与哈希及恢复状态 |

从完整旧设计恢复时，`education_exclusion_reconsidered` 记录之前仅因教育被排除的人；原始重建路径没有教育筛选，因此该列不会把这些人标为旧教育恢复，实际脑区恢复以 recovery log 为准。没有恢复需求时恢复文件不生成。

## 两套结果分别是什么

`scheme=joint`：一次 PLSRegression 同时拟合年龄和成功，两列经过各训练集的协变量调整和标准化。R² 为两个结局的 pooled 外层 R² 的等权平均，不是分别做年龄和成功再比较两个样本。

`scheme=shared`：共同目标为每折训练集全脑参考 PLS 第一 Y 方向得分。该折四个模型的 `shared_true` 相同；预测的二维 Y 投影到同一方向得到 `shared_pred`。它是脑信息参与定义的行为方向，不能视为独立纯行为因子。

R²、MAE、RMSE 的主要报告空间均为**调整后的标准化残差/分数**；M0 在这个空间预测零。原始单位年龄和成功预测值在 OOF 文件中保留，包含训练协变量模型的预测。因而本目录 R² 不能直接与旧版原始单位单目标总预测 R² 混用。

## 主要表：先看两行，再看八行

`report/primary_comparisons.csv`：**2 行**，每套方案各一行 `distributed_minus_pfc`。建议正式汇报：方案、共同 N、站点、PFC R²、Distributed R²、ΔR²、置换 P、BH q。两个单模型 R² 从八行表现表提取。

| 字段 | 含义 |
|---|---|
| `delta_r2` | 全脑减 PFC，正值为观察到的预测增益；本身不是显著性 |
| `p_permutation` | 对应约简模型零分布的单侧近似置换 P |
| `q_primary` | 两套主比较一起做 BH；仅主比较有 q |
| `permutations_requested/valid` | 正式应为 1000/1000 |
| `inference` | 未运行 `not_run`、未完成 `incomplete`、或完整近似条件置换状态 |
| `supported_primary` | 完整主比较且 ΔR²>0、q<0.05 才 true；推断未就绪时留空 |

`model_performance.csv` / `report/model_performance_full_precision.csv`：**8 行**（2 方案×4 模型），字段含 N、站点数、脑区数、R²、MAE、RMSE。负 R² 保留，表示不如评价空间中的常数均值参照，不得截成零。

`model_comparisons.csv` / `report/model_comparisons_full_precision.csv`：**10 行**（2×5 预设比较）：全脑−PFC、PFC−协变量、全脑−协变量、非PFC−协变量、全脑−非PFC。负 ΔMAE/ΔRMSE 表示误差降低。只跑默认 `--null pfc` 时，其他比较 P 留空；它们不是“全部不显著”，而是尚未执行相应零分布。

## 成分与泛化诊断

- `reference_axes.csv`：每个外层训练集的共同参考 q、年龄/成功载荷、残差相关、是否同向。年龄系数统一定向为正。若成功权重接近零或反向，应如实解释成分内容，不预先称之为正向共同发展。
- `model_first_components.csv`：PFC、全脑、非PFC 各自第一成分的行为权重、行为方向载荷、脑得分与两行为的相关。是**全样本描述性**拟合，不能当作外部预测显著性。
- `model_first_component_weights.csv`：上述三模型逐脑区第一成分权重，脑区编号对应 432 分区；不同模型学习的是不同方向，权重不能直接当作共同空间中的效应量。
- `outer_folds.csv`：各站点、各方案、各模型的测试表现、选中的 k、实际特征数。
- `inner_scores.csv`：各 k 的内层 pooled R²、负 R² 选择损失及标准化 MSE。选中的 k 可以为 0，表示脑活动未带来内层可泛化增益。
- `oof_predictions.csv`：每人×4模型×2方案，共 8N 行。含年龄/成功真实值、原始单位预测值、调整后的标准化真实值与预测值、共同分数真实值与预测值。所有正式表现取这些外层留出预测，不能取全样本拟合相关。

## 置换与 bootstrap 记录

`resampling_status.csv` 检查请求、完成、有效重复数。未运行不产生正式 P；失败不会自动作为 false 或被丢弃后重新计算较好看的 P。

`permutations/pfc/` 是主零分布；可选 `covariates/`、`nonpfc/` 分别对应 README 中说明的探索性比较。`manifest.json` 固定参数、零假设与源哈希，`rep_*.csv` 记录每次状态与 ΔR²，`indices_*.npy` 保存确切站点内置换索引。

`bootstrap/independent/` 与置换互不嵌套。每次从各站点有放回抽原站点 N_s 人，总数固定 N。输出：

- `bootstrap_intervals.csv`：ΔR² 的 2.5%–97.5% 百分位区间及完成状态。它用于稳定性评估，不生成置换 P。部分完成时不得当作预设 1000 次最终区间。
- `brain_first_component_stability.csv`：三模型全样本第一成分脑权重、bootstrap SE 与 ratio，符号与原始权重对齐；它不代表多成分最终预测模型的全部贡献，也不是逐脑区校正显著性。
- `reference_axis_stability.csv`：共同参考方向的年龄/成功载荷及 bootstrap 变化；用于辅助解释共同方向是否稳定。

## 可直接查看的文件

`report/response_tables.docx` 为两页三线表：每套方案一页，四模型表现＋五项模型比较＋英文结果段落。`report/result_summary.md` 含相同结果和更完整的方法边界；`joint_model_performance.svg`、`shared_model_performance.svg` 为各自表现图，仅输出 SVG。

正式结果段落只有在主置换完成后才会写“支持/未支持正向增量”；未运行或未完成会明确写推断待完成。即使 P/q 不显著，也不能据此证明 PFC 与非 PFC 完全等效或发展相互独立。

英文结果必须限定为：在指定约简线性均值、站点内残差可交换和嵌套验证流程下得到的条件信息增量证据；不能写成机制、因果或不受模型假设限制的显著预测优势。
