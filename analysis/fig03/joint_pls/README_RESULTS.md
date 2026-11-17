# 两套联合 PLS：结果读取说明

执行与方法见 [README](README.md)。修复迭代上限后的 meta 正式根目录为 `outputs/fig03/joint_pls_432_FD05_meta_pfc_both_iter10000/`；旧的 80 区结果在 `outputs/fig03/joint_pls_432_FD05_both/`，两者不能混用。`mock` 路径下全部数值均为模拟，不能用于论文。

## 先核实样本

`precheck/sample_metadata.json` 的 `n`、`sites` 是重新预检后的共同样本量和站点数，两个结局、四个模型及两套方案完全一致。`parent_education_used` 必须为 false。不要把旧分析的 N 直接抄到新表。

其中 `pfc_definition.source` 应为 `Buhle2014_reappraisal_PFC`，`n_pfc=10`、`n_remaining=422`。旧 80 区版本为 `Donahue2018_liberal_bilateral`。`precheck/pfc_parcel_membership.csv` 是本次冻结使用的具体成员，报告开头也会标注定义和脑区数。

Meta 报告中的 `Remaining regions` 是其余 422 个分区，可能包括未被元分析定义选中的解剖 PFC；不能解释为严格的解剖非前额叶。CSV 内部模型名 `M3_nonPFC` 为兼容性保留。

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
| `q_primary` | 两套 M2−M1 主比较一起做 BH，仅主比较填此列 |
| `q_secondary` | 两套 M3−M1 次要比较一起做 BH，仅次要比较填此列 |
| `q_exploratory` | M1−M0、M2−M0、M3−M0 × 两套方案，共六项一起做 BH |
| `q_bh` | 该行所属检验族的 BH q，供统一读取和绘图；不代表10项统一校正 |
| `bh_family_size` | 对应族的预定检验数：主比较2、次要比较2、探索性比较6；描述性比较0 |
| `permutations_requested/valid` | 正式应为 1000/1000 |
| `inference` | 未运行 `not_run`、未完成 `incomplete`、或完整近似条件置换状态 |
| `supported_primary` | 完整主比较且 ΔR²>0、q<0.05 才 true；推断未就绪时留空 |

`model_performance.csv` / `report/model_performance_full_precision.csv`：**8 行**（2 方案×4 模型），字段含 N、站点数、脑区数、R²、MAE、RMSE。负 R² 保留，表示不如评价空间中的常数均值参照，不得截成零。

`report/pfc_comparisons.csv`：**4 行**，两套方案各有全脑−PFC（M2−M1）与剩余脑区−PFC（M3−M1）。这是本次直接比较 PFC 时优先读取的表。`q_primary` 只校正原来的两项 M2−M1；`q_secondary` 只校正新加的两项 M3−M1。两个检验族分别 BH，不是四项统一校正。`supported_secondary` 为新比较的 ΔR²>0 且次要 q<0.05 的状态，推断待完成时为空。

`report/model_comparisons_full_precision.csv`：**10 行**（2×5），每套按 M1−M0、M2−M0、M3−M0、M2−M1、M3−M1 排列。负 ΔMAE/ΔRMSE 表示误差降低。`--null pfc` 可恢复 M2−M1 和 M3−M1 的各自 P；三种相对 M0 的比较只有另跑 `--null covariates` 才有 P。M3−M2 / M2−M3 不是分析目标，不出现在新报告中。根目录 `model_comparisons.csv` 与 `observed_comparisons.csv` 为旧计算原始记录；正式汇报使用 `report/` 新表。

相对 M0 的六项比较现在一起输出 `q_exploratory`；两套方案共同进入同一个探索性六检验族，不是每张图各自校正三项。`q_bh` 汇集各行相应的 `q_primary/q_secondary/q_exploratory`，图和三线表统一读取此列。原主/次比较族保持不变，三个族分别校正，不提供全部10项比较的统一 FDR 保证。`supported_exploratory` 仅在探索性 q 就绪且 ΔR²>0、q<0.05 时为 true。若族内任一 P 未就绪，整族 q 留空，不把未完成当作不显著。补 q 只需重新 report，无需重跑已完成的置换。

新增统计量不是沿用 M2−M1 的 P。新导出器逐次计算 `(M3−M0)−(M1−M0)`，并使用 PFC 约简零分布下该差值计算近似条件 P。其解释受既有线性约简均值与可交换性假设限制，不是两个非嵌套算法普遍相等的检验。M3−M1 也不是单个脑区显著性；它比较剩余脑区集合和 PFC 集合。真正的合并增量问题仍看 M2−M1。

## 成分与泛化诊断

- `reference_axes.csv`：每个外层训练集的共同参考 q、年龄/成功载荷、残差相关、是否同向。年龄系数统一定向为正。若成功权重接近零或反向，应如实解释成分内容，不预先称之为正向共同发展。
- `model_first_components.csv`：PFC、全脑、非PFC 各自第一成分的行为权重、行为方向载荷、脑得分与两行为的相关。是**全样本描述性**拟合，不能当作外部预测显著性。
- `model_first_component_weights.csv`：上述三模型逐脑区第一成分权重，脑区编号对应 432 分区；不同模型学习的是不同方向，权重不能直接当作共同空间中的效应量。
- `outer_folds.csv`：各站点、各方案、各模型的测试表现、选中的 k、实际特征数。
- `inner_scores.csv`：各 k 的内层 pooled R²、负 R² 选择损失及标准化 MSE。选中的 k 可以为 0，表示脑活动未带来内层可泛化增益。
- `oof_predictions.csv`：每人×4模型×2方案，共 8N 行。含年龄/成功真实值、原始单位预测值、调整后的标准化真实值与预测值、共同分数真实值与预测值。所有正式表现取这些外层留出预测，不能取全样本拟合相关。

## 置换与 bootstrap 记录

`resampling_status.csv` 检查请求、完成、有效重复数。未运行不产生正式 P；失败不会自动作为 false 或被丢弃后重新计算较好看的 P。

`completed=1000` 只表示 1000 次都返回了状态，不能代替 `valid=1000`。例如置换 `completed=1000, valid=989` 表示 11 次失败，P/q 留空；bootstrap `valid=994` 表示 6 次失败，最终权重稳定性汇总未就绪。进度日志的最后一行 `1000/1000 ok` 仅指编号 1000 的那一次成功，不代表前面每一次都成功。

已确认的旧批次错误是 `ConvergenceWarning: Maximum number of iterations reached`。新版将 PLS 上限从默认 500 增至 10000，容差仍为 `1e-6`；读取 `analysis_manifest.json` 的 `solver` 可核对。修复后用新的 meta 输出目录重新运行，先核对 `PFC=10`。不能把 80 区结果改名当作 meta，也不能忽略失败记录直接用 989 次生成正式 P。正式有效数仍以服务器运行后实际输出为准。

`permutations/pfc/` 用于 M2−M1 和新 M3−M1；可选 `covariates/` 对应三个脑模型相对 M0 的探索性比较。本次不使用 `nonpfc/`。`manifest.json` 固定参数、零假设与源哈希，`rep_*.csv` 记录每次状态与 ΔR²，`indices_*.npy` 保存确切站点内置换索引。新报告直接核对这些原始逐次记录，不以旧聚合文件中的“complete”代替检查。

`bootstrap/independent/` 与置换互不嵌套。每次从各站点有放回抽原站点 N_s 人，总数固定 N。输出：

- `bootstrap_intervals.csv`：ΔR² 的 2.5%–97.5% 百分位区间及完成状态。它用于稳定性评估，不生成置换 P。部分完成时不得当作预设 1000 次最终区间。
- `brain_first_component_stability.csv`：三模型全样本第一成分脑权重、bootstrap SE 与 ratio，符号与原始权重对齐；它不代表多成分最终预测模型的全部贡献，也不是逐脑区校正显著性。
- `reference_axis_stability.csv`：共同参考方向的年龄/成功载荷及 bootstrap 变化；用于辅助解释共同方向是否稳定。

## 可直接查看的文件

`report/response_tables.docx` 为每套方案的三线表：四模型表现＋五项模型比较＋英文结果段落。`report/result_summary.md` 含相同结果和完整的方法边界。

默认输出 **joint 和 shared 各自独立的一张图，每张五行比较**：`report/joint_all_model_comparisons.svg`、`report/shared_all_model_comparisons.svg`。模型身份固定使用高饱和度配色：M0 蓝（#0072FF）、M1 绿（#00A651）、M2 紫（#8B2CFF）、M3 橙（#FF7A00），跨所有比较和方案一致。左图为两个模型 R² 点估计和灰色连接线；**单模型端点CI默认关闭**。右图保留可用的 ΔR² 配对95%CI，直接在对应差值点上方**统一标注 BH q**，读取 `q_bh`，以 `BH-adjusted q` 标明其性质；无就绪 q 时不标注，不回退显示 P。CSV 和 Word 同时保留原始 P 与相应族的校正后 q。正式图不再单设统计列或重复列出估计和区间数值，完整统计与方法说明保留在表格中。显著性不由单模型区间重叠与否自动判定。

只有需要单模型端点CI时，才执行 README 中 `bootstrap_endpoints` 恢复命令，并在 report 命令加 **`--model-ci`**。默认不读取已有端点区间。恢复按原索引重建外层 M0 评估目标，再用三个 model−M0 差值恢复绝对 R²；shared 仅重建外层参考第一 Y 方向，不重新做内层调参或候选模型拟合。`bootstrap_model_r2_intervals.csv` 含8行区间，另有逐次绝对 R² 及来源哈希文件。端点CI和差值CI使用同一批抽样，不能直接相减两个端点区间来替代配对差值区间。

每套图的 `*_all_model_comparisons_data.csv` 含该图五行全精度端点、端点区间、差值、差值区间和 P/q。旧版合并图只是历史产物，默认不再生成或刷新，请使用新的两张独立图。

`figure_comparison_data.csv` 为主图逐项全精度数据；`paired_permutation_statistics.csv` 与 `paired_bootstrap_statistics.csv` 为逐次配对恢复后的统计量；`comparison_provenance.json` 记录派生公式、校正范围、请求/完成/有效数及失败记录。失败或未完成时不输出最终 P 或最终区间。正式图默认仅 SVG，`--preview` 额外生成 PNG 供预览。

正式结果段落只有在主置换完成后才会写“支持/未支持正向增量”；未运行或未完成会明确写推断待完成。即使 P/q 不显著，也不能据此证明 PFC 与非 PFC 完全等效或发展相互独立。

英文结果必须限定为：在指定约简线性均值、站点内残差可交换和嵌套验证流程下得到的条件信息增量证据；不能写成机制、因果或不受模型假设限制的显著预测优势。
