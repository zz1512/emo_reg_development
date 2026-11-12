# 年龄 BF 与 ERQ 三模型：结果读取和解释指南

适用版本：2026-11-10 新增的年龄贝叶斯样条比较，以及 `erq_age_smooth_pygam_v3` 的 ERQ 三模型分析。本文件说明读取顺序、字段、判定和论文报告方式；运行命令见 [README_BAYES.md](README_BAYES.md) 和 [README_ERQ.md](README_ERQ.md)。阅读本说明不需要重跑分析。

`outputs/fig03/inference_mock/` 中的结果仅为模拟数据，Word 标有 SIMULATED，不能作为真实发现。下文中的数值例子也仅用于解释。

## 1. 先区分两套模型

| 分析 | 模型 | 含义 |
|---|---|---|
| 原年龄分析及新增年龄 BF | M1：截距 + C | 无年龄项 |
| 同上 | M2：截距 + 线性年龄 + C | 线性年龄关联 |
| 同上 | M3：截距 + 年龄样条 + C | 允许非线性年龄关联，也包含线性情况 |
| ERQ 三模型 | E1 / baseline：截距 + s(age) + C | 控制年龄与协变量 |
| 同上 | E2 / additive：E1 + β_ERQ × ERQ_c | 增加恒定 ERQ 斜率 |
| 同上 | E3 / full：E2 + ERQ_c × g(age) | 允许 ERQ 斜率随年龄变化 |

四类协变量为性别、站点、父母教育和 meanFD。年龄 BF 保留历史年龄分析的数值站点编码；ERQ 三模型统一使用分类站点编码。ERQ 的 E1 在 ERQ 完整样本上重新拟合，不等于直接复制历史 M3 的结果。两套分析的样本数也可能不同。

432 版包括 parcel432、Yeo7、system4、Juelich6，合计 449 单元；232 版只分析 parcel232。Tian 杏仁核属于 parcel 家族，独立 Juelich6 属于 amygdala 家族，不要混为同一套分区。

## 2. 推荐阅读顺序和文件位置

以下目录相对于工程根目录。

| 内容 | 目录 |
|---|---|
| 原年龄 F/p/q、partial R² | `outputs/fig03/gam_pygam_volume_432_models3_FD05/statistics/` |
| 432 年龄 BF | `outputs/fig03/gam_pygam_volume_432_models3_FD05/bayes/` |
| 232 年龄 BF | `outputs/fig03/gam_pygam_volume_232_models3_FD05/bayes/` |
| ERQ 三模型 | `outputs/fig03/erq_age_pygam_432_models3_FD05/` |

实际运行时若指定了其他 `--output`，以该目录为准。

1. **确认来源和版本**：读 `run_metadata.json`，核对分析类型、源缓存、FD≤0.5、样本及状态计数；ERQ 必须是 v3。读取统计主表的 `status` 和 `n_obs`。
2. **确认 ERQ 样本整合**：读 `diagnostics/erq_sample_summary.json` 和 `erq_filter_summary.csv`；需要定位被试时再查逐行、逐被试清单。
3. **看全量比较 CSV**：年龄读 `bayes/age_bayes_comparisons.csv`；ERQ 读 `statistics/erq_sequential_comparisons.csv`。不要只从显著列表开始。
4. **看完整宽表和敏感性**：需要三模型拟合指标、E3 系数、原年龄统计或先验敏感性时，读相应宽表和明细。
5. **看图及 Word**：图用于理解形状和方向；Word 用于阅读、摘选和论文排版。精确数值以 CSV 为准。

`name` 在 parcel 家族中是图谱编号，在其他家族中是组名；跨表连接使用 `dimension + name`，不要只按编号、名称或行序连接。

## 3. 年龄 BF：先看比较，再看大小

### 文件怎么选

| 文件 | 用途 |
|---|---|
| `age_bayes_statistics.csv` | 每单元一行的完整宽表；主先验 BF、敏感性摘要及可合并的原 F/p/q |
| `age_bayes_comparisons.csv` | 每单元三行，适合筛选和摘表 |
| `age_bayes_prior_sensitivity.csv` | 每有效单元九行，查看九套先验下的证据和积分精度 |
| `age_bayes_manuscript_table.docx` | 全量 Word 三线表 |
| `unit_response_exclusions.csv` | 各单元共同排除的被试 |

432 版比较长表为 449×3=1347 行；232 版为 232×3=696 行，均不含表头。九套先验明细最多分别为 4041 和 2088 行；不可计算单元会保留在主表和比较表中，但没有有效积分明细。

### 三种比较各自回答什么

| `comparison` | 简称 | BF10 > 1 时相对支持 | BF10 < 1 时相对支持 |
|---|---|---|---|
| `m3_vs_m1` | BF31 | 含年龄函数的 M3 | 无年龄项的 M1 |
| `m3_vs_m2` | BF32 | 允许非线性的 M3 | 线性年龄模型 M2 |
| `m2_vs_m1` | BF21 | 线性年龄模型 M2 | 无年龄项的 M1 |

**BF32 小只表示数据更支持线性模型，不能直接写成没有年龄关联。** 要结合 BF31、BF21 及原年龄 FDR 结果解读。

例如 BF31=12、BF32=0.2，则 BF21=60：相对于无年龄模型，数据支持年龄关联；相对于更灵活的样条模型，数据更支持线性模型。这个例子不表示年龄曲线向上还是向下，方向仍需看曲线或斜率。

### 字段和证据标签

| 比较长表字段 | 含义 |
|---|---|
| `bf10` | `comparison` 中第一个模型相对于第二个模型的边际证据比 |
| `bf01` | 反方向 BF，数值允许时等于 `1 / bf10` |
| `log_bf10` | BF10 的自然对数；0 对应 BF=1 |
| `evidence` | 根据本工程 3、10 阈值生成的描述标签 |
| `sensitivity_log_bf10_min/max` | 九套先验下 log BF 的最小/最大值；不是置信区间 |
| `sensitivity_direction_stable` | 九套先验下是否始终在 BF=1 的同一侧；不表示都达到 BF≥3 或≤1/3 |
| `sensitivity_evidence_labels` | 九套先验中出现的证据标签集合 |
| `frequentist_*` | 若存在匹配的原统计表，则带入原年龄 F/df/p/q/partial R²；不是从 BF 换算 |

宽表使用带比较后缀的列，例如 `bf10_m3_vs_m1`、`log_bf10_m3_vs_m2`；长表把后缀放在 `comparison` 列，因此字段缩短为 `bf10`、`log_bf10`。

| BF10 范围 | `evidence` | 按本工程约定的解释 |
|---|---|---|
| ≥10 | `strong_support_full` | 较强支持比较中的第一个模型 |
| ≥3 且 <10 | `support_full` | 支持比较中的第一个模型 |
| >1/3 且 <3 | `inconclusive` | 区分两模型的证据不足 |
| >0.1 且 ≤1/3 | `support_reduced` | 支持比较中的第二个模型 |
| ≤0.1 | `strong_support_reduced` | 较强支持比较中的第二个模型 |

这里的 `full` 指当前比较的第一个模型，**不是固定指 ERQ full，也不总是 M3**。这些分级是描述约定，不是显著性阈值。BF 是证据比，不是模型为真的概率；后验模型赔率还取决于先验模型赔率，见 [JASP 贝叶斯分析报告指南](https://pmc.ncbi.nlm.nih.gov/articles/PMC8219590/)。

### 先验敏感性与计算状态

主结果使用线性和非线性先验尺度均为 0.5；另计算 `{0.25, 0.5, 1}` 的九套组合，定义见 [贝叶斯方法说明](README_BAYES.md)。

- 敏感性 log BF 范围跨过 0：比较方向随先验改变，应报告这一点，不能只挑最大的 BF。
- 范围始终大于 0 或小于 0：方向一致，但证据强度仍可能在不同标签之间变化。
- 明细中的 `integration_converged` 应为 True，`quadrature_relative_error` 和 `quadrature_check_log_difference` 均需≤`1e-6`。这表示积分数值检查通过，不表示数据支持某模型，也不表示模型拟合充分。
- 线性尺度 BF 出现 `inf` 或 0 时先检查 log BF，前者可能上溢，后者可能下溢。不能据此写成确切的无穷证据或零概率。

BF 不做 BH，本次使用它补充模型证据；原频率学发现仍按对应范围/比较的 q 值报告。不能把 BF≥3 称为“全脑校正后显著”。q 不显著也不自动等于支持无年龄模型：BF 接近 1 与 BF 明确支持 M1 是不同结果。

## 4. ERQ：两个主要比较都要读

### 优先打开这几张表

| `statistics/` 下文件 | 读取目的 |
|---|---|
| `erq_sequential_comparisons.csv` | 首选阅读表：449 单元×3比较=1347行，含两项主要比较和一项补充比较 |
| `erq_age_statistics.csv` | 完整宽表：三模型拟合指标、全部检验、E2 与 E3 的系数 |
| `erq_manuscript_table.docx` | 同一批单元的三线表；非显著和不可计算单元也保留 |
| `erq_table_definitions.json` | 字段来源、单位、每范围/每比较有效 BH 数 |
| `additive_vs_baseline_significant.csv` | 第一项比较 q<0.05 的单元 |
| `full_vs_additive_significant.csv` | 第二项比较 q<0.05 的单元 |
| `sequential_significant.csv` | 任一主要比较通过的单元，用于描述性随访绘图 |

旧名称 `erq_full_vs_baseline.csv`、`joint_erq_significant.csv` 仍存在，但现在属于补充联合比较，不包含两项逐步比较的全部结论。不要继续只靠该表判断本次 ERQ 结果。

### 对应列和判断

| 长表 `comparison` | 宽表首要列 | 问题 |
|---|---|---|
| `additive_vs_baseline` | `q_additive_vs_baseline` | 不加入交互时，ERQ 是否在年龄和协变量之外提供关联信息？ |
| `full_vs_additive` | `q_full_vs_additive` | 加入 ERQ 主效应后，ERQ 与激活的关联是否随年龄变化？ |
| `full_vs_baseline` | `q_full_vs_baseline` | ERQ 主效应与交互合起来，是否提供额外信息？补充比较 |

`q_full_vs_additive` 与 `q_interaction` 相同；`q_full_vs_baseline` 与 `q_joint_erq` 相同。它们是别名，不是重复的独立证据。

两项主要比较分别在 parcel432、Yeo7、system4、Juelich6 的全部有效单元中做 BH，以 q<0.05 判定。没有先后显著性门槛，不在摘选后的脑区内重新 BH，也不声称控制两项比较及四个范围合起来的总体 FDR。

| E2–E1 的 q<0.05 | E3–E2 的 q<0.05 | 可得结论 |
|---|---|---|
| 是 | 否 | 发现 E2 中的 ERQ 关联；未发现 ERQ 斜率随年龄变化的校正后证据，不能证明斜率恒定 |
| 否 | 是 | 发现 ERQ 关联随年龄变化；恒定斜率不显著不能否定这种交互，应读随龄 ERQ 斜率 |
| 是 | 是 | 两项比较均有证据；E2 斜率不能独自概括 E3 中各年龄的关联 |
| 否 | 否 | 两项预设比较均未达到校正阈值；不能证明无关联、无交互或两种发育过程彼此独立 |

补充联合比较与两项逐步检验的检验维数、分母模型不同，可能给出不同的显著性结果。若仅联合比较通过，就报告“整体 ERQ 相关项有证据”，不能据此指定主效应或交互显著。

### β、F、partial R² 分别怎么读

| 长表/Word 中的内容 | 含义 |
|---|---|
| `ERQ main effect` 的 `beta`、`se`、`ci_low/high` | **E2** 的 ERQ 恒定斜率、标准误和95% CI；单位为 contrast / 原始 ERQ 一分 |
| 主效应的 `t` / `statistic` | E2 的 t 检验；相应一自由度比较满足 t²=F |
| `Age by ERQ` 的 F、df1、df2 | E3 对 E2 的整体交互检验；F 不提供交互正负方向 |
| `Joint ERQ` 的 F、df1、df2 | E3 对 E1 的补充联合检验 |
| `p` | 对应比较未校正 p |
| `q`、`reject` | 对应范围和比较的 BH 结果及是否 q<0.05 |
| `partial_r2_reference` | 该步骤解释的受限模型剩余变异比例，`(RSS_reduced − RSS_target) / RSS_reduced` |

例如 E2 β=0.04 表示在 E2 模型控制年龄和协变量后，ERQ 每增加一分，contrast 平均增加0.04个原始单位；它不表示激活增加4%。β 的符号指向本项目的重评相对观看 contrast，不能单凭正负称为“调节能力更好/更差”。

非线性交互有多个基函数系数，主表交互 β 栏的“—”表示不适用，不是0。要看不同年龄的正负关系，读 E3 的 ERQ 斜率，而不是试图从整体 F 推方向。

宽表 `beta_erq_reference` 及相应 SE/CI/t 是 **E3 参考拟合在观测年龄分布上的平均 ERQ 斜率**，辅助校正列为 `q_erq`。它不同于主表使用的 E2 系数 `beta_erq_additive_reference`，也不等于某一指定年龄的斜率。

三线表 β/CI 和正式 F/p/q 均来自与 pyGAM 相同固定样条空间内的无惩罚 OLS 参考拟合。`rss_*`、`r2_*`、`edf_*`、`gcv_*`、lambda 和 `delta_r2_pygam` 等则描述实际 GCV 惩罚拟合；不要跨两套拟合混用 RSS 或自由度。参考 F 的适用假设包括独立、同方差误差，见 [statsmodels 官方说明](https://www.statsmodels.org/stable/generated/statsmodels.regression.linear_model.RegressionResults.compare_f_test.html)。

`partial_r2_reference=0.01` 表示解释了该受限模型剩余变异的1%，不是整个大脑变异的1%，也不是 ERQ 中介了1%的年龄效应。独立 GCV 选择后，实际 pyGAM 的 `delta_r2_*` 可能为负，不能直接把它当作参考检验的 partial R²。

95% CI 是未作全脑多重比较调整的区间，所以可能出现 CI 不含0、未校正 p<0.05，但 q≥0.05；最终发现判定仍按预设 q。Word 显示 `<.001` 不是精确等于0，全精度数字查 CSV。

## 5. ERQ 图和指定年龄斜率怎么读

| 文件 | 表示什么 |
|---|---|
| `statistics/erq_reference_slopes_at_ages.csv` | 全部有效单元在实际年龄范围内的6/8/10/12/14/16/18岁 E3 ERQ 斜率，来自 OLS 参考拟合，含 SE/95% CI |
| `statistics/conditional_age_curves.csv`、`figures/age_by_erq/` | 任一主要比较通过单元，在 ERQ P16/P50/P84 处的同一个 full pyGAM 条件年龄曲线 |
| `statistics/erq_slopes_by_age.csv`、`figures/erq_slope_by_age/` | 上述单元中，full pyGAM 的 `d(contrast)/d(ERQ)` 随龄变化及95%点态 CI |
| `statistics/figure_manifest.csv` | 单元、两项 q 与 SVG 文件的对应关系 |

三条条件曲线不是三个 ERQ 分组分别拟合的结果；共同平移中心化后保留 ERQ 水平之间的差异。斜率图纵轴是 ERQ 斜率，不是年龄的一阶导数。

指定年龄的参考斜率与绘图斜率分别来自 OLS 参考拟合和惩罚 pyGAM，数值不要求完全一致。它们用于理解方向和形状；逐年龄 CI 不含0不能直接当作又通过了全脑和年龄搜索校正的发现。绘图 CI 以选定平滑参数为条件，既不是个体预测区间，也不是同时置信带。没有任一步显著单元时，图清单和曲线表只有表头属于正常输出。

## 6. ERQ 数据整合和被试排除怎么追溯

1. `diagnostics/erq_id_crosswalk.csv`：按源数据行号核对匿名 ID → 原表前6位真实 ID → 规范化影像 ID。ID 来自原表，分数来自新表 `ERQ_CR`。
2. `diagnostics/erq_rows.csv`：查看每一条 ERQ 源行的原始分数、数值分数、转换状态、是否匹配影像、原来和本次是否纳入。包含不在影像名册中的源行。
3. `diagnostics/original_exclusion_summary.csv` 和 `original_exclusions_by_filter.csv`：追溯原年龄缓存已经排除的人和记录的原因。
4. `diagnostics/erq_filter_summary.csv`：按 `n_newly_excluded` 读人数流水；某人可能同时不满足多项条件，不能把各条件不合格总人数直接相加。
5. `diagnostics/erq_exclusions_by_filter.csv`：定位被试和各项不合格条件，`is_first_exclusion` 标记实际首次排除。
6. `statistics/unit_response_exclusions.csv`：进一步查看每个神经单元 contrast 非有限造成的排除，因此各单元 N 可能不同。

完整拟合还把样本审计保存到 `inputs/` 作为快照。ERQ 分数按一人一行整合，不对重复 ID 平均、不重新计算问卷条目总分。行号还原依赖两张源表已经确认的逐行对应，不能分别排序后再合并。

同一个单元应有 `n_obs = n_obs_baseline = n_obs_additive = n_obs_full`。E1 即使不使用 ERQ，也不能另外纳入 ERQ 缺失者，否则比较不再使用同一响应样本。

## 7. 空值、False 和失败状态

| `status` | 含义 |
|---|---|
| `ok` | 该入口的拟合或证据计算有效；不等于显著，也不代表全部统计假设已验证 |
| `insufficient_observations` | 完整样本太少或相对于模型秩缺少必要残差维度 |
| `constant_age_or_response` | 年龄 BF：年龄或响应没有变异 |
| `age_confounded_with_nuisance` | 年龄 BF：年龄空间与协变量混淆，不能辨识目标比较 |
| `zero_reference_residual_variance` | 年龄 BF：参考残差变异为零或数值上过小 |
| `rank_deficient_targets` | ERQ：年龄/ERQ/交互目标项不满足可辨识的秩要求 |
| `zero_residual_variance` | ERQ：响应无变异，或 full 参考模型残差变异数值上为零 |

q 为空且 `reject=False` 可能表示未能检验，而不是已经检验且不显著。先看 `status`，再看有效 p/q。失败单元仍占全量结果表的行，但不进入有效 p 值的 BH 家族；实际有效数查 `erq_table_definitions.json`。

年龄积分若未通过精度检查会明确报错并中止，不会把失败悄悄替换成 BF=0。不要把中止运行遗留的旧文件或不完整文件当作本次成功输出。

## 8. 论文结果怎么摘选和表述

### 年龄分析

从同一范围、同一单元摘取：原年龄检验的 F/df/q、partial R²、N；新增 BF31、BF32（可加 BF21）、对应 log BF、主先验尺度及敏感性范围。若未找到匹配的原频率学表，BF 长表的 `frequentist_*` 会空缺，应回查对应年龄结果目录，不要补0。

建议表述模板：

> 在控制协变量后，〔脑区〕的年龄关联检验为〔F、df、q、partial R²、N〕。预先规定先验下，M3 相对 M1 的 BF31=〔值〕，M3 相对线性模型 M2 的 BF32=〔值〕；九套先验下 log BF31 范围为〔最小值，最大值〕，证据方向〔一致/随先验改变〕。

若 q 不显著而 BF 也没有明确区分模型，应写“证据不足”；只有 BF 相对支持指定无效应模型时，才按该比较及先验范围描述“相对支持无年龄项模型”，仍不等于证明年龄效应绝对为零。

### ERQ 主效应和交互

每个脑区可在同一三线表中保留三行。主效应行报告 E2 β、SE、95% CI、t(df)、p/q、partial R²、N；交互行报告 E3–E2 的 F(df1,df2)、p/q、partial R²、N；整体行报告 E3–E1 的相应指标，并标明为补充比较。方法中注明分范围、分比较 BH 和参考检验空间。

建议表述模板：

> 控制平滑年龄效应及协变量后，〔脑区〕的恒定 ERQ 斜率为 β=〔值〕，95% CI〔区间〕，t(〔df〕)=〔值〕，q=〔值〕。加入 ERQ×年龄平滑交互后，模型比较 F(〔df1〕,〔df2〕)=〔值〕，q=〔值〕，partial R²=〔值〕。交互方向结合预定年龄处的模型斜率作描述。

若所有主要 q≥0.05，可写：

> 在各预设分析范围内分别进行 BH 校正后，未发现超出年龄及协变量的恒定 ERQ 关联或 ERQ×年龄交互证据。

这不能改写为“ERQ 与神经激活无关”或“二者是独立发育过程”。本轮 ERQ 没有 BF，不能拿年龄 BF 代替 ERQ 无效应证据。当前模型回答关联和年龄调节关系，不能直接证明 ERQ 中介了年龄效应或给出年龄效应被解释的比例。

## 9. 文档更新与执行

本次只新增阅读说明，并在三个 README 中增加入口，不改变模型或结果，不需要重跑。需要云端也能阅读时，同步以下文档到同一路径即可：

- `analysis/fig03/gam_pygam/README_RESULTS.md`
- `analysis/fig03/gam_pygam/README.md`
- `analysis/fig03/gam_pygam/README_BAYES.md`
- `analysis/fig03/gam_pygam/README_ERQ.md`

需要运行或重新导表时，使用对应方法 README 中已有的完整命令；仅阅读说明无需执行分析命令。
