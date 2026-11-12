# ERQ 与年龄的神经激活模型

结果出来后如何阅读：见 [结果读取和解释指南](README_RESULTS.md)，包括两个主要比较的判定矩阵、E2/E3 系数区别、三线表字段、ERQ 样本整合追溯及论文报告模板。

本分析在原 Tian432 的 m3 年龄模型上加入习惯性认知重评使用（ERQ_CR）及其与年龄的平滑交互。2026-11-10 升级为三模型：baseline、additive、full，在完全相同的有效样本上分别拟合。主要比较为 `additive_vs_baseline` 和 `full_vs_additive`；两项分别在四个范围内做 BH，保留全部 449 单元及 Word 三线表/全精度 CSV。ERQ 暂不计算 BF。

## 模型

```text
baseline (m3 structure): Y = intercept + f(age) + C
additive:               Y = intercept + f(age) + beta_ERQ * ERQ_c + C
full:                   Y = intercept + f(age) + beta_ERQ * ERQ_c + ERQ_c * g(age) + C
C = gender + categorical(site_id) + parent_edu + meanFD
ERQ_c = ERQ_CR - mean(ERQ_CR)
```

`Y` 与原分析相同：重评相对观看的神经激活差值。`f(age)` 与原 m3 一样使用 5 个三次 B 样条；`g(age)` 使用相同形式的年龄基函数，但有独立系数，描述 ERQ 与激活的关联如何随年龄改变。这里的交互对应用户指定的 `s(age) * l(erq)`，不是先分别拟合两条曲线再相乘。

baseline 的 pyGAM 特征按 `age_c, nuisance...` 排列，使用 `s(0, n_splines=5)` 加不惩罚的协变量和截距。full 的特征按 `age_c, ERQ_c, nuisance...` 排列，模型项为：

```python
s(0, n_splines=5) + l(1, penalties=None) + s(0, n_splines=5, by=1)
# 再加入不惩罚的协变量和截距
```

`by=1` 使年龄样条乘以 ERQ_c，见 [pyGAM SplineTerm 文档](https://pygam.readthedocs.io/en/latest/reference/_autosummary/pygam.terms.SplineTerm.html)。ERQ 保持线性；三个模型都用 GCV 选择平滑参数，搜索范围为 `1e-5` 至 `1e5` 的 12 个值。baseline 和 additive 分别搜索 **12 个年龄惩罚参数**；full 的年龄项和交互项独立搜索 **144 种组合**。三个模型不强制使用相同最佳年龄惩罚值；ERQ、截距和协变量不加惩罚。additive 项为 `s(0,n_splines=5) + l(1,penalties=None)` 加协变量。

年龄、ERQ 和协变量设计列在每个单元的共同完整样本内中心化，保留原始单位。三个模型使用相同年龄样条空间和中心化参考。站点统一使用分类固定效应，即去掉一个参考站点的哑变量，不把站点编号之间的数值距离当作连续量。

这里的 baseline 沿用原 m3 的模型结构，但会在 ERQ 有效样本上重新拟合，不直接读取历史 m3 的 RSS、EDF 或统计结果。旧年龄入口按数值列处理站点，本次三个模型均按分类站点编码，因此 baseline 不应被视为历史 m3 输出的逐值复现。历史年龄分析文件保持不变。

## 数据和样本

ERQ 读取接口参考 `analysis/fig02/run.py`。本次使用新表的分数和用户已确认逐行对应的原表 ID：

- 分数表为 `/public/home/dingrui/fmri_analysis/data/beh/er_profiles_new.csv`（CLI 默认值），使用 `sub-0001` 等匿名顺序编号。其 ID 不能直接与影像的字母数字 ID 匹配。
- 云端实际文件的分数列为 **`ERQ_CR`**，使用 `--erq-column ERQ_CR`（也是 CLI 默认值）。不要截短为 `CR`，也不要用 `relative_reappraisal_success` 或 `reappraisal_success` 任务表现列替代 ERQ。其他文件若采用不同的量表列名，应通过 `--erq-column` 显式指定。
- 云端实际文件的 ID 列为 **`sub_ids`**，命令显式指定 `--erq-id-column sub_ids`。不指定时依次识别 `sub_id`、`sub_ids`、`subject`。列名不匹配时会分别说明缺少的是 ID 还是量表列，并显示实际可用列及对应参数。

通过 `--erq-id-source "/public/home/dingrui/fmri_analysis/data/beh/er_profiles.csv"` 指定原表。**用户已确认这两张表的数据行逐行对应，原表每条数据行的前 6 位是真实 `sub_id`。** 脚本跳过原表的一行表头，在任何排序或筛选前逐行还原 ID，再通过真实 ID 与影像表一对一合并。两表非空数据行数必须一致；还原出的 ID 必须是 6 位字母数字且规范化后唯一。原表仅提供 ID，分数始终取新表的 `ERQ_CR`，不会改读原表的 `ERQ-CR`。两张源表保持只读。

此参数明确依赖已经核实的行对应关系；脚本可检查行数和 ID，但不能单凭数据证明两表行序正确。因此，不要先独立排序、去重或删掉任一源表的行。原表允许制表符、逗号或整行连接的内容，只提取每条非空数据行开头的 6 位 ID。分数表自动识别逗号/制表符。不传该参数时，仍按分数表自身的 ID 匹配。

原分析的 `design.csv` 决定影像被试范围、年龄及四类协变量。先保留原 `in_analysis` 样本，再匹配 ERQ、剔除缺失值。ERQ 文件内的 age/sex 不覆盖影像设计表。

被试 ID 按字符串匹配：去首尾空格、统一大小写、兼容 `sub-` 前缀和纯数字 ID 的 `.0` 后缀，保留前导零。规范化后重复的 ID 会报错，不自动平均。原来排除的被试不会重新纳入；是否采用 FD≤0.5 由所选缓存的原分析纳入标记决定。

本次“ERQ 聚合”是被试级数据整合：每个被试取新表的一行量表总分，保留原始分数字符串，再转换为数值并合并影像表；不重新计算问卷条目总分、不对重复被试取平均。匿名 ID 或还原 ID 重复时直接报错，防止多对多合并扩大样本。

样本审计同时保存以下内容：

- 原筛选：保留 `design.csv` 全名册、`original_in_analysis` 和原 `excluded_by`，按原原因分别列出被试和人数。原缓存只记录首个原因时，本次如实继承，不能反推出缓存中未记录的其他历史原因；原名册形成前已被删掉的人也不在本次审计范围。
- 本次筛选依次为：无匹配 ERQ ID → ERQ 分数缺失/非数值/非有限 → 年龄缺失/非有限 → 性别缺失/非有限 → 父母教育缺失/非有限 → meanFD 缺失/非有限 → 性别不为 0/1 → 站点缺失。原排除被试不重新纳入，FD 阈值沿用原缓存。
- `erq_filter_summary.csv` 记录每步筛选前人数、实际新增排除人数、剩余人数，以及原纳入样本中该条件不合格的总人数。`erq_exclusions_by_filter.csv` 列出每项不合格条件对应的真实被试 ID；同一人可出现多次，`is_first_exclusion` 标记首个实际排除步骤。因此人数流水用 `n_newly_excluded`，不要将所有条件不合格人数直接相加。
- `erq_rows.csv` 保留**全部 ERQ 源数据行**，包括不在影像名册中的人：数据行号（不含表头，从 1 开始）、匿名 ID、还原 ID、规范化 ID、原始分数、数值分数、分数转换状态、影像匹配及原/本次纳入状态。`erq_score_status` 区分 `valid`、`missing`、`not_numeric` 和 `nonfinite`。不在影像名册中的源行标为 `not_in_imaging_roster`；影像名册里没有 ERQ 的人列在全名册及排除清单中。
- `erq_sample_summary.json` 记录来源路径、还原方式、源行数、各分数转换状态人数、影像内/外源行数、原纳入/排除的源行数和最终保留数，说明采用一人一行、不平均的整合规则。

每个神经单元再剔除自身激活缺失的行，逐单元、逐被试写入 `statistics/unit_response_exclusions.csv`；baseline、full 及全部参考比较均使用同一批有效行。即使 baseline 本身不含 ERQ，也不额外纳入 ERQ 缺失的被试。样本不足、无残差变异或目标项无法辨识的单元保留在统计主表中，记录 `status`，拟合及检验指标留空；`n_obs`、`n_obs_baseline`、`n_obs_full` 记录同一个共同有效样本的行数。

需要完整的 Tian432 缓存，包含 `design.csv`、`parcel_response.csv`、`parcel_info.csv`、`group_response.csv`、`amygdala_response.csv`、`tian_combined.json` 和 `combined_atlas.nii.gz`。不重新读取原始被试 effectSize 图像。本入口不支持 232 版。

## 检验与校正

两项预先指定的逐步比较都计算和校正；不以第一项显著或整体联合检验显著作为第二项的门槛。不能只凭包含年龄的总体模型显著判断 ERQ 有作用。

| 字段后缀 | full 对应的受限模型/零假设 | 通常分子自由度 | 用途 |
|---|---|---:|---|
| `additive_vs_baseline` | additive 与 baseline 比较 | 1 | 主要比较一，控制年龄/协变量后的恒定 ERQ 斜率；分母自由度来自 additive |
| `full_vs_additive` | full 与 additive 比较 | 4 | 主要比较二，ERQ 斜率随年龄变化；分母自由度来自 full；与 `interaction` 完全相同 |
| `full_vs_baseline` | 与 baseline：`intercept + f(age) + C` 比较 | 5 | 补充联合比较；`joint_erq` 是同一检验的兼容名称 |
| `joint_age_erq` | 与 `intercept + C` 比较 | 9 | 年龄、ERQ 和交互的整体检验；可由年龄单独驱动 |
| `age` | 平均 ERQ 处的年龄函数为零，保留交互 | 4 | full 内的年龄项 |
| `erq` | 观测年龄分布上平均的 ERQ 斜率为零，保留交互 | 1 | full 内的 ERQ 主效应 |
| `interaction` | 与 `intercept + f(age) + ERQ_c + C` 比较 | 4 | ERQ 斜率是否随年龄变化 |

为避免 GCV 有效自由度差接近零和主效应/交互基函数重叠引起的伪显著，全部检验使用 **相同固定样条空间内的无惩罚 OLS 参考拟合**。将年龄基函数在实际观测年龄上中心化，去掉最后一列，得到 4 个独立年龄列；与 ERQ_c 的乘积提供 4 个交互列。这样 ERQ 主效应与交互常数部分不会混淆。三套参考空间均检查与实际 pyGAM 空间一致。

三个模型的实际惩罚拟合指标用于描述拟合结果；正式 F/p/q 使用对应参考模型的 RSS 和矩阵秩，不把独立 GCV 后的 RSS、EDF 差直接代入普通 F。代码验证三套实际模型与参考模型张成相同空间。additive 的 beta/SE/CI/t 来自 additive 自己的残差方差，不能用 full 的分母自由度替代；其 t² 等于 additive vs baseline 的 F。

嵌套 F 由参考拟合的 RSS 和矩阵秩计算，沿用独立、同方差高斯误差假设，见 [statsmodels 检验说明](https://www.statsmodels.org/stable/generated/statsmodels.regression.linear_model.RegressionResults.compare_f_test.html)。它检验预设样条空间中的效应，不等同于 mgcv 的 EDF 检验。也不直接采用 pyGAM 自带的单项 p 值；其估计平滑参数后的显著性局限见 [pyGAM Quick Start](https://pygam.readthedocs.io/en/latest/notebooks/quick_start.html)。

四个分析范围分别为：**parcel432、Yeo7、system4、Juelich6**，合计 449 个单元。Tian 杏仁核分区仍属于 parcel432；Juelich6 为 LB/CM/SF × 左右。每个独立检验在每个范围的全部有效单元中做 BH，阈值 `q < 0.05`；不合并范围，不在筛选后的单元里重新校正，不代表跨范围/比较的总体 FDR 控制。

两项主要比较 **分别 BH**，不将它们合并为 864 个脑区检验，也不声称控制两项比较合起来的整体 FDR。联合比较和 full 内辅助项另外按范围校正；别名不重复计数。全部单元进入完整结果表，无显著性筛选。任一主要比较 `q < 0.05` 的单元可生成同一 full 模型的条件年龄曲线和 ERQ 斜率图，不另拟合低/中/高 ERQ 子组；它们是描述性随访，不增加确认性的分项结论。

`q_joint_erq` 显著不代表交互显著，仍须看 `q_interaction`。联合检验不显著表示本次分析未发现超出年龄和协变量的 ERQ 相关证据，不能据此证明两个发育过程独立。本模型可回答关联及调节关系，不能单独证明 ERQ 中介或解释了多少年龄效应；若要回答后者，需要后续同样本的调整前后比较。

## 输出

| 工程输出目录内的文件 | 内容 |
|---|---|
| `inputs/analysis_design.csv` | ERQ 匹配后的完整预测变量样本，保留原始单位 |
| `inputs/erq_sample_roster.csv` | 原被试名册、ERQ、纳入标记和首个排除原因 |
| `diagnostics/erq_sample_summary.json` | 拟合前样本诊断：原纳入数、ID 匹配数、有效 ERQ 匹配数、逐步筛选人数、首个排除原因及 ID 示例 |
| `diagnostics/erq_sample_roster.csv` | 拟合前逐被试诊断，含规范化 ID、匹配标记和 ERQ 源 ID；即使最终 0 人也保存 |
| `diagnostics/erq_rows.csv` | 全部 ERQ 源行、匿名/真实 ID、原始/数值分数及合并和纳入状态 |
| `diagnostics/erq_id_crosswalk.csv` | 源数据行号、匿名 ID、真实 ID 与影像 ID 的对应表 |
| `diagnostics/erq_filter_summary.csv` | 各步骤实际新增排除人数、剩余人数及该条件不合格人数 |
| `diagnostics/erq_exclusions_by_filter.csv` | 每项本次过滤不合格的被试名单，含首个实际排除标记 |
| `diagnostics/original_exclusion_summary.csv` | 按原年龄分析排除原因汇总的人数 |
| `diagnostics/original_exclusions_by_filter.csv` | 原分析已排除的被试，逐人列出原原因 |
| `statistics/unit_response_exclusions.csv` | 各单元 contrast 缺失/非有限排除 ID；三个模型同步排除 |
| `statistics/erq_full_vs_baseline.csv` | 保留的补充联合比较，含三个模型指标和整体 F/p/q |
| `statistics/full_vs_baseline_significant.csv` | 补充联合比较 `q_full_vs_baseline < 0.05` 的单元及配对指标 |
| `statistics/erq_age_statistics.csv` | 全部 449 单元宽表：三模型指标、两项逐步比较、联合/辅助检验、两种 ERQ 系数及兼容列 |
| `statistics/erq_sequential_comparisons.csv` | 全精度三线表数据，每单元三行，合计 1347 行；主效应/交互/整体贡献 |
| `statistics/erq_manuscript_table.docx` | 全部 449 单元的可编辑 Word 三线表，显示舍入，不筛掉非显著或失败单元 |
| `statistics/erq_table_definitions.json` | 字段来源、效应单位、每范围/每比较有效校正数及解释边界 |
| `statistics/additive_vs_baseline_significant.csv` | 第一步 ERQ 主效应显著单元 |
| `statistics/full_vs_additive_significant.csv` | 第二步交互显著单元 |
| `statistics/sequential_significant.csv` | 任一逐步比较通过的单元 |
| `statistics/erq_reference_slopes_at_ages.csv` | 所有有效单元在实际年龄范围内的 6/8/10/12/14/16/18 岁 ERQ 斜率、SE、95%点态 CI；描述性，不作新的校正显著性判定 |
| `statistics/joint_erq_significant.csv` | 与补充联合比较显著列表相同，附完整统计，保留旧名 |
| `statistics/conditional_age_curves.csv` | 显著单元的 ERQ P16/P50/P84 条件年龄曲线和 95% CI |
| `statistics/erq_slopes_by_age.csv` | 显著单元的 ERQ 斜率随龄曲线和 95% CI |
| `statistics/figure_manifest.csv` | 图与单元的对应关系，路径相对于 `figures/` |
| `figures/age_by_erq/` | 同一个连续模型在三个 ERQ 参考值处的条件预测 SVG |
| `figures/erq_slope_by_age/` | `d(activation)/d(ERQ)` 随年龄变化的 SVG，不是年龄一阶导数 |
| `run_metadata.json` | 模型、平滑参数搜索、样本来源、协变量编码、检验和校正定义 |

`--step check` 和完整拟合都会先写 `diagnostics/` 的上述样本审计文件；完整拟合成功后，相同审计还保存到 `inputs/` 作为本次模型输入快照。预检不改写已有模型统计及正式输入快照。单元 contrast 排除名单仅在完整拟合时产生。

完整结果表中的指标分为以下几类（模型后缀均包含 baseline、additive、full）：

| 列名 | 含义 |
|---|---|
| `rss_baseline/additive/full` | 三个独立 GCV 模型在共同样本上的残差平方和 |
| `edf_baseline/additive/full` | 三模型的有效自由度 |
| `r2_baseline/additive/full` | 同一 contrast 总变异下的训练 R² |
| `gcv_baseline/additive/full` | 选定参数处 GCV |
| `lambda_age_baseline`、`lambda_age_full`、`lambda_interaction_full` | 实际选定的惩罚参数；完整表保留 `lambda_age`、`lambda_interaction` 作为 full 的旧列名 |
| `delta_rss_full_vs_baseline` | `RSS_baseline - RSS_full`，正值表示 full 的训练 RSS 更小 |
| `delta_r2_full_vs_baseline` | `R²_full - R²_baseline`，正值表示 full 的训练 R² 更高 |
| `rss_reference_baseline/additive/full` | 三套正式检验参考 RSS |
| `rank_baseline/additive/full_reference` | 参考矩阵秩 |
| `F/df1/df2/p/q_additive_vs_baseline`、`F/df1/df2/p/q_full_vs_additive` | 两项主要比较及各自 BH |
| `F/df1/df2/p/q_full_vs_baseline` | 补充联合比较及其 BH |
| `partial_r2_reference_full_vs_baseline` | 参考空间内增加 ERQ 相关项的 partial R² |

三个 GCV 模型的各步拟合增益保留真实符号、不截零：独立选择惩罚时，更复杂模型的训练 RSS 不保证更小。这些描述指标不是检验 p 值，也不是 ERQ 中介比例。所有步骤均有 `delta_rss_*`、`delta_r2_*`。

`beta/se/ci_low/ci_high/t/df/p_erq_additive_reference` 是 E2 的恒定 ERQ 斜率，单位为 contrast / 原始 ERQ 分；`beta/se/ci_low/ci_high/t/df/p_erq_reference` 是 E3 的年龄平均 ERQ 斜率，后者的分项校正看 `q_erq`。两者不能互相替代。非线性交互没有单个 beta，用整体 F、自由度、p/q、partial R² 报告；三线表交互 beta 栏留空线。

`rss_reference_*`、`partial_r2_reference_*`、`beta_erq_reference` 及其标准误/区间描述用于检验的 OLS 参考拟合，不能与实际惩罚拟合的 RSS 混用。参考 partial R² 为 `max(0, 1 - RSS_reference_full / RSS_reference_baseline)`。完整表中的 `*_joint_erq` 与相应 `*_full_vs_baseline` 比较列保持一致。

条件年龄图使用同一个中心化参考：实际年龄分布上、平均 ERQ 和平均协变量处的拟合预测均值。三个 ERQ 参考值的曲线共同减去这一数值，保留它们之间的 ERQ 主效应差异。95% CI 来自完整 pyGAM 预测协方差，曲线与区间整体平移，保留总体均值不确定性。ERQ 斜率图用模型矩阵对 ERQ 的差分及对应协方差计算。

两类区间都是以选定平滑参数为条件的 95% 点态区间，不是个体预测区间或同时置信带，也未包含显著后选择的不确定性。绘图区宽高比 3:4，横轴明确显示 6、8、10、12、14、16、18 岁，纵轴按曲线及真实 CI 自适应。ERQ 斜率图沿用浅黄色置信带、无描边；SVG 保留可编辑文字。图文件使用 `parcel_0156.svg`、`yeo_Default.svg` 等范围前缀区分，实际名称以 manifest 为准。

没有任一步显著单元时，显著列表、曲线表、图清单仍保留表头；全量 Word/CSV 照常输出。重画/重跑仅清理旧清单拥有的失效显著图，不删除其他用户文件。

## 云端执行

将以下文件按工程中的相同路径同步到云端，不需要打包：

1. `analysis/fig03/gam_pygam/run_erq.py`（新增分析入口）
2. `analysis/fig03/gam_pygam/plot_erq.py`（新增绘图模块）
3. `analysis/fig03/gam_pygam/README_ERQ.md`（本文档）
4. `analysis/fig03/gam_pygam/README.md`（补充新分析入口）
5. `analysis/fig03/gam_pygam/result_tables.py`（全精度比较 CSV 与 Word 三线表）
6. `analysis/fig03/gam_pygam/README_RESULTS.md`（结果读取和解释指南）

入口复用云端已有的最新 `run_pygam.py`、`plot_region_age.py` 和 `gam_nested/amygdala.py`。不修改原 m3 年龄模型入口。下面使用既有 FD≤0.5 的 Tian432 缓存，结果写入新的 ERQ 目录；在工程根目录执行：

可先做样本预检，查看上述筛选名单和 ERQ 逐行整合记录；命令只检查和保存诊断，不拟合模型：

```bash
python -m analysis.fig03.gam_pygam.run_erq \
  --step check \
  --inputs "outputs/fig03/gam_pygam_volume_432_models3_FD05/inputs" \
  --erq "/public/home/dingrui/fmri_analysis/data/beh/er_profiles_new.csv" \
  --erq-id-column sub_ids \
  --erq-column ERQ_CR \
  --erq-id-source "/public/home/dingrui/fmri_analysis/data/beh/er_profiles.csv" \
  --output "outputs/fig03/erq_age_pygam_432_models3_FD05"
```

完整命令自动预检、拟合、校正、导出三线表及绘图。沿用原 conda 环境，新增 Word 依赖，在工程根目录执行：

```bash
python -m pip install python-docx

python -m analysis.fig03.gam_pygam.run_erq \
  --inputs "outputs/fig03/gam_pygam_volume_432_models3_FD05/inputs" \
  --erq "/public/home/dingrui/fmri_analysis/data/beh/er_profiles_new.csv" \
  --erq-id-column sub_ids \
  --erq-column ERQ_CR \
  --erq-id-source "/public/home/dingrui/fmri_analysis/data/beh/er_profiles.csv" \
  --output "outputs/fig03/erq_age_pygam_432_models3_FD05"
```

不传 `--step` 即执行三模型完整流程。缓存不同只调整 `--inputs`；沿用原筛选标记，不重复传 FD 阈值。需要 pyGAM、NumPy、pandas、SciPy、Matplotlib、python-docx 及原年龄入口依赖。可用 `--skip-docx` 先输出全部 CSV，之后独立导出 Word。

如果报错 `ERQ score column 'CR' not found` 且可用列包含 `ERQ_CR`，应将参数改为 `--erq-column ERQ_CR`，同时保留 `--erq-id-column sub_ids`。程序在拟合前会打印实际选用的 ID 和量表列，当前文件应显示 `ID='sub_ids'; score='ERQ_CR'`。不需要重命名原始文件或删除输出目录。

如果报错 `ERQ subject ID overlap is 0`，表示原纳入影像样本没有匹配到任何 ERQ ID。本数据应检查是否传入上述 `--erq-id-source` 及其路径，并查看 `erq_id_crosswalk.csv`；不能用匿名顺序号直接合并影像，也不能降低样本数门槛掩盖问题。若原表文件名不同，将该参数改为实际文件路径。若 ID 匹配数大于 0、完整样本仍不足 12，错误会列出缺失分数或协变量等具体排除原因。预检写入 `diagnostics/`，不改写已有模型统计和正式 `inputs/` 快照；完整拟合也先执行同样的预检。

**旧 ERQ 结果必须完整重跑，不能靠重画/导表补中间模型。** 本次命令推荐新目录保留旧结果；也允许在原 ERQ 目录升级 v1/v2，metadata 变为 `erq_age_smooth_pygam_v3`。原年龄目录只读。旧 schema 的 `--step figures/tables` 会要求重新拟合。

拟合完成后，只调整样式并重画：

```bash
python -m analysis.fig03.gam_pygam.run_erq \
  --step figures \
  --output "outputs/fig03/erq_age_pygam_432_models3_FD05"
```

单独重画读取该 ERQ 目录中的统计和曲线表，不需要 ERQ 原始文件，也不会重新拟合。变更模型、样本或区间计算方法时应重新执行完整命令。

只重新导出三线表，不重新拟合：

```bash
python -m analysis.fig03.gam_pygam.run_erq \
  --step tables \
  --output "outputs/fig03/erq_age_pygam_432_models3_FD05"
```

也可以在仅有 pandas、python-docx 的环境中独立导表：

```bash
python -m analysis.fig03.gam_pygam.result_tables \
  --kind erq \
  --statistics "outputs/fig03/erq_age_pygam_432_models3_FD05/statistics/erq_age_statistics.csv" \
  --output "outputs/fig03/erq_age_pygam_432_models3_FD05/statistics"
```

## 本地验证

`test_erq.py` 使用独立 statsmodels OLS 验证逐步 F/自由度、E2 系数/SE/CI、参考 partial R²、三模型共享样本、比较别名、分范围分比较 BH、拟合指标和兼容升级。端到端验证 449 单元、1347 行全精度表、失效图清理及源文件只读；保留逐行 ID 还原、全部 ERQ 源行及逐被试排除审计。测试使用小平滑网格，不改变检验空间。

```bash
python -m unittest analysis.fig03.gam_pygam.test_erq -v
```

测试额外需要 statsmodels。这些合成数据验证不代表真实样本的显著性结果，真实结果需执行上面的云端命令。

完整的 449 单元模拟验证保存在 `outputs/fig03/inference_mock/`，其中 ERQ 比较长表有 1347 行。可重跑：

```bash
python -m analysis.fig03.gam_pygam.mock_inference
python -m analysis.fig03.gam_pygam.result_tables \
  --kind erq \
  --statistics "outputs/fig03/inference_mock/erq/erq_age_statistics.csv" \
  --output "outputs/fig03/inference_mock/erq" --mock
```

模拟 Word 明确标记为 SIMULATED，不是实际分析结果。已检查三线边框、表头重复、全量行数和替代渲染的可读性。本机没有 LibreOffice，标准 `render_docx.py` 无法执行；替代渲染不能保证与 Word 原生分页一致。生成 Word 本身不需要 LibreOffice。
