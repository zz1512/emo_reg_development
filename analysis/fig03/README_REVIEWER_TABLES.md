# BF10和PLS的审稿汇报表

这两套入口沿用ERQ汇总的方式：读取已经完成的分析，导出简洁汇报表、全精度CSV、可复制的Markdown和Word三线表。不会重新拟合GAM/PLS、重新积分BF、重做bootstrap或改变校正。

真实数据及完整结果在服务器。本地的`outputs/fig03/reviewer_tables_mock/`仅作模拟验证，不能用于论文填写。

## 一、BF表：比较方向和证据数量

入口：`analysis/fig03/gam_pygam/summarize_bayes.py`。

- **432版主表8行**：4范围×BF31/BF32两项比较；BF21补充表4行。
- **232版主表2行**：仅parcel232×BF31/BF32；BF21补充表1行。
- 主表每行报告范围、N/范围、比较、有效单元数、BF≥3数量、1/3<BF<3数量、BF≤1/3数量。
- 全精度摘要另保留BF及log BF的最小/最大值、BF≥10和≤0.1数量、先验方向稳定数量及无效单元数。不能平均各脑区BF来制造一个全脑BF。

| 比较 | BF10分子/分母 | 解读 |
|---|---|---|
| BF31 | M3/M1 | 允许年龄函数的模型相对无年龄模型 |
| BF32 | M3/M2 | 允许非线性相对线性年龄模型 |
| BF21 | M2/M1 | 线性年龄相对无年龄模型，单独补充 |

主先验线性和非线性尺度均为0.5。BF≥3支持分子模型、BF≤1/3支持分母模型、中间为证据不足；3和10仅为原工程的描述性约定。BF32≤1/3支持线性模型相对M3，不能写成没有年龄关联。BF不做BH，也不能把这些计数称为FDR显著数。

脚本核对BF和log BF、三种BF的恒等关系、主结果积分精度及保存的长表。存在九套先验明细时，还核对完整组合、积分精度、主先验及敏感性摘要。失败单元保留在明细和无效清单，不以BF=0替代失败；数值溢出优先读log BF。敏感性范围不是置信区间。

### 服务器命令

在服务器工程根目录执行；如果结果目录名不同，调整`--analysis-dir`：

```bash
python -m analysis.fig03.gam_pygam.summarize_bayes \
  --analysis-dir "outputs/fig03/gam_pygam_volume_432_models3_FD05/bayes"

python -m analysis.fig03.gam_pygam.summarize_bayes \
  --analysis-dir "outputs/fig03/gam_pygam_volume_232_models3_FD05/bayes"
```

输入是包含`age_bayes_statistics.csv`的**bayes目录**。默认输出到该目录下`reviewer_summary/`。

| 输出 | 432版数据行数 | 232版数据行数 |
|---|---:|---:|
| `bf_main_summary.csv` | 8 | 2 |
| `bf_linear_summary.csv` | 4 | 1 |
| `bf_all_summary.csv` | 12 | 3 |
| `bf_primary_details.csv` | 898 | 464 |
| `bf_linear_details.csv` | 449 | 232 |
| `bf_invalid_tests.csv` | 按实际无效比较数；每无效单元3行 | 同左 |

明细保留`dimension/name/label`、`n_obs/status`、`comparison`、`bf10/bf01/log_bf10/evidence`、先验敏感性范围与标记，以及源表中已有的`frequentist_*`。原频率学统计未匹配时保持空值，不用BF反推p、q或效应方向。

例如`n_support_first`表示BF≥3的单元数，`n_support_second`表示BF≤1/3的单元数；对应模型由本行`comparison`决定。`n_strong_first/second`已包含在前述两组数量中，不能重复相加。`n_prior_direction_stable`只说明先验网格中BF始终位于1的同一侧，不保证始终达到3或1/3。

### 只汇报关注的6个脑区

若论文只需汇报此前绘图的152、154、182、185，以及另外关注的156、157，直接从完整432版结果摘选，不再使用全脑计数总结：

```bash
python -m analysis.fig03.gam_pygam.summarize_bayes \
  --analysis-dir "outputs/fig03/gam_pygam_volume_432_models3_FD05/bayes" \
  --parcels 152 154 182 185 156 157
```

本次仅需替换工程中的`analysis/fig03/gam_pygam/summarize_bayes.py`；共同依赖`analysis/fig03/reporting.py`和`analysis/fig03/gam_pygam/result_tables.py`沿用此前已同步版本。命令只读取完整结果，不重新拟合、积分或校正。

对应Schaefer400编号：152=LH Default Temp 4，154=LH Default Temp 6，182=LH Default PFC 17，185=LH Default PFC 20，156=LH Default Temp 8，157=LH Default Temp 9。232图谱的同号分区不等价，该选区入口会阻止混用。

默认输出到`bayes/reviewer_selected_parcels/`，保留原`reviewer_summary/`全量报告：

- `selected_region_summary.csv`：6行，每脑区一行，包含N、状态及三种模型比较的BF、log BF和先验敏感性范围。
- `selected_region_comparisons.csv`：18行，六区×三比较，保留源全精度字段与原频率学统计，并增加直观比较名称。
- `selected_invalid_tests.csv`：不可计算的比较；不将失败写成BF=0。
- `response_tables.docx`：一张6行Word三线表；三个BF列直接标为`Age spline vs. no age`、`Age spline vs. linear age`、`Linear age vs. no age`，不使用BF31/32/21作展示标题。
- `result_summary.md`：同一张可复制表格和按这6区真实BF生成的英文结果段落，逐比较列出对应脑区与证据方向。
- `summary_metadata.json`：所选编号、来源哈希及核验状态。

这6区是依据此前区域分析与绘图选取的关注区，应作为后续证据描述，不称为独立确认；尤其不能仅因被选中就预设BF支持年龄效应，也不将156/157自动称为显著或边缘显著。实际支持方向由每区BF决定，旧全脑计数无法反推出具体6区数值。

## 二、PLS表：预测性能和增量

入口：`analysis/fig03/predictive_pls/summarize_pls.py`。

| 表 | 数据行数 | 填写字段 |
|---|---:|---|
| `pls_main_summary.csv` | 2 | success和age各1行；N、站点数、四模型R²、M2−M1的ΔR²及95%CI、有效bootstrap数、报告状态、实际PFC定义 |
| `pls_model_performance.csv` | 8 | 两任务×四模型；特征数、成分数摘要、R²及95%CI、MAE、RMSE、Pearson r、相对M0误差改善 |
| `pls_all_comparisons.csv` | 10 | 两任务×五项预设比较；ΔR²及95%CI、ΔMAE、ΔRMSE、主/辅助角色、原推断标记及新增报告状态 |

**优先填2行主表，再配8行模型表现表。** 完整10行比较供补充材料。M2包含432分区；本PLS没有另把7网络、4系统和6杏仁核作为额外模型特征。不能套用ERQ的四范围汇总结构。

主比较：`task=success`且`comparison=M2_distributed432-M1_PFC`；年龄同一比较是次要分析。ΔR²为前一模型减后一模型，正数表示分布式模型R²更高；ΔMAE/ΔRMSE负数表示前一模型误差更小。负R²保留，不能改成0。R²来自拼接外层测试预测，不是训练R²或各站点R²的简单平均。

本脚本核对模型表现、比较差值和外层预测是否一致；读取`frozen/pfc_parcel_membership.csv`识别实际PFC集合，不依据文件夹名称猜测。Buhle元分析PFC与Donahue宽范围解剖PFC必须分别运行、分别解释。

### bootstrap状态必须同时报告

| `report_status` | 含义 |
|---|---|
| `bootstrap_not_run` | 只有外层点估计，尚无可用bootstrap推断 |
| `bootstrap_completion_unverified` | 有保存区间和状态，但缺少`bootstrap_metrics.csv`，无法核实完整重复记录 |
| `bootstrap_incomplete` | 请求的重复编号尚未全部完成；区间只能作暂定描述 |
| `bootstrap_quality_insufficient` | 重复记录完整，但未达到原分析的有效重复数/比例要求 |
| `inference_available` | 重复完整且满足原分析质量要求，可按预设角色解释 |

原判据为请求至少2000次、有效至少1000次且比例至少95%；只有success的M2−M1、95%CI下界>0才允许`支持主比较`。本脚本还检查保存的重复记录是否完整，防止提前合并的部分重复被当作最终推断。

`supported`和`inference_valid`原值保留；新增`confirmation_applicable`说明是否属于确认性主比较，`supported_for_report`同时要求重复完整。其他比较和年龄即使区间为正，`supported=False`也不等于“没有差异”，因为它们并非确认性主假设。原标记、区间、重复计数相互矛盾时明确报错，不修改源统计文件。

未运行或未完成bootstrap时，不能把False解释成无增量。有效推断未支持主比较时，也不能写成PFC与全脑等效。若M2优于M1，还应结合M2−M0与四模型绝对性能判断是否超过协变量。年龄预测改善不能单独证明神经活动解释了行为发展或网络重组。

### 服务器命令

```bash
python -m analysis.fig03.predictive_pls.summarize_pls \
  --analysis-dir "outputs/fig03/predictive_pls_432_FD05"
```

如果要汇总元分析PFC版本，单独运行：

```bash
python -m analysis.fig03.predictive_pls.summarize_pls \
  --analysis-dir "outputs/fig03/predictive_pls_432_FD05_meta_pfc"
```

输入是包含`model_performance.csv`、`model_comparisons.csv`、`oof_predictions.csv`和`frozen/`的预测结果根目录。bootstrap已经合并时还会读取`bootstrap_metrics.csv`及`frozen/bootstrap_status.json`，核对百分位区间。汇总不运行bootstrap，不从区间推算p、q或BF。

默认输出到各自结果目录下`reviewer_summary/`。不同定义不会合并；分数保持原方向，在评分公式未核实前不自动写成“能力更强”。

## 三、共同输出与云端同步文件

两套脚本各自输出：

- `response_tables.docx`：英文Word三线表。
- `result_summary.md`：可复制的表格、中文数值总结及解释说明。
- `summary_metadata.json`：输入SHA256、输出行数、模拟标记和来源核验状态。
- 上述全精度CSV：17位有效数字；Word和Markdown只作显示舍入。

共同参数：`--output`指定独立汇总目录；`--skip-docx`只导出CSV/Markdown；`--mock`明确标记模拟验证。原结果文件保持只读，输出目录不覆盖原分析目录，只更新该脚本自己创建的报告目录。

**服务器需同步以下工程文件，不需要压缩包：**

1. 新增`analysis/fig03/reporting.py`：共同导出与Word排版工具。
2. 新增`analysis/fig03/gam_pygam/summarize_bayes.py`。
3. 新增`analysis/fig03/predictive_pls/summarize_pls.py`。
4. 确认已有最新版`analysis/fig03/gam_pygam/result_tables.py`；若服务器版本较旧，一并替换。
5. 阅读说明为本文件；拟合、bootstrap入口无需替换或重跑。

依赖仅NumPy、pandas、python-docx，不要求导入pyGAM或scikit-learn。原分析环境一般已具备；缺少Word依赖可执行`python -m pip install python-docx`或使用`--skip-docx`。

## 四、本地验证

```bash
python -m unittest analysis.fig03.test_reviewer_summaries -v
```

测试文件同时提供确定性的模拟结果结构，测试不重新拟合分析。另使用已有449单元、4041行先验敏感性明细的BF拟合模拟结果核对兼容性。所有模拟数字不能用于真实研究结果。

Word经过三线表结构、表头和行数检查。本机运行环境未提供LibreOffice，标准DOCX渲染器无法核验原生分页；正式使用前仍需在Word/LibreOffice查看版式。该限制不影响CSV提取或Word导出。
