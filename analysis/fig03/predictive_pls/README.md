# 跨站点 PLS：宽范围解剖 PFC

结果目录：`outputs/fig03/predictive_pls_432_FD05`。设计全文见 [实现方案](PLAN_MAIN.md)，PFC 文献核查见 [PFC 定义核查](README_PFC_LITERATURE.md)。本分析不改 GAM、ERQ 或贝叶斯因子。

主分析只用 432 分区、meanFD ≤ 0.5。四个模型是协变量、协变量加 PFC 分区、协变量加全部 432 分区、协变量加非 PFC 分区。外层和内层都是留一站点。确认性终点是重评成功的 M2−M1 外层 R² 增量。年龄用同一套模型，作为次要分析。

所有命令都在含 `analysis/` 的工程根目录执行：

```bash
cd /public/home/dingrui/fmri_analysis/zz_analysis
conda activate brain
```

不要在 `outputs/fig03/predictive_pls_432_FD05` 里运行，也不要使用 `(base)`。相对路径 `outputs/fig03/...` 只在工程根目录下有效。

## 1. 服务器需要的文件

预测代码：

- `analysis/fig03/predictive_pls/run_pls.py`
- `analysis/fig03/predictive_pls/models.py`
- `analysis/fig03/predictive_pls/bootstrap.py`
- `analysis/fig03/predictive_pls/report_results.py`
- `analysis/fig03/predictive_pls/precheck.py`
- `analysis/fig03/predictive_pls/config_main.json`
- `analysis/fig03/predictive_pls/build_buhle_pfc.py`：只在本地重做元分析 PFC 成员表时需要；服务器拟合不运行它

成员表只上传这一个 CSV，放到与本地相同的相对路径：

`outputs/fig03/predictive_pls_432_FD05/pfc_definition/pfc_parcel_membership.csv`

`pfc_definition.json` 和 `pfc_definition_preview.svg` 不参与拟合，可以不上传。服务器不要运行 `make_pfc_definition`。`build_local_pfc.py` 和 `resources/atlases/donahue2018/` 只在本地重做成员表时使用。

输入沿用已经完成的 432、FD≤0.5 年龄分析缓存：

`outputs/fig03/gam_pygam_volume_432_models3_FD05/inputs`

其中需要 `design.csv`、`parcel_response.csv`、`parcel_info.csv`、`tian_combined.json`。行为分数读 `er_profiles_new.csv` 的 `reappraisal_success`，真实编号用同目录 `er_profiles.csv` 每一对应行的前 6 位。不使用按 ERQ 完整性筛过的样本，不重读 NIfTI。

## 2. PFC 成员表

liberal PFC 取 Donahue 等（2018，PNAS，115:E5183–E5192）补充表 S2 的全部 HCP-MMP 分区。顶点标签是公开的组平均 HCP_MMP1.0，面积是 TemplateFlow fsLR 32k midthickness 三角形面积的三分之一。这不是需要登录的 BALSA 场景文件。

皮层分区覆盖率 ≥ 0.50 归 PFC，Tian 401–432 覆盖率为 0，归非 PFC。本地核对结果：

| 项目 | 结果 |
| --- | --- |
| PFC 分区 | 80，左 40、右 40 |
| 跨界分区 `boundary` | 47，覆盖率在 0 和 1 之间 |
| 重点核查 `review_band` | 11，覆盖率在 0.20 到 0.80 |
| 名称含 PFC 或 OFC | 83 个里 80 个纳入 |
| 名称含 PFC 或 OFC 但未纳入 | 355 `RH_Cont_PFCl_15` 0.074；340 `RH_Cont_PFCv_1` 0.068；143 `LH_Cont_PFCv_1` 0.015 |
| FEF | 未纳入，最高覆盖率 0.223 |
| 视觉分区、Tian | 覆盖率 0 |

左右掩膜分别有 6139 和 5970 个顶点，只与同侧 Schaefer 分区重叠。

`pfc_definition.json` 记录来源、规则、上述计数和源文件 SHA256。哈希用于确认文件没有被替换，不能代替空间核对。`pfc_definition_preview.svg` 是覆盖率直方图，图中写明没有 inflated 表面，所以它不是左右脑位置图。

## 2b. 元分析 PFC：回应审稿人的方案

这一版回答审稿人举的例子：用交叉验证比较前额叶模型和分布式模型对重评成功、年龄的预测。它不覆盖第 2 节的 Donahue 宽范围 80 区，结果写到单独目录 `outputs/fig03/predictive_pls_432_FD05_meta_pfc`。

四个模型、两个因变量、留一站点和 PLS 规则与第 4 节相同。变化只在 `M1_PFC` 使用哪些分区：

| 模型 | 这一版的脑区 |
| --- | --- |
| `M0_covariates` | 无 |
| `M1_PFC` | 下面 10 个元分析 PFC 分区，各自保留 |
| `M2_distributed432` | 全部 432 分区 |
| `M3_nonPFC` | 其余 422 个分区 |

主比较仍是重评成功的 `M2_distributed432-M1_PFC`。`delta_r2` 是外层 R² 的差。大于 0 表示全部 432 区的点估计更高。站点重抽样区间是补充，不是审稿人要求的标准。

PFC 范围来自 Buhle、Silvers、Wager 等（2014，Cerebral Cortex，24:2981–2990）认知重评元分析的公开共识 z 图，文件是 `resources/atlases/buhle2014/Buhle_Silvers_2014_Emotion_Regulation_Meta.nii.gz`。Neurosynth 的 `reappraisal` 页面不再提供可下载 nii，所以不用那个自动词条图。

对应规则：

1. 保留 z ≥ 2.3 的正簇。这个阈值对上原文 Table 2 的前额叶、顶叶和颞叶峰。
2. 只保留 Harvard-Oxford 的额极、额上回、额中回、额下回、旁扣带和前扣带。顶叶、颞叶、岛叶，以及原文未作为重评执行网络的眶额和腹内侧皮层，不进入 PFC。
3. 把掩膜投到 Schaefer 400 的 2 mm 体积网格。分区体素有至少 25% 落在掩膜内才纳入。50% 只会留下 3 个分区，右外侧前额叶的峰会被丢掉。
4. Tian 401–432 覆盖率记为 0，不是 PFC。杏仁核会出现在重评调制结果里，但不算这次的 PFC 范围。

本地重做成员表：

```bash
python -m analysis.fig03.predictive_pls.build_buhle_pfc
```

纳入的 10 个分区：

| 编号 | 名称 | 覆盖率 |
| --- | --- | --- |
| 107 | 7Networks_LH_SalVentAttn_Med_1 | 0.571 |
| 148 | 7Networks_LH_Cont_PFCmp_1 | 0.412 |
| 175 | 7Networks_LH_Default_PFC_10 | 0.655 |
| 182 | 7Networks_LH_Default_PFC_17 | 0.278 |
| 185 | 7Networks_LH_Default_PFC_20 | 0.563 |
| 189 | 7Networks_LH_Default_PFC_24 | 0.448 |
| 352 | 7Networks_RH_Cont_PFCl_12 | 0.319 |
| 353 | 7Networks_RH_Cont_PFCl_13 | 0.253 |
| 361 | 7Networks_RH_Cont_PFCmp_2 | 0.274 |
| 378 | 7Networks_RH_Default_PFCv_4 | 0.255 |

上传 `outputs/fig03/predictive_pls_432_FD05_meta_pfc/pfc_definition/pfc_parcel_membership.csv`，并同步新的 `run_pls.py`。旧脚本只接受 Donahue 成员表。

在工程根目录、`brain` 环境下执行。外层预测：

```bash
python -m analysis.fig03.predictive_pls.run_pls \
  --inputs outputs/fig03/gam_pygam_volume_432_models3_FD05/inputs \
  --behavior /public/home/dingrui/fmri_analysis/data/beh/er_profiles_new.csv \
  --behavior-id-column sub_ids \
  --score-column reappraisal_success \
  --behavior-id-source /public/home/dingrui/fmri_analysis/data/beh/er_profiles.csv \
  --membership outputs/fig03/predictive_pls_432_FD05_meta_pfc/pfc_definition/pfc_parcel_membership.csv \
  --output outputs/fig03/predictive_pls_432_FD05_meta_pfc
```

站点 bootstrap。跑完 0–1999 会自动合并，不必再单独加 `--merge-only`：

```bash
python -m analysis.fig03.predictive_pls.bootstrap \
  --analysis outputs/fig03/predictive_pls_432_FD05_meta_pfc \
  --start 0 --count 2000 \
  --jobs 8
```

图表：

```bash
python -m analysis.fig03.predictive_pls.report_results \
  --analysis outputs/fig03/predictive_pls_432_FD05_meta_pfc
```

结果字段与第 5–8 节相同，路径换成 `outputs/fig03/predictive_pls_432_FD05_meta_pfc/`。回复审稿人时先看 `model_comparisons.csv` 里 `task=success`、`comparison=M2_distributed432-M1_PFC` 的 `delta_r2`，再看 `task=age` 的同一比较，并从 `model_performance.csv` 读取 `M0_covariates`、`M1_PFC`、`M2_distributed432` 的 `r2`。

## 3. 服务器命令

先跑外层预测。这一步内部会写 `precheck/`，再拟合四个模型。文件名 `pfc_parcel_membership.csv` 必须保持在同一行。

```bash
python -m analysis.fig03.predictive_pls.run_pls \
  --inputs outputs/fig03/gam_pygam_volume_432_models3_FD05/inputs \
  --behavior /public/home/dingrui/fmri_analysis/data/beh/er_profiles_new.csv \
  --behavior-id-column sub_ids \
  --score-column reappraisal_success \
  --behavior-id-source /public/home/dingrui/fmri_analysis/data/beh/er_profiles.csv \
  --membership outputs/fig03/predictive_pls_432_FD05/pfc_definition/pfc_parcel_membership.csv \
  --output outputs/fig03/predictive_pls_432_FD05
```

此时已有点估计，还没有区间。只想先看样本、不拟合时，可以把同一套输入交给预检，输出改成 `outputs/fig03/predictive_pls_432_FD05/precheck`：

```bash
python -m analysis.fig03.predictive_pls.precheck \
  --inputs outputs/fig03/gam_pygam_volume_432_models3_FD05/inputs \
  --behavior /public/home/dingrui/fmri_analysis/data/beh/er_profiles_new.csv \
  --behavior-id-column sub_ids \
  --score-column reappraisal_success \
  --behavior-id-source /public/home/dingrui/fmri_analysis/data/beh/er_profiles.csv \
  --output outputs/fig03/predictive_pls_432_FD05/precheck
```

站点 bootstrap 默认 2000 次，种子 20261111。每次完整重拟合。它读取上一步已经写好的 `frozen/*_arrays.npz`，不必重跑 `run_pls`。每个训练折只拟合一次 PLS，再取出各个成分数的系数。`--jobs` 是并行进程数，默认用当前机器的全部 CPU；每个进程只用一个计算线程。一次跑完：

```bash
python -m analysis.fig03.predictive_pls.bootstrap \
  --analysis outputs/fig03/predictive_pls_432_FD05 \
  --start 0 --count 2000 \
  --jobs 8
```

把 `--jobs` 改成服务器上打算占用的核数。进程会打印每完成 20 个任务的进度。时间长时仍可按 100 次分块，`--start` 依次取 0、100、200，直到 1900。分块编号不重叠，抽样由种子和重复编号决定。0–1999 到齐后自动合并。已有全部分块、只补合并时加 `--merge-only`。

```bash
python -m analysis.fig03.predictive_pls.bootstrap \
  --analysis outputs/fig03/predictive_pls_432_FD05 \
  --merge-only
```

最后出表和图：

```bash
python -m analysis.fig03.predictive_pls.report_results \
  --analysis outputs/fig03/predictive_pls_432_FD05
```

有效 bootstrap 少于 1000 次，或有效比例低于 95% 时，区间仍可阅读，`supported` 保持为假。无效重复保留原因，不改记成 0，也不为了凑满次数而重新抽样。

## 4. 样本和模型

两个任务分开建样本。年龄保留原 `in_analysis=1`，并且年龄、性别、父母教育、meanFD 有效，性别为 0/1，meanFD ≤ 0.5，站点非空，432 列全部有限。重评成功在此之上还要求行为编号匹配、`reappraisal_success` 有效。年龄任务不因行为或 ERQ 缺失删人。两个任务的人数可以不同。同一任务内四个模型使用同一批被试和同一套站点划分。

重评成功的协变量是截距、年龄的 4 列三次样条、性别、父母教育和 meanFD。年龄预测的协变量是截距、性别、父母教育和 meanFD，不再把年龄放进协变量。样条结点、连续协变量缩放、近零方差列和 PLS 成分数都在每一层训练折内拟合。成分数从 0 到 10 中选择，内层均方误差并列时取更小的 k。k = 0 表示只用协变量。

| 模型 | 神经特征 |
| --- | --- |
| `M0_covariates` | 无 |
| `M1_PFC` | `is_pfc` 为真的分区，各自保留，不先合成一个 PFC 均值 |
| `M2_distributed432` | 全部 432 分区 |
| `M3_nonPFC` | 其余分区 |

特征不含小脑，也不另加 7 网络、4 系统或 6 杏仁核聚合值。站点只用于留站点分组。缓存里没有家系列；同一被试出现在两个站点会停止。

`reappraisal_success` 保持原分数和原方向。公式、所用试次和“越高是否表示越好”尚未核实，正文在核实前只称原分数。

## 5. 结果阅读顺序

路径都相对于 `outputs/fig03/predictive_pls_432_FD05/`。精确数字以 CSV 为准。`tables_main.docx` 使用同一套 `.12g` 格式。跑完第 3 节最后一条命令后，输出目录里还会有一份较短的 `README_RESULTS.md`；字段解释以本文件为准。

1. `frozen/config_snapshot.json`：输入哈希、PFC 区数、实际 scikit-learn 版本、分数方向尚未改写。
2. `precheck/`：两个任务的人数、站点和逐项排除。
3. `pfc_definition/pfc_parcel_membership.csv`：432 区的覆盖率和 PFC 标记。拟合时使用的副本在 `frozen/pfc_parcel_membership.csv`。
4. `oof_predictions.csv`：每名被试、每个任务一次外层测试预测。
5. `model_performance.csv`、`model_comparisons.csv`：四个模型的点估计；bootstrap 合并后才有区间。
6. `site_performance.csv`：各站点。目标恒定或人数少于 2 时，该站点 R² 为 NA，MAE 和 RMSE 仍保留。
7. `inner_cv_scores.csv`、`fold_models.csv`：内层候选 MSE 和选定成分数。
8. `bootstrap_metrics.csv` 与 `frozen/bootstrap_status.json`：每次重复是否有效。
9. `figures/` 与 `tables_main.docx`。

预检里的 `pfc_named_candidates.csv`、`pfc_plus_ofc_candidates.csv` 和 `pfc_boundary_review.csv` 只是名称候选。本次模型用的是成员表里的 `is_pfc`，不是这三份名称表。

## 6. 预检文件

这些文件在 `precheck/`。

| 文件 | 内容 |
| --- | --- |
| `REPORT.md`、`precheck_summary.json` | 两任务人数、站点数、可行性提示、输入哈希 |
| `site_summary.csv` | 各站点人数、年龄范围、行为分布、头动、性别计数 |
| `loso_feasibility.csv` | 外层测试人数、训练站点数、内层是否可做、年龄是否超出训练范围、测试目标是否恒定 |
| `site_age_bins.csv` | 各站点年龄分箱人数 |
| `filter_summary.csv` | 每项筛选前后人数；`n_newly_excluded` 可以相加 |
| `exclusions_by_filter.csv` | 每项条件失败的被试，以及是否为首次排除 |
| `sample_roster.csv` | 原排除原因、行为匹配、两个任务的新排除原因 |
| `behavior_rows.csv` | 每一行为源行、匿名编号、真实编号、原分数、是否进入两个任务 |
| `age_sample.csv`、`success_sample.csv` | 两个任务的有效名单和协变量 |
| `parcel_quality.csv` | 432 列的非有限值个数；只作诊断，不据此删特征 |

同一被试可以同时不满足多项条件，所以各项失败人数不能直接相加。不足 3 个站点时不能做嵌套留一站点。这些提示不自动删站点，这里也不计算预测性能。

## 7. 预测结果字段

`oof_predictions.csv` 每个任务每名被试一行。

| 字段 | 含义 |
| --- | --- |
| `task` | `success` 或 `age` |
| `sub_id`、`site_id`、`outer_site` | 被试和被留出的站点；后两列相同 |
| `y_true` | 原尺度的重评成功或年龄 |
| `pred_<model>`、`error_<model>` | 该模型的外层测试预测和 `y_true` 减预测 |
| `k_<model>` | 该外层折选定的成分数；M0 为 0 |

`model_performance.csv` 每个任务每个模型一行。

| 字段 | 含义 |
| --- | --- |
| `n`、`n_sites`、`n_features` | 样本量、站点数、定义中的特征数；M0 为 0 |
| `median_k` | 各被试所在外层折选定成分数的中位数 |
| `r2`、`mae`、`rmse`、`pearson_r` | 全部外层测试拼在一起后的指标 |
| `prediction_error_improvement_vs_m0` | `1 - SSE_M / SSE_M0`，是误差改善，不是因果解释比例 |
| `ci_low_r2`、`ci_high_r2` | 合并 bootstrap 后，该模型 R² 的 2.5% 和 97.5% 分位数 |
| `n_bootstrap_valid`、`bootstrap_fraction_valid` | 该任务的有效重复次数和比例 |
| `confirmatory` | 此列保持为假；是否支持看比较表的 `supported` |

`model_comparisons.csv` 每个任务五行。主行是 `task=success` 且 `comparison=M2_distributed432-M1_PFC`，`role` 为 `primary`。其余比较和年龄任务的 `role` 为 `support`。

| 比较 | 问题 |
| --- | --- |
| `M2_distributed432-M1_PFC` | 全部 432 区是否高于 PFC 区，主比较 |
| `M1_PFC-M0_covariates` | PFC 区是否高于协变量 |
| `M2_distributed432-M0_covariates` | 全部 432 区是否高于协变量 |
| `M3_nonPFC-M0_covariates` | 非 PFC 区是否高于协变量 |
| `M2_distributed432-M3_nonPFC` | 全部 432 区是否高于非 PFC 区 |

`delta_r2` 大于 0 表示前一个模型的外层 R² 更高。`delta_mae` 和 `delta_rmse` 用同一顺序相减，负数表示前一个模型的绝对误差更小。R² 可以是负数，表示外层预测比用测试集总均值更差，保留负数。Pearson r 只描述预测和实测的线性关系，不用它的普通相关 p 值代替模型比较。

`inference_valid` 为真需要请求满 2000 次、有效不少于 1000 次、有效比例不少于 95%。`supported` 只在重评成功的 M2−M1 上，并且 `inference_valid` 为真、区间下界大于 0 时才为真。年龄任务和其他比较可以有区间，不作为确认性发现。

`site_performance.csv` 是每个任务、每个站点、每个模型一行，并附有该站点的 `delta_r2_M2_minus_M1`、`delta_mae_M2_minus_M1`、`delta_rmse_M2_minus_M1`。`y_constant` 为真时该站点 R² 和 Pearson r 为 NA。

`inner_cv_scores.csv` 的 `pooled_mse` 是某个外层训练折里、某个模型和某个 k 的内层留站点均方误差。`fold_models.csv` 记录最终选定的 `k`、训练人数、训练站点数、协变量秩、特征秩、用到的特征数、被丢掉的近零方差分区、共同候选上界 `common_k_max`，以及选定 k 是否在迭代上限内收敛。

`bootstrap_metrics.csv` 每个任务每次重复一行。`valid` 为假时读 `reason`。`site_counts` 记录每个原始站点被抽到的次数；重复抽到的站点仍使用原来的站点编号，整组进入训练或测试。

## 8. 图和 Word

| 文件 | 内容 |
| --- | --- |
| `figures/model_r2.svg` | 两个任务、四个模型的外层 R²；已有区间时带分位数误差线 |
| `figures/site_delta_success.svg`、`site_delta_age.svg` | 各站点 M2−M1 的 R² 差 |
| `figures/oof_scatter_success.svg`、`oof_scatter_age.svg` | 实测值对 M2 的外层测试预测，按站点着色 |
| `tables_main.docx` | 表现表和比较表的三线表 |

SVG 文字可编辑。散点是外层留出站点的预测，不是全样本拟合值。图不隐藏表现差的站点。

## 9. 报告边界

- 未发现 M2−M1 增量时，写成这一套预定义分区和外层验证下没有支持增量，不写成 PFC 与全脑等效。
- 年龄预测变好，只说明这些脑活动带有可泛化的年龄信息，不写成脑活动解释了重评能力的发育，也不写成网络重组。
- M2 高于一个很差的 M1、但没有高于 M0 时，不写成神经指标提供了可靠的行为增量。
- 跨界分区已经有的均值不能拆成 PFC 与非 PFC 体素。结论写成超出预定义 PFC 分区集合的增量预测信息。
- 重评成功在方向核实前只称原分数。

## 10. 本地检查

```bash
python -m unittest analysis.fig03.predictive_pls.test_precheck analysis.fig03.predictive_pls.test_pipeline -v
```

测试使用合成数据，检查编号还原、排除追溯、训练折不读取测试目标、成分数并列时取较小值、无信号与仅有 PFC 信号、非 PFC 信号，以及 bootstrap 分组。这些测试不代替服务器上的真实样本结果。
