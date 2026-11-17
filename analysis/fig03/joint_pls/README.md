# 年龄＋重评成功：两套 PLSRegression 模型比较

入口：`python -m analysis.fig03.joint_pls.run`。两套分析在一次运行中完成，旧版单目标 `predictive_pls` 代码和输出保留。结果字段见 [结果读取说明](README_RESULTS.md)。

## 两套方案的具体含义

共同输入为 432 个 `Reappraise − LookNegative` 脑区均值，Y 为同一批人的 `[age, reappraisal_success]`。`reappraisal_success` 直接取行为表，不把脑活动 contrast 当作行为成功，不自行反转分数方向。

| 方案 | 学习与评价 | 主要终点 |
|---|---|---|
| `joint` 联合预测 | 每个脑区集合独立学习一个多输出 PLSRegression，同时预测两列 Y；不拆成两个单目标分析 | 两个结局的交叉验证 R² 等权平均 |
| `shared` 共有成分预测 | 每折训练集的全脑参考 PLS 提取第一 Y 方向，四个模型预测同一行为分数 | 同一参考分数的交叉验证 R² |

`joint` 的均值 R² 衡量联合预测，不能单凭此值声称提取了年龄和成功的纯共有方差。

`shared` 中，先在训练集用全脑 X 与标准化 Y 拟合参考 PLS2，取归一化第一 Y 权重 q，定义 `u = Y_z @ q / SD_train(Y_z @ q)`。q 的年龄系数定向为正。四个候选模型都先预测二维 `Y_z`，再以同一个 q 投影预测值 `u_hat = Y_z_hat @ q / SD_train`。因此共有成分方案仍由多输出 PLSRegression 实现，没有再拟合独立的 PCA 或单目标 PLS1。

同一折四个模型的真实 u 完全相同，不能各自提取 Y 分数后直接比较其 R²。q 在每个内层训练折、外层训练折、置换和 bootstrap 中重新学习；测试集 Y 仅用于计分。不同外层折的 q 可以不同，最终指标评价的是这套训练和预测流程，不是一个预先冻结的全体被试心理测量分数。

参考方向使用了训练集全脑信息，是一个**脑信息参与定义的行为方向**；参考选择可能影响模型比较的结论。它不保证两个行为载荷相等、同向或都很大。必须检查 `reference_axes.csv`、`model_first_components.csv`，不能预先命名为共同发展因子。

## 样本、协变量与交叉验证

- 只做 432 分区、meanFD≤0.5。通过预检的 `--membership` 选择已验证定义：`Donahue2018_liberal_bilateral`（80 区解剖 PFC）或 `Buhle2014_reappraisal_PFC`（10 区元分析 PFC），不根据本次显著结果挑选脑区。下面服务器命令使用当前要求的 **meta 版本**。
- 年龄、成功、性别、meanFD、站点和全部 432 脑区完整的共同样本，两套方案和所有模型完全一致。
- **父母教育既不是协变量，也不是纳入条件。** 推荐从原始人口学表＋有效名单＋二级排除名单重建样本，避免沿用旧样本隐含的教育筛选。
- 性别和 meanFD 为预测协变量。仅在训练数据拟合 `X~C`、`Y~C`，再对残差标准化；应用到验证和测试数据。年龄属于 Y，不作为协变量。Y 的训练残差各自单位方差，避免年龄单位主导拟合。
- 站点用于内外层 LOSO 和限制重抽样。没有给新站点拟合截距或用测试站点数据做回归调整。有效站点指每个被试有真实可对齐的站点标签且有至少三个站点；不能用行号构造站点。
- 外层留一站点测试，内层留一训练站点选成分数。候选 k 为 0–10（样本量/特征数不足时自动缩减为所有模型共同可行范围），k=0 允许神经模型回退到协变量基线。并列优先较少成分。内层最大化与该方案一致的 pooled R²。

四模型：M0 协变量；M1 协变量＋PFC；M2 协变量＋全部 432；M3 协变量＋剩余脑区。每套方案按 **M1−M0、M2−M0、M3−M0、M2−M1、M3−M1** 五项顺序输出，joint 和 shared 各自单独绘图。M2−M1 检验合并剩余脑区后相对 PFC 的增益；M3−M1 比较两个不同特征集合的预测表现，本身不等于“控制 PFC 后的增益”。两套方案的 M2−M1 仍为原主检验族，两个检验一起 BH；两套 M3−M1 为次要比较，单独两个一起 BH。三种脑模型各自相对 M0 的比较仍为探索性，新增以“三项比较×两套方案”的六个 P 一起做 BH。三族大小分别为2、2、6，不是全部10项推断比较的统一校正。M3−M2 / M2−M3 不是分析目标，不进入比较表和图。

Meta 版本的 M1 为 10 个分区，M3 为剩余 422 个分区；M3 是元分析 PFC 集合的补集，其中可能包含未纳入元分析掩膜的解剖 PFC。因此报告中用 `Meta-PFC` 和 `Remaining regions` 区分，内部模型 ID `M1_PFC`、`M3_nonPFC` 保持一致。预检冻结成员表到 `precheck/pfc_parcel_membership.csv`，并在 `sample_metadata.json` 和报告中记录定义、规则及脑区数。

已有 80 区结果不能只改输出标签或只重跑 report 变成 meta 结果；需要用 meta 成员表重新预检，在新目录重跑 fit、置换和 bootstrap。两套联合方案仍会一起输出。

## PLS 收敛设置与旧批次失败

所有观察拟合、内外层验证、置换和 bootstrap 统一使用 `PLSRegression(scale=False, max_iter=10000, tol=1e-6)`，设置冻结在 `analysis_manifest.json` 的 `solver` 中。原脚本采用默认 500 次上限，在服务器 80 区批次出现置换 989/1000、bootstrap 994/1000 有效；失败原因为 `ConvergenceWarning: Maximum number of iterations reached`。

新版提高迭代上限，保持原收敛容差。上限是最多允许次数，达到容差即提前结束，不是每次都跑 10000 步。仍未收敛会保留失败记录，附训练样本数、特征数、成分数和求解设置，不把失败拟合纳入推断，也不放宽 P/q 的完整性要求。本地慢收敛回归案例在 500 步失败，在 637 步收敛。

本次同时切换到 10 区 meta 定义，应按下面完整命令从预检开始运行到**新的 meta 输出目录** `outputs/fig03/joint_pls_432_FD05_meta_pfc_both_iter10000`，不能只补旧 80 区批次的 11 次置换。旧结果保留。样本、核心代码或求解设置变化时使用新目录，不要更改旧 manifest 中的哈希来混用结果。

`fit` 完成会打印 `PFC=10, remaining=422`，应在开始长时间重抽样之前核对。若仍为 80，使用的成员表未切换成功。

## 置换和 bootstrap：两个独立分支

**置换 1000 次，单侧检验正向 ΔR²。** 主检验不是把 Y 整体乱序后声称检验了 M2−M1。它使用近似 Freedman–Lane：在共同样本上拟合约简线性均值 `Y~性别+meanFD+站点哑变量+全部PFC`；保留拟合值，把二维残差按被试整行在站点内置换，构造 Y*，每次完整重跑两套嵌套 LOSO（包括 q 与成分数选择）。

这保留了约简模型解释的 PFC/协变量关系，检验是否有超出该约简模型的条件信息；**不是无假设的两个预测算法相等检验**。有效性依赖约简条件均值的线性设定和站点内残差行可交换性，P 是近似置换 P。完整样本约简拟合仅用于生成零分布，不用于外层测试预测。配对的年龄＋成功残差使用同一索引，不分别打乱两列。

`P=(1+置换统计量≥观察统计量的次数)/(1000+1)`。只有请求的 1000 次全部完成且有效才给出 P；失败记录和不足次数不会自动被当作不显著。`q_primary` 为两套方案 M2−M1 的 BH q。

**Bootstrap 1000 次，独立运行，不嵌套在置换里。** 每个站点有 N_s 人，就在该站点有放回抽 N_s 人；总量仍为共同 N。每次重新执行嵌套 LOSO，重复抽中的人始终属于原站点，不跨训练/测试。它估计固定这些站点条件下的模型差异和第一成分权重稳定性，不代表抽取新站点的总体不确定性。第一成分的 bootstrap ratio 不是脑区逐一校正后的显著性。

`--null pfc` 的现有重复可同时为 M2−M1 与新 M3−M1 计算各自的 P；后者仍是在指定 PFC 约简线性均值成立时评估该差值，不是任意非嵌套算法相等的无条件检验。`--null covariates` 为三个脑模型各自相对 M0 的比较提供原始探索性 P，报告再把这三项×两套方案的六个 P 做 BH，输出 `q_exploratory`。当前报告不再使用 `--null nonpfc`，不需要执行它。Bootstrap 不生成 P。

## 已完成 meta 分析后只更新比较与图表

本次不修改冻结哈希覆盖的 `models.py`、`pipeline.py`、`resampling.py`、`io.py`。每次旧重复已保存 M3−M0 与 M1−M0；新报告严格按同一方案、同一 replicate 计算 `(M3−M0)−(M1−M0)=M3−M1`，再用该比较自己的零分布计算 P。Bootstrap 区间也逐次恢复后取百分位，绝不相减两个区间端点。原始计算记录保留，面向阅读的表在 `report/` 中替换为新比较。

若自提出 M3−M1 比较后还没有同步代码，且已经完成 meta-PFC 10区批次的 1000 次 PFC 置换与 1000 次 bootstrap，请把以下文件传到服务器相同的 `analysis/fig03/joint_pls/` 目录：

| 文件 | 用途 |
|---|---|
| `comparison_summary.py` | **新增且必需**；逐次恢复五项配对比较，计算 M3−M1 的 P、BH q 和区间 |
| `report.py` | 更新 CSV、Word、结果说明并调用新版绘图 |
| `visualize.py` | 两套独立 SVG、高饱和度固定配色、两栏精简布局，P/q 直接标在差值图上 |
| `run.py` | 更新命令入口并支持 `--model-ci` |
| `bootstrap_endpoints.py` | 建议一并上传，供以后显式开启单模型端点CI；当前默认不运行 |
| `README.md`、`README_RESULTS.md` | 执行及结果读取说明 |

本次结果更新保留服务器上的 `models.py`、`pipeline.py`、`resampling.py` 和 `io.py`，以便继续使用原冻结批次的计算记录。无需上传本地 mock、缓存或整个 outputs 目录。以下命令只针对已经算好的 **meta-PFC 10区** 批次；先核对结果目录中的 `precheck/sample_metadata.json`，应为 `source=Buhle2014_reappraisal_PFC`、`n_pfc=10`、`n_remaining=422`。目录名本身不决定实际脑区定义。

重新生成报告和图表：

```bash
cd /public/home/dingrui/fmri_analysis/zz_analysis
conda activate brain
python -m analysis.fig03.joint_pls.report \
  --output outputs/fig03/joint_pls_432_FD05_meta_pfc_both_iter10000
```

**无需重跑 fit、置换或 bootstrap，也不需要默认执行端点区间恢复。** 默认关闭单模型 R² 的CI，仅保留各比较 ΔR² 的配对95%CI。即使目录已经有恢复的端点区间，默认也不会读取或显示。

上述命令只导出表和图，不拟合模型；默认同时更新 Word 三线表。仅需 CSV 和 SVG 时加 `--skip-docx`；加 `--preview` 可额外生成 PNG 预览。

若需要前三项比较 M1−M0、M2−M0、M3−M0 也显示 P，且此前尚未完成相应零假设置换，在报告命令之前补跑：

```bash
python -m analysis.fig03.joint_pls.run permutation \
  --output outputs/fig03/joint_pls_432_FD05_meta_pfc_both_iter10000 \
  --null covariates --total 1000 --workers 4
```

这会新增协变量零假设的置换分支，不替换已有的 PFC 零假设或 bootstrap；相同配置下已有重复自动跳过。当前无需运行 `--null nonpfc`。全部比较应读取新 `report/` 子目录的 CSV，根目录旧比较表不包含本次汇总得到的新比较。

### 已有全部 P 后补齐探索性 BH q

若 PFC 和 covariates 两个零假设均已完成，只需同步新版 `comparison_summary.py`、`report.py`、`visualize.py`，再执行上述 report 命令。无需新增独立补算脚本，无需重跑置换或 bootstrap。报告新增 `q_exploratory`（相对 M0 的六项比较共同 BH）、`q_bh`（每行对应检验族的 q）、`bh_family_size` 和 `supported_exploratory`。原 `q_primary`、`q_secondary` 保持原校正范围；每个检验族只有全部 P 就绪才输出 q，不用部分检验临时缩小校正范围。图统一标 `q_bh`，不回退显示 P。M3−M2 / M2−M3 已从比较表和图中去掉。

主要输出：`report/joint_all_model_comparisons.svg` 和 `report/shared_all_model_comparisons.svg`，每张图各五行，分别展示该方案的所有比较；`model_comparisons_full_precision.csv` 合计10行。`pfc_comparisons.csv` 仍保留原来四行 PFC 比较。区间只使用全部有效的请求重复；P/q 仍来自相应条件置换和既定校正族。

### 杠铃图与置信区间

每个模型固定一种高饱和度颜色，跨比较和两套方案始终一致：**M0 蓝色（#0072FF）、M1 绿色（#00A651）、M2 紫色（#8B2CFF）、M3 橙色（#FF7A00）**。颜色按模型身份分配，不按 reference/compared 角色分配。左图两个端点轻微上下错开并连接灰线；默认只画 R² 点估计。右图保留 ΔR² 和配对95%CI，在对应差值点上方**统一显示 BH q**，读取 `q_bh`，并注明 `BH-adjusted q`。没有就绪的 q 则不标注，不回退显示 P。CSV 和 Word 同时保留原始 P 与校正后 q；检验族范围为主比较2项、次要比较2项、探索性比较6项，分别校正。删除独立统计列、重复数值和长说明，完整信息保留在 CSV 和结果说明中。CI 不能从两个端点区间相减得到。

若以后需要单模型端点区间，先恢复（只需一次），再显式加 `--model-ci`：

```bash
python -m analysis.fig03.joint_pls.bootstrap_endpoints \
  --output outputs/fig03/joint_pls_432_FD05_meta_pfc_both_iter10000 --workers 4
python -m analysis.fig03.joint_pls.report \
  --output outputs/fig03/joint_pls_432_FD05_meta_pfc_both_iter10000 --skip-docx --model-ci
```

可选恢复使用原抽样索引，每折只重建协变量调整和全脑参考第一 Y 方向，再用原配对差值恢复绝对 R²，不重新做内层调参或候选模型拟合。需要原预检 `data.npz`、全部有效 bootstrap CSV/索引与原 sklearn 环境。缓存位于 `report_cache/model_r2_bootstrap/`。显式要求 `--model-ci` 但区间缺失或原 bootstrap 未完成时会报错，不自动计算或静默跳过。`run report` 也支持该开关。

可选恢复输出 `bootstrap_model_r2_intervals.csv`、`paired_bootstrap_model_r2.csv` 和 `bootstrap_model_r2_provenance.json`。作图数据中的 `endpoint_interval=not_shown_optional` 表示本次主动不显示端点CI；`reference_color/compared_color` 记录端点实际模型的固定颜色。旧版合并两套方案的图不再由默认 report 生成，历史文件保留；请读取上述两个新文件名。

## 上传哪些工程文件

上传整个 `analysis/fig03/joint_pls/` 目录中的源文件和两份 README，无需打包，无需上传 `__pycache__` 或本地 mock 输出。具体文件：

`__init__.py`、`io.py`、`models.py`、`pipeline.py`、`sample.py`、`resampling.py`、`comparison_summary.py`、`bootstrap_endpoints.py`、`visualize.py`、`report.py`、`run.py`、`test_joint_pls.py`、`test_comparison_summary.py`、`test_bootstrap_endpoints.py`、`README.md`、`README_RESULTS.md`。

沿用服务器工程已有的 `predictive_pls/precheck.py`、`predictive_pls/run_pls.py`、`gam_supplement/prepare_inputs.py` 和 `gam_pygam/result_tables.py`；这些计算依赖本次没有修改。不需要运行旧版 PLS 的 fit。需要 numpy、pandas、scikit-learn、scipy、nibabel、threadpoolctl、matplotlib、python-docx。

Meta 版本还需上传已有的 `outputs/fig03/predictive_pls_432_FD05_meta_pfc/pfc_definition/pfc_parcel_membership.csv`。服务器的 `predictive_pls/run_pls.py` 必须是已经支持 Buhle 成员表的版本。若此前已运行过 meta 分析，该依赖和成员表通常已在服务器。本次同步新版 `models.py`、`run.py`，以及此前 meta 支持所需的 `sample.py`、`report.py`、`test_joint_pls.py` 和两份 README；建议同步整个 `analysis/fig03/joint_pls/` 源码目录以避免遗漏。联合预测和共有成分仍一起输出。

## 服务器执行命令

以下路径沿用既有工程。先确认冻结的 PFC 成员表在服务器；`--volume-dir` 是**原始 effectSize 数据目录**。如果旧缓存少了因父母教育缺失而被排除的人，预检从这里恢复脑区均值，并要求至少一个缓存被试的提取结果可复现，防止混用图谱或提取算法。默认恢复图谱为 `--inputs/combined_atlas.nii.gz`，不是单独的 Schaefer400。

```bash
cd /public/home/dingrui/fmri_analysis/zz_analysis
conda activate brain

python -m analysis.fig03.joint_pls.run precheck \
  --inputs outputs/fig03/gam_pygam_volume_432_models3_FD05/inputs \
  --demographics /public/home/dingrui/fmri_analysis/data/beh/participants_in_tfmri_demographics.csv \
  --subjects /public/home/dingrui/fmri_analysis/data/beh/sub_list_valid_4_ER.txt \
  --excluded /public/home/dingrui/fmri_analysis/data/fmri/subs_removed_from_2ndlevel.txt \
  --meanfd /public/home/dingrui/fmri_analysis/data/fmri/qc_res_subs_runs_fmri.csv \
  --behavior /public/home/dingrui/fmri_analysis/data/beh/er_profiles_new.csv \
  --id-source /public/home/dingrui/fmri_analysis/data/beh/er_profiles.csv \
  --membership outputs/fig03/predictive_pls_432_FD05_meta_pfc/pfc_definition/pfc_parcel_membership.csv \
  --volume-dir /public/home/dingrui/fmri_analysis/data/fmri/res_1stlevel \
  --output outputs/fig03/joint_pls_432_FD05_meta_pfc_both_iter10000

python -m analysis.fig03.joint_pls.run fit --output outputs/fig03/joint_pls_432_FD05_meta_pfc_both_iter10000 --k-max 10
python -m analysis.fig03.joint_pls.run permutation --output outputs/fig03/joint_pls_432_FD05_meta_pfc_both_iter10000 --null pfc --total 1000 --workers 4
python -m analysis.fig03.joint_pls.run bootstrap --output outputs/fig03/joint_pls_432_FD05_meta_pfc_both_iter10000 --total 1000 --workers 4
python -m analysis.fig03.joint_pls.run report --output outputs/fig03/joint_pls_432_FD05_meta_pfc_both_iter10000
```

`--id-source` 使用已核实的行对应关系，按 `er_profiles.csv` 各数据行前六位取得真实 ID，保留完整对齐审计。行为表默认字段 `sub_ids` 与 `reappraisal_success`；需要其他成功指标时显式传 `--score-column`。预检从不默认采用旧的 1169 或 991；实际共同 N 以新输出为准。

若已有一份确实**未按父母教育筛选**的完整设计表，也可以用 `--cohort-design /path/design.csv`，替代 `--demographics/--subjects/--excluded` 三个参数。其他参数不变。只有旧缓存完整保存排除审计时才可不传这四种样本来源参数。

若希望所有探索性模型比较也有 P，在 report 之前执行：

```bash
python -m analysis.fig03.joint_pls.run permutation --output outputs/fig03/joint_pls_432_FD05_meta_pfc_both_iter10000 --null covariates --total 1000 --workers 4
```

同一条置换或 bootstrap 命令可续跑，会跳过已有重复。可用 `--start 0 --count 100 --total 1000` 先跑前 100 次；后续省略 start/count 即补齐。不要同时运行写同一个输出目录的任务。样本、核心代码、k、seed、total 不一致会拒绝混用；发现失败后先诊断，修复代码或改变配置后使用新的输出目录。观察拟合用时写在 `timing.json`，服务器耗时应据实际预检和小批重复估计。

## 本地验证（模拟数据不能写进论文）

```bash
python -m unittest analysis.fig03.joint_pls.test_joint_pls analysis.fig03.joint_pls.test_comparison_summary analysis.fig03.joint_pls.test_bootstrap_endpoints -v
python -m analysis.fig03.joint_pls.run mock --output outputs/fig03/joint_pls_mock_final --permutations 3 --bootstraps 3
```

测试覆盖 sklearn 多输出 PLS 系数复现、测试站点 Y 隔离、共有目标一致、年龄单位不影响结果、成对站点内置换、bootstrap 固定 N、失败/未完成不输出 P、真实 ID 行映射与教育排除恢复。Mock 只是工程贯通，三次重复不能用于正式推断。

方法依据：[PLSRegression 官方定义](https://scikit-learn.org/stable/modules/generated/sklearn.cross_decomposition.PLSRegression.html)、[多输出 R²](https://scikit-learn.org/stable/modules/generated/sklearn.metrics.r2_score.html)、[Winkler 等关于约简模型与可交换块的置换方法](https://pmc.ncbi.nlm.nih.gov/articles/PMC4010955/)。
