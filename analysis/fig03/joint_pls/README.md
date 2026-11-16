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

- 只做 432 分区、meanFD≤0.5、宽范围解剖 PFC 主分析。冻结成员表必须是 `Donahue2018_liberal_bilateral`，不根据本次显著结果挑选脑区。
- 年龄、成功、性别、meanFD、站点和全部 432 脑区完整的共同样本，两套方案和所有模型完全一致。
- **父母教育既不是协变量，也不是纳入条件。** 推荐从原始人口学表＋有效名单＋二级排除名单重建样本，避免沿用旧样本隐含的教育筛选。
- 性别和 meanFD 为预测协变量。仅在训练数据拟合 `X~C`、`Y~C`，再对残差标准化；应用到验证和测试数据。年龄属于 Y，不作为协变量。Y 的训练残差各自单位方差，避免年龄单位主导拟合。
- 站点用于内外层 LOSO 和限制重抽样。没有给新站点拟合截距或用测试站点数据做回归调整。有效站点指每个被试有真实可对齐的站点标签且有至少三个站点；不能用行号构造站点。
- 外层留一站点测试，内层留一训练站点选成分数。候选 k 为 0–10（样本量/特征数不足时自动缩减为所有模型共同可行范围），k=0 允许神经模型回退到协变量基线。并列优先较少成分。内层最大化与该方案一致的 pooled R²。

四模型：M0 协变量；M1 协变量＋PFC；M2 协变量＋全部 432；M3 协变量＋非 PFC。主比较是 M2−M1。两套方案各有一个主比较，共两个检验，放在同一族做 BH。其他比较预设输出，属于探索性。

## 置换和 bootstrap：两个独立分支

**置换 1000 次，单侧检验正向 ΔR²。** 主检验不是把 Y 整体乱序后声称检验了 M2−M1。它使用近似 Freedman–Lane：在共同样本上拟合约简线性均值 `Y~性别+meanFD+站点哑变量+全部PFC`；保留拟合值，把二维残差按被试整行在站点内置换，构造 Y*，每次完整重跑两套嵌套 LOSO（包括 q 与成分数选择）。

这保留了约简模型解释的 PFC/协变量关系，检验是否有超出该约简模型的条件信息；**不是无假设的两个预测算法相等检验**。有效性依赖约简条件均值的线性设定和站点内残差行可交换性，P 是近似置换 P。完整样本约简拟合仅用于生成零分布，不用于外层测试预测。配对的年龄＋成功残差使用同一索引，不分别打乱两列。

`P=(1+置换统计量≥观察统计量的次数)/(1000+1)`。只有请求的 1000 次全部完成且有效才给出 P；失败记录和不足次数不会自动被当作不显著。`q_primary` 为两套方案 M2−M1 的 BH q。

**Bootstrap 1000 次，独立运行，不嵌套在置换里。** 每个站点有 N_s 人，就在该站点有放回抽 N_s 人；总量仍为共同 N。每次重新执行嵌套 LOSO，重复抽中的人始终属于原站点，不跨训练/测试。它估计固定这些站点条件下的模型差异和第一成分权重稳定性，不代表抽取新站点的总体不确定性。第一成分的 bootstrap ratio 不是脑区逐一校正后的显著性。

`--null covariates` 为三个脑模型各自相对 M0 的比较提供 P；`--null nonpfc` 为 M2−M3 提供 P。这两种可选零分布与主检验分开生成，探索性 P 暂不做额外校正。Bootstrap 不生成 P。

## 上传哪些工程文件

上传整个 `analysis/fig03/joint_pls/` 目录中的源文件和两份 README，无需打包，无需上传 `__pycache__` 或本地 mock 输出。具体文件：

`__init__.py`、`io.py`、`models.py`、`pipeline.py`、`sample.py`、`resampling.py`、`report.py`、`run.py`、`test_joint_pls.py`、`README.md`、`README_RESULTS.md`。

沿用服务器工程已有的 `predictive_pls/precheck.py`、`predictive_pls/run_pls.py`、`gam_supplement/prepare_inputs.py` 和 `gam_pygam/result_tables.py`；这些计算依赖本次没有修改。不需要运行旧版 PLS 的 fit。需要 numpy、pandas、scikit-learn、scipy、nibabel、threadpoolctl、matplotlib、python-docx。

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
  --membership outputs/fig03/predictive_pls_432_FD05/pfc_definition/pfc_parcel_membership.csv \
  --volume-dir /public/home/dingrui/fmri_analysis/data/fmri/res_1stlevel \
  --output outputs/fig03/joint_pls_432_FD05_both

python -m analysis.fig03.joint_pls.run fit --output outputs/fig03/joint_pls_432_FD05_both --k-max 10
python -m analysis.fig03.joint_pls.run permutation --output outputs/fig03/joint_pls_432_FD05_both --null pfc --total 1000 --workers 4
python -m analysis.fig03.joint_pls.run bootstrap --output outputs/fig03/joint_pls_432_FD05_both --total 1000 --workers 4
python -m analysis.fig03.joint_pls.run report --output outputs/fig03/joint_pls_432_FD05_both
```

`--id-source` 使用已核实的行对应关系，按 `er_profiles.csv` 各数据行前六位取得真实 ID，保留完整对齐审计。行为表默认字段 `sub_ids` 与 `reappraisal_success`；需要其他成功指标时显式传 `--score-column`。预检从不默认采用旧的 1169 或 991；实际共同 N 以新输出为准。

若已有一份确实**未按父母教育筛选**的完整设计表，也可以用 `--cohort-design /path/design.csv`，替代 `--demographics/--subjects/--excluded` 三个参数。其他参数不变。只有旧缓存完整保存排除审计时才可不传这四种样本来源参数。

若希望所有探索性模型比较也有 P，在 report 之前执行：

```bash
python -m analysis.fig03.joint_pls.run permutation --output outputs/fig03/joint_pls_432_FD05_both --null covariates --total 1000 --workers 4
python -m analysis.fig03.joint_pls.run permutation --output outputs/fig03/joint_pls_432_FD05_both --null nonpfc --total 1000 --workers 4
```

同一条置换或 bootstrap 命令可续跑，会跳过已有重复。可用 `--start 0 --count 100 --total 1000` 先跑前 100 次；后续省略 start/count 即补齐。不要同时运行写同一个输出目录的任务。样本、核心代码、k、seed、total 不一致会拒绝混用；发现失败后先诊断，修复代码或改变配置后使用新的输出目录。观察拟合用时写在 `timing.json`，服务器耗时应据实际预检和小批重复估计。

## 本地验证（模拟数据不能写进论文）

```bash
python -m unittest analysis.fig03.joint_pls.test_joint_pls -v
python -m analysis.fig03.joint_pls.run mock --output outputs/fig03/joint_pls_mock_final --permutations 3 --bootstraps 3
```

测试覆盖 sklearn 多输出 PLS 系数复现、测试站点 Y 隔离、共有目标一致、年龄单位不影响结果、成对站点内置换、bootstrap 固定 N、失败/未完成不输出 P、真实 ID 行映射与教育排除恢复。Mock 只是工程贯通，三次重复不能用于正式推断。

方法依据：[PLSRegression 官方定义](https://scikit-learn.org/stable/modules/generated/sklearn.cross_decomposition.PLSRegression.html)、[多输出 R²](https://scikit-learn.org/stable/modules/generated/sklearn.metrics.r2_score.html)、[Winkler 等关于约简模型与可交换块的置换方法](https://pmc.ncbi.nlm.nih.gov/articles/PMC4010955/)。
