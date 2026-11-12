# 年龄贝叶斯样条模型比较

结果出来后如何阅读：见 [结果读取和解释指南](README_RESULTS.md)，包括 BF31/BF32/BF21 的方向、证据标签、先验敏感性、原年龄 q 与 BF 的合并解释及论文报告模板。

本入口只读取当前 Tian432/232、FD≤0.5 的年龄分析缓存，增加正式贝叶斯模型比较。不重新提取 contrast、不重拟合原 pyGAM、不改变原 F/p/q、曲线或地图。ERQ 本轮不计算 BF。

## 模型与先验

同一单元三个模型使用相同完整被试和原年龄分析的四个数值协变量编码：

```text
M1 = intercept + C
M2 = intercept + C + beta_age * age_z
M3 = intercept + C + beta_age * age_z + Z(age) * u
C = gender + site_id + parent_edu + meanFD
```

这里保留年龄入口现有的数值站点编码，避免补 BF 时悄悄改变原模型；ERQ 入口仍使用它已经采用的分类站点编码，两者不混用。

年龄使用每单元有效样本的总体标准差标准化。构造与原模型相同的五个三次 B 样条，使用原二阶差分惩罚 `D2.T @ D2`。将正特征值方向按惩罚白化，去掉截距和线性成分，得到三列非线性基函数 Z。验证 `{1, age_z, Z}` 完全张成原五维样条空间，随后使 `tr(Z Z') / n = 1`。

主分析先验：

```text
beta_age | sigma ~ Normal(0, (r_linear * sigma)^2)
u | sigma,tau ~ Normal(0, sigma^2 * tau^2 * I_3)
tau ~ HalfNormal(0, r_nonlinear)
r_linear = r_nonlinear = 0.5
```

M2/M3 共享相同线性年龄先验。Z 的归一化使 tau 描述非线性函数在观测年龄分布上的先验 RMS，相对于残差标准差 sigma，而不是任意基函数系数尺度。0.5 是本项目预先确定的候选先验，不宣称是该领域通用默认值。

共同截距和协变量采用可识别 nuisance 列空间中的相同 Lebesgue 参考测度，共同方差采用 `d(sigma²)/sigma²`。这些参考先验不正规，但只用于所有模型共同存在的参数；同一数据和参考测度的任意常数在 BF 中抵消。新增年龄参数和非线性强度的先验全部正规。代码只报告相对 M1 的对数证据，不报告具有任意常数的绝对边际密度。秩亏的重复 nuisance 列只保留其可识别列空间；年龄被 nuisance 完全混淆的单元不计算 BF。

该分解与贝叶斯平滑先验的思路见 [Sabanés Bové 等的样条模型比较方法](https://arxiv.org/html/1108.3520v2)。本工程使用上述已明确约定的先验，不是直接复现该论文的 hyper-g 默认先验。

## 计算和验证

先投影消去共同 nuisance 空间，在同一尺度先验下解析积分年龄系数、非线性系数及残差方差。M3 再对半正态非线性强度进行自适应数值积分；额外采用不同分段和更严格容差独立复核，要求相对积分误差及 log 证据差均不超过 `1e-6`。不使用 GCV 最优 lambda 代替积分，也不采用 BIC 近似。

```text
log BF31 = log evidence(M3) - log evidence(M1)
log BF32 = log evidence(M3) - log evidence(M2)
log BF21 = log evidence(M2) - log evidence(M1)
```

每个有效单元计算 `r_linear,r_nonlinear ∈ {0.25,0.5,1}` 的九套组合，保存全部结果及每种比较的敏感性范围。M2 不依赖非线性先验尺度。三种 BF 使用同一套模型证据，满足 `log BF31 - log BF32 = log BF21`。

BF10 > 1 支持比较中的第一个模型，BF10 < 1 支持第二个模型。3 和 10 的标签仅用于描述证据强度，不是显著性或多重比较校正。BF32 很小支持线性模型相对于非线性模型，不能表述为没有年龄关联。BF 不做 BH；若将全脑 BF 用于发现判定，需要另行制定模型先验/多重比较策略。

线性尺度 BF 极大或极小时可出现浮点上溢/下溢，log BF 仍为主要记录；不人为截断 log BF 或把数值极限当作确切无穷证据。常量响应、年龄与协变量混淆、样本不足等会保留 status，不能以 BF=0 代替失败。

## 云端替换和命令

按原工程路径同步以下文件，不需要打包：

- `analysis/fig03/gam_pygam/age_bayes.py`
- `analysis/fig03/gam_pygam/run_age_bayes.py`
- `analysis/fig03/gam_pygam/result_tables.py`
- `analysis/fig03/gam_pygam/README_BAYES.md`
- `analysis/fig03/gam_pygam/README_RESULTS.md`
- `analysis/fig03/gam_pygam/README.md`

复用云端当前 `run_pygam.py` 和其输入准备依赖。环境沿用原 pyGAM 环境；Word 导出额外需要 python-docx：

```bash
python -m pip install python-docx

python -m analysis.fig03.gam_pygam.run_age_bayes \
  --inputs "outputs/fig03/gam_pygam_volume_432_models3_FD05/inputs" \
  --output "outputs/fig03/gam_pygam_volume_432_models3_FD05/bayes"

python -m analysis.fig03.gam_pygam.run_age_bayes \
  --inputs "outputs/fig03/gam_pygam_volume_232_models3_FD05/inputs" \
  --output "outputs/fig03/gam_pygam_volume_232_models3_FD05/bayes"
```

若已有 `statistics/pygam_age_comparison.csv`，会按范围和名称合并原频率学统计，并验证每个单元的 N 与 BF 一致；原表不改写。432 输出 449 单元，232 只输出 232 脑区。实际缓存目录不同只调整命令路径，入口会核验保留样本的 meanFD 全部≤0.5。

| 输出文件 | 内容 |
|---|---|
| `age_bayes_statistics.csv` | 主先验的全量宽表、三种 BF/log BF、九套先验敏感性摘要及原统计（若存在） |
| `age_bayes_prior_sensitivity.csv` | 每有效单元九套先验组合的完整证据、BF 和积分精度 |
| `age_bayes_comparisons.csv` | 每单元三行的全精度比较长表 |
| `age_bayes_manuscript_table.docx` | 全量可编辑 Word 三线表，按显示需要舍入 |
| `unit_response_exclusions.csv` | 各单元共同排除的真实被试 ID |
| `analysis_design.csv` | 本次读取的原年龄分析纳入样本 |
| `run_metadata.json` | 完整先验、尺度、积分方法、源数据哈希、范围与状态计数 |

仅重新导出 Word/CSV、不重新积分：

```bash
python -m analysis.fig03.gam_pygam.run_age_bayes \
  --step tables \
  --output "outputs/fig03/gam_pygam_volume_432_models3_FD05/bayes"
```

如果分析环境暂时没有 python-docx，可以拟合时加 `--skip-docx`，之后在有 pandas、python-docx 的环境中使用独立 `result_tables` 入口导出，不需要 pyGAM 或 SciPy：

```bash
python -m analysis.fig03.gam_pygam.result_tables \
  --kind age \
  --statistics "outputs/fig03/gam_pygam_volume_432_models3_FD05/bayes/age_bayes_statistics.csv" \
  --output "outputs/fig03/gam_pygam_volume_432_models3_FD05/bayes"
```

## 测试

```bash
python -m unittest analysis.fig03.gam_pygam.test_age_bayes -v
```

验证包括独立的完整 n×n 协方差边际积分、无效应/线性/非线性模拟、响应单位和 nuisance 平移不变性、九套先验的一致性、缺失/不可识别情况、432/232 范围及源缓存不变。测试依赖 statsmodels（合成缓存辅助函数使用），真实统计需执行上述云端命令。

全量模拟结果保存在 `outputs/fig03/inference_mock/age/`：449 单元、4041 行先验敏感性结果，以及 1347 行 BF 比较长表。模拟数据由 `python -m analysis.fig03.gam_pygam.mock_inference` 生成。模拟 Word 使用独立导表入口并加 `--mock`，明确标记 SIMULATED：

```bash
python -m analysis.fig03.gam_pygam.result_tables \
  --kind age \
  --statistics "outputs/fig03/inference_mock/age/age_bayes_statistics.csv" \
  --output "outputs/fig03/inference_mock/age" --mock
```

Word 已检查三线边框、表头重复、全量行数和替代渲染的可读性。本机缺少 LibreOffice，标准 DOCX 渲染器无法执行；替代渲染不能保证与 Word 原生分页一致。生成 Word 本身不需要 LibreOffice。
