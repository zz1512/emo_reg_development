# 体积 GAM 补充分析

这份分析回答审稿意见里的两个问题。

1. 脑区、Yeo 7 网络和四个情绪调节系统的年龄轨迹是否偏离直线。
2. 脑区与脑区、网络与网络、系统与系统的年龄轨迹是否不同。

它不替换正文原来的 Python `pygam` GAM，也不写进 `analysis/fig03/run.py`。模型在 R 的 `mgcv` 里拟合。Python 只准备被试均值、设计矩阵和作图。

已完成的体积分析样本是 1171 人，结果目录是 `outputs/fig03/gam_supplement_volume`。确切名单在该目录的 `inputs/design.csv`。这张表只含进入模型的人，没有 `in_analysis` 列。后面的嵌套模型和 pygam 把同一次筛选写成包含全部被试的 `design.csv`，用 `in_analysis` 和 `excluded_by` 标出是否进入分析、排除原因，说明在 `analysis/fig03/gam_nested/README.md`。

## 方法

### 样本和协变量

样本从原来的情绪调节有效名单出发，去掉二阶分析剔除的被试，再依次接上协变量和体积图。

1. `participants_in_tfmri_demographics.csv` 提供 `sub_id`、`age`、`gender`、`site_id`。`gender` 保持 0/1。`site_id` 保持数值，不转成因子。
2. 父母教育来自 `env_family_student_cbcl.tsv`。被试编号列是 `被试编号`，取值列是 `父母受教育水平_父母受教育水平`。编号去掉 `sub-`、忽略大小写后对齐。空白、表里没有、同一被试多行数值不一致的人不进入模型。去掉的人写在 `inputs/dropped_parent_edu.csv`。
3. 体积图必须同时有再评价和被动看负性两张 `effectSize`。缺文件的编号写在 `inputs/dropped_no_volume.txt`。
4. 头动来自 `qc_res_subs_runs_fmri.csv`。每个被试取 `run_name` 为 EMO 和 SOC 的 `mean_fd` 的平均，写入 `meanFD`。分析样本里还有人对不上时，脚本停止，不删人。

四个协变量都是线性项：`gender`、`site_id`、`parent_edu`、`meanFD`。

### 因变量

因变量是体积图上的再评价减被动看负性。每个被试读：

- `res_1stlevel/<sub>/effectSize-reg_neg_Observe_across_runs.nii.gz`
- `res_1stlevel/<sub>/effectSize-pas_neg_Observe_across_runs.nii.gz`

两张图相减。被试图留在自己的网格上，这次运行的网格是 `(76, 91, 76)`。Schaefer 200 的 2 mm 图谱是 `(91, 109, 91)`，只做一次最近邻，投到效应图网格上，用来给体素标脑区号。脑区值是该标签下有限体素的平均，NaN 不进入平均。200 个脑区都有体素，所以 200 个都进入 GAM。

Yeo 7 网络不是另一张图。每个网络是它包含的脑区均值，按这些脑区在效应图网格上的体素数加权平均。顺序是 Vis、SomMot、DorsAttn、SalVentAttn、Limbic、Cont、Default。

四个系统用 `res/res_bayes/mask_<system>.nii.gz`。掩膜在效应图网格上，系统值是掩膜内有限体素的对比均值。这次读到的体素数是：

| 系统 | 掩膜体素 |
|---|---|
| `reappraisal_only` | 954 |
| `common_appraisal` | 177 |
| `nonmodifiable_emotion` | 2648 |
| `modifiable_emotion` | 177 |

20% 规则只给脑区贴系统标签，供扇形图上色。系统 GAM 用的是整个掩膜的均值。一个脑区只有在某一系统里的体素占比达到 0.20 时，`parcel_info.csv` 的 `system` 才写成该系统。这次达到标准的是再评价特异 3 个、不可调节情绪 9 个；共同评估和可调节情绪是 0 个。

### 年龄模型

每个脑区、网络和系统各拟合一次，式子相同：

```text
y ~ age + h + gender + site_id + parent_edu + meanFD
```

写成

\[
y=\beta_0+\beta_{\text{age}}\cdot\text{age}+\gamma\cdot h(\text{age})+\beta_g\cdot\text{gender}+\beta_s\cdot\text{site}+\beta_p\cdot\text{parent\_edu}+\beta_m\cdot\text{meanFD}+\varepsilon
\]

`h(age)` 来自 `s(age, k=3, bs="cr")`。这个三次回归基里已经包含截距和直线。脚本先去掉与截距和直线年龄平行的部分，留下一个和它们正交的方向，再乘上系数 \(\gamma\)。这个方向带自己的三次回归惩罚，平滑参数由 REML 选择。拟合函数是 `mgcv::gam()`。脚本不调用 `anova()`。

\(\gamma=0\) 就是这条弯曲函数处处为 0。`p_smooth` 是 \(\gamma\) 的 Wald 检验，`p_linear` 是 \(\beta_{\text{age}}\) 的 Wald 检验。p 值比较的是系数除以标准误，不是系数的绝对值。残差大小和平滑参数在每个模型里单独估计，所以 `gamma_h` 更大时，`p_smooth` 可以更大。标准误在 `se_h`、`se_age` 列，不在 `linearity_deviation.csv` 里。

BH 校正按下面的范围分开做。

- 200 个脑区的 `p_linear`、`p_smooth` 在 200 个脑区里校正。
- 7 个网络在 7 个网络里校正。
- 4 个系统在 4 个系统里校正。
- 第二步里，网络对、系统对、脑区对各自校正。7 个网络有 42 个有序对，4 个系统有 12 个有序对，200 个脑区有 39800 个有序对。A 减 B 和 B 减 A 方向相反、p 值相同，两格都进入校正。

`q_smooth < 0.05` 表示偏离直线。`q_linear < 0.05` 表示直线斜率不为 0。

带符号的偏 \(R^2\) 只描述年龄项在四个协变量之外多解释了多少变异。符号是中心化曲线从最小年龄到最大年龄的净变化方向。它不是显著性。

作图用的曲线在样本平均年龄处减成 0。曲线在平均年龄附近相交，是这个中心化的结果。

### 轨迹比较

第二步对差值拟合同一个模型：

\[
d(\text{age})=f_A(\text{age})-f_B(\text{age})
\]

被试水平先做 A 的均值减 B 的均值，再拟合 `y ~ age + h + gender + site_id + parent_edu + meanFD`。

- 这一行的 `beta_age` 是差值的直线斜率，也就是两条斜率之差。`q_linear` 检验它是否为 0。
- 这一行的 `gamma_h` 是差值的弯曲系数。`q_smooth` 检验两条轨迹的弯曲是否不同。

比较做三组：Yeo 网络的有序对、四个系统的有序对、200 个 Schaefer 脑区的有序对。脑区对只保存系数和方阵，不保存每对的曲线，也不画成轨迹格子。

## 运行

在工程根目录执行。第一步写线性和弯曲项，第二步读第一步的输入做比较。

```bash
python -m analysis.fig03.gam_supplement.run_supplement --step terms --input "/public/home/dingrui/fmri_analysis/data/beh/participants_in_tfmri_demographics.csv" --subjects "/public/home/dingrui/fmri_analysis/data/beh/sub_list_valid_4_ER.txt" --excluded "/public/home/dingrui/fmri_analysis/data/fmri/subs_removed_from_2ndlevel.txt" --parent-edu "/public/home/dingrui/fmri_analysis/data/beh/env_family_student_cbcl.tsv" --meanfd "/public/home/dingrui/fmri_analysis/data/fmri/qc_res_subs_runs_fmri.csv" --fig3-volume-dir "/public/home/dingrui/fmri_analysis/data/fmri/res_1stlevel" --fig3-volume-atlas "resources/atlases/schaefer2018/Schaefer2018_200Parcels_7Networks_order_FSLMNI152_2mm.nii.gz" --schaefer-atlas "resources/atlases/schaefer2018/Schaefer2018_200Parcels_7Networks_order.dlabel.nii" --system-mask-dir "/public/home/dingrui/fmri_analysis/res/res_bayes" --output "outputs/fig03/gam_supplement_volume"
```

```bash
python -m analysis.fig03.gam_supplement.run_supplement --step matrix --schaefer-atlas "resources/atlases/schaefer2018/Schaefer2018_200Parcels_7Networks_order.dlabel.nii" --output "outputs/fig03/gam_supplement_volume"
```

脚本在 `analysis/fig03/gam_supplement/`：

| 文件 | 作用 |
|---|---|
| `prepare_inputs.py` | 对齐样本，计算脑区、网络和系统均值 |
| `mgcv_age_decomposition.R` | 拟合 `gam()`，写统计表和曲线 |
| `plot_suites.py` | 写线性偏离表、方阵和图 |
| `run_supplement.py` | 两步入口 |

## 输出文件

下面的路径都相对于 `outputs/fig03/gam_supplement_volume`。

### 输入

| 文件 | 内容 |
|---|---|
| `inputs/design.csv` | 进入模型的被试和协变量。这一目录没有排除行 |
| `inputs/parcel_response.csv` | 每个被试、每个脑区的 Reappraise − Look 均值 |
| `inputs/parcel_info.csv` | 脑区编号、名称、Yeo 网络、系统标签和掩膜占比 |
| `inputs/group_response.csv` | 每个被试的网络均值和系统均值，长表 |
| `inputs/dropped_parent_edu.csv` | 因父母教育去掉的被试和原因 |
| `inputs/dropped_no_volume.txt` | 缺两张效应图的被试 |
| `inputs/parent_edu_source.txt` | 父母教育的文件和列名 |
| `inputs/meanfd_source.txt` | 头动的文件和平均方式 |

### 第一步：是否偏离直线

| 文件 | 内容 |
|---|---|
| `statistics/parcel_gam_statistics.csv` | 200 个脑区各自的模型 |
| `statistics/parcel_curves.csv` | 200 条中心化年龄曲线，每个脑区 100 个年龄点 |
| `statistics/group_gam_statistics.csv` | 7 个网络和 4 个系统各自的模型 |
| `statistics/group_curves.csv` | 这 11 条中心化曲线 |
| `statistics/linearity_deviation.csv` | 脑区、网络、系统合在一张表里的偏离直线结果 |
| `statistics/model_reference.csv` | 平均年龄、`k=3`、`cr`、当前步骤名 |
| `figures/linearity_deviation_groups.png` | 网络和系统的曲线。金色或紫色表示 `q_smooth < 0.05`，灰色表示弯曲项不显著 |
| `figures/linearity_deviation_parcels.png` | `q_smooth < 0.05` 的脑区。金色是 `gamma_h > 0`，紫色是 `gamma_h < 0` |

给老师看偏离直线时，先看 `linearity_deviation_groups.png` 和 `linearity_deviation.csv`。脑区水平再看 `linearity_deviation_parcels.png`。

### 第二步：轨迹是否不同

| 文件 | 内容 |
|---|---|
| `statistics/difference_statistics.csv` | 网络对、系统对、39800 个脑区对的模型。每行是 A 减 B |
| `statistics/difference_curves.csv` | 只有网络对和系统对的差值曲线。脑区对的曲线不保存 |
| `statistics/yeo_q_linear_matrix.csv` | 网络对直线斜率差的 q。行减列 |
| `statistics/yeo_q_smooth_matrix.csv` | 网络对弯曲差的 q |
| `statistics/system_q_linear_matrix.csv` | 系统对直线斜率差的 q |
| `statistics/system_q_smooth_matrix.csv` | 系统对弯曲差的 q |
| `statistics/yeo_p_linear_matrix.csv`、`yeo_p_smooth_matrix.csv` | 校正前的 p |
| `statistics/system_p_linear_matrix.csv`、`system_p_smooth_matrix.csv` | 校正前的 p |
| `statistics/yeo_signed_partial_r2_matrix.csv`、`system_signed_partial_r2_matrix.csv` | 差值的带符号偏 \(R^2\) |
| `statistics/yeo_shape_matrix.csv`、`system_shape_matrix.csv` | 差值曲线的形状 |
| `statistics/parcel_q_linear_matrix.csv`、`parcel_q_smooth_matrix.csv` | 200×200 脑区对的 q |
| `statistics/parcel_signed_partial_r2_matrix.csv`、`parcel_net_change_matrix.csv`、`parcel_shape_matrix.csv` | 脑区对的描述统计 |
| `statistics/parcel_p_linear_matrix.csv`、`parcel_p_smooth_matrix.csv` | 脑区对校正前的 p |
| `figures/yeo_difference_matrix.png` | 7×7 的差值曲线。颜色是带符号偏 \(R^2\) |
| `figures/system_difference_matrix.png` | 4×4 的差值曲线 |
| `figures/yeo_parcel_trajectories.png` | 200 个脑区各自的中心化曲线，按方向分上升和下降 |
| `figures/system_parcel_trajectories.png` | 只画 `parcel_info.csv` 里有系统标签的脑区 |
| `figures/parcel_pair_net_change.png` | 200×200 净变化热图 |
| `figures/parcel_pair_q_smooth.png` | 200×200 的 `log10(q_smooth)` 热图 |

方阵的对角线是空的，因为一组不和自己比较。下三角是上三角的相反方向。

给老师看轨迹差异时，看两张比较图和四张 q 方阵：`yeo_difference_matrix.png`、`system_difference_matrix.png`、`yeo_q_linear_matrix.csv`、`yeo_q_smooth_matrix.csv`、`system_q_linear_matrix.csv`、`system_q_smooth_matrix.csv`。某个网络或系统自己的斜率在 `group_gam_statistics.csv` 的 `beta_age`。某对的斜率差在 `difference_statistics.csv` 的 `beta_age`。

`yeo_parcel_trajectories.png` 表示每个脑区曲线的方向，不是两两比较，也不是显著性。`parcel_pair_*.png` 留到被问到 200 对时再打开。

`system_parcel_trajectories.png` 依赖系统标签过滤。当前代码会去掉空标签。若服务器上的这张图标题仍是 `n=200`，它是过滤更新前画的，和 `yeo_parcel_trajectories.png` 是同一批脑区，需要用更新后的 `plot_suites.py` 重画。重画只读已有统计表，不必重跑 39800 次拟合。

## 重要字段

### `inputs/design.csv`

| 列 | 含义 |
|---|---|
| `sub_id` | 成像被试编号，保持原名单里的写法 |
| `age` | 年龄 |
| `gender` | 0/1 |
| `site_id` | 数值型站点 |
| `parent_edu` | `父母受教育水平_父母受教育水平` |
| `meanFD` | EMO 与 SOC 的 `mean_fd` 平均，单位毫米 |

这一目录的 `design.csv` 每一行都进入了模型。嵌套模型和 pygam 的 `design.csv` 另外有 `in_analysis`、`excluded_by`，并把未进入分析的被试留在同一张表里。

### `inputs/parcel_info.csv`

| 列 | 含义 |
|---|---|
| `parcel` | Schaefer 脑区编号，1–200 |
| `name` | 图谱标签名 |
| `network` | Yeo 网络名 |
| `system` | 掩膜占比达到 20% 的系统名。达不到时为空 |
| `system_fraction` | 该脑区落在占比最高的那个系统里的体素比例 |

### `statistics/parcel_gam_statistics.csv` 和 `group_gam_statistics.csv`

两张表的统计列相同。脑区表用 `parcel` 标识。组表用 `dimension`（`yeo` 或 `system`）和 `group` 标识。

| 列 | 含义 |
|---|---|
| `beta_age` | 直线项 `age` 的系数，即这条轨迹的斜率。`p_linear` 检验的就是它是否为 0 |
| `se_age` | `beta_age` 的标准误 |
| `p_linear` | `beta_age = 0` 的 Wald p 值，未校正 |
| `q_linear` | `p_linear` 的 BH 校正。脑区、网络、系统分开校正 |
| `gamma_h` | 弯曲项 `h(age)` 的系数 |
| `se_h` | `gamma_h` 的标准误。`p_smooth` 由 `gamma_h / se_h` 决定 |
| `p_smooth` | `gamma_h = 0` 的 Wald p 值，未校正 |
| `q_smooth` | `p_smooth` 的 BH 校正。小于 0.05 表示偏离直线 |
| `edf_h` | `h(age)` 的有效自由度。惩罚把它压向 0 时小于 1 |
| `edf_model` | 整个模型的有效自由度，含截距、直线、四个协变量和 `h` |
| `sp_h` | `h(age)` 的 REML 平滑参数。越大，弯曲越接近 0 |
| `partial_r2` | 年龄项相对“只有四个协变量”的模型多解释的变异比例 |
| `signed_partial_r2` | `partial_r2` 乘净变化的符号 |
| `net_change` | 中心化曲线在最大年龄处的值减去最小年龄处的值 |
| `shape` | `increasing`、`decreasing`、`nonmonotone` 或 `flat` |

### `statistics/linearity_deviation.csv`

由上面两张统计表拼成。`level` 是 `parcel`、`yeo` 或 `system`，`name` 是脑区编号或组名。它保留 `gamma_h`、`p_smooth`、`q_smooth`、`p_linear`、`q_linear`、`shape`、`signed_partial_r2`。`departs_from_linear` 在 `q_smooth < 0.05` 时为真。这张表没有 `beta_age` 和标准误。

### `statistics/difference_statistics.csv`

列和组统计表相同，另外有：

| 列 | 含义 |
|---|---|
| `dimension` | `yeo`、`system` 或 `parcel` |
| `group_a` | 被减的一组。脑区对里是脑区编号 |
| `group_b` | 减去的一组 |
| `beta_age` | A 减 B 之后，差值轨迹的直线斜率 |
| `q_linear` | 这个斜率差是否为 0。校正范围是同一 `dimension` 里的全部有序对 |
| `gamma_h` | 差值轨迹的弯曲系数 |
| `q_smooth` | 两条轨迹的弯曲是否不同 |

网络和系统自己的斜率在 `group_gam_statistics.csv`，不在 q 方阵里。q 方阵只保存校正后的检验值。

### 曲线表

`parcel_curves.csv`、`group_curves.csv`、`difference_curves.csv` 的 `age` 是从样本最小年龄到最大年龄的 100 个点。`centered_effect` 是

\[
\beta_{\text{age}}(\text{age}-\bar{\text{age}})+\gamma\left(h(\text{age})-h(\bar{\text{age}})\right)
\]

脑区对的曲线可以由 `difference_statistics.csv` 里的 `beta_age` 和 `gamma_h` 按这个式子再算。

## 本轮结果里已经能直接说的部分

网络和系统的弯曲项没有过校正。`linearity_deviation_groups.png` 上七个网络的 `q_smooth` 都是 1，四个系统大约从 0.67 到 0.999，曲线都是灰色。这一层的年龄轨迹没有偏离直线。

轨迹图上的方向是 Reappraise − Look 的年龄变化，不是激活高低。默认网络相对其余网络，差值随年龄更正。视觉网络相对其余网络，差值随年龄更负。可调节情绪系统相对另外三个系统，差值随年龄更负。再评价特异减不可调节情绪的净变化为正，中段较高。共同评估与再评价特异、不可调节情绪的差值接近 0。这些方向的校正后显著性以四张 q 方阵为准。偏 \(R^2\) 大约在 0.006 这一量级，年龄项在协变量之外增加的解释量很小。
