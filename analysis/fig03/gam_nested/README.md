# 嵌套年龄模型 m1、m2、m3

这份分析用三个嵌套模型比较直线年龄和平滑年龄。它单独写在 `analysis/fig03/gam_nested/`，结果目录是 `outputs/fig03/gam_nested_volume`。已经汇报过的 `outputs/fig03/gam_supplement_volume` 不在这里重算。

要回答的是两层年龄问题。

1. 加上年龄之后，拟合是否比只有协变量更好。m2 对 m1 看一条直线，m3 对 m1 看平滑曲线。
2. 平滑曲线是否比一条直线多解释了年龄。这一步是 m3 对 m2。m3 里 `s(age)` 的汇总 p 值包含直线部分，脚本不使用它。

## 方法

### 样本和协变量

样本规则与体积 GAM 补充分析相同。

1. 有效情绪调节名单去掉二阶分析剔除的被试。
2. `participants_in_tfmri_demographics.csv` 提供 `sub_id`、`age`、`gender`、`site_id`。`gender` 是 0/1，`site_id` 保持数值。
3. 父母教育来自 `env_family_student_cbcl.tsv` 的 `父母受教育水平_父母受教育水平`，用 `被试编号` 对齐。空白、对不上、同一被试多行数值不一致的人不进入模型。去掉的人写在 `inputs/dropped_parent_edu.csv`。
4. 每个被试必须同时有再评价和被动看负性两张 `effectSize`。缺文件的编号写在 `inputs/dropped_no_volume.txt`。
5. 头动来自 `qc_res_subs_runs_fmri.csv`。`meanFD` 是该被试 EMO 和 SOC 两次 `mean_fd` 的平均。分析样本里还有人对不上时，脚本停止。

全部被试和协变量都在同一张 `inputs/design.csv`，不再另写过滤前的名单。`in_analysis` 是是否进入分析，`excluded_by` 是排除原因。性别和站点不用来删行。数值不一致、家庭表匹配不足一半，或分析样本里还有人缺 `meanFD` 时，脚本仍会停止；停止前这张表已经写好。没有单独排除原因的人写成 `sample_stop:` 加上停止原因。各原因的写法见文末的字段表。

进入模型的人数是 `in_analysis` 为 1 的行数。模型只使用这些行。带 `--meanfd-max` 时，共享缓存里的 `design.csv` 还没按头动幅度排除。`_FD05` 和 `_FD02` 仍是同一张名单，只把已经进入分析、且 `meanFD` 大于阈值的人改成 `in_analysis` 为 0，排除原因写成 `meanFD_gt_0.5` 或 `meanFD_gt_0.2`。等于阈值的人保留。旧缓存若只有进入分析的人，下次沿用体积均值时会把排除的人补回这张表，不重读 nii。

### 因变量

因变量是体积图上的再评价减被动看负性：

- `res_1stlevel/<sub>/effectSize-reg_neg_Observe_across_runs.nii.gz`
- `res_1stlevel/<sub>/effectSize-pas_neg_Observe_across_runs.nii.gz`

被试图留在自己的网格上。Schaefer 200 的 2 mm 图谱只做一次最近邻，投到效应图网格上。脑区值是该标签下有限体素的平均。

七个 Yeo 网络是所属脑区均值按体素数加权的平均。顺序是 Vis、SomMot、DorsAttn、SalVentAttn、Limbic、Cont、Default。

四个系统用 `res/res_bayes/mask_<system>.nii.gz`。系统值是掩膜内有限体素的对比均值。顺序是 `reappraisal_only`、`common_appraisal`、`nonmodifiable_emotion`、`modifiable_emotion`。

### 六个杏仁核 ROI

图谱是 Amunts 等（2005）的三个亚核，用 FSL 的 `Juelich-maxprob-thr25-2mm`。一个体素只有在某个结构概率最高、而且概率至少 25% 时才得到该标签。三个标签是 `Amygdala_LB`、`Amygdala_CM`、`Amygdala_SF`。

左右由 MNI 坐标分开：x < 0 为左，x > 0 为右。标签名里已经写了左右时，直接用该半球，不再按坐标切。六个 ROI 的顺序是 `LB_L`、`LB_R`、`CM_L`、`CM_R`、`SF_L`、`SF_R`。

六张掩膜只最近邻投到效应图网格一次，和系统掩膜使用同一张参考图。ROI 的值是掩膜内有限体素的平均。投到效应图网格后，任一 ROI 少于 5 个体素，或者一个体素落入两个 ROI，脚本停止，并留下 `inputs/amygdala_voxel_counts.csv`。

脚本按这个顺序找图谱：`--juelich-atlas` 加 `--juelich-labels`，否则 `FSLDIR` 下的 `Juelich-maxprob-thr25-2mm` 和 `Juelich.xml`，否则用 nilearn 下载 FSL 5.0 的同一份 `maxprob-thr25-2mm`。

这份分析不计算杏仁核与四个系统掩膜的重叠。

### 三个模型

每个脑区、网络、系统和杏仁核 ROI 各拟合三次，式子相同，只是 `y` 不同。拟合函数是 `mgcv::gam()`。`--method gcv` 使用 mgcv 的 `GCV.Cp`，`--method reml` 使用 `REML`。不写这个参数时是 `REML`。`--basis cr` 是三次回归样条，`--basis ps` 是 P 样条。不写时是 `cr`。`--basis-k` 是基维数，默认 3，最小是 3。系统和网络的成对 m3 斜率用同一次运行的估计方法、基和 `k`。

```text
m1: y ~ gender + site_id + parent_edu + meanFD
m2: y ~ age + gender + site_id + parent_edu + meanFD
m3: y ~ s(age, k = 3, bs = "cr") + gender + site_id + parent_edu + meanFD
```

`s(age, bs = "cr")` 是三次回归样条。默认 `k = 3`：基有 3 个函数，截距吸收常数项之后，平滑项最多约 2 个自由度，其中直线不惩罚，弯曲方向只有 1 个。`k = 4` 时弯曲方向约有 2 个。比较使用 `anova(..., test = "F")`。

`k = 4` 的三次回归样条和 pygam 的 `n_splines = 5` 不是同一套曲线。pygam 用的是 5 个 P 样条。mgcv 里靠近它的设置是 `--basis ps --basis-k 5`：5 个系数，二阶差分惩罚。平滑参数仍由 mgcv 自己搜索，不是 pygam 那 12 个 λ。截距的处理也不相同，所以 q 值不会逐格一样。对照时看 `q_m3_vs_m1` 和 pygam 的 `q_age`。

| 比较 | 含义 |
|---|---|
| m2 对 m1 | 直线年龄是否有贡献 |
| m3 对 m1 | 平滑年龄是否有贡献 |
| m3 对 m2 | 平滑是否比直线多出弯曲 |

BH 校正分四组进行：200 个脑区、7 个网络、4 个系统、6 个杏仁核。三种比较各自校正，不把不同维度或不同比较放在同一次校正里。`q < 0.05` 表示该比较在所属维度内校正后显著。

`beta_age` 是 m2 里直线年龄的系数，也就是 Reappraise − Look 随年龄的斜率。它的检验是 m2 对 m1，对应列是 `p_m2_vs_m1` 和 `q_m2_vs_m1`。

曲线在样本平均年龄处对齐到 0。m2 的曲线是 `beta_age * (age - 平均年龄)`。m3 的曲线是协变量固定在各自均值时，预测值减去平均年龄处的预测值。图里的相交只表示这个零点。

## 运行

在工程根目录执行。需要更新 `analysis/fig03/gam_supplement/prepare_inputs.py`，并放上 `analysis/fig03/gam_nested/`。

```bash
python -m analysis.fig03.gam_nested.run_nested --input "/public/home/dingrui/fmri_analysis/data/beh/participants_in_tfmri_demographics.csv" --subjects "/public/home/dingrui/fmri_analysis/data/beh/sub_list_valid_4_ER.txt" --excluded "/public/home/dingrui/fmri_analysis/data/fmri/subs_removed_from_2ndlevel.txt" --parent-edu "/public/home/dingrui/fmri_analysis/data/beh/env_family_student_cbcl.tsv" --meanfd "/public/home/dingrui/fmri_analysis/data/fmri/qc_res_subs_runs_fmri.csv" --fig3-volume-dir "/public/home/dingrui/fmri_analysis/data/fmri/res_1stlevel" --fig3-volume-atlas "resources/atlases/schaefer2018/Schaefer2018_200Parcels_7Networks_order_FSLMNI152_2mm.nii.gz" --schaefer-atlas "resources/atlases/schaefer2018/Schaefer2018_200Parcels_7Networks_order.dlabel.nii" --system-mask-dir "/public/home/dingrui/fmri_analysis/res/res_bayes" --output "outputs/fig03/gam_nested_volume"
```

上面这条不按头动幅度删人，结果仍在 `outputs/fig03/gam_nested_volume`。再跑一版不过滤的 200 脑区时，不要用这个 `--output`。

## 200 脑区，不按头动幅度过滤

范围与下面的头动版本相同：Schaefer 200 × Yeo 7，7 个网络、4 个系统、6 个杏仁核，四个协变量都保留。不写 `--method`、`--basis`、`--basis-k`，分别是 `REML`、`cr` 和 `k = 3`。不写 `--meanfd-max`，不按 0.5 或 0.2 删人。`meanFD` 仍是协变量。

结果直接写在 `outputs/fig03/gam_nested_volume_200`。不写入 `gam_nested_volume`，也不写入 `gam_nested_volume_200_FD05`、`gam_nested_volume_200_FD02` 和 `gam_nested_volume_200_shared_inputs`。

```bash
python -m analysis.fig03.gam_nested.run_nested --input "/public/home/dingrui/fmri_analysis/data/beh/participants_in_tfmri_demographics.csv" --subjects "/public/home/dingrui/fmri_analysis/data/beh/sub_list_valid_4_ER.txt" --excluded "/public/home/dingrui/fmri_analysis/data/fmri/subs_removed_from_2ndlevel.txt" --parent-edu "/public/home/dingrui/fmri_analysis/data/beh/env_family_student_cbcl.tsv" --meanfd "/public/home/dingrui/fmri_analysis/data/fmri/qc_res_subs_runs_fmri.csv" --fig3-volume-dir "/public/home/dingrui/fmri_analysis/data/fmri/res_1stlevel" --fig3-volume-atlas "resources/atlases/schaefer2018/Schaefer2018_200Parcels_7Networks_order_FSLMNI152_2mm.nii.gz" --schaefer-atlas "resources/atlases/schaefer2018/Schaefer2018_200Parcels_7Networks_order.dlabel.nii" --system-mask-dir "/public/home/dingrui/fmri_analysis/res/res_bayes" --output "outputs/fig03/gam_nested_volume_200"
```

## 200 脑区，头动 0.5 和 0.2

这一版是 Schaefer 200 × Yeo 7，同时拟合 7 个网络、4 个系统和 6 个杏仁核。四个协变量都保留。不写 `--method`、`--basis`、`--basis-k` 时，分别是 `REML`、`cr` 和 `k = 3`。`meanFD` 大于 0.5 和大于 0.2 的被试分别去掉。

体积图只读一次，放在 `outputs/fig03/gam_nested_volume_200_shared_inputs`。结果在 `outputs/fig03/gam_nested_volume_200_FD05` 和 `outputs/fig03/gam_nested_volume_200_FD02`。不写入 `gam_nested_volume`，也不写入 400 脑区的目录。

```bash
python -m analysis.fig03.gam_nested.run_nested --meanfd-max 0.5 0.2 --input "/public/home/dingrui/fmri_analysis/data/beh/participants_in_tfmri_demographics.csv" --subjects "/public/home/dingrui/fmri_analysis/data/beh/sub_list_valid_4_ER.txt" --excluded "/public/home/dingrui/fmri_analysis/data/fmri/subs_removed_from_2ndlevel.txt" --parent-edu "/public/home/dingrui/fmri_analysis/data/beh/env_family_student_cbcl.tsv" --meanfd "/public/home/dingrui/fmri_analysis/data/fmri/qc_res_subs_runs_fmri.csv" --fig3-volume-dir "/public/home/dingrui/fmri_analysis/data/fmri/res_1stlevel" --fig3-volume-atlas "resources/atlases/schaefer2018/Schaefer2018_200Parcels_7Networks_order_FSLMNI152_2mm.nii.gz" --schaefer-atlas "resources/atlases/schaefer2018/Schaefer2018_200Parcels_7Networks_order.dlabel.nii" --system-mask-dir "/public/home/dingrui/fmri_analysis/res/res_bayes" --output "outputs/fig03/gam_nested_volume_200"
```

## 400 脑区，Yeo 7，不按头动幅度过滤

范围与下面的头动版本相同：Schaefer 400 × Yeo 7，7 个网络、4 个系统、6 个杏仁核，四个协变量都保留。不写 `--method`、`--basis`、`--basis-k`，分别是 `REML`、`cr` 和 `k = 3`。不写 `--meanfd-max`，不按 0.5 或 0.2 删人。`meanFD` 仍是协变量。

结果直接写在 `outputs/fig03/gam_nested_volume_400`。不写入 `gam_nested_volume_400_FD05`、`gam_nested_volume_400_FD02`、`gam_nested_volume_400_shared_inputs`，也不写入 `gam_nested_volume_400_gcv`、`gam_nested_volume_400_gcv_ps5`、`gam_nested_volume_400_reml` 和 `gam_nested_volume_400_17net`。

```bash
python -m analysis.fig03.gam_nested.run_nested --input "/public/home/dingrui/fmri_analysis/data/beh/participants_in_tfmri_demographics.csv" --subjects "/public/home/dingrui/fmri_analysis/data/beh/sub_list_valid_4_ER.txt" --excluded "/public/home/dingrui/fmri_analysis/data/fmri/subs_removed_from_2ndlevel.txt" --parent-edu "/public/home/dingrui/fmri_analysis/data/beh/env_family_student_cbcl.tsv" --meanfd "/public/home/dingrui/fmri_analysis/data/fmri/qc_res_subs_runs_fmri.csv" --fig3-volume-dir "/public/home/dingrui/fmri_analysis/data/fmri/res_1stlevel" --fig3-volume-atlas "resources/atlases/schaefer2018/Schaefer2018_400Parcels_7Networks_order_FSLMNI152_2mm.nii.gz" --schaefer-atlas "resources/atlases/schaefer2018/Schaefer2018_400Parcels_7Networks_order.dlabel.nii" --system-mask-dir "/public/home/dingrui/fmri_analysis/res/res_bayes" --output "outputs/fig03/gam_nested_volume_400"
```

## 400 脑区，头动 0.5 和 0.2

这一版仍是 m1、m2、m3，四个协变量都保留，包括 `parent_edu`。分区换成 Schaefer 400 × Yeo 7，网络、4 个系统和 6 个杏仁核照旧。`meanFD` 大于阈值的被试被去掉。0.5 和 0.2 各写一个目录，不覆盖 `outputs/fig03/gam_nested_volume`。

体积图只读一次，放在 `outputs/fig03/gam_nested_volume_400_shared_inputs`。两个结果目录是 `outputs/fig03/gam_nested_volume_400_FD05` 和 `outputs/fig03/gam_nested_volume_400_FD02`。表面标签是 `Schaefer2018_400Parcels_7Networks_order.dlabel.nii`，编号和名称与同目录的 `order_info.txt` 一致。

图都是 svg，文字保持为文字。`figures/significant_m3_vs_m1_parcels.svg` 只画 q(m3 vs m1) < 0.05 的脑区，每个脑区一种颜色。`figures/age_parent_edu.svg` 是 design 表里年龄和 `parent_edu` 的散点。

```bash
python -m analysis.fig03.gam_nested.run_nested --meanfd-max 0.5 0.2 --input "/public/home/dingrui/fmri_analysis/data/beh/participants_in_tfmri_demographics.csv" --subjects "/public/home/dingrui/fmri_analysis/data/beh/sub_list_valid_4_ER.txt" --excluded "/public/home/dingrui/fmri_analysis/data/fmri/subs_removed_from_2ndlevel.txt" --parent-edu "/public/home/dingrui/fmri_analysis/data/beh/env_family_student_cbcl.tsv" --meanfd "/public/home/dingrui/fmri_analysis/data/fmri/qc_res_subs_runs_fmri.csv" --fig3-volume-dir "/public/home/dingrui/fmri_analysis/data/fmri/res_1stlevel" --fig3-volume-atlas "resources/atlases/schaefer2018/Schaefer2018_400Parcels_7Networks_order_FSLMNI152_2mm.nii.gz" --schaefer-atlas "resources/atlases/schaefer2018/Schaefer2018_400Parcels_7Networks_order.dlabel.nii" --system-mask-dir "/public/home/dingrui/fmri_analysis/res/res_bayes" --output "outputs/fig03/gam_nested_volume_400"
```

服务器上已有 FSL 时，上面这条 Yeo 7 命令读取 `$FSLDIR/data/atlases/Juelich/Juelich-maxprob-thr25-2mm.nii.gz`。图谱在别的路径时加上：

```bash
--juelich-atlas "/path/to/Juelich-maxprob-thr25-2mm.nii.gz" --juelich-labels "/path/to/Juelich.xml"
```

## 400 × Yeo 7，平滑估计用 GCV

范围与上面的 Yeo 7 命令相同：400 个脑区、7 个网络、4 个系统、6 个杏仁核，四个协变量，头动 0.5 和 0.2。差别是 `--method gcv`，mgcv 内部用 `GCV.Cp`。输出是 `outputs/fig03/gam_nested_volume_400_gcv_FD05` 和 `outputs/fig03/gam_nested_volume_400_gcv_FD02`，体积均值缓存在 `outputs/fig03/gam_nested_volume_400_gcv_shared_inputs`。不写入 `gam_nested_volume_400_FD05`、`gam_nested_volume_400_FD02` 和 `gam_nested_volume`。

要在另一个新目录里用 REML 再拟合同一范围，把 `--method gcv` 换成 `--method reml`，`--output` 换成 `outputs/fig03/gam_nested_volume_400_reml`。

```bash
python -m analysis.fig03.gam_nested.run_nested --method gcv --meanfd-max 0.5 0.2 --input "/public/home/dingrui/fmri_analysis/data/beh/participants_in_tfmri_demographics.csv" --subjects "/public/home/dingrui/fmri_analysis/data/beh/sub_list_valid_4_ER.txt" --excluded "/public/home/dingrui/fmri_analysis/data/fmri/subs_removed_from_2ndlevel.txt" --parent-edu "/public/home/dingrui/fmri_analysis/data/beh/env_family_student_cbcl.tsv" --meanfd "/public/home/dingrui/fmri_analysis/data/fmri/qc_res_subs_runs_fmri.csv" --fig3-volume-dir "/public/home/dingrui/fmri_analysis/data/fmri/res_1stlevel" --fig3-volume-atlas "resources/atlases/schaefer2018/Schaefer2018_400Parcels_7Networks_order_FSLMNI152_2mm.nii.gz" --schaefer-atlas "resources/atlases/schaefer2018/Schaefer2018_400Parcels_7Networks_order.dlabel.nii" --system-mask-dir "/public/home/dingrui/fmri_analysis/res/res_bayes" --output "outputs/fig03/gam_nested_volume_400_gcv"
```

## 400 × Yeo 7，P 样条 k = 5，GCV

年龄项是 `s(age, k = 5, bs = "ps")`，平滑参数是 `GCV.Cp`。范围仍是 400 个脑区、7 个网络、4 个系统、6 个杏仁核，四个协变量，头动 0.5 和 0.2。输出是 `outputs/fig03/gam_nested_volume_400_gcv_ps5_FD05` 和 `outputs/fig03/gam_nested_volume_400_gcv_ps5_FD02`，体积均值缓存在 `outputs/fig03/gam_nested_volume_400_gcv_ps5_shared_inputs`。不写入 `gam_nested_volume_400`、`gam_nested_volume_400_gcv` 和 `gam_nested_volume_400_gcv_k4`。

```bash
python -m analysis.fig03.gam_nested.run_nested --method gcv --basis ps --basis-k 5 --meanfd-max 0.5 0.2 --input "/public/home/dingrui/fmri_analysis/data/beh/participants_in_tfmri_demographics.csv" --subjects "/public/home/dingrui/fmri_analysis/data/beh/sub_list_valid_4_ER.txt" --excluded "/public/home/dingrui/fmri_analysis/data/fmri/subs_removed_from_2ndlevel.txt" --parent-edu "/public/home/dingrui/fmri_analysis/data/beh/env_family_student_cbcl.tsv" --meanfd "/public/home/dingrui/fmri_analysis/data/fmri/qc_res_subs_runs_fmri.csv" --fig3-volume-dir "/public/home/dingrui/fmri_analysis/data/fmri/res_1stlevel" --fig3-volume-atlas "resources/atlases/schaefer2018/Schaefer2018_400Parcels_7Networks_order_FSLMNI152_2mm.nii.gz" --schaefer-atlas "resources/atlases/schaefer2018/Schaefer2018_400Parcels_7Networks_order.dlabel.nii" --system-mask-dir "/public/home/dingrui/fmri_analysis/res/res_bayes" --output "outputs/fig03/gam_nested_volume_400_gcv_ps5"
```

## Kong 2022 的 17 网络，400 脑区

这一版只拟合 400 个脑区。17 个网络名只用来标脑区，不对网络做平均，也不拟合 4 个系统和 6 个杏仁核。模型仍是 m1、m2、m3，`method = "REML"`，四个协变量都保留。`meanFD` 大于 0.5 和大于 0.2 的被试分别去掉。输出是 `outputs/fig03/gam_nested_volume_400_17net_FD05` 和 `outputs/fig03/gam_nested_volume_400_17net_FD02`，体积均值缓存在 `outputs/fig03/gam_nested_volume_400_17net_shared_inputs`。不覆盖 Yeo 7 的 `gam_nested_volume_400`，也不覆盖 `gam_pygam_volume_400`。

脑区名称来自 `Schaefer2018_400Parcels_Kong2022_17Networks_order_info.txt`。皮层图用同目录的 `Schaefer2018_400Parcels_Kong2022_17Networks_order.dlabel.nii`。显著脑区曲线和年龄与 `parent_edu` 的散点仍是 svg。这条命令不做斜率比较。

```bash
python -m analysis.fig03.gam_nested.run_nested --parcels-only --meanfd-max 0.5 0.2 --input "/public/home/dingrui/fmri_analysis/data/beh/participants_in_tfmri_demographics.csv" --subjects "/public/home/dingrui/fmri_analysis/data/beh/sub_list_valid_4_ER.txt" --excluded "/public/home/dingrui/fmri_analysis/data/fmri/subs_removed_from_2ndlevel.txt" --parent-edu "/public/home/dingrui/fmri_analysis/data/beh/env_family_student_cbcl.tsv" --meanfd "/public/home/dingrui/fmri_analysis/data/fmri/qc_res_subs_runs_fmri.csv" --fig3-volume-dir "/public/home/dingrui/fmri_analysis/data/fmri/res_1stlevel" --fig3-volume-atlas "resources/atlases/schaefer2018/Schaefer2018_400Parcels_Kong2022_17Networks_order_FSLMNI152_2mm.nii.gz" --parcel-info "resources/atlases/schaefer2018/Schaefer2018_400Parcels_Kong2022_17Networks_order_info.txt" --schaefer-atlas "resources/atlases/schaefer2018/Schaefer2018_400Parcels_Kong2022_17Networks_order.dlabel.nii" --output "outputs/fig03/gam_nested_volume_400_17net"
```

| 文件 | 作用 |
|---|---|
| `amygdala.py` | 从 Juelich 最大概率图生成 6 个 ROI |
| `run_nested.py` | 准备设计表、体积均值，调用 R 并作图 |
| `mgcv_nested_models.R` | 拟合 m1、m2、m3，做嵌套 F 检验和 BH 校正 |
| `plot_nested.py` | 画网络、系统和杏仁核的 m2、m3 曲线 |
| `../gam_supplement/prepare_inputs.py` | 读取效应图，计算脑区、网络和系统均值 |

## 输出文件

路径相对于 `outputs/fig03/gam_nested_volume`。

### 输入

| 文件 | 内容 |
|---|---|
| `inputs/design.csv` | 全部被试、四个协变量、是否进入分析和排除原因 |
| `inputs/parcel_response.csv` | 每个被试、每个脑区的 Reappraise − Look 均值 |
| `inputs/parcel_info.csv` | 脑区编号、名称和 Yeo 网络 |
| `inputs/group_response.csv` | 网络和系统的长表均值 |
| `inputs/yeo_order.csv`、`inputs/system_order.csv` | 网络和系统的顺序 |
| `inputs/amygdala_response.csv` | 6 个杏仁核 ROI 的被试均值 |
| `inputs/amygdala_order.csv` | 6 个 ROI 的顺序 |
| `inputs/amygdala_voxel_counts.csv` | 每个 ROI 投到效应图网格后的体素数 |
| `inputs/amygdala_masks/` | 投到效应图网格后的 6 张二值掩膜 |
| `inputs/amygdala_source.txt` | 图谱路径、左右划分和对比定义 |
| `inputs/dropped_parent_edu.csv` | 因父母教育去掉的被试。完整原因仍以 `design.csv` 的 `excluded_by` 为准 |
| `inputs/dropped_no_volume.txt` | 缺两张效应图的被试 |
| `inputs/dropped_meanfd.csv` | 这一头动阈值额外去掉的被试。只出现在 `_FD05`、`_FD02` 目录 |
| `inputs/parent_edu_source.txt`、`inputs/meanfd_source.txt` | 两个协变量的来源 |

### 结果

| 文件 | 内容 |
|---|---|
| `statistics/nested_model_comparison.csv` | 全部脑区、网络、系统和杏仁核的三次拟合与三种比较 |
| `statistics/nested_curves.csv` | m2 和 m3 的中心化年龄曲线，每个对象 100 个年龄点 |
| `statistics/model_reference.csv` | 平均年龄、本次的 `k`、`cr` 或 `ps`、`GCV.Cp` 或 `REML`，以及比较和校正的范围 |
| `figures/yeo_m2_m3_curves.svg` | 7 个网络。灰色虚线是 m2，棕色实线是 m3 |
| `figures/system_m2_m3_curves.svg` | 4 个系统 |
| `figures/amygdala_m2_m3_curves.svg` | 6 个杏仁核 ROI |
| `run_metadata.json` | 模型公式、校正范围和样本量 |

200 个脑区的结果在 `nested_model_comparison.csv` 里，`dimension` 为 `parcel`。脚本不为 200 个脑区画曲线格子。

给老师看时，主表用 `statistics/nested_model_comparison.csv`。偏离直线看 `q_m3_vs_m2`。直线年龄看 `q_m2_vs_m1` 和 `beta_age`。有没有年龄看 `q_m3_vs_m1`。网络、系统和杏仁核的曲线看 `figures/` 里的三张图。

## 系统和网络的 m3 斜率比较

每个系统、每个网络自己的 m3 拟合不比较两组之间的年龄斜率。比较在已有的 `inputs/group_response.csv` 上另做，不重读体积图。

每一对有序组合先做被试水平的差值，行减列，再拟合和单区域相同的 m3。斜率是这条 m3 年龄曲线从最小年龄到最大年龄的平均变化，单位是每年。协变量固定在各自的均值上，所以截距和协变量不进入这个差。p 值是该斜率是否为 0 的 Wald 检验，协方差用惩罚后的贝叶斯协方差。BH 校正在 42 个网络有序对内部做一次，在 12 个系统有序对内部做一次。A 减 B 和 B 减 A 方向相反、p 值相同，两格都进入校正。

在已经跑完嵌套模型的目录上执行：

```bash
python -m analysis.fig03.gam_nested.run_nested --step slopes --output "outputs/fig03/gam_nested_volume"
```

| 文件 | 内容 |
|---|---|
| `statistics/m3_pair_slopes.csv` | 每一对的斜率、标准误、p 和 q |
| `statistics/m3_pair_curves.csv` | 每一对差值的 m3 中心化年龄曲线 |
| `statistics/yeo_m3_slope_matrix.csv` | 7×7 斜率。行减列，对角线为空 |
| `statistics/yeo_m3_p_slope_matrix.csv`、`yeo_m3_q_slope_matrix.csv` | 对应的 p 和 q |
| `statistics/system_m3_slope_matrix.csv` | 4×4 斜率 |
| `statistics/system_m3_p_slope_matrix.csv`、`system_m3_q_slope_matrix.csv` | 对应的 p 和 q |
| `figures/yeo_m3_slope_matrix.svg`、`figures/system_m3_slope_matrix.svg` | 差值曲线。颜色是每年的 m3 斜率 |

`m3_pair_slopes.csv` 的 `slope` 为正，表示行的 m3 年龄曲线比列更向正方向变化。`q_slope < 0.05` 表示这对斜率差在该维度内校正后不为 0。图上的曲线是整条 `d(age)`，不只是一条直线；检验的数是这条曲线两端的平均斜率。

## 重要字段

### `inputs/design.csv`

每一行是有效名单去掉二阶剔除之后的一个被试。二阶剔除的人不在这张表里。协变量能对齐就写上，对不齐就空着。

| 列 | 含义 |
|---|---|
| `sub_id` | 成像被试编号 |
| `age` | 年龄 |
| `gender` | 0/1。不用于删行 |
| `site_id` | 数值型站点。不用于删行 |
| `parent_edu` | `父母受教育水平_父母受教育水平`。没有单一数值时为空 |
| `meanFD` | EMO 与 SOC 的 `mean_fd` 平均，单位毫米。没有数值时为空 |
| `in_analysis` | 是否进入分析。1 进入，0 不进入 |
| `excluded_by` | 排除原因。进入分析的人为空 |

`excluded_by` 按流水线记下第一条原因：`parent_edu_not_in_table`、`parent_edu_missing_value`、`parent_edu_conflicting_values`、`no_effect_volume`、`meanFD_missing`。家庭表或头动表本身读不到时，原因是 `parent_edu_file_missing`、`parent_edu_columns_missing`、`parent_edu_unreadable`、`meanFD_file_missing` 或 `meanFD_unreadable`。

整批样本因此停止、而这个人自己没有上面的原因时，写成 `sample_stop:` 加停止原因，例如 `sample_stop:parent_edu_conflicting_values`、`sample_stop:parent_edu_match_rate`、`sample_stop:meanFD_missing`。本人的原因仍优先，例如缺头动的那一行是 `meanFD_missing`，其余本来可以进入的人是 `sample_stop:meanFD_missing`。

头动幅度阈值写在同一个 `excluded_by` 里，不另加列。已经进入分析且 `meanFD` 大于阈值时改为 `meanFD_gt_0.5` 或 `meanFD_gt_0.2`，`in_analysis` 改为 0。更早被排除的人保持原来的原因。

### `statistics/nested_model_comparison.csv`

每一行是一个脑区、网络、系统或杏仁核 ROI。

| 列 | 含义 |
|---|---|
| `dimension` | `parcel`、`yeo`、`system` 或 `amygdala` |
| `name` | 脑区编号、网络名、系统名或 ROI 名 |
| `n` | 该模型的被试数 |
| `beta_age` | m2 中直线年龄的系数。正值表示 Reappraise − Look 随年龄升高 |
| `se_age` | `beta_age` 的标准误 |
| `edf_m1`、`edf_m2`、`edf_m3` | 三个模型的有效自由度 |
| `deviance_m1`、`deviance_m2`、`deviance_m3` | 三个模型的残差偏差。越小拟合越好 |
| `F_m2_vs_m1`、`df_m2_vs_m1`、`p_m2_vs_m1` | m2 对 m1 的嵌套 F、自由度差和未校正 p |
| `q_m2_vs_m1` | `p_m2_vs_m1` 在该 `dimension` 内的 BH 校正 |
| `F_m3_vs_m1`、`df_m3_vs_m1`、`p_m3_vs_m1` | m3 对 m1，即平滑年龄是否有贡献 |
| `q_m3_vs_m1` | `p_m3_vs_m1` 在该维度内的 BH 校正 |
| `F_m3_vs_m2`、`df_m3_vs_m2`、`p_m3_vs_m2` | m3 对 m2，即是否偏离直线 |
| `q_m3_vs_m2` | `p_m3_vs_m2` 在该维度内的 BH 校正。小于 0.05 表示偏离直线 |

`q_m2_vs_m1` 只在同一种比较、同一个维度里校正。脑区的 q 不进入网络、系统或杏仁核的校正。

### `statistics/nested_curves.csv`

| 列 | 含义 |
|---|---|
| `dimension`、`name` | 与比较表相同 |
| `model` | `m2` 或 `m3` |
| `age` | 从样本最小年龄到最大年龄的 100 个点 |
| `centered_effect` | 该年龄点的中心化年龄效应。平均年龄处为 0 |

### `statistics/m3_pair_slopes.csv`

每一行是一对有序比较。因变量是该被试的 `group_a` 减 `group_b`，模型是 m3。网络 42 行，系统 12 行。

| 列 | 含义 |
|---|---|
| `dimension` | `yeo` 或 `system` |
| `group_a` | 被减的一组。矩阵里的行 |
| `group_b` | 减去的一组。矩阵里的列 |
| `slope` | m3 年龄曲线两端之差除以年龄跨度，单位是每年。正值表示 `group_a` 的 Reappraise − Look 比 `group_b` 更向正方向变化 |
| `se_slope` | `slope` 的标准误，来自惩罚后协方差上的 Wald 对比 |
| `z_slope` | `slope / se_slope` |
| `p_slope` | `slope` 是否为 0 的双侧 Wald p，未校正 |
| `q_slope` | `p_slope` 在该 `dimension` 内的 BH 校正。小于 0.05 表示这对斜率差校正后不为 0 |

`group_a` 与 `group_b` 对调后，`slope` 和 `z_slope` 变号，`se_slope` 与 `p_slope` 不变。两行都进入 BH。

### `statistics/m3_pair_curves.csv`

每一对有 100 行，对应图中的一条差值曲线。

| 列 | 含义 |
|---|---|
| `dimension`、`group_a`、`group_b` | 与 `m3_pair_slopes.csv` 相同 |
| `age` | 从样本最小年龄到最大年龄的 100 个点 |
| `centered_effect` | 该年龄点的 m3 预测减去平均年龄处的预测。协变量固定在均值。平均年龄处为 0 |

`slope` 是这条曲线两端的平均变化，不是曲线上每一个点的局部斜率。

### 斜率矩阵

`yeo_m3_slope_matrix.csv`、`yeo_m3_p_slope_matrix.csv`、`yeo_m3_q_slope_matrix.csv` 是 7×7。`system_m3_slope_matrix.csv`、`system_m3_p_slope_matrix.csv`、`system_m3_q_slope_matrix.csv` 是 4×4。数值与长表相同。

| 位置 | 含义 |
|---|---|
| 第一列 `group` | 行名，即 `group_a` |
| 其余列名 | 列名，即 `group_b`，顺序与 `inputs/yeo_order.csv` 或 `inputs/system_order.csv` 相同 |
| 单元格 | 行减列的 `slope`、`p_slope` 或 `q_slope` |
| 对角线 | 空。一组不与自己比较 |

### `inputs/amygdala_voxel_counts.csv`

| 列 | 含义 |
|---|---|
| `roi` | `LB_L`、`LB_R`、`CM_L`、`CM_R`、`SF_L`、`SF_R` 之一 |
| `n_voxels` | 该 ROI 在效应图网格上的体素数 |
