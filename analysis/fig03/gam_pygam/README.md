# pygam 年龄模型

全模型是截距 + `s(age)` + `gender` + `site_id` + `parent_edu` + `meanFD`。`s(age)` 是 5 个 B 样条，平滑参数在 `1e-5` 到 `1e5` 的 12 个值上用 GCV 选择，只惩罚年龄项。小模型是截距加上同样的四个协变量。`p_age` 是这两个模型的嵌套 F 检验，问的是年龄曲线是否有贡献。它对应 mgcv 的 m3 对 m1，不检验相对直线的弯曲。

BH 在脑区、Yeo 7、4 个系统、6 个杏仁核内部各做一次。

## 正式运行清单

在工程根目录、能导入 pygam 的环境里按顺序执行。四条都调用 `python -m analysis.fig03.gam_pygam.run_pygam`，都不加 `--parcels-only`。每条同时拟合脑区、Yeo 7、4 个系统和 6 个杏仁核。年龄项是 5 个 B 样条，平滑参数在 `1e-5` 到 `1e5` 的 12 个值上用 GCV 选择。四个协变量都保留。

头动 0.5 和 0.2 写在同一条命令里。不过滤要单独一条，因为不写 `--meanfd-max` 时结果直接进 `--output`。这四条不写入 `gam_pygam_volume` 和 Kong 2022 的 `gam_pygam_volume_400`。不过滤不会沿用头动目录里的体积缓存，体积图会再读一遍。

| 顺序 | 范围 | 头动 | 结果目录 |
|---|---|---|---|
| 1 | Schaefer 200 × Yeo 7 | 不过滤 | `outputs/fig03/gam_pygam_volume_200` |
| 2 | Schaefer 200 × Yeo 7 | 0.5 和 0.2 | `outputs/fig03/gam_pygam_volume_200_FD05`、`_FD02`，缓存 `_shared_inputs` |
| 3 | Schaefer 400 × Yeo 7 | 不过滤 | `outputs/fig03/gam_pygam_volume_400_7net` |
| 4 | Schaefer 400 × Yeo 7 | 0.5 和 0.2 | `outputs/fig03/gam_pygam_volume_400_7net_FD05`、`_FD02`，缓存 `_shared_inputs` |

```bash
python -m analysis.fig03.gam_pygam.run_pygam --input "/public/home/dingrui/fmri_analysis/data/beh/participants_in_tfmri_demographics.csv" --subjects "/public/home/dingrui/fmri_analysis/data/beh/sub_list_valid_4_ER.txt" --excluded "/public/home/dingrui/fmri_analysis/data/fmri/subs_removed_from_2ndlevel.txt" --parent-edu "/public/home/dingrui/fmri_analysis/data/beh/env_family_student_cbcl.tsv" --meanfd "/public/home/dingrui/fmri_analysis/data/fmri/qc_res_subs_runs_fmri.csv" --fig3-volume-dir "/public/home/dingrui/fmri_analysis/data/fmri/res_1stlevel" --fig3-volume-atlas "resources/atlases/schaefer2018/Schaefer2018_200Parcels_7Networks_order_FSLMNI152_2mm.nii.gz" --schaefer-atlas "resources/atlases/schaefer2018/Schaefer2018_200Parcels_7Networks_order.dlabel.nii" --system-mask-dir "/public/home/dingrui/fmri_analysis/res/res_bayes" --output "outputs/fig03/gam_pygam_volume_200"
```

```bash
python -m analysis.fig03.gam_pygam.run_pygam --meanfd-max 0.5 0.2 --input "/public/home/dingrui/fmri_analysis/data/beh/participants_in_tfmri_demographics.csv" --subjects "/public/home/dingrui/fmri_analysis/data/beh/sub_list_valid_4_ER.txt" --excluded "/public/home/dingrui/fmri_analysis/data/fmri/subs_removed_from_2ndlevel.txt" --parent-edu "/public/home/dingrui/fmri_analysis/data/beh/env_family_student_cbcl.tsv" --meanfd "/public/home/dingrui/fmri_analysis/data/fmri/qc_res_subs_runs_fmri.csv" --fig3-volume-dir "/public/home/dingrui/fmri_analysis/data/fmri/res_1stlevel" --fig3-volume-atlas "resources/atlases/schaefer2018/Schaefer2018_200Parcels_7Networks_order_FSLMNI152_2mm.nii.gz" --schaefer-atlas "resources/atlases/schaefer2018/Schaefer2018_200Parcels_7Networks_order.dlabel.nii" --system-mask-dir "/public/home/dingrui/fmri_analysis/res/res_bayes" --output "outputs/fig03/gam_pygam_volume_200"
```

```bash
python -m analysis.fig03.gam_pygam.run_pygam --input "/public/home/dingrui/fmri_analysis/data/beh/participants_in_tfmri_demographics.csv" --subjects "/public/home/dingrui/fmri_analysis/data/beh/sub_list_valid_4_ER.txt" --excluded "/public/home/dingrui/fmri_analysis/data/fmri/subs_removed_from_2ndlevel.txt" --parent-edu "/public/home/dingrui/fmri_analysis/data/beh/env_family_student_cbcl.tsv" --meanfd "/public/home/dingrui/fmri_analysis/data/fmri/qc_res_subs_runs_fmri.csv" --fig3-volume-dir "/public/home/dingrui/fmri_analysis/data/fmri/res_1stlevel" --fig3-volume-atlas "resources/atlases/schaefer2018/Schaefer2018_400Parcels_7Networks_order_FSLMNI152_2mm.nii.gz" --schaefer-atlas "resources/atlases/schaefer2018/Schaefer2018_400Parcels_7Networks_order.dlabel.nii" --system-mask-dir "/public/home/dingrui/fmri_analysis/res/res_bayes" --output "outputs/fig03/gam_pygam_volume_400_7net"
```

```bash
python -m analysis.fig03.gam_pygam.run_pygam --meanfd-max 0.5 0.2 --input "/public/home/dingrui/fmri_analysis/data/beh/participants_in_tfmri_demographics.csv" --subjects "/public/home/dingrui/fmri_analysis/data/beh/sub_list_valid_4_ER.txt" --excluded "/public/home/dingrui/fmri_analysis/data/fmri/subs_removed_from_2ndlevel.txt" --parent-edu "/public/home/dingrui/fmri_analysis/data/beh/env_family_student_cbcl.tsv" --meanfd "/public/home/dingrui/fmri_analysis/data/fmri/qc_res_subs_runs_fmri.csv" --fig3-volume-dir "/public/home/dingrui/fmri_analysis/data/fmri/res_1stlevel" --fig3-volume-atlas "resources/atlases/schaefer2018/Schaefer2018_400Parcels_7Networks_order_FSLMNI152_2mm.nii.gz" --schaefer-atlas "resources/atlases/schaefer2018/Schaefer2018_400Parcels_7Networks_order.dlabel.nii" --system-mask-dir "/public/home/dingrui/fmri_analysis/res/res_bayes" --output "outputs/fig03/gam_pygam_volume_400_7net"
```

## 432 和 232，Tian S2，只过滤头动 0.5

这两条在 Schaefer 皮层后面接上 Tian Subcortex S2 的 32 个皮下区。432 是 Schaefer 400 加上编号 401–432。232 是 Schaefer 200 加上编号 201–232。Yeo 7 只由皮层 parcel 合成。4 个系统仍用贝叶斯 mask。BH 分三族：全部脑区一起校正，7 个网络一起校正，4 个系统一起校正。

杏仁核用 Tian S2 3T 的正式名称，不读 NIfTI 头扩展。1–16 是右侧，17–32 是左侧。四个杏仁核是 3 `lAMY-rh`、4 `mAMY-rh`、19 `lAMY-lh`、20 `mAMY-lh`，留在脑区族里。fig06 的 `TIAN_S2_LABELS` 把它们收成 V_3 和 V_19 两个 AMY，这里不用那张表。432 版把 Schaefer 编号 156 和 157 标成边缘显著。232 版没有这两个编号。图谱旁边若有完整的 32 行名称表，用那张表代替内置顺序。

只去掉 `meanFD > 0.5` 的被试。体积缓存分别在 `gam_pygam_volume_432_shared_inputs` 和 `gam_pygam_volume_232_shared_inputs`。拟合结果在 `gam_pygam_volume_432_FD05` 和 `gam_pygam_volume_232_FD05`。不写入 `gam_pygam_volume`、`gam_pygam_volume_200` 和 `gam_pygam_volume_400_7net`。

全脑 `maps/partial_r2.nii.gz` 填 fig3 同一个 partial R²。`maps/significance_mask.nii.gz` 里 1 是 `q_age < 0.05`，2 是 156 和 157。这两张 volume 再转到 fsLR 32k，写出同名的左右 `func.gii` 和 `dscalar.nii`。代码不到 Workbench 里出图。

随龄曲线是一个脑区一张 SVG，写在 `figures/age_curves/`。纳入 `q_age < 0.05` 的脑区、432 版的 156 和 157、以及 4 个杏仁核。散点是去掉截距和四个协变量后的残差。曲线和 95% 置信区间用同一套年龄效应。不画 0 线，去掉上边框和右边框，文字可编辑。7 个网络和 4 个系统只写统计，不画曲线。

432：

```bash
python -m analysis.fig03.gam_pygam.run_pygam --meanfd-max 0.5 --tian-atlas "/public/home/dingrui/tools/masks_atlas/Tian_Subcortex_S2_3T_2009cAsym.nii.gz" --input "/public/home/dingrui/fmri_analysis/data/beh/participants_in_tfmri_demographics.csv" --subjects "/public/home/dingrui/fmri_analysis/data/beh/sub_list_valid_4_ER.txt" --excluded "/public/home/dingrui/fmri_analysis/data/fmri/subs_removed_from_2ndlevel.txt" --parent-edu "/public/home/dingrui/fmri_analysis/data/beh/env_family_student_cbcl.tsv" --meanfd "/public/home/dingrui/fmri_analysis/data/fmri/qc_res_subs_runs_fmri.csv" --fig3-volume-dir "/public/home/dingrui/fmri_analysis/data/fmri/res_1stlevel" --fig3-volume-atlas "resources/atlases/schaefer2018/Schaefer2018_400Parcels_7Networks_order_FSLMNI152_2mm.nii.gz" --schaefer-atlas "resources/atlases/schaefer2018/Schaefer2018_400Parcels_7Networks_order.dlabel.nii" --system-mask-dir "/public/home/dingrui/fmri_analysis/res/res_bayes" --output "outputs/fig03/gam_pygam_volume_432"
```

432 的拟合如果已经写完统计表，只是 fsLR 的 dscalar 没写成，用 `--step maps` 重出脑图。这一步不再拟合，也不重读被试体积。

```bash
python -m analysis.fig03.gam_pygam.run_pygam --step maps --meanfd-max 0.5 --schaefer-atlas "resources/atlases/schaefer2018/Schaefer2018_400Parcels_7Networks_order.dlabel.nii" --output "outputs/fig03/gam_pygam_volume_432"
```

232：

```bash
python -m analysis.fig03.gam_pygam.run_pygam --meanfd-max 0.5 --tian-atlas "/public/home/dingrui/tools/masks_atlas/Tian_Subcortex_S2_3T_2009cAsym.nii.gz" --input "/public/home/dingrui/fmri_analysis/data/beh/participants_in_tfmri_demographics.csv" --subjects "/public/home/dingrui/fmri_analysis/data/beh/sub_list_valid_4_ER.txt" --excluded "/public/home/dingrui/fmri_analysis/data/fmri/subs_removed_from_2ndlevel.txt" --parent-edu "/public/home/dingrui/fmri_analysis/data/beh/env_family_student_cbcl.tsv" --meanfd "/public/home/dingrui/fmri_analysis/data/fmri/qc_res_subs_runs_fmri.csv" --fig3-volume-dir "/public/home/dingrui/fmri_analysis/data/fmri/res_1stlevel" --fig3-volume-atlas "resources/atlases/schaefer2018/Schaefer2018_200Parcels_7Networks_order_FSLMNI152_2mm.nii.gz" --schaefer-atlas "resources/atlases/schaefer2018/Schaefer2018_200Parcels_7Networks_order.dlabel.nii" --system-mask-dir "/public/home/dingrui/fmri_analysis/res/res_bayes" --output "outputs/fig03/gam_pygam_volume_232"
```

## 200 脑区，Yeo 7，不按头动幅度过滤

范围与下面的头动版本相同：Schaefer 200 × Yeo 7，7 个网络、4 个系统、6 个杏仁核。年龄项仍是 5 个 B 样条，平滑参数仍在 `1e-5` 到 `1e5` 的 12 个值上用 GCV 选择。四个协变量都保留，包括 `meanFD`。不写 `--meanfd-max`，不按 0.5 或 0.2 删人。

结果直接写在 `outputs/fig03/gam_pygam_volume_200`。不覆盖 `gam_pygam_volume`，也不写入 `gam_pygam_volume_200_FD05`、`gam_pygam_volume_200_FD02` 和 `gam_pygam_volume_200_shared_inputs`。

```bash
python -m analysis.fig03.gam_pygam.run_pygam --input "/public/home/dingrui/fmri_analysis/data/beh/participants_in_tfmri_demographics.csv" --subjects "/public/home/dingrui/fmri_analysis/data/beh/sub_list_valid_4_ER.txt" --excluded "/public/home/dingrui/fmri_analysis/data/fmri/subs_removed_from_2ndlevel.txt" --parent-edu "/public/home/dingrui/fmri_analysis/data/beh/env_family_student_cbcl.tsv" --meanfd "/public/home/dingrui/fmri_analysis/data/fmri/qc_res_subs_runs_fmri.csv" --fig3-volume-dir "/public/home/dingrui/fmri_analysis/data/fmri/res_1stlevel" --fig3-volume-atlas "resources/atlases/schaefer2018/Schaefer2018_200Parcels_7Networks_order_FSLMNI152_2mm.nii.gz" --schaefer-atlas "resources/atlases/schaefer2018/Schaefer2018_200Parcels_7Networks_order.dlabel.nii" --system-mask-dir "/public/home/dingrui/fmri_analysis/res/res_bayes" --output "outputs/fig03/gam_pygam_volume_200"
```

## 200 脑区，Yeo 7，头动 0.5 和 0.2

分区是 Schaefer 200 × Yeo 7。同时拟合 7 个网络、4 个系统和 6 个杏仁核。年龄项仍是 5 个 B 样条，平滑参数仍在 `1e-5` 到 `1e5` 的 12 个值上用 GCV 选择。`meanFD` 大于 0.5 和大于 0.2 的被试分别去掉。

体积图只读一次，放在 `outputs/fig03/gam_pygam_volume_200_shared_inputs`。结果在 `outputs/fig03/gam_pygam_volume_200_FD05` 和 `outputs/fig03/gam_pygam_volume_200_FD02`。不覆盖 `gam_pygam_volume` 和 `gam_pygam_volume_400_7net`。

`inputs/design.csv` 的列与 mgcv 的 200 脑区头动版本相同。对照时看同一头动目录里的 `q_age` 和 `q_m3_vs_m1`。mgcv 这一版是 `k = 3` 的三次回归样条加 `REML`，不是 5 个 B 样条。

```bash
python -m analysis.fig03.gam_pygam.run_pygam --meanfd-max 0.5 0.2 --input "/public/home/dingrui/fmri_analysis/data/beh/participants_in_tfmri_demographics.csv" --subjects "/public/home/dingrui/fmri_analysis/data/beh/sub_list_valid_4_ER.txt" --excluded "/public/home/dingrui/fmri_analysis/data/fmri/subs_removed_from_2ndlevel.txt" --parent-edu "/public/home/dingrui/fmri_analysis/data/beh/env_family_student_cbcl.tsv" --meanfd "/public/home/dingrui/fmri_analysis/data/fmri/qc_res_subs_runs_fmri.csv" --fig3-volume-dir "/public/home/dingrui/fmri_analysis/data/fmri/res_1stlevel" --fig3-volume-atlas "resources/atlases/schaefer2018/Schaefer2018_200Parcels_7Networks_order_FSLMNI152_2mm.nii.gz" --schaefer-atlas "resources/atlases/schaefer2018/Schaefer2018_200Parcels_7Networks_order.dlabel.nii" --system-mask-dir "/public/home/dingrui/fmri_analysis/res/res_bayes" --output "outputs/fig03/gam_pygam_volume_200"
```

## 400 脑区，Yeo 7，不按头动幅度过滤

范围与下面的头动版本相同：Schaefer 400 × Yeo 7，7 个网络、4 个系统、6 个杏仁核。年龄项仍是 5 个 B 样条，平滑参数仍在 `1e-5` 到 `1e5` 的 12 个值上用 GCV 选择。四个协变量都保留，包括 `meanFD`。不写 `--meanfd-max`，不按 0.5 或 0.2 删人。

结果直接写在 `outputs/fig03/gam_pygam_volume_400_7net`。不覆盖 `gam_pygam_volume`、`gam_pygam_volume_400`，也不写入 `gam_pygam_volume_400_7net_FD05`、`gam_pygam_volume_400_7net_FD02` 和 `gam_pygam_volume_400_7net_shared_inputs`。

```bash
python -m analysis.fig03.gam_pygam.run_pygam --input "/public/home/dingrui/fmri_analysis/data/beh/participants_in_tfmri_demographics.csv" --subjects "/public/home/dingrui/fmri_analysis/data/beh/sub_list_valid_4_ER.txt" --excluded "/public/home/dingrui/fmri_analysis/data/fmri/subs_removed_from_2ndlevel.txt" --parent-edu "/public/home/dingrui/fmri_analysis/data/beh/env_family_student_cbcl.tsv" --meanfd "/public/home/dingrui/fmri_analysis/data/fmri/qc_res_subs_runs_fmri.csv" --fig3-volume-dir "/public/home/dingrui/fmri_analysis/data/fmri/res_1stlevel" --fig3-volume-atlas "resources/atlases/schaefer2018/Schaefer2018_400Parcels_7Networks_order_FSLMNI152_2mm.nii.gz" --schaefer-atlas "resources/atlases/schaefer2018/Schaefer2018_400Parcels_7Networks_order.dlabel.nii" --system-mask-dir "/public/home/dingrui/fmri_analysis/res/res_bayes" --output "outputs/fig03/gam_pygam_volume_400_7net"
```

## 400 脑区，Yeo 7，头动 0.5 和 0.2

分区是 Schaefer 400 × Yeo 7。同时拟合 7 个网络、4 个系统和 6 个杏仁核。体积图只读一次，放在 `outputs/fig03/gam_pygam_volume_400_7net_shared_inputs`。`meanFD` 大于 0.5 和大于 0.2 的被试分别去掉，结果在 `outputs/fig03/gam_pygam_volume_400_7net_FD05` 和 `outputs/fig03/gam_pygam_volume_400_7net_FD02`。

这个目录不覆盖 `gam_pygam_volume`、`gam_pygam_volume_400` 和 `gam_nested_volume_400`。

图都是 svg，文字保持为文字。

| 文件 | 内容 |
|---|---|
| `figures/significant_parcel_trajectories.svg` | `q_age < 0.05` 的脑区。每个脑区一种颜色。零点是这条曲线自身的平均值 |
| `figures/yeo_age_curves.svg` | 7 个网络 |
| `figures/system_age_curves.svg` | 4 个系统 |
| `figures/amygdala_age_curves.svg` | 6 个杏仁核 ROI |
| `figures/age_parent_edu.svg` | design 表里年龄和 `parent_edu` 的散点 |
| `statistics/pygam_age_comparison.csv` | 全部单元的 `p_age` 和 `q_age` |
| `statistics/significant_parcels.csv` | `q_age < 0.05` 的脑区。没有显著脑区时表是空的，不画轨迹图 |
| `inputs/design.csv` | 全部被试和协变量，外加是否进入分析、排除原因。与 mgcv 嵌套分析是同一张名单 |
| `inputs/dropped_parent_edu.csv` | 因父母教育去掉的被试 |
| `inputs/dropped_no_volume.txt` | 缺两张效应图的被试 |
| `inputs/dropped_meanfd.csv` | 这一头动阈值额外去掉的被试。只出现在 `_FD05`、`_FD02` |

样本从有效名单去掉二阶剔除开始。每个人一行，写在 `inputs/design.csv`。父母教育、两张效应图和 `meanFD` 对得上就填进协变量，对不上就空着。`in_analysis` 为 1 表示进入分析，`excluded_by` 是排除原因。模型只使用进入分析的行。性别和站点不用于删行。

| 列 | 含义 |
|---|---|
| `sub_id`、`age`、`gender`、`site_id` | 与人口学表一致 |
| `parent_edu` | 父母教育。没有单一数值时为空 |
| `meanFD` | EMO 与 SOC 的 `mean_fd` 平均。没有数值时为空 |
| `in_analysis` | 是否进入分析。1 进入，0 不进入 |
| `excluded_by` | 排除原因。进入分析的人为空 |

`excluded_by` 只记第一条：`parent_edu_not_in_table`、`parent_edu_missing_value`、`parent_edu_conflicting_values`、`no_effect_volume`、`meanFD_missing`。表本身读不到时是 `parent_edu_file_missing`、`parent_edu_columns_missing`、`parent_edu_unreadable`、`meanFD_file_missing` 或 `meanFD_unreadable`。整批因此停止、而这个人没有单独原因时，写成 `sample_stop:` 加上停止原因。冲突、匹配率不足一半或分析样本里还有人缺头动时，脚本仍会停止，停止前 `design.csv` 已经写好。

`--meanfd-max 0.5 0.2` 不另写一张表。共享缓存 `gam_pygam_volume_400_7net_shared_inputs/inputs/design.csv` 还没按头动幅度排除。`_FD05` 和 `_FD02` 的 `design.csv` 仍列出所有人，只把已经进入分析且超过该阈值的人改为 `in_analysis` 为 0，原因是 `meanFD_gt_0.5` 或 `meanFD_gt_0.2`。等于阈值的人保留。更早被排除的人保持原来的原因。旧缓存若只有进入分析的人，下次沿用体积均值时会补成这张全名单，不重读 nii。

在工程根目录、能导入 pygam 的环境里执行：

```bash
python -m analysis.fig03.gam_pygam.run_pygam --meanfd-max 0.5 0.2 --input "/public/home/dingrui/fmri_analysis/data/beh/participants_in_tfmri_demographics.csv" --subjects "/public/home/dingrui/fmri_analysis/data/beh/sub_list_valid_4_ER.txt" --excluded "/public/home/dingrui/fmri_analysis/data/fmri/subs_removed_from_2ndlevel.txt" --parent-edu "/public/home/dingrui/fmri_analysis/data/beh/env_family_student_cbcl.tsv" --meanfd "/public/home/dingrui/fmri_analysis/data/fmri/qc_res_subs_runs_fmri.csv" --fig3-volume-dir "/public/home/dingrui/fmri_analysis/data/fmri/res_1stlevel" --fig3-volume-atlas "resources/atlases/schaefer2018/Schaefer2018_400Parcels_7Networks_order_FSLMNI152_2mm.nii.gz" --schaefer-atlas "resources/atlases/schaefer2018/Schaefer2018_400Parcels_7Networks_order.dlabel.nii" --system-mask-dir "/public/home/dingrui/fmri_analysis/res/res_bayes" --output "outputs/fig03/gam_pygam_volume_400_7net"
```

服务器上已有 FSL 时，6 个杏仁核读取 `$FSLDIR/data/atlases/Juelich/Juelich-maxprob-thr25-2mm.nii.gz`。图谱在别的路径时加上 `--juelich-atlas` 和 `--juelich-labels`。

和 mgcv 的 400 × Yeo 7 对照时，看同一头动目录里的 `q_age` 和 `q_m3_vs_m1`，以及两边的 `inputs/design.csv`。被试、图谱、四个协变量和杏仁核定义一致。mgcv 使用 `--method gcv --basis ps --basis-k 5` 时，年龄基最接近这里的 5 个 B 样条。平滑参数搜索仍不同：这里是 12 个取值，mgcv 是连续的 `GCV.Cp`。
