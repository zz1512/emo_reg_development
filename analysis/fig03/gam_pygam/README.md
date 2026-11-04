# pygam 年龄模型

四个协变量 `C` 为 `gender + site_id + parent_edu + meanFD`。当前拟合三个模型：

| 模型 | 公式与估计 |
|---|---|
| m1 | 截距 + C，OLS |
| m2 | 截距 + age + C，OLS |
| m3 | 截距 + s(age) + C，pyGAM；5 个 B 样条，只惩罚年龄项 |

m3 的平滑参数继续在 `1e-5` 到 `1e5` 的 12 个对数等距值上用 GCV 选择。m3 对 m1 检验整体年龄效应，保留历史的 RSS/EDF 近似 F 检验。

m3 对 m2 检验年龄函数是否偏离直线：在m3相同的固定B样条基函数空间内额外做一个无惩罚OLS参考拟合，与m2做嵌套F检验。其分子自由度为 `rank(B_m3) - rank(X_m2)`，分母自由度为 `n - rank(B_m3)`，F使用这个参考拟合的RSS。这样避免平滑曲线接近直线时，直接用极小EDF差导致虚假的显著结果。参考拟合仅用于非线性检验；绘图、三个模型的RSS/EDF及两种partial R²仍描述原来的平滑m3。该比较检验的是给定样条空间内的非线性，结果不等同于mgcv的 `anova.gam`。

两种检验沿用独立、同方差高斯误差假设。mgcv也提醒多模型EDF检验对可被完全惩罚掉的部分可能严重高估显著性，见[mgcv说明](https://stat.ethz.ch/R-manual/R-devel/library/mgcv/html/anova.gam.html)；普通嵌套F的假设见[statsmodels说明](https://www.statsmodels.org/stable/generated/statsmodels.regression.linear_model.RegressionResults.compare_f_test.html)。

BH 对每种比较、每个分析家族分别进行。当前 Tian432 有 parcel432、Yeo7、system4、Juelich6 四个家族；Tian232 仅有 parcel232。旧200/400皮层与Kong方案继续保留，命令和适用范围见后文。

## 历史运行清单：200/400皮层方案

本节保留已有200/400皮层方案的执行方式。当前Tian232/432方案使用下一节的独立 `models3` 输出目录；重复执行相同输出名会更新该目录自己的结果。

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

## 当前方案：Tian432/232，三个模型，只过滤头动0.5

两版均在Schaefer皮层后接Tian Subcortex S2的32个皮下区。432的皮下编号为401–432；232的皮下编号为201–232。Tian名称读取沿用S2 3T正式顺序：1–16右侧、17–32左侧，图谱旁若有完整名称表则采用该表。Tian四个杏仁核是3 `lAMY-rh`、4 `mAMY-rh`、19 `lAMY-lh`、20 `mAMY-lh`，继续属于parcel家族。

| 版本 | 拟合及曲线范围 | 每种比较的BH家族 |
|---|---|---|
| 432 | 432个parcel、7个Yeo网络、4个系统、6个Juelich杏仁核 | 432、7、4、6分别校正 |
| 232 | 232个parcel | 232一起校正 |

432的Yeo7只使用皮层parcel的体素加权均值，四系统仍用贝叶斯mask均值。独立的Juelich6是LB、CM、SF × 左右，与Tian四杏仁核分区分别记录；Juelich6单独校正，Tian四区参加432家族校正。232仅拟合parcel，忽略旧缓存中可能保留的组表及独立杏仁核表。

两个比较 `m3_vs_m1`、`m3_vs_m2` 各自按上述家族执行BH，阈值 `q < 0.05`。两种partial R²分别是 `max(0, 1 - RSS_m3 / RSS_m1)` 和 `max(0, 1 - RSS_m3 / RSS_m2)`。统计主表仍为 `statistics/pygam_age_comparison.csv`，为每个比较写出 `F_<comparison>`、`df1_<comparison>`、`df2_<comparison>`、`p_<comparison>`、`q_<comparison>`、`reject_<comparison>_fdr`、`partial_r2_<comparison>`，并保存三个模型的RSS、EDF及m2年龄系数/标准误。另存 `rss_m3_basis`、`rank_m3_basis`，用于复核m3对m2检验的OLS参考拟合；它们与平滑模型的 `rss_m3`、`edf_m3` 分开。旧列 `p_age`、`q_age`、`partial_r2` 保留为m3对m1的兼容别名。

只去掉 `meanFD > 0.5` 的被试，等于0.5保留，meanFD继续作为协变量。新命令使用独立目录：

| 版本 | 输出base | 拟合结果 | 全量提取时的共享均值缓存 |
|---|---|---|---|
| 432 | `outputs/fig03/gam_pygam_volume_432_models3` | `outputs/fig03/gam_pygam_volume_432_models3_FD05` | `outputs/fig03/gam_pygam_volume_432_models3_shared_inputs/inputs` |
| 232 | `outputs/fig03/gam_pygam_volume_232_models3` | `outputs/fig03/gam_pygam_volume_232_models3_FD05` | `outputs/fig03/gam_pygam_volume_232_models3_shared_inputs/inputs` |

### 曲线和脑图

所有纳入绘图的单元各输出两张可编辑文字的SVG：固定纵轴的散点总览图，以及仅含同一曲线和95%置信区间的局部放大图。parcel取任一比较显著的并集，再加入Tian四杏仁核及432版的人工边缘156/157，同一parcel只入选一次。432另为全部7个网络、4个系统、6个Juelich杏仁核各输出这两张图；232仅画parcel。

散点、曲线及置信区间统一采用m3的中心化年龄效应，不再混用m1残差与m3曲线。每个单元以有效被试的实际年龄分布为中心化参考，协变量固定在这些被试的均值处。设 `B(a, C均值)` 为m3设计矩阵的一行，`B均值` 为该单元有效被试年龄对应设计行的平均值，则：

- 年龄曲线：`L(a) = B(a, C均值) - B均值`，`effect(a) = L(a) × beta_m3`。
- 部分残差散点：`y_i - m3(age_i, C_i) + effect(age_i)`。
- 95%点态置信区间：`effect(a) ± t(0.975, n-edf_m3) × sqrt(L(a) × Cov(beta_m3) × L(a)转置)`；已知误差尺度时使用正态分位数。

这里包含了中心化参考与曲线之间的协方差，去除了共同截距和协变量基线的不确定性。参考来自实际被试年龄，不随绘图网格密度改变。该区间描述中心化平均年龄效应，仍以GCV选定的平滑参数为条件，未包含平滑参数选择的不确定性，也不是个体预测区间。即使散点很分散，区间仍可能很窄；接近线性的中心化曲线在平均年龄附近也可能明显收窄。不会为了图形效果放大区间数值。长表 `statistics/covariate_residuals.csv` 的列仍为 `sub_id, age, dimension, name, residual`，其中 `residual` 现在保存上述m3部分残差；Tian432/232的 `pygam_curves.csv` 中 `effect` 与 `centered_effect` 均保存同一中心化曲线。三个模型、比较检验、BH校正和partial R²的计算保持原有定义。

总览图纵轴按以下规则固定，绘图区宽高比仍为3:4，曲线与置信带使用亮金色：

| 绘图单元 | 纵轴范围 |
|---|---|
| 432版任一比较显著的parcel，以及156、157 | `[-1.5, 1.5]` |
| 其余单元：432版仅因Tian杏仁核标注而绘制的非显著parcel、Yeo7、四系统、Juelich6，以及232版全部parcel | `[-1, 1]` |

显著指 `q_m3_vs_m1 < 0.05` 或 `q_m3_vs_m2 < 0.05`。散点较淡，置信带位于散点之上；去掉0参考线、上边框及右边框。范围外的残差和曲线/置信区间保留在数据表中，清单的 `y_min`、`y_max` 记录总览图实际纵轴，`n_outside_y_limits` 按各图自身范围统计画面外的部分残差散点数。`statistics/trajectory_units.csv` 记录各类别绘图单元；`trajectory_parcels.csv` 仅记录parcel。清单保留两个比较的q值。

总览图继续写入 `figures/age_curves/`，文件名保持不变，例如 `parcel_0154.svg`、`yeo_Vis.svg`、`system_reappraisal_only.svg`、`amygdala_Juelich_LB_L.svg`。放大图写入独立的 `figures/age_curves_zoom/`，对应文件名增加 `_zoom`，例如 `parcel_0154_zoom.svg`。放大图不叠加散点，使用与总览图完全相同的年龄曲线、CI和横轴范围；纵轴包含完整曲线和CI包络，并在两端各留出包络跨度的15%。极小或退化包络采用最小留白，避免相同上下限。放大图同样使用3:4绘图区；各单元纵轴会不同，跨单元比较效应幅度时应看总览图及数值表。清单增加 `zoom_figure_file`、`zoom_y_min`、`zoom_y_max` 以记录放大图文件及纵轴。

统一m3绘图口径并增加放大图时需要本README、`run_pygam.py`、`plot_region_age.py`；当前版本还应同步 `brain_maps.py`，以支持下文新增的边界导出。**首次升级绘图口径必须从现有均值重新拟合一次**，旧CSV未保存重算部分残差和中心化CI所需的模型信息，不能仅重画。在工程根目录、`er-manuscript` 环境执行下面两条命令，不重读原始被试effectSize；它们会更新现有结果目录内的统计表、两套曲线图及后续脑图和边界：

```bash
python -m analysis.fig03.gam_pygam.run_pygam --step fit --schaefer-atlas "resources/atlases/schaefer2018/Schaefer2018_400Parcels_7Networks_order.dlabel.nii" --output "outputs/fig03/gam_pygam_volume_432_models3_FD05"
python -m analysis.fig03.gam_pygam.run_pygam --step fit --schaefer-atlas "resources/atlases/schaefer2018/Schaefer2018_200Parcels_7Networks_order.dlabel.nii" --output "outputs/fig03/gam_pygam_volume_232_models3_FD05"
```

这两条命令沿用目标目录已经完成的FD0.5样本筛选。升级后的 `run_metadata.json` 含 `age_effect_definition: m3_empirical_centered_age_v1`。之后只调整样式可运行 `--step curves --output "..._FD05"`，同时重画两套图并刷新metadata中的绘图样式，保留原统计信息。未带此版本标记的旧结果会明确提示先执行fit，防止把旧残差和旧区间当成新版结果。

脑图按两种比较分别输出，在 `maps/` 下保存全脑partial R²和显著mask：

| 文件模式 | 内容 |
|---|---|
| `partial_r2_<comparison>_volume_MNI152.nii.gz` | 对应比较的全脑parcel partial R² |
| `significance_mask_<comparison>_volume_MNI152.nii.gz` | 0背景、1对应比较q<0.05、2人工指定156/157 |
| `partial_r2_<comparison>_fslr32k_lh.func.gii`、`_rh.func.gii` | fsLR32k皮层partial R² |
| `significance_mask_<comparison>_fslr32k_lh.func.gii`、`_rh.func.gii` | fsLR32k皮层显著mask |
| `partial_r2_<comparison>_fslr32k.dscalar.nii` | 皮层partial R² CIFTI |
| `significance_mask_<comparison>_fslr32k.dscalar.nii` | 皮层显著mask CIFTI |

其中 `<comparison>` 为 `m3_vs_m1` 或 `m3_vs_m2`。156、157仅适用于432版，两个mask中均人工标成2，并优先于显著标签1；232没有这个人工边缘集合。volume保留皮下区，fsLR文件存皮层。输出文件供手动导入Workbench出图。

左右GIFTI在文件级metadata中明确写入 `AnatomicalStructurePrimary`：`lh`为 `CortexLeft`，`rh`为 `CortexRight`，覆盖两种比较的partial R²和显著mask。Workbench据此识别半球；CIFTI本身已有左右皮层的BrainModelAxis。此处仅补充结构标记，保留原有投影metadata和数值。

2026-11-04修复Workbench半球识别时，云端需替换 `analysis/fig03/gam_pygam/brain_maps.py` 和本README。以后全量运行及 `--step maps` 导出均自动写入标记。

已经生成的文件可以用[Workbench官方的 `-set-structure` 命令](https://www.humanconnectome.org/software/workbench-command/-set-structure)直接修复。在工程根目录的Bash中执行下面的循环，覆盖432/232的两种比较及两类GIFTI；缺失的目录或文件自动跳过。只就地补充文件内部半球信息，无需重新拟合、投影或重画曲线；修复后在Workbench中重新加载GIFTI：

```bash
for maps_dir in outputs/fig03/gam_pygam_volume_432_models3_FD05/maps outputs/fig03/gam_pygam_volume_232_models3_FD05/maps; do
  for file in "$maps_dir"/*_fslr32k_lh.func.gii; do
    [ -f "$file" ] || continue
    wb_command -set-structure "$file" CORTEX_LEFT
  done
  for file in "$maps_dir"/*_fslr32k_rh.func.gii; do
    [ -f "$file" ] || continue
    wb_command -set-structure "$file" CORTEX_RIGHT
  done
done
```

### 独立FDR边界与156/157红框

当前Tian方案的全量运行、fit和maps步骤会自动调用Workbench的 [`-metric-rois-to-border`](https://www.humanconnectome.org/software/workbench-command/-metric-rois-to-border)，在 `maps/borders/` 输出独立ROI及 `.border`。已有统计结果只补边界时，可在工程根目录直接执行，无需重新拟合、重提取被试均值或重写partial R²图：

```bash
conda activate er-manuscript
python -m analysis.fig03.gam_pygam.run_pygam --step borders --output "outputs/fig03/gam_pygam_volume_432_models3_FD05"
python -m analysis.fig03.gam_pygam.run_pygam --step borders --output "outputs/fig03/gam_pygam_volume_232_models3_FD05"
```

该步骤读取目标目录的 `statistics/pygam_age_comparison.csv`、`inputs/tian_combined.json` 和其中记录的 `combined_atlas.nii.gz`，将组合图谱编号以nearest方法投到fsLR32k一次，再根据parcel编号生成0/1 ROI。所用MNI到fsLR映射与现有mask一致；不会对二值ROI进行平滑或线性插值。FDR边界直接取相应比较的有限 `q < 0.05`，不是partial R²阈值，也不把旧0/1/2混合mask中的所有正值视为显著。156、157是独立人工标注；若其中某区本身通过FDR，它也会保留在对应FDR边界中。

| `maps/borders/` 中的文件 | 内容 |
|---|---|
| `fdr_significant_m3_vs_m1_fslr32k_lh.border`、`_rh.border` | m3对m1的FDR边界 |
| `fdr_significant_m3_vs_m2_fslr32k_lh.border`、`_rh.border` | m3对m2的FDR边界 |
| `marginal_156_157_fslr32k_lh.border`、`_rh.border` | 仅432版，156与157各有独立名称及轮廓，存于一份文件 |
| 同名 `.func.gii` | 边界的二值ROI来源；FDR文件1列，156/157文件2列，每列独立0/1 |
| `border_manifest.json` | 文件、半球、使用的surface、每个ROI顶点数及边界段数 |

432版共6份 `.border`，232版只有两个比较×左右半球的4份，不添加156/157。某半球没有相应ROI时仍输出带半球信息的有效空边界，清单中的顶点数及边界段数为0；空文件不代表转换失败。`.border`用于皮层表面显示，Tian皮下区结果仍在volume中查看。

默认使用工程已有的 `resources/templates/fslr_32k/fsaverage.L.inflated.32k_fs_LR.surf.gii` 及对应R文件，两者实际采用fsLR32k编号、每侧32492顶点。此命令通过surface的网格邻接关系追踪边界，使用同网格的inflated表面即可，不必另行下载midthickness。导出前检查左右半球、顶点数和带编号的三角网格；自定义surface须与这套fsLR32k网格一致，不能只凭相同顶点数判断。

若云端资源目录没有这两份surface，将工程中的对应文件复制到相同位置。Workbench需在PATH中；否则在命令后追加 `--wb-command "/实际路径/wb_command"`。已有匹配的midthickness时可追加 `--border-surfaces "/实际路径/L.midthickness.surf.gii" "/实际路径/R.midthickness.surf.gii"`，顺序固定为左、右。缺少Workbench或网格不匹配会明确报错。

输出文件的Name和Class颜色均设为红色。加载partial R²图作为底图，再加载所需FDR及156/157 `.border`；在Workbench的 **Features Toolbox → Borders → Attributes** 设置 **Display Borders**、**Draw As = Polylines**、**Coloring = Standard Color**、**Standard Color = Red**，**Line Diameter先设2–3 px**。显示开关、绘制方式及线宽由Workbench界面/scene保存，`.border`本身不固定这些设置；Name/Class着色时文件内的颜色也为红色。相关控件见[Workbench官方界面实现](https://github.com/Washington-University/workbench/blob/master/src/GuiQt/BorderSelectionViewController.cxx)。FDR和人工边界可以独立加载及显隐；partial R²底图保持连续颜色，不用显著性ROI替换。

### 云端需替换的文件

每次修改分析脚本时同步本README、可运行命令及替换清单。本次独立红框/border更新，只需将以下三个文件同步到服务器工程相同相对位置：

- `analysis/fig03/gam_pygam/README.md`
- `analysis/fig03/gam_pygam/run_pygam.py`
- `analysis/fig03/gam_pygam/brain_maps.py`

若服务器尚未同步此前的m3部分残差/CI和两套随龄图更新，还需同步 `analysis/fig03/gam_pygam/plot_region_age.py`。此前三模型、四家族的提取修复另涉及 `analysis/fig03/gam_supplement/prepare_inputs.py`；已同步则无需再次替换。

默认使用前述两份工程已有的fsLR32k surface，不额外生成或打包模板。

继续依赖工程已有的 `tian_combined.py`、`gam_nested/amygdala.py` 和 `src/prj_analysis/src/mylib/stat/stats.py`。在含 `analysis/` 的工程根目录、`er-manuscript` 环境执行下列命令。

本地回归测试另见 `analysis/fig03/gam_pygam/test_models.py`、`test_borders.py`，不参与云端分析；如需运行，使用含pyGAM的环境执行 `python -m unittest analysis.fig03.gam_pygam.test_models analysis.fig03.gam_pygam.test_borders -v`。边界集成测试实际调用Workbench生成并回读红色边界；未安装Workbench时跳过这些集成项。

下文直接给出432/232的完整全量命令；需要指定独立Juelich图谱时，可在432命令后追加 `--juelich-atlas ... --juelich-labels ...`。已有完整均值时，本次采用前述fit升级命令；完成升级后才可只重画。半球标记可用前述补标记命令单独修复。

### 432全量命令

从原始effectSize提取432、Yeo7、四系统和Juelich6，再按FD0.5筛选及拟合。Juelich默认从 `$FSLDIR` 查找，再尝试nilearn下载；另存图谱时传 `--juelich-atlas` 及 `--juelich-labels`。

```bash
python -m analysis.fig03.gam_pygam.run_pygam --meanfd-max 0.5 --tian-atlas "/public/home/dingrui/tools/masks_atlas/Tian_Subcortex_S2_3T_2009cAsym.nii.gz" --input "/public/home/dingrui/fmri_analysis/data/beh/participants_in_tfmri_demographics.csv" --subjects "/public/home/dingrui/fmri_analysis/data/beh/sub_list_valid_4_ER.txt" --excluded "/public/home/dingrui/fmri_analysis/data/fmri/subs_removed_from_2ndlevel.txt" --parent-edu "/public/home/dingrui/fmri_analysis/data/beh/env_family_student_cbcl.tsv" --meanfd "/public/home/dingrui/fmri_analysis/data/fmri/qc_res_subs_runs_fmri.csv" --fig3-volume-dir "/public/home/dingrui/fmri_analysis/data/fmri/res_1stlevel" --fig3-volume-atlas "resources/atlases/schaefer2018/Schaefer2018_400Parcels_7Networks_order_FSLMNI152_2mm.nii.gz" --schaefer-atlas "resources/atlases/schaefer2018/Schaefer2018_400Parcels_7Networks_order.dlabel.nii" --system-mask-dir "/public/home/dingrui/fmri_analysis/res/res_bayes" --output "outputs/fig03/gam_pygam_volume_432_models3"
```

### 232全量命令

仅提取232个parcel并拟合。此命令不需要系统mask及Juelich图谱。

```bash
python -m analysis.fig03.gam_pygam.run_pygam --meanfd-max 0.5 --tian-atlas "/public/home/dingrui/tools/masks_atlas/Tian_Subcortex_S2_3T_2009cAsym.nii.gz" --input "/public/home/dingrui/fmri_analysis/data/beh/participants_in_tfmri_demographics.csv" --subjects "/public/home/dingrui/fmri_analysis/data/beh/sub_list_valid_4_ER.txt" --excluded "/public/home/dingrui/fmri_analysis/data/fmri/subs_removed_from_2ndlevel.txt" --parent-edu "/public/home/dingrui/fmri_analysis/data/beh/env_family_student_cbcl.tsv" --meanfd "/public/home/dingrui/fmri_analysis/data/fmri/qc_res_subs_runs_fmri.csv" --fig3-volume-dir "/public/home/dingrui/fmri_analysis/data/fmri/res_1stlevel" --fig3-volume-atlas "resources/atlases/schaefer2018/Schaefer2018_200Parcels_7Networks_order_FSLMNI152_2mm.nii.gz" --schaefer-atlas "resources/atlases/schaefer2018/Schaefer2018_200Parcels_7Networks_order.dlabel.nii" --output "outputs/fig03/gam_pygam_volume_232_models3"
```

### 复用已保存的均值

`--reuse-inputs PATH` 仅在默认 `--step all` 时使用，PATH指向已保存的 `inputs/` 目录。源目录须有完整的 `design.csv`、`parcel_info.csv`、`parcel_response.csv`、`group_response.csv`、`tian_combined.json`、`combined_atlas.nii.gz`；232的group表可以为空。源均值复制到新输出，再根据 `--meanfd-max` 筛选，源目录保持不变。

已有旧232共享缓存时，可以直接运行：

```bash
python -m analysis.fig03.gam_pygam.run_pygam --reuse-inputs "outputs/fig03/gam_pygam_volume_232_shared_inputs/inputs" --meanfd-max 0.5 --schaefer-atlas "resources/atlases/schaefer2018/Schaefer2018_200Parcels_7Networks_order.dlabel.nii" --output "outputs/fig03/gam_pygam_volume_232_models3"
```

432复用还要求源目录包含完整的Yeo7、四系统响应及Juelich6响应，六核按固定LB_L、LB_R、CM_L、CM_R、SF_L、SF_R顺序拟合。旧432缓存只有Tian四杏仁核时，脚本明确拒绝复用；使用上面的432全量命令补齐提取。对已经做过FD筛选的inputs再筛选时，先前排除者继续保留原排除状态；要从完整候选样本重筛，应选择 `_shared_inputs/inputs`。复用表示沿用该缓存的样本和协变量；数据源有变化时，换一个新的 `--output` base执行全量命令，确保重新提取。相同base的全量命令仍可能命中原共享均值缓存。

### 拟合及绘图恢复命令

`--step fit` 从目标输出中的完整均值重新拟合三个模型并更新后续图，不重读被试effectSize。输入仅有旧两模型统计，或绘图仍使用旧的m1残差/平移CI定义时，应先运行fit生成新比较和曲线表。432需已有完整的Yeo7、四系统及Juelich6输入，否则使用全量命令。

```bash
python -m analysis.fig03.gam_pygam.run_pygam --step fit --meanfd-max 0.5 --schaefer-atlas "resources/atlases/schaefer2018/Schaefer2018_400Parcels_7Networks_order.dlabel.nii" --output "outputs/fig03/gam_pygam_volume_432_models3"
python -m analysis.fig03.gam_pygam.run_pygam --step fit --meanfd-max 0.5 --schaefer-atlas "resources/atlases/schaefer2018/Schaefer2018_200Parcels_7Networks_order.dlabel.nii" --output "outputs/fig03/gam_pygam_volume_232_models3"
```

已有本次fit生成的 `pygam_age_comparison.csv`、`pygam_curves.csv`、长格式 `covariate_residuals.csv`，且 `run_metadata.json` 中的 `age_effect_definition` 为 `m3_empirical_centered_age_v1` 时，`--step curves` 只重画总览图及放大图：

```bash
python -m analysis.fig03.gam_pygam.run_pygam --step curves --meanfd-max 0.5 --output "outputs/fig03/gam_pygam_volume_432_models3"
python -m analysis.fig03.gam_pygam.run_pygam --step curves --meanfd-max 0.5 --output "outputs/fig03/gam_pygam_volume_232_models3"
```

已有新版统计表和组合图谱时，`--step maps` 重出脑图并生成独立边界，不重新拟合：

```bash
python -m analysis.fig03.gam_pygam.run_pygam --step maps --meanfd-max 0.5 --schaefer-atlas "resources/atlases/schaefer2018/Schaefer2018_400Parcels_7Networks_order.dlabel.nii" --output "outputs/fig03/gam_pygam_volume_432_models3"
python -m analysis.fig03.gam_pygam.run_pygam --step maps --meanfd-max 0.5 --schaefer-atlas "resources/atlases/schaefer2018/Schaefer2018_200Parcels_7Networks_order.dlabel.nii" --output "outputs/fig03/gam_pygam_volume_232_models3"
```

fit、curves、maps、borders均从目标inputs中的Tian标记识别模式，无需再传 `--tian-atlas`。fit/maps须指定与结果匹配的400或200 dlabel；curves和borders不需要dlabel参数。只补边界优先用前述 `--step borders` 命令。

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

## 历史图谱能力：Kong 2022的400 × 17网络标签

Kong图谱沿用 `--parcels-only` 模式：只拟合400个parcel，不合成17个网络，不提取系统或独立杏仁核。使用配套 `order_info.txt` 提供名称。此模式与Tian232/432分别运行，结果写入独立的旧17net目录；重新执行会更新该目录自己的结果。

```bash
python -m analysis.fig03.gam_pygam.run_pygam --parcels-only --meanfd-max 0.5 0.2 --input "/public/home/dingrui/fmri_analysis/data/beh/participants_in_tfmri_demographics.csv" --subjects "/public/home/dingrui/fmri_analysis/data/beh/sub_list_valid_4_ER.txt" --excluded "/public/home/dingrui/fmri_analysis/data/fmri/subs_removed_from_2ndlevel.txt" --parent-edu "/public/home/dingrui/fmri_analysis/data/beh/env_family_student_cbcl.tsv" --meanfd "/public/home/dingrui/fmri_analysis/data/fmri/qc_res_subs_runs_fmri.csv" --fig3-volume-dir "/public/home/dingrui/fmri_analysis/data/fmri/res_1stlevel" --fig3-volume-atlas "resources/atlases/schaefer2018/Schaefer2018_400Parcels_Kong2022_17Networks_order_FSLMNI152_2mm.nii.gz" --schaefer-atlas "resources/atlases/schaefer2018/Schaefer2018_400Parcels_Kong2022_17Networks_order.dlabel.nii" --parcel-info "resources/atlases/schaefer2018/Schaefer2018_400Parcels_Kong2022_17Networks_order_info.txt" --output "outputs/fig03/gam_pygam_volume_400_17net"
```
