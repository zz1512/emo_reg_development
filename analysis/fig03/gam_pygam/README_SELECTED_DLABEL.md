# 与原生 Schaefer surface 一致的地图与六脑区 dlabel

本地已生成的主要文件：

`E:/emo_reg_development/outputs/fig03/selected_roi_dlabel_native/Schaefer400_7Networks_selected6_nativeSurface_fslr32k.dlabel.nii`

| 原图谱编号 | 脑区 | 标签颜色 |
|---|---|---|
| 182 | 7Networks_LH_Default_PFC_17 | 绿色 `#00FF00` |
| 185 | 7Networks_LH_Default_PFC_20 | 绿色 `#00FF00` |
| 152 | 7Networks_LH_Default_Temp_4 | 红色 `#FF0000` |
| 154 | 7Networks_LH_Default_Temp_6 | 红色 `#FF0000` |
| 156 | 7Networks_LH_Default_Temp_8 | 沿用原 border 亮青色 `#00E5FF` |
| 157 | 7Networks_LH_Default_Temp_9 | 沿用原 border 亮青色 `#00E5FF` |

六区保留独立标签编号。相邻脑区即使颜色相同，也没有合并成一个 ROI。其余位置为透明的背景标签 0，右半球全为背景。该文件是用户指定脑区的手动标注，不代表六区都达到统计显著。

## 边界来源

2026-11-08 按用户确认的目标改为**原生 Schaefer surface 标准边界**。顶点赋值直接来自：

`resources/atlases/schaefer2018/Schaefer2018_400Parcels_7Networks_order.dlabel.nii`

六区 dlabel 保留原始编号对应的全部原生顶点，其余顶点设为透明背景 0。`brain_maps.py` 的 partial R²、显著 mask、FDR border 和人工边缘 border 同样直接按这份表面图谱的编号赋值；因此各文件使用相同脑区边界。不做体积投影、插值、平滑、膨胀或腐蚀。无需 neuromaps 或注册融合数据下载。

指定的 `Schaefer2018_400Parcels_7Networks_order_FSLMNI152_2mm.nii.gz` 继续用于生成独立的体积参考掩膜和核对编号，不再决定 surface 标签。432 分析仍从原来的体积 ROI 提取数据；模型、partial R²、p/q 和显著结论均不变，只改变表面可视化几何。232 自动使用对应 200 区原生 surface；若 200/400 图谱与组合体积编号不匹配，会明确报错。

旧目录 `outputs/fig03/selected_roi_dlabel/` 中的 `fromMNI` 文件保持原样，属于旧投影版本。它与原生表面标签不相同；升级后请从 Workbench 移除旧 dlabel/border，并重新加载新文件及更新后的统计图，不要混用旧 scene 中缓存的图层。

## Workbench 使用

1. 通过 `File → Open File` 加载主 `.dlabel.nii`；配合 fsLR32k 左半球表面，例如工程中的 `resources/templates/fslr_32k/fsaverage.L.inflated.32k_fs_LR.surf.gii`。`Q1-Q6_R440.L.very_inflated.32k_fs_LR.surf.gii` 也是工程中可用的表面。
2. 在 Overlay Toolbox 的图层中启用该 dlabel，并将其放在 partial R² 统计图层上方。
3. 打开该图层的扳手设置（Overlay and Map Settings / Display and Map Settings），进入 `Labels`。
4. 将 `Drawing Type` 设置为 **`Outline Label Color`**，使用文件内预先写好的绿、红、亮青色。`Outline Color` 模式使用统一轮廓色，无法显示这里的三组标签颜色。
5. 保持叠加层不透明度为 1。标签背景已经透明，轮廓模式下脑区内部由下面的统计图层显示。
6. 需要保留叠加顺序和轮廓模式时，保存 Workbench scene。颜色保存在 dlabel 标签表中；Drawing Type 是查看器/scene 设置，单独打开 dlabel 时可能默认显示实心填色。[Workbench 绘图模式定义](https://github.com/Washington-University/workbench/blob/master/src/Files/LabelDrawingTypeEnum.cxx)、[查看器属性与 scene 保存](https://github.com/Washington-University/workbench/blob/master/src/Files/LabelDrawingProperties.cxx)

旧版本也可能在 `Features Toolbox → Labels → Attributes` 提供相应的 `Draw as` 选项。这里没有把 border 的 `Line Diameter` 写入 dlabel；需要可调粗细的折线时，仍应使用 `.border`。

## 同目录辅助文件

- `Schaefer400_7Networks_selected6_nativeSurface_fslr32k.L.label.gii` / `.R.label.gii`：相同顶点标签、相同颜色，显式包含左右脑结构信息。通常只加载主 dlabel 即可。
- `Schaefer400_7Networks_selected6_source_volume.nii.gz`：原体积空间中的六区掩膜，保留原编号、体素位置和 affine，供核对。
- `selected_rois.tsv`：六区名称、RGB、体素数及左右顶点数。
- `manifest.json`：原生 surface 图谱及参考体积的路径、SHA256、无投影的赋值方法。
- `workbench_label_table.txt`：由本地 Workbench 从实际生成的 dlabel 导出的标签表，用于核对颜色。
- `native_selected6_preview.svg`：真实左半球侧面/内侧六区范围和轮廓预览；颜色表示手动选择的 ROI 身份，不表示统计效应或显著性。这是图谱核对图，不是 Workbench scene。
- `native_vs_volume_alignment.svg`：旧投影填色与新原生填色的对照，两幅都叠加相同的原生轮廓；仅当旧 dlabel 存在时输出。
- `preview_validation.json`：新版六区标签与原生标签逐顶点完全一致的检查结果，以及新旧顶点数对照。

预览仅输出 SVG；其中密集三维网格嵌入为 300 dpi 渲染图以控制文件大小，标题、图例等文字保留可编辑 SVG。此设置不改变 dlabel、GIFTI 或边界几何。

## 重建命令

入口为 `analysis/fig03/gam_pygam/make_selected_dlabel.py`，不读取被试数据或重新拟合统计模型。复用现有 nibabel、NumPy、pandas 依赖；预览另需 matplotlib、nilearn。

在工程根目录执行：

```bash
python -m analysis.fig03.gam_pygam.make_selected_dlabel
python -m analysis.fig03.gam_pygam.preview_native_surface
```

本地 PowerShell 命令：

```powershell
Set-Location E:/emo_reg_development
& 'C:/Users/dell/anaconda3/envs/nilearn/python.exe' -m analysis.fig03.gam_pygam.make_selected_dlabel
& 'C:/Users/dell/anaconda3/envs/nilearn/python.exe' -m analysis.fig03.gam_pygam.preview_native_surface
```

主入口生成 dlabel、左右 GIFTI、体积参考掩膜、清单及 manifest；预览由第二个入口生成。Workbench 标签表可另行导出：

```bash
wb_command -cifti-label-export-table "outputs/fig03/selected_roi_dlabel_native/Schaefer400_7Networks_selected6_nativeSurface_fslr32k.dlabel.nii" 1 "outputs/fig03/selected_roi_dlabel_native/workbench_label_table.txt"
```

## 云端替换文件与执行命令

在工程中直接复制以下文件到云端相同相对路径，不需要打包：

- `analysis/fig03/gam_pygam/brain_maps.py`：所有表面地图和 border 改用原生标签；保存来源校验信息。
- `analysis/fig03/gam_pygam/run_pygam.py`：`--step borders` 与 maps 使用同一个指定的表面图谱。
- `analysis/fig03/gam_pygam/make_selected_dlabel.py`：原生六区 dlabel。
- `analysis/fig03/gam_pygam/README.md`、`README_SELECTED_DLABEL.md`：本次说明。
- 可选 `analysis/fig03/gam_pygam/preview_native_surface.py`：云端也生成预览时复制。
- 可选 `analysis/fig03/gam_pygam/test_borders.py`、`test_models.py`：运行验证时复制。

复用云端已有 Schaefer400 原生 dlabel 和左右 fsLR32k surface。Workbench 仍需在 PATH；缺失时使用 `--wb-command /实际路径/wb_command`。在云端工程根目录、原分析环境中执行：

```bash
python -m analysis.fig03.gam_pygam.run_pygam \
  --step maps \
  --schaefer-atlas "resources/atlases/schaefer2018/Schaefer2018_400Parcels_7Networks_order.dlabel.nii" \
  --output "outputs/fig03/gam_pygam_volume_432_models3_FD05"

python -m analysis.fig03.gam_pygam.make_selected_dlabel
```

第一条同时覆盖目标目录 `maps/` 下两种比较的 surface partial R²、mask 和 border，并按原方法重写体积文件；不修改 `statistics/` 或 `inputs/`，不拟合、不重新提取、不计算 FDR、不改随龄图。第二条生成独立目录中的新版六区 dlabel。**从旧版升级必须使用 maps，不能只重建 border 而继续使用旧的线性投影底图。**

本地没有云端的真实 432 统计表，所以此次本地交付的是实际图谱的新版六区文件与几何效果图；真实 partial R² 地图在云端运行第一条命令即可获得。

验证命令：`python -m unittest analysis.fig03.gam_pygam.test_borders analysis.fig03.gam_pygam.test_models -q`。验证覆盖逐顶点的地图/dlabel/掩膜赋值、原生标签与体积分辨率匹配、实际 Workbench border 转换与颜色、232/432 导出及无需重新拟合的重建路径。
