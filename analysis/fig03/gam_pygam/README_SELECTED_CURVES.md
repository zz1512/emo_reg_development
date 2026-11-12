# 四脑区随龄图与 mock 样式预览

本次只重画指定 Schaefer400/Tian432 中的 152、154、182、185：

| 编号 | 脑区 | 散点/CI 基色 | 拟合线 |
|---|---|---|---|
| 152、154 | LH Default Temp 4、6 | 红色 `#FF0000` | 深红 `#CC2020` |
| 182、185 | LH Default PFC 17、20 | 绿色 `#00FF00` | 深绿 `#008A20` |

散点基色与本次六区 dlabel 的红/绿标签一致。拟合线用较深的同色系，以便从散点中辨认；95% 均值置信带用对应红/绿色，透明度 0.18，不加边框，不改变区间数值。散点面积从原图 8 pt² 增至 22 pt²，透明度 0.65，增加深灰色 0.4 pt 描边，对齐参考图的圆点风格。

字体统一声明 Arial，进一步增大为：坐标轴名称 24 pt、刻度数字 18 pt、脑区标题 16 pt、左上角注释 14 pt。横轴为 `Age (years)`，明确显示 6、8、10、12、14、16、18；纵轴保留 `Partial residual (m3)`，说明散点是模型部分残差。绘图区保持宽:高 = 3:4，去掉上、右边框。轴名称和刻度均为黑色，轴线为灰色。

纵轴按脑区固定：152、154 为 `[-1.5, 1.5]`，182、185 为 `[-1, 1]`，不再根据数据自动扩展或四图共用同一范围。两组纵轴刻度间隔均为 0.5。范围外的散点仅在画面中裁切，不从数据中删除；`trajectory_parcels_selected.csv` 记录实际纵轴、总散点数及范围外散点数，mock 的对应数量写入 `mock_metadata.json`。曲线与 CI 的数值保持原样。

正式图及 mock 图均只输出 SVG。重画时会移除该绘图函数先前生成的同名 PNG；不清理其他用户文件。

## 已生成的本地 mock 示例

- `E:/emo_reg_development/outputs/fig03/age_curves_mock/mock_comparison.svg`：152 红色与 182 绿色并排矢量预览。
- `mock_parcel_0152.svg`、`mock_parcel_0182.svg`：独立的可编辑矢量图，不生成 PNG。
- `MOCK_pygam_curves.csv`、`MOCK_partial_residuals.csv`：合成数据拟合得到的曲线、CI 与部分残差。
- `mock_metadata.json`：随机种子 20261108、每区 991 个合成观测、实际纵轴和样式参数。

mock 数据不来自真实被试，正负斜率仅用于评判显示风格，不代表脑区的真实发育方向。示例用真实 pyGAM 拟合流程计算中心化曲线和纳入均值不确定性的 95% 点态 CI，调用与真实数据相同的绘图函数。mock 图左上角明确标注 `MOCK DATA`；真实图改为原始全家族校正的两个 q 值，不计算或仿造参考图片中的 Pearson r。

在可导入 pyGAM 的环境、工程根目录执行：

```bash
python -m analysis.fig03.gam_pygam.mock_age_curves
```

本次 Windows 本地实际使用：

```powershell
Set-Location E:/emo_reg_development
$env:PYTHONPATH = 'C:/Users/dell/AppData/Local/Temp/fig03-age-effect-3kjh9arl/packages;E:/emo_reg_development'
& 'C:/Users/dell/anaconda3/envs/nilearn/python.exe' -m analysis.fig03.gam_pygam.mock_age_curves
```

临时目录提供独立安装的 pyGAM；已有 pyGAM 的环境无须此临时路径。需要生成全部四区的 mock 时增加 `--parcels 152 154 182 185`。

## 云端真实数据：只重画四区

本次样式更新需要替换绘图脚本并同步说明；此前已同步过的 `run_pygam.py` 本次没有改动：

1. `analysis/fig03/gam_pygam/plot_region_age.py`：固定各区纵轴、更大字号、只输出 SVG。
2. `analysis/fig03/gam_pygam/README.md` 和本说明 `README_SELECTED_CURVES.md`。

如果还未同步上一版四区绘图入口，也需替换 `analysis/fig03/gam_pygam/run_pygam.py`，它提供 `--curve-parcels` 参数。

`mock_age_curves.py` 仅用于样式预览，云端真实重画不需要；`test_models.py` 为验证文件，也不需要上传。

在云端工程根目录、原来可运行分析的 Python 环境执行：

```bash
python -m analysis.fig03.gam_pygam.run_pygam \
  --step curves \
  --curve-parcels 152 154 182 185 \
  --output "outputs/fig03/gam_pygam_volume_432_models3_FD05"
```

只读取目标目录已经保存的统计、曲线和部分残差，不读取原始影像、不重新拟合、不重做 FDR。需要已有当前版本的中心化年龄效应及纳入均值不确定性的 CI 标记；如果缓存尚未升级，程序仍按既有规则提示先拟合，不会偷偷转换旧区间。

四区按指定清单绘制，不以 q 值决定是否绘制。真实图左上角显示原表 `q_m3_vs_m1` 和 `q_m3_vs_m2`，不能因手动绘图就称它们均显著。选择只作用于绘图，不缩小原统计检验的校正家族。

新文件写入：

```text
<output>/figures/age_curves_selected/parcel_0152.svg
<output>/figures/age_curves_selected/parcel_0154.svg
<output>/figures/age_curves_selected/parcel_0182.svg
<output>/figures/age_curves_selected/parcel_0185.svg
<output>/statistics/trajectory_parcels_selected.csv
```

本模式只生成这四区的带散点总览图，不生成网络、系统、杏仁核图或 zoom 图。原 `figures/age_curves/`、`age_curves_zoom/`、原轨迹清单及模型统计不被覆盖；`run_metadata.json` 单独增加 `selected_parcel_curve_style`，保留原绘图样式记录。不加 `--curve-parcels` 时，原来的全部选区/放大图绘制流程保持原样。

Arial 在本地已确认可用。SVG 保留可编辑文字和 Arial 字体声明；查看/编辑时安装 Arial 可使字形与本地一致。缺少 Arial 的环境会退回可用的 Liberation Sans / DejaVu Sans。

## 验证

本次重画检查使用包含四区及范围外散点的合成缓存，验证：非显著指定区也能输出、四区分别使用指定固定纵轴、范围外散点计数正确、只输出 SVG 并移除旧同名 PNG、原曲线/CI/统计文件不变、原图不被覆盖、重画过程不调用拟合。

```bash
python -m unittest analysis.fig03.gam_pygam.test_models
```
