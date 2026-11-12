# PFC 范围的原文核查与本项目建议

核查日期：2026-11-11。目的：为“PFC vs PFC+非PFC”的跨站点预测固定独立ROI定义。
以下区分作者实际使用的掩膜与本项目建议；没有依据当前样本年龄/行为显著性选区，也没有生成最终PFC名单。

## 1. 用户提供的 Nature 论文并未提取全部解剖 PFC

Sicorello et al. (2026), *Emotion regulation capacity lacks a shared neural signature*.
[论文](https://www.nature.com/articles/s41586-026-11203-1)，DOI：10.1038/s41586-026-11203-1。
用户提供的 `s41586-026-11203-1.pdf`：PDF第10页，Methods → Statistical analyses；
PDF第24页，Reporting summary → Anatomical location(s)。已核对第10页原文及页面。

- 理论表述讨论前额叶调节环路，但实际先验分析使用 Morawetz et al. (2020) 的两个情绪调节网络。
- 网络1、2主要包含前额叶、顶叶及颞叶；它们属于**分布式功能网络**，不能直接当作解剖PFC掩膜。
- 作者先分别取两个网络平均值，再在合并网络内做体素分析，最后扩展到全脑。
- 杏仁核ROI另行定义：SPM Anatomy Toolbox v2.2b解剖杏仁核与功能性情绪生成网络的交集。
- 正文Fig.3部分还使用Yeo网络汇总；那也不是全部解剖PFC的定义。

Reporting summary给出的[作者掩膜目录](https://github.com/MaurizioSicorello/NIDER_project/tree/main/data/masks)已核查。
其中 `morawetzMA/` 包含四个网络：

```text
EmotionRegulation_Cluster_1_of_4_cFWE05_001_104.nii
EmotionRegulation_Cluster_2_of_4_cFWE05_001_205.nii
EmotionRegulation_Cluster_3_of_4_cFWE05_001_118.nii
EmotionRegulation_Cluster_4_of_4_cFWE05_001_114.nii
```

父目录还有 `morawetz_oneTwoComb.nii`。文件名/目录来自作者仓库，尚未对本工程的参考图逐体素验证，不能据名称宣称两者完全一致。

## 2. Buhle (2014) 与 Doré (2017)：重评相关 PFC 的功能定义

Buhle et al. (2014), *Cognitive Reappraisal of Emotion: A Meta-Analysis of Human Neuroimaging Studies*.
[作者提供PDF](https://canlab.yale.edu/sites/default/files/Buhle_2014_Emo_Reg_Meta_Analysis.pdf)，
DOI：10.1093/cercor/bht154；48项研究、116个对比，使用MKDA识别重评相关一致激活区。
主要前额叶结果包括dlPFC、vlPFC、dmPFC，同时存在后部非前额叶结果。
这是一张功能性元分析图，不能把整张图叫PFC，也不能把未通过元分析的解剖PFC排除后称剩余为非PFC。

Doré, Weber & Ochsner (2017), *Neural Predictors of Decisions to Cognitively Control Emotion*.
[全文](https://pmc.ncbi.nlm.nih.gov/articles/PMC5354318/)，DOI：10.1523/JNEUROSCI.2526-16.2016。
Methods → Regions of interest明确写明：

- PFC由**双侧vlPFC、双侧dlPFC及跨双侧dmPFC**三个ROI组成。
- 作者从Buhle元分析原统计图，选取重评相对自然反应达到全脑显著性的连续体素簇；不是自己数据中显著的区域。
- 这些ROI不是整个额上/中/下回，也不是单纯围绕一个峰坐标画球。
- 后续将三个PFC ROI的活动均值合成为一个PFC指标，以降低模型复杂度。
- 所预测的是后续是否选择重评，不能直接当成本项目重评成功、年龄预测的同任务验证。

对本项目的启示：可用作“经典重评相关PFC”的独立敏感性定义；PLS主模型应保留各分区为独立特征，
不能照搬其一个合成PFC均值与432维全脑模型直接比较。

## 3. Morawetz (2020)：四个功能网络，不能等同四个 PFC 子区

Morawetz et al. (2020), *Multiple large-scale neural networks underlying emotion regulation*.
[PubMed](https://pubmed.ncbi.nlm.nih.gov/32659287/)，
[作者实验室说明](https://nbclab.github.io/papers/morawetz-multiple-largescale-neural)，
DOI：10.1016/j.neubiorev.2020.07.001。

基于既有实验的元分析数据用k-means分出四个大尺度网络。前两个偏调节过程，分别更突出反应抑制/执行控制、
评价/语言过程；后两个偏情绪生成、评价和生理过程。它为上述Nature文章的网络ROI提供依据。
本核查通过原文摘要、作者实验室和Nature方法确认其用途；没有从未取到的补充表推断每个网络的完整解剖清单。

## 4. Donahue (2018)：可复现的解剖 PFC，存在保守与宽范围两套边界

Donahue et al. (2018), *Quantitative assessment of prefrontal cortex in humans relative to nonhuman primates*.
[全文](https://pmc.ncbi.nlm.nih.gov/articles/PMC5984508/)，DOI：10.1073/pnas.1721653115。
Results → Cytoarchitectural Designations…，Fig.1及其图注。

- **Conservative PFC**：按细胞构筑等依据划定运动、前运动区前方的颗粒/异颗粒额叶皮层，
  包括相应外侧、内侧及眶额区域。
- **Liberal PFC**：在保守边界上纳入部分功能上非运动/非前运动的额叶区域。
  人类新增区包括外侧 IFJa/p、i/s6-8、55b；内侧 a/p24、25、a32pr、d/p/s32、SFL；以及 pOFC。
- 该文两套边界均排除FEF，作者同时说明其他文献可能把FEF计入PFC。
- 因而不能一概规定“所有ACC都算PFC”或“所有ACC都不算PFC”，也不能用整个额叶代替PFC。
- 本文不是重评任务研究；引用它的用途是为解剖边界提供独立、公开的操作定义。

[作者BALSA人类Fig.1资源](https://balsa.wustl.edu/GrK7)基于HCP_MMP1.0/表面空间；
网页可见文件清单，但数据下载要求登录并同意相应访问条款。本次没有下载掩膜或完成与Schaefer的映射。

## 5. 映射到本工程432的判断

432 = Schaefer400皮层 + Tian32皮层下；没有小脑。PFC来源应为400皮层中的一部分。
按当前本地 `Schaefer2018_400Parcels_7Networks_order_info.txt` 的实际名称统计：

| 名称筛选 | LH | RH | 合计 |
| --- | ---: | ---: | ---: |
| 含PFC | 35 | 36 | 71 |
| 含OFC | 6 | 6 | 12 |
| PFC或OFC | 41 | 42 | 83 |

**83只是名称候选数，不能据此宣布已确定83个解剖PFC。** 额盖/岛叶混合分区、扣带、FEF及其他边界区仍须核查。
名称本身也有版本问题；不能下载新版名称表后未经编号核验覆盖既有统计结果的名称。
[CBIG作者说明](https://groups.google.com/g/cbig_users/c/6ra-LLlogUw)明确提醒：component名称较粗，
有时不能完整描述一个高分辨率parcel的空间范围；并记录过名称修正及部分版本的编号调整。

## 6. 本项目建议（属于分析设计建议，不是上述作者的方法）

主分析采用**独立图谱定义的双侧、宽范围解剖PFC**；可优先以Donahue liberal定义作为依据。
这样核心问题为：在已包含较全面PFC信息时，其他脑区是否仍提供跨站点增量预测价值。
保守解剖定义、Buhle/Doré功能PFC可作预先指定的敏感性分析，不依据哪一个结果更显著决定主定义。

实际使用前必须：

1. 取得可复现的独立PFC空间掩膜，固定来源、版本及包含的皮层区域。
2. 与本工程Schaefer400放到相同表面网格/配准空间。不能直接把MMP编号当Schaefer编号。
3. 输出每个parcel的PFC覆盖比例、半球、纳入理由和跨界标记，先看空间分布，再固定映射阈值/跨界处理。
   只拥有parcel均值时，跨界parcel不能精确拆成其PFC及非PFC部分；最终推论需保留分区层面的限制。
4. 锁定一份全432分区的PFC/非PFC清单和模型特征列；四个模型同一任务使用相同被试、分割及处理流程。
5. “PFC合集”保留所有所选parcel作为多个特征，而非一个全PFC均值。M2使用完整432特征。

不能直接拿Nature的两个完整调节网络作PFC模型，否则其内部已包含顶叶/颞叶，无法清楚回答前额叶之外的增量价值。
不能先看本样本年龄效应、行为效应或PLS权重再定义主PFC集合。
当前结论是**定义依据已核查，最终空间掩膜及Schaefer成员表尚待制作**；没有声称已完成解剖分组。
