# ERQ审稿回复结果汇总

**SIMULATED DATA ONLY：以下均为模拟结果，不能用于论文填写。**

## 回复信主表

建议4行数据，不含表头；每范围一行，同时列两项主要比较。

| Analysis scope | N | ERQ main q range | Sig/valid | Interaction q range | Sig/valid |
|---|---|---|---|---|---|
| Whole-brain parcels (432) | 998 to 1000 | 0.600 to 0.800 | 0/432 | 0.600 to 0.800 | 0/432 |
| Yeo networks (7) | 998 to 1000 | 0.600 to 0.800 | 0/7 | 0.600 to 0.800 | 0/7 |
| Systems (4) | 998 to 1000 | 0.600 to 0.800 | 0/4 | 0.600 to 0.800 | 0/4 |
| Amygdala subdivisions (6) | 998 to 1000 | 0.600 to 0.800 | 0/6 | 0.600 to 0.800 | 0/6 |

联合比较另见erq_joint_summary.csv（4行），不能称为单独ERQ主效应。
主效应和交互明细898行；联合明细449行。失败单元保留，不能填成0效应。

## 数值总结

恒定ERQ关联的449项有效检验中，0项通过BH校正（最小q=0.600）。
年龄×ERQ交互的449项有效检验中，0项通过BH校正（最小q=0.600）。
ERQ与交互的联合贡献的449项有效检验中，0项通过BH校正（最小q=0.600）。

## 英文结果段落

BH correction was performed separately for each comparison within the parcel432, Yeo7, system4 and Juelich6 families. For constant ERQ associations, none of the 449 valid tests survived BH correction (minimum q = 0.600). For age by ERQ interactions, none of the 449 valid tests survived BH correction (minimum q = 0.600). For joint ERQ and interaction contributions, none of the 449 valid tests survived BH correction (minimum q = 0.600). Nonsignificant results do not establish absence of associations, developmental independence, or mediation of age-related neural differences.

## 填写与解释

不显著仍报告有效单元数、N、q范围和显著数；不将全部449单元平均成一个beta、F或p。
详细beta、SE、95%CI、t/F、自由度、p、q及partial R²保存在两份明细表，不因不显著而删除。
主效应beta来自E2；平滑交互没有单个beta。q区间是跨单元的最小到最大值，不是置信区间。
若有无效检验，仅能说有效单元中未发现校正后显著结果；不对无效单元作阴性结论。
只读已有结果，不重新拟合、不重新BH；未校正p<0.05不等于校正后显著。

E1 includes smooth age and gender, categorical site, parental education and mean FD. E2 adds centered ERQ; E3 adds the smooth ERQ by age interaction. The main comparison is E2 versus E1; the interaction comparison is E3 versus E2. Each comparison was BH corrected separately within each family. Sig/valid gives q < .05 tests divided by valid tests; invalid results are not nonsignificant. N is the common sample size or its range across valid units. Ranges summarize units, not confidence intervals. Display values are rounded; counts use unrounded q. Coefficients and inference use the fixed-spline OLS reference, not the independently GCV-selected penalized fits.

完整来源、哈希、行数和校验状态见summary_metadata.json。
