**SIMULATED DATA ONLY：仅为模拟，不可填写论文。**

# 关注脑区的逐区贝叶斯结果

| Parcel | Region | N | Age spline vs. no age | Age spline vs. linear age | Linear age vs. no age |
|---|---|---|---|---|---|
| 152 | SIMULATED parcel 152 | 120 | 6.317 | 0.365 | 17.314 |
| 154 | SIMULATED parcel 154 | 120 | 265.711 | 0.327 | 812.185 |
| 182 | SIMULATED parcel 182 | 120 | 2.036 | 0.458 | 4.448 |
| 185 | SIMULATED parcel 185 | 120 | 41.086 | 0.382 | 107.485 |
| 156 | SIMULATED parcel 156 | 120 | 1.632 | 2.804 | 0.582 |
| 157 | SIMULATED parcel 157 | 120 | 0.706 | 0.483 | 1.461 |

## 英文结果段落

We examined Bayesian evidence for age associations in 6 parcels highlighted by the preceding regional analyses (parcel IDs: 152, 154, 182, 185, 156, 157; N = 120 across valid parcels). When comparing the age spline model with the model without age terms, evidence favored the age spline model in parcels 152 (BF = 6.317), 154 (BF = 265.711), 185 (BF = 41.086); evidence was inconclusive in parcels 182 (BF = 2.036), 156 (BF = 1.632), 157 (BF = 0.706). When comparing the age spline model with the linear age model, evidence favored the linear age model in parcels 154 (BF = 0.327); evidence was inconclusive in parcels 152 (BF = 0.365), 182 (BF = 0.458), 185 (BF = 0.382), 156 (BF = 2.804), 157 (BF = 0.483). When comparing the linear age model with the model without age terms, evidence favored the linear age model in parcels 152 (BF = 17.314), 154 (BF = 812.185), 182 (BF = 4.448), 185 (BF = 107.485); evidence was inconclusive in parcels 156 (BF = 0.582), 157 (BF = 1.461). These results are follow-up characterizations of selected parcels, rather than an independent whole-brain confirmation.

## 读取说明

Entries are Bayes factors for the first model relative to the second. All models include the same covariates. The spline age model allows both linear and nonlinear associations. BF >= 3 favors the first model; BF <= 1/3 favors the second; intermediate values are inconclusive. Main prior scales are 0.5/0.5. Selected reporting does not change priors or recompute any BF. These parcels were highlighted in earlier analyses; selection is not an independent confirmation. Invalid fits are displayed as NA. The sensitivity ranges in CSV summarize priors, not confidence intervals. See log BF if a linear-scale value overflows.

本表为既有432版结果的摘选，不重新计算BF，不在6区内重新校正。不把156、157自动称为显著或边缘显著，不把其余4区预设为BF支持年龄关联。具体支持方向只由各区真实BF决定；BF支持线性相对样条，不等于没有年龄关联。
