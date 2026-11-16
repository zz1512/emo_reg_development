**SIMULATED DATA ONLY — do not use these values in a manuscript.**

Common N = 96; sites = 4. Both outcomes are modeled in one PLS2 fit. Gender and meanFD adjustment and X/Y scaling are fitted only on training rows. Parent education is neither a covariate nor an eligibility criterion. Site labels group inner/outer LOSO and restrict permutation; no test-site regression is fitted.


## Joint age and success prediction

| Model | N | Regions | R2 | MAE | RMSE |
|---|---|---|---|---|---|
| Covariates | 96 | 0 | -0.000 | 0.861 | 1.108 |
| PFC | 96 | 80 | 0.344 | 0.687 | 0.898 |
| Distributed | 96 | 432 | 0.605 | 0.512 | 0.697 |
| Non-PFC | 96 | 352 | 0.605 | 0.511 | 0.697 |

R2 is the equally weighted mean of the two R2 values for cross-fitted, covariate-adjusted, training-standardized outcomes. MAE/RMSE pool both standardized outcomes. This is joint prediction, not pure shared variance.

| Comparison | Delta R2 | P | Primary q | Inference status |
|---|---|---|---|---|
| Distributed vs PFC | 0.261 | 0.667 | 0.667 | complete_approximate_reduced_model_permutation |
| PFC vs Covariates | 0.344 | NA | NA | not_run |
| Distributed vs Covariates | 0.606 | NA | NA | not_run |
| Non-PFC vs Covariates | 0.605 | NA | NA | not_run |
| Distributed vs Non-PFC | 0.000 | NA | NA | not_run |

The distributed-minus-PFC difference was 0.261 (permutation P = 0.667, q = 0.667); the conditional permutation test did not support a positive incremental contribution after the prespecified correction.


## Common reference behavior score prediction

| Model | N | Regions | R2 | MAE | RMSE |
|---|---|---|---|---|---|
| Covariates | 96 | 0 | -0.000 | 0.913 | 1.115 |
| PFC | 96 | 80 | 0.446 | 0.633 | 0.830 |
| Distributed | 96 | 432 | 0.779 | 0.411 | 0.525 |
| Non-PFC | 96 | 352 | 0.778 | 0.411 | 0.525 |

The target is the first Y direction of a full432 reference PLS fitted within each training fold. Every model uses the same held-out target in that fold. Target construction is repeated inside inner CV. The score is training-unit-variance and age-oriented; it is brain-informed and not an independent behavioral factor.

| Comparison | Delta R2 | P | Primary q | Inference status |
|---|---|---|---|---|
| Distributed vs PFC | 0.333 | 0.333 | 0.667 | complete_approximate_reduced_model_permutation |
| PFC vs Covariates | 0.446 | NA | NA | not_run |
| Distributed vs Covariates | 0.779 | NA | NA | not_run |
| Non-PFC vs Covariates | 0.779 | NA | NA | not_run |
| Distributed vs Non-PFC | 0.000 | NA | NA | not_run |

The distributed-minus-PFC difference was 0.333 (permutation P = 0.333, q = 0.667); the conditional permutation test did not support a positive incremental contribution after the prespecified correction.

Permutation inference uses an approximate Freedman-Lane reduced linear model with site-restricted paired bivariate residual permutations. The PFC null includes all PFC features, gender, meanFD and categorical site. It preserves modeled PFC associations and tests additional conditional information; it is not a generic test of equality of prediction algorithms. Two primary route comparisons are jointly BH-corrected. Other comparisons are descriptive unless their corresponding null is separately run. Bootstrap is an independent optional branch and is never nested in permutation. Negative R2 is retained. A shared score may be dominated by a single behavior; inspect reference_axes.csv.