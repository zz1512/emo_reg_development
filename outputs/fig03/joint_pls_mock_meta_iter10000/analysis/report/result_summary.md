**SIMULATED DATA ONLY — do not use these values in a manuscript.**

PFC definition: Buhle2014_reappraisal_PFC; 10 PFC regions, 422 remaining regions.

Membership rule: `schaefer2mm_voxel_fraction>=0.25_inside_buhle_z2.3_and_harvard_oxford_pfc`.

Remaining regions are the complement of the selected meta-PFC set; they may include anatomical prefrontal regions outside the meta-analysis mask.

Common N = 96; sites = 4. Both outcomes are modeled in one PLS2 fit. Gender and meanFD adjustment and X/Y scaling are fitted only on training rows. Parent education is neither a covariate nor an eligibility criterion. Site labels group inner/outer LOSO and restrict permutation; no test-site regression is fitted.

All six signed pairwise comparisons are displayed in this order: M1-M0, M2-M0, M3-M0, M2-M1, M3-M1 and M3-M2. Joint and shared routes have separate figures. The PFC comparisons are M2 minus M1 (all regions versus PFC) and M3 minus M1 (remaining regions versus PFC). The first comparison evaluates adding remaining regions to PFC; the second contrasts two separate feature sets and does not itself establish incremental information after combining them. The new contrast is recovered within each original paired repeat as (M3-M0)-(M1-M0); neither model fits nor P values from different contrasts are reused as if equivalent.


## Joint age and success prediction

| Model | N | Regions | R2 | MAE | RMSE |
|---|---|---|---|---|---|
| Covariates | 96 | 0 | 0.000 | 0.861 | 1.108 |
| Meta-PFC | 96 | 10 | 0.096 | 0.820 | 1.053 |
| Distributed | 96 | 432 | 0.605 | 0.515 | 0.697 |
| Remaining regions | 96 | 422 | 0.605 | 0.515 | 0.697 |

R2 is the equally weighted mean of the two R2 values for cross-fitted, covariate-adjusted, training-standardized outcomes. MAE/RMSE pool both standardized outcomes. This is joint prediction, not pure shared variance.

| Comparison | Delta R2 | P | BH q | Inference status |
|---|---|---|---|---|
| Meta-PFC vs Covariates | 0.096 | 0.25 | 0.25 | complete_approximate_reduced_model_permutation |
| Distributed vs Covariates | 0.605 | 0.25 | 0.25 | complete_approximate_reduced_model_permutation |
| Remaining regions vs Covariates | 0.605 | 0.25 | 0.25 | complete_approximate_reduced_model_permutation |
| Distributed vs Meta-PFC | 0.509 | 0.25 | 0.25 | complete_approximate_reduced_model_permutation |
| Remaining regions vs Meta-PFC | 0.508 | 0.25 | 0.25 | complete_approximate_reduced_model_permutation |
| Remaining regions vs Distributed | 0.000 | pending | — | descriptive_only_no_matching_null |

The distributed-minus-PFC difference was 0.509 (permutation P = 0.25, q = 0.25); the conditional permutation test did not support a positive incremental contribution after the prespecified correction. The remaining-regions-minus-PFC difference was 0.508 (conditional permutation P = 0.25, secondary BH q = 0.25).


## Common reference behavior score prediction

| Model | N | Regions | R2 | MAE | RMSE |
|---|---|---|---|---|---|
| Covariates | 96 | 0 | 0.000 | 0.912 | 1.116 |
| Meta-PFC | 96 | 10 | 0.125 | 0.841 | 1.044 |
| Distributed | 96 | 432 | 0.779 | 0.411 | 0.524 |
| Remaining regions | 96 | 422 | 0.778 | 0.410 | 0.525 |

The target is the first Y direction of a full432 reference PLS fitted within each training fold. Every model uses the same held-out target in that fold. Target construction is repeated inside inner CV. The score is training-unit-variance and age-oriented; it is brain-informed and not an independent behavioral factor.

| Comparison | Delta R2 | P | BH q | Inference status |
|---|---|---|---|---|
| Meta-PFC vs Covariates | 0.125 | 0.25 | 0.25 | complete_approximate_reduced_model_permutation |
| Distributed vs Covariates | 0.779 | 0.25 | 0.25 | complete_approximate_reduced_model_permutation |
| Remaining regions vs Covariates | 0.779 | 0.25 | 0.25 | complete_approximate_reduced_model_permutation |
| Distributed vs Meta-PFC | 0.654 | 0.25 | 0.25 | complete_approximate_reduced_model_permutation |
| Remaining regions vs Meta-PFC | 0.654 | 0.25 | 0.25 | complete_approximate_reduced_model_permutation |
| Remaining regions vs Distributed | 0.000 | pending | — | descriptive_only_no_matching_null |

The distributed-minus-PFC difference was 0.654 (permutation P = 0.25, q = 0.25); the conditional permutation test did not support a positive incremental contribution after the prespecified correction. The remaining-regions-minus-PFC difference was 0.654 (conditional permutation P = 0.25, secondary BH q = 0.25).

Permutation inference uses an approximate Freedman-Lane reduced linear model with site-restricted paired bivariate residual permutations. The PFC null includes all PFC features, gender, meanFD and categorical site. It preserves modeled PFC associations and tests additional conditional information; it is not a generic test of equality of prediction algorithms. Two primary M2-minus-M1 route comparisons are jointly BH-corrected. The two newly requested M3-minus-M1 route comparisons form a separate secondary BH family. The three model-minus-M0 contrasts across both routes form an additional exploratory BH family of six tests; these require the covariates-null permutations. This is not a combined ten-test correction. The q_bh column reports the appropriate within-family adjusted P. M3-M2 is descriptive, retaining its signed difference and paired bootstrap interval without borrowing a reverse-direction P value. Bootstrap is an independent optional branch and is never nested in permutation. Negative R2 is retained. A shared score may be dominated by a single behavior; inspect reference_axes.csv.