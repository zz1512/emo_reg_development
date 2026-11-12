"""Independent statistical and cached-input checks for the smooth ERQ model.

Run from the project root with:
    python -m unittest analysis.fig03.gam_pygam.test_erq -v

The statistical oracle is statsmodels OLS, rather than pyGAM's approximate
term tests. All input tables and responses used here are synthetic.
"""

from __future__ import annotations

import contextlib
import importlib.util
import io
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

import numpy as np
import pandas as pd
from scipy.stats import false_discovery_control
import statsmodels.api as sm

from analysis.fig03.gam_nested.amygdala import AMYGDALA_ORDER
from analysis.fig03.gam_pygam import run_erq as runner
from analysis.fig03.gam_pygam.run_pygam import SYSTEM_ORDER, YEO_ORDER


TESTS = ("joint_erq", "joint_age_erq", "age", "erq", "interaction", "additive_vs_baseline")


def _oracle(y, age, erq, covariates):
    """Use an independent regression library on exactly the complete rows."""
    complete = (np.isfinite(y) & np.isfinite(age) & np.isfinite(erq)
                & np.isfinite(covariates).all(axis=1))
    y, age, erq, covariates = (values[complete] for values in (y, age, erq, covariates))
    # Five cubic B-splines span cubics plus one truncated cubic, with the
    # sole internal knot halfway between observed age extremes. This is an
    # independent parameterization of the production inference reference.
    x = (age - (age.min() + age.max()) / 2) / ((age.max() - age.min()) / 2)
    basis = np.column_stack([x, x ** 2, x ** 3, np.maximum(x, 0) ** 3])
    basis -= basis.mean(axis=0)
    erq = erq - erq.mean()
    covariates = covariates - covariates.mean(axis=0)
    full_x = sm.add_constant(np.column_stack([basis, erq, erq[:, None] * basis, covariates]),
                             has_constant="add")
    age_x = sm.add_constant(np.column_stack([basis, covariates]), has_constant="add")
    covariate_x = sm.add_constant(covariates, has_constant="add")
    full = sm.OLS(y, full_x).fit()
    reduced_age = sm.OLS(y, age_x).fit()
    reduced_covariates = sm.OLS(y, covariate_x).fit()
    return complete, full, reduced_age, reduced_covariates


def _write_cache(input_dir: Path, n_subjects: int = 48, real_ids: bool = False) -> tuple[pd.DataFrame, pd.DataFrame]:
    """A genuine 449-unit table contract; only parcel 1 has an ERQ signal."""
    input_dir.mkdir(parents=True)
    rng = np.random.default_rng(20261104)
    identifiers = [f"r{index:05d}" if real_ids else f"sub-{index:03d}" for index in range(n_subjects)]
    age = np.linspace(6.2, 18.7, n_subjects)
    erq = 3 + .03 * (age - age.mean()) + rng.normal(0, .6, n_subjects)
    design = pd.DataFrame({
        "sub_id": identifiers, "age": age, "gender": np.arange(n_subjects) % 2,
        "site_id": np.array(["1", "9", "101"])[np.arange(n_subjects) % 3],
        "parent_edu": rng.integers(1, 6, n_subjects),
        "meanFD": rng.uniform(.02, .18, n_subjects),
        "in_analysis": [1] * (n_subjects - 1) + [0],
        "excluded_by": [""] * (n_subjects - 1) + ["previous_exclusion"],
    })
    design.to_csv(input_dir / "design.csv", index=False)
    signal = .01 * (age - age.mean()) + 1.5 * (erq - erq.mean()) + rng.normal(0, .06, n_subjects)
    parcel_values = {str(index): np.full(n_subjects, float(index)) for index in range(1, 433)}
    parcel_values["1"] = signal
    pd.DataFrame({"sub_id": identifiers, **parcel_values}).to_csv(
        input_dir / "parcel_response.csv", index=False)
    pd.DataFrame({"parcel": range(1, 433), "name": [f"synthetic_{index}" for index in range(1, 433)]}).to_csv(
        input_dir / "parcel_info.csv", index=False)
    rows = []
    for dimension, names in (("yeo", YEO_ORDER), ("system", SYSTEM_ORDER)):
        for index, name in enumerate(names):
            rows.extend({"sub_id": identifier, "dimension": dimension, "group": name, "y": 20. + index}
                        for identifier in identifiers)
    pd.DataFrame(rows).to_csv(input_dir / "group_response.csv", index=False)
    pd.DataFrame({"sub_id": identifiers,
                  **{name: np.full(n_subjects, 30. + index) for index, name in enumerate(AMYGDALA_ORDER)}}).to_csv(
        input_dir / "amygdala_response.csv", index=False)
    pd.DataFrame({"roi": AMYGDALA_ORDER}).to_csv(input_dir / "amygdala_order.csv", index=False)
    (input_dir / "combined_atlas.nii.gz").write_bytes(b"not loaded by ERQ statistics")
    (input_dir / "tian_combined.json").write_text(json.dumps({
        "n_cortical": 400, "n_parcels": 432, "combined_atlas": "combined_atlas.nii.gz",
        "marginal_parcels": [156, 157], "amygdala_parcels": [403, 404, 419, 420],
    }), encoding="utf-8")
    # Deliberately conflicting profile age/sex must never replace imaging covariates.
    profiles = pd.DataFrame({
        "subject": [f"{index:03d}.0" for index in range(n_subjects)],
        "ERQ_CR": erq, "age": np.full(n_subjects, 99.), "sex": np.full(n_subjects, 7),
    })
    return design, profiles


@unittest.skipUnless(importlib.util.find_spec("pygam") is not None, "pygam is not installed")
class SmoothInferenceTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        rng = np.random.default_rng(54625)
        cls.age = np.linspace(6.2, 18.8, 240)
        cls.erq = 3 + .06 * (cls.age - cls.age.mean()) + rng.normal(0, .8, len(cls.age))
        cls.covariates = rng.normal(size=(len(cls.age), 5))
        cls.noise = rng.normal(0, .25, len(cls.age))

    def fit(self, y, age=None, erq=None, covariates=None):
        return runner.fit_unit(y, self.age if age is None else age,
                               self.erq if erq is None else erq,
                               self.covariates if covariates is None else covariates,
                               lam_grid=np.array([1e-4, 1e4]))

    def test_coefficients_intervals_and_joint_f_match_independent_ols(self):
        age_c, erq_c = self.age - self.age.mean(), self.erq - self.erq.mean()
        y = (.4 + .03 * age_c + .09 * erq_c + .02 * age_c * erq_c
             + self.covariates @ np.array([.2, -.1, .1, .05, -.03]) + self.noise)
        result, details = self.fit(y)
        complete, full, age_model, covariate_model = _oracle(y, self.age, self.erq, self.covariates)
        self.assertEqual(result["status"], "ok")
        self.assertIsNotNone(details)
        self.assertEqual(result["n_obs"], int(complete.sum()))
        self.assertAlmostEqual(result["rss_reference_full"], full.ssr, places=10)
        self.assertGreaterEqual(result["rss_full"] + 1e-8, full.ssr)
        self.assertAlmostEqual(result["beta_erq_reference"], full.params[5], places=10)
        self.assertAlmostEqual(result["se_erq_reference"], full.bse[5], places=10)
        self.assertAlmostEqual(result["p_erq"], full.pvalues[5], places=10)
        interval = full.conf_int()[5]
        self.assertAlmostEqual(result["ci_low_erq_reference"], interval[0], places=10)
        self.assertAlmostEqual(result["ci_high_erq_reference"], interval[1], places=10)
        for name, reduced in (("joint_erq", age_model), ("joint_age_erq", covariate_model)):
            with self.subTest(comparison=name):
                f_value, p_value, df1 = full.compare_f_test(reduced)
                self.assertAlmostEqual(result[f"F_{name}"], f_value, places=8)
                self.assertAlmostEqual(result[f"p_{name}"], p_value, places=10)
                self.assertEqual(result[f"df1_{name}"], df1)
                self.assertEqual(result[f"df2_{name}"], full.df_resid)
                self.assertAlmostEqual(result[f"partial_r2_reference_{name}"], 1 - full.ssr / reduced.ssr)
        additive = sm.OLS(y, np.delete(full.model.exog, np.s_[6:10], axis=1)).fit()
        additive_f, additive_p, additive_df = additive.compare_f_test(age_model)
        self.assertAlmostEqual(result["F_additive_vs_baseline"], additive_f, places=8)
        self.assertAlmostEqual(result["p_additive_vs_baseline"], additive_p, places=10)
        self.assertEqual(result["df1_additive_vs_baseline"], additive_df)
        self.assertEqual(result["df2_additive_vs_baseline"], additive.df_resid)
        self.assertAlmostEqual(result["beta_erq_additive_reference"], additive.params[5], places=10)
        self.assertAlmostEqual(result["se_erq_additive_reference"], additive.bse[5], places=10)
        np.testing.assert_allclose([result["ci_low_erq_additive_reference"], result["ci_high_erq_additive_reference"]],
                                   additive.conf_int()[5], atol=1e-10)
        self.assertAlmostEqual(result["partial_r2_reference_additive_vs_baseline"], 1 - additive.ssr / age_model.ssr)
        self.assertEqual(result["n_obs_additive"], result["n_obs_full"])
        self.assertAlmostEqual(result["rss_reference_additive"], additive.ssr, places=10)
        for field in runner.COMPARISON_FIELDS:
            self.assertEqual(result[f"{field}_full_vs_additive"], result[f"{field}_interaction"])
        slopes = runner.reference_slopes("parcel", "1", details)
        self.assertTrue(slopes["age"].between(self.age.min(), self.age.max()).all())
        observed_x = (self.age - (self.age.min() + self.age.max()) / 2) / ((self.age.max() - self.age.min()) / 2)
        original_mean = np.column_stack([observed_x, observed_x ** 2, observed_x ** 3, np.maximum(observed_x, 0) ** 3]).mean(axis=0)
        for row in slopes.itertuples():
            x = (row.age - (self.age.min() + self.age.max()) / 2) / ((self.age.max() - self.age.min()) / 2)
            c = np.zeros(len(full.params))
            c[5] = 1
            c[6:10] = np.array([x, x ** 2, x ** 3, max(x, 0) ** 3]) - original_mean
            expected = full.t_test(c)
            self.assertAlmostEqual(row.beta_erq, float(np.asarray(expected.effect).item()), places=9)
            self.assertAlmostEqual(row.se, float(np.asarray(expected.sd).item()), places=9)
            np.testing.assert_allclose([row.ci_low, row.ci_high], expected.conf_int()[0], atol=1e-9)
        no_age = sm.OLS(y, np.delete(full.model.exog, np.s_[1:5], axis=1)).fit()
        for name, reduced in (("interaction", additive), ("age", no_age)):
            expected_f, expected_p, df1 = full.compare_f_test(reduced)
            self.assertAlmostEqual(result[f"F_{name}"], expected_f, places=8)
            self.assertAlmostEqual(result[f"p_{name}"], expected_p, places=10)
            self.assertEqual(result[f"df1_{name}"], df1)
        self.assertEqual(result["df1_joint_erq"], 5)
        self.assertEqual(result["df1_joint_age_erq"], 9)
        self.assertEqual(result["df1_interaction"], 4)
        self.assertAlmostEqual(result["rss_reference_baseline"], age_model.ssr, places=10)
        self.assertEqual(result["rank_baseline_reference"], np.linalg.matrix_rank(age_model.model.exog))
        self.assertEqual(result["rank_full_reference"], np.linalg.matrix_rank(full.model.exog))
        for field in runner.COMPARISON_FIELDS:
            self.assertEqual(result[f"{field}_full_vs_baseline"], result[f"{field}_joint_erq"])
        self.assertAlmostEqual(result["delta_rss_full_vs_baseline"], result["rss_baseline"] - result["rss_full"])
        self.assertAlmostEqual(result["delta_r2_full_vs_baseline"], result["r2_full"] - result["r2_baseline"])

    def test_baseline_matches_independent_gcv_fit_on_the_common_sample(self):
        from pygam import LinearGAM, l, s

        age, erq, covariates = (value.copy() for value in (self.age, self.erq, self.covariates))
        y = .3 * age + .5 * erq + .3 * covariates[:, 0] + self.noise
        y[2], erq[4], covariates[6, 0] = np.nan, np.nan, np.nan
        result, _ = self.fit(y, age, erq, covariates)
        complete = np.isfinite(y) & np.isfinite(erq) & np.isfinite(covariates).all(axis=1)
        ac = age[complete] - age[complete].mean()
        nuisance = covariates[complete] - covariates[complete].mean(axis=0)
        baseline_x = np.column_stack([ac, nuisance])
        terms = s(0, n_splines=5)
        for feature in range(1, baseline_x.shape[1]):
            terms += l(feature, penalties=None)
        expected = LinearGAM(terms, max_iter=1000).gridsearch(
            baseline_x, y[complete], progress=False, objective="GCV",
            lam=np.array([[penalty, *([0] * nuisance.shape[1])] for penalty in (1e-4, 1e4)]))
        expected_rss = np.sum((y[complete] - expected.predict(baseline_x)) ** 2)
        self.assertEqual(result["n_obs_baseline"], int(complete.sum()))
        self.assertEqual(result["n_obs_baseline"], result["n_obs_full"])
        self.assertAlmostEqual(result["rss_baseline"], expected_rss, places=10)
        self.assertAlmostEqual(result["edf_baseline"], expected.statistics_["edof"], places=10)
        self.assertAlmostEqual(result["gcv_baseline"], expected.statistics_["GCV"], places=10)
        self.assertEqual(result["lambda_age_baseline"], expected.lam[0][0])

    def test_baseline_is_unchanged_when_only_observed_erq_values_change(self):
        y = .05 * self.age + .8 * self.erq + self.noise
        first, _ = self.fit(y)
        changed_erq = np.random.default_rng(19).permutation(self.erq)
        second, _ = self.fit(y, erq=changed_erq)
        for column in ("rss_baseline", "r2_baseline", "edf_baseline", "gcv_baseline", "lambda_age_baseline"):
            self.assertAlmostEqual(first[column], second[column], places=12)
        self.assertGreater(abs(first["rss_full"] - second["rss_full"]), 1)

    def _orthogonal_noise(self):
        _, full, _, _ = _oracle(self.noise, self.age, self.erq, self.covariates)
        return full.resid

    def test_erq_main_effect_can_pass_joint_test_without_interaction(self):
        y = (.1 * (self.age - self.age.mean()) + 1.2 * (self.erq - self.erq.mean())
             + self.covariates[:, 0] * .2 + self._orthogonal_noise())
        result, _ = self.fit(y)
        self.assertLess(result["p_joint_erq"], 1e-20)
        self.assertLess(result["p_erq"], 1e-20)
        self.assertGreater(result["p_interaction"], .9)
        self.assertLess(result["p_additive_vs_baseline"], 1e-20)

    def test_pure_age_signal_does_not_open_the_erq_followup_gate(self):
        y = .8 * (self.age - self.age.mean()) + self._orthogonal_noise()
        result, _ = self.fit(y)
        self.assertLess(result["p_joint_age_erq"], 1e-20)
        self.assertGreater(result["p_joint_erq"], .5)

    def test_interaction_is_tested_and_selected_when_the_additive_step_is_null(self):
        _, full, _, _ = _oracle(self.noise, self.age, self.erq, self.covariates)
        additive_x = np.delete(full.model.exog, np.s_[6:10], axis=1)
        signal = full.model.exog[:, 7]
        signal -= additive_x @ np.linalg.lstsq(additive_x, signal, rcond=None)[0]
        y = signal * 3 + self._orthogonal_noise() * .1
        result, _ = self.fit(y)
        self.assertGreater(result["p_additive_vs_baseline"], .9)
        self.assertLess(result["p_full_vs_additive"], 1e-20)
        result.update(dimension="parcel", name="1", label="synthetic")
        corrected = runner.correct_families(pd.DataFrame([result]))
        self.assertEqual(len(runner.selected_units(corrected)), 1)

    def test_all_comparisons_use_the_same_complete_rows(self):
        y = .2 * self.age + .4 * self.erq + self.noise
        y, age, erq, covariates = (values.copy() for values in (y, self.age, self.erq, self.covariates))
        y[2], age[5], erq[9], covariates[12, 0] = np.nan, np.nan, np.inf, np.nan
        result, _ = self.fit(y, age, erq, covariates)
        complete, full, reduced, _ = _oracle(y, age, erq, covariates)
        self.assertEqual(result["n_obs"], int(complete.sum()))
        self.assertEqual(result["n_obs_baseline"], int(complete.sum()))
        self.assertEqual(result["n_obs_full"], int(complete.sum()))
        self.assertAlmostEqual(result["beta_erq_reference"], full.params[5], places=10)
        self.assertAlmostEqual(result["F_joint_erq"], full.compare_f_test(reduced)[0], places=8)

    def test_nonidentifiable_or_constant_units_do_not_receive_significance(self):
        for y, age, erq in ((self.noise, self.age, 2 * self.age),
                            (self.noise, self.age, np.ones(len(self.age))),
                            (self.noise, np.ones(len(self.age)), self.erq),
                            (np.ones(len(self.age)), self.age, self.erq)):
            with self.subTest(constant_response=np.ptp(y) == 0, constant_erq=np.ptp(erq) == 0):
                result, details = self.fit(y, age=age, erq=erq)
                self.assertNotEqual(result["status"], "ok")
                self.assertTrue(np.isnan(result["p_joint_erq"]))
                self.assertIsNone(details)

    def test_varying_erq_slope_can_curve_and_predictions_use_one_shared_reference(self):
        x = (self.age - 12.5) / 6.3
        ec = self.erq - self.erq.mean()
        y = .1 * x + .1 * x ** 2 + ec * (.7 + .8 * x ** 2) + self.noise * .15
        result, details = self.fit(y)
        self.assertLess(result["p_interaction"], 1e-20)
        curves, slopes = runner.conditional_predictions("parcel", "1", details)
        self.assertGreater(slopes["erq_slope"].iloc[0], slopes["erq_slope"].iloc[len(slopes) // 2] + .4)
        model = details["model"]
        for row in curves.iloc[::35].itertuples():
            at = np.zeros((1, 2 + details["n_covariates"]))
            at[0, :2] = row.age - details["age_center"], row.erq_value - details["erq_center"]
            reference_rows = np.zeros((len(details["age"]), at.shape[1]))
            reference_rows[:, 0] = details["age"] - details["age_center"]
            reference = model.predict(reference_rows).mean()
            expected_mean = model.predict(at)[0] - reference
            expected_ci = model.confidence_intervals(at)[0] - reference
            self.assertAlmostEqual(row.effect, expected_mean, places=10)
            np.testing.assert_allclose([row.ci_low, row.ci_high], expected_ci, atol=1e-10)
        low = curves[curves["erq_level"] == "low"].reset_index(drop=True)
        high = curves[curves["erq_level"] == "high"].reset_index(drop=True)
        difference = (high["effect"] - low["effect"]) / (high["erq_value"] - low["erq_value"])
        np.testing.assert_allclose(difference, slopes["erq_slope"], atol=1e-10)

    def test_fixed_space_tests_are_invariant_to_gcv_penalty_grid(self):
        y = .2 * self.age + .3 * self.erq + self.noise
        results = [runner.fit_unit(y, self.age, self.erq, self.covariates, lam_grid=grid)[0]
                   for grid in (np.array([1e-5, 1e-4]), np.array([1e4, 1e5]))]
        for name in TESTS:
            self.assertAlmostEqual(results[0][f"p_{name}"], results[1][f"p_{name}"], places=12)


class CorrectionAndSampleTests(unittest.TestCase):
    def test_bh_is_separate_for_each_dimension_and_each_full_family(self):
        rows = []
        for family, count in (("parcel", 4), ("yeo", 3), ("system", 2), ("amygdala", 3)):
            for index in range(count):
                row = {"dimension": family, "name": str(index)}
                row.update({f"p_{name}": (.003, .04, .2, np.nan)[(index + offset) % 4]
                            for offset, name in enumerate(TESTS)})
                rows.append(row)
        source = pd.DataFrame(rows)
        corrected = runner.correct_families(source.copy())
        np.testing.assert_allclose(corrected["p_full_vs_baseline"], corrected["p_joint_erq"], equal_nan=True)
        np.testing.assert_allclose(corrected["q_full_vs_baseline"], corrected["q_joint_erq"], equal_nan=True)
        np.testing.assert_array_equal(corrected["reject_full_vs_baseline"], corrected["reject_joint_erq"])
        for family in source["dimension"].unique():
            for name in TESTS:
                with self.subTest(family=family, test=name):
                    select = source["dimension"] == family
                    p_values = source.loc[select, f"p_{name}"].to_numpy()
                    finite = np.isfinite(p_values)
                    expected = np.full(len(p_values), np.nan)
                    expected[finite] = false_discovery_control(p_values[finite], method="bh")
                    np.testing.assert_allclose(corrected.loc[select, f"q_{name}"], expected, equal_nan=True)

    def test_erq_matching_preserves_imaging_sample_and_covariates(self):
        with tempfile.TemporaryDirectory(prefix="erq-sample-") as temporary:
            root = Path(temporary)
            design, profiles = _write_cache(root / "inputs", n_subjects=24)
            profiles = profiles.iloc[1:].copy()  # sub-000 has no ERQ profile.
            profiles.loc[profiles.index[0], "ERQ_CR"] = np.nan
            design.loc[2, "parent_edu"] = np.nan
            design.to_csv(root / "inputs/design.csv", index=False)
            profile_path = root / "profiles.csv"
            profiles.to_csv(profile_path, index=False)
            sample, roster, _, _ = runner.prepare_sample(root / "inputs", profile_path)
            self.assertEqual(len(sample), 20)
            self.assertEqual(len(roster), 24)
            self.assertNotIn("sub-023", set(sample["sub_id"]))
            self.assertNotIn("sub-000", set(sample["sub_id"]))
            self.assertLess(sample["age"].max(), 20)
            self.assertEqual(set(sample["gender"]), {0, 1})
            self.assertTrue(np.isfinite(sample["erq"]).all())
            expected = profiles.set_index(profiles["subject"].str.replace(".0", "", regex=False))["ERQ_CR"]
            for row in sample.itertuples():
                self.assertAlmostEqual(row.erq, expected.loc[row.sub_id.removeprefix("sub-")])

    def test_normalized_duplicate_erq_ids_fail_instead_of_being_averaged(self):
        with tempfile.TemporaryDirectory(prefix="erq-duplicates-") as temporary:
            root = Path(temporary)
            _, profiles = _write_cache(root / "inputs")
            profiles = pd.concat([profiles, pd.DataFrame({"subject": ["sub-001"], "ERQ_CR": [6.]})],
                                 ignore_index=True)
            profile_path = root / "profiles.csv"
            profiles.to_csv(profile_path, index=False)
            with self.assertRaisesRegex(ValueError, "(?i)duplicat|conflict"):
                runner.prepare_sample(root / "inputs", profile_path)

    def test_verified_row_ids_are_restored_before_filtering_and_scores_stay_in_new_table(self):
        with tempfile.TemporaryDirectory(prefix="erq-row-ids-") as temporary:
            root = Path(temporary)
            design, profiles = _write_cache(root / "inputs", n_subjects=24, real_ids=True)
            expected = profiles["ERQ_CR"].iloc[:-1].to_numpy()
            profiles = profiles.rename(columns={"subject": "sub_ids"})
            profiles["sub_ids"] = [f"sub-{number + 1:04d}" for number in range(24)]
            # Permute both source tables together; the final join must follow imaging IDs.
            order = np.random.default_rng(5).permutation(24)
            profile_path, id_path = root / "er_profiles_new.csv", root / "er_profiles.csv"
            profiles.iloc[order].to_csv(profile_path, index=False)
            id_path.write_text("subject\tage\tERQ-CR\n" + "\n".join(
                f'{design.iloc[index]["sub_id"]}\t9\t999' for index in order), encoding="utf-8-sig")
            sample, roster, metadata, rows = runner.prepare_sample(
                root / "inputs", profile_path, "ERQ_CR", "sub_ids", id_path)
            np.testing.assert_allclose(sample["erq"], expected)
            self.assertEqual(sample["sub_id"].tolist(), design["sub_id"].iloc[:-1].tolist())
            self.assertEqual(len(roster), 24)
            self.assertEqual(rows["erq_real_sub_id"].tolist(), design.iloc[order]["sub_id"].tolist())
            self.assertEqual(rows["erq_source_sub_id"].tolist(), profiles.iloc[order]["sub_ids"].tolist())
            self.assertEqual(rows["erq_source_row"].tolist(), list(range(1, 25)))
            self.assertTrue(metadata["id_restoration"]["enabled"])
            self.assertEqual(metadata["erq_aggregation"]["n_source_rows_retained"], 23)
            self.assertFalse(np.isclose(sample["erq"], 999).any())

    def test_row_id_restoration_rejects_count_mismatch_duplicates_and_invalid_prefix(self):
        with tempfile.TemporaryDirectory(prefix="erq-invalid-ids-") as temporary:
            path = Path(temporary) / "er_profiles.csv"
            for records, count, message in [
                (["r00001\t9", "r00002\t8"], 3, "row count mismatch"),
                (["r00001\t9", "r00001\t8"], 2, "duplicate subject IDs"),
                (["bad!\t9", "r00002\t8"], 2, "six alphanumeric"),
            ]:
                with self.subTest(message=message):
                    path.write_text("subject\tage\n" + "\n".join(records), encoding="utf-8")
                    with self.assertRaisesRegex(ValueError, message):
                        runner._rowwise_real_ids(path, count)

    def test_audit_lists_each_failed_filter_and_source_row_without_double_counting(self):
        with tempfile.TemporaryDirectory(prefix="erq-filter-audit-") as temporary:
            root = Path(temporary)
            design, profiles = _write_cache(root / "inputs", n_subjects=24)
            profiles = profiles.iloc[1:].copy()  # First imaging ID has no ERQ row.
            profiles["ERQ_CR"] = profiles["ERQ_CR"].astype(object)
            profiles.loc[1, "ERQ_CR"] = "NA"
            profiles = pd.concat([profiles, pd.DataFrame({"subject": ["outside"], "ERQ_CR": ["bad"]})], ignore_index=True)
            design.loc[[1, 2], "parent_edu"] = np.nan
            design.loc[23, "excluded_by"] = "meanFD_gt_0.5"
            design.to_csv(root / "inputs/design.csv", index=False)
            profile_path = root / "profiles.csv"
            profiles.to_csv(profile_path, index=False)
            sample, roster, metadata, rows = runner.prepare_sample(root / "inputs", profile_path)
            runner._write_sample_audits(root / "audit", roster, rows, metadata)
            self.assertEqual(len(sample), 20)
            summary = pd.read_csv(root / "audit/erq_filter_summary.csv").set_index("filter")
            self.assertEqual(summary.loc["ERQ_subject_not_found", "n_newly_excluded"], 1)
            self.assertEqual(summary.loc["erq_missing_or_nonfinite", "n_newly_excluded"], 1)
            self.assertEqual(summary.loc["parent_edu_missing_or_nonfinite", "n_newly_excluded"], 1)
            self.assertEqual(summary.loc["parent_edu_missing_or_nonfinite", "n_failed_criterion"], 2)
            self.assertEqual(summary["n_newly_excluded"].sum(), 4)
            exclusions = pd.read_csv(root / "audit/erq_exclusions_by_filter.csv")
            twice = exclusions.loc[exclusions["sub_id"] == "sub-001"]
            self.assertEqual(set(twice["filter"]), {"erq_missing_or_nonfinite", "parent_edu_missing_or_nonfinite"})
            self.assertEqual(twice["is_first_exclusion"].sum(), 1)
            audited_rows = pd.read_csv(root / "audit/erq_rows.csv", keep_default_na=False)
            missing = audited_rows.loc[audited_rows["erq_source_sub_id"] == "001.0"].iloc[0]
            self.assertEqual(missing["erq_score_raw"], "NA")
            self.assertEqual(missing["erq_score_status"], "missing")
            outside = audited_rows.loc[audited_rows["erq_source_sub_id"] == "outside"].iloc[0]
            self.assertEqual(outside["erq_excluded_by"], "not_in_imaging_roster")
            self.assertEqual(outside["erq_score_status"], "not_numeric")
            self.assertEqual(len(audited_rows), len(profiles))
            original = pd.read_csv(root / "audit/original_exclusions_by_filter.csv")
            self.assertEqual(original["sub_id"].tolist(), ["sub-023"])
            self.assertEqual(original["filter"].tolist(), ["meanFD_gt_0.5"])
            self.assertEqual(metadata["erq_aggregation"]["n_source_rows_outside_imaging"], 1)

    def test_cloud_sub_ids_and_erq_cr_columns_and_wrong_cr_hint(self):
        with tempfile.TemporaryDirectory(prefix="erq-cloud-columns-") as temporary:
            root = Path(temporary)
            _, profiles = _write_cache(root / "inputs")
            profiles = profiles.rename(columns={"subject": "sub_ids"})
            profiles["relative_reappraisal_success"] = 99.0
            profiles["reappraisal_success"] = -99.0
            profile_path = root / "er_profiles_new.csv"
            profiles.to_csv(profile_path, index=False, encoding="utf-8-sig")
            sample, _, metadata, _ = runner.prepare_sample(root / "inputs", profile_path)
            explicit, _, _, _ = runner.prepare_sample(root / "inputs", profile_path, "ERQ_CR", "sub_ids")
            pd.testing.assert_frame_equal(sample, explicit)
            self.assertEqual(metadata["erq_id_column"], "sub_ids")
            self.assertEqual(metadata["erq_column"], "ERQ_CR")
            np.testing.assert_allclose(sample["erq"], profiles["ERQ_CR"].iloc[:-1])
            with self.assertRaisesRegex(ValueError, r"score column 'CR'.*--erq-column ERQ_CR"):
                runner.prepare_sample(root / "inputs", profile_path, "CR", "sub_ids")

    def test_explicit_alternative_score_column_is_matched_without_changing_the_sample(self):
        with tempfile.TemporaryDirectory(prefix="erq-cloud-columns-") as temporary:
            root = Path(temporary)
            _, profiles = _write_cache(root / "inputs")
            cloud_profiles = profiles.rename(columns={"subject": "sub_ids", "ERQ_CR": "CR"})
            profile_path = root / "er_profiles_new.csv"
            cloud_profiles.to_csv(profile_path, index=False, encoding="utf-8-sig")
            sample, _, metadata, _ = runner.prepare_sample(root / "inputs", profile_path, score_column="CR")
            explicit, _, _, _ = runner.prepare_sample(root / "inputs", profile_path, "CR", "sub_ids")
            pd.testing.assert_frame_equal(sample, explicit)
            self.assertEqual(len(sample), 47)
            self.assertEqual(metadata["erq_id_column"], "sub_ids")
            self.assertEqual(metadata["erq_column"], "CR")
            np.testing.assert_allclose(sample["erq"], cloud_profiles["CR"].iloc[:-1])
            with self.assertRaisesRegex(ValueError, r"ERQ score column 'ERQ_CR'.*Available columns:.*'CR'.*--erq-column"):
                runner.prepare_sample(root / "inputs", profile_path)
            with self.assertRaisesRegex(ValueError, r"subject ID column 'wrong_id'.*Available columns:.*sub_ids.*--erq-id-column"):
                runner.prepare_sample(root / "inputs", profile_path, "CR", "wrong_id")

    def test_station_labels_are_categorical_not_numeric_distances(self):
        sample = pd.DataFrame({
            "gender": [0, 1, 0, 1, 0, 1], "parent_edu": [1, 3, 2, 4, 5, 1],
            "meanFD": [.03, .06, .04, .11, .07, .09],
            "site_id": ["1", "9", "101", "1", "9", "101"],
        })
        matrix, names, _ = runner.build_covariates(sample)
        self.assertEqual(matrix.shape, (6, 5))  # three numeric covariates plus two site indicators.
        self.assertEqual(len(names), 5)
        alternative = sample.copy()
        alternative["site_id"] = alternative["site_id"].map({"1": "200", "9": "2", "101": "1"})
        recoded, _, _ = runner.build_covariates(alternative)
        original_x = sm.add_constant(matrix, has_constant="add")
        recoded_x = sm.add_constant(recoded, has_constant="add")
        np.testing.assert_allclose(original_x @ np.linalg.pinv(original_x),
                                   recoded_x @ np.linalg.pinv(recoded_x), atol=1e-10)

    def test_zero_id_overlap_is_audited_and_rejected_before_fitting(self):
        with tempfile.TemporaryDirectory(prefix="erq-zero-overlap-") as temporary:
            root = Path(temporary)
            _, profiles = _write_cache(root / "source/inputs", n_subjects=24)
            profiles["subject"] = [f"sub-anonymous-{number:04d}" for number in range(len(profiles))]
            profile_path = root / "profiles.csv"
            profiles.to_csv(profile_path, index=False)
            destination = root / "result"
            with patch.object(runner, "fit_unit", side_effect=AssertionError("Must not fit")), \
                    contextlib.redirect_stdout(io.StringIO()), self.assertRaisesRegex(ValueError, "ID overlap is 0"):
                runner.run_analysis(root / "source/inputs", profile_path, destination)
            metadata = json.loads((destination / "diagnostics/erq_sample_summary.json").read_text())
            self.assertEqual(metadata["n_original_in_analysis"], 23)
            self.assertEqual(metadata["n_erq_rows"], 24)
            self.assertEqual(metadata["n_original_with_erq_id"], 0)
            self.assertEqual(metadata["exclusions"]["ERQ_subject_not_found"], 23)
            audit = pd.read_csv(destination / "diagnostics/erq_sample_roster.csv")
            self.assertEqual(len(audit), 24)
            self.assertEqual(audit["erq_subject_matched"].sum(), 0)
            self.assertEqual(audit["sub_id"].iloc[0], "sub-000")
            self.assertIn("erq_source_sub_id", audit)
            self.assertFalse((destination / "statistics").exists())

    def test_zero_complete_score_is_distinguished_from_id_mismatch(self):
        with tempfile.TemporaryDirectory(prefix="erq-empty-score-") as temporary:
            root = Path(temporary)
            _, profiles = _write_cache(root / "source/inputs", n_subjects=24)
            profiles["ERQ_CR"] = "not recorded"
            profile_path = root / "profiles.csv"
            profiles.to_csv(profile_path, index=False)
            destination = root / "result"
            with contextlib.redirect_stdout(io.StringIO()), self.assertRaisesRegex(ValueError, "Only 0 complete subjects"):
                runner.prepare_checked_sample(root / "source/inputs", profile_path, destination)
            metadata = json.loads((destination / "diagnostics/erq_sample_summary.json").read_text())
            self.assertEqual(metadata["n_original_with_erq_id"], 23)
            self.assertEqual(metadata["n_original_with_finite_erq"], 0)
            self.assertEqual(metadata["exclusions"]["erq_missing_or_nonfinite"], 23)

    def test_tab_delimited_original_profiles_pass_cli_check_in_id_order(self):
        with tempfile.TemporaryDirectory(prefix="erq-tab-check-") as temporary:
            root = Path(temporary)
            _, profiles = _write_cache(root / "source/inputs", n_subjects=24)
            expected = profiles["ERQ_CR"].iloc[:-1].to_numpy()
            profiles = profiles.rename(columns={"ERQ_CR": "ERQ-CR"}).sample(frac=1, random_state=9)
            profile_path = root / "er_profiles.csv"
            profiles.to_csv(profile_path, sep="\t", index=False)
            destination = root / "result"
            arguments = ["run_erq", "--step", "check", "--inputs", str(root / "source/inputs"),
                         "--erq", str(profile_path), "--erq-id-column", "subject", "--erq-column", "ERQ-CR",
                         "--output", str(destination)]
            with patch.object(sys, "argv", arguments), patch.object(
                    runner, "fit_unit", side_effect=AssertionError("Check must not fit")), \
                    contextlib.redirect_stdout(io.StringIO()):
                runner.main()
            audit = pd.read_csv(destination / "diagnostics/erq_sample_roster.csv")
            np.testing.assert_allclose(audit.loc[audit["in_erq_analysis"] == 1, "erq"], expected)
            self.assertFalse((destination / "statistics").exists())


@unittest.skipUnless(importlib.util.find_spec("pygam") is not None, "pygam is not installed")
class CachedPipelineTests(unittest.TestCase):
    def test_cli_covers_449_units_and_only_plots_joint_erq_discoveries(self):
        with tempfile.TemporaryDirectory(prefix="erq-pipeline-") as temporary:
            root = Path(temporary)
            design, profiles = _write_cache(root / "source/inputs", real_ids=True)
            profiles = profiles.rename(columns={"subject": "sub_ids"})
            profiles["sub_ids"] = [f"sub-{number + 1:04d}" for number in range(len(profiles))]
            profile_path = root / "profiles.csv"
            profiles.to_csv(profile_path, index=False)
            id_path = root / "er_profiles.csv"
            id_path.write_text("subject\tage\tERQ-CR\n" + "\n".join(
                f"{identifier}\t9\t999" for identifier in design["sub_id"]), encoding="utf-8")
            responses = pd.read_csv(root / "source/inputs/parcel_response.csv")
            responses.loc[[2, 9], "1"] = np.nan
            responses.to_csv(root / "source/inputs/parcel_response.csv", index=False)
            before = {path.name: path.read_bytes() for path in (root / "source/inputs").iterdir()}
            profile_before = profile_path.read_bytes()
            id_before = id_path.read_bytes()
            destination = root / "erq_result"
            destination.mkdir()
            (destination / "run_metadata.json").write_text(json.dumps({
                "analysis_kind": runner.PREVIOUS_ANALYSIS_KIND}), encoding="utf-8")
            redraw = ["run_erq", "--step", "figures", "--output", str(destination)]
            with patch.object(sys, "argv", redraw), self.assertRaisesRegex(ValueError, "rerun --step fit"):
                runner.main()
            arguments = ["run_erq", "--inputs", str(root / "source/inputs"),
                         "--erq", str(profile_path), "--erq-id-column", "sub_ids", "--erq-column", "ERQ_CR",
                         "--erq-id-source", str(id_path),
                         "--skip-docx",
                         "--output", str(destination)]
            with patch.object(sys, "argv", arguments), patch.object(runner, "LAM_GRID", np.array([1e-4, 1e4])), \
                    contextlib.redirect_stdout(io.StringIO()):
                runner.main()
            statistics = pd.read_csv(destination / "statistics/erq_age_statistics.csv")
            self.assertEqual(statistics.groupby("dimension").size().to_dict(),
                             {"parcel": 432, "yeo": 7, "system": 4, "amygdala": 6})
            comparison = pd.read_csv(destination / "statistics/erq_full_vs_baseline.csv")
            self.assertEqual(len(comparison), 449)
            manuscript = pd.read_csv(destination / "statistics/erq_sequential_comparisons.csv")
            self.assertEqual(len(manuscript), 449 * 3)
            self.assertEqual(set(manuscript["comparison"]), {"additive_vs_baseline", "full_vs_additive", "full_vs_baseline"})
            self.assertTrue(manuscript.loc[manuscript["comparison"] == "full_vs_additive", "beta"].isna().all())
            self.assertEqual(set(comparison["comparison"]), {"full_vs_baseline"})
            np.testing.assert_array_equal(comparison["n_obs_baseline"], comparison["n_obs_full"])
            np.testing.assert_allclose(comparison["q_full_vs_baseline"], statistics["q_joint_erq"], equal_nan=True)
            self.assertTrue(np.isfinite(comparison.loc[comparison["status"] == "ok", list(runner.PAIRED_MODEL_COLUMNS)]).all().all())
            primary_selected = pd.read_csv(destination / "statistics/full_vs_baseline_significant.csv")
            self.assertEqual(len(primary_selected), 1)
            selected = pd.read_csv(destination / "statistics/joint_erq_significant.csv")
            self.assertEqual(len(selected), 1)
            self.assertEqual(selected.iloc[0]["dimension"], "parcel")
            self.assertEqual(str(selected.iloc[0]["name"]), "1")
            self.assertTrue((selected["q_joint_erq"] < .05).all())
            self.assertEqual(len(list((destination / "figures").rglob("*.svg"))), 2)
            self.assertEqual(len(pd.read_csv(destination / "inputs/erq_sample_roster.csv")), 48)
            self.assertEqual(len(pd.read_csv(destination / "inputs/analysis_design.csv")), 47)
            crosswalk = pd.read_csv(destination / "inputs/erq_id_crosswalk.csv")
            self.assertEqual(crosswalk["erq_real_sub_id"].tolist(), design["sub_id"].tolist())
            self.assertEqual(crosswalk["erq_source_sub_id"].tolist(), profiles["sub_ids"].tolist())
            response_exclusions = pd.read_csv(destination / "statistics/unit_response_exclusions.csv")
            self.assertEqual(response_exclusions["sub_id"].tolist(), ["r00002", "r00009"])
            self.assertEqual(set(response_exclusions["name"].astype(str)), {"1"})
            self.assertEqual(comparison.loc[comparison["status"] == "ok", "n_obs_full"].tolist(), [45])
            self.assertTrue((destination / "run_metadata.json").is_file())
            metadata = json.loads((destination / "run_metadata.json").read_text(encoding="utf-8"))
            self.assertEqual(metadata["analysis_kind"], runner.ANALYSIS_KIND)
            self.assertEqual(metadata["primary_comparison"], "full_vs_baseline")
            self.assertEqual(metadata["sample"]["erq_id_column"], "sub_ids")
            self.assertEqual(metadata["sample"]["erq_column"], "ERQ_CR")
            self.assertTrue(metadata["sample"]["id_restoration"]["enabled"])
            self.assertEqual(metadata["sample"]["erq_aggregation"]["n_source_rows_retained"], 47)
            self.assertEqual(metadata["smoothing"]["baseline_combinations"], 2)
            self.assertEqual(metadata["smoothing"]["full_combinations"], 4)
            self.assertEqual(before, {path.name: path.read_bytes() for path in (root / "source/inputs").iterdir()})
            self.assertEqual(profile_before, profile_path.read_bytes())
            self.assertEqual(id_before, id_path.read_bytes())
            # Redrawing is self-contained: the original ERQ table need not exist.
            profile_path.unlink()
            redraw = ["run_erq", "--step", "figures", "--output", str(destination)]
            with patch.object(sys, "argv", redraw), patch.object(
                    runner, "fit_unit", side_effect=AssertionError("Redraw must not fit models")):
                runner.main()
            self.assertEqual(len(list((destination / "figures").rglob("*.svg"))), 2)
            # A changed selection removes only figures owned by the old manifest.
            unrelated = destination / "figures/user_annotation.svg"
            unrelated.write_text("user file", encoding="utf-8")
            statistics["q_joint_erq"] = 1.0
            statistics["q_full_vs_baseline"] = 1.0
            statistics["q_additive_vs_baseline"] = 1.0
            statistics["q_full_vs_additive"] = 1.0
            statistics.to_csv(destination / "statistics/erq_age_statistics.csv", index=False)
            with patch.object(sys, "argv", redraw):
                runner.main()
            self.assertEqual(list((destination / "figures").rglob("*.svg")), [unrelated])
            self.assertTrue(pd.read_csv(destination / "statistics/figure_manifest.csv").empty)
            with self.assertRaisesRegex(ValueError, "separate"):
                runner.run_analysis(root / "source/inputs", profile_path, root / "source")


if __name__ == "__main__":
    unittest.main()
