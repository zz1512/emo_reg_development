"""Independent dense-covariance oracle and synthetic age-BF pipeline checks."""

import contextlib
import io
import json
from pathlib import Path
import tempfile
import unittest

import numpy as np
import pandas as pd
from scipy.integrate import quad

from analysis.fig03.gam_pygam.age_bayes import fit_age_bayes, spline_components
from analysis.fig03.gam_pygam import run_age_bayes as runner
from analysis.fig03.gam_pygam.test_erq import _write_cache


class BayesianSplineTests(unittest.TestCase):
    def setUp(self):
        rng = np.random.default_rng(46502)
        self.age = np.linspace(6.2, 18.8, 160)
        self.covariates = rng.normal(size=(160, 2))
        self.x, self.z = spline_components(self.age)
        design = np.column_stack([np.ones(160), self.covariates, self.x, self.z])
        noise = rng.normal(size=160)
        self.noise = noise - design @ np.linalg.lstsq(design, noise, rcond=None)[0]

    def test_penalty_basis_is_centered_normalized_and_spans_original_space(self):
        np.testing.assert_allclose(self.z.mean(axis=0), 0, atol=1e-14)
        np.testing.assert_allclose(self.x @ self.z, 0, atol=1e-12)
        self.assertAlmostEqual(np.sum(self.z ** 2) / len(self.age), 1)

    def test_dense_marginal_integral_matches_low_rank_calculation(self):
        age, cov = self.age[:28], self.covariates[:28]
        x, z = spline_components(age)
        y = .3 * x + .2 * z[:, 0] + self.noise[:28]
        nuisance = np.column_stack([np.ones(len(age)), cov])
        nu = len(age) - nuisance.shape[1]

        def dense(tau, linear=True):
            v = np.eye(len(age)) + tau ** 2 * z @ z.T
            if linear:
                v += .5 ** 2 * np.outer(x, x)
            inverse = np.linalg.inv(v)
            gram = nuisance.T @ inverse @ nuisance
            residual = y - nuisance @ np.linalg.solve(gram, nuisance.T @ inverse @ y)
            return -.5 * np.linalg.slogdet(v)[1] - .5 * np.linalg.slogdet(gram)[1] - nu / 2 * np.log(residual @ inverse @ residual)

        null = dense(0, False)
        expected_m2 = dense(0) - null
        value = quad(lambda tau: np.exp(dense(tau) - null) * np.sqrt(2 / np.pi) / .5 * np.exp(-tau ** 2 / (.5 ** 2 * 2)),
                     0, np.inf, epsabs=1e-9, epsrel=1e-9)[0]
        result, _ = fit_age_bayes(y, age, cov, scales=(.5,))
        self.assertAlmostEqual(result["log_bf10_m2_vs_m1"], expected_m2, places=9)
        self.assertAlmostEqual(result["log_bf10_m3_vs_m1"], np.log(value), places=8)

    def test_null_linear_and_nonlinear_signals_have_correct_evidence_direction(self):
        null, _ = fit_age_bayes(self.noise, self.age, self.covariates, scales=(.5,))
        linear, _ = fit_age_bayes(1.5 * self.x + self.noise, self.age, self.covariates, scales=(.5,))
        curved, _ = fit_age_bayes(1.5 * self.z[:, 0] + self.noise, self.age, self.covariates, scales=(.5,))
        self.assertLess(null["bf10_m2_vs_m1"], 1 / 3)
        self.assertLess(null["bf10_m3_vs_m1"], 1 / 3)
        self.assertGreater(linear["bf10_m2_vs_m1"], 10)
        self.assertLess(linear["bf10_m3_vs_m2"], 1 / 3)
        self.assertGreater(curved["bf10_m3_vs_m2"], 10)

    def test_all_nine_priors_are_coherent_and_invariant_to_response_units_and_nuisance_shift(self):
        y = .4 * self.x + .2 * self.z[:, 0] + self.noise
        first, sensitivity = fit_age_bayes(y, self.age, self.covariates)
        second, _ = fit_age_bayes(7 * y + 10 + self.covariates @ [2, -4], self.age * 3 + 20, self.covariates)
        self.assertEqual(len(sensitivity), 9)
        for row in sensitivity:
            self.assertAlmostEqual(row["log_bf10_m3_vs_m1"] - row["log_bf10_m3_vs_m2"], row["log_bf10_m2_vs_m1"], places=12)
            self.assertAlmostEqual(row["bf10_m3_vs_m2"] * row["bf01_m3_vs_m2"], 1, places=12)
            self.assertLess(row["quadrature_check_log_difference"], 1e-6)
        for comparison in ("m3_vs_m1", "m3_vs_m2", "m2_vs_m1"):
            self.assertAlmostEqual(first[f"log_bf10_{comparison}"], second[f"log_bf10_{comparison}"], places=9)

    def test_missing_rows_and_nonidentifiability_do_not_fabricate_evidence(self):
        y = self.noise.copy()
        y[1] = np.nan
        result, _ = fit_age_bayes(y, self.age, self.covariates, scales=(.5,))
        self.assertEqual(result["n_obs"], 159)
        for response, age, cov in ((np.ones(160), self.age, self.covariates),
                                    (self.noise, np.ones(160), self.covariates),
                                    (self.noise, self.age, np.column_stack([self.age, self.covariates]))):
            invalid, sensitivity = fit_age_bayes(response, age, cov)
            self.assertNotEqual(invalid["status"], "ok")
            self.assertEqual(sensitivity, [])


class CachedAgePipelineTests(unittest.TestCase):
    def test_432_and_232_use_only_authorized_families_and_preserve_source_inputs(self):
        for cortical, expected in ((400, {"parcel": 432, "yeo": 7, "system": 4, "amygdala": 6}),
                                  (200, {"parcel": 232})):
            with self.subTest(cortical=cortical), tempfile.TemporaryDirectory() as temporary:
                root = Path(temporary)
                input_dir = root / "source/inputs"
                _write_cache(input_dir, n_subjects=48)
                if cortical == 200:
                    parcels = pd.read_csv(input_dir / "parcel_response.csv", dtype={"sub_id": str})
                    parcels[["sub_id", *map(str, range(1, 233))]].to_csv(input_dir / "parcel_response.csv", index=False)
                    note = json.loads((input_dir / "tian_combined.json").read_text())
                    note.update(n_cortical=200, n_parcels=232)
                    (input_dir / "tian_combined.json").write_text(json.dumps(note))
                before = {path.name: path.read_bytes() for path in input_dir.iterdir()}
                with contextlib.redirect_stdout(io.StringIO()):
                    statistics = runner.run_analysis(input_dir, root / "bayes", skip_docx=True)
                self.assertEqual(statistics.groupby("dimension").size().to_dict(), expected)
                self.assertEqual(before, {path.name: path.read_bytes() for path in input_dir.iterdir()})
                self.assertEqual(len(pd.read_csv(root / "bayes/age_bayes_comparisons.csv")), sum(expected.values()) * 3)
                self.assertEqual(len(pd.read_csv(root / "bayes/age_bayes_prior_sensitivity.csv")), 9)
                design = pd.read_csv(input_dir / "design.csv")
                design.loc[0, "meanFD"] = .6
                design.to_csv(input_dir / "design.csv", index=False)
                with self.assertRaisesRegex(ValueError, "FD <= 0.5"):
                    runner.run_analysis(input_dir, root / "unauthorized", skip_docx=True)


if __name__ == "__main__":
    unittest.main()
