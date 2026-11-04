"""Regression checks for the three-model and multi-dimension pyGAM contract.

Run from the repository root with:
    python -m unittest analysis.fig03.gam_pygam.test_models

The saved-input checks use small synthetic CSVs. The statistical checks use
the real repository pyGAM helper when pygam is installed.
"""

from __future__ import annotations

import contextlib
import importlib.util
import io
import json
from pathlib import Path
import sys
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch

import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix
from scipy.stats import f as f_distribution, false_discovery_control, t

from analysis.fig03.gam_nested.amygdala import AMYGDALA_ORDER
from analysis.fig03.gam_pygam import run_pygam as runner


class _ScipyBH:
    """Independent BH oracle for tests that replace the expensive unit fit."""

    def __init__(self) -> None:
        self.calls: list[np.ndarray] = []

    def fdr_bh(self, p_values, alpha=0.05):
        values = np.asarray(p_values, dtype=float)
        self.calls.append(values.copy())
        q_values = false_discovery_control(values, method="bh")
        return {"qvals": q_values, "reject": q_values < alpha}


def _synthetic_fit(p_values: dict[int, tuple[float, float]] | None = None):
    """Give each response a distinguishable fitted result and residual."""
    def fit(stat_utils, y, age, covariates, grid, lam_grid, store_partial=False):
        tag = int(y[0])
        first, second = (p_values or {}).get(tag, (0.4, 0.8))
        statistics = {
            "n_obs": len(y), "p_m3_vs_m1": first, "p_m3_vs_m2": second,
            "partial_r2_m3_vs_m1": 0.1, "partial_r2_m3_vs_m2": 0.02,
        }
        effect = np.full(len(grid), tag, dtype=float)
        partial = None
        if store_partial:
            partial = {
                "effect": effect, "ci_low": effect - 0.1, "ci_high": effect + 0.1,
                "residual": np.full(len(y), tag, dtype=float),
            }
        return statistics, effect, partial
    return fit


def _write_inputs(input_dir: Path, n_parcels: int = 2,
                  n_cortical: int | None = None,
                  group_order: tuple[tuple[str, tuple[str, ...]], ...] | None = None,
                  amygdala_order: tuple[str, ...] = ()) -> None:
    input_dir.mkdir(parents=True, exist_ok=True)
    subject_ids = [f"sub-{index:03d}" for index in range(16)]
    design = pd.DataFrame({
        "sub_id": subject_ids, "age": np.linspace(6, 18, 16),
        "gender": np.arange(16) % 2, "site_id": np.arange(16) % 4 + 1,
        "parent_edu": np.arange(16) % 3 + 1, "meanFD": np.linspace(0.05, 0.2, 16),
        "in_analysis": [1] * 15 + [0], "excluded_by": [""] * 15 + ["synthetic_exclusion"],
    })
    design.to_csv(input_dir / "design.csv", index=False)
    values = {"sub_id": subject_ids}
    values.update({str(parcel): np.full(16, parcel) for parcel in range(1, n_parcels + 1)})
    pd.DataFrame(values).to_csv(input_dir / "parcel_response.csv", index=False)
    pd.DataFrame({"parcel": range(1, n_parcels + 1),
                  "name": [f"parcel_{parcel}" for parcel in range(1, n_parcels + 1)]}).to_csv(
        input_dir / "parcel_info.csv", index=False)
    if group_order is None:
        group_order = (("yeo", runner.YEO_ORDER), ("system", runner.SYSTEM_ORDER))
    rows = []
    tag = n_parcels
    for dimension, names in group_order:
        for name in names:
            tag += 1
            rows.extend({"sub_id": subject, "dimension": dimension,
                         "group": name, "y": tag} for subject in subject_ids)
    pd.DataFrame(rows, columns=["sub_id", "dimension", "group", "y"]).to_csv(
        input_dir / "group_response.csv", index=False)
    if amygdala_order:
        amygdala = {"sub_id": subject_ids}
        for name in amygdala_order:
            tag += 1
            amygdala[name] = np.full(16, tag)
        pd.DataFrame(amygdala).to_csv(input_dir / "amygdala_response.csv", index=False)
        pd.DataFrame({"roi": amygdala_order}).to_csv(input_dir / "amygdala_order.csv", index=False)
    if n_cortical is not None:
        # validate_saved_inputs checks presence; no spatial operation runs here.
        (input_dir / "combined_atlas.nii.gz").write_bytes(b"synthetic atlas sentinel")
        (input_dir / "tian_combined.json").write_text(json.dumps({
            "n_cortical": n_cortical, "n_parcels": n_cortical + 32,
            "combined_atlas": "combined_atlas.nii.gz",
            "amygdala_parcels": [n_cortical + label for label in (3, 4, 19, 20)],
            "marginal_parcels": [156, 157] if n_cortical == 400 else [],
        }), encoding="utf-8")


@unittest.skipUnless(importlib.util.find_spec("pygam") is not None, "pygam is not installed")
class ThreeModelStatisticsTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.stats = runner.load_historical_stats()
        rng = np.random.default_rng(8403)
        cls.age = np.linspace(6, 18, 180)
        cls.covariates = rng.normal(size=(180, 4))
        cls.covariate_effect = cls.covariates @ np.array([0.3, -0.2, 0.1, 0.4])
        cls.noise = rng.normal(0, 0.08, 180)
        cls.grid = np.linspace(6.2, 17.8, 41)
        cls.lam_grid = 10.0 ** np.linspace(*runner.LAM_LOG_RANGE)

    def fit(self, y, age=None, covariates=None):
        return runner._fit_one(
            self.stats, y, self.age if age is None else age,
            self.covariates if covariates is None else covariates,
            self.grid, self.lam_grid, store_partial=True,
        )

    def test_linear_signal_has_overall_age_without_evidence_of_bending(self):
        y = 0.5 * (self.age - 12) + self.covariate_effect + self.noise
        result, centered, partial = self.fit(y)
        self.assertLess(result["p_m3_vs_m1"], 1e-6)
        self.assertGreater(result["p_m3_vs_m2"], 0.05)
        self.assertLess(result["rss_m2"], result["rss_m1"])
        self.assertEqual(result["edf_m1"], 5)
        self.assertEqual(result["edf_m2"], 6)
        self.assertAlmostEqual(result["p_age"], result["p_m3_vs_m1"])
        self.assertAlmostEqual(result["partial_r2"], result["partial_r2_m3_vs_m1"])
        self.assertEqual(partial["effect"].shape, self.grid.shape)
        np.testing.assert_array_equal(centered, partial["effect"])
        self.assertTrue(np.all(partial["ci_high"] >= partial["ci_low"]))

    def test_curved_signal_adds_explanation_beyond_linear_age(self):
        y = 0.6 * ((self.age - 12) / 3) ** 2 + self.covariate_effect + self.noise
        result, _, _ = self.fit(y)
        self.assertLess(result["p_m3_vs_m1"], 1e-6)
        self.assertLess(result["p_m3_vs_m2"], 1e-6)
        self.assertLess(result["rss_m3"], result["rss_m2"] * 0.1)
        self.assertGreater(result["edf_m3"], result["edf_m2"])
        self.assertAlmostEqual(result["partial_r2_m3_vs_m2"],
                               1 - result["rss_m3"] / result["rss_m2"])

    def test_nonlinearity_matches_fixed_cubic_spline_space_across_penalties(self):
        # Five cubic B-splines span cubics plus a truncated cubic at the sole
        # internal knot. This independent parameterization is the F oracle.
        x = (self.age - 12) / 6
        y = 0.2 * x ** 2 + self.covariate_effect + self.noise
        reduced = np.column_stack([np.ones(len(x)), x, self.covariates])
        full = np.column_stack([reduced, x ** 2, x ** 3, np.maximum(x, 0) ** 3])
        residual_reduced = y - reduced @ np.linalg.lstsq(reduced, y, rcond=None)[0]
        residual_full = y - full @ np.linalg.lstsq(full, y, rcond=None)[0]
        rss_reduced, rss_full = residual_reduced @ residual_reduced, residual_full @ residual_full
        expected_f = ((rss_reduced - rss_full) / 3) / (rss_full / (len(x) - 9))
        expected_p = f_distribution.sf(expected_f, 3, len(x) - 9)
        for penalties in (np.array([1e-5, 1e-4]), np.array([1e4, 1e5])):
            with self.subTest(penalties=penalties):
                result, _, _ = runner._fit_one(
                    self.stats, y, self.age, self.covariates, self.grid, penalties)
                self.assertEqual(result["rank_m3_basis"], 9)
                self.assertEqual(result["df1_m3_vs_m2"], 3)
                self.assertEqual(result["df2_m3_vs_m2"], len(x) - 9)
                self.assertAlmostEqual(result["rss_m3_basis"], rss_full, places=10)
                self.assertAlmostEqual(result["F_m3_vs_m2"], expected_f, places=8)
                self.assertAlmostEqual(result["p_m3_vs_m2"], expected_p, places=12)

    def test_all_models_use_the_same_complete_observations(self):
        y = 0.4 * (self.age - 12) + self.covariate_effect + self.noise
        y = y.copy()
        age = self.age.copy()
        covariates = self.covariates.copy()
        y[2], age[5], covariates[9, 1] = np.nan, np.nan, np.nan
        result, _, partial = self.fit(y, age, covariates)
        keep = np.ones(len(y), dtype=bool)
        keep[[2, 5, 9]] = False
        expected_design = np.column_stack([np.ones(keep.sum()), age[keep], covariates[keep]])
        expected_beta = np.linalg.lstsq(expected_design, y[keep], rcond=None)[0]
        expected_residual = y[keep] - expected_design @ expected_beta
        self.assertEqual(result["n_obs"], int(keep.sum()))
        self.assertAlmostEqual(result["rss_m2"], float(expected_residual @ expected_residual))
        self.assertAlmostEqual(result["beta_age_m2"], expected_beta[1])
        self.assertTrue(np.isnan(partial["residual"][[2, 5, 9]]).all())
        self.assertTrue(np.isfinite(partial["residual"][keep]).all())

    def test_partial_residuals_use_m3_when_age_correlates_with_covariates(self):
        covariates = self.covariates.copy()
        covariates[:, 0] += (self.age - self.age.mean()) / 2
        y = 0.3 * (self.age - 12) + covariates @ np.array([0.7, -.2, .1, .4]) + self.noise
        _, _, _, _, model = self.stats.fit_gam_full_and_reduced(
            y, self.age[:, None], covariates, k_splines=runner.K_SPLINES,
            lam_grid=self.lam_grid, verbose=False)
        observed = np.column_stack([self.age, covariates])
        mean_covariates = covariates.mean(axis=0)
        grid_design = np.column_stack([self.grid, np.tile(mean_covariates, (len(self.grid), 1))])
        original_coefficients = model.coef_.copy()
        original_covariance = model.statistics_["cov"].copy()
        effect, _, _, points = runner.age_effect_on_residual_scale(
            model, y, self.age, covariates, np.ones(len(y), dtype=bool), self.grid, mean_covariates)
        np.testing.assert_array_equal(model.coef_, original_coefficients)
        np.testing.assert_array_equal(model.statistics_["cov"], original_covariance)
        # pyGAM's separate age-term interface is an independent oracle for
        # the full-model contrast used by the plotting path.
        observed_term = model.partial_dependence(term=0, X=observed)
        expected_age = observed_term - observed_term.mean()
        full_residual = y - model.predict(observed)
        np.testing.assert_allclose(points, full_residual + expected_age, atol=1e-10)
        np.testing.assert_allclose(effect, model.partial_dependence(term=0, X=grid_design)
                                   - observed_term.mean(), atol=1e-10)
        self.assertAlmostEqual(float((points - full_residual).mean()), 0, places=10)
        no_age = np.column_stack([np.ones(len(y)), covariates])
        old_residual = y - no_age @ np.linalg.lstsq(no_age, y, rcond=None)[0]
        self.assertGreater(float(np.std(old_residual - points)), 0.1)


class CenteredAgeIntervalTests(unittest.TestCase):
    def test_linear_effect_ci_removes_shared_baseline_uncertainty(self):
        # In this simple model the exact contrast is beta_age * (age-mean_age).
        # A large intercept/covariate variance must cancel, rather than being
        # retained by shifting an ordinary full-prediction interval.
        age = np.array([6., 8., 10., 12., 14., 16., 18.])
        covariates = np.column_stack([np.linspace(-2, 2, len(age))])
        coefficients = np.array([10., .1, .5])
        model = SimpleNamespace(
            coef_=coefficients,
            statistics_={"edof": 3., "cov": np.diag([100., .01, 50.])},
            distribution=SimpleNamespace(_known_scale=False),
            _modelmat=lambda design: csr_matrix(np.column_stack([np.ones(len(design)), design])),
            predict=lambda design: np.column_stack([np.ones(len(design)), design]) @ coefficients,
        )
        y = model.predict(np.column_stack([age, covariates])) + np.array([.1, -.1, 0, 0, 0, -.1, .1])
        grid = np.array([6., 12., 18.])
        effect, low, high, points = runner.age_effect_on_residual_scale(
            model, y, age, covariates, np.ones(len(y), dtype=bool), grid, covariates.mean(axis=0))
        expected = .1 * (grid - age.mean())
        half_width = t.ppf(.975, len(y) - 3) * .1 * np.abs(grid - age.mean())
        np.testing.assert_allclose(effect, expected, atol=1e-12)
        np.testing.assert_allclose(low, expected - half_width, atol=1e-12)
        np.testing.assert_allclose(high, expected + half_width, atol=1e-12)
        self.assertAlmostEqual(low[1], high[1], places=12)
        # Changing the intercept shifts y but leaves centered curves, points,
        # and their interval untouched.
        model.coef_[0] += 100
        shifted = runner.age_effect_on_residual_scale(
            model, y + 100, age, covariates, np.ones(len(y), dtype=bool), grid, covariates.mean(axis=0))
        for original, updated in zip((effect, low, high, points), shifted):
            np.testing.assert_allclose(original, updated, atol=1e-12)


class SavedInputContractTests(unittest.TestCase):
    def test_each_dimension_and_comparison_has_its_own_bh_family(self):
        p_values = {
            1: (0.01, 0.9), 2: (0.04, 0.9),
            3: (0.02, 0.002), 4: (0.8, 0.01),
            5: (0.06, 0.02), 6: (0.07, 0.9),
            7: (0.003, 0.005), 8: (0.6, 0.8),
        }
        oracle = _ScipyBH()
        with tempfile.TemporaryDirectory(prefix="pygam-model-check-") as temporary:
            input_dir = Path(temporary)
            order = (("yeo", ("shared", "other")), ("system", ("shared", "other")))
            _write_inputs(input_dir, group_order=order, amygdala_order=("shared", "other"))
            with patch.object(runner, "YEO_ORDER", ("shared", "other")), \
                    patch.object(runner, "SYSTEM_ORDER", ("shared", "other")), \
                    patch.object(runner, "_fit_one", side_effect=_synthetic_fit(p_values)), \
                    contextlib.redirect_stdout(io.StringIO()):
                statistics, curves, residuals = runner.fit_saved_means(oracle, input_dir, store_partial=True)
        self.assertEqual(len(oracle.calls), 8)
        expected_tags = {"parcel": (1, 2), "yeo": (3, 4), "system": (5, 6), "amygdala": (7, 8)}
        for dimension, tags in expected_tags.items():
            rows = statistics[statistics["dimension"] == dimension]
            for index, comparison in enumerate(runner.COMPARISONS):
                expected = false_discovery_control([p_values[tag][index] for tag in tags], method="bh")
                np.testing.assert_allclose(rows[f"q_{comparison}"], expected)
        np.testing.assert_allclose(statistics["q_age"], statistics["q_m3_vs_m1"])
        self.assertEqual(set(curves["dimension"]), set(expected_tags))
        self.assertFalse(residuals.duplicated(["sub_id", "dimension", "name"]).any())
        self.assertEqual(residuals["sub_id"].nunique(), 15)
        for dimension, tag in (("yeo", 3), ("system", 5), ("amygdala", 7)):
            selected = residuals[(residuals["dimension"] == dimension) & (residuals["name"] == "shared")]
            self.assertEqual(len(selected), 15)
            np.testing.assert_array_equal(selected["residual"], np.full(15, tag))

    def test_tian_232_ignores_cached_network_system_and_juelich_responses(self):
        with tempfile.TemporaryDirectory(prefix="pygam-232-check-") as temporary:
            input_dir = Path(temporary)
            _write_inputs(input_dir, n_parcels=232, n_cortical=200, amygdala_order=AMYGDALA_ORDER)
            with patch.object(runner, "_fit_one", side_effect=_synthetic_fit()) as fit, \
                    contextlib.redirect_stdout(io.StringIO()):
                statistics, curves, residuals = runner.fit_saved_means(_ScipyBH(), input_dir, store_partial=True)
            self.assertEqual(fit.call_count, 232)
            for frame in (statistics, curves, residuals):
                self.assertEqual(set(frame["dimension"]), {"parcel"})
            counts = runner.saved_counts(input_dir)
            self.assertTrue(counts["parcels_only"])
            self.assertEqual(counts["n_amygdala"], 0)

    def test_tian_432_requires_six_juelich_responses_before_fitting(self):
        with tempfile.TemporaryDirectory(prefix="pygam-432-check-") as temporary:
            input_dir = Path(temporary)
            _write_inputs(input_dir, n_parcels=432, n_cortical=400)
            with patch.object(runner, "_fit_one") as fit:
                with self.assertRaisesRegex(ValueError, "Juelich six-ROI"):
                    runner.fit_saved_means(_ScipyBH(), input_dir)
                fit.assert_not_called()
            pd.DataFrame({"sub_id": ["sub-000"], **{roi: [0] for roi in AMYGDALA_ORDER[:-1]}}).to_csv(
                input_dir / "amygdala_response.csv", index=False)
            with self.assertRaisesRegex(ValueError, "Juelich six-ROI"):
                runner.validate_saved_inputs(input_dir)

    def test_tian_432_complete_cache_fits_all_four_scales(self):
        with tempfile.TemporaryDirectory(prefix="pygam-432-complete-") as temporary:
            input_dir = Path(temporary)
            _write_inputs(input_dir, n_parcels=432, n_cortical=400, amygdala_order=AMYGDALA_ORDER)
            oracle = _ScipyBH()
            with patch.object(runner, "_fit_one", side_effect=_synthetic_fit()) as fit, \
                    contextlib.redirect_stdout(io.StringIO()):
                statistics, curves, residuals = runner.fit_saved_means(oracle, input_dir, store_partial=True)
            self.assertEqual(fit.call_count, 449)
            self.assertEqual(statistics.groupby("dimension").size().to_dict(),
                             {"parcel": 432, "yeo": 7, "system": 4, "amygdala": 6})
            self.assertEqual(sorted(map(len, oracle.calls)), [4, 4, 6, 6, 7, 7, 432, 432])
            self.assertEqual(set(curves["dimension"]), {"parcel", "yeo", "system", "amygdala"})
            self.assertEqual(residuals.groupby(["dimension", "name"]).size().nunique(), 1)
            self.assertFalse(residuals.duplicated(["sub_id", "dimension", "name"]).any())

    def test_tian_432_rejects_incomplete_network_family(self):
        with tempfile.TemporaryDirectory(prefix="pygam-432-groups-") as temporary:
            input_dir = Path(temporary)
            _write_inputs(input_dir, n_parcels=432, n_cortical=400, amygdala_order=AMYGDALA_ORDER)
            path = input_dir / "group_response.csv"
            groups = pd.read_csv(path)
            groups = groups[~((groups["dimension"] == "yeo") & (groups["group"] == runner.YEO_ORDER[0]))]
            groups.to_csv(path, index=False)
            with self.assertRaisesRegex(ValueError, "complete yeo responses"):
                runner.validate_saved_inputs(input_dir)


class PipelineOutputTests(unittest.TestCase):
    def test_reuse_csv_figures_maps_and_redraw_for_both_versions(self):
        import nibabel as nib
        from nibabel.cifti2.cifti2_axes import BrainModelAxis, LabelAxis
        from analysis.fig03.gam_pygam import brain_maps

        def project(image, method):
            # The spatial transform is deliberately replaced: this test checks
            # export wiring and file formats without downloading registration.
            self.assertIn(method, ("linear", "nearest"))
            data = np.full(8, 1.0 if method == "nearest" else 0.123, dtype=np.float32)
            return tuple(nib.GiftiImage(darrays=[nib.gifti.GiftiDataArray(data)]) for _ in range(2))

        with tempfile.TemporaryDirectory(prefix="pygam-pipeline-") as temporary:
            root = Path(temporary)
            brain = (BrainModelAxis.from_surface(np.arange(8), 8, "CortexLeft")
                     + BrainModelAxis.from_surface(np.arange(8), 8, "CortexRight"))
            template = root / "synthetic.dlabel.nii"
            nib.save(nib.Cifti2Image(np.ones((1, 16), dtype=np.int32),
                                   (LabelAxis(["synthetic"], [{1: ("region", (1, 0, 0, 1))}]), brain)),
                     str(template))
            for n_cortical in (200, 400):
                with self.subTest(n_cortical=n_cortical):
                    n_parcels = n_cortical + 32
                    source = root / f"source_{n_parcels}"
                    _write_inputs(source, n_parcels=n_parcels, n_cortical=n_cortical,
                                  amygdala_order=AMYGDALA_ORDER)
                    volume = np.arange(1, n_parcels + 1, dtype=np.int16).reshape(n_parcels, 1, 1)
                    nib.save(nib.Nifti1Image(volume, np.eye(4)), str(source / "combined_atlas.nii.gz"))
                    before = {path.name: path.read_bytes() for path in source.iterdir()}
                    output = root / f"output_{n_parcels}"
                    arguments = ["run_pygam", "--reuse-inputs", str(source), "--meanfd-max", "0.5",
                                 "--schaefer-atlas", str(template), "--output", str(output)]
                    with patch.object(sys, "argv", arguments), \
                            patch.object(runner, "load_historical_stats", return_value=_ScipyBH()), \
                            patch.object(runner, "_fit_one", side_effect=_synthetic_fit()), \
                            patch.object(brain_maps, "_project_volume", side_effect=project), \
                            patch.object(brain_maps, "write_brain_borders") as border_export, \
                            contextlib.redirect_stdout(io.StringIO()):
                        runner.main()
                    border_export.assert_called_once()
                    self.assertEqual(border_export.call_args.args[2], [156, 157] if n_cortical == 400 else [])
                    destination = Path(str(output) + "_FD05")
                    statistics = pd.read_csv(destination / "statistics/pygam_age_comparison.csv")
                    units = pd.read_csv(destination / "statistics/trajectory_units.csv")
                    expected = {"parcel": n_parcels}
                    if n_cortical == 400:
                        expected.update(yeo=7, system=4, amygdala=6)
                    self.assertEqual(statistics.groupby("dimension").size().to_dict(), expected)
                    self.assertEqual(len(units), 4 if n_cortical == 200 else 23)
                    self.assertEqual(len(list((destination / "figures/age_curves").glob("*.svg"))), len(units))
                    self.assertEqual(len(list((destination / "figures/age_curves_zoom").glob("*_zoom.svg"))), len(units))
                    self.assertTrue((units["zoom_y_min"] < units["zoom_y_max"]).all())
                    for filename in units["zoom_figure_file"]:
                        self.assertTrue((destination / "figures" / filename).is_file())
                    self.assertEqual(len(list((destination / "maps").iterdir())), 17)
                    metadata = json.loads((destination / "run_metadata.json").read_text())
                    self.assertEqual(metadata["counts"]["n_subjects"], 15)
                    self.assertEqual(metadata["counts"]["meanfd_max"], 0.5)
                    self.assertEqual(metadata["age_effect_definition"], runner.AGE_EFFECT_DEFINITION)
                    self.assertEqual(before, {path.name: path.read_bytes() for path in source.iterdir()})
                    for comparison in runner.COMPARISONS:
                        mapped = nib.load(str(destination / "maps" / f"partial_r2_{comparison}_volume_MNI152.nii.gz"))
                        self.assertEqual(mapped.get_data_dtype(), np.dtype("float32"))
                        scalar = nib.load(str(destination / "maps" / f"partial_r2_{comparison}_fslr32k.dscalar.nii"))
                        self.assertEqual(scalar.shape, (1, 16))
                    # Redraw commands must work without fitting any model.
                    for step in ("curves", "maps", "borders"):
                        arguments = ["run_pygam", "--step", step, "--output", str(destination),
                                     "--schaefer-atlas", str(template)]
                        with patch.object(sys, "argv", arguments), \
                                patch.object(runner, "load_historical_stats", side_effect=AssertionError("unexpected refit")), \
                                patch.object(brain_maps, "_project_volume", side_effect=project), \
                                patch.object(brain_maps, "write_brain_borders") as border_export, \
                                contextlib.redirect_stdout(io.StringIO()):
                            runner.main()
                        self.assertEqual(border_export.call_count, int(step != "curves"))

    def test_reuse_rejects_missing_tian_note(self):
        with tempfile.TemporaryDirectory(prefix="pygam-reuse-invalid-") as temporary:
            source = Path(temporary) / "source"
            _write_inputs(source)
            arguments = ["run_pygam", "--reuse-inputs", str(source), "--output", str(source.parent / "out")]
            with patch.object(sys, "argv", arguments), \
                    patch.object(runner, "load_historical_stats", return_value=_ScipyBH()), \
                    self.assertRaisesRegex(ValueError, "requires Tian saved inputs"):
                runner.main()

    def test_old_curve_definition_requires_fit_before_redraw(self):
        with tempfile.TemporaryDirectory(prefix="pygam-old-curves-") as temporary:
            output = Path(temporary)
            _write_inputs(output / "inputs", n_parcels=232, n_cortical=200)
            (output / "statistics").mkdir()
            for name in ("pygam_age_comparison.csv", "pygam_curves.csv", "covariate_residuals.csv"):
                (output / "statistics" / name).write_text("old-data-sentinel", encoding="utf-8")
            (output / "run_metadata.json").write_text("{}", encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "Run --step fit"):
                runner.rewrite_saved_curves(output)
            self.assertFalse((output / "figures").exists())


if __name__ == "__main__":
    unittest.main()
