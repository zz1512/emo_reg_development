"""Scientific invariants: shared targets, train-only processing and paired nulls."""
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import numpy as np
import pandas as pd
from sklearn.cross_decomposition import PLSRegression
from sklearn.exceptions import ConvergenceWarning
from threadpoolctl import threadpool_limits

from .models import fit_curve, prepare, PLS_MAX_ITER, PLS_TOL
from .pipeline import Data, nested_loso
from .resampling import reduced_null, within_site_permutation, bootstrap_index, merge
from .sample import precheck
from .run import mock


def fixture(seed=8, n=72):
    rng = np.random.default_rng(seed)
    x = rng.normal(size=(n, 12))
    x[:, 6:] += x[:, :6] * .3
    cov = np.column_stack([np.arange(n) % 2, rng.uniform(.03, .4, n)])
    y = np.column_stack([12 + 1.7*x[:, 7] + .4*x[:, 0] + rng.normal(size=n),
                         .5*x[:, 7] + .2*x[:, 0] + rng.normal(scale=.6, size=n)])
    y += cov @ np.array([[.1, .1], [1., .3]])
    return Data(x, y, cov, np.asarray([str(i % 4) for i in range(n)]),
                np.asarray([f"s{i:05d}" for i in range(n)]), np.arange(12) < 6)


class JointPLSTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.threads = threadpool_limits(limits=1)

    @classmethod
    def tearDownClass(cls):
        cls.threads.restore_original_limits()

    def test_prefix_curve_matches_separately_fitted_multivariate_pls(self):
        d = fixture()
        p = prepare(d.x[:50], d.y[:50], d.cov[:50])
        xt, _, _ = p.transform(d.x[50:], None, d.cov[50:])
        curve = fit_curve(p.x, p.y, 4, p.rank_c)
        for k in range(1, 5):
            reference = PLSRegression(n_components=k, scale=False).fit(p.x, p.y)
            np.testing.assert_allclose(curve.predict(xt, k), reference.predict(xt), atol=1e-10)

    def test_slow_pls_converges_with_extended_budget_without_relaxing_tolerance(self):
        rng = np.random.default_rng(421)
        x = rng.normal(size=(60, 2))
        x = np.linalg.qr(x - x.mean(0))[0]
        angle = 1.2
        rotation = np.array([[np.cos(angle), -np.sin(angle)],
                             [np.sin(angle), np.cos(angle)]])
        y = x @ np.diag([1., np.sqrt(.997)]) @ rotation
        # A genuine numerical regression: the old 500-step fit emits a warning.
        with self.assertWarns(ConvergenceWarning):
            PLSRegression(n_components=1, scale=False, max_iter=500, tol=PLS_TOL).fit(x, y)
        curve = fit_curve(x, y, 1)
        self.assertGreater(curve.n_iter[0], 500)
        self.assertLess(curve.n_iter[0], PLS_MAX_ITER)
        reference = PLSRegression(n_components=1, scale=False,
                                  max_iter=PLS_MAX_ITER, tol=PLS_TOL).fit(x, y)
        np.testing.assert_allclose(curve.predict(x, 1), reference.predict(x), atol=1e-12)
        # Hitting even the enlarged budget must still fail, never silently pass.
        with patch("analysis.fig03.joint_pls.models.PLS_MAX_ITER", 1):
            with self.assertRaisesRegex(ConvergenceWarning, "max_iter=1"):
                fit_curve(x, y, 1)

    def test_heldout_y_cannot_change_axis_tuning_or_predictions_for_that_site(self):
        d = fixture()
        original = nested_loso(d, 3)
        y = d.y.copy()
        y[d.sites == "0"] += np.array([1000., -2000.])
        changed = nested_loso(d.with_y(y), 3)
        a = original["predictions"].query("site_id == '0'").reset_index(drop=True)
        b = changed["predictions"].query("site_id == '0'").reset_index(drop=True)
        np.testing.assert_allclose(a[["age_pred", "success_pred", "shared_pred", "k"]],
                                   b[["age_pred", "success_pred", "shared_pred", "k"]], atol=1e-10)
        np.testing.assert_allclose(original["axes"].query("outer_site == '0'")[["age_weight", "success_weight"]],
                                   changed["axes"].query("outer_site == '0'")[["age_weight", "success_weight"]])

    def test_same_common_target_in_every_model_and_both_routes(self):
        result = nested_loso(fixture(), 3)
        rows = result["predictions"]
        self.assertEqual(len(rows), 72*4*2)
        self.assertTrue(rows.groupby("sub_id").shared_true.nunique().eq(1).all())
        self.assertTrue(rows.groupby("sub_id").age_z.nunique().eq(1).all())
        self.assertEqual(len(result["performance"]), 8)
        self.assertEqual(len(result["comparisons"]), 10)

    def test_change_age_units_does_not_dominate_joint_training(self):
        d = fixture()
        a = nested_loso(d, 2, details=False)
        y = d.y.copy()
        y[:, 0] *= 1000.
        b = nested_loso(d.with_y(y), 2, details=False)
        np.testing.assert_allclose(a["performance"][["r2", "mae", "rmse"]],
                                   b["performance"][["r2", "mae", "rmse"]], atol=1e-9)

    def test_paired_permutation_and_stratified_bootstrap(self):
        d = fixture()
        index = within_site_permutation(d.sites, 88, 2)
        np.testing.assert_array_equal(np.sort(index), np.arange(72))
        np.testing.assert_array_equal(d.sites[index], d.sites)
        np.testing.assert_array_equal(index, within_site_permutation(d.sites, 88, 2))
        fit, residual, meta = reduced_null(d, "pfc")
        np.testing.assert_allclose(fit+residual, d.y)
        np.testing.assert_allclose(np.sort(residual[index, 0]), np.sort(residual[:, 0]))
        # Both targets must move together, not independently shuffled columns.
        null_y = fit + residual[index]
        np.testing.assert_allclose(null_y-fit, residual[index])
        self.assertGreater(meta["reduced_residual_df"], 5)
        for site in np.unique(d.sites):
            np.testing.assert_allclose(residual[d.sites == site].mean(0), 0, atol=1e-12)
        draw = bootstrap_index(d.sites, 88, 2)
        self.assertEqual(len(draw), len(d.y))
        for site in np.unique(d.sites):
            self.assertEqual(sum(d.sites[draw] == site), sum(d.sites == site))
        self.assertLess(len(np.unique(draw)), len(draw))

    def test_failed_or_partial_permutations_do_not_produce_p_values(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            observed = nested_loso(fixture(), 1, details=False)["comparisons"]
            observed.to_csv(root / "observed_comparisons.csv", index=False)
            dest = root / "permutations" / "pfc"
            dest.mkdir(parents=True)
            (dest / "manifest.json").write_text(json.dumps(dict(total=2, null="pfc")))
            table = observed.copy()
            table["replicate"], table["status"] = 0, "ok"
            table.to_csv(dest / "rep_000000.csv", index=False)
            partial = merge(root)
            self.assertTrue(partial.p_permutation.isna().all())
            pd.DataFrame([dict(replicate=1, status="failed", reason="deliberate")]).to_csv(dest / "rep_000001.csv", index=False)
            failed = merge(root)
            self.assertTrue(failed.p_permutation.isna().all())
            table["replicate"] = 1
            table.to_csv(dest / "rep_000001.csv", index=False)
            complete = merge(root)
            primary = complete.query("comparison == 'distributed_minus_pfc'")
            np.testing.assert_allclose(primary.p_permutation, 1.)
            np.testing.assert_allclose(primary.q_primary, 1.)
            self.assertTrue(complete.query("comparison != 'distributed_minus_pfc'").p_permutation.isna().all())

    def test_complete_mock_recovers_education_exclusions_and_real_row_ids(self):
        membership = Path(__file__).resolve().parents[3] / "outputs/fig03/predictive_pls_432_FD05/pfc_definition/pfc_parcel_membership.csv"
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder) / "mock"
            analysis = mock(root, membership, n=96, permutations=0, bootstraps=0, docx=False)
            meta = json.loads((analysis / "precheck/sample_metadata.json").read_text(encoding="utf-8"))
            self.assertEqual(meta["n"], 96)
            self.assertEqual(meta["restored_education_excluded_n"], 2)
            self.assertFalse(meta["parent_education_used"])
            sample = pd.read_csv(analysis / "precheck/common_sample.csv")
            self.assertTrue(sample.sub_id.str.startswith("s").all())
            with self.assertRaisesRegex(ValueError, "lack cached means"):
                precheck(root / "source", root / "source/behavior.csv", membership, root / "missing_recovery",
                         id_source=root / "source/row_ids.csv")
            # Reconstruct the source cohort rather than inheriting education-filtered IDs.
            source = root / "source"
            people = pd.read_csv(source / "design.csv", dtype={"sub_id": str})
            people[["sub_id", "age", "gender", "site_id"]].to_csv(source / "demographics.csv", index=False)
            (source / "valid.txt").write_text("\n".join(people.sub_id.iloc[:-1]), encoding="utf-8")
            (source / "excluded.txt").write_text(people.sub_id.iloc[-2], encoding="utf-8")
            pd.DataFrame(dict(sub_id=people.sub_id, mean_fd=people.meanFD,
                              run_name="EMO")).to_csv(source / "motion.csv", index=False)
            rebuilt = root / "rebuilt"
            meta = precheck(source, source / "behavior.csv", membership, rebuilt,
                            id_source=source / "row_ids.csv", volume_dir=source / "volumes",
                            demographics=source / "demographics.csv", subjects=source / "valid.txt",
                            excluded=source / "excluded.txt", meanfd=source / "motion.csv")
            self.assertEqual(meta["n"], 94)
            audit = pd.read_csv(rebuilt / "upstream_cohort_audit.csv")
            self.assertEqual(audit.upstream_exclusion.eq("not_in_valid_subject_list").sum(), 1)
            self.assertEqual(audit.upstream_exclusion.eq("second_level_removed").sum(), 1)
            roster = pd.read_csv(rebuilt / "sample_roster.csv")
            self.assertEqual(roster.first_exclusion.eq("independent_original_qc").sum(), 1)
            self.assertTrue(set(people.sub_id.iloc[:2]).issubset(set(pd.read_csv(rebuilt / "common_sample.csv").sub_id)))

    def test_meta_pfc_membership_is_frozen_and_identified_in_report(self):
        membership = Path(__file__).resolve().parents[3] / "outputs/fig03/predictive_pls_432_FD05_meta_pfc/pfc_definition/pfc_parcel_membership.csv"
        with tempfile.TemporaryDirectory() as folder:
            analysis = mock(Path(folder) / "meta", membership, n=96, permutations=0, bootstraps=0, docx=False)
            meta = json.loads((analysis / "precheck/sample_metadata.json").read_text(encoding="utf-8"))
            self.assertEqual(meta["pfc_definition"]["source"], "Buhle2014_reappraisal_PFC")
            self.assertEqual(meta["pfc_definition"]["n_pfc"], 10)
            self.assertEqual(meta["pfc_definition"]["n_remaining"], 422)
            frozen = pd.read_csv(analysis / "precheck/pfc_parcel_membership.csv")
            with np.load(analysis / "precheck/data.npz", allow_pickle=False) as archive:
                np.testing.assert_array_equal(frozen.is_pfc, archive["pfc"])
            performance = pd.read_csv(analysis / "model_performance.csv")
            self.assertTrue(performance.loc[performance.model.eq("M1_PFC"), "n_features"].eq(10).all())
            self.assertTrue(performance.loc[performance.model.eq("M3_nonPFC"), "n_features"].eq(422).all())
            text = (analysis / "report/result_summary.md").read_text(encoding="utf-8")
            self.assertIn("Buhle2014_reappraisal_PFC", text)
            self.assertIn("Remaining regions", text)


if __name__ == "__main__":
    unittest.main()
