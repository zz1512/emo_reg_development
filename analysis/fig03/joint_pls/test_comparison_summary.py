"""Results-only recovery must retain pairing, nulls and incomplete-run gates."""
import json
from pathlib import Path
import tempfile
import unittest

import numpy as np
import pandas as pd

from .comparison_summary import paired_comparisons, summarize
from .io import MODELS, SCHEMES


def performance():
    rows = []
    for scheme, scores in zip(SCHEMES, [(-.125, .125, .375, .25), (-.1, .05, .4, .35)]):
        for model, r2, count in zip(MODELS, scores, [0, 10, 432, 422]):
            rows.append(dict(scheme=scheme, model=model, n=96, sites=4, n_features=count,
                             r2=r2, mae=1-r2/2, rmse=1-r2/3))
    return pd.DataFrame(rows)


def archived(perf):
    rows = []
    pairs = [("pfc_minus_covariates", MODELS[1], MODELS[0]),
             ("distributed_minus_covariates", MODELS[2], MODELS[0]),
             ("nonpfc_minus_covariates", MODELS[3], MODELS[0]),
             ("distributed_minus_pfc", MODELS[2], MODELS[1]),
             ("distributed_minus_nonpfc", MODELS[2], MODELS[3])]
    for scheme in SCHEMES:
        p = perf.loc[perf.scheme.eq(scheme)].set_index("model")
        for name, a, b in pairs:
            rows.append(dict(scheme=scheme, comparison=name,
                             **{f"delta_{m}": p.at[a, m]-p.at[b, m] for m in ("r2", "mae", "rmse")}))
    return pd.DataFrame(rows)


class ComparisonSummaryTests(unittest.TestCase):
    def fixture(self, root, repeats=3):
        (root / "analysis_manifest.json").write_text(json.dumps(dict(data_sha256="frozen", code_sha256={"pipeline": "old"}, k_max=2)))
        performance().to_csv(root / "model_performance.csv", index=False)
        for kind, folder, null in [("permutation", "permutations/pfc", "pfc"),
                                   ("bootstrap", "bootstrap/independent", "independent")]:
            dest = root / folder
            dest.mkdir(parents=True)
            meta = dict(kind=kind, null=null, total=repeats, data_sha256="frozen", code_sha256={"pipeline": "old"}, k_max=2)
            (dest / "manifest.json").write_text(json.dumps(meta))
            for i in range(repeats):
                p = performance()
                if kind == "permutation":
                    p.loc[p.model.eq(MODELS[2]), "r2"] = [0., .2, .5][i]
                    p.loc[p.model.eq(MODELS[3]), "r2"] = [0., .4, .6][i]
                else:
                    p.loc[p.model.eq(MODELS[3]), "r2"] += (i-1)*.01
                table = archived(p)
                table["replicate"], table["status"] = i, "ok"
                # Shuffle row order so subtraction cannot accidentally rely on order.
                table.sample(frac=1, random_state=i).to_csv(dest / f"rep_{i:06d}.csv", index=False)

    def covariates_fixture(self, root):
        dest = root / "permutations/covariates"
        dest.mkdir(parents=True)
        meta = json.loads((root / "permutations/pfc/manifest.json").read_text())
        meta["null"] = "covariates"
        (dest / "manifest.json").write_text(json.dumps(meta))
        # Exceedance counts give raw P = .25, .5, .75, 1, .25, .5.
        counts = {"joint": (0, 1, 2), "shared": (3, 0, 1)}
        for i in range(3):
            perf = performance()
            for scheme in SCHEMES:
                for model, exceedances in zip(MODELS[1:], counts[scheme]):
                    selected = perf.scheme.eq(scheme) & perf.model.eq(model)
                    perf.loc[selected, "r2"] += 1 if i < exceedances else -1
            table = archived(perf)
            table["replicate"], table["status"] = i, "ok"
            table.to_csv(dest / f"rep_{i:06d}.csv", index=False)

    def test_paired_recovery_in_all_three_metrics_and_shuffled_rows(self):
        result = paired_comparisons(archived(performance()).sample(frac=1, random_state=1))
        self.assertEqual(len(result), 10)
        self.assertNotIn("distributed_minus_nonpfc", set(result.comparison))
        self.assertNotIn("nonpfc_minus_distributed", set(result.comparison))
        for scheme, delta in [("joint", .125), ("shared", .3)]:
            row = result.loc[result.scheme.eq(scheme) & result.comparison.eq("nonpfc_minus_pfc")].iloc[0]
            np.testing.assert_allclose([row.delta_r2, row.delta_mae, row.delta_rmse], [delta, -delta/2, -delta/3])

    def test_each_contrast_uses_its_own_null_and_secondary_bh(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            self.fixture(root)
            result, nulls = summarize(root, root / "report")
            self.assertEqual(len(nulls), 30)
            p = result.set_index(["scheme", "comparison"])
            self.assertEqual(p.at[("joint", "distributed_minus_pfc"), "p_permutation"], .5)
            self.assertEqual(p.at[("joint", "nonpfc_minus_pfc"), "p_permutation"], .75)
            secondary = result.query("family == 'secondary'")
            self.assertTrue(secondary.q_secondary.notna().all())
            self.assertTrue(secondary.q_primary.isna().all())
            self.assertTrue(result.query("family == 'exploratory'").p_permutation.isna().all())
            self.assertTrue(result.ci_low.notna().all())
            row = p.loc[("joint", "nonpfc_minus_pfc")]
            np.testing.assert_allclose([row.ci_low, row.ci_high], [.1155, .1345])
            self.assertEqual(len(pd.read_csv(root / "report/pfc_comparisons.csv")), 4)
            self.assertEqual(len(result), 10)
            self.assertNotIn("nonpfc_minus_distributed", set(result.comparison))
            self.assertNotIn("distributed_minus_nonpfc", set(result.comparison))

    def test_missing_or_failed_repeat_keeps_p_and_final_intervals_empty(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            self.fixture(root)
            (root / "permutations/pfc/rep_000002.csv").unlink()
            pd.DataFrame([dict(replicate=1, status="failed", reason="nonconvergence")]).to_csv(
                root / "bootstrap/independent/rep_000001.csv", index=False)
            result, _ = summarize(root, root / "report")
            self.assertTrue(result.p_permutation.isna().all())
            self.assertTrue(result.q_secondary.isna().all())
            self.assertTrue(result.ci_low.isna().all())
            self.assertTrue(result.query("family == 'secondary'").inference.eq("incomplete").all())

    def test_baseline_bh_pools_six_tests_and_preserves_existing_families(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            self.fixture(root)
            original, _ = summarize(root, root / "report")
            self.covariates_fixture(root)
            result, _ = summarize(root, root / "report")
            baseline = result.query("family == 'exploratory'")
            np.testing.assert_allclose(baseline.p_permutation, [.25, .5, .75, 1., .25, .5])
            np.testing.assert_allclose(baseline.q_exploratory, [.75, .75, .9, 1., .75, .75])
            np.testing.assert_allclose(baseline.q_bh, baseline.q_exploratory)
            self.assertTrue(baseline.bh_family_size.eq(6).all())
            for family, column in (("primary", "q_primary"), ("secondary", "q_secondary")):
                selected = result.family.eq(family)
                np.testing.assert_array_equal(result.loc[selected, column], original.loc[selected, column])
                np.testing.assert_array_equal(result.loc[selected, "q_bh"], result.loc[selected, column])
                self.assertTrue(result.loc[selected, "bh_family_size"].eq(2).all())
            self.assertFalse(result.family.eq("descriptive").any())
            self.assertEqual(len(result), 10)
            provenance = json.loads((root / "report/comparison_provenance.json").read_text())
            self.assertIn("6 exploratory", provenance["exploratory_family"])

    def test_incomplete_baseline_null_does_not_adjust_a_partial_family(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            self.fixture(root)
            self.covariates_fixture(root)
            (root / "permutations/covariates/rep_000002.csv").unlink()
            result, _ = summarize(root, root / "report")
            baseline = result.query("family == 'exploratory'")
            self.assertTrue(baseline.q_exploratory.isna().all())
            self.assertTrue(baseline.q_bh.isna().all())
            self.assertTrue(baseline.supported_exploratory.isna().all())
            self.assertTrue(result.query("family in ['primary', 'secondary']").q_bh.notna().all())

    def test_mixed_source_or_duplicate_records_are_rejected(self):
        table = archived(performance())
        with self.assertRaisesRegex(ValueError, "Duplicate"):
            paired_comparisons(pd.concat([table, table.iloc[:1]]))
        with self.assertRaisesRegex(ValueError, "Missing paired"):
            paired_comparisons(table.iloc[1:])
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            self.fixture(root)
            path = root / "permutations/pfc/manifest.json"
            meta = json.loads(path.read_text())
            meta["data_sha256"] = "different sample"
            path.write_text(json.dumps(meta))
            with self.assertRaisesRegex(ValueError, "data_sha256 differs"):
                summarize(root, root / "report")


if __name__ == "__main__":
    unittest.main()
