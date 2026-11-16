"""Synthetic saved-result fixtures and tests; never real study results."""
from __future__ import annotations

import contextlib
import io as stream
import json
from pathlib import Path
import tempfile
import unittest

import numpy as np
import pandas as pd

from analysis.fig03 import reporting as io
from analysis.fig03.gam_pygam import result_tables, summarize_bayes as bf
from analysis.fig03.predictive_pls import summarize_pls as pls


def mock_bayes(root, atlas=432):
    root.mkdir(parents=True)
    rows = []
    families = {"parcel": atlas, **({"yeo": 7, "system": 4, "amygdala": 6} if atlas == 432 else {})}
    for dimension, count in families.items():
        for unit in range(1, count + 1):
            # Three known categories and an exact log BF identity.
            logs = (4., 0., -4.)[(unit - 1) % 3], 0.
            row = dict(dimension=dimension, name=str(unit), label=f"SIMULATED {dimension} {unit}", status="ok", n_obs=991,
                       integration_converged=True, quadrature_relative_error=1e-10, quadrature_check_log_difference=1e-12)
            for comparison, log in zip(bf.COMPARISONS, (*logs, logs[0] - logs[1])):
                evidence = "strong_support_full" if log >= np.log(10) else "strong_support_reduced" if log <= -np.log(10) else "inconclusive"
                row.update({f"{field}_{comparison}": value for field, value in {
                    "log_bf10": log, "bf10": np.exp(log), "bf01": np.exp(-log), "evidence": evidence,
                    "sensitivity_log_bf10_min": log - .1, "sensitivity_log_bf10_max": log + .1,
                    "sensitivity_direction_stable": log != 0}.items()})
            rows.append(row)
    wide = pd.DataFrame(rows)
    wide.to_csv(root / "age_bayes_statistics.csv", index=False, float_format="%.17g")
    result_tables.age_long_table(wide).to_csv(root / "age_bayes_comparisons.csv", index=False, float_format="%.17g")
    (root / "run_metadata.json").write_text(json.dumps(dict(analysis_kind="age_bayesian_spline_v1", n_cortical=atlas - 32,
                                                         max_meanFD=.5, prior={"main_scales": [.5, .5]}, simulated=True)), encoding="utf-8")


def mock_pls(root, bootstrap=True, complete=True, n_valid=1950):
    root.mkdir(parents=True)
    frozen = root / "frozen"
    frozen.mkdir()
    membership = pd.DataFrame({"parcel": range(1, 433), "is_pfc": [i <= 100 for i in range(1, 433)],
                               "definition_source": "SIMULATED anatomical PFC", "rule": "SIMULATED rule"})
    membership.to_csv(frozen / "pfc_parcel_membership.csv", index=False)
    (frozen / "config_snapshot.json").write_text(json.dumps(dict(n_pfc=100, n_non_pfc=332, primary_comparison=pls.PRIMARY,
                                                                simulated=True)), encoding="utf-8")
    prediction_rows, p_rows, c_rows, boot_rows = [], [], [], []
    rng = np.random.default_rng(1120)
    for task in pls.TASKS:
        y = rng.normal(size=24)
        site = np.repeat(["A", "B", "C"], 8)
        row = pd.DataFrame(dict(task=task, sub_id=[f"SIMULATED{i}" for i in range(24)], site_id=site, outer_site=site, y_true=y))
        for model, weight in zip(pls.MODELS, (0., .1, .35, .2)):
            row[f"pred_{model}"] = weight * y + .8 * rng.normal(size=24)
        prediction_rows.append(row)
        metrics = {}
        for model in pls.MODELS:
            error = y - row[f"pred_{model}"]
            metrics[model] = dict(r2=1 - np.sum(error**2) / np.sum((y - y.mean())**2), mae=np.mean(abs(error)),
                                  rmse=np.sqrt(np.mean(error**2)), pearson_r=np.corrcoef(y, row[f"pred_{model}"])[0, 1])
        # This is a schema fixture, not a refitted bootstrap or a scientific estimate.
        reps = 2000 if complete else 1900
        draws = {model: metrics[model]["r2"] + rng.normal(0, .03, reps) for model in pls.MODELS}
        valid = np.arange(reps) < min(n_valid, reps)
        for i in range(reps):
            b = dict(task=task, replicate=i, valid=bool(valid[i]))
            b.update({f"r2_{model}": values[i] for model, values in draws.items()})
            for comparison in pls.COMPARISONS:
                left, right = comparison.split("-")
                b[f"delta_r2_{comparison}"] = draws[left][i] - draws[right][i]
            boot_rows.append(b)
        count = int(valid.sum()) if bootstrap else 0
        quality = bootstrap and count >= 1000 and count / 2000 >= .95
        for model, features in zip(pls.MODELS, (0, 100, 432, 332)):
            low, high = np.percentile(draws[model][valid], [2.5, 97.5]) if bootstrap and count else (np.nan, np.nan)
            p_rows.append(dict(task=task, model=model, n=24, n_sites=3, n_features=features, median_k=2,
                               **metrics[model], ci_low_r2=low, ci_high_r2=high, n_bootstrap_valid=count,
                               prediction_error_improvement_vs_m0=1 - (1 - metrics[model]["r2"]) / (1 - metrics[pls.MODELS[0]]["r2"])))
        for comparison in pls.COMPARISONS:
            left, right = comparison.split("-")
            low, high = np.percentile((draws[left] - draws[right])[valid], [2.5, 97.5]) if bootstrap and count else (np.nan, np.nan)
            applicable = task == "success" and comparison == pls.PRIMARY
            c_rows.append(dict(task=task, comparison=comparison, role="primary" if applicable else "support",
                               **{f"delta_{key}": metrics[left][key] - metrics[right][key] for key in ("r2", "mae", "rmse")},
                               ci_low_delta_r2=low, ci_high_delta_r2=high, n_bootstrap_valid=count,
                               inference_valid=bool(quality), supported=bool(quality and applicable and low > 0)))
    pd.concat(prediction_rows).to_csv(root / "oof_predictions.csv", index=False, float_format="%.17g")
    pd.DataFrame(p_rows).to_csv(root / "model_performance.csv", index=False, float_format="%.17g")
    pd.DataFrame(c_rows).to_csv(root / "model_comparisons.csv", index=False, float_format="%.17g")
    if bootstrap:
        pd.DataFrame(boot_rows).to_csv(root / "bootstrap_metrics.csv", index=False, float_format="%.17g")
        (frozen / "bootstrap_status.json").write_text(json.dumps(dict(requested_replicates=2000,
                                                                   valid_by_task={t: count for t in pls.TASKS})))


class ReviewerTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)

    def tearDown(self):
        self.temp.cleanup()

    def run_quiet(self, function, source, **kwargs):
        with contextlib.redirect_stdout(stream.StringIO()):
            return function(source, mock=True, docx=False, **kwargs)

    def test_bayes_432_counts_source_untouched(self):
        source = self.root / "bf"
        mock_bayes(source)
        before = {p: p.read_bytes() for p in source.iterdir()}
        result = self.run_quiet(bf.run_summary, source)
        self.assertEqual(result["row_counts"]["bf_main_summary.csv"], 8)
        self.assertEqual(result["row_counts"]["bf_primary_details.csv"], 898)
        self.assertEqual(result["row_counts"]["bf_linear_details.csv"], 449)
        frame = io.read_csv(source / "reviewer_summary/bf_main_summary.csv")
        row = frame.iloc[0]
        self.assertEqual((row.n_support_first, row.n_inconclusive, row.n_support_second), (144, 144, 144))
        for path, data in before.items():
            self.assertEqual(path.read_bytes(), data)

    def test_bayes_232_and_failed_unit(self):
        source = self.root / "bf"
        mock_bayes(source, 232)
        path = source / "age_bayes_statistics.csv"
        wide = io.read_csv(path)
        wide.loc[0, "status"] = "constant_age_or_response"
        wide.to_csv(path, index=False, float_format="%.17g")
        (source / "age_bayes_comparisons.csv").unlink()
        result = self.run_quiet(bf.run_summary, source)
        self.assertEqual(result["row_counts"]["bf_main_summary.csv"], 2)
        self.assertEqual(result["row_counts"]["bf_primary_details.csv"], 464)
        self.assertEqual(result["n_invalid_units"], 1)

    def test_bayes_overflow_and_identity(self):
        source = self.root / "bf"
        mock_bayes(source, 232)
        path = source / "age_bayes_statistics.csv"
        wide = io.read_csv(path)
        for comparison in ("m3_vs_m1", "m2_vs_m1"):
            wide.loc[0, [f"log_bf10_{comparison}", f"bf10_{comparison}", f"bf01_{comparison}",
                         f"sensitivity_log_bf10_min_{comparison}", f"sensitivity_log_bf10_max_{comparison}"]] = [1000., np.inf, 0., 999., 1001.]
        wide.to_csv(path, index=False, float_format="%.17g")
        (source / "age_bayes_comparisons.csv").unlink()
        self.run_quiet(bf.run_summary, source)
        wide.loc[0, "log_bf10_m2_vs_m1"] = 999.
        wide.to_csv(path, index=False, float_format="%.17g")
        with self.assertRaisesRegex(ValueError, "log BF31"):
            self.run_quiet(bf.run_summary, source)

    def test_bayes_integration_and_mock_guard(self):
        source = self.root / "bf"
        mock_bayes(source, 232)
        with self.assertRaisesRegex(ValueError, "Simulated input"):
            bf.load_results(source)
        path = source / "age_bayes_statistics.csv"
        wide = io.read_csv(path)
        wide.loc[0, "quadrature_relative_error"] = .01
        wide.to_csv(path, index=False)
        with self.assertRaisesRegex(ValueError, "integration failed"):
            bf.load_results(source, mock=True)

    def test_selected_bayes_six_individual_parcels(self):
        source = self.root / "bf"
        mock_bayes(source)
        parcels = [152, 154, 182, 185, 156, 157]
        result = self.run_quiet(lambda root, **kwargs: bf.run_selected_summary(root, parcels, **kwargs), source)
        self.assertEqual(result["row_counts"]["selected_region_summary.csv"], 6)
        self.assertEqual(result["row_counts"]["selected_region_comparisons.csv"], 18)
        output = source / "reviewer_selected_parcels"
        compact = io.read_csv(output / "selected_region_summary.csv")
        self.assertEqual(compact.parcel.tolist(), parcels)
        original = io.read_csv(source / "age_bayes_statistics.csv").set_index(["dimension", "name"])
        for row in compact.itertuples():
            self.assertEqual(row.bf10_age_spline_vs_no_age, original.at[("parcel", str(row.parcel)), "bf10_m3_vs_m1"])
        text = (output / "result_summary.md").read_text(encoding="utf-8")
        self.assertIn("6 parcels highlighted", text)
        self.assertNotIn("BF31", text)
        self.assertNotIn("BF32", text)
        self.assertNotIn("BF21", text)
        self.assertIn("rather than an independent whole-brain confirmation", text)

    def test_selected_bayes_rejects_other_atlas_and_duplicate_ids(self):
        source = self.root / "bf"
        mock_bayes(source, 232)
        with self.assertRaisesRegex(ValueError, "do not reuse"):
            bf.run_selected_summary(source, [152, 154], mock=True, docx=False)
        with self.assertRaisesRegex(ValueError, "duplicate IDs"):
            bf.run_selected_summary(source, [152, 152], mock=True, docx=False)

    def test_pls_complete_tables_and_originals(self):
        source = self.root / "pls"
        mock_pls(source)
        before = {p: p.read_bytes() for p in source.rglob("*") if p.is_file()}
        result = self.run_quiet(pls.run_summary, source)
        self.assertEqual(result["row_counts"], {"pls_main_summary.csv": 2, "pls_model_performance.csv": 8, "pls_all_comparisons.csv": 10})
        self.assertEqual(set(result["bootstrap_report_status"].values()), {"inference_available"})
        for path, data in before.items():
            self.assertEqual(path.read_bytes(), data)
        frame = io.read_csv(source / "reviewer_summary/pls_main_summary.csv")
        self.assertFalse(frame.loc[frame.task.eq("age"), "confirmation_applicable"].iloc[0])
        self.assertTrue(io.read_csv(source / "reviewer_summary/pls_model_performance.csv").r2.lt(0).any())

    def test_pls_missing_incomplete_and_poor_bootstrap(self):
        for name, kwargs, state in (("none", {"bootstrap": False}, "bootstrap_not_run"),
                                    ("partial", {"complete": False}, "bootstrap_incomplete"),
                                    ("poor", {"n_valid": 1800}, "bootstrap_quality_insufficient")):
            with self.subTest(name=name):
                source = self.root / name
                mock_pls(source, **kwargs)
                result = self.run_quiet(pls.run_summary, source)
                self.assertEqual(set(result["bootstrap_report_status"].values()), {state})
                frame = io.read_csv(source / "reviewer_summary/pls_all_comparisons.csv")
                self.assertFalse(frame.supported_for_report.any())

    def test_pls_stale_predictions_and_flags(self):
        source = self.root / "pls"
        mock_pls(source)
        path = source / "model_performance.csv"
        frame = io.read_csv(path)
        frame.loc[0, "r2"] += .1
        frame.to_csv(path, index=False)
        with self.assertRaisesRegex(ValueError, "OOF vs saved"):
            self.run_quiet(pls.run_summary, source)
        mock_path = self.root / "flags"
        mock_pls(mock_path, bootstrap=False)
        path = mock_path / "model_comparisons.csv"
        frame = io.read_csv(path)
        frame.loc[0, "supported"] = True
        frame.to_csv(path, index=False)
        with self.assertRaisesRegex(ValueError, "Inference flags"):
            self.run_quiet(pls.run_summary, mock_path)

    def test_pls_stale_interval(self):
        source = self.root / "pls"
        mock_pls(source)
        path = source / "model_comparisons.csv"
        frame = io.read_csv(path)
        frame.loc[0, "ci_high_delta_r2"] += .1
        frame.to_csv(path, index=False)
        with self.assertRaisesRegex(ValueError, "bootstrap comparison interval"):
            self.run_quiet(pls.run_summary, source)

    def test_pls_zero_valid_repetitions_and_unverified_completion(self):
        source = self.root / "zero"
        mock_pls(source, n_valid=0)
        path = source / "bootstrap_metrics.csv"
        frame = io.read_csv(path)[["task", "replicate", "valid"]]
        frame.to_csv(path, index=False)
        result = self.run_quiet(pls.run_summary, source)
        self.assertEqual(set(result["bootstrap_report_status"].values()), {"bootstrap_quality_insufficient"})
        main = io.read_csv(source / "reviewer_summary/pls_main_summary.csv")
        self.assertTrue(main.ci_low_delta_r2.isna().all())
        other = self.root / "unverified"
        mock_pls(other)
        (other / "bootstrap_metrics.csv").unlink()
        result = self.run_quiet(pls.run_summary, other)
        self.assertEqual(set(result["bootstrap_report_status"].values()), {"bootstrap_completion_unverified"})
        self.assertFalse(io.read_csv(other / "reviewer_summary/pls_main_summary.csv").supported_for_report.any())

    def test_pls_duplicate_outer_predictions_rejected(self):
        source = self.root / "duplicate"
        mock_pls(source, bootstrap=False)
        path = source / "oof_predictions.csv"
        frame = io.read_csv(path)
        frame.loc[1, "sub_id"] = frame.loc[0, "sub_id"]
        frame.to_csv(path, index=False)
        with self.assertRaisesRegex(ValueError, "exactly one outer-test"):
            self.run_quiet(pls.run_summary, source)

    def test_word_table_row_counts(self):
        from docx import Document
        for name, builder, reporter, expected in (("bf", mock_bayes, bf.run_summary, [9, 5]),
                                                   ("pls", mock_pls, pls.run_summary, [3, 9])):
            source = self.root / name
            builder(source)
            with contextlib.redirect_stdout(stream.StringIO()):
                reporter(source, mock=True)
            document = Document(source / "reviewer_summary/response_tables.docx")
            self.assertEqual([len(t.rows) for t in document.tables], expected)
            self.assertTrue(any("SIMULATED" in p.text for p in document.paragraphs))


if __name__ == "__main__":
    unittest.main()
