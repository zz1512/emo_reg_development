"""Sample alignment/exclusion checks; synthetic data never represent real results."""

import json
from pathlib import Path
import tempfile
import unittest

import numpy as np
import pandas as pd

from .precheck import PARCELS, PROJECT_ROOT, pfc_candidates, read_behavior, run_precheck, subject_keys


class PrecheckTests(unittest.TestCase):
    def fixture(self, root):
        inputs = root / "inputs"
        inputs.mkdir()
        ids = [f"s{i:05}" for i in range(12)]
        design = pd.DataFrame({"sub_id": ids, "age": np.arange(12) + 7.,
                               "gender": np.arange(12) % 2, "site_id": np.arange(12) % 3,
                               "parent_edu": 14., "meanFD": .2, "in_analysis": 1, "excluded_by": ""})
        design.loc[0, ["in_analysis", "excluded_by"]] = [0, "original_image_missing"]
        design.loc[[1, 2], "meanFD"] = .8
        design.loc[2, "age"] = np.nan
        design.to_csv(inputs / "design.csv", index=False)
        brain = pd.DataFrame(np.random.default_rng(4).normal(size=(12, 432)), columns=PARCELS)
        brain.insert(0, "sub_id", ids)
        brain.loc[3, "432"] = np.nan
        brain.to_csv(inputs / "parcel_response.csv", index=False)
        pd.DataFrame({"parcel": range(1, 433), "name": [f"parcel_{p}" for p in PARCELS]}).to_csv(
            inputs / "parcel_info.csv", index=False)
        (inputs / "tian_combined.json").write_text(json.dumps({"n_cortical": 400}))
        matched_ids = [subject for i, subject in enumerate(ids) if i != 5]
        behavior = root / "behavior.csv"
        pd.DataFrame({"sub_ids": [f"sub-{i:04}" for i in range(len(matched_ids))],
                      "reappraisal_success": ["" if s == ids[4] else str(i / 10) for i, s in enumerate(matched_ids)],
                      "ERQ_CR": ""}).to_csv(behavior, index=False)
        id_source = root / "original.csv"
        id_source.write_text("subject,other_fields\n" + "\n".join(f"{s},unused" for s in matched_ids))
        return inputs, behavior, id_source

    def test_alignment_and_filter_audit_without_erq_requirement(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            inputs, behavior, source = self.fixture(root)
            original = {p: p.read_bytes() for p in inputs.iterdir()}
            output = root / "report"
            result = run_precheck(inputs, behavior, output, id_source=source)
            self.assertEqual(result["tasks"], {"age": {"n": 8, "n_sites": 3}, "success": {"n": 6, "n_sites": 3}})
            roster = pd.read_csv(output / "sample_roster.csv").set_index("sub_id")
            self.assertEqual(roster.loc["s00000", "original_excluded_by"], "original_image_missing")
            self.assertEqual(roster.loc["s00002", "excluded_by_success"], "age_missing_or_nonfinite")
            self.assertEqual(roster.loc["s00004", "excluded_by_success"], "success_missing_or_invalid")
            self.assertEqual(roster.loc["s00005", "excluded_by_success"], "behavior_subject_not_found")
            failures = pd.read_csv(output / "exclusions_by_filter.csv")
            multiple = failures.loc[(failures["task"] == "success") & (failures["sub_id"] == "s00002")]
            self.assertEqual(set(multiple["filter"]), {"age_missing_or_nonfinite", "meanFD_above_0_5"})
            self.assertEqual(multiple["is_first_exclusion"].sum(), 1)
            summary = pd.read_csv(output / "filter_summary.csv")
            for task, n in (("age", 8), ("success", 6)):
                self.assertEqual(summary.loc[summary["task"] == task, "n_newly_excluded"].sum(), 12 - n)
            splits = pd.read_csv(output / "loso_feasibility.csv")
            self.assertTrue(splits["inner_loso_possible"].all())
            self.assertEqual(original, {p: p.read_bytes() for p in inputs.iterdir()})
            with self.assertRaisesRegex(ValueError, "separate precheck"):
                run_precheck(inputs, behavior, inputs, id_source=source)
            # Empty cohorts remain readable and cannot silently proceed to nested LOSO.
            design = pd.read_csv(inputs / "design.csv")
            design["in_analysis"] = 0
            design.to_csv(inputs / "design.csv", index=False)
            empty = run_precheck(inputs, behavior, output, id_source=source)
            self.assertEqual(empty["tasks"]["age"]["n"], 0)
            self.assertEqual(pd.read_csv(output / "site_summary.csv").shape[0], 0)

    def test_invalid_or_duplicate_ids_fail_instead_of_aggregation(self):
        with self.assertRaisesRegex(ValueError, "duplicate"):
            subject_keys(pd.Series(["sub-ABC123", "abc123"]), "test")
        self.assertEqual(subject_keys(pd.Series(["001234.0", "000001"]), "test").tolist(), ["001234", "000001"])
        with tempfile.TemporaryDirectory() as directory:
            _, behavior, source = self.fixture(Path(directory))
            source.write_text("header\nabc123\n")
            with self.assertRaisesRegex(ValueError, "row counts"):
                read_behavior(behavior, "sub_ids", "reappraisal_success", source)

    def test_actual_atlas_labels_are_candidates_not_anatomical_ground_truth(self):
        path = PROJECT_ROOT / "resources/atlases/schaefer2018/Schaefer2018_400Parcels_7Networks_order_info.txt"
        lines = path.read_text().splitlines()
        names = lines[::2] + [f"Tian_{i}" for i in range(401, 433)]
        candidates = pfc_candidates(pd.DataFrame({"parcel": range(1, 433), "name": names}))
        self.assertEqual(candidates["candidate_pfc_named"].sum(), 71)
        self.assertEqual(candidates["candidate_ofc_named"].sum(), 12)
        self.assertEqual(candidates["candidate_pfc_plus_ofc"].sum(), 83)
        self.assertFalse(candidates.loc[candidates["parcel"] > 400, "candidate_pfc_plus_ofc"].any())
        self.assertTrue(candidates["definition_status"].str.contains("not_final").all())


if __name__ == "__main__":
    unittest.main()
