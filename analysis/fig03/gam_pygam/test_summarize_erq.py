"""Validate saved-result extraction; these synthetic results are not study data."""

from __future__ import annotations

import contextlib
import io
import json
from pathlib import Path
import tempfile
import unittest

import pandas as pd

from analysis.fig03.gam_pygam import result_tables, summarize_erq as summary


def mock_statistics():
    """Deterministic all-nonsignificant fixture, never a fitted scientific result."""
    rows = []
    for dimension, (_, count) in summary.FAMILIES.items():
        for unit in range(1, count + 1):
            beta = .0012345678901234567 * (-1 if unit % 2 else 1)
            row = dict(dimension=dimension, name=str(unit), label=f"SIMULATED {dimension} {unit}",
                       status="ok", n_obs=1000 - unit % 3,
                       beta_erq_additive_reference=beta, se_erq_additive_reference=.01,
                       ci_low_erq_additive_reference=beta - .02, ci_high_erq_additive_reference=beta + .02,
                       t_erq_additive_reference=beta / .01,
                       # E3's ERQ coefficient deliberately differs: the extractor must use E2.
                       beta_erq_reference=999., rss_reference_additive=99., rss_reference_full=98.)
            for i, comparison in enumerate(summary.COMPARISONS):
                row.update({f"{field}_{comparison}": value for field, value in {
                    "F": (beta / .01) ** 2 if i == 0 else 1.01,
                    "df1": (1, 4, 5)[i], "df2": row["n_obs"] - (16 if i == 0 else 20),
                    "p": .2 + unit % 5 / 10, "q": .6 + unit % 3 / 10,
                    "reject": False, "partial_r2_reference": .00012345678901234567,
                    "rss_reference_reduced": 100.}.items()})
            rows.append(row)
    return pd.DataFrame(rows)


def write_mock_source(root: Path):
    directory = root / "statistics"
    directory.mkdir(parents=True)
    frame = mock_statistics()
    frame.to_csv(directory / "erq_age_statistics.csv", index=False, float_format="%.17g")
    result_tables.erq_long_table(frame).to_csv(directory / "erq_sequential_comparisons.csv",
                                             index=False, float_format="%.17g")
    (directory / "erq_table_definitions.json").write_text('{"simulated": true}', encoding="utf-8")
    (root / "run_metadata.json").write_text(json.dumps({"analysis_kind": "erq_age_smooth_pygam_v3",
                                                       "simulated": True}), encoding="utf-8")
    return directory


class SummaryTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary.name) / "source"
        self.statistics = write_mock_source(self.root)

    def tearDown(self):
        self.temporary.cleanup()

    def source_frame(self):
        return summary.load_results(self.root, mock=True)[0]

    def run_export(self, **kwargs):
        with contextlib.redirect_stdout(io.StringIO()):
            return summary.run_summary(self.root, mock=True, include_docx=False, **kwargs)

    def test_all_false_counts_precision_and_inputs_untouched(self):
        before = {path: path.read_bytes() for path in self.root.rglob("*") if path.is_file()}
        result = self.run_export(expect_all_nonsignificant=True)
        self.assertTrue(result["all_tests_valid_and_nonsignificant"])
        self.assertEqual(result["row_counts"], {
            "erq_response_summary.csv": 4, "erq_comparison_summary.csv": 12,
            "erq_joint_summary.csv": 4, "erq_primary_details.csv": 898,
            "erq_joint_details.csv": 449, "erq_invalid_tests.csv": 0})
        output = self.root / "reviewer_R1_6"
        compact = summary._read(output / "erq_response_summary.csv")
        self.assertEqual(compact["erq_main_n_valid"].tolist(), [432, 7, 4, 6])
        self.assertEqual(compact["interaction_n_significant"].sum(), 0)
        detailed = summary._read(output / "erq_primary_details.csv")
        main = detailed[detailed.comparison == summary.PRIMARY[0]]
        expected = mock_statistics().iloc[0].beta_erq_additive_reference
        self.assertEqual(main.iloc[0].beta, expected)
        self.assertEqual(main["name"].iloc[:3].tolist(), ["1", "2", "3"])
        self.assertTrue(detailed.loc[detailed.comparison == summary.PRIMARY[1], "beta"].isna().all())
        text = (output / "erq_result_summary.md").read_text(encoding="utf-8")
        self.assertIn("SIMULATED", text)
        self.assertIn("none of the 449 valid tests", text)
        for path, content in before.items():
            self.assertEqual(path.read_bytes(), content)
        self.assertEqual(set(result["source_sha256"]), {str(path) for path in before})

    def test_invalid_not_counted_as_nonsignificant(self):
        frame = self.source_frame()
        index = frame.index[(frame.dimension == "parcel") & (frame.comparison == summary.PRIMARY[0])][0]
        frame.loc[index, "q"] = float("nan")
        validated = summary.validate_results(frame)
        compact, full = summary.aggregate(validated)
        self.assertEqual(compact.iloc[0].erq_main_n_valid, 431)
        self.assertEqual(compact.iloc[0].erq_main_n_invalid, 1)
        english, chinese = summary.result_text(full)
        self.assertIn("none of the 448 valid tests", english)
        self.assertIn("未通过结果校验", chinese)
        self.assertEqual(validated.significant_from_q.sum(), 0)

    def test_failed_status_and_missing_flag_are_invalid(self):
        frame = self.source_frame()
        frame.loc[0, "status"] = "rank_deficient_targets"
        frame["reject"] = frame["reject"].astype(object)
        frame.loc[1, "reject"] = None
        validated = summary.validate_results(frame)
        self.assertEqual((~validated.valid_test).sum(), 2)
        self.assertIn("reject_missing", ";".join(validated.invalid_reason))

    def test_no_valid_tests_does_not_claim_null_evidence(self):
        frame = self.source_frame()
        frame["status"] = "insufficient_observations"
        _, full = summary.aggregate(summary.validate_results(frame))
        english, _ = summary.result_text(full)
        self.assertEqual(english.count("No valid tests were available"), 3)
        self.assertNotIn("none of the", english)

    def test_string_false_and_flag_q_mismatch(self):
        frame = self.source_frame()
        frame["reject"] = "False"
        self.assertTrue(summary.validate_results(frame).valid_test.all())
        frame.loc[0, "q"] = .01
        with self.assertRaisesRegex(ValueError, "reject/q inconsistency"):
            summary.validate_results(frame)

    def test_missing_unit_and_duplicate_rejected(self):
        frame = self.source_frame()
        with self.assertRaisesRegex(ValueError, "expected the same"):
            summary.validate_results(frame.iloc[1:])
        with self.assertRaisesRegex(ValueError, "Duplicate"):
            summary.validate_results(pd.concat([frame, frame.iloc[:1]]))

    def test_stale_long_table_rejected(self):
        path = self.statistics / "erq_sequential_comparisons.csv"
        frame = summary._read(path)
        frame.loc[0, "q"] = .91
        frame.to_csv(path, index=False, float_format="%.17g")
        with self.assertRaisesRegex(ValueError, "Wide/long result mismatch"):
            self.source_frame()

    def test_all_false_assertion_catches_actual_discovery(self):
        path = self.statistics / "erq_age_statistics.csv"
        frame = summary._read(path)
        frame.loc[0, "q_additive_vs_baseline"] = .01
        frame.loc[0, "reject_additive_vs_baseline"] = True
        frame.to_csv(path, index=False, float_format="%.17g")
        (self.statistics / "erq_sequential_comparisons.csv").unlink()
        with self.assertRaisesRegex(ValueError, "contradicted by q"):
            self.run_export(expect_all_nonsignificant=True)
        self.assertFalse((self.root / "reviewer_R1_6").exists())

    def test_source_version_and_output_guard(self):
        with self.assertRaisesRegex(ValueError, "separate reviewer"):
            self.run_export(output=self.statistics)
        with self.assertRaisesRegex(ValueError, "explicitly simulated"):
            summary.load_results(self.root)
        (self.root / "run_metadata.json").write_text('{"analysis_kind":"erq_age_smooth_pygam_v2"}')
        with self.assertRaisesRegex(ValueError, "three-model v3"):
            self.source_frame()

    def test_metadata_alone_marks_simulated_source(self):
        for filename in ("erq_age_statistics.csv", "erq_sequential_comparisons.csv"):
            path = self.statistics / filename
            frame = summary._read(path)
            frame["label"] = frame["label"].str.replace("SIMULATED ", "", regex=False)
            frame.to_csv(path, index=False, float_format="%.17g")
        (self.statistics / "erq_table_definitions.json").unlink()
        with self.assertRaisesRegex(ValueError, "explicitly simulated"):
            summary.load_results(self.root)

    def test_word_three_line_tables(self):
        from docx import Document
        from docx.oxml.ns import qn

        compact, full = summary.aggregate(summary.validate_results(self.source_frame()))
        output = Path(self.temporary.name)
        summary.write_word(compact, full, output, mock=True)
        document = Document(output / "erq_response_table.docx")
        self.assertEqual(len(document.tables), 2)
        self.assertTrue(any("SIMULATED" in p.text for p in document.paragraphs))
        for table in document.tables:
            self.assertEqual(len(table.rows), 5)
            borders = table._tbl.tblPr.find(qn("w:tblBorders"))
            for edge in ("top", "bottom"):
                self.assertEqual(borders.find(qn("w:" + edge)).get(qn("w:val")), "single")
            for edge in ("left", "right", "insideV", "insideH"):
                self.assertEqual(borders.find(qn("w:" + edge)).get(qn("w:val")), "nil")
            self.assertIsNotNone(table.rows[0]._tr.trPr.find(qn("w:tblHeader")))


if __name__ == "__main__":
    unittest.main()
