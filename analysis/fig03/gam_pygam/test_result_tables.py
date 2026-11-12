"""Document mapping/OOXML checks; no model dependencies required."""

from pathlib import Path
import tempfile
import unittest
from zipfile import ZipFile
import xml.etree.ElementTree as ET

import pandas as pd

from analysis.fig03.gam_pygam import result_tables as tables

NS = {"w": "http://schemas.openxmlformats.org/wordprocessingml/2006/main"}


class ResultTableTests(unittest.TestCase):
    def test_three_effect_rows_keep_main_coefficient_and_own_test_df(self):
        row = {"dimension": "parcel", "name": "152", "label": "Example", "n_obs": 120, "status": "ok",
               "beta_erq_additive_reference": .12345678901234567, "se_erq_additive_reference": .02,
               "ci_low_erq_additive_reference": .08, "ci_high_erq_additive_reference": .16,
               "t_erq_additive_reference": 6.172839450617284}
        for comparison, _, _ in tables.EFFECTS:
            row.update({f"p_{comparison}": 1.234567890123456e-7, f"q_{comparison}": .001,
                        f"df1_{comparison}": 1 if comparison == "additive_vs_baseline" else 4,
                        f"df2_{comparison}": 109 if comparison == "additive_vs_baseline" else 105,
                        f"F_{comparison}": 3.1, f"partial_r2_reference_{comparison}": .02})
        source = pd.DataFrame([row])
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            output = tables.write_erq_tables(source, root, mock=True)
            self.assertEqual(len(output), 3)
            self.assertEqual(output.iloc[0]["df2"], 109)
            self.assertEqual(output.iloc[1]["df2"], 105)
            self.assertEqual(output.iloc[0]["beta"], row["beta_erq_additive_reference"])
            self.assertTrue(output.iloc[1:]["beta"].isna().all())
            saved = pd.read_csv(root / "erq_sequential_comparisons.csv", float_precision="round_trip")
            self.assertEqual(saved.iloc[0]["p"], row["p_additive_vs_baseline"])
            with ZipFile(root / "erq_manuscript_table.docx") as archive:
                document = ET.fromstring(archive.read("word/document.xml"))
                styles = ET.fromstring(archive.read("word/styles.xml"))
            self.assertEqual(len(document.findall(".//w:tbl/w:tr", NS)), 4)
            self.assertIsNotNone(document.find(".//w:tblHeader", NS))
            borders = document.find(".//w:tblBorders", NS)
            for name in ("left", "right", "insideH", "insideV"):
                self.assertEqual(borders.find(f"w:{name}", NS).get(f"{{{NS['w']}}}val"), "nil")
            title = next(element for element in styles if element.get(f"{{{NS['w']}}}styleId") == "Title")
            self.assertIsNone(title.find("w:pPr/w:pBdr", NS))

    def test_failed_bayesian_units_are_retained_with_missing_evidence(self):
        source = pd.DataFrame([{"dimension": "parcel", "name": "1", "label": "Constant synthetic response",
                                "n_obs": 20, "status": "constant_age_or_response"}])
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            output = tables.write_age_tables(source, root, mock=True)
            self.assertEqual(len(output), 3)
            self.assertTrue(output["bf10"].isna().all())
            with ZipFile(root / "age_bayes_manuscript_table.docx") as archive:
                xml = archive.read("word/document.xml").decode()
                self.assertIn("constant_age_or_response", xml)
                self.assertNotIn(">nan<", xml)


if __name__ == "__main__":
    unittest.main()
