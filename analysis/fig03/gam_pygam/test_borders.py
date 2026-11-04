"""Binary border definitions and real Workbench conversion on the fsLR32k mesh.

Only the MNI-to-surface transform is replaced by known parcel labels. The
border integration tests invoke wb_command and read its generated files back.
"""

from __future__ import annotations

import contextlib
import io
import json
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest
from unittest.mock import patch
import xml.etree.ElementTree as ET

import nibabel as nib
import numpy as np
import pandas as pd

from analysis.fig03.gam_pygam import brain_maps as maps


class BorderInputTests(unittest.TestCase):
    def test_swapped_hemispheres_are_rejected(self):
        with self.assertRaisesRegex(ValueError, "requires CortexLeft"):
            maps._border_surfaces(tuple(reversed(maps.DEFAULT_BORDER_SURFACES)))

    def test_equal_vertex_count_does_not_accept_a_different_mesh(self):
        surface = nib.load(str(maps.DEFAULT_BORDER_SURFACES[0]))
        triangles = surface.get_arrays_from_intent("NIFTI_INTENT_TRIANGLE")[0].data
        triangles[0, 0] = (triangles[0, 0] + 10) % 32492
        with tempfile.TemporaryDirectory(prefix="fig03-bad-surface-") as temporary:
            path = Path(temporary) / "wrong_mesh.surf.gii"
            nib.save(surface, str(path))
            with self.assertRaisesRegex(ValueError, "indexed mesh"):
                maps._border_surfaces((path, maps.DEFAULT_BORDER_SURFACES[1]))

    def test_missing_workbench_gives_an_actionable_error(self):
        with patch.object(maps.shutil, "which", return_value=None), \
                patch.object(maps, "_project_volume") as project, \
                self.assertRaisesRegex(FileNotFoundError, "--wb-command"):
            maps.write_brain_borders(Path("unused"), pd.DataFrame(), [], Path("unused"))
        project.assert_not_called()


@unittest.skipUnless(shutil.which("wb_command"), "Workbench is not installed")
class WorkbenchBorderTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        # Known surface parcel IDs; no real subject data or registration download.
        template = maps.SURFACE_DIR.parents[1] / "atlases" / "schaefer2018" / "Schaefer2018_400Parcels_7Networks_order.dlabel.nii"
        image = nib.load(str(template))
        values = np.asarray(image.dataobj)[0]
        cls.labels = []
        for _, section, axis in image.header.get_axis(1).iter_structures():
            data = np.zeros(32492, dtype=np.float32)
            data[axis.vertex] = values[section]
            cls.labels.append(data)

    def _write(self, root, marginal):
        atlas = root / "atlas.nii.gz"
        nib.save(nib.Nifti1Image(np.array([1, 156, 157, 201], dtype=np.int16).reshape(4, 1, 1), np.eye(4)), str(atlas))
        stats = pd.DataFrame({
            "dimension": "parcel", "name": [1, 156, 157, 201],
            "q_m3_vs_m1": [.01, .01, .3, .9],
            "q_m3_vs_m2": [.9, .3, .01, .01],
            "partial_r2_m3_vs_m1": .1, "partial_r2_m3_vs_m2": .02,
        })
        projected = tuple(nib.GiftiImage(darrays=[nib.gifti.GiftiDataArray(data)]) for data in self.labels)
        with patch.object(maps, "_project_volume", return_value=projected) as project, \
                contextlib.redirect_stdout(io.StringIO()):
            maps.write_brain_borders(atlas, stats, marginal, root / "maps")
        self.assertEqual(project.call_count, 1)
        self.assertEqual(project.call_args.args[1], "nearest")
        return root / "maps" / "borders"

    def test_fdr_and_manual_rois_remain_independent_and_workbench_reads_distinct_colors(self):
        with tempfile.TemporaryDirectory(prefix="fig03-borders-wb-") as temporary:
            output = self._write(Path(temporary), [156, 157])
            self.assertEqual(len(list(output.glob("*.border"))), 6)
            fdr = nib.load(str(output / "fdr_significant_m3_vs_m1_fslr32k_lh.func.gii"))
            np.testing.assert_array_equal(fdr.darrays[0].data, np.isin(self.labels[0], [1, 156]))
            other = nib.load(str(output / "fdr_significant_m3_vs_m2_fslr32k_lh.func.gii"))
            np.testing.assert_array_equal(other.darrays[0].data, self.labels[0] == 157)
            manual = nib.load(str(output / "marginal_156_157_fslr32k_lh.func.gii"))
            self.assertEqual(len(manual.darrays), 2)
            for column, parcel in zip(manual.darrays, (156, 157)):
                np.testing.assert_array_equal(column.data, self.labels[0] == parcel)
            root = ET.parse(output / "marginal_156_157_fslr32k_lh.border").getroot()
            self.assertEqual({element.get("Name") for element in root.findall(".//Border")},
                             {"Manual_marginal_parcel_156", "Manual_marginal_parcel_157"})
            for path in output.glob("*.border"):
                root = ET.parse(path).getroot()
                self.assertEqual(root.get("Structure"), "CORTEX_LEFT" if "_lh." in path.name else "CORTEX_RIGHT")
                self.assertEqual(int(root.get("SurfaceNumberOfVertices")), 32492)
                expected_color = (0., 229. / 255., 1.) if path.name.startswith("marginal_") else (1., 0., 0.)
                for element in (*root.findall("Class"), *root.findall(".//Border")):
                    np.testing.assert_allclose(
                        [float(element.get(color)) for color in ("Red", "Green", "Blue")],
                        expected_color, atol=1e-8, rtol=0)
                result = subprocess.run([shutil.which("wb_command"), "-file-information", str(path)],
                                        check=True, capture_output=True, text=True)
                self.assertIn("Border", result.stdout)
            manifest = json.loads((output / "border_manifest.json").read_text())
            self.assertEqual(manifest["colors"], {"fdr_significant": "#FF0000", "manual_marginal": "#00E5FF"})
            self.assertIn("Coloring = Class", manifest["workbench_display"])
            for record in manifest["files"]:
                expected_color = (0., 229. / 255., 1.) if record["border_file"].startswith("marginal_") else (1., 0., 0.)
                np.testing.assert_allclose(record["color_rgb"], expected_color)
            empty = next(row for row in manifest["files"] if row["border_file"] == "marginal_156_157_fslr32k_rh.border")
            self.assertEqual(empty["border_parts"], 0)
            self.assertEqual(sum(empty["vertices_per_roi"].values()), 0)

    def test_no_marginal_input_produces_only_comparison_borders(self):
        with tempfile.TemporaryDirectory(prefix="fig03-borders-parcel-only-") as temporary:
            output = self._write(Path(temporary), [])
            self.assertEqual(len(list(output.glob("*.border"))), 4)
            self.assertEqual(len(list(output.glob("marginal*"))), 0)

    def test_mismatched_projected_grid_is_rejected_before_conversion(self):
        with tempfile.TemporaryDirectory(prefix="fig03-border-grid-") as temporary:
            root = Path(temporary)
            atlas = root / "atlas.nii.gz"
            nib.save(nib.Nifti1Image(np.zeros((2, 2, 2), dtype=np.int16), np.eye(4)), str(atlas))
            stats = pd.DataFrame({"dimension": ["parcel"], "name": [1], "q_age": [.01], "partial_r2": [.1]})
            projected = tuple(nib.GiftiImage(darrays=[nib.gifti.GiftiDataArray(np.ones(8, dtype=np.float32))]) for _ in range(2))
            with patch.object(maps, "_project_volume", return_value=projected), \
                    patch.object(maps.subprocess, "run") as command, \
                    self.assertRaisesRegex(ValueError, "matching fsLR32k"):
                maps.write_brain_borders(atlas, stats, [], root / "maps")
            command.assert_not_called()


if __name__ == "__main__":
    unittest.main()
