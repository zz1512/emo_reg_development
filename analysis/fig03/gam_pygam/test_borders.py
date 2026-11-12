"""Binary border definitions and real Workbench conversion on the fsLR32k mesh.

Native atlas vertex assignments are checked directly. Border integration
tests invoke wb_command and read its generated files back.
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

ATLAS_DIR = maps.SURFACE_DIR.parents[1] / "atlases" / "schaefer2018"
NATIVE = ATLAS_DIR / "Schaefer2018_400Parcels_7Networks_order.dlabel.nii"


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
                patch.object(maps, "load_surface_atlas") as load, \
                self.assertRaisesRegex(FileNotFoundError, "--wb-command"):
            maps.write_brain_borders(Path("unused"), pd.DataFrame(), [], Path("unused"), NATIVE)
        load.assert_not_called()


class NativeSurfaceTests(unittest.TestCase):
    def test_surface_maps_masks_and_selected_labels_match_native_at_every_vertex(self):
        from analysis.fig03.gam_pygam.make_selected_dlabel import make_selected_dlabel, SELECTED

        native, labels = maps.load_surface_atlas(NATIVE)
        stats = pd.DataFrame({"dimension": "parcel", "name": np.arange(1, 433),
                              "partial_r2_m3_vs_m1": np.arange(1, 433) / 10000.,
                              "partial_r2_m3_vs_m2": np.arange(1, 433) / 20000.,
                              "q_m3_vs_m1": .8, "q_m3_vs_m2": .9})
        stats.loc[stats.name.isin([152, 154, 156]), "q_m3_vs_m1"] = .01
        stats.loc[stats.name == 185, "q_m3_vs_m2"] = .02
        original = stats.copy(deep=True)
        with tempfile.TemporaryDirectory(prefix="fig03-native-maps-") as temporary:
            root = Path(temporary)
            atlas = root / "combined.nii.gz"
            voxels = np.arange(433, dtype=np.int16).reshape(433, 1, 1)
            nib.save(nib.Nifti1Image(voxels, np.eye(4)), str(atlas))
            with patch.object(maps, "write_brain_borders") as borders, contextlib.redirect_stdout(io.StringIO()):
                maps.write_brain_maps(atlas, stats, [156, 157], root / "maps", NATIVE)
                selected_path = make_selected_dlabel(atlas, NATIVE, root / "selected")
            self.assertEqual(borders.call_args.args[4], NATIVE)
            for comparison, divisor, significant in (("m3_vs_m1", 10000., [152, 154, 156]),
                                                     ("m3_vs_m2", 20000., [185])):
                expected_scalar = np.asarray(native.dataobj, dtype=np.float32) / divisor
                scalar = nib.load(root / f"maps/partial_r2_{comparison}_fslr32k.dscalar.nii")
                self.assertEqual(int(scalar.nifti_header['intent_code']), 3006)
                np.testing.assert_array_equal(np.asarray(scalar.dataobj), expected_scalar)
                volume = nib.load(root / f"maps/partial_r2_{comparison}_volume_MNI152.nii.gz")
                np.testing.assert_array_equal(np.asarray(volume.dataobj), (voxels / divisor).astype(np.float32))
                for hemi, native_labels, structure in zip(("lh", "rh"), labels, ("CortexLeft", "CortexRight")):
                    partial = nib.load(root / f"maps/partial_r2_{comparison}_fslr32k_{hemi}.func.gii")
                    np.testing.assert_array_equal(partial.darrays[0].data, (native_labels / divisor).astype(np.float32))
                    self.assertEqual(partial.meta['AnatomicalStructurePrimary'], structure)
                    mask = nib.load(root / f"maps/significance_mask_{comparison}_fslr32k_{hemi}.func.gii")
                    expected_mask = np.isin(native_labels, significant).astype(int)
                    expected_mask[np.isin(native_labels, [156, 157])] = 2
                    np.testing.assert_array_equal(mask.darrays[0].data, expected_mask)
            selected = nib.load(selected_path)
            source = np.asarray(native.dataobj)
            np.testing.assert_array_equal(np.asarray(selected.dataobj), np.where(np.isin(source, list(SELECTED)), source, 0))
            self.assertEqual(int(selected.nifti_header['intent_code']), 3007)
            self.assertEqual(selected.header.get_axis(1), native.header.get_axis(1))
            for key, (_, rgb) in SELECTED.items():
                np.testing.assert_allclose(selected.header.get_axis(0).label[0][key][1], [*[c/255. for c in rgb], 1.])
            manifest = json.loads((root / "maps/surface_manifest.json").read_text())
            self.assertIn("no volume projection", manifest["surface_assignment"])
        pd.testing.assert_frame_equal(stats, original)

    def test_wrong_cortical_resolution_is_rejected(self):
        with tempfile.TemporaryDirectory(prefix="fig03-native-wrong-resolution-") as temporary:
            atlas = Path(temporary) / "combined432.nii.gz"
            nib.save(nib.Nifti1Image(np.arange(433, dtype=np.int16).reshape(433, 1, 1), np.eye(4)), str(atlas))
            other = ATLAS_DIR / "Schaefer2018_200Parcels_7Networks_order.dlabel.nii"
            with self.assertRaisesRegex(ValueError, "--schaefer-atlas"):
                maps._native_labels_for_volume(other, atlas)


@unittest.skipUnless(shutil.which("wb_command"), "Workbench is not installed")
class WorkbenchBorderTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        _, cls.labels = maps.load_surface_atlas(NATIVE)

    def _write(self, root, marginal):
        atlas = root / "atlas.nii.gz"
        nib.save(nib.Nifti1Image(np.arange(433, dtype=np.int16).reshape(433, 1, 1), np.eye(4)), str(atlas))
        stats = pd.DataFrame({
            "dimension": "parcel", "name": [1, 156, 157, 201],
            "q_m3_vs_m1": [.01, .01, .3, .9],
            "q_m3_vs_m2": [.9, .3, .01, .01],
            "partial_r2_m3_vs_m1": .1, "partial_r2_m3_vs_m2": .02,
        })
        with contextlib.redirect_stdout(io.StringIO()):
            maps.write_brain_borders(atlas, stats, marginal, root / "maps", NATIVE)
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

    def test_mismatched_native_grid_is_rejected_before_conversion(self):
        with tempfile.TemporaryDirectory(prefix="fig03-border-grid-") as temporary:
            root = Path(temporary)
            atlas = root / "atlas.nii.gz"
            nib.save(nib.Nifti1Image(np.zeros((2, 2, 2), dtype=np.int16), np.eye(4)), str(atlas))
            stats = pd.DataFrame({"dimension": ["parcel"], "name": [1], "q_age": [.01], "partial_r2": [.1]})
            with patch.object(maps, "_native_labels_for_volume", return_value=(None, (np.ones(8), np.ones(8)))), \
                    patch.object(maps.subprocess, "run") as command, \
                    self.assertRaisesRegex(ValueError, "matching fsLR32k"):
                maps.write_brain_borders(atlas, stats, [], root / "maps", NATIVE)
            command.assert_not_called()


if __name__ == "__main__":
    unittest.main()
