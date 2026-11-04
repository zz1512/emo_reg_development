"""Comparison-specific partial R² volumes, masks, and fsLR cortical copies.

The volumes keep cortex and Tian subcortex. Each comparison's mask uses 1
for its BH q < 0.05 and 2 for manually selected marginal parcels. Marginal
labels take priority. File names identify the comparison and spatial domain.
GIFTI metadata explicitly identifies the hemisphere for Workbench.
Independent binary FDR ROIs and marginal-parcel ROIs produce red Workbench
borders in maps/borders; the legacy two-class mask is not their source.
"""

from __future__ import annotations

from pathlib import Path
import json
import shutil
import subprocess
import xml.etree.ElementTree as ET

import nibabel as nib
import numpy as np
import pandas as pd
from nibabel.cifti2 import Cifti2Image
from nibabel.cifti2.cifti2_axes import BrainModelAxis, ScalarAxis

MASK_SIGNIFICANT = 1
MASK_MARGINAL = 2
SURFACE_DIR = Path(__file__).resolve().parents[3] / "resources" / "templates" / "fslr_32k"
DEFAULT_BORDER_SURFACES = tuple(
    SURFACE_DIR / f"fsaverage.{hemi}.inflated.32k_fs_LR.surf.gii" for hemi in ("L", "R"))


def paint_statistic_and_mask(labels: np.ndarray, partial_r2: dict[int, float],
                             significant: set[int], marginal: set[int],
                             ) -> tuple[np.ndarray, np.ndarray]:
    """Paint partial R² and the two-class mask onto atlas labels."""
    labels = np.asarray(labels)
    volume = np.zeros(labels.shape, dtype=np.float32)
    mask = np.zeros(labels.shape, dtype=np.int16)
    for parcel, value in partial_r2.items():
        if np.isfinite(value):
            volume[labels == int(parcel)] = np.float32(value)
    for parcel in significant:
        if int(parcel) in marginal:
            continue
        mask[labels == int(parcel)] = MASK_SIGNIFICANT
    for parcel in marginal:
        mask[labels == int(parcel)] = MASK_MARGINAL
    return volume, mask


def _project_volume(image: nib.spatialimages.SpatialImage, method: str):
    from neuromaps.transforms import mni152_to_fslr

    try:
        return mni152_to_fslr(image, fslr_density="32k", method=method)
    except TypeError as exc:
        raise RuntimeError(
            f"fsLR projection failed with interpolation={method!r}. "
            "Use a neuromaps mni152_to_fslr implementation accepting the method argument. "
            "The requested interpolation was not replaced by a default."
        ) from exc


def _gifti_matrix(path: Path) -> np.ndarray:
    """One row per GIFTI data array, one column per vertex."""
    image = nib.load(str(path))
    if not isinstance(image, nib.GiftiImage) or len(image.darrays) == 0:
        raise ValueError(f"{path} is not a GIFTI metric with data arrays")
    rows = [np.asarray(array.data).reshape(-1) for array in image.darrays]
    width = rows[0].size
    if any(row.size != width for row in rows):
        raise ValueError(f"{path} data arrays do not share one vertex count")
    return np.vstack(rows)


def _write_cifti(left_path: Path, right_path: Path, output_path: Path,
                 template: Path, map_name: str) -> None:
    """Place fsLR metrics on the Schaefer dlabel grayordinate axis.

    Subcortical grayordinates stay 0. The dscalar is written with nibabel, so
    this step does not import neuro_utils.
    """
    template_image = nib.load(str(template))
    brain_axis = None
    for index in range(template_image.ndim):
        axis = template_image.header.get_axis(index)
        if isinstance(axis, BrainModelAxis):
            brain_axis = axis
            break
    if brain_axis is None:
        raise ValueError(f"{template} has no CIFTI brain model axis")
    matrices = {
        "CIFTI_STRUCTURE_CORTEX_LEFT": _gifti_matrix(left_path),
        "CIFTI_STRUCTURE_CORTEX_RIGHT": _gifti_matrix(right_path),
    }
    n_maps = next(iter(matrices.values())).shape[0]
    data = np.zeros((n_maps, len(brain_axis)), dtype=np.float32)
    for structure_name, cifti_slice, sub_axis in brain_axis.iter_structures():
        matrix = matrices.get(structure_name)
        if matrix is None:
            continue
        vertices = np.asarray(sub_axis.vertex)
        if matrix.shape[1] != int(brain_axis.nvertices[structure_name]):
            raise ValueError(
                f"{structure_name} has {matrix.shape[1]} GIFTI vertices and "
                f"{int(brain_axis.nvertices[structure_name])} template vertices"
            )
        data[:, cifti_slice] = matrix[:, vertices]
    image = Cifti2Image(data, (ScalarAxis([map_name] * n_maps), brain_axis))
    nib.save(image, str(output_path))


def _parcel_comparisons(statistics: pd.DataFrame) -> tuple[pd.DataFrame, list[tuple[str, str, str]]]:
    parcels = statistics[statistics["dimension"].astype(str) == "parcel"].copy()
    parcels["parcel"] = parcels["name"].astype(int)
    if parcels["parcel"].duplicated().any():
        raise ValueError("Brain maps require one statistics row per parcel")
    comparisons = []
    for comparison in ("m3_vs_m1", "m3_vs_m2"):
        q_column = f"q_{comparison}"
        partial_column = f"partial_r2_{comparison}"
        if comparison == "m3_vs_m1":
            if q_column not in parcels.columns:
                q_column = "q_age"
            if partial_column not in parcels.columns:
                partial_column = "partial_r2"
        available = (q_column in parcels.columns, partial_column in parcels.columns)
        if comparison == "m3_vs_m2" and not any(available):
            print("Saved statistics have no m3_vs_m2 comparison; refit to write its brain maps.")
            continue
        if not all(available):
            missing = [column for column in (q_column, partial_column) if column not in parcels.columns]
            raise ValueError(f"Brain maps for {comparison} require columns: {', '.join(missing)}")
        comparisons.append((comparison, q_column, partial_column))
    return parcels, comparisons


def _surface_mesh(path: Path, hemisphere: str) -> tuple[np.ndarray, np.ndarray]:
    image = nib.load(str(path))
    if not isinstance(image, nib.GiftiImage):
        raise ValueError(f"{path} is not a GIFTI surface")
    points = image.get_arrays_from_intent("NIFTI_INTENT_POINTSET")
    triangles = image.get_arrays_from_intent("NIFTI_INTENT_TRIANGLE")
    if len(points) != 1 or len(triangles) != 1:
        raise ValueError(f"{path} requires one pointset and one triangle array")
    expected = "CortexLeft" if hemisphere == "lh" else "CortexRight"
    structure = image.meta.get("AnatomicalStructurePrimary", points[0].meta.get("AnatomicalStructurePrimary"))
    if structure != expected:
        raise ValueError(f"{path} has structure {structure!r}; {hemisphere} requires {expected}")
    return np.asarray(points[0].data), np.asarray(triangles[0].data)


def _canonical_triangles(triangles: np.ndarray) -> np.ndarray:
    triangles = np.sort(triangles, axis=1)
    return triangles[np.lexsort(triangles.T)]


def _border_surfaces(surfaces: tuple[Path, Path] | None) -> tuple[Path, Path]:
    """Check hemisphere, vertex count and indexed topology against fsLR32k."""
    paths = DEFAULT_BORDER_SURFACES if surfaces is None else tuple(Path(path) for path in surfaces)
    if len(paths) != 2:
        raise ValueError("Border surfaces must be supplied as LEFT RIGHT")
    for hemi, path, template in zip(("lh", "rh"), paths, DEFAULT_BORDER_SURFACES):
        points, triangles = _surface_mesh(path, hemi)
        canonical_points, canonical_mesh = _surface_mesh(template, hemi)
        if points.shape != canonical_points.shape or not np.array_equal(
                _canonical_triangles(triangles), _canonical_triangles(canonical_mesh)):
            raise ValueError(f"{path} does not use the bundled fsLR32k indexed mesh")
    return paths


def _color_border_red(path: Path) -> int:
    """Set Workbench v3 name/class colors; drawing width belongs to the viewer."""
    tree = ET.parse(path)
    root = tree.getroot()
    if root.tag != "BorderFile" or root.get("Version") != "3":
        raise ValueError(f"{path}: expected a Workbench version 3 BorderFile")
    for element in (*root.findall("Class"), *root.findall(".//Border")):
        element.set("Red", "1")
        element.set("Green", "0")
        element.set("Blue", "0")
    tree.write(path, encoding="utf-8", xml_declaration=True)
    return len(root.findall(".//BorderPart"))


def write_brain_borders(atlas_path: Path, statistics: pd.DataFrame, marginal: list[int],
                        output_dir: Path, border_surfaces: tuple[Path, Path] | None = None,
                        wb_command: str = "wb_command") -> None:
    """Write independent FDR and marginal ROIs/borders without refitting.

    Nearest-neighbor projection of atlas IDs is shared across all binary ROIs:
    applying a parcel-membership mask commutes with this discrete sampling.
    Each marginal parcel gets its own column so adjacent 156/157 retain their
    individual outlines, in a file separate from the two FDR comparisons.
    """
    executable = shutil.which(str(wb_command))
    if executable is None:
        raise FileNotFoundError("Workbench wb_command was not found. Add it to PATH or use --wb-command PATH.")
    surfaces = _border_surfaces(border_surfaces)
    parcels, comparisons = _parcel_comparisons(statistics)
    parcel_ids = parcels["parcel"].to_numpy(dtype=int)
    marginal_ids = sorted(set(map(int, marginal)))
    if not set(marginal_ids) <= set(parcel_ids):
        raise ValueError("Marginal border parcels are missing from the saved parcel statistics")
    rois = []
    for comparison, q_column, _ in comparisons:
        q = parcels[q_column].to_numpy(dtype=float)
        significant = set(parcel_ids[np.isfinite(q) & (q < 0.05)])
        rois.append((f"fdr_significant_{comparison}", f"FDR_{comparison}",
                     [(f"FDR_q_lt_0.05_{comparison}", significant)]))
    if marginal_ids:
        suffix = "_".join(map(str, marginal_ids))
        rois.append((f"marginal_{suffix}", f"Manual_marginal_{suffix}",
                     [(f"Manual_marginal_parcel_{parcel}", {parcel}) for parcel in marginal_ids]))
    atlas = nib.load(str(atlas_path))
    projected = tuple(_project_volume(atlas, "nearest"))
    if len(projected) != 2:
        raise ValueError("Border projection requires left and right fsLR metrics")
    labels = []
    for hemi, metric, surface in zip(("lh", "rh"), projected, surfaces):
        if len(metric.darrays) != 1:
            raise ValueError("The projected parcel atlas must have exactly one metric column")
        values = np.asarray(metric.darrays[0].data).reshape(-1)
        points, _ = _surface_mesh(surface, hemi)
        if len(values) != len(points) or not np.isfinite(values).all() or not np.allclose(values, np.rint(values)):
            raise ValueError(f"Invalid nearest-neighbor parcel projection for {hemi}; requires matching fsLR32k vertices and integer labels")
        labels.append(np.rint(values).astype(np.int32))
    border_dir = output_dir / "borders"
    border_dir.mkdir(parents=True, exist_ok=True)
    records = []
    for stem, class_name, columns in rois:
        for hemi, surface, vertex_labels in zip(("lh", "rh"), surfaces, labels):
            metric = nib.GiftiImage()
            metric.meta["AnatomicalStructurePrimary"] = "CortexLeft" if hemi == "lh" else "CortexRight"
            counts = {}
            for name, ids in columns:
                data = np.isin(vertex_labels, list(ids)).astype(np.float32)
                array = nib.gifti.GiftiDataArray(data, intent="NIFTI_INTENT_SHAPE")
                array.meta["Name"] = name
                metric.add_gifti_data_array(array)
                counts[name] = int(data.sum())
            metric_path = border_dir / f"{stem}_fslr32k_{hemi}.func.gii"
            border_path = border_dir / f"{stem}_fslr32k_{hemi}.border"
            nib.save(metric, str(metric_path))
            command = [executable, "-metric-rois-to-border", str(surface), str(metric_path),
                       class_name, str(border_path)]
            if len(columns) == 1:
                command += ["-column", "1"]
            try:
                subprocess.run(command, check=True, capture_output=True, text=True)
            except subprocess.CalledProcessError as exc:
                raise RuntimeError(f"Workbench border conversion failed for {metric_path}: {exc.stderr or exc.stdout}") from exc
            parts = _color_border_red(border_path)
            records.append({"border_file": border_path.name, "roi_file": metric_path.name,
                            "hemisphere": hemi, "surface": str(surface.resolve()),
                            "vertices_per_roi": counts, "border_parts": parts})
            print(f"Brain border: {border_path} ({parts} parts)")
    (border_dir / "border_manifest.json").write_text(json.dumps({
        "definition": "FDR borders use q < 0.05 directly, including any significant marginal parcel; manual marginal borders are independent annotations.",
        "source_atlas": str(atlas_path.resolve()),
        "projection": "nearest-neighbor parcel IDs, then binary membership; same MNI-to-fsLR32k mapping as the masks",
        "color": "red (name and class colors)",
        "workbench_display": "Display Borders; Polylines; Standard Color = Red; Line Diameter 2-3 px. Width is a viewer/scene setting.",
        "empty_rois": "A hemisphere with no ROI vertices produces a valid empty border file; see vertices_per_roi and border_parts.",
        "files": records,
    }, indent=2), encoding="utf-8")


def write_brain_maps(atlas_path: Path, statistics: pd.DataFrame, marginal: list[int],
                     output_dir: Path, dlabel_path: Path,
                     border_surfaces: tuple[Path, Path] | None = None,
                     wb_command: str = "wb_command") -> None:
    """Write comparison-specific maps plus independent cortical ROI borders.

    Old q_age/partial_r2 tables can still produce m3_vs_m1 maps. A saved fit
    must contain both new m3_vs_m2 columns to produce that comparison.
    """
    parcels, comparisons = _parcel_comparisons(statistics)
    marginal_ids = {int(parcel) for parcel in marginal}
    atlas = nib.load(str(atlas_path))
    labels = np.rint(np.asarray(atlas.get_fdata())).astype(np.int32)
    output_dir.mkdir(parents=True, exist_ok=True)
    projections = []
    for comparison, q_column, partial_column in comparisons:
        partial_values = parcels[partial_column].to_numpy(dtype=float)
        q_values = parcels[q_column].to_numpy(dtype=float)
        parcel_ids = parcels["parcel"].to_numpy(dtype=int)
        partial = {int(parcel): float(value) for parcel, value in zip(parcel_ids, partial_values)
                   if np.isfinite(value)}
        significant = {int(parcel) for parcel, q_value in zip(parcel_ids, q_values)
                       if np.isfinite(q_value) and q_value < 0.05}
        partial_volume, mask_volume = paint_statistic_and_mask(
            labels, partial, significant, marginal_ids)
        for kind, volume, method, description in (
            ("partial_r2", partial_volume, "linear", "partial R2"),
            ("significance_mask", mask_volume, "nearest", "significance mask"),
        ):
            stem = f"{kind}_{comparison}"
            volume_path = output_dir / f"{stem}_volume_MNI152.nii.gz"
            image = nib.Nifti1Image(volume, atlas.affine, atlas.header)
            image.set_data_dtype(volume.dtype)
            nib.save(image, str(volume_path))
            print(f"Brain volume: {volume_path}")
            projections.append((stem, volume_path, method, f"{comparison} {description}"))
    (output_dir / "mask_values.txt").write_text(
        "\n".join([
            "comparisons written: " + ", ".join(row[0] for row in comparisons),
            "partial_r2_m3_vs_m1: max(0, 1 - RSS_m3 / RSS_m1), the fig3 overall-age partial R2.",
            "partial_r2_m3_vs_m2: max(0, 1 - RSS_m3 / RSS_m2), improvement beyond linear age.",
            "significance_mask_<comparison>: 1 = q_<comparison> < 0.05; 2 = manually selected marginal parcel.",
            "Value 2 is an annotation, not a new significant finding from either comparison.",
            "Manual marginal annotation takes priority over value 1 in both comparisons.",
            "manually selected marginal parcels: " + (", ".join(str(parcel) for parcel in marginal) or "none"),
            "*_volume_MNI152.nii.gz keeps every parcel, including Tian subcortex, on the combined atlas grid.",
            "*_fslr32k_lh/rh.func.gii and *_fslr32k.dscalar.nii are cortical projections.",
            "GIFTI AnatomicalStructurePrimary: lh = CortexLeft; rh = CortexRight.",
            "Noncortical grayordinates in the dscalar stay 0; subcortical results remain in the volumes.",
            "Projection interpolation: partial_r2 = linear; significance_mask = nearest.",
            "borders/ has independent binary FDR ROIs (q < 0.05 directly) and manual marginal ROIs, plus red .border files.",
            "The combined 0/1/2 mask is not used as a binary border ROI; adjacent marginal parcels keep separate columns/borders.",
        ]) + "\n",
        encoding="utf-8",
    )
    # Write every volume before projection so a surface dependency failure
    # still leaves both comparisons available for volume-space inspection.
    for stem, volume_path, method, map_name in projections:
        projected = _project_volume(nib.load(str(volume_path)), method)
        giftis = []
        for hemi, metric in zip(("lh", "rh"), projected):
            path = output_dir / f"{stem}_fslr32k_{hemi}.func.gii"
            # Workbench reads the GIFTI file-level structure, not the lh/rh
            # filename suffix. Preserve other metadata from the projection.
            metric.meta["AnatomicalStructurePrimary"] = (
                "CortexLeft" if hemi == "lh" else "CortexRight")
            nib.save(metric, str(path))
            giftis.append(path)
        dscalar = output_dir / f"{stem}_fslr32k.dscalar.nii"
        _write_cifti(giftis[0], giftis[1], dscalar, dlabel_path, map_name)
        print(f"Brain map: {dscalar}")
    write_brain_borders(atlas_path, statistics, marginal, output_dir, border_surfaces, wb_command)
