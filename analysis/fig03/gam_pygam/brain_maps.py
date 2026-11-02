"""Partial R² volume, significance mask, and fsLR cortical copies.

The volume keeps every parcel, cortex and Tian subcortex. The mask uses 1 for
q_age < 0.05 and 2 for the marginal parcels. Marginal wins when a parcel is
both. Workbench screenshots are left to the user.
"""

from __future__ import annotations

from pathlib import Path

import nibabel as nib
import numpy as np
import pandas as pd
from nibabel.cifti2 import Cifti2Image
from nibabel.cifti2.cifti2_axes import BrainModelAxis, ScalarAxis

MASK_SIGNIFICANT = 1
MASK_MARGINAL = 2


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
    except TypeError:
        return mni152_to_fslr(image, fslr_density="32k")


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


def write_brain_maps(atlas_path: Path, statistics: pd.DataFrame, marginal: list[int],
                     output_dir: Path, dlabel_path: Path) -> None:
    """Write the MNI volumes and their fsLR cortical GIFTI and dscalar copies."""
    parcels = statistics[statistics["dimension"].astype(str) == "parcel"].copy()
    parcels["parcel"] = parcels["name"].astype(int)
    partial = {
        int(row.parcel): float(row.partial_r2)
        for row in parcels.itertuples()
        if np.isfinite(row.partial_r2)
    }
    significant = {
        int(row.parcel)
        for row in parcels.itertuples()
        if np.isfinite(row.q_age) and float(row.q_age) < 0.05
    }
    marginal_ids = {int(parcel) for parcel in marginal}
    atlas = nib.load(str(atlas_path))
    labels = np.rint(np.asarray(atlas.get_fdata())).astype(np.int32)
    partial_volume, mask_volume = paint_statistic_and_mask(
        labels, partial, significant, marginal_ids)
    output_dir.mkdir(parents=True, exist_ok=True)
    partial_path = output_dir / "partial_r2.nii.gz"
    mask_path = output_dir / "significance_mask.nii.gz"
    nib.save(nib.Nifti1Image(partial_volume, atlas.affine, atlas.header), str(partial_path))
    nib.save(nib.Nifti1Image(mask_volume, atlas.affine, atlas.header), str(mask_path))
    print(f"Brain volumes: {partial_path}")
    print(f"Brain mask: {mask_path}")
    (output_dir / "mask_values.txt").write_text(
        "\n".join([
            "partial_r2.nii.gz: max(0, 1 - RSS_full / RSS_reduced), the fig3 partial R2.",
            "significance_mask.nii.gz: 1 = q_age < 0.05, 2 = marginal parcel.",
            "A parcel that is both significant and marginal is stored as 2.",
            "marginal parcels: " + (", ".join(str(parcel) for parcel in marginal) or "none"),
            "fsLR files are the cortical projection of these two volumes.",
        ]) + "\n",
        encoding="utf-8",
    )
    partial_image = nib.load(str(partial_path))
    mask_image = nib.load(str(mask_path))
    for stem, image, method, map_name in (
        ("partial_r2", partial_image, "linear", "partial R2"),
        ("significance_mask", mask_image, "nearest", "significance mask"),
    ):
        projected = _project_volume(image, method)
        giftis = []
        for hemi, metric in zip(("lh", "rh"), projected):
            path = output_dir / f"{stem}_{hemi}.func.gii"
            nib.save(metric, str(path))
            giftis.append(path)
        dscalar = output_dir / f"{stem}.dscalar.nii"
        _write_cifti(giftis[0], giftis[1], dscalar, dlabel_path, map_name)
        print(f"Brain map: {dscalar}")
