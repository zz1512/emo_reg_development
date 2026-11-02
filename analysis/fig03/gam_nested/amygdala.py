"""Six amygdala ROIs from the Juelich maximum-probability atlas.

The atlas is Amunts et al. (2005) as packaged by FSL: Juelich-maxprob-thr25-2mm.
Each of Amygdala_LB, Amygdala_CM, and Amygdala_SF is one bilateral label. Left
and right are split by the sign of the MNI x coordinate. The resulting masks
stay in atlas space; the volume reader resamples them once onto the effect-size
grid.
"""

from __future__ import annotations

import os
import re
import xml.etree.ElementTree as ET
from pathlib import Path

import nibabel as nib
import numpy as np

NUCLEI = ("LB", "CM", "SF")
AMYGDALA_ORDER = ("LB_L", "LB_R", "CM_L", "CM_R", "SF_L", "SF_R")
MIN_AMYGDALA_VOXELS = 5

_NUCLEUS_PATTERNS = {
    "LB": re.compile(r"\blb\b|latero ?basal"),
    "CM": re.compile(r"\bcm\b|centro ?medi"),
    "SF": re.compile(r"\bsf\b|superficial"),
}


def nucleus_of(name: str) -> str | None:
    """Return LB, CM, or SF when the atlas label is one of those amygdala groups."""
    text = name.lower().replace("_", " ").replace("-", " ")
    text = " ".join(text.split())
    if "amygdala" not in text:
        return None
    for nucleus, pattern in _NUCLEUS_PATTERNS.items():
        if pattern.search(text):
            return nucleus
    return None


def side_of(name: str) -> str | None:
    """Return L or R when the label name already names a hemisphere."""
    text = name.lower().replace("_", " ").replace("-", " ")
    if re.search(r"\bleft\b|\blh\b", text):
        return "L"
    if re.search(r"\bright\b|\brh\b", text):
        return "R"
    return None


def parse_juelich_xml(path: Path) -> list[tuple[int, str]]:
    """Read FSL atlas XML. The label index is the integer stored in the NIfTI."""
    if not path.is_file():
        raise FileNotFoundError(f"Juelich label file not found: {path}")
    root = ET.parse(path).getroot()
    labels = []
    for node in root.iter("label"):
        if "index" not in node.attrib:
            continue
        name = " ".join((node.text or "").split())
        if not name:
            continue
        labels.append((int(node.attrib["index"]), name))
    if not labels:
        raise ValueError(f"No label entries in {path}")
    return labels


def _fsl_atlas_paths() -> tuple[Path, Path] | None:
    environment = os.environ.get("FSLDIR")
    if not environment:
        return None
    root = Path(environment)
    images = [
        root / "data" / "atlases" / "Juelich" / "Juelich-maxprob-thr25-2mm.nii.gz",
        root / "data" / "atlases" / "Juelich" / "Juelich-maxprob-thr25-2mm.nii",
    ]
    labels = [
        root / "data" / "atlases" / "Juelich.xml",
        root / "data" / "atlases" / "Juelich" / "Juelich-maxprob-thr25-2mm.xml",
    ]
    image = next((path for path in images if path.is_file()), None)
    label = next((path for path in labels if path.is_file()), None)
    if image is not None and label is not None:
        return image, label
    return None


def load_juelich_atlas(atlas_path: Path | None, label_path: Path | None) -> tuple[nib.spatialimages.SpatialImage, list[tuple[int, str]], str]:
    """Load the thr25 2 mm maximum-probability map and its label names."""
    if atlas_path is not None:
        if label_path is None:
            raise ValueError("Pass --juelich-labels with the FSL Juelich XML alongside --juelich-atlas")
        if not atlas_path.is_file():
            raise FileNotFoundError(f"Juelich atlas not found: {atlas_path}")
        return nib.load(str(atlas_path)), parse_juelich_xml(label_path), str(atlas_path)
    found = _fsl_atlas_paths()
    if found is not None:
        image_path, xml_path = found
        print(f"Juelich atlas: {image_path}")
        return nib.load(str(image_path)), parse_juelich_xml(xml_path), str(image_path)
    try:
        from nilearn.datasets import fetch_atlas_juelich
    except ImportError as exc:
        raise FileNotFoundError(
            "The Juelich thr25 2 mm atlas was not found. Set FSLDIR, or pass "
            "--juelich-atlas and --juelich-labels. nilearn is not installed, so it "
            "cannot be downloaded."
        ) from exc
    print("FSLDIR is unset. Downloading the FSL 5.0 Juelich maxprob-thr25-2mm atlas with nilearn.")
    data = fetch_atlas_juelich("maxprob-thr25-2mm")
    maps = data.maps
    image = maps if hasattr(maps, "shape") else nib.load(str(maps))
    labels = [(index, str(name)) for index, name in enumerate(data.labels)]
    source = "nilearn:maxprob-thr25-2mm"
    filename = getattr(image, "get_filename", lambda: None)()
    if filename:
        source = str(filename)
    return image, labels, source


def _mask_from_indices(data: np.ndarray, indices: list[int]) -> np.ndarray:
    return np.isin(data, np.asarray(indices, dtype=np.int32))


def _split_bilateral(data: np.ndarray, image: nib.spatialimages.SpatialImage, nucleus: str) -> tuple[np.ndarray, np.ndarray]:
    coordinates = np.argwhere(data)
    world = nib.affines.apply_affine(image.affine, coordinates)
    left = np.zeros(data.shape, dtype=bool)
    right = np.zeros(data.shape, dtype=bool)
    left_rows = world[:, 0] < 0
    right_rows = world[:, 0] > 0
    midline = int(np.sum(~left_rows & ~right_rows))
    if midline:
        print(f"{nucleus}: {midline} midline voxels were not assigned to either hemisphere")
    if np.any(left_rows):
        left[tuple(coordinates[left_rows].T)] = True
    if np.any(right_rows):
        right[tuple(coordinates[right_rows].T)] = True
    return left, right


def amygdala_masks(image: nib.spatialimages.SpatialImage,
                   labels: list[tuple[int, str]]) -> dict[str, nib.spatialimages.SpatialImage]:
    """Build six binary masks in the atlas grid.

    A label that already says left or right is used as that hemisphere. A
    bilateral label is split by MNI x. The six masks must all be non-empty
    before resampling.
    """
    data = np.rint(np.asarray(image.get_fdata(dtype=np.float32))).astype(np.int32)
    grouped: dict[str, dict[str, list[tuple[int, str]]]] = {
        nucleus: {"L": [], "R": [], "both": []} for nucleus in NUCLEI
    }
    for index, name in labels:
        nucleus = nucleus_of(name)
        if nucleus is None:
            continue
        side = side_of(name) or "both"
        grouped[nucleus][side].append((index, name))
    present = {int(value) for value in np.unique(data)}
    masks: dict[str, np.ndarray] = {}
    for nucleus in NUCLEI:
        found = grouped[nucleus]
        names = [name for side in found.values() for _index, name in side]
        if not names:
            matched = [name for _index, name in labels if "amygdala" in name.lower()]
            raise ValueError(
                f"No Juelich label matched amygdala {nucleus}. "
                f"Amygdala labels in the file: {matched[:12]}"
            )
        missing_from_image = [
            f"{name} (index {index})"
            for side in found.values()
            for index, name in side
            if index not in present
        ]
        if missing_from_image and not any(index in present for side in found.values() for index, _name in side):
            raise ValueError(
                f"Amygdala {nucleus} labels are absent from the atlas image: {missing_from_image[:6]}"
            )
        left_ids = [index for index, _name in found["L"] if index in present]
        right_ids = [index for index, _name in found["R"] if index in present]
        both_ids = [index for index, _name in found["both"] if index in present]
        if left_ids and right_ids and not both_ids:
            masks[f"{nucleus}_L"] = _mask_from_indices(data, left_ids)
            masks[f"{nucleus}_R"] = _mask_from_indices(data, right_ids)
            continue
        if both_ids and not left_ids and not right_ids:
            region = _mask_from_indices(data, both_ids)
            left, right = _split_bilateral(region, image, nucleus)
            masks[f"{nucleus}_L"] = left
            masks[f"{nucleus}_R"] = right
            continue
        raise ValueError(
            f"Amygdala {nucleus} labels mix bilateral and hemispheric names, so the six ROIs "
            f"are ambiguous: {names}"
        )
    empty = [name for name in AMYGDALA_ORDER if int(masks[name].sum()) == 0]
    if empty:
        raise ValueError(f"Amygdala masks are empty in atlas space: {empty}")
    images = {}
    for name in AMYGDALA_ORDER:
        images[name] = nib.Nifti1Image(masks[name].astype(np.uint8), image.affine, image.header)
        print(f"{name}: {int(masks[name].sum())} voxels in atlas space")
    return images
