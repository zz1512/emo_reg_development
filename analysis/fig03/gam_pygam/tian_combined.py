"""Schaefer cortex plus the 32 Tian subcortical parcels.

Schaefer ids stay as they are. Tian labels 1..32 are shifted by the cortical
count, so 200 cortical parcels become 201..232 and 400 cortical parcels become
401..432. Voxel labels are read the same way as fig06: integer labels only,
nearest-neighbor resampling, and no NIfTI header extension. The built-in names
are the Tian Subcortex S2 3T order. Labels 1..16 are right and 17..32 are left.
lamy is lAMY at 3 and 19; mamy is mAMY at 4 and 20. They stay inside the parcel
family. fig06's TIAN_S2_LABELS collapses those four to AMY at V_3 and V_19.
"""

from __future__ import annotations

import re
from pathlib import Path

import nibabel as nib
import numpy as np
import pandas as pd

from analysis.fig03.gam_supplement.prepare_inputs import _resample_to_atlas

TIAN_LABELS = tuple(range(1, 33))
AMYGDALA_KINDS = ("lamy", "mamy")
HEMISPHERES = ("L", "R")
_LABEL_SUFFIXES = {".xml", ".txt", ".lut", ".csv", ".tsv", ".label"}
# Tian_Subcortex_S2_3T_label.txt. Right hemisphere is 1..16, left is 17..32.
TIAN_S2_NAMES = (
    "aHIP-rh", "pHIP-rh", "lAMY-rh", "mAMY-rh",
    "THA-DP-rh", "THA-VP-rh", "THA-VA-rh", "THA-DA-rh",
    "NAc-shell-rh", "NAc-core-rh", "pGP-rh", "aGP-rh",
    "aPUT-rh", "pPUT-rh", "aCAU-rh", "pCAU-rh",
    "aHIP-lh", "pHIP-lh", "lAMY-lh", "mAMY-lh",
    "THA-DP-lh", "THA-VP-lh", "THA-VA-lh", "THA-DA-lh",
    "NAc-shell-lh", "NAc-core-lh", "pGP-lh", "aGP-lh",
    "aPUT-lh", "pPUT-lh", "aCAU-lh", "pCAU-lh",
)
_BARE_NAME = re.compile(r"^[A-Za-z][A-Za-z0-9_+\-]*$")
_INDEX_NAME = re.compile(
    r"(?:label\s+)?index\s*=\s*['\"](\d+)['\"][^>]*>([^<]+)",
    flags=re.IGNORECASE,
)
_HEMI_LEFT = re.compile(r"(^|[^a-z])(left|lh|l)([^a-z]|$)")
_HEMI_RIGHT = re.compile(r"(^|[^a-z])(right|rh|r)([^a-z]|$)")


def marginal_parcels_for(n_cortical: int, cortical_ids: set[int]) -> list[int]:
    """156 and 157 are Schaefer 400 ids. The 200-parcel run has no marginal pair."""
    if int(n_cortical) != 400:
        return []
    marginal = [156, 157]
    missing = [parcel for parcel in marginal if parcel not in cortical_ids]
    if missing:
        raise ValueError(
            "Marginal parcels are missing from the Schaefer dlabel: "
            + ", ".join(str(parcel) for parcel in missing)
        )
    return marginal


def _atlas_stem(path: Path) -> str:
    name = path.name
    if name.endswith(".nii.gz"):
        return name[:-7]
    if name.endswith(".nii"):
        return name[:-4]
    return path.stem


def _clean_name(parts: list[str]) -> str:
    name_parts = list(parts)
    while len(name_parts) > 1 and name_parts[-1].isdigit():
        name_parts.pop()
    return " ".join(name_parts)


def _bare_name_table(text: str) -> dict[int, str]:
    """Official Tian label files are 32 names, one per line, in label order."""
    names = []
    for raw_line in text.splitlines():
        line = raw_line.strip().strip("\x00")
        if not line or line.startswith("#"):
            continue
        if not _BARE_NAME.fullmatch(line):
            return {}
        names.append(line)
    if len(names) != len(TIAN_LABELS):
        return {}
    return {label: name for label, name in zip(TIAN_LABELS, names)}


def _parse_label_text(text: str) -> dict[int, str]:
    """Read integer labels from XML, CSV, 'index name' lines, or 32 bare names."""
    labels: dict[int, str] = {}
    for match in _INDEX_NAME.finditer(text):
        labels[int(match.group(1))] = " ".join(match.group(2).split())
    if set(TIAN_LABELS) <= set(labels):
        return labels
    for raw_line in text.splitlines():
        line = raw_line.strip()
        if not line or line.startswith("#"):
            continue
        lowered = line.lower().replace(" ", "")
        if lowered in {"index,name", "index\tname", "label,name", "id,name", "roi,name"}:
            continue
        parts = re.split(r"[\s,;\t]+", line)
        if len(parts) < 2:
            continue
        if parts[0].isdigit():
            labels[int(parts[0])] = _clean_name(parts[1:])
        elif parts[-1].isdigit() and not parts[0].isdigit():
            labels[int(parts[-1])] = _clean_name(parts[:-1])
    if set(TIAN_LABELS) <= set(labels):
        return labels
    bare = _bare_name_table(text)
    return bare or labels


def _sidecar_paths(atlas_path: Path) -> list[Path]:
    """Label tables for this Tian S2 file only.

    masks_atlas also holds other atlases whose names contain 'label'. Those
    files are not Tian S2 names and must not replace the built-in order.
    """
    stem = _atlas_stem(atlas_path).lower()
    found = []
    if not atlas_path.parent.is_dir():
        return found
    for path in atlas_path.parent.iterdir():
        if not path.is_file() or path.suffix.lower() not in _LABEL_SUFFIXES:
            continue
        if path.stat().st_size > 2_000_000:
            continue
        name = path.name.lower()
        if stem in name or "tian_subcortex_s2_3t" in name:
            found.append(path)
    return sorted(found, key=lambda path: (0 if stem in path.name.lower() else 1, path.name.lower()))


def builtin_tian_names() -> dict[int, str]:
    """Tian S2 3T names. Label 3 is lAMY-rh and label 4 is mAMY-rh."""
    return {label: name for label, name in zip(TIAN_LABELS, TIAN_S2_NAMES)}


def read_tian_label_names(atlas_path: Path) -> dict[int, str]:
    """Names for Tian labels 1..32.

    fig06 loads this atlas as integer voxels and does not read header
    extensions. nibabel's get_content() raises NotImplementedError on the
    extension stored in Tian_Subcortex_S2_3T_2009cAsym.nii.gz, so this function
    does not open extensions either. A complete sidecar replaces the built-in
    Melbourne order; otherwise labels 1..16 are right and 17..32 are left.
    """
    for path in _sidecar_paths(atlas_path):
        table = _parse_label_text(path.read_text(encoding="utf-8", errors="ignore"))
        named = {int(label): str(table[label]) for label in TIAN_LABELS if int(label) in set(map(int, table))}
        if set(named) != set(TIAN_LABELS):
            continue
        if not _has_four_amygdala(named):
            print(f"Skipping Tian label file without lAMY and mAMY on both sides: {path}")
            continue
        print(f"Tian label names: {path}")
        return named
    print(
        "Tian label names: built-in Tian Subcortex S2 3T order. "
        "Labels 1-16 are right, 17-32 are left. "
        "Amygdala labels are 3 lAMY-rh, 4 mAMY-rh, 19 lAMY-lh, 20 mAMY-lh."
    )
    return builtin_tian_names()


def name_kind(name: str) -> str | None:
    """Return lamy or mamy when the label is one of those amygdala parcels."""
    compact = re.sub(r"[^a-z]", "", str(name).lower())
    if "lamy" in compact or "laterobasal" in compact:
        return "lamy"
    if "mamy" in compact:
        return "mamy"
    return None


def _has_four_amygdala(names: dict[int, str]) -> bool:
    found = set()
    for label, name in names.items():
        kind = name_kind(name)
        if kind is None:
            continue
        found.add((kind, name_hemisphere(name, int(label))))
    return found == {(kind, hemi) for kind in AMYGDALA_KINDS for hemi in HEMISPHERES}


def name_hemisphere(name: str, label: int) -> str:
    """L or R from rh/lh in the name. Without a side token, labels 1..16 are left."""
    text = str(name).lower().replace("_", "-")
    left = _HEMI_LEFT.search(text) is not None
    right = _HEMI_RIGHT.search(text) is not None
    if left and not right:
        return "L"
    if right and not left:
        return "R"
    return "L" if int(label) <= 16 else "R"


def amygdala_from_names(names: dict[int, str], n_cortical: int) -> pd.DataFrame:
    """The four lamy/mamy parcels, with ids shifted onto the combined atlas."""
    rows = []
    for label, name in sorted(names.items()):
        kind = name_kind(name)
        if kind is None:
            continue
        rows.append({
            "tian_label": int(label),
            "parcel": int(n_cortical) + int(label),
            "structure": kind,
            "hemi": name_hemisphere(name, int(label)),
            "source_name": str(name),
        })
    table = pd.DataFrame(rows)
    if table.empty:
        preview = ", ".join(f"{label}:{name}" for label, name in sorted(names.items())[:12])
        raise ValueError(f"Tian label names do not include lamy or mamy. First labels: {preview}")
    duplicate = table.duplicated(["structure", "hemi"])
    if duplicate.any():
        collisions = table.loc[duplicate, "source_name"].tolist()
        raise ValueError(f"More than one Tian label matches an amygdala side: {collisions}")
    expected = {(kind, hemi) for kind in AMYGDALA_KINDS for hemi in HEMISPHERES}
    found = set(zip(table["structure"], table["hemi"]))
    missing = sorted(expected - found)
    if missing:
        found_text = ", ".join(f"{row.structure}_{row.hemi}={row.source_name}" for row in table.itertuples())
        missing_text = ", ".join(f"{kind}_{hemi}" for kind, hemi in missing)
        raise ValueError(
            f"Tian amygdala labels are {found_text}. Missing {missing_text}. "
            "lamy and mamy each need a left and a right parcel."
        )
    order = {(kind, hemi): index for index, (kind, hemi) in enumerate(
        (kind, hemi) for kind in AMYGDALA_KINDS for hemi in HEMISPHERES
    )}
    table["_order"] = [order[(row.structure, row.hemi)] for row in table.itertuples()]
    return table.sort_values("_order").drop(columns="_order").reset_index(drop=True)


def _positive_labels(data: np.ndarray) -> set[int]:
    return {int(value) for value in np.unique(data) if int(value) > 0}


def build_combined_atlas(schaefer_path: Path, tian_path: Path, n_cortical: int,
                         output_path: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Write Schaefer labels plus Tian labels shifted by n_cortical.

    Cortical labels win where the resampled Tian parcel lands on cortex.
    Returns the Tian rows for parcel_info and the four amygdala rows.
    """
    if int(n_cortical) <= 0:
        raise ValueError(f"Cortical parcel count must be positive, got {n_cortical}")
    schaefer_image = nib.load(str(schaefer_path))
    tian_image = nib.load(str(tian_path))
    names = read_tian_label_names(tian_path)
    tian_on_cortex = _resample_to_atlas(
        tian_image, schaefer_image, "nearest", tian_path, 0.0, "Schaefer grid")
    schaefer = np.rint(np.asarray(schaefer_image.get_fdata())).astype(np.int32)
    tian = np.rint(np.asarray(tian_on_cortex.get_fdata())).astype(np.int32)
    found = _positive_labels(tian)
    if found != set(TIAN_LABELS):
        missing = sorted(set(TIAN_LABELS) - found)
        extra = sorted(found - set(TIAN_LABELS))
        raise ValueError(
            "Tian atlas must contain labels 1..32 after resampling onto the Schaefer grid. "
            f"Missing {missing[:8]}. Unexpected {extra[:8]}."
        )
    unknown_names = sorted(set(TIAN_LABELS) - set(int(label) for label in names))
    if unknown_names:
        raise ValueError(
            "Tian label table has no name for "
            + ", ".join(str(label) for label in unknown_names[:8])
        )
    covered = int(np.sum((schaefer > 0) & (tian > 0)))
    combined = np.zeros(schaefer.shape, dtype=np.int16)
    combined[schaefer > 0] = schaefer[schaefer > 0]
    subcortical = (schaefer == 0) & (tian > 0)
    combined[subcortical] = tian[subcortical] + int(n_cortical)
    surviving = {int(value) - int(n_cortical) for value in np.unique(combined) if int(value) > int(n_cortical)}
    lost = sorted(set(TIAN_LABELS) - surviving)
    if lost:
        raise ValueError(
            "Tian labels have no voxels left after cortex is given priority: "
            + ", ".join(str(label) for label in lost)
        )
    output_path.parent.mkdir(parents=True, exist_ok=True)
    nib.save(nib.Nifti1Image(combined, schaefer_image.affine, schaefer_image.header), str(output_path))
    print(
        f"Combined atlas: {int(np.sum(schaefer > 0))} cortical voxels, "
        f"{int(subcortical.sum())} Tian voxels, {covered} Tian voxels covered by cortex. "
        f"Wrote {output_path.name}."
    )
    amygdala = amygdala_from_names(names, int(n_cortical))
    rows = []
    for label in TIAN_LABELS:
        source = str(names[label])
        matched = amygdala[amygdala["tian_label"] == label]
        if matched.empty:
            hemi = name_hemisphere(source, label)
            structure = source
            parcel_name = source
            is_amygdala = False
        else:
            row = matched.iloc[0]
            hemi = str(row["hemi"])
            structure = str(row["structure"])
            parcel_name = f"{structure}_{hemi}"
            is_amygdala = True
        rows.append({
            "parcel": int(n_cortical) + int(label),
            "name": parcel_name,
            "hemi": hemi,
            "network": "",
            "structure": structure,
            "region_type": "subcortical",
            "is_amygdala": bool(is_amygdala),
            "source_name": source,
        })
    return pd.DataFrame(rows), amygdala
