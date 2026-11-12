"""Map the Buhle et al. 2014 reappraisal meta-analysis onto the 432 parcels.

The source map is the published consensus z image. Positive clusters at
z >= 2.3 recover the paper's Table 2 peaks. Voxels are then limited to
Harvard-Oxford prefrontal labels, so parietal, temporal, and insular clusters
in the same map are not called PFC. A Schaefer parcel is included when at
least 25% of its 2 mm voxels fall in that intersection. Tian parcels stay
non-PFC.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import nibabel as nib
import numpy as np
import pandas as pd
from nilearn.datasets import fetch_atlas_harvard_oxford
from nilearn.image import resample_to_img

from analysis.fig03.gam_pygam.tian_combined import TIAN_S2_NAMES

PROJECT = Path(__file__).resolve().parents[3]
SOURCE = PROJECT / "resources" / "atlases" / "buhle2014" / "Buhle_Silvers_2014_Emotion_Regulation_Meta.nii.gz"
SCHAEFER = PROJECT / "resources" / "atlases" / "schaefer2018" / "Schaefer2018_400Parcels_7Networks_order_FSLMNI152_2mm.nii.gz"
INFO = PROJECT / "resources" / "atlases" / "schaefer2018" / "Schaefer2018_400Parcels_7Networks_order_info.txt"
OUTPUT = PROJECT / "outputs" / "fig03" / "predictive_pls_432_FD05_meta_pfc" / "pfc_definition"

SOURCE_NAME = "Buhle2014_reappraisal_PFC"
RULE = "schaefer2mm_voxel_fraction>=0.25_inside_buhle_z2.3_and_harvard_oxford_pfc"
Z_THRESHOLD = 2.3
COVERAGE_THRESHOLD = 0.25
# Harvard-Oxford cort-maxprob-thr25-2mm indices: frontal pole, SFG, MFG,
# IFG triangularis, IFG opercularis, paracingulate, anterior cingulate.
PFC_LABELS = {1, 3, 4, 5, 6, 28, 29}


def schaefer_names(path: Path) -> dict[int, str]:
    lines = [line.strip() for line in path.read_text(encoding="utf-8").splitlines() if line.strip()]
    names = {}
    for name, coded in zip(lines[0::2], lines[1::2]):
        names[int(coded.split()[0])] = name
    if set(names) != set(range(1, 401)):
        raise ValueError("Schaefer info file does not name parcels 1..400")
    return names


def pfc_mask(schaefer_image) -> np.ndarray:
    z_image = nib.load(str(SOURCE))
    meta = nib.Nifti1Image((np.asarray(z_image.dataobj) >= Z_THRESHOLD).astype(np.float32), z_image.affine)
    meta_on_schaefer = np.asarray(resample_to_img(meta, schaefer_image, interpolation="nearest").dataobj) > 0
    harvard = fetch_atlas_harvard_oxford("cort-maxprob-thr25-2mm")
    cortical = np.asarray(resample_to_img(harvard.maps, schaefer_image, interpolation="nearest").dataobj)
    return meta_on_schaefer & np.isin(cortical, list(PFC_LABELS))


def build_table() -> pd.DataFrame:
    if not SOURCE.is_file():
        raise FileNotFoundError(f"Missing Buhle meta-analysis map: {SOURCE}")
    schaefer_image = nib.load(str(SCHAEFER))
    labels = np.asarray(schaefer_image.dataobj).astype(int)
    mask = pfc_mask(schaefer_image)
    names = schaefer_names(INFO)
    rows = []
    for parcel, name in names.items():
        voxels = labels == parcel
        fraction = float(mask[voxels].mean()) if voxels.any() else 0.0
        hemisphere = "LH" if "_LH_" in name else "RH" if "_RH_" in name else ""
        rows.append({"parcel": parcel, "name": name, "hemisphere": hemisphere,
                     "coverage_fraction": fraction, "is_pfc": fraction >= COVERAGE_THRESHOLD,
                     "boundary": 0 < fraction < 1, "review_band": 0.15 <= fraction < COVERAGE_THRESHOLD,
                     "definition_source": SOURCE_NAME, "rule": RULE})
    for index, name in enumerate(TIAN_S2_NAMES, start=1):
        rows.append({"parcel": 400 + index, "name": name, "hemisphere": "subcortical",
                     "coverage_fraction": 0.0, "is_pfc": False, "boundary": False, "review_band": False,
                     "definition_source": SOURCE_NAME, "rule": RULE})
    table = pd.DataFrame(rows).sort_values("parcel").reset_index(drop=True)
    if set(table["parcel"]) != set(range(1, 433)) or not table["is_pfc"].any():
        raise ValueError("Meta PFC membership does not partition parcels 1..432")
    return table


def main():
    table = build_table()
    OUTPUT.mkdir(parents=True, exist_ok=True)
    table.to_csv(OUTPUT / "pfc_parcel_membership.csv", index=False)
    included = table.loc[table["is_pfc"], ["parcel", "name", "coverage_fraction"]]
    note = {
        "definition_source": SOURCE_NAME,
        "rule": RULE,
        "citation": "Buhle, Silvers, Wager et al., 2014, Cerebral Cortex 24:2981-2990",
        "map": "CANlab public consensus z map for that meta-analysis",
        "z_threshold": Z_THRESHOLD,
        "anatomical_gate": "Harvard-Oxford frontal pole, SFG, MFG, IFG, paracingulate, anterior cingulate",
        "excluded_from_same_map": "parietal, lateral temporal, insula, and the vmPFC/orbital cortex that the paper did not treat as the reappraisal implementation network",
        "n_pfc": int(table["is_pfc"].sum()),
        "parcels": included.to_dict(orient="records"),
        "sha256": hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
    }
    (OUTPUT / "pfc_definition.json").write_text(json.dumps(note, indent=2, ensure_ascii=False), encoding="utf-8")
    print(f"PFC parcels: {int(table['is_pfc'].sum())}")
    print(included.to_string(index=False))
    print(f"Membership: {OUTPUT / 'pfc_parcel_membership.csv'}")


if __name__ == "__main__":
    main()
