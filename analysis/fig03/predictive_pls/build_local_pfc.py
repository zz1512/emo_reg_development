"""Build the Donahue 2018 liberal PFC membership on this machine.

Liberal PFC is every HCP-MMP area named in Donahue et al. (2018) Table S2.
That table is the published areal list; red rows are the conservative subset
and blue rows are the liberal-only additions. The vertex labels are the
public group-average HCP_MMP1.0 map on fsLR 32k, not the BALSA scene file,
which requires a login. Area uses the TemplateFlow fsLR 32k midthickness.
"""

from __future__ import annotations

import shutil
from pathlib import Path

import nibabel as nib
import numpy as np
import pandas as pd

from analysis.fig03.gam_pygam.tian_combined import TIAN_S2_NAMES
from analysis.fig03.predictive_pls.make_pfc_definition import build_definition

PROJECT = Path(__file__).resolve().parents[3]
SOURCE = PROJECT / "resources" / "atlases" / "donahue2018"
SCHAEFER_INFO = PROJECT / "resources" / "atlases" / "schaefer2018" / "Schaefer2018_400Parcels_7Networks_order_info.txt"
SCHAEFER_DLABEL = PROJECT / "resources" / "atlases" / "schaefer2018" / "Schaefer2018_400Parcels_7Networks_order.dlabel.nii"
OUTPUT = PROJECT / "outputs" / "fig03" / "predictive_pls_432_FD05" / "pfc_definition"

# Table S2, Donahue et al. 2018 SI. Liberal PFC is the whole table.
LIBERAL_AREAS = {
    "i6-8", "s6-8", "8Ad", "8Av", "8BL", "8BM", "8C",
    "9-46d", "a9-46v", "p9-46v", "9a", "9m", "9p",
    "10d", "10pp", "10r", "10v", "a10p", "p10p",
    "11l", "13l", "a24", "p24", "25", "a32pr", "d32", "p32", "s32",
    "44", "45", "46", "47l", "47m", "47s", "a47r", "p47r",
    "55b", "IFJa", "IFJp", "IFSa", "IFSp", "OFC", "pOFC", "SFL",
}
EXCLUDED = {"FEF", "PEF"}


def schaefer_names(path: Path) -> dict[int, str]:
    lines = [line.strip() for line in path.read_text(encoding="utf-8").splitlines() if line.strip()]
    names = {}
    for name, coded in zip(lines[0::2], lines[1::2]):
        names[int(coded.split()[0])] = name
    if set(names) != set(range(1, 401)):
        raise ValueError("Schaefer info file does not name parcels 1..400")
    return names


def parcel_info(path: Path) -> None:
    names = schaefer_names(SCHAEFER_INFO)
    rows = [{"parcel": parcel, "name": names[parcel]} for parcel in range(1, 401)]
    rows += [{"parcel": 400 + index, "name": name} for index, name in enumerate(TIAN_S2_NAMES, start=1)]
    pd.DataFrame(rows).to_csv(path, index=False)


def short_name(label: str) -> str:
    if label.endswith("_ROI") and label[1:2] == "_":
        return label[2:-4].replace("–", "-")
    return label.replace("–", "-")


def hemisphere_masks(glasser_path: Path) -> tuple[np.ndarray, np.ndarray]:
    image = nib.load(str(glasser_path))
    values = np.asarray(image.get_fdata()).reshape(-1)
    labels = image.header.get_axis(0).label[0]
    present = {short_name(text) for text, _ in labels.values()}
    missing = sorted(LIBERAL_AREAS - present)
    if missing or EXCLUDED & LIBERAL_AREAS:
        raise ValueError(f"Liberal area names missing from the Glasser label table: {missing}")
    inside_keys = {key for key, (text, _) in labels.items() if short_name(text) in LIBERAL_AREAS}
    structures = {
        "CIFTI_STRUCTURE_CORTEX_LEFT": np.zeros(32492, dtype=np.float32),
        "CIFTI_STRUCTURE_CORTEX_RIGHT": np.zeros(32492, dtype=np.float32),
    }
    brain_axis = image.header.get_axis(1)
    for name, section, axis in brain_axis.iter_structures():
        structures[name][np.asarray(axis.vertex)] = np.isin(values[section], list(inside_keys)).astype(np.float32)
    return structures["CIFTI_STRUCTURE_CORTEX_LEFT"], structures["CIFTI_STRUCTURE_CORTEX_RIGHT"]


def write_metric(path: Path, values: np.ndarray) -> None:
    array = nib.gifti.GiftiDataArray(np.asarray(values, dtype=np.float32), intent="NIFTI_INTENT_SHAPE")
    nib.save(nib.gifti.GiftiImage(darrays=[array]), str(path))


def main():
    source_names = {
        "glasser.dlabel.nii": SOURCE / "Q1-Q6_RelatedValidation210.CorticalAreas_dil_Final_Final_Areas_Group_Colors.32k_fs_LR.dlabel.nii",
        "L.midthickness.surf.gii": SOURCE / "tpl-fsLR_den-32k_hemi-L_midthickness.surf.gii",
        "R.midthickness.surf.gii": SOURCE / "tpl-fsLR_den-32k_hemi-R_midthickness.surf.gii",
    }
    staged = Path(r"C:\Users\dell\AppData\Local\Temp\pfc_src")
    SOURCE.mkdir(parents=True, exist_ok=True)
    for staged_name, destination in source_names.items():
        if not destination.is_file():
            shutil.copy2(staged / staged_name, destination)
    left, right = hemisphere_masks(source_names["glasser.dlabel.nii"])
    mask_left = SOURCE / "donahue2018_liberal_pfc_lh.func.gii"
    mask_right = SOURCE / "donahue2018_liberal_pfc_rh.func.gii"
    write_metric(mask_left, left)
    write_metric(mask_right, right)
    info = SOURCE / "parcel_info_432.csv"
    parcel_info(info)
    table = build_definition(
        SCHAEFER_DLABEL, info, mask_left, mask_right,
        source_names["L.midthickness.surf.gii"], source_names["R.midthickness.surf.gii"], OUTPUT)
    print(f"PFC parcels: {int(table['is_pfc'].sum())}")
    print(f"Boundary parcels: {int(table['boundary'].sum())}")
    print(f"Review-band parcels: {int(table['review_band'].sum())}")


if __name__ == "__main__":
    main()
