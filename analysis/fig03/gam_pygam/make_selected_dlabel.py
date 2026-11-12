"""Make six colored labels using the standard native Schaefer surface atlas.

The MNI volume is retained only for a separate volume-space reference mask.
All surface outputs use the original dlabel's exact vertex assignments.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path

import nibabel as nib
import numpy as np
from nibabel.cifti2.cifti2_axes import LabelAxis

from analysis.fig03.gam_pygam.brain_maps import MARGINAL_BORDER_RGB, _native_labels_for_volume, _surface_provenance

PROJECT_ROOT = Path(__file__).resolve().parents[3]
ATLAS_DIR = PROJECT_ROOT / "resources/atlases/schaefer2018"
STEM = "Schaefer400_7Networks_selected6_nativeSurface_fslr32k"
SELECTED = {
    152: ("7Networks_LH_Default_Temp_4", (255, 0, 0)),
    154: ("7Networks_LH_Default_Temp_6", (255, 0, 0)),
    156: ("7Networks_LH_Default_Temp_8", tuple(round(c * 255) for c in MARGINAL_BORDER_RGB)),
    157: ("7Networks_LH_Default_Temp_9", tuple(round(c * 255) for c in MARGINAL_BORDER_RGB)),
    182: ("7Networks_LH_Default_PFC_17", (0, 255, 0)),
    185: ("7Networks_LH_Default_PFC_20", (0, 255, 0)),
}


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def make_selected_dlabel(atlas_path: Path, template_path: Path, output_dir: Path) -> Path:
    atlas_path, template_path, output_dir = (path.resolve() for path in (atlas_path, template_path, output_dir))
    source = nib.load(str(atlas_path))
    volume = np.asarray(source.dataobj)
    if volume.ndim != 3 or not np.isfinite(volume).all() or not np.array_equal(volume, np.rint(volume)):
        raise ValueError("Source must be a finite, integer-valued 3D label volume")
    template, labels = _native_labels_for_volume(template_path, atlas_path)
    label_axis, brain_axis = (template.header.get_axis(index) for index in (0, 1))
    structures = ("CIFTI_STRUCTURE_CORTEX_LEFT", "CIFTI_STRUCTURE_CORTEX_RIGHT")
    for key, (expected, _) in SELECTED.items():
        if label_axis.label[0].get(key, (None,))[0] != expected or not np.any(volume == key):
            raise ValueError(f"Parcel {key} does not match the expected Schaefer400 label {expected}")

    vertex_labels = dict(zip(structures, labels))
    for key in SELECTED:
        if not np.any(vertex_labels[structures[0]] == key) or np.any(vertex_labels[structures[1]] == key):
            raise ValueError(f"Parcel {key} must be present in the left hemisphere and absent from the right")

    table = {0: ("???", (0., 0., 0., 0.))}
    table.update({key: (name, (*[channel / 255. for channel in rgb], 1.))
                  for key, (name, rgb) in SELECTED.items()})
    selected_keys = list(SELECTED)
    data = np.zeros((1, len(brain_axis)), dtype=np.int32)
    for name, indices, axis in brain_axis.iter_structures():
        values = vertex_labels[name][axis.vertex]
        data[0, indices] = np.where(np.isin(values, selected_keys), values, 0)
    if set(np.unique(data)) != {0, *selected_keys}:
        raise ValueError("The target brain-model axis lost one of the selected parcels")
    image = nib.Cifti2Image(data, (LabelAxis(["Selected six left-hemisphere ROIs"], [table]), brain_axis))
    image.nifti_header.set_intent(3007)  # NIFTI_INTENT_CONNECTIVITY_DENSE_LABEL
    image.header.matrix.metadata = nib.cifti2.Cifti2MetaData({
        "SourceAtlas": str(template_path), "SurfaceAssignment": "Native Schaefer labels; no projection or smoothing",
        "SelectedParcels": ",".join(map(str, selected_keys)),
        "DisplayInstructions": "Drawing Type = Outline Label Color; zero is transparent; labels retain independent IDs",
    })
    output_dir.mkdir(parents=True, exist_ok=True)
    output_path = output_dir / f"{STEM}.dlabel.nii"
    nib.save(image, str(output_path))

    # Separate GIFTI files make the explicit hemisphere metadata inspectable.
    for name, hemi, structure in zip(structures, ("L", "R"), ("CortexLeft", "CortexRight")):
        values = vertex_labels[name]
        selected = np.where(np.isin(values, selected_keys), values, 0).astype(np.int32)
        gifti = nib.GiftiImage(darrays=[nib.gifti.GiftiDataArray(selected, intent="NIFTI_INTENT_LABEL")])
        gifti.meta["AnatomicalStructurePrimary"] = structure
        gifti.darrays[0].meta["Name"] = "Selected six ROIs"
        for key, (label, rgba) in table.items():
            entry = nib.gifti.GiftiLabel(key, *rgba)
            entry.label = label
            gifti.labeltable.labels.append(entry)
        nib.save(gifti, str(output_dir / f"{STEM}.{hemi}.label.gii"))

    selected_volume = np.where(np.isin(volume, selected_keys), volume, 0).astype(np.int16)
    header = source.header.copy()
    header.set_data_dtype(np.int16)
    nib.save(nib.Nifti1Image(selected_volume, source.affine, header),
             str(output_dir / "Schaefer400_7Networks_selected6_source_volume.nii.gz"))
    records = [{"parcel": key, "name": name, "color_hex": "#" + "".join(f"{c:02X}" for c in rgb),
                "red": rgb[0], "green": rgb[1], "blue": rgb[2],
                "source_voxels": int((volume == key).sum()),
                "left_vertices": int((vertex_labels[structures[0]] == key).sum()),
                "right_vertices": int((vertex_labels[structures[1]] == key).sum())}
               for key, (name, rgb) in SELECTED.items()]
    with (output_dir / "selected_rois.tsv").open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(records[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(records)
    manifest = {
        **_surface_provenance(template_path),
        "reference_volume_only": str(atlas_path), "reference_volume_sha256": _sha256(atlas_path),
        "projection": "none; select the original surface parcel IDs directly",
        "output_dlabel": output_path.name, "nifti_intent": 3007,
        "background": "key 0, transparent RGBA 0/0/0/0",
        "geometry": "exact native surface vertex labels; no smoothing, dilation, erosion or relabeling; distinct atlas keys retained",
        "display": "Outline Label Color uses embedded per-label colors. Drawing mode and width are viewer/scene properties, not dlabel label-table properties.",
        "rois": records,
    }
    (output_dir / "manifest.json").write_text(json.dumps(manifest, indent=2, ensure_ascii=False), encoding="utf-8")
    print(output_path)
    for record in records:
        print(f"{record['parcel']}: {record['name']} {record['color_hex']} L={record['left_vertices']} R={record['right_vertices']}")
    return output_path


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--atlas", type=Path, default=ATLAS_DIR / "Schaefer2018_400Parcels_7Networks_order_FSLMNI152_2mm.nii.gz")
    parser.add_argument("--template", type=Path, default=ATLAS_DIR / "Schaefer2018_400Parcels_7Networks_order.dlabel.nii")
    parser.add_argument("--output", type=Path, default=PROJECT_ROOT / "outputs/fig03/selected_roi_dlabel_native")
    args = parser.parse_args()
    make_selected_dlabel(args.atlas, args.template, args.output)


if __name__ == "__main__":
    main()
