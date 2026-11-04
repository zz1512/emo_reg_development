"""Export the design, parcel responses, and group means used by the mgcv supplement."""

from __future__ import annotations

import inspect
from pathlib import Path

import nibabel as nib
import numpy as np
import pandas as pd
from nilearn.image import resample_to_img

CONDITIONS = ("rpsl", "lkng", "lknt")

YEO_ORDER = ("Vis", "SomMot", "DorsAttn", "SalVentAttn", "Limbic", "Cont", "Default")
SYSTEM_ORDER = (
    "reappraisal_only",
    "common_appraisal",
    "nonmodifiable_emotion",
    "modifiable_emotion",
)
SYSTEM_MASK_NAMES = {name: f"mask_{name}.nii.gz" for name in SYSTEM_ORDER}
# Paths written in original_code/brain_act_uni/brain_act_er_age.py and
# original_code/er_systems/er_systems_identification.py.
ORIGINAL_DEMOGRAPHICS = Path("/public/home/dingrui/fmri_analysis/data/beh/participants_in_tfmri_demographics.csv")
ORIGINAL_SUBJECTS = Path("/public/home/dingrui/fmri_analysis/data/beh/sub_list_valid_4_ER.txt")
ORIGINAL_EXCLUDED = Path("/public/home/dingrui/fmri_analysis/data/fmri/subs_removed_from_2ndlevel.txt")
ORIGINAL_PARENT_EDU = Path("/public/home/dingrui/fmri_analysis/data/beh/env_family_student_cbcl.tsv")
PARENT_EDU_ID_COLUMN = "被试编号"
PARENT_EDU_VALUE_COLUMN = "父母受教育水平_父母受教育水平"
ORIGINAL_MEANFD = Path("/public/home/dingrui/fmri_analysis/data/fmri/qc_res_subs_runs_fmri.csv")
ORIGINAL_VOLUME_DIR = Path("/public/home/dingrui/fmri_analysis/data/fmri/res_1stlevel")
ORIGINAL_SURFACE_DIR = Path("/public/home/dingrui/fmri_analysis/data/fmri/res_1stlevel_surf")
ORIGINAL_SYSTEM_MASKS = Path("/public/home/dingrui/fmri_analysis/res/res_bayes")
# emo+soc branch of brain_act_er_age.py: one GIFTI per hemisphere.
SURFACE_METRICS = {
    "rpsl": ("reg_neg_Observe_beta_lh.func.gii", "reg_neg_Observe_beta_rh.func.gii"),
    "lkng": ("pas_neg_Observe_beta_lh.func.gii", "pas_neg_Observe_beta_rh.func.gii"),
}
VOLUME_FILENAMES = {
    "rpsl": ("beta_0001.nii.gz", "effectSize-reg_neg_Observe_across_runs.nii.gz"),
    "lkng": ("beta_0002.nii.gz", "effectSize-pas_neg_Observe_across_runs.nii.gz"),
    "lknt": ("beta_0003.nii.gz", "effectSize-pas_neu_Observe_across_runs.nii.gz"),
}
# Subject-level volumes read by er_systems_identification.py.
EFFECT_VOLUME_FILES = {
    "rpsl": "effectSize-reg_neg_Observe_across_runs.nii.gz",
    "lkng": "effectSize-pas_neg_Observe_across_runs.nii.gz",
}
# Toy masks only exercise the file layout. They are not the paper's systems.
TOY_SYSTEM_NETWORKS = {
    "reappraisal_only": ("DorsAttn", "SalVentAttn"),
    "common_appraisal": ("Cont", "Default"),
    "nonmodifiable_emotion": ("Vis", "SomMot"),
    "modifiable_emotion": ("Limbic",),
}
MIN_SYSTEM_FRACTION = 0.20
_RESAMPLE_NOTE = {"printed": False}


def _resample_to_atlas(image: nib.spatialimages.SpatialImage, atlas_img: nib.spatialimages.SpatialImage,
                       interpolation: str, source: Path, fill_value: float,
                       grid_name: str = "Schaefer grid") -> nib.spatialimages.SpatialImage:
    """Return an image on the target grid. Matching geometry is left unchanged."""
    same_grid = image.shape == atlas_img.shape and np.allclose(image.affine, atlas_img.affine)
    if same_grid:
        return image
    if not _RESAMPLE_NOTE["printed"]:
        print(
            f"Resampling maps onto the {grid_name} {tuple(int(v) for v in atlas_img.shape)} "
            f"with {interpolation} interpolation. Example {source.name} has shape {tuple(int(v) for v in image.shape)}."
        )
        _RESAMPLE_NOTE["printed"] = True
    parameters = inspect.signature(resample_to_img).parameters
    kwargs = {"interpolation": interpolation}
    if "force_resample" in parameters:
        kwargs["force_resample"] = True
    if "copy_header" in parameters:
        kwargs["copy_header"] = True
    if "fill_value" in parameters:
        kwargs["fill_value"] = fill_value
    return resample_to_img(image, atlas_img, **kwargs)


def _delta_on_atlas(rpsl_path: Path, lkng_path: Path, atlas_img: nib.spatialimages.SpatialImage) -> np.ndarray:
    """Reappraise minus Look on the Schaefer voxel grid."""
    rpsl = _resample_to_atlas(nib.load(str(rpsl_path)), atlas_img, "continuous", rpsl_path, np.nan)
    lkng = _resample_to_atlas(nib.load(str(lkng_path)), atlas_img, "continuous", lkng_path, np.nan)
    return rpsl.get_fdata(dtype=np.float32) - lkng.get_fdata(dtype=np.float32)


def schaefer_parcel_info(dlabel_path: Path) -> pd.DataFrame:
    """Read Schaefer 200 x 7 parcel names from the CIFTI label table."""
    image = nib.load(str(dlabel_path))
    label_axis = image.header.get_axis(0)
    if not hasattr(label_axis, "label"):
        fallback = Path(__file__).resolve().parents[3] / (
            "resources/atlases/schaefer2018/Schaefer2018_200Parcels_7Networks_order.dlabel.nii"
        )
        if dlabel_path.resolve() != fallback.resolve() and fallback.is_file():
            return schaefer_parcel_info(fallback)
        raise ValueError(f"Schaefer CIFTI has no parcel-name table: {dlabel_path}")
    rows = []
    for key, value in label_axis.label[0].items():
        parcel = int(key)
        if parcel <= 0:
            continue
        name = value[0] if isinstance(value, (tuple, list)) else str(value)
        parts = str(name).split("_")
        if len(parts) < 3:
            raise ValueError(f"Unrecognized Schaefer parcel name: {name}")
        rows.append({"parcel": parcel, "name": str(name), "hemi": parts[1], "network": parts[2]})
    info = pd.DataFrame(rows).sort_values("parcel").reset_index(drop=True)
    unknown = sorted(set(info["network"]) - set(YEO_ORDER))
    if unknown:
        raise ValueError(f"Schaefer names contain networks outside Yeo 7: {unknown}")
    return info


def read_schaefer_order_info(path: Path) -> pd.DataFrame:
    """Read a Schaefer order_info.txt table. Networks are not restricted to Yeo 7."""
    if not path.is_file():
        raise FileNotFoundError(f"Schaefer order info not found: {path}")
    rows = []
    pending_name = ""
    for line_number, raw_line in enumerate(path.read_text(encoding="utf-8").splitlines(), start=1):
        text = raw_line.strip()
        if not text:
            continue
        parts = text.split()
        if parts[0].isdigit():
            if not pending_name:
                raise ValueError(f"{path} line {line_number} has a parcel id before its name")
            name_parts = pending_name.split("_")
            if len(name_parts) < 4:
                raise ValueError(f"Unrecognized Schaefer parcel name: {pending_name}")
            rows.append({
                "parcel": int(parts[0]),
                "name": pending_name,
                "hemi": name_parts[1],
                "network": name_parts[2],
            })
            pending_name = ""
        else:
            pending_name = text
    if pending_name:
        raise ValueError(f"{path} ends on a parcel name with no id line: {pending_name}")
    info = pd.DataFrame(rows)
    if info.empty:
        raise ValueError(f"Schaefer order info has no parcels: {path}")
    if info["parcel"].duplicated().any():
        raise ValueError(f"Schaefer order info has duplicate parcel ids: {path}")
    return info.sort_values("parcel").reset_index(drop=True)


def _parcel_table_for_labels(table: pd.DataFrame, roi_ids: np.ndarray) -> pd.DataFrame:
    required = {"parcel", "name", "hemi", "network"}
    missing_columns = sorted(required - set(table.columns))
    if missing_columns:
        raise ValueError(f"Parcel info table is missing columns: {', '.join(missing_columns)}")
    info = table.copy()
    info["parcel"] = info["parcel"].astype(int)
    if info["parcel"].duplicated().any():
        raise ValueError("Parcel info table has duplicate parcel ids")
    missing = sorted(set(int(parcel) for parcel in roi_ids) - set(info["parcel"].astype(int)))
    if missing:
        preview = ", ".join(str(parcel) for parcel in missing[:8])
        raise ValueError(f"Volume parcel ids are absent from the parcel info table: {preview}")
    info = info[info["parcel"].isin(roi_ids)].copy()
    info["system"] = ""
    return info.sort_values("parcel").reset_index(drop=True)


def _hemisphere_vertex_labels(cifti_path: Path, hemi: str) -> np.ndarray:
    """Schaefer labels on one fsLR hemisphere, in the same order as a metric GIFTI."""
    image = nib.load(str(cifti_path))
    data = np.asarray(image.get_fdata()).reshape(-1)
    structure = "CIFTI_STRUCTURE_CORTEX_LEFT" if hemi == "lh" else "CIFTI_STRUCTURE_CORTEX_RIGHT"
    brain_axis = image.header.get_axis(1)
    if not hasattr(brain_axis, "iter_structures"):
        brain_axis = image.header.get_axis(0)
        data = np.asarray(image.get_fdata()).reshape(-1)
    for name, row_slice, _sub_axis in brain_axis.iter_structures():
        if str(name) == structure:
            return np.rint(data[row_slice]).astype(np.int32)
    raise ValueError(f"{cifti_path} has no {structure} vertices")


def _gifti_vertices(path: Path) -> np.ndarray:
    image = nib.load(str(path))
    return np.asarray(image.darrays[0].data, dtype=np.float32).reshape(-1)


def _surface_subject_dir(surface_dir: Path, subject: str) -> Path | None:
    for name in (str(subject), f"sub-{_subject_key(subject)}"):
        folder = surface_dir / name
        if folder.is_dir():
            return folder
    return None


def surface_metric_paths(surface_dir: Path, metadata: pd.DataFrame) -> tuple[pd.DataFrame, list[dict[str, tuple[Path, Path]]], list[str]]:
    """Keep subjects who have the same four GIFTI files the original GAM reads."""
    kept = []
    found = []
    missing = []
    for row in metadata.itertuples(index=False):
        folder = _surface_subject_dir(surface_dir, row.sub_id)
        paths = None
        if folder is not None:
            paths = {}
            for condition, names in SURFACE_METRICS.items():
                pair = (folder / names[0], folder / names[1])
                if not pair[0].is_file() or not pair[1].is_file():
                    paths = None
                    break
                paths[condition] = pair
        if paths is None:
            missing.append(str(row.sub_id))
        else:
            kept.append(row._asdict() if hasattr(row, "_asdict") else row)
            found.append(paths)
    if not found:
        example = surface_dir / str(metadata["sub_id"].iloc[0]) / SURFACE_METRICS["rpsl"][0]
        raise FileNotFoundError(
            "No subject has the surface maps used by the original GAM. "
            f"Expected files such as {example}"
        )
    if missing:
        preview = ", ".join(missing[:30])
        print(
            f"Left out {len(missing)} subjects with no surface maps; {len(found)} subjects remain. "
            f"Ids: {preview}"
        )
    design = pd.DataFrame(kept)
    if "Index" in design.columns:
        design = design.drop(columns=["Index"])
    return design.reset_index(drop=True), found, missing


def _parcel_vertex_means(labels_by_hemi: dict[str, np.ndarray],
                         deltas_by_hemi: dict[str, np.ndarray],
                         roi_ids: np.ndarray) -> np.ndarray:
    """Mean of finite vertices in each parcel, pooling the two hemispheres."""
    index = np.full(int(roi_ids.max()) + 1, -1, dtype=np.int32)
    index[roi_ids] = np.arange(roi_ids.size)
    sums = np.zeros(roi_ids.size, dtype=np.float64)
    counts = np.zeros(roi_ids.size, dtype=np.float64)
    for hemi, labels in labels_by_hemi.items():
        delta = deltas_by_hemi[hemi]
        if delta.shape != labels.shape:
            raise ValueError(
                f"{hemi} metric has {delta.size} vertices and the Schaefer labels have {labels.size}"
            )
        finite = np.isfinite(delta) & (labels > 0)
        columns = index[labels[finite]]
        valid = columns >= 0
        sums += np.bincount(columns[valid], weights=delta[finite][valid], minlength=roi_ids.size)
        counts += np.bincount(columns[valid], minlength=roi_ids.size)
    if np.any(counts == 0):
        raise ValueError("A Schaefer parcel has no finite surface vertices in this contrast")
    return sums / counts


def parcel_contrast_matrix(metadata: pd.DataFrame, rpsl_paths: list[Path], lkng_paths: list[Path],
                           volume_atlas: Path, system_mask_dir: Path,
                           ) -> tuple[pd.DataFrame, np.ndarray, nib.Nifti1Image, dict[str, np.ndarray], dict[str, np.ndarray]]:
    """Subject-by-parcel means, plus system-mask means from the same resampled contrast."""
    atlas_img = nib.load(str(volume_atlas))
    atlas = np.rint(atlas_img.get_fdata()).astype(np.int32)
    roi_ids = np.unique(atlas[atlas > 0])
    if roi_ids.size == 0:
        raise ValueError(f"Volume atlas contains no positive parcel labels: {volume_atlas}")
    flat_atlas = atlas.ravel()
    present = flat_atlas > 0
    roi_index = np.searchsorted(roi_ids, flat_atlas[present])
    roi_counts = np.bincount(roi_index, minlength=roi_ids.size).astype(float)
    masks = _load_system_masks(system_mask_dir, atlas_img)
    flat_masks = {}
    for system, mask in masks.items():
        flat_mask = mask.ravel()
        if int(flat_mask.sum()) == 0:
            raise ValueError(f"System mask {system} has no voxels on the Schaefer grid")
        flat_masks[system] = flat_mask
    n_subjects = len(metadata)
    values = np.empty((n_subjects, roi_ids.size), dtype=np.float64)
    system_values = {system: np.empty(n_subjects, dtype=np.float64) for system in masks}
    for subject_idx, (rpsl_path, lkng_path) in enumerate(zip(rpsl_paths, lkng_paths)):
        if subject_idx % 100 == 0:
            print(f"Contrast {subject_idx + 1} / {n_subjects}")
        delta = _delta_on_atlas(rpsl_path, lkng_path, atlas_img).ravel()
        weights = delta[present]
        finite = np.isfinite(weights)
        counts = np.bincount(roi_index[finite], minlength=roi_ids.size).astype(float)
        if np.any(counts == 0):
            raise ValueError(f"A Schaefer parcel has no finite contrast voxels: {rpsl_path}")
        sums = np.bincount(roi_index[finite], weights=weights[finite], minlength=roi_ids.size)
        values[subject_idx, :] = sums / counts
        for system, flat_mask in flat_masks.items():
            observed = delta[flat_mask]
            observed = observed[np.isfinite(observed)]
            if observed.size == 0:
                raise ValueError(f"System mask {system} has no finite contrast voxels: {rpsl_path}")
            system_values[system][subject_idx] = float(np.mean(observed))
    if not np.isfinite(values).all():
        raise ValueError("Parcel Reappraise-minus-Look means contain non-finite values")
    frame = pd.DataFrame(values, columns=[str(int(parcel)) for parcel in roi_ids])
    frame.insert(0, "sub_id", metadata["sub_id"].astype(str).to_numpy())
    return frame, roi_counts, atlas_img, masks, system_values


def _load_system_masks(mask_dir: Path, atlas_img: nib.Nifti1Image) -> dict[str, np.ndarray]:
    masks = {}
    for name, filename in SYSTEM_MASK_NAMES.items():
        path = mask_dir / filename
        if not path.is_file():
            raise FileNotFoundError(
                f"Missing system mask {path}. Expected mask_<system>.nii.gz for {', '.join(SYSTEM_ORDER)}."
            )
        image = _resample_to_atlas(nib.load(str(path)), atlas_img, "nearest", path, 0.0)
        masks[name] = np.rint(image.get_fdata()).astype(np.int32) > 0
    return masks


def write_toy_system_masks(atlas_img: nib.Nifti1Image, parcel_info: pd.DataFrame, mask_dir: Path) -> None:
    """Paint disjoint parcel groups into the four system-mask filenames."""
    mask_dir.mkdir(parents=True, exist_ok=True)
    atlas = np.rint(atlas_img.get_fdata()).astype(np.int32)
    network_by_parcel = parcel_info.set_index("parcel")["network"].to_dict()
    for system, networks in TOY_SYSTEM_NETWORKS.items():
        member = np.array([
            parcel for parcel, network in network_by_parcel.items() if network in networks
        ], dtype=np.int32)
        data = np.isin(atlas, member).astype(np.uint8)
        nib.save(nib.Nifti1Image(data, atlas_img.affine, atlas_img.header),
                 str(mask_dir / SYSTEM_MASK_NAMES[system]))


def _load_native_system_masks(mask_dir: Path) -> tuple[dict[str, np.ndarray], nib.Nifti1Image]:
    """Read the four system masks on their own grid.

    These masks were saved from the effect-size volume. Overlap with Schaefer
    is counted after the atlas is moved onto that grid.
    """
    loaded = []
    for name in SYSTEM_ORDER:
        path = mask_dir / SYSTEM_MASK_NAMES[name]
        if not path.is_file():
            raise FileNotFoundError(
                f"Missing system mask {path}. Expected mask_<system>.nii.gz for {', '.join(SYSTEM_ORDER)}."
            )
        loaded.append((name, path, nib.load(str(path))))
    reference = loaded[0][2]
    masks = {}
    for name, path, image in loaded:
        if image.shape != reference.shape or not np.allclose(image.affine, reference.affine):
            image = _resample_to_atlas(image, reference, "nearest", path, 0.0, "system-mask grid")
        mask = np.rint(np.asarray(image.get_fdata())).astype(np.int32) > 0
        print(f"{name}: {int(mask.sum())} voxels")
        if int(mask.sum()) == 0:
            raise ValueError(f"System mask {name} has no positive voxels: {path}")
        masks[name] = mask
    return masks, reference


def parcel_overlap_on_mask_grid(atlas_img: nib.Nifti1Image, atlas_path: Path,
                                masks: dict[str, np.ndarray],
                                reference: nib.Nifti1Image) -> tuple[pd.DataFrame, dict[str, dict[int, int]]]:
    """Parcel overlap measured on the system-mask grid, not the full Schaefer grid."""
    atlas_on_mask = _resample_to_atlas(atlas_img, reference, "nearest", atlas_path, 0.0, "system-mask grid")
    labels = np.rint(np.asarray(atlas_on_mask.get_fdata())).astype(np.int32)
    flat = labels.ravel()
    parcel_sizes = np.bincount(np.clip(flat, 0, None))
    hit_counts = {}
    overlap: dict[str, dict[int, int]] = {}
    for system, mask in masks.items():
        selected = flat[mask.ravel()]
        counts = np.bincount(np.clip(selected, 0, None), minlength=parcel_sizes.size)
        hit_counts[system] = counts
        overlap[system] = {
            int(parcel): int(counts[parcel])
            for parcel in np.flatnonzero(counts)
            if int(parcel) > 0
        }
        fractions = [
            count / parcel_sizes[parcel]
            for parcel, count in overlap[system].items()
            if parcel_sizes[parcel]
        ]
        n_labeled = sum(1 for fraction in fractions if fraction >= MIN_SYSTEM_FRACTION)
        best = max(fractions) if fractions else 0.0
        print(
            f"{system}: {sum(overlap[system].values())} Schaefer voxels, "
            f"{len(overlap[system])} parcels overlap, {n_labeled} parcels reach "
            f"{MIN_SYSTEM_FRACTION:.0%}, max parcel fraction {best:.3f}"
        )
    assignment = []
    parcels = sorted(int(value) for value in np.unique(flat) if value > 0)
    for parcel in parcels:
        fractions = {
            system: float(counts[parcel] / parcel_sizes[parcel]) if parcel < counts.size and parcel_sizes[parcel] else 0.0
            for system, counts in hit_counts.items()
        }
        best = max(SYSTEM_ORDER, key=lambda system: (fractions[system], -SYSTEM_ORDER.index(system)))
        assignment.append({
            "parcel": parcel,
            "system": best if fractions[best] >= MIN_SYSTEM_FRACTION else "",
            "system_fraction": fractions[best],
        })
    return pd.DataFrame(assignment), overlap


def parcel_system_assignment(atlas_img: nib.Nifti1Image, masks: dict[str, np.ndarray]) -> pd.DataFrame:
    """Modal system of each Schaefer parcel from masks already on the atlas grid."""
    atlas = np.rint(atlas_img.get_fdata()).astype(np.int32)
    flat_atlas = atlas.ravel()
    parcel_sizes = np.bincount(flat_atlas)
    hit_counts = {}
    for system, mask in masks.items():
        flat_mask = mask.ravel()
        if flat_mask.shape != flat_atlas.shape:
            raise ValueError(f"System mask {system} does not match the Schaefer volume grid")
        hit_counts[system] = np.bincount(flat_atlas[flat_mask], minlength=parcel_sizes.size)
    assignment = []
    for parcel in sorted(int(value) for value in np.unique(flat_atlas) if value > 0):
        fractions = {
            system: float(counts[parcel] / parcel_sizes[parcel]) if parcel_sizes[parcel] else 0.0
            for system, counts in hit_counts.items()
        }
        best = max(SYSTEM_ORDER, key=lambda system: (fractions[system], -SYSTEM_ORDER.index(system)))
        assignment.append({
            "parcel": parcel,
            "system": best if fractions[best] >= MIN_SYSTEM_FRACTION else "",
            "system_fraction": fractions[best],
        })
    return pd.DataFrame(assignment)


def system_means_from_parcels(parcel_response: pd.DataFrame,
                              overlap: dict[str, dict[int, int]]) -> pd.DataFrame:
    """Mean of surface parcels, weighted by how many system-mask voxels each parcel contains."""
    rows = []
    subject_ids = parcel_response["sub_id"].astype(str)
    available = {int(column) for column in parcel_response.columns if column != "sub_id"}
    for system in SYSTEM_ORDER:
        parcels = [parcel for parcel, count in overlap[system].items() if count > 0 and parcel in available]
        if not parcels:
            raise ValueError(f"System {system} shares no Schaefer parcels on the mask grid")
        columns = [str(parcel) for parcel in parcels]
        parcel_weights = np.array([overlap[system][parcel] for parcel in parcels], dtype=float)
        means = np.average(parcel_response[columns].to_numpy(dtype=float), axis=1, weights=parcel_weights)
        for sub_id, value in zip(subject_ids, means):
            rows.append({"sub_id": sub_id, "dimension": "system", "group": system, "y": float(value)})
    return pd.DataFrame(rows)


def system_responses_and_assignment(metadata: pd.DataFrame, atlas_img: nib.Nifti1Image,
                                    masks: dict[str, np.ndarray],
                                    system_values: dict[str, np.ndarray]) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Record the mask means and the modal system of each parcel."""
    assignment = parcel_system_assignment(atlas_img, masks)
    subject_ids = metadata["sub_id"].astype(str).to_numpy()
    response_rows = []
    for system in masks:
        for sub_id, value in zip(subject_ids, system_values[system]):
            response_rows.append({"sub_id": sub_id, "dimension": "system", "group": system, "y": float(value)})
    return pd.DataFrame(response_rows), assignment


def yeo_group_responses(parcel_response: pd.DataFrame, parcel_info: pd.DataFrame,
                        voxel_counts: np.ndarray, roi_ids: np.ndarray) -> pd.DataFrame:
    """Voxel-weighted mean of parcel contrasts within each Yeo network."""
    weights = {int(parcel): float(count) for parcel, count in zip(roi_ids, voxel_counts)}
    rows = []
    subject_ids = parcel_response["sub_id"].astype(str)
    for network in YEO_ORDER:
        parcels = parcel_info.loc[parcel_info["network"] == network, "parcel"].astype(int).tolist()
        if not parcels:
            raise ValueError(f"Yeo network {network} has no Schaefer parcels")
        columns = [str(parcel) for parcel in parcels]
        parcel_weights = np.array([weights[parcel] for parcel in parcels], dtype=float)
        means = np.average(parcel_response[columns].to_numpy(dtype=float), axis=1, weights=parcel_weights)
        for sub_id, value in zip(subject_ids, means):
            rows.append({"sub_id": sub_id, "dimension": "yeo", "group": network, "y": float(value)})
    return pd.DataFrame(rows)


def _read_table(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep=None, engine="python")


def _read_id_list(path: Path) -> list[str]:
    if not path.is_file():
        raise FileNotFoundError(f"Subject list not found: {path}")
    identifiers = []
    for line in path.read_text(encoding="utf-8").splitlines():
        token = line.strip().split()[0] if line.strip() else ""
        if token and token.lower() not in {"sub_id", "subject_id"}:
            identifiers.append(token)
    if not identifiers:
        raise ValueError(f"Subject list is empty: {path}")
    return identifiers


def _subject_key(value: object) -> str:
    """Match ids with or without the BIDS sub- prefix."""
    text = str(value).strip()
    if text.lower().startswith("sub-"):
        return text[4:]
    return text


def _parent_edu_key(value: object) -> str:
    """Key used only to align 被试编号 with the imaging sub_id."""
    if value is None or (isinstance(value, float) and not np.isfinite(value)):
        return ""
    text = str(value).strip()
    if text.lower() in {"", "nan", "none", "na"}:
        return ""
    if text.endswith(".0") and text[:-2].replace("-", "", 1).isdigit():
        text = text[:-2]
    return _subject_key(text).casefold()


def _read_text_table(path: Path) -> pd.DataFrame:
    last_error: Exception | None = None
    for encoding in ("utf-8-sig", "utf-8", "gb18030"):
        try:
            table = pd.read_csv(path, sep="\t", encoding=encoding, dtype=str)
            break
        except UnicodeDecodeError as exc:
            last_error = exc
    else:
        raise ValueError(f"Could not decode {path}") from last_error
    table.columns = [str(column).replace("\ufeff", "").strip() for column in table.columns]
    return table


OBSOLETE_BEFORE_DESIGN = "design_before_covariate_filters.csv"
ROSTER_COLUMNS = (
    "sub_id", "age", "gender", "site_id", "parent_edu", "meanFD", "in_analysis", "excluded_by",
)


def _load_parent_education(family_path: Path) -> tuple[dict[str, float], list[tuple[str, list[float]]], pd.DataFrame]:
    """Read 父母受教育水平_父母受教育水平 keyed by 被试编号.

    Conflicting duplicate rows are returned and not resolved. Callers that fit
    a model still stop on those rows. The before-filter list records them.
    """
    if not family_path.is_file():
        raise FileNotFoundError(f"Parental-education table not found: {family_path}")
    family = _read_text_table(family_path)
    required = (PARENT_EDU_ID_COLUMN, PARENT_EDU_VALUE_COLUMN)
    missing = [column for column in required if column not in family.columns]
    if missing:
        education_columns = [column for column in family.columns if "教育" in column or "被试" in column]
        raise ValueError(
            f"{family_path.name} is missing {missing}. Columns containing 被试 or 教育: {education_columns}"
        )
    family = family.loc[:, list(required)].copy()
    family["_key"] = family[PARENT_EDU_ID_COLUMN].map(_parent_edu_key)
    family = family[family["_key"] != ""]
    values = pd.to_numeric(family[PARENT_EDU_VALUE_COLUMN], errors="coerce")
    conflicts = []
    resolved: dict[str, float] = {}
    for key, numeric in values.groupby(family["_key"], sort=False):
        finite = pd.unique(numeric[np.isfinite(numeric.to_numpy(dtype=float))])
        if len(finite) > 1:
            conflicts.append((str(key), [float(item) for item in finite]))
        elif len(finite) == 1:
            resolved[str(key)] = float(finite[0])
        else:
            resolved[str(key)] = np.nan
    return resolved, conflicts, family


def _parent_education_status(metadata: pd.DataFrame, family_path: Path) -> tuple[pd.DataFrame, list[tuple[str, list[float]]], str]:
    """Per-subject parent_edu status for the cohort. Does not drop rows.

    A missing or unreadable table is a status, so the before-filter list can
    be written before the model sample stops. Conflicts stay in the returned
    list even when that id is not in this cohort.
    """
    n_rows = len(metadata)
    parent_edu = np.full(n_rows, np.nan)
    reason = np.full(n_rows, "", dtype=object)
    value_text = np.full(n_rows, "", dtype=object)
    if not family_path.is_file():
        reason[:] = "file_missing"
        return _parent_status_frame(parent_edu, reason, value_text), [], "file_missing"
    try:
        resolved, conflicts, _family = _load_parent_education(family_path)
    except ValueError as exc:
        token = "columns_missing" if "is missing" in str(exc) else "unreadable"
        reason[:] = token
        return _parent_status_frame(parent_edu, reason, value_text), [], token
    conflict_map = {key: levels for key, levels in conflicts}
    keys = metadata["sub_id"].map(_parent_edu_key)
    parent_edu = keys.map(resolved).to_numpy(dtype=float)
    is_conflict = keys.isin(set(conflict_map)).to_numpy()
    in_table = keys.isin(set(resolved)).to_numpy()
    has_value = np.isfinite(parent_edu)
    reason[is_conflict] = "conflicting_values"
    reason[~is_conflict & ~in_table] = "not_in_table"
    reason[~is_conflict & in_table & ~has_value] = "missing_value"
    for index, key in enumerate(keys):
        if key in conflict_map:
            value_text[index] = ";".join(f"{value:g}" for value in conflict_map[key])
    return _parent_status_frame(parent_edu, reason, value_text), conflicts, ""


def _parent_status_frame(parent_edu: np.ndarray, reason: np.ndarray, value_text: np.ndarray) -> pd.DataFrame:
    return pd.DataFrame({
        "parent_edu": parent_edu,
        "parent_edu_reason": reason,
        "parent_edu_values": value_text,
    })


def attach_parent_education(metadata: pd.DataFrame, family_path: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Add parent_edu from 父母受教育水平_父母受教育水平, aligned on 被试编号.

    Subjects with no row, a blank value, or disagreeing duplicate rows do not
    stay in the design. The returned drop table records which of those happened.
    A low id-match rate stops the run instead of fitting a mismatched sample.
    """
    resolved, conflicts, family = _load_parent_education(family_path)
    if conflicts:
        preview = "; ".join(f"{key} -> {levels}" for key, levels in conflicts[:8])
        raise ValueError(
            f"{len(conflicts)} subjects have conflicting values in {PARENT_EDU_VALUE_COLUMN}. Examples: {preview}"
        )
    design = metadata.copy()
    design["_key"] = design["sub_id"].map(_parent_edu_key)
    design["parent_edu"] = design["_key"].map(resolved)
    in_table = design["_key"].isin(set(resolved))
    has_value = np.isfinite(design["parent_edu"].to_numpy(dtype=float))
    n_design = len(design)
    n_in_table = int(in_table.sum())
    if n_design == 0 or n_in_table < 0.5 * n_design:
        design_examples = design.loc[~in_table, "sub_id"].astype(str).head(8).tolist()
        family_examples = family[PARENT_EDU_ID_COLUMN].astype(str).head(8).tolist()
        raise ValueError(
            f"Only {n_in_table} of {n_design} imaging subjects match {PARENT_EDU_ID_COLUMN} "
            f"in {family_path.name}. Imaging ids: {design_examples}. "
            f"Family-table ids: {family_examples}."
        )
    reasons = np.full(n_design, "", dtype=object)
    reasons[~in_table.to_numpy()] = "not_in_table"
    reasons[in_table.to_numpy() & ~has_value] = "missing_value"
    dropped = design.loc[reasons != "", ["sub_id"]].copy()
    dropped["reason"] = reasons[reasons != ""]
    kept = design.loc[has_value].drop(columns="_key").reset_index(drop=True)
    print(
        f"Parent education aligned {len(kept)} of {n_design} subjects "
        f"from {family_path.name} column {PARENT_EDU_VALUE_COLUMN}. "
        f"Dropped {int((reasons == 'missing_value').sum())} blank values and "
        f"{int((reasons == 'not_in_table').sum())} ids absent from the table."
    )
    return kept, dropped.reset_index(drop=True)


def _read_meanfd_table(path: Path, quiet: bool = False) -> pd.DataFrame:
    """Read run-level mean_fd and return one averaged value per subject.

    qc_res_subs_runs_fmri.csv has one row per run. The emotion-regulation
    scans are EMO and SOC. A subject's meanFD is the mean of those runs.
    """
    last_error: Exception | None = None
    table = None
    for encoding in ("utf-8-sig", "utf-8", "gb18030"):
        try:
            table = pd.read_csv(path, encoding=encoding, dtype=str)
            break
        except UnicodeDecodeError as exc:
            last_error = exc
    if table is None:
        raise ValueError(f"Could not decode {path}") from last_error
    table.columns = [str(column).replace("\ufeff", "").strip() for column in table.columns]
    id_column = next((name for name in ("sub_id", "subject_id") if name in table.columns), None)
    value_column = next((name for name in ("mean_fd", "mean_FD", "meanFD") if name in table.columns), None)
    if id_column is None or value_column is None:
        raise ValueError(
            f"{path.name} must contain sub_id and mean_fd. Columns: {list(table.columns)}"
        )
    runs = table.loc[:, [id_column, value_column] + (["run_name"] if "run_name" in table.columns else [])].copy()
    if "run_name" in runs.columns:
        task = runs["run_name"].astype(str).str.strip().str.upper()
        runs = runs[task.isin({"EMO", "SOC"})]
    if runs.empty:
        raise ValueError(f"{path.name} has no EMO or SOC rows with mean_fd")
    runs["_key"] = runs[id_column].map(_parent_edu_key)
    runs["mean_fd"] = pd.to_numeric(runs[value_column], errors="coerce")
    runs = runs[runs["_key"] != ""]
    grouped = runs.groupby("_key", sort=False)["mean_fd"]
    resolved = grouped.mean()
    n_runs = grouped.size()
    summary = pd.DataFrame({
        "sub_id": runs.groupby("_key", sort=False)[id_column].first(),
        "mean_FD": resolved,
        "n_runs": n_runs,
    }).reset_index(drop=True)
    if not quiet:
        print(
            f"mean_fd averaged over EMO and SOC: {len(summary)} subjects from {len(runs)} run rows in {path.name}."
        )
    return summary


def attach_mean_fd(metadata: pd.DataFrame, meanfd_path: Path) -> pd.DataFrame:
    """Add meanFD for every remaining subject from qc_res_subs_runs_fmri.csv.

    Run-level mean_fd values for EMO and SOC are averaged within sub_id.
    A subject still in the design who is absent from the file stops the run.
    """
    if not meanfd_path.is_file():
        raise FileNotFoundError(f"meanFD table not found: {meanfd_path}")
    motion = _read_meanfd_table(meanfd_path)
    motion["_key"] = motion["sub_id"].map(_parent_edu_key)
    resolved = {
        key: float(value)
        for key, value in zip(motion["_key"], motion["mean_FD"])
        if np.isfinite(value)
    }
    design = metadata.copy()
    design["_key"] = design["sub_id"].map(_parent_edu_key)
    design["meanFD"] = design["_key"].map(resolved)
    has_value = np.isfinite(design["meanFD"].to_numpy(dtype=float))
    if not bool(has_value.all()):
        missing_ids = design.loc[~has_value, "sub_id"].astype(str).head(8).tolist()
        file_ids = motion["sub_id"].astype(str).head(8).tolist()
        raise ValueError(
            f"{int((~has_value).sum())} of {len(design)} analysis subjects have no mean FD "
            f"in {meanfd_path.name}. The file has {motion['_key'].nunique()} subjects. "
            f"Imaging ids: {missing_ids}. File ids: {file_ids}."
        )
    kept = design.drop(columns="_key").reset_index(drop=True)
    print(
        f"meanFD aligned all {len(kept)} analysis subjects from {meanfd_path.name}."
    )
    return kept


def _meanfd_status(metadata: pd.DataFrame, meanfd_path: Path) -> tuple[np.ndarray, np.ndarray, str]:
    """Per-subject meanFD for the whole cohort. A missing value does not drop the row."""
    n_rows = len(metadata)
    meanfd = np.full(n_rows, np.nan)
    reason = np.full(n_rows, "", dtype=object)
    if not meanfd_path.is_file():
        reason[:] = "file_missing"
        return meanfd, reason, "file_missing"
    try:
        motion = _read_meanfd_table(meanfd_path, quiet=True)
    except ValueError:
        reason[:] = "unreadable"
        return meanfd, reason, "unreadable"
    motion["_key"] = motion["sub_id"].map(_parent_edu_key)
    resolved = {
        key: float(value)
        for key, value in zip(motion["_key"], motion["mean_FD"])
        if np.isfinite(value)
    }
    mapped = metadata["sub_id"].map(_parent_edu_key).map(resolved).to_numpy(dtype=float)
    reason[~np.isfinite(mapped)] = "missing"
    return mapped, reason, ""


def _has_effect_volumes(volume_dir: Path, metadata: pd.DataFrame) -> np.ndarray:
    """True when both Reappraise and Look effectSize files exist. Does not open them."""
    found = np.zeros(len(metadata), dtype=bool)
    if not volume_dir.is_dir():
        return found
    for index, subject in enumerate(metadata["sub_id"].astype(str)):
        folder = _volume_subject_dir(volume_dir, subject)
        if folder is None:
            continue
        rpsl = folder / EFFECT_VOLUME_FILES["rpsl"]
        lkng = folder / EFFECT_VOLUME_FILES["lkng"]
        found[index] = rpsl.is_file() and lkng.is_file()
    return found


def _first_exclusion(parent_reason: np.ndarray, has_volume: np.ndarray, meanfd_reason: np.ndarray) -> np.ndarray:
    """First filter that removes the subject, in the order the model sample is built."""
    excluded = np.full(len(parent_reason), "", dtype=object)
    for index, reason in enumerate(parent_reason):
        if reason:
            excluded[index] = f"parent_edu_{reason}"
        elif not bool(has_volume[index]):
            excluded[index] = "no_effect_volume"
        elif meanfd_reason[index]:
            excluded[index] = f"meanFD_{meanfd_reason[index]}"
    return excluded


def _ids_preview(identifiers: list[str], limit: int = 30) -> str:
    shown = identifiers[:limit]
    text = ", ".join(shown)
    if len(identifiers) > limit:
        return f"{text} ({len(identifiers)} subjects)"
    return text


def build_subject_roster(metadata: pd.DataFrame, parent_edu_path: Path,
                         volume_dir: Path, meanfd_path: Path) -> pd.DataFrame:
    """One design table for every subject after the valid-list and second-level cut.

    Covariates are filled when a single value exists and left blank otherwise.
    in_analysis is 1 only for subjects who enter the model. excluded_by is the
    first reason a subject is left out. A sample-wide stop is written on the
    rows that had no earlier reason, as sample_stop:<reason>. Conflicts, a
    match rate below one half, and a missing meanFD still stop the fit after
    this table is written. gender and site_id are not row filters.
    """
    parent, conflicts, parent_problem = _parent_education_status(metadata, parent_edu_path)
    has_volume = _has_effect_volumes(volume_dir, metadata)
    meanfd, meanfd_reason, meanfd_problem = _meanfd_status(metadata, meanfd_path)
    parent_reason = parent["parent_edu_reason"].to_numpy(dtype=object)
    excluded_by = _first_exclusion(parent_reason, has_volume, meanfd_reason)
    sample_stop, _sample_stop_detail = _sample_stop(
        metadata, parent_reason, has_volume, meanfd_reason, parent_problem, meanfd_problem,
        conflicts, parent_edu_path, volume_dir, meanfd_path,
    )
    if sample_stop:
        excluded_by = excluded_by.copy()
        excluded_by[excluded_by == ""] = f"sample_stop:{sample_stop}"
    in_analysis = (excluded_by == "").astype(int)
    table = pd.DataFrame({
        "sub_id": metadata["sub_id"].astype(str).to_numpy(),
        "age": pd.to_numeric(metadata["age"], errors="coerce").to_numpy(),
        "gender": pd.to_numeric(metadata["gender"], errors="coerce").to_numpy(),
        "site_id": pd.to_numeric(metadata["site_id"], errors="coerce").to_numpy(),
        "parent_edu": parent["parent_edu"].to_numpy(dtype=float),
        "meanFD": meanfd,
        "in_analysis": in_analysis,
        "excluded_by": excluded_by,
    })
    return table.loc[:, list(ROSTER_COLUMNS)]


def _sample_stop(metadata: pd.DataFrame, parent_reason: np.ndarray, has_volume: np.ndarray,
                 meanfd_reason: np.ndarray, parent_problem: str, meanfd_problem: str,
                 conflicts: list[tuple[str, list[float]]], parent_edu_path: Path,
                 volume_dir: Path, meanfd_path: Path) -> tuple[str, str]:
    """Reason the fit stops, matching the order of the checks that raise."""
    if parent_problem == "file_missing":
        return "parent_edu_file_missing", str(parent_edu_path)
    if parent_problem in {"columns_missing", "unreadable"}:
        return f"parent_edu_{parent_problem}", parent_edu_path.name
    if conflicts:
        preview = "; ".join(f"{key} -> {levels}" for key, levels in conflicts[:8])
        return "parent_edu_conflicting_values", preview
    n_design = len(metadata)
    in_table = np.isin(parent_reason, ["", "missing_value"])
    n_in_table = int(in_table.sum())
    if n_design == 0 or n_in_table < 0.5 * n_design:
        return "parent_edu_match_rate", f"{n_in_table} of {n_design} imaging ids match the family table"
    if not volume_dir.is_dir():
        return "volume_dir_missing", str(volume_dir)
    passed_parent = parent_reason == ""
    if passed_parent.any() and not bool(has_volume[passed_parent].any()):
        return "no_effect_volume", "no subject who passed parent education has both effectSize files"
    if meanfd_problem == "file_missing":
        return "meanFD_file_missing", str(meanfd_path)
    if meanfd_problem == "unreadable":
        return "meanFD_unreadable", meanfd_path.name
    still_in = passed_parent & has_volume
    missing_motion = still_in & (meanfd_reason != "")
    if bool(missing_motion.any()):
        missing_ids = metadata.loc[missing_motion, "sub_id"].astype(str).tolist()
        return "meanFD_missing", _ids_preview(missing_ids)
    return "", ""


def analysis_design(table: pd.DataFrame) -> pd.DataFrame:
    """Subjects who enter the model, in design.csv order.

    A table without in_analysis is already the analysis sample.
    """
    columns = [column for column in ("sub_id", "age", "gender", "site_id", "parent_edu", "meanFD") if column in table.columns]
    if "in_analysis" not in table.columns:
        return table.loc[:, columns].reset_index(drop=True)
    entered = pd.to_numeric(table["in_analysis"], errors="coerce").fillna(0).to_numpy() == 1
    return table.loc[entered, columns].reset_index(drop=True)


def n_analysis_subjects(table: pd.DataFrame) -> int:
    """How many design.csv rows enter the model."""
    if "in_analysis" not in table.columns:
        return int(len(table))
    return int((pd.to_numeric(table["in_analysis"], errors="coerce").fillna(0) == 1).sum())


def _remove_obsolete_before_design(input_dir: Path) -> None:
    obsolete = input_dir / OBSOLETE_BEFORE_DESIGN
    if obsolete.is_file():
        obsolete.unlink()


def write_subject_roster(input_dir: Path, metadata: pd.DataFrame, parent_edu_path: Path,
                         volume_dir: Path, meanfd_path: Path) -> pd.DataFrame:
    """Write design.csv for the whole cohort before parent_edu or meanFD can stop the run."""
    input_dir.mkdir(parents=True, exist_ok=True)
    table = build_subject_roster(metadata, parent_edu_path, volume_dir, meanfd_path)
    path = input_dir / "design.csv"
    table.to_csv(path, index=False)
    _remove_obsolete_before_design(input_dir)
    n_in = int(table["in_analysis"].sum())
    stopped = table.loc[table["excluded_by"].astype(str).str.startswith("sample_stop:"), "excluded_by"]
    stop = str(stopped.iloc[0]) if len(stopped) else ""
    print(
        f"Design roster: {len(table)} subjects, {n_in} enter the analysis"
        + (f", sample stops at {stop}" if stop else "")
        + f". Wrote {path}."
    )
    reasons = table["excluded_by"].replace("", "(entered)")
    summary = ", ".join(f"{name}={int(count)}" for name, count in reasons.value_counts().items())
    print(f"Design roster reasons: {summary}")
    return table


def ensure_subject_roster(input_dir: Path, demographics_path: Path, subjects_path: Path,
                          excluded_path: Path, parent_edu_path: Path, meanfd_path: Path,
                          volume_dir: Path) -> None:
    """Upgrade a saved design.csv that lists only the analysis sample.

    The saved analysis ids must match the current source tables. EffectSize
    images are not opened. A design.csv that already has both roster columns
    is left unchanged.
    """
    path = input_dir / "design.csv"
    if path.is_file():
        header = pd.read_csv(path, nrows=0)
        if {"in_analysis", "excluded_by"}.issubset(set(header.columns)):
            _remove_obsolete_before_design(input_dir)
            return
    if not (demographics_path.is_file() and subjects_path.is_file() and excluded_path.is_file()):
        print(f"{input_dir} has no full design roster, and the subject lists are unavailable, so it was not rebuilt.")
        return
    metadata = build_design_from_sources(demographics_path, subjects_path, excluded_path)
    roster = build_subject_roster(metadata, parent_edu_path, volume_dir, meanfd_path)
    if path.is_file():
        saved = pd.read_csv(path, dtype={"sub_id": str})
        saved_ids = saved["sub_id"].astype(str).tolist()
        roster_ids = analysis_design(roster)["sub_id"].astype(str).tolist()
        if saved_ids != roster_ids:
            print(
                "Warning: the saved design.csv does not match the current source tables. "
                "The full roster was not written over it."
            )
            return
    input_dir.mkdir(parents=True, exist_ok=True)
    roster.to_csv(path, index=False)
    _remove_obsolete_before_design(input_dir)
    print(f"Rebuilt design.csv as the full subject roster without reading volume images: {path}")


def apply_meanfd_exclusion(design: pd.DataFrame, max_meanfd: float) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Apply a meanFD ceiling without removing rows from a full roster.

    People already excluded keep their reason. People in the analysis with
    meanFD above the threshold get excluded_by meanFD_gt_<threshold> and
    in_analysis 0. A design without those columns is still returned as the
    kept rows only.
    """
    table = design.copy()
    table["sub_id"] = table["sub_id"].astype(str)
    threshold = float(max_meanfd)
    if "in_analysis" not in table.columns or "excluded_by" not in table.columns:
        motion = pd.to_numeric(table["meanFD"], errors="raise")
        keep = motion.to_numpy(dtype=float) <= threshold
        dropped = table.loc[~keep, ["sub_id", "meanFD"]].copy()
        kept = table.loc[keep].reset_index(drop=True)
        return kept, dropped.reset_index(drop=True)
    entered = pd.to_numeric(table["in_analysis"], errors="coerce").fillna(0).to_numpy() == 1
    motion = pd.to_numeric(table["meanFD"], errors="coerce").to_numpy(dtype=float)
    over = entered & np.isfinite(motion) & (motion > threshold)
    table.loc[over, "in_analysis"] = 0
    table.loc[over, "excluded_by"] = f"meanFD_gt_{threshold:g}"
    dropped = table.loc[over, ["sub_id", "meanFD"]].copy()
    return table.reset_index(drop=True), dropped.reset_index(drop=True)


def _finalize_design(metadata: pd.DataFrame) -> pd.DataFrame:
    required = ["sub_id", "age", "gender", "site_id", "parent_edu", "meanFD"]
    missing = [column for column in required if column not in metadata.columns]
    if missing:
        raise ValueError(f"Design is missing columns: {missing}")
    metadata = metadata.loc[:, required].copy()
    metadata["sub_id"] = metadata["sub_id"].astype(str)
    for column in ("age", "gender", "site_id", "parent_edu", "meanFD"):
        metadata[column] = pd.to_numeric(metadata[column], errors="raise")
    if metadata["sub_id"].duplicated().any():
        raise ValueError("Design must contain one row per subject")
    if set(metadata["gender"].unique()) != {0, 1}:
        raise ValueError("gender codes must be 0 and 1")
    numeric = ["age", "gender", "site_id", "parent_edu", "meanFD"]
    if not np.isfinite(metadata[numeric].to_numpy(dtype=float)).all():
        raise ValueError("Design columns age, gender, site_id, parent_edu, and meanFD must be finite")
    return metadata.reset_index(drop=True)


def read_design_tsv(input_path: Path) -> pd.DataFrame:
    """Read a finished design table. It must already contain parent_edu and meanFD."""
    table = _read_table(input_path)
    if "site" in table.columns and "site_id" not in table.columns:
        table = table.rename(columns={"site": "site_id"})
    return _finalize_design(table)


def build_design_from_sources(demographics_path: Path, subjects_path: Path,
                              excluded_path: Path) -> pd.DataFrame:
    """Keep the original GAM sample: valid list minus the second-level removals.

    Parental education is attached afterwards. meanFD is attached after the
    volume-file filter and must match every remaining subject.
    """
    for path in (demographics_path, subjects_path, excluded_path):
        if not path.is_file():
            raise FileNotFoundError(f"Design source not found: {path}")
    people = _read_table(demographics_path)
    if "site" in people.columns and "site_id" not in people.columns:
        people = people.rename(columns={"site": "site_id"})
    included = _read_id_list(subjects_path)
    excluded = {_subject_key(subject) for subject in _read_id_list(excluded_path)}
    order = [subject for subject in included if _subject_key(subject) not in excluded]
    people = people.copy()
    people["_key"] = people["sub_id"].map(_subject_key)
    people = people.drop_duplicates("_key").set_index("_key")
    absent = [subject for subject in order if _subject_key(subject) not in people.index]
    if absent:
        raise ValueError(
            f"{len(absent)} subjects are missing from the demographics table. Examples: {absent[:8]}"
        )
    keys = [_subject_key(subject) for subject in order]
    design = people.loc[keys, ["age", "gender", "site_id"]].copy()
    design.insert(0, "sub_id", order)
    design = design.reset_index(drop=True)
    design["sub_id"] = design["sub_id"].astype(str)
    for column in ("age", "gender", "site_id"):
        design[column] = pd.to_numeric(design[column], errors="raise")
    if set(design["gender"].unique()) != {0, 1}:
        raise ValueError("gender codes must be 0 and 1")
    return design


def _volume_subject_dir(volume_dir: Path, subject: str) -> Path | None:
    for name in (str(subject), f"sub-{_subject_key(subject)}"):
        folder = volume_dir / name
        if folder.is_dir():
            return folder
    return None


def _same_grid(image: nib.spatialimages.SpatialImage, reference: nib.spatialimages.SpatialImage) -> bool:
    return image.shape == reference.shape and np.allclose(image.affine, reference.affine)


def effect_volume_paths(volume_dir: Path, metadata: pd.DataFrame) -> tuple[pd.DataFrame, list[Path], list[Path], list[str]]:
    """Keep subjects who have the two effectSize volumes used for Reappraise minus Look."""
    if not volume_dir.is_dir():
        raise FileNotFoundError(f"Volume directory not found: {volume_dir}")
    kept = []
    rpsl_paths = []
    lkng_paths = []
    missing = []
    for row in metadata.itertuples(index=False):
        folder = _volume_subject_dir(volume_dir, row.sub_id)
        rpsl = folder / EFFECT_VOLUME_FILES["rpsl"] if folder is not None else None
        lkng = folder / EFFECT_VOLUME_FILES["lkng"] if folder is not None else None
        if rpsl is not None and rpsl.is_file() and lkng is not None and lkng.is_file():
            kept.append(row._asdict() if hasattr(row, "_asdict") else row)
            rpsl_paths.append(rpsl)
            lkng_paths.append(lkng)
        else:
            missing.append(str(row.sub_id))
    if not rpsl_paths:
        example = volume_dir / str(metadata["sub_id"].iloc[0]) / EFFECT_VOLUME_FILES["rpsl"]
        raise FileNotFoundError(f"No subject has the effectSize volumes. Expected files such as {example}")
    if missing:
        preview = ", ".join(missing[:30])
        print(
            f"Left out {len(missing)} subjects with no effectSize volumes; {len(rpsl_paths)} subjects remain. "
            f"Ids: {preview}"
        )
    design = pd.DataFrame(kept)
    if "Index" in design.columns:
        design = design.drop(columns=["Index"])
    return design.reset_index(drop=True), rpsl_paths, lkng_paths, missing


def volume_paths(metadata: pd.DataFrame, volume_dir: Path) -> dict[str, list[Path]]:
    """Resolve either Fig. 3 beta_0001 names or the original effectSize names."""
    if not volume_dir.is_dir():
        raise FileNotFoundError(f"Volume directory not found: {volume_dir}")
    subjects = metadata["sub_id"].astype(str).tolist()
    chosen = {}
    first_dir = volume_dir / subjects[0]
    for condition, names in VOLUME_FILENAMES.items():
        found = next((name for name in names if (first_dir / name).is_file()), None)
        if found is None:
            raise FileNotFoundError(
                f"Missing {condition} volume for {subjects[0]}. "
                f"Looked for {', '.join(names)} in {first_dir}"
            )
        chosen[condition] = found
    paths = {condition: [] for condition in CONDITIONS}
    for sub_id in subjects:
        for condition, filename in chosen.items():
            path = volume_dir / sub_id / filename
            if not path.is_file():
                raise FileNotFoundError(f"Missing {condition} volume for {sub_id}: {path}")
            paths[condition].append(path)
    return paths


def make_toy_volumes(directory: Path, volume_atlas: Path, n_subjects: int = 36) -> tuple[pd.DataFrame, dict[str, list[Path]]]:
    """Small deterministic cohort on the Schaefer volume grid."""
    rng = np.random.default_rng(20260925)
    metadata = pd.DataFrame({
        "sub_id": [f"toy-{i:03d}" for i in range(n_subjects)],
        "age": rng.uniform(6.2, 18.8, n_subjects),
        "gender": np.tile([0, 1], n_subjects // 2),
        "site_id": np.tile([1, 2, 3, 4], n_subjects // 4),
        "meanFD": rng.uniform(0.04, 0.24, n_subjects),
        "parent_edu": np.tile([1, 2, 3, 4, 5, 6, 7], n_subjects // 7 + 1)[:n_subjects],
    })
    atlas_img = nib.load(str(volume_atlas))
    atlas = np.rint(atlas_img.get_fdata()).astype(np.int32)
    roi_ids = np.unique(atlas[atlas > 0])
    roi_code = np.zeros(int(atlas.max()) + 1, dtype=np.float32)
    roi_code[roi_ids] = rng.normal(0.3, 0.2, roi_ids.size)
    base_pattern = roi_code[np.clip(atlas, 0, roi_code.size - 1)]
    paths = {condition: [] for condition in CONDITIONS}
    for _index, row in metadata.iterrows():
        subject_dir = directory / str(row.sub_id)
        subject_dir.mkdir(parents=True, exist_ok=True)
        sampled = rng.normal(0, 0.03, roi_code.size)[np.clip(atlas, 0, roi_code.size - 1)]
        baseline = base_pattern + sampled
        age_effect = 0.2 * np.sin((row.age - 12.0) / 3.0 + roi_code[np.clip(atlas, 0, roi_code.size - 1)] * 5)
        condition_maps = {
            "rpsl": baseline + 0.35 * base_pattern + age_effect,
            "lkng": baseline,
            "lknt": baseline - 0.1 * base_pattern,
        }
        for condition, beta_index in zip(CONDITIONS, (1, 2, 3)):
            path = subject_dir / f"beta_{beta_index:04d}.nii.gz"
            nib.save(nib.Nifti1Image(condition_maps[condition].astype(np.float32),
                                     atlas_img.affine, atlas_img.header), str(path))
            paths[condition].append(path)
    return metadata, paths


def write_inputs(output_dir: Path, metadata: pd.DataFrame, parcel_response: pd.DataFrame,
                 parcel_info: pd.DataFrame, group_response: pd.DataFrame) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    design_path = output_dir / "design.csv"
    header = pd.read_csv(design_path, nrows=0) if design_path.is_file() else pd.DataFrame()
    roster_present = {"in_analysis", "excluded_by"}.issubset(set(header.columns))
    if roster_present:
        roster = pd.read_csv(design_path, dtype={"sub_id": str})
        entered = analysis_design(roster)["sub_id"].astype(str).tolist()
        model_ids = metadata["sub_id"].astype(str).tolist()
        if entered != model_ids:
            raise ValueError(
                "design.csv rows with in_analysis = 1 do not match the subjects whose volumes were read"
            )
    else:
        design_columns = ["sub_id", "age", "gender", "site_id", "parent_edu", "meanFD"]
        design = metadata.loc[:, design_columns].copy()
        design["sub_id"] = design["sub_id"].astype(str)
        if not np.isfinite(design[["age", "gender", "site_id", "parent_edu", "meanFD"]].to_numpy(dtype=float)).all():
            raise ValueError("Design columns age, gender, site_id, parent_edu, and meanFD must be finite")
        design.to_csv(design_path, index=False)
    parcel_response.to_csv(output_dir / "parcel_response.csv", index=False)
    parcel_info.to_csv(output_dir / "parcel_info.csv", index=False)
    group_response.to_csv(output_dir / "group_response.csv", index=False)
    if not group_response.empty:
        pd.DataFrame({"group": list(YEO_ORDER)}).to_csv(output_dir / "yeo_order.csv", index=False)
        pd.DataFrame({"group": list(SYSTEM_ORDER)}).to_csv(output_dir / "system_order.csv", index=False)


def _masks_on_image_grid(mask_dir: Path, reference: nib.spatialimages.SpatialImage) -> dict[str, np.ndarray]:
    """System masks on the effect-size grid. Each mask is resampled at most once."""
    masks = {}
    for name in SYSTEM_ORDER:
        path = mask_dir / SYSTEM_MASK_NAMES[name]
        if not path.is_file():
            raise FileNotFoundError(
                f"Missing system mask {path}. Expected mask_<system>.nii.gz for {', '.join(SYSTEM_ORDER)}."
            )
        image = _resample_to_atlas(nib.load(str(path)), reference, "nearest", path, 0.0, "effect-size grid")
        mask = np.rint(np.asarray(image.get_fdata())).astype(np.int32) > 0
        print(f"{name}: {int(mask.sum())} voxels on the effect-size grid")
        if int(mask.sum()) == 0:
            raise ValueError(f"System mask {name} has no positive voxels on the effect-size grid: {path}")
        masks[name] = mask
    return masks


def prepare_from_effect_volumes(metadata: pd.DataFrame, rpsl_paths: list[Path], lkng_paths: list[Path],
                                volume_atlas: Path, dlabel_path: Path, system_mask_dir: Path,
                                output_dir: Path,
                                extra_masks: dict[str, nib.spatialimages.SpatialImage] | None = None,
                                min_extra_voxels: int = 5,
                                parcel_info_table: pd.DataFrame | None = None,
                                include_groups: bool = True,
                                extra_parcel_info: pd.DataFrame | None = None) -> dict[str, int]:
    """Parcel and system means from res_1stlevel effectSize volumes.

    Subject maps stay on their own grid. The Schaefer volume atlas is moved
    onto that grid once with nearest-neighbor interpolation, then each parcel
    mean uses the voxels that received its label. System means are the mean of
    finite contrast voxels inside each mask. extra_masks, when supplied, are
    resampled once onto the same grid and averaged the same way. Added parcel
    metadata is retained when include_groups disables network and system means.
    """
    if len(rpsl_paths) != len(metadata) or len(lkng_paths) != len(metadata):
        raise ValueError("Effect-size paths and the design must have one row per subject")
    reference = nib.load(str(rpsl_paths[0]))
    atlas_img = nib.load(str(volume_atlas))
    print(
        f"Effect-size grid {tuple(int(v) for v in reference.shape)}. "
        f"Schaefer atlas {tuple(int(v) for v in atlas_img.shape)} is labeled onto it once."
    )
    labeled = _resample_to_atlas(atlas_img, reference, "nearest", volume_atlas, 0.0, "effect-size grid")
    labels = np.rint(np.asarray(labeled.get_fdata())).astype(np.int32)
    roi_ids = np.unique(labels[labels > 0]).astype(np.int32)
    if roi_ids.size == 0:
        raise ValueError(f"Schaefer atlas has no positive labels on the effect-size grid: {volume_atlas}")
    flat_labels = labels.ravel()
    present = flat_labels > 0
    roi_index = np.searchsorted(roi_ids, flat_labels[present])
    voxel_counts = np.bincount(roi_index, minlength=roi_ids.size).astype(float)
    if include_groups:
        masks = _masks_on_image_grid(system_mask_dir, reference)
    else:
        if parcel_info_table is None and extra_parcel_info is None:
            raise ValueError("Parcel-only volume analysis needs a parcel info table")
        masks = {}
    flat_masks = {system: mask.ravel() for system, mask in masks.items()}
    extra_flat: dict[str, np.ndarray] = {}
    if extra_masks:
        output_dir.mkdir(parents=True, exist_ok=True)
        mask_dir = output_dir / "amygdala_masks"
        mask_dir.mkdir(parents=True, exist_ok=True)
        for name, image in extra_masks.items():
            resampled = _resample_to_atlas(
                image, reference, "nearest", Path(name), 0.0, "effect-size grid")
            mask = np.rint(np.asarray(resampled.get_fdata())).astype(np.int32) > 0
            print(f"{name}: {int(mask.sum())} voxels on the effect-size grid")
            nib.save(nib.Nifti1Image(mask.astype(np.uint8), reference.affine, reference.header),
                     str(mask_dir / f"{name}.nii.gz"))
            extra_flat[name] = mask.ravel()
        count_rows = [{"roi": name, "n_voxels": int(extra_flat[name].sum())} for name in extra_masks]
        pd.DataFrame(count_rows).to_csv(output_dir / "amygdala_voxel_counts.csv", index=False)
        overlap = int(np.sum(np.vstack(list(extra_flat.values())).sum(axis=0) > 1))
        if overlap:
            raise ValueError(
                f"{overlap} effect-size voxels belong to more than one amygdala ROI "
                "after nearest-neighbor resampling"
            )
        small = [row for row in count_rows if row["n_voxels"] < min_extra_voxels]
        if small:
            preview = ", ".join(f"{row['roi']}={row['n_voxels']}" for row in small)
            raise ValueError(
                "Amygdala ROI voxel count is below "
                f"{min_extra_voxels} after resampling onto the effect-size grid: {preview}"
            )
    n_subjects = len(metadata)
    values = np.empty((n_subjects, roi_ids.size), dtype=np.float64)
    system_values = {system: np.empty(n_subjects, dtype=np.float64) for system in SYSTEM_ORDER}
    extra_values = {
        name: np.empty(n_subjects, dtype=np.float64) for name in extra_flat
    }
    for subject_idx, (rpsl_path, lkng_path) in enumerate(zip(rpsl_paths, lkng_paths)):
        if subject_idx % 100 == 0:
            print(f"Volume contrast {subject_idx + 1} / {n_subjects}")
        rpsl_img = nib.load(str(rpsl_path))
        lkng_img = nib.load(str(lkng_path))
        if not _same_grid(rpsl_img, reference) or not _same_grid(lkng_img, reference):
            raise ValueError(
                f"Effect-size grid differs from {rpsl_paths[0].name}: {rpsl_path} has shape {rpsl_img.shape}"
            )
        delta = rpsl_img.get_fdata(dtype=np.float32).ravel() - lkng_img.get_fdata(dtype=np.float32).ravel()
        weights = delta[present]
        finite = np.isfinite(weights)
        counts = np.bincount(roi_index[finite], minlength=roi_ids.size).astype(float)
        empty = np.flatnonzero(counts == 0)
        if empty.size:
            raise ValueError(f"Parcel {int(roi_ids[empty[0]])} has no finite contrast voxels: {rpsl_path}")
        sums = np.bincount(roi_index[finite], weights=weights[finite], minlength=roi_ids.size)
        values[subject_idx, :] = sums / counts
        for system, flat_mask in flat_masks.items():
            observed = delta[flat_mask]
            observed = observed[np.isfinite(observed)]
            if observed.size == 0:
                raise ValueError(f"System mask {system} has no finite contrast voxels: {rpsl_path}")
            system_values[system][subject_idx] = float(np.mean(observed))
        for name, flat_mask in extra_flat.items():
            observed = delta[flat_mask]
            observed = observed[np.isfinite(observed)]
            if observed.size == 0:
                raise ValueError(f"Amygdala ROI {name} has no finite contrast voxels: {rpsl_path}")
            extra_values[name][subject_idx] = float(np.mean(observed))
    parcel_response = pd.DataFrame(values, columns=[str(int(parcel)) for parcel in roi_ids])
    parcel_response.insert(0, "sub_id", metadata["sub_id"].astype(str).to_numpy())
    if include_groups or extra_parcel_info is not None:
        parcel_info = schaefer_parcel_info(dlabel_path)
        missing_names = sorted(set(int(parcel) for parcel in roi_ids) - set(parcel_info["parcel"].astype(int)))
        if extra_parcel_info is not None and missing_names:
            extra = extra_parcel_info.copy()
            extra["parcel"] = extra["parcel"].astype(int)
            covered = set(extra["parcel"].astype(int))
            still_missing = [parcel for parcel in missing_names if parcel not in covered]
            if still_missing:
                raise ValueError(
                    "Volume parcel ids are absent from the Schaefer dlabel and the added parcel table: "
                    f"{still_missing[:8]}"
                )
            overlap_ids = set(parcel_info["parcel"].astype(int)) & covered
            if overlap_ids:
                preview = ", ".join(str(parcel) for parcel in sorted(overlap_ids)[:8])
                raise ValueError(f"Added parcel ids collide with Schaefer ids: {preview}")
            if "network" not in extra.columns:
                extra["network"] = ""
            bad_network = sorted(set(extra["network"].fillna("").astype(str)) & set(YEO_ORDER))
            if bad_network:
                raise ValueError(f"Added parcel networks collide with Yeo 7: {bad_network}")
            parcel_info = parcel_info.copy()
            parcel_info["region_type"] = "cortical"
            parcel_info["structure"] = ""
            parcel_info["is_amygdala"] = False
            extra_rows = extra[extra["parcel"].isin(missing_names)].copy()
            parcel_info = pd.concat([parcel_info, extra_rows], ignore_index=True)
        elif missing_names:
            raise ValueError(f"Volume parcel ids are absent from the Schaefer dlabel: {missing_names[:8]}")
        parcel_info = parcel_info[parcel_info["parcel"].isin(roi_ids)].copy()
    else:
        parcel_info = _parcel_table_for_labels(parcel_info_table, roi_ids)
    if include_groups:
        assignment, _overlap = parcel_overlap_on_mask_grid(atlas_img, volume_atlas, masks, reference)
        parcel_info = parcel_info.merge(assignment, on="parcel", how="left")
        parcel_info["system"] = parcel_info["system"].fillna("")
        yeo_response = yeo_group_responses(parcel_response, parcel_info, voxel_counts, roi_ids)
        subject_ids = metadata["sub_id"].astype(str).to_numpy()
        system_rows = []
        for system in SYSTEM_ORDER:
            for sub_id, value in zip(subject_ids, system_values[system]):
                system_rows.append({"sub_id": sub_id, "dimension": "system", "group": system, "y": float(value)})
        group_response = pd.concat([yeo_response, pd.DataFrame(system_rows)], ignore_index=True)
    else:
        parcel_info["system"] = ""
        group_response = pd.DataFrame(columns=["sub_id", "dimension", "group", "y"])
    write_inputs(output_dir, metadata, parcel_response, parcel_info, group_response)
    if extra_values:
        amygdala_response = pd.DataFrame(extra_values)
        amygdala_response.insert(0, "sub_id", metadata["sub_id"].astype(str).to_numpy())
        amygdala_response.to_csv(output_dir / "amygdala_response.csv", index=False)
        pd.DataFrame({"roi": list(extra_values)}).to_csv(output_dir / "amygdala_order.csv", index=False)
    parcel_word = "parcels" if extra_parcel_info is not None else "Schaefer parcels"
    print(f"Volume parcel means use {roi_ids.size} {parcel_word} on the effect-size grid.")
    return {
        "n_subjects": int(len(metadata)),
        "n_parcels": int(len(roi_ids)),
        "n_system_parcels": int((parcel_info["system"].astype(str) != "").sum()),
        "n_amygdala": int(len(extra_values)),
        "input": "volume effectSize",
    }


def prepare_from_surfaces(metadata: pd.DataFrame, metric_paths: list[dict[str, tuple[Path, Path]]],
                          cifti_atlas: Path, volume_atlas: Path, system_mask_dir: Path,
                          output_dir: Path) -> dict[str, int]:
    """Parcel means from the surface GIFTI files used by gam_AmB_surf_roi."""
    parcel_info = schaefer_parcel_info(cifti_atlas)
    labels = {
        "lh": _hemisphere_vertex_labels(cifti_atlas, "lh"),
        "rh": _hemisphere_vertex_labels(cifti_atlas, "rh"),
    }
    roi_ids = np.array(sorted(set(parcel_info["parcel"].astype(int))), dtype=np.int32)
    present = {int(parcel) for hemi in labels.values() for parcel in np.unique(hemi) if parcel > 0}
    missing = sorted(set(int(parcel) for parcel in roi_ids) - present)
    if missing:
        raise ValueError(f"Schaefer names have no surface vertices: {missing[:8]}")
    roi_ids = np.array(sorted(present), dtype=np.int32)
    parcel_info = parcel_info[parcel_info["parcel"].isin(roi_ids)].copy()
    vertex_counts = np.array([
        int(np.sum(labels["lh"] == parcel) + np.sum(labels["rh"] == parcel)) for parcel in roi_ids
    ], dtype=float)
    values = np.empty((len(metadata), roi_ids.size), dtype=np.float64)
    for subject_idx, paths in enumerate(metric_paths):
        if subject_idx % 100 == 0:
            print(f"Surface contrast {subject_idx + 1} / {len(metadata)}")
        deltas = {}
        for hemi, side in (("lh", 0), ("rh", 1)):
            rpsl = _gifti_vertices(paths["rpsl"][side])
            lkng = _gifti_vertices(paths["lkng"][side])
            if rpsl.shape != lkng.shape:
                raise ValueError(f"Reappraise and Look vertex counts differ: {paths['rpsl'][side]}")
            deltas[hemi] = rpsl - lkng
        values[subject_idx, :] = _parcel_vertex_means(labels, deltas, roi_ids)
    if not np.isfinite(values).all():
        raise ValueError("Parcel Reappraise-minus-Look means contain non-finite values")
    parcel_response = pd.DataFrame(values, columns=[str(int(parcel)) for parcel in roi_ids])
    parcel_response.insert(0, "sub_id", metadata["sub_id"].astype(str).to_numpy())
    atlas_img = nib.load(str(volume_atlas))
    masks, mask_reference = _load_native_system_masks(system_mask_dir)
    assignment, overlap = parcel_overlap_on_mask_grid(atlas_img, volume_atlas, masks, mask_reference)
    parcel_info = parcel_info.merge(assignment, on="parcel", how="left")
    parcel_info["system"] = parcel_info["system"].fillna("")
    yeo_response = yeo_group_responses(parcel_response, parcel_info, vertex_counts, roi_ids)
    system_response = system_means_from_parcels(parcel_response, overlap)
    group_response = pd.concat([yeo_response, system_response], ignore_index=True)
    write_inputs(output_dir, metadata, parcel_response, parcel_info, group_response)
    return {
        "n_subjects": int(len(metadata)),
        "n_parcels": int(len(roi_ids)),
        "n_system_parcels": int((parcel_info["system"].astype(str) != "").sum()),
        "input": "surface gifti",
    }


def prepare_from_volumes(metadata: pd.DataFrame, rpsl_paths: list[Path], lkng_paths: list[Path],
                         volume_atlas: Path, dlabel_path: Path, system_mask_dir: Path,
                         output_dir: Path) -> dict[str, int]:
    parcel_info = schaefer_parcel_info(dlabel_path)
    parcel_response, voxel_counts, atlas_img, masks, system_values = parcel_contrast_matrix(
        metadata, rpsl_paths, lkng_paths, volume_atlas, system_mask_dir)
    roi_ids = parcel_response.columns[1:].astype(int).to_numpy()
    missing = sorted(set(int(parcel) for parcel in roi_ids) - set(parcel_info["parcel"].astype(int)))
    if missing:
        raise ValueError(f"Volume atlas parcel ids are absent from the Schaefer dlabel: {missing[:8]}")
    parcel_info = parcel_info[parcel_info["parcel"].isin(roi_ids)].copy()
    system_response, assignment = system_responses_and_assignment(
        metadata, atlas_img, masks, system_values)
    parcel_info = parcel_info.merge(assignment, on="parcel", how="left")
    yeo_response = yeo_group_responses(parcel_response, parcel_info, voxel_counts, roi_ids)
    group_response = pd.concat([yeo_response, system_response], ignore_index=True)
    write_inputs(output_dir, metadata, parcel_response, parcel_info, group_response)
    return {
        "n_subjects": int(len(metadata)),
        "n_parcels": int(len(roi_ids)),
        "n_system_parcels": int((parcel_info["system"].astype(str) != "").sum()),
    }
