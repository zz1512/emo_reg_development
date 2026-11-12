"""Map Donahue 2018 liberal PFC onto Schaefer 400 by midthickness area.

The step stops when the author mask or the midthickness surfaces are absent.
It does not fall back to parcel-name search or to vertex counts.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

from analysis.fig03.predictive_pls.models import PredictionFailure
from analysis.fig03.predictive_pls.precheck import project_path
from analysis.fig03.predictive_pls.run_pls import MEMBERSHIP_RULE, MEMBERSHIP_SOURCE

N_VERTICES = 32492


def vertex_areas(coordinates: np.ndarray, faces: np.ndarray) -> np.ndarray:
    """One-third of each incident face area, on the supplied mesh."""
    coordinates = np.asarray(coordinates, dtype=float)
    faces = np.asarray(faces, dtype=int)
    triangles = coordinates[faces]
    cross = np.cross(triangles[:, 1] - triangles[:, 0], triangles[:, 2] - triangles[:, 0])
    face_area = 0.5 * np.linalg.norm(cross, axis=1)
    area = np.zeros(coordinates.shape[0], dtype=float)
    for corner in range(3):
        np.add.at(area, faces[:, corner], face_area / 3.0)
    return area


def coverage_by_parcel(labels: np.ndarray, inside: np.ndarray, area: np.ndarray) -> dict[int, tuple[float, float]]:
    """Return parcel -> (pfc area, total area) for positive parcel labels."""
    labels = np.asarray(labels, dtype=int)
    inside = np.asarray(inside, dtype=bool)
    area = np.asarray(area, dtype=float)
    if not (labels.shape == inside.shape == area.shape):
        raise PredictionFailure("Parcel labels, PFC mask, and vertex area must share one vertex axis")
    if not np.isfinite(area).all() or np.any(area < 0):
        raise PredictionFailure("Vertex area must be finite and nonnegative")
    totals: dict[int, list[float]] = {}
    for parcel in np.unique(labels):
        if int(parcel) <= 0:
            continue
        selected = labels == parcel
        if not selected.any():
            continue
        totals[int(parcel)] = [float(np.sum(area[selected & inside])), float(np.sum(area[selected]))]
    return {parcel: (values[0], values[1]) for parcel, values in totals.items()}


def membership_from_coverage(coverage: dict[int, tuple[float, float]], info: pd.DataFrame) -> pd.DataFrame:
    info = info.copy()
    info["parcel"] = pd.to_numeric(info["parcel"], errors="raise").astype(int)
    if set(info["parcel"]) != set(range(1, 433)):
        raise PredictionFailure("parcel_info.csv must name parcels 1..432")
    rows = []
    for record in info.sort_values("parcel").itertuples(index=False):
        parcel = int(record.parcel)
        name = str(record.name)
        hemisphere = "subcortical"
        if parcel <= 400:
            matched = pd.Series([name]).str.extract(r"_(LH|RH)_", expand=False).iloc[0]
            hemisphere = "" if pd.isna(matched) else str(matched)
            if parcel not in coverage:
                raise PredictionFailure(f"Schaefer parcel {parcel} has no vertices in the native dlabel")
            pfc_area, total_area = coverage[parcel]
            if total_area <= 0:
                raise PredictionFailure(f"Schaefer parcel {parcel} has no positive midthickness area")
            fraction = pfc_area / total_area
        else:
            fraction = 0.0
        rows.append({
            "parcel": parcel, "name": name, "hemisphere": hemisphere,
            "coverage_fraction": fraction, "is_pfc": bool(parcel <= 400 and fraction >= 0.50),
            "boundary": bool(parcel <= 400 and 0 < fraction < 1),
            "review_band": bool(parcel <= 400 and 0.20 <= fraction <= 0.80),
            "definition_source": MEMBERSHIP_SOURCE, "rule": MEMBERSHIP_RULE,
        })
    table = pd.DataFrame(rows)
    pfc = set(table.loc[table["is_pfc"], "parcel"])
    other = set(table.loc[~table["is_pfc"], "parcel"])
    if pfc & other or pfc | other != set(range(1, 433)):
        raise PredictionFailure("PFC and non-PFC labels do not partition parcels 1..432")
    return table


def _gifti_surface(path: Path):
    import nibabel as nib

    image = nib.load(str(path))
    if not hasattr(image, "agg_data"):
        raise PredictionFailure(f"{path} is not a GIFTI surface")
    coordinates, faces = image.agg_data(("pointset", "triangle"))
    coordinates = np.asarray(coordinates)
    faces = np.asarray(faces)
    if coordinates.shape != (N_VERTICES, 3):
        raise PredictionFailure(f"{path} must be an fsLR32k surface with {N_VERTICES} vertices")
    return coordinates, faces


def _gifti_mask(path: Path) -> np.ndarray:
    import nibabel as nib

    image = nib.load(str(path))
    arrays = getattr(image, "darrays", None)
    if not arrays:
        raise PredictionFailure(f"{path} is not a GIFTI metric")
    values = np.asarray(arrays[0].data, dtype=float).reshape(-1)
    if values.size != N_VERTICES or not np.isfinite(values).all():
        raise PredictionFailure(f"{path} must contain {N_VERTICES} finite PFC-mask values")
    return values > 0


def build_definition(schaefer_dlabel: Path, parcel_info: Path, pfc_left: Path, pfc_right: Path,
                     midthickness_left: Path, midthickness_right: Path, output: Path,
                     inflated_left: Path | None = None, inflated_right: Path | None = None) -> pd.DataFrame:
    from analysis.fig03.gam_pygam.brain_maps import load_surface_atlas

    paths = [schaefer_dlabel, parcel_info, pfc_left, pfc_right, midthickness_left, midthickness_right]
    missing = [str(path) for path in paths if not path.is_file()]
    if missing:
        raise PredictionFailure(
            "Liberal PFC mapping stopped. Required files are absent: "
            + "; ".join(missing)
            + ". Provide the Donahue 2018 bilateral liberal PFC mask on fsLR32k and matching midthickness surfaces. "
            "Parcel-name lists are not a substitute.")
    _, (left_labels, right_labels) = load_surface_atlas(schaefer_dlabel)
    info = pd.read_csv(parcel_info)
    coverage = {}
    for labels, mask_path, surface_path in (
            (left_labels, pfc_left, midthickness_left),
            (right_labels, pfc_right, midthickness_right)):
        area = vertex_areas(*_gifti_surface(surface_path))
        coverage.update(coverage_by_parcel(labels, _gifti_mask(mask_path), area))
    table = membership_from_coverage(coverage, info)
    output.mkdir(parents=True, exist_ok=True)
    table.to_csv(output / "pfc_parcel_membership.csv", index=False)
    provenance = {
        "definition_source": MEMBERSHIP_SOURCE,
        "rule": MEMBERSHIP_RULE,
        "citation": "Donahue et al., 2018, PNAS 115:E5183-E5192; liberal PFC is every area in supplementary Table S2",
        "source_page": "https://pmc.ncbi.nlm.nih.gov/articles/PMC5984508/",
        "n_pfc": int(table["is_pfc"].sum()),
        "n_boundary": int(table["boundary"].sum()),
        "n_review_band": int(table["review_band"].sum()),
        "sha256": {str(path): hashlib.sha256(path.read_bytes()).hexdigest() for path in paths},
        "area": "one-third of incident midthickness face area; inflated surfaces are not used for the fraction",
    }
    (output / "pfc_definition.json").write_text(json.dumps(provenance, indent=2), encoding="utf-8")
    _write_preview(table, output, inflated_left, inflated_right, left_labels, right_labels)
    print(f"PFC membership: {output / 'pfc_parcel_membership.csv'}", flush=True)
    return table


def _write_preview(table, output, inflated_left, inflated_right, left_labels, right_labels):
    import matplotlib
    matplotlib.use("Agg")
    matplotlib.rcParams["svg.fonttype"] = "none"
    import matplotlib.pyplot as plt

    figure, axes = plt.subplots(1, 2, figsize=(8, 3.6))
    if inflated_left and inflated_right and inflated_left.is_file() and inflated_right.is_file():
        import nibabel as nib
        for axis, surface_path, labels, title in (
                (axes[0], inflated_left, left_labels, "Left"),
                (axes[1], inflated_right, right_labels, "Right")):
            coordinates = np.asarray(nib.load(str(surface_path)).agg_data("pointset"))
            color = np.where(np.isin(labels, table.loc[table["is_pfc"], "parcel"]), 1.0, 0.0)
            color[labels == 0] = np.nan
            axis.scatter(coordinates[:, 0], coordinates[:, 2], c=color, s=0.2, cmap="cividis", linewidths=0)
            axis.set_title(f"{title} inflated display")
            axis.set_aspect("equal")
            axis.axis("off")
    else:
        cortical = table.loc[table["parcel"].le(400)]
        axes[0].hist(cortical["coverage_fraction"], bins=20, color="#4C78A8")
        axes[0].set_xlabel("Midthickness PFC fraction")
        axes[0].set_ylabel("Parcels")
        axes[1].axis("off")
        axes[1].text(0.0, 0.5, "Inflated display surfaces were not supplied.\nArea fractions still use midthickness.",
                     transform=axes[1].transAxes)
    figure.tight_layout()
    figure.savefig(output / "pfc_definition_preview.svg", format="svg")
    plt.close(figure)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--schaefer-dlabel", type=Path, required=True)
    parser.add_argument("--parcel-info", type=Path, required=True)
    parser.add_argument("--pfc-left", type=Path, required=True)
    parser.add_argument("--pfc-right", type=Path, required=True)
    parser.add_argument("--midthickness-left", type=Path, required=True)
    parser.add_argument("--midthickness-right", type=Path, required=True)
    parser.add_argument("--inflated-left", type=Path)
    parser.add_argument("--inflated-right", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    build_definition(*(project_path(path) for path in (
        args.schaefer_dlabel, args.parcel_info, args.pfc_left, args.pfc_right,
        args.midthickness_left, args.midthickness_right, args.output)),
        project_path(args.inflated_left) if args.inflated_left else None,
        project_path(args.inflated_right) if args.inflated_right else None)


if __name__ == "__main__":
    main()
