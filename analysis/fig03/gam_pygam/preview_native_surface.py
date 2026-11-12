"""Preview real Schaefer parcel geometry; colors identify ROIs, not effects."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib.patches import Patch
import nibabel as nib
import numpy as np
from nilearn import plotting

from analysis.fig03.gam_pygam.brain_maps import DEFAULT_BORDER_SURFACES, load_surface_atlas
from analysis.fig03.gam_pygam.make_selected_dlabel import ATLAS_DIR, PROJECT_ROOT, SELECTED, STEM


def _draw(figure, position, mesh, labels, native_labels, view, title):
    axes = figure.add_subplot(position, projection="3d")
    # The gray parcel pattern is the original native surface atlas. Six colors
    # encode ROI identity only; there are no simulated statistical values.
    color_table = np.ones((401, 4))
    color_table[:, :3] = np.array([.78 + .055 * (key % 3) for key in range(401)])[:, None]
    color_table[0, :3] = .94
    for key, (_, rgb) in SELECTED.items():
        color_table[key, :3] = .58 + .42 * np.asarray(rgb) / 255.
    plotting.plot_surf_roi(mesh, labels, hemi="left", view=view, cmap=ListedColormap(color_table),
                           threshold=None, vmin=0, vmax=400, alpha=1., colorbar=False,
                           axes=axes, figure=figure)
    # These six parcels lie on the lateral surface. The matplotlib contour
    # artist does not depth-test back-facing lines, so omit it in the medial
    # panel rather than showing ghost outlines through the brain.
    if view == "lateral":
        plotting.plot_surf_contours(mesh, native_labels, levels=list(SELECTED),
                                    colors=[tuple(np.asarray(rgb)/255.) for _, rgb in SELECTED.values()],
                                    axes=axes, figure=figure, linewidths=1.8)
    axes.set_box_aspect(axes.get_box_aspect(), zoom=1.20)
    # Keep the preview SVG compact: rasterize the dense 3D mesh at export
    # resolution, while titles, legend and other text remain editable SVG.
    for collection in axes.collections:
        collection.set_rasterized(True)
    axes.set_title(title, fontsize=14, pad=0)
    return axes


def make_preview(output_dir: Path, template_path: Path, legacy_path: Path | None = None) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    _, labels = load_surface_atlas(template_path)
    native = labels[0]
    selected_path = output_dir / f"{STEM}.dlabel.nii"
    _, selected_labels = load_surface_atlas(selected_path)
    expected = np.where(np.isin(native, list(SELECTED)), native, 0)
    np.testing.assert_array_equal(selected_labels[0], expected)
    np.testing.assert_array_equal(selected_labels[1], 0)
    mesh = str(DEFAULT_BORDER_SURFACES[0])
    with plt.rc_context({"font.family": "Arial", "svg.fonttype": "none"}):
        figure = plt.figure(figsize=(12, 5.8), facecolor="white")
        _draw(figure, 121, mesh, native, native, "lateral", "Left hemisphere | lateral")
        _draw(figure, 122, mesh, native, native, "medial", "Left hemisphere | medial")
        figure.suptitle("Native Schaefer400 fsLR32k: selected six ROIs", fontsize=18, y=.97)
        handles = [Patch(facecolor=.58 + .42*np.asarray(rgb)/255., edgecolor=np.asarray(rgb)/255.,
                         label=f"{key}  {name.split('_Default_')[1]}", linewidth=1.8)
                   for key, (name, rgb) in SELECTED.items()]
        figure.legend(handles=handles, loc="lower center", ncol=3, frameon=False, fontsize=12,
                      bbox_to_anchor=(.5, .035))
        figure.text(.5, .015, "Real atlas geometry; colors denote ROI identity, not partial R² or statistical significance.",
                    ha="center", fontsize=10, color="#555555")
        figure.subplots_adjust(left=0, right=1, bottom=.17, top=.90, wspace=-.1)
        figure.savefig(output_dir / "native_selected6_preview.svg", dpi=300)
        plt.close(figure)

        comparisons = []
        if legacy_path is not None and legacy_path.is_file():
            _, old_labels = load_surface_atlas(legacy_path)
            # Keep the same gray atlas background in both panels; only the
            # selected colored regions change. Outlines always use native IDs.
            old_map = native.copy()
            old_map[np.isin(old_map, list(SELECTED))] = 0
            present = old_labels[0] > 0
            old_map[present] = old_labels[0][present]
            figure = plt.figure(figsize=(12, 5.2), facecolor="white")
            _draw(figure, 121, mesh, old_map, native, "lateral", "Previous: volume-derived colored ROIs")
            _draw(figure, 122, mesh, native, native, "lateral", "Updated: native surface colored ROIs")
            figure.suptitle("Same native outlines in both panels", fontsize=18, y=.98)
            figure.text(.5, .06, "Left: color and native outlines may disagree.   Right: colors and outlines share the same native parcel IDs.",
                        ha="center", fontsize=11)
            figure.subplots_adjust(left=0, right=1, bottom=.09, top=.88, wspace=-.1)
            figure.savefig(output_dir / "native_vs_volume_alignment.svg", dpi=300)
            plt.close(figure)
            for key in SELECTED:
                old, new = old_labels[0] == key, native == key
                comparisons.append({"parcel": key, "old_volume_projected_vertices": int(old.sum()),
                                    "native_vertices": int(new.sum()),
                                    "old_native_overlap": int((old & new).sum()),
                                    "new_native_mismatches": int(np.count_nonzero((selected_labels[0] == key) != new))})
    (output_dir / "preview_validation.json").write_text(json.dumps({
        "source_surface": str(template_path.resolve()), "dlabel": str(selected_path.resolve()),
        "selected_matches_native_at_every_vertex": True,
        "note": "Real anatomical atlas preview; no subject data or simulated effect sizes; not a Workbench scene.",
        "comparisons": comparisons,
    }, indent=2), encoding="utf-8")
    print(output_dir / "native_selected6_preview.svg")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--template", type=Path, default=ATLAS_DIR / "Schaefer2018_400Parcels_7Networks_order.dlabel.nii")
    parser.add_argument("--output", type=Path, default=PROJECT_ROOT / "outputs/fig03/selected_roi_dlabel_native")
    parser.add_argument("--legacy-dlabel", type=Path, default=PROJECT_ROOT / "outputs/fig03/selected_roi_dlabel/Schaefer400_7Networks_selected6_fromMNI_fslr32k.dlabel.nii")
    args = parser.parse_args()
    make_preview(args.output, args.template, args.legacy_dlabel)


if __name__ == "__main__":
    main()
