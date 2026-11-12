"""Render clearly labeled synthetic previews with the real age-curve renderer."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import re
import xml.etree.ElementTree as ET

import numpy as np
import pandas as pd

from analysis.fig03.gam_pygam import run_pygam as runner
from analysis.fig03.gam_pygam.plot_region_age import (
    PAPER_PARCELS, PAPER_Y_LIMITS, curve_style_metadata, plot_one_region,
)


def _comparison_svg(output_dir: Path, parcels: list[int]) -> None:
    """Compose vector panels with unique element IDs; no raster output."""
    namespace = "http://www.w3.org/2000/svg"
    ET.register_namespace("", namespace)
    ET.register_namespace("xlink", "http://www.w3.org/1999/xlink")
    panels = [ET.parse(output_dir / f"mock_parcel_{parcel:04d}.svg").getroot() for parcel in parcels]
    sizes = [list(map(float, panel.attrib["viewBox"].split()))[2:] for panel in panels]
    width, height = sum(size[0] for size in sizes), max(size[1] for size in sizes)
    root = ET.Element(f"{{{namespace}}}svg", {"width": f"{width}pt", "height": f"{height}pt",
                                           "viewBox": f"0 0 {width} {height}", "version": "1.1"})
    offset = 0.0
    for parcel, panel, (panel_width, _) in zip(parcels, panels, sizes):
        prefix = f"roi{parcel}_"
        for element in panel.iter():
            for key, value in list(element.attrib.items()):
                if key == "id":
                    element.set(key, prefix + value)
                elif key.endswith("href") and value.startswith("#"):
                    element.set(key, "#" + prefix + value[1:])
                else:
                    element.set(key, re.sub(r"url\(#([^)]+)\)", lambda match: f"url(#{prefix}{match.group(1)})", value))
        group = ET.SubElement(root, f"{{{namespace}}}g", {"transform": f"translate({offset}, 0)"})
        group.extend(child for child in panel if child.tag != f"{{{namespace}}}metadata")
        offset += panel_width
    ET.ElementTree(root).write(output_dir / "mock_comparison.svg", encoding="utf-8", xml_declaration=True)


def render_mock(output_dir: Path, parcels: list[int], n_subjects: int = 991,
                seed: int = 20261108) -> None:
    """Fit mock contrasts, including true mean CI; no real participant data."""
    rng = np.random.default_rng(seed)
    age = []
    while len(age) < n_subjects:
        candidates = rng.normal(12.5, 2.65, n_subjects)
        age.extend(candidates[(candidates >= 6.) & (candidates <= 18.9)].tolist())
    age = np.asarray(age[:n_subjects])
    covariates = np.column_stack([rng.integers(0, 2, n_subjects), rng.integers(1, 8, n_subjects),
                                 rng.integers(1, 6, n_subjects), rng.uniform(.02, .3, n_subjects)])
    nuisance = covariates @ np.array([.08, -.015, .035, -.2])
    grid = runner._age_grid(age)
    lam_grid = 10. ** np.linspace(*runner.LAM_LOG_RANGE)
    stats = runner.load_historical_stats()
    parameters = {152: (.040, .0007), 154: (.034, .0003),
                  182: (-.038, .0005), 185: (-.032, -.0003)}
    names = {152: "LH Default Temp 4", 154: "LH Default Temp 6",
             182: "LH Default PFC 17", 185: "LH Default PFC 20"}
    curve_frames, residual_frames, details = [], [], {}
    for parcel in parcels:
        slope, curvature = parameters[parcel]
        ac = age - age.mean()
        response = slope * ac + curvature * (ac ** 2 - np.mean(ac ** 2)) + nuisance + rng.normal(0, .42, n_subjects)
        _, _, partial = runner._fit_one(stats, response, age, covariates, grid, lam_grid, store_partial=True)
        curve_frames.append(pd.DataFrame({"dimension": "parcel", "name": str(parcel), "age": grid,
                                          "effect": partial["effect"], "ci_low": partial["ci_low"], "ci_high": partial["ci_high"]}))
        residual_frames.append(pd.DataFrame({"dimension": "parcel", "name": str(parcel),
                                             "age": age, "residual": partial["residual"]}))
        details[parcel] = partial
    curves, residuals = pd.concat(curve_frames, ignore_index=True), pd.concat(residual_frames, ignore_index=True)
    output_dir.mkdir(parents=True, exist_ok=True)
    for parcel in parcels:
        partial = details[parcel]
        plot_one_region(age, partial["residual"], grid, partial["effect"], partial["ci_low"], partial["ci_high"],
                        title=f"{parcel}  {names[parcel]}", output_path=output_dir / f"mock_parcel_{parcel:04d}.svg",
                        y_limits=PAPER_Y_LIMITS[parcel], paper_parcel=parcel, annotation="MOCK DATA\n95% mean CI")
    curves.to_csv(output_dir / "MOCK_pygam_curves.csv", index=False)
    residuals.to_csv(output_dir / "MOCK_partial_residuals.csv", index=False)
    metadata = {"data_kind": "SYNTHETIC ONLY; not participant data or research results", "seed": seed,
                "n_subjects_per_parcel": n_subjects, "parcels": parcels,
                "y_limits": {str(parcel): list(PAPER_Y_LIMITS[parcel]) for parcel in parcels},
                "n_outside_y_limits": {str(parcel): int((~residuals.loc[residuals["name"] == str(parcel), "residual"].between(*PAPER_Y_LIMITS[parcel])).sum()) for parcel in parcels},
                "style": curve_style_metadata(parcels),
                "ci": "Real GCV-selected pyGAM mean-inclusive 95% pointwise CI fitted to synthetic responses",
                "real_data_redraw": "python -m analysis.fig03.gam_pygam.run_pygam --step curves --curve-parcels 152 154 182 185 --output outputs/fig03/gam_pygam_volume_432_models3_FD05"}
    (output_dir / "mock_metadata.json").write_text(json.dumps(metadata, indent=2), encoding="utf-8")
    _comparison_svg(output_dir, parcels)
    legacy_png = output_dir / "mock_comparison.png"
    if legacy_png.is_file():
        legacy_png.unlink()
    print(f"SVG-only mock previews: {output_dir.resolve()} (y limits: {metadata['y_limits']})")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=runner.PROJECT_ROOT / "outputs/fig03/age_curves_mock")
    parser.add_argument("--parcels", type=int, nargs="+", choices=PAPER_PARCELS, default=[152, 182])
    args = parser.parse_args()
    render_mock(args.output, args.parcels)


if __name__ == "__main__":
    main()
