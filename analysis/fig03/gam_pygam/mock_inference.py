"""Deterministic, explicitly synthetic 449-unit inference/table demonstration.

This creates CSVs only. Export Word through result_tables in a document
environment; no real participant data or real brain results are used.
"""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from analysis.fig03.gam_pygam.age_bayes import fit_age_bayes
from analysis.fig03.gam_pygam import run_erq
from analysis.fig03.gam_pygam.result_tables import write_erq_tables, write_age_tables


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=Path("outputs/fig03/inference_mock"))
    args = parser.parse_args()
    output = args.output
    output.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(20261110)
    n = 120
    age = np.linspace(6.2, 18.8, n)
    ac = age - age.mean()
    erq = 3 + .035 * ac + rng.normal(0, .7, n)
    ec = erq - erq.mean()
    sample = pd.DataFrame({"gender": np.arange(n) % 2, "parent_edu": rng.normal(4, 1, n),
                           "meanFD": rng.uniform(.02, .2, n),
                           "site_id": np.array(["1", "2", "3"])[np.arange(n) % 3]})
    covariates, _, _ = run_erq.build_covariates(sample)
    units = [("parcel", str(number)) for number in range(1, 433)]
    units += [("yeo", name) for name in run_erq.YEO_ORDER]
    units += [("system", name) for name in run_erq.SYSTEM_ORDER]
    units += [("amygdala", name) for name in run_erq.AMYGDALA_ORDER]
    erq_rows, age_rows, sensitivities, slope_frames = [], [], [], []
    for index, (dimension, name) in enumerate(units, start=1):
        y = .015 * ac + rng.normal(0, .3, n)
        if dimension == "parcel" and name in ("152", "154"):
            y += .8 * ec
        if dimension == "parcel" and name in ("182", "185"):
            y += ec * (.08 * ac + .02 * ac ** 2)
        if dimension == "parcel" and name == "403":
            y += .04 * ac ** 2
        key = {"dimension": dimension, "name": name, "label": f"SIMULATED {dimension} {name}"}
        erq_result, details = run_erq.fit_unit(y, age, erq, covariates, lam_grid=np.array([1e-3, 1e3]))
        erq_rows.append({**key, **erq_result})
        if details is not None:
            slope_frames.append(run_erq.reference_slopes(dimension, name, details))
        age_result, sensitivity = fit_age_bayes(y, age, covariates)
        age_rows.append({**key, **age_result})
        sensitivities.extend({**key, **row} for row in sensitivity)
        if index % 50 == 0 or index == len(units):
            print(f"Synthetic inference units: {index}/{len(units)}", flush=True)
    erq_statistics = run_erq.correct_families(pd.DataFrame(erq_rows))
    age_statistics = pd.DataFrame(age_rows)
    for kind, statistics in (("erq", erq_statistics), ("age", age_statistics)):
        directory = output / kind
        directory.mkdir(parents=True, exist_ok=True)
        statistics.to_csv(directory / ("erq_age_statistics.csv" if kind == "erq" else "age_bayes_statistics.csv"), index=False)
    pd.DataFrame(sensitivities).to_csv(output / "age/age_bayes_prior_sensitivity.csv", index=False)
    pd.concat(slope_frames, ignore_index=True).to_csv(output / "erq/erq_reference_slopes_at_ages.csv", index=False)
    write_erq_tables(erq_statistics, output / "erq", include_docx=False, mock=True)
    write_age_tables(age_statistics, output / "age", include_docx=False, mock=True)
    # A small, deliberately labelled preview allows all pages to be visually
    # checked; the full tables retain all 449 units for later selection.
    for kind, statistics in (("erq", erq_statistics), ("age", age_statistics)):
        preview = statistics.loc[((statistics["dimension"] == "parcel") & statistics["name"].isin(["152", "182", "403"]))
                                 | statistics["dimension"].ne("parcel")].groupby("dimension", sort=False).head(3)
        preview.to_csv(output / kind / "preview_statistics.csv", index=False)
    (output / "README.txt").write_text(
        "SIMULATED DATA ONLY. 449 units, 120 synthetic participants, two-value GCV grids.\n"
        "These are implementation/layout demonstrations, not real neural findings.\n"
        "Full CSVs retain all 449 units; preview_statistics.csv is a small selection for visual QA.\n", encoding="utf-8")


if __name__ == "__main__":
    main()
