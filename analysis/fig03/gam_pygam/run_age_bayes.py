"""Add formal age-spline Bayes factors to saved Tian432/232 FD <= 0.5 inputs."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

from analysis.fig03.gam_pygam.age_bayes import MAIN_SCALE, PRIOR_SCALES, fit_age_bayes
from analysis.fig03.gam_pygam.run_pygam import (
    COVARIATES, SYSTEM_ORDER, YEO_ORDER, analysis_design, resolve_project_path, validate_saved_inputs,
)


def load_age_units(input_dir: Path, design: pd.DataFrame, n_cortical: int):
    """Read each cached response exactly once, joining by string subject IDs."""
    info = pd.read_csv(input_dir / "parcel_info.csv").set_index("parcel")["name"].astype(str).to_dict()

    def align(table, columns):
        if table["sub_id"].duplicated().any():
            raise ValueError("Saved response has duplicate subject IDs")
        return design[["sub_id"]].merge(table[["sub_id", *columns]], on="sub_id", how="left", validate="one_to_one")

    parcels = pd.read_csv(input_dir / "parcel_response.csv", dtype={"sub_id": str})
    names = [str(number) for number in range(1, n_cortical + 33)]
    matched = align(parcels, names)
    for name in names:
        yield "parcel", name, f"{name} {info.get(int(name), '')}".strip(), matched[name].to_numpy(dtype=float)
    if n_cortical == 200:
        return
    groups = pd.read_csv(input_dir / "group_response.csv", dtype={"sub_id": str})
    for dimension, names in (("yeo", YEO_ORDER), ("system", SYSTEM_ORDER)):
        for name in names:
            matched = align(groups.loc[(groups["dimension"] == dimension) & (groups["group"] == name)], ["y"])
            yield dimension, str(name), f"{dimension} {name}", matched["y"].to_numpy(dtype=float)
    from analysis.fig03.gam_nested.amygdala import AMYGDALA_ORDER

    amygdala = pd.read_csv(input_dir / "amygdala_response.csv", dtype={"sub_id": str})
    matched = align(amygdala, list(AMYGDALA_ORDER))
    for name in AMYGDALA_ORDER:
        yield "amygdala", name, f"Juelich {name}", matched[name].to_numpy(dtype=float)


def run_analysis(input_dir: Path, output_dir: Path, skip_docx: bool = False) -> pd.DataFrame:
    from analysis.fig03.gam_pygam.result_tables import check_docx_dependency, write_age_tables

    input_dir, output_dir = input_dir.resolve(), output_dir.resolve()
    if output_dir in (input_dir, input_dir.parent) or input_dir in output_dir.parents:
        raise ValueError("Bayesian output must not overwrite or be nested inside source inputs")
    metadata_path = output_dir / "run_metadata.json"
    if metadata_path.is_file() and json.loads(metadata_path.read_text(encoding="utf-8")).get("analysis_kind") != "age_bayesian_spline_v1":
        raise ValueError("Bayesian output contains another analysis; choose a separate bayes directory")
    note = validate_saved_inputs(input_dir)
    if note is None or note.get("n_cortical") not in (200, 400):
        raise ValueError("Age Bayes factors require Tian232 or Tian432 saved inputs")
    design = analysis_design(pd.read_csv(input_dir / "design.csv", dtype={"sub_id": str})).copy()
    if design["sub_id"].isna().any() or design["sub_id"].duplicated().any():
        raise ValueError("Age design must have unique nonmissing subject IDs")
    if design.empty or not np.isfinite(design["meanFD"]).all() or (design["meanFD"] > .5).any():
        raise ValueError("Only the existing FD <= 0.5 retained sample is authorized; check the source cache")
    if not skip_docx:
        check_docx_dependency()
    output_dir.mkdir(parents=True, exist_ok=True)
    age = design["age"].to_numpy(dtype=float)
    covariates = design[list(COVARIATES)].to_numpy(dtype=float)
    rows, sensitivities, exclusions = [], [], []
    units = list(load_age_units(input_dir, design, int(note["n_cortical"])))
    for index, (dimension, name, label, response) in enumerate(units, start=1):
        finite = np.isfinite(response) & np.isfinite(age) & np.isfinite(covariates).all(axis=1)
        exclusions.extend({"dimension": dimension, "name": name, "sub_id": subject,
                           "reason": "nonfinite_response_age_or_covariate", "models": "m1_m2_m3"}
                          for subject in design.loc[~finite, "sub_id"])
        row, sensitivity = fit_age_bayes(response, age, covariates)
        key = {"dimension": dimension, "name": name, "label": label}
        rows.append({**key, **row})
        sensitivities.extend({**key, **item} for item in sensitivity)
        if index % 25 == 0 or index == len(units):
            print(f"Age Bayesian spline comparisons: {index}/{len(units)}", flush=True)
    statistics = pd.DataFrame(rows)
    original_path = input_dir.parent / "statistics" / "pygam_age_comparison.csv"
    if original_path.is_file():
        original = pd.read_csv(original_path, dtype={"name": str})
        joined = statistics.merge(original, on=["dimension", "name"], how="left", validate="one_to_one", suffixes=("", "_frequentist"))
        if joined["n_obs_frequentist"].isna().any() or not np.array_equal(joined["n_obs"], joined["n_obs_frequentist"]):
            raise ValueError("Bayesian and saved frequentist samples differ; verify the age cache/results pairing")
        statistics = joined
    statistics.to_csv(output_dir / "age_bayes_statistics.csv", index=False)
    sensitivity_table = pd.DataFrame(sensitivities) if sensitivities else pd.DataFrame(columns=[
        "dimension", "name", "label", "linear_prior_scale", "nonlinear_prior_scale",
        "log_evidence_m1_relative", "log_evidence_m2_relative", "log_evidence_m3_relative",
        "quadrature_relative_error", "quadrature_check_log_difference", "integration_converged"])
    sensitivity_table.to_csv(output_dir / "age_bayes_prior_sensitivity.csv", index=False)
    pd.DataFrame(exclusions, columns=["dimension", "name", "sub_id", "reason", "models"]).to_csv(
        output_dir / "unit_response_exclusions.csv", index=False)
    design.to_csv(output_dir / "analysis_design.csv", index=False)
    metadata = {
        "analysis_kind": "age_bayesian_spline_v1", "source_inputs": str(input_dir),
        "n_cortical": int(note["n_cortical"]), "max_meanFD": .5,
        "models": {"m1": "common nuisance only", "m2": "nuisance + linear age",
                   "m3": "nuisance + linear age + second-difference spline nonlinearity"},
        "nuisance_encoding": "same numeric four covariates as the source age pyGAM analysis; no silent site recoding",
        "basis": "same five cubic pyGAM B splines; original second-difference penalty; remove intercept/linear components; normalize nonlinear prior RMS to 1",
        "prior": {"linear": "beta_age | sigma ~ Normal(0, (r_linear*sigma)^2); age has population SD 1",
                  "nonlinear": "penalty-whitened u | tau,sigma ~ Normal(0, sigma^2*tau^2 I); tr(Z Z')/n=1",
                  "strength": "tau ~ HalfNormal(0, r_nonlinear)",
                  "common_reference": "same Lebesgue nuisance measure and d(sigma^2)/sigma^2 in every model; constants cancel; only relative evidences are defined",
                  "main_scales": [MAIN_SCALE, MAIN_SCALE], "sensitivity_scales": list(PRIOR_SCALES)},
        "integration": "coefficients and variance analytic; half-normal strength adaptive quadrature with independent partition and tolerance check; not GCV plug-in or BIC",
        "bayes_factor_reporting": "raw pairwise BF and log BF; evidence labels at 3 and 10 are descriptive, not multiplicity-corrected discoveries; no BH applied to BF",
        "source_hashes": {path.name: hashlib.sha256(path.read_bytes()).hexdigest() for path in (
            input_dir / "design.csv", input_dir / "parcel_response.csv")},
        "n_units": statistics.groupby("dimension").size().to_dict(), "statuses": statistics["status"].value_counts().to_dict(),
    }
    (output_dir / "run_metadata.json").write_text(json.dumps(metadata, indent=2, ensure_ascii=False), encoding="utf-8")
    write_age_tables(statistics, output_dir, include_docx=not skip_docx)
    return statistics


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--inputs", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--skip-docx", action="store_true")
    parser.add_argument("--step", choices=("fit", "tables"), default="fit")
    args = parser.parse_args()
    output = resolve_project_path(args.output)
    if args.step == "tables":
        from analysis.fig03.gam_pygam.result_tables import write_age_tables

        write_age_tables(pd.read_csv(output / "age_bayes_statistics.csv", dtype={"name": str}), output, not args.skip_docx)
    else:
        if args.inputs is None:
            parser.error("--inputs is required for fitting")
        run_analysis(resolve_project_path(args.inputs), output, args.skip_docx)


if __name__ == "__main__":
    main()
