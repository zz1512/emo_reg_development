"""CLI: precheck, fit both routes, permutation, bootstrap, report, or mock."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import time
import platform
import numpy as np
import pandas as pd
import sklearn
from threadpoolctl import threadpool_limits

from analysis.fig03.predictive_pls.precheck import project_path
from .sample import precheck, load_data, own_output, sha256, write_csv
from .pipeline import nested_loso, component_tables
from .models import PLS_MAX_ITER, PLS_TOL
from .report import report
from . import resampling


def fit(output, k_max=10, simulated=False):
    output = Path(output).resolve()
    own_output(output, "joint_pls_analysis_v1")
    if k_max < 1:
        raise ValueError("k-max must be positive; k=0 is automatically included as a baseline option")
    manifest = dict(kind="joint_pls_both_routes_v1", k_max=k_max, simulated=bool(simulated),
                    data_sha256=sha256(output / "precheck" / "data.npz"),
                    code_sha256={name: sha256(Path(__file__).with_name(name))
                                 for name in ("models.py", "pipeline.py", "resampling.py", "io.py")},
                    versions=dict(python=platform.python_version(), numpy=np.__version__,
                                  pandas=pd.__version__, sklearn=sklearn.__version__),
                    reference="training-only all432 PLS first Y direction, common to all candidate models",
                    solver=dict(estimator="sklearn.PLSRegression", scale=False,
                                max_iter=PLS_MAX_ITER, tol=PLS_TOL,
                                convergence_warning="failure; no unconverged result accepted"),
                    prediction_covariates=["intercept", "gender", "meanFD"],
                    primary_family=["joint:distributed_minus_pfc", "shared:distributed_minus_pfc"],
                    default_permutations=1000, default_bootstraps=1000, resampling_nested=False)
    manifest_path = output / "analysis_manifest.json"
    if manifest_path.is_file() and json.loads(manifest_path.read_text(encoding="utf-8")) != manifest:
        raise ValueError("Sample/configuration/code differs from existing analysis. Use a new output directory")
    started = time.perf_counter()
    with threadpool_limits(limits=1):
        data = load_data(output / "precheck")
        result = nested_loso(data, k_max)
        component_summary, brain_weights = component_tables(data, k_max)
    write_csv(component_summary, output / "model_first_components.csv")
    write_csv(brain_weights, output / "model_first_component_weights.csv")
    for key, name in (("performance", "model_performance"), ("comparisons", "observed_comparisons"),
                      ("predictions", "oof_predictions"), ("folds", "outer_folds"),
                      ("inner", "inner_scores"), ("axes", "reference_axes")):
        write_csv(result[key], output / f"{name}.csv")
    manifest_path.write_text(json.dumps(manifest, ensure_ascii=False, indent=2), encoding="utf-8")
    resampling.merge(output)
    timing = dict(observed_fit_seconds=time.perf_counter()-started,
                  note="One full both-route nested LOSO pass; resampling wall time must be benchmarked on the server")
    (output / "timing.json").write_text(json.dumps(timing, indent=2), encoding="utf-8")
    print(f"Both routes fitted: {output}; N={len(data.y)}, PFC={int(data.pfc.sum())}, "
          f"remaining={int((~data.pfc).sum())}; {timing['observed_fit_seconds']:.2f}s", flush=True)
    return result


def mock(output, membership, n=96, seed=5, k_max=2, permutations=3, bootstraps=3, docx=True):
    from analysis.fig03.predictive_pls.run_pls import load_membership
    import nibabel as nib
    output = output.resolve()
    own_output(output, "joint_pls_mock_v1")
    source, analysis = output / "source", output / "analysis"
    source.mkdir(exist_ok=True)
    membership = membership.resolve()
    member = load_membership(membership)
    rng = np.random.default_rng(seed)
    sites = np.asarray([str(i % 4 + 1) for i in range(n)])
    ids = np.asarray([f"s{i:05d}" for i in range(n)])
    gender, fd = np.arange(n) % 2, rng.uniform(.03, .4, n)
    common = rng.normal(size=n)
    x = rng.normal(size=(n, 432))
    pfc = member.is_pfc.to_numpy(bool)
    x[:, pfc] += .15 * common[:, None]
    x[:, ~pfc] += .5 * common[:, None]
    age = 12 + 2*common + rng.normal(size=n)
    success = .35*common + .1*gender + .2*fd + rng.normal(scale=.5, size=n)
    design = pd.DataFrame(dict(sub_id=ids, age=age, gender=gender, site_id=sites, meanFD=fd,
                               parent_edu=np.nan, in_analysis=1, excluded_by=""))
    design.loc[:1, ["in_analysis", "excluded_by"]] = [0, "parent_edu_missing_value"]
    write_csv(design, source / "design.csv")
    brain = pd.DataFrame(x, columns=[str(i) for i in range(1, 433)])
    brain.insert(0, "sub_id", ids)
    write_csv(brain.iloc[2:], source / "parcel_response.csv")
    write_csv(member[["parcel", "name"]], source / "parcel_info.csv")
    (source / "tian_combined.json").write_text(json.dumps(dict(n_cortical=400, simulated=True)), encoding="utf-8")
    # Wrong behavior IDs + independently supplied rowwise real IDs exercise the actual mapping route.
    write_csv(pd.DataFrame(dict(sub_ids=[f"wrong{i}" for i in range(n)], reappraisal_success=success)), source / "behavior.csv")
    (source / "row_ids.csv").write_text("subject,other\n" + "\n".join(f"{s},unused" for s in ids) + "\n", encoding="utf-8")
    shape = (8, 9, 6)
    atlas = np.arange(1, 433).reshape(shape).astype(np.int16)
    nib.save(nib.Nifti1Image(atlas, np.eye(4)), str(source / "combined_atlas.nii.gz"))
    for index in (0, 1, 2):
        folder = source / "volumes" / ids[index]
        folder.mkdir(parents=True, exist_ok=True)
        for name, value in (("reg_neg", x[index].reshape(shape)), ("pas_neg", np.zeros(shape))):
            nib.save(nib.Nifti1Image(value.astype(np.float32), np.eye(4)),
                     str(folder / f"effectSize-{name}_Observe_across_runs.nii.gz"))
    own_output(analysis, "joint_pls_analysis_v1")
    precheck(source, source / "behavior.csv", membership, analysis / "precheck",
             id_source=source / "row_ids.csv", volume_dir=source / "volumes")
    fit(analysis, k_max, simulated=True)
    if permutations:
        resampling.run(analysis, total=permutations)
    if bootstraps:
        resampling.run(analysis, kind="bootstrap", total=bootstraps)
    report(analysis, docx=docx)
    return analysis


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    p = sub.add_parser("precheck")
    for name in ("inputs", "behavior", "membership", "output"):
        p.add_argument(f"--{name}", type=Path, required=True)
    for name in ("cohort-design", "id-source", "volume-dir", "atlas", "meanfd", "demographics", "subjects", "excluded"):
        p.add_argument(f"--{name}", type=Path)
    p.add_argument("--id-column", default="sub_ids")
    p.add_argument("--score-column", default="reappraisal_success")
    p = sub.add_parser("fit")
    p.add_argument("--output", type=Path, required=True)
    p.add_argument("--k-max", type=int, default=10)
    for name in ("permutation", "bootstrap"):
        p = sub.add_parser(name)
        p.add_argument("--output", type=Path, required=True)
        p.add_argument("--total", type=int, default=1000)
        p.add_argument("--start", type=int, default=0)
        p.add_argument("--count", type=int)
        p.add_argument("--workers", type=int, default=1)
        p.add_argument("--seed", type=int, default=20261117)
        if name == "permutation":
            p.add_argument("--null", choices=["pfc", "covariates", "nonpfc"], default="pfc")
    p = sub.add_parser("report")
    p.add_argument("--output", type=Path, required=True)
    p.add_argument("--skip-docx", action="store_true")
    p.add_argument("--preview", action="store_true")
    p.add_argument("--model-ci", action="store_true", help="Show optional recovered single-model R2 intervals")
    p = sub.add_parser("mock")
    p.add_argument("--output", type=Path, default=Path("outputs/fig03/joint_pls_mock"))
    p.add_argument("--membership", type=Path, default=Path("outputs/fig03/predictive_pls_432_FD05/pfc_definition/pfc_parcel_membership.csv"))
    p.add_argument("--n", type=int, default=96)
    p.add_argument("--k-max", type=int, default=2)
    p.add_argument("--permutations", type=int, default=3)
    p.add_argument("--bootstraps", type=int, default=3)
    p.add_argument("--skip-docx", action="store_true")
    args = parser.parse_args(argv)
    for name, value in vars(args).items():
        if isinstance(value, Path):
            setattr(args, name, project_path(value))
    if args.command == "precheck":
        own_output(args.output, "joint_pls_analysis_v1")
        print(precheck(args.inputs, args.behavior, args.membership, args.output / "precheck",
                       args.id_column, args.score_column, args.id_source, args.cohort_design,
                       args.volume_dir, args.atlas, args.meanfd, args.demographics,
                       args.subjects, args.excluded))
    elif args.command == "fit":
        fit(args.output, args.k_max)
    elif args.command in ("permutation", "bootstrap"):
        resampling.run(args.output, kind=args.command, null=getattr(args, "null", "pfc"),
                       total=args.total, start=args.start, count=args.count, workers=args.workers, seed=args.seed)
    elif args.command == "report":
        # Results-only recovery preserves the frozen fits and repeat archives.
        print(report(args.output, not args.skip_docx, preview=args.preview, model_ci=args.model_ci))
    elif args.command == "mock":
        print(mock(args.output, args.membership, args.n, k_max=args.k_max,
                   permutations=args.permutations, bootstraps=args.bootstraps, docx=not args.skip_docx))


if __name__ == "__main__":
    main()
