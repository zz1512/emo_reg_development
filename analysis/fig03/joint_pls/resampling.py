"""Independent bootstrap and conditional reduced-model permutation branches.

Freedman-Lane-style permutation is approximate: the reduced conditional mean
is linear and residual *rows* must be exchangeable within sites. It tests
incremental conditional information, not equality of two arbitrary algorithms.
"""
from __future__ import annotations

from concurrent.futures import ProcessPoolExecutor
import json
from pathlib import Path
import numpy as np
import pandas as pd
from threadpoolctl import threadpool_limits

from .models import MODELS, covariates
from .pipeline import nested_loso, _fit
from .sample import load_data, sha256, write_csv


def reduced_null(data, null):
    c, _, _ = covariates(data.cov)
    sites = np.unique(data.sites)
    c = np.column_stack([c, *[(data.sites == s).astype(float) for s in sites[1:]]])
    if null != "covariates":
        mask = data.pfc if null == "pfc" else ~data.pfc if null == "nonpfc" else None
        if mask is None:
            raise ValueError("Null must be pfc, nonpfc or covariates")
        x = data.x[:, mask]
        scale = x.std(0)
        keep = scale > 1e-12
        c = np.column_stack([c, (x[:, keep] - x[:, keep].mean(0)) / scale[keep]])
    beta, _, rank, _ = np.linalg.lstsq(c, data.y, rcond=None)
    fitted = c @ beta
    residual = data.y - fitted
    df = len(c) - int(rank)
    if df < 5 or np.any(residual.std(0) < 1e-12):
        raise ValueError(f"Reduced null has insufficient residual information: df={df}")
    return fitted, residual, dict(reduced_rank=int(rank), reduced_residual_df=df,
                                 null=null, nuisance_sites=list(map(str, sites)),
                                 assumption="linear reduced conditional mean; bivariate residual rows exchangeable within site; approximate Freedman-Lane")


def within_site_permutation(sites, seed, replicate):
    rng = np.random.default_rng(np.random.SeedSequence([int(seed), int(replicate), 811]))
    index = np.arange(len(sites))
    for site in np.unique(sites):
        rows = np.flatnonzero(sites == site)
        index[rows] = rng.permutation(rows)
    return index


def bootstrap_index(sites, seed, replicate):
    rng = np.random.default_rng(np.random.SeedSequence([int(seed), int(replicate), 977]))
    return np.concatenate([rng.choice(rows, size=len(rows), replace=True)
                           for site in np.unique(sites)
                           for rows in [np.flatnonzero(sites == site)]])


def full_patterns(data, k_max):
    fitted = _fit(data, max(1, k_max))
    _, masks, curves, q, _, loadings = fitted
    weights = np.zeros((3, data.x.shape[1]))
    for j, name in enumerate(MODELS[1:]):
        curve = curves[name]
        if curve.weights_x is not None:
            flip = 1 if curve.weights_y[0] >= 0 else -1
            weights[j, masks[name]] = curve.weights_x * flip
    return weights, np.asarray(q), np.asarray(loadings)


_STATE = None


def _initialize(root, kind, null, seed):
    global _STATE
    root = Path(root)
    config = json.loads((root / "analysis_manifest.json").read_text(encoding="utf-8"))
    data = load_data(root / "precheck")
    null_fit = reduced_null(data, null) if kind == "permutation" else None
    pattern = full_patterns(data, config["k_max"])[0] if kind == "bootstrap" else None
    _STATE = root, data, config["k_max"], kind, null, seed, null_fit, pattern


def _replicate(replicate):
    root, data, k_max, kind, null, seed, null_fit, pattern = _STATE
    folder = root / ("permutations" if kind == "permutation" else "bootstrap") / null
    path = folder / f"rep_{replicate:06d}.csv"
    if path.is_file():
        return replicate, "cached"
    try:
        with threadpool_limits(limits=1):
            if kind == "permutation":
                index = within_site_permutation(data.sites, seed, replicate)
                current = data.with_y(null_fit[0] + null_fit[1][index])
            else:
                index = bootstrap_index(data.sites, seed, replicate)
                current = data.take(index)
            result = nested_loso(current, k_max, details=False)
            table = result["comparisons"].copy()
            table["replicate"], table["status"], table["reason"] = replicate, "ok", ""
            if not np.isfinite(table.delta_r2).all():
                raise ValueError("Nonfinite comparison statistic")
            if kind == "bootstrap":
                weights, q, loadings = full_patterns(current, k_max)
                for j in range(3):
                    if np.dot(weights[j], pattern[j]) < 0:
                        weights[j] *= -1
                np.savez_compressed(folder / f"weights_{replicate:06d}.npz", weights=weights, q=q, loadings=loadings)
            # Save exact row indices for traceable paired resampling.
            np.save(folder / f"indices_{replicate:06d}.npy", index, allow_pickle=False)
    except Exception as error:
        table = pd.DataFrame([dict(replicate=replicate, status="failed", reason=f"{type(error).__name__}: {error}")])
    temporary = path.with_suffix(".tmp")
    write_csv(table, temporary)
    temporary.replace(path)
    return replicate, str(table.status.iloc[0])


def run(root, kind="permutation", null="pfc", total=1000, start=0, count=None, workers=1, seed=20261117):
    root = Path(root).resolve()
    manifest = json.loads((root / "analysis_manifest.json").read_text(encoding="utf-8"))
    if total < 1 or workers < 1 or start < 0:
        raise ValueError("total/workers must be positive; start must be nonnegative")
    code_hashes = {name: sha256(Path(__file__).with_name(name))
                   for name in ("models.py", "pipeline.py", "resampling.py", "io.py")}
    if code_hashes != manifest["code_sha256"] or sha256(root / "precheck" / "data.npz") != manifest["data_sha256"]:
        raise ValueError("Code or sample changed since observed fit. Use a new output and refit before resampling")
    if kind == "bootstrap":
        null = "independent"
    data = load_data(root / "precheck")
    metadata = dict(kind=kind, null=null, total=total, seed=seed, data_sha256=manifest["data_sha256"],
                    code_sha256=code_hashes, k_max=manifest["k_max"], bootstrap_nested=False)
    if kind == "permutation":
        _, _, description = reduced_null(data, null)
        metadata.update(description)
    folder = root / ("permutations" if kind == "permutation" else "bootstrap") / null
    folder.mkdir(parents=True, exist_ok=True)
    saved = folder / "manifest.json"
    if saved.is_file() and json.loads(saved.read_text(encoding="utf-8")) != metadata:
        raise ValueError("Resampling parameters differ from existing records; use a new analysis output")
    saved.write_text(json.dumps(metadata, ensure_ascii=False, indent=2), encoding="utf-8")
    stop = total if count is None else min(total, start + count)
    jobs = [i for i in range(start, stop) if not (folder / f"rep_{i:06d}.csv").is_file()]
    if workers == 1:
        _initialize(str(root), kind, null, seed)
        iterator = map(_replicate, jobs)
        for i, status in iterator:
            print(f"{kind} {null}: {i+1}/{total} {status}", flush=True)
    else:
        with ProcessPoolExecutor(workers, initializer=_initialize, initargs=(str(root), kind, null, seed)) as pool:
            for i, status in pool.map(_replicate, jobs):
                print(f"{kind} {null}: {i+1}/{total} {status}", flush=True)
    return merge(root)


def bh(p):
    p = np.asarray(p, float)
    order = np.argsort(p)
    adjusted = np.minimum.accumulate((p[order] * len(p) / np.arange(1, len(p)+1))[::-1])[::-1]
    result = np.empty_like(p)
    result[order] = np.minimum(adjusted, 1.)
    return result


NULL_COMPARISONS = {"pfc": ["distributed_minus_pfc"], "nonpfc": ["distributed_minus_nonpfc"],
                    "covariates": ["pfc_minus_covariates", "distributed_minus_covariates", "nonpfc_minus_covariates"]}


def merge(root):
    root = Path(root)
    observed = pd.read_csv(root / "observed_comparisons.csv")
    merged = observed.copy()
    quality = []
    for manifest_path in sorted((root / "permutations").glob("*/manifest.json")):
        meta = json.loads(manifest_path.read_text(encoding="utf-8"))
        folder = manifest_path.parent
        frames = [pd.read_csv(folder / f"rep_{i:06d}.csv") for i in range(meta["total"])
                  if (folder / f"rep_{i:06d}.csv").is_file()]
        if not frames:
            continue
        all_rows = pd.concat(frames, ignore_index=True)
        write_csv(all_rows, folder / "replicates.csv")
        successful = all_rows.loc[all_rows.status.eq("ok")]
        n_done = all_rows.replicate.nunique()
        n_valid = successful.replicate.nunique()
        # A failed replicate is not silently removed to obtain a favorable p.
        complete = n_done == meta["total"] and n_valid == meta["total"]
        quality.append(dict(kind="permutation", null=meta["null"], requested=meta["total"],
                            completed=n_done, valid=n_valid, inference_complete=complete))
        if successful.empty:
            continue
        for index, row in merged.iterrows():
            if row.comparison not in NULL_COMPARISONS[meta["null"]]:
                continue
            values = successful.loc[successful.scheme.eq(row.scheme) & successful.comparison.eq(row.comparison), "delta_r2"]
            if len(values) != n_valid:
                raise ValueError("Duplicate or incomplete scheme/comparison records")
            merged.loc[index, "null"] = meta["null"]
            merged.loc[index, "permutations_valid"] = n_valid
            merged.loc[index, "permutations_requested"] = meta["total"]
            merged.loc[index, "null_mean_delta_r2"] = values.mean()
            merged.loc[index, "p_permutation"] = ((1 + (values >= row.delta_r2).sum()) / (meta["total"] + 1)
                                                  if complete else np.nan)
            merged.loc[index, "inference"] = "complete_approximate_reduced_model_permutation" if complete else "incomplete"
    primary = merged.comparison.eq("distributed_minus_pfc") & merged.p_permutation.notna()
    # The two prespecified analysis routes form the primary family of two tests.
    if primary.sum() == 2:
        merged.loc[primary, "q_primary"] = bh(merged.loc[primary, "p_permutation"])
    merged["supported_primary"] = pd.Series(pd.NA, index=merged.index, dtype="boolean")
    merged.loc[primary, "supported_primary"] = (merged.loc[primary, "q_primary"].lt(.05) &
                                                merged.loc[primary, "delta_r2"].gt(0))
    write_csv(merged, root / "model_comparisons.csv")
    for manifest_path in (root / "bootstrap").glob("*/manifest.json"):
        meta = json.loads(manifest_path.read_text(encoding="utf-8"))
        folder = manifest_path.parent
        frames = [pd.read_csv(folder / f"rep_{i:06d}.csv") for i in range(meta["total"])
                  if (folder / f"rep_{i:06d}.csv").is_file()]
        if not frames:
            continue
        table = pd.concat(frames, ignore_index=True)
        write_csv(table, folder / "replicates.csv")
        valid = table.loc[table.status.eq("ok")]
        n_valid = valid.replicate.nunique()
        complete = len(frames) == meta["total"] and n_valid == meta["total"]
        quality.append(dict(kind="bootstrap", null="independent", requested=meta["total"],
                            completed=len(frames), valid=n_valid, inference_complete=complete))
        if valid.empty:
            continue
        interval_rows = []
        for (scheme, comparison), group in valid.groupby(["scheme", "comparison"]):
            low, high = np.percentile(group.delta_r2, [2.5, 97.5])
            interval_rows.append(dict(scheme=scheme, comparison=comparison, ci_low=low, ci_high=high,
                                      n_bootstrap_valid=n_valid, complete=complete,
                                      use="stability only; no permutation p inferred"))
        write_csv(pd.DataFrame(interval_rows), root / "bootstrap_intervals.csv")
        if complete:
            _merge_patterns(root, folder, valid.replicate.unique(), meta["k_max"])
    write_csv(pd.DataFrame(quality), root / "resampling_status.csv")
    return merged


def _merge_patterns(root, folder, replicates, k_max):
    data = load_data(root / "precheck")
    original, q, loadings = full_patterns(data, k_max)
    paths = [folder / f"weights_{int(i):06d}.npz" for i in replicates]
    if not all(p.is_file() for p in paths):
        raise ValueError("A successful bootstrap is missing its weight archive")
    weights, qs, ls = [], [], []
    for path in paths:
        with np.load(path) as a:
            weights.append(a["weights"])
            qs.append(a["q"])
            ls.append(a["loadings"])
    weights, qs, ls = np.stack(weights), np.stack(qs), np.stack(ls)
    se = weights.std(0, ddof=1)
    ratios = np.divide(original, se, out=np.full_like(original, np.nan), where=se > 1e-12)
    rows = []
    for m, name in enumerate(MODELS[1:]):
        mask = data.pfc if m == 0 else np.ones(len(data.pfc), bool) if m == 1 else ~data.pfc
        for parcel in np.flatnonzero(mask):
            rows.append(dict(model=name, parcel=parcel+1, first_component_weight=original[m, parcel],
                             bootstrap_se=se[m, parcel], bootstrap_ratio=ratios[m, parcel],
                             n_valid=len(replicates), interpretation="first full-sample component stability; not regionwise corrected significance"))
    write_csv(pd.DataFrame(rows), root / "brain_first_component_stability.csv")
    axes = []
    for j, name in enumerate(("age", "success")):
        low, high = np.percentile(ls[:, j], [2.5, 97.5])
        axes.append(dict(variable=name, weight=q[j], loading=loadings[j], loading_ci_low=low,
                         loading_ci_high=high, weight_bootstrap_se=qs[:, j].std(ddof=1), n_valid=len(replicates)))
    write_csv(pd.DataFrame(axes), root / "reference_axis_stability.csv")
