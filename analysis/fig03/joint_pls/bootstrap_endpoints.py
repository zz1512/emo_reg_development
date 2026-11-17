"""Recover absolute model R² for archived paired bootstrap draws.

Only M0's evaluation targets are rebuilt. No candidate model, inner CV or
permutation is rerun. Shared targets require the outer-training reference's
first component, with exactly the original preprocessing and solver.
"""
from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor
import json
from pathlib import Path

import numpy as np
import pandas as pd

from .comparison_summary import BASELINE, _records
from .io import MODELS, SCHEMES, own_output, sha256, write_csv


def baseline_r2(data):
    from .models import prepare, fit_curve, reference_axis, r2
    truth = {"joint": np.empty((len(data.y), 2)), "shared": np.empty(len(data.y))}
    for held in np.unique(data.sites):
        train = data.take(np.flatnonzero(data.sites != held))
        indices = np.flatnonzero(data.sites == held)
        test = data.take(indices)
        p = prepare(train.x, train.y, train.cov)
        # Later components cannot change the reference's first Y direction.
        curve = fit_curve(p.x[:, p.keep], p.y, 1, p.rank_c)
        q, scale, _ = reference_axis(p.y, curve)
        _, yz, _ = p.transform(test.x, test.y, test.cov)
        truth["joint"][indices], truth["shared"][indices] = yz, yz @ q / scale
    return {s: r2(y, np.zeros_like(y)) for s, y in truth.items()}


_STATE = None


def _initialize(root, cache):
    from .sample import load_data
    global _STATE
    root = Path(root)
    _STATE = root, Path(cache), load_data(root / "precheck")


def _recover(replicate):
    from threadpoolctl import threadpool_limits
    root, cache, data = _STATE
    source = root / "bootstrap/independent" / f"rep_{replicate:06d}.csv"
    index_path = source.with_name(f"indices_{replicate:06d}.npy")
    hashes = dict(source_sha256=sha256(source), indices_sha256=sha256(index_path))
    destination = cache / source.name
    if destination.is_file():
        previous = pd.read_csv(destination)
        if (len(previous) == 8 and all(previous[k].eq(v).all() for k, v in hashes.items()) and
                np.isfinite(previous.r2).all() and not previous.duplicated(["scheme", "model"]).any()):
            return replicate, "cached"
    indices = np.load(index_path, allow_pickle=False)
    if (indices.shape != (len(data.y),) or not np.issubdtype(indices.dtype, np.integer) or
            np.any(indices < 0) or np.any(indices >= len(data.y))):
        raise ValueError(f"Invalid stored bootstrap indices: {index_path}")
    draw = data.take(indices)
    for site in np.unique(data.sites):
        if np.sum(draw.sites == site) != np.sum(data.sites == site):
            raise ValueError("Stored bootstrap changes the original within-site sample sizes")
    table = pd.read_csv(source)
    if not table.status.eq("ok").all():
        raise ValueError(f"Cannot recover a failed original replicate: {source}")
    with threadpool_limits(limits=1):
        baseline = baseline_r2(draw)
    rows = []
    for scheme in SCHEMES:
        part = table.loc[table.scheme.eq(scheme)].set_index("comparison")
        for model in MODELS:
            delta = 0. if model == MODELS[0] else float(part.at[BASELINE[model], "delta_r2"])
            rows.append(dict(scheme=scheme, model=model, replicate=replicate,
                             r2=baseline[scheme]+delta, m0_r2=baseline[scheme], delta_from_m0=delta, **hashes))
    frame = pd.DataFrame(rows)
    if not np.isfinite(frame.r2).all():
        raise ValueError("Recovered nonfinite model R²")
    temporary = destination.with_suffix(".tmp")
    write_csv(frame, temporary)
    temporary.replace(destination)
    return replicate, "ok"


def recover(root, workers=1):
    from .sample import load_data
    from threadpoolctl import threadpool_limits
    root = Path(root).resolve()
    if workers < 1:
        raise ValueError("workers must be positive")
    original = json.loads((root / "analysis_manifest.json").read_text(encoding="utf-8"))
    current = {name: sha256(Path(__file__).with_name(name)) for name in original["code_sha256"]}
    if current != original["code_sha256"] or sha256(root / "precheck/data.npz") != original["data_sha256"]:
        raise ValueError("Frozen fitting code or sample changed; cannot reconstruct old reference targets")
    _, quality = _records(root, Path("bootstrap/independent"), "bootstrap", "independent")
    if not quality["complete"] or quality["valid"] < 2:
        raise ValueError("Endpoint CIs require every requested original bootstrap repeat to be valid (at least two)")
    bootstrap_manifest = root / "bootstrap/independent/manifest.json"
    config = dict(method="rebuild_M0_outer_targets_then_add_original_paired_differences_v1",
                  data_sha256=original["data_sha256"], code_sha256=current,
                  recovery_code_sha256=sha256(Path(__file__)),
                  bootstrap_manifest_sha256=sha256(bootstrap_manifest), total=quality["requested"])
    cache = root / "report_cache/model_r2_bootstrap"
    own_output(cache, "joint_pls_model_r2_recovery_v1")
    manifest = cache / "manifest.json"
    if manifest.is_file() and json.loads(manifest.read_text(encoding="utf-8")) != config:
        raise ValueError("Endpoint recovery configuration changed; preserve the old report_cache and use a fresh recovery cache")
    manifest.write_text(json.dumps(config, indent=2), encoding="utf-8")
    # Before using any repeat, verify reconstruction against the observed M0.
    with threadpool_limits(limits=1):
        baseline = baseline_r2(load_data(root / "precheck"))
    observed = pd.read_csv(root / "model_performance.csv").set_index(["scheme", "model"])
    for scheme in SCHEMES:
        np.testing.assert_allclose(baseline[scheme], observed.at[(scheme, MODELS[0]), "r2"], atol=1e-10, rtol=1e-9,
                                   err_msg="Recovered reference evaluation does not reproduce observed M0")
    total = quality["requested"]
    if workers == 1:
        _initialize(root, cache)
        iterator = map(_recover, range(total))
        for i, status in iterator:
            if (i+1) % 25 == 0 or i+1 == total:
                print(f"Recover bootstrap model R2: {i+1}/{total} {status}", flush=True)
    else:
        with ProcessPoolExecutor(workers, initializer=_initialize, initargs=(root, cache)) as pool:
            for i, status in pool.map(_recover, range(total)):
                if (i+1) % 25 == 0 or i+1 == total:
                    print(f"Recover bootstrap model R2: {i+1}/{total} {status}", flush=True)
    records = pd.concat([pd.read_csv(cache / f"rep_{i:06d}.csv") for i in range(total)], ignore_index=True)
    rows = []
    for (scheme, model), part in records.groupby(["scheme", "model"], sort=False):
        low, high = np.percentile(part.r2, [2.5, 97.5])
        rows.append(dict(scheme=scheme, model=model, ci_low=low, ci_high=high,
                         bootstrap_valid=part.replicate.nunique(), bootstrap_requested=total,
                         complete=True, data_sha256=original["data_sha256"],
                         interpretation="Absolute model R2 CI from the same original paired bootstrap draws"))
    output = root / "report"
    own_output(output, "joint_pls_report_v1")
    write_csv(records, output / "paired_bootstrap_model_r2.csv")
    write_csv(pd.DataFrame(rows), output / "bootstrap_model_r2_intervals.csv")
    provenance = dict(config, observed_m0_verified=True, baseline_r2=baseline,
                      intervals_sha256=sha256(output / "bootstrap_model_r2_intervals.csv"),
                      note="No inner tuning or candidate model refitting; first outer reference direction only")
    (output / "bootstrap_model_r2_provenance.json").write_text(json.dumps(provenance, indent=2), encoding="utf-8")
    return output / "bootstrap_model_r2_intervals.csv"


def read_intervals(root):
    """Do not silently draw stale endpoint intervals on a newer analysis."""
    root = Path(root)
    path = root / "report/bootstrap_model_r2_intervals.csv"
    if not path.is_file():
        return None
    provenance = json.loads((path.parent / "bootstrap_model_r2_provenance.json").read_text(encoding="utf-8"))
    original = json.loads((root / "analysis_manifest.json").read_text(encoding="utf-8"))
    if (provenance["data_sha256"] != original["data_sha256"] or
            provenance["code_sha256"] != original["code_sha256"] or
            provenance["bootstrap_manifest_sha256"] != sha256(root / "bootstrap/independent/manifest.json") or
            provenance["intervals_sha256"] != sha256(path)):
        raise ValueError("Endpoint intervals no longer match the frozen analysis and bootstrap configuration")
    return pd.read_csv(path)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--workers", type=int, default=1)
    args = parser.parse_args()
    print(recover(args.output, args.workers))
