"""Site-level bootstrap that refits the nested LOSO on each replicate.

A site drawn more than once is repeated as extra rows with its original site
id, so copies stay in the same train or test fold. Replicate draws depend on
the master seed and the replicate number, which makes chunks restartable.
"""

from __future__ import annotations

import argparse
import json
import os
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd

from analysis.fig03.predictive_pls.models import (
    MODELS,
    PRIMARY_COMPARISON,
    SUPPORT_COMPARISONS,
    PredictionFailure,
    comparison_rows,
    point_performance,
    stored_text,
)
from analysis.fig03.predictive_pls.precheck import project_path
from analysis.fig03.predictive_pls.run_pls import load_config, nested_loso

VALID_FRACTION_MIN = 0.95
VALID_N_MIN = 1000
REQUESTED_B = 2000


def site_draw(seed: int, replicate: int, site_order: np.ndarray) -> np.ndarray:
    generator = np.random.default_rng(np.random.SeedSequence([int(seed), int(replicate)]))
    order = np.asarray(site_order).astype(str)
    return generator.choice(order, size=order.size, replace=True)


def expand_index(sites: np.ndarray, draw: np.ndarray) -> np.ndarray:
    sites = np.asarray(sites).astype(str)
    pieces = [np.flatnonzero(sites == site) for site in draw]
    if any(piece.size == 0 for piece in pieces):
        raise PredictionFailure("Bootstrap draw names a site that is absent from the frozen sample")
    return np.concatenate(pieces)


def confirmatory_flags(n_valid: int, n_requested: int, lower: float, task: str, comparison: str) -> dict:
    fraction = n_valid / n_requested if n_requested else 0.0
    inference_valid = bool(n_requested >= REQUESTED_B and n_valid >= VALID_N_MIN and fraction >= VALID_FRACTION_MIN)
    supported = bool(inference_valid and task == "success" and comparison == PRIMARY_COMPARISON
                     and np.isfinite(lower) and lower > 0)
    return {"n_bootstrap_requested": n_requested, "n_bootstrap_valid": n_valid,
            "bootstrap_fraction_valid": fraction, "inference_valid": inference_valid, "supported": supported}


def _load_task(analysis: Path, task: str) -> dict:
    archive = np.load(analysis / "frozen" / f"{task}_arrays.npz", allow_pickle=True)
    sites = stored_text(archive["sites"])
    mask = archive["pfc_mask"].astype(bool)
    matrix = archive["X"]
    parcels = archive["parcel_ids"].astype(int)
    features = {
        "M1_PFC": (matrix[:, mask], parcels[mask]),
        "M2_distributed432": (matrix, parcels),
        "M3_nonPFC": (matrix[:, ~mask], parcels[~mask]),
    }
    raw = {key: archive[key] for key in ("age", "gender", "parent_edu", "meanFD")}
    sub_id = stored_text(archive["sub_id"])
    return {"raw": raw, "features": features, "y": archive["y"], "sites": sites,
            "sub_id": sub_id, "site_order": pd.unique(sites)}


def evaluate_replicate(task_data: dict, task: str, seed: int, replicate: int) -> dict:
    draw = site_draw(seed, replicate, task_data["site_order"])
    index = expand_index(task_data["sites"], draw)
    counts = pd.Series(draw).value_counts()
    record = {"task": task, "replicate": int(replicate),
              "site_counts": ";".join(f"{site}:{int(counts.get(site, 0))}" for site in task_data["site_order"]),
              "n_unique_sites": int(pd.unique(task_data["sites"][index]).size), "valid": False, "reason": ""}
    y = task_data["y"][index]
    if record["n_unique_sites"] < 3:
        record["reason"] = "fewer_than_3_sites"
        return record
    if np.std(y) <= 1e-12:
        record["reason"] = "target_without_variation"
        return record
    raw = {key: values[index] for key, values in task_data["raw"].items()}
    features = {model: (matrix[index], parcels) for model, (matrix, parcels) in task_data["features"].items()}
    try:
        predictions, _, _ = nested_loso(raw, features, y, task_data["sites"][index], task,
                                        task_data["sub_id"][index], allow_repeated_rows=True)
    except PredictionFailure as error:
        record["reason"] = str(error)
        return record
    performance = point_performance(predictions).set_index("model")
    for model in MODELS:
        record[f"r2_{model}"] = performance.at[model, "r2"]
        record[f"mae_{model}"] = performance.at[model, "mae"]
        record[f"rmse_{model}"] = performance.at[model, "rmse"]
    for name in (PRIMARY_COMPARISON, *SUPPORT_COMPARISONS):
        left, right = name.split("-")
        record[f"delta_r2_{name}"] = performance.at[left, "r2"] - performance.at[right, "r2"]
    record["valid"] = True
    return record


_WORKER: dict = {}


def _worker_init(analysis: str, seed: int) -> None:
    os.environ["OMP_NUM_THREADS"] = "1"
    os.environ["MKL_NUM_THREADS"] = "1"
    os.environ["OPENBLAS_NUM_THREADS"] = "1"
    root = Path(analysis)
    _WORKER["seed"] = seed
    _WORKER["data"] = {task: _load_task(root, task) for task in ("success", "age")}


def _worker_job(item: tuple[int, str]) -> dict:
    from threadpoolctl import threadpool_limits

    replicate, task = item
    with threadpool_limits(limits=1):
        return evaluate_replicate(_WORKER["data"][task], task, _WORKER["seed"], replicate)


def write_chunk(analysis: Path, start: int, count: int, seed: int, jobs: int = 1) -> Path:
    analysis = analysis.resolve()
    chunk_dir = analysis / "bootstrap_chunks"
    chunk_dir.mkdir(exist_ok=True)
    items = [(replicate, task) for replicate in range(start, start + count) for task in ("success", "age")]
    if jobs <= 1:
        loaded = {task: _load_task(analysis, task) for task in ("success", "age")}
        rows = [evaluate_replicate(loaded[task], task, seed, replicate) for replicate, task in items]
    else:
        with ProcessPoolExecutor(max_workers=jobs, initializer=_worker_init,
                                 initargs=(str(analysis), seed)) as pool:
            rows = []
            for row in pool.map(_worker_job, items, chunksize=1):
                rows.append(row)
                if len(rows) % 20 == 0:
                    print(f"Bootstrap finished {len(rows)} of {len(items)} task-replicates", flush=True)
    rows = sorted(rows, key=lambda row: (row["task"], row["replicate"]))
    path = chunk_dir / f"replicates_{start:05d}_{start + count:05d}.csv"
    pd.DataFrame(rows).to_csv(path, index=False)
    print(f"Bootstrap chunk: {path}", flush=True)
    return path


def merge_chunks(analysis: Path, requested: int) -> pd.DataFrame:
    paths = sorted((analysis / "bootstrap_chunks").glob("replicates_*.csv"))
    if not paths:
        raise PredictionFailure("No bootstrap chunks are available to merge")
    table = pd.concat([pd.read_csv(path) for path in paths], ignore_index=True)
    if table.duplicated(["task", "replicate"]).any():
        raise PredictionFailure("Bootstrap chunks overlap on task and replicate")
    table = table.sort_values(["task", "replicate"]).reset_index(drop=True)
    table.to_csv(analysis / "bootstrap_metrics.csv", index=False)
    _write_intervals(analysis, table, requested)
    return table


def _write_intervals(analysis: Path, replicates: pd.DataFrame, requested: int) -> None:
    predictions = pd.read_csv(analysis / "oof_predictions.csv")
    performance = point_performance(predictions)
    feature_path = analysis / "model_performance.csv"
    if feature_path.is_file():
        previous = pd.read_csv(feature_path)
        performance["n_features"] = previous["n_features"].to_numpy()
    intervals = {}
    for task, frame in replicates.groupby("task"):
        valid = frame.loc[frame["valid"].astype(str).str.lower().isin(("true", "1"))]
        fraction = len(valid) / requested if requested else 0.0
        for model in MODELS:
            column = f"r2_{model}"
            low, high = (np.nan, np.nan) if valid.empty else np.nanpercentile(valid[column], [2.5, 97.5])
            matched = performance["task"].eq(task) & performance["model"].eq(model)
            performance.loc[matched, ["ci_low_r2", "ci_high_r2", "n_bootstrap_valid", "bootstrap_fraction_valid"]] = (
                low, high, len(valid), fraction)
        for name in (PRIMARY_COMPARISON, *SUPPORT_COMPARISONS):
            column = f"delta_r2_{name}"
            low, high = (np.nan, np.nan) if valid.empty or column not in valid else np.nanpercentile(valid[column], [2.5, 97.5])
            flags = confirmatory_flags(len(valid), requested, low, task, name)
            intervals[(task, name)] = {"ci_low_delta_r2": low, "ci_high_delta_r2": high, **flags}
    comparisons = comparison_rows(predictions, intervals)
    performance.to_csv(analysis / "model_performance.csv", index=False)
    comparisons.to_csv(analysis / "model_comparisons.csv", index=False)
    status = {"requested_replicates": requested,
              "valid_by_task": {task: int(frame["valid"].astype(str).str.lower().isin(("true", "1")).sum())
                                for task, frame in replicates.groupby("task")}}
    (analysis / "frozen" / "bootstrap_status.json").write_text(json.dumps(status, indent=2), encoding="utf-8")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--analysis", type=Path, required=True)
    parser.add_argument("--start", type=int, default=0)
    parser.add_argument("--count", type=int, default=REQUESTED_B)
    parser.add_argument("--seed", type=int, default=None)
    parser.add_argument("--merge-only", action="store_true")
    parser.add_argument("--jobs", type=int, default=os.cpu_count() or 1)
    args = parser.parse_args()
    config = load_config()
    seed = config["bootstrap"]["seed"] if args.seed is None else args.seed
    requested = config["bootstrap"]["B"]
    analysis = project_path(args.analysis)
    if not args.merge_only:
        if args.start < 0 or args.count <= 0 or args.start + args.count > requested:
            raise PredictionFailure(f"Bootstrap chunk must lie inside 0..{requested}")
        write_chunk(analysis, args.start, args.count, seed, max(1, args.jobs))
    present = 0
    chunk_dir = analysis / "bootstrap_chunks"
    if chunk_dir.is_dir():
        present = sum(pd.read_csv(path)["replicate"].nunique() for path in chunk_dir.glob("replicates_*.csv"))
    if present >= requested:
        merge_chunks(analysis, requested)
        print(f"Bootstrap metrics: {analysis / 'bootstrap_metrics.csv'}", flush=True)


if __name__ == "__main__":
    main()
