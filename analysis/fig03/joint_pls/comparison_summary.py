"""Results-only PFC comparisons, including paired recovery of M3 minus M1.

No fitting code is imported or changed. Historical repeat archives contain
M3-M0 and M1-M0 on the same resample, so their difference is exactly M3-M1.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

from .io import SCHEMES, MODELS, write_csv


COMPARISONS = (
    ("pfc_minus_covariates", "M1_PFC", "M0_covariates", "exploratory"),
    ("distributed_minus_covariates", "M2_distributed432", "M0_covariates", "exploratory"),
    ("nonpfc_minus_covariates", "M3_nonPFC", "M0_covariates", "exploratory"),
    ("distributed_minus_pfc", "M2_distributed432", "M1_PFC", "primary"),
    ("nonpfc_minus_pfc", "M3_nonPFC", "M1_PFC", "secondary"),
)
METRICS = ("delta_r2", "delta_mae", "delta_rmse")
BH_FAMILIES = (("primary", "q_primary", 2), ("secondary", "q_secondary", 2),
               ("exploratory", "q_exploratory", 6))
BASELINE = {"M1_PFC": "pfc_minus_covariates", "M2_distributed432": "distributed_minus_covariates",
            "M3_nonPFC": "nonpfc_minus_covariates"}


def _bh(values):
    p = np.asarray(values, float)
    order = np.argsort(p)
    adjusted = np.minimum.accumulate((p[order] * len(p) / np.arange(1, len(p)+1))[::-1])[::-1]
    result = np.empty_like(p)
    result[order] = np.minimum(adjusted, 1.)
    return result


def paired_comparisons(table):
    """Derive from paired differences; never subtract P values or interval ends."""
    if table.duplicated(["scheme", "comparison"]).any():
        raise ValueError("Duplicate scheme/comparison in a repeat archive")
    rows = []
    for scheme in SCHEMES:
        part = table.loc[table.scheme.eq(scheme)].set_index("comparison")
        if not set(BASELINE.values()).issubset(part.index):
            raise ValueError(f"Missing paired model-minus-baseline statistics for {scheme}")
        for name, full, reduced, family in COMPARISONS:
            a = part.loc[BASELINE[full], list(METRICS)].to_numpy(float)
            b = (np.zeros(3) if reduced == MODELS[0] else
                 part.loc[BASELINE[reduced], list(METRICS)].to_numpy(float))
            difference = a - b
            if name in part.index:
                saved = part.loc[name, list(METRICS)].to_numpy(float)
                np.testing.assert_allclose(difference, saved,
                                           atol=1e-10, rtol=1e-9,
                                           err_msg="Archived paired differences are inconsistent")
                difference = saved
            if not np.isfinite(difference).all():
                raise ValueError("Nonfinite paired comparison statistics")
            rows.append(dict(scheme=scheme, comparison=name, full_model=full, reduced_model=reduced,
                             family=family, **dict(zip(METRICS, difference))))
    return pd.DataFrame(rows)


def _records(root, relative, kind, null):
    folder = root / relative
    path = folder / "manifest.json"
    empty = pd.DataFrame(columns=["scheme", "comparison", "replicate", *METRICS])
    if not path.is_file():
        return empty, dict(kind=kind, null=null, requested=0, completed=0, valid=0, complete=False)
    meta = json.loads(path.read_text(encoding="utf-8"))
    original = json.loads((root / "analysis_manifest.json").read_text(encoding="utf-8"))
    for key in ("data_sha256", "code_sha256", "k_max"):
        if key in original and meta.get(key) != original[key]:
            raise ValueError(f"Resampling {key} differs from the frozen observed analysis: {folder}")
    if meta.get("kind", kind) != kind or meta.get("null", null) != null:
        raise ValueError(f"Unexpected resampling kind/null in {path}")
    total = int(meta["total"])
    if total < 1:
        raise ValueError("Requested repeat count must be positive")
    frames, completed, failures = [], 0, []
    for i in range(total):
        file = folder / f"rep_{i:06d}.csv"
        if not file.is_file():
            continue
        completed += 1
        table = pd.read_csv(file)
        if table.empty or not table.replicate.eq(i).all():
            raise ValueError(f"Empty or mismatched replicate ID: {file}")
        if not table.status.eq("ok").all():
            failures.append(dict(replicate=i, reason="; ".join(table.get("reason", pd.Series(dtype=str)).dropna().astype(str))))
            continue
        paired = paired_comparisons(table)
        paired["replicate"] = i
        frames.append(paired)
    valid = len(frames)
    quality = dict(kind=kind, null=null, requested=total, completed=completed, valid=valid,
                   complete=completed == valid == total, failures=failures)
    return pd.concat(frames, ignore_index=True) if frames else empty, quality


def summarize(root: Path, output: Path):
    """Refresh report outputs from immutable observed fits and raw repeat files."""
    root, output = Path(root), Path(output)
    performance = pd.read_csv(root / "model_performance.csv")
    if (len(performance) != 8 or performance.duplicated(["scheme", "model"]).any() or
            set(performance.scheme) != set(SCHEMES) or
            not np.isfinite(performance[["r2", "mae", "rmse"]].to_numpy(float)).all()):
        raise ValueError("Expected eight unique finite scheme/model performance records")
    if performance["n"].nunique() != 1 or performance.sites.nunique() != 1:
        raise ValueError("All model comparisons require the same participants and held-out sites")
    differences = []
    for scheme in SCHEMES:
        p = performance.loc[performance.scheme.eq(scheme)].set_index("model")
        for name, full, reduced, family in COMPARISONS:
            differences.append(dict(scheme=scheme, comparison=name, full_model=full, reduced_model=reduced,
                                    n=int(p.at[full, "n"]), sites=int(p.at[full, "sites"]), family=family,
                                    **{f"delta_{m}": p.at[full, m]-p.at[reduced, m] for m in ("r2", "mae", "rmse")}))
    comparisons = pd.DataFrame(differences)
    for column in ("p_permutation", "q_primary", "q_secondary", "q_exploratory", "q_bh",
                   "null_mean_delta_r2", "ci_low", "ci_high"):
        comparisons[column] = np.nan
    comparisons["bh_family_size"] = 0
    comparisons["inference"] = "not_run"
    comparisons["null"] = comparisons.family.map({"primary": "pfc", "secondary": "pfc", "exploratory": "covariates"})
    comparisons["permutations_valid"] = comparisons["permutations_requested"] = 0
    comparisons["bootstrap_valid"] = comparisons["bootstrap_requested"] = 0
    comparisons["bootstrap_complete"] = False
    repeats, quality = [], []
    for null in ("pfc", "covariates"):
        frame, status = _records(root, Path("permutations") / null, "permutation", null)
        quality.append(status)
        frame["null"] = null
        repeats.append(frame)
        for index, row in comparisons.loc[comparisons["null"].eq(null)].iterrows():
            values = frame.loc[frame.scheme.eq(row.scheme) & frame.comparison.eq(row.comparison), "delta_r2"]
            comparisons.loc[index, "permutations_valid"] = len(values)
            comparisons.loc[index, "permutations_requested"] = status["requested"]
            comparisons.loc[index, "null_mean_delta_r2"] = values.mean()
            complete = status["complete"] and len(values) == status["requested"]
            if complete:
                tolerance = 100 * np.finfo(float).eps * max(1., abs(row.delta_r2))
                comparisons.loc[index, "p_permutation"] = (1 + (values >= row.delta_r2-tolerance).sum()) / (len(values)+1)
            comparisons.loc[index, "inference"] = ("complete_approximate_reduced_model_permutation" if complete else
                                                    "incomplete" if status["requested"] else "not_run")
    # Preserve the two original PFC comparison families. The three baseline
    # contrasts across both routes form one additional exploratory family.
    for family, column, size in BH_FAMILIES:
        selected = comparisons.family.eq(family)
        if int(selected.sum()) != size:
            raise ValueError(f"Expected {size} tests in the {family} BH family")
        comparisons.loc[selected, "bh_family_size"] = size
        if comparisons.loc[selected, "p_permutation"].notna().all():
            comparisons.loc[selected, column] = _bh(comparisons.loc[selected, "p_permutation"])
            comparisons.loc[selected, "q_bh"] = comparisons.loc[selected, column]
    for family, column, _ in BH_FAMILIES:
        comparisons[f"supported_{family}"] = pd.Series(pd.NA, index=comparisons.index, dtype="boolean")
        ready = comparisons.family.eq(family) & comparisons[column].notna()
        comparisons.loc[ready, f"supported_{family}"] = comparisons.loc[ready, column].lt(.05) & comparisons.loc[ready, "delta_r2"].gt(0)
    boots, status = _records(root, Path("bootstrap") / "independent", "bootstrap", "independent")
    quality.append(status)
    comparisons["bootstrap_requested"], comparisons["bootstrap_valid"] = status["requested"], status["valid"]
    comparisons["bootstrap_complete"] = status["complete"]
    if status["complete"] and status["valid"] >= 2:
        for index, row in comparisons.iterrows():
            values = boots.loc[boots.scheme.eq(row.scheme) & boots.comparison.eq(row.comparison), "delta_r2"]
            comparisons.loc[index, ["ci_low", "ci_high"]] = np.percentile(values, [2.5, 97.5])
    output.mkdir(parents=True, exist_ok=True)
    write_csv(comparisons, output / "model_comparisons_full_precision.csv")
    write_csv(comparisons.loc[comparisons.family.eq("primary")], output / "primary_comparisons.csv")
    write_csv(comparisons.loc[comparisons.family.isin(["primary", "secondary"])], output / "pfc_comparisons.csv")
    perm = pd.concat(repeats, ignore_index=True)
    write_csv(perm, output / "paired_permutation_statistics.csv")
    write_csv(boots, output / "paired_bootstrap_statistics.csv")
    provenance = dict(derivation="Within each original replicate, all five pairwise differences derive from the three model-minus-M0 values",
                      primary_family="joint/shared distributed_minus_pfc, BH over 2 tests (unchanged)",
                      secondary_family="joint/shared nonpfc_minus_pfc, BH over 2 secondary tests",
                      exploratory_family="Three model-minus-M0 contrasts across joint/shared, BH over 6 exploratory tests",
                      q_bh="Within-family adjusted P, copied from q_primary, q_secondary or q_exploratory",
                      correction_scope="Three separate BH families of 2, 2 and 6 tests; not a combined ten-test correction",
                      omitted="M3-M2 / M2-M3 is not an analysis target and is omitted from report tables and figures",
                      null="Approximate PFC-reduced conditional null; not generic equality of nonnested algorithms",
                      quality=quality)
    (output / "comparison_provenance.json").write_text(json.dumps(provenance, indent=2), encoding="utf-8")
    return comparisons, perm
