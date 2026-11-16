"""Auditable common cohort without any education eligibility requirement."""
from __future__ import annotations

import json
from pathlib import Path
import numpy as np
import pandas as pd

from analysis.fig03.predictive_pls.precheck import PARCELS, SCORE_MISSING, read_behavior, subject_keys
from analysis.fig03.predictive_pls.run_pls import load_membership
from .pipeline import Data
from .io import sha256, write_csv, own_output


def _cohort_from_sources(demographics, subjects, excluded, output):
    """Rebuild the original eligible roster before any education filter."""
    from analysis.fig03.gam_supplement.prepare_inputs import _read_id_list
    people = pd.read_csv(demographics, sep=None, engine="python", dtype=str, keep_default_na=False)
    if "site" in people and "site_id" not in people:
        people = people.rename(columns={"site": "site_id"})
    required = ["sub_id", "age", "gender", "site_id"]
    if set(required) - set(people):
        raise ValueError(f"Demographics missing {sorted(set(required)-set(people))}")
    people = people[required].copy()
    people["subject_key"] = subject_keys(people.sub_id, "Demographics")
    valid = pd.DataFrame({"sub_id": _read_id_list(Path(subjects))})
    valid["subject_key"] = subject_keys(valid.sub_id, "Valid subject list")
    removed = subject_keys(pd.Series(_read_id_list(Path(excluded))), "Second-level exclusions")
    audit = people.merge(valid[["subject_key"]].assign(in_valid_list=True),
                         on="subject_key", how="outer", validate="one_to_one")
    audit["in_valid_list"] = audit.in_valid_list.fillna(False).astype(bool)
    audit["second_level_removed"] = audit.subject_key.isin(removed)
    audit["upstream_exclusion"] = np.select(
        [~audit.in_valid_list, audit.second_level_removed],
        ["not_in_valid_subject_list", "second_level_removed"], default="")
    audit["sub_id"] = audit.sub_id.fillna(audit.subject_key.map(valid.set_index("subject_key").sub_id))
    write_csv(audit, output / "upstream_cohort_audit.csv")
    design = valid.merge(people.drop(columns="sub_id"), on="subject_key", how="left", validate="one_to_one")
    design["in_analysis"] = (~design.subject_key.isin(removed)).astype(int)
    design["excluded_by"] = np.where(design.in_analysis.eq(1), "", "second_level_removed")
    design["meanFD"] = np.nan  # Filled from the original motion table below.
    design = design.drop(columns="subject_key").fillna("")
    path = output / "cohort_without_education_filter.csv"
    write_csv(design, path)
    return path


def _recover(subjects, volume_dir, atlas_path, cached, output):
    """Match the existing effect-grid, finite-voxel mean algorithm exactly."""
    import nibabel as nib
    from nibabel.processing import resample_from_to
    names = ("effectSize-reg_neg_Observe_across_runs.nii.gz", "effectSize-pas_neg_Observe_across_runs.nii.gz")
    atlas = nib.load(str(atlas_path))
    labels_by_grid = {}
    provenance, rows = [], []

    def extract(subject):
        key = str(subject).removeprefix("sub-")
        folder = next((p for p in (volume_dir / str(subject), volume_dir / key,
                                  volume_dir / f"sub-{key}") if p.is_dir()), None)
        paths = [] if folder is None else [folder / name for name in names]
        if not paths or not all(p.is_file() for p in paths):
            return None, []
        a, b = [nib.load(str(p)) for p in paths]
        if a.shape != b.shape or not np.allclose(a.affine, b.affine):
            raise ValueError(f"Contrast image grids differ: {subject}")
        grid = (a.shape, tuple(a.affine.ravel()))
        if grid not in labels_by_grid:
            image = (atlas if atlas.shape == a.shape and np.allclose(atlas.affine, a.affine)
                     else resample_from_to(atlas, (a.shape, a.affine), order=0))
            labels_by_grid[grid] = np.rint(image.get_fdata()).astype(int).ravel()
        labels = labels_by_grid[grid]
        delta = a.get_fdata(dtype=np.float32).ravel() - b.get_fdata(dtype=np.float32).ravel()
        use = (labels > 0) & (labels <= 432) & np.isfinite(delta)
        counts = np.bincount(labels[use], minlength=433)[1:433]
        sums = np.bincount(labels[use], weights=delta[use], minlength=433)[1:433]
        values = np.divide(sums, counts, out=np.full(432, np.nan), where=counts > 0)
        return values, paths

    # A recovered row must not be mixed with a differently extracted cache.
    verified = False
    for cached_row in cached.itertuples(index=False, name=None):
        values, _ = extract(cached_row[0])
        if values is not None:
            if not np.allclose(values, np.asarray(cached_row[1:], float), rtol=1e-5, atol=1e-6, equal_nan=True):
                raise ValueError("Recovery algorithm/atlas does not reproduce cached parcel means; rebuild a consistent cache")
            verified = True
            break
    if not verified:
        raise ValueError("Cannot verify recovery against any cached subject with available original volumes")
    for subject in subjects:
        values, paths = extract(subject)
        provenance.append({"sub_id": subject, "status": "no_effect_volume" if values is None else "recovered",
                           "atlas_sha256": sha256(atlas_path),
                           "reappraise_path": str(paths[0]) if paths else "",
                           "looknegative_path": str(paths[1]) if paths else "",
                           "reappraise_sha256": sha256(paths[0]) if paths else "",
                           "looknegative_sha256": sha256(paths[1]) if paths else ""})
        if values is not None:
            rows.append([subject, *values])
    write_csv(pd.DataFrame(provenance), output / "brain_recovery_log.csv")
    recovered = pd.DataFrame(rows, columns=["sub_id", *PARCELS])
    write_csv(recovered, output / "recovered_brain.csv")
    return recovered


def precheck(inputs: Path, behavior: Path, membership: Path, output: Path,
             id_column="sub_ids", score_column="reappraisal_success", id_source=None,
             cohort_design=None, volume_dir=None, atlas=None, meanfd=None,
             demographics=None, subjects=None, excluded=None):
    inputs, output = inputs.resolve(), output.resolve()
    if output == inputs or output in inputs.parents or inputs in output.parents:
        raise ValueError("Precheck output must be separate from source inputs")
    own_output(output, "joint_pls_precheck_v1")
    upstream_sources = []
    if any(p is not None for p in (demographics, subjects, excluded)):
        if any(p is None for p in (demographics, subjects, excluded, meanfd)) or cohort_design:
            raise ValueError("Raw cohort requires demographics, subjects, excluded and meanfd together; do not also pass cohort-design")
        upstream_sources = [demographics, subjects, excluded]
        cohort_design = _cohort_from_sources(demographics, subjects, excluded, output)
    design_path = Path(cohort_design) if cohort_design else inputs / "design.csv"
    design = pd.read_csv(design_path, dtype=str, keep_default_na=False)
    required = {"sub_id", "age", "gender", "site_id", "meanFD"}
    if required - set(design):
        raise ValueError(f"Cohort design missing {sorted(required-set(design))}")
    if cohort_design is None and not {"in_analysis", "excluded_by"}.issubset(design):
        raise ValueError("Cache design lacks a full-cohort audit. Supply --cohort-design with a cohort unfiltered by education")
    note_path = inputs / "tian_combined.json"
    if int(json.loads(note_path.read_text())["n_cortical"]) != 400:
        raise ValueError("Only Schaefer400 + Tian32 (432) is supported")
    design["subject_key"] = subject_keys(design["sub_id"], "Cohort")
    design["site_id"] = design["site_id"].str.strip().str.replace(r"^(\d+)\.0$", r"\1", regex=True)
    for name in ("age", "gender", "meanFD"):
        design[name] = pd.to_numeric(design[name], errors="coerce")
    if meanfd:
        from analysis.fig03.gam_supplement.prepare_inputs import _read_meanfd_table
        fd = _read_meanfd_table(Path(meanfd))
        fd["subject_key"] = subject_keys(fd["sub_id"].astype(str), "Motion")
        design["meanFD"] = design["subject_key"].map(fd.set_index("subject_key")["mean_FD"])
    source_in = pd.to_numeric(design.get("in_analysis", pd.Series(1, index=design.index)), errors="coerce").eq(1)
    reason = design.get("excluded_by", pd.Series("", index=design.index)).str.strip()
    design["original_in_analysis"] = source_in
    design["original_excluded_by"] = reason
    education = reason.str.startswith(("parent_edu_", "sample_stop:parent_edu_"))
    # Availability/motion rules are re-evaluated; unknown QC exclusions remain binding.
    reevaluable = education | reason.str.match(r"^(meanFD_|no_effect_volume$)")
    design["independent_qc_failed"] = ~source_in & ~reevaluable
    design["education_exclusion_reconsidered"] = ~source_in & education
    behavior_rows = read_behavior(behavior, id_column, score_column, id_source)
    roster = design.merge(behavior_rows.drop(columns="source_sub_id"), on="subject_key", how="left",
                          validate="one_to_one", sort=False)
    rules = [("independent_original_qc", roster["independent_qc_failed"])]
    rules += [(f"{col}_missing_or_nonfinite", ~np.isfinite(roster[col])) for col in ("age", "gender", "meanFD")]
    rules += [("gender_not_0_or_1", ~roster.gender.isin([0, 1])),
              ("meanFD_above_0_5", roster.meanFD.gt(.5)),
              ("meanFD_negative", roster.meanFD.lt(0)),
              ("site_id_missing", roster.site_id.str.casefold().isin(SCORE_MISSING)),
              ("behavior_subject_not_found", roster.source_row.isna()),
              ("success_missing_or_invalid", ~np.isfinite(roster.success))]
    eligible = ~np.column_stack([np.asarray(failed, bool) for _, failed in rules]).any(1)
    brain_path = inputs / "parcel_response.csv"
    brain = pd.read_csv(brain_path, dtype={"sub_id": str}, keep_default_na=False)
    if set(brain) != {"sub_id", *PARCELS}:
        raise ValueError("Response table must contain exactly sub_id and parcels 1..432")
    brain = pd.concat([brain[["sub_id"]], brain[PARCELS].apply(pd.to_numeric, errors="coerce")], axis=1).copy()
    brain["subject_key"] = subject_keys(brain.sub_id, "Brain cache")
    missing = roster.loc[eligible & ~roster.subject_key.isin(brain.subject_key), "sub_id"].tolist()
    if missing:
        write_csv(roster.loc[roster.sub_id.isin(missing)], output / "requires_brain_recovery.csv")
        if volume_dir is None:
            write_csv(roster, output / "sample_roster.csv")
            raise ValueError(f"{len(missing)} eligible subjects lack cached means. Supply --volume-dir and --atlas; education exclusions must not silently persist")
        recovered = _recover(missing, Path(volume_dir), Path(atlas) if atlas else inputs / "combined_atlas.nii.gz",
                             brain[["sub_id", *PARCELS]], output)
        recovered["subject_key"] = subject_keys(recovered.sub_id, "Recovered brain") if len(recovered) else pd.Series(dtype=str)
        brain = pd.concat([brain, recovered], ignore_index=True)
    roster = roster.merge(brain.drop(columns="sub_id"), on="subject_key", how="left", validate="one_to_one", indicator="brain_match").copy()
    rules += [("brain_subject_not_found", roster.brain_match.ne("both")),
              ("brain_has_nonfinite_parcels", ~np.isfinite(roster[PARCELS].to_numpy(float)).all(1))]
    first, all_reasons = np.full(len(roster), "", object), [[] for _ in range(len(roster))]
    summaries, exclusions = [], []
    for name, failed in rules:
        failed = np.asarray(failed, bool)
        before = first == ""
        new = before & failed
        first[new] = name
        summaries.append(dict(filter=name, n_failed_total=int(failed.sum()), n_before=int(before.sum()),
                              n_newly_excluded=int(new.sum()), n_remaining=int((first == "").sum())))
        for index in np.flatnonzero(failed):
            all_reasons[index].append(name)
            exclusions.append(dict(sub_id=roster.at[index, "sub_id"], site_id=roster.at[index, "site_id"],
                                   filter=name, first_exclusion=bool(new[index])))
    roster["first_exclusion"], roster["all_exclusions"] = first, [";".join(v) for v in all_reasons]
    roster["in_joint_analysis"] = first == ""
    sample = roster.loc[roster.in_joint_analysis].copy().reset_index(drop=True)
    write_csv(roster.drop(columns=PARCELS), output / "sample_roster.csv")
    write_csv(pd.DataFrame(summaries), output / "filter_summary.csv")
    write_csv(pd.DataFrame(exclusions, columns=["sub_id", "site_id", "filter", "first_exclusion"]), output / "exclusions_by_filter.csv")
    alignment = behavior_rows.merge(roster[["subject_key", "sub_id", "in_joint_analysis", "all_exclusions"]],
                                    on="subject_key", how="left", validate="one_to_one")
    write_csv(alignment, output / "behavior_alignment.csv")
    write_csv(sample.drop(columns=PARCELS), output / "common_sample.csv")
    write_csv(sample.groupby("site_id").agg(n=("sub_id", "size"), age_mean=("age", "mean"),
              age_min=("age", "min"), age_max=("age", "max"), success_mean=("success", "mean"),
              success_sd=("success", "std")).reset_index(), output / "site_summary.csv")
    members = load_membership(membership)
    if members.definition_source.iloc[0] != "Donahue2018_liberal_bilateral":
        raise ValueError("This implementation is the broad anatomical PFC main analysis")
    data = Data(sample[PARCELS].to_numpy(float), sample[["age", "success"]].to_numpy(float),
                sample[["gender", "meanFD"]].to_numpy(float), sample.site_id.to_numpy(str),
                sample.sub_id.to_numpy(str), members.is_pfc.to_numpy(bool))
    data.validate()
    np.savez_compressed(output / "data.npz", x=data.x, y=data.y, cov=data.cov, sites=data.sites, ids=data.ids, pfc=data.pfc)
    sources = [design_path, brain_path, note_path, inputs / "parcel_info.csv", behavior, membership]
    sources += [p for p in (id_source, meanfd) if p is not None] + upstream_sources
    metadata = dict(kind="joint_pls_common_sample_v1", n=len(sample), sites=len(np.unique(data.sites)),
                    parent_education_used=False, age_is_outcome=True, score_column=score_column,
                    score_direction="source score; no sign reversal; verify source formula before interpretation",
                    reference_axis="training-only full432 PLS first Y direction; brain-informed, not a pure behavioral factor",
                    restored_education_excluded_n=int(sample.education_exclusion_reconsidered.sum()),
                    sources={str(Path(p).resolve()): sha256(p) for p in sources})
    (output / "sample_metadata.json").write_text(json.dumps(metadata, indent=2, ensure_ascii=False), encoding="utf-8")
    return metadata


def load_data(path):
    with np.load(Path(path) / "data.npz", allow_pickle=False) as f:
        return Data(**{key: f[key] for key in ("x", "y", "cov", "sites", "ids", "pfc")})
