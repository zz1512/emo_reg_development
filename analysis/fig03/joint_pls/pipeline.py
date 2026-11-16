"""One nested LOSO pass, producing both joint and common-score predictions."""
from __future__ import annotations

from dataclasses import dataclass, replace
import numpy as np
import pandas as pd

from .models import MODELS, SCHEMES, prepare, fit_curve, reference_axis, metric, select_k, r2


@dataclass
class Data:
    x: np.ndarray
    y: np.ndarray
    cov: np.ndarray
    sites: np.ndarray
    ids: np.ndarray
    pfc: np.ndarray

    def take(self, indices):
        return Data(self.x[indices], self.y[indices], self.cov[indices], self.sites[indices],
                    self.ids[indices], self.pfc)

    def with_y(self, y):
        return replace(self, y=y)

    def validate(self):
        n = len(self.y)
        if (self.x.ndim != 2 or self.y.shape != (n, 2) or self.cov.shape != (n, 2)
                or self.x.shape[0] != n or len(self.sites) != n or len(self.ids) != n
                or self.pfc.shape != (self.x.shape[1],)):
            raise ValueError("Data arrays have inconsistent shapes")
        if not all(np.isfinite(a).all() for a in (self.x, self.y, self.cov)):
            raise ValueError("Input arrays contain nonfinite values")
        if len(np.unique(self.sites)) < 3 or not self.pfc.any() or self.pfc.all():
            raise ValueError("Need >=3 sites and nonempty PFC / non-PFC feature sets")


def _fit(train: Data, k_max):
    p = prepare(train.x, train.y, train.cov)
    masks = {"M1_PFC": p.keep & train.pfc,
             "M2_distributed432": p.keep,
             "M3_nonPFC": p.keep & ~train.pfc}
    curves = {m: fit_curve(p.x[:, mask], p.y, k_max, p.rank_c) for m, mask in masks.items()}
    q, scale, loadings = reference_axis(p.y, curves["M2_distributed432"])
    return p, masks, curves, q, scale, loadings


def _predictions(fitted, validation: Data, model, k):
    p, masks, curves, q, scale, _ = fitted
    x, y, base = p.transform(validation.x, validation.y, validation.cov)
    predicted = (np.zeros_like(y) if model == "M0_covariates"
                 else curves[model].predict(x[:, masks[model]], k))
    return y, predicted, y @ q / scale, predicted @ q / scale, base + predicted * p.scale_y


def component_tables(data: Data, k_max=10):
    """Descriptive first components of each full-sample feature-set fit."""
    fitted = _fit(data, k_max)
    p, masks, curves, _, _, _ = fitted
    summaries, weights = [], []
    for model in MODELS[1:]:
        curve, mask = curves[model], masks[model]
        q, scale, loadings = reference_axis(p.y, curve)
        flip = 1 if np.dot(q, curve.weights_y) >= 0 else -1
        w = curve.weights_x * flip
        t = (p.x[:, mask] - curve.mean_x) @ w
        summaries.append(dict(model=model, n=len(data.y), n_features=int(mask.sum()),
                              age_weight=q[0], success_weight=q[1], score_scale=scale,
                              age_loading=loadings[0], success_loading=loadings[1],
                              brain_score_age_correlation=float(np.corrcoef(t, p.y[:, 0])[0, 1]),
                              brain_score_success_correlation=float(np.corrcoef(t, p.y[:, 1])[0, 1]),
                              same_direction_weights=bool(q[0]*q[1] > 0),
                              interpretation="full-sample descriptive first component; not held-out performance or significance"))
        weights.extend(dict(model=model, parcel=int(j+1), first_component_weight=value)
                       for j, value in zip(np.flatnonzero(mask), w))
    return pd.DataFrame(summaries), pd.DataFrame(weights)


def nested_loso(data: Data, k_max=10, details=True):
    data.validate()
    sites = np.unique(data.sites)
    n = len(data.y)
    predictions = {(s, m): np.full((n, 2) if s == "joint" else n, np.nan)
                   for s in SCHEMES for m in MODELS}
    truth = {"joint": np.full((n, 2), np.nan), "shared": np.full(n, np.nan)}
    records, folds, inner_rows, axes = [], [], [], []
    for held in sites:
        train_idx, test_idx = np.flatnonzero(data.sites != held), np.flatnonzero(data.sites == held)
        train, test = data.take(train_idx), data.take(test_idx)
        inners = []
        for inner in np.unique(train.sites):
            a, b = np.flatnonzero(train.sites != inner), np.flatnonzero(train.sites == inner)
            fitted = _fit(train.take(a), k_max)
            inners.append((inner, fitted, train.take(b)))
        common_cap = min(curve.max_k for _, f, _ in inners for curve in f[2].values())
        chosen = {}
        for scheme in SCHEMES:
            for model in MODELS[1:]:
                losses = []
                for k in range(common_cap + 1):
                    sse, count = 0., 0
                    targets, estimates = [], []
                    for _, f, valid in inners:
                        yz, pred, u, up, _ = _predictions(f, valid, model, k)
                        target, estimate = (yz, pred) if scheme == "joint" else (u, up)
                        err = target - estimate
                        sse += float(np.sum(err ** 2))
                        count += err.size
                        targets.append(target)
                        estimates.append(estimate)
                    pooled_r2 = r2(np.concatenate(targets), np.concatenate(estimates))
                    losses.append(-pooled_r2)
                    if details:
                        inner_rows.append(dict(scheme=scheme, outer_site=held, model=model, k=k,
                                               standardized_mse=sse/count, pooled_r2=pooled_r2,
                                               selection_loss=losses[-1], n_values=count))
                chosen[scheme, model] = select_k(losses)
        fitted = _fit(train, k_max)
        p, masks, curves, q, scale, loadings = fitted
        if details:
            axes.append(dict(outer_site=held, reference="all432_first_Y_direction",
                             n_train=len(train.y), age_weight=q[0], success_weight=q[1],
                             score_scale=scale, age_loading=loadings[0], success_loading=loadings[1],
                             behavior_residual_correlation=float(np.corrcoef(p.y.T)[0, 1]),
                             same_direction_weights=bool(q[0] * q[1] > 0)))
        for scheme in SCHEMES:
            for model in MODELS:
                k = 0 if model == MODELS[0] else chosen[scheme, model]
                if model != MODELS[0] and k > curves[model].max_k:
                    raise ValueError("Outer training cannot fit the inner-selected component count")
                yz, pred, u, up, raw_pred = _predictions(fitted, test, model, k)
                target, prediction = (yz, pred) if scheme == "joint" else (u, up)
                truth[scheme][test_idx] = target
                predictions[scheme, model][test_idx] = prediction
                if details:
                    folds.append(dict(scheme=scheme, outer_site=held, model=model, k=k,
                                      common_k_max=common_cap, n_train=len(train.y), n_test=len(test.y),
                                      n_features=0 if model == MODELS[0] else int(masks[model].sum()),
                                      **metric(target, prediction)))
                    for j, row in enumerate(test_idx):
                        records.append(dict(scheme=scheme, model=model, sub_id=data.ids[row],
                                            site_id=held, k=k, age_true=data.y[row, 0],
                                            success_true=data.y[row, 1], age_z=yz[j, 0], success_z=yz[j, 1],
                                            age_pred_z=pred[j, 0], success_pred_z=pred[j, 1],
                                            age_pred=raw_pred[j, 0], success_pred=raw_pred[j, 1],
                                            shared_true=u[j], shared_pred=up[j]))
    performance = []
    for scheme in SCHEMES:
        for model in MODELS:
            prediction = predictions[scheme, model]
            if not np.isfinite(prediction).all():
                raise ValueError("Incomplete outer predictions")
            performance.append(dict(scheme=scheme, model=model, n=n, sites=len(sites),
                                    n_features=0 if model == MODELS[0] else
                                    int(data.pfc.sum()) if model == MODELS[1] else
                                    data.x.shape[1] if model == MODELS[2] else int((~data.pfc).sum()),
                                    **metric(truth[scheme], prediction)))
    perf = pd.DataFrame(performance)
    comparisons = comparison_table(perf)
    return dict(performance=perf, comparisons=comparisons, predictions=pd.DataFrame(records),
                folds=pd.DataFrame(folds), inner=pd.DataFrame(inner_rows), axes=pd.DataFrame(axes))


COMPARISONS = (("distributed_minus_pfc", "M2_distributed432", "M1_PFC"),
               ("pfc_minus_covariates", "M1_PFC", "M0_covariates"),
               ("distributed_minus_covariates", "M2_distributed432", "M0_covariates"),
               ("nonpfc_minus_covariates", "M3_nonPFC", "M0_covariates"),
               ("distributed_minus_nonpfc", "M2_distributed432", "M3_nonPFC"))


def comparison_table(performance):
    rows = []
    for scheme, frame in performance.groupby("scheme", sort=False):
        table = frame.set_index("model")
        for name, full, reduced in COMPARISONS:
            a, b = table.loc[full], table.loc[reduced]
            rows.append(dict(scheme=scheme, comparison=name, full_model=full, reduced_model=reduced,
                             n=int(a["n"]), delta_r2=a.r2-b.r2, delta_mae=a.mae-b.mae,
                             delta_rmse=a.rmse-b.rmse, p_permutation=np.nan, q_primary=np.nan,
                             inference="not_run"))
    return pd.DataFrame(rows)
