"""Train-fold covariate adjustment and PLS prediction.

Every coder, residual regression, variance screen, and scaler is fit on the
training rows of that call. k=0 returns the covariate prediction alone.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from sklearn.cross_decomposition import PLSRegression
from sklearn.preprocessing import SplineTransformer

K_MAX = 10
TIE_TOLERANCE = 1e-10
NEAR_ZERO_STD = 1e-12
MODELS = ("M0_covariates", "M1_PFC", "M2_distributed432", "M3_nonPFC")
NEURAL_MODELS = MODELS[1:]
PRIMARY_COMPARISON = "M2_distributed432-M1_PFC"
SUPPORT_COMPARISONS = (
    "M1_PFC-M0_covariates",
    "M2_distributed432-M0_covariates",
    "M3_nonPFC-M0_covariates",
    "M2_distributed432-M3_nonPFC",
)


class PredictionFailure(RuntimeError):
    """A fold or replicate cannot be completed under the locked rules."""


def stored_text(values) -> np.ndarray:
    """Fixed-width text. Object arrays of IDs are pickled by np.savez and then fail to load."""
    items = ["" if value is None else str(value) for value in np.asarray(values, dtype=object).reshape(-1)]
    width = max((len(item) for item in items), default=1)
    return np.asarray(items, dtype=f"U{width}")


@dataclass
class CovariateCoder:
    task: str
    spline: SplineTransformer | None = None
    continuous_center: dict = field(default_factory=dict)
    continuous_scale: dict = field(default_factory=dict)
    spline_knots: np.ndarray | None = None

    def fit(self, age, gender, parent_edu, meanfd):
        age = np.asarray(age, dtype=float)
        if self.task == "success":
            self.spline = SplineTransformer(
                n_knots=3, degree=3, knots="quantile", include_bias=False, extrapolation="linear")
            self.spline.fit(age.reshape(-1, 1))
            self.spline_knots = np.asarray(self.spline.bsplines_[0].t, dtype=float)
        for name, values in (("parent_edu", parent_edu), ("meanFD", meanfd)):
            values = np.asarray(values, dtype=float)
            center = float(np.mean(values))
            scale = float(np.std(values))
            self.continuous_center[name] = center
            self.continuous_scale[name] = scale
        self._check_gender(gender)
        return self

    @staticmethod
    def _check_gender(gender):
        gender = np.asarray(gender, dtype=float)
        if not np.isin(gender, (0.0, 1.0)).all():
            raise PredictionFailure("Gender must be coded 0 or 1 before covariate fitting")

    def transform(self, age, gender, parent_edu, meanfd) -> np.ndarray:
        age = np.asarray(age, dtype=float)
        gender = np.asarray(gender, dtype=float).reshape(-1, 1)
        columns = [np.ones((age.size, 1))]
        if self.task == "success":
            columns.append(self.spline.transform(age.reshape(-1, 1)))
        columns.append(gender)
        for name, values in (("parent_edu", parent_edu), ("meanFD", meanfd)):
            values = np.asarray(values, dtype=float)
            scale = self.continuous_scale[name]
            centered = values - self.continuous_center[name]
            columns.append((np.zeros(age.size) if not np.isfinite(scale) or scale <= NEAR_ZERO_STD
                            else centered / scale).reshape(-1, 1))
        design = np.column_stack(columns)
        if self.task == "success" and design.shape[1] != 8:
            raise PredictionFailure(f"Success covariate design has {design.shape[1]} columns, expected 8")
        if self.task == "age" and design.shape[1] != 4:
            raise PredictionFailure(f"Age covariate design has {design.shape[1]} columns, expected 4")
        return design


@dataclass
class FoldState:
    coder: CovariateCoder
    a: np.ndarray
    B: np.ndarray | None
    keep: np.ndarray
    x_mean: np.ndarray
    x_scale: np.ndarray
    parcel_ids: np.ndarray
    pls: dict
    rank_c: int
    rank_x: int
    n_train: int
    max_k: int
    limit_name: str
    converged: dict = field(default_factory=dict)
    coef_by_k: dict = field(default_factory=dict)
    pls_x_mean: np.ndarray | None = None
    pls_intercept: float = 0.0

    @property
    def dropped_parcels(self) -> list[int]:
        if self.parcel_ids.size == 0:
            return []
        return [int(parcel) for parcel in self.parcel_ids[~self.keep]]


def _lstsq(design: np.ndarray, target: np.ndarray):
    coef, _, rank, _ = np.linalg.lstsq(design, target, rcond=None)
    return coef, int(rank)


def _component_limit(n_train, rank_c, n_features, rank_x) -> tuple[int, str]:
    caps = {
        "grid": K_MAX,
        "n_features": int(n_features),
        "rank_X": int(rank_x),
        "residual_df": int(max(0, n_train - rank_c - 1)),
    }
    name = min(caps, key=caps.get)
    return int(caps[name]), name


def fit_fold(age, gender, parent_edu, meanfd, y, task: str,
             X: np.ndarray | None = None, parcel_ids: np.ndarray | None = None,
             components_through: int | None = None) -> FoldState:
    """Fit one training fold. X=None builds the covariate-only state."""
    y = np.asarray(y, dtype=float).reshape(-1)
    coder = CovariateCoder(task).fit(age, gender, parent_edu, meanfd)
    design = coder.transform(age, gender, parent_edu, meanfd)
    if design.shape[0] != y.size:
        raise PredictionFailure("Covariate rows and target rows differ")
    a, rank_c = _lstsq(design, y)
    y_residual = y - design @ a
    if X is not None:
        X = np.asarray(X, dtype=float)
        parcel_ids = np.asarray(parcel_ids, dtype=int)
        if X.ndim != 2 or X.shape[0] != y.size or X.shape[1] != parcel_ids.size:
            raise PredictionFailure("Neural feature matrix does not match the training rows or parcel ids")
        if not np.isfinite(X).all():
            raise PredictionFailure("Training neural features contain nonfinite values")
    if np.std(y_residual) <= NEAR_ZERO_STD:
        empty = np.zeros(0 if X is None else X.shape[1], dtype=bool)
        parcels = np.zeros(0, dtype=int) if parcel_ids is None else np.asarray(parcel_ids, dtype=int)
        return FoldState(coder, a, None, empty, np.zeros(0), np.zeros(0), parcels, {},
                         rank_c, 0, y.size, 0, "target_residual_constant")
    if X is None:
        limit, limit_name = 0, "no_neural_features"
        return FoldState(coder, a, None, np.zeros(0, dtype=bool), np.zeros(0), np.zeros(0),
                         np.zeros(0, dtype=int), {}, rank_c, 0, y.size, limit, limit_name)
    X = np.asarray(X, dtype=float)
    parcel_ids = np.asarray(parcel_ids, dtype=int)
    if X.ndim != 2 or X.shape[0] != y.size or X.shape[1] != parcel_ids.size:
        raise PredictionFailure("Neural feature matrix does not match the training rows or parcel ids")
    if not np.isfinite(X).all():
        raise PredictionFailure("Training neural features contain nonfinite values")
    B, _ = _lstsq(design, X)
    residual_X = X - design @ B
    scale = residual_X.std(axis=0)
    keep = scale > NEAR_ZERO_STD
    if not keep.any():
        limit, limit_name = 0, "all_features_near_zero_variance"
        return FoldState(coder, a, B, keep, np.zeros(0), np.zeros(0), parcel_ids, {},
                         rank_c, 0, y.size, limit, limit_name)
    kept = residual_X[:, keep]
    center = kept.mean(axis=0)
    kept_scale = kept.std(axis=0)
    standardized = (kept - center) / kept_scale
    rank_x = int(np.linalg.matrix_rank(standardized))
    limit, limit_name = _component_limit(y.size, rank_c, standardized.shape[1], rank_x)
    through = limit if components_through is None else min(limit, int(components_through))
    if through == 0:
        return FoldState(coder, a, B, keep, center, kept_scale, parcel_ids, {},
                         rank_c, rank_x, y.size, limit, limit_name)
    try:
        pls = PLSRegression(n_components=through, scale=False)
        pls.fit(standardized, y_residual)
    except np.linalg.LinAlgError as error:
        raise PredictionFailure(f"PLS failed at k={through}: {error}") from error
    obtained = len(pls.n_iter_)
    coefs = _prefix_coefs(pls, obtained)
    converged = {component: bool(pls.n_iter_[component - 1] < pls.max_iter) for component in range(1, obtained + 1)}
    for component in range(obtained + 1, through + 1):
        # A stopped fit pads unused components with zeros, so those requested widths are fit alone.
        one = PLSRegression(n_components=component, scale=False)
        try:
            one.fit(standardized, y_residual)
        except np.linalg.LinAlgError as error:
            raise PredictionFailure(f"PLS failed at k={component}: {error}") from error
        coefs[component] = np.asarray(one.coef_, dtype=float).reshape(-1)
        iterations = np.asarray(getattr(one, "n_iter_", [0]), dtype=float)
        converged[component] = bool(np.max(iterations) < one.max_iter)
    state = FoldState(coder, a, B, keep, center, kept_scale, parcel_ids, {},
                      rank_c, rank_x, y.size, limit, limit_name, converged)
    state.coef_by_k = coefs
    state.pls_x_mean = np.asarray(pls._x_mean, dtype=float).reshape(-1)
    state.pls_intercept = float(np.asarray(pls.intercept_).reshape(-1)[0])
    return state


def _prefix_coefs(pls: PLSRegression, k_max: int) -> dict[int, np.ndarray]:
    """Coefficients for the first k components of one completed PLS fit.

    sklearn builds rotations from the components that were requested, so each
    prefix is inverted on its own. That matches a model fit with n_components=k
    when every requested component was actually computed.
    """
    from scipy.linalg import pinv

    coefs = {}
    for component in range(1, k_max + 1):
        weights = pls.x_weights_[:, :component]
        loadings = pls.x_loadings_[:, :component]
        y_loadings = pls.y_loadings_[:, :component]
        rotation = weights @ pinv(loadings.T @ weights, check_finite=False)
        coef = rotation @ y_loadings.T
        coef = (coef * pls._y_std).T / pls._x_std
        coefs[component] = np.asarray(coef, dtype=float).reshape(-1)
    return coefs


def predict_fold(state: FoldState, age, gender, parent_edu, meanfd,
                 X: np.ndarray | None, k: int) -> np.ndarray:
    """Apply a stored training fold. The test target is not an input."""
    design = state.coder.transform(age, gender, parent_edu, meanfd)
    baseline = design @ state.a
    if int(k) == 0:
        return baseline
    if int(k) not in state.coef_by_k and int(k) not in state.pls:
        raise PredictionFailure(f"Requested k={k} was not fit on this training fold")
    if X is None:
        raise PredictionFailure("Neural prediction requires the feature matrix")
    residual_X = np.asarray(X, dtype=float) - design @ state.B
    standardized = (residual_X[:, state.keep] - state.x_mean) / state.x_scale
    if int(k) in state.coef_by_k:
        centered = standardized - state.pls_x_mean
        neural = centered @ state.coef_by_k[int(k)] + state.pls_intercept
        return baseline + neural
    neural = np.asarray(state.pls[int(k)].predict(standardized), dtype=float).reshape(-1)
    return baseline + neural


def select_k(mse_by_k: list[float]) -> int:
    """Smallest k whose pooled MSE is within the absolute tie tolerance of the best."""
    scores = np.asarray(mse_by_k, dtype=float)
    if scores.ndim != 1 or not np.isfinite(scores).all():
        raise PredictionFailure("Component selection requires one finite MSE per k")
    best = float(np.min(scores))
    tied = np.flatnonzero(scores <= best + TIE_TOLERANCE)
    return int(tied[0])


def point_performance(predictions) -> "pd.DataFrame":
    import pandas as pd

    rows = []
    for task, frame in predictions.groupby("task", sort=True):
        y = frame["y_true"].to_numpy(dtype=float)
        for model in MODELS:
            metrics = regression_metrics(y, frame[f"pred_{model}"])
            rows.append({
                "task": task, "model": model, "n": metrics["n"],
                "n_sites": int(frame["site_id"].nunique()),
                "n_features": np.nan,
                "median_k": float(np.median(frame[f"k_{model}"])),
                "r2": metrics["r2"], "mae": metrics["mae"], "rmse": metrics["rmse"],
                "pearson_r": metrics["pearson_r"], "sse": metrics["sse"],
                "ci_low_r2": np.nan, "ci_high_r2": np.nan, "n_bootstrap_valid": 0,
                "bootstrap_fraction_valid": np.nan, "confirmatory": False,
            })
    table = pd.DataFrame(rows)
    sse_m0 = table.loc[table["model"].eq("M0_covariates"), ["task", "sse"]].rename(columns={"sse": "sse_m0"})
    table = table.merge(sse_m0, on="task", how="left")
    table["prediction_error_improvement_vs_m0"] = np.where(
        table["sse_m0"] > NEAR_ZERO_STD, 1.0 - table["sse"] / table["sse_m0"], np.nan)
    return table.drop(columns=["sse", "sse_m0"])


def comparison_rows(predictions, intervals=None) -> "pd.DataFrame":
    import pandas as pd

    point = point_performance(predictions).set_index(["task", "model"])
    rows = []
    for task in point.index.get_level_values(0).unique():
        for name in (PRIMARY_COMPARISON, *SUPPORT_COMPARISONS):
            left, right = name.split("-")
            record = {"task": task, "comparison": name,
                      "role": "primary" if task == "success" and name == PRIMARY_COMPARISON else "support",
                      "delta_r2": point.at[(task, left), "r2"] - point.at[(task, right), "r2"],
                      "delta_mae": point.at[(task, left), "mae"] - point.at[(task, right), "mae"],
                      "delta_rmse": point.at[(task, left), "rmse"] - point.at[(task, right), "rmse"]}
            record.update({"ci_low_delta_r2": np.nan, "ci_high_delta_r2": np.nan,
                           "n_bootstrap_valid": 0, "inference_valid": False, "supported": False})
            if intervals is not None:
                record.update(intervals.get((task, name), {}))
            rows.append(record)
    return pd.DataFrame(rows)


def regression_metrics(y, prediction) -> dict:
    y = np.asarray(y, dtype=float).reshape(-1)
    prediction = np.asarray(prediction, dtype=float).reshape(-1)
    error = y - prediction
    sse = float(np.sum(error ** 2))
    center = float(np.sum((y - np.mean(y)) ** 2))
    r2 = np.nan if center <= NEAR_ZERO_STD else 1.0 - sse / center
    if np.std(y) <= NEAR_ZERO_STD or np.std(prediction) <= NEAR_ZERO_STD:
        correlation = np.nan
    else:
        correlation = float(np.corrcoef(y, prediction)[0, 1])
    return {"n": int(y.size), "sse": sse, "mae": float(np.mean(np.abs(error))),
            "rmse": float(np.sqrt(np.mean(error ** 2))), "r2": float(r2) if np.isfinite(r2) else np.nan,
            "pearson_r": correlation}
