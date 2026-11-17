"""Training-only two-output PLS curves, shared reference score, and metrics."""
from __future__ import annotations

from dataclasses import dataclass
import warnings

import numpy as np
from sklearn.cross_decomposition import PLSRegression
from sklearn.exceptions import ConvergenceWarning
from .io import MODELS, SCHEMES

EPS = 1e-12
PLS_MAX_ITER = 10000
PLS_TOL = 1e-6


def covariates(cov: np.ndarray, center=None, scale=None):
    """Only gender and meanFD. Sites group validation, never enter prediction."""
    cov = np.asarray(cov, float)
    if cov.ndim != 2 or cov.shape[1] != 2 or not np.isfinite(cov).all():
        raise ValueError("Covariates must be finite [gender, meanFD], without age or parent education")
    if not np.isin(cov[:, 0], [0, 1]).all():
        raise ValueError("Gender must be coded 0/1")
    center = cov.mean(0) if center is None else np.asarray(center)
    scale = cov.std(0) if scale is None else np.asarray(scale)
    scale = np.where(scale > EPS, scale, 1.)
    return np.column_stack([np.ones(len(cov)), (cov - center) / scale]), center, scale


@dataclass
class Prepared:
    center_c: np.ndarray
    scale_c: np.ndarray
    beta_x: np.ndarray
    beta_y: np.ndarray
    scale_x: np.ndarray
    scale_y: np.ndarray
    keep: np.ndarray
    x: np.ndarray
    y: np.ndarray
    rank_c: int

    def transform(self, x, y, cov):
        c, _, _ = covariates(cov, self.center_c, self.scale_c)
        xz = (np.asarray(x) - c @ self.beta_x) / self.scale_x
        yz = None if y is None else (np.asarray(y) - c @ self.beta_y) / self.scale_y
        return xz, yz, c @ self.beta_y


def prepare(x, y, cov) -> Prepared:
    x, y = np.asarray(x, float), np.asarray(y, float)
    if y.shape != (len(x), 2) or not np.isfinite(x).all() or not np.isfinite(y).all():
        raise ValueError("PLS2 requires finite X and an aligned N x 2 Y [age, success]")
    c, center, scale = covariates(cov)
    bx, _, rank, _ = np.linalg.lstsq(c, x, rcond=None)
    by = np.linalg.lstsq(c, y, rcond=None)[0]
    rx, ry = x - c @ bx, y - c @ by
    sx, sy = rx.std(0, ddof=1), ry.std(0, ddof=1)
    if np.any(sy <= EPS):
        raise ValueError("A behavior has constant training residuals")
    keep = sx > EPS
    sx = np.where(keep, sx, 1.)
    return Prepared(center, scale, bx, by, sx, sy, keep, rx / sx, ry / sy, int(rank))


@dataclass
class Curve:
    mean_x: np.ndarray
    mean_y: np.ndarray
    coefs: dict
    weights_y: np.ndarray | None
    weights_x: np.ndarray | None
    n_iter: tuple[int, ...] = ()

    @property
    def max_k(self):
        return max(self.coefs, default=0)

    def predict(self, x, k):
        if k == 0:
            return np.zeros((len(x), 2))
        return (x - self.mean_x) @ self.coefs[k] + self.mean_y


def fit_curve(x, y, k_max=10, rank_c=1) -> Curve:
    cap = min(int(k_max), x.shape[1], len(x) - rank_c - 1)
    if cap <= 0:
        return Curve(np.zeros(x.shape[1]), np.zeros(2), {}, None, None)
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", message="(?i)y residual is constant.*", category=UserWarning)
        warnings.simplefilter("error", ConvergenceWarning)
        try:
            model = PLSRegression(n_components=cap, scale=False,
                                  max_iter=PLS_MAX_ITER, tol=PLS_TOL).fit(x, y)
        except ConvergenceWarning as error:
            raise ConvergenceWarning(
                f"PLS did not converge: max_iter={PLS_MAX_ITER}, tol={PLS_TOL}, "
                f"n_train={len(x)}, n_features={x.shape[1]}, n_components={cap}; {error}"
            ) from error
    obtained = len(model.n_iter_)
    coefs = {}
    for k in range(1, obtained + 1):
        w, p, q = model.x_weights_[:, :k], model.x_loadings_[:, :k], model.y_loadings_[:, :k]
        coefs[k] = w @ np.linalg.pinv(p.T @ w) @ q.T
    return Curve(x.mean(0), y.mean(0), coefs,
                 model.y_weights_[:, 0].copy() if obtained else None,
                 model.x_weights_[:, 0].copy() if obtained else None,
                 tuple(int(v) for v in model.n_iter_))


def reference_axis(y, curve: Curve):
    """First Y direction of the *single full-X reference*, shared by all models."""
    if curve.weights_y is None or np.linalg.norm(curve.weights_y) <= EPS:
        raise ValueError("The full-X reference PLS has no nonzero first Y direction")
    q = curve.weights_y / np.linalg.norm(curve.weights_y)
    pivot = 0 if abs(q[0]) > EPS else int(np.argmax(abs(q)))
    if q[pivot] < 0:
        q = -q
    u = y @ q
    scale = float(u.std(ddof=1))
    if scale <= EPS:
        raise ValueError("Reference behavior score has zero variance")
    loadings = [float(np.corrcoef(y[:, j], u)[0, 1]) for j in range(2)]
    return q, scale, loadings


def r2(y, pred):
    y, pred = np.asarray(y, float), np.asarray(pred, float)
    ss = np.sum((y - y.mean(axis=0)) ** 2, axis=0)
    if np.any(ss <= EPS):
        return float("nan")
    return float(np.mean(1 - np.sum((y - pred) ** 2, axis=0) / ss))


def metric(y, pred):
    error = np.asarray(y) - np.asarray(pred)
    return {"r2": r2(y, pred), "mae": float(np.mean(abs(error))),
            "rmse": float(np.sqrt(np.mean(error ** 2)))}


def select_k(losses):
    values = np.asarray(losses)
    if not np.isfinite(values).all():
        raise ValueError("Nonfinite inner validation loss")
    return int(np.flatnonzero(values <= values.min() + 1e-10)[0])
