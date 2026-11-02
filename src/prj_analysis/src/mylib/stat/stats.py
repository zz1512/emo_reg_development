import numpy as np
from typing import Tuple, Optional, Sequence, Dict
from scipy.stats import f as f_dist
from scipy.stats import norm
from pygam import LinearGAM, s, l
from dataclasses import dataclass

@dataclass
class GAMConfig:
    smooth_var: str = "age"
    k_splines: int = 4
    lam_grid: Tuple[float, float, int] = (-4, 4, 9)  # log10 range -> 10**linspace
    center_smooth: bool = True
    grid_kind: str = "quantile"
    grid_size: int = 101
    lam_vals: Optional[np.ndarray] = None
    bootstrap_deriv_se: bool = False
    n_boot: int = 200
    verbose: bool = False

def nested_anova_gaussian(
    rss_full: float,
    rss_reduced: float,
    n_obs: int,
    edf_full: float,
    edf_reduced: float,
    eps: float = 1e-12
) -> Tuple[float, float, float, float, float]:
    """
    Generic nested-model ANOVA (Gaussian). Works for GAM/GLM if you can supply RSS and edf.
    Returns:
      F, pval, df1, df2, partial_R2

    Formulas:
      df1 = edf_full - edf_reduced          (effective added dof by smooth term)
      df2 = n_obs  - edf_full               (residual dof of full model)
      F   = ((RSS_reduced - RSS_full)/df1) / (RSS_full/df2)
      partial R^2 = max(0, 1 - RSS_full/RSS_reduced)
    """
    rss_full = float(rss_full)
    rss_reduced = float(rss_reduced)

    # Guard: reduced RSS should be >= full RSS; otherwise clamp at equality
    num = max(rss_reduced - rss_full, 0.0)

    df1 = max(float(edf_full) - float(edf_reduced), eps)
    df2 = max(float(n_obs) - float(edf_full), eps)

    # If no additional dof or degenerate denominator, fall back to p=1
    if df1 <= eps or df2 <= eps or rss_full <= eps:
        return 0.0, 1.0, df1, df2, 0.0

    F = (num / df1) / (rss_full / df2)
    # Upper tail p-value
    p = float(f_dist.sf(F, df1, df2))

    # Partial R^2 under the Gaussian model
    pr2 = 0.0
    if rss_reduced > eps:
        pr2 = max(0.0, 1.0 - (rss_full / rss_reduced))

    return float(F), p, float(df1), float(df2), float(pr2)

def fdr_bh(
    pvals: Sequence[float],
    alpha: float = 0.05
) -> Dict[str, np.ndarray]:
    """
    Benjamini–Hochberg FDR control.
    Returns a dict with:
      'reject': boolean array
      'qvals' : BH-adjusted p-values (a.k.a. FDR q-values)
      'crit_p': critical p-value threshold (scalar)
      'order' : argsort indices used internally
    """
    p = np.asarray(pvals, dtype=float)
    m = p.size
    order = np.argsort(p)
    p_sorted = p[order]
    # BH critical line
    bh_line = (np.arange(1, m+1) / m) * alpha
    below = p_sorted <= bh_line
    k = np.max(np.where(below)[0]) + 1 if np.any(below) else 0
    crit_p = p_sorted[k-1] if k > 0 else 0.0

    # q-values (step-up)
    q = np.empty_like(p_sorted)
    cmin = 1.0
    for i in range(m-1, -1, -1):
        cmin = min(cmin, (m / (i+1)) * p_sorted[i])
        q[i] = cmin
    # Undo sorting
    qvals = np.empty_like(q)
    qvals[order] = q
    reject = p <= crit_p
    return {"reject": reject, "qvals": qvals, "crit_p": np.array(crit_p), "order": order}

def rss_from_predictions(y_true: np.ndarray, y_hat: np.ndarray) -> float:
    """Residual Sum of Squares."""
    y = np.asarray(y_true, float).ravel()
    yhat = np.asarray(y_hat, float).ravel()
    mask = np.isfinite(y) & np.isfinite(yhat)
    if not np.any(mask):
        return np.nan
    r = y[mask] - yhat[mask]
    return float(np.dot(r, r))

def edf_from_pygam(model) -> Optional[float]:
    try:
        stats = getattr(model, "statistics_", None)
        if stats is not None and "edof" in stats:
            return float(stats["edof"])
        # older pygam may store 'edof_per_term' and 'edof'
        if stats is not None and "edof_per_term" in stats and np.size(stats["edof_per_term"])>0:
            return float(np.sum(stats["edof_per_term"]))
    except Exception:
        pass
    return None

def _reduced_design(n_obs: int, Xlin: Optional[np.ndarray]) -> np.ndarray:
    """Intercept column, then the linear covariates. This is the no-age model."""
    intercept = np.ones((n_obs, 1), dtype=float)
    if Xlin is None or np.size(Xlin) == 0:
        return intercept
    covariates = np.asarray(Xlin, dtype=float)
    if covariates.ndim == 1:
        covariates = covariates.reshape(-1, 1)
    if covariates.shape[0] != n_obs:
        raise ValueError(
            f"Reduced-model covariates have {covariates.shape[0]} rows; expected {n_obs}"
        )
    return np.hstack([intercept, covariates])


def fit_gam_full_and_reduced(
    ysub: np.ndarray,           # ROI y for valid subjects
    Xs: np.ndarray,             # shape (n, 1), the smooth var (age)
    Xlin: Optional[np.ndarray], # shape (n, p_lin) or None
    *,
    k_splines: int = 4,         # closer to the paper (k=3)
    lam_grid = 10.0**np.linspace(-4, 4, 9),  # wide search
    max_iter: int = 1000,
    verbose: bool = False
):
    """
    Fit FULL (smooth + linear) and REDUCED (intercept + linear) models.
    The reduced design always contains a column of ones, so it is nested in
    the full model, which uses fit_intercept=True. Return
    (rss_full, edf_full, rss_red, edf_red, model_full).
    """
    # Build design for FULL. Only the age smooth is penalized. A 1-d lam
    # grid would be copied onto every term, which ridge-shrinks the covariates
    # and makes the full model unable to contain the reduced model.
    if Xlin is not None and Xlin.size > 0:
        Z_full = np.hstack([Xs, Xlin])  # col 0 will be smoothed
        terms = s(0, n_splines=k_splines)
        for j in range(1, Z_full.shape[1]):
            terms += l(j, penalties=None)
        n_linear = Z_full.shape[1] - 1
    else:
        Z_full = Xs
        terms = s(0, n_splines=k_splines)
        n_linear = 0

    mf = LinearGAM(terms, max_iter=max_iter, fit_intercept=True, verbose=verbose)
    smooth_lam = np.asarray(lam_grid, dtype=float).reshape(-1)
    if n_linear:
        lam_search = np.column_stack([
            smooth_lam,
            np.zeros((smooth_lam.size, n_linear), dtype=float),
        ])
    else:
        lam_search = smooth_lam
    try:
        mf.gridsearch(Z_full, ysub, lam=lam_search, progress=False)
    except Exception:
        mf.fit(Z_full, ysub)

    yhat_full = mf.predict(Z_full)
    rss_full = rss_from_predictions(ysub, yhat_full)
    edf_full = edf_from_pygam(mf)

    # REDUCED: intercept plus the same linear covariates, without the age smooth.
    Xr = _reduced_design(ysub.shape[0], Xlin)
    XtX = Xr.T @ Xr
    try:
        beta = np.linalg.solve(XtX, Xr.T @ ysub)
    except np.linalg.LinAlgError:
        beta = np.linalg.pinv(XtX, rcond=1e-12) @ (Xr.T @ ysub)
    yhat_red = Xr @ beta
    edf_red = float(np.linalg.matrix_rank(Xr))

    rss_red = rss_from_predictions(ysub, yhat_red)

    # If pygam didn't provide edf, approximate: edf_full ≈ number of columns in model matrix
    if edf_full is None:
        # crude fallback: intercept + smooth basis + linear columns
        p_lin = 0 if Xlin is None or np.size(Xlin) == 0 else np.asarray(Xlin).reshape(ysub.shape[0], -1).shape[1]
        edf_full = float(min(ysub.size - 1, 1 + k_splines + p_lin))

    return rss_full, edf_full, rss_red, edf_red, mf

def simul_band_from_pointwise(
    mean: np.ndarray,
    se: np.ndarray,
    alpha: float = 0.05,
    method: str = "sidak"
) -> Tuple[np.ndarray, np.ndarray, float]:
    """
    Build a simultaneous (family-wise) band from pointwise SEs (approximate).
    	- mean, se: arrays over the evaluation grid.
    	- method: 'sidak' or 'bonferroni'.
    Returns (lower, upper, crit_z)
    """
    mean = np.asarray(mean, float)
    se   = np.asarray(se, float)
    m = np.sum(np.isfinite(se) & (se > 0))
    if m == 0:
        return np.full_like(mean, np.nan), np.full_like(mean, np.nan), np.nan

    if method == "sidak":
        # 1 - (1 - alpha)^(1/m) ≈ alpha/m for small alpha
        alpha_pt = 1.0 - (1.0 - alpha) ** (1.0 / m)
    else:  # bonferroni
        alpha_pt = alpha / m
    zc = float(norm.ppf(1.0 - alpha_pt / 2.0))
    lw = mean - zc * se
    up = mean + zc * se
    return lw, up, zc
