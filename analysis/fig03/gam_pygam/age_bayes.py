"""Gaussian spline Bayes factors with proper age priors and shared references.

The common nuisance coefficients have Lebesgue reference measure and the
common variance prior is d(sigma^2)/sigma^2. Their arbitrary constants cancel
between models on identical rows. Model-specific age priors are proper.
Evidence is reported relative to M1, never as an absolute marginal density.
"""

from __future__ import annotations

import warnings

import numpy as np
from scipy.integrate import IntegrationWarning, quad
from scipy.linalg import cho_factor, cho_solve
from scipy.optimize import minimize_scalar

PRIOR_SCALES = (0.25, 0.5, 1.0)
MAIN_SCALE = 0.5
COMPARISONS = {"m3_vs_m1": (3, 1), "m3_vs_m2": (3, 2), "m2_vs_m1": (2, 1)}


def spline_components(age: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Preserve pyGAM's five-spline space and second-difference penalty.

    Z has no intercept/linear component. Whiten the positive penalty
    eigenvalues, then scale tr(Z Z')/n to one, so tau measures prior RMS
    nonlinearity in residual-SD units, independent of arbitrary basis units.
    """
    from pygam import s

    age = np.asarray(age, dtype=float)
    x = (age - age.mean()) / age.std(ddof=0)
    term = s(0, n_splines=5)
    term.compile(age[:, None])
    basis = term.build_columns(age[:, None]).toarray()
    difference = np.diff(np.eye(5), n=2, axis=0)
    values, vectors = np.linalg.eigh(difference.T @ difference)
    positive = values > values.max() * 1e-10
    z = basis @ (vectors[:, positive] / np.sqrt(values[positive]))
    linear = np.column_stack([np.ones(len(age)), x])
    z -= linear @ np.linalg.lstsq(linear, z, rcond=None)[0]
    z /= np.sqrt(np.sum(z * z) / len(age))
    target = np.column_stack([linear, z])
    if (np.linalg.matrix_rank(target) != 5 or
            not np.allclose(target @ np.linalg.lstsq(target, basis, rcond=None)[0],
                            basis, rtol=1e-8, atol=1e-9)):
        raise ValueError("The Bayesian age decomposition does not span the original five-spline space")
    return x, z


def _safe_exp(value: float) -> float:
    # log BF remains authoritative when a linear-scale BF over/underflows.
    return float(np.exp(value)) if value <= np.log(np.finfo(float).max) else float("inf")


def evidence_label(log_bf: float) -> str:
    if log_bf >= np.log(10):
        return "strong_support_full"
    if log_bf >= np.log(3):
        return "support_full"
    if log_bf <= -np.log(10):
        return "strong_support_reduced"
    if log_bf <= -np.log(3):
        return "support_reduced"
    return "inconclusive"


def fit_age_bayes(y, age, covariates, scales=PRIOR_SCALES) -> tuple[dict, list[dict]]:
    """Integrate coefficients/variance analytically, tau by adaptive quadrature."""
    y, age = np.asarray(y, dtype=float), np.asarray(age, dtype=float)
    covariates = np.asarray(covariates, dtype=float)
    if covariates.ndim == 1:
        covariates = covariates[:, None]
    finite = np.isfinite(y) & np.isfinite(age) & np.isfinite(covariates).all(axis=1)
    n = int(finite.sum())
    result = {"status": "insufficient_observations", "n_obs": n}
    if n < 12:
        return result, []
    y, age, covariates = y[finite], age[finite], covariates[finite]
    if age.std() == 0 or np.ptp(y) == 0:
        result["status"] = "constant_age_or_response"
        return result, []
    nuisance = np.column_stack([np.ones(n), covariates - covariates.mean(axis=0)])
    q0, singular, _ = np.linalg.svd(nuisance, full_matrices=False)
    rank = int(np.sum(singular > singular[0] * max(nuisance.shape) * np.finfo(float).eps))
    q0 = q0[:, :rank]
    nu = n - rank
    result.update(rank_nuisance=rank, residual_dimension=nu, age_mean=float(age.mean()),
                  age_sd=float(age.std()), age_min=float(age.min()), age_max=float(age.max()))
    if nu <= 4:
        return result, []
    x, z = spline_components(age)
    w = np.column_stack([x, z])
    w -= q0 @ (q0.T @ w)
    yr = y - q0 @ (q0.T @ y)
    u, singular, vt = np.linalg.svd(w, full_matrices=False)
    if singular[-1] <= singular[0] * max(w.shape) * np.finfo(float).eps:
        result["status"] = "age_confounded_with_nuisance"
        return result, []
    projected = u.T @ yr
    remainder = yr - u @ projected
    rss_remainder, rss0 = float(remainder @ remainder), float(yr @ yr)
    if rss_remainder <= np.finfo(float).eps * rss0:
        result["status"] = "zero_reference_residual_variance"
        return result, []
    r = singular[:, None] * vt

    def log_evidence(linear_scale, tau):
        scaled = r * np.array([linear_scale, tau, tau, tau])
        factor, lower = cho_factor(np.eye(4) + scaled @ scaled.T, lower=True)
        rss = rss_remainder + float(projected @ cho_solve((factor, lower), projected))
        return -float(np.log(np.diag(factor)).sum()) - nu / 2 * np.log(rss / rss0)

    def integrate(linear_scale, nonlinear_scale):
        constant = .5 * np.log(2 / np.pi)

        def log_integrand(tau_standardized):
            return log_evidence(linear_scale, nonlinear_scale * tau_standardized) + constant - tau_standardized ** 2 / 2

        upper = np.sqrt(2 * nu) + 10
        optimum = minimize_scalar(lambda value: -log_integrand(value), bounds=(0, upper), method="bounded")
        mode = float(optimum.x)
        offset = max(log_integrand(0), log_integrand(mode))

        def integrand(value):
            return np.exp(log_integrand(value) - offset)

        with warnings.catch_warnings():
            warnings.simplefilter("error", IntegrationWarning)
            left, error_left = quad(integrand, 0, mode, epsabs=1e-10, epsrel=1e-8, limit=180)
            right, error_right = quad(integrand, mode, np.inf, epsabs=1e-10, epsrel=1e-8, limit=180)
            total = left + right
            # Independent partition/tolerance check prevents silently accepting
            # a missed narrow mode or unreliable numerical evidence.
            check_left, _ = quad(integrand, 0, mode / 2, epsabs=1e-12, epsrel=1e-10, limit=220)
            check_mid, _ = quad(integrand, mode / 2, mode + 1, epsabs=1e-12, epsrel=1e-10, limit=220)
            check_right, _ = quad(integrand, mode + 1, np.inf, epsabs=1e-12, epsrel=1e-10, limit=220)
        relative_error = (error_left + error_right) / total
        log_difference = abs(np.log((check_left + check_mid + check_right) / total))
        if not np.isfinite(total) or total <= 0 or relative_error > 1e-6 or log_difference > 1e-6:
            raise ArithmeticError("Bayesian spline evidence failed quadrature convergence checks")
        return offset + np.log(total), relative_error, log_difference

    scales = tuple(float(value) for value in scales)
    if MAIN_SCALE not in scales or any(not np.isfinite(value) or value <= 0 for value in scales):
        raise ValueError("Prior scales must be positive and include the main scale 0.5")
    rows = []
    for linear_scale in scales:
        log_m2 = log_evidence(linear_scale, 0)
        for nonlinear_scale in scales:
            log_m3, error, difference = integrate(linear_scale, nonlinear_scale)
            row = {"linear_prior_scale": linear_scale, "nonlinear_prior_scale": nonlinear_scale,
                   "log_evidence_m1_relative": 0.0, "log_evidence_m2_relative": log_m2,
                   "log_evidence_m3_relative": log_m3, "quadrature_relative_error": error,
                   "quadrature_check_log_difference": difference, "integration_converged": True}
            logs = {1: 0.0, 2: log_m2, 3: log_m3}
            for comparison, (full, reduced) in COMPARISONS.items():
                value = float(logs[full] - logs[reduced])
                row.update({f"log_bf10_{comparison}": value, f"bf10_{comparison}": _safe_exp(value),
                            f"bf01_{comparison}": _safe_exp(-value), f"evidence_{comparison}": evidence_label(value)})
            rows.append(row)
            if linear_scale == nonlinear_scale == MAIN_SCALE:
                result.update(row)
    for comparison in COMPARISONS:
        logs = [row[f"log_bf10_{comparison}"] for row in rows]
        result.update({f"sensitivity_log_bf10_min_{comparison}": min(logs),
                       f"sensitivity_log_bf10_max_{comparison}": max(logs),
                       f"sensitivity_direction_stable_{comparison}": min(logs) > 0 or max(logs) < 0,
                       f"sensitivity_evidence_labels_{comparison}": ";".join(sorted({row[f"evidence_{comparison}"] for row in rows}))})
    result["status"] = "ok"
    return result, rows
