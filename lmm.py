"""
Minimal linear mixed model (REML/ML) for a single grouping factor,
fit via exact closed-form per-cluster marginal likelihood. This is
valid (not an approximation) here because each patient contributes
only 2-3 observations, so per-cluster covariance matrices are tiny
and computed exactly.

Why a hand-written implementation: no internet access was available
to install statsmodels/lme4 while this was developed, so this was
built and validated against synthetic data with known ground-truth
parameters (see validate_lmm.py) before being used on the real
patient data. If you have statsmodels or lme4 available, you can
cross-check the results against `smf.mixedlm(...)` / `lmer(...)` --
they should agree, since this implements the same underlying model.

No external dependencies beyond numpy/scipy.
"""
import numpy as np
from scipy.optimize import minimize
from scipy.stats import chi2


def _L_from_theta(theta, q):
    """Unconstrained theta -> lower-triangular Cholesky factor L
    (random-effects covariance Psi = L L') and residual variance
    sigma2. Diagonal entries of L are exp-transformed so Psi is
    always guaranteed positive semi-definite during optimization."""
    L = np.zeros((q, q))
    idx = 0
    for i in range(q):
        for j in range(i + 1):
            if i == j:
                L[i, j] = np.exp(theta[idx])
            else:
                L[i, j] = theta[idx]
            idx += 1
    log_sigma = theta[idx]
    sigma2 = np.exp(log_sigma) ** 2
    return L, sigma2


def _n_theta(q):
    return q * (q + 1) // 2 + 1  # lower-triangular entries + log_sigma


def _neg_loglik(theta, y_list, X_list, Z_list, reml=True):
    q = Z_list[0].shape[1]
    L, sigma2 = _L_from_theta(theta, q)
    Psi = L @ L.T
    p = X_list[0].shape[1]
    XtVinvX = np.zeros((p, p))
    XtVinvy = np.zeros(p)
    logdetV_sum = 0.0
    Vinv_list = []
    for y_i, X_i, Z_i in zip(y_list, X_list, Z_list):
        n_i = len(y_i)
        V_i = Z_i @ Psi @ Z_i.T + sigma2 * np.eye(n_i)
        sign, logdet = np.linalg.slogdet(V_i)
        if sign <= 0:
            return 1e10
        logdetV_sum += logdet
        Vinv = np.linalg.inv(V_i)
        Vinv_list.append(Vinv)
        XtVinvX += X_i.T @ Vinv @ X_i
        XtVinvy += X_i.T @ Vinv @ y_i
    try:
        beta = np.linalg.solve(XtVinvX, XtVinvy)
    except np.linalg.LinAlgError:
        return 1e10
    rss = 0.0
    for y_i, X_i, Vinv in zip(y_list, X_list, Vinv_list):
        r = y_i - X_i @ beta
        rss += r @ Vinv @ r
    ll = -0.5 * (logdetV_sum + rss)
    if reml:
        sign, logdet_XtVinvX = np.linalg.slogdet(XtVinvX)
        ll -= 0.5 * logdet_XtVinvX
    return -ll


def fit_lmm(y_list, X_list, Z_list, reml=True, n_restarts=10, seed=0):
    """Fit a linear mixed model: y_i = X_i*beta + Z_i*b_i + eps_i,
    b_i ~ N(0, Psi), eps_i ~ N(0, sigma2*I). Random-restart
    Nelder-Mead is used because the covariance parameterization is
    non-convex; multiple restarts guard against local optima."""
    q = Z_list[0].shape[1]
    n_th = _n_theta(q)
    rng = np.random.default_rng(seed)
    best = None
    for _ in range(n_restarts):
        theta0 = rng.normal(size=n_th) * 0.5
        res = minimize(_neg_loglik, theta0, args=(y_list, X_list, Z_list, reml),
                        method="Nelder-Mead",
                        options=dict(xatol=1e-9, fatol=1e-9, maxiter=10000))
        if best is None or res.fun < best.fun:
            best = res
    L, sigma2 = _L_from_theta(best.x, q)
    Psi = L @ L.T
    p = X_list[0].shape[1]
    XtVinvX = np.zeros((p, p))
    XtVinvy = np.zeros(p)
    for y_i, X_i, Z_i in zip(y_list, X_list, Z_list):
        n_i = len(y_i)
        V_i = Z_i @ Psi @ Z_i.T + sigma2 * np.eye(n_i)
        Vinv = np.linalg.inv(V_i)
        XtVinvX += X_i.T @ Vinv @ X_i
        XtVinvy += X_i.T @ Vinv @ y_i
    beta = np.linalg.solve(XtVinvX, XtVinvy)
    beta_cov = np.linalg.inv(XtVinvX)
    n_obs = sum(len(y_i) for y_i in y_list)
    n_params = p + n_th
    aic = 2 * best.fun + 2 * n_params
    bic = 2 * best.fun + n_params * np.log(n_obs)
    return dict(beta=beta, beta_se=np.sqrt(np.diag(beta_cov)), beta_cov=beta_cov,
                Psi=Psi, sigma2=sigma2, negloglik=best.fun, aic=aic, bic=bic,
                converged=best.success, n_params=n_params, reml=reml)


def lrt(fit_reduced, fit_full, df):
    """Likelihood ratio test between nested models.
    - Comparing different RANDOM-effects structures with the SAME
      fixed effects: fit both with reml=True (valid REML LRT).
    - Comparing different FIXED-effects structures: fit both with
      reml=False (REML LRT is invalid here; must use ML).
    Note: LRTs for variance components tested at a boundary (e.g.
    variance = 0) are conservative under the standard chi2
    reference distribution -- treat borderline p-values on random
    effects tests with that caveat in mind."""
    stat = 2 * (fit_reduced["negloglik"] - fit_full["negloglik"])
    pval = chi2.sf(stat, df)
    return stat, pval


def blup(y_list, X_list, Z_list, beta, Psi, sigma2):
    """Per-cluster BLUPs (empirical-Bayes / shrinkage estimates of
    the random effects): b_i_hat = Psi Z_i' V_i^-1 (y_i - X_i beta)"""
    out = []
    for y_i, X_i, Z_i in zip(y_list, X_list, Z_list):
        n_i = len(y_i)
        V_i = Z_i @ Psi @ Z_i.T + sigma2 * np.eye(n_i)
        Vinv = np.linalg.inv(V_i)
        r = y_i - X_i @ beta
        b_i = Psi @ Z_i.T @ Vinv @ r
        out.append(b_i)
    return np.array(out)
