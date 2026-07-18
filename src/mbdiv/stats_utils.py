"""Shared statistical utility functions."""

import numpy as np
from scipy.stats import f_oneway


def bh_fdr(pvals):
    """Benjamini-Hochberg FDR correction. Handles NaN p-values."""
    p = np.array(pvals, dtype=float)
    n = len(p)
    if n == 0:
        return np.array([])

    nan_mask = np.isnan(p)
    p_clean = p.copy()
    p_clean[nan_mask] = 1.0

    order = np.argsort(p_clean)
    ranked = p_clean[order]
    fdr = ranked * n / (np.arange(n) + 1)
    fdr = np.minimum.accumulate(fdr[::-1])[::-1]
    fdr = np.clip(fdr, 0, 1)
    result = np.empty(n)
    result[order] = fdr
    result[nan_mask] = np.nan
    return result


def permdisp(D, groups, nperm=999):
    """Multivariate dispersion test (Anderson 2006).

    Tests whether within-group dispersions are equal — a key
    assumption of PERMANOVA. Returns (F, p).
    """
    n = D.shape[0]
    grp = np.array(groups)
    ug = np.unique(grp)

    H = np.eye(n) - np.ones((n, n)) / n
    B = -0.5 * H @ D @ H
    ev, evec = np.linalg.eigh(B)
    idx = np.argsort(ev)[::-1]
    ev, evec = ev[idx], evec[:, idx]
    pos = ev > 0
    coords = evec[:, pos] * np.sqrt(ev[pos])

    def disp_F(g):
        grouped_dists = [
            np.sqrt(((coords[g == u] - coords[g == u].mean(axis=0)) ** 2).sum(axis=1))
            for u in ug if (g == u).sum() >= 2
        ]
        if len(grouped_dists) < 2:
            return 0.0
        F, _ = f_oneway(*grouped_dists)
        return F

    F_obs = disp_F(grp)
    if np.isnan(F_obs) or F_obs == 0.0:
        return 0.0, 1.0

    np.random.seed(42)
    perm_Fs = [disp_F(np.random.permutation(grp)) for _ in range(nperm)]
    perm_Fs = [f for f in perm_Fs if not np.isnan(f)]
    perm_p = np.mean(np.array(perm_Fs) >= F_obs) if perm_Fs else 1.0
    return round(F_obs, 2), round(perm_p, 4)


def permanova_manual(D, groups, nperm=999):
    """PERMANOVA (Anderson 2001). Returns (F, p, R2).

    R2 = SS_between / SS_total — proportion of variance explained.
    """
    n = D.shape[0]
    grp = np.array(groups)
    ug = np.unique(grp)
    total_ss = np.sum(D ** 2) / n

    def f_stat(g):
        w = sum(
            np.sum(D[np.ix_(g == u, g == u)] ** 2) / max((g == u).sum(), 1)
            for u in ug
        )
        b = total_ss - w
        F = (b / (len(ug) - 1)) / (w / (n - len(ug))) if w > 0 else 0
        R2 = b / total_ss if total_ss > 0 else 0
        return F, R2

    F_obs, R2_obs = f_stat(grp)
    np.random.seed(42)
    perm_results = [f_stat(np.random.permutation(grp)) for _ in range(nperm)]
    F_perm = [r[0] for r in perm_results]
    p_val = np.mean(np.array(F_perm) >= F_obs)
    return F_obs, p_val, R2_obs
