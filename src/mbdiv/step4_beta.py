"""Step 4: beta diversity - Bray-Curtis distance, PCoA, PERMANOVA, PERMDISP."""

import os
import numpy as np
import pandas as pd
from scipy.spatial.distance import pdist, squareform

from .config import PipelineConfig
from .theme import set_theme, plot_pcoa
from .stats_utils import bh_fdr, permdisp, permanova_manual


def run_step4(cfg: PipelineConfig, rel_path: str = None, meta_path: str = None):
    print("\n" + "=" * 60)
    print("Step 4: Beta diversity analysis")
    print("=" * 60)

    if rel_path is None:
        rel_path = os.path.join(cfg.output_dir, "data", "relative_abundance.xlsx")
    if meta_path is None:
        meta_path = cfg.meta_data

    base_dir = os.path.join(cfg.output_dir, "result", "step4_beta")
    dist_dir = os.path.join(base_dir, "distance")
    pcoa_dir = os.path.join(base_dir, "pcoa")
    stat_dir = os.path.join(base_dir, "statistics")
    fig_dir = os.path.join(base_dir, "figures")
    for d in [dist_dir, pcoa_dir, stat_dir, fig_dir]:
        os.makedirs(d, exist_ok=True)

    otu = pd.read_excel(rel_path, index_col=0)
    meta = pd.read_excel(meta_path)

    otu_t = otu.T
    print(f"  Distance matrix input: {otu_t.shape[0]} samples x {otu_t.shape[1]} species")

    distance = pd.DataFrame(
        squareform(pdist(otu_t, metric=cfg.beta_distance)),
        index=otu_t.index,
        columns=otu_t.index,
    )
    dist_path = os.path.join(dist_dir, f"{cfg.beta_distance}_distance.xlsx")
    distance.to_excel(dist_path)
    print(f"  Distance matrix saved: {dist_path}")

    # PCoA - prefer scikit-bio, fall back to manual
    try:
        from skbio.stats.ordination import pcoa as skbio_pcoa
        from skbio.stats.distance import DistanceMatrix

        dm = DistanceMatrix(distance.values, ids=distance.index.tolist())
        pcoa_result = skbio_pcoa(dm)

        coords = pcoa_result.samples.iloc[:, 0:2].copy()
        coords.columns = ["PC1", "PC2"]
        coords[cfg.meta_sample_col] = coords.index
        coords = coords[[cfg.meta_sample_col, "PC1", "PC2"]]

        variance = pd.DataFrame({
            "Axis": [str(i) for i in pcoa_result.proportion_explained.index],
            "Explained_variance": pcoa_result.proportion_explained.values,
        }).iloc[:5]

        print("  PCoA computed via scikit-bio")
        print(f"  Variance explained: PC1={variance.iloc[0]['Explained_variance']*100:.1f}%, "
              f"PC2={variance.iloc[1]['Explained_variance']*100:.1f}%")

    except ImportError:
        print("  scikit-bio not available, using manual PCoA (SVD-based)")
        coords, variance = _manual_pcoa(distance, cfg)
        print("  Manual PCoA computed")

    coords_path = os.path.join(pcoa_dir, "pcoa_coordinates.xlsx")
    coords.to_excel(coords_path, index=False)
    print(f"  PCoA coordinates saved: {coords_path}")

    var_path = os.path.join(pcoa_dir, "pcoa_variance.xlsx")
    variance.to_excel(var_path, index=False)
    print(f"  PCoA variance saved: {var_path}")

    # Build group labels aligned with distance matrix
    meta_indexed = meta.set_index(cfg.meta_sample_col)
    samples = list(distance.index)
    groups = meta_indexed.loc[samples, cfg.meta_group_col].values
    D = distance.values.astype(float)

    # Handle NaN/Inf in distance matrix
    if np.any(np.isnan(D)) or np.any(np.isinf(D)):
        print("  WARNING: NaN/Inf in distance matrix — replacing with 0")
        D = np.nan_to_num(D, nan=0.0, posinf=1.0, neginf=0.0)

    # PERMANOVA (manual implementation for R²)
    nperm = cfg.beta_permanova_permutations
    perm_F, perm_p, perm_R2 = permanova_manual(D, groups, nperm=nperm)
    print(f"  PERMANOVA: F={perm_F:.2f}, p={perm_p:.4f}, R²={perm_R2:.4f}")

    # PERMDISP (betadisper)
    disp_F, disp_p = permdisp(D, groups, nperm=nperm)
    disp_sig = "YES" if disp_p > 0.05 else "NO — PERMANOVA may be confounded"
    print(f"  PERMDISP:  F={disp_F:.2f}, p={disp_p:.4f}  (equal dispersion: {disp_sig})")

    # BH FDR correction across PERMANOVA + PERMDISP
    raw_pvals = {"PERMANOVA": perm_p, "PERMDISP": disp_p}
    test_names = list(raw_pvals.keys())
    p_list = [raw_pvals[k] for k in test_names]
    fdr_vals = bh_fdr(p_list)
    fdr_map = {k: round(float(f), 4) for k, f in zip(test_names, fdr_vals)}

    print("\n  --- FDR-corrected p-values (Benjamini-Hochberg) ---")
    for k in test_names:
        p_raw = raw_pvals[k]
        p_adj = fdr_map[k]
        sig = "***" if p_adj < 0.001 else "**" if p_adj < 0.01 else "*" if p_adj < 0.05 else "ns"
        print(f"    {k:12s}: p_raw={p_raw:.4f}  p_adj={p_adj:.4f} {sig}")

    # Save statistics
    stat_df = pd.DataFrame([
        {
            "Test": "PERMANOVA",
            "F_statistic": round(perm_F, 2),
            "R_squared": round(perm_R2, 4),
            "p_value": round(perm_p, 4),
            "p_fdr": fdr_map["PERMANOVA"],
            "permutations": nperm,
        },
        {
            "Test": "PERMDISP",
            "F_statistic": disp_F,
            "R_squared": "",
            "p_value": disp_p,
            "p_fdr": fdr_map["PERMDISP"],
            "permutations": nperm,
        },
    ])
    stat_path = os.path.join(stat_dir, "beta_statistics.xlsx")
    stat_df.to_excel(stat_path, index=False)
    print(f"  Statistics saved: {stat_path}")

    permanova_stats = {
        "F": perm_F,
        "p": perm_p,
        "p_fdr": fdr_map["PERMANOVA"],
        "R2": perm_R2,
    }
    permdisp_stats = {
        "F": disp_F,
        "p": disp_p,
        "p_fdr": fdr_map["PERMDISP"],
    }

    # PCoA plot
    set_theme(cfg)
    group_order = cfg.get_group_order()
    if not group_order:
        group_order = sorted(meta[cfg.meta_group_col].dropna().unique().tolist())
    colors = cfg.get_group_colors()

    plot_pcoa(
        coords, variance, meta,
        cfg.meta_group_col, cfg.meta_sample_col,
        group_order, colors,
        permanova_stats, permdisp_stats,
        fig_dir, cfg,
    )
    print(f"  PCoA plot saved to: {fig_dir}")

    print("  Beta diversity analysis complete.")
    return dist_path, coords_path, perm_p


def _manual_pcoa(distance: pd.DataFrame, cfg: PipelineConfig):
    n = distance.shape[0]
    D = distance.values.astype(float)

    D_sq = D ** 2
    J = np.eye(n) - np.ones((n, n)) / n
    G = -0.5 * J @ D_sq @ J

    eigenvalues, eigenvectors = np.linalg.eigh(G)

    idx = np.argsort(eigenvalues)[::-1]
    eigenvalues = eigenvalues[idx]
    eigenvectors = eigenvectors[:, idx]

    coords = eigenvectors[:, :2] * np.sqrt(np.maximum(eigenvalues[:2], 0))

    total = np.sum(eigenvalues[eigenvalues > 0])
    prop_explained = eigenvalues / total if total > 0 else eigenvalues * 0

    coords_df = pd.DataFrame({
        cfg.meta_sample_col: distance.index,
        "PC1": coords[:, 0],
        "PC2": coords[:, 1],
    })

    variance_df = pd.DataFrame({
        "Axis": [f"PC{i+1}" for i in range(min(5, len(eigenvalues)))],
        "Explained_variance": prop_explained[:5],
    })

    return coords_df, variance_df
