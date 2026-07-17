"""
Step 3: Alpha diversity analysis.

- Calculate Observed, Shannon, Simpson, Chao1 indices
- Kruskal-Wallis (or ANOVA) test across groups
- Optional Spearman correlation with numeric metadata
- Boxplot + jitter visualization

Input:  relative_abundance.xlsx, metadata.xlsx
Output: result/step3_alpha/alpha_diversity.xlsx, alpha_statistics.xlsx,
        alpha_spearman.xlsx, figures/*.pdf
"""

import os
import numpy as np
import pandas as pd
from scipy.stats import kruskal, f_oneway, spearmanr

from .config import PipelineConfig
from .theme import set_theme, plot_alpha_boxplot


def shannon_index(x):
    """Shannon diversity index (natural log)."""
    x = x[x > 0]
    if len(x) == 0:
        return 0.0
    p = x / x.sum()
    return -(p * np.log(p)).sum()


def simpson_index(x):
    """Gini-Simpson diversity index (1 - D)."""
    x = x[x > 0]
    if len(x) == 0:
        return 0.0
    p = x / x.sum()
    return 1 - (p * p).sum()


def chao1_index(x):
    """Chao1 richness estimator."""
    x = x[x > 0]
    obs = len(x)
    if obs == 0:
        return 0.0
    # Singletons and doubletons
    singles = np.sum(x == 1)
    doubles = np.sum(x == 2)
    return obs + (singles * (singles - 1)) / (2 * (doubles + 1)) if doubles > 0 else obs + (singles * (singles - 1)) / 2


METRIC_FUNCS = {
    "Observed": lambda x: int(np.sum(x > 0)),
    "Shannon": shannon_index,
    "Simpson": simpson_index,
    "Chao1": chao1_index,
}


def run_step3(cfg: PipelineConfig, rel_path: str = None, meta_path: str = None) -> str:
    """
    Execute Step 3: alpha diversity.
    Returns path to alpha_diversity.xlsx.
    """
    print("\n" + "=" * 60)
    print("Step 3: Alpha diversity analysis")
    print("=" * 60)

    if rel_path is None:
        rel_path = os.path.join(cfg.output_dir, "data", "relative_abundance.xlsx")
    if meta_path is None:
        meta_path = cfg.meta_data

    otu = pd.read_excel(rel_path)
    meta = pd.read_excel(meta_path)

    species_col = "Species_name"
    sample_cols = cfg._sample_cols_detected
    if not sample_cols:
        non_meta = [c for c in otu.columns if c != species_col and c != cfg._tax_col_detected]
        sample_cols = non_meta

    # Abundance matrix: rows=samples, cols=species
    otu_matrix = otu[sample_cols].T
    otu_matrix.index.name = cfg.meta_sample_col

    # Calculate alpha metrics
    results = []
    for sample, row in otu_matrix.iterrows():
        abundance = row.values.astype(float)
        entry = {cfg.meta_sample_col: sample}
        for metric in cfg.alpha_metrics:
            if metric in METRIC_FUNCS:
                entry[metric] = METRIC_FUNCS[metric](abundance)
        results.append(entry)

    alpha = pd.DataFrame(results)
    alpha = alpha.merge(meta, on=cfg.meta_sample_col, how="left")

    print(f"  Alpha diversity calculated for {len(alpha)} samples")
    print(f"  Metrics: {cfg.alpha_metrics}")
    print(alpha.head().to_string())

    # Save alpha diversity
    outdir = os.path.join(cfg.output_dir, "result", "step3_alpha")
    os.makedirs(outdir, exist_ok=True)

    alpha_path = os.path.join(outdir, "alpha_diversity.xlsx")
    alpha.to_excel(alpha_path, index=False)
    print(f"  Saved: {alpha_path}")

    # --- Group comparison statistics ---
    group_col = cfg.meta_group_col
    group_order = sorted(alpha[group_col].dropna().unique())
    cfg._group_order_detected = group_order
    # Use configured order if provided
    if cfg.group_order:
        group_order = [g for g in cfg.group_order if g in group_order]

    stat_results = []
    for metric in cfg.alpha_metrics:
        data = []
        for g in group_order:
            vals = alpha.loc[alpha[group_col] == g, metric].dropna()
            data.append(vals.values)

        if len(data) >= 2:
            if cfg.alpha_test == "anova":
                stat, p = f_oneway(*data)
                stat_name = "F_statistic"
            else:
                stat, p = kruskal(*data)
                stat_name = "Kruskal_H"
            stat_results.append({
                "Metric": metric,
                stat_name: stat,
                "p_value": p,
                "n_groups": len(data),
            })

    stat_df = pd.DataFrame(stat_results)
    stat_path = os.path.join(outdir, "alpha_statistics.xlsx")
    stat_df.to_excel(stat_path, index=False)
    print(f"  Statistics saved: {stat_path}")
    print(stat_df.to_string())

    # --- Spearman correlation (optional) ---
    if cfg.alpha_spearman and cfg.meta_numeric_col and cfg.meta_numeric_col in alpha.columns:
        spearman_results = []
        for metric in cfg.alpha_metrics:
            valid = alpha[[cfg.meta_numeric_col, metric]].dropna()
            if len(valid) >= 3:
                rho, p = spearmanr(valid[cfg.meta_numeric_col], valid[metric])
                spearman_results.append({
                    "Metric": metric,
                    "Spearman_rho": rho,
                    "p_value": p,
                    "n": len(valid),
                })
        if spearman_results:
            spearman_df = pd.DataFrame(spearman_results)
            spearman_path = os.path.join(outdir, "alpha_spearman.xlsx")
            spearman_df.to_excel(spearman_path, index=False)
            print(f"  Spearman correlation saved: {spearman_path}")
            print(spearman_df.to_string())

    # --- Plots ---
    set_theme(cfg)
    fig_dir = os.path.join(outdir, "figures")
    os.makedirs(fig_dir, exist_ok=True)

    colors = cfg.get_group_colors()
    for metric in cfg.alpha_metrics:
        plot_alpha_boxplot(
            alpha, metric, group_col,
            group_order, colors,
            fig_dir, cfg,
        )
        print(f"  Plot saved: {metric}_boxplot")

    print("  Alpha diversity analysis complete.")
    return alpha_path
