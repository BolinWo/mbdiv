"""Step 3: alpha diversity (Observed, Shannon, Simpson, Chao1) + stats + boxplots."""

import os
import numpy as np
import pandas as pd
from scipy.stats import kruskal, f_oneway, spearmanr

from .config import PipelineConfig
from .theme import set_theme, plot_alpha_boxplot
from .stats_utils import bh_fdr


def shannon_index(x):
    x = x[x > 0]
    if len(x) == 0:
        return 0.0
    p = x / x.sum()
    return -(p * np.log(p)).sum()


def simpson_index(x):
    x = x[x > 0]
    if len(x) == 0:
        return 0.0
    p = x / x.sum()
    return 1 - (p * p).sum()


def chao1_index(x):
    x = x[x > 0]
    obs = len(x)
    if obs == 0:
        return 0.0
    singles = np.sum(x == 1)
    doubles = np.sum(x == 2)
    if doubles > 0:
        return obs + (singles * (singles - 1)) / (2 * (doubles + 1))
    return obs + (singles * (singles - 1)) / 2


def is_count_data(matrix):
    vals = matrix.values.astype(float).ravel()
    vals = vals[vals > 0]
    if len(vals) == 0:
        return False
    return np.allclose(vals, np.round(vals), atol=1e-6)


def detect_group_order(alpha_df, meta_df, group_col, numeric_col, cfg):
    groups = alpha_df[group_col].dropna().unique().tolist()

    if cfg.group_order:
        ordered = [g for g in cfg.group_order if g in groups]
        remaining = [g for g in groups if g not in ordered]
        return ordered + remaining

    if numeric_col and numeric_col in alpha_df.columns:
        group_means = alpha_df.groupby(group_col)[numeric_col].mean()
        group_means = group_means.reindex(groups)
        ordered = group_means.sort_values().index.tolist()
        print(f"  Group order (by {numeric_col} mean, ascending): {ordered}")
        return ordered

    if meta_df is not None and group_col in meta_df.columns:
        meta_order = meta_df[group_col].dropna().drop_duplicates().tolist()
        ordered = [g for g in meta_order if g in groups]
        remaining = [g for g in groups if g not in ordered]
        result = ordered + remaining
        print(f"  Group order (by metadata appearance): {result}")
        return result

    return sorted(groups)


METRIC_FUNCS = {
    "Observed": lambda x: int(np.sum(x > 0)),
    "Shannon": shannon_index,
    "Simpson": simpson_index,
    "Chao1": chao1_index,
}


def run_step3(cfg: PipelineConfig, merged_path: str = None, meta_path: str = None) -> str:
    print("\n" + "=" * 60)
    print("Step 3: Alpha diversity analysis")
    print("=" * 60)

    if merged_path is None:
        merged_path = os.path.join(cfg.output_dir, "data", "merged_species_clean.xlsx")
    if meta_path is None:
        meta_path = cfg.meta_data

    otu = pd.read_excel(merged_path)
    meta = pd.read_excel(meta_path)

    species_col = "Species_name"
    sample_cols = cfg._sample_cols_detected
    if not sample_cols:
        sample_cols = [c for c in otu.columns if c != species_col and c != cfg._tax_col_detected]

    abundance_matrix = otu[sample_cols]
    count_data = is_count_data(abundance_matrix)

    metrics = list(cfg.alpha_metrics)
    if "Chao1" in metrics and not count_data:
        print("  WARNING: Input data is not integer count data (RPKM or relative abundance).")
        print("           Chao1 requires integer read counts to identify singletons/doubletons.")
        print("           Chao1 will be excluded from this analysis.")
        metrics = [m for m in metrics if m != "Chao1"]
    elif "Chao1" in metrics and count_data:
        print("  Input detected as integer count data -> Chao1 is valid.")

    otu_matrix = otu[sample_cols].T
    otu_matrix.index.name = cfg.meta_sample_col

    results = []
    for sample, row in otu_matrix.iterrows():
        abundance = row.values.astype(float)
        entry = {cfg.meta_sample_col: sample}
        for metric in metrics:
            if metric in METRIC_FUNCS:
                entry[metric] = METRIC_FUNCS[metric](abundance)
        results.append(entry)

    alpha = pd.DataFrame(results)
    alpha = alpha.merge(meta, on=cfg.meta_sample_col, how="left")

    if not cfg.meta_numeric_col:
        numeric_cols = meta.select_dtypes(include=[np.number]).columns.tolist()
        numeric_cols = [c for c in numeric_cols if c != cfg.meta_sample_col]
        if numeric_cols:
            cfg.meta_numeric_col = numeric_cols[0]
            print(f"  Auto-detected numeric column: {cfg.meta_numeric_col}")

    print(f"  Alpha diversity calculated for {len(alpha)} samples")
    print(f"  Metrics: {metrics}")
    print(alpha.head().to_string())

    outdir = os.path.join(cfg.output_dir, "result", "step3_alpha")
    os.makedirs(outdir, exist_ok=True)

    alpha_path = os.path.join(outdir, "alpha_diversity.xlsx")
    alpha.to_excel(alpha_path, index=False)
    print(f"  Saved: {alpha_path}")

    group_col = cfg.meta_group_col
    group_order = detect_group_order(alpha, meta, group_col, cfg.meta_numeric_col, cfg)
    cfg._group_order_detected = group_order

    n_total = len(alpha)

    # Kruskal-Wallis tests + effect size (epsilon-squared)
    raw_pvals = {}
    stat_results = []
    for metric in metrics:
        data = []
        for g in group_order:
            vals = alpha.loc[alpha[group_col] == g, metric].dropna()
            data.append(vals.values)

        if len(data) >= 2:
            if cfg.alpha_test == "anova":
                stat, p = f_oneway(*data)
                stat_name = "F_statistic"
                effect_size = stat
                effect_name = "eta_squared"
                n_groups = len(data)
                df_between = n_groups - 1
                df_within = n_total - n_groups
                if df_within > 0 and df_between > 0:
                    effect_size = stat * df_between / (stat * df_between + df_within)
            else:
                stat, p = kruskal(*data)
                stat_name = "Kruskal_H"
                effect_name = "epsilon_squared"
                effect_size = stat / (n_total - 1) if n_total > 1 else 0.0

            raw_pvals[metric] = p
            stat_results.append({
                "Metric": metric,
                stat_name: stat,
                "p_value": p,
                effect_name: effect_size,
                "n_groups": len(data),
            })

            sig = "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else "ns"
            print(f"  K-W {metric}: H={stat:.2f}, p={p:.4f} {sig}, {effect_name}={effect_size:.4f}")

    # BH FDR correction across all alpha tests
    if raw_pvals:
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

        for row in stat_results:
            row["p_fdr"] = fdr_map.get(row["Metric"], row["p_value"])

    stat_df = pd.DataFrame(stat_results)
    stat_path = os.path.join(outdir, "alpha_statistics.xlsx")
    stat_df.to_excel(stat_path, index=False)
    print(f"  Statistics saved: {stat_path}")
    print(stat_df.to_string())

    # Build stats dict for plotting
    plot_stats = {}
    for row in stat_results:
        metric = row["Metric"]
        if cfg.alpha_test == "anova":
            plot_stats[metric] = {
                "stat": row.get("F_statistic", 0),
                "p_raw": row["p_value"],
                "effect_size": row.get("eta_squared", 0),
                "effect_name": "η²",
                "p_fdr": row.get("p_fdr", row["p_value"]),
            }
        else:
            plot_stats[metric] = {
                "stat": row.get("Kruskal_H", 0),
                "p_raw": row["p_value"],
                "effect_size": row.get("epsilon_squared", 0),
                "effect_name": "ε²",
                "p_fdr": row.get("p_fdr", row["p_value"]),
            }

    if cfg.alpha_spearman and cfg.meta_numeric_col and cfg.meta_numeric_col in alpha.columns:
        spearman_results = []
        for metric in metrics:
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
            # FDR correction for Spearman tests too
            s_pvals = [r["p_value"] for r in spearman_results]
            s_fdr = bh_fdr(s_pvals)
            for r, f in zip(spearman_results, s_fdr):
                r["p_fdr"] = round(float(f), 4)

            spearman_df = pd.DataFrame(spearman_results)
            spearman_path = os.path.join(outdir, "alpha_spearman.xlsx")
            spearman_df.to_excel(spearman_path, index=False)
            print(f"  Spearman correlation saved: {spearman_path}")
            print(spearman_df.to_string())

    set_theme(cfg)
    fig_dir = os.path.join(outdir, "figures")
    os.makedirs(fig_dir, exist_ok=True)

    colors = cfg.get_group_colors()
    for metric in metrics:
        plot_alpha_boxplot(
            alpha, metric, group_col,
            group_order, colors,
            fig_dir, cfg,
            stats_info=plot_stats.get(metric),
        )
        print(f"  Plot saved: {metric}_boxplot")

    print("  Alpha diversity analysis complete.")
    return alpha_path
