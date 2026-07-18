"""Step 5: filter viruses/fungi, pick top-N species, draw stacked bar + heatmap."""

import os
import pandas as pd

from .config import PipelineConfig
from .theme import set_theme, plot_top10_stacked


def run_step5(cfg: PipelineConfig, merged_path: str = None, meta_path: str = None):
    print("\n" + "=" * 60)
    print("Step 5: Filter viruses/fungi + Top-N species stacked bar")
    print("=" * 60)

    if merged_path is None:
        merged_path = os.path.join(cfg.output_dir, "data", "merged_species_clean.xlsx")
    if meta_path is None:
        meta_path = cfg.meta_data

    base_dir = os.path.join(cfg.output_dir, "result", "step5")
    fig_dir = os.path.join(base_dir, "figures")
    os.makedirs(fig_dir, exist_ok=True)

    df = pd.read_excel(merged_path)
    meta = pd.read_excel(meta_path)
    print(f"  Input shape: {df.shape}")

    species_col = "Species_name"
    tax_col = cfg._tax_col_detected or "#Taxonomy"
    sample_cols = cfg._sample_cols_detected
    if not sample_cols:
        sample_cols = [c for c in df.columns if c != species_col and c != tax_col]

    # filter out viruses and fungi by taxonomy keywords
    if cfg.filter_viruses_fungi:
        taxonomy = df[tax_col].fillna("").astype(str)
        remove_mask = pd.Series(False, index=df.index)

        for key in cfg.virus_keywords + cfg.fungi_keywords:
            remove_mask |= taxonomy.str.contains(key, case=False, regex=False)

        for domain in cfg.filter_domains:
            remove_mask |= taxonomy.str.contains(domain, case=False, regex=False)

        df_clean = df.loc[~remove_mask].copy()
        print(f"  Removed {remove_mask.sum()} viral/fungal taxa")
        print(f"  Remaining: {df_clean.shape[0]} species")
    else:
        df_clean = df.copy()
        print("  Virus/fungi filtering disabled")

    zero_mask = df_clean[sample_cols].sum(axis=1) == 0
    df_clean = df_clean[~zero_mask].copy()
    print(f"  After zero removal: {df_clean.shape[0]} species")

    bacteria_path = os.path.join(base_dir, "bacteria_species.xlsx")
    df_clean.to_excel(bacteria_path, index=False)
    print(f"  Saved: {bacteria_path}")

    # recalculate relative abundance on filtered data
    abundance = df_clean[sample_cols].copy()
    col_sums = abundance.sum(axis=0)

    zero_cols = col_sums[col_sums == 0].index.tolist()
    if zero_cols:
        print(f"  WARNING: {len(zero_cols)} samples have zero total after filtering, excluded.")
        sample_cols = [s for s in sample_cols if s not in zero_cols]
        abundance = abundance[sample_cols]
        col_sums = col_sums[sample_cols]

    relative = abundance.div(col_sums, axis=1) * 100

    result = pd.concat(
        [df_clean[[species_col]], relative, df_clean[[tax_col]]],
        axis=1,
    )

    bacteria_rel_path = os.path.join(base_dir, "bacteria_relative_abundance.xlsx")
    result.to_excel(bacteria_rel_path, index=False)
    print(f"  Saved: {bacteria_rel_path}")

    # top-N selection
    result["Mean_abundance"] = result[sample_cols].mean(axis=1)

    group_order = cfg.get_group_order()
    if not group_order:
        group_order = sorted(meta[cfg.meta_group_col].dropna().unique().tolist())

    if cfg.top_n_mode == "group_union":
        # Each group selects its own Top-N, then take the union
        print(f"\n  Top-N mode: group_union (each group Top-{cfg.top_n}, then union)")
        meta_indexed = meta.set_index(cfg.meta_sample_col)
        sample_to_group = meta_indexed[cfg.meta_group_col].to_dict()

        group_samples = {}
        for s in sample_cols:
            g = sample_to_group.get(s)
            if g is not None:
                group_samples.setdefault(g, []).append(s)

        all_top = set()
        for g in group_order:
            g_samples = group_samples.get(g, [])
            if g_samples:
                g_means = result[g_samples].mean(axis=1)
                g_top = g_means.sort_values(ascending=False).head(cfg.top_n)
                top_names = result.loc[g_top.index, species_col].tolist()
                all_top.update(top_names)
                print(f"    {g}: {len(top_names)} species selected (n={len(g_samples)} samples)")

        # Sort union by overall mean abundance
        top_species = (
            result[result[species_col].isin(all_top)]
            .sort_values("Mean_abundance", ascending=False)[species_col]
            .tolist()
        )
        print(f"  Union total: {len(top_species)} species")
    else:
        # Overall mean Top-N (default)
        print(f"\n  Top-N mode: overall_mean (Top-{cfg.top_n} by mean across all samples)")
        top_n_df = result.sort_values("Mean_abundance", ascending=False).head(cfg.top_n).copy()
        top_species = top_n_df[species_col].tolist()

    top_n_df = result[result[species_col].isin(top_species)].sort_values("Mean_abundance", ascending=False).copy()

    print(f"\n  Selected species ({len(top_species)}):")
    for i, s in enumerate(top_species, 1):
        row = top_n_df[top_n_df[species_col] == s].iloc[0]
        print(f"    {i:2d}. {s} (mean={row['Mean_abundance']:.2f}%)")

    top10_path = os.path.join(base_dir, f"top{cfg.top_n}_species.xlsx")
    top_n_df.to_excel(top10_path, index=False)
    print(f"  Saved: {top10_path}")

    # group-level composition
    composition = result[[species_col] + sample_cols].set_index(species_col).T
    composition.index.name = cfg.meta_sample_col

    comp_top = composition[top_species].copy()
    comp_top["Others"] = composition.drop(columns=top_species).sum(axis=1)

    meta_indexed = meta.set_index(cfg.meta_sample_col)
    comp_top = comp_top.join(meta_indexed[[cfg.meta_group_col]])

    group_summary = comp_top.groupby(cfg.meta_group_col).mean()

    group_order = [g for g in group_order if g in group_summary.index]
    group_summary = group_summary.loc[group_order]

    pct_path = os.path.join(base_dir, f"top{cfg.top_n}_group_percentage.xlsx")
    group_summary.to_excel(pct_path)
    print(f"  Group percentages saved: {pct_path}")

    # stacked bar
    set_theme(cfg)

    species_order = (
        group_summary.drop(columns="Others")
        .mean()
        .sort_values(ascending=False)
        .index.tolist()
    )
    species_order.append("Others")

    plot_top10_stacked(
        group_summary, species_order,
        fig_dir, cfg,
        filename=f"top{cfg.top_n}_group_bar",
    )
    print(f"  Stacked bar plot saved to: {fig_dir}")

    # per-sample heatmap
    try:
        import matplotlib.pyplot as plt
        import seaborn as sns

        individual_data = comp_top.drop(columns=cfg.meta_group_col)[species_order]

        fig, ax = plt.subplots(figsize=(10, 6), dpi=cfg.fig_dpi)
        sns.heatmap(
            individual_data.T,
            cmap="viridis",
            ax=ax,
            cbar_kws={"label": "Relative abundance (%)"},
            linewidths=0.3,
            linecolor="white",
        )
        ax.set_xlabel("")
        ax.set_ylabel("")
        ax.set_title("Top-{} species composition (individual samples)".format(cfg.top_n),
                      fontsize=13, weight="bold")
        plt.xticks(rotation=45, ha="right", fontsize=8)
        plt.yticks(fontsize=9)
        plt.tight_layout()

        from .theme import save_figure
        save_figure(fig, os.path.join(fig_dir, f"top{cfg.top_n}_individual_heatmap"), cfg)
        plt.close(fig)
        print("  Individual heatmap saved")
    except Exception as e:
        print(f"  Heatmap skipped: {e}")

    print("  Step 5 complete.")
    return bacteria_rel_path, top10_path
