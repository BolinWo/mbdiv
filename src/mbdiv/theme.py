"""Shared plotting functions."""

import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from matplotlib.patches import Ellipse, Patch
from scipy.stats import chi2
from typing import Dict, List, Optional
from .config import PipelineConfig


def set_theme(cfg: PipelineConfig):
    plt.rcParams.update({
        "font.family": "Arial",
        "font.size": 12,
        "axes.linewidth": 1,
        "axes.grid": False,
        "figure.dpi": cfg.fig_dpi,
        "savefig.dpi": cfg.savefig_dpi,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    })
    sns.set_theme(style="white", context="paper")


def save_figure(fig, filepath_no_ext: str, cfg: PipelineConfig):
    fmt = cfg.fig_format.lower()
    if fmt in ("pdf", "both"):
        fig.savefig(f"{filepath_no_ext}.pdf", dpi=cfg.savefig_dpi, bbox_inches="tight")
    if fmt in ("png", "both"):
        fig.savefig(f"{filepath_no_ext}.png", dpi=cfg.savefig_dpi, bbox_inches="tight")


def confidence_ellipse(x, y, ax, color, alpha=0.12):
    if len(x) < 3:
        return
    cov = np.cov(x, y)
    mean_x, mean_y = np.mean(x), np.mean(y)
    eigenvalues, eigenvectors = np.linalg.eigh(cov)
    order = eigenvalues.argsort()[::-1]
    eigenvalues = eigenvalues[order]
    eigenvectors = eigenvectors[:, order]
    theta = np.degrees(np.arctan2(eigenvectors[1, 0], eigenvectors[0, 0]))
    width, height = 2 * np.sqrt(chi2.ppf(0.95, 2) * eigenvalues)
    ellipse = Ellipse(
        xy=(mean_x, mean_y),
        width=width, height=height, angle=theta,
        facecolor=color, edgecolor="none",
        alpha=alpha, linewidth=1.5,
    )
    ax.add_patch(ellipse)


def plot_alpha_boxplot(
    alpha_df, metric: str, group_col: str,
    group_order: List[str], colors: Dict[str, str],
    output_dir: str, cfg: PipelineConfig,
    stats_info: Optional[dict] = None,
):
    fig, ax = plt.subplots(figsize=(5.5, 4.5), dpi=cfg.fig_dpi)

    sns.boxplot(
        data=alpha_df, x=group_col, y=metric,
        order=group_order, hue=group_col, hue_order=group_order,
        palette=colors, legend=False,
        width=0.55, linewidth=1.2, fliersize=0,
        boxprops={"alpha": 0.75}, ax=ax,
    )
    sns.stripplot(
        data=alpha_df, x=group_col, y=metric,
        order=group_order, color="black",
        size=5, jitter=0.15, alpha=0.8, ax=ax,
    )

    ax.set_xlabel("")
    ax.set_ylabel(metric, fontsize=13, weight="bold")
    ax.set_title(metric, fontsize=15, weight="bold")
    ax.tick_params(labelsize=11)
    sns.despine()

    if stats_info:
        stat_val = stats_info.get("stat", 0)
        p_raw = stats_info.get("p_raw", 1)
        effect_size = stats_info.get("effect_size", 0)
        effect_name = stats_info.get("effect_name", "")
        p_fdr = stats_info.get("p_fdr", p_raw)

        sig = "***" if p_fdr < 0.001 else "**" if p_fdr < 0.01 else "*" if p_fdr < 0.05 else "ns"
        clr = "#C0392B" if p_fdr < 0.05 else "#95A5A6"

        stat_label = "H" if effect_name == "ε²" else "F"
        text_str = (
            f"{stat_label} = {stat_val:.2f}\n"
            f"{effect_name} = {effect_size:.4f}\n"
            f"p = {p_raw:.4f}\n"
            f"p_adj = {p_fdr:.4f} {sig}"
        )
        ax.text(
            0.97, 0.97, text_str,
            transform=ax.transAxes, ha="right", va="top", fontsize=8.5,
            fontweight="bold", color=clr,
            bbox=dict(boxstyle="round,pad=0.35", fc="white", ec=clr, alpha=0.92, lw=1),
        )

    plt.tight_layout()

    os.makedirs(output_dir, exist_ok=True)
    save_figure(fig, os.path.join(output_dir, f"{metric}_boxplot"), cfg)
    plt.close(fig)


def plot_pcoa(
    pcoa_df, variance_df, meta_df,
    group_col: str, sample_col: str,
    group_order: List[str], colors: Dict[str, str],
    permanova_stats: Optional[dict],
    permdisp_stats: Optional[dict],
    output_dir: str, cfg: PipelineConfig,
):
    df = pcoa_df.merge(meta_df, on=sample_col, how="left")

    var1 = variance_df.iloc[0]["Explained_variance"] if len(variance_df) > 0 else 0
    var2 = variance_df.iloc[1]["Explained_variance"] if len(variance_df) > 1 else 0
    pc1_label = f"PC1 ({var1*100:.1f}%)"
    pc2_label = f"PC2 ({var2*100:.1f}%)"

    fig, ax = plt.subplots(figsize=(8, 7), dpi=cfg.fig_dpi)

    for g in group_order:
        sub = df[df[group_col] == g]
        if len(sub) == 0:
            continue
        if cfg.pcoa_ellipse and len(sub) >= 3:
            confidence_ellipse(sub["PC1"].values, sub["PC2"].values, ax, colors.get(g, "#888888"))
        ax.scatter(
            sub["PC1"], sub["PC2"],
            s=70, color=colors.get(g, "#888888"),
            edgecolor="black", linewidth=0.7,
            alpha=0.9, label=g,
        )

    ax.set_xlabel(pc1_label, fontsize=14, weight="bold")
    ax.set_ylabel(pc2_label, fontsize=14, weight="bold")
    ax.legend(frameon=False, fontsize=12)

    if permanova_stats:
        perm_F = permanova_stats.get("F", 0)
        perm_p = permanova_stats.get("p", 1)
        perm_p_fdr = permanova_stats.get("p_fdr", perm_p)
        perm_R2 = permanova_stats.get("R2", 0)
        perm_sig = "***" if perm_p_fdr < 0.001 else "**" if perm_p_fdr < 0.01 else "*" if perm_p_fdr < 0.05 else "ns"
        perm_clr = "#C0392B" if perm_p_fdr < 0.05 else "#95A5A6"

        text_lines = [
            f"PERMANOVA  F = {perm_F:.2f}  R² = {perm_R2:.4f}",
            f"p = {perm_p:.4f}  p_adj = {perm_p_fdr:.4f} {perm_sig}",
        ]

        if permdisp_stats:
            disp_F = permdisp_stats.get("F", 0)
            disp_p = permdisp_stats.get("p", 1)
            disp_p_fdr = permdisp_stats.get("p_fdr", disp_p)
            disp_sig = "***" if disp_p_fdr < 0.001 else "**" if disp_p_fdr < 0.01 else "*" if disp_p_fdr < 0.05 else "ns"
            text_lines.append(f"PERMDISP   F = {disp_F:.2f}  p = {disp_p:.4f}  p_adj = {disp_p_fdr:.4f} {disp_sig}")

        ax.text(
            0.95, 0.95,
            "\n".join(text_lines),
            transform=ax.transAxes,
            ha="right", va="top", fontsize=9,
            fontweight="bold", color=perm_clr,
            bbox=dict(boxstyle="round,pad=0.35", fc="white", ec=perm_clr, alpha=0.92, lw=1),
        )

    sns.despine(top=True, right=True)
    ax.tick_params(axis="both", which="major", labelsize=12, width=1)
    plt.tight_layout()

    os.makedirs(output_dir, exist_ok=True)
    save_figure(fig, os.path.join(output_dir, "PCoA_BrayCurtis"), cfg)
    plt.close(fig)


def plot_top10_stacked(
    group_data, species_order: List[str],
    output_dir: str, cfg: PipelineConfig,
    filename: str = "top10_group_bar",
):
    n_species = len(species_order)
    palette = cfg.species_palette
    colors = [palette[i % len(palette)] for i in range(n_species)]

    fig, ax = plt.subplots(figsize=(6.5, 5.5), dpi=cfg.fig_dpi)

    plot_data = group_data[species_order]
    bottom = np.zeros(len(plot_data))

    for i, sp in enumerate(species_order):
        ax.bar(
            plot_data.index, plot_data[sp],
            bottom=bottom, width=0.55,
            color=colors[i], linewidth=0.0,
        )
        bottom += plot_data[sp].values

    ax.set_ylim(0, 100)
    ax.set_ylabel("Relative abundance (%)", fontsize=13, fontweight="bold")
    ax.tick_params(axis="x", labelsize=12)
    ax.tick_params(axis="y", labelsize=11)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    legend_elements = [
        Patch(facecolor=colors[i], label=sp)
        for i, sp in enumerate(species_order)
    ]
    ax.legend(
        handles=legend_elements,
        bbox_to_anchor=(1.02, 1), loc="upper left",
        frameon=False, fontsize=8,
    )

    plt.tight_layout()
    os.makedirs(output_dir, exist_ok=True)
    save_figure(fig, os.path.join(output_dir, filename), cfg)
    plt.close(fig)
