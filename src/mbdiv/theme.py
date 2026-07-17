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
    plt.tight_layout()

    os.makedirs(output_dir, exist_ok=True)
    save_figure(fig, os.path.join(output_dir, f"{metric}_boxplot"), cfg)
    plt.close(fig)


def plot_pcoa(
    pcoa_df, variance_df, meta_df,
    group_col: str, sample_col: str,
    group_order: List[str], colors: Dict[str, str],
    permanova_p: Optional[float],
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

    if permanova_p is not None:
        ax.text(
            0.95, 0.95,
            f"PERMANOVA\np = {permanova_p:.3f}",
            transform=ax.transAxes,
            ha="right", va="top", fontsize=11,
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
