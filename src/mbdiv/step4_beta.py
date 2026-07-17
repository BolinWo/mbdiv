"""
Step 4: Beta diversity analysis.

- Bray-Curtis distance matrix
- PCoA ordination
- PERMANOVA test
- PCoA scatter plot with confidence ellipses

Input:  relative_abundance.xlsx, metadata.xlsx
Output: result/step4_beta/distance/bray_curtis_distance.xlsx,
        result/step4_beta/pcoa/pcoa_coordinates.xlsx, pcoa_variance.xlsx,
        result/step4_beta/statistics/PERMANOVA_result.txt,
        result/step4_beta/figures/PCoA_BrayCurtis.{pdf,png}
"""

import os
import numpy as np
import pandas as pd
from scipy.spatial.distance import pdist, squareform

from .config import PipelineConfig
from .theme import set_theme, plot_pcoa


def run_step4(cfg: PipelineConfig, rel_path: str = None, meta_path: str = None):
    """
    Execute Step 4: beta diversity.
    Returns (distance_path, pcoa_coords_path, permanova_p_value).
    """
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

    # --- Load data ---
    otu = pd.read_excel(rel_path, index_col=0)
    meta = pd.read_excel(meta_path)

    # Transpose: samples x species
    otu_t = otu.T
    print(f"  Distance matrix input: {otu_t.shape[0]} samples x {otu_t.shape[1]} species")

    # --- Bray-Curtis distance ---
    distance = pd.DataFrame(
        squareform(pdist(otu_t, metric=cfg.beta_distance)),
        index=otu_t.index,
        columns=otu_t.index,
    )
    dist_path = os.path.join(dist_dir, f"{cfg.beta_distance}_distance.xlsx")
    distance.to_excel(dist_path)
    print(f"  Distance matrix saved: {dist_path}")

    # --- PCoA ---
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

    # Save PCoA results
    coords_path = os.path.join(pcoa_dir, "pcoa_coordinates.xlsx")
    coords.to_excel(coords_path, index=False)
    print(f"  PCoA coordinates saved: {coords_path}")

    var_path = os.path.join(pcoa_dir, "pcoa_variance.xlsx")
    variance.to_excel(var_path, index=False)
    print(f"  PCoA variance saved: {var_path}")

    # --- PERMANOVA ---
    permanova_p = None
    try:
        from skbio.stats.distance import permanova

        meta_indexed = meta.set_index(cfg.meta_sample_col)
        samples = list(distance.index)
        meta_matched = meta_indexed.loc[samples]

        dm = DistanceMatrix(
            np.ascontiguousarray(distance.values, dtype=float),
            ids=samples,
        )

        permanova_result = permanova(
            dm, meta_matched,
            column=cfg.meta_group_col,
            permutations=cfg.beta_permanova_permutations,
        )
        permanova_p = permanova_result.get("p-value", None)

        result_text = str(permanova_result)
        with open(os.path.join(stat_dir, "PERMANOVA_result.txt"), "w") as f:
            f.write(result_text)
        print(f"  PERMANOVA p-value: {permanova_p}")
        print(f"  PERMANOVA result saved: {os.path.join(stat_dir, 'PERMANOVA_result.txt')}")

    except ImportError:
        print("  scikit-bio not available, skipping PERMANOVA")
        with open(os.path.join(stat_dir, "PERMANOVA_result.txt"), "w") as f:
            f.write("PERMANOVA: scikit-bio not available\n")

    except Exception as e:
        print(f"  PERMANOVA error: {e}")
        with open(os.path.join(stat_dir, "PERMANOVA_result.txt"), "w") as f:
            f.write(f"PERMANOVA error: {e}\n")

    # --- PCoA plot ---
    set_theme(cfg)
    group_order = cfg.get_group_order()
    colors = cfg.get_group_colors()

    plot_pcoa(
        coords, variance, meta,
        cfg.meta_group_col, cfg.meta_sample_col,
        group_order, colors,
        permanova_p,
        fig_dir, cfg,
    )
    print(f"  PCoA plot saved to: {fig_dir}")

    print("  Beta diversity analysis complete.")
    return dist_path, coords_path, permanova_p


def _manual_pcoa(distance: pd.DataFrame, cfg: PipelineConfig):
    """Manual PCoA via eigendecomposition of the distance matrix."""
    n = distance.shape[0]
    D = distance.values.astype(float)

    # Double-centering: G = -0.5 * J * D^2 * J
    D_sq = D ** 2
    J = np.eye(n) - np.ones((n, n)) / n
    G = -0.5 * J @ D_sq @ J

    # Eigendecomposition
    eigenvalues, eigenvectors = np.linalg.eigh(G)

    # Sort descending
    idx = np.argsort(eigenvalues)[::-1]
    eigenvalues = eigenvalues[idx]
    eigenvectors = eigenvectors[:, idx]

    # Coordinates
    coords = eigenvectors[:, :2] * np.sqrt(np.maximum(eigenvalues[:2], 0))

    # Proportion explained
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
