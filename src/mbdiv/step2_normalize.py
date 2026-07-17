"""
Step 2: Normalize abundance to relative abundance.

Input:  merged_species_clean.xlsx
Output: relative_abundance.xlsx, relative_abundance_percent.xlsx, sample_rpkm_summary.xlsx
"""

import os
import numpy as np
import pandas as pd
from .config import PipelineConfig


def run_step2(cfg: PipelineConfig, merged_path: str = None) -> str:
    """
    Execute Step 2: relative abundance normalization.
    Returns path to relative_abundance.xlsx (or _percent version).
    """
    print("\n" + "=" * 60)
    print("Step 2: Normalization (relative abundance)")
    print("=" * 60)

    if merged_path is None:
        merged_path = os.path.join(cfg.output_dir, "data", "merged_species_clean.xlsx")

    df = pd.read_excel(merged_path)
    print(f"  Input shape: {df.shape}")

    species_col = "Species_name"
    sample_cols = cfg._sample_cols_detected

    if not sample_cols:
        # Re-detect
        from .step1_merge import detect_columns
        detect_columns(df, cfg)
        sample_cols = cfg._sample_cols_detected

    # Abundance matrix
    abundance = df[[species_col] + sample_cols].copy()
    sample_sum = abundance[sample_cols].sum(axis=0)

    # Guard: samples with zero total abundance
    zero_samples = sample_sum[sample_sum == 0].index.tolist()
    if zero_samples:
        print(f"  WARNING: {len(zero_samples)} samples have zero total abundance: {zero_samples[:5]}")
        print("  These will be excluded from normalization.")
        sample_cols = [s for s in sample_cols if s not in zero_samples]
        abundance = abundance[sample_cols]
        sample_sum = sample_sum[sample_cols]

    print("  Total reads per sample (first 5):")
    for s, v in sample_sum.head(5).items():
        print(f"    {s}: {v:.2f}")

    # Relative abundance
    relative = abundance.copy()
    relative[sample_cols] = relative[sample_cols] / sample_sum

    # Verify normalization
    check = relative[sample_cols].sum(axis=0)
    assert np.allclose(check, 1.0), f"Normalization check failed: {check}"
    print("  Normalization check passed (all columns sum to 1.0)")

    outdir = os.path.join(cfg.output_dir, "data")
    os.makedirs(outdir, exist_ok=True)

    # Save fractional relative abundance
    frac_path = os.path.join(outdir, "relative_abundance.xlsx")
    relative.to_excel(frac_path, index=False)
    print(f"  Saved: {frac_path}")

    # Save percentage version
    if cfg.as_percent:
        percent = relative.copy()
        percent[sample_cols] = percent[sample_cols] * 100
        pct_path = os.path.join(outdir, "relative_abundance_percent.xlsx")
        percent.to_excel(pct_path, index=False)
        print(f"  Saved: {pct_path}")

    # Save sequencing depth summary
    summary = pd.DataFrame({
        "Sample": sample_sum.index,
        "Total_RPKM": sample_sum.values,
    })
    summary_path = os.path.join(outdir, "sample_rpkm_summary.xlsx")
    summary.to_excel(summary_path, index=False)
    print(f"  Saved: {summary_path}")

    return frac_path
