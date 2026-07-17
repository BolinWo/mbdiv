"""
Step 1: Merge same-name species and remove all-zero abundance taxa.

Input:  raw abundance Excel with #Species / Species_name and #Taxonomy columns
Output: merged_species_clean.xlsx, zero_abundance_taxa.xlsx
"""

import os
import pandas as pd
from .config import PipelineConfig


def detect_columns(df: pd.DataFrame, cfg: PipelineConfig):
    """Auto-detect species, taxonomy, and sample columns from the dataframe."""
    cols = df.columns.tolist()

    # --- Species column ---
    if cfg.species_col and cfg.species_col in cols:
        species_col = cfg.species_col
    else:
        candidates = [c for c in cols if "species" in c.lower() or c.startswith("#")]
        species_col = candidates[0] if candidates else cols[0]
    cfg._species_col_detected = species_col

    # --- Taxonomy column ---
    if cfg.tax_col and cfg.tax_col in cols:
        tax_col = cfg.tax_col
    else:
        candidates = [c for c in cols if "tax" in c.lower()]
        tax_col = candidates[0] if candidates else None
    cfg._tax_col_detected = tax_col or ""

    # --- Sample columns ---
    non_sample = {species_col}
    if tax_col:
        non_sample.add(tax_col)

    if cfg.sample_prefix:
        sample_cols = [c for c in cols if c.startswith(cfg.sample_prefix) and c not in non_sample]
    elif cfg.sample_col_keyword:
        sample_cols = [c for c in cols if cfg.sample_col_keyword in c and c not in non_sample]
    else:
        # Auto-detect: numeric columns that are not species/taxonomy
        sample_cols = []
        for c in cols:
            if c in non_sample:
                continue
            if pd.api.types.is_numeric_dtype(df[c]):
                sample_cols.append(c)

    if not sample_cols:
        # Fallback: all non-meta columns
        sample_cols = [c for c in cols if c not in non_sample]

    cfg._sample_cols_detected = sample_cols
    return species_col, tax_col, sample_cols


def extract_species_name(raw_name: str, split_char: str) -> str:
    """Extract clean species name (e.g. 'Bacteriophage_sp. [TAX_009734245.1]' -> 'Bacteriophage_sp.')."""
    return str(raw_name).split(split_char)[0].strip()


def run_step1(cfg: PipelineConfig) -> str:
    """
    Execute Step 1: merge species + remove zeros.
    Returns path to merged_species_clean.xlsx.
    """
    print("\n" + "=" * 60)
    print("Step 1: Merge same-name species & remove zero-abundance taxa")
    print("=" * 60)

    df = pd.read_excel(cfg.raw_data)
    print(f"  Input shape: {df.shape}")

    species_col, tax_col, sample_cols = detect_columns(df, cfg)
    print(f"  Species column: {species_col}")
    print(f"  Taxonomy column: {tax_col}")
    print(f"  Sample columns ({len(sample_cols)}): {sample_cols[:5]}...")

    # Extract clean species name
    df["Species_name"] = df[species_col].apply(
        lambda x: extract_species_name(x, cfg.species_name_split)
    )

    # Merge by species name (sum abundance)
    merged = df.groupby("Species_name")[sample_cols].sum().reset_index()

    # Attach taxonomy (first occurrence)
    if tax_col:
        taxonomy = df[["Species_name", tax_col]].drop_duplicates("Species_name")
        merged = merged.merge(taxonomy, on="Species_name", how="left")
    else:
        tax_col = "#Taxonomy"
        merged["#Taxonomy"] = ""

    # Reorder columns: Species_name, samples..., Taxonomy
    merged = merged[["Species_name"] + sample_cols + [tax_col]]

    # Identify and save zero-abundance taxa
    zero_mask = merged[sample_cols].sum(axis=1) == 0
    zero_taxa = merged[zero_mask].copy()

    outdir = os.path.join(cfg.output_dir, "data")
    os.makedirs(outdir, exist_ok=True)

    zero_path = os.path.join(outdir, "zero_abundance_taxa.xlsx")
    zero_taxa.to_excel(zero_path, index=False)

    before = len(merged)
    merged = merged[~zero_mask].copy()
    after = len(merged)
    print(f"  Zero-abundance taxa removed: {before - after}")
    print(f"  Remaining species: {after}")

    # Save merged result
    merged_path = os.path.join(outdir, "merged_species_clean.xlsx")
    merged.to_excel(merged_path, index=False)
    print(f"  Saved: {merged_path}")

    return merged_path
