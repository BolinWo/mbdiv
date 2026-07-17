"""
Pipeline configuration module.
All configurable parameters live here with sensible defaults.
"""

from dataclasses import dataclass, field
from typing import List, Dict


@dataclass
class PipelineConfig:
    """Configuration for the microbiome diversity analysis pipeline."""

    # ---- Input files ----
    raw_data: str = ""          # Path to raw species abundance Excel
    meta_data: str = ""         # Path to metadata Excel
    output_dir: str = "result"  # Base output directory

    # ---- Column names (auto-detected if empty) ----
    species_col: str = ""       # e.g. "#Species"
    tax_col: str = ""           # e.g. "#Taxonomy"
    sample_prefix: str = ""     # e.g. "P" — if set, sample cols = those starting with prefix
    sample_col_keyword: str = ""  # alternative: any keyword to match sample columns

    # ---- Metadata ----
    meta_sample_col: str = "Sample"    # column name for sample IDs in metadata
    meta_group_col: str = "Group"      # column name for grouping variable
    meta_numeric_col: str = ""         # optional numeric column for Spearman correlation

    # ---- Step 1: Merge ----
    # Regex pattern to extract clean species name (remove [TAX_xxx] suffix etc.)
    species_name_split: str = "["      # species name = part before this char
    # Whether to keep taxonomy column
    keep_taxonomy: bool = True

    # ---- Step 2: Normalization ----
    # Output relative abundance as fraction (False) or percentage (True)
    as_percent: bool = True

    # ---- Step 3: Alpha diversity ----
    alpha_metrics: List[str] = field(default_factory=lambda: [
        "Observed", "Shannon", "Simpson", "Chao1"
    ])
    # Statistical test for group comparison: "kruskal" or "anova"
    alpha_test: str = "kruskal"
    # Whether to compute Spearman correlation with numeric metadata column
    alpha_spearman: bool = True

    # ---- Step 4: Beta diversity ----
    beta_distance: str = "braycurtis"   # scipy pdist metric name
    beta_permanova_permutations: int = 999
    # Draw confidence ellipses on PCoA (requires >= 3 samples per group)
    pcoa_ellipse: bool = True

    # ---- Step 5: Top10 species ----
    # Keywords to filter OUT (case-insensitive)
    virus_keywords: List[str] = field(default_factory=lambda: [
        "Viruses", "Duplodnaviria", "Caudoviricetes",
        "phage", "bacteriophage",
    ])
    fungi_keywords: List[str] = field(default_factory=lambda: [
        "Fungi", "Ascomycota", "Basidiomycota",
        "Saccharomycetes", "Zygomycota", "Chytridiomycota",
    ])
    # Also filter by domain prefix in taxonomy string
    filter_domains: List[str] = field(default_factory=lambda: [
        "d__Viruses", "k__Viruses",
    ])
    top_n: int = 10

    # ---- Plotting ----
    # Group display order (auto-detected from metadata if empty)
    group_order: List[str] = field(default_factory=list)
    # Group colors (auto-assigned if group not in dict)
    group_colors: Dict[str, str] = field(default_factory=lambda: {
        "Control": "#55A868",
        "Treatment": "#4C72B0",
        "Group_C": "#DD8452",
        "Group_D": "#C44E52",
    })
    # Species color palette for stacked bar
    species_palette: List[str] = field(default_factory=lambda: [
        "#4C78A8", "#F58518", "#54A24B", "#E45756", "#72B7B2",
        "#B279A2", "#FF9DA6", "#9D755D", "#BAB0AC", "#79706E",
        "#D9D9D9",
    ])
    # Figure format: "pdf", "png", or "both"
    fig_format: str = "both"
    fig_dpi: int = 300
    savefig_dpi: int = 600

    # ---- Internal (set during runtime) ----
    _species_col_detected: str = ""
    _tax_col_detected: str = ""
    _sample_cols_detected: List[str] = field(default_factory=list)
    _group_order_detected: List[str] = field(default_factory=list)

    def get_group_order(self) -> List[str]:
        """Return configured or detected group order."""
        if self.group_order:
            return self.group_order
        return self._group_order_detected

    def get_group_colors(self) -> Dict[str, str]:
        """Return colors, auto-assigning missing groups."""
        colors = dict(self.group_colors)
        palette = [
            "#55A868", "#4C72B0", "#DD8452", "#C44E52",
            "#8172B3", "#937860", "#DA8BC3", "#8C8C8C",
            "#CCB974", "#64B5CD",
        ]
        order = self.get_group_order()
        for i, g in enumerate(order):
            if g not in colors:
                colors[g] = palette[i % len(palette)]
        return colors
