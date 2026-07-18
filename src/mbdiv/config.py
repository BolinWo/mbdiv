"""Pipeline configuration with sensible defaults."""

from dataclasses import dataclass, field
from typing import List, Dict


@dataclass
class PipelineConfig:
    """All configurable parameters for the analysis pipeline."""

    # input files
    raw_data: str = ""
    meta_data: str = ""
    output_dir: str = "result"

    # column names - left empty for auto-detection
    species_col: str = ""
    tax_col: str = ""
    sample_prefix: str = ""
    sample_col_keyword: str = ""

    # metadata
    meta_sample_col: str = "Sample"
    meta_group_col: str = "Group"
    meta_numeric_col: str = ""  # for Spearman correlation

    # step 1: merge
    species_name_split: str = "["  # species name = text before this char
    keep_taxonomy: bool = True

    # step 2: normalize
    as_percent: bool = True

    # step 3: alpha diversity
    alpha_metrics: List[str] = field(default_factory=lambda: [
        "Observed", "Shannon", "Simpson", "Chao1"
    ])
    alpha_test: str = "kruskal"  # "kruskal" or "anova"
    alpha_spearman: bool = True

    # step 4: beta diversity
    beta_distance: str = "braycurtis"
    beta_permanova_permutations: int = 999
    pcoa_engine: str = "manual"  # "manual" or "skbio"
    pcoa_ellipse: bool = True  # needs >= 3 samples per group

    # step 5: top-N
    filter_viruses_fungi: bool = True
    virus_keywords: List[str] = field(default_factory=lambda: [
        "Viruses", "Duplodnaviria", "Caudoviricetes",
        "phage", "bacteriophage",
    ])
    fungi_keywords: List[str] = field(default_factory=lambda: [
        "Fungi", "Ascomycota", "Basidiomycota",
        "Saccharomycetes", "Zygomycota", "Chytridiomycota",
    ])
    filter_domains: List[str] = field(default_factory=lambda: [
        "d__Viruses", "k__Viruses",
    ])
    top_n: int = 10
    top_n_mode: str = "group_union"  # "overall_mean" or "group_union"

    # plotting
    group_order: List[str] = field(default_factory=list)
    group_colors: Dict[str, str] = field(default_factory=lambda: {
        "Control": "#55A868",
        "Treatment": "#4C72B0",
        "Group_C": "#DD8452",
        "Group_D": "#C44E52",
    })
    species_palette: List[str] = field(default_factory=lambda: [
        "#4C78A8", "#F58518", "#54A24B", "#E45756", "#72B7B2",
        "#B279A2", "#FF9DA6", "#9D755D", "#BAB0AC", "#79706E",
        "#D9D9D9",
    ])
    fig_format: str = "both"  # "pdf", "png", or "both"
    fig_dpi: int = 300
    savefig_dpi: int = 600

    # runtime detection results
    _species_col_detected: str = ""
    _tax_col_detected: str = ""
    _sample_cols_detected: List[str] = field(default_factory=list)
    _group_order_detected: List[str] = field(default_factory=list)

    def get_group_order(self) -> List[str]:
        return self.group_order or self._group_order_detected

    def get_group_colors(self) -> Dict[str, str]:
        colors = dict(self.group_colors)
        palette = [
            "#55A868", "#4C72B0", "#DD8452", "#C44E52",
            "#8172B3", "#937860", "#DA8BC3", "#8C8C8C",
            "#CCB974", "#64B5CD",
        ]
        for i, g in enumerate(self.get_group_order()):
            if g not in colors:
                colors[g] = palette[i % len(palette)]
        return colors
