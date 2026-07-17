# Pipeline Design Notes

## Original Project Decomposition

The source project (`CSVD_analysis/`) was a step-by-step manual pipeline
broken into 5 folders with 14 individual scripts. Each step read and wrote
Excel files from disk, with hardcoded paths and group names throughout.

### Issues in the original

1. **Hardcoded paths in step 5** — `r"E:/metadata/CSVD_analysis/..."` made
   scripts non-portable.
2. **Hardcoded group names** — `"Normal", "Low", "Medium", "Heavy"` repeated
   across step3 and step4.
3. **Hardcoded PCoA axis labels** — `"PC1 (25.1%)"`, `"PC2 (20.4%)"` were
   written as text instead of being read from the analysis output.
4. **Hardcoded PERMANOVA p-value** — `"p = 0.435"` was typed into the plot
   code rather than read from `PERMANOVA_result.txt`.
5. **Fragile sample detection** — `c.startswith("P")` failed when sample
   names used different prefixes.
6. **No config system** — every parameter had to be edited in source code.
7. **No main runner** — 14 separate scripts had to be invoked manually.
8. **No error handling** — failed mid-pipeline left partial results with
   no report.

### Refactoring goals

- Single command runs the whole pipeline end-to-end.
- All paths derived from `--raw_data`, `--meta_data`, `--output_dir`.
- All column names, group names, and colors configurable.
- All hardcoded numbers in plots replaced with dynamic values.
- Pipeline report generated for reproducibility.
- Modest enhancements added: Chao1 alpha index, individual sample heatmap,
  manual PCoA fallback when `scikit-bio` is unavailable.

## Output structure

```
<output_dir>/
├── data/                              # Step 1-2 intermediate tables
│   ├── merged_species_clean.xlsx
│   ├── zero_abundance_taxa.xlsx
│   ├── relative_abundance.xlsx
│   ├── relative_abundance_percent.xlsx
│   └── sample_rpkm_summary.xlsx
├── result/
│   ├── step3_alpha/
│   │   ├── alpha_diversity.xlsx
│   │   ├── alpha_statistics.xlsx      # Kruskal-Wallis / ANOVA
│   │   ├── alpha_spearman.xlsx        # if numeric_col provided
│   │   └── figures/
│   │       ├── Observed_boxplot.{pdf,png}
│   │       ├── Shannon_boxplot.{pdf,png}
│   │       ├── Simpson_boxplot.{pdf,png}
│   │       └── Chao1_boxplot.{pdf,png}
│   ├── step4_beta/
│   │   ├── distance/braycurtis_distance.xlsx
│   │   ├── pcoa/pcoa_coordinates.xlsx
│   │   ├── pcoa/pcoa_variance.xlsx
│   │   ├── statistics/PERMANOVA_result.txt
│   │   └── figures/PCoA_BrayCurtis.{pdf,png}
│   └── step5/
│       ├── bacteria_species.xlsx
│       ├── bacteria_relative_abundance.xlsx
│       ├── topN_species.xlsx
│       ├── topN_group_percentage.xlsx
│       └── figures/
│           ├── topN_group_bar.{pdf,png}
│           └── topN_individual_heatmap.{pdf,png}
└── pipeline_report.txt
```

## Algorithm choices

### Step 1: Merge & clean
- Extract clean species name by splitting on `[` (removes `[TAX_xxx]` suffix).
- Sum abundances for duplicate names.
- Drop rows where all sample columns are 0.

### Step 2: Normalize
- Divide each column by its column-sum (sample total read count).
- Equivalent to relative abundance normalization.
- Output both fraction (0–1) and percent (0–100) versions.

### Step 3: Alpha diversity
- **Observed** = count of species with abundance > 0
- **Shannon** = `-Σ p·ln(p)` for p>0
- **Simpson** = `1 - Σ p²` (Gini-Simpson)
- **Chao1** = `Obs + f1·(f1-1) / (2·(f2+1))` with f1=1, f2=2 correction
- Group comparison: Kruskal-Wallis H-test (default) or one-way ANOVA
- Optional: Spearman correlation with a numeric metadata column

### Step 4: Beta diversity
- Distance: Bray-Curtis (default) or any `scipy.spatial.distance.pdist` metric
- Ordination: PCoA via `scikit-bio` (falls back to manual eigendecomposition
  if scikit-bio is unavailable)
- Significance: PERMANOVA (999 permutations) on the Group column
- Visualization: scatter with 95% confidence ellipses (chi² based)

### Step 5: Filter & Top-N
- Filter out viral and fungal taxa using keyword matching on the
  `#Taxonomy` string (case-insensitive). Keywords are configurable.
- Recalculate relative abundance on the filtered table.
- Rank species by mean relative abundance across all samples.
- Take top-N (default 10), collapse the rest into "Others".
- Group-level stacked bar + per-sample heatmap.

## Why a single pipeline runner?

The original 14-script structure forced the user to:
1. Remember which script to run next.
2. Check the output of one step before running the next.
3. Manually update paths/parameters when sharing with collaborators.

The new design:
- A single `run_pipeline.py` entry point.
- Configuration via CLI flags OR a YAML config file.
- Validation at startup.
- A summary report at the end.
