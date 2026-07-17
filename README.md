# mbdiv — Microbiome Diversity Analysis Pipeline

[![Python 3.9+](https://img.shields.io/badge/python-3.9+-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![pyflakes clean](https://img.shields.io/badge/pyflakes-clean-green.svg)](https://github.com/PyCQA/pyflakes)

One-command pipeline for species-level metagenomic diversity analysis.

```
raw_data  ──►  merge species  ──►  normalize  ──►  alpha diversity (4 indices + boxplots)
                                      │
              beta diversity  ◄───────┘  (Bray-Curtis + PCoA + PERMANOVA)
                                      │
              filter viruses/fungi  ──►  Top-N stacked bar + heatmap
```

## Quick Start

```bash
# Install
pip install mbdiv-1.0.0-py3-none-any.whl[full]

# Run — that's it
mbdiv species.xlsx metadata.xlsx
```

Output appears in `./result/` with all tables, statistics, and figures.

## Installation

### Option 1: From wheel (recommended)

```bash
# Core features only (alpha, beta, top-N without scikit-bio)
pip install mbdiv-1.0.0-py3-none-any.whl

# Full features (PCoA via scikit-bio + YAML config support)
pip install mbdiv-1.0.0-py3-none-any.whl[full]
```

### Option 2: From source

```bash
git clone https://github.com/BolinWo/mbdiv.git
cd mbdiv
pip install .[full]
```

### Option 3: Portable zip (offline machines)

```bash
# Unzip mbdiv-1.0.0-portable.zip
# Windows: double-click setup.bat
# Linux/Mac: bash setup.sh
```

## Usage

### Basic

```bash
# Minimal — auto-detects everything (column names, sample prefix, groups)
mbdiv species.xlsx metadata.xlsx

# Specify output directory and top-N
mbdiv species.xlsx metadata.xlsx -o my_results -t 15
```

### All Parameters (single-letter flags)

| Flag | Long form | Default | Description |
|------|-----------|---------|-------------|
| *(positional)* | — | — | Species abundance Excel file |
| *(positional)* | — | — | Sample metadata Excel file |
| `-r` | `--raw` | — | Species abundance Excel (alternative to positional) |
| `-m` | `--meta` | — | Sample metadata Excel (alternative to positional) |
| `-o` | `--output` | `result` | Output directory |
| `-t` | `--top` | `10` | Top-N species count for stacked bar |
| `-c` | `--config` | — | YAML config file for advanced options |
| `-h` | `--help` | — | Show help message |

### Two invocation styles

```bash
# Style 1: positional (concise)
mbdiv species.xlsx metadata.xlsx -o results -t 10

# Style 2: explicit flags (readable)
mbdiv -r species.xlsx -m metadata.xlsx -o results -t 10
```

### Advanced: YAML config

For reproducible analyses, save your settings as `config.yaml`:

```bash
mbdiv species.xlsx metadata.xlsx -c config.yaml
```

See `assets/config_template.yaml` for all configurable options including:
- Column name overrides (species, taxonomy, sample prefix)
- Group order and colors
- Alpha diversity metrics and statistical test
- Beta distance metric and PERMANOVA permutations
- Virus/fungi filter keywords
- Figure format (PDF/PNG/both) and DPI

## Input Format

### raw_data.xlsx — Species abundance table

| #Species | P001 | P002 | ... | #Taxonomy |
|----------|------|------|-----|-----------|
| E.coli [TAX_001] | 152.3 | 89.1 | ... | d__Bacteria;p__Proteobacteria;... |
| B.fragilis [TAX_002] | 203.5 | 0.0 | ... | d__Bacteria;p__Bacteroidota;... |

- **First column**: species name (auto-detected: `#Species`, `Species`, or any column starting with `#`)
- **Middle columns**: sample abundance values (numeric, auto-detected)
- **Last column**: taxonomy string (auto-detected: any column containing "tax")

### meta_data.xlsx — Sample metadata

| Sample | Group | CSVD |
|--------|-------|------|
| P001 | Normal | 0 |
| P002 | Heavy | 3 |

- **Sample** column: sample IDs matching raw_data column headers
- **Group** column: group assignment for comparison
- **Numeric column** (optional): for Spearman correlation with alpha diversity

**All column names are auto-detected.** No need to specify them unless your format is unusual.

## Output Structure

```
<output_dir>/
├── pipeline_report.txt              # Full audit report
├── data/                            # Intermediate data tables
│   ├── merged_species_clean.xlsx    #   Step 1: merged + zero-removed
│   ├── zero_abundance_taxa.xlsx     #   Step 1: removed taxa log
│   ├── relative_abundance.xlsx      #   Step 2: normalized (fraction)
│   ├── relative_abundance_percent.xlsx  #  Step 2: normalized (%)
│   └── sample_rpkm_summary.xlsx     #   Step 2: sequencing depth
└── result/
    ├── step3_alpha/                 # Alpha diversity
    │   ├── alpha_diversity.xlsx     #   4 indices × N samples
    │   ├── alpha_statistics.xlsx    #   Kruskal-Wallis p-values
    │   ├── alpha_spearman.xlsx      #   Spearman correlation (if numeric col)
    │   └── figures/
    │       ├── Observed_boxplot.{pdf,png}
    │       ├── Shannon_boxplot.{pdf,png}
    │       ├── Simpson_boxplot.{pdf,png}
    │       └── Chao1_boxplot.{pdf,png}
    ├── step4_beta/                  # Beta diversity
    │   ├── distance/
    │   │   └── braycurtis_distance.xlsx
    │   ├── pcoa/
    │   │   ├── pcoa_coordinates.xlsx
    │   │   └── pcoa_variance.xlsx
    │   ├── statistics/
    │   │   └── PERMANOVA_result.txt
    │   └── figures/
    │       └── PCoA_BrayCurtis.{pdf,png}
    └── step5/                       # Top-N composition
        ├── bacteria_species.xlsx
        ├── bacteria_relative_abundance.xlsx
        ├── top10_species.xlsx
        ├── top10_group_percentage.xlsx
        └── figures/
            ├── top10_group_bar.{pdf,png}
            └── top10_individual_heatmap.{pdf,png}
```

## Pipeline Steps

### Step 1: Merge & Clean
- Merge same-name species (sum abundance)
- Remove all-zero abundance taxa
- Extract clean species names

### Step 2: Normalize
- Per-sample relative abundance (column sum → 1.0)
- Save both fractional and percentage versions
- Sequencing depth summary

### Step 3: Alpha Diversity
- **Observed** — species richness
- **Shannon** — H' index (natural log)
- **Simpson** — Gini-Simpson (1-D)
- **Chao1** — richness estimator
- Kruskal-Wallis test across groups
- Optional Spearman correlation with numeric metadata
- Boxplot + jitter visualization (PDF + PNG)

### Step 4: Beta Diversity
- Bray-Curtis distance matrix
- PCoA ordination (scikit-bio or SVD fallback)
- PERMANOVA test (999 permutations)
- PCoA scatter plot with 95% confidence ellipses
- Dynamic axis labels (variance %) and PERMANOVA p-value annotation

### Step 5: Top-N Composition
- Filter out viruses and fungi by taxonomy keywords
- Recalculate relative abundance on bacteria-only data
- Select Top-N species by mean abundance
- Group-level stacked bar chart
- Individual sample heatmap

## Dependencies

| Package | Required | Purpose |
|---------|----------|---------|
| pandas | Yes | Data manipulation |
| numpy | Yes | Numerical computation |
| scipy | Yes | Statistics + distance |
| matplotlib | Yes | Plotting |
| seaborn | Yes | Statistical visualization |
| openpyxl | Yes | Excel I/O |
| scikit-bio | Optional | PCoA + PERMANOVA (falls back to manual SVD) |
| pyyaml | Optional | YAML config file support |

## Reproducibility

The pipeline generates a `pipeline_report.txt` in the output directory with:
- Input file paths
- All output file paths
- Key statistical results (PERMANOVA p-value)
- Total elapsed time

For exact reproducibility, save your YAML config and use the same input files.

## License

MIT — see [LICENSE](LICENSE)

## Citation

If you use mbdiv in your research, please cite:

```
mbdiv: A single-command microbiome diversity analysis pipeline
https://github.com/BolinWo/mbdiv
```
