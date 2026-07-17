# Microbiome Diversity Pipeline — Test Examples

## Test 1: Original CSVD data

```bash
python scripts/run_pipeline.py \
  --raw_data "E:/metadata/CSVD_analysis_extracted/CSVD_analysis/data/CSVD_faeces_merged.xlsx" \
  --meta_data "E:/metadata/CSVD_analysis_extracted/CSVD_analysis/data/metadata.xlsx" \
  --output_dir "result" \
  --sample_prefix "P" \
  --numeric_col "CSVD" \
  --group_order "Normal,Low,Medium,Heavy" \
  --top_n 10
```

Expected:
- Input: 5235 species, 31 samples
- Step 1: removes 666 zero-abundance → 2354 species
- Step 3: Kruskal-Wallis p>0.05 (no significant group difference)
- Step 4: PCoA PC1=25.1%, PC2=20.4%, PERMANOVA p≈0.45
- Step 5: removes 734 viral taxa → 1620 bacteria, top10 dominated by E. coli

## Test 2: Synthetic data (different format)

```bash
python scripts/run_pipeline.py \
  --raw_data "synthetic_raw.xlsx" \
  --meta_data "synthetic_meta.xlsx" \
  --output_dir "synthetic_test_output" \
  --species_col "Organism" \
  --sample_prefix "S" \
  --numeric_col "Score" \
  --group_order "Control,Treatment" \
  --top_n 5
```

Expected:
- Input: 14 species, 20 samples
- Step 1: 14 species remain (no zero-abundance)
- Step 5: removes 4 viral/fungal → 10 bacteria
- Top-5 identified: K.pneumoniae, F.prausnitzii, B.longum, B.theta, L.rhamnosus

## Test 3: Config-file mode

```bash
# Use a copy of the template
cp assets/config_template.yaml my_config.yaml
# Edit my_config.yaml
python scripts/run_pipeline.py \
  --config my_config.yaml \
  --raw_data input.xlsx \
  --meta_data meta.xlsx
```

## Verification checklist

After running, verify these files exist:
- `result/data/merged_species_clean.xlsx`
- `result/data/relative_abundance.xlsx`
- `result/result/step3_alpha/alpha_diversity.xlsx`
- `result/result/step3_alpha/figures/{Observed,Shannon,Simpson,Chao1}_boxplot.{pdf,png}`
- `result/result/step4_beta/distance/braycurtis_distance.xlsx`
- `result/result/step4_beta/figures/PCoA_BrayCurtis.{pdf,png}`
- `result/result/step5/bacteria_relative_abundance.xlsx`
- `result/result/step5/figures/topN_group_bar.{pdf,png}`
- `result/pipeline_report.txt`
