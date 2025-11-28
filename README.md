# RNA-seq exploratory analysis (airway dataset)

Small project to practice RNA-seq QC and differential expression in R using the airway dataset.

## Data source

The data come from the airway RNA-seq experiment (Himes et al., 2014):
human airway smooth muscle cells treated with dexamethasone vs. control.
The processed counts are provided by the Bioconductor `airway` package.

## Contents

- `data/`:
  - `airway_scaledcounts.csv` (RNA-seq counts exported from the airway object)
  - `airway_metadata.csv` (sample metadata)

- `scripts/01_qc_airway.R`: QC and exploratory analysis
  - load counts and metadata
  - library size QC plot (`library_sizes_airway.png`)
  - log2 transformation
  - PCA coloured by treatment and cell type (`pca_airway.png`)

- `scripts/02_deseq2_airway.R`: differential expression analysis with DESeq2
  - build `DESeqDataSet` from the `airway` object
  - run DESeq2 for dex (trt vs untrt)
  - results table: `results/deseq2_airway_results.csv`
  - volcano plot: `results/volcano_airway_dex_trt_vs_untrt.png`
  - heatmap of top DE genes: `results/heatmap_top30_airway.png`

## Future work

- Additional QC (vst transformation, sample–sample distances)
- Gene set enrichment / pathway analysis
