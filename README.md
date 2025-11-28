# RNA-seq exploratory analysis (airway dataset)

Small project to practice RNA-seq QC and PCA in R using the airway dataset.

## Contents

-   `data/`: `airway_scaledcounts.csv` (RNA-seq counts) and `airway_metadata.csv` (sample metadata)
-   `scripts/01_qc_airway.R`: R script for:
    -   loading counts and metadata
    -   library size QC
    -   log2 transformation
    -   PCA coloured by treatment and cell type
-   `results/`:
    -   `library_sizes_airway.png`
    -   `pca_airway.png`
    -   `airway_qc_objects.rds`

## Future work

-   Differential expression analysis (DESeq2)
-   Gene set enrichment / pathway analysis
-   More QC plots (e.g. sample–sample distances)
