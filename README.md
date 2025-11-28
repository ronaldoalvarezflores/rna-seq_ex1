# RNA-seq exploratory analysis (airway dataset)

Small project to practice RNA-seq QC and differential expression in R using the airway dataset.

## Data source

The data come from the airway RNA-seq experiment (Himes et al., 2014): human airway smooth muscle cells treated with dexamethasone vs. control. The processed counts are provided by the Bioconductor `airway` package (original GEO accession: GSE52778).

## Contents

-   `data/`:
    -   `airway_scaledcounts.csv` (RNA-seq counts exported from the airway object)
    -   `airway_metadata.csv` (sample metadata)
-   `scripts/01_qc_airway.R`: QC and exploratory analysis
    -   load counts and metadata
    -   library size QC plot (`results/library_sizes_airway.png`)
    -   log2 transformation
    -   PCA coloured by treatment and cell type (`results/pca_airway.png`)
-   `scripts/02_deseq2_airway.R`: differential expression analysis with DESeq2
    -   build `DESeqDataSet` from the `airway` object
    -   run DESeq2 for dex (trt vs untrt)
    -   raw results table: `results/deseq2_airway_results.csv`
    -   **annotated results table:** `results/deseq2_airway_results_annotated.csv`
    -   **top 10 DE genes:** `results/deseq2_airway_top10_genes.csv`
    -   volcano plot: `results/volcano_airway_dex_trt_vs_untrt.png`
    -   heatmap of top DE genes: `results/heatmap_top30_airway.png`

## Summary of results

### Quality control

-   Library sizes are between \~18 and \~32 million reads per sample, with no obvious outliers.
-   The PCA on log2-transformed counts shows that dex-treated samples separate clearly from controls along the first principal component, while samples from the same cell line cluster together. This suggests a strong and consistent transcriptional response to dexamethasone on top of donor-specific baseline differences.

### Differential expression (DESeq2)

Using DESeq2 with the design `~ dex` (dex-treated vs untreated controls):

-   **2,694 genes** show an adjusted p-value \< 0.05.

-   With an additional effect-size threshold (\|log2FC\| \> 1), **475 genes** are up-regulated and **397 genes** are down-regulated in dex-treated samples.

-   The top 10 differentially expressed genes (sorted by adjusted p-value) are:

    1.  **SPARCL1 (SPARC like 1)** – log2FC ≈ **4.60**\
        → \~**24-fold higher** expression in dex-treated samples; padj ≈ **1.7 × 10⁻¹⁰⁰**.

    2.  **STOM (stomatin)** – log2FC ≈ **1.45**\
        → \~**2.7-fold higher** expression; padj ≈ **1.4 × 10⁻⁶¹**.

    3.  **PER1 (period circadian regulator 1)** – log2FC ≈ **3.18**\
        → \~**9.1-fold higher** expression; padj ≈ **1.9 × 10⁻⁵²**.

    4.  **PHC2 (polyhomeotic homolog 2)** – log2FC ≈ **1.39**\
        → \~**2.6-fold higher** expression; padj ≈ **4.4 × 10⁻⁴⁸**.

    5.  **MT2A (metallothionein 2A)** – log2FC ≈ **2.20**\
        → \~**4.6-fold higher** expression; padj ≈ **6.0 × 10⁻⁴⁷**.

    6.  **DUSP1 (dual specificity phosphatase 1)** – log2FC ≈ **2.95**\
        → \~**7.7-fold higher** expression; padj ≈ **5.4 × 10⁻⁴⁵**.

    7.  **MAOA (monoamine oxidase A)** – log2FC ≈ **3.31**\
        → \~**9.9-fold higher** expression; padj ≈ **3.4 × 10⁻⁴⁴**.

    8.  **ZBTB16 (zinc finger and BTB domain containing 16)** – log2FC ≈ **7.18**\
        → \~**145-fold higher** expression (very strong induction); padj ≈ **5.2 × 10⁻⁴⁴**.

    9.  **KCTD12 (potassium channel tetramerization domain containing 12)** – log2FC ≈ **–2.50**\
        → \~**5–6-fold lower** expression in dex-treated samples; padj ≈ **4.6 × 10⁻⁴¹**.

    10. **SAMHD1 (SAM and HD domain containing deoxynucleoside triphosphate triphosphohydrolase 1)** – log2FC ≈ **3.86**\
        → \~**14.5-fold higher** expression; padj ≈ **4.0 × 10⁻⁴⁰**.

Overall, the top DE genes show very large fold changes and extremely small adjusted p-values, consistent with a strong transcriptional response to dexamethasone in airway smooth muscle cells. The volcano plot highlights a large number of significantly regulated genes with moderate-to-large effect sizes, and the heatmap of the top 30 DE genes shows a clear separation between dex-treated and control samples, with good clustering of biological replicates.

## How to run

-   Open the project in RStudio (`rna-seq_ex1.Rproj`).
-   Make sure the following packages are installed: `tidyverse`, `DESeq2`, `airway`, `pheatmap`, `AnnotationDbi`, `org.Hs.eg.db`.
-   Run:
    -   `scripts/01_qc_airway.R` for QC and PCA.
    -   `scripts/02_deseq2_airway.R` for differential expression, annotated results and plots.

## Future work

-   Additional QC (vst transformation, sample–sample distances)
-   Gene set enrichment / pathway analysis
