# 02_deseq2_airway.R
# Differential expression analysis for the airway dataset (dex vs untrt)

# 1. Libraries ----
library(DESeq2)
library(airway)
library(tidyverse)
library(pheatmap)
library(AnnotationDbi)
library(org.Hs.eg.db)

# 2. Load airway data ----
data("airway")
se <- airway

# Relevel dex so that "untrt" is the reference level
colData(se)$dex <- relevel(colData(se)$dex, ref = "untrt")

# 3. Create DESeqDataSet ----
dds <- DESeqDataSet(se, design = ~ dex)

# Filter out very lowly expressed genes
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep, ]

# 4. Run DESeq2 ----
dds <- DESeq(dds)

res <- results(dds, contrast = c("dex", "trt", "untrt"))
res_ordered <- res[order(res$padj), ]

res_df <- as.data.frame(res_ordered) %>%
  rownames_to_column(var = "gene_id")

# Create results folder if needed
if (!dir.exists("results")) dir.create("results", showWarnings = FALSE)

# Save raw DESeq2 results
write.csv(
  res_df,
  file = "results/deseq2_airway_results.csv",
  row.names = FALSE
)

# 5. Annotate genes with SYMBOL and GENENAME ----
annot <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys    = res_df$gene_id,
  keytype = "ENSEMBL",
  columns = c("SYMBOL", "GENENAME")
)

res_annot <- res_df %>%
  dplyr::left_join(annot, by = c("gene_id" = "ENSEMBL"))

# Save annotated results
write.csv(
  res_annot,
  file = "results/deseq2_airway_results_annotated.csv",
  row.names = FALSE
)

# 6. Volcano plot ----

volcano_df <- res_df %>%
  mutate(
    neg_log10_padj = -log10(padj),
    sig = case_when(
      padj < 0.05 & log2FoldChange >  1 ~ "Up (padj<0.05, LFC>1)",
      padj < 0.05 & log2FoldChange < -1 ~ "Down (padj<0.05, LFC<-1)",
      TRUE ~ "NS"
    )
  ) %>%
  filter(!is.na(padj))

p_volcano <- ggplot(volcano_df,
                    aes(x = log2FoldChange, y = neg_log10_padj, color = sig)) +
  geom_point(alpha = 0.7, size = 1.5) +
  theme_minimal() +
  labs(
    title = "DESeq2 airway: dex trt vs untrt",
    x = "log2 fold change",
    y = "-log10(padj)",
    color = "Status"
  )

ggsave(
  filename = "results/volcano_airway_dex_trt_vs_untrt.png",
  plot     = p_volcano,
  width    = 6,
  height   = 4,
  dpi      = 300
)

# 7. Heatmap of top 30 DE genes ----

top_genes <- head(res_df$gene_id[order(res_df$padj)], 30)
mat_top <- counts(dds, normalized = TRUE)[top_genes, ]

# log2-transform and scale by gene (row)
mat_top_scaled <- t(scale(t(log2(mat_top + 1))))

sample_anno <- as.data.frame(colData(dds)[, c("dex", "cell")])
sample_anno$dex  <- as.factor(sample_anno$dex)
sample_anno$cell <- as.factor(sample_anno$cell)

pheatmap(
  mat_top_scaled,
  annotation_col = sample_anno,
  show_rownames = TRUE,
  show_colnames = TRUE,
  filename = "results/heatmap_top30_airway.png",
  width = 6,
  height = 6
)
# 8. Save top 10 DE genes ----
top10 <- res_annot[order(res_annot$padj),
                   c("gene_id", "SYMBOL", "GENENAME", "log2FoldChange", "padj")][1:10, ]

write.csv(
  top10,
  file = "results/deseq2_airway_top10_genes.csv",
  row.names = FALSE
)

print(top10)  # Optional: To also see them in the console when running the script
