# 02_deseq2_airway.R
# Expresión diferencial dex (tratado) vs untrt (control)
# usando el dataset airway original (Bioconductor)

# 1. Librerías ----
library(DESeq2)
library(airway)
library(tidyverse)

# 2. Cargar datos airway ----
data("airway")
se <- airway

# Echamos un vistazo rápido a la metadata
colData(se)

# Asegurar niveles de la variable de interés: dex (untrt = control)
colData(se)$dex <- relevel(colData(se)$dex, ref = "untrt")

# 3. Crear objeto DESeqDataSet ----
dds <- DESeqDataSet(se, design = ~ dex)

# Filtrar genes muy poco expresados (ruido)
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep, ]

# 4. Ejecutar DESeq2 ----
dds <- DESeq(dds)

# Resultado principal: tratado (trt) vs control (untrt)
res <- results(dds, contrast = c("dex", "trt", "untrt"))

# Ordenar por p-ajustada
res_ordered <- res[order(res$padj), ]

# Convertimos a data.frame y añadimos columna de gen
res_df <- as.data.frame(res_ordered) %>%
  rownames_to_column(var = "gene_id")

# 5. Guardar tabla de resultados ----
if (!dir.exists("results")) dir.create("results", showWarnings = FALSE)

write.csv(
  res_df,
  file = "results/deseq2_airway_results.csv",
  row.names = FALSE
)

# 6. Volcano plot ----

# Añadimos columna para significancia
volcano_df <- res_df %>%
  mutate(
    neg_log10_padj = -log10(padj),
    sig = case_when(
      padj < 0.05 & log2FoldChange >  1 ~ "Up (padj<0.05, LFC>1)",
      padj < 0.05 & log2FoldChange < -1 ~ "Down (padj<0.05, LFC<-1)",
      TRUE ~ "NS"
    )
  )

# Pequeño filtro para no hacer infinitos con NA
volcano_df <- volcano_df %>% filter(!is.na(padj))

p_volcano <- ggplot(volcano_df,
                    aes(x = log2FoldChange, y = neg_log10_padj, color = sig)) +
  geom_point(alpha = 0.7, size = 1.5) +
  theme_minimal() +
  labs(
    title = "DESeq2 airway: dex trt vs untrt",
    x = "log2 fold change",
    y = "-log10(padj)",
    color = "Estado"
  )

ggsave(
  filename = "results/volcano_airway_dex_trt_vs_untrt.png",
  plot     = p_volcano,
  width    = 6,
  height   = 4,
  dpi      = 300
)

# 7. Heatmap simple de genes top ----

# Vamos a coger los 30 genes con menor padj
top_genes <- head(res_df$gene_id[order(res_df$padj)], 30)

mat_top <- counts(dds, normalized = TRUE)[top_genes, ]

# Escalamos por gen (fila)
mat_top_scaled <- t(scale(t(log2(mat_top + 1))))

# Para el heatmap usaremos pheatmap
if (!requireNamespace("pheatmap", quietly = TRUE)) {
  install.packages("pheatmap")
}
library(pheatmap)

# Construimos anotación de columnas con la condición
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
