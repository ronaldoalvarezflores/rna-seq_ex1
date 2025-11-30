# 01_qc_airway.R
# Basic QC and PCA from airway dataset (RNA-seq)

# 1. Librarys ----
library(tidyverse)

# 2. Read data ----
counts_raw <- read_csv("data/airway_scaledcounts.csv")
metadata   <- read_csv("data/airway_metadata.csv")

# counts_raw: column 'ensgene' = gene ID, rest = samples
counts_mat <- counts_raw %>%
  column_to_rownames("ensgene") %>%
  as.matrix()

# Check that columns and metadata are aligned
stopifnot(all(colnames(counts_mat) == metadata$id))

# 3. QC: Library sizes ----
library_sizes <- colSums(counts_mat)

libsizes_df <- tibble(
  sample   = names(library_sizes),
  lib_size = library_sizes
)

p_libsize <- libsizes_df %>%
  ggplot(aes(x = sample, y = lib_size, fill = sample)) +
  geom_col() +
  theme_minimal() +
  labs(
    title = "Library size per sample (airway)",
    x = "Sample",
    y = "Total readings (counts)"
  ) +
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    axis.title = element_text(size = 13),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y = element_text(size = 11),
    legend.position = "none"   # Remove the legend; you can comment if you want it.
  )

if (!dir.exists("results")) dir.create("results")

ggsave(
  filename = "results/library_sizes_airway.png",
  plot     = p_libsize,
  width    = 6,
  height   = 4,
  dpi      = 300
)

# 4. log2 transformation and PCA ----
log_counts <- log2(counts_mat + 1)

# 4.1. Eliminate genes with no variability (variance = 0)
gene_var <- apply(log_counts, 1, var)
sum(gene_var == 0)  # Just out of curiosity

log_counts_filt <- log_counts[gene_var > 0, ]

# 4.2. PCA (samples in rows, genes in columns)
pca <- prcomp(t(log_counts_filt), scale. = TRUE)

pca_df <- as.data.frame(pca$x)
pca_df$id <- rownames(pca_df)
pca_df <- dplyr::left_join(pca_df, metadata, by = "id")

p_pca <- ggplot(pca_df, aes(x = PC1, y = PC2,
                            color = dex,
                            shape = celltype,
                            label = id)) +
  geom_point(size = 3) +
  geom_text(vjust = -1.1, size = 3) +
  theme_minimal() +
  labs(
    title = "Gene expression PCA (airway)",
    x = "PC1",
    y = "PC2",
    color = "Treatment (dex)",
    shape = "Cell type"
  )

ggsave(
  filename = "results/pca_airway.png",
  plot     = p_pca,
  width    = 6,
  height   = 4,
  dpi      = 300
)

# 5. Save objects for future analysis ----
saveRDS(
  object = list(
    counts_mat     = counts_mat,
    metadata       = metadata,
    log_counts_filt = log_counts_filt,
    pca            = pca
  ),
  file = "results/airway_qc_objects.rds"
)
