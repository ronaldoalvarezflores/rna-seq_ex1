# 01_qc_airway.R
# QC básico y PCA del dataset airway (RNA-seq)

# 1. Librerías ----
library(tidyverse)

# 2. Leer datos ----
counts_raw <- read_csv("data/airway_scaledcounts.csv")
metadata   <- read_csv("data/airway_metadata.csv")

# counts_raw: columna 'ensgene' = ID de gen, resto = muestras
counts_mat <- counts_raw %>%
  column_to_rownames("ensgene") %>%
  as.matrix()

# Comprobar que columnas y metadata están alineadas
stopifnot(all(colnames(counts_mat) == metadata$id))

# 3. QC: tamaños de librería ----
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
    title = "Tamaño de librería por muestra (airway)",
    x = "Muestra",
    y = "Total de lecturas (conteos)"
  )

if (!dir.exists("results")) dir.create("results")

ggsave(
  filename = "results/library_sizes_airway.png",
  plot     = p_libsize,
  width    = 6,
  height   = 4,
  dpi      = 300
)

# 4. Transformación log2 y PCA ----
log_counts <- log2(counts_mat + 1)

# 4.1. Eliminar genes sin variabilidad (varianza = 0)
gene_var <- apply(log_counts, 1, var)
sum(gene_var == 0)  # solo por curiosidad

log_counts_filt <- log_counts[gene_var > 0, ]

# 4.2. PCA (muestras en filas, genes en columnas)
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
    title = "PCA de expresión génica (airway)",
    x = "PC1",
    y = "PC2",
    color = "Tratamiento (dex)",
    shape = "Tipo celular"
  )

ggsave(
  filename = "results/pca_airway.png",
  plot     = p_pca,
  width    = 6,
  height   = 4,
  dpi      = 300
)

# 5. Guardar objetos para futuros análisis ----
saveRDS(
  object = list(
    counts_mat     = counts_mat,
    metadata       = metadata,
    log_counts_filt = log_counts_filt,
    pca            = pca
  ),
  file = "results/airway_qc_objects.rds"
)
