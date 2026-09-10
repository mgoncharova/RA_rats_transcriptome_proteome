library(tidyverse) #load required packages
library(ggplot2) 


#read knn-imputed matrix
proteome_knn_imputed <- read.csv(
  "data/processed/02_proteome_lfq_log2_knn_imputed.csv",
  row.names = 1,
  check.names = FALSE  
)

#read sample's metadata
metadata <- read.delim(
  "data/metadata/proteome_metadata.tsv",
  stringsAsFactors = FALSE
)

 #put metadata rows in the same order as matrix columns
metadata <- metadata[match(colnames(proteome_knn_imputed), metadata$sample), ]

#preparing for pca
pca <- prcomp(t(proteome_knn_imputed))
pca_data <- data.frame(
  sample = metadata$sample,
  group = metadata$group,
  PC1 = pca$x[, 1],
  PC2 = pca$x[, 2]
)

#plot
plot_pca <- ggplot(pca_data, aes(x = PC1, y = PC2, color = group)) +
  geom_point(size = 4) +
  geom_text(aes(label = sample), vjust = -1) +
  theme_minimal()

ggsave(
  "results/figures/03_proteome_PCA.jpeg",
  plot_pca,
  width = 8,
  height = 6,
  dpi = 600
)