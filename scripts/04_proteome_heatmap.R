library(pheatmap)

#Read protein abundance matrix and metadata
proteome <- read.csv(
  "data/processed/02_proteome_lfq_log2_knn_imputed.csv",
  row.names = 1,
  check.names = FALSE
)
metadata <- read.delim(
  "data/metadata/proteome_metadata.tsv",
  stringsAsFactors = FALSE
)

#Put metadata in the same order as matrix columns
metadata <- metadata[match(colnames(proteome), metadata$sample), ]

stopifnot(identical(colnames(proteome), metadata$sample))

#row names
sample_annotation <- data.frame(
  group = metadata$group
)
rownames(sample_annotation) <- metadata$sample

#Pearson-correlation distances
sample_distance <- as.dist(1 - cor(proteome))
protein_distance <- as.dist(1 - cor(t(proteome)))

#create and save heatmap
jpeg(
  "results/figures/04_proteome_heatmap.jpeg",
  width = 8,
  height = 10,
  units = "in",
  res = 600
)

pheatmap(
  proteome,
  scale = "row",
  show_rownames = FALSE,
  annotation_col = sample_annotation,
  clustering_distance_rows = protein_distance,
  clustering_distance_cols = sample_distance,
  clustering_method = "ward.D2",
  main = "Proteome heatmap"
)

dev.off()