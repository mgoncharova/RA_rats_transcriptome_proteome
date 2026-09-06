library(impute)

# read filtered matrix from file
proteome_filtered <- read.csv(
  "data/processed/01_proteome_lfq_log2_filtered.csv",
  row.names = 1,
  check.names = FALSE
)

message("Proteins before KNN: ", nrow(proteome_filtered))
message("Missing values before KNN: ", sum(is.na(proteome_filtered)))

# K-nearest-neighbour imputation
proteome_knn <- impute.knn(
  as.matrix(proteome_filtered),
  k = 3
)$data

# save new imputed matrix
knn_file <- "data/processed/02_proteome_lfq_log2_knn_imputed.csv"

write.csv(
  as.data.frame(proteome_knn),
  file = knn_file,
  row.names = TRUE
)

message("Missing values after KNN: ", sum(is.na(proteome_knn)))
message("Saved KNN-imputed matrix: ", normalizePath(knn_file))