library(impute)

# reading of preprocessed log2 LFQ matrix
proteome <- read.csv(
  "data/processed/00_proteome_lfq_log2_preprocessed.csv",
  row.names = 1,
  check.names = FALSE
)

# define sample groups
sba <- c("P17410_B_10", "P17410_B_8", "P17410_B_9")

normal <- c("P17410_B_25","P17410_B_28","P17410_B_29")

# keep proteins detected in at least 2 of 3 samples
# in either the SBA or normal group (explanation in mark down file))
proteome_filtered <- proteome[
  rowSums(!is.na(proteome[, sba])) >= 2 |
  rowSums(!is.na(proteome[, normal])) >= 2,
]

# save filtered matrix
filtered_file <- "data/processed/01_proteome_lfq_log2_filtered.csv"

write.csv(
  proteome_filtered,
  file = filtered_file,
  row.names = TRUE
)

message("Proteins before filtering: ", nrow(proteome))
message("Proteins after filtering: ", nrow(proteome_filtered))
message("Saved filtered matrix: ", normalizePath(filtered_file))