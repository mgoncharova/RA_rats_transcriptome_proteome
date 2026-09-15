library(limma)

# LFQ matrix after KNN imputation
proteome <- read.csv(
  "data/processed/02_proteome_lfq_log2_knn_imputed.csv",
  row.names = 1,
  check.names = FALSE
)

#metadata
metadata <- read.delim(
  "data/metadata/proteome_metadata.tsv",
  stringsAsFactors = FALSE
)

# Put metadata in the same order as proteome columns
metadata <- metadata[match(colnames(proteome), metadata$sample), ]
stopifnot(identical(colnames(proteome), metadata$sample))

#define experimental groups
group <- factor(metadata$group, levels = c("Normal", "SBA"))

# one mean abundance value per group
design <- model.matrix(~ 0 + group)
colnames(design) <- c("Normal", "SBA")

#fit linear model for every protein
fit <- lmFit(proteome, design)

# Test SBA minus Normal
contrast <- makeContrasts(SBA - Normal, levels = design)
fit <- contrasts.fit(fit, contrast)

#Bayes moderation
fit <- eBayes(fit)

#save results for all proteins
results <- topTable(
  fit,
  number = Inf,
  adjust.method = "BH",
  sort.by = "P"
)

write.csv(
  results,
  "results/tables/05_proteome_limma_SBA_vs_Normal.csv"
)