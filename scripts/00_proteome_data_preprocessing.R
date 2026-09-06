#1 libraries activation first

library(limma) #different abundance analysis using linear model
library(impute)  #missing value imputation (KNN included)
library(tidyverse)   #data import, cleaning, transformation, building plots


#2 file reading 

protein_groups <- read.delim(
  "data/raw/proteinGroups.txt",
  check.names = FALSE #dont change column's names
)

dim(protein_groups)  #table dimensions
colnames(protein_groups)[1:20]




#3 find LFQ intensity columns
lfq_cols <- grep("^LFQ intensity", colnames(protein_groups), value = TRUE) #find
#columns that start with LFQ intensity and show their name

#4 create sample annotation
sample_annotation <- data.frame(
  label = sub("^LFQ intensity ", "", lfq_cols),
  condition = ifelse(
    sub("^LFQ intensity ", "", lfq_cols) %in% c("P17410_B_10", "P17410_B_8", "P17410_B_9"),
    "SBA",
    "normal"
  )
)




#5 remove low-confidence and technical protein groups
protein_groups_clean <- protein_groups %>%
  filter(
    is.na(`Potential contaminant`) | `Potential contaminant` != "+",#keep lines
    #if this line does'nt have +
    is.na(Reverse) | Reverse != "+",
    is.na(`Only identified by site`) | `Only identified by site` != "+"
  ) # 

dim(protein_groups)
dim(protein_groups_clean)



#6 create LFQ abundance matrix
lfq_matrix <- protein_groups_clean %>%
  select(all_of(lfq_cols)) %>%
  as.matrix()

colnames(lfq_matrix) <- sample_annotation$label
rownames(lfq_matrix) <- protein_groups_clean$`Gene names`

dim(lfq_matrix)
lfq_matrix[1:5, ]



#7 prepare matrix for analysis
lfq_matrix[lfq_matrix == 0] <- NA #we should'n have 0 in lines, only NA

sum(is.na(lfq_matrix))
colSums(is.na(lfq_matrix))

lfq_log2 <- log2(lfq_matrix)#we convert LFQ intensities of proteins to log2 
#instad of 0.0000 we get log2

dim(lfq_log2)
lfq_log2[1:5, ]

# 8 save cleaned log2-transformed LFQ matrix
output_file <- "data/processed/00_proteome_lfq_log2_preprocessed.csv"

write.csv(
  as.data.frame(lfq_log2),
  file = output_file,
  row.names = TRUE
)



