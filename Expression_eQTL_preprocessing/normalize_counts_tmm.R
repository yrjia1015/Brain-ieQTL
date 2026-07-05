#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(edgeR)
})

usage <- function() {
  cat(
    "Usage:\n",
    "  Rscript normalize_counts_tmm.R <count_matrix> <output_matrix>\n\n",
    "Run edgeR TMM normalization and output log2 CPM values.\n\n",
    "Input format:\n",
    "  Tab-delimited matrix with gene ID in column 1 and sample IDs in the header.\n\n",
    "Output format:\n",
    "  Tab-delimited matrix with ID in column 1 and sample IDs in remaining columns.\n\n",
    "Example:\n",
    "  Rscript normalize_counts_tmm.R genes.filtered.count.txt genes.TMM.txt\n",
    sep = ""
  )
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  usage()
  quit(status = 1)
}

input_file <- args[[1]]
output_file <- args[[2]]

if (!file.exists(input_file)) {
  stop("Input count matrix not found: ", input_file, call. = FALSE)
}

counts_df <- read.table(
  file = input_file,
  sep = "\t",
  header = TRUE,
  check.names = FALSE,
  comment.char = "",
  row.names = 1
)

if (nrow(counts_df) == 0 || ncol(counts_df) == 0) {
  stop("Input matrix must contain at least one gene and one sample.", call. = FALSE)
}

counts_mat <- as.matrix(counts_df)
storage.mode(counts_mat) <- "numeric"

dge <- DGEList(counts = counts_mat)
dge <- calcNormFactors(dge)
log_cpm <- cpm(dge, log = TRUE)

output_df <- data.frame(ID = rownames(log_cpm), log_cpm, check.names = FALSE)
write.table(
  output_df,
  file = output_file,
  row.names = FALSE,
  col.names = TRUE,
  quote = FALSE,
  sep = "\t"
)
