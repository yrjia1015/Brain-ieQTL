#!/usr/bin/env Rscript

usage <- function() {
  cat(
    "Usage:\n",
    "  Rscript inverse_normal_transform.R <expression_matrix> <output_matrix>\n\n",
    "Apply rank-based inverse normal transformation to each gene/row.\n\n",
    "Input format:\n",
    "  Tab-delimited matrix with gene ID in column 1 and sample IDs in the header.\n\n",
    "Output format:\n",
    "  Tab-delimited matrix with ID in column 1 and transformed sample values.\n\n",
    "Example:\n",
    "  Rscript inverse_normal_transform.R genes.TMM.txt genes.INT.txt\n",
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
  stop("Input expression matrix not found: ", input_file, call. = FALSE)
}

rank_inverse_normal <- function(x) {
  qnorm((rank(x, na.last = "keep", ties.method = "random") - 3 / 8) / sum(!is.na(x)))
}

rank_inverse_normal_matrix <- function(x) {
  t(apply(x, 1, rank_inverse_normal))
}

expr_df <- read.table(
  file = input_file,
  sep = "\t",
  header = TRUE,
  check.names = FALSE,
  comment.char = "",
  row.names = 1
)

if (nrow(expr_df) == 0 || ncol(expr_df) == 0) {
  stop("Input matrix must contain at least one gene and one sample.", call. = FALSE)
}

expr_mat <- as.matrix(expr_df)
storage.mode(expr_mat) <- "numeric"

set.seed(20251028)
int_mat <- rank_inverse_normal_matrix(expr_mat)
rownames(int_mat) <- rownames(expr_mat)
colnames(int_mat) <- colnames(expr_mat)

output_df <- data.frame(ID = rownames(int_mat), int_mat, check.names = FALSE)
write.table(
  output_df,
  file = output_file,
  col.names = TRUE,
  row.names = FALSE,
  quote = FALSE,
  sep = "\t"
)
