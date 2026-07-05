#!/usr/bin/env Rscript

# Estimate PEER factors from an expression matrix or phenotype BED.
# The number of factors commonly follows GTEx practice:
#   15 factors for N < 150
#   30 factors for 150 <= N < 250
#   45 factors for 250 <= N < 350
#   60 factors for N >= 350

suppressPackageStartupMessages({
  library(peer)
  library(argparser)
})

write_table <- function(data, filename, index_name) {
  datafile <- file(filename, open = "wt")
  on.exit(close(datafile))

  header <- c(index_name, colnames(data))
  writeLines(paste0(header, collapse = "\t"), con = datafile, sep = "\n")
  write.table(data, datafile, sep = "\t", col.names = FALSE, quote = FALSE)
}

read_expression <- function(expr_file) {
  if (!file.exists(expr_file)) {
    stop("Expression file not found: ", expr_file, call. = FALSE)
  }

  if (grepl("\\.bed(\\.gz)?$", expr_file)) {
    expr_df <- read.table(
      expr_file,
      sep = "\t",
      header = TRUE,
      check.names = FALSE,
      comment.char = ""
    )
    if (ncol(expr_df) < 5) {
      stop("BED expression file must contain chrom/start/end/id plus sample columns.", call. = FALSE)
    }
    row.names(expr_df) <- expr_df[[4]]
    expr_df <- expr_df[, 5:ncol(expr_df), drop = FALSE]
  } else {
    expr_df <- read.table(
      expr_file,
      sep = "\t",
      header = TRUE,
      check.names = FALSE,
      comment.char = "",
      row.names = 1
    )
  }

  if (nrow(expr_df) == 0 || ncol(expr_df) == 0) {
    stop("Expression data must contain at least one gene and one sample.", call. = FALSE)
  }

  expr_mat <- as.matrix(expr_df)
  storage.mode(expr_mat) <- "numeric"
  expr_mat
}

p <- arg_parser("Run PEER factor estimation")
p <- add_argument(p, "expr.file", help = "Expression matrix or phenotype BED/BED.GZ")
p <- add_argument(p, "prefix", help = "Output file prefix")
p <- add_argument(p, "n", help = "Number of PEER factors to estimate")
p <- add_argument(p, "--covariates", help = "Observed covariates: samples x covariates table with sample IDs in column 1")
p <- add_argument(p, "--alphaprior_a", help = "PEER alpha prior a", default = 0.001)
p <- add_argument(p, "--alphaprior_b", help = "PEER alpha prior b", default = 0.01)
p <- add_argument(p, "--epsprior_a", help = "PEER eps prior a", default = 0.1)
p <- add_argument(p, "--epsprior_b", help = "PEER eps prior b", default = 10)
p <- add_argument(p, "--max_iter", help = "Maximum PEER iterations", default = 1000)
p <- add_argument(p, "--output_dir", short = "-o", help = "Output directory", default = ".")
argv <- parse_args(p)

n_factors <- as.integer(argv$n)
if (!is.finite(n_factors) || n_factors < 1) {
  stop("n must be a positive integer.", call. = FALSE)
}

dir.create(argv$output_dir, recursive = TRUE, showWarnings = FALSE)

cat("PEER: loading expression data ... ")
expr_mat <- read_expression(argv$expr.file)
M <- t(expr_mat)
cat("done.\n")

cat(paste0("PEER: estimating hidden confounders (", n_factors, ")\n"))
model <- PEER()
invisible(PEER_setNk(model, n_factors))
invisible(PEER_setPhenoMean(model, M))
invisible(PEER_setPriorAlpha(model, as.numeric(argv$alphaprior_a), as.numeric(argv$alphaprior_b)))
invisible(PEER_setPriorEps(model, as.numeric(argv$epsprior_a), as.numeric(argv$epsprior_b)))
invisible(PEER_setNmax_iterations(model, as.integer(argv$max_iter)))

has_covariates <- !is.null(argv$covariates) && !is.na(argv$covariates)
if (has_covariates) {
  if (!file.exists(argv$covariates)) {
    stop("Covariate file not found: ", argv$covariates, call. = FALSE)
  }

  covar_df <- read.table(argv$covariates, sep = "\t", header = TRUE, row.names = 1, as.is = TRUE)
  missing_samples <- setdiff(rownames(M), rownames(covar_df))
  if (length(missing_samples) > 0) {
    stop("Covariate file is missing samples: ", paste(head(missing_samples, 10), collapse = ", "), call. = FALSE)
  }
  covar_df <- covar_df[rownames(M), , drop = FALSE]
  covar_df[] <- lapply(covar_df, as.numeric)

  cat(paste0("  * including ", ncol(covar_df), " observed covariates\n"))
  invisible(PEER_setCovariates(model, as.matrix(covar_df)))
}

invisible(system.time(PEER_update(model)))

X <- PEER_getX(model)
A <- PEER_getAlpha(model)
R <- t(PEER_getResiduals(model))

if (has_covariates) {
  inferred_count <- ncol(X) - ncol(covar_df)
  inferred_names <- if (inferred_count > 0) paste0("InferredCov", seq_len(inferred_count)) else character(0)
  cols <- c(colnames(covar_df), inferred_names)
} else {
  cols <- paste0("InferredCov", seq_len(ncol(X)))
}

rownames(X) <- rownames(M)
colnames(X) <- cols
rownames(A) <- cols
colnames(A) <- "Alpha"
A <- as.data.frame(A)
A$Relevance <- 1.0 / A$Alpha
rownames(R) <- colnames(M)
colnames(R) <- rownames(M)

cat("PEER: writing results ... ")
write_table(
  t(X),
  file.path(argv$output_dir, paste0(argv$prefix, ".PEER_covariates.txt")),
  "ID"
)
write_table(
  A,
  file.path(argv$output_dir, paste0(argv$prefix, ".PEER_alpha.txt")),
  "ID"
)
write_table(
  R,
  file.path(argv$output_dir, paste0(argv$prefix, ".PEER_residuals.txt")),
  "ID"
)
cat("done.\n")
