#!/usr/bin/env Rscript

# ============================================================
# Run mashr on protein-coding "strong set" using fitted prior
# learned from a protein-coding "random set".
#
# Command-line arguments (all required unless stated):
#   --random_set   <path>   (TSV/TXT readable by data.table::fread)
#   --strong_set   <path>
#   --out_rds                     <path>   (output .rds)
# Optional:
#   --lfsr_threshold              <float>  (default: 0.05)
#   --pair_id_col                 <string> (default: Pair)
#   --beta_cols                   <csv>    (default: Ast_b_gi,Ex_b_gi,In_b_gi,Oli_b_gi)
#   --se_cols                     <csv>    (default: Ast_b_gi_se,Ex_b_gi_se,In_b_gi_se,Oli_b_gi_se)
#   --pca_components              <int>    (default: 4)
#
# Example:
# Rscript MashR_analysis.R \
#   --random_set /path/to/random.txt \
#   --strong_set /path/to/strong.txt \
#   --out_rds /path/to/results/strong_mash.rds
# ============================================================

suppressPackageStartupMessages({
  library(ashr)
  library(mashr)
  library(data.table)
})

# -----------------------------
# Minimal CLI parser (base R)
# -----------------------------
parse_args <- function(argv) {
  if (length(argv) == 0) return(list())

  args <- list()
  i <- 1
  while (i <= length(argv)) {
    key <- argv[i]
    if (!startsWith(key, "--")) {
      stop(sprintf("Unexpected token '%s'. Arguments must be in '--key value' form.", key), call. = FALSE)
    }
    if (i == length(argv)) {
      stop(sprintf("Missing value for '%s'.", key), call. = FALSE)
    }
    value <- argv[i + 1]
    args[[sub("^--", "", key)]] <- value
    i <- i + 2
  }
  args
}

require_arg <- function(args, name) {
  if (is.null(args[[name]]) || args[[name]] == "") {
    stop(sprintf("Missing required argument: --%s", name), call. = FALSE)
  }
  args[[name]]
}

file_must_exist <- function(path) {
  if (!file.exists(path)) stop(sprintf("File not found: %s", path), call. = FALSE)
  invisible(TRUE)
}

parse_csv <- function(x) {
  trimws(strsplit(x, ",", fixed = TRUE)[[1]])
}

# -----------------------------
# Extract matrices for mashr
# -----------------------------
extract_effect_and_se <- function(dt, effect_cols, se_cols, row_id_col) {
  required_cols <- c(effect_cols, se_cols, row_id_col)
  missing_cols <- setdiff(required_cols, names(dt))
  if (length(missing_cols) > 0) {
    stop(sprintf("Missing required columns: %s", paste(missing_cols, collapse = ", ")), call. = FALSE)
  }

  effect_matrix <- as.matrix(dt[, ..effect_cols])
  se_matrix     <- as.matrix(dt[, ..se_cols])

  rownames(effect_matrix) <- dt[[row_id_col]]
  rownames(se_matrix)     <- dt[[row_id_col]]

  if (!all(dim(effect_matrix) == dim(se_matrix))) {
    stop("Effect and SE matrices have mismatched dimensions.", call. = FALSE)
  }
  list(effect = effect_matrix, se = se_matrix)
}

# -----------------------------
# Main
# -----------------------------
args <- parse_args(commandArgs(trailingOnly = TRUE))

random_path <- require_arg(args, "random_set")
strong_path <- require_arg(args, "strong_set")
output_rds_path            <- require_arg(args, "out_rds")

lfsr_threshold <- as.numeric(ifelse(is.null(args[["lfsr_threshold"]]), "0.05", args[["lfsr_threshold"]]))
pair_id_col    <- ifelse(is.null(args[["pair_id_col"]]), "Pair", args[["pair_id_col"]])

effect_cols <- if (is.null(args[["beta_cols"]])) {
  c("Ast_b_gi", "Ex_b_gi", "In_b_gi", "Oli_b_gi")
} else {
  parse_csv(args[["beta_cols"]])
}

se_cols <- if (is.null(args[["se_cols"]])) {
  c("Ast_b_gi_se", "Ex_b_gi_se", "In_b_gi_se", "Oli_b_gi_se")
} else {
  parse_csv(args[["se_cols"]])
}

pca_components <- as.integer(ifelse(is.null(args[["pca_components"]]), "4", args[["pca_components"]]))

# Validate inputs
file_must_exist(random_path)
file_must_exist(strong_path)
dir.create(dirname(output_rds_path), recursive = TRUE, showWarnings = FALSE)

message("Loading protein-coding random set: ", random_path)
random_set <- fread(random_path)

message("Loading protein-coding strong set: ", strong_path)
strong_set <- fread(strong_path)

message("Extracting effect sizes and standard errors...")
random_matrices <- extract_effect_and_se(
  dt = random_set,
  effect_cols = effect_cols,
  se_cols = se_cols,
  row_id_col = pair_id_col
)
strong_matrices <- extract_effect_and_se(
  dt = strong_set,
  effect_cols = effect_cols,
  se_cols = se_cols,
  row_id_col = pair_id_col
)

random_effect <- random_matrices$effect
random_se     <- random_matrices$se
strong_effect <- strong_matrices$effect
strong_se     <- strong_matrices$se

# Basic sanity checks
if (ncol(random_effect) != length(effect_cols)) stop("Unexpected number of conditions in random effect matrix.", call. = FALSE)
if (ncol(strong_effect) != length(effect_cols)) stop("Unexpected number of conditions in strong effect matrix.", call. = FALSE)

# 1) Estimate null correlation (Vhat) from random set
message("Estimating null correlation matrix (Vhat) from random set...")
temporary_data_object <- mash_set_data(random_effect, random_se)
null_correlation_matrix <- estimate_null_correlation_simple(temporary_data_object)
rm(temporary_data_object)

message("Null correlation matrix (Vhat):")
print(null_correlation_matrix)

# 2) Build mash data objects with Vhat
random_data_object <- mash_set_data(random_effect, random_se, V = null_correlation_matrix)
strong_data_object <- mash_set_data(strong_effect, strong_se, V = null_correlation_matrix)

# 3) Construct covariance matrices
message("Constructing covariance matrices for mashr...")
pca_covariances <- cov_pca(strong_data_object, npc = pca_components)
refined_covariances <- cov_ed(strong_data_object, pca_covariances)
canonical_covariances <- cov_canonical(random_data_object)

covariance_list <- c(refined_covariances, canonical_covariances)

# 4) Fit mash on random set to learn prior g
message("Fitting mash model on random set (learning prior g)...")
random_fit <- mash(random_data_object, Ulist = covariance_list, outputlevel = 1)

# 5) Fit strong set with fixed prior g
message("Fitting mash model on strong set (fixing g learned from random set)...")
strong_fit <- mash(
  strong_data_object,
  g = get_fitted_g(random_fit),
  fixg = TRUE,
  verbose = TRUE
)

# 6) Identify significant results and cell-type-specific ones
message("Extracting significant results using LFSR threshold: ", lfsr_threshold)

significant_indices <- get_significant_results(
  strong_fit,
  thresh = lfsr_threshold,
  condition = NULL,
  sig_fn = get_lfsr
)

lfsr_matrix <- get_lfsr(strong_fit)
significant_lfsr <- lfsr_matrix[significant_indices, , drop = FALSE]

cell_type_specific_lfsr <- significant_lfsr[
  apply(significant_lfsr, 1, function(row) sum(row < lfsr_threshold) == 1),
  ,
  drop = FALSE
]

message(sprintf("Number of significant pairs (any condition, LFSR < %.3g): %d",
                lfsr_threshold, nrow(significant_lfsr)))
message(sprintf("Number of cell-type-specific pairs (exactly one condition passes): %d",
                nrow(cell_type_specific_lfsr)))

message("Preview of cell-type-specific pairs:")
print(utils::head(cell_type_specific_lfsr))

# 7) Save results
saveRDS(strong_fit, file = output_rds_path)
message("Saved mash fit object to: ", output_rds_path)