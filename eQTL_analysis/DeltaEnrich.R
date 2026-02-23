#!/usr/bin/env Rscript

# ============================================================
# Peak enrichment test by MAF-TSS matched sampling
#
# This script:
#   1) Loads ieQTL table and selects one top SNP per gene (min p)
#   2) Joins SNP annotations (maf_tss_index, in_peak) by SNP_probe
#   3) Computes observed peak overlap proportion for query set
#   4) Builds a null distribution by repeated sampling from background
#      with matched maf_tss_index composition
#   5) Saves null hits, and a summary log including fold enrichment and CI
#
# Required inputs:
#   --ieqtl_file        path to ieQTL table
#   --anno_file         path to SNP annotation table (must include SNP_probe, maf_tss_index, in_peak)
#   --out_null          output file for null distribution (one column)
#   --out_log           output file for summary (one row TSV)
#
# Query filters:
#   --cell_type         e.g. Oli
#   --gene_type         e.g. LncRNA
#   --direction         e.g. Positive
#
# Optional:
#   --rep_num           number of repeats (default: 1000)
#   --seed              random seed (default: 20251028)
#   --pair_sep          separator for SNP_probe (default: "_")
#   --allow_replace     whether to allow sampling with replacement if bin is insufficient
#                       (default: TRUE)
#
# Example:
# Rscript peak_enrichment_matched_sampling.R \
#   --ieqtl_file /path/ieQTL.txt \
#   --anno_file  /path/Oli_Glia_SNP_Annotation.txt \
#   --cell_type Oli --gene_type LncRNA --direction Positive \
#   --rep_num 1000 --seed 20251028 \
#   --out_null /path/null_hits.txt \
#   --out_log  /path/summary_log.tsv
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
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
      stop(sprintf("Unexpected token '%s'. Use '--key value' pairs.", key), call. = FALSE)
    }
    if (i == length(argv)) stop(sprintf("Missing value for '%s'.", key), call. = FALSE)
    args[[sub("^--", "", key)]] <- argv[i + 1]
    i <- i + 2
  }
  args
}

require_arg <- function(args, name) {
  val <- args[[name]]
  if (is.null(val) || val == "") stop(sprintf("Missing required argument: --%s", name), call. = FALSE)
  val
}

as_int <- function(x, name) {
  v <- suppressWarnings(as.integer(x))
  if (!is.finite(v)) stop(sprintf("Argument --%s must be an integer. Got: %s", name, x), call. = FALSE)
  v
}

as_num <- function(x, name) {
  v <- suppressWarnings(as.numeric(x))
  if (!is.finite(v)) stop(sprintf("Argument --%s must be numeric. Got: %s", name, x), call. = FALSE)
  v
}

as_bool <- function(x, name) {
  if (tolower(x) %in% c("true", "t", "1", "yes", "y")) return(TRUE)
  if (tolower(x) %in% c("false", "f", "0", "no", "n")) return(FALSE)
  stop(sprintf("Argument --%s must be boolean (true/false). Got: %s", name, x), call. = FALSE)
}

file_must_exist <- function(path) {
  if (!file.exists(path)) stop(sprintf("File not found: %s", path), call. = FALSE)
  invisible(TRUE)
}

dir_must_exist <- function(path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  invisible(TRUE)
}

# -----------------------------
# Delta method for ratio x / ybar
# -----------------------------
# We treat x (observed proportion) as fixed and ybar = mean(null) as random.
# Var(ybar) = Var(y) / m.
# Then Var(x / ybar) ≈ (x^2 / ybar^4) * Var(ybar).
# If you want to incorporate Var(x) as well, you can extend this function.
ratio_se_delta <- function(x, ybar, var_y, m) {
  var_ybar <- var_y / m
  var_ratio <- (x^2 / (ybar^4)) * var_ybar
  sqrt(var_ratio)
}

# -----------------------------
# Main
# -----------------------------
args <- parse_args(commandArgs(trailingOnly = TRUE))

ieqtl_file  <- require_arg(args, "ieqtl_file")
anno_file   <- require_arg(args, "anno_file")
out_null    <- require_arg(args, "out_null")
out_log     <- require_arg(args, "out_log")

cell_type   <- require_arg(args, "cell_type")
gene_type   <- require_arg(args, "gene_type")
direction   <- require_arg(args, "direction")

rep_num     <- if (is.null(args[["rep_num"]])) 1000L else as_int(args[["rep_num"]], "rep_num")
seed        <- if (is.null(args[["seed"]])) 20251028L else as_int(args[["seed"]], "seed")
pair_sep    <- if (is.null(args[["pair_sep"]])) "_" else args[["pair_sep"]]
allow_replace <- if (is.null(args[["allow_replace"]])) TRUE else as_bool(args[["allow_replace"]], "allow_replace")

set.seed(seed)

file_must_exist(ieqtl_file)
file_must_exist(anno_file)
dir_must_exist(out_null)
dir_must_exist(out_log)

message(">> Loading ieQTL table: ", ieqtl_file)
ieqtl_raw <- fread(ieqtl_file)

# Validate ieQTL required columns
required_ieqtl_cols <- c(
  "CellType", "phenotype_id", "variant_id",
  "pval_gi", "b_gi", "b_gi_se", "pval_adj_bh",
  "eGeneType", "SourceGene", "Direction"
)
missing_ieqtl <- setdiff(required_ieqtl_cols, names(ieqtl_raw))
if (length(missing_ieqtl) > 0) {
  stop("ieQTL file is missing required columns: ", paste(missing_ieqtl, collapse = ", "), call. = FALSE)
}

message(">> Filtering query set: cell_type=", cell_type,
        ", gene_type=", gene_type, ", direction=", direction)

query_ieqtl <- ieqtl_raw %>%
  select(CellType, phenotype_id, variant_id, start_distance, af,
         pval_gi, b_gi, b_gi_se, pval_adj_bh, eGeneType, SourceGene, Direction) %>%
  rename(
    Gene   = phenotype_id,
    TopSNP = variant_id,
    p      = pval_gi,
    b      = b_gi,
    se     = b_gi_se,
    FDR    = pval_adj_bh
  ) %>%
  filter(
    CellType == cell_type,
    eGeneType == gene_type,
    Direction == direction
  ) %>%
  as.data.table()

if (nrow(query_ieqtl) == 0) {
  stop("No ieQTL records matched the query filters. Please check cell_type/gene_type/direction.", call. = FALSE)
}

# Keep one SNP per gene (smallest p-value)
query_ieqtl <- query_ieqtl[, .SD[which.min(p)], by = Gene]
query_ieqtl[, SNP_probe := paste0(Gene, pair_sep, TopSNP)]

message(">> Query genes: ", uniqueN(query_ieqtl$Gene),
        " | Query pairs: ", nrow(query_ieqtl))

message(">> Loading annotation table: ", anno_file)
annotation_table <- fread(anno_file)

required_anno_cols <- c("SNP_probe", "maf_tss_index", "in_peak")
missing_anno <- setdiff(required_anno_cols, names(annotation_table))
if (length(missing_anno) > 0) {
  stop("Annotation file is missing required columns: ", paste(missing_anno, collapse = ", "), call. = FALSE)
}

# Join annotation into query set
query_ieqtl <- query_ieqtl %>%
  left_join(annotation_table[, .(SNP_probe, maf_tss_index, in_peak)], by = "SNP_probe") %>%
  as.data.table()

# Drop query SNPs without bin info (cannot be matched)
query_ieqtl <- query_ieqtl[!is.na(maf_tss_index)]
if (nrow(query_ieqtl) == 0) {
  stop("All query SNPs are missing maf_tss_index after join. Cannot proceed.", call. = FALSE)
}

observed_hit_rate <- mean(query_ieqtl$in_peak, na.rm = TRUE)
query_bin_counts <- table(query_ieqtl$maf_tss_index)
query_bin_counts <- query_bin_counts[query_bin_counts > 0]

message(">> Observed peak overlap proportion: ", round(observed_hit_rate, 6))
message(">> Number of bins in query: ", length(query_bin_counts))

# Background pool: all annotated SNPs excluding query SNP_probe
background_pool <- annotation_table[!(SNP_probe %in% query_ieqtl$SNP_probe)]
background_pool <- background_pool[!is.na(maf_tss_index)]

# Prepare bin requirements
bin_requirement <- data.table(
  maf_tss_index = names(query_bin_counts),
  need = as.integer(query_bin_counts)
)

bin_available <- background_pool[, .(available = .N), by = maf_tss_index]
bin_info <- merge(bin_requirement, bin_available, by = "maf_tss_index", all.x = TRUE)
bin_info[is.na(available), available := 0L]

# Warn on bins missing from background
missing_bins <- bin_info[available == 0L, maf_tss_index]
if (length(missing_bins) > 0) {
  message("!! Warning: The following bins exist in query but have zero background SNPs: ",
          paste(missing_bins, collapse = ", "),
          ". These bins will be skipped in sampling and reduce effective matched size.")
}

# Matched sampling
random_hit_rates <- numeric(rep_num)

for (i in seq_len(rep_num)) {
  sampled_list <- vector("list", nrow(bin_info))

  for (k in seq_len(nrow(bin_info))) {
    b <- bin_info$maf_tss_index[k]
    n <- bin_info$need[k]
    a <- bin_info$available[k]

    if (a <= 0L) next

    candidates <- background_pool[maf_tss_index == b]
    if (a < n && !allow_replace) {
      # Not enough candidates and replacement is not allowed -> skip this bin
      next
    }

    sampled_list[[k]] <- candidates[sample(.N, n, replace = (a < n))]
  }

  sampled_all <- rbindlist(Filter(Negate(is.null), sampled_list))

  # If everything got skipped (rare but possible), mark as NA
  random_hit_rates[i] <- if (nrow(sampled_all) == 0) NA_real_ else mean(sampled_all$in_peak, na.rm = TRUE)
}

# Save null distribution
fwrite(data.table(random_hits = random_hit_rates), file = out_null, sep = "\t", col.names = FALSE)
message(">> Null distribution saved: ", out_null)

# Summarize null
valid_null <- random_hit_rates[is.finite(random_hit_rates)]
if (length(valid_null) < 10) {
  stop("Too few valid null samples were generated. Check background bins / allow_replace.", call. = FALSE)
}

null_mean <- mean(valid_null)
null_var  <- var(valid_null)
m_eff     <- length(valid_null)

fold_enrichment <- observed_hit_rate / null_mean
se_fold <- ratio_se_delta(observed_hit_rate, null_mean, null_var, m_eff)
ci_low  <- fold_enrichment - 1.96 * se_fold
ci_high <- fold_enrichment + 1.96 * se_fold

# Print summary
cat("==== Peak Enrichment (Matched Sampling) ====\n")
cat("Cell type   :", cell_type, "\n")
cat("Peak file   :", anno_file, "\n")
cat("Gene type   :", gene_type, "\n")
cat("Direction   :", direction, "\n")
cat("Repeats     :", rep_num, " (valid:", m_eff, ")\n")
cat("Observed hit rate     :", round(observed_hit_rate, 6), "\n")
cat("Null mean hit rate    :", round(null_mean, 6), "\n")
cat("Fold enrichment       :", round(fold_enrichment, 6), "\n")
cat("SE (delta method)     :", round(se_fold, 6), "\n")
cat("95% CI                : [", round(ci_low, 6), ", ", round(ci_high, 6), "]\n", sep = "")

# Save log table
log_table <- data.table(
  cell_type = cell_type,
  gene_type = gene_type,
  direction = direction,
  rep_num   = rep_num,
  seed      = seed,
  allow_replace = allow_replace,
  query_n_pairs = nrow(query_ieqtl),
  query_n_genes = uniqueN(query_ieqtl$Gene),
  observed_hit_rate = observed_hit_rate,
  null_mean_hit_rate = null_mean,
  null_var_hit_rate  = null_var,
  valid_repeats      = m_eff,
  fold_enrichment    = fold_enrichment,
  se_fold            = se_fold,
  ci_low             = ci_low,
  ci_high            = ci_high
)

fwrite(log_table, file = out_log, sep = "\t", col.names = TRUE)
message(">> Summary log saved: ", out_log)