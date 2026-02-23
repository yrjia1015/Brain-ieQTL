#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(coloc)
  library(foreach)
  library(doParallel)
})

# ============================================================
# Generic coloc pipeline (GWAS vs eQTL) with chr-level parallelization
#
# Required:
#   --trait             Trait name (used only for output naming)
#   --cell              Cell name (used only for output naming)
#   --trait_type        "cc" or "quant"
#   --gwas_file         GWAS summary file (single file with all chr)
#   --eqtl_pattern      eQTL file pattern with {cell} and {chr}
#                       e.g. "/path/eqtl/{cell}/{cell}_chr{chr}.txt"
#   --out_dir           Output root directory
#
# Trait meta:
#   If --trait_type cc:
#       --cases --controls  (required)
#   If --trait_type quant:
#       --gwas_N            (required if GWAS file lacks N column)
#
# Optional:
#   --eqtl_N              eQTL sample size (default: 2443)
#   --n_cores             parallel cores (default: 8)
#   --tss_window          abs(tss_distance) filter (default: 1e6)
#   --min_snps            minimum overlapping SNPs per gene to run coloc (default: 50)
#   --pp4_threshold       threshold to flag/keep results (default: 0.0; keep all)
#   --gwas_N_col          column name for per-SNP N in GWAS (default: "N")
#   --gwas_freq_col       column name for GWAS A1 frequency (default: auto)
#   --log_dir             log directory (default: <out_dir>/logs/<cell>/<trait>)
#
# Expected columns:
#   GWAS: SNP, beta, varbeta, A1, A2, p (or P), se (or SE), N (optional), frq (optional)
#   eQTL: SNP, Probe, beta, varbeta, A1, A2, p, tss_distance, Chr, Probe_bp, Freq
#
# Notes:
#   - This script does a simple allele flip if A1/A2 are swapped between datasets.
#   - Strand complements are NOT handled here (add if needed).
# ============================================================

# -----------------------------
# CLI parsing (optparse-like minimal)
# -----------------------------
parse_args <- function(argv) {
  if (length(argv) == 0) return(list())
  args <- list()
  i <- 1
  while (i <= length(argv)) {
    key <- argv[i]
    if (!startsWith(key, "--")) stop("Arguments must be provided as --key value pairs.", call. = FALSE)
    if (i == length(argv)) stop(sprintf("Missing value for %s", key), call. = FALSE)
    args[[sub("^--", "", key)]] <- argv[i + 1]
    i <- i + 2
  }
  args
}

require_arg <- function(args, name) {
  v <- args[[name]]
  if (is.null(v) || v == "") stop(sprintf("Missing required argument: --%s", name), call. = FALSE)
  v
}

as_int <- function(x, name) {
  v <- suppressWarnings(as.integer(x))
  if (!is.finite(v)) stop(sprintf("--%s must be an integer. Got: %s", name, x), call. = FALSE)
  v
}

as_num <- function(x, name) {
  v <- suppressWarnings(as.numeric(x))
  if (!is.finite(v)) stop(sprintf("--%s must be numeric. Got: %s", name, x), call. = FALSE)
  v
}

as_bool <- function(x, name) {
  if (tolower(x) %in% c("true", "t", "1", "yes", "y")) return(TRUE)
  if (tolower(x) %in% c("false", "f", "0", "no", "n")) return(FALSE)
  stop(sprintf("--%s must be boolean. Got: %s", name, x), call. = FALSE)
}

stop_if_missing_file <- function(path) {
  if (!file.exists(path)) stop("File not found: ", path, call. = FALSE)
}

mkdir_parent <- function(path) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
}

# -----------------------------
# Helpers: column standardization
# -----------------------------
standardize_gwas_cols <- function(gwas, gwas_freq_col = NULL, gwas_N_col = "N") {
  # Required core columns
  required <- c("SNP", "beta", "varbeta", "A1", "A2")
  missing <- setdiff(required, names(gwas))
  if (length(missing) > 0) stop("GWAS missing required columns: ", paste(missing, collapse = ", "), call. = FALSE)

  # p-value column
  if (!("p" %in% names(gwas))) {
    if ("P" %in% names(gwas)) setnames(gwas, "P", "p")
  }
  if (!("p" %in% names(gwas))) stop("GWAS missing p/P column.", call. = FALSE)

  # SE column is not required if varbeta exists; keep compatibility anyway
  if (!("se" %in% names(gwas))) {
    if ("SE" %in% names(gwas)) setnames(gwas, "SE", "se")
  }

  # Frequency column (A1 frequency): optional
  if (is.null(gwas_freq_col) || gwas_freq_col == "" || is.na(gwas_freq_col)) {
    # Try common names
    candidates <- c("frq_A1", "FRQ", "Freq", "freq", "EAF", "eaf", "AF", "af", "A1FREQ", "freqA1")
    freq_hit <- candidates[candidates %in% names(gwas)]
    if (length(freq_hit) > 0) {
      setnames(gwas, freq_hit[1], "frq_A1")
    } else {
      # keep absent; only needed for quant case if you want MAF
      # coloc.abf can work without MAF if beta/varbeta/p are consistent, but MAF is recommended
      # We'll not force it, just warn later if needed.
      gwas[, frq_A1 := NA_real_]
    }
  } else {
    if (!(gwas_freq_col %in% names(gwas))) stop("GWAS freq column not found: ", gwas_freq_col, call. = FALSE)
    setnames(gwas, gwas_freq_col, "frq_A1")
  }

  # N column: optional; if present, standardize to N
  if (!(gwas_N_col %in% names(gwas))) {
    # ok to be absent; handled by --gwas_N
  } else if (gwas_N_col != "N") {
    setnames(gwas, gwas_N_col, "N")
  }

  gwas[]
}

standardize_eqtl_cols <- function(eqtl) {
  required <- c("SNP", "Probe", "beta", "varbeta", "A1", "A2", "p", "tss_distance", "Chr", "Probe_bp", "Freq")
  missing <- setdiff(required, names(eqtl))
  if (length(missing) > 0) stop("eQTL missing required columns: ", paste(missing, collapse = ", "), call. = FALSE)
  eqtl[]
}

# -----------------------------
# Core coloc per gene
# -----------------------------
run_coloc_for_gene <- function(merged, eqtl_N, trait_type, cases = NA, controls = NA, gwas_N_fallback = NA_real_) {
  # Filter duplicates and window already applied outside; assume merged contains:
  # p_eqtl, beta_eqtl, varbeta_eqtl, Freq, p_gwas, beta_gwas, varbeta_gwas, frq_A1, SNP, etc.

  # Allele flip if A1/A2 swapped between datasets
  flip <- merged$A1_eqtl == merged$A2_gwas & merged$A2_eqtl == merged$A1_gwas
  if (any(flip, na.rm = TRUE)) {
    merged$frq_A1[flip] <- 1 - merged$frq_A1[flip]
    merged$beta_gwas[flip] <- -merged$beta_gwas[flip]
  }

  eqtl_list <- list(
    pvalues = merged$p_eqtl,
    beta    = merged$beta_eqtl,
    varbeta = merged$varbeta_eqtl,
    snp     = merged$SNP,
    MAF     = merged$Freq,
    N       = eqtl_N,
    type    = "quant"
  )

  if (trait_type == "cc") {
    s <- cases / (cases + controls)
    gwas_N <- if ("N" %in% names(merged)) max(merged$N, na.rm = TRUE) else gwas_N_fallback
    if (!is.finite(gwas_N)) stop("GWAS N is missing (no N column and --gwas_N not provided).", call. = FALSE)

    gwas_list <- list(
      pvalues = merged$p_gwas,
      beta    = merged$beta_gwas,
      varbeta = merged$varbeta_gwas,
      s       = s,
      N       = gwas_N,
      snp     = merged$SNP,
      type    = "cc"
    )
  } else {
    gwas_N <- if ("N" %in% names(merged)) max(merged$N, na.rm = TRUE) else gwas_N_fallback
    if (!is.finite(gwas_N)) stop("GWAS N is missing (no N column and --gwas_N not provided).", call. = FALSE)

    gwas_list <- list(
      pvalues = merged$p_gwas,
      beta    = merged$beta_gwas,
      varbeta = merged$varbeta_gwas,
      MAF     = merged$frq_A1,
      N       = gwas_N,
      snp     = merged$SNP,
      type    = "quant"
    )
  }

  suppressWarnings({
    res <- coloc.abf(dataset1 = eqtl_list, dataset2 = gwas_list)
  })

  as.data.frame(t(res$summary))
}

# -----------------------------
# Main parallel runner
# -----------------------------
run_coloc_parallel <- function(
  trait, cell,
  trait_type,
  gwas_file, eqtl_pattern,
  out_dir,
  eqtl_N = 2443,
  gwas_N = NA_real_,
  cases = NA, controls = NA,
  tss_window = 1e6,
  min_snps = 50,
  pp4_threshold = 0.0,
  allow_replace = TRUE, # placeholder (not used here, kept for symmetry)
  n_cores = 8,
  gwas_N_col = "N",
  gwas_freq_col = NULL,
  log_dir = NULL
) {
  stop_if_missing_file(gwas_file)

  if (!(trait_type %in% c("cc", "quant"))) stop("--trait_type must be 'cc' or 'quant'.", call. = FALSE)
  if (trait_type == "cc") {
    if (!is.finite(cases) || !is.finite(controls)) {
      stop("For cc traits, you must provide --cases and --controls.", call. = FALSE)
    }
  } else {
    # quant
    # gwas_N can be missing only if GWAS file includes per-SNP N
    # We'll validate later after reading GWAS.
  }

  if (is.null(log_dir) || is.na(log_dir) || log_dir == "") {
    log_dir <- file.path(out_dir, "logs", cell, trait)
  }
  mkdir_parent(log_dir)

  mkdir_parent(out_dir)

  # Load GWAS once in master (but each worker can also read; we broadcast to workers to reduce IO)
  message(">> Loading GWAS: ", gwas_file)
  gwas <- fread(gwas_file)
  gwas <- standardize_gwas_cols(gwas, gwas_freq_col = gwas_freq_col, gwas_N_col = gwas_N_col)

  # If quant and no N column, need --gwas_N
  if (trait_type == "quant" && !("N" %in% names(gwas)) && !is.finite(gwas_N)) {
    stop("Quant trait requires GWAS N: provide --gwas_N if GWAS file does not have an N column.", call. = FALSE)
  }
  if (trait_type == "cc" && !("N" %in% names(gwas)) && !is.finite(gwas_N)) {
    stop("CC trait requires GWAS N: provide --gwas_N if GWAS file does not have an N column.", call. = FALSE)
  }

  # Setup cluster
  cl <- makeCluster(n_cores)
  on.exit(stopCluster(cl), add = TRUE)
  registerDoParallel(cl)

  foreach(chr = 1:22,
          .packages = c("data.table", "coloc"),
          .export = c("standardize_eqtl_cols", "run_coloc_for_gene")) %dopar% {

    logf <- file.path(log_dir, sprintf("%s_%s_chr%s.log", trait, cell, chr))
    cat(sprintf("[%s] Start chr%s\n", Sys.time(), chr), file = logf, append = TRUE)

    eqtl_file <- gsub("\\{cell\\}", cell, eqtl_pattern)
    eqtl_file <- gsub("\\{chr\\}", as.character(chr), eqtl_file)

    if (!file.exists(eqtl_file)) {
      cat(sprintf("[%s] Missing eQTL file: %s\n", Sys.time(), eqtl_file), file = logf, append = TRUE)
      return(NULL)
    }

    eqtl <- fread(eqtl_file)
    eqtl <- standardize_eqtl_cols(eqtl)

    genes <- unique(eqtl$Probe)
    cat(sprintf("[%s] Loaded eQTL chr%s: %d genes\n", Sys.time(), chr, length(genes)),
        file = logf, append = TRUE)

    # Pre-filter by tss window early to reduce work
    eqtl <- eqtl[abs(tss_distance) <= tss_window]

    # Results table
    chr_results <- data.table(
      Gene = character(),
      eGene_chr = integer(),
      eGene_pos = numeric(),
      eGene_topSNP = character(),
      eGene_topSNP_p = numeric(),
      GWAS_topSNP = character(),
      GWAS_topSNP_p = numeric(),
      nsnp = integer(),
      PP.H0 = numeric(), PP.H1 = numeric(), PP.H2 = numeric(), PP.H3 = numeric(), PP.H4 = numeric()
    )

    # For speed: key GWAS by SNP once
    setkey(gwas, SNP)

    for (i in seq_along(genes)) {
      gene <- genes[i]
      if (i %% 50 == 0) {
        cat(sprintf("[%s] chr%s progress: %d/%d genes\n", Sys.time(), chr, i, length(genes)),
            file = logf, append = TRUE)
      }

      eqtl_g <- eqtl[Probe == gene]
      if (nrow(eqtl_g) == 0) next

      merged <- merge(eqtl_g, gwas, by = "SNP", suffixes = c("_eqtl", "_gwas"))
      if (nrow(merged) < min_snps) next

      # Remove duplicated SNP if any
      merged <- merged[!duplicated(SNP)]

      # coloc requires complete vectors; drop NA in key stats
      merged <- merged[is.finite(p_eqtl) & is.finite(beta_eqtl) & is.finite(varbeta_eqtl) &
                         is.finite(p_gwas) & is.finite(beta_gwas) & is.finite(varbeta_gwas)]

      if (nrow(merged) < min_snps) next

      # Run coloc
      coloc_sum <- tryCatch(
        run_coloc_for_gene(
          merged = merged,
          eqtl_N = eqtl_N,
          trait_type = trait_type,
          cases = cases,
          controls = controls,
          gwas_N_fallback = gwas_N
        ),
        error = function(e) {
          cat(sprintf("[%s] ERROR gene=%s chr%s: %s\n", Sys.time(), gene, chr, e$message),
              file = logf, append = TRUE)
          return(NULL)
        }
      )
      if (is.null(coloc_sum)) next

      # Extract top SNPs
      eqtl_top_idx <- which.min(merged$p_eqtl)
      gwas_top_idx <- which.min(merged$p_gwas)

      new_row <- data.table(
        Gene = merged$Probe[1],
        eGene_chr = as.integer(merged$Chr[1]),
        eGene_pos = as.numeric(merged$Probe_bp[1]),
        eGene_topSNP = merged$SNP[eqtl_top_idx],
        eGene_topSNP_p = min(merged$p_eqtl, na.rm = TRUE),
        GWAS_topSNP = merged$SNP[gwas_top_idx],
        GWAS_topSNP_p = min(merged$p_gwas, na.rm = TRUE),
        nsnp = as.integer(coloc_sum$nsnps),
        PP.H0 = as.numeric(coloc_sum$PP.H0.abf),
        PP.H1 = as.numeric(coloc_sum$PP.H1.abf),
        PP.H2 = as.numeric(coloc_sum$PP.H2.abf),
        PP.H3 = as.numeric(coloc_sum$PP.H3.abf),
        PP.H4 = as.numeric(coloc_sum$PP.H4.abf)
      )

      if (pp4_threshold <= 0 || (is.finite(new_row$PP.H4) && new_row$PP.H4 >= pp4_threshold)) {
        chr_results <- rbind(chr_results, new_row)
      }
    }

    # Save chr results
    out_subdir <- file.path(out_dir, cell, trait)
    dir.create(out_subdir, recursive = TRUE, showWarnings = FALSE)
    out_file <- file.path(out_subdir, sprintf("%s_%s_chr%s.csv", trait, cell, chr))
    fwrite(chr_results, file = out_file)

    cat(sprintf("[%s] Saved: %s (%d rows)\n", Sys.time(), out_file, nrow(chr_results)),
        file = logf, append = TRUE)

    out_file
  }
}

# -----------------------------
# Entry point
# -----------------------------
args <- parse_args(commandArgs(trailingOnly = TRUE))

trait        <- require_arg(args, "trait")
cell         <- require_arg(args, "cell")
trait_type   <- require_arg(args, "trait_type")

gwas_file    <- require_arg(args, "gwas_file")
eqtl_pattern <- require_arg(args, "eqtl_pattern")
out_dir      <- require_arg(args, "out_dir")

eqtl_N       <- if (is.null(args[["eqtl_N"]])) 2443 else as_int(args[["eqtl_N"]], "eqtl_N")
n_cores      <- if (is.null(args[["n_cores"]])) 8 else as_int(args[["n_cores"]], "n_cores")
tss_window   <- if (is.null(args[["tss_window"]])) 1e6 else as_num(args[["tss_window"]], "tss_window")
min_snps     <- if (is.null(args[["min_snps"]])) 50 else as_int(args[["min_snps"]], "min_snps")
pp4_thr      <- if (is.null(args[["pp4_threshold"]])) 0.0 else as_num(args[["pp4_threshold"]], "pp4_threshold")
seed         <- if (is.null(args[["seed"]])) 20251028 else as_int(args[["seed"]], "seed")

gwas_N       <- if (is.null(args[["gwas_N"]])) NA_real_ else as_num(args[["gwas_N"]], "gwas_N")
cases        <- if (is.null(args[["cases"]])) NA_real_ else as_num(args[["cases"]], "cases")
controls     <- if (is.null(args[["controls"]])) NA_real_ else as_num(args[["controls"]], "controls")

gwas_N_col   <- if (is.null(args[["gwas_N_col"]])) "N" else args[["gwas_N_col"]]
gwas_freq_col<- if (is.null(args[["gwas_freq_col"]])) NULL else args[["gwas_freq_col"]]
log_dir      <- if (is.null(args[["log_dir"]])) NULL else args[["log_dir"]]

set.seed(seed)

message(">> Start coloc")
message("   trait      = ", trait)
message("   cell       = ", cell)
message("   trait_type = ", trait_type)
message("   gwas_file  = ", gwas_file)
message("   eqtl_pattern = ", eqtl_pattern)
message("   out_dir    = ", out_dir)
message("   n_cores    = ", n_cores)

run_coloc_parallel(
  trait = trait,
  cell = cell,
  trait_type = trait_type,
  gwas_file = gwas_file,
  eqtl_pattern = eqtl_pattern,
  out_dir = out_dir,
  eqtl_N = eqtl_N,
  gwas_N = gwas_N,
  cases = cases,
  controls = controls,
  tss_window = tss_window,
  min_snps = min_snps,
  pp4_threshold = pp4_thr,
  n_cores = n_cores,
  gwas_N_col = gwas_N_col,
  gwas_freq_col = gwas_freq_col,
  log_dir = log_dir
)

message(">> Done")