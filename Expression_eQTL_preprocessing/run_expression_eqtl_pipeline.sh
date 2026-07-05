#!/usr/bin/env bash
set -euo pipefail

# ============================================================
# eQTL expression preprocessing and tensorQTL command template
# ============================================================
#
# This file is a cleaned, reusable version of the historical commands used for
# HBCC/new lncRNA and lncRNAKB analyses. Edit the variables in the "User
# settings" block, then run sections step by step or submit the whole script.
#
# Main workflow:
#   1. Filter lowly expressed genes using TPM and raw count matrices.
#   2. Run edgeR TMM normalization and output log2 CPM.
#   3. Convert gene annotation GTF to BED.
#   4. Subset/relabel samples if needed.
#   5. Apply rank inverse normal transformation (INT/RINT).
#   6. Convert expression matrix to tensorQTL phenotype BED.
#   7. bgzip/tabix the phenotype BED.
#   8. Estimate PEER factors and residual expression.
#   9. Build covariates and run tensorQTL.
#
# Tools expected in PATH:
#   perl, Rscript, sort, bgzip, tabix, zcat, python
#
# R packages:
#   edgeR, peer, argparser
#
# Python packages for tensorQTL:
#   pandas, numpy, torch, tensorqtl, pyarrow
# ============================================================

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# -----------------------------
# User settings
# -----------------------------
PROJECT_DIR="${PROJECT_DIR:-$PWD}"
PREFIX="${PREFIX:-newlncRNA_exp}"

TPM_MATRIX="${TPM_MATRIX:-/path/to/gencode.genes.collapsed.lncRNA.TPM.txt}"
COUNT_MATRIX="${COUNT_MATRIX:-/path/to/gencode.genes.collapsed.lncRNA.fragments.txt}"
GTF_FILE="${GTF_FILE:-/path/to/gencode.genes.collapsed.lncRNA.gtf}"

# Optional sample files. SAMPLE_KEEP contains the expression-matrix sample IDs
# to keep. SAMPLE_RENAME contains replacement sample IDs in the same order.
SAMPLE_KEEP="${SAMPLE_KEEP:-}"
SAMPLE_RENAME="${SAMPLE_RENAME:-}"

# PLINK prefix used by tensorQTL, without .bed/.bim/.fam suffix.
PLINK_PREFIX="${PLINK_PREFIX:-/path/to/plink/genotypes_prefix}"

# Optional covariates for PEER. Format: samples x covariates, header present,
# sample IDs in column 1. If empty, PEER is run without observed covariates.
PEER_COVARIATES="${PEER_COVARIATES:-}"
N_PEER="${N_PEER:-15}"
MAX_PEER_ITER="${MAX_PEER_ITER:-1000}"

# Covariates for tensorQTL nominal/permutation scans. Format: covariates x
# samples, first column ID, matching tensorQTL conventions.
TENSORQTL_COVARIATES="${TENSORQTL_COVARIATES:-covar3_transpose.txt}"

# Historical filters used in the paper commands.
MIN_TPM="${MIN_TPM:-0.1}"
MIN_TPM_FRACTION="${MIN_TPM_FRACTION:-0.2}"
MIN_COUNT="${MIN_COUNT:-6}"
MIN_COUNT_FRACTION="${MIN_COUNT_FRACTION:-0.2}"

cd "$PROJECT_DIR"

require_file() {
  local file="$1"
  [[ -n "$file" && -e "$file" ]] || {
    echo "ERROR: missing file: $file" >&2
    exit 1
  }
}

subset_samples() {
  local matrix="$1"
  local keep_file="$2"
  local output="$3"

  if [[ -z "$keep_file" ]]; then
    cp "$matrix" "$output"
    return
  fi

  require_file "$keep_file"
  perl -Mstrict -Mwarnings -e '
    my ($matrix, $keep_file) = @ARGV;
    open my $kf, "<", $keep_file or die "Cannot open $keep_file: $!\n";
    my %keep = map { chomp; $_ => 1 } <$kf>;
    close $kf;

    open my $mf, "<", $matrix or die "Cannot open $matrix: $!\n";
    my $header = <$mf>;
    chomp $header;
    my @header = split /\t/, $header;
    my @idx = (0);
    for my $i (1..$#header) {
      push @idx, $i if exists $keep{$header[$i]};
    }
    print join("\t", @header[@idx]), "\n";
    while (my $line = <$mf>) {
      chomp $line;
      my @fields = split /\t/, $line;
      print join("\t", @fields[@idx]), "\n";
    }
  ' "$matrix" > "$output"
}

rename_header_samples() {
  local matrix="$1"
  local rename_file="$2"
  local output="$3"

  if [[ -z "$rename_file" ]]; then
    cp "$matrix" "$output"
    return
  fi

  require_file "$rename_file"
  perl -Mstrict -Mwarnings -e '
    my ($rename_file, $matrix) = @ARGV;
    open my $rf, "<", $rename_file or die "Cannot open $rename_file: $!\n";
    my @names = map { chomp; $_ } <$rf>;
    close $rf;

    open my $mf, "<", $matrix or die "Cannot open $matrix: $!\n";
    my $header = <$mf>;
    chomp $header;
    my @header = split /\t/, $header;
    die "Rename list has ", scalar(@names), " samples but matrix has ", scalar(@header) - 1, "\n"
      if @names != @header - 1;
    print join("\t", $header[0], @names), "\n";
    while (my $line = <$mf>) {
      print $line;
    }
  ' "$rename_file" "$matrix" > "$output"
}

prepare_expression_bed() {
  require_file "$TPM_MATRIX"
  require_file "$COUNT_MATRIX"
  require_file "$GTF_FILE"

  echo "[1/8] Filtering lowly expressed genes"
  perl "$ROOT_DIR/filter_low_expressed_genes.pl" \
    "$TPM_MATRIX" "$COUNT_MATRIX" \
    "$MIN_TPM" "$MIN_TPM_FRACTION" "$MIN_COUNT" "$MIN_COUNT_FRACTION" 2 \
    > "${PREFIX}.filtered.count.txt"

  echo "[2/8] TMM normalization"
  Rscript "$ROOT_DIR/normalize_counts_tmm.R" "${PREFIX}.filtered.count.txt" "${PREFIX}.TMM.txt"

  echo "[3/8] GTF to BED"
  perl "$ROOT_DIR/gtf_to_gene_bed.pl" "$GTF_FILE" gene > "${PREFIX}.genes.bed"
  sed -i -e 's/^chr//' "${PREFIX}.genes.bed"

  echo "[4/8] Sample subset/rename"
  subset_samples "${PREFIX}.TMM.txt" "$SAMPLE_KEEP" "${PREFIX}.sampled.txt"
  rename_header_samples "${PREFIX}.sampled.txt" "$SAMPLE_RENAME" "${PREFIX}.matrix.txt"

  echo "[5/8] Rank inverse normal transformation"
  Rscript "$ROOT_DIR/inverse_normal_transform.R" "${PREFIX}.matrix.txt" "${PREFIX}.INT.txt"

  echo "[6/8] Matrix to phenotype BED"
  perl "$ROOT_DIR/expression_matrix_to_phenotype_bed.pl" "${PREFIX}.genes.bed" "${PREFIX}.INT.txt" \
    > "${PREFIX}.INT.bed"

  echo "[7/8] Sort/index phenotype BED"
  sort -k1,1 -k2,2n "${PREFIX}.INT.bed" > "${PREFIX}.INT.sort.bed"
  bgzip -f "${PREFIX}.INT.sort.bed"
  tabix -f -p bed "${PREFIX}.INT.sort.bed.gz"
}

run_peer_and_residual_bed() {
  echo "[8/8] PEER factor estimation"
  if [[ -n "$PEER_COVARIATES" ]]; then
    require_file "$PEER_COVARIATES"
    Rscript "$ROOT_DIR/estimate_peer_factors.R" \
      --max_iter "$MAX_PEER_ITER" \
      --covariates "$PEER_COVARIATES" \
      "${PREFIX}.INT.sort.bed.gz" "${PREFIX}.PEER" "$N_PEER"
  else
    Rscript "$ROOT_DIR/estimate_peer_factors.R" \
      --max_iter "$MAX_PEER_ITER" \
      "${PREFIX}.INT.sort.bed.gz" "${PREFIX}.PEER" "$N_PEER"
  fi

  echo "[9/9] Convert PEER residuals to indexed phenotype BED"
  Rscript "$ROOT_DIR/inverse_normal_transform.R" "${PREFIX}.PEER.PEER_residuals.txt" "${PREFIX}.PEER_residuals.INT.txt"
  perl "$ROOT_DIR/expression_matrix_to_phenotype_bed.pl" "${PREFIX}.genes.bed" "${PREFIX}.PEER_residuals.INT.txt" \
    > "${PREFIX}.PEER_residuals.INT.bed"
  sort -k1,1 -k2,2n "${PREFIX}.PEER_residuals.INT.bed" > "${PREFIX}.PEER_residuals.INT.sort.bed"
  bgzip -f "${PREFIX}.PEER_residuals.INT.sort.bed"
  tabix -f -p bed "${PREFIX}.PEER_residuals.INT.sort.bed.gz"
}

run_tensorqtl_examples() {
  require_file "${PLINK_PREFIX}.bed"
  require_file "${PREFIX}.PEER_residuals.INT.sort.bed.gz"

  cat <<EOF

Run tensorQTL with one of your project-specific wrappers, for example:

  python /path/to/tensorQTL_CMC.py \\
    "$PLINK_PREFIX" \\
    "${PREFIX}.PEER_residuals.INT.sort.bed.gz" \\
    "$PREFIX"

  python /path/to/tensorQTL_CMC_withcovar.py \\
    "$PLINK_PREFIX" \\
    "${PREFIX}.PEER_residuals.INT.sort.bed.gz" \\
    "$TENSORQTL_COVARIATES" \\
    "$PREFIX"

For the cell-type interaction workflow in this repository, see:
  Cell_type_interaction_eQTL/TensorQTL_cisQTL.py
EOF
}

main() {
  local step="${1:-all}"
  case "$step" in
    expression)
      prepare_expression_bed
      ;;
    peer)
      run_peer_and_residual_bed
      ;;
    tensorqtl)
      run_tensorqtl_examples
      ;;
    all)
      prepare_expression_bed
      run_peer_and_residual_bed
      run_tensorqtl_examples
      ;;
    *)
      echo "Usage: $0 [expression|peer|tensorqtl|all]" >&2
      exit 1
      ;;
  esac
}

main "$@"
