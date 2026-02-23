#!/bin/bash
#SBATCH -p q07
#SBATCH -J SMR
#SBATCH -c 20
#SBATCH --mem 10G
#SBATCH -o smr.%A.%a.log
#SBATCH -a 1-1

set -euo pipefail

# ============================================================
# Generic SMR runner
#
# Fixed cell types: Ast, Ex, In, Oli
# GWAS is provided from command line: --gwas <path>
#
# Required arguments:
#   --smr_bin        path to smr executable
#   --bfile_prefix   plink bfile prefix template with {chr}
#                   e.g. "/path/ref/hrs_1kg_hwe1e-6_chr{chr}"
#   --eqtl_dir       directory containing cell subfolders and chr prefixes
#                   e.g. "/path/BESD/brainmeta_lncRNA"
#                   expects: <eqtl_dir>/<cell>/<cell>_chr{chr}
#   --gwas           gwas summary path/prefix (passed to --gwas-summary)
#   --out_dir        output root directory
#
# Optional:
#   --threads        threads for SMR (default: SLURM_CPUS_PER_TASK or 1)
#   --chr_start      default 1
#   --chr_end        default 22
#
# Example:
# sbatch -a 1-1 smr_by_gwas.sh -- \
#   --smr_bin /path/smr-1.3.1 \
#   --bfile_prefix "/path/ref/hrs_1kg_hwe1e-6_chr{chr}" \
#   --eqtl_dir "/path/BESD/brainmeta_lncRNA" \
#   --gwas "/path/GWAS/ADHD" \
#   --out_dir "/path/SMR/Result/Second" \
#   --threads 20
# ============================================================

CELL_TYPES=("Ast" "Ex" "In" "Oli")

die(){ echo "ERROR: $*" >&2; exit 1; }

# ----------------------------
# Minimal CLI parsing
# ----------------------------
SMR_BIN=""
BFILE_PREFIX=""
EQTL_DIR=""
GWAS_FILE=""
OUT_DIR=""
THREADS="${SLURM_CPUS_PER_TASK:-1}"
CHR_START=1
CHR_END=22

# Allow passing args after "--"
while [[ $# -gt 0 ]]; do
  case "$1" in
    --smr_bin)      SMR_BIN="$2"; shift 2;;
    --bfile_prefix) BFILE_PREFIX="$2"; shift 2;;
    --eqtl_dir)     EQTL_DIR="$2"; shift 2;;
    --gwas)         GWAS_FILE="$2"; shift 2;;
    --out_dir)      OUT_DIR="$2"; shift 2;;
    --threads)      THREADS="$2"; shift 2;;
    --chr_start)    CHR_START="$2"; shift 2;;
    --chr_end)      CHR_END="$2"; shift 2;;
    --) shift; break;;
    *) die "Unknown argument: $1";;
  esac
done

# ----------------------------
# Validate required args
# ----------------------------
[[ -n "$SMR_BIN" ]]      || die "Missing --smr_bin"
[[ -n "$BFILE_PREFIX" ]] || die "Missing --bfile_prefix"
[[ -n "$EQTL_DIR" ]]     || die "Missing --eqtl_dir"
[[ -n "$GWAS_FILE" ]]    || die "Missing --gwas"
[[ -n "$OUT_DIR" ]]      || die "Missing --out_dir"

[[ -x "$SMR_BIN" ]] || die "SMR binary not executable: $SMR_BIN"
[[ -d "$EQTL_DIR" ]] || die "eqtl_dir not found: $EQTL_DIR"
[[ -e "$GWAS_FILE" ]] || die "GWAS not found: $GWAS_FILE"

# Name outputs by GWAS basename
GWAS_NAME="$(basename "$GWAS_FILE")"

echo "[INFO] GWAS=$GWAS_NAME"
echo "[INFO] Threads=$THREADS chr=${CHR_START}-${CHR_END}"
echo "[INFO] EQTL_DIR=$EQTL_DIR"
echo "[INFO] OUT_DIR=$OUT_DIR"
echo "[INFO] BFILE_PREFIX template=$BFILE_PREFIX"

mkdir -p "$OUT_DIR"

# ----------------------------
# Run SMR
# ----------------------------
for cell in "${CELL_TYPES[@]}"; do
  out_subdir="${OUT_DIR}/${cell}/${GWAS_NAME}"
  mkdir -p "$out_subdir"

  for chr in $(seq "$CHR_START" "$CHR_END"); do
    bfile="${BFILE_PREFIX/\{chr\}/$chr}"
    eqtl_prefix="${EQTL_DIR}/${cell}/${cell}_chr${chr}"
    out_prefix="${out_subdir}/${cell}_chr${chr}"

    [[ -e "${bfile}.bed" ]] || die "Missing PLINK bed: ${bfile}.bed"
    # eQTL prefix may correspond to multiple files; check common BESD extensions loosely
    if [[ ! -e "${eqtl_prefix}.besd" && ! -e "${eqtl_prefix}" ]]; then
      echo "[WARN] eQTL prefix not found for cell=$cell chr=$chr: $eqtl_prefix (still running; SMR may fail)"
    fi

    echo "[INFO] Running cell=$cell chr=$chr"
    "$SMR_BIN" \
      --bfile "$bfile" \
      --gwas-summary "$GWAS_FILE" \
      --beqtl-summary "$eqtl_prefix" \
      --out "$out_prefix" \
      --thread-num "$THREADS"
  done
done

echo "[INFO] Done. GWAS=$GWAS_NAME"