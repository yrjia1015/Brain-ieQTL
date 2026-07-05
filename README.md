# Brain lncRNA-ieQTL Analysis Code

This folder contains the code used for expression preprocessing, cell-type interaction eQTL mapping, COLOC, SMR, enrichment, and mashr analyses.

## Data Access

ieQTL summary statistics and the transcriptome model are available from the
[Brain-ieQTL Atlas](https://brain.hitxqtl.org.cn/Brain-ieQTLAtlas/#/).

## Folder Structure

```text
Brain-ieQTL/
|-- Expression_eQTL_preprocessing/   # gene expression matrix -> tensorQTL phenotype BED
|-- Cell_type_interaction_eQTL/      # cell-type interaction eQTL
|-- COLOC/                           # GWAS-eQTL colocalization
|-- SMR/                             # SMR analysis
|-- eQTL_analysis/                   # enrichment and mashr analysis
`-- README.md
```

## 1. Expression Matrix and eQTL Preprocessing

Run from your analysis working directory:

```bash
PROJECT_DIR="$PWD" \
PREFIX="analysis_name" \
TPM_MATRIX="/path/to/gene_tpm_matrix.txt" \
COUNT_MATRIX="/path/to/gene_count_matrix.txt" \
GTF_FILE="/path/to/gene_annotation.gtf" \
PLINK_PREFIX="/path/to/plink_prefix" \
N_PEER=15 \
bash /path/to/Brain-ieQTL/Expression_eQTL_preprocessing/run_expression_eqtl_pipeline.sh all
```

Fill in:

```text
PROJECT_DIR     output/work directory
PREFIX          output prefix, for example newlncRNA_exp
TPM_MATRIX      gene TPM matrix, first column is gene ID
COUNT_MATRIX    gene raw count/fragments matrix, first column is gene ID
GTF_FILE        gene annotation file
PLINK_PREFIX    genotype PLINK prefix, without .bed/.bim/.fam
N_PEER          number of PEER factors, usually 15/30/45/60 by sample size
```

Optional:

```text
SAMPLE_KEEP        one sample ID per line, used to keep selected samples
SAMPLE_RENAME      one replacement sample ID per line, same order as kept samples
PEER_COVARIATES    observed covariates for PEER, samples x covariates
MAX_PEER_ITER      PEER max iterations
```

Run only one part if needed:

```bash
bash /path/to/Brain-ieQTL/Expression_eQTL_preprocessing/run_expression_eqtl_pipeline.sh expression
bash /path/to/Brain-ieQTL/Expression_eQTL_preprocessing/run_expression_eqtl_pipeline.sh peer
bash /path/to/Brain-ieQTL/Expression_eQTL_preprocessing/run_expression_eqtl_pipeline.sh tensorqtl
```

Main outputs:

```text
PREFIX.filtered.count.txt
PREFIX.TMM.txt
PREFIX.INT.sort.bed.gz
PREFIX.PEER.PEER_covariates.txt
PREFIX.PEER.PEER_residuals.txt
PREFIX.PEER_residuals.INT.sort.bed.gz
```

Use the final indexed phenotype BED and covariate file for tensorQTL.

## 2. Cell-Type Interaction eQTL

Edit the values in the cell-type interaction script:

```text
PLINK prefix
phenotype BED
cell proportion matrix
covariate matrix
cohort name
cell type, such as Ex/Oli/In/Ast
```

Then run:

```bash
python /path/to/Brain-ieQTL/Cell_type_interaction_eQTL/TensorQTL_cisQTL.py
```

## 3. COLOC

Run:

```bash
Rscript /path/to/Brain-ieQTL/COLOC/COLOC_analysis.R \
  --trait TRAIT_NAME \
  --cell CELL_NAME \
  --trait_type cc \
  --gwas_file /path/to/gwas_summary.txt \
  --eqtl_pattern "/path/to/eqtl/{cell}/{cell}_chr{chr}.txt" \
  --out_dir /path/to/coloc_output \
  --cases CASE_NUMBER \
  --controls CONTROL_NUMBER \
  --gwas_N GWAS_SAMPLE_SIZE \
  --eqtl_N EQTL_SAMPLE_SIZE \
  --n_cores 8
```

Fill in:

```text
trait/cell       names used for output
trait_type       cc for case-control, quant for quantitative trait
gwas_file        GWAS summary statistics
eqtl_pattern     eQTL file pattern; keep {cell} and {chr}
cases/controls   required for case-control traits
gwas_N/eqtl_N    sample sizes
```

## 4. SMR

Run:

```bash
sbatch /path/to/Brain-ieQTL/SMR/SMR_analysis.bash -- \
  --smr_bin /path/to/smr \
  --bfile_prefix "/path/to/ref_chr{chr}" \
  --eqtl_dir /path/to/eqtl_besd_dir \
  --gwas /path/to/gwas_summary \
  --out_dir /path/to/smr_output \
  --threads 20
```

Fill in:

```text
smr_bin        SMR executable
bfile_prefix   reference genotype prefix, keep {chr}
eqtl_dir       eQTL BESD directory
gwas           GWAS summary file
out_dir        output directory
```

## 5. Enrichment and mashr

Peak enrichment:

```bash
Rscript /path/to/Brain-ieQTL/eQTL_analysis/DeltaEnrich.R \
  --ieqtl_file /path/to/ieqtl_result.txt \
  --anno_file /path/to/snp_annotation.txt \
  --cell_type CELL_NAME \
  --gene_type GENE_TYPE \
  --direction Positive \
  --out_null /path/to/null_output.txt \
  --out_log /path/to/summary_output.txt
```

mashr:

```bash
Rscript /path/to/Brain-ieQTL/eQTL_analysis/MashR_analysis.R \
  --random_set /path/to/random_set.txt \
  --strong_set /path/to/strong_set.txt \
  --out_rds /path/to/output.rds
```

## Notes

Make sure chromosome names match between phenotype BED and genotype files. If genotype chromosomes are `1-22`, remove `chr` from annotation-derived BED files.

Sample IDs must match exactly across expression, genotype, PEER covariates, and tensorQTL covariates.

The old scattered root scripts have been moved into the expression preprocessing folder and renamed by function.
