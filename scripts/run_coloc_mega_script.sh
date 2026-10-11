#!/bin/bash
#SBATCH --job-name=coloc_gene_locus
#SBATCH --output=logs/coloc_%j.out
#SBATCH --error=logs/coloc_%j.err
#SBATCH --time=04:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=1
#SBATCH --partition=day

# =============================================================================
# PATHS  –  edit these
# =============================================================================

GWAS="nightmare_MAF_sumstat.coloc_ready.txt"
EQTL="Whole_Blood.v8.signif_variant_gene_pairs_annotated.txt.gz"
LOCI="/gpfs/gibbs/pi/gelernter/jh3583/LAVA/test.loci"
LD_PREFIX="/gpfs/gibbs/pi/gelernter/jh3583/LAVA/g1000_test"
OUTDIR="coloc_results"
SCRIPT="coloc_mega_script_v2.sh"

# GTEx v8 Brain Cerebellar Hemisphere N = 175  (Amygdala = 163)
EQTL_N=175

# =============================================================================
# GENE FILTER  –  choose ONE of the two options below, comment out the other
# =============================================================================

# Option A: run on ALL genes in the eQTL file (no filtering)
#GENE_LIST="ALL"

# Option B: supply a file with one gene name per line (e.g. genes.txt)
GENE_LIST="genes.txt"

# =============================================================================
# RUN  –  no need to edit below this line
# =============================================================================

mkdir -p logs "${OUTDIR}"

module load R/4.3.1   # adjust to your cluster's module name if needed

# Build the --gene_list argument dynamically
if [[ "${GENE_LIST}" == "ALL" ]]; then
    GENE_ARG="--gene_list ALL"
else
    GENE_ARG="--gene_list ${GENE_LIST}"
fi

Rscript "${SCRIPT}"   \
    --gwas      "${GWAS}"        \
    --eqtl      "${EQTL}"        \
    --loci      "${LOCI}"        \
    --ld_prefix "${LD_PREFIX}"   \
    --outdir    "${OUTDIR}"      \
    --eqtl_N    "${EQTL_N}"     \
    --gwas_type quant            \
    --window    0                \
    --p1        1e-4             \
    --p2        1e-4             \
    --p12       5e-6             \
    --min_snps  50               \
    --h4_thresh 0.8              \
    ${GENE_ARG}
