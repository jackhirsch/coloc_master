#!/bin/bash
#SBATCH --job-name=eqtl_array
#SBATCH --output=eqtl_array_%A_%a.out
#SBATCH --error=eqtl_array_%A_%a.err
#SBATCH --time=04:00:00
#SBATCH --mem=12G
#SBATCH --cpus-per-task=1
#SBATCH --array=1-15

set -euo pipefail

# =========================================
# USER INPUTS
# =========================================
REF_FILE="/gpfs/gibbs/pi/gelernter/beckse/for_Jack/hg38_snps_sortbychromposrefalt_LC_ALL.txt"

OUTDIR="annotated_eqtl_outputs"
TMP_DIR="tmp_eqtl_array"

mkdir -p "$OUTDIR" "$TMP_DIR"

# =========================================
# BUILD FILE LIST
# =========================================
mapfile -t FILES < <(ls *.signif_variant_gene_pairs.txt.gz | sort)

N=${#FILES[@]}
echo "Found $N files"

if [[ -z "${SLURM_ARRAY_TASK_ID:-}" ]]; then
    echo "ERROR: submit with sbatch"
    exit 1
fi

INPUT_FILE="${FILES[$SLURM_ARRAY_TASK_ID-1]}"
echo "Processing: $INPUT_FILE"

BASE=$(basename "$INPUT_FILE" .txt.gz)

TMP_IDS="$TMP_DIR/${BASE}_ids.txt"
TMP_EQTL="$TMP_DIR/${BASE}_eqtl.tsv"
TMP_EQTL_BODY="$TMP_DIR/${BASE}_eqtl_body.tsv"

TMP_SORT_EQTL="$TMP_DIR/${BASE}_eqtl_sorted.tsv"
TMP_SORT_REF="$TMP_DIR/${BASE}_ref_sorted.tsv"

TMP_JOIN="$TMP_DIR/${BASE}_joined.tsv"

TMP_R="$TMP_DIR/${BASE}_map.R"

OUTPUT_FILE="$OUTDIR/${BASE}_annotated.txt"

# =========================================
# STEP 1: EXTRACT GENE IDS
# =========================================
zcat "$INPUT_FILE" \
| awk 'NR>1 {print $2}' \
| sed 's/\..*//' \
| sort -u > "$TMP_IDS"

export INPUT_FILE TMP_IDS TMP_EQTL_BODY

# =========================================
# STEP 2: R PROCESSING
# =========================================
cat << 'EOF' > "$TMP_R"
library(data.table)
library(httr)
library(jsonlite)

input_file <- Sys.getenv("INPUT_FILE")
ids_file <- Sys.getenv("TMP_IDS")
tmp_eqtl <- Sys.getenv("TMP_EQTL_BODY")

ids <- readLines(ids_file)

chunks <- split(ids, ceiling(seq_along(ids)/500))

mapping <- c()

for (i in seq_along(chunks)) {

  cat("chunk", i, "\n")

  r <- POST(
    "https://biotools.fr/human/ensembl_symbol_converter/",
    body=list(
      api=1,
      ids=toJSON(chunks[[i]])
    )
  )

  res <- fromJSON(content(r, "text"), simplifyVector=TRUE)

  mapping <- c(mapping, res)
}

map_dt <- data.table(
  gene_id_clean = names(mapping),
  gene = as.character(unname(mapping))
)

dt <- fread(input_file)

# clean gene IDs
dt[, gene_id_clean := sub("\\..*", "", gene_id)]

# variant conversion
# chr1_845402_A_G_b38 -> 1:845402:A:G
dt[, variant_id := gsub("_b38", "", variant_id)]
dt[, variant_id := gsub("_", ":", variant_id)]
dt[, variant_id := sub("^chr", "", variant_id)]

# merge gene symbols
dt <- merge(dt, map_dt, by="gene_id_clean", all.x=TRUE)

# explicit column order
desired_order <- c(
  "variant_id",
  "gene_id",
  "gene",
  "gene_id_clean",
  "tss_distance",
  "ma_samples",
  "ma_count",
  "maf",
  "pval_nominal",
  "slope",
  "slope_se",
  "pval_nominal_threshold",
  "min_pval_nominal",
  "pval_beta"
)

desired_order <- desired_order[desired_order %in% names(dt)]

setcolorder(dt, desired_order)

fwrite(
  dt,
  file=tmp_eqtl,
  sep="\t",
  quote=FALSE,
  col.names=FALSE
)
EOF

Rscript "$TMP_R"

# =========================================
# STEP 3: SORT FILES BEFORE JOIN
# =========================================
LC_ALL=C sort -t $'\t' -k1,1 "$TMP_EQTL_BODY" > "$TMP_SORT_EQTL"

LC_ALL=C sort -t $'\t' -k1,1 "$REF_FILE" > "$TMP_SORT_REF"

# =========================================
# STEP 4: JOIN rsIDs
# =========================================
LC_ALL=C join \
  -t $'\t' \
  -1 1 \
  -2 1 \
  "$TMP_SORT_EQTL" \
  "$TMP_SORT_REF" \
  > "$TMP_JOIN"

# =========================================
# STEP 5: ADD HEADER
# =========================================
{
echo -e "variant_id\tgene_id\tgene\tgene_id_clean\ttss_distance\tma_samples\tma_count\tmaf\tpval_nominal\tslope\tslope_se\tpval_nominal_threshold\tmin_pval_nominal\tpval_beta\trs_id"

cat "$TMP_JOIN"

} > "$OUTPUT_FILE"

# =========================================
# STEP 6: COMPRESS
# =========================================
gzip -f "$OUTPUT_FILE"

# =========================================
# STEP 7: CLEANUP
# =========================================
rm -f \
"$TMP_IDS" \
"$TMP_EQTL_BODY" \
"$TMP_SORT_EQTL" \
"$TMP_SORT_REF" \
"$TMP_JOIN" \
"$TMP_R"

echo "Done $INPUT_FILE"
