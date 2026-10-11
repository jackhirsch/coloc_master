#!/bin/bash
# Run from the folder containing prep_sumstats.sh, AFTER setup_coloc_repo.sh
set -euo pipefail

PREP_SRC="$(pwd)/prep_sumstats.sh"
REPO="$HOME/coloc_master"
cd "$REPO"

# -----------------------------------------------------------------------------
# Commit 1: rsID reference builder + GWAS prep script
# -----------------------------------------------------------------------------
cat > scripts/00_build_rsid_reference.sh <<'EOF'
#!/bin/bash
#SBATCH --job-name=build_rsid_ref
#SBATCH --time=02:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=1
# Combine per-chromosome dbSNP files (chr:pos:ref:alt <TAB> rsID, GRCh38)
# into one file sorted for `join` in 01_preprocess_eqtl.sh.
set -euo pipefail

IN_DIR="/path/to/rsids_38"                  # build157_allmultiallelic_chr*.txt
OUT="/path/to/hg38_b157_rsid_sorted.txt"

cat "$IN_DIR"/build157_allmultiallelic_chr{1..23}.txt \
  | LC_ALL=C sort -t $'\t' -k1,1 -S 12G \
  > "$OUT"
echo "Wrote $(wc -l < "$OUT") variants to $OUT"
EOF

cp "$PREP_SRC" scripts/02_prep_gwas_sumstats.sh
sed -i '1i #!/bin/bash' scripts/02_prep_gwas_sumstats.sh
sed -i 's|^    "nightmare_MAF_sumstat"$|    "/path/to/your_sumstats.txt"|' scripts/02_prep_gwas_sumstats.sh
sed -i '/^#SORTED_LOOKUP=/d' scripts/02_prep_gwas_sumstats.sh
sed -i 's|^SORTED_LOOKUP=.*|SORTED_LOOKUP="${SORTED_LOOKUP:-/path/to/rsid_lookup_sorted_by_rsid.txt}"|' \
  scripts/02_prep_gwas_sumstats.sh
perl -0pi -e 's/echo "\[prep\] Next step:.*?coordinates\."\n/echo "[prep] Next step: set GWAS= in scripts\/submit_coloc.sh to a .coloc_ready.txt file above."\n/s' \
  scripts/02_prep_gwas_sumstats.sh
chmod +x scripts/*.sh

# README: point step 2 at the new script, add reference notes
perl -0pi -e 's/2\. \*\*GWAS prep\*\* – sumstats must contain/2. **`scripts\/02_prep_gwas_sumstats.sh`** – joins sumstats (`RSID A1 A2 BETA SE Z P N`, optional MAF\/EAF) to an rsID→chr\/pos lookup and writes `.coloc_ready.txt`. Output must contain/' README.md

cat >> README.md <<'EOF'

## Reference files

- **eQTL rsID reference (GRCh38)** – built with `scripts/00_build_rsid_reference.sh`
  from per-chromosome dbSNP b157 files (`chr:pos:ref:alt <TAB> rsID`, no header).
  Must be GRCh38 to match GTEx v8.
- **GWAS rsID lookup** – `CHR_POS CHR POS RSID`, no header, `LC_ALL=C` sorted by RSID.
  Its build must match the loci file.

eQTL variants are assigned to loci by rsID overlap with the GWAS SNPs in each
locus, so the eQTL build does not need to match the loci build.
EOF

# protect the new reference outputs from being committed
printf '\n# reference files\n*rsid*.txt\n' >> .gitignore

git add .
git status --short
git commit -q -m "Add rsID reference builder and GWAS sumstat prep"

# -----------------------------------------------------------------------------
# Commit 2: select eQTL variants per locus by rsID (build-agnostic)
# -----------------------------------------------------------------------------
perl -0pi -e 's/  eqtl_loc       <- eqtl\[eqtl\$chr_eqtl == loc_chr &\s+eqtl\$pos_eqtl >= loc_s\s+&\s+eqtl\$pos_eqtl <= loc_e, \]/  # match eQTL to locus via GWAS rsIDs (GTEx is hg38; loci may be hg19)\n  eqtl_loc       <- eqtl[eqtl\$rs_id %in% gwas_loc\$snp, ]/' \
  scripts/run_coloc_per_gene_per_locus.R

grep -n 'eqtl_loc       <-' scripts/run_coloc_per_gene_per_locus.R
git commit -q -am "Assign eQTL variants to loci by rsID instead of position"

git log --oneline
