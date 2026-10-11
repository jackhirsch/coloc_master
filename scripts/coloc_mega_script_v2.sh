#!/usr/bin/env Rscript
# =============================================================================
# run_coloc_per_gene_per_locus.R
#
# Colocalization between nightmare GWAS and GTEx brain eQTL sumstats,
# run per gene per locus.
#
# Outputs:
#   <outdir>/
#     per_pair/              locus{LOC}_{gene}.coloc.txt
#     coloc_summary_all.txt  all results sorted by PP.H4
#     coloc_locus_report.txt loci file + best gene + all PP values
#     plots/
#       pp4_dotplot.pdf
#       pp_barplot_by_locus.pdf
#       pp_heatmap.pdf
#       pp4_bin_summary.pdf
#       locus{LOC}_{gene}_scatter.pdf   (coloc pairs only)
#
# Usage:
#   Rscript run_coloc_per_gene_per_locus.R \
#       --gwas        nightmare_MAF_sumstat.coloc_ready.txt \
#       --eqtl        Brain_Amygdala.v8.signif_variant_gene_pairs_annotated.txt.gz \
#       --loci        test.loci \
#       --outdir      coloc_results \
#       --gene_list   ALL              # run all genes
#       --gene_list   genes.txt        # or filter to genes listed in this file
# =============================================================================

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(coloc)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(scales)
  library(ggrepel)
})

# ---------------------------------------------------------------------------
# 1.  CLI arguments
# ---------------------------------------------------------------------------
option_list <- list(
  make_option("--gwas",       type = "character", help = "Path to GWAS sumstat file"),
  make_option("--eqtl",       type = "character", help = "Path to eQTL file (can be .gz)"),
  make_option("--loci",       type = "character", help = "Path to loci file (LOC CHR START STOP)"),
  make_option("--ld_prefix",  type = "character", default = NULL,
              help = "Prefix for PLINK LD reference (bed/bim/fam)."),
  make_option("--outdir",     type = "character", default = "coloc_results",
              help = "Output directory [default: coloc_results]"),
  make_option("--gene_list",  type = "character", default = "ALL",
              help = paste0("Gene filter. Either:\n",
                            "  'ALL'         – run every gene in the eQTL file (default)\n",
                            "  path/genes.txt – plain text file, one gene name per line;\n",
                            "                   gene names must match column 3 (gene) of eQTL file")),
  make_option("--window",     type = "integer",   default = 0,
              help = "Extra bp window around each locus [default: 0]"),
  make_option("--p1",         type = "double",    default = 1e-4,  help = "coloc prior p1"),
  make_option("--p2",         type = "double",    default = 1e-4,  help = "coloc prior p2"),
  make_option("--p12",        type = "double",    default = 5e-6,  help = "coloc prior p12"),
  make_option("--gwas_type",  type = "character", default = "quant",
              help = "GWAS trait type: 'quant' or 'cc' [default: quant]"),
  make_option("--gwas_N",     type = "integer",   default = NULL,
              help = "Override GWAS N (uses column N in sumstats if not set)"),
  make_option("--eqtl_N",     type = "integer",   default = 163,
              help = "eQTL sample size [default: 163]"),
  make_option("--min_snps",   type = "integer",   default = 50,
              help = "Min overlapping SNPs required to run coloc [default: 50]"),
  make_option("--h4_thresh",  type = "double",    default = 0.8,
              help = "PP.H4 threshold to flag as colocalising [default: 0.8]")
)

opt <- parse_args(OptionParser(option_list = option_list))
for (arg in c("gwas", "eqtl", "loci")) {
  if (is.null(opt[[arg]])) stop(sprintf("Required argument --%s is missing.", arg))
}

dir_main  <- opt$outdir
dir_pairs <- file.path(dir_main, "per_pair")
dir_plots <- file.path(dir_main, "plots")
for (d in c(dir_main, dir_pairs, dir_plots))
  dir.create(d, showWarnings = FALSE, recursive = TRUE)

# ---------------------------------------------------------------------------
# 2.  Parse gene filter
# ---------------------------------------------------------------------------
if (toupper(opt$gene_list) == "ALL") {
  gene_filter <- NULL     # NULL means "keep everything"
  message("Gene filter   : ALL genes")
} else {
  if (!file.exists(opt$gene_list))
    stop(sprintf("Gene list file not found: %s", opt$gene_list))
  gene_filter <- trimws(readLines(opt$gene_list))
  gene_filter <- gene_filter[nchar(gene_filter) > 0]   # drop blank lines
  message(sprintf("Gene filter   : %d genes loaded from %s",
                  length(gene_filter), opt$gene_list))
  message("  First 5: ", paste(head(gene_filter, 5), collapse = ", "))
}

message("=== coloc per gene per locus ===")
message("GWAS          : ", opt$gwas)
message("eQTL          : ", opt$eqtl)
message("Loci          : ", opt$loci)
message("Out dir       : ", opt$outdir)
message("Window        : ", opt$window, " bp")
message("eQTL N        : ", opt$eqtl_N)
message("Priors        : p1=", opt$p1, "  p2=", opt$p2, "  p12=", opt$p12)
message("H4 threshold  : ", opt$h4_thresh)

# ---------------------------------------------------------------------------
# Colour palette (colourblind-friendly)
# ---------------------------------------------------------------------------
PP_COLS <- c(
  H0 = "#4E4E4E",
  H1 = "#4292C6",
  H2 = "#41AB5D",
  H3 = "#FD8D3C",
  H4 = "#CB181D"
)

# ---------------------------------------------------------------------------
# 3.  Load data
# ---------------------------------------------------------------------------
message("\n[1/4] Loading GWAS sumstats...")
gwas <- fread(opt$gwas, data.table = FALSE)
colnames(gwas) <- tolower(colnames(gwas))
if (!"varbeta" %in% colnames(gwas)) {
  message("  varbeta not found – computing as se^2")
  gwas$varbeta <- gwas$se^2
}
gwas$chr <- as.integer(gwas$chr)
gwas$pos <- as.integer(gwas$pos)
message("  Loaded ", nrow(gwas), " GWAS SNPs")

message("\n[2/4] Loading eQTL data...")
eqtl <- fread(opt$eqtl, data.table = FALSE)
colnames(eqtl) <- tolower(colnames(eqtl))

# Parse chr/pos from variant_id (format chr:pos:REF:ALT)
vid_parts         <- strsplit(eqtl$variant_id, ":", fixed = TRUE)
eqtl$chr_eqtl     <- as.integer(sapply(vid_parts, `[`, 1))
eqtl$pos_eqtl     <- as.integer(sapply(vid_parts, `[`, 2))
eqtl$varbeta_eqtl <- eqtl$slope_se^2

message("  Loaded ", nrow(eqtl), " eQTL records | ",
        length(unique(eqtl$gene)), " unique genes before filtering")

# ── Apply gene filter ──────────────────────────────────────────────────────
if (!is.null(gene_filter)) {

  # Report any requested genes that are completely absent from the eQTL file
  present_in_eqtl <- unique(eqtl$gene)
  missing_genes   <- setdiff(gene_filter, present_in_eqtl)
  if (length(missing_genes) > 0) {
    message(sprintf("  WARNING: %d requested gene(s) not found in eQTL file:",
                    length(missing_genes)))
    message("    ", paste(missing_genes, collapse = ", "))
  }

  eqtl <- eqtl[eqtl$gene %in% gene_filter, ]
  message(sprintf("  After gene filter: %d eQTL records | %d unique genes",
                  nrow(eqtl), length(unique(eqtl$gene))))

  if (nrow(eqtl) == 0)
    stop("No eQTL records remain after gene filtering. Check gene names match column 3 of the eQTL file.")
}

message("\n[3/4] Loading loci...")
loci <- fread(opt$loci, data.table = FALSE)
colnames(loci) <- toupper(colnames(loci))
loci$CHR   <- as.integer(loci$CHR)
loci$START <- as.integer(loci$START)
loci$STOP  <- as.integer(loci$STOP)
message("  Loaded ", nrow(loci), " loci")

# ---------------------------------------------------------------------------
# 4.  coloc helper
# ---------------------------------------------------------------------------
run_coloc_one <- function(gwas_sub, eqtl_sub, gene_name, locus_id,
                          gwas_type, gwas_N_override, eqtl_N, priors, min_snps) {

  EMPTY <- data.frame(locus = locus_id, gene = gene_name, nsnps = 0L,
                      PP.H0 = NA_real_, PP.H1 = NA_real_,
                      PP.H2 = NA_real_, PP.H3 = NA_real_, PP.H4 = NA_real_,
                      status = NA_character_, stringsAsFactors = FALSE)

  merged <- merge(
    gwas_sub[, c("snp","maf","beta","varbeta","n","pvalue")],
    eqtl_sub[, c("rs_id","maf","slope","varbeta_eqtl","pval_nominal")],
    by.x = "snp", by.y = "rs_id", suffixes = c("_gwas","_eqtl")
  )

  if (nrow(merged) < min_snps) {
    EMPTY$nsnps  <- nrow(merged)
    EMPTY$status <- paste0("SKIPPED: only ", nrow(merged), " overlapping SNPs")
    return(list(result = EMPTY, merged = merged))
  }
  merged <- merged[!duplicated(merged$snp), ]

  N_gwas <- if (!is.null(gwas_N_override)) gwas_N_override else
               round(mean(merged$n, na.rm = TRUE))

  D1 <- list(snp = merged$snp, beta = merged$beta, varbeta = merged$varbeta,
             MAF = merged$maf_gwas, N = N_gwas, type = gwas_type)
  D2 <- list(snp = merged$snp, beta = merged$slope, varbeta = merged$varbeta_eqtl,
             MAF = merged$maf_eqtl, N = eqtl_N, type = "quant")

  chk1 <- tryCatch(check_dataset(D1, warn.minp = FALSE), error = function(e) e)
  chk2 <- tryCatch(check_dataset(D2, warn.minp = FALSE), error = function(e) e)
  if (inherits(chk1, "error") || inherits(chk2, "error")) {
    EMPTY$nsnps  <- nrow(merged)
    EMPTY$status <- paste("Dataset check failed:",
                          if (inherits(chk1,"error")) chk1$message else chk2$message)
    return(list(result = EMPTY, merged = merged))
  }

  res <- tryCatch(
    coloc.abf(D1, D2, p1 = priors$p1, p2 = priors$p2, p12 = priors$p12),
    error = function(e) e
  )
  if (inherits(res, "error")) {
    EMPTY$nsnps  <- nrow(merged)
    EMPTY$status <- paste("coloc error:", res$message)
    return(list(result = EMPTY, merged = merged))
  }

  pp <- res$summary
  result <- data.frame(
    locus  = locus_id, gene = gene_name,
    nsnps  = as.integer(pp["nsnps"]),
    PP.H0  = as.numeric(pp["PP.H0.abf"]),
    PP.H1  = as.numeric(pp["PP.H1.abf"]),
    PP.H2  = as.numeric(pp["PP.H2.abf"]),
    PP.H3  = as.numeric(pp["PP.H3.abf"]),
    PP.H4  = as.numeric(pp["PP.H4.abf"]),
    status = "OK", stringsAsFactors = FALSE, row.names = NULL
  )
  list(result = result, merged = merged)
}

# ---------------------------------------------------------------------------
# 5.  Main loop
# ---------------------------------------------------------------------------
message("\n[4/4] Running coloc...\n")

all_results  <- list()
scatter_data <- list()
priors       <- list(p1 = opt$p1, p2 = opt$p2, p12 = opt$p12)

for (i in seq_len(nrow(loci))) {

  loc_id  <- loci$LOC[i]
  loc_chr <- loci$CHR[i]
  loc_s   <- loci$START[i] - opt$window
  loc_e   <- loci$STOP[i]  + opt$window

  message(sprintf("Locus %-6s  chr%s:%s-%s", loc_id, loc_chr,
                  format(loc_s, big.mark=","), format(loc_e, big.mark=",")))

  gwas_loc <- gwas[gwas$chr == loc_chr & gwas$pos >= loc_s & gwas$pos <= loc_e, ]
  if (nrow(gwas_loc) == 0) { message("  No GWAS SNPs – skipping"); next }

  eqtl_loc       <- eqtl[eqtl$chr_eqtl == loc_chr &
                          eqtl$pos_eqtl >= loc_s   &
                          eqtl$pos_eqtl <= loc_e, ]
  genes_in_locus <- unique(eqtl_loc$gene)

  # Drop NA/NULL gene names that sometimes appear in annotation files
  genes_in_locus <- genes_in_locus[!is.na(genes_in_locus) &
                                   genes_in_locus != "NULL" &
                                   genes_in_locus != "NA"  &
                                   nchar(genes_in_locus) > 0]

  if (length(genes_in_locus) == 0) { message("  No eQTL genes – skipping"); next }

  message(sprintf("  %d GWAS SNPs | %d eQTL gene(s): %s",
                  nrow(gwas_loc), length(genes_in_locus),
                  paste(genes_in_locus, collapse = ", ")))

  for (gene in genes_in_locus) {

    eqtl_gene <- eqtl_loc[eqtl_loc$gene == gene, ]
    eqtl_gene <- eqtl_gene[order(eqtl_gene$pval_nominal), ]
    eqtl_gene <- eqtl_gene[!duplicated(eqtl_gene$rs_id), ]

    out    <- run_coloc_one(gwas_loc, eqtl_gene, gene, loc_id,
                            opt$gwas_type, opt$gwas_N, opt$eqtl_N,
                            priors, opt$min_snps)
    result <- out$result

    pp4_str <- if (!is.na(result$PP.H4)) sprintf("%.3f", result$PP.H4) else "NA"
    message(sprintf("    %-20s  nSNPs=%-5d  PP.H4=%s  [%s]",
                    gene, result$nsnps, pp4_str, result$status))

    write.table(result,
                file.path(dir_pairs, sprintf("locus%s_%s.coloc.txt", loc_id, gene)),
                sep = "\t", quote = FALSE, row.names = FALSE)

    all_results[[length(all_results) + 1]] <- result

    if (result$status == "OK" && !is.null(out$merged) && nrow(out$merged) > 0)
      scatter_data[[length(scatter_data) + 1]] <- list(
        locus = loc_id, gene = gene, chr = loc_chr,
        PP.H4 = result$PP.H4, merged = out$merged
      )
  }
}

# ---------------------------------------------------------------------------
# 6.  Summary table
# ---------------------------------------------------------------------------
if (length(all_results) == 0) {
  message("\nNo coloc results generated. Exiting.")
  quit(save = "no", status = 0)
}

summary_df <- do.call(rbind, all_results)
rownames(summary_df) <- NULL
summary_df <- summary_df[order(-summary_df$PP.H4, na.last = TRUE), ]
pp_cols    <- c("PP.H0","PP.H1","PP.H2","PP.H3","PP.H4")
summary_df[pp_cols] <- lapply(summary_df[pp_cols],
                               function(x) round(as.numeric(x), 4))

write.table(summary_df,
            file.path(dir_main, "coloc_summary_all.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)
message("\nWrote coloc_summary_all.txt  (", nrow(summary_df), " rows)")

# ---------------------------------------------------------------------------
# 7.  Locus report
# ---------------------------------------------------------------------------
message("Building locus report...")

ok_df <- summary_df[summary_df$status == "OK", ]

if (nrow(ok_df) > 0) {

  best_df <- ok_df[!duplicated(ok_df$locus), ]
  colnames(best_df)[colnames(best_df) == "gene"]  <- "best_gene"
  colnames(best_df)[colnames(best_df) == "nsnps"] <- "best_gene_nsnps"
  colnames(best_df)[colnames(best_df) %in% pp_cols] <-
    paste0("best_", colnames(best_df)[colnames(best_df) %in% pp_cols])

  genes_per_locus <- aggregate(gene ~ locus, data = ok_df,
    FUN = function(x) paste(sort(unique(x)), collapse = "|"))
  colnames(genes_per_locus)[2] <- "all_tested_genes"

  ok_df$gene_pp4 <- sprintf("%s=%.3f", ok_df$gene, ok_df$PP.H4)
  pp4_per_locus  <- aggregate(gene_pp4 ~ locus, data = ok_df,
    FUN = function(x) paste(x, collapse = "; "))
  colnames(pp4_per_locus)[2] <- "all_gene_PP.H4"

  locus_report <- loci
  locus_report <- merge(locus_report, genes_per_locus,
                        by.x = "LOC", by.y = "locus", all.x = TRUE)
  locus_report <- merge(locus_report, pp4_per_locus,
                        by.x = "LOC", by.y = "locus", all.x = TRUE)
  locus_report <- merge(locus_report,
                        best_df[, c("locus","best_gene","best_gene_nsnps",
                                    paste0("best_", pp_cols), "status")],
                        by.x = "LOC", by.y = "locus", all.x = TRUE)
  locus_report$colocalises <- ifelse(
    !is.na(locus_report$best_PP.H4) & locus_report$best_PP.H4 >= opt$h4_thresh,
    "YES", "NO"
  )
  locus_report <- locus_report[match(loci$LOC, locus_report$LOC), ]

} else {
  locus_report <- loci
  locus_report[, c("all_tested_genes","all_gene_PP.H4","best_gene",
                   "best_gene_nsnps", paste0("best_", pp_cols),
                   "status","colocalises")] <- NA
}

# If a gene_list was used, record it in the report
locus_report$gene_filter_applied <- if (is.null(gene_filter)) "ALL" else
  paste(gene_filter, collapse = "|")

write.table(locus_report,
            file.path(dir_main, "coloc_locus_report.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
message("Wrote coloc_locus_report.txt")

# ---------------------------------------------------------------------------
# 8.  Plots
# ---------------------------------------------------------------------------
message("Generating plots...")

theme_coloc <- theme_bw(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "#F0F0F0", colour = "grey70"),
    plot.title       = element_text(face = "bold", size = 13),
    plot.subtitle    = element_text(size = 10, colour = "grey40"),
    legend.position  = "right"
  )

filter_label <- if (is.null(gene_filter)) "all genes" else
  sprintf("%d target genes", length(gene_filter))

if (nrow(ok_df) > 0) {

  # ── 8a. PP.H4 dot-plot ──────────────────────────────────────────────────
  ok_df$label   <- with(ok_df, ifelse(PP.H4 >= opt$h4_thresh,
                                      paste0(gene, "\n(L", locus, ")"), NA))
  ok_df$locus_f <- factor(paste0("L", ok_df$locus))

  p_dot <- ggplot(ok_df,
                  aes(x = locus_f, y = PP.H4,
                      colour = PP.H3, size = nsnps, label = label)) +
    geom_hline(yintercept = opt$h4_thresh, linetype = "dashed",
               colour = "grey50", linewidth = 0.6) +
    geom_point(alpha = 0.85) +
    geom_text_repel(size = 3, colour = "black", na.rm = TRUE,
                    min.segment.length = 0.2, max.overlaps = 20) +
    scale_colour_gradient2(low = "#4292C6", mid = "#FDD0A2", high = "#FD8D3C",
                           midpoint = 0.5, name = "PP.H3\n(distinct signal)") +
    scale_size_continuous(range = c(2, 8), name = "nSNPs") +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
    labs(
      title    = "Colocalization probability (PP.H4) per gene x locus",
      subtitle = sprintf("Nightmare GWAS x Brain eQTL | %s | dashed = H4 threshold (%.2f)",
                         filter_label, opt$h4_thresh),
      x = "Locus", y = "PP.H4 (shared causal variant)"
    ) +
    theme_coloc

  ggsave(file.path(dir_plots, "pp4_dotplot.pdf"), p_dot,
         width = max(8, length(unique(ok_df$locus)) * 1.2), height = 6,
         limitsize = FALSE)
  message("  pp4_dotplot.pdf")

  # ── 8b. Stacked PP bars ─────────────────────────────────────────────────
  ok_df$pair <- paste0(ok_df$gene, "\nL", ok_df$locus)
  ok_df$pair <- factor(ok_df$pair, levels = ok_df$pair[order(ok_df$PP.H4)])

  long_df <- ok_df %>%
    select(pair, PP.H0, PP.H1, PP.H2, PP.H3, PP.H4) %>%
    pivot_longer(cols = starts_with("PP."), names_to = "Hypothesis",
                 values_to = "Probability") %>%
    mutate(Hypothesis = factor(Hypothesis,
                               levels = c("PP.H0","PP.H1","PP.H2","PP.H3","PP.H4")))

  p_bar <- ggplot(long_df, aes(x = pair, y = Probability, fill = Hypothesis)) +
    geom_col(width = 0.75, colour = "white", linewidth = 0.2) +
    scale_fill_manual(
      values = PP_COLS,
      labels = c(H0 = "H0: no assoc.",
                 H1 = "H1: GWAS only",
                 H2 = "H2: eQTL only",
                 H3 = "H3: both, diff. SNP",
                 H4 = "H4: shared SNP  *")
    ) +
    coord_flip() +
    scale_y_continuous(expand = c(0, 0), limits = c(0, 1.01)) +
    labs(
      title    = "Posterior probability breakdown per gene x locus",
      subtitle = sprintf("Nightmare GWAS x Brain eQTL | %s | sorted by PP.H4",
                         filter_label),
      x = NULL, y = "Posterior probability", fill = "Hypothesis"
    ) +
    theme_coloc +
    theme(legend.position = "bottom", axis.text.y = element_text(size = 8))

  ggsave(file.path(dir_plots, "pp_barplot_by_locus.pdf"), p_bar,
         width = 9, height = max(5, nrow(ok_df) * 0.45), limitsize = FALSE)
  message("  pp_barplot_by_locus.pdf")

  # ── 8c. Heatmap ──────────────────────────────────────────────────────────
  if (nrow(ok_df) > 1) {
    heat_df <- ok_df %>%
      select(gene, locus, PP.H4) %>%
      mutate(locus = paste0("L", locus))

    p_heat <- ggplot(heat_df, aes(x = locus, y = gene, fill = PP.H4)) +
      geom_tile(colour = "white", linewidth = 0.5) +
      geom_text(aes(label = sprintf("%.2f", PP.H4)),
                size = 3.2, colour = "white",
                data = filter(heat_df, PP.H4 >= 0.5)) +
      scale_fill_gradientn(
        colours = c("#EFF3FF","#BDD7E7","#6BAED6","#2171B5","#CB181D"),
        values  = rescale(c(0, 0.2, 0.5, 0.8, 1)),
        limits  = c(0, 1), name = "PP.H4"
      ) +
      labs(
        title    = "PP.H4 heatmap: gene x locus",
        subtitle = sprintf("Nightmare GWAS x Brain eQTL | %s", filter_label),
        x = "Locus", y = "Gene"
      ) +
      theme_coloc +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))

    ggsave(file.path(dir_plots, "pp_heatmap.pdf"), p_heat,
           width  = max(6, length(unique(heat_df$locus)) * 0.9),
           height = max(4, length(unique(heat_df$gene))  * 0.5),
           limitsize = FALSE)
    message("  pp_heatmap.pdf")
  }

  # ── 8d. Bin summary bar ──────────────────────────────────────────────────
  bin_df <- ok_df %>%
    mutate(PP4_bin = cut(PP.H4,
                         breaks = c(-Inf, 0.2, 0.5, 0.8, 1.0),
                         labels = c("< 0.2","0.2-0.5","0.5-0.8",">= 0.8"),
                         include.lowest = TRUE)) %>%
    count(PP4_bin, name = "n_pairs")

  p_bin <- ggplot(bin_df, aes(x = PP4_bin, y = n_pairs, fill = PP4_bin)) +
    geom_col(width = 0.6, show.legend = FALSE) +
    geom_text(aes(label = n_pairs), vjust = -0.4, size = 4.5, fontface = "bold") +
    scale_fill_manual(values = c("< 0.2"   = "#DEEBF7",
                                 "0.2-0.5" = "#9ECAE1",
                                 "0.5-0.8" = "#FD8D3C",
                                 ">= 0.8"  = "#CB181D")) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.18))) +
    labs(
      title    = "Gene x locus pairs by PP.H4 bin",
      subtitle = sprintf("Nightmare GWAS x Brain eQTL | %s", filter_label),
      x = "PP.H4 bin", y = "Number of gene x locus pairs"
    ) +
    theme_coloc

  ggsave(file.path(dir_plots, "pp4_bin_summary.pdf"), p_bin, width = 6, height = 5)
  message("  pp4_bin_summary.pdf")
}

# ── 8e. GWAS vs eQTL scatter for colocalising pairs ─────────────────────────
coloc_pairs <- Filter(function(x) !is.na(x$PP.H4) && x$PP.H4 >= opt$h4_thresh,
                      scatter_data)
message(sprintf("  Generating %d scatter plot(s) for PP.H4 >= %.2f pairs...",
                length(coloc_pairs), opt$h4_thresh))

for (sp in coloc_pairs) {
  df            <- sp$merged
  df$gwas_logp  <- -log10(pmax(df$pvalue,       1e-300))
  df$eqtl_logp  <- -log10(pmax(df$pval_nominal, 1e-300))
  top5          <- head(df[order(df$pvalue), ], 5)

  p_scat <- ggplot(df, aes(x = gwas_logp, y = eqtl_logp)) +
    geom_point(aes(colour = gwas_logp), alpha = 0.65, size = 1.8) +
    geom_point(data = top5, colour = "#CB181D", size = 3.5, shape = 18) +
    geom_text_repel(data = top5, aes(label = snp), size = 2.8, colour = "#CB181D",
                    min.segment.length = 0.1, max.overlaps = 10) +
    scale_colour_gradientn(colours = c("#DEEBF7","#9ECAE1","#3182BD","#08519C"),
                           name = "-log10(p)\nGWAS") +
    labs(
      title    = sprintf("Locus %s  |  Gene: %s", sp$locus, sp$gene),
      subtitle = sprintf("PP.H4 = %.3f  |  chr%s  |  %d shared SNPs",
                         sp$PP.H4, sp$chr, nrow(df)),
      x = "-log10(p)  Nightmare GWAS",
      y = "-log10(p)  Brain eQTL"
    ) +
    theme_coloc

  ggsave(file.path(dir_plots,
                   sprintf("locus%s_%s_scatter.pdf", sp$locus, sp$gene)),
         p_scat, width = 6, height = 5)
}

# ---------------------------------------------------------------------------
# 9.  Console summary
# ---------------------------------------------------------------------------
n_ok   <- sum(summary_df$status == "OK",  na.rm = TRUE)
n_skip <- sum(summary_df$status != "OK",  na.rm = TRUE)
n_h4   <- sum(!is.na(summary_df$PP.H4) & summary_df$PP.H4 >= opt$h4_thresh)

message("\n", strrep("=", 58))
message("DONE")
message(strrep("=", 58))
message(sprintf("  Gene filter             : %s", filter_label))
message(sprintf("  Total gene x locus tests: %d", nrow(summary_df)))
message(sprintf("  Ran OK                  : %d", n_ok))
message(sprintf("  Skipped (too few SNPs)  : %d", n_skip))
message(sprintf("  PP.H4 >= %-5s          : %d", opt$h4_thresh, n_h4))
message("")
message("Output files:")
message("  ", file.path(dir_main, "coloc_summary_all.txt"),
        "     <- all results")
message("  ", file.path(dir_main, "coloc_locus_report.txt"),
        "  <- loci file + PP values")
message("  ", dir_pairs, "/          <- per gene x locus flat files")
message("  ", dir_plots, "/          <- PDF plots")

if (n_h4 > 0) {
  message(sprintf("\nColocalising pairs (PP.H4 >= %s):", opt$h4_thresh))
  hits <- summary_df[!is.na(summary_df$PP.H4) & summary_df$PP.H4 >= opt$h4_thresh,
                     c("locus","gene","nsnps","PP.H3","PP.H4")]
  print(hits, row.names = FALSE)
}
