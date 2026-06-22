################################################################################
# GENOMIC HALLMARKS OF DEPOT MEDROXYPROGESTERONE ACETATE-ASSOCIATED MENINGIOMAS
#
# RNA-seq analysis: TRAF7-mutant vs TRAF7-wildtype within the DMPA cohort.
#
# Self-contained script. The only required input is the raw featureCounts matrix
# (counts_file, set in Section 1). TRAF7 mutation status is taken from the
# institutional GlioSeq targeted panel (TRAF7_MUT / TRAF7_WT, defined below).
#
# FIGURES PRODUCED (saved to fig_dir):
#   volcano_TRAF7.pdf             — Volcano plot, nominal p-values
#   hallmark_barplot.pdf          — MSigDB Hallmark pathway enrichment
#   boxplot_PGR_supp.pdf          — PGR gene expression
#   boxplot_GO0032570_*.pdf       — Composite progesterone response pathway score
#   boxplot_GO0050847_*.pdf       — Composite PGR signaling pathway score
#
# OUTPUT TABLES (saved to results_dir):
#   Supplementary_Table_DGE_TRAF7mut_vs_WT.csv  — Full differential-expression results
#   fgsea_hallmark_TRAF7mut_vs_WT.csv            — Hallmark pathway enrichment results
################################################################################

##==============================================================================
###-----1. PACKAGES AND PATHS-----###
##==============================================================================

suppressPackageStartupMessages({
  pkgs_cran <- c("ggplot2", "dplyr", "tidyr", "stringr", "forcats",
                 "readr", "scales")
  pkgs_bioc <- c("edgeR", "limma", "fgsea", "msigdbr", "org.Hs.eg.db",
                 "AnnotationDbi")
  all_pkgs  <- c(pkgs_cran, pkgs_bioc)
  
  need_cran <- pkgs_cran[!sapply(pkgs_cran, requireNamespace, quietly = TRUE)]
  need_bioc <- pkgs_bioc[!sapply(pkgs_bioc, requireNamespace, quietly = TRUE)]
  
  if (length(need_cran)) install.packages(need_cran, repos = "https://cloud.r-project.org")
  if (length(need_bioc)) {
    if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
    BiocManager::install(need_bioc, ask = FALSE, update = FALSE)
  }
  
  lapply(all_pkgs, library, character.only = TRUE)
})

#dplyr verbs take precedence
select    <- dplyr::select
filter    <- dplyr::filter
mutate    <- dplyr::mutate
arrange   <- dplyr::arrange
left_join <- dplyr::left_join

###------Paths (EDIT these two)------###
rnaseq_dir  <- path.expand("~/Depo_Meningiomas/RNAseq")  # <-- output directory for figures/tables
counts_file <- path.expand("~/Depo_Meningiomas/RNAseq/DMPA_featureCounts_raw.txt")  # <-- raw featureCounts matrix
fig_dir     <- file.path(rnaseq_dir, "figures")
results_dir <- file.path(rnaseq_dir, "results")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)


###------Group membership (GlioSeq targeted panel)------###
# TRAF7-mutant: confirmed missense mutation on the GlioSeq panel.
# TRAF7-wildtype: no detected TRAF7 mutation.
# All 10 patients are NF2-wildtype, a central finding of the study.
TRAF7_MUT <- c("DMPA-2", "DMPA-4", "DMPA-7", "DMPA-9", "DMPA-10")
TRAF7_WT  <- c("DMPA-1", "DMPA-3", "DMPA-5", "DMPA-6", "DMPA-8")
DMPA_ALL  <- c(TRAF7_WT, TRAF7_MUT)

#Verification
stopifnot(length(TRAF7_MUT) == 5)
stopifnot(length(TRAF7_WT) == 5)
stopifnot(length(intersect(TRAF7_MUT, TRAF7_WT)) == 0)
stopifnot("DMPA-8" %in% TRAF7_WT)  # PIK3CA/FGFR1 patient must be in WT group

###------Shared aesthetics------###
GROUP_COLORS <- c("TRAF7-Mut" = "#d95f0e", "TRAF7-WT" = "#2c7fb8")
GROUP_LEVELS <- c("TRAF7-Mut", "TRAF7-WT")


##==============================================================================
###-----2. LOAD AND NORMALIZE EXPRESSION DATA-----###
##==============================================================================

cat("Loading count matrix...\n")

#Read raw featureCounts output (tab-delimited, first 6 cols are annotation)
stopifnot(file.exists(counts_file))

fc_raw <- read.delim(counts_file, comment.char = "#", check.names = FALSE)

#Count columns start at column 7 (1-6 are Geneid, Chr, Start, End, Strand, Length)
counts_mat <- as.matrix(fc_raw[, 7:ncol(fc_raw)])
rownames(counts_mat) <- fc_raw$Geneid

#Reduce featureCounts column headers (the original per-sample input file paths)
#to their DMPA-N sample identifiers, independent of upstream naming.
colnames(counts_mat) <- sub(".*(DMPA[-_][0-9]+).*", "\\1",
                            basename(colnames(counts_mat)))

#Subset to the 10 DMPA samples
missing_samples <- setdiff(DMPA_ALL, colnames(counts_mat))
if (length(missing_samples) > 0) {
  stop("Missing DMPA samples in count matrix: ", paste(missing_samples, collapse = ", "))
}
counts_dmpa <- counts_mat[, DMPA_ALL]

cat(sprintf("Count matrix: %d genes x %d samples\n", nrow(counts_dmpa), ncol(counts_dmpa)))

###-----TMM normalization via edgeR-----###
group_factor <- factor(
  ifelse(colnames(counts_dmpa) %in% TRAF7_MUT, "TRAF7-Mut", "TRAF7-WT"),
  levels = GROUP_LEVELS
)

dge <- edgeR::DGEList(counts = counts_dmpa, group = group_factor)

#Low-expression filter (edgeR::filterByExpr); permissive settings for the n=10 cohort
keep <- edgeR::filterByExpr(dge, min.count = 5, min.prop = 0.3)
dge  <- dge[keep, ]
dge  <- edgeR::calcNormFactors(dge, method = "TMM")
cat(sprintf("After filtering: %d genes\n", nrow(dge)))

#Log2 CPM matrix for plotting
logcpm <- edgeR::cpm(dge, log = TRUE, prior.count = 2)

###-----Sample metadata-----###
sample_meta <- data.frame(
  Sample  = colnames(logcpm),
  DMPA_ID = colnames(logcpm),
  Group   = factor(
    ifelse(colnames(logcpm) %in% TRAF7_MUT, "TRAF7-Mut", "TRAF7-WT"),
    levels = GROUP_LEVELS
  ),
  stringsAsFactors = FALSE
)

cat("\n===== TRAF7 GROUP ASSIGNMENTS =====\n")
print(table(sample_meta$Group))
cat("TRAF7-Mut:", paste(sample_meta$DMPA_ID[sample_meta$Group == "TRAF7-Mut"], collapse = ", "), "\n")
cat("TRAF7-WT: ", paste(sample_meta$DMPA_ID[sample_meta$Group == "TRAF7-WT"], collapse = ", "), "\n")
cat("====================================\n\n")


##==============================================================================
###-----3. GENE ANNOTATION (Ensembl → Symbol)-----###
##==============================================================================

#Strip Gencode version suffixes for annotation lookup
strip_version <- function(x) sub("\\.\\d+$", "", x)
gene_ids_stripped <- strip_version(rownames(dge))

#Map Ensembl IDs to gene symbols
ens2sym <- AnnotationDbi::mapIds(
  org.Hs.eg.db,
  keys    = gene_ids_stripped,
  keytype = "ENSEMBL",
  column  = "SYMBOL",
  multiVals = "first"
)

gene_symbols <- ens2sym[gene_ids_stripped]
names(gene_symbols) <- rownames(dge)

cat(sprintf("Mapped %d / %d genes to symbols\n",
            sum(!is.na(gene_symbols)), length(gene_symbols)))

###-----Build symbol-indexed logcpm matrix (deduplicate by AveExpr)-----###
has_symbol <- !is.na(gene_symbols)
logcpm_annotated <- logcpm[has_symbol, , drop = FALSE]
symbols_annotated <- gene_symbols[has_symbol]

#For duplicate symbols, keep the row with highest mean expression
avg_expr <- rowMeans(logcpm_annotated)
dedup_df <- data.frame(
  ensembl = rownames(logcpm_annotated),
  symbol  = symbols_annotated,
  avg     = avg_expr,
  stringsAsFactors = FALSE
) %>%
  group_by(symbol) %>%
  slice_max(avg, n = 1, with_ties = FALSE) %>%
  ungroup()

logcpm_sym <- logcpm_annotated[dedup_df$ensembl, , drop = FALSE]
rownames(logcpm_sym) <- dedup_df$symbol
cat(sprintf("Symbol-indexed matrix: %d unique genes\n", nrow(logcpm_sym)))


##==============================================================================
###-----4. LIMMA-VOOM DIFFERENTIAL EXPRESSION-----###
##==============================================================================

#Design matrix — use syntactically valid column names (limma requires this
#for makeContrasts); display labels elsewhere keep the hyphenated form.
design <- model.matrix(~ 0 + group_factor)
colnames(design) <- c("TRAF7_Mut", "TRAF7_WT")

#Voom transformation
v <- voom(dge, design, plot = FALSE)

#Fit and contrast
fit <- lmFit(v, design)
contrast_mat <- makeContrasts(
  TRAF7_effect = TRAF7_Mut - TRAF7_WT,
  levels = design
)
fit2 <- contrasts.fit(fit, contrast_mat)
fit2 <- eBayes(fit2)

#Extract full results
limma_res <- topTable(fit2, coef = "TRAF7_effect", number = Inf, sort.by = "none")
limma_res$gene_id <- rownames(limma_res)
limma_res$gene_id_stripped <- strip_version(limma_res$gene_id)
limma_res$SYMBOL <- gene_symbols[limma_res$gene_id]

n_up   <- sum(limma_res$P.Value < 0.05 & limma_res$logFC >= 1, na.rm = TRUE)
n_down <- sum(limma_res$P.Value < 0.05 & limma_res$logFC <= -1, na.rm = TRUE)
cat(sprintf("\nDifferential expression (nominal p<0.05, |LFC|>=1):\n"))
cat(sprintf("  Up in TRAF7-Mut:   %d\n", n_up))
cat(sprintf("  Down in TRAF7-Mut: %d\n", n_down))
cat(sprintf("  FDR < 0.05:        %d\n", sum(limma_res$adj.P.Val < 0.05, na.rm = TRUE)))


##==============================================================================
###-----5. fgsea — HALLMARK PATHWAYS-----###
##==============================================================================

#Ranked gene list (t-statistic, named by symbol)
ranks_traf7 <- setNames(limma_res$t, limma_res$SYMBOL)
ranks_traf7 <- ranks_traf7[!is.na(names(ranks_traf7)) & names(ranks_traf7) != ""]
if (any(duplicated(names(ranks_traf7)))) {
  ranks_traf7 <- ranks_traf7[!duplicated(names(ranks_traf7))]
}
ranks_traf7 <- sort(ranks_traf7, decreasing = TRUE)

#MSigDB Hallmark gene sets
hallmark_sets <- msigdbr(species = "Homo sapiens", collection = "H") %>%
  split(x = .$gene_symbol, f = .$gs_name)

set.seed(42)
fgsea_hall_traf7 <- fgsea(
  pathways    = hallmark_sets,
  stats       = ranks_traf7,
  minSize     = 15,
  maxSize     = 500,
  nPermSimple = 10000
)

hallmark <- as.data.frame(fgsea_hall_traf7)
cat(sprintf("\nHallmark fgsea: %d pathways at padj < 0.25\n",
            sum(hallmark$padj < 0.25, na.rm = TRUE)))


##==============================================================================
###-----6. fgsea — PROGESTERONE PATHWAYS-----###
##==============================================================================

#Pull progesterone-pathway gene sets directly from org.Hs.eg.db via GOALL
#(includes child terms through GO DAG propagation — matches methylation script)
get_go_genes <- function(go_id) {
  df <- AnnotationDbi::select(
    org.Hs.eg.db,
    keys    = go_id,
    keytype = "GOALL",
    columns = "SYMBOL"
  )
  unique(df$SYMBOL[!is.na(df$SYMBOL)])
}

go_0032570_genes <- get_go_genes("GO:0032570")
go_0050847_genes <- get_go_genes("GO:0050847")

cat(sprintf("GO:0032570 (response to progesterone): %d genes\n", length(go_0032570_genes)))
cat(sprintf("GO:0050847 (PR signaling pathway): %d genes\n", length(go_0050847_genes)))

fgsea_prog_traf7 <- fgsea(
  pathways = list(
    "GO:0032570_Response_to_progesterone" = go_0032570_genes,
    "GO:0050847_PR_signaling"             = go_0050847_genes
  ),
  stats       = ranks_traf7,
  minSize     = 3,
  maxSize     = 500,
  nPermSimple = 10000
)

cat("\n--- Progesterone pathway fgsea ---\n")
print(fgsea_prog_traf7[, c("pathway", "NES", "pval", "padj", "size")])


##==============================================================================
###-----7. SUPPLEMENTARY TABLES-----###
##==============================================================================

#Full DGE results
write.csv(limma_res, file.path(results_dir, "Supplementary_Table_DGE_TRAF7mut_vs_WT.csv"),
          row.names = FALSE)
cat(sprintf("\nSaved DGE table: %d genes\n", nrow(limma_res)))

#Hallmark fgsea results — flatten leadingEdge list column for CSV export
hallmark_out <- hallmark %>%
  mutate(leadingEdge = sapply(leadingEdge, paste, collapse = ";"))
write.csv(hallmark_out, file.path(results_dir, "fgsea_hallmark_TRAF7mut_vs_WT.csv"),
          row.names = FALSE)


##==============================================================================
###-----8. VOLCANO PLOT-----###
##==============================================================================

LOGFC_CUT <- 1
PVAL_CUT  <- 0.05

volcano_df <- limma_res %>%
  mutate(
    neg_log10_p = -log10(pmax(P.Value, .Machine$double.xmin)),
    sig = case_when(
      P.Value < PVAL_CUT & logFC >=  LOGFC_CUT ~ "Up in TRAF7-Mut",
      P.Value < PVAL_CUT & logFC <= -LOGFC_CUT ~ "Down in TRAF7-Mut",
      TRUE ~ "NS"
    )
  )

n_up   <- sum(volcano_df$sig == "Up in TRAF7-Mut")
n_down <- sum(volcano_df$sig == "Down in TRAF7-Mut")
n_ns   <- sum(volcano_df$sig == "NS")

volcano_df <- volcano_df %>% arrange(sig == "NS")

p_volcano <- ggplot(volcano_df, aes(x = logFC, y = neg_log10_p)) +
  geom_point(aes(color = sig), alpha = 0.85, size = 1.8) +
  geom_vline(xintercept = c(-LOGFC_CUT, LOGFC_CUT),
             linetype = "dashed", color = "grey55", linewidth = 0.4) +
  geom_hline(yintercept = -log10(PVAL_CUT),
             linetype = "dashed", color = "grey55", linewidth = 0.4) +
  scale_color_manual(
    values = c("Up in TRAF7-Mut"   = "#d95f0e",
               "NS"                = "grey55",
               "Down in TRAF7-Mut" = "#2c7fb8"),
    labels = c(
      "Up in TRAF7-Mut"   = sprintf("Up in TRAF7-Mut (n=%d)", n_up),
      "Down in TRAF7-Mut" = sprintf("Down in TRAF7-Mut (n=%d)", n_down),
      "NS"                = sprintf("Not significant (n=%s)", formatC(n_ns, big.mark = ","))
    ),
    breaks = c("Up in TRAF7-Mut", "Down in TRAF7-Mut", "NS")
  ) +
  labs(
    title = "TRAF7-Mut vs. TRAF7-WT differential gene expression",
    x     = expression(Log[2]*" fold change"),
    y     = expression("-"*Log[10]*" p-value (nominal)"),
    color = NULL
  ) +
  theme_classic(base_size = 16) +
  theme(
    plot.title       = element_text(hjust = 0.5, face = "bold"),
    plot.title.position = "plot", 
    axis.line.x      = element_line(linewidth = 1.1, color = "black"),
    axis.line.y      = element_line(linewidth = 1.1, color = "black"),
    axis.ticks        = element_line(linewidth = 0.8, color = "black"),
    axis.ticks.length = unit(5, "pt"),
    axis.text.x       = element_text(size = 13),
    axis.text.y       = element_text(size = 13),
    axis.title.x      = element_text(size = 15, margin = ggplot2::margin(t = 10)),
    axis.title.y      = element_text(size = 15, margin = ggplot2::margin(r = 16)),
    legend.position      = c(0.99, 0.03),
    legend.justification = c("right", "bottom"),
    legend.background    = element_rect(fill = "transparent", color = NA),
    legend.key           = element_rect(fill = "transparent", color = NA),
    legend.text          = element_text(size = 10),
    legend.key.size      = unit(10, "pt"),
    legend.spacing.y     = unit(1, "pt"),
    panel.grid           = element_blank()
  ) +
  coord_cartesian(clip = "off")

ggsave(file.path(fig_dir, "volcano_TRAF7.pdf"),
       p_volcano, width = 8, height = 7)
ggsave(file.path(fig_dir, "volcano_TRAF7.png"),
       p_volcano, width = 8, height = 7, dpi = 300)
message("\u2705 Saved: volcano_TRAF7.pdf / .png")

ggsave(
  file.path(fig_dir, "volcano_TRAF7_600dpi.tiff"),
  p_volcano,
  width = 8,
  height = 7,
  dpi = 600,
  compression = "lzw"
)

##==============================================================================
###-----9. PATHWAY ENRICHMENT BARPLOT-----###
##==============================================================================

hallmark_plot <- hallmark %>%
  filter(padj < 0.25) %>%
  mutate(
    pathway_label = pathway %>%
      str_replace("^HALLMARK_", "") %>%
      str_replace_all("_", " ") %>%
      str_to_title() %>%
      str_replace("Tnfa",  "TNFα") %>%
      str_replace("Nfkb",  "NFκB") %>%
      str_replace("Il2",   "IL2") %>%
      str_replace("Il6",   "IL6") %>%
      str_replace("Stat3", "STAT3") %>%
      str_replace("Stat5", "Stat5") %>%
      str_replace("Jak",   "Jak") %>%
      str_replace("Kras",  "KRAS") %>%
      str_replace("Myc",   "MYC") %>%
      str_replace("E2f",   "E2F") %>%
      str_replace("Dna",   "DNA") %>%
      str_replace("G2m",   "G2M") %>%
      str_replace("Uv",    "UV"),
    Direction = ifelse(NES > 0, "Up in TRAF7-Mut", "Down in TRAF7-Mut")
  ) %>%
  arrange(NES) %>%
  mutate(pathway_label = fct_inorder(pathway_label))

p_hallmark <- ggplot(hallmark_plot, aes(x = NES, y = pathway_label, fill = Direction)) +
  geom_bar(stat = "identity", width = 0.72) +
  geom_vline(xintercept = 0, color = "black", linewidth = 0.5) +
  scale_fill_manual(
    values = c("Up in TRAF7-Mut" = "#d95f0e", "Down in TRAF7-Mut" = "#2c7fb8"),
    breaks = c("Up in TRAF7-Mut", "Down in TRAF7-Mut")
  ) +
  scale_x_continuous(expand = expansion(mult = 0.05)) +
  labs(
    title = "TRAF7-Mut vs. TRAF7-WT pathway enrichment analysis",
    x     = "Normalized Enrichment Score",
    y     = NULL,
    fill  = NULL
  ) +
  theme_classic(base_size = 16) +
  theme(
    plot.title          = element_text(hjust = 0.5, face = "bold"),
    plot.title.position = "plot",
    axis.line.x         = element_line(linewidth = 1.1, color = "black"),
    axis.line.y         = element_line(linewidth = 1.1, color = "black"),
    axis.ticks          = element_line(linewidth = 0.8, color = "black"),
    axis.ticks.length   = unit(5, "pt"),
    axis.text.x         = element_text(size = 13),
    axis.text.y         = element_text(size = 12),
    axis.title.x        = element_text(size = 15, margin = ggplot2::margin(t = 10)),
    legend.position      = c(0.99, 0.02),
    legend.justification = c("right", "bottom"),
    legend.background    = element_rect(fill = alpha("white", 0.85), color = NA),
    legend.key           = element_rect(fill = "transparent", color = NA),
    legend.text          = element_text(size = 11),
    legend.key.size      = unit(12, "pt"),
    panel.grid           = element_blank()
  )

ggsave(file.path(fig_dir, "hallmark_barplot.pdf"),
       p_hallmark, width = 8, height = 7, device = cairo_pdf)
ggsave(file.path(fig_dir, "hallmark_barplot.png"),
       p_hallmark, width = 8, height = 7, dpi = 300)
message("\u2705 Saved: hallmark_barplot.pdf / .png")


##==============================================================================
###-----10. PGR EXPRESSION BOXPLOT-----###
##==============================================================================

if ("PGR" %in% rownames(logcpm_sym)) {
  
  df_pgr <- data.frame(
    Expression = logcpm_sym["PGR", ],
    Group      = sample_meta$Group,
    stringsAsFactors = FALSE
  )
  
  gene_row <- limma_res %>% filter(SYMBOL == "PGR")
  if (nrow(gene_row) > 1) gene_row <- gene_row %>% slice_max(AveExpr, n = 1)
  if (nrow(gene_row) == 1) {
    cat(sprintf("\nPGR: p = %.3f, log2FC = %.2f\n", gene_row$P.Value, gene_row$logFC))
  }
  
  p_pgr <- ggplot(df_pgr, aes(x = Group, y = Expression, fill = Group)) +
    stat_boxplot(geom = "errorbar", width = 0.25, linewidth = 0.8) +
    geom_boxplot(width = 0.5, outlier.shape = NA, alpha = 0.7, color = NA) +
    geom_boxplot(width = 0.5, outlier.shape = NA, fill = NA,
                 color = "black", linewidth = 0.6) +
    geom_jitter(width = 0.12, size = 2.5, alpha = 0.8, aes(color = Group)) +
    scale_fill_manual(values  = GROUP_COLORS) +
    scale_color_manual(values = GROUP_COLORS) +
    scale_x_discrete(labels = c("TRAF7-Mut" = "TRAF7\nMut",
                                "TRAF7-WT"  = "TRAF7\nWT")) +
    labs(
      title = "PGR gene expression",
      x     = NULL,
      y     = expression(log[2]*" CPM")
    ) +
    theme_classic(base_size = 16) +
    theme(
      plot.title       = element_text(hjust = 0.5, face = "bold"),
      axis.line.x      = element_line(linewidth = 1.1, color = "black"),
      axis.line.y      = element_line(linewidth = 1.1, color = "black"),
      axis.ticks        = element_line(linewidth = 0.8, color = "black"),
      axis.ticks.length = unit(5, "pt"),
      axis.text.x       = element_text(size = 14, face = "bold", lineheight = 0.9),
      axis.text.y       = element_text(size = 13),
      axis.title.y      = element_text(size = 15, margin = ggplot2::margin(r = 10)),
      legend.position   = "none",
      panel.grid        = element_blank(),
      plot.margin       = ggplot2::margin(t = 25, r = 70, b = 10, l = 70)
    )
  
  ggsave(file.path(fig_dir, "boxplot_PGR_supp.pdf"),
         p_pgr, width = 6, height = 6.5)
  ggsave(file.path(fig_dir, "boxplot_PGR_supp.png"),
         p_pgr, width = 6, height = 6.5, dpi = 300)
  message("\u2705 Saved: boxplot_PGR_supp.pdf / .png")
  
} else {
  message("\u26A0\uFE0F  PGR not found in logcpm_sym — skipping PGR boxplot.")
}


##==============================================================================
###-----11. COMPOSITE PROGESTERONE PATHWAY EXPRESSION SCORES-----###
##==============================================================================

#Define pathway info
pathway_info <- list(
  list(
    genes    = go_0032570_genes,
    go_id    = "GO:0032570",
    title    = "Composite progesterone response\npathway expression score",
    filename = "boxplot_GO0032570_response_to_progesterone"
  ),
  list(
    genes    = go_0050847_genes,
    go_id    = "GO:0050847",
    title    = "Composite PGR signaling\npathway expression score",
    filename = "boxplot_GO0050847_PR_signaling"
  )
)

for (pw in pathway_info) {
  
  pw_genes_avail <- intersect(pw$genes, rownames(logcpm_sym))
  cat(sprintf("\n%s: %d / %d genes present in expression data\n",
              pw$go_id, length(pw_genes_avail), length(pw$genes)))
  
  if (length(pw_genes_avail) < 3) {
    message("\u26A0\uFE0F  Fewer than 3 genes available for ", pw$go_id, " — skipping.")
    next
  }
  
  #Composite score = mean log2 CPM across pathway genes per sample
  pw_matrix <- logcpm_sym[pw_genes_avail, , drop = FALSE]
  composite_score <- colMeans(pw_matrix, na.rm = TRUE)
  
  df_pw <- data.frame(
    Score = composite_score,
    Group = sample_meta$Group,
    stringsAsFactors = FALSE
  )
  
  #Wilcoxon test
  wt <- wilcox.test(Score ~ Group, data = df_pw, exact = FALSE)
  pval_pw <- wt$p.value
  plab_pw <- if (pval_pw < 0.001) sprintf("p = %.2e", pval_pw) else sprintf("p = %.3f", pval_pw)
  cat(sprintf("  Wilcoxon %s, n genes = %d\n", plab_pw, length(pw_genes_avail)))
  
  p_pw <- ggplot(df_pw, aes(x = Group, y = Score, fill = Group)) +
    stat_boxplot(geom = "errorbar", width = 0.25, linewidth = 0.8) +
    geom_boxplot(width = 0.5, outlier.shape = NA, alpha = 0.7, color = NA) +
    geom_boxplot(width = 0.5, outlier.shape = NA, fill = NA,
                 color = "black", linewidth = 0.6) +
    geom_jitter(width = 0.12, size = 2.5, alpha = 0.8, aes(color = Group)) +
    scale_fill_manual(values  = GROUP_COLORS) +
    scale_color_manual(values = GROUP_COLORS) +
    scale_x_discrete(labels = c("TRAF7-Mut" = "TRAF7\nMut",
                                "TRAF7-WT"  = "TRAF7\nWT")) +
    labs(
      title = pw$title,
      x     = NULL,
      y     = expression("Mean "*log[2]*" CPM")
    ) +
    theme_classic(base_size = 16) +
    theme(
      plot.title        = element_text(hjust = 0.5, face = "bold"),
      axis.line.x       = element_line(linewidth = 1.1, color = "black"),
      axis.line.y       = element_line(linewidth = 1.1, color = "black"),
      axis.ticks        = element_line(linewidth = 0.8, color = "black"),
      axis.ticks.length = unit(5, "pt"),
      axis.text.x       = element_text(size = 14, face = "bold", lineheight = 0.9),
      axis.text.y       = element_text(size = 13),
      axis.title.y      = element_text(size = 15, margin = ggplot2::margin(r = 10)),
      legend.position   = "none",
      panel.grid        = element_blank(),
      plot.margin       = ggplot2::margin(t = 25, r = 70, b = 10, l = 70)
    )
  
  ggsave(file.path(fig_dir, sprintf("%s.pdf", pw$filename)),
         p_pw, width = 6, height = 6.5)
  ggsave(file.path(fig_dir, sprintf("%s.png", pw$filename)),
         p_pw, width = 6, height = 6.5, dpi = 300)
  message(sprintf("\u2705 Saved: %s.pdf / .png", pw$filename))
}


##==============================================================================
###-----12. COMPLETION SUMMARY-----###
##==============================================================================

cat("\n\n")
cat("=========================================================\n")
cat("RNA-seq analysis complete\n")
cat("=========================================================\n")
cat("\nTRAF7-Mut: DMPA-2, DMPA-4, DMPA-7, DMPA-9, DMPA-10\n")
cat("TRAF7-WT:  DMPA-1, DMPA-3, DMPA-5, DMPA-6, DMPA-8\n")
cat(sprintf("\nDifferential expression: %d up, %d down (nominal p<0.05, |LFC|>=1)\n",
            sum(limma_res$P.Value < 0.05 & limma_res$logFC >= 1, na.rm = TRUE),
            sum(limma_res$P.Value < 0.05 & limma_res$logFC <= -1, na.rm = TRUE)))
cat(sprintf("Hallmark pathways at padj < 0.25: %d\n",
            sum(hallmark$padj < 0.25, na.rm = TRUE)))
cat(sprintf("\nFigures saved to: %s\n", normalizePath(fig_dir)))
cat(sprintf("Tables saved to:  %s\n", normalizePath(results_dir)))
cat("=========================================================\n")

sessionInfo()



