################################################################################
# TCGA COAD/READ – MYC-Centered Integrative Analysis
# Author  : Nitish Kumar R
#
# Description:
#   End-to-end reproducible pipeline for MYC-driven colorectal cancer
#   transcriptomics using TCGA COAD + READ RNA-seq data (STAR counts).
#
# Outputs
#   figures/   – 15 publication-ready TIFF figures (600 DPI)
#   tables/    – 12 CSV result tables
#
# Sections:
#   00. Setup & packages
#   01. TCGA download & preparation
#   02. Sample annotation & QC
#   03. Gene filtering
#   04. DESeq2  – Tumor vs Normal
#   05. PCA
#   06. Heatmap – Tumor vs Normal (top 50 DEGs)
#   07. MYC stratification
#   08. DESeq2  – MYC-High vs MYC-Low (+ LFC shrinkage)
#   09. Heatmap – MYC-High vs MYC-Low (top 30 DEGs)
#   10. Volcano plots – Tumor vs Normal & MYC High/Low
#   11. GO enrichment – MYC-High upregulated genes
#   12. WGCNA  – Tumor samples
#   13. MYC-correlated module & annotated hub genes
#   14. GO enrichment – ME4 module hub genes
#   15. GSEA   – MYC + extended Hallmark gene sets
#   16. Survival – ME4 score KM + continuous & stage-adjusted Cox
#   17. ML dataset – leakage-free train/test split (gene symbols)
#   18. ceRNA network – miRTarBase + ENCORI
#
# External files needed only for section 18:
#   miRTarBase_MTI.csv       https://mirtarbase.cuhk.edu.cn
#   ENCORI_miRNA_lncRNA.tsv  https://rnasysu.com/encori
################################################################################


################################################################################
## 00. SETUP
################################################################################

set.seed(123)
options(stringsAsFactors = FALSE)

dir.create("figures", showWarnings = FALSE)
dir.create("tables",  showWarnings = FALSE)

# ── Auto-install missing packages (runs once; skips already-installed) ─────────

# CRAN packages
cran_pkgs <- c(
  "ggplot2", "ggrepel", "dplyr", "stringr",
  "data.table", "tibble", "survival", "survminer", "scales"
)
missing_cran <- cran_pkgs[!cran_pkgs %in% rownames(installed.packages())]
if (length(missing_cran)) {
  message("Installing missing CRAN packages: ",
          paste(missing_cran, collapse = ", "))
  install.packages(missing_cran, dependencies = TRUE)
}

# Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

bioc_pkgs <- c(
  "TCGAbiolinks", "SummarizedExperiment", "DESeq2",
  "clusterProfiler", "org.Hs.eg.db", "AnnotationDbi",
  "msigdbr", "enrichplot", "WGCNA", "pheatmap"
)
missing_bioc <- bioc_pkgs[!bioc_pkgs %in% rownames(installed.packages())]
if (length(missing_bioc)) {
  message("Installing missing Bioconductor packages: ",
          paste(missing_bioc, collapse = ", "))
  BiocManager::install(missing_bioc, ask = FALSE, update = FALSE)
}

message("All packages present.")

suppressPackageStartupMessages({
  library(TCGAbiolinks)
  library(SummarizedExperiment)
  library(DESeq2)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
  library(msigdbr)
  library(enrichplot)
  library(WGCNA)
  library(pheatmap)
  library(ggplot2)
  library(ggrepel)
  library(dplyr)
  library(stringr)
  library(data.table)
  library(tibble)
  library(survival)
  library(survminer)
  library(scales)
})

enableWGCNAThreads()

# CRITICAL: WGCNA::blockwiseModules internally dispatches WGCNA::cor(), which
# accepts extra arguments (weights.x, weights.y, cosine) that stats::cor()
# does not. The workaround is to swap cor in the global environment:
#   • set cor <- WGCNA::cor   immediately BEFORE each blockwiseModules() call
#   • set cor <- stats::cor   immediately AFTER  each blockwiseModules() call
# All other correlation calls in this script (moduleTraitCor, geneSignificance,
# modCor_train) use stats::cor via the restored binding.
cor <- stats::cor

# ── Shared helper functions ────────────────────────────────────────────────────

# Save a ggsurvplot (requires explicit print inside tiff device)
save_km <- function(fit, data, filename, title,
                    legend_title, legend_labs, palette,
                    break_by = 1000, width = 7, height = 7) {
  tiff(filename, width = width, height = height,
       units = "in", res = 600, compression = "lzw")
  print(
    ggsurvplot(
      fit, data         = data,
      risk.table        = TRUE,
      pval              = TRUE,
      conf.int          = FALSE,
      censor            = TRUE,
      censor.size       = 2.5,
      censor.shape      = 124,
      risk.table.height = 0.28,
      break.time.by     = break_by,
      xlab              = "Time (days)",
      ylab              = "Overall Survival Probability",
      legend.title      = legend_title,
      legend.labs       = legend_labs,
      palette           = palette,
      ggtheme           = theme_classic(base_size = 14),
      title             = title
    )
  )
  dev.off()
}

# Convert Ensembl rownames (with optional version suffix) to HGNC symbols.
# Removes rows with no symbol and duplicate symbols.
ensembl_to_symbols <- function(expr_matrix) {
  ens_clean <- sub("\\..*$", "", rownames(expr_matrix))
  syms <- mapIds(org.Hs.eg.db,
                 keys      = ens_clean,
                 keytype   = "ENSEMBL",
                 column    = "SYMBOL",
                 multiVals = "first")
  rownames(expr_matrix) <- syms
  expr_matrix <- expr_matrix[!is.na(rownames(expr_matrix)), , drop = FALSE]
  expr_matrix <- expr_matrix[!duplicated(rownames(expr_matrix)), , drop = FALSE]
  expr_matrix
}

# Build a publication-quality volcano plot with top-N gene labels
make_volcano <- function(df, lfc_col, padj_col, label_col,
                         color_col, color_pal, up_level, down_level, title) {
  top_up <- df %>%
    filter(.data[[color_col]] == up_level) %>%
    arrange(desc(abs(.data[[lfc_col]]))) %>%
    slice_head(n = 5)
  top_dn <- df %>%
    filter(.data[[color_col]] == down_level) %>%
    arrange(.data[[lfc_col]]) %>%
    slice_head(n = 5)
  labels <- bind_rows(top_up, top_dn)

  ggplot(df, aes(.data[[lfc_col]], -log10(.data[[padj_col]]))) +
    geom_point(aes(color = .data[[color_col]]),
               alpha = 0.65, size = 1.6, show.legend = TRUE) +
    geom_text_repel(
      data          = labels,
      aes(label     = .data[[label_col]]),
      size          = 3.2,
      box.padding   = 0.35,
      point.padding = 0.25,
      segment.color = "grey50",
      max.overlaps  = Inf
    ) +
    scale_color_manual(values = color_pal) +
    geom_vline(xintercept = c(-1, 1),
               linetype = "dashed", linewidth = 0.4) +
    geom_hline(yintercept = -log10(0.05),
               linetype = "dashed", linewidth = 0.4) +
    coord_cartesian(clip = "off") +
    labs(
      title = title,
      x     = expression(Log[2] ~ Fold ~ Change),
      y     = expression(-Log[10] ~ Adjusted ~ P-value),
      color = NULL
    ) +
    theme_bw(base_size = 11) +
    theme(
      plot.title  = element_text(hjust = 0.5, face = "bold"),
      plot.margin = margin(10, 25, 10, 10)
    )
}

# Render a GO dotplot and save as 600 DPI TIFF
save_go_dotplot <- function(ego_obj, filename, title,
                             show_n = 10, wrap_width = 38,
                             width = 8, height = 6) {
  p <- dotplot(ego_obj, showCategory = show_n, font.size = 11) +
    scale_y_discrete(labels = function(x) str_wrap(x, width = wrap_width)) +
    scale_color_gradient(low = "#2166AC", high = "#B2182B",
                         name = "Adjusted\np-value") +
    labs(title = title, x = "Gene Ratio", y = NULL) +
    theme_bw() +
    theme(
      plot.title         = element_text(hjust = 0.5, face = "bold", size = 13),
      plot.margin        = margin(20, 30, 15, 20),
      panel.grid.major.y = element_blank(),
      legend.position    = "right"
    ) +
    coord_cartesian(clip = "off")
  ggsave(filename, p, dpi = 600, width = width, height = height,
         device = "tiff", compression = "lzw")
}


################################################################################
## 01. TCGA DOWNLOAD & PREPARATION
################################################################################

query <- GDCquery(
  project       = c("TCGA-COAD", "TCGA-READ"),
  data.category = "Transcriptome Profiling",
  data.type     = "Gene Expression Quantification",
  workflow.type = "STAR - Counts"
)

GDCdownload(query)
rna_se  <- GDCprepare(query)
counts  <- assay(rna_se)
coldata <- as.data.frame(colData(rna_se))


################################################################################
## 02. SAMPLE ANNOTATION & QC
################################################################################

message("Raw matrix: ", nrow(counts), " genes x ", ncol(counts), " samples")
print(table(coldata$shortLetterCode))

coldata$condition <- dplyr::case_when(
  coldata$shortLetterCode == "TP" ~ "Tumor",
  coldata$shortLetterCode == "NT" ~ "Normal",
  TRUE                            ~ NA_character_
)

keep    <- !is.na(coldata$condition)
counts  <- counts[,  keep]
coldata <- coldata[keep, ]
coldata$condition <- factor(coldata$condition, levels = c("Normal", "Tumor"))

message("Samples kept  – Tumor: ", sum(coldata$condition == "Tumor"),
        "  Normal: ", sum(coldata$condition == "Normal"))


################################################################################
## 03. GENE FILTERING  (≥10 counts in ≥20% of samples)
################################################################################

keep_genes  <- rowSums(counts >= 10) >= 0.2 * ncol(counts)
counts_filt <- counts[keep_genes, ]

message("Genes after filtering: ", nrow(counts_filt),
        "  (removed ", sum(!keep_genes), ")")


################################################################################
## 04. DESEQ2 – TUMOR vs NORMAL
################################################################################

dds <- DESeqDataSetFromMatrix(
  countData = counts_filt,
  colData   = coldata,
  design    = ~ condition
)
dds    <- DESeq(dds)
res_TN <- results(dds, contrast = c("condition", "Tumor", "Normal"))

# VST expression matrix — used in every downstream step
vsd       <- vst(dds, blind = TRUE)
norm_expr <- assay(vsd)

res_TN_df <- as.data.frame(res_TN) %>%
  rownames_to_column("ensembl") %>%
  mutate(
    gene_symbol  = rowData(rna_se)$gene_name[
      match(ensembl, rownames(rowData(rna_se)))
    ],
    Significance = case_when(
      padj < 0.05 & log2FoldChange >  1 ~ "Upregulated",
      padj < 0.05 & log2FoldChange < -1 ~ "Downregulated",
      TRUE                              ~ "Not Significant"
    )
  )

write.csv(res_TN_df, "tables/01_DGE_Tumor_vs_Normal.csv", row.names = FALSE)

message("Tumor vs Normal – up: ",
        sum(res_TN_df$Significance == "Upregulated",   na.rm = TRUE),
        "  down: ",
        sum(res_TN_df$Significance == "Downregulated", na.rm = TRUE))


################################################################################
## 05. PCA  (600 DPI)
################################################################################

pca_df     <- plotPCA(vsd, intgroup = "condition", ntop = 500, returnData = TRUE)
percentVar <- round(100 * attr(pca_df, "percentVar"))

p_pca <- ggplot(pca_df, aes(PC1, PC2, color = condition)) +
  geom_point(size = 3, alpha = 0.8) +
  scale_color_manual(values = c(Normal = "#00BFC4", Tumor = "#F8766D")) +
  xlab(paste0("PC1 (", percentVar[1], "%)")) +
  ylab(paste0("PC2 (", percentVar[2], "%)")) +
  ggtitle("PCA – TCGA Colorectal Cancer: Tumor vs Normal") +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

ggsave("figures/Fig_01_PCA_Tumor_vs_Normal.tiff",
       p_pca, dpi = 600, width = 6, height = 5,
       device = "tiff", compression = "lzw")


################################################################################
## 06. HEATMAP – TUMOR vs NORMAL  (top 50 DEGs, 600 DPI)
################################################################################

top50_ens <- res_TN_df %>%
  filter(!is.na(padj), !is.na(gene_symbol)) %>%
  arrange(padj) %>%
  slice_head(n = 50) %>%
  pull(ensembl)

expr_heat_TN <- norm_expr[top50_ens, ]
rownames(expr_heat_TN) <- res_TN_df$gene_symbol[
  match(top50_ens, res_TN_df$ensembl)
]

# Row-wise Z-score, clipped to [-2, 2]
z_TN <- t(scale(t(expr_heat_TN)))
z_TN[z_TN >  2] <-  2
z_TN[z_TN < -2] <- -2

sample_order_TN <- rownames(coldata)[order(coldata$condition)]
ann_col_TN      <- data.frame(Condition = coldata$condition,
                               row.names = rownames(coldata))
ann_col_TN      <- ann_col_TN[sample_order_TN, , drop = FALSE]
ann_colors_TN   <- list(Condition = c(Normal = "#00BFC4", Tumor = "#F8766D"))

tiff("figures/Fig_02_Heatmap_Tumor_vs_Normal.tiff",
     width = 8, height = 9, units = "in", res = 600, compression = "lzw")
pheatmap(
  z_TN[, sample_order_TN],
  annotation_col    = ann_col_TN,
  annotation_colors = ann_colors_TN,
  show_colnames     = FALSE,
  show_rownames     = TRUE,
  fontsize_row      = 7,
  cluster_cols      = FALSE,
  cluster_rows      = TRUE,
  color             = colorRampPalette(c("navy", "white", "firebrick3"))(100),
  breaks            = seq(-2, 2, length.out = 101),
  border_color      = NA,
  main              = "Top 50 DEGs: Tumor vs Normal (CRC)"
)
dev.off()


################################################################################
## 07. MYC STRATIFICATION
################################################################################

myc_row <- which(rowData(rna_se)$gene_name == "MYC")
stopifnot("MYC not found in rowData$gene_name" = length(myc_row) == 1)

myc_exp <- norm_expr[myc_row, ]      # named numeric vector, all samples

coldata$MYC_group <- factor(
  ifelse(myc_exp >= median(myc_exp, na.rm = TRUE), "MYC_high", "MYC_low"),
  levels = c("MYC_low", "MYC_high")
)

message("MYC stratification – high: ",
        sum(coldata$MYC_group == "MYC_high"),
        "  low: ",
        sum(coldata$MYC_group == "MYC_low"))


################################################################################
## 08. DESEQ2 – MYC-HIGH vs MYC-LOW  (+ LFC shrinkage)
################################################################################

dds_myc <- DESeqDataSetFromMatrix(
  countData = counts_filt,
  colData   = coldata,
  design    = ~ MYC_group
)
dds_myc$MYC_group <- factor(dds_myc$MYC_group,
                             levels = c("MYC_low", "MYC_high"))
dds_myc <- DESeq(dds_myc)

res_myc <- results(dds_myc,
                   contrast = c("MYC_group", "MYC_high", "MYC_low"))

# apeglm LFC shrinkage: stabilises estimates for GSEA and volcano ranking
res_myc <- lfcShrink(
  dds_myc,
  coef = "MYC_group_MYC_high_vs_MYC_low",
  res  = res_myc
)

res_myc_df <- as.data.frame(res_myc) %>%
  rownames_to_column("ensembl") %>%
  mutate(
    gene_symbol  = rowData(rna_se)$gene_name[
      match(ensembl, rownames(rowData(rna_se)))
    ],
    Significance = case_when(
      padj < 0.05 & log2FoldChange >  1 ~ "MYC-High",
      padj < 0.05 & log2FoldChange < -1 ~ "MYC-Low",
      TRUE                              ~ "Not Significant"
    )
  ) %>%
  filter(!is.na(gene_symbol))

write.csv(res_myc_df, "tables/02_DGE_MYC_High_vs_Low.csv", row.names = FALSE)

message("MYC High vs Low – up: ",
        sum(res_myc_df$Significance == "MYC-High", na.rm = TRUE),
        "  down: ",
        sum(res_myc_df$Significance == "MYC-Low",  na.rm = TRUE))


################################################################################
## 09. HEATMAP – MYC-HIGH vs MYC-LOW  (top 30 DEGs, 600 DPI)
################################################################################

tumor_mask <- coldata$condition == "Tumor"

top30_myc_ens <- res_myc_df %>%
  filter(!is.na(padj), padj < 0.05) %>%
  arrange(desc(abs(log2FoldChange))) %>%
  slice_head(n = 30) %>%
  pull(ensembl)

expr_heat_myc <- norm_expr[top30_myc_ens, tumor_mask]
rownames(expr_heat_myc) <- res_myc_df$gene_symbol[
  match(top30_myc_ens, res_myc_df$ensembl)
]

z_myc <- t(scale(t(expr_heat_myc)))
z_myc[z_myc >  2] <-  2
z_myc[z_myc < -2] <- -2

myc_sample_order <- rownames(coldata)[tumor_mask][
  order(coldata$MYC_group[tumor_mask])
]
ann_col_myc <- data.frame(
  MYC_Status = coldata$MYC_group[tumor_mask],
  row.names  = rownames(coldata)[tumor_mask]
)[myc_sample_order, , drop = FALSE]

ann_colors_myc <- list(MYC_Status = c(MYC_low  = "#4575B4",
                                       MYC_high = "#D73027"))

tiff("figures/Fig_03_Heatmap_MYC_High_vs_Low.tiff",
     width = 8, height = 8, units = "in", res = 600, compression = "lzw")
pheatmap(
  z_myc[, myc_sample_order],
  annotation_col    = ann_col_myc,
  annotation_colors = ann_colors_myc,
  show_colnames     = FALSE,
  show_rownames     = TRUE,
  fontsize_row      = 8,
  cluster_cols      = FALSE,
  cluster_rows      = TRUE,
  color             = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
  breaks            = seq(-2, 2, length.out = 101),
  border_color      = NA,
  main              = "Top 30 DEGs: MYC-High vs MYC-Low (CRC Tumors)"
)
dev.off()


################################################################################
## 10. VOLCANO PLOTS  (600 DPI)
################################################################################

# Tumor vs Normal
p_vol_TN <- make_volcano(
  df         = res_TN_df %>% filter(!is.na(padj)),
  lfc_col    = "log2FoldChange",
  padj_col   = "padj",
  label_col  = "gene_symbol",
  color_col  = "Significance",
  color_pal  = c(Upregulated       = "#D73027",
                 Downregulated     = "#4575B4",
                 `Not Significant` = "grey70"),
  up_level   = "Upregulated",
  down_level = "Downregulated",
  title      = "Differential Expression: Tumor vs Normal (CRC)"
)
ggsave("figures/Fig_04a_Volcano_Tumor_vs_Normal.tiff",
       p_vol_TN, width = 7.5, height = 6, dpi = 600,
       device = "tiff", compression = "lzw")

# MYC High vs Low
p_vol_myc <- make_volcano(
  df         = res_myc_df %>% filter(!is.na(padj)),
  lfc_col    = "log2FoldChange",
  padj_col   = "padj",
  label_col  = "gene_symbol",
  color_col  = "Significance",
  color_pal  = c(`MYC-High`        = "#D73027",
                 `MYC-Low`         = "#4575B4",
                 `Not Significant` = "grey70"),
  up_level   = "MYC-High",
  down_level = "MYC-Low",
  title      = "Differential Expression: MYC-High vs MYC-Low (CRC)"
)
ggsave("figures/Fig_04b_Volcano_MYC_High_vs_Low.tiff",
       p_vol_myc, width = 7.5, height = 6, dpi = 600,
       device = "tiff", compression = "lzw")


################################################################################
## 11. GO ENRICHMENT – MYC-HIGH UPREGULATED GENES
################################################################################

myc_up_genes <- res_myc_df %>%
  filter(padj < 0.05, log2FoldChange > 1) %>%
  pull(gene_symbol)

ego_MYC <- enrichGO(
  gene          = myc_up_genes,
  OrgDb         = org.Hs.eg.db,
  keyType       = "SYMBOL",
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  readable      = TRUE
)

ego_MYC_simp <- simplify(ego_MYC, cutoff = 0.7,
                          by = "p.adjust", select_fun = min)

save_go_dotplot(
  ego_MYC_simp,
  filename = "figures/Fig_05_GO_MYC_High.tiff",
  title    = "GO Biological Processes – MYC-High CRC"
)

write.csv(as.data.frame(ego_MYC_simp),
          "tables/03_GO_MYC_High_simplified.csv",
          row.names = FALSE)


################################################################################
## 12. WGCNA – TUMOR SAMPLES
################################################################################

expr_tumor <- norm_expr[, tumor_mask]

# Top 5 000 most variable genes by MAD
mad_vals <- apply(expr_tumor, 1, mad)
top5k    <- names(sort(mad_vals, decreasing = TRUE))[1:5000]
datExpr  <- t(expr_tumor[top5k, ])

# Remove bad samples / genes (rarely needed but good practice)
gsg <- goodSamplesGenes(datExpr, verbose = 0)
if (!gsg$allOK) {
  datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes]
  message("WGCNA QC removed ",
          sum(!gsg$goodSamples), " samples, ",
          sum(!gsg$goodGenes),   " genes")
}

# Soft-threshold selection — examine Fig_06 before running blockwiseModules
powers <- 1:20
sft    <- pickSoftThreshold(datExpr, powerVector = powers, verbose = 0)

tiff("figures/Fig_06_WGCNA_SoftThreshold.tiff",
     width = 8, height = 4, units = "in", res = 600, compression = "lzw")
par(mfrow = c(1, 2), mar = c(4, 4, 2, 1))
# Scale-free fit
plot(sft$fitIndices[, 1],
     -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     xlab = "Soft Threshold (power)",
     ylab = expression(R^2 ~ "Scale-Free Fit"),
     type = "n", main = "Scale Independence")
text(sft$fitIndices[, 1],
     -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     labels = powers, col = "red")
abline(h = 0.85, col = "blue", lty = 2)
# Mean connectivity
plot(sft$fitIndices[, 1], sft$fitIndices[, 5],
     xlab = "Soft Threshold (power)",
     ylab = "Mean Connectivity",
     type = "n", main = "Mean Connectivity")
text(sft$fitIndices[, 1], sft$fitIndices[, 5],
     labels = powers, col = "red")
dev.off()

# ── Set softPower after inspecting Fig_06 (aim for R² ≥ 0.85) ────────────────
softPower <- 6

cor <- WGCNA::cor          # WGCNA::blockwiseModules requires its own cor()
net          <- blockwiseModules(
  datExpr,
  power             = softPower,
  TOMType           = "unsigned",
  minModuleSize     = 50,
  reassignThreshold = 0,
  mergeCutHeight    = 0.25,
  numericLabels     = TRUE,
  pamRespectsDendro = FALSE,
  saveTOMs          = FALSE,
  verbose           = 0
)
cor <- stats::cor          # restore immediately after

moduleColors <- labels2colors(net$colors)
message("WGCNA modules identified: ", length(unique(moduleColors)))

tiff("figures/Fig_07_WGCNA_Dendrogram.tiff",
     width = 10, height = 5, units = "in", res = 600, compression = "lzw")
plotDendroAndColors(
  net$dendrograms[[1]],
  moduleColors[net$blockGenes[[1]]],
  "Module colours",
  dendroLabels = FALSE, hang = 0.03,
  addGuide = TRUE, guideHang = 0.05
)
dev.off()


################################################################################
## 13. MYC-CORRELATED MODULE & ANNOTATED HUB GENES
################################################################################

# Align MYC expression vector to datExpr row order
# (rows may have been pruned by goodSamplesGenes)
myc_expr_tumor <- myc_exp[tumor_mask][rownames(datExpr)]

MEs            <- moduleEigengenes(datExpr, moduleColors)$eigengenes
moduleTraitCor <- WGCNA::cor(MEs, myc_expr_tumor, use = "p")
moduleTraitP   <- corPvalueStudent(moduleTraitCor, nrow(datExpr))

# Module–MYC labelled heatmap
tiff("figures/Fig_08_WGCNA_Module_MYC_Heatmap.tiff",
     width = 4, height = 7, units = "in", res = 600, compression = "lzw")
labeledHeatmap(
  Matrix      = matrix(moduleTraitCor, ncol = 1),
  xLabels     = "MYC expression",
  yLabels     = rownames(moduleTraitCor),
  ySymbols    = rownames(moduleTraitCor),
  colorLabels = FALSE,
  colors      = blueWhiteRed(50),
  textMatrix  = signif(matrix(moduleTraitCor, ncol = 1), 2),
  main        = "Module–MYC Correlation"
)
dev.off()

# Best module (highest positive correlation with MYC)
module_of_interest <- rownames(moduleTraitCor)[which.max(moduleTraitCor[, 1])]
target_color       <- gsub("^ME", "", module_of_interest)
ME4_genes          <- colnames(datExpr)[moduleColors == target_color]

message("MYC-correlated module : ", module_of_interest,
        "  colour = ", target_color,
        "  genes  = ", length(ME4_genes))

# Hub gene scoring
kME              <- signedKME(datExpr, MEs)
geneSignificance <- cor(datExpr, myc_expr_tumor, use = "p")  # samples × genes

ME4_col <- grep(paste0("kME", target_color), colnames(kME), value = TRUE)
stopifnot("ME4 kME column not found" = length(ME4_col) == 1)

hub_ME4 <- data.frame(
  Gene   = ME4_genes,
  kME    = kME[ME4_genes, ME4_col],
  GS_MYC = geneSignificance[ME4_genes, 1]
)

stopifnot(
  nrow(hub_ME4) == length(ME4_genes),
  !any(is.na(hub_ME4$kME)),
  !any(is.na(hub_ME4$GS_MYC))
)

# Annotate: strip Ensembl version suffix, map to HGNC symbol
hub_ME4$ensembl_clean <- sub("\\..*$", "", hub_ME4$Gene)
hub_ME4$GeneSymbol    <- mapIds(
  org.Hs.eg.db,
  keys      = hub_ME4$ensembl_clean,
  keytype   = "ENSEMBL",
  column    = "SYMBOL",
  multiVals = "first"
)
hub_ME4 <- hub_ME4 %>%
  filter(!is.na(GeneSymbol)) %>%
  arrange(desc(abs(kME)))

message("Annotated hub genes: ", nrow(hub_ME4))
message("Top 15:")
print(head(hub_ME4[, c("GeneSymbol", "kME", "GS_MYC")], 15))

write.csv(data.frame(ensembl = ME4_genes),
          "tables/04_ME4_gene_list.csv",
          row.names = FALSE)
write.csv(hub_ME4,
          "tables/05_ME4_hub_genes_annotated.csv",
          row.names = FALSE)

# Eigengene score for every sample (used in survival section 16)
ME4_eigengene        <- MEs[, module_of_interest]
names(ME4_eigengene) <- rownames(datExpr)


################################################################################
## 14. GO ENRICHMENT – ME4 MODULE HUB GENES
################################################################################

ego_ME4 <- enrichGO(
  gene          = hub_ME4$GeneSymbol,
  OrgDb         = org.Hs.eg.db,
  keyType       = "SYMBOL",
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  readable      = TRUE
)

ego_ME4_simp <- simplify(ego_ME4, cutoff = 0.7,
                          by = "p.adjust", select_fun = min)

save_go_dotplot(
  ego_ME4_simp,
  filename = "figures/Fig_09_GO_ME4_Module.tiff",
  title    = "GO Biological Processes – MYC-Associated ME4 Module"
)

write.csv(as.data.frame(ego_ME4_simp),
          "tables/06_GO_ME4_module_simplified.csv",
          row.names = FALSE)


################################################################################
## 15. GSEA – MYC + EXTENDED HALLMARK GENE SETS
################################################################################

# Ranked gene list: LFC-ranked, symbol-keyed, no duplicates
gene_list <- res_myc_df %>%
  filter(!is.na(log2FoldChange)) %>%
  distinct(gene_symbol, .keep_all = TRUE) %>%
  arrange(desc(log2FoldChange)) %>%
  pull(log2FoldChange, name = gene_symbol)

msig_H <- msigdbr(species = "Homo sapiens", category = "H")

# MYC hallmark sets (V1 + V2)
myc_sets <- msig_H %>%
  filter(gs_name %in% c("HALLMARK_MYC_TARGETS_V1",
                        "HALLMARK_MYC_TARGETS_V2")) %>%
  select(gs_name, gene_symbol)

gsea_myc <- GSEA(geneList     = gene_list,
                 TERM2GENE    = myc_sets,
                 pvalueCutoff = 0.05,
                 verbose      = FALSE)

# Extended sets: MYC + related proliferation/repair pathways
hallmark_ext <- msig_H %>%
  filter(gs_name %in% c("HALLMARK_MYC_TARGETS_V1",
                        "HALLMARK_MYC_TARGETS_V2",
                        "HALLMARK_E2F_TARGETS",
                        "HALLMARK_G2M_CHECKPOINT",
                        "HALLMARK_DNA_REPAIR")) %>%
  select(gs_name, gene_symbol)

gsea_hallmark <- GSEA(geneList     = gene_list,
                      TERM2GENE    = hallmark_ext,
                      pvalueCutoff = 0.05,
                      verbose      = FALSE)

# Per-geneset enrichment plots for MYC V1 and V2
for (gs in c("HALLMARK_MYC_TARGETS_V1", "HALLMARK_MYC_TARGETS_V2")) {
  if (gs %in% gsea_myc@result$ID) {
    tiff(paste0("figures/Fig_10_GSEA_", gs, ".tiff"),
         width = 7, height = 5, units = "in", res = 600, compression = "lzw")
    print(gseaplot2(gsea_myc, geneSetID = gs,
                   title = paste("GSEA:", gsub("_", " ", gs))))
    dev.off()
  }
}

# Summary dotplot for extended sets
p_gsea_dot <- dotplot(gsea_hallmark, showCategory = 5) +
  coord_cartesian(clip = "off") +
  theme(plot.margin = margin(10, 30, 10, 10)) +
  ggtitle("GSEA – MYC-Associated Hallmark Pathways (CRC)")

ggsave("figures/Fig_11_GSEA_Hallmark_Summary.tiff",
       p_gsea_dot, dpi = 600, width = 7, height = 5,
       device = "tiff", compression = "lzw")

write.csv(as.data.frame(gsea_hallmark),
          "tables/07_GSEA_Hallmark_results.csv",
          row.names = FALSE)


################################################################################
## 16. SURVIVAL – ME4 SCORE KM + CONTINUOUS & STAGE-ADJUSTED COX
################################################################################

# Build clinical_surv: one row per tumor sample with ME4 eigengene score,
# survival endpoints, MYC group, and stage — all drawn from coldata.
# ME4_eigengene is keyed by datExpr rownames (tumor sample barcodes).
tumor_barcodes <- rownames(coldata)[tumor_mask]

clinical_surv <- coldata[tumor_mask, ] %>%
  mutate(
    patient_id = rownames(.),
    ME4_score  = ME4_eigengene[patient_id],          # NA if sample was pruned by WGCNA QC
    OS_time    = ifelse(!is.na(days_to_death),
                        days_to_death,
                        days_to_last_follow_up),
    OS_event   = ifelse(vital_status == "Dead", 1L, 0L)
  ) %>%
  filter(!is.na(OS_time), OS_time > 0)

message("Survival cohort: ", nrow(clinical_surv),
        " tumor samples  |  events: ", sum(clinical_surv$OS_event, na.rm = TRUE))

clinical_tumor_surv <- clinical_surv %>%
  filter(!is.na(ME4_score))

stopifnot(
  "Insufficient tumor patients for ME4 survival analysis" =
    nrow(clinical_tumor_surv) > 10
)

# Median-split KM
clinical_tumor_surv <- clinical_tumor_surv %>%
  mutate(ME4_group = factor(
    ifelse(ME4_score >= median(ME4_score, na.rm = TRUE),
           "ME4_high", "ME4_low"),
    levels = c("ME4_low", "ME4_high")
  ))

fit_ME4 <- survfit(Surv(OS_time, OS_event) ~ ME4_group,
                   data = clinical_tumor_surv)

save_km(
  fit          = fit_ME4,
  data         = clinical_tumor_surv,
  filename     = "figures/Fig_12_KM_ME4_Group.tiff",
  title        = "Overall Survival – ME4 High vs Low (CRC)",
  legend_title = "ME4 Score",
  legend_labs  = c("ME4 Low", "ME4 High"),
  palette      = c("#4575B4", "#D73027")
)

# Continuous Cox
cox_ME4_cont <- coxph(Surv(OS_time, OS_event) ~ ME4_score,
                      data = clinical_tumor_surv)
message("── Cox PH  (continuous ME4 score) ──")
print(summary(cox_ME4_cont))

# Stage-adjusted Cox
clinical_tumor_surv <- clinical_tumor_surv %>%
  mutate(stage_simple = factor(substr(ajcc_pathologic_stage, 1, 7)))

df_stage <- clinical_tumor_surv %>%
  filter(!is.na(stage_simple), stage_simple != "")

if (nrow(df_stage) > 20 && nlevels(droplevels(df_stage$stage_simple)) > 1) {
  cox_ME4_stage <- coxph(
    Surv(OS_time, OS_event) ~ ME4_score + stage_simple,
    data = df_stage
  )
  message("── Cox PH  (ME4 score + AJCC stage) ──")
  print(summary(cox_ME4_stage))
} else {
  warning("Insufficient or single-level stage data – skipping adjusted Cox.")
}

write.csv(
  clinical_tumor_surv %>%
    select(patient_id, MYC_group, ME4_score, ME4_group,
           OS_time, OS_event, stage_simple),
  "tables/08_Survival_ME4_clinical.csv",
  row.names = FALSE
)

################################################################################
## 17. ML DATASET – ME4 GENES + SURVIVAL (FULL COHORT + LEAKAGE-FREE SPLIT)
#
#  Design:
#   A. FULL COHORT export  →  ML_dataset_ME4_GeneSymbols_Survival.csv
#      All tumor samples, ME4 genes (from section 13 WGCNA), HGNC symbols,
#      MYC_group, MYC_label, OS_time, OS_event.
#      This is the file read directly by the Python XGBoost / SHAP script.
#
#   B. LEAKAGE-FREE SPLIT  →  09_ML_train_ME4.csv / 10_ML_test_ME4.csv
#      Stratified 70/30 split by MYC group; WGCNA run on training set only
#      so no test-set information leaks into feature selection.
#      Use these for rigorous cross-validation experiments.
################################################################################

# ── Helper: build a samples-as-rows ML dataframe ──────────────────────────────
build_ml_df <- function(expr_mat, sample_idx) {
  df           <- as.data.frame(t(expr_mat))
  df$MYC_group <- as.character(coldata$MYC_group[sample_idx])
  df$MYC_label <- ifelse(df$MYC_group == "MYC_high", 1L, 0L)
  df$OS_time   <- ifelse(!is.na(coldata$days_to_death[sample_idx]),
                          coldata$days_to_death[sample_idx],
                          coldata$days_to_last_follow_up[sample_idx])
  df$OS_event  <- ifelse(coldata$vital_status[sample_idx] == "Dead", 1L, 0L)
  df
}

# ── A. FULL COHORT EXPORT ─────────────────────────────────────────────────────
# ME4_genes comes from section 13 (full-cohort WGCNA on all tumor samples).
# Convert those Ensembl IDs to HGNC symbols for Python-friendly column names.

tumor_idx_int <- which(coldata$condition == "Tumor")

expr_full_ME4 <- ensembl_to_symbols(norm_expr[ME4_genes, tumor_idx_int])

# Drop any samples with missing survival before export
full_df <- build_ml_df(expr_full_ME4, tumor_idx_int)
full_df <- full_df[!is.na(full_df$OS_time) & full_df$OS_time > 0, ]

message("Full ME4 ML dataset: ", nrow(full_df), " samples x ",
        ncol(full_df) - 4, " genes  (+MYC_group, MYC_label, OS_time, OS_event)")

# This filename matches exactly what the Python script reads
write.csv(full_df,
          "ML_dataset_ME4_GeneSymbols_Survival.csv",
          row.names = TRUE)

message("Exported: ML_dataset_ME4_GeneSymbols_Survival.csv")

# ── B. LEAKAGE-FREE STRATIFIED 70/30 SPLIT ───────────────────────────────────
set.seed(123)
train_idx <- unlist(lapply(
  split(tumor_idx_int, coldata$MYC_group[tumor_idx_int]),
  function(idx) sample(idx, floor(0.7 * length(idx)))
))
test_idx <- setdiff(tumor_idx_int, train_idx)

message("ML split – train: ", length(train_idx),
        "  test: ",           length(test_idx))

# WGCNA on training set only (no test information used in feature selection)
expr_train_raw <- norm_expr[, train_idx]
mad_train      <- apply(expr_train_raw, 1, mad)
top5k_train    <- names(sort(mad_train, decreasing = TRUE))[1:5000]
datExpr_train  <- t(expr_train_raw[top5k_train, ])

cor <- WGCNA::cor          # WGCNA::blockwiseModules requires its own cor()
net_train <- blockwiseModules(
  datExpr_train,
  power             = softPower,
  TOMType           = "unsigned",
  minModuleSize     = 50,
  mergeCutHeight    = 0.25,
  numericLabels     = TRUE,
  pamRespectsDendro = FALSE,
  saveTOMs          = FALSE,
  verbose           = 0
)
cor <- stats::cor          # restore immediately after

moduleColors_train <- labels2colors(net_train$colors)
MEs_train          <- moduleEigengenes(datExpr_train,
                                       moduleColors_train)$eigengenes
MYC_num_train      <- as.numeric(coldata$MYC_group[train_idx] == "MYC_high")

modCor_train    <- WGCNA::cor(MEs_train, MYC_num_train, use = "p")
mod_train       <- rownames(modCor_train)[which.max(abs(modCor_train[, 1]))]
color_train     <- gsub("^ME", "", mod_train)
ME4_genes_train <- colnames(datExpr_train)[moduleColors_train == color_train]

message("Train ME4 module (", mod_train, "): ",
        length(ME4_genes_train), " genes")

# Apply training gene list to both splits; convert to HGNC symbols
expr_train_ME4 <- ensembl_to_symbols(norm_expr[ME4_genes_train, train_idx])
expr_test_ME4  <- ensembl_to_symbols(norm_expr[ME4_genes_train, test_idx])

# Align to the same symbol set across both splits
common_syms    <- intersect(rownames(expr_train_ME4), rownames(expr_test_ME4))
expr_train_ME4 <- expr_train_ME4[common_syms, ]
expr_test_ME4  <- expr_test_ME4[common_syms, ]

message("Shared gene symbols after annotation: ", length(common_syms))

train_df <- build_ml_df(expr_train_ME4, train_idx)
test_df  <- build_ml_df(expr_test_ME4,  test_idx)

message("Train shape: ", nrow(train_df), " samples x ",
        ncol(train_df) - 4, " genes  (+MYC_group, MYC_label, OS_time, OS_event)")
message("Test  shape: ", nrow(test_df),  " samples x ",
        ncol(test_df)  - 4, " genes")

write.csv(train_df, "tables/09_ML_train_ME4.csv",  row.names = TRUE)
write.csv(test_df,  "tables/10_ML_test_ME4.csv",   row.names = TRUE)

# Figure 13: class balance barplot
bal_df <- bind_rows(
  data.frame(Split = "Train", MYC_group = train_df$MYC_group),
  data.frame(Split = "Test",  MYC_group = test_df$MYC_group)
)
p_bal <- ggplot(bal_df, aes(Split, fill = MYC_group)) +
  geom_bar(position = "fill", width = 0.5) +
  scale_fill_manual(values = c(MYC_low = "#4575B4", MYC_high = "#D73027")) +
  scale_y_continuous(labels = percent_format()) +
  labs(title = "ML Split – Class Balance (ME4 Genes)",
       y = "Proportion", x = NULL, fill = "MYC Group") +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

ggsave("figures/Fig_13_ML_ClassBalance.tiff",
       p_bal, dpi = 600, width = 5, height = 4,
       device = "tiff", compression = "lzw")


################################################################################
## 18. ceRNA NETWORK – miRTarBase + ENCORI
#
#  Place these files in your working directory before running:
#    miRTarBase_MTI.csv       https://mirtarbase.cuhk.edu.cn
#    ENCORI_miRNA_lncRNA.tsv  https://rnasysu.com/encori
#
#  Output tables are ready for direct import into Cytoscape.
################################################################################

MIRTARBASE_FILE <- "miRTarBase_MTI.csv"
ENCORI_FILE     <- "ENCORI_miRNA_lncRNA.tsv"

if (!file.exists(MIRTARBASE_FILE)) {
  warning("miRTarBase file not found – skipping ceRNA network.\n",
          "Download: https://mirtarbase.cuhk.edu.cn")
} else if (!file.exists(ENCORI_FILE)) {
  warning("ENCORI file not found – skipping ceRNA network.\n",
          "Download: https://rnasysu.com/encori")
} else {

  # Seed mRNAs: top 40 MYC-high upregulated genes
  myc_mRNAs <- res_myc_df %>%
    filter(padj < 0.05, log2FoldChange > 1, !is.na(gene_symbol)) %>%
    distinct(gene_symbol) %>%
    slice_head(n = 40) %>%
    pull(gene_symbol)

  # miRNA–mRNA: human miRNAs with functional MTI evidence
  mirtar   <- read.csv(MIRTARBASE_FILE, stringsAsFactors = FALSE)
  mir_mrna <- mirtar %>%
    dplyr::filter(grepl("^hsa", miRNA),
                  Target.Gene %in% myc_mRNAs,
                  Support.Type == "Functional MTI") %>%
    dplyr::select(miRNA, mRNA = Target.Gene) %>%
    dplyr::distinct()

  message("miRNA–mRNA interactions: ", nrow(mir_mrna))

  # lncRNA–miRNA: lncRNAs sponging the selected miRNAs
  mir_lnc_raw   <- fread(ENCORI_FILE, data.table = FALSE)
  mir_lnc_clean <- mir_lnc_raw %>%
    dplyr::filter(geneType == "lncRNA",
                  miRNAname %in% mir_mrna$miRNA) %>%
    dplyr::select(lncRNA = geneName, miRNA = miRNAname) %>%
    dplyr::distinct()

  message("lncRNA–miRNA interactions: ", nrow(mir_lnc_clean))

  # ceRNA triplets: lncRNA --|miRNA|--> mRNA
  ceRNA_triplets <- mir_lnc_clean %>%
    inner_join(mir_mrna, by = "miRNA") %>%
    distinct()

  message("ceRNA triplets: ", nrow(ceRNA_triplets))

  # Cytoscape-ready edge and node tables
  edges <- bind_rows(
    ceRNA_triplets %>%
      transmute(source = lncRNA, target = miRNA,
                interaction = "lncRNA-miRNA"),
    ceRNA_triplets %>%
      transmute(source = miRNA,  target = mRNA,
                interaction = "miRNA-mRNA")
  )

  nodes <- data.frame(
    id = unique(c(edges$source, edges$target)),
    stringsAsFactors = FALSE
  ) %>%
    mutate(type = case_when(
      grepl("^hsa",                              id) ~ "miRNA",
      grepl("^LINC|^MALAT|^NEAT|^HOTAIR|^OIP5", id) ~ "lncRNA",
      TRUE                                           ~ "mRNA"
    ))

  write.csv(edges, "tables/11_ceRNA_network_edges.csv", row.names = FALSE)
  write.csv(nodes, "tables/12_ceRNA_network_nodes.csv", row.names = FALSE)

  message("ceRNA network exported to tables/ (",
          nrow(edges), " edges, ", nrow(nodes), " nodes)")
}


################################################################################
message("── All outputs written to figures/ and tables/ ──")
message("── Pipeline complete ──")
## END OF PIPELINE
################################################################################
