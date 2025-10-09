# ============================================================
# BEESWARM PLOT SCRIPT — EXPLAINED -- xx KH
# TRPV1 RNA-seq dataset: visualization of top 4 DEGs (up/down)
# Notes:
#  • Object names, filenames my own, logic same as in Kls scripts
#  • Code reorganized and commented for clarity to best of my knowledge - no guarantee.
#  • Each contrast processed sequentially with identical logic.
#  • Visual output: beeswarm points = sample-level log2 expression per group; semi-transparent boxplots overlay group medians/IQR.
# ============================================================


# ============================================================
# 1) PACKAGES
# ============================================================
library(beeswarm)


# ============================================================
# 2) INTRO / WORKING DIRECTORY
# ============================================================
# Working directory -- adjust to your system !! this is absolute path (for now)
setwd("/Users/katharinahollmann/GitHubDesktop/Wagner_comparative_analyses/TRPV1/data")
dir()  # verify available files; choose the newest normalized (batch-corrected) file below


# ============================================================
# 3) LOAD DATA FILES
#    dat: normalized expression table with gene annotations + sample columns
#    ta:  target/metadata table with sample_ID and group
# ============================================================
dat <- read.table("norm_Wagner_CLP_TRPV1_2023_080525.txt", header = TRUE)
names(dat); dim(dat)
dat[1:5, c(2,16:18)]  # preview a few genes across selected samples

ta <- read.table("target_Wagner_CLP_TRPV1_2023_240124.txt")
head(ta); unique(ta$group)


# ============================================================
# 4) PLOT PARAMETERS + REUSABLE FUNCTION
#    Function plot_gene(gene_name):
#      • Extracts one gene’s log2 expression across samples
#      • Groups samples by ta$group (fixed level order)
#      • Draws beeswarm + overlayed boxplot
#      • Performs one-way ANOVA (u ~ f) and reports p-value in console
# ============================================================
par(mfrow = c(2,2), font = 2, font.axis = 2, font.lab = 3, mar = c(8,5,2,2))  # default 2x2 layout

plot_gene <- function(gene_name) {
  # target gene
  query <- gene_name
  gene_row <- dat[dat$gene_symbol %in% query, ]
  selected_row <- as.numeric(rownames(gene_row))  # row index of the selected gene
  
  # expression matrix = only sample columns (9:ncol)
  test1 <- as.matrix(dat[, 9:ncol(dat)])
  
  # load and check target/metadata
  target_data <- read.table("target_Wagner_CLP_TRPV1_2023_240124.txt", header = TRUE, na.strings = c("", "na"))
  if (!identical(colnames(test1), as.character(target_data[, 1]))) {
    stop("Column names do not match between data and target files!")
  }
  identical(colnames(test1), as.character(target_data[, 1]))  # should be TRUE
  
  # group factor in desired plotting order
  f <- factor(target_data$group, levels = c("WT_CLP", "TRPV1_KO_CLP", "TRPV1_KO_Sham", "WT_Sham"))
  labels <- unique(levels(f))
  
  # single-gene vector across samples
  u <- test1[selected_row, ]
  
  # point color setup (all black by default)
  all <- target_data[, 1]
  all.2 <- rep("black", length(all))
  cc1 <- all.2
  
  # basic stats in console for QC
  tapply(u, f, mean)
  aov1 <- summary(aov(u ~ f))
  pval1 <- aov1[[1]][1, 5]
  pval2 <- round(pval1, digits = 6)
  pval3 <- paste("ANOVA:", pval2)
  
  # beeswarm + semi-transparent boxplot overlay
  beeswarm(u ~ f,
           data   = test1,
           corral = "random",    # jitter method to reduce overlap
           method = "center",    # center per group
           pch    = 16,
           pwcol  = cc1,
           xlab   = "",
           ylab   = "log2 expression",
           labels = labels,
           main   = paste(dat[selected_row, c(2)]),  # gene symbol as title
           las    = 2)
  abline(v = 2.5, lty = 5, col = "black", lwd = 3)
  boxplot(u ~ f, data = test1, names = c("", "", "", ""), col = "#0000ff22", add = TRUE, yaxt = "n")
  # legend can be enabled if desired:
  # legend("topleft", legend = pval3, pch = 16, col = "white", bty = "n", cex = 0.5)
}

# Quick search helper (optional)
test.query <- grep("Il17", dat$gene_symbol, value = TRUE); test.query; length(test.query)


# ============================================================
# 5) DISCOVER DESEQ2 RESULT FILES FOR TRPV1
#    all_deseq: list of 5 contrast files
# ============================================================
setwd("/Users/katharinahollmann/GitHubDesktop/Projects/Wagner_comparative_analyses/TRPV1/data")
dir()
all_data  <- dir(); all_data
all_deseq <- grep("DESeq2_TRPV1", all_data, value = TRUE); all_deseq

# (Optional) remove specific files if needed:
# to_remove <- c("pathways_DEG_DESeq2_comb_sita_CLP_vs_Sham_260325.txt",
#                "pathways_down_DEG_DESeq2_comb_sita_CLP_vs_Sham_260325.txt")
# all_deseq <- setdiff(all_deseq, to_remove)


# ============================================================
# 6) CONTRAST 1 — KO CLP vs KO Sham
#    Visual output: 2 PDFs
#      • Top four upregulated genes KO CLP vs KO Sham
#      • Top four downregulated genes KO CLP vs KO Sham
# ============================================================

# ---- Load DE table and filter by thresholds ----
comp_file <- all_deseq[1]; comp_file
title <- gsub("_010525.txt", "", comp_file); title

comp_file_deg <- read.table(comp_file, header = TRUE)
names(comp_file_deg); head(comp_file_deg)

# gene_symbol QC
if (any(is.na(comp_file_deg$gene_symbol))) cat("Warning: NA in gene_symbol.\n") else cat("No NA in gene_symbol.\n")
if (any(duplicated(comp_file_deg$gene_symbol))) comp_file_deg$gene_symbol <- make.unique(comp_file_deg$gene_symbol) else cat("No duplicates in gene_symbol.\n")

lim.P  <- 0.05
lim.FC <- log2(2)
filtered_deg <- subset(comp_file_deg, comp_file_deg$padj < lim.P & abs(comp_file_deg$log2FoldChange) > lim.FC)
dim(comp_file_deg); dim(filtered_deg)

# ---- UP: choose top 4 by log2FC and plot with plot_gene() ----
sorted_deg <- filtered_deg[order(filtered_deg$log2FoldChange, decreasing = TRUE), ]
names(sorted_deg); sorted_deg[1:20, c(1,2,28,32)]
all_genes <- sorted_deg$gene_symbol
genes_for_plot <- as.character(all_genes[1:4]); genes_for_plot; length(genes_for_plot)

# PDF device name uses comp_file variable as in your script (kept)
comp_file
comp_file <- "Top4_upregulated_genesKOCLPvsKOSham.pdf"
pdf(comp_file, width = 10, height = 7)
par(mfrow = c(2,2), mar = c(8,5,2,2), oma = c(0,0,4,0))
for (i in 1:length(genes_for_plot)) {
  gen <- genes_for_plot[i]
  filtered_data <- filtered_deg[filtered_deg$gene_symbol == gen, ]
  if (nrow(filtered_data) == 0) { message("Skipping gene ", gen, " (no data found)"); next }
  if ("group" %in% colnames(filtered_data)) {
    if (length(unique(filtered_data$group)) < 2) { message("Skipping gene ", gen, " (not enough factor levels)"); next }
  }
  print(paste("Plotting gene:", gen))
  plot_gene(gen)
}
mtext("Top four upregulated genes KO CLP vs KO Sham", outer = TRUE, line = 1, cex = 1, font = 2)
dev.off(); cat("Saved PDF:", comp_file, "\n")

# ---- DOWN: choose top 4 by lowest log2FC and plot ----
sorted_deg <- filtered_deg[order(filtered_deg$log2FoldChange, decreasing = FALSE), ]
names(sorted_deg); sorted_deg[1:20, c(1,2,28,32)]
all_genes <- sorted_deg$gene_symbol
genes_for_plot <- as.character(all_genes[1:4]); genes_for_plot

comp_file
comp_file <- "Top4_downregulated_genesKOCLPvsKOSham.pdf"
pdf(comp_file, width = 10, height = 7)
par(mfrow = c(2,2), mar = c(8,5,2,2), oma = c(0,0,4,0))
for (i in 1:length(genes_for_plot)) {
  gen <- genes_for_plot[i]
  plot_gene(gen)
}
mtext("Top four downregulated genes KO CLP vs KO Sham", outer = TRUE, line = 1, cex = 1, font = 2)
dev.off(); cat("Saved PDF:", comp_file, "\n")


# ============================================================
# 7) CONTRAST 2 — KO CLP vs WT CLP
#    Visual output: 2 PDFs
# ============================================================

all_deseq
comp_file <- all_deseq[2]; comp_file
title <- gsub("_010525.txt", "", comp_file); title

comp_file_deg <- read.table(comp_file, header = TRUE)
names(comp_file_deg); head(comp_file_deg)
if (any(is.na(comp_file_deg$gene_symbol))) cat("Warning: NA in gene_symbol.\n") else cat("No NA in gene_symbol.\n")
if (any(duplicated(comp_file_deg$gene_symbol))) comp_file_deg$gene_symbol <- make.unique(comp_file_deg$gene_symbol) else cat("No duplicates in gene_symbol.\n")

lim.P  <- 0.05
lim.FC <- log2(2)
filtered_deg <- subset(comp_file_deg, comp_file_deg$padj < lim.P & abs(comp_file_deg$log2FoldChange) > lim.FC)
dim(comp_file_deg); dim(filtered_deg)

# UP
sorted_deg <- filtered_deg[order(filtered_deg$log2FoldChange, decreasing = TRUE), ]
names(sorted_deg); sorted_deg[1:20, c(1,2,28,32)]
all_genes <- sorted_deg$gene_symbol; genes_for_plot <- as.character(all_genes[1:4]); genes_for_plot; length(genes_for_plot)

comp_file
comp_file <- "Top4_upregulated_genesKOCLPvsWTCLP.pdf"
pdf(comp_file, width = 10, height = 7)
par(mfrow = c(2,2), mar = c(8,5,2,2), oma = c(0,0,4,0))
for (i in 1:length(genes_for_plot)) {
  gen <- genes_for_plot[i]
  filtered_data <- filtered_deg[filtered_deg$gene_symbol == gen, ]
  if (nrow(filtered_data) == 0) { message("Skipping gene ", gen); next }
  if ("group" %in% colnames(filtered_data)) {
    if (length(unique(filtered_data$group)) < 2) { message("Skipping gene ", gen); next }
  }
  print(paste("Plotting gene:", gen))
  plot_gene(gen)
}
mtext("Top four upregulated genes KO CLP vs WT CLP", outer = TRUE, line = 1, cex = 1, font = 2)
dev.off(); cat("Saved PDF:", comp_file, "\n")

# DOWN
sorted_deg <- filtered_deg[order(filtered_deg$log2FoldChange, decreasing = FALSE), ]
names(sorted_deg); sorted_deg[1:20, c(1,2,28,32)]
all_genes <- sorted_deg$gene_symbol; genes_for_plot <- as.character(all_genes[1:4]); genes_for_plot

comp_file
comp_file <- "Top4_downregulated_genesKOCLPvsWTCLP.pdf"
pdf(comp_file, width = 10, height = 7)
par(mfrow = c(2,2), mar = c(8,5,2,2), oma = c(0,0,4,0))
for (i in 1:length(genes_for_plot)) {
  gen <- genes_for_plot[i]
  plot_gene(gen)
}
mtext("Top four downregulated genes KO CLP vs WT CLP", outer = TRUE, line = 1, cex = 1, font = 2)
dev.off(); cat("Saved PDF:", comp_file, "\n")


# ============================================================
# 8) CONTRAST 3 — KO CLP vs WT Sham
#    Visual output: 2 PDFs
# ============================================================

all_deseq
comp_file <- all_deseq[3]; comp_file
title <- gsub("_010525.txt", "", comp_file); title

comp_file_deg <- read.table(comp_file, header = TRUE)
names(comp_file_deg); head(comp_file_deg)
if (any(is.na(comp_file_deg$gene_symbol))) cat("Warning: NA in gene_symbol.\n") else cat("No NA in gene_symbol.\n")
if (any(duplicated(comp_file_deg$gene_symbol))) comp_file_deg$gene_symbol <- make.unique(comp_file_deg$gene_symbol) else cat("No duplicates in gene_symbol.\n")

lim.P  <- 0.05
lim.FC <- log2(2)
filtered_deg <- subset(comp_file_deg, comp_file_deg$padj < lim.P & abs(comp_file_deg$log2FoldChange) > lim.FC)
dim(comp_file_deg); dim(filtered_deg)

# UP
sorted_deg <- filtered_deg[order(filtered_deg$log2FoldChange, decreasing = TRUE), ]
names(sorted_deg); sorted_deg[1:20, c(1,2,28,32)]
all_genes <- sorted_deg$gene_symbol; genes_for_plot <- as.character(all_genes[1:4]); genes_for_plot; length(genes_for_plot)

comp_file
comp_file <- "Top4_upregulated_genesKOCLPvsWTSham.pdf"
pdf(comp_file, width = 10, height = 7)
par(mfrow = c(2,2), mar = c(8,5,2,2), oma = c(0,0,4,0))
for (i in 1:length(genes_for_plot)) {
  gen <- genes_for_plot[i]
  filtered_data <- filtered_deg[filtered_deg$gene_symbol == gen, ]
  if (nrow(filtered_data) == 0) { message("Skipping gene ", gen); next }
  if ("group" %in% colnames(filtered_data)) {
    if (length(unique(filtered_data$group)) < 2) { message("Skipping gene ", gen); next }
  }
  print(paste("Plotting gene:", gen))
  plot_gene(gen)
}
mtext("Top four upregulated genes KO CLP vs WT Sham", outer = TRUE, line = 1, cex = 1, font = 2)
dev.off(); cat("Saved PDF:", comp_file, "\n")

# DOWN
sorted_deg <- filtered_deg[order(filtered_deg$log2FoldChange, decreasing = FALSE), ]
names(sorted_deg); sorted_deg[1:20, c(1,2,28,32)]
all_genes <- sorted_deg$gene_symbol; genes_for_plot <- as.character(all_genes[1:4]); genes_for_plot

comp_file
comp_file <- "Top4_downregulated_genesKOCLPvsWTSham.pdf"
pdf(comp_file, width = 10, height = 7)
par(mfrow = c(2,2), mar = c(8,5,2,2), oma = c(0,0,4,0))
for (i in 1:length(genes_for_plot)) {
  gen <- genes_for_plot[i]
  plot_gene(gen)
}
mtext("Top four downregulated genes KO CLP vs WT Sham", outer = TRUE, line = 1, cex = 1, font = 2)
dev.off(); cat("Saved PDF:", comp_file, "\n")


# ============================================================
# 9) CONTRAST 4 — KO Sham vs WT Sham
#    Visual output: 2 PDFs
# ============================================================

all_deseq
comp_file <- all_deseq[4]; comp_file
title <- gsub("_010525.txt", "", comp_file); title

comp_file_deg <- read.table(comp_file, header = TRUE)
names(comp_file_deg); head(comp_file_deg)
if (any(is.na(comp_file_deg$gene_symbol))) cat("Warning: NA in gene_symbol.\n") else cat("No NA in gene_symbol.\n")
if (any(duplicated(comp_file_deg$gene_symbol))) comp_file_deg$gene_symbol <- make.unique(comp_file_deg$gene_symbol) else cat("No duplicates in gene_symbol.\n")

lim.P  <- 0.05
lim.FC <- log2(2)
filtered_deg <- subset(comp_file_deg, comp_file_deg$padj < lim.P & abs(comp_file_deg$log2FoldChange) > lim.FC)
dim(comp_file_deg); dim(filtered_deg)

# UP
sorted_deg <- filtered_deg[order(filtered_deg$log2FoldChange, decreasing = TRUE), ]
names(sorted_deg); sorted_deg[1:20, c(1,2,28,32)]
all_genes <- sorted_deg$gene_symbol; genes_for_plot <- as.character(all_genes[1:4]); genes_for_plot; length(genes_for_plot)

comp_file
comp_file <- "Top4_upregulated_genesKOShamvsWTSham.pdf"
pdf(comp_file, width = 10, height = 7)
par(mfrow = c(2,2), mar = c(8,5,2,2), oma = c(0,0,4,0))
for (i in 1:length(genes_for_plot)) {
  gen <- genes_for_plot[i]
  filtered_data <- filtered_deg[filtered_deg$gene_symbol == gen, ]
  if (nrow(filtered_data) == 0) { message("Skipping gene ", gen); next }
  if ("group" %in% colnames(filtered_data)) {
    if (length(unique(filtered_data$group)) < 2) { message("Skipping gene ", gen); next }
  }
  print(paste("Plotting gene:", gen))
  plot_gene(gen)
}
mtext("Top four upregulated genes KO Sham vs WT Sham", outer = TRUE, line = 1, cex = 1, font = 2)
dev.off(); cat("Saved PDF:", comp_file, "\n")

# DOWN
sorted_deg <- filtered_deg[order(filtered_deg$log2FoldChange, decreasing = FALSE), ]
names(sorted_deg); sorted_deg[1:20, c(1,2,28,32)]
all_genes <- sorted_deg$gene_symbol; genes_for_plot <- as.character(all_genes[1:4]); genes_for_plot

comp_file
comp_file <- "Top4_downregulated_genesKOShamvsWTSham.pdf"
pdf(comp_file, width = 10, height = 7)
par(mfrow = c(2,2), mar = c(8,5,2,2), oma = c(0,0,4,0))
for (i in 1:length(genes_for_plot)) {
  gen <- genes_for_plot[i]
  plot_gene(gen)
}
mtext("Top four downregulated genes KO Sham vs WT Sham", outer = TRUE, line = 1, cex = 1, font = 2)
dev.off(); cat("Saved PDF:", comp_file, "\n")


# ============================================================
# 10) CONTRAST 5 — WT CLP vs WT Sham
#     Visual output: 2 PDFs
# ============================================================

all_deseq
comp_file <- all_deseq[5]; comp_file
title <- gsub("_010525.txt", "", comp_file); title

comp_file_deg <- read.table(comp_file, header = TRUE)
names(comp_file_deg); head(comp_file_deg)
if (any(is.na(comp_file_deg$gene_symbol))) cat("Warning: NA in gene_symbol.\n") else cat("No NA in gene_symbol.\n")
if (any(duplicated(comp_file_deg$gene_symbol))) comp_file_deg$gene_symbol <- make.unique(comp_file_deg$gene_symbol) else cat("No duplicates in gene_symbol.\n")

lim.P  <- 0.05
lim.FC <- log2(2)
filtered_deg <- subset(comp_file_deg, comp_file_deg$padj < lim.P & abs(comp_file_deg$log2FoldChange) > lim.FC)
dim(comp_file_deg); dim(filtered_deg)

# UP
sorted_deg <- filtered_deg[order(filtered_deg$log2FoldChange, decreasing = TRUE), ]
names(sorted_deg); sorted_deg[1:20, c(1,2,28,32)]
all_genes <- sorted_deg$gene_symbol; genes_for_plot <- as.character(all_genes[1:4]); genes_for_plot; length(genes_for_plot)

comp_file
comp_file <- "Top4_upregulated_genesWTCLPvsWTSham.pdf"
pdf(comp_file, width = 10, height = 7)
par(mfrow = c(2,2), mar = c(8,5,2,2), oma = c(0,0,4,0))
for (i in 1:length(genes_for_plot)) {
  gen <- genes_for_plot[i]
  filtered_data <- filtered_deg[filtered_deg$gene_symbol == gen, ]
  if (nrow(filtered_data) == 0) { message("Skipping gene ", gen); next }
  if ("group" %in% colnames(filtered_data)) {
    if (length(unique(filtered_data$group)) < 2) { message("Skipping gene ", gen); next }
  }
  print(paste("Plotting gene:", gen))
  plot_gene(gen)
}
mtext("Top four upregulated genes WT CLP vs WT Sham", outer = TRUE, line = 1, cex = 1, font = 2)
dev.off(); cat("Saved PDF:", comp_file, "\n")

# DOWN
sorted_deg <- filtered_deg[order(filtered_deg$log2FoldChange, decreasing = FALSE), ]
names(sorted_deg); sorted_deg[1:20, c(1,2,28,32)]
all_genes <- sorted_deg$gene_symbol; genes_for_plot <- as.character(all_genes[1:4]); genes_for_plot

comp_file
comp_file <- "Top4_downregulated_genesWTCLPvsWTSham.pdf"
pdf(comp_file, width = 10, height = 7)
par(mfrow = c(2,2), mar = c(8,5,2,2), oma = c(0,0,4,0))
for (i in 1:length(genes_for_plot)) {
  gen <- genes_for_plot[i]
  plot_gene(gen)
}
mtext("Top four downregulated genes WT CLP vs WT Sham", outer = TRUE, line = 1, cex = 1, font = 2)
dev.off(); cat("Saved PDF:", comp_file, "\n")


# ============================================================
# 11) VISUAL OUTPUT NOTES 
# ============================================================
# • Each PDF contains 4 panels (2x2). Each panel is one high-priority gene from the selected contrast.
# • Beeswarm points: single samples. Y-axis = log2 expression from normalized matrix.
# • X-axis groups: WT_CLP, TRPV1_KO_CLP, TRPV1_KO_Sham, WT_Sham (fixed order for comparability).
# • Boxplots summarize median and IQR per group to aid visual comparison.
# • The internal ANOVA (u ~ f) is printed in console for QC (not printed on the plot by default).
# • Use these plots to visually confirm direction and dispersion of top DEGs before downstream validation.
# ============================================================