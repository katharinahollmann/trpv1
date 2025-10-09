# ============================================================
# VOLCANO PLOT SCRIPT — EXPLAINED -- xx KH
# TRPV1 RNA-seq dataset: DEG visualization and top-20 gene extraction
# Notes:
#  • Object names, filenames my own, logic unchanged to first skript from Kls
#  • Code reorganized and commented for clarity to best of my knowledge - no guarantee.
#  • Each contrast processed sequentially with identical logic.
# ============================================================


# ============================================================
# 1) SETUP AND INPUT DISCOVERY
# ============================================================

# Packages
# if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
# BiocManager::install("EnhancedVolcano")
library(EnhancedVolcano)

# Working directory -- adjust to your system !! this is absolute path (for now)
setwd("/Users/katharinahollmann/GitHubDesktop/Wagner_comparative_analyses/TRPV1/data")
dir()

# Discover DESeq2 result files for TRPV1
all_data  <- dir(); all_data
all_deseq <- grep("DESeq2_TRPV1", all_data, value = TRUE); all_deseq


# ============================================================
# 2) CONTRAST 1 — FROM all_deseq[1]
#    KO CLP vs KO Sham (per your naming)
# ============================================================

# ---- Select file and derive title ----
comp_file <- all_deseq[1]        # adjust index if needed
comp_file                         # confirm selection
title <- gsub("_080525.txt", "", comp_file); title

# ---- Load DE table ----
comp_file_deg <- read.table(comp_file, header = TRUE)
names(comp_file_deg)
head(comp_file_deg)

# ---- Thresholds ----
lim.P  <- 0.05
lim.FC <- log2(2)

# ---- Filter DEGs ----
filtered_deg <- subset(comp_file_deg, comp_file_deg$padj < lim.P & abs(comp_file_deg$log2FoldChange) > lim.FC)
dim(comp_file_deg); dim(filtered_deg)

# ---- Order by log2FC (increasing) ----
sorted_deg <- filtered_deg[order(filtered_deg$log2FoldChange, decreasing = FALSE),]
dim(sorted_deg); names(sorted_deg)

# Preview first 20 (ID, symbol, log2FC, padj columns as in your indices)
sorted_deg[1:20, c(1, 2, 28, 32)]

# ---- Split up/down after ordering ----
up_reg   <- which(sorted_deg$log2FoldChange > 0)
up       <- sorted_deg[up_reg,]
down_reg <- which(sorted_deg$log2FoldChange < 0)
down     <- sorted_deg[down_reg,]

# Counts for slides
dim(sorted_deg); dim(up); dim(down)

# ============================================================
# Volcano Plot Preparation — build labels and exports
# ============================================================

# Sort full table by decreasing log2FC for plotting
all_genes <- comp_file_deg[order(comp_file_deg$log2FoldChange, decreasing = TRUE),]
row.names(all_genes) <- all_genes$gene_symbol
names(all_genes)
head(row.names(all_genes), 10)

# -------- UP set (export + top 20 selection) --------
top_up_deg <- subset(comp_file_deg, comp_file_deg$log2FoldChange > lim.FC & comp_file_deg$padj < lim.P)
row.names(top_up_deg) <- as.character(top_up_deg$gene_symbol)
top_up_deg_sorted <- top_up_deg[order(top_up_deg$log2FoldChange, decreasing = TRUE),]

# Export table (unchanged filename)
deg_export_up <- top_up_deg_sorted[, c("entrezgene_id", "gene_symbol", "log2FoldChange", "padj")]
head(deg_export_up); dim(deg_export_up)
comp_file
write.table(deg_export_up, file = "TRPV1SortedDegKOCLPvsKOSham_UP_020625.txt")

# Show top 20 for slide
top_up_deg_sorted[1:20, c(1, 2, 28, 32)]

# Keep only top 20 rows (all columns kept), filter NA symbols, build labels
top_up_deg_sorted <- top_up_deg_sorted[1:20,]
top_up_deg_sorted[, c(1,2,28,32)]
head(top_up_deg_sorted)
top_up_final <- top_up_deg_sorted[complete.cases(top_up_deg_sorted$gene_symbol),]
row.names(top_up_final) <- as.character(top_up_final$gene_symbol)
names_top_up <- row.names(top_up_final)[1:20]
names_top_up <- names_top_up[complete.cases(names_top_up)]
names_top_up

# -------- DOWN set (export + top 20 selection) --------
top_down_deg <- subset(comp_file_deg, comp_file_deg$log2FoldChange < (-log2(2)) & comp_file_deg$padj < lim.P)
row.names(top_down_deg) <- as.character(top_down_deg$gene_symbol)
top_down_deg_sorted <- top_down_deg[order(top_down_deg$log2FoldChange, decreasing = FALSE),]

# Export table (unchanged filename commented in your original)
class(top_down_deg_sorted); head(top_down_deg_sorted)
deg_export_down <- top_down_deg_sorted[, c("entrezgene_id", "gene_symbol", "log2FoldChange", "padj")]
head(deg_export_down); dim(deg_export_down)
comp_file
# write.table(deg_export_down, file = "TRPV1SortedDegKOCLPvsKOSham_DOWN_020625.txt")

# Show top 20 for slide
top_down_deg_sorted[1:20, c(1, 2, 28, 32)]

# Keep only top 20 rows, filter NA symbols, build labels
top_down_deg_sorted <- top_down_deg_sorted[1:20,]
top_down_deg_sorted[, c(1,2,28,32)]
head(top_down_deg_sorted)
top_down_final <- top_down_deg_sorted[complete.cases(top_down_deg_sorted$gene_symbol),]
row.names(top_down_final) <- as.character(top_down_final$gene_symbol)
names_top_down <- row.names(top_down_final)[1:20]
names_top_down <- names_top_down[complete.cases(names_top_down)]
names_top_down

# Combined labels for plot
all_names <- c(names_top_up, names_top_down)

# Legend labels and cutoffs
pCut1  <- lim.P
LFC1   <- paste0("absLFC>", round(lim.FC, 2))
PP1    <- paste0("adjP<", pCut1)
leg_lab <- c("NS", LFC1, PP1, "DEGs")

# Ranges for safety
range(all_genes$padj, na.rm = TRUE)
range(all_genes$log2FoldChange, na.rm = TRUE)
head(row.names(all_genes), 20)

# ============================================================
# Volcano Plot — Contrast 1
# ============================================================

par(mfrow = c(1,1), font = 2, font.axis = 2, font.lab = 3, mar = c(4,4,3,2), xpd = FALSE)
first <- EnhancedVolcano(
  all_genes,
  lab = row.names(all_genes),
  x = "log2FoldChange",
  y = "padj",
  xlab = bquote(~Log[2]~ "fold change"),
  col = c("grey", "orange", "lightblue", "red"),
  colAlpha = 4/5,
  pCutoff = pCut1,
  FCcutoff = lim.FC,
  pointSize = 3.0,
  labSize = 4,
  max.overlaps = 50,
  selectLab = all_names,
  xlim = c(-10, 10),
  # title = paste(title, "TOP 20UP + 20DN by LFC"),
  legendLabels = leg_lab,
  drawConnectors = TRUE,
  widthConnectors = 0.75
)
first
ggsave("Volcano_KOCLPvsKOSham.pdf", plot = first, width = 10, height = 8)


# ============================================================
# 3) CONTRAST 2 — FROM all_deseq[2]
#    KO CLP vs WT CLP
# ============================================================

all_deseq
comp_file <- all_deseq[2]
comp_file
title <- gsub("_080525.txt", "", comp_file); title

comp_file_deg <- read.table(comp_file, header = TRUE)
names(comp_file_deg); head(comp_file_deg)

lim.P  <- 0.05
lim.FC <- log2(2)

filtered_deg <- subset(comp_file_deg, comp_file_deg$padj < lim.P & abs(comp_file_deg$log2FoldChange) > lim.FC)
dim(comp_file_deg); dim(filtered_deg)

sorted_deg <- filtered_deg[order(filtered_deg$log2FoldChange, decreasing = FALSE),]
dim(sorted_deg); names(sorted_deg)
sorted_deg[1:20, c(1, 2, 28, 32)]

up_reg   <- which(sorted_deg$log2FoldChange > 0);   up   <- sorted_deg[up_reg,]
down_reg <- which(sorted_deg$log2FoldChange < 0);   down <- sorted_deg[down_reg,]
dim(sorted_deg); dim(up); dim(down)

# Prep
library(EnhancedVolcano)
names(comp_file_deg)
all_genes <- comp_file_deg[order(comp_file_deg$log2FoldChange, decreasing = TRUE),]
row.names(all_genes) <- all_genes$gene_symbol
names(all_genes); head(row.names(all_genes), 10)

# UP
comp_file
top_up_deg <- subset(comp_file_deg, comp_file_deg$log2FoldChange > lim.FC & comp_file_deg$padj < lim.P)
row.names(top_up_deg) <- as.character(top_up_deg$gene_symbol)
top_up_deg_sorted <- top_up_deg[order(top_up_deg$log2FoldChange, decreasing = TRUE),]
class(top_up_deg_sorted); head(top_up_deg_sorted)
deg_export_up <- top_up_deg_sorted[, c("entrezgene_id", "gene_symbol", "log2FoldChange", "padj")]
head(deg_export_up); dim(deg_export_up)
comp_file
write.table(deg_export_up, file = "TRPV1SortedDegKOCLPvsCLP_UP_020625.txt")
top_up_deg_sorted[1:20, c(1,2, 28, 32)]
top_up_deg_sorted <- top_up_deg_sorted[1:20,]
top_up_deg_sorted[, c(1,2,28,32)]; head(top_up_deg_sorted)
top_up_final <- top_up_deg_sorted[complete.cases(top_up_deg_sorted$gene_symbol),]
row.names(top_up_final) <- as.character(top_up_final$gene_symbol)
names_top_up <- row.names(top_up_final)[1:20]
names_top_up <- names_top_up[complete.cases(names_top_up)]
names_top_up

# DOWN
top_down_deg <- subset(comp_file_deg, comp_file_deg$log2FoldChange < (-log2(2)) & comp_file_deg$padj < lim.P)
row.names(top_down_deg) <- as.character(top_down_deg$gene_symbol)
top_down_deg_sorted <- top_down_deg[order(top_down_deg$log2FoldChange, decreasing = FALSE),]
class(top_down_deg_sorted); head(top_down_deg_sorted)
deg_export_down <- top_down_deg_sorted[, c("entrezgene_id", "gene_symbol", "log2FoldChange", "padj")]
head(deg_export_down); dim(deg_export_down)
comp_file
write.table(deg_export_down, file = "TRPV1SortedDegKOCLPvsCLP_DOWN_020625.txt")
top_down_deg_sorted[1:20, c(1,2, 28, 32)]
top_down_deg_sorted <- top_down_deg_sorted[1:20,]
top_down_deg_sorted[, c(1,2,28,32)]; head(top_down_deg_sorted)
top_down_final <- top_down_deg_sorted[complete.cases(top_down_deg_sorted$gene_symbol),]
row.names(top_down_final) <- as.character(top_down_final$gene_symbol)
names_top_down <- row.names(top_down_final)[1:20]
names_top_down <- names_top_down[complete.cases(names_top_down)]
names_top_down

# Labels + legend
all_names <- c(names_top_up, names_top_down)
pCut1  <- lim.P
LFC1   <- paste0("absLFC>", round(lim.FC, 2))
PP1    <- paste0("adjP<", pCut1)
leg_lab <- c("NS", LFC1, PP1, "DEGs")
range(all_genes$padj, na.rm = TRUE)
range(all_genes$log2FoldChange, na.rm = TRUE)
head(row.names(all_genes), 20)

# Plot
par(mfrow = c(1,1), font = 2, font.axis = 2, font.lab = 3, mar = c(4,4,3,2), xpd = FALSE)
second <- EnhancedVolcano(
  all_genes,
  lab = row.names(all_genes),
  x = "log2FoldChange",
  y = "padj",
  xlab = bquote(~Log[2]~ "fold change"),
  col = c("grey", "orange", "lightblue", "red"),
  colAlpha = 4/5,
  pCutoff = pCut1,
  FCcutoff = lim.FC,
  pointSize = 3.0,
  labSize = 4,
  max.overlaps = 50,
  selectLab = all_names,
  xlim = c(-10, 10),
  legendLabels = leg_lab,
  drawConnectors = TRUE,
  widthConnectors = 0.75
)
second
ggsave("Volcano_KOCLPWTCLP.pdf", plot = second, width = 10, height = 8)


# ============================================================
# 4) CONTRAST 3 — FROM all_deseq[3]
#    KO CLP vs WT Sham
# ============================================================

comp_file <- all_deseq[3]
comp_file
title <- gsub("_080525.txt", "", comp_file); title

comp_file_deg <- read.table(comp_file, header = TRUE)
names(comp_file_deg); head(comp_file_deg)

lim.P  <- 0.05
lim.FC <- log2(2)

filtered_deg <- subset(comp_file_deg, comp_file_deg$padj < lim.P & abs(comp_file_deg$log2FoldChange) > lim.FC)
dim(comp_file_deg); dim(filtered_deg)

sorted_deg <- filtered_deg[order(filtered_deg$log2FoldChange, decreasing = FALSE),]
dim(sorted_deg); names(sorted_deg)
sorted_deg[1:20, c(1, 2, 28, 32)]

up_reg   <- which(sorted_deg$log2FoldChange > 0);   up   <- sorted_deg[up_reg,]
down_reg <- which(sorted_deg$log2FoldChange < 0);   down <- sorted_deg[down_reg,]
dim(sorted_deg); dim(up); dim(down)

# Prep
library(EnhancedVolcano)
names(comp_file_deg)
all_genes <- comp_file_deg[order(comp_file_deg$log2FoldChange, decreasing = TRUE),]
row.names(all_genes) <- all_genes$gene_symbol
names(all_genes); head(row.names(all_genes), 10)

# UP
comp_file
top_up_deg <- subset(comp_file_deg, comp_file_deg$log2FoldChange > lim.FC & comp_file_deg$padj < lim.P)
row.names(top_up_deg) <- as.character(top_up_deg$gene_symbol)
top_up_deg_sorted <- top_up_deg[order(top_up_deg$log2FoldChange, decreasing = TRUE),]
class(top_up_deg_sorted); head(top_up_deg_sorted)
deg_export_up <- top_up_deg_sorted[, c("entrezgene_id", "gene_symbol", "log2FoldChange", "padj")]
head(deg_export_up); dim(deg_export_up)
comp_file
write.table(deg_export_up, file = "TRPV1SortedDegKOCLPvsSham_UP_020625.txt")
top_up_deg_sorted[1:20, c(1,2, 28, 32)]
top_up_deg_sorted <- top_up_deg_sorted[1:20,]
top_up_deg_sorted[, c(1,2,28,32)]; head(top_up_deg_sorted)
top_up_final <- top_up_deg_sorted[complete.cases(top_up_deg_sorted$gene_symbol),]
row.names(top_up_final) <- as.character(top_up_final$gene_symbol)
names_top_up <- row.names(top_up_final)[1:20]
names_top_up <- names_top_up[complete.cases(names_top_up)]
names_top_up

# DOWN
top_down_deg <- subset(comp_file_deg, comp_file_deg$log2FoldChange < (-log2(2)) & comp_file_deg$padj < lim.P)
row.names(top_down_deg) <- as.character(top_down_deg$gene_symbol)
top_down_deg_sorted <- top_down_deg[order(top_down_deg$log2FoldChange, decreasing = FALSE),]
class(top_down_deg_sorted); head(top_down_deg_sorted)
deg_export_down <- top_down_deg_sorted[, c("entrezgene_id", "gene_symbol", "log2FoldChange", "padj")]
head(deg_export_down); dim(deg_export_down)
comp_file
write.table(deg_export_down, file = "TRPV1SortedDegKOCLPvsSham_DOWN_020625.txt")
top_down_deg_sorted[1:20, c(1,2, 28, 32)]
top_down_deg_sorted <- top_down_deg_sorted[1:20,]
top_down_deg_sorted[, c(1,2,28,32)]; head(top_down_deg_sorted)
top_down_final <- top_down_deg_sorted[complete.cases(top_down_deg_sorted$gene_symbol),]
row.names(top_down_final) <- as.character(top_down_final$gene_symbol)
names_top_down <- row.names(top_down_final)[1:20]
names_top_down <- names_top_down[complete.cases(names_top_down)]
names_top_down

# Labels + legend
all_names <- c(names_top_up, names_top_down)
pCut1  <- lim.P
LFC1   <- paste0("absLFC>", round(lim.FC, 2))
PP1    <- paste0("adjP<", pCut1)
leg_lab <- c("NS", LFC1, PP1, "DEGs")
range(all_genes$padj, na.rm = TRUE)
range(all_genes$log2FoldChange, na.rm = TRUE)
head(row.names(all_genes), 20)

# Plot
par(mfrow = c(1,1), font = 2, font.axis = 2, font.lab = 3, mar = c(4,4,3,2), xpd = FALSE)
third <- EnhancedVolcano(
  all_genes,
  lab = row.names(all_genes),
  x = "log2FoldChange",
  y = "padj",
  xlab = bquote(~Log[2]~ "fold change"),
  col = c("grey", "orange", "lightblue", "red"),
  colAlpha = 4/5,
  pCutoff = pCut1,
  FCcutoff = lim.FC,
  pointSize = 3.0,
  labSize = 4,
  max.overlaps = 50,
  selectLab = all_names,
  xlim = c(-10, 10),
  legendLabels = leg_lab,
  drawConnectors = TRUE,
  widthConnectors = 0.75
)
ggsave("Volcano_KOCLPWTSham.pdf", plot = third, width = 10, height = 8)


# ============================================================
# 5) CONTRAST 4 — FROM all_deseq[4]
#    KO Sham vs WT Sham
# ============================================================

comp_file <- all_deseq[4]
comp_file
title <- gsub("_080525.txt", "", comp_file); title

comp_file_deg <- read.table(comp_file, header = TRUE)
names(comp_file_deg); head(comp_file_deg)

lim.P  <- 0.05
lim.FC <- log2(2)

filtered_deg <- subset(comp_file_deg, comp_file_deg$padj < lim.P & abs(comp_file_deg$log2FoldChange) > lim.FC)
dim(comp_file_deg); dim(filtered_deg)

sorted_deg <- filtered_deg[order(filtered_deg$log2FoldChange, decreasing = FALSE),]
dim(sorted_deg); names(sorted_deg)
sorted_deg[1:20, c(1, 2, 28, 32)]

up_reg   <- which(sorted_deg$log2FoldChange > 0);   up   <- sorted_deg[up_reg,]
down_reg <- which(sorted_deg$log2FoldChange < 0);   down <- sorted_deg[down_reg,]
dim(sorted_deg); dim(up); dim(down)

# Prep
library(EnhancedVolcano)
names(comp_file_deg)
all_genes <- comp_file_deg[order(comp_file_deg$log2FoldChange, decreasing = TRUE),]
row.names(all_genes) <- all_genes$gene_symbol
names(all_genes); head(row.names(all_genes), 10)

# UP
comp_file
top_up_deg <- subset(comp_file_deg, comp_file_deg$log2FoldChange > lim.FC & comp_file_deg$padj < lim.P)
row.names(top_up_deg) <- as.character(top_up_deg$gene_symbol)
top_up_deg_sorted <- top_up_deg[order(top_up_deg$log2FoldChange, decreasing = TRUE),]
class(top_up_deg_sorted); head(top_up_deg_sorted)
deg_export_up <- top_up_deg_sorted[, c("entrezgene_id", "gene_symbol", "log2FoldChange", "padj")]
head(deg_export_up); dim(deg_export_up)
comp_file
write.table(deg_export_up, file = "TRPV1SortedDegKOShamvsWTSham_UP_020625.txt")
top_up_deg_sorted[1:20, c(1,2, 28, 32)]
top_up_deg_sorted <- top_up_deg_sorted[1:20,]
top_up_deg_sorted[, c(1,2,28,32)]; head(top_up_deg_sorted)
top_up_final <- top_up_deg_sorted[complete.cases(top_up_deg_sorted$gene_symbol),]
row.names(top_up_final) <- as.character(top_up_final$gene_symbol)
names_top_up <- row.names(top_up_final)[1:20]
names_top_up <- names_top_up[complete.cases(names_top_up)]
names_top_up

# DOWN
top_down_deg <- subset(comp_file_deg, comp_file_deg$log2FoldChange < (-log2(2)) & comp_file_deg$padj < lim.P)
row.names(top_down_deg) <- as.character(top_down_deg$gene_symbol)
top_down_deg_sorted <- top_down_deg[order(top_down_deg$log2FoldChange, decreasing = FALSE),]
class(top_down_deg_sorted); head(top_down_deg_sorted)
deg_export_down <- top_down_deg_sorted[, c("entrezgene_id", "gene_symbol", "log2FoldChange", "padj")]
head(deg_export_down); dim(deg_export_down)
comp_file
write.table(deg_export_down, file = "TRPV1SortedDegKOShamvsWTSham_DOWN_020625.txt")
top_down_deg_sorted[1:20, c(1,2, 28, 32)]
top_down_deg_sorted <- top_down_deg_sorted[1:20,]
top_down_deg_sorted[, c(1,2,28,32)]; head(top_down_deg_sorted)
top_down_final <- top_down_deg_sorted[complete.cases(top_down_deg_sorted$gene_symbol),]
row.names(top_down_final) <- as.character(top_down_final$gene_symbol)
names_top_down <- row.names(top_down_final)[1:20]
names_top_down <- names_top_down[complete.cases(names_top_down)]
names_top_down

# Labels + legend
all_names <- c(names_top_up, names_top_down)
pCut1  <- lim.P
LFC1   <- paste0("absLFC>", round(lim.FC, 2))
PP1    <- paste0("adjP<", pCut1)
leg_lab <- c("NS", LFC1, PP1, "DEGs")
range(all_genes$padj, na.rm = TRUE)
range(all_genes$log2FoldChange, na.rm = TRUE)
head(row.names(all_genes), 20)

# Plot
par(mfrow = c(1,1), font = 2, font.axis = 2, font.lab = 3, mar = c(4,4,3,2), xpd = FALSE)
forth <- EnhancedVolcano(
  all_genes,
  lab = row.names(all_genes),
  x = "log2FoldChange",
  y = "padj",
  xlab = bquote(~Log[2]~ "fold change"),
  col = c("grey", "orange", "lightblue", "red"),
  colAlpha = 4/5,
  pCutoff = pCut1,
  FCcutoff = lim.FC,
  pointSize = 3.0,
  labSize = 4,
  max.overlaps = 50,
  selectLab = all_names,
  xlim = c(-10, 10),
  legendLabels = leg_lab,
  drawConnectors = TRUE,
  widthConnectors = 0.75
)
forth
ggsave("Volcano_KOWTSham.pdf", plot = forth, width = 10, height = 8)

comp_file


# ============================================================
# 6) CONTRAST 5 — FROM all_deseq[5]
#    WT CLP vs WT Sham
# ============================================================

comp_file <- all_deseq[5]
comp_file
# Note: 5th title used a different suffix in your original
title <- gsub("_010525.txt", "", comp_file); title

comp_file_deg <- read.table(comp_file, header = TRUE)
names(comp_file_deg); head(comp_file_deg)

lim.P  <- 0.05
lim.FC <- log2(2)

filtered_deg <- subset(comp_file_deg, comp_file_deg$padj < lim.P & abs(comp_file_deg$log2FoldChange) > lim.FC)
dim(comp_file_deg); dim(filtered_deg)

sorted_deg <- filtered_deg[order(filtered_deg$log2FoldChange, decreasing = FALSE),]
dim(sorted_deg); names(sorted_deg)
sorted_deg[1:20, c(1, 2, 28, 32)]

up_reg   <- which(sorted_deg$log2FoldChange > 0);   up   <- sorted_deg[up_reg,]
down_reg <- which(sorted_deg$log2FoldChange < 0);   down <- sorted_deg[down_reg,]
dim(sorted_deg); dim(up); dim(down)

# Prep
library(EnhancedVolcano)
names(comp_file_deg)
all_genes <- comp_file_deg[order(comp_file_deg$log2FoldChange, decreasing = TRUE),]
row.names(all_genes) <- all_genes$gene_symbol
names(all_genes); head(row.names(all_genes), 10)

# UP
comp_file
top_up_deg <- subset(comp_file_deg, comp_file_deg$log2FoldChange > lim.FC & comp_file_deg$padj < lim.P)
row.names(top_up_deg) <- as.character(top_up_deg$gene_symbol)
top_up_deg_sorted <- top_up_deg[order(top_up_deg$log2FoldChange, decreasing = TRUE),]
class(top_up_deg_sorted); head(top_up_deg_sorted)
deg_export_up <- top_up_deg_sorted[, c("entrezgene_id", "gene_symbol", "log2FoldChange", "padj")]
head(deg_export_up); dim(deg_export_up)
comp_file
write.table(deg_export_up, file = "TRPV1SortedDegCLPvsSham_UP_020625.txt")
top_up_deg_sorted[1:20, c(1,2, 28, 32)]
top_up_deg_sorted <- top_up_deg_sorted[1:20,]
top_up_deg_sorted[, c(1,2,28,32)]; head(top_up_deg_sorted)
top_up_final <- top_up_deg_sorted[complete.cases(top_up_deg_sorted$gene_symbol),]
row.names(top_up_final) <- as.character(top_up_final$gene_symbol)
names_top_up <- row.names(top_up_final)[1:20]
names_top_up <- names_top_up[complete.cases(names_top_up)]
names_top_up

# DOWN
top_down_deg <- subset(comp_file_deg, comp_file_deg$log2FoldChange < (-log2(2)) & comp_file_deg$padj < lim.P)
row.names(top_down_deg) <- as.character(top_down_deg$gene_symbol)
top_down_deg_sorted <- top_down_deg[order(top_down_deg$log2FoldChange, decreasing = FALSE),]
class(top_down_deg_sorted); head(top_down_deg_sorted)
deg_export_down <- top_down_deg_sorted[, c("entrezgene_id", "gene_symbol", "log2FoldChange", "padj")]
head(deg_export_down); dim(deg_export_down)
comp_file
write.table(deg_export_down, file = "TRPV1SortedDegCLPvsSham_DOWN_020625.txt")
top_down_deg_sorted[1:20, c(1,2, 28, 32)]
top_down_deg_sorted <- top_down_deg_sorted[1:20,]
top_down_deg_sorted[, c(1,2,28,32)]; head(top_down_deg_sorted)
top_down_final <- top_down_deg_sorted[complete.cases(top_down_deg_sorted$gene_symbol),]
row.names(top_down_final) <- as.character(top_down_final$gene_symbol)
names_top_down <- row.names(top_down_final)[1:20]
names_top_down <- names_top_down[complete.cases(names_top_down)]
names_top_down

# Labels + legend
all_names <- c(names_top_up, names_top_down)
pCut1  <- lim.P
LFC1   <- paste0("absLFC>", round(lim.FC, 2))
PP1    <- paste0("adjP<", pCut1)
leg_lab <- c("NS", LFC1, PP1, "DEGs")
range(all_genes$padj, na.rm = TRUE)
range(all_genes$log2FoldChange, na.rm = TRUE)
head(row.names(all_genes), 20)

# Plot
par(mfrow = c(1,1), font = 2, font.axis = 2, font.lab = 3, mar = c(4,4,3,2), xpd = FALSE)
fifth <- EnhancedVolcano(
  all_genes,
  lab = row.names(all_genes),
  x = "log2FoldChange",
  y = "padj",
  xlab = bquote(~Log[2]~ "fold change"),
  col = c("grey", "orange", "lightblue", "red"),
  colAlpha = 4/5,
  pCutoff = pCut1,
  FCcutoff = lim.FC,
  pointSize = 3.0,
  labSize = 4,
  max.overlaps = 50,
  selectLab = all_names,
  xlim = c(-10, 10),
  legendLabels = leg_lab,
  drawConnectors = TRUE,
  widthConnectors = 0.75
)
fifth
ggsave("Volcano_CLPvsShamTRPV1.pdf", plot = fifth, width = 10, height = 8)


# ============================================================
# 7) NOTES
# ============================================================
# • Top-20 extraction is performed per direction after thresholding.
# • Export to your filenames for reproducibility and reuse.
# • adjust pdf size//xlim/ylim for your preference
# • Label set = top 20 up + top 20 down. Adjust N by editing the [1:20] slices.
# • Kept EnhancedVolcano parameters consistent across contrasts for comparability.
# ============================================================