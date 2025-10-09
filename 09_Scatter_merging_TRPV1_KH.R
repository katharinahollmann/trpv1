################################################################################################
# DIFFERENTIAL EXPRESSION COMPARISON SCRIPT — EXPLAINED — xx KH
# TRPV1 RNA-seq | DEG selection, export lists, and CLP vs TRPV1 KO CLP comparisons
#
# Notes on preserved objects (unchanged):
# - all_data, all_deseq : file discovery
# - clpKO, clp : chosen DESeq2 results files (indices preserved)
# - title : derived label for plots
# - clpKO_deg, clp_deg : full DE tables for each contrast
# - lim.P, lim.FC : thresholds for padj and |log2FC|
# - filtered_deg, sorted_deg : DEGs after filtering, then ordered by log2FC
# - up, down : data frames of direction-specific DEGs
# - top_up_deg(_sorted), top_down_deg(_sorted) : directional subsets for export
# - deg_export_* : compact export tables for merging / plotting
# - deg_*_merged(_log2FC) : merged tables by gene_symbol for cross-contrast scatter plots
#
# Brief interpretive notes are added inline (biological/analytical context) without changing logic.
################################################################################################


##############################################
# 1) SETUP: working directory & file inventory
##############################################

setwd("/Users/katharinahollmann/GitHubDesktop/Projects/Wagner_comparative_analyses/TRPV1/data")
dir()

# Load all Files ----------------------------------------------------------
all_data  <- dir(); all_data                          # list all files
all_deseq <- grep("DESeq2", all_data, value = TRUE);  # DESeq2 result files only
all_deseq

# (Interpretation) These are the per-contrast DESeq2 tables you generated earlier. 
# We’ll pick contrasts by their index as you’ve been doing, preserving your selection logic.


##############################################
# 2) CLP TRPV1 KO (SECOND in your list)
##############################################

# Choose comparison (index preserved) -------------------------------------
all_deseq
clpKO <- all_deseq[6]             # Select the DESeq2 results file (adjust index to switch contrast)
clpKO                             # sanity check: is this the intended group?

# Title for plots ---------------------------------------------------------
title <- gsub("_080525.txt", "", clpKO)
title

# Load differential expression table -------------------------------------
clpKO_deg <- read.table(clpKO, header = TRUE)
names(clpKO_deg)
head(clpKO_deg)

# Thresholds --------------------------------------------------------------
lim.P  <- 0.05           # adjusted p-value
lim.FC <- log2(2)        # absolute log2FC > 1

# Filter DEGs -------------------------------------------------------------
filtered_deg <- subset(
  clpKO_deg,
  clpKO_deg$padj < lim.P & abs(clpKO_deg$log2FoldChange) > lim.FC
)

# Counts before/after -----------------------------------------------------
dim(clpKO_deg)       # all genes tested
dim(filtered_deg)    # DEGs at thresholds

# Order by log2FC (ascending) --------------------------------------------
sorted_deg <- filtered_deg[order(filtered_deg$log2FoldChange, decreasing = FALSE), ]
dim(sorted_deg)
names(sorted_deg)

# Split by direction ------------------------------------------------------
up_reg   <- which(sorted_deg$log2FoldChange > 0)
up       <- sorted_deg[up_reg,]

down_reg <- which(sorted_deg$log2FoldChange < 0)
down     <- sorted_deg[down_reg,]

# In PowerPoint (sanity counts) -------------------------------------------
dim(sorted_deg)   # all DEGs at thresholds (ordered)
dim(up)           # upregulated DEGs
dim(down)         # downregulated DEGs

# (Interpretation) Upregulated genes here show higher expression in the CLP TRPV1 KO setting
# versus its control; downregulated genes show lower expression. Magnitude reflects effect size.


# ------------------------ UP (CLP TRPV1 KO) -------------------------

clpKO

top_up_deg <- subset(
  clpKO_deg,
  clpKO_deg$log2FoldChange > log2(2) & clpKO_deg$padj < lim.P
)
names(top_up_deg)

# Ensure gene symbols act as row labels (useful in Excel/export)
row.names(top_up_deg) <- as.character(top_up_deg$gene_symbol)

# Order decreasing by log2FC ----------------------------------------------
top_up_deg_sorted <- top_up_deg[order(top_up_deg$log2FoldChange, decreasing = TRUE), ]

# Prepare export table ----------------------------------------------------
class(top_up_deg_sorted)
head(top_up_deg_sorted)

deg_export_up_clpKO <- top_up_deg_sorted[, c("entrezgene_id", "gene_symbol", "log2FoldChange", "padj")]
head(deg_export_up_clpKO)
dim(deg_export_up_clpKO)

clpKO
# write.table(deg_export_up_clpKO, file = "SortedDegclpKOvsShamAKDATA_UP_020625.txt",
#             sep = "\t", quote = FALSE, row.names = FALSE)

# (Interpretation) This UP list is ideal for GO/pathway input, heatmaps, and 
# candidate prioritization (high log2FC & significant padj).


# ------------------------ DOWN (CLP TRPV1 KO) -----------------------

top_down_deg <- subset(
  clpKO_deg,
  clpKO_deg$log2FoldChange < (-log2(2)) & clpKO_deg$padj < lim.P
)
names(top_down_deg)

row.names(top_down_deg) <- as.character(top_down_deg$gene_symbol)

# Order increasing by log2FC (more negative first) ------------------------
top_down_deg_sorted <- top_down_deg[order(top_down_deg$log2FoldChange, decreasing = FALSE), ]

# Prepare export table ----------------------------------------------------
class(top_down_deg_sorted)
head(top_down_deg_sorted)

deg_export_down_clpKO <- top_down_deg_sorted[, c("entrezgene_id", "gene_symbol", "log2FoldChange", "padj")]
head(deg_export_down_clpKO)
dim(deg_export_down_clpKO)

clpKO
# write.table(deg_export_down_clpKO, file = "SortedDegclpKOvsShamAKDATA_DOWN_020625.txt",
#             sep = "\t", quote = FALSE, row.names = FALSE)

# (Interpretation) This DOWN list captures suppressed programs; check for pathways indicating
# loss of function or dampened inflammatory/metabolic response in TRPV1 KO 


##############################################
# 3) CLP only (THIRD in list)
##############################################

all_deseq
clp <- all_deseq[8]               # Select the DESeq2 results file (index preserved)
clp                               # sanity check

# Title -------------------------------------------------------------------
title <- gsub("_080525.txt", "", clp)
title

# Load table --------------------------------------------------------------
clp_deg <- read.table(clp, header = TRUE)
names(clp_deg)
head(clp_deg)

# Thresholds (same) -------------------------------------------------------
lim.P  <- 0.05
lim.FC <- log2(2)

# Filter DEGs -------------------------------------------------------------
filtered_deg <- subset(
  clp_deg,
  clp_deg$padj < lim.P & abs(clp_deg$log2FoldChange) > lim.FC
)

# Counts before/after -----------------------------------------------------
dim(clp_deg)
dim(filtered_deg)

# Order by log2FC (ascending) --------------------------------------------
sorted_deg <- filtered_deg[order(filtered_deg$log2FoldChange, decreasing = FALSE), ]
dim(sorted_deg)
names(sorted_deg)

# Split by direction ------------------------------------------------------
up_reg <- which(sorted_deg$log2FoldChange > 0)
up     <- sorted_deg[up_reg,]

down_reg <- which(sorted_deg$log2FoldChange < 0)
down     <- sorted_deg[down_reg,]

# In PowerPoint (sanity counts) -------------------------------------------
dim(sorted_deg)
dim(up)
dim(down)

# (Interpretation) This contrast reflects WT CLP vs sham (or your control per file),
# providing the baseline septic response without TRPV1 knockout. Useful for KO vs WT comparison.


# ------------------------ UP (CLP) ---------------------------------------

clp

top_up_deg <- subset(
  clp_deg,
  clp_deg$log2FoldChange > log2(2) & clp_deg$padj < lim.P
)
names(top_up_deg)
# top_up_deg[1:20, c(1, 2, 33, 37)] # quick peek if needed

row.names(top_up_deg) <- as.character(top_up_deg$gene_symbol)

# Order decreasing by log2FC ----------------------------------------------
top_up_deg_sorted <- top_up_deg[order(top_up_deg$log2FoldChange, decreasing = TRUE), ]

# Prepare export table ----------------------------------------------------
class(top_up_deg_sorted)
head(top_up_deg_sorted)

deg_export_up_CLP <- top_up_deg_sorted[, c("entrezgene_id", "gene_symbol", "log2FoldChange", "padj")]
head(deg_export_up_CLP)
dim(deg_export_up_CLP)

clp
# write.table(deg_export_up_CLP, file = "SortedDegCLPvsSham_UP_020625.txt",
#             sep = "\t", quote = FALSE, row.names = FALSE)


# ------------------------ DOWN (CLP) -------------------------------------

top_down_deg <- subset(
  clp_deg,
  clp_deg$log2FoldChange < (-log2(2)) & clp_deg$padj < lim.P
)
names(top_down_deg)
# top_down_deg[1:20, c(1, 2, 33, 37)] # quick peek if needed

row.names(top_down_deg) <- as.character(top_down_deg$gene_symbol)

# Order increasing by log2FC ----------------------------------------------
top_down_deg_sorted <- top_down_deg[order(top_down_deg$log2FoldChange, decreasing = FALSE), ]

# Prepare export table ----------------------------------------------------
class(top_down_deg_sorted)
head(top_down_deg_sorted)

deg_export_down_CLP <- top_down_deg_sorted[, c("entrezgene_id", "gene_symbol", "log2FoldChange", "padj")]
head(deg_export_down_CLP)
dim(deg_export_down_CLP)

clp
# write.table(deg_export_down_CLP, file = "SortedDegCLPvsShamAKDATA_DOWN_020625.txt",
#             sep = "\t", quote = FALSE, row.names = FALSE)


##############################################
# 4) MERGE & SCATTER — UP lists and DOWN lists (CLP vs TRPV1 KO CLP)
##############################################

# (Biology) By merging UP (or DOWN) genes that pass thresholds in both contrasts,
# we focus on the shared responders and directly compare their effect sizes (log2FC).
# Points far from the diagonal (y = x) highlight genes with KO-sensitive modulation.

library(dplyr)

# ------------------------ UP: merge & plot --------------------------------

deg_up_merged <- merge(
  deg_export_up_CLP,
  deg_export_up_clpKO,
  by = "gene_symbol",
  suffixes = c("_clp", "_clpKO")
)
head(deg_up_merged)

deg_up_merged_log2FC <- deg_up_merged |>
  select(gene_symbol, log2FoldChange_clp, log2FoldChange_clpKO)

head(deg_up_merged_log2FC)

# Sort by WT CLP log2FC (descending) --------------------------------------
deg_up_merged_log2FC <- deg_up_merged_log2FC[
  order(deg_up_merged_log2FC$log2FoldChange_clp, decreasing = TRUE), 
]
head(deg_up_merged_log2FC)
dim(deg_up_merged_log2FC)

write.table(
  deg_up_merged_log2FC,
  file = "DEG_UP_merged_CLPandCLPKOSham_log2FC_050625.txt",
  sep = "\t", quote = FALSE, row.names = FALSE
)

# Scatter: UP genes -------------------------------------------------------
library(ggplot2)
library(ggrepel)

ggplot(deg_up_merged_log2FC, aes(x = log2FoldChange_clp, y = log2FoldChange_clpKO)) +
  geom_point(color = "steelblue", size = 1.5) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  labs(
    title    = "Comparison of log2FC of merged genes",
    subtitle = "upregulated in CLP and TRPV1 KO CLP compared to Sham",
    x        = "log2FC (CLP)",
    y        = "log2FC (TRPV1 KO CLP)"
  ) +
  theme_minimal()

# Label genes with large divergence between conditions --------------------
deg_up_merged_log2FC$label <- ifelse(
  abs(deg_up_merged_log2FC$log2FoldChange_clp - deg_up_merged_log2FC$log2FoldChange_clpKO) > 2 |
    abs(deg_up_merged_log2FC$log2FoldChange_clpKO - deg_up_merged_log2FC$log2FoldChange_clp) > 2,
  deg_up_merged_log2FC$gene_symbol, NA
)

ggplot(deg_up_merged_log2FC, aes(x = log2FoldChange_clp, y = log2FoldChange_clpKO)) +
  geom_point(color = "steelblue", size = 2) +
  geom_text_repel(aes(label = label), size = 5, na.rm = TRUE, max.overlaps = 20) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  labs(
    title    = "Comparison of upregulated DEGs",
    subtitle = "upregulated in CLP and TRPV1 KO CLP compared to Sham",
    x        = "log2FC (WT CLP)",
    y        = "log2FC (TRPV1 KO CLP)"
  ) +
  theme_minimal()

ggsave("scatterUPKO_labeledover2diff.pdf", width = 6, height = 4.5)

# (Interpretation) Points above the diagonal: stronger upregulation in KO than WT (KO-boosted).
# Points below: weaker upregulation in KO (KO-attenuated). Labeled genes show |Δlog2FC| > 2.


# ------------------------ DOWN: merge & plot ------------------------------

deg_down_merged <- merge(
  deg_export_down_CLP,
  deg_export_down_clpKO,
  by = "gene_symbol",
  suffixes = c("_clp", "_clpKO")
)
head(deg_down_merged)

deg_down_merged_log2FC <- deg_down_merged |>
  select(gene_symbol, log2FoldChange_clp, log2FoldChange_clpKO)

head(deg_down_merged_log2FC)

# Sort by WT CLP log2FC (increasing; more negative first) -----------------
deg_down_merged_log2FC <- deg_down_merged_log2FC[
  order(deg_down_merged_log2FC$log2FoldChange_clp, decreasing = FALSE),
]
head(deg_down_merged_log2FC)
dim(deg_down_merged_log2FC)

# If you want a table output here as well, uncomment:
# write.table(
#   deg_down_merged_log2FC,
#   file = "DEG_DOWN_merged_CLPSham_and_clpKOSham_log2FC_050625.txt",
#   sep = "\t", quote = FALSE, row.names = FALSE
# )

# Scatter: DOWN genes (plotting negatives with sign preserved) ------------
library(ggplot2)
library(ggrepel)

deg_down_merged_log2FC$label <- ifelse(
  abs(deg_down_merged_log2FC$log2FoldChange_clp - deg_down_merged_log2FC$log2FoldChange_clpKO) > 2 |
    abs(deg_down_merged_log2FC$log2FoldChange_clpKO - deg_down_merged_log2FC$log2FoldChange_clp) > 2,
  deg_down_merged_log2FC$gene_symbol, NA
)

ggplot(deg_down_merged_log2FC, aes(x = -log2FoldChange_clp, y = -log2FoldChange_clpKO)) +
  geom_point(color = "steelblue", size = 2) +
  geom_text_repel(aes(label = label), size = 5, na.rm = TRUE, max.overlaps = 20) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  labs(
    title    = "Comparison of downregulated DEGs",
    subtitle = "downregulated in CLP and TRPV1 KO CLP compared to Sham",
    x        = "- log2FC (WT CLP)",
    y        = "- log2FC (TRPV1 KO CLP)"
  ) +
  theme_minimal()

ggsave("scatterDOWNKO_labeledover2diff.pdf", width = 6, height = 4.5)

# (Interpretation) Because both axes are negated, points above the diagonal correspond to 
# stronger downregulation in KO than WT; below the diagonal indicate weaker downregulation in KO.
# Labels highlight genes with |Δlog2FC| > 2—prime candidates for TRPV1-dependent effects.


################################################################################################
# END — Script preserves your logic, object names, and file indices.
# Added explanatory headers.
################################################################################################