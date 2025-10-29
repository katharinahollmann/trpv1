################################################################################
# HEATMAP SCRIPT — EXPLAINED — KH
# TRPV1 RNA-seq dataset: heatmaps of DEGs (per DESeq2 contrast) and an
# intersect-gene set; uses gplots::heatmap.2
#
# Notes:
# • Object names, filenames, and overall logic kept the same.
# • Reorganized and commented for clarity; no guarantees.
# • Key variable legend:
#     dat1 -> gene_expr_data
#     ta   -> metadata
#     ta2  -> sorted_metadata
#     sm   -> comp_file (for “set” / intersect section)
#     gg1  -> selected_genes
#     dd2  -> filtered_gene_expr
#     sel1 -> sample_id_sorted_metadata
#     dd3  -> ordered_gene_expr
#     x4   -> ordered_gene_expr_matrix
#     xx5  -> scaled_matrix
#     xx6  -> final_scaled_matrix
################################################################################


## ---------------------------------------------------------------------------
## Setup: paths, input tables
## ---------------------------------------------------------------------------

## Working directory with TRPV1 data -- adjust to your system !! this is absolute path (for now)
setwd("/Users/katharinahollmann/GitHubDesktop/Wagner_comparative_analyses/TRPV1/data")
dir()

## Expression matrix (batch-corrected normalized)
gene_expr_data  <- read.table("norm_Wagner_CLP_TRPV1_2023_080525.txt", header = TRUE)
names(gene_expr_data); dim(gene_expr_data)
gene_expr_data[1:5, c(3,16:18)]  # quick preview

## Sample metadata (target)
ta <- read.table("target_Wagner_CLP_TRPV1_2023_240124.txt", header = TRUE, na.strings = c("", "na"))
names(ta); head(ta); table(ta$group)

## Sanity check: column names of expression match sample order in metadata
identical(colnames(gene_expr_data)[9:26], as.character(ta$sample_ID))

## Sort metadata by group for consistent column order in heatmaps
sorted_metadata <- ta[order(ta$group), ]
sorted_metadata$group


## ---------------------------------------------------------------------------
## Helper: list DESeq2 result files (TRPV1 contrasts)
## ---------------------------------------------------------------------------

all_data  <- dir(); all_data
all_deseq <- grep("DESeq2_TRPV1", all_data, value = TRUE); all_deseq
## Expecting 5 TRPV1 DESeq2 contrast files in all_deseq[1:5]


## ---------------------------------------------------------------------------
## Common thresholds for DEG selection
## ---------------------------------------------------------------------------

lim.P  <- 0.05    # padj cutoff
lim.FC <- log2(2) # |log2FC| > 1


## ===========================================================================
## Contrast 1 (all_deseq[1]) — build heatmap from filtered DEGs
## ===========================================================================

comp_file <- all_deseq[1]; comp_file
title     <- gsub("_080525.txt", "", comp_file); title

all_genes <- read.table(comp_file, header = TRUE)
names(all_genes); head(all_genes)

## Filter by thresholds
filtered_genes <- subset(all_genes, all_genes$padj < lim.P & abs(all_genes$log2FoldChange) > lim.FC)
dim(all_genes); dim(filtered_genes)

## Sort by log2FC (descending -> top up at top of table)
sorted_genes <- filtered_genes[order(filtered_genes$log2FoldChange, decreasing = TRUE), ]
names(sorted_genes)
sorted_genes[1:20, c(1,2,28,32)]  # preview: top 20 upregulated (gene_symbol, etc.)

## Count up/down just for reporting sanity
up.reg   <- which(sorted_genes$log2FoldChange > 0)
down.reg <- which(sorted_genes$log2FoldChange < 0)
dim(sorted_genes); dim(sorted_genes[up.reg,]); dim(sorted_genes[down.reg,])

## Keep gene symbols to fetch expression from gene_expr_data
gg1 <- sorted_genes$gene_symbol; head(gg1, 10); length(gg1)

## Filter expression matrix to selected gene symbols
dd2 <- gene_expr_data[gene_expr_data$gene_symbol %in% gg1, ]
dim(dd2); names(dd2)

## Reorder columns by sorted metadata sample list
sel1 <- sorted_metadata$sample_ID; sel1
dd3  <- dd2[, sel1]; names(dd3)
identical(colnames(dd3), as.character(sorted_metadata[,1]))

## Matrix for heatmap
x4 <- as.matrix(dd3); dim(x4); head(x4); range(x4)

## ---------------------------------------------------------------------------
## Heatmap-2 (gplots) pipeline (scaling + caps + clustering + colors)
## ---------------------------------------------------------------------------

library(gplots)

row.names <- dd2$gene_symbol
col.names <- sorted_metadata$group
table(col.names)

## Z-score per gene (row), cap to [-2, 2]
xx5 <- t(scale(t(x4))); range(xx5); hist(xx5)
xx5[xx5 >  2] <-  2
xx6 <- xx5
xx6[xx6 < -2] <- -2
x4 <- xx6; range(x4); hist(x4)

## Color palette
my_palette <- bluered(75)

## Distances and clustering
mydatascale <- t(scale(t(x4)))
hr <- hclust(as.dist(1 - cor(t(mydatascale), method = "pearson")), method = "complete")
hc <- hclust(as.dist(1 - cor(mydatascale,        method = "spearman")), method = "complete")

## Optional row-cluster color bar (not shown by default)
cut_val1 <- 1.5
mycl <- cutree(hr, h = max(hr$height / cut_val1))
clusterCols        <- rainbow(length(unique(mycl)))
myClusterSideBar   <- clusterCols[mycl]

## Column-side group colors (match your preferred mapping)
cc1  <- as.numeric(as.factor(col.names))
from <- unique(cc1)
to   <- c("red", "orange", "darkgreen", "black")
map  <- setNames(to, from)
cc1[] <- map[cc1]
col_col1 <- cc1

## Draw heatmap
par(mfrow = c(1,1), font = 2, font.axis = 2, font.lab = 3, mar = c(4,3,3,2))
hv <- heatmap.2(
  x4, col = my_palette,
  margins = c(9,4),
  trace = "none",
  Colv  = NA,               # keep fixed column order (sorted by group)
  # Rowv = NA,              # (optionally disable row clustering)
  colsep = c(5,10,14), sepcolor = "black", sepwidth = c(0.1, 0.1),
  ColSideColors = col_col1,
  # RowSideColors = myClusterSideBar,
  dendrogram = "row",
  labRow = row.names,
  labCol = col.names,
  cexRow = 0.2, cexCol = 0.8,
  main = paste(title)
)

## Legend (unique group labels and their colors)
col2  <- unique(col_col1)
nam10 <- unique(col.names)
legend("topright",
       legend = nam10,
       pch = 15, cex = 1.1, pt.cex = 1.6, bty = "n", ncol = 2,
       col = col2, inset = c(0.001, 0.03))


## ===========================================================================
## Contrast 2 (all_deseq[2])
## ===========================================================================

all_deseq
comp_file <- all_deseq[2]; comp_file
title     <- gsub("_080525.txt", "", comp_file); title

all_genes <- read.table(comp_file, header = TRUE)
names(all_genes); head(all_genes)

filtered_genes <- subset(all_genes, all_genes$padj < lim.P & abs(all_genes$log2FoldChange) > lim.FC)
dim(all_genes); dim(filtered_genes)

sorted_genes <- filtered_genes[order(filtered_genes$log2FoldChange, decreasing = TRUE), ]
sorted_genes[1:20, c(1,2,28,32)]

up.reg   <- which(sorted_genes$log2FoldChange > 0)
down.reg <- which(sorted_genes$log2FoldChange < 0)
dim(sorted_genes); dim(sorted_genes[up.reg,]); dim(sorted_genes[down.reg,])

gg1 <- sorted_genes$gene_symbol; head(gg1,10); length(gg1)

dd2 <- gene_expr_data[gene_expr_data$gene_symbol %in% gg1, ]
dim(dd2); names(dd2)

sel1 <- sorted_metadata$sample_ID
dd3  <- dd2[, sel1]; names(dd3)
identical(colnames(dd3), as.character(sorted_metadata[,1]))

x4 <- as.matrix(dd3); dim(x4); range(x4)

## Reuse the same heatmap pipeline
library(gplots)

row.names <- dd2$gene_symbol
col.names <- sorted_metadata$group

xx5 <- t(scale(t(x4))); xx5[xx5 > 2] <- 2; xx6 <- xx5; xx6[xx6 < -2] <- -2; x4 <- xx6

my_palette  <- bluered(75)
mydatascale <- t(scale(t(x4)))
hr <- hclust(as.dist(1 - cor(t(mydatascale), method = "pearson")), method = "complete")
hc <- hclust(as.dist(1 - cor(mydatascale,        method = "spearman")), method = "complete")

cut_val1 <- 1.5
mycl <- cutree(hr, h = max(hr$height / cut_val1))
clusterCols      <- rainbow(length(unique(mycl)))
myClusterSideBar <- clusterCols[mycl]

cc1  <- as.numeric(as.factor(col.names))
from <- unique(cc1); to <- c("red", "orange", "darkgreen", "black")
map  <- setNames(to, from)
cc1[] <- map[cc1]
col_col1 <- cc1

par(mfrow = c(1,1), mar = c(4,3,3,2))
hv <- heatmap.2(
  x4, col = my_palette,
  margins = c(9,4),
  trace = "none",
  Colv  = NA,
  colsep = c(5,10,14), sepcolor = "black", sepwidth = c(0.1, 0.1),
  ColSideColors = col_col1,
  dendrogram = "row",
  labRow = row.names,
  labCol = col.names,
  cexRow = 0.2, cexCol = 0.8,
  main = paste(title)
)

col2  <- unique(col_col1)
nam10 <- unique(col.names)
legend("topright", legend = nam10, pch = 15, cex = 1.1, pt.cex = 1.6, bty = "n",
       ncol = 2, col = col2, inset = c(0.001, 0.03))


## ===========================================================================
## Contrast 3 (all_deseq[3])
## ===========================================================================

all_deseq
comp_file <- all_deseq[3]; comp_file
title     <- gsub("_080525.txt", "", comp_file); title

all_genes <- read.table(comp_file, header = TRUE)
names(all_genes); head(all_genes)

filtered_genes <- subset(all_genes, all_genes$padj < lim.P & abs(all_genes$log2FoldChange) > lim.FC)
dim(all_genes); dim(filtered_genes)

sorted_genes <- filtered_genes[order(filtered_genes$log2FoldChange, decreasing = TRUE), ]
sorted_genes[1:20, c(1,2,28,32)]
up.reg   <- which(sorted_genes$log2FoldChange > 0)
down.reg <- which(sorted_genes$log2FoldChange < 0)

gg1 <- sorted_genes$gene_symbol
dd2 <- gene_expr_data[gene_expr_data$gene_symbol %in% gg1, ]
sel1 <- sorted_metadata$sample_ID
dd3  <- dd2[, sel1]
x4   <- as.matrix(dd3)

## Heatmap
library(gplots)
row.names <- dd2$gene_symbol
col.names <- sorted_metadata$group

xx5 <- t(scale(t(x4))); xx5[xx5 > 2] <- 2; xx6 <- xx5; xx6[xx6 < -2] <- -2; x4 <- xx6

my_palette  <- bluered(75)
mydatascale <- t(scale(t(x4)))
hr <- hclust(as.dist(1 - cor(t(mydatascale), method = "pearson")), method = "complete")
hc <- hclust(as.dist(1 - cor(mydatascale,        method = "spearman")), method = "complete")

cut_val1 <- 1.5
mycl <- cutree(hr, h = max(hr$height / cut_val1))
clusterCols      <- rainbow(length(unique(mycl)))
myClusterSideBar <- clusterCols[mycl]

cc1  <- as.numeric(as.factor(col.names))
from <- unique(cc1); to <- c("red", "orange", "darkgreen", "black")
map  <- setNames(to, from)
cc1[] <- map[cc1]
col_col1 <- cc1

par(mfrow = c(1,1), mar = c(4,3,3,2))
hv <- heatmap.2(
  x4, col = my_palette,
  margins = c(9,4),
  trace = "none",
  Colv  = NA,
  colsep = c(5,10,14), sepcolor = "black", sepwidth = c(0.1, 0.1),
  ColSideColors = col_col1,
  dendrogram = "row",
  labRow = row.names,
  labCol = col.names,
  cexRow = 0.2, cexCol = 0.8,
  main = paste(title)
)

col2  <- unique(col_col1)
nam10 <- unique(col.names)
legend("topright", legend = nam10, pch = 15, cex = 1.1, pt.cex = 1.6, bty = "n",
       ncol = 2, col = col2, inset = c(0.001, 0.03))


## ===========================================================================
## Contrast 4 (all_deseq[4])
## ===========================================================================

all_deseq
comp_file <- all_deseq[4]; comp_file
title     <- gsub("_080525.txt", "", comp_file); title

all_genes <- read.table(comp_file, header = TRUE)
filtered_genes <- subset(all_genes, all_genes$padj < lim.P & abs(all_genes$log2FoldChange) > lim.FC)
sorted_genes <- filtered_genes[order(filtered_genes$log2FoldChange, decreasing = TRUE), ]

gg1 <- sorted_genes$gene_symbol
dd2 <- gene_expr_data[gene_expr_data$gene_symbol %in% gg1, ]
sel1 <- sorted_metadata$sample_ID
dd3  <- dd2[, sel1]
x4   <- as.matrix(dd3)

library(gplots)
row.names <- dd2$gene_symbol
col.names <- sorted_metadata$group

xx5 <- t(scale(t(x4))); xx5[xx5 > 2] <- 2; xx6 <- xx5; xx6[xx6 < -2] <- -2; x4 <- xx6

my_palette  <- bluered(75)
mydatascale <- t(scale(t(x4)))
hr <- hclust(as.dist(1 - cor(t(mydatascale), method = "pearson")), method = "complete")
hc <- hclust(as.dist(1 - cor(mydatascale,        method = "spearman")), method = "complete")

cut_val1 <- 1.5
mycl <- cutree(hr, h = max(hr$height / cut_val1))
clusterCols      <- rainbow(length(unique(mycl)))
myClusterSideBar <- clusterCols[mycl]

cc1  <- as.numeric(as.factor(col.names))
from <- unique(cc1); to <- c("red", "orange", "darkgreen", "black")
map  <- setNames(to, from)
cc1[] <- map[cc1]
col_col1 <- cc1

par(mfrow = c(1,1), mar = c(4,3,3,2))
hv <- heatmap.2(
  x4, col = my_palette,
  margins = c(9,4),
  trace = "none",
  Colv  = NA,
  colsep = c(5,10,14), sepcolor = "black", sepwidth = c(0.1, 0.1),
  ColSideColors = col_col1,
  dendrogram = "row",
  labRow = row.names,
  labCol = col.names,
  cexRow = 0.2, cexCol = 0.8,
  main = paste(title)
)

col2  <- unique(col_col1)
nam10 <- unique(col.names)
legend("topright", legend = nam10, pch = 15, cex = 1.1, pt.cex = 1.6, bty = "n",
       ncol = 2, col = col2, inset = c(0.001, 0.03))


## ===========================================================================
## Contrast 5 (all_deseq[5]) — CLP vs Sham (WT)
## ===========================================================================

all_deseq
comp_file <- all_deseq[5]; comp_file
title     <- gsub("_080525.txt", "", comp_file); title

all_genes <- read.table(comp_file, header = TRUE)
filtered_genes <- subset(all_genes, all_genes$padj < lim.P & abs(all_genes$log2FoldChange) > lim.FC)
sorted_genes <- filtered_genes[order(filtered_genes$log2FoldChange, decreasing = TRUE), ]

gg1 <- sorted_genes$gene_symbol
dd2 <- gene_expr_data[gene_expr_data$gene_symbol %in% gg1, ]
sel1 <- sorted_metadata$sample_ID
dd3  <- dd2[, sel1]
x4   <- as.matrix(dd3)

library(gplots)
row.names <- dd2$gene_symbol
col.names <- sorted_metadata$group

xx5 <- t(scale(t(x4))); xx5[xx5 > 2] <- 2; xx6 <- xx5; xx6[xx6 < -2] <- -2; x4 <- xx6

my_palette  <- bluered(75)
mydatascale <- t(scale(t(x4)))
hr <- hclust(as.dist(1 - cor(t(mydatascale), method = "pearson")), method = "complete")
hc <- hclust(as.dist(1 - cor(mydatascale,        method = "spearman")), method = "complete")

cut_val1 <- 1.5
mycl <- cutree(hr, h = max(hr$height / cut_val1))
clusterCols      <- rainbow(length(unique(mycl)))
myClusterSideBar <- clusterCols[mycl]

cc1  <- as.numeric(as.factor(col.names))
from <- unique(cc1); to <- c("red", "orange", "darkgreen", "black")
map  <- setNames(to, from)
cc1[] <- map[cc1]
col_col1 <- cc1

par(mfrow = c(1,1), mar = c(4,3,3,2))
hv <- heatmap.2(
  x4, col = my_palette,
  margins = c(9,4),
  trace = "none",
  Colv  = NA,
  colsep = c(5,10,14), sepcolor = "black", sepwidth = c(0.1, 0.1),
  ColSideColors = col_col1,
  dendrogram = "row",
  labRow = row.names,
  labCol = col.names,
  cexRow = 0.2, cexCol = 0.8,
  main = paste(title)
)

col2  <- unique(col_col1)
nam10 <- unique(col.names)
legend("topright", legend = nam10, pch = 15, cex = 1.1, pt.cex = 1.6, bty = "n",
       ncol = 2, col = col2, inset = c(0.001, 0.03))


## ---------------------------------------------------------------------------
## “Nicer” fixed-order heatmap examples (PDF outputs)
## ---------------------------------------------------------------------------

library(gplots); library(RColorBrewer)

pdf("heatmapAbCLPvsShamTRPV1_plot.pdf", width = 10, height = 8)

my_palette <- colorRampPalette(c("blue", "white", "red"))(255)

## Desired left-to-right group order
new_order   <- c("WT_CLP", "TRPV1_KO_CLP", "TRPV1_KO_Sham", "WT_Sham")
reorder_idx <- order(factor(col.names, levels = new_order))

x4        <- x4[, reorder_idx, drop = FALSE]
col.names <- col.names[reorder_idx]

## Group colors (match order)
group_colors <- c("red", "orange", "darkgreen", "black")
color_map    <- setNames(group_colors, new_order)
col_col1     <- color_map[col.names]

## Optional row cluster colors
cut_val1 <- 1.5
mycl <- cutree(hr, h = max(hr$height / cut_val1))
clusterCols      <- rainbow(length(unique(mycl)))
myClusterSideBar <- clusterCols[mycl]

heatmap.2(
  x4,
  col = my_palette,
  margins = c(9, 6),
  trace = "none",
  Colv = NA,
  dendrogram = "row",
  labRow = row.names,
  labCol = col.names,
  cexRow = 0.5,
  cexCol = 1.0,
  ColSideColors = col_col1,
  # RowSideColors = myClusterSideBar,
  colsep = c(4, 9, 14),
  sepcolor = "black",
  sepwidth = c(0.1, 0.1),
  main = "DEGs in sepsis - WT CLP vs WT Sham"
)

legend_labels <- c(
  "WT_CLP"        = "WT Sepsis",
  "TRPV1_KO_CLP"  = "KO Sepsis",
  "TRPV1_KO_Sham" = "KO Sham",
  "WT_Sham"       = "WT Sham"
)

legend("topright",
       legend = legend_labels,
       title = "groups",
       pch = 15,
       col = color_map,
       pt.cex = 2,
       cex = 1.0,
       bty = "n",
       inset = c(0.01, 0.0001),
       ncol = 2
)
dev.off()


## Same idea with automatic separators computed from counts
pdf("heatmapAbCLPvsShamTRPV1NEWORDER_plot.pdf", width = 10, height = 8)

my_palette <- colorRampPalette(c("blue", "white", "red"))(255)

new_order   <- c("WT_CLP", "TRPV1_KO_CLP", "TRPV1_KO_Sham", "WT_Sham")
reorder_idx <- order(factor(col.names, levels = new_order))
x4        <- x4[, reorder_idx, drop = FALSE]
col.names <- col.names[reorder_idx]

group_colors <- c("red", "orange", "darkgreen", "black")
color_map    <- setNames(group_colors, new_order)
col_col1     <- color_map[col.names]

cut_val1 <- 1.5
mycl <- cutree(hr, h = max(hr$height / cut_val1))
clusterCols      <- rainbow(length(unique(mycl)))
myClusterSideBar <- clusterCols[mycl]

## compute col separators between groups
group_counts <- as.integer(table(factor(col.names, levels = new_order)))
colsep <- cumsum(group_counts); if (length(colsep) > 1) colsep <- colsep[-length(colsep)] else colsep <- integer(0)

heatmap.2(
  x4,
  col = my_palette,
  margins = c(9, 6),
  trace = "none",
  Colv = NA,
  dendrogram = "row",
  labRow = row.names,
  labCol = col.names,
  cexRow = 0.5,
  cexCol = 1.0,
  ColSideColors = col_col1,
  # RowSideColors = myClusterSideBar,
  colsep = colsep,
  sepcolor = "black",
  sepwidth = c(0.1, 0.1),
  main = "DEGs in sepsis - WT CLP vs WT Sham"
)

legend("topright",
       legend = legend_labels[new_order],
       title  = "groups",
       pch = 15,
       col = color_map[new_order],
       pt.cex = 2,
       cex = 1.0,
       bty = "n",
       inset = c(0.01, 0.0001),
       ncol = 2
)
dev.off()


################################################################################
## Intersect/Custom gene set heatmap (from a saved list)
################################################################################

## sm points to a file with a gene list (e.g., intersect from a Venn step)
sm <- "intersect2_CLP_vs_CLP_sitag_260325.txt"; sm
nam_sm <- gsub(".txt", "", sm)

## Read gene list (first column)
ss1 <- read.table(sm[1], header = TRUE); head(ss1)
gg1 <- as.character(ss1[,1]); class(gg1); length(gg1); head(gg1)

## Subset expression matrix
dd2 <- gene_expr_data[gene_expr_data$gene_symbol %in% gg1, ]
dim(dd2); names(dd2)

## Order columns by sorted metadata
sel1 <- sorted_metadata$sample_ID; sel1
ordered_gene_expr <- dd2[, sel1]
identical(colnames(ordered_gene_expr), as.character(sorted_metadata[,1]))

ordered_gene_expr_matrix <- as.matrix(ordered_gene_expr); dim(ordered_gene_expr_matrix)
x4 <- ordered_gene_expr_matrix

## Standard gplots heatmap pipeline
library(gplots)
row.names <- dd2$gene_symbol
col.names <- sorted_metadata$group
table(col.names)

## Scale per gene and cap to [-2, 2]
scaled_matrix <- t(scale(t(x4))); range(scaled_matrix); hist(scaled_matrix)
scaled_matrix[scaled_matrix >  2] <-  2
final_scaled_matrix <- scaled_matrix
final_scaled_matrix[final_scaled_matrix < -2] <- -2
x4 <- final_scaled_matrix; range(x4); hist(x4)

## Colors & clustering
my_palette  <- bluered(75)
mydatascale <- t(scale(t(x4)))

hr <- hclust(as.dist(1 - cor(t(mydatascale), method = "pearson")),  method = "complete")
hc <- hclust(as.dist(1 - cor(   mydatascale, method = "spearman")), method = "complete")

cut_val1 <- 1.5
mycl <- cutree(hr, h = max(hr$height / cut_val1))
clusterCols      <- rainbow(length(unique(mycl)))
myClusterSideBar <- clusterCols[mycl]

## Group color bar
cc1  <- as.numeric(as.factor(col.names))
from <- unique(cc1)
to   <- c("red", "orange", "darkgreen", "black")
map  <- setNames(to, from)
cc1[] <- map[cc1]
col_col1 <- cc1

## Plot
par(mfrow = c(1,1), font = 2, font.axis = 2, font.lab = 3, mar = c(4,3,3,2))
hv <- heatmap.2(
  x4, col = my_palette,
  margins = c(9,4),
  trace = "none",
  Colv  = NA,
  colsep = c(5,10,15), sepcolor = "black", sepwidth = c(0.1, 0.1),
  ColSideColors = col_col1,
  # RowSideColors = myClusterSideBar,
  dendrogram = "row",
  labRow = row.names,
  labCol = col.names,
  cexRow = 0.2, cexCol = 0.8,
  main = paste(nam_sm)
)

col2  <- unique(col_col1)
nam10 <- unique(col.names)
legend("topright",
       legend = nam10,
       pch = 15, cex = 1, pt.cex = 1.5, bty = "n", ncol = 2,
       col = col2, inset = 0.01)

################################################################################
# End
################################################################################
