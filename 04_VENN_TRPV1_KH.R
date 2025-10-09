# ============================================================
# VENN DIAGRAM ANALYSIS SCRIPT — EXPLAINED -- xx KH
# TRPV1 RNA-seq dataset: DEG overlap visualization and extraction
# Notes:
#  • Object names, filenames, and logic unchanged.
#  • Code reorganized and commented for clarity to best of my knowledge - no guarantee
# ============================================================


# ============================================================
# 1) SETUP AND FILE DISCOVERY
# ============================================================

# Working directory: adjust to your machine if needed
setwd("/Users/katharinahollmann/GitHubDesktop/Wagner_comparative_analyses/TRPV1/data")
dir()

# List DESeq2 result files for TRPV1 contrasts
all_files <- dir(); all_files
deseq_files <- grep("DEG_DESeq2_TRPV1", all_files, value = TRUE); deseq_files  # list of DEG files
# Optional: drop non-DEG pathway files if present
# deseq_files <- deseq_files[-5]; deseq_files

# Expected order reference (info only):
# "DEG_DESeq2_TRPV1_TRPV1_KO_CLP_vs_TRPV1_KO_Sham_080525.txt"
# "DEG_DESeq2_TRPV1_TRPV1_KO_CLP_vs_WT_CLP_080525.txt"
# "DEG_DESeq2_TRPV1_TRPV1_KO_CLP_vs_WT_Sham_080525.txt"
# "DEG_DESeq2_TRPV1_TRPV1_KO_Sham_vs_WT_Sham_080525.txt"
# "DEG_DESeq2_TRPV1_WT_CLP_vs_WT_Sham_080525.txt"


# ============================================================
# 2) PER-CONTRAST FILTERING (TRPV1 ONLY)
#    Thresholds: padj < 0.05 and |log2FC| > log2(2)
#    Collect gene_symbol sets: genes1..genes5
# ============================================================

# 2.1 First contrast — TRPV1_KO_CLP vs TRPV1_KO_Sham
deseq_files  # check files
comp_file <- read.table("DEG_DESeq2_TRPV1_TRPV1_KO_CLP_vs_TRPV1_KO_Sham_080525.txt", header = TRUE)
names(comp_file); dim(comp_file); comp_file$log2FoldChange
range(abs(comp_file$log2FoldChange), na.rm = TRUE); range(comp_file$padj, na.rm = TRUE)

lim.P <- 0.05; lim.P
lim.FC <- log2(2); lim.FC

filtered_data <- subset(comp_file, comp_file$padj < lim.P & abs(comp_file$log2FoldChange) > lim.FC)
range(filtered_data$log2FoldChange, na.rm = TRUE); range(filtered_data$padj, na.rm = TRUE)
dim(filtered_data); dim(comp_file)

up.reg1   <- which(filtered_data$log2FoldChange > 0)
down.reg1 <- which(filtered_data$log2FoldChange < 0)
dim(filtered_data); length(up.reg1); length(down.reg1); names(filtered_data)

genes1 <- filtered_data$gene_symbol
head(genes1); length(genes1)


# 2.2 Second contrast — TRPV1_KO_CLP vs WT_CLP
deseq_files
comp_file <- read.table("DEG_DESeq2_TRPV1_TRPV1_KO_CLP_vs_WT_CLP_080525.txt", header = TRUE)
names(comp_file); dim(comp_file); comp_file$log2FoldChange
range(abs(comp_file$log2FoldChange), na.rm = TRUE); range(comp_file$padj, na.rm = TRUE)

lim.P <- 0.05; lim.P
lim.FC <- log2(2); lim.FC

filtered_data <- subset(comp_file, comp_file$padj < lim.P & abs(comp_file$log2FoldChange) > lim.FC)
range(filtered_data$log2FoldChange, na.rm = TRUE); range(filtered_data$padj, na.rm = TRUE)
dim(filtered_data); dim(comp_file)

up.reg2   <- which(filtered_data$log2FoldChange > 0)
down.reg2 <- which(filtered_data$log2FoldChange < 0)
dim(filtered_data); length(up.reg2); length(down.reg2); names(filtered_data)

genes2 <- filtered_data$gene_symbol
head(genes2); length(genes2)


# 2.3 Third contrast — TRPV1_KO_CLP vs WT_Sham
deseq_files
comp_file <- read.table("DEG_DESeq2_TRPV1_TRPV1_KO_CLP_vs_WT_Sham_080525.txt", header = TRUE)
names(comp_file); dim(comp_file); comp_file$log2FoldChange
range(abs(comp_file$log2FoldChange), na.rm = TRUE); range(comp_file$padj, na.rm = TRUE)

lim.P <- 0.05; lim.P
lim.FC <- log2(2); lim.FC

filtered_data3 <- subset(comp_file, comp_file$padj < lim.P & abs(comp_file$log2FoldChange) > lim.FC)
range(filtered_data3$log2FoldChange, na.rm = TRUE); range(filtered_data3$padj, na.rm = TRUE)
dim(filtered_data3); dim(comp_file)

up.reg3   <- which(filtered_data3$log2FoldChange > 0)
down.reg3 <- which(filtered_data3$log2FoldChange < 0)
dim(filtered_data3); length(up.reg3); length(down.reg3); names(filtered_data3)

genes3 <- filtered_data3$gene_symbol
head(genes3); length(genes3)


# 2.4 Fourth contrast — TRPV1_KO_Sham vs WT_Sham
deseq_files
comp_file <- read.table("DEG_DESeq2_TRPV1_TRPV1_KO_Sham_vs_WT_Sham_080525.txt", header = TRUE)
names(comp_file); dim(comp_file); head(comp_file$log2FoldChange, 10)
range(abs(comp_file$log2FoldChange), na.rm = TRUE); range(comp_file$padj, na.rm = TRUE)

lim.P <- 0.05; lim.P
lim.FC <- log2(2); lim.FC

filtered_data <- subset(comp_file, comp_file$padj < lim.P & abs(comp_file$log2FoldChange) > lim.FC)
range(filtered_data$log2FoldChange, na.rm = TRUE); range(filtered_data$padj, na.rm = TRUE)
dim(filtered_data); dim(comp_file)

up.reg4   <- which(filtered_data$log2FoldChange > 0)
down.reg4 <- which(filtered_data$log2FoldChange < 0)
dim(filtered_data); length(up.reg4); length(down.reg4); names(filtered_data)

genes4 <- filtered_data$gene_symbol
head(genes4); length(genes4)


# 2.5 Fifth contrast — WT_CLP vs WT_Sham
deseq_files
comp_file <- read.table("DEG_DESeq2_TRPV1_WT_CLP_vs_WT_Sham_080525.txt", header = TRUE)
names(comp_file); dim(comp_file); head(comp_file$log2FoldChange, 10)
range(abs(comp_file$log2FoldChange), na.rm = TRUE); range(comp_file$padj, na.rm = TRUE)

lim.P <- 0.05; lim.P
lim.FC <- log2(2); lim.FC

filtered_data <- subset(comp_file, comp_file$padj < lim.P & abs(comp_file$log2FoldChange) > lim.FC)
range(filtered_data$log2FoldChange, na.rm = TRUE); range(filtered_data$padj, na.rm = TRUE)
dim(filtered_data); dim(comp_file)

up.reg5   <- which(filtered_data$log2FoldChange > 0)
down.reg5 <- which(filtered_data$log2FoldChange < 0)
dim(filtered_data); length(up.reg5); length(down.reg5); names(filtered_data)

genes5 <- filtered_data$gene_symbol
head(genes5); length(genes5)


# ============================================================
# 3) VENN DIAGRAMS — BASE R overLapper + vennPlot (TRPV1)
#    First Venn selection: choose contrasts of interest
#    Below example: 1, 2, 3 (KO_CLP_vs_KO_Sham, KO_CLP_vs_WT_CLP, WT_CLP_vs_WT_Sham)
# ============================================================

# Combine and de-duplicate genes used in first Venn
all1 <- c(as.character(genes1), as.character(genes2), as.character(genes3))
length(all1)
all2 <- unique(all1)
length(all2)

# Keep counts for sanity
length(genes1); length(genes2); length(genes3)
length(all1); length(all2)

# Load overLapper/vennPlot utility
source("/Users/katharinahollmann/GitHubDesktop/Projects/Wagner_comparative_analyses/TRPV1/function/script_VENN.R")

# Set list for Venn (labelled by human-readable comparison)
deseq_files
setlist <- list(
  "KO CLP vs KO Sham" = genes1,
  "KO CLP vs WT CLP"  = genes2,
  "WT CLP vs WT Sham" = genes3
)

# Quick counts per set for reporting
cat("TRPV1_KO_CLP_vs_TRPV1_KO_Sham:", length(genes1), "UP:", length(up.reg1), "DOWN:", length(down.reg1), "\n",
    "TRPV1_KO_CLP_vs_WT_CLP:",      length(genes2), "UP:", length(up.reg2), "DOWN:", length(down.reg2), "\n",
    "TRPV1_WT_CLP_vs_WT_Sham:",     length(genes3), "UP:", length(up.reg3), "DOWN:", length(down.reg3), "\n")

# Compute venn sets and plot
str(setlist)
OLlist <- overLapper(setlist = setlist, sep = "", type = "vennsets")
OLlist; names(OLlist)
OLlist$Venn_List
counts <- sapply(OLlist$Venn_List, length); counts

par(mfrow = c(1,1), font = 2, font.axis = 2, font.lab = 3, mar = c(4,4,3,6))
vennPlot(counts = counts, mymain = "Venn DEG", lcex = 0.8, ccex = 1)
# PDF guidance: landscape 7x6
str(OLlist); names(OLlist$Venn_List)


# ============================================================
# 4) VENN DIAGRAMS — ggVennDiagram (TRPV1)
#    Alternate visualization with fill gradients and cleaner styling
#    Example selection: genes3, genes2, genes5
# ============================================================

library(ggplot2)
library(ggVennDiagram)

setlist <- list(
  "KO CLP vs WT Sham" = genes3,
  "KO CLP vs WT CLP"  = genes2,
  "WT CLP vs WT Sham" = genes5
)

p <- ggVennDiagram(setlist,
                   label = "count",
                   label_alpha = 0) +
  scale_fill_gradient(low = "#F4FAFE", high = "#7CAEDC") +
  theme_void() +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, size = 12),
    plot.margin = margin(t = 60, r = 160, b = 60, l = 250)
  ) +
  coord_fixed(ratio = 1, clip = "off") +
  ggtitle("Venn DEG")

p
ggsave("VennplotKOmargin.pdf", plot = p, width = 9.5, height = 6)
ggsave("venn_deg.png", p, width = 8, height = 6, units = "in", dpi = 300)


# ============================================================
# 5) OVERLAP LISTS AND ANNOTATION MERGE (OPTIONAL REPORTING)
#    Example: collect all CLP vs Sham overlaps across KO and WT
# ============================================================

# Example two-set overlap using OLlist from section 3:
overlapCLPSham <- OLlist$Venn_List$`KO CLP vs KO ShamWT CLP vs WT Sham`
overlapCLPSham

allCLPSham <- c(
  OLlist$Venn_List$`KO CLP vs KO Sham`,
  OLlist$Venn_List$`WT CLP vs WT Sham`,
  OLlist$Venn_List$`KO CLP vs KO ShamWT CLP vs WT Sham`
)
length(allCLPSham)

# Convert to data.frame and deduplicate gene symbols
overlap_df <- data.frame(gene_symbol = allCLPSham)
head(overlap_df); dim(overlap_df)

# Use dplyr for distinct and merges
library(dplyr)

overlap_df_unique <- overlap_df %>% distinct(gene_symbol)
dim(overlap_df_unique)

# Merge KO_CLP_vs_KO_Sham and WT_CLP_vs_WT_Sham filtered tables to bring back entrez and stats
filteredKOKO <- filtered_data %>% select(gene_symbol, entrezgene_id, log2FoldChange, padj)  # from 2.4 block (last filtered_data)
dim(filteredKOKO)

filteredWTWT <- filtered_data3 %>% select(gene_symbol, entrezgene_id, log2FoldChange, padj) # from 2.3 block
dim(filteredWTWT)

# Union and bind_rows as in original script
mergefiltered <- union(filteredKOKO, filteredWTWT)
mergefiltered <- bind_rows(filteredKOKO, filteredWTWT)
head(mergefiltered); dim(mergefiltered)

# Join overlap set with annotations
matched_genesALL <- inner_join(
  overlap_df, mergefiltered[, c("gene_symbol", "entrezgene_id", "log2FoldChange", "padj")],
  by = "gene_symbol"
)
head(matched_genesALL); dim(matched_genesALL); names(matched_genesALL)

# Remove duplicates by gene_symbol
matched_genesALL$gene_symbol <- as.character(matched_genesALL$gene_symbol)
matchedgenesALL_nodup <- distinct(matched_genesALL, gene_symbol, .keep_all = TRUE)
dim(matchedgenesALL_nodup); head(matchedgenesALL_nodup); names(matchedgenesALL_nodup); class(matchedgenesALL_nodup)

# Save annotated overlap table
write.table(matchedgenesALL_nodup, file = "TRPV1DEGCLPShamWTandKO_260625.txt")


# ============================================================
# 6) SECOND VENN: KO CLP vs WT CLP vs CLP vs Sham (KO|WT pooled)
#    Build combined set (CLP vs Sham KO + WT) and compare with KO vs WT under CLP
# ============================================================

setlist <- list(
  "KO CLP vs WT CLP"   = genes2,
  "CLP vs Sham KO|WT"  = allCLPSham # genes modified in sepsis independent of genotype 
)

str(setlist)
OLlist <- overLapper(setlist = setlist, sep = "", type = "vennsets")
OLlist; names(OLlist)
OLlist$Venn_List
counts <- sapply(OLlist$Venn_List, length); counts

par(mfrow = c(1,1), font = 2, font.axis = 2, font.lab = 3, mar = c(4,4,3,6))
vennPlot(counts = counts, mymain = "Venn DEG", lcex = 0.8, ccex = 1)

# Extract overlap genes across these two sets
DEGregulated <- OLlist$Venn_List$`KO CLP vs WT CLPCLP vs Sham KO|WT`
length(DEGregulated)

overlap_deg <- data.frame(gene_symbol = DEGregulated)
head(overlap_deg); dim(overlap_deg)

# Merge back effect sizes and IDs
overlap_deg$gene_symbol <- as.character(overlap_deg$gene_symbol)
matchedgenesALL_nodup$gene_symbol <- as.character(matchedgenesALL_nodup$gene_symbol)

sum(overlap_deg$gene_symbol %in% matchedgenesALL_nodup$gene_symbol)  # Should be 375

DEGafterinfection <- merge(
  overlap_deg,
  matchedgenesALL_nodup[, c("gene_symbol", "entrezgene_id", "log2FoldChange", "padj")],
  by = "gene_symbol"
)
dim(DEGafterinfection); head(DEGafterinfection) 

# i called them DEG after infection -> those are the genes which are modified by trpv1 knockout AND altered in sepsis 
# for promotion and data visualization i recreated venn just using smart art in microsoft office 

# Save downstream list for profiler use
write.table(DEGafterinfection, file = "TRPV1DEGafterinfection_260625.txt")


# ============================================================
# 7) VENN BY DIRECTION (UP/DOWN) — TRPV1 SET A (1,3,5)
#    Up/Down Venn for KO_CLP_vs_KO_Sham, KO_CLP_vs_WT_Sham, WT_CLP_vs_WT_Sham
# ============================================================

# Confirm counts by direction
length(up.reg1); length(up.reg3); length(up.reg5)
length(down.reg1); length(down.reg3); length(down.reg5)
deseq_files

# Upregulated sets
setlist_up <- list(
  "TRPV1_KO_CLP_vs_TRPV1_KO_Sham" = genes1[up.reg1],
  "TRPV1_KO_CLP_vs_WT_Sham"       = genes3[up.reg3],
  "TRPV1_WT_CLP_vs_WT_Sham"       = genes5[up.reg5]
)
OLlist_up  <- overLapper(setlist = setlist_up, sep = "", type = "vennsets")
counts_up  <- sapply(OLlist_up$Venn_List, length)
vennPlot(counts = counts_up, mymain = "Venn DEG - Upregulated", lcex = 0.8, ccex = 1)

# Downregulated sets
setlist_down <- list(
  "TRPV1_KO_CLP_vs_TRPV1_KO_Sham" = genes1[down.reg1],
  "TRPV1_KO_CLP_vs_WT_Sham"       = genes3[down.reg3],
  "TRPV1_WT_CLP_vs_WT_Sham"       = genes5[down.reg5]
)
OLlist_down <- overLapper(setlist = setlist_down, sep = "", type = "vennsets")
counts_down <- sapply(OLlist_down$Venn_List, length)
vennPlot(counts = counts_down, mymain = "Venn DEG - Downregulated", lcex = 0.8, ccex = 1)


# ============================================================
# 8) VENN BY DIRECTION (UP/DOWN) — TRPV1 SET B (2,3,4)
#    Up/Down Venn for TRPV1_KO_CLP_vs_WT_CLP, TRPV1_KO_CLP_vs_WT_Sham, TRPV1_KO_Sham_vs_WT_Sham
# ============================================================

length(up.reg2); length(up.reg3); length(up.reg4)
length(down.reg2); length(down.reg3); length(down.reg4)
deseq_files

# Upregulated sets
setlist_up <- list(
  "TRPV1_KO_CLP_vs_WT_CLP"   = genes1[up.reg2],
  "TRPV1_KO_CLP_vs_WT_Sham"  = genes3[up.reg3],
  "TRPV1_KO_Sham_vs_WT_Sham" = genes5[up.reg4]
)
OLlist_up  <- overLapper(setlist = setlist_up, sep = "", type = "vennsets")
counts_up  <- sapply(OLlist_up$Venn_List, length)
vennPlot(counts = counts_up, mymain = "Venn DEG - Upregulated", lcex = 0.8, ccex = 1)

# Downregulated sets
setlist_down <- list(
  "TRPV1_KO_CLP_vs_WT_CLP"   = genes1[down.reg2],
  "TRPV1_KO_CLP_vs_WT_Sham"  = genes3[down.reg3],
  "TRPV1_KO_Sham_vs_WT_Sham" = genes5[down.reg4]
)
OLlist_down <- overLapper(setlist = setlist_down, sep = "", type = "vennsets")
counts_down <- sapply(OLlist_down$Venn_List, length)
vennPlot(counts = counts_down, mymain = "Venn DEG - Downregulated", lcex = 0.8, ccex = 1)


# ============================================================
# END — TRPV1 VENN WORKFLOW ONLY (antibody/“more analyses” removed)
# ============================================================