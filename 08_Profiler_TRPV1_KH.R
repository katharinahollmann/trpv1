################################################################################################
# FUNCTIONAL ANALYSIS SCRIPT — EXPLAINED — KH
# TRPV1 RNA-seq | GO enrichment with clusterProfiler
#
# Notes on key objects (preserved):
# - file_list, deseq_files, deseq_txt : file discovery
# - comp_file : the chosen DESeq2 results file for a contrast
# - title : human-readable title derived from comp_file
# - dat : full DESeq2 table for a contrast
# - lim.P, lim.FC : thresholds for padj and |log2FC|
# - filtered_dat : DEGs (padj < lim.P & |log2FC| > lim.FC)
# - filtered_datUP / filtered_datDN : split by direction
# - up / dn : index vectors (sanity check sizes only)
# - background : full gene list (as Entrez) for universe (where used)
# - geneList1 : background Entrez IDs
# - genes_for_pathways_all / gene3 : query Entrez IDs for enrichment
# - ego : enrichGO result
# - dp : dotplot ggplot object
# - res1 : data.frame of enrichment results (for optional write.table)
################################################################################################


##############################################
# 0) SETUP: working directory & packages
##############################################

# Intro, Working Directory ------------------------------------------------
# -- adjust to your system !! this is absolute path (for now)
setwd("/Users/katharinahollmann/GitHubDesktop/Wagner_comparative_analyses/TRPV1/data")
dir()

# Libraries ---------------------------------------------------------------
library(org.Mm.eg.db)
library(clusterProfiler)
library(IRanges)
# ATTENTION: KEGG.db is old (now proprietary). Using GO via org.Mm.eg.db.
# sessionInfo()


##############################################
# 1) LOAD ALL FILES & BUILD CANDIDATE LISTS
##############################################

# Load All Files ----------------------------------------------------------
file_list <- dir()
file_list

deseq_files <- grep("TRPV1", file_list, value = TRUE)
deseq_files

deseq_txt <- grep(".txt", deseq_files, value = TRUE)
deseq_txt   # list with all DESeq2 txt files
# (Removed "DESeq2" from first grep — overlap files are at indices 9 and 10 in your current folder)


################################################################################################
# 2) GLOBAL ENRICHMENT BLOCK — FIRST SELECTION (comp_file <- deseq_txt[9])
################################################################################################

#---------------------------------------------
# Choose file (FIRST)
comp_file <- deseq_txt[9]
comp_file

title <- gsub("_260625.txt", "", comp_file)
title

# Read full table & define thresholds -------------------------------------
dat <- read.table(comp_file, header = TRUE)
dim(dat)
names(dat)

lim.P  <- 0.05
lim.FC <- log2(2)

# Define DEG sets ----------------------------------------------------------
filtered_dat <- subset(
  dat,
  dat$padj < lim.P & abs(dat$log2FoldChange) > lim.FC
)
dim(filtered_dat)

filtered_datUP <- subset(
  dat,
  dat$log2FoldChange > lim.FC & dat$padj < lim.P
)
filtered_datUP <- droplevels(filtered_datUP)

filtered_datDN <- subset(
  dat,
  dat$log2FoldChange < -lim.FC & dat$padj < lim.P
)
filtered_datDN <- droplevels(filtered_datDN)

up <- which(filtered_datUP$log2FoldChange > 0)
dn <- which(filtered_datDN$log2FoldChange < 0)

dim(filtered_dat)
length(up)
length(dn)

# Background / Universe (from another file) --------------------------------
background <- deseq_txt[3]
background <- read.table(background, header = TRUE)
dim(background)
names(background)

geneList1 <- as.character(background$entrezgene_id)
head(geneList1, 10)   # ALL genes as background

# ALL DEGs (both directions) as query --------------------------------------
genes_for_pathways_all <- as.character(filtered_dat$entrezgene_id)
length(genes_for_pathways_all)
head(genes_for_pathways_all)

# GO enrichment (BP) & simple barplot --------------------------------------
ego <- enrichGO(
  gene          = genes_for_pathways_all,
  universe      = geneList1,
  OrgDb         = "org.Mm.eg.db",
  ont           = "BP",
  pvalueCutoff  = 0.05,
  readable      = TRUE
)
comp_file

barplot(ego, showCategory = 10)

# Dotplot (ggplot) ---------------------------------------------------------
library(ggplot2)

title_text <- "TRPV1 // "

dp <- dotplot(
  ego,
  showCategory = 10,
  title       = paste(title_text, "DEG after infection"),
  font.size   = 9
)

dp + theme(
  axis.text.y = element_text(size = 16),
  plot.title  = element_text(size = 14, face = "bold")
)

# Collect results table (optional write) -----------------------------------
dim(filtered_dat)
names(filtered_dat)

res1 <- as.data.frame(ego)
head(res1)
names(res1)
comp_file
# write.table(
#   res1,
#   file = "pathwaysTRPV1DEGafterinfection_260625.txt",
#   sep = "\t", quote = FALSE, row.names = FALSE
# )


# UP only ------------------------------------------------------------------
names(dat)
names(filtered_datUP)
length(up)

geneList1 <- as.character(background$entrezgene_id)

gene3 <- as.character(filtered_datUP$entrezgene_id)
head(gene3, 10)
length(gene3)

ego <- enrichGO(
  gene          = gene3,
  universe      = geneList1,
  OrgDb         = "org.Mm.eg.db",
  ont           = "BP",
  pvalueCutoff  = 0.05,
  readable      = TRUE
)
comp_file

library(ggplot2)

title_text <- "TRPV1// "

dp <- dotplot(
  ego,
  showCategory = 10,
  title       = paste(title_text, "DEG after infection UP"),
  font.size   = 9
)

dp + theme(
  axis.text.y = element_text(size = 16),
  plot.title  = element_text(size = 14, face = "bold")
)

dotplot(
  ego,
  showCategory = 30,
  title       = paste("EnrichGO //", title),
  font.size   = 9
)

dim(filtered_datUP)
names(filtered_datUP)

res1 <- as.data.frame(ego)
head(res1)
names(res1)
comp_file
# write.table(
#   res1,
#   file = "pathwaysTRPV1_UP_DEGafterinfection_260625.txt",
#   sep = "\t", quote = FALSE, row.names = FALSE
# )


# DOWN only ----------------------------------------------------------------
library(org.Mm.eg.db)
library(clusterProfiler)

names(dat)
names(filtered_datDN)
length(dn)

geneList1 <- as.character(background$entrezgene_id)

gene3 <- as.character(filtered_datDN$entrezgene_id)
head(gene3, 10)
length(gene3)

ego <- enrichGO(
  gene          = gene3,
  universe      = geneList1,
  OrgDb         = "org.Mm.eg.db",
  ont           = "BP",
  pvalueCutoff  = 0.05,
  readable      = TRUE
)
comp_file

dotplot(
  ego,
  showCategory = 10,
  title       = paste("EnrichGO //", title),
  font.size   = 9
)

library(ggplot2)

title_text <- "TRPV1// "

dp <- dotplot(
  ego,
  showCategory = 10,
  title       = paste(title_text, "DEG after infection DN"),
  font.size   = 9
)

dp + theme(
  axis.text.y = element_text(size = 16),
  plot.title  = element_text(size = 14, face = "bold")
)

dim(filtered_datDN)
names(filtered_datDN)

res1 <- as.data.frame(ego)
head(res1)
names(res1)
comp_file
# write.table(
#   res1,
#   file = "pathways_down_KO_CLP_vs_KO_Sham.txt",
#   sep = "\t", quote = FALSE, row.names = FALSE
# )



################################################################################################
# 3) INDIVIDUAL CONTRASTS — FIVE BLOCKS (FIRST ... FIFTH)
#    Each block repeats the same logic with comp_file <- deseq_txt[i].
################################################################################################

deseq_txt   # view available indices


# ============================================================================================
# FIRST CONTRAST (example): comp_file <- deseq_txt[2]
# ============================================================================================

comp_file <- deseq_txt[2]
comp_file

title <- gsub("_080525.txt", "", comp_file)
title

dat <- read.table(comp_file, header = TRUE)
dim(dat)
names(dat)

lim.P  <- 0.05
lim.FC <- log2(2)

filtered_dat <- subset(
  dat,
  dat$padj < lim.P & abs(dat$log2FoldChange) > lim.FC
)
dim(filtered_dat)

filtered_datUP <- subset(
  dat,
  dat$log2FoldChange > lim.FC & dat$padj < lim.P
)
filtered_datUP <- droplevels(filtered_datUP)

filtered_datDN <- subset(
  dat,
  dat$log2FoldChange < -lim.FC & dat$padj < lim.P
)
filtered_datDN <- droplevels(filtered_datDN)

up <- which(filtered_datUP$log2FoldChange > 0)
dn <- which(filtered_datDN$log2FoldChange < 0)

dim(filtered_dat)
length(up)
length(dn)

# ALL (universe = all genes in this 'dat') ---------------------------------
geneList1 <- as.character(dat$entrezgene_id)

genes_for_pathways_all <- as.character(filtered_dat$entrezgene_id)
length(genes_for_pathways_all)
head(genes_for_pathways_all)

ego <- enrichGO(
  gene          = genes_for_pathways_all,
  universe      = geneList1,
  OrgDb         = "org.Mm.eg.db",
  ont           = "BP",
  pvalueCutoff  = 0.05,
  readable      = TRUE
)
comp_file

library(ggplot2)

title_text <- "EnrichGO // "

dp <- dotplot(
  ego,
  showCategory = 10,
  title       = paste(title_text, "KO CLP vs KO Sham"),
  font.size   = 9
)

dp + theme(
  axis.text.y = element_text(size = 16),
  plot.title  = element_text(size = 14, face = "bold")
)

dim(filtered_dat)
names(filtered_dat)

res1 <- as.data.frame(ego)
head(res1)
names(res1)
comp_file
# write.table(
#   res1,
#   file = "pathwaysALL_KO_CLP_vs_Sham_160525.txt",
#   sep = "\t", quote = FALSE, row.names = FALSE
# )


# UP ======================================================================
geneList1 <- as.character(dat$entrezgene_id)

gene3 <- as.character(filtered_datUP$entrezgene_id)
head(gene3, 10)
length(gene3)

ego <- enrichGO(
  gene          = gene3,
  universe      = geneList1,
  OrgDb         = "org.Mm.eg.db",
  ont           = "BP",
  pvalueCutoff  = 0.05,
  readable      = TRUE
)
comp_file

dim(filtered_datUP)
names(filtered_datUP)

res1 <- as.data.frame(ego)
head(res1)
names(res1)
comp_file

### prettier ---------------------------------------------------------------
library(ggplot2)

dp <- dotplot(
  ego,
  showCategory = 10,
  font.size   = 40,
  title       = "KO CLP vs KO Sham up"
)

dp <- dp + 
  scale_size(range = c(6, 20)) +
  theme(
    axis.text.y  = element_text(size = 50),
    legend.text  = element_text(size = 20),
    legend.title = element_text(size = 22),
    plot.title   = element_text(size = 30, face = "bold")
  )

comp_file

ggsave(
  "enrichGO_dotplotUPKOCLPvsKOSham.pdf",
  plot = dp, width = 16.5, height = 19
)

# write.table(
#   res1,
#   file = "pathways_up_TRPV1_KO_CLP_vs_TRPV1_KO_Sham.txt",
#   sep = "\t", quote = FALSE, row.names = FALSE
# )

# cnet (UP) ---------------------------------------------------------------
gl1 <- filtered_datUP$log2FoldChange
names(gl1) <- filtered_datUP$entrezgene_id

geneList10 <- gl1

cnetplot(
  ego,
  categorySize       = "pvalue",
  foldChange         = geneList10,
  showCategory       = 10,
  max.overlaps       = 800,
  cex_category       = 0.8,
  cex_gene           = 0.5,
  cex_label_category = 0.8,
  cex_label_gene     = 0.6
)

# additional views --------------------------------------------------------
edox <- ego

p1 <- cnetplot(
  edox,
  node_label        = "category",
  showCategory      = 20,
  max.overlaps      = 100,
  cex_label_category= 0.8,
  foldChange        = geneList10
)

p2 <- cnetplot(
  edox,
  node_label    = "gene",
  cex_label_gene= 1,
  showCategory  = 20,
  max.overlaps  = 100
)

cowplot::plot_grid(p1, p2, ncol = 2, labels = LETTERS[1:2])


# DOWN  ===================================================================
library(org.Mm.eg.db)
library(clusterProfiler)

geneList1 <- as.character(dat$entrezgene_id)

gene3 <- as.character(filtered_datDN$entrezgene_id)
head(gene3, 10)
length(gene3)

ego <- enrichGO(
  gene          = gene3,
  universe      = geneList1,
  OrgDb         = "org.Mm.eg.db",
  ont           = "BP",
  pvalueCutoff  = 0.05,
  readable      = TRUE
)
comp_file

dim(filtered_datDN)
names(filtered_datDN)

res1 <- as.data.frame(ego)
head(res1)
names(res1)
comp_file

### prettier ---------------------------------------------------------------
dp <- dotplot(
  ego,
  showCategory = 10,
  font.size   = 40,
  title       = "KO CLP vs KO Sham down"
)

dp <- dp +
  scale_size(range = c(6, 20)) +
  theme(
    axis.text.y  = element_text(size = 50),
    legend.text  = element_text(size = 20),
    legend.title = element_text(size = 22),
    plot.title   = element_text(size = 30, face = "bold")
  )

ggsave(
  "enrichGO_dotplotDOWNKOCLPvsKOSham.pdf",
  plot = dp, width = 16.5, height = 19
)

# write.table(
#   res1,
#   file = "pathways_down_KO_CLP_vs_KO_Sham.txt",
#   sep = "\t", quote = FALSE, row.names = FALSE
# )

gl1 <- filtered_datDN$log2FoldChange
names(gl1) <- filtered_datDN$entrezgene_id

geneList10 <- gl1

cnetplot(
  ego,
  categorySize       = "pvalue",
  foldChange         = geneList10,
  showCategory       = 10,
  max.overlaps       = 800,
  cex_category       = 0.8,
  cex_gene           = 0.5,
  cex_label_category = 0.8,
  cex_label_gene     = 0.6
)

edox <- ego

p1 <- cnetplot(
  edox,
  node_label        = "category",
  showCategory      = 20,
  max.overlaps      = 100,
  cex_label_category= 0.8,
  foldChange        = geneList10
)

p2 <- cnetplot(
  edox,
  node_label    = "gene",
  cex_label_gene= 1,
  showCategory  = 20,
  max.overlaps  = 100
)

cowplot::plot_grid(p1, p2, ncol = 2, labels = LETTERS[1:2])

p1 <- heatplot(edox, showCategory = 5)
p1

p2 <- heatplot(edox, foldChange = geneList10, showCategory = 10)
p2



# ============================================================================================
# SECOND CONTRAST: comp_file <- deseq_txt[3]
# ============================================================================================

comp_file <- deseq_txt[3]
comp_file

title <- gsub("_080525.txt", "", comp_file)
title

dat <- read.table(comp_file, header = TRUE)
dim(dat)
names(dat)

lim.P  <- 0.05
lim.FC <- log2(2)

filtered_dat <- subset(
  dat,
  dat$padj < lim.P & abs(dat$log2FoldChange) > lim.FC
)
dim(filtered_dat)

filtered_datUP <- subset(
  dat,
  dat$log2FoldChange > lim.FC & dat$padj < lim.P
)
filtered_datUP <- droplevels(filtered_datUP)

filtered_datDN <- subset(
  dat,
  dat$log2FoldChange < -lim.FC & dat$padj < lim.P
)
filtered_datDN <- droplevels(filtered_datDN)

up <- which(filtered_datUP$log2FoldChange > 0)
dn <- which(filtered_datDN$log2FoldChange < 0)

dim(filtered_dat)
length(up)
length(dn)

# ALL ---------------------------------------------------------------------
geneList1 <- as.character(dat$entrezgene_id)

genes_for_pathways_all <- as.character(filtered_dat$entrezgene_id)

ego <- enrichGO(
  gene          = genes_for_pathways_all,
  universe      = geneList1,
  OrgDb         = "org.Mm.eg.db",
  ont           = "BP",
  pvalueCutoff  = 0.05,
  readable      = TRUE
)

res1 <- as.data.frame(ego)
head(res1)

library(ggplot2)

dp <- dotplot(
  ego,
  showCategory = 10,
  title       = "EnrichGO // KO CLP vs WT CLP",
  font.size   = 9
)

dp + theme(
  axis.text.y = element_text(size = 16),
  plot.title  = element_text(size = 14, face = "bold")
)

# write.table(
#   res1,
#   file = "pathwaysALL_KO_CLP_vs_WT_CLP_160525.txt",
#   sep = "\t", quote = FALSE, row.names = FALSE
# )

# UP ======================================================================
gene3 <- as.character(filtered_datUP$entrezgene_id)

ego <- enrichGO(
  gene          = gene3,
  universe      = geneList1,
  OrgDb         = "org.Mm.eg.db",
  ont           = "BP",
  pvalueCutoff  = 0.05,
  readable      = TRUE
)

res1 <- as.data.frame(ego)

### prettier ---------------------------------------------------------------
dp <- dotplot(
  ego,
  showCategory = 10,
  font.size   = 40,
  title       = "KO CLP vs WT CLP up"
)

dp <- dp +
  scale_size(range = c(6, 20)) +
  theme(
    axis.text.y  = element_text(size = 50),
    legend.text  = element_text(size = 20),
    legend.title = element_text(size = 22),
    plot.title   = element_text(size = 30, face = "bold")
  )

ggsave(
  "enrichGO_dotplotUPKOCLPvsWTCLP.pdf",
  plot = dp, width = 16, height = 19
)

write.table(
  res1,
  file = "pathways_up_TRPV1_KO_CLP_vs_WT_CLP.txt",
  sep = "\t", quote = FALSE, row.names = FALSE
)

dotplot(
  ego,
  showCategory = 30,
  title       = paste("EnrichGO //", title),
  font.size   = 9
)

gl1 <- filtered_datUP$log2FoldChange
names(gl1) <- filtered_datUP$entrezgene_id

geneList10 <- gl1

edox <- ego

p1 <- cnetplot(
  edox,
  node_label        = "category",
  showCategory      = 20,
  max.overlaps      = 100,
  cex_label_category= 0.8,
  foldChange        = geneList10
)

p2 <- cnetplot(
  edox,
  node_label    = "gene",
  cex_label_gene= 1,
  showCategory  = 20,
  max.overlaps  = 100
)

cowplot::plot_grid(p1, p2, ncol = 2, labels = LETTERS[1:2])

# DOWN ====================================================================
gene3 <- as.character(filtered_datDN$entrezgene_id)

ego <- enrichGO(
  gene          = gene3,
  universe      = geneList1,
  OrgDb         = "org.Mm.eg.db",
  ont           = "BP",
  pvalueCutoff  = 0.05,
  readable      = TRUE
)

### prettier ---------------------------------------------------------------
dp <- dotplot(
  ego,
  showCategory = 10,
  font.size   = 40,
  title       = "KO CLP vs WT Sham down"
)

dp <- dp +
  scale_size(range = c(6, 20)) +
  theme(
    axis.text.y  = element_text(size = 50),
    legend.text  = element_text(size = 20),
    legend.title = element_text(size = 22),
    plot.title   = element_text(size = 30, face = "bold")
  )

ggsave(
  "enrichGO_dotplotDOWNKOCLPvsWTSham.pdf",
  plot = dp, width = 13, height = 9
)

dp <- dotplot(
  ego,
  showCategory = 10,
  title       = "EnrichGO // KO CLP vs WT Sham - down",
  font.size   = 9
)

dp + theme(
  axis.text.y = element_text(size = 16),
  plot.title  = element_text(size = 14, face = "bold")
)
# (optional write suppressed)



# ============================================================================================
# THIRD CONTRAST: comp_file <- deseq_txt[4]
# ============================================================================================

comp_file <- deseq_txt[4]
comp_file

title <- gsub("_080525.txt", "", comp_file)
title

dat <- read.table(comp_file, header = TRUE)
dim(dat)
names(dat)

lim.P  <- 0.05
lim.FC <- log2(2)

filtered_dat <- subset(
  dat,
  dat$padj < lim.P & abs(dat$log2FoldChange) > lim.FC
)
dim(filtered_dat)

filtered_datUP <- subset(
  dat,
  dat$log2FoldChange > lim.FC & dat$padj < lim.P
)
filtered_datUP <- droplevels(filtered_datUP)

filtered_datDN <- subset(
  dat,
  dat$log2FoldChange < -lim.FC & dat$padj < lim.P
)
filtered_datDN <- droplevels(filtered_datDN)

up <- which(filtered_datUP$log2FoldChange > 0)
dn <- which(filtered_datDN$log2FoldChange < 0)

# ALL ---------------------------------------------------------------------
geneList1 <- as.character(dat$entrezgene_id)

genes_for_pathways_all <- as.character(filtered_dat$entrezgene_id)

ego <- enrichGO(
  gene          = genes_for_pathways_all,
  universe      = geneList1,
  OrgDb         = "org.Mm.eg.db",
  ont           = "BP",
  pvalueCutoff  = 0.05,
  readable      = TRUE
)

library(ggplot2)

dp <- dotplot(
  ego,
  showCategory = 10,
  title       = "EnrichGO // KO CLP vs WT Sham",
  font.size   = 9
)

dp + theme(
  axis.text.y = element_text(size = 16),
  plot.title  = element_text(size = 14, face = "bold")
)

res1 <- as.data.frame(ego)
head(res1)

# write.table(
#   res1,
#   file = "pathwaysALL_KO_CLP_WT_Sham_160525.txt",
#   sep = "\t", quote = FALSE, row.names = FALSE
# )

# UP ----------------------------------------------------------------------
gene3 <- as.character(filtered_datUP$entrezgene_id)

ego <- enrichGO(
  gene          = gene3,
  universe      = geneList1,
  OrgDb         = "org.Mm.eg.db",
  ont           = "BP",
  pvalueCutoff  = 0.05,
  readable      = TRUE
)

### prettier ---------------------------------------------------------------
dp <- dotplot(
  ego,
  showCategory = 10,
  font.size   = 40,
  title       = "KO CLP vs WT Sham up"
)

dp <- dp +
  scale_size(range = c(6, 20)) +
  theme(
    axis.text.y  = element_text(size = 50),
    legend.text  = element_text(size = 20),
    legend.title = element_text(size = 22),
    plot.title   = element_text(size = 30, face = "bold")
  )

ggsave(
  "enrichGO_dotplotUPKOCLPvsWTSham.pdf",
  plot = dp, width = 16.7, height = 19
)

dotplot(
  ego,
  showCategory = 30,
  title       = paste("EnrichGO //", title),
  font.size   = 9
)

res1 <- as.data.frame(ego)
head(res1)

write.table(
  res1,
  file = "pathways_up_TRPV1_KO_CLP_vs_WT_Sham.txt",
  sep = "\t", quote = FALSE, row.names = FALSE
)

gl1 <- filtered_datUP$log2FoldChange
names(gl1) <- filtered_datUP$entrezgene_id

geneList10 <- gl1

edox <- ego

p1 <- cnetplot(
  edox,
  node_label        = "category",
  showCategory      = 20,
  max.overlaps      = 100,
  cex_label_category= 0.8,
  foldChange        = geneList10
)

p2 <- cnetplot(
  edox,
  node_label    = "gene",
  cex_label_gene= 1,
  showCategory  = 20,
  max.overlaps  = 100
)

cowplot::plot_grid(p1, p2, ncol = 2, labels = LETTERS[1:2])

# DOWN --------------------------------------------------------------------
gene3 <- as.character(filtered_datDN$entrezgene_id)

ego <- enrichGO(
  gene          = gene3,
  universe      = geneList1,
  OrgDb         = "org.Mm.eg.db",
  ont           = "BP",
  pvalueCutoff  = 0.05,
  readable      = TRUE
)
comp_file

### prettier ---------------------------------------------------------------
dp <- dotplot(
  ego,
  showCategory = 10,
  font.size   = 40,
  title       = "KO CLP vs WT Sham down"
)

dp <- dp +
  scale_size(range = c(6, 20)) +
  theme(
    axis.text.y  = element_text(size = 50),
    legend.text  = element_text(size = 20),
    legend.title = element_text(size = 22),
    plot.title   = element_text(size = 30, face = "bold")
  )

ggsave(
  "enrichGO_dotplotDOWNKOCLPvsWTSham.pdf",
  plot = dp, width = 13, height = 9
)

dp <- dotplot(
  ego,
  showCategory = 10,
  title       = "EnrichGO // KO CLP vs WT Sham - down",
  font.size   = 9
)
dp + theme(
  axis.text.y = element_text(size = 16),
  plot.title  = element_text(size = 14, face = "bold")
)
# (optional write suppressed)



# ============================================================================================
# FORTH CONTRAST: comp_file <- deseq_txt[5]
# ============================================================================================

comp_file <- deseq_txt[5]
comp_file

title <- gsub("_080525.txt", "", comp_file)
title

dat <- read.table(comp_file, header = TRUE)
dim(dat)
names(dat)

lim.P  <- 0.05
# comment in source says “change to 0”, but preserved code keeps:
lim.FC <- log2(2)

filtered_dat <- subset(
  dat,
  dat$padj < lim.P & abs(dat$log2FoldChange) > lim.FC
)
dim(filtered_dat)

filtered_datUP <- subset(
  dat,
  dat$log2FoldChange > lim.FC & dat$padj < lim.P
)
filtered_datUP <- droplevels(filtered_datUP)

filtered_datDN <- subset(
  dat,
  dat$log2FoldChange < -lim.FC & dat$padj < lim.P
)
filtered_datDN <- droplevels(filtered_datDN)

up <- which(filtered_datUP$log2FoldChange > 0)
dn <- which(filtered_datDN$log2FoldChange < 0)

# ALL ---------------------------------------------------------------------
geneList1 <- as.character(dat$entrezgene_id)

genes_for_pathways_all <- as.character(filtered_dat$entrezgene_id)

ego <- enrichGO(
  gene          = genes_for_pathways_all,
  universe      = geneList1,
  OrgDb         = "org.Mm.eg.db",
  ont           = "BP",
  pvalueCutoff  = 0.05,
  readable      = TRUE
)

library(ggplot2)

dp <- dotplot(
  ego,
  showCategory = 10,
  title       = "EnrichGO // KO vs WT Sham",
  font.size   = 9
)

dp + theme(
  axis.text.y = element_text(size = 16),
  plot.title  = element_text(size = 14, face = "bold")
)

res1 <- as.data.frame(ego)
head(res1)

write.table(
  res1,
  file = "pathwaysALL_KO_Sham_WT_Sham_160525.txt",
  sep = "\t", quote = FALSE, row.names = FALSE
)

# UP ----------------------------------------------------------------------
gene3 <- as.character(filtered_datUP$entrezgene_id)

ego <- enrichGO(
  gene          = gene3,
  universe      = geneList1,
  OrgDb         = "org.Mm.eg.db",
  ont           = "BP",
  pvalueCutoff  = 0.05,
  readable      = TRUE
)

### prettier ---------------------------------------------------------------
dp <- dotplot(
  ego,
  showCategory = 10,
  font.size   = 40,
  title       = "KO Sham vs WT Sham"
)

dp <- dp +
  scale_size(range = c(6, 20)) +
  theme(
    axis.text.y  = element_text(size = 50),
    legend.text  = element_text(size = 20),
    legend.title = element_text(size = 22),
    plot.title   = element_text(size = 30, face = "bold")
  )

ggsave(
  "enrichGO_dotplotUPKOShamWTSham.pdf",
  plot = dp, width = 15.5, height = 18
)

dp <- dotplot(
  ego,
  showCategory = 10,
  title       = "EnrichGO // KO vs WT Sham - up",
  font.size   = 9
)

dp + theme(
  axis.text.y = element_text(size = 16),
  plot.title  = element_text(size = 14, face = "bold")
)

dotplot(
  ego,
  showCategory = 30,
  title       = paste("EnrichGO //", title),
  font.size   = 9
)

res1 <- as.data.frame(ego)
head(res1)
comp_file
# write.table(
#   res1,
#   file = "pathways_up_TRPV1_KO_Sham_vs_WT_Sham.txt",
#   sep = "\t", quote = FALSE, row.names = FALSE
# )

gl1 <- filtered_datUP$log2FoldChange
names(gl1) <- filtered_datUP$entrezgene_id

geneList10 <- gl1

edox <- ego

p1 <- cnetplot(
  edox,
  node_label        = "category",
  showCategory      = 15,
  max.overlaps      = 100,
  cex_label_category= 0.8,
  foldChange        = geneList10
)

p2 <- cnetplot(
  edox,
  node_label    = "gene",
  cex_label_gene= 1,
  showCategory  = 15,
  max.overlaps  = 100
)

cowplot::plot_grid(p1, p2, ncol = 2, labels = LETTERS[1:2])

# DOWN --------------------------------------------------------------------
gene_symbols <- filtered_datDN$gene_symbol

gene3 <- bitr(
  gene_symbols,
  fromType = "SYMBOL",
  toType   = "ENTREZID",
  OrgDb    = "org.Mm.eg.db"
)$ENTREZID

gene3 <- as.character(filtered_datDN$entrezgene_id)

ego <- enrichGO(
  gene          = gene3,
  universe      = geneList1,
  OrgDb         = "org.Mm.eg.db",
  ont           = "BP",
  pvalueCutoff  = 0.05,
  readable      = TRUE
)

res1 <- as.data.frame(ego)
head(res1)

write.table(
  res1,
  file = "pathways_down_CLP_PCT_AB_vs_Sham_PCT_AB_280425.txt",
  sep = "\t", quote = FALSE, row.names = FALSE
)

dotplot(
  ego,
  showCategory = 30,
  title       = paste("EnrichGO //", title),
  font.size   = 9
)

gl1 <- filtered_datDN$log2FoldChange
names(gl1) <- filtered_datDN$entrezgene_id

geneList10 <- gl1

cnetplot(
  ego,
  categorySize       = "pvalue",
  foldChange         = geneList10,
  showCategory       = 15,
  max.overlaps       = 500,
  cex_category       = 0.8,
  cex_gene           = 0.5,
  cex_label_category = 0.8,
  cex_label_gene     = 0.7
)

edox <- ego

p1 <- cnetplot(
  edox,
  node_label        = "category",
  showCategory      = 30,
  max.overlaps      = 100,
  cex_label_category= 0.8,
  foldChange        = geneList10
)

p2 <- cnetplot(
  edox,
  node_label    = "gene",
  cex_label_gene= 1,
  showCategory  = 30,
  max.overlaps  = 100
)

cowplot::plot_grid(p1, p2, ncol = 2, labels = LETTERS[1:2])



# ============================================================================================
# FIFTH CONTRAST: comp_file <- deseq_txt[6]
# ============================================================================================

comp_file <- deseq_txt[6]
comp_file

title <- gsub("_080525.txt", "", comp_file)
title

dat <- read.table(comp_file, header = TRUE)
dim(dat)
names(dat)

lim.P  <- 0.05
lim.FC <- log2(2)

filtered_dat <- subset(
  dat,
  dat$padj < lim.P & abs(dat$log2FoldChange) > lim.FC
)
dim(filtered_dat)

filtered_datUP <- subset(
  dat,
  dat$log2FoldChange > lim.FC & dat$padj < lim.P
)
filtered_datUP <- droplevels(filtered_datUP)

filtered_datDN <- subset(
  dat,
  dat$log2FoldChange < -lim.FC & dat$padj < lim.P
)
filtered_datDN <- droplevels(filtered_datDN)

up <- which(filtered_datUP$log2FoldChange > 0)
dn <- which(filtered_datDN$log2FoldChange < 0)

# ALL ---------------------------------------------------------------------
geneList1 <- as.character(dat$entrezgene_id)

genes_for_pathways_all <- as.character(filtered_dat$entrezgene_id)

ego <- enrichGO(
  gene          = genes_for_pathways_all,
  universe      = geneList1,
  OrgDb         = "org.Mm.eg.db",
  ont           = "BP",
  pvalueCutoff  = 0.05,
  readable      = TRUE
)

library(ggplot2)

dp <- dotplot(
  ego,
  showCategory = 10,
  title       = "EnrichGO // WT CLP vs Sham",
  font.size   = 9
)

dp + theme(
  axis.text.y = element_text(size = 16),
  plot.title  = element_text(size = 14, face = "bold")
)

res1 <- as.data.frame(ego)
head(res1)

write.table(
  res1,
  file = "pathwaysALL_WT_CLP_vs_WT_Sham_160525.txt",
  sep = "\t", quote = FALSE, row.names = FALSE
)

# UP ----------------------------------------------------------------------
gene3 <- as.character(filtered_datUP$entrezgene_id)

ego <- enrichGO(
  gene          = gene3,
  universe      = geneList1,
  OrgDb         = "org.Mm.eg.db",
  ont           = "BP",
  pvalueCutoff  = 0.05,
  readable      = TRUE
)

### prettier ---------------------------------------------------------------
dp <- dotplot(
  ego,
  showCategory = 10,
  font.size   = 40,
  title       = "WT CLP vs WT Sham up"
)

dp <- dp +
  scale_size(range = c(6, 20)) +
  theme(
    axis.text.y  = element_text(size = 50),
    legend.text  = element_text(size = 20),
    legend.title = element_text(size = 22),
    plot.title   = element_text(size = 30, face = "bold")
  )

ggsave(
  "enrichGO_dotplotUPCLPvsShamTRPV1.pdf",
  plot = dp, width = 15.5, height = 20
)

dp <- dotplot(
  ego,
  showCategory = 10,
  title       = "EnrichGO // WT CLP vs Sham - up",
  font.size   = 9
)

dp + theme(
  axis.text.y = element_text(size = 16),
  plot.title  = element_text(size = 14, face = "bold")
)

dotplot(
  ego,
  showCategory = 30,
  title       = paste("EnrichGO //", title),
  font.size   = 9
)

res1 <- as.data.frame(ego)
head(res1)

write.table(
  res1,
  file = "pathways_up_WT_CLP_vs_WT_Sham.txt",
  sep = "\t", quote = FALSE, row.names = FALSE
)

gl1 <- filtered_datUP$log2FoldChange
names(gl1) <- filtered_datUP$entrezgene_id

geneList10 <- gl1

edox <- ego

p1 <- cnetplot(
  edox,
  node_label        = "category",
  showCategory      = 20,
  max.overlaps      = 100,
  cex_label_category= 0.8,
  foldChange        = geneList10
)

p2 <- cnetplot(
  edox,
  node_label    = "gene",
  cex_label_gene= 1,
  showCategory  = 20,
  max.overlaps  = 100
)

cowplot::plot_grid(p1, p2, ncol = 2, labels = LETTERS[1:2])

# DOWN --------------------------------------------------------------------
gene_symbols <- filtered_datDN$gene_symbol

gene3 <- bitr(
  gene_symbols,
  fromType = "SYMBOL",
  toType   = "ENTREZID",
  OrgDb    = "org.Mm.eg.db"
)$ENTREZID

gene3 <- as.character(filtered_datDN$entrezgene_id)

ego <- enrichGO(
  gene          = gene3,
  universe      = geneList1,
  OrgDb         = "org.Mm.eg.db",
  ont           = "BP",
  pvalueCutoff  = 0.05,
  readable      = TRUE
)

res1 <- as.data.frame(ego)
head(res1)

write.table(
  res1,
  file = "pathways_down_CLP_PCT_AB_vs_Sham_PCT_AB_280425.txt",
  sep = "\t", quote = FALSE, row.names = FALSE
)

dotplot(
  ego,
  showCategory = 30,
  title       = paste("EnrichGO //", title),
  font.size   = 9
)

gl1 <- filtered_datDN$log2FoldChange
names(gl1) <- filtered_datDN$entrezgene_id

geneList10 <- gl1

cnetplot(
  ego,
  categorySize       = "pvalue",
  foldChange         = geneList10,
  showCategory       = 15,
  max.overlaps       = 500,
  cex_category       = 0.8,
  cex_gene           = 0.5,
  cex_label_category = 0.8,
  cex_label_gene     = 0.7
)

edox <- ego

p1 <- cnetplot(
  edox,
  node_label        = "category",
  showCategory      = 15,
  max.overlaps      = 100,
  cex_label_category= 0.8,
  foldChange        = geneList10
)

p2 <- cnetplot(
  edox,
  node_label    = "gene",
  cex_label_gene= 1,
  showCategory  = 15,
  max.overlaps  = 100
)

cowplot::plot_grid(p1, p2, ncol = 2, labels = LETTERS[1:2])


################################################################################################
# END — All five contrasts processed with identical logic as Kls. 
# Headers normalized; object names and code flow preserved.
################################################################################################
