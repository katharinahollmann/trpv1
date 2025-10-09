################################################################################################
# ENRICHR SCRIPT — EXPLAINED — xx KH
# TRPV1 RNA-seq | Enrichr-based functional enrichment
#
# Notes on key objects (preserved):
# - ll1 : directory listing
# - sm  : chosen DEG file (character path)
# - nam1: title derived from sm
# - dd01: full DEG table from sm
# - dd1 : filtered DEGs (|log2FC| > 1 & padj < 0.05)
# - up / down : counts of up- and down-regulated genes
# - genALL / genUP / genDN : gene symbol vectors for enrichment
# - dbs : Enrichr libraries to query
# - enriched / res1 : Enrichr results containers
################################################################################################


######################################################################
# 0) SETUP: working directory & packages // adjust to your system !!
######################################################################

setwd("/Users/katharinahollmann/GitHubDesktop/Wagner_comparative_analyses/TRPV1/data")
dir()

# Packages ----------------------------------------------------------------
# BiocManager::install("enrichR")
library(enrichR)
library(dplyr)


##############################################
# 1) INPUT DEGs: choose file & basic filtering
##############################################

# Get files & choose DEG table --------------------------------------------
ll1 <- dir(); ll1

# choose the file -- deg files are 12-16 (original note)
sm <- ll1[18]; sm                       # << adjust index as needed
nam1 <- gsub("_020625.txt", "", sm); nam1

# Read DEG table -----------------------------------------------------------
dd01 <- read.table(sm, header = TRUE)
names(dd01)
head(dd01)

# Filter DEGs (|log2FC| > 1 & padj < 0.05) --------------------------------
dd1 <- subset(
  dd01,
  abs(dd01$log2FoldChange) > log2(2) & dd01$padj < 0.05
)

dim(dd01)    # total rows
dim(dd1)     # passing thresholds

# Counts: Up vs Down -------------------------------------------------------
up   <- length(which(dd1$log2FoldChange > 0))
down <- length(which(dd1$log2FoldChange < 0))

dim(dd1)
up
down


##############################################
# 2) QUICK BARPLOT: Up vs Down (base R)
##############################################

counts_df <- dd1 |>
  mutate(Direction = ifelse(log2FoldChange > 0, "Up", "Down")) |>
  count(Direction) |>
  arrange(Direction)

# Stacked-bar matrix
bar_matrix <- matrix(counts_df$n, ncol = 1)
rownames(bar_matrix) <- counts_df$Direction

# Legend labels with counts
legend_labels <- paste0(counts_df$Direction, " (", counts_df$n, ")")

# Plot (default margins) ---------------------------------------------------
barplot(
  bar_matrix,
  col = c("blue", "red"),
  ylab = "Number of Genes",
  main = paste("Up vs Down -", nam1),
  names.arg = "DEGs",
  legend.text = legend_labels,
  args.legend = list(x = "topright", fill = c("blue", "red"))
)

# Optional: reset margins & alternative layout -----------------------------
par(mar = c(5, 5, 4, 8))  # bottom, left, top, right

barplot(
  bar_matrix,
  col = c("blue", "red"),
  ylab = "Number of Genes",
  main = paste("Up vs Down -", nam1),
  names.arg = "DEGs",
  legend.text = legend_labels,
  args.legend = list(x = "topright", inset = c(-0.2, 0), fill = c("blue", "red"), bty = "n"),
  width = 0.1,    # narrower bar
  cex.names = 0.8,
  cex.main  = 1
)


##############################################
# 3) BUILD GENE LISTS (ALL / UP / DOWN)
##############################################

# ALL ---------------------------------------------------------------------
names(dd1)
dim(dd1)

genALL <- as.character(dd1$gene_symbol)
genALL
head(genALL)
length(genALL)

# UP ----------------------------------------------------------------------
UP1 <- subset(dd01, dd01$log2FoldChange >  log2(2) & dd01$padj < 0.05)
names(UP1)
dim(UP1)

genUP <- as.character(UP1$gene_symbol)
genUP
head(genUP)
length(genUP)

# DOWN --------------------------------------------------------------------
DN1 <- subset(dd01, dd01$log2FoldChange < -log2(2) & dd01$padj < 0.05)
names(DN1)
dim(DN1)

genDN <- as.character(DN1$gene_symbol)
genDN
head(genDN)
length(genDN)


##############################################################
# 4) ENRICHR: database selection & quick plots 
##############################################################

listEnrichrSites()
setEnrichrSite("Enrichr")   # Human genes endpoint
websiteLive <- TRUE

# Available libraries ------------------------------------------------------
dbs <- listEnrichrDbs()
if (is.null(dbs)) websiteLive <- FALSE
if (websiteLive) head(dbs)
names(dbs)
dbs[3]    # original peek

# Choose working set of libraries (as in your script) ---------------------
dbs <- c(
  "Transcription_Factor_PPIs",
  "Reactome_2016",
  "GO_Biological_Process_2021",
  "KEGG_2021_Human",
  "WikiPathway_2021_Human",
  "MSigDB_Hallmark_2020"
)

# ################### for short fast analysis
# dbs <- c("MSigDB_Hallmark_2020","GO_Biological_Process_2021")
# dbs <- c("MSigDB_Hallmark_2020")


# === UP genes: Enrichr call + base plot ==================================
if (websiteLive) {
  enriched <- enrichr(genUP, dbs)  # input gene list
}

names(enriched)
sel1 <- 3
names(enriched)[sel1]     # "GO_Biological_Process_2021" (in this setup)

res1 <- if (websiteLive) enriched[[sel1]]
res1
names(res1)
head(res1)
dim(res1)

nam2 <- paste(nam1, "// DB:", names(enriched)[sel1]); nam2

plotEnrich(
  enriched[[sel1]],
  showTerms = 10,
  numChar   = 50,
  y         = "Ratio",
  orderBy   = "P.value",
  title     = nam2
)   # y = "Ratio" or "Count"


##############################################
# 5) PRETTIER — Tidy helpers & publication plots
##############################################

# ---- packages (as in your prettier block) --------------------------------
library(enrichR)
library(dplyr)
library(stringr)
library(ggplot2)

# Dependencies (optional install) -----------------------------------------
pkgs <- c("enrichR", "dplyr", "stringr", "ggplot2")
to_install <- pkgs[!pkgs %in% rownames(installed.packages())]
if (length(to_install)) install.packages(to_install)
lapply(pkgs, library, character.only = TRUE)

# Config -------------------------------------------------------------------
websiteLive <- TRUE   # required by enrichR

# Pick Enrichr DB robustly -------------------------------------------------
pick_enrichr_db <- function(preferred = "KEGG_2021_Mouse") {
  all <- enrichR::listEnrichrDbs()$libraryName
  if (preferred %in% all) return(preferred)
  kegg_like <- grep("^KEGG", all, value = TRUE)
  if (length(kegg_like) == 0) stop("No KEGG library found in Enrichr.")
  message("Using KEGG db: ", kegg_like[1])
  kegg_like[1]
}

# Tidy Enrichr result ------------------------------------------------------
tidy_enrichr <- function(genes, db_name, websiteLive = TRUE) {
  stopifnot(length(genes) > 0)
  if (!websiteLive) stop("Set websiteLive = TRUE for Enrichr queries.")
  enr <- enrichR::enrichr(genes, db_name)
  if (!db_name %in% names(enr)) {
    stop("Returned Enrichr object does not contain db: ", db_name)
  }
  res <- enr[[db_name]]
  if (is.null(res) || nrow(res) == 0) {
    warning("No enrichment terms returned.")
    return(res)
  }
  res |>
    mutate(
      Term        = stringr::str_replace_all(Term, "_", " "),
      Count       = as.numeric(stringr::str_extract(Overlap, "^[0-9]+")),
      SetSize     = as.numeric(stringr::str_extract(Overlap, "(?<=/)[0-9]+")),
      GeneRatio   = Count / SetSize,
      negLog10FDR = -log10(Adjusted.P.value)
    )
}

# Plotting helper (dot or bar) --------------------------------------------
plot_enrichr <- function(res,
                         top_n       = 10,
                         title       = "KEGG enrichment",
                         type        = c("dot", "bar"),
                         font_axis_y = 50,
                         font_title  = 30,
                         font_legend = 22,
                         dot_range   = c(6, 20)) {
  type <- match.arg(type)
  if (is.null(res) || nrow(res) == 0) stop("Nothing to plot: result is empty.")
  
  df <- res |>
    arrange(desc(Count)) |>
    slice_head(n = top_n) |>
    mutate(Term = factor(Term, levels = rev(unique(Term))))
  
  if (type == "dot") {
    p <- ggplot(df, aes(x = GeneRatio, y = Term, size = Count, color = negLog10FDR)) +
      geom_point() +
      scale_size(range = dot_range) +
      labs(
        title = title, x = "Gene Ratio", y = NULL,
        size  = "Count", color = "-log10(FDR)"
      ) +
      theme_minimal(base_size = 16) +
      theme(
        axis.text.y  = element_text(size = font_axis_y),
        legend.text  = element_text(size = font_legend - 2),
        legend.title = element_text(size = font_legend),
        plot.title   = element_text(size = font_title, face = "bold")
      )
  } else {
    p <- ggplot(df, aes(x = Term, y = negLog10FDR)) +
      geom_col(width = 0.5, fill = "steelblue") +
      coord_flip() +
      labs(title = title, x = NULL, y = "-log10(FDR)") +
      theme_minimal(base_size = 30) +
      theme(
        axis.text.y = element_text(size = font_axis_y),
        plot.title  = element_text(size = font_title, face = "bold")
      )
  }
  p
}

# Save PDF helper ----------------------------------------------------------
save_pdf <- function(plot, file, width = 12, height = 9, use_cairo = FALSE) {
  dev <- if (use_cairo) cairo_pdf else "pdf"
  ggplot2::ggsave(filename = file, plot = plot, width = width, height = height, device = dev)
}

# Choose DB for pretty plots ----------------------------------------------
db <- pick_enrichr_db("GO_Biological_Process_2021")   # change to "KEGG_2021_Human" for human as needed
db

# UP — pretty dot plot ----------------------------------------------------
res_up <- tidy_enrichr(genUP, db)

p_up_dot <- plot_enrichr(
  res_up,
  top_n      = 10,
  title      = "WT CLP vs WT Sham (UP) — KEGG pathways",
  type       = "dot",
  font_axis_y = 50,
  font_title  = 30,
  font_legend = 22,
  dot_range   = c(6, 20)
)
print(p_up_dot)

save_pdf(p_up_dot, "KEGG_CLPvsShamTRPV1_UP_dot.pdf", width = 28, height = 20)

# DOWN — pretty dot plot --------------------------------------------------
res_dn <- tidy_enrichr(genDN, db)

p_dn_dot <- plot_enrichr(
  res_dn,
  top_n      = 10,
  title      = "WT CLP vs WT Sham (DOWN) — KEGG pathways",
  type       = "dot",
  font_axis_y = 50,
  font_title  = 30,
  font_legend = 22,
  dot_range   = c(6, 20)
)
print(p_dn_dot)

save_pdf(p_dn_dot, "KEGG_CLPvsShamTRPV1_DN_dot.pdf", width = 28, height = 20)

# (just to keep original variable visible)
sm


##############################################
# 6) OPTIONAL: DOWN bar-plot version
##############################################

res_dn <- tidy_enrichr(genDN, db)

p_dn_bar <- plot_enrichr(
  res_dn,
  top_n       = 10,
  title       = "CLP vs WT Sham (DOWN) — KEGG pathways",
  type        = "bar",
  font_axis_y = 60,
  font_title  = 50,
  font_legend = 50
)

print(p_dn_bar)

save_pdf(p_dn_bar, "KOCLP_vs_WTSham_DOWN_KEGG_bar.pdf", width = 48, height = 30)


##############################################
# 7) ORIGINAL ENRICHR CALLS (DOWN / ALL) — base plotEnrich
##############################################

# for DOWN ---------------------------------------------------------------
if (websiteLive) {
  enriched <- enrichr(genDN, dbs)
}

names(enriched)
sel1 <- 3
names(enriched)[sel1]   # "GO_Biological_Process_2021"

res1 <- if (websiteLive) enriched[[sel1]]
res1
names(res1)
head(res1)
dim(res1)

nam2 <- paste(nam1, "// DB:", names(enriched)[sel1]); nam2

plotEnrich(
  enriched[[sel1]],
  showTerms = 10,
  numChar   = 50,
  y         = "Ratio",
  orderBy   = "P.value",
  title     = nam2
)


# for ALL ----------------------------------------------------------------
if (websiteLive) {
  enriched <- enrichr(genALL, dbs)
}

names(enriched)
sel1 <- 3
names(enriched)[sel1]   # "GO_Biological_Process_2021"

res1 <- if (websiteLive) enriched[[sel1]]
res1
names(res1)
head(res1)
dim(res1)

nam2 <- paste(nam1, "// DB:", names(enriched)[sel1]); nam2

plotEnrich(
  enriched[[sel1]],
  showTerms = 10,
  numChar   = 50,
  y         = "Ratio",
  orderBy   = "P.value",
  title     = nam2
)

################################################################################################
# END — All original logic preserved; sections & prettier helpers integrated
################################################################################################