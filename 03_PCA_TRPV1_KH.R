# ============================================================
# PCA ANALYSIS SCRIPT — EXPLAINED -- xx KH
# TRPV1 RNA-seq dataset: PCA visualization of normalized counts
# Notes:
#  • Object names, filenames, and logic unchanged.
#  • Code reorganized and commented for clarity to best of my knowledge - no guarantee
# ============================================================


# ============================================================
# 1) SETUP AND DATA IMPORT
# ============================================================

# --- Set working directory ---
# Adjust this path to your system. All plots and outputs will be written here.
# To find your absolute path in Terminal: navigate to your folder, run `pwd`, paste below.
# Generated outputs are written to this folder unless a different path is set in write commands.

setwd("/Users/katharinahollmann/GitHubDesktop/Wagner_comparative_analyses/TRPV1/data")
dir()  # confirm files are visible

# --- Load normalized data (rlog output from DESeq2) ---
# Use batch-corrected file if available/need (SITAG Dataset only for now); otherwise use this normalization.
data <- read.table("norm_Wagner_CLP_TRPV1_2023_080525.txt", header = TRUE, sep = "")
names(data)
dim(data)
class(data)
data[1:5, 1:8]  # preview first annotation columns

# --- Check for duplicated gene symbols (sanity check) ---
table(duplicated(data$gene_symbol))  # duplicates should be 0


# ============================================================
# 2) PREPARE MATRIX FOR PCA
# ============================================================

# DESeq2 normalized table: first 8 columns are annotations, counts start at column 9
countdata <- as.matrix(data[, 9:ncol(data)])
colnames(countdata)
dim(countdata)
head(countdata)

# Quick QC plot to inspect data distribution
par(mfrow = c(1,1), font = 2, font.axis = 2, font.lab = 3, mar = c(12,4,2,4))
boxplot(countdata, 
        las = 3, 
        col = "red", 
        main = "Normalized counts per sample")


# ============================================================
# 3) LOAD SAMPLE METADATA (TARGET FILE)
# ============================================================

dir()
ta <- read.table("target_Wagner_CLP_TRPV1_2023_240124.txt", header = TRUE)
names(ta)
head(ta)
dim(ta)

# Confirm that column names (samples) match between count matrix and metadata
identical(colnames(countdata), as.character(ta$sample_ID))
which(colnames(countdata) != as.character(ta$sample_ID))  # should return integer(0)
length(colnames(countdata))
length(as.character(ta$sample_ID))

# --- Quick summary for PowerPoint or methods section ---
dim(ta)
table(ta$group)

# Ensure consistent group order across datasets
ta$group <- factor(ta$group,
                   levels = c("WT_CLP", "TRPV1_KO_CLP", "WT_Sham", "TRPV1_KO_Sham"))
table(ta$group)

group <- ta$group  # shorthand for plotting


# ============================================================
# 4) PCA COMPUTATION
# ============================================================

# Perform PCA
# t(countdata) transposes so that samples = rows and genes = columns (required)
# scale = FALSE because DESeq2-normalized data is already variance-stabilized
pca <- prcomp(t(countdata), scale = FALSE)
names(pca)
pca$x[1:3, 1:3]  # preview first principal components
summary(pca)

# Scree plot: visual check of explained variance
par(font = 2, font.axis = 2, font.lab = 3, mar = c(4,4,2,2))
plot(prcomp(t(countdata)), main = "Scree plot of PCA")

# Percent variance explained by first PCs
summary(pca)$importance[, 1:5]
PC1 <- round(summary(pca)$importance[2, 1] * 100)
PC2 <- round(summary(pca)$importance[2, 2] * 100)
PC3 <- round(summary(pca)$importance[2, 3] * 100)
PC4 <- round(summary(pca)$importance[2, 4] * 100)

# Combine PCA scores with group metadata
scores <- data.frame(group, pca$x[, 1:5])
head(scores)


# ============================================================
# 5) PCA VISUALIZATION (GGPLOT WITH ELLIPSES)
# ============================================================

library(tidyverse)
library(viridis)
library(ggplot2)

sample_ID <- ta$sample_ID

# Combine PCA scores and metadata
pca_df <- data.frame(
  PC1 = scores$PC1,
  PC2 = scores$PC2,
  group = as.factor(group),
  sample_ID = sample_ID
)

# PCA plot with ellipses (“clouds”) per group
# viridis for colorblind-safe option
pca_plot <- ggplot(pca_df, aes(PC1, PC2, color = group, fill = group)) +
  geom_point(size = 3, shape = 21, stroke = 1.2) +
  stat_ellipse(type = "norm", alpha = 0.1, geom = "polygon") +
  labs(title = "PCA plot by group",
       x = paste("PC1:", PC1, "%"),
       y = paste("PC2:", PC2, "%")) +
  scale_color_viridis(
    discrete = TRUE, option = "H", name = "Groups",
    breaks = c("WT_CLP","TRPV1_KO_CLP","WT_Sham","TRPV1_KO_Sham"),
    labels = c("CLP","KO CLP","Sham","KO Sham")
  ) +
  scale_fill_viridis(
    discrete = TRUE, option = "H", name = "Groups",
    breaks = c("WT_CLP","TRPV1_KO_CLP","WT_Sham","TRPV1_KO_Sham"),
    labels = c("CLP","KO CLP","Sham","KO Sham")
  ) +
  theme_minimal(base_size = 10) +
  theme(legend.position = c(0.85, 1),
        legend.justification = c("top"))

pca_plot
ggsave("PCA_plotTRPV1.pdf", plot = pca_plot, width = 6, height = 6)


# ============================================================
# 6) ALTERNATIVE VISUALIZATION (MANUAL COLORS) -- NO CLOUDS 
# ============================================================

# Define colors manually to match desired palette
group_levels <- levels(pca_df$group)
col2 <- c("darkred","magenta","darkblue","darkgreen")[1:length(group_levels)]
names(col2) <- group_levels
shape2 <- rep(19, length(group_levels))
names(shape2) <- group_levels

p <- ggplot(pca_df, aes(x = PC1, y = PC2, color = group, shape = group)) +
  geom_point(size = 3) +
  scale_color_manual(values = col2) +
  scale_shape_manual(values = shape2) +
  labs(
    title = "PCA plot by group",
    x = paste0("PC1: ", PC1, "%"),
    y = paste0("PC2: ", PC2, "%"),
    color = "Group",
    shape = "Group"
  ) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        legend.position = "right",
        legend.box.margin = margin(0, 10, 0, 0))

print(p)
ggsave("PCA_plot_ggplot2.pdf", plot = p, width = 8, height = 6)


# ============================================================
# 7) NOTES
# ============================================================
# • countdata was already normalized (rlog-transformed) → no scaling.
# • Group order consistent across scripts for reproducibility.
# • Viridis and manual color versions provide colorblind-safe options.
# • Ellipse type = "norm" shows multivariate normal distribution boundaries.
# • Output: PCA_plotTRPV1.pdf and PCA_plot_ggplot2.pdf (vector, high-quality)
# ============================================================