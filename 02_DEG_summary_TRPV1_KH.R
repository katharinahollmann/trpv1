# ============================================================
# TRPV1 DEG SUMMARY SCRIPT — EXPLAINED -- KH
# Purpose: Summarize DESeq2 results for all TRPV1 contrasts,
#           count up-/downregulated genes, and visualize.
# ============================================================

### 1. SETUP AND FILE SELECTION ---------------------------------------------

# Set working directory to your analysis folder
# Adjust the path to your system if needed
# To find your absolute path in Terminal: navigate to your folder, run `pwd`, paste below.
# Generated outputs are written to this folder unless a different path is set in write commands.
setwd("/Users/katharinahollmann/GitHubDesktop/Wagner_comparative_analyses/TRPV1/data")

# List all files to confirm available datasets
dir()

# Identify and select DESeq2 result files
allfiles <- dir()
deseq_files <- grep("DESeq2", allfiles, value = TRUE)   # list files containing "DESeq2"
deseq_files
deseq_files <- deseq_files[c(5:9)]                      # select only DEG result files
deseq_files

#what is the log2fc of trpv1 here?
wt_clp <- read.table(deseq_files[5], header = TRUE)

dim(wt_clp)
names(wt_clp)

trpv1_wt_clp <- wt_clp[wt_clp$gene_symbol == "Trpv1", ]
trpv1_wt_clp #2.629654


### 2. CREATE EMPTY MATRIX TO STORE RESULTS -------------------------------

# Prepare an empty 2-column matrix (UP, DOWN) for each file
mat1 <- matrix(ncol = 2, nrow = length(deseq_files))
mat2 <- as.data.frame(mat1) ## this is empty
names(mat2) <- c("UP", "DOWN")
row.names(mat2) <- as.character(deseq_files)  # use filenames as row names
mat2  # check setup

# Preview one DESeq2 file to understand its structure
dd1 <- read.table(deseq_files[1], header = TRUE)
names(dd1)
head(dd1)

### 3. DEFINE FILTER THRESHOLDS -------------------------------------------

# Set significance thresholds for differential expression
lim.P <- 0.05               # adjusted p-value cutoff
lim.FC <- log2(2)           # log2 fold change cutoff = ±1
lim.FC

### 4. LOOP THROUGH FILES AND COUNT DEGS ----------------------------------

# Loop over each DESeq2 file and extract DEG counts
for (i in 1:length(deseq_files)) {
  file1 <- as.character(deseq_files[i])
  y2 <- read.table(file1, header = TRUE)  # read file
  head(y2)
  
  # Filter for genes meeting the significance thresholds
  y3 <- subset(y2, y2$padj < lim.P & abs(y2$log2FoldChange) > lim.FC)
  y4 <- y3[order(y3$log2FoldChange, decreasing = TRUE), ]  # order by fold change
  
  # Identify up- and downregulated genes
  up.reg   <- which(y4$log2FoldChange > 0)
  up1      <- y4[up.reg, ]
  down.reg <- which(y4$log2FoldChange < 0)
  down1    <- y4[down.reg, ]
  
  # Store results in the summary matrix
  mat2[i, 1:2] <- c(as.numeric(dim(up1)[1]), as.numeric(dim(down1)[1]))
}

# Inspect DEG counts per comparison
mat2 

# Add a total column (UP + DOWN)
sum1 <- apply(mat2, 1, sum)
mat2$Total <- sum1
mat2 #this is now filled with data & total --------- CHECK 

### 5. REORDER CONTRASTS AND EXPORT TABLE ---------------------------------

row.names(mat2)
#ordered <- c(5,3,2,4,1)

# Define preferred contrast order for readability
ordered <- c(
  "DEG_DESeq2_TRPV1_WT_CLP_vs_WT_Sham_080525.txt",
  "DEG_DESeq2_TRPV1_TRPV1_KO_CLP_vs_WT_Sham_080525.txt",
  "DEG_DESeq2_TRPV1_TRPV1_KO_CLP_vs_TRPV1_KO_Sham_080525.txt",
  "DEG_DESeq2_TRPV1_TRPV1_KO_CLP_vs_WT_CLP_080525.txt",
  "DEG_DESeq2_TRPV1_TRPV1_KO_Sham_vs_WT_Sham_080525.txt"
)

mat2 <- mat2[ordered, ]
mat2

# Save summarized DEG counts
write.table(mat2, file = "DEG_DESeq2_summary_table_080925.txt")

### 6. PLOT DEG SUMMARY — GGPLOT VERSION ----------------------------------

# Reload summary table
mat2 <- read.table("DEG_DESeq2_summary_table_080925.txt", header = TRUE)

# Rename rows for readability in the plot
new_names <- c(
  "WT CLP vs WT Sham",
  "KO CLP vs WT Sham",
  "KO CLP vs KO Sham",
  "KO CLP vs WT CLP",
  "KO Sham vs WT Sham"
)
row.names(mat2) <- new_names
lab1 <- as.character(row.names(mat2))
mat2

# (Optional) overwrite matrix manually if needed (e.g., summarized values)
######## CHECK matrix and put in numbers accordingly 
mat2 <- data.frame(
  UP   = c(1450, 1175, 583, 327, 411),
  DOWN = c(568, 509, 396, 467, 10),
  row.names = new_names
)
mat2$Comparison <- rownames(mat2)

# Convert to long format for ggplot
library(ggplot2)
library(tidyr)
library(forcats)
library(tidyverse)

df_long <- pivot_longer(mat2, cols = c("UP", "DOWN"),
                        names_to = "Category", values_to = "Count")

# Preserve order of input comparisons
df_long$Comparison <- fct_inorder(df_long$Comparison)
df_long$Category <- factor(df_long$Category, levels = c("UP", "DOWN"))

# Bar plot of up/downregulated genes per comparison
b <- ggplot(df_long, aes(x = Comparison, y = Count, fill = Category)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.8) +
  labs(
    x = "Comparison of Groups",
    y = "Number of DEGs",
    title = "DEG summary TRPV1",
    subtitle = "up- and downregulated DEGs per group comparison"
  ) +
  scale_fill_manual(values = c("UP" = "red", "DOWN" = "blue")) +
  labs(fill = "Regulation") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = c(0.85, 0.85),
    legend.justification = c("right", "top")
  )

b  # display plot

# Export high-quality figures for publication or reports
ggsave("2610barplotDEGsummaryKOneworder.pdf", plot = b,
       device = "pdf", width = 4.3, height = 6, units = "in")

ggsave("barplotDEGsummaryKOneworder2.png", plot = b,
       device = "png", width = 3.0, height = 7, units = "in", dpi = 300)


### 7. ALTERNATIVE BASIC R BARPLOT FROM KLS -------------------------------

# Convert matrix for base R plotting
mat4 <- as.matrix(mat2[, c("UP", "DOWN")])
mat5 <- t(mat4)  # transpose for stacked display

# Base R plot — older version (from Klaus’ script)
par(mfrow = c(1, 1), font = 2, font.axis = 2, font.lab = 3, mar = c(16, 4, 3, 2), xpd = FALSE)
col1 <- c("red", "blue")

barplot(mat5, col = col1, las = 3, names.arg = lab1,
        main = "DEGs per contrast")
box(bty = "l")

# Add legend
groups.2 <- as.factor(row.names((mat5)))
shape2 <- rep(15, length(unique(groups.2)))
legend("topright", legend = groups.2, pch = shape2, cex = 1, pt.cex = 2, bty = "n", ncol = 1, col = col1)

# ============================================================
# END OF SCRIPT
# ============================================================
