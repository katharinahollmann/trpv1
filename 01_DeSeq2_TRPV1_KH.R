# =========================================================================
# THIS IS KATS DESEQ SKRIPT —— EXPLAINED -- KH 
# TRPV1 RNA-seq workflow: normalization, DEG calling, QC, export
# Notes:
#  • Object names, file names, and commands preserved exactly.
#  • Comments added for clarity, and to best of my knowledge, no guarantee
# =========================================================================

# ---------------------------
#  ENVIRONMENT INFORMATION
# ---------------------------
# Useful to record R session, packages, and platform for reproducibility
sessionInfo()

# check if you need to install packages

# if (!require("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
# BiocManager::install("DESeq2")
#
# ========================================================================
#  DATA LOADING AND BASICS
# ============================================================

# --- Working directory ---
# Absolute path used here. Adjust to your machine if needed.
# To find your absolute path in Terminal: navigate to your folder, run `pwd`, paste below.
# Generated outputs are written to this folder unless a different path is set in write commands.
setwd("/Users/katharinahollmann/GitHubDesktop/Wagner_comparative_analyses/TRPV1/data")
getwd()

# --- Libraries ---
# If not installed, use the commented BiocManager lines in your original header.
library(DESeq2)
library(genefilter)
library(RColorBrewer)
library(gplots)

# load your file data and check the file 
dir() # look at directory
raw_counts <- read.table("annot_raw_counts_Wagner_CLP_TRPV1_2023_070525.txt", header = TRUE, check.names = FALSE)
# header = TRUE means first lines contain the column names, check.names würde mit = TRUE die Namen auf "valid" variables anpassen 
names(raw_counts)
dim(raw_counts) ### note your number of columns to later choose the ones with log2FC
head(raw_counts)

## TEST ##

# control for genes and values -> is data correct, e.g. ENSMUSG00000005952 should be value = Tryp1
test1 <- grep("ENSMUSG00000005952", raw_counts$ENS_gene_ID, value = TRUE) #not necessary
test2 <- grep("ENSMUSG00000005952", raw_counts$ENS_gene_ID, value = FALSE) 
# value documentation: value
# if FALSE, a vector containing the (integer) indices of the matches determined by grep is returned, 
# and if TRUE, a vector containing the matching elements themselves is returned.
raw_counts[test2, ]
###################------------------------------------------------------------------

# prepare matrix for DESeq2 - it can only work with matrix as data type --------

# cut off the columns containing names et. -> only the columns with the countdata 
names(raw_counts) # look were (columnnr) count data starts
countdata <- raw_counts[,9:ncol(raw_counts)] ## from 9 (including) onwards (in 9 starts the first sample)
names(countdata)
dim(countdata)
head(countdata)

# visualize countdata before normalization
par(mfrow=c(1,1), mar=c(12,4,3,2)) #graphical settings for plot
boxplot(log2(countdata + 1), 
        col="blue", #colour of bars
        las=3, # x axis legend at 90 degree
        main="before normalization") # title
#boxplot for every column in countdata, transformed in log2 format with +1 (to avoid log(0)

# load design, match countdata matrix with metadata file (target), which includes groupnames, batch etc.
dir() #list of files again -> copy name of metadata file (target)

coldata <- read.table("target_Wagner_CLP_TRPV1_2023_240124.txt", header = TRUE) #coldata is metadata file 
head(coldata)
dim(coldata)

coldata$sample_ID
row.names(coldata) <- coldata$sample_ID # rows should be named like samples (sample_ID), which are on the first row at target file 

# check for you matrix prep, if you have the same samples in countdata(raw counts) and coldata (targetfile)
if (identical(colnames(countdata), as.character(coldata$sample_ID))) {
  #if they are identical, you want to continue with your analysis
  # as.character() function is used to coerce (convert) an object to a character vector
  # a character in R is a string in python
  # i need it here because sample_ID also contains numbers -> convert to character vector (to not be seen as numbers)
  cat("Column names match sample IDs. Proceeding ... \n")
} else {
  # if not identical, stop and print a warning
  stop("Column names of countdata do not match sample_ID in coldata.")
}

# check if nameing scheme is correct 
names(countdata)
names(coldata)
coldata$group


# running Deseq2 ----------------------------------------------------------

# RUN DESEQ ----- (you have to correct for batch in Sitag -> design for group and batch --> no need for batch in other datasets)
coldata$group <- as.factor(coldata$group)
#coldata$batch <- as.factor(coldata$batch)
#### A factor is how R stores categorical data — like labels or group, converts a vector into a factor

# create DeSeq object - copy design in pptx/method description, design ~ ... is important 
dds <- DESeqDataSetFromMatrix(countData=countdata, colData=coldata, design=  ~ group)
# DESeqDataSetFromMatrix(...) creates a DESeq2 dataset object from a matrix of counts and metadata
# countData = countdata	Your raw gene expression count matrix (rows = genes, columns = samples)
# colData = coldata	A metadata data frame (with info about each sample, like treatment group, condition, etc.)
# design = ~ group	A design formula telling DESeq2 which variable to use for modeling differences (in this case, a column called group in coldata)
head(dds)

dds$group <- as.factor(dds$group) # factor() is necessary, because group is a categorical variable (not continouus like age e.g.)
## !!!!!! By default, if your group column is a character or even numeric, R might not treat it as a proper categorical variable unless you explicitly convert it to factor  
table(dds$group) # how many in each group -> for pptx 
design(dds) <- ~ group

# set the levels so WT_Sham is the base/reference -- checked, no need
#dds$group <- factor(dds$group, levels=c("WT_Sham", "WT_CLP", "TRPV1_KO_Sham", "TRPV1_KO_CLP"))

levels(dds$group)

############ this is the process/calculation
# --> run DESeq2 ----- can take a while --------
dds <- DESeq(dds) ##### in console you can view what DESeq does --- estimating size factors ....--- fitting model and testing
# ----------------------------------------------
############ Deseq was run, now we check results 

resultsNames(dds) # to check what DESeq did 

# gather results without filter/threshold 
res <- results(dds)
res <- res[order(res$padj), ]
#  genes ranked by the smallest adjusted p-value, so you can focus on the most significant genes
# before the comma order generates a vector of row indices that sorts res by padj column and returns a list of numbers
# the sorted rows will be kept and put baack in new dataframe, the part after the comma is the column -> we keep all columns (empty)
# adjusted p-value in differential expression analysis to prioritize genes with the most evidence of being differentially expressed after accounting for multiple testing
head(res)
list(res) # the res data frame now becomes a list, which only contains one element, which is the res data frame
dim(res)
summary(res)

# done running and checking DESEQ -----------------------------------------

# Data prep and filtering (thresholds) ------------------------------------

# prepare for data filtering
table(is.na(res$log2FoldChange)) # frequency of NA, TRUE means NA (missing value)

# filter the data, filtering performed by DESeq
library("genefilter")
attr(results(dds), "filterThreshold") 
# I want to filter for the thresholds from the DESeq analyses result
# The attr() function in R is used to retrieve or set the attributes of an object. 
# In R, many objects (like vectors, matrices, and data frames) have additional attributes that store metadata, such as:
# Names of the elements, Dimensions of matrices or arrays, Class of the object (e.g., data frame, matrix)
# The attr() function allows you to access or modify these attributes.
# 
#  Syntax:
#   attr(x, which)
# x: The object from which you want to retrieve or modify the attribute.
# which: The name of the attribute you want to retrieve or modify (e.g., "names", "class", "dim").
# 🚨 Setting attributes:
#   You can also use attr() to set an attribute by passing a second argument:
#   
#   attr(x, which) <- value
plot(attr(res,"filterNumRej"),type="b",
     ylab="number of rejections")
# attr retrieves the numbers of rejection from res, the ordered results from DESeq2
#"filterNumRej" is an attribute within res that stores information about the number of rejections (possibly in some statistical test or filtering step).
# type="b":
#   The type argument controls the type of plot.
# "b" stands for both points and lines. So, this will plot the data with points at each value and connect those points with lines.
# Other options for type include:
#   "p" for points only
# "l" for lines only
# "o" for both points and lines but with different style
# ylab="number of rejections":
#   This sets the label for the y-axis of the plot.

resNoFilt <- results(dds, independentFiltering=FALSE)
# independentFiltering=FALSE:
#   By default, DESeq2 applies independent filtering when generating results -> RES HAS FILTERED APPLIED 
# Independent filtering is a method to remove genes with very low counts 
# or those that are less likely to show differential expression, to improve the overall power of the analysis.
# It adjusts the p-values of genes based on their count information.
# Setting independentFiltering=FALSE means you want to skip this step. 
# !! You are requesting the results without the filtering applied, meaning 
# !! all genes will be considered, even those with low counts or low significance.
#### --- !! we do this here and compare with filtering applied !! #########

table(filtering = (res$padj < 0.1), noFiltering = (resNoFilt$padj < 0.1))
# compair res (with independent filtering applied) and resNoFilt (without independent filtering applied)


# outlier correction ------------------------------------------------------

# replace NA values -> to not delete all genes from all samples if one has NA 
dds_outliercorrection <- replaceOutliersWithTrimmedMean(dds)
res_outliercorrection <- results(dds_outliercorrection)
res_outliercorrection <- res_outliercorrection[order(res_outliercorrection$padj), ] #order due to padj
head(res_outliercorrection)
dim(res_outliercorrection)

# compair data before and after outliercorrection 

tab <- table(initial = results(dds)$padj < 0.1,
             cleaned = results(dds_outliercorrection)$padj < 0.1) # how many before and after outlier correction were significant 
addmargins(tab) # to campair 
# addmargins() will show totals per row, column, and overall 
# addmargins() = add row and column totals (or other functions like mean, etc.) to a table or matrix.
# “Margins” means the edges — just like margins on a page — and in R, it refers to the sums (or other summaries) along rows/columns.
# -> adds sum columns of row totald and columns total 

dds <- dds_outliercorrection #after compaired, can be named dds again 


# outlier corrrection included in dds object -> done ----------------------


# visualization -----------------------------------------------------------
# now we prepare the plot 
# show values after outlier correction
range(res_outliercorrection$log2FoldChange,na.rm=T) #for log2FC
range(res_outliercorrection$padj,na.rm=T) #for padj
# anyNA(x)  # returns TRUE if any NA exists in x
# na.rm = FALSE (default) → the function returns NA if there's any missing value.
# na.rm = TRUE → removes all NAs before computing
# na.rm = TRUE in R, the NAs are temporarily ignored for that calculation — they’re not deleted or altered in your original data.

## PLOTS 
par(mfrow=c(1,1), mar=c(4,4,3,2)) #bottom, left, top, right 
plotMA(dds, 
       ylim=c(-5, 7.5), 
       main="// DESeq2 // red: adj.p <0.1")## dispersion plot (to see if extreme values pull the curve out)
# A warning will be printed that one should use plotDispEsts to check the quality of the fit
## whether the curve is pulled dramatically by a few outlier points.
plotDispEsts(dds)
## if this is the case perform local fit
# ddsLocal <- estimateDispersions(dds, fitType="local")

### look at plots, check data quality --- curve should look nicely with everything deseq does for you --> dispersion shrinkage etc.. ----


# contrast generation -----------------------------------------------------

#### generate tables with contrasts: 
####contrasts refers to the specific pairwise comparison between two groups (or conditions) that you want to test for differential expression #####
table(dds$group)


# contrast as function, which will automatically run the contrast for two groups
# in the following:
### you set the thresholds for DEGs -> log2(2) and padj < 0.05 here 
####### you could also define the thresholds as seperate objects e.g. LFC1 <- 1, padj1 <- 0.05 -> than padj < padj1 etc.
### you set the names for the generated contrast files -> adjust date etc.

run_contrast <- function(dds, raw_counts, cntr1, cntr2, log2FC_thresh=log2(2)) {  #set your log threshold; -- minimum log2FC threshold is default = 1.5
  nam <- paste("DEG_DESeq2_TRPV1", cntr1, "vs", cntr2, "080525.txt", sep="_") #name for generated data with date
  res_contrast <- results(dds, contrast=c("group", cntr1, cntr2)) # contrasts
  deg <- subset(res_contrast, abs(log2FoldChange) > log2FC_thresh & padj < 0.05) #filter, set thresholds
  cntr_annot <- cbind(raw_counts, res_contrast)
  write.table(cntr_annot, file=nam) ##this is the command for file with name nam abd data cntr_annot 
  return(dim(deg)) #gives back how many genes were significant
}

##### list of contrasts you want to run !!!!!!!!!!!!! ############
contrasts <- list(
  c("TRPV1_KO_CLP", "WT_CLP"),
  c("WT_CLP", "WT_Sham"),
  c("TRPV1_KO_CLP", "WT_Sham"),
  c("TRPV1_KO_Sham", "WT_Sham"),
  c("TRPV1_KO_CLP", "TRPV1_KO_Sham")
)

lapply(contrasts, function(x) run_contrast(dds, raw_counts, x[1], x[2]))
############## generate norm matrix -- takes MANY minutes !!! -----------

results(dds)

norms <- counts(dds, normalized = TRUE)
summary(norms)

head(norms)
row.names(norms) <- raw_counts$gene_symbol
head(norms)

trpv1row <- norms["Trpv1", ]

trpv1row

######### i want to visualize normalized counts dds #########
library(ggplot2)
library(tidyr)

# # ensure it's a data frame
# trpv1row <- as.data.frame(trpv1row)

# transpose to make samples a column
df <- data.frame(
  sample = names(trpv1row),
  count = as.numeric(trpv1row)
)

df# 



# bar plot
n <- ggplot(df, aes(x = sample, y = count)) +
  geom_col(fill = "steelblue", width = 0.6) +
  labs(title = "Norm. Trpv1 expression counts (before rlog)",
       x = "Sample",
       y = "Normalized count") +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(size = 14, angle = 90, vjust = 0.5, hjust = 1))

ggsave(filename = "trpv1countsNORM.pdf", plot = n, width = 20, height = 10)


#########
rld <- rlogTransformation(dds, blind = TRUE)

# blind = TRUE means:
# The transformation does not use information about the experimental design (i.e., group labels).
# So it's "blind" to the condition (like control vs. treatment).
# This is useful when you want unbiased exploration of your data, like visualizing overall structure or detecting batch effects.
# ➤ It's blind to groupings, so it's safe for exploratory data analysis like PCA or clustering

# The rlog transformation is used to:
# Stabilize variance across samples — especially for low-count genes
# Make the data more normally distributed (important for visualization and clustering)
# Reduce the influence of extreme outliers compared to using just log(count + 1)
# This is especially useful for things like PCA plots, Heatmaps, Sample distance calculations, Hierarchical clustering

datafromRld <- as.data.frame(assays(rld)[[1]]) 
# [ ] returns a sublist (still a list).
# to extract the first element from that list: 
# [[ ]] returns the actual element inside — the object itself.

# assays(rld) retrieves the assay data from the rld object.
# assay is a matrix like object
# rld is a type of SummarizedExperiment, and it stores data layers (like counts, normalized counts, etc.) in a list-like structure called assays.

#assays(rld)[[1]]  # extracts the first assay from the list — in this case, it's the rlog-transformed count matrix.
# It’s a numeric matrix: genes (rows) × samples (columns).

# as.data.frame(...)
# Converts that matrix into a data frame — which is easier to manipulate, filter, or plot using functions like ggplot2, dplyr, etc.

# make all positive --- add small value ############
### add contast if downstream steps forbids non-positive values, e.g. 
# •	you plan to take another log() on the rlog matrix,
# •	a plotting layer enforces a strictly positive scale,
# •	an external tool rejects negatives
## adding a global constant is acceptable (1e-5) 
## !! note: the shift does not change pairwise difference, but it moves color scales and means by the global constant
## i took this from Kls -- shouldn't bias much since for heatmaps we redefine zero center
add1 <- -min(datafromRld, na.rm = TRUE) + (1e-5) ##1e-5 = 0.00001 (10^-5)
datafromRld <- datafromRld + add1
range(datafromRld)

# adjust row and column names from original data with all columns and annotations
rownames(datafromRld) <- rownames(raw_counts)
colnames(datafromRld) <- colnames(raw_counts)[9:ncol(raw_counts)]

# control
head(datafromRld)

# visualization of normalization ------------------------------------------
#### BOX PLOTS
range1 <- range(log2(countdata + 1), na.rm = TRUE) #+1 to avoid log(0)
range2 <- range(datafromRld, na.rm = TRUE)
common_range <- range(c(range1, range2))

# boxplot before normalization
par(mfrow=c(1,1), mar=c(12,4,3,2)) #graphical settings for plot
boxplot(log2(countdata + 1), 
        col="blue", 
        las=3, 
        main="before normalization")
#boxplot for every column in countdata, transformed in log2 format with +1 (to avoid log(0), las=3 is turning the x-axis legend at 90 degrees, main is title

# boxplot after normalization
par(mfrow = c(1,1), font = 2, font.axis = 2, font.lab = 3, mar = c(12,4,3,2))
boxplot(datafromRld, 
        col = "red", 
        las = 3, 
        main = "after normalization (positive values)")


# done visualization of normalization -------------------------------------

# retrieve annotation and metadata  ---------------------------------------
# retrieve annotations from raw_counts table that match the DEG results and merge results
### look at your raw count table and add the columns back which you cut off for count data matrix generation 
names(raw_counts)
ann1 <- raw_counts[, c(1:8)] #first 8 columns (which were cut off for the matrix in the beginning)
lg_dat2 <- merge(ann1, datafromRld, by.x = "row.names", by.y = "row.names", all = FALSE)
lg_dat2[1:3,] #has merge worked?
row.names(lg_dat2) <- 1:nrow(lg_dat2) #row indizes from 1 onwards
cat("Verschmolzene Tabelle:", dim(lg_dat2), "\n")
head(lg_dat2)


# safe data -> norm_file --------------------------------------------------
final_table <- lg_dat2[, -1]  # Erste Spalte (Row.names) entfernen, could also be written like <- lg_dat2[,2:ncol(lg_dat2)]
head(final_table)
dir()
# write.table with final_table, file = "CHOOSE YOUR NAME, put date in the end.txt"
write.table(final_table, file = "norm_Wagner_CLP_TRPV1_2023_080525.txt", row.names = FALSE)
# done file saving --------------------------------------------------------

## look at norm genes 

norm_counts <- read.table("norm_Wagner_CLP_TRPV1_2023_080525.txt", header = TRUE, check.names = FALSE)

dim(norm_counts)
summary(norm_counts)
names(norm_counts)

head(norm_counts)

trpv1 <- norm_counts[norm_counts$gene_symbol == "Trpv1", ]

trpv1

# ============================================================
# SAMPLE-LEVEL HEATMAP AND PCA (EXPLORATORY)
# ============================================================

# Select top genes by mean normalized counts for demonstration heatmap
select <- order(rowMeans(counts(dds, normalized = TRUE)), decreasing = TRUE)[1:30]
hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)

# Euclidean distance heatmap between samples
distsRL <- dist(t(assay(rld)))
mat <- as.matrix(distsRL)
rownames(mat) <- colnames(mat) <- coldata$group
heatmap.2(mat, trace = "none", col = rev(hmcol), margins = c(6,8), cexRow = 0.7, cexCol = 0.7)

# PCA for quick structure check (blind = TRUE rlog above)
head(coldata)
list(rld)
names(colData(dds))
print(plotPCA(rld, intgroup = c("group")))
# print(plotPCA(rld, intgroup = c("batch")))  # not for presentation in this dataset

# ===========================================================================================
# NOTES
# • This script mirrors Kls' original with added explanations and new object names.
# • Logic, object names, filenames, and thresholds remain identical.
# • Global positive shift is minimal and acceptable for downstream tools that require >0 inputs.
# ===========================================================================================
