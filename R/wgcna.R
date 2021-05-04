################################################################################
# wgcna.R
################################################################################

## Weighted gene co-expression analysis (WGCNA)
## Author: Sonja Grath
## Input:
### count data: ALL_GCounts_M9sampleremoved.csv
### metadata: metadata.csv
## Last update: 2021-04-29

################################################################################

##############
# Preparation
##############

# General preparation (libraries, scripts, folders)

# Load necessary libraries

# DESeq2 analysis
library(DESeq2)

# WGCNA analysis
library(WGCNA)

# Data preparation
library(tidyverse)

# Functional annotation
# BiocManager::install("ReactomePA")
# BiocManager::install("org.Hs.eg.db")
# BiocManager::install("org.Dm.eg.db")
# ReactomePA appears to need org.Hs.eg.db as well
library(ReactomePA)
library(org.Hs.eg.db)
library(org.Dm.eg.db)

# Visualization
library(RColorBrewer)
library(flashClust)
library(dendextend)
library(cowplot)

# Load self-defined functions

source("R/own_functions.R")

# Data preparation

# Load count data
rawcounts <- read.csv("data/ALL_GCounts_M9sampleremoved.csv")

# Set rownames to first column (column has title and contains Flybase ids)
rownames(rawcounts) <- rawcounts[, 1]
rawcounts <- rawcounts[,-1]

# Load meta data
metadata <- read.csv("data/metadata.csv", row.names = 1)

# Use the match() function to reorder the columns of the raw counts
reorder_idx <- match(rownames(metadata), colnames(rawcounts[4:50, ]))
# Reorder the columns of the count data
reordered_rawcounts <- rawcounts[ , reorder_idx]

#######################################################
# Weighted Gene Co-Expression Network Analysis (WGCNA)
#######################################################

# Set parameters for WGCNA

# Define network type (parameter 'NetworkType' or 'type')
# current possibilities: "unsigned", "signed"
type = "signed"

# Create folders for output
dir.create("results")
dir.create("results/figures")
dir.create("results/files")
# used to save figures and files
figures <- "results/figures"
files <- "results/files"

# Data preparation

# Create DESeq2 object
dds <- DESeqDataSetFromMatrix(countData = reordered_rawcounts,
                              colData = metadata,
                              design = ~ treatment)

# Run the DESeq2 analysis
dds <- DESeq(dds)
resultsNames(dds)

# Create DESeq result
res <- results(dds,
               #specifiy contrast as c(covariate, val1, val2)
               contrast = c("treatment", "S", "C"),
               alpha = 0.1)

# Perform variance stabilization
vsd <- varianceStabilizingTransformation(dds)

# Extract a matrix of normalized counts from the S4 object vsd
vsdMat <- assay(vsd)

# Further, we want genes with enough information so that the variance can be
# reliably estimated and we do not get spurious correlations. Here, we can
# re-use the DESeq2 results: p-values of genes with too little information for
# DE-analysis were set to NA in the results object
nonNAgenes <- rownames(res[is.na(res$padj) == FALSE, ])
vsdMat <- vsdMat[nonNAgenes, ]

# Cleaning up the data: WGCNA wants samples as rows
mat <- t(vsdMat)
# We can use WGCNA functions that check for genes and samples with too little information
gsg <- goodSamplesGenes(mat)
# should be 'TRUE'
gsg$allOK
# TRUE

# Next, we cluster the samples (in contrast to clustering genes that will come later)
# to see if there are any obvious outliers.
sampleTree = hclust(dist(mat),
                    method = "average")
# Plot the sample tree: Open a graphic output window of size 12 by 9 inches
# The user should change the dimensions if the window is too large or too small.sizeGrWindow(12,9)
pdf(file = "sampleClustering.pdf", width = 12, height = 9);
par(cex = 0.6)
par(mar = c(0, 4, 2, 0))
plot(sampleTree, main = "Sample clustering to detect outliers",
     sub = "",
     xlab = "",
     cex.lab = 1.5,
     cex.axis = 1.5,
     cex.main = 2)

# Measure co-Expression of genes

# Choose a set of soft-thresholding powers (exponent in equation)
powers <- c(c(1:15), seq(from = 15, to = 21, by = 2))

# Allowing multi-threading with up to 20 threads
# Warning in allowWGCNAThreads: Requested number of threads is higher than number
# of available processors (or cores). Using too many threads may degrade code
# performance. It is recommended that the number of threads is no more than number
# of available processors.
# -> Not useful/necessary on laptop, relevant on servers
# allowWGCNAThreads(20)

# Note: 'type' was defined before
sft <- pickSoftThreshold(mat, powerVector = powers, networkType = type)

## Scale-free topology fit index as a function of the soft-thresholding power
rsq <- ggplot(sft$fitIndices,
              aes(x = Power,
                  y = -sign(slope)*SFT.R.sq,
                  label = Power)) +
    geom_text(col = "red", size = 3) +
    ylab(expression("Scale Free Topology Model Fit, signed"~R^{2})) +
    xlab("Soft Threshold (power)") +
    ggtitle("Scale independence") +
    geom_vline(xintercept = sft$powerEstimate, col = "red", alpha = 0.5) +
    geom_hline(yintercept = 0.90, alpha = 0.5) +
    theme_minimal() +
    theme(axis.title = element_text(size = 10)) +
    theme(axis.text = element_text(size = 10)) +
    # square plot area
    theme(aspect.ratio = 1)
# save figure as png, pdf and svg
ggsave(paste(figures, "/", "rsq.png", sep = ""),
       plot = rsq, width = 10, height = 10, units = "cm")
ggsave(paste(figures, "/", "rsq.pdf", sep = ""),
       plot = rsq, width = 10, height = 10, units = "cm")
svg(paste(figures, "/", "rsq.svg", sep = ""), width = 5, height = 5, pointsize = 10)
rsq
dev.off()

# Mean connectivity as a function of the soft-thresholding power
df <- sft$fitIndices %>% gather("type", "Connectivity",
                                median.k., mean.k., max.k.)
conne <- ggplot(df, aes(x = Power,
                        y = Connectivity,
                        label = Power,
                        col = type)) +
    geom_point() +
    scale_y_log10() +
    xlab("Soft Threshold (power)") +
    ggtitle("Connectivity") +
    theme_minimal() +
    theme(axis.title = element_text(size = 10)) +
    theme(axis.text = element_text(size = 10)) +
    # square plot area
    theme(aspect.ratio = 1)
# save figure as png, pdf and svg
ggsave(paste(figures, "/", "conne.png", sep = ""),
       plot = conne, width = 10, height = 10, units = "cm")
ggsave(paste(figures, "/", "conne.pdf", sep = ""),
       plot = conne, width = 10, height = 10, units = "cm")
svg(paste(figures, "/", "conne.svg", sep = ""), width = 5, height = 5, pointsize = 10)
conne
dev.off()

# Display both plots next to each other (package cowplot)
rsq_conne <- cowplot::plot_grid(rsq, conne, nrow = 1,
                   labels = LETTERS[1:2],
                   rel_widths = c(0.5, 0.5))
# save figure as png, pdf and svg
ggsave(paste(figures, "/", "rsq_conne.png", sep = ""),
       plot = rsq_conne, width = 20, height = 10, units = "cm")
ggsave(paste(figures, "/", "rsq_conne.pdf", sep = ""),
       plot = rsq_conne, width = 20, height = 10, units = "cm")
svg(paste(figures, "/", "rsq_conne.svg", sep = ""), width = 10, height = 5, pointsize = 10)
rsq_conne
dev.off()

### TO CHECK:
### Default for soft-thresholding power?
### DEFAULT: softPower = 5
### Good value?!

# Get mean, median and max connectivity for sft
# mean
sft$fitIndices$mean.k.[sft$powerEstimate]
# [1] 784.5154
# max
sft$fitIndices$max.k.[sft$powerEstimate]
# [1] 1629.21
# median
sft$fitIndices$median.k.[sft$powerEstimate]
# [1] 714.5975
# Set soft power to the estimated value
softPower <- sft$powerEstimate

# ALTERNATIVE: Set soft power to 12 according to guidelines by Langfelder and Horvath (2017, FAQ)
softPower <- 12

# Get mean, median and max connectivity for power of 12
# mean
sft$fitIndices$mean.k.[softPower]
# [1] 73.1647
# max
sft$fitIndices$max.k.[softPower]
# [1] 499.9192
# median
sft$fitIndices$median.k.[softPower]
# [1] 41.37947


# calculate connectivity
# Note: NetworkType is called type here and was predefined above
k <- softConnectivity(mat, power = softPower, type = type)
names(k) <- colnames(mat)

# look at connectivity in histogram
connectivityHist <- ggplot(data.frame(k = k), aes(x = k)) +
    geom_histogram(fill = NA, col = 1) +
    theme_minimal() +
    theme(axis.title = element_text(size = 10)) +
    theme(axis.text = element_text(size = 10)) +
    # square plot area
    theme(aspect.ratio = 1)
ggsave(paste(figures, "/", "connectivityHist.png", sep = ""),
       plot = connectivityHist, width = 20, height = 10, units = "cm")
ggsave(paste(figures, "/", "connectivityHist.pdf", sep = ""),
       plot = connectivityHist, width = 20, height = 10, units = "cm")
svg(paste(figures, "/", "connectivityHist.svg", sep = ""), width = 10, height = 5, pointsize = 10)
connectivityHist
dev.off()

# plot the fit to the scale free topology
png(paste(figures, "/", "scaleFreePlot.png", sep = ""))
scaleFreePlot(k)
dev.off()
pdf(paste(figures, "/", "scaleFreePlot.pdf", sep = ""))
scaleFreePlot(k)
dev.off()
svg(paste(figures, "/", "scaleFreePlot.svg", sep = ""))
scaleFreePlot(k)
dev.off()

# Power transformation to get adjacency matrix
# The adjacency matrix is a metric to determine which genes j and i are connected
# Note: NetworkType (type) is defined as type above
adjacency <- adjacency(mat, power = softPower, type = type)

# Moreover, we can estimate how many connections gene j and i have in common
# This is called topological overlap matrix (TOM), a different measure that is more robust than adjacency
# Calculate the dissimilarity matrix based on TOM
TOM <- TOMsimilarity(adjacency)
colnames(TOM) <- colnames(mat)
rownames(TOM) <- colnames(mat)
saveRDS(TOM, file = paste(files, "/", "TOM.rds", sep = ""))
saveRDS(adjacency, file = paste(files, "/", "adjacency.rds", sep = ""))
collectGarbage()

# Identify co-expression modules
require(flashClust)
disTOM <- 1-TOM
geneTree <- flashClust(as.dist(disTOM), method = "average")
# Module identification using dynamic tree cut
# each module should have a meaningful size, we require 30 genes
minModuleSize <- 30
# assign a cluster id to each gene
# ids are ususally ordered according to size, starting with the biggest cluster
dynamicMods <- cutreeDynamic(dendro = geneTree,
                             distM = disTOM,
                             cutHeight = 0.97,
                             minClusterSize = minModuleSize)
dynamicColors <- labels2colors(dynamicMods)
table(dynamicColors)

# Simply cutting the TOM-tree tends to over-split similar clusters. Therefore,
# we now identify clusters that have a very similar average expression pattern
# and merge them. This is achieved using principal component analysis.

# Calculate eigengenes
MEList <- moduleEigengenes(mat, colors = dynamicColors)
MEs <- MEList$eigengenes

# Calculate dissimilarity of module eigengenes
MEDiss <- 1-cor(MEs)

# Cluster module eigengenes
METree<- flashClust(as.dist(MEDiss), method = "average")

# Plot the result
require(dendextend)
d <- as.dendrogram(METree)
lab <- METree$labels[order.dendrogram(d)]
d <- color_labels(d ,labels=lab, col = sub("ME","",lab) )

png(paste(figures, "/", "clusteringModuleEigengenes.png", sep = ""))
par(mar = c(5,4,4,6)) # make figure margins bigger
plot(d,
     main = "Clustering of module eigengenes",
     horiz = TRUE,
     xlab = "",
     sub = "",
     cex = 0.5)
abline(v = 0.2, col = 3)
abline(v = 0.15, col = 2)
dev.off()

pdf(paste(figures, "/", "clusteringModuleEigengenes.pdf", sep = ""))
par(mar = c(5,4,4,6)) # make figure margins bigger
plot(d,
     main = "Clustering of module eigengenes",
     horiz = TRUE,
     xlab = "",
     sub = "",
     cex = 0.5)
abline(v = 0.2, col = 3)
abline(v = 0.15, col = 2)
dev.off()

# Based on the hierarchical clustering of the TOM-module eigengenes, we decide
# to merge modules for which the eigengenes have a distance smaller than 0.2,
# which roughly corresponds to a correlation of 0.8.
MEDissThres = 0.2

# Call automatic merging function
merge <- mergeCloseModules(mat,
                           dynamicColors,
                           cutHeight = MEDissThres,
                           verbose = 0)
# mergeCloseModules: Merging modules whose distance is less than 0.2
# Calculating new MEs...
mergedColors <- merge$colors

# Plot gene dendrogram with new & old module association
# WGCNA function
dendroAndColors <- plotDendroAndColors(geneTree,
                    cbind(dynamicColors, mergedColors),
                    c("DynamicTree Cut", "Merged dynamic"),
                    dendroLabels = FALSE,
                    hang = 0.03,
                    addGuide = TRUE,
                    guideHang = 0.05)

png(paste(figures, "/", "clusterDendrogram.png", sep = ""))
dendroAndColors
dev.off()

# associate finale module colors
moduleColors <- mergedColors
write.csv(moduleColors, file = paste(files, "/", "moduleColors.csv", sep = ""))

colorOrder <- c("grey", standardColors(50))
moduleLabels <- match(moduleColors, colorOrder)-1
MEs <- merge$newMEs

# We draw a heatmap displaying the relationship of the identified modules using
# the topological overlap dissimilarity matrix (TOM). Darker shades of yellow
# indicate a closer expression neighborhood.
# This figure becomes too big as a pdf, therefore, we save it as jpeg.
jpeg(filename = paste(figures, "/", "tom.jpg", sep = ""), quality = 100)
TOMplot(disTOM, geneTree, as.character(mergedColors))
dev.off()

# How many genes are in the final merged modules?
df <- data.frame(dc = moduleColors) %>% dplyr::count(dc)
mergedModuleSizes <- ggplot(df , aes(x = reorder(dc, n),
                y = n,
                fill = dc)) +
    geom_bar(stat = "identity") +
    coord_flip() +
    scale_fill_manual(values = as.character(df$dc)) +
    ylim(c(0, 7500)) +
    theme_minimal() +
    ggtitle("Merged Module sizes") +
    ylab("Number of genes") +
    xlab("") +
    geom_text(aes(label=n), color="black", size=3, hjust = -0.1, vjust = 0.5) +
    theme(legend.position = 'None',
          axis.title = element_text(size = 10),
          axis.text = element_text(size = 10),
          aspect.ratio = 1) # square plot area
ggsave(paste(figures, "/", "mergedModuleSizes.png", sep = ""),
       plot = mergedModuleSizes, width = 10, height = 10, units = "cm")
ggsave(paste(figures, "/", "mergedModuleSizes.pdf", sep = ""),
       plot = mergedModuleSizes, width = 10, height = 10, units = "cm")
# save svg version
svg(paste(figures, "/", "mergedModuleSizes.svg", sep = ""), width = 10, height = 10)
mergedModuleSizes
dev.off()

# QUESTION:
# Which modules are actually interesting?
# Associate modules with sample information
# We want to know which modules are associated with interesting sample characteristics.
# This is achieved by correlating the sample property, i.e. stress condition or control,
# to to module eigene
# The function to do this (self-written, loaded with own_functions.R): Module2SampleFeat

# Define numbers of genes and samples
nGenes = ncol(mat)
nSamples = nrow(mat)
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(mat, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
cat.inf <- data.frame(colData(dds)) %>% dplyr::select(treatment)

# Associate categorical varibales with module eigengenes of modules
moduleWithSampleInfo <- Module2SampleFeat(inf = cat.inf,
                                          MEs = MEs,
                                          type = "categorical",
                                          lab = "")
ggsave(paste(figures, "/", "moduleWithSampleInfo.png", sep = ""),
       plot = moduleWithSampleInfo, width = 10, height = 10, units = "cm")
ggsave(paste(figures, "/", "moduleWithSampleInfo.pdf", sep = ""),
       plot = moduleWithSampleInfo, width = 10, height = 10, units = "cm")
svg(paste(figures, "/", "moduleWithSampleInfo.svg", sep = ""), width = 10, height = 5, pointsize = 10)
moduleWithSampleInfo
dev.off()

# Colors indicate the slope (correlation coefficient) for the sample variable.
# The printed values are the p-values.
# TO CHECK:
# Genes of the yellow, turquoise, and greenyellow modules appear to be down-regulated
# under stressed condition, while genes of the pink and blue modules are up-regulated.

# Associate modules with gene information
kim <- intramodularConnectivity(adjacency, moduleColors, scaleByMax = TRUE)

# save as csv (-> needed for functional annotation analysis)
write.csv(kim, file = paste(files, "/", "kim.csv", sep = ""))


# QUESTION:
# Are significantly down-regulated genes enriched in one of the modules?
# To answer this question, we use Fisher's exact test
deseq2res <- data.frame(gene = rownames(res), res)
module2Gene <- data.frame(module = moduleColors, gene = colnames(mat))
# save as csv file (-> needed for functional annotation analysis)
write.csv(module2Gene, file = paste(files, "/", "module2Gene.csv", sep = ""))

## Are significantly down-regulated genes enriched in any of the modules?
de.enrich.down <- deseq2res %>%
    inner_join(module2Gene, by = "gene") %>%
    mutate(sigDown = ifelse(log2FoldChange < 0 & padj < 0.1, "down", "ns")) %>%
    dplyr::count(sigDown, module) %>%
    group_by(sigDown) %>%
    mutate(N = sum(n)) %>%
    ungroup() %>%
    mutate(N = N-n) %>%
    gather(stype, nn, N, n) %>%
    group_by(module) %>%
    do(tidy(fisher.test(matrix(.$nn, 2, 2), alternative = "greater")))
de.enrich.down %>% filter(p.value < 0.01)

## Are significantly up-regulated genes enriched in any of the modules?

de.enrich.up <- deseq2res %>%
    inner_join(module2Gene, by = "gene") %>%
    mutate(sigUp = ifelse(log2FoldChange > 0 & padj < 0.1, "up", "ns")) %>%
    dplyr::count(sigUp, module) %>%
    group_by(sigUp) %>%
    mutate(N = sum(n)) %>%
    ungroup() %>%
    mutate(N = N-n) %>%
    gather(stype, nn, N, n) %>%
    group_by(module) %>%
    do(tidy(fisher.test(matrix(.$nn, 2, 2), alternative = "greater")))
de.enrich.up %>% filter(p.value < 0.01)


################################################################################
# Functional Module Annotation
################################################################################

# Load data from WGCNA

module2Gene <- read.table(paste(files, "module2Gene.csv", sep = "/"),
                          header = TRUE,
                          row.names = 1,
                          sep = ",")
kim <- read.table(paste(files, "kim.csv", sep = "/"),
                  header = TRUE,
                  row.names = 1,
                  sep = ",")

# Annotate data

g4reac <- AnnotationDbi::select(org.Dm.eg.db,
                                keys = colnames(mat),
                                keytype = "ENSEMBL",
                                columns = c("ENSEMBL", "ENTREZID")) %>%
    dplyr::left_join(module2Gene, by = c("ENSEMBL" = "gene"))
# save as csv file (-> needed for functional annotation analysis)
write.csv(g4reac, file = paste("../", files, "/", "g4reac.csv", sep = ""))

## Create gene list (necessary for plot annotations with fold change of genes)
## For gene set enrichment analysis, we need a ranked list of genes which we
## will call geneList.
## The geneList contains three features:
### 1. numeric vector: fold change or other type of numerical variable
### 2. named vector: every number was named by the corresponding gene ID
### 3. sorted vector: number should be sorted in decreasing order

## the rownames are IDs, res is the results data from the DESeq2 analysis above
id <- rownames(res)
res <- cbind(res, id)
data <- AnnotationDbi::select(org.Dm.eg.db,
                              keys = id,
                              keytype = "ENSEMBL",
                              columns = c("ENSEMBL", "ENTREZID"))
resultTable <- merge(res, data, by.x="id", by.y="ENSEMBL")
## feature 1: numeric vector (fold change in column 3)
geneList <- resultTable[ , 3]
## feature 2: named vector (ENTREZ id in column 8)
names(geneList) <- as.character(resultTable[ , 8])
## feature 3: decreasing order
geneList <- sort(geneList, decreasing = TRUE)

################################
# Functional module annotation
################################

# Complete visualization for interesting modules
# We defined the modules pink, blue, turquoise, magenta and
# green as interesting

# Set path for subfolders that contain figures for each module
my_path <- paste(getwd(), "/", figures, sep = "")

# Create folders for these modules
create_color_folder(my_path, "pink")
create_color_folder(my_path, "blue")
create_color_folder(my_path, "turquoise")
create_color_folder(my_path, "magenta")
create_color_folder(my_path, "green")

# Create plots (self-defined functions have to be loaded)

# yellow
complete_visualization(color = "pink",
                       nr_categories = 15,
                       type = "dotplot",
                       path = my_path)
complete_visualization(color = "pink",
                       nr_categories = 15,
                       type = "emapplot",
                       path = my_path)
complete_visualization(color = "pink",
                       nr_categories = 5,
                       type = "cnetplot",
                       path = my_path)

# blue
complete_visualization(color = "blue",
                       nr_categories = 15,
                       path = my_path)
complete_visualization(color = "blue",
                       nr_categories = 15,
                       type = "emapplot",
                       path = my_path)
complete_visualization("blue",
                       nr_categories = 5,
                       type = "cnetplot",
                       path = my_path)


# turquoise
complete_visualization(color = "turquoise",
                       nr_categories = 15,
                       type = "dotplot",
                       path = my_path)
complete_visualization(color = "turquoise",
                       nr_categories = 15,
                       type = "emapplot",
                       path = my_path)
complete_visualization(color = "turquoise",
                       nr_categories = 5,
                       type = "cnetplot",
                       path = my_path)

# magenta
complete_visualization(color = "magenta",
                       nr_categories = 15,
                       type = "dotplot",
                       path = my_path)
complete_visualization(color = "magenta",
                       nr_categories = 15,
                       type = "emapplot",
                       path = my_path)
complete_visualization(color = "magenta",
                       nr_categories = 5,
                       type = "cnetplot",
                       path = my_path)

# green
complete_visualization(color = "green",
                       nr_categories = 15,
                       type = "dotplot",
                       path = my_path)
complete_visualization(color = "green",
                       nr_categories = 15,
                       type = "emapplot",
                       path = my_path)
complete_visualization(color = "green",
                       nr_categories = 5,
                       type = "cnetplot",
                       path = my_path)

#################
# Session Info
#################

# session info is saved in a .txt file into the current folder
sink("session_info.txt")
devtools::session_info()
sink()

