################################################################################
# Figure_S1.R
################################################################################

## Principal Component Analysis - all samples
## Authors: Timothy Ramnarine, Sonja Grath
## Input:
### count data: ALL_GCounts.csv
### metadata: metadata_all.csv
## Last update: 2021-05-12

################################################################################

##############
# Preparation
##############

# General preparation (libraries, scripts, folders)

# Load necessary libraries

# DESeq2 analysis
library(DESeq2)

# Data preparation
library(tidyverse)

# Visualization
library(ggplot2)

# Data preparation

# Load count data
rawcounts <- read.csv("../data/ALL_GCounts.csv")

# Set rownames to first column (column has title and contains Flybase ids)
rownames(rawcounts) <- rawcounts[, 1]
rawcounts <- rawcounts[,-1]

# Load meta data
metadata <- read.csv("../data/metadata_all.csv", row.names = 1)

# Use the match() function to reorder the columns of the raw counts
reorder_idx <- match(rownames(metadata), colnames(rawcounts[4:51, ]))
# Reorder the columns of the count data
reordered_rawcounts <- rawcounts[ , reorder_idx]

# Create DESeq2 object
dds <- DESeqDataSetFromMatrix(countData = reordered_rawcounts,
                              colData = metadata,
                              design = ~ treatment)

# Run the DESeq2 analysis
dds <- DESeq(dds, betaPrior = TRUE)
resultsNames(dds)

# Perform variance stabilization
vsd <- vst(dds, blind = FALSE) ## variance stabilizing transformation, better because variance is the same across high and low read counts
pcaData <- plotPCA(vsd, 
                   intgroup = c("line.treat", "population", "pop.line"), 
                   returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

pcaData$population <- factor(pcaData$population, levels = c("M", "C2", "KL", "ZI"))
pcaData$line.treat <- factor(pcaData$line.treat, levels = c("MC", "C2C", "KLC","ZIC",  "MS", "C2S",  "KLS",  "ZIS"))

pca <- ggplot(pcaData, aes(PC1, PC2 )) +
  geom_point(aes(fill = line.treat, 
                 shape = line.treat, 
                 color = line.treat),
             size = 1.7, 
             stroke = 1) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  ggtitle("") +
  coord_fixed() +
  theme_classic() +
  scale_color_manual(values = c("#440154FF", "#1B9E77", "#66A61E", "#D95F02", "black","black","black","black")) +
  scale_fill_manual(values = c("white", "white", "white", "white","#440154FF", "#1B9E77", "#66A61E", "#D95F02"), 
                    name = "Control                  Stress",
                    labels = c("     - Munich -",  "     - Cyprus -","- Kuala Lumpur -", "     - Zambia -",
                               "",  "","", ""))+ 
  scale_shape_manual(name = "Control                  Stress", 
                     labels = c("     - Munich -",  "     - Cyprus -","- Kuala Lumpur -", "     - Zambia -",
                                "",  "","", ""),
                     values = c(22, 24, 21, 23, 22, 24, 21, 23))+
  theme(text = element_text(size = 10))+
  theme(axis.text.x = element_text(size = 10))+
  theme(axis.text.y = element_text(size = 10))+
  theme(legend.position = "none")+
  guides(color=guide_legend(ncol =2)) +
  theme(panel.background = element_blank())+
  theme(plot.title = element_text(hjust = 0.5))
pca

