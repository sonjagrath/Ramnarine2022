#Supplemental Table 1 & 2 (DEG tables) Outlier removed (treatment only)
#Supplemental Table 3 & 4 - Lists from DEG tables as input to --> Flymine (GO terms) + REVIGO (GO term summary)

library("DESeq2")
library(readxl)
library(tidyr)

allcounts <- read_excel("../data/ALL_GCounts_M9sampleremoved.xlsx")


allcounts$newhead <- paste(allcounts$FBgn, allcounts$CGnum, allcounts$Chrom, allcounts$Length, sep=",")
names(allcounts)

counts<-data.frame(allcounts[c(52, 5:51)], row.names=1)
head(counts) #check table loaded correctly

dim(counts)
#[1] 13957    47

samples <- read.csv("../data/metadata_M9_excluded.csv")


### using one term ###
dds4 <- DESeqDataSetFromMatrix(countData = counts, 
                               colData = samples, 
                               design = ~treatment)
dds_4 <- DESeq(dds4, betaPrior = TRUE)


#contrasting all stress samples against control samples 

all_s_unorder <- results(dds_4, contrast = c("treatment","s", "c"), alpha = 0.05)
all_s <- all_s_unorder[order(all_s_unorder$padj),]
summary(all_s)

#out of 13860 with nonzero total read count
#adjusted p-value < 0.05
#LFC > 0 (up)       : 3434, 25%
#LFC < 0 (down)     : 2839, 20%
#outliers [1]       : 0, 0%
#low counts [2]     : 0, 0%
#(mean count < 0)
#[1] see 'cooksCutoff' argument of ?results
#[2] see 'independentFiltering' argument of ?results


all_s_noNA <- na.omit(all_s)
head(all_s)
summary(all_s_noNA)
#out of 13860 with nonzero total read count
#adjusted p-value < 0.05
#LFC > 0 (up)       : 3434, 25%
#LFC < 0 (down)     : 2839, 20%
#outliers [1]       : 0, 0%
#low counts [2]     : 0, 0%
#(mean count < 0)
#[1] see 'cooksCutoff' argument of ?results
#[2] see 'independentFiltering' argument of ?results

all_s_unorder <- as.data.frame(cbind(Gene_info = rownames(all_s_unorder), all_s_unorder))
all_s_unorder <- separate(data=all_s_unorder, col = "Gene_info", into = c("FBgn", "CGnum", "Chrom", "Length"), sep = ",")
all_s_unorder <- na.omit(all_s_unorder)
head(all_s_unorder)

all_s <- as.data.frame(cbind(Gene_info = rownames(all_s), all_s))
all_s <- separate(data=all_s, col = "Gene_info", into = c("FBgn", "CGnum", "Chrom", "Length"), sep = ",")
all_s <- na.omit(all_s)
head(all_s)

# save total DEG results file
write.csv(all_s_unorder, file = "../results/Stress_v_Control_ALL.csv", row.names = FALSE, quote = FALSE)

# generating file for up-regulated genes in all stress samples compared to all control
alls_up <- all_s[all_s[,6] > 0,]
alls_up <- alls_up[alls_up[,10] <= 0.05,]                                        
head (alls_up)
dim(alls_up)
#[1] 3434   10

write.csv(alls_up, file = "../results/ALL_svc_up.csv", row.names = FALSE, quote = FALSE)

# verify Gene IDs have been saved
all_up_list <- alls_up[[2]]
head(all_up_list)
#[1]"CG10383" "CG7059"  "CG4181"  "CG2914"  "CG11086" "CG4423"  

# generating file for down-regulated genes in all stress samples compared to all control
alls_down <-all_s[all_s[,6] < 0,]
alls_down <- alls_down[alls_down[,10] <= 0.05,]                                        
head (alls_down)
dim(alls_down)
#[1]  2839   10

write.csv(alls_down, file = "../results/ALL_svc_down.csv", row.names=FALSE, quote=FALSE)

all_down_list <- alls_down[[2]]
head(all_down_list)
