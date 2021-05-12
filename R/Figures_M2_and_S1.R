library("DESeq2")
library(readxl)
library(ggplot2)

##Supplemental Figure (S1)- PCA (full dataset, using three factors)##

allcounts <- read_excel("../data/ALL_GCounts_excel_import.xlsx")
allcounts$newhead <- paste(allcounts$FBgn, allcounts$CGnum, allcounts$Chrom, allcounts$Length, sep=",")
names(allcounts)
counts<-data.frame(allcounts[c(53, 5:52)], row.names=1)
head(counts) #check table loaded correctly
dim(counts)
#[1] 13957    48

#create sample information table#
samples <- data.frame(row.names=c("C2D3c_H3","C2D3c_H4","C2D3s_H1","C2D3s_H2",
                                  "C2N1c_G3","C2N1c_G4","C2N1s_G1","C2N1s_G2",
                                  "KL01c_E3","KL01c_E4","KL01s_E1","KL01s_E2", 
                                  "KL02c_F3","KL02c_F4","KL02s_F1", "KL02s_F2",
                                  "M12D3c_C3","M12D3c_C4","M12D3s_C1","M12D3s_C2", 
                                  "M12N3c_D3","M12N3c_D4","M12N3s_D1","M12N3s_D2", 
                                  "M9D3c_B3", "M9D3c_B4","M9D3s_B1", "M9D3s_B2",  
                                  "M9N1c_A3","M9N1c_A4", " M9N1s_A1", "M9N1s_A2", 
                                  "ZI197c_I2","ZI197c_I3","ZI197s_I1",
                                  "ZI254c_J3","ZI254c_J4","ZI254s_J1","ZI254s_J2","ZI254s_J5", 
                                  "ZI273c_K3","ZI273c_K4","ZI273s_K1","ZI273s_K2","ZI273s_K5",
                                  "ZI418c_L2","ZI418c_L3","ZI418s_L1"), 
                      sample_names=as.factor(c("H3", "H4", "H1", "H2",
                                               "G3", "G4", "G1", "G2",
                                               "E3", "E4", "E1", "E2",
                                               "F3", "F4", "F1", "F2",
                                               "C3", "C4", "C1", "C2",
                                               "D3", "D4", "D1", "D2",
                                               "B3", "B4", "B1", "B2",
                                               "A3", "A4", "A1", "A2",
                                               "I2", "I3", "I1",
                                               "J3", "J4", "J1", "J2", "J5",
                                               "K3", "K4", "K1", "K2", "K5",
                                               "L2", "L3", "L1")),
                      pop.geno.treat=as.factor(c("C2DelC", "C2DelC", "C2DelS", "C2DelS", 
                                                 "C2NonC", "C2NonC", "C2NonS", "C2NonS",
                                                 "KL1DelC", "KL1DelC",  "KL1DelS", "KL1DelS",
                                                 "KL2NonC",  "KL2NonC","KL2NonS", "KL2NonS",
                                                 "M12DelC", "M12DelC", "M12DelS", "M12DelS", 
                                                 "M12NonC","M12NonC", "M12NonS", "M12NonS",
                                                 "M9DelC", "M9DelC", "M9DelS", "M9DelS",
                                                 "M9NonC", "M9NonC", "M9NonS","M9NonS",
                                                 "ZINonC", "ZINonC", "ZINonS",
                                                 "ZINonC", "ZINonC", "ZINonS", "ZINonS", "ZINonS",
                                                 "ZINonC","ZINonC", "ZINonS", "ZINonS", "ZINonS",
                                                 "ZINonC", "ZINonC", "ZINonS")),
                      population=as.factor(c("Cyprus","Cyprus","Cyprus","Cyprus",
                                             "Cyprus","Cyprus","Cyprus","Cyprus" ,
                                             "Kuala Lumpur","Kuala Lumpur", "Kuala Lumpur","Kuala Lumpur",
                                             "Kuala Lumpur", "Kuala Lumpur", "Kuala Lumpur", "Kuala Lumpur",
                                             "Munich", "Munich", "Munich", "Munich",
                                             "Munich", "Munich", "Munich", "Munich",
                                             "Munich", "Munich", "Munich", "Munich",
                                             "Munich", "Munich", "Munich", "Munich",
                                             "Zambia","Zambia","Zambia",
                                             "Zambia","Zambia","Zambia","Zambia","Zambia",
                                             "Zambia","Zambia","Zambia","Zambia","Zambia",
                                             "Zambia","Zambia","Zambia")), 
                      genotype=as.factor(c("del", "del", "del", "del",
                                           "non", "non", "non", "non",
                                           "del", "del", "del", "del",
                                           "non", "non", "non", "non",
                                           "del", "del", "del", "del",
                                           "non", "non", "non", "non",
                                           "del", "del", "del", "del",
                                           "non", "non", "non", "non",
                                           "non", "non", "non",
                                           "non", "non", "non", "non", "non",
                                           "non", "non", "non", "non", "non",
                                           "non", "non", "non")),
                      treatment=as.factor(c("control","control","stress","stress",
                                            "control","control","stress","stress",
                                            "control","control","stress","stress",
                                            "control","control","stress","stress",
                                            "control","control","stress","stress",
                                            "control","control","stress","stress",
                                            "control","control","stress","stress",
                                            "control","control","stress","stress",
                                            "control","control","stress",
                                            "control","control","stress","stress","stress", 
                                            "control","control","stress","stress","stress", 
                                            "control","control","stress")),
                      geno.treat=as.factor(c("DelC", "DelC", "DelS", "DelS", 
                                             "NonC", "NonC", "NonS", "NonS",
                                             "DelC", "DelC", "DelS", "DelS",
                                             "NonC", "NonC", "NonS", "NonS",
                                             "DelC", "DelC", "DelS", "DelS",
                                             "NonC", "NonC", "NonS", "NonS",
                                             "DelC", "DelC", "DelS", "DelS",
                                             "NonC", "NonC", "NonS", "NonS",
                                             "NonC", "NonC", "NonS",
                                             "NonC", "NonC", "NonS", "NonS", "NonS",
                                             "NonC", "NonC", "NonS", "NonS", "NonS",
                                             "NonC", "NonC", "NonS")),
                      line.treat=as.factor(c("C2C", "C2C", "C2S", "C2S", 
                                             "C2C", "C2C", "C2S", "C2S",
                                             "KLC", "KLC", "KLS", "KLS",
                                             "KLC", "KLC", "KLS", "KLS",
                                             "MC", "MC", "MS", "MS",
                                             "MC", "MC", "MS", "MS",
                                             "MC", "MC", "MS", "MS",
                                             "MC", "MC", "MS", "MS",
                                             "ZIC", "ZIC", "ZIS",
                                             "ZIC", "ZIC", "ZIS", "ZIS", "ZIS",
                                             "ZIC", "ZIC", "ZIS", "ZIS", "ZIS",
                                             "ZIC", "ZIC", "ZIS")),
                      pop.line=as.factor(c("C2Del", "C2Del", "C2Del", "C2Del", 
                                           "C2Non", "C2Non", "C2Non", "C2Non",
                                           "KL1Del", "KL1Del",  "KL1Del", "KL1Del",
                                           "KL2Non",  "KL2Non","KL2Non", "KL2Non",
                                           "M12Del", "M12Del", "M12Del", "M12Del", 
                                           "M12Non","M12Non", "M12Non", "M12Non",
                                           "M9Del", "M9Del", "M9Del", "M9Del",
                                           "M9Non", "M9Non", "M9Non", "M9Non",
                                           "ZI197Non", "ZI197Non", "ZI197Non",
                                           "ZI254Non", "ZI254Non", "ZI254Non", "ZI254Non", "ZI254Non",
                                           "ZI273Non","ZI273Non", "ZI273Non", "ZI273Non", "ZI273Non",
                                           "ZI418Non", "ZI418Non", "ZI418Non")))

# must be in same order as column of count matrix!
samples



dds_full <- DESeqDataSetFromMatrix(countData = counts, 
                               colData = samples, 
                               design = ~population + genotype + treatment)
##note: levels of factors in the design contain characters other than
#letters, numbers, '_' and '.'. It is recommended (but not required) to use
#only letters, numbers, and delimiters '_' or '.', as these are safe characters
#for column names in R. [This is a message, not a warning or an error]

#this is due to some labels in "population" having spaces


dds_full 
#class: DESeqDataSet 
#dim: 13957 48 
#metadata(1): version
#assays(1): counts
#rownames(13957): FBgn0032191,CG5734,2L,2205 FBgn0053476,CG33476,2R,498 ... FBgn0053814,CG33814,2L,493
#FBgn0034310,CG5733,2R,2146
#rowData names(0):
#  colnames(48): C2D3c_H3 C2D3c_H4 ... ZI418c_L3 ZI418s_L1
#colData names(7): pop.geno.treat population ... line.treat pop.line

deseq_full <- DESeq(dds_full, betaPrior = TRUE)

resultsNames(deseq_full)
#[1] "Intercept"              "populationCyprus"       "populationKuala.Lumpur" "populationMunich"      
#[5] "populationZambia"       "genotypedel"            "genotypenon"            "treatmentcontrol"      
#[9] "treatmentstress"  

vsd_full <- vst(deseq_full, blind = FALSE)

pcaData1 <- plotPCA(vsd_full, intgroup = c("line.treat", "population", "genotype", "sample_names"), returnData = TRUE)
percentVar1 <- round(100 * attr(pcaData1, "percentVar"))

pcaData1$population <- factor(pcaData1$population,levels = c("M", "C2", "KL", "ZI"))
pcaData1$line.treat <- factor(pcaData1$line.treat,levels = c("MC", "C2C", "KLC","ZIC",  "MS", "C2S",  "KLS",  "ZIS"))
pcaData1$genotype <- factor(pcaData1$genotype,levels = c("del", "non"))

pca_out <- ggplot(pcaData1, aes(PC1, PC2)) +
  geom_point(aes(fill= line.treat, shape = line.treat, color = line.treat),size=1.7, stroke = 1) +
  xlab(paste0("PC1: ",percentVar1[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar1[2],"% variance")) +
  ggtitle("") +
  coord_fixed() +
  theme_classic() +
  scale_color_manual(values = c("#440154FF", "#1B9E77", "#66A61E", "#D95F02", "black","black","black","black"))+
  scale_fill_manual(values = c("white", "white", "white", "white","#440154FF", "#1B9E77", "#66A61E", "#D95F02"))+ 
  scale_shape_manual(values = c(22, 24, 21, 23, 22, 24, 21, 23))+
  theme(text = element_text(size = 10))+
  theme(axis.text.x = element_text(size = 10))+
  theme(axis.text.y = element_text(size = 10))+
  theme(legend.position = "none") +
  guides(color=guide_legend(ncol=2)) +
  theme(panel.background = element_blank())+
  theme(plot.title = element_text(hjust = 0.5))
pca_out 




##Main Figure (M2) - PCA##

allcounts <- read_excel("../data/ALL_GCounts_M9sampleremoved.xlsx")

allcounts$newhead <- paste(allcounts$FBgn, allcounts$CGnum, allcounts$Chrom, allcounts$Length, sep = ",")
names(allcounts)

counts <- data.frame(allcounts[c(52, 5:51)], row.names = 1)
head(counts) #check table loaded correctly

dim(counts)
#[1] 13957    47

samples <- data.frame(row.names=c("C2D3c_H3","C2D3c_H4","C2D3s_H1","C2D3s_H2",
                                  "C2N1c_G3","C2N1c_G4","C2N1s_G1","C2N1s_G2",
                                  "KL01c_E3","KL01c_E4","KL01s_E1","KL01s_E2", 
                                  "KL02c_F3","KL02c_F4","KL02s_F1", "KL02s_F2",
                                  "M12D3c_C3","M12D3c_C4","M12D3s_C1","M12D3s_C2", 
                                  "M12N3c_D3","M12N3c_D4","M12N3s_D1","M12N3s_D2", 
                                  "M9D3c_B3", "M9D3c_B4","M9D3s_B1", "M9D3s_B2",  
                                  "M9N1c_A3","M9N1c_A4", "M9N1s_A2", 
                                  "ZI197c_I2","ZI197c_I3","ZI197s_I1",
                                  "ZI254c_J3","ZI254c_J4","ZI254s_J1","ZI254s_J2","ZI254s_J5", 
                                  "ZI273c_K3","ZI273c_K4","ZI273s_K1","ZI273s_K2","ZI273s_K5",
                                  "ZI418c_L2","ZI418c_L3","ZI418s_L1"), 
                      pop.geno.treat=as.factor(c("C2DelC", "C2DelC", "C2DelS", "C2DelS", 
                                                 "C2NonC", "C2NonC", "C2NonS", "C2NonS",
                                                 "KL1DelC", "KL1DelC",  "KL1DelS", "KL1DelS",
                                                 "KL2NonC",  "KL2NonC","KL2NonS", "KL2NonS",
                                                 "M12DelC", "M12DelC", "M12DelS", "M12DelS", 
                                                 "M12NonC","M12NonC", "M12NonS", "M12NonS",
                                                 "M9DelC", "M9DelC", "M9DelS", "M9DelS",
                                                 "M9NonC", "M9NonC", "M9NonS",
                                                 "ZINonC", "ZINonC", "ZINonS",
                                                 "ZINonC", "ZINonC", "ZINonS", "ZINonS", "ZINonS",
                                                 "ZINonC","ZINonC", "ZINonS", "ZINonS", "ZINonS",
                                                 "ZINonC", "ZINonC", "ZINonS")),
                      population=as.factor(c("C2","C2","C2","C2",
                                             "C2","C2","C2","C2" ,
                                             "KL","KL", "KL","KL",
                                             "KL", "KL", "KL", "KL",
                                             "M", "M", "M", "M",
                                             "M", "M", "M", "M",
                                             "M", "M", "M", "M",
                                             "M", "M", "M",
                                             "ZI","ZI","ZI",
                                             "ZI","ZI","ZI","ZI","ZI",
                                             "ZI","ZI","ZI","ZI", "ZI",
                                             "ZI","ZI","ZI" )), 
                      genotype=as.factor(c("del", "del", "del", "del",
                                           "non", "non", "non", "non",
                                           "del", "del", "del", "del",
                                           "non", "non", "non", "non",
                                           "del", "del", "del", "del",
                                           "non", "non", "non", "non",
                                           "del", "del", "del", "del",
                                           "non", "non", "non",
                                           "non", "non", "non",
                                           "non", "non", "non", "non", "non",
                                           "non", "non", "non", "non", "non",
                                           "non", "non", "non")),
                      treatment=as.factor(c("c","c","s","s",
                                            "c","c","s","s",
                                            "c","c","s","s", 
                                            "c","c","s","s",
                                            "c","c","s","s", 
                                            "c","c","s","s", 
                                            "c","c","s","s",  
                                            "c","c","s", 
                                            "c","c","s",
                                            "c","c","s","s","s", 
                                            "c","c","s","s","s",
                                            "c","c","s")), 
                      geno.treat=as.factor(c("DelC", "DelC", "DelS", "DelS", 
                                             "NonC", "NonC", "NonS", "NonS",
                                             "DelC", "DelC", "DelS", "DelS",
                                             "NonC", "NonC", "NonS", "NonS",
                                             "DelC", "DelC", "DelS", "DelS",
                                             "NonC", "NonC", "NonS", "NonS",
                                             "DelC", "DelC", "DelS", "DelS",
                                             "NonC", "NonC", "NonS",
                                             "NonC", "NonC", "NonS",
                                             "NonC", "NonC", "NonS", "NonS", "NonS",
                                             "NonC", "NonC", "NonS", "NonS", "NonS",
                                             "NonC", "NonC", "NonS")),
                      line.treat=as.factor(c("C2C", "C2C", "C2S", "C2S", 
                                             "C2C", "C2C", "C2S", "C2S",
                                             "KLC", "KLC", "KLS", "KLS",
                                             "KLC", "KLC", "KLS", "KLS",
                                             "MC", "MC", "MS", "MS",
                                             "MC", "MC", "MS", "MS",
                                             "MC", "MC", "MS", "MS",
                                             "MC", "MC", "MS",
                                             "ZIC", "ZIC", "ZIS",
                                             "ZIC", "ZIC", "ZIS", "ZIS", "ZIS",
                                             "ZIC", "ZIC", "ZIS", "ZIS", "ZIS",
                                             "ZIC", "ZIC", "ZIS")),
                      pop.line = as.factor(c("C2Del", "C2Del", "C2Del", "C2Del", 
                                             "C2Non", "C2Non", "C2Non", "C2Non",
                                             "KL1Del", "KL1Del",  "KL1Del", "KL1Del",
                                             "KL2Non",  "KL2Non","KL2Non", "KL2Non",
                                             "M12Del", "M12Del", "M12Del", "M12Del", 
                                             "M12Non","M12Non", "M12Non", "M12Non",
                                             "M9Del", "M9Del", "M9Del", "M9Del",
                                             "M9Non", "M9Non", "M9Non",
                                             "ZI197Non", "ZI197Non", "ZI197Non",
                                             "ZI254Non", "ZI254Non", "ZI254Non", "ZI254Non", "ZI254Non",
                                             "ZI273Non","ZI273Non", "ZI273Non", "ZI273Non", "ZI273Non",
                                             "ZI418Non", "ZI418Non", "ZI418Non")))

#Make table with sample information, must be in same order as column of count matrix!
samples


dds <- DESeqDataSetFromMatrix(countData = counts, 
                              colData = samples, 
                              design = ~treatment)

dds
#class: DESeqDataSet 
#dim: 13957 47 
#metadata(1): version
#assays(1): counts
#rownames(13957): FBgn0032191,CG5734,2L,2205 FBgn0053476,CG33476,2R,498 ...
#FBgn0053814,CG33814,2L,493 FBgn0034310,CG5733,2R,2146
#rowData names(0):
#  colnames(47): C2D3c_H3 C2D3c_H4 ... ZI418c_L3 ZI418s_L1
#colData names(5): pop.geno.treat population genotype treatment geno.treat

deseq1 <- DESeq(dds, betaPrior = TRUE)
#Shrink results: only need if betprior set to false in deseq function!!! as this does not shrink log2 values

resultsNames(deseq1)
# [1] "Intercept"  "treatmentc" "treatments"

as.data.frame(mcols(mcols(deseq1), use.names = TRUE))

# rdsfile <- paste(work_direct, "/dds3.rds", sep = "")
# rdsfile
# saveRDS(deseq1, file=rdsfile)


vsd3 <- vst(deseq1, blind = FALSE) ## variance stabilizing transformation, better because variance is the same across high and low read counts
pcaData3 <- plotPCA(vsd3, intgroup=c("line.treat", "population", "pop.line"), returnData = TRUE)
percentVar3 <- round(100 * attr(pcaData3, "percentVar"))

pcaData3$population <- factor(pcaData3$population, levels = c("M", "C2", "KL", "ZI"))
pcaData3$line.treat <- factor(pcaData3$line.treat, levels = c("MC", "C2C", "KLC","ZIC",  "MS", "C2S",  "KLS",  "ZIS"))

pca1 <- ggplot(pcaData3, aes(PC1, PC2 )) +
  geom_point(aes(fill = line.treat, 
                 shape = line.treat, 
                 color = line.treat),
             size = 1.7, 
             stroke = 1) +
  xlab(paste0("PC1: ",percentVar3[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar3[2],"% variance")) +
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
pca1


