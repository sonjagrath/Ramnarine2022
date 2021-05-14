################################################################################
# Figure_M3_and_SupFigures_2-5.R
################################################################################

## 
## Authors: Timothy Ramnarine
## Input:
### 
### 
## Last update: 2021-05-14

################################################################################


library(readxl)
library(ggplot2)


#Supplemental Figure 2 - comparison of Genotype-environment interaction genes (von Heckel et al 2016)
#von Heckel, Korbinian, Wolfgang Stephan, and Stephan Hutter. 
#"Canalization of gene expression is a major signature of regulatory cold adaptation in temperate Drosophila melanogaster." 
#BMC genomics 17, no. 1 (2016): 1-14.
# MSB vs rec90 log2 fold change for Genotype-environment interaction candidate genes

#using simplified p-value to indicate y/n significance! - only reflects if the p-value is less than 0.05


gei <- read_excel("../data/L2FC_GEI_rec90_v_MSB.xlsx")
geiplot <-  ggplot(gei, aes(x = Name , y = L2FC)) +
  geom_linerange(
    aes(x = Name, ymin = 0, ymax = L2FC, group = Dataset), 
    color = "lightgray", size = 1.5,
    position = position_dodge(0.55)
  )+
  geom_point(
    aes(color = Dataset),
    position = position_dodge(0.55), size = 3
  )+
  scale_y_continuous(limits = c(-3, 5),breaks=c(-3,-2,-1,0,1,2,3,4,5))+
  geom_text(label = gei$psimple_msb, y = 4.7, angle = 0,nudge_x=-0.15,  size = 4)+
  geom_text(label = gei$psimple_rec90, y = 4.7, angle = 0,nudge_x =0.15, size = 4)+
  scale_color_manual(values = c("#EB9622", "#4F17C6"))+
  labs(x = "", y = "Log2(stress/control)")+
  theme(text = element_text(size=10))+
  theme(axis.text.x = element_text(size=10,angle = 0, hjust=0.5))+
  theme(axis.text.y = element_text(size=10))+
  theme(panel.background = element_blank())+ 
  theme(legend.position="bottom")+
  theme(legend.title = element_blank())+
  theme(axis.line.y = element_line(color="black"))+
  geom_hline(yintercept = 0)
geiplot1 <- geiplot + coord_flip() 
geiplot1
########################################################

#Supplemental Figure 3 - comparison of cold tolerance candidate gene list (von Heckel et al 2016)
#von Heckel, Korbinian, Wolfgang Stephan, and Stephan Hutter. 
#"Canalization of gene expression is a major signature of regulatory cold adaptation in temperate Drosophila melanogaster." 
#BMC genomics 17, no. 1 (2016): 1-14.
#MSB vs rec90 log2 fold change for cold tolerance candidate gene list 

#using simplified p-value to indicate y/n significance! - only reflects if the p-value is less than 0.05

l2fc_rec90MSB <- read_excel("../data/L2FC_rec90_v_MSB.xlsx")

#subsetting for easier plotting 
subv <- subset(l2fc_rec90MSB, Viability == "viable") #viable genes
subnv <- subset(l2fc_rec90MSB, Viability == "not viable" | Viability == "malformed") #not viable + malformed

p1 <-  ggplot(subv, aes(x = Gene_name , y = L2FC)) +
  geom_linerange(
    aes(x = Gene_name, ymin = 0, ymax = L2FC, group = Dataset), 
    color = "lightgray", size = 1.5,
    position = position_dodge(0.55)
  )+
  geom_point(
    aes(color = Dataset),
    position = position_dodge(0.55), size = 3
  )+
  scale_y_continuous(limits = c(-2, 8.2),breaks=c(-2,-1,0,1,2,3,4,5,6,7,8))+
  geom_text(label = subv$psimple_msb, y = 8, angle = 0,nudge_x=-0.15,  size = 4)+
  geom_text(label = subv$psimple_rec90, y = 8, angle = 0,nudge_x =0.15, size = 4)+
  scale_color_manual(values = c("#EB9622", "#4F17C6"))+
  ggpubr::theme_pubclean()+
  labs(x = "", y = "Log2(stress/control)")+
  theme(text = element_text(size=10))+
  theme(axis.text.x = element_text(size=10,angle = 0, hjust=1))+
  theme(axis.text.y = element_text(size=10))+
  theme(panel.background = element_blank())+ 
  guides(color = FALSE, size = FALSE, fill = FALSE)+
  theme(axis.line.y = element_line(color="black"))+
  geom_hline(yintercept = 0)

p1a <- p1 + coord_flip()
p1a


p2<- ggplot(subnv, aes(x = Gene_name, y = L2FC)) +
  geom_linerange(
    aes(x = Gene_name, ymin = 0, ymax = L2FC, group = Dataset), 
    color = "lightgray", size = 1.5,
    position = position_dodge(0.55)
  )+
  geom_point(
    aes(color = Dataset),
    position = position_dodge(0.55), size = 3
  )+
  scale_y_continuous(limits = c(-2, 8.2),breaks=c(-2,-1,0,1,2,3,4,5,6,7,8))+
  geom_text(label = subnv$psimple_msb, y = 8, angle = 0,nudge_x=-0.15,  size = 4)+
  geom_text(label = subnv$psimple_rec90, y = 8, angle = 0,nudge_x =0.15, size = 4)+
  scale_color_manual(values = c("#EB9622", "#4F17C6"))+
  ggpubr::theme_pubclean()+
  labs(x = "", y = "Log2(stress/control)")+
  theme(text = element_text(size=10))+
  theme(axis.text.x = element_text(size=10,angle = 0, hjust=1))+
  theme(axis.text.y = element_text(size=10))+
  theme(panel.background = element_blank())+ 
  guides(color = FALSE, size = FALSE, fill = FALSE)+
  theme(axis.line.y = element_line(color="black"))+
  geom_hline(yintercept = 0)

p2a <- p2 + coord_flip()
p2a

#for legend#
plotall <-  ggplot(l2fc_rec90MSB ,aes(factor(Gene_name), L2FC)) + 
  geom_bar(aes(x = as.factor(Gene_name), y = L2FC, fill = as.factor(Dataset)), stat = "identity", width = 0.9) +
  theme(panel.background = element_blank())+ 
  scale_fill_manual(values = c ("#EB9622", "#4F17C6"), labels = c("MSB", "rec90"), name = "")+
  labs(x = "Gene Name", y = "Log2(stress/control)")+
  theme(legend.position="bottom", legend.direction="horizontal")

legend <- ggpubr::get_legend(plotall)

#arrange plots as single figure with legend
plot1_2 <- ggpubr::ggarrange(p1a,p2a, nrow =1, ncol = 2, legend.grob = legend, legend = "bottom", labels = c("A", "B"))

####################################################


#Main Figure 4 & Supplemental Figure 4 and 5

cand <- read_excel("../data/GeneGroups.xlsx")

mtn <- subset(cand, Group == "1") #Metallothionein group - Main Figure 3

web <- subset(cand, Group == "3") #Weber candidates - Supp Figure 4

bari_jheh <- subset(cand, Group == "2") #Bari-Jheh Epoxide hydrolase - Supp Figure 5
########################################################################################################

#Main Figure 3
#using actual p-values as *P < 0.05, **P < 0.01, ***P < 0.001

mtnplot <- ggplot(mtn ,aes(factor(Name), MSB_LFC)) + 
  geom_bar(aes(x = as.factor(Name), y = MSB_LFC, fill = as.factor(Group)), stat = "identity", width = 0.5) +
  theme(panel.background = element_blank())+ 
  scale_y_continuous(limits = c(0, 3.7),breaks=c(0,0.5,1,1.5,2,2.5,3,3.5))+
  geom_text(label = mtn$psig_msb, y = 3.5, angle = 0, size =7)+
  scale_fill_manual(values = c ("#3c764f"), labels = c("Metallothionein Group"))+
  labs(x = "", y = "Log2(stress/control)")+
  theme(text = element_text(size=10))+
  theme(axis.text.x = element_text(size=10,angle = 45, hjust=0.5))+
  theme(axis.text.y = element_text(size=10))+
  theme(panel.background = element_blank())+ 
  guides(color = FALSE, size = FALSE, fill = FALSE)+
  theme(axis.line.y = element_line(color="black"))+
  geom_hline(yintercept = 0)
mtnplot


#Supplemental Figure 4 - Weber et al. 2012
#Weber, Allison L., George F. Khan, Michael M. Magwire, Crystal L. Tabor, Trudy FC Mackay, and Robert RH Anholt. 
#"Genome-wide association analysis of oxidative stress resistance in Drosophila melanogaster." 
#PloS one 7, no. 4 (2012): e34745.
#CG9650, Ecodysone-induced protein 75 B (Eip75B), ena, fog, homeobrain (hbn), nACR??-30D, and rg) 
#associated with oxidative stress
#using simplified p-value to indicate y/n significance! - only reflects if the p-value is less than 0.05

webplot <- ggplot(web ,aes(factor(Name), MSB_LFC)) + 
  geom_bar(aes(x = as.factor(Name), y = MSB_LFC, fill = as.factor(Group)), stat = "identity", width = 0.5) +
  theme(panel.background = element_blank())+ 
  scale_y_continuous(limits = c(-0.2, 0.4),breaks=c(-0.2 ,0, 0.2,0.4))+
  geom_text(label = web$psimple_msb, y = 0.3, angle = 0, size = 5)+
  scale_fill_manual(values = c ("#42ab83"), labels = c("Weber 2012 candidates"))+
  labs(x = "", y = "Log2(stress/control)")+
  theme(text = element_text(size=10))+
  theme(axis.text.x = element_text(size=10,angle = 45, hjust=0.5))+
  theme(axis.text.y = element_text(size=10))+
  theme(panel.background = element_blank())+ 
  guides(color = FALSE, size = FALSE, fill = FALSE)+
  theme(axis.line.y = element_line(color="black"))+
  geom_hline(yintercept = 0)
webplot


#Supplemental Figure 5 - Bari Jheh
#using simplified p-value to indicate y/n significance! - only reflects if the p-value is less than 0.05

bariplot <- ggplot(bari_jheh ,aes(factor(Name), MSB_LFC)) + 
  geom_bar(aes(x = as.factor(Name), y = MSB_LFC, fill = as.factor(Group)), stat = "identity", width = 0.5) +
  theme(panel.background = element_blank())+ 
  scale_y_continuous(limits = c(-1, 2),breaks=c(-1, -0.5,0,0.5,1,1.5,2))+
  geom_text(label = bari_jheh$psimple_msb, y = 1.7, angle = 0, size = 5)+
  scale_fill_manual(values = c ("#549422"), labels = c("Epoxide Hydrolase"))+
  labs(x = "", y = "Log2(stress/control)")+
  theme(text = element_text(size=10))+
  theme(axis.text.x = element_text(size=10,angle = 45, hjust=0.5))+
  theme(axis.text.y = element_text(size=10))+
  theme(panel.background = element_blank())+ 
  guides(color = FALSE, size = FALSE, fill = FALSE)+
  theme(axis.line.y = element_line(color="black"))+
  geom_hline(yintercept = 0)
bariplot
