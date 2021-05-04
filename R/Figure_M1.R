################ Main Figure 1 #######################
#Mortality MSB at 48 hours 
#checking first indel comparisons to get p-values for merged plot 

#tolerance plot of different lines at 48 hours
work_direct <- "C:/Users/Tim/Desktop/R_working_directory/MSB_stress_data"
work_direct
setwd(work_direct)

mun <- read_excel("MSB_Munich.xlsx")
cyp <- read_excel("MSB_Cyprus.xlsx")
kl <- read_excel("MSB_KL.xlsx")
zi <- read_excel("MSB_ZI.xlsx")


#plotting stress vials only 
#controls recorded (See separate excel file) but no death in these vials during experiment
#--> does not affect the mortality (%) in stress vials

############Mortality at 48 hours
#M9
s1 <- subset(mun, strain == "M9N1" | strain == "M9D3")
s1_48 <- subset(s1, Time =="48")
m9_plot <- ggplot(s1_48, aes(x = as.factor(Time), y = mortality, fill=genotype)) +
  geom_boxplot(width =0.5) +
  xlab("")+
  ylab("Mortality")+
  ggtitle("MUNICH 9, GERMANY")+
  theme(text = element_text(size=15))+
  stat_compare_means(label = "p.format", label.y = 1.05, size = 5, method = "wilcox.test" )+
  scale_fill_manual(values=c("#74ADD1", "#5AAE61"),
                    name="Genotype",
                    labels=c("Deletion", "Non-deletion"))+
  theme(panel.background = element_blank())+
  scale_y_continuous(limits = c(0, 1.1), breaks = c(0,0.2,0.4,0.6,0.8,1.00))+
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(face="bold", 
                                   size=15))


#M12
s2 <- subset(mun, strain == "M12D3" | strain == "M12N3")
s2_48 <- subset(s2, Time =="48")
m12_plot <- ggplot(s2_48, aes(x = as.factor(Time), y = mortality, fill=genotype)) +
  geom_boxplot(width =0.5) +
  xlab("") +
  ylab("")+
  ggtitle("MUNICH 12, GERMANY")+
  theme(text = element_text(size=15))+ 
  stat_compare_means(label = "p.format", label.y = 1.05, size = 5, method = "wilcox.test" )+
  scale_fill_manual(values=c("#74ADD1", "#5AAE61"),
                    name="Munich 12",
                    labels=c("Deletion", "Non-deletion"))+
  theme(panel.background = element_blank())+
  scale_y_continuous(limits = c(0, 1.1), breaks = c(0,0.2,0.4,0.6,0.8,1.00))+
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(face="bold", 
                                   size=15))


#KL

s3_48 <- subset(kl, Time =="48")
kl_plot <- ggplot(s3_48, aes(x = as.factor(strain), y = mortality, fill=genotype)) +
  geom_boxplot(width =0.5) +
  xlab("") +
  ylab("Mortality")+
  ggtitle("KUALA LUMPUR, MALAYSIA")+
  theme(text = element_text(size=15))+
  stat_compare_means(label = "p.format", label.y = 1.05,label.x = 1.5, size = 5, method = "wilcox.test" )+
  scale_fill_manual(values=c("#74ADD1", "#5AAE61"),
                    name="Genotype",
                    labels=c("Deletion", "Non-deletion"))+
  theme(panel.background = element_blank())+
  scale_y_continuous(limits = c(0, 1.1), breaks = c(0,0.2,0.4,0.6,0.8,1.00))+ 
  theme(axis.text.x = element_text(face="bold", 
                                   size=15, angle=45),
        axis.text.y = element_text(face="bold", 
                                   size=15))
#C2

s4_48 <- subset(cyp, Time =="48")
c2_plot <- ggplot(s4_48, aes(x = as.factor(Time), y = mortality, fill=genotype)) +
  geom_boxplot(width =0.5) +
  xlab("") +
  ylab("")+
  ggtitle("NICOSIA, CYPRUS")+
  theme(text = element_text(size=15))+ 
  stat_compare_means(label = "p.format", label.y = 1.05, size = 5, method = "wilcox.test" )+
  scale_fill_manual(values=c("#74ADD1", "#5AAE61"),
                    name="Genotype",
                    labels=c("Deletion", "Non-deletion"))+
  theme(panel.background = element_blank())+
  scale_y_continuous(limits = c(0, 1.1), breaks = c(0,0.2,0.4,0.6,0.8,1.00))+
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(face="bold", 
                                   size=15))


#ZI
s5_48 <- subset(zi, Time =="48")
zi_plot <- ggplot(s5_48, aes(factor(strain), mortality)) +
  geom_boxplot(aes(x = as.factor(strain), y = mortality, fill=genotype)) +
  xlab("") +
  ylab("Mortality")+
  ggtitle("SIAVONGA, ZAMBIA")+
  theme(text = element_text(size=15))+
  scale_fill_manual(name = "",values = c("#5AAE61"), 
                    labels =c("Non-deletion"))+
  theme(panel.background = element_blank())+
  scale_y_continuous(limits = c(0, 1.6), breaks = c(0,0.2,0.4,0.6,0.8,1.00))+ 
  theme(axis.text.x = element_text(face="bold", 
                                   size=15, angle=45),
        axis.text.y = element_text(face="bold", 
                                   size=15))


##########combining all lines into on plot####

#computed the p-values between indel lines for each group as above, 
#and manually added the p value symbols with annotate when 
#re-plotting the separate plots into one figure (below)

#"#7570B3", Germany
#"#1B9E77", Cyprus
#"#66A61E", Malaysia
#"#D95F02", Zambia

merge1 <- rbind(s1_48, s2_48, s3_48, s4_48, s5_48)
merge1
merge1$strain <- factor(merge1$strain,levels = c("M12N3", "M12D3", "M9N1", "M9D3", "C2N1", "C2D3",
                                                 "KL2", "KL1", "ZI197", "ZI254", "ZI273", "ZI418"))
merge1$population <- factor(merge1$population,levels = c("Germany", "Cyprus", "Malaysia", "Zambia"))

merge1_plot <- ggplot(merge1, aes(factor(strain), mortality)) +
  geom_boxplot(aes(x = as.factor(strain), y = mortality, fill = population), width =0.5, outlier.color = "white") +
  xlab("Fly Line") +
  ylab("Mortality")+
  theme(text = element_text(size=10))+
  scale_fill_manual(values = c("#7570B3", "#1B9E77", "#66A61E", "#D95F02")) +
  theme(axis.text.x = element_text(size=10,angle = 35, hjust=1))+
  scale_y_continuous(limits = c(0, 1.1), breaks = c(0,0.2,0.4,0.6,0.8,1.00))+
  theme(axis.text.y = element_text(size=10))+
  theme(panel.background = element_blank())+ 
  scale_x_discrete(labels = c(expression("M12", "M12"*Delta*"", "M9", "M9"*Delta*"", "C2", "C2"*Delta*"", "KL2", "KL1"*Delta*"", "ZI197", "ZI254", "ZI273", "ZI418")))+
  theme(axis.line = element_line(color="black"))+
  theme(legend.position = "none")
merge1_plot

#adding p-values
merge2 <- merge1_plot + annotate("text", x = 1.5, y = 1.06, label = "***", size = 7)+
  annotate("text", x = 3.5, y = 1.06, label = "**", size = 7)+
  annotate("text", x = 7.5, y = 1.06, label = "***", size = 7)+
  annotate("text", x = 1.5, y = 1.08, label = "___", size = 5)+
  annotate("text", x = 3.5, y = 1.08, label = "___", size = 5)+
  annotate("text", x = 7.5, y = 1.08, label = "___", size = 5)
merge2 

ggsave(file = "C:/Users/Tim/Desktop/2020Publication/Figures_2020pub/msb_mortality.svg", 
       plot = merge2, device = "svg",
       scale = 1,
       dpi = 200,
       limitsize = TRUE)


ggsave(file = "C:/Users/Tim/Desktop/2020Publication/Figures_2020pub/msb_mortality.pdf", 
       plot = merge2, device = "pdf",
       scale = 1,
       dpi = 200,
       limitsize = TRUE) 

ggsave(file = "C:/Users/Tim/Desktop/2020Publication/Figures_2020pub/msb_mortality.eps", 
       plot = merge2, device = "eps",
       scale = 1,
       dpi = 200,
       limitsize = TRUE)

