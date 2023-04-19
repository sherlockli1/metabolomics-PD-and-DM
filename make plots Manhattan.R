#Create plots
library(stringr)
library(tibble)
library(mixOmics)
library(qqman)
library(ggplot2)
library(ggpubr)
library(dplyr)

#################################################################################
#MWAS Plots
load("mwasResults/hilic_medcombat_diabetes_MWAS_new_sd.RData")
load("mwasResults/c18_medcombat_diabetes_MWAS_new_sd.RData")

hilic_medcombat_diabetes_MWAS$direction<-ifelse(hilic_medcombat_diabetes_MWAS$beta>0&hilic_medcombat_diabetes_MWAS$FDR_pvalue<0.05,"positive",
                                                ifelse(hilic_medcombat_diabetes_MWAS$beta<0&hilic_medcombat_diabetes_MWAS$FDR_pvalue<0.05,"negative","FDR-adjusted p-value<0.05"))
c18_medcombat_diabetes_MWAS$direction<-ifelse(c18_medcombat_diabetes_MWAS$beta>0&c18_medcombat_diabetes_MWAS$FDR_pvalue<0.05,"positive",
                                              ifelse(c18_medcombat_diabetes_MWAS$beta<0&c18_medcombat_diabetes_MWAS$FDR_pvalue<0.05,"negative","FDR-adjusted p-value<0.05"))
hilic_medcombat_diabetes_MWAS$direction <- factor(hilic_medcombat_diabetes_MWAS$direction, levels = c("positive", "negative", "FDR-adjusted p-value<0.05"))
c18_medcombat_diabetes_MWAS$direction <- factor(c18_medcombat_diabetes_MWAS$direction, levels = c("positive", "negative", "FDR-adjusted p-value<0.05"))


hilic_medcombat_diabetes_MWAS[c('drop1', 'drop2','mz','rt')] <- str_split_fixed(hilic_medcombat_diabetes_MWAS$metabo, '_', 4)
c18_medcombat_diabetes_MWAS[c('drop1', 'drop2','mz','rt')] <- str_split_fixed(c18_medcombat_diabetes_MWAS$metabo, '_', 4)

hilic_medcombat_diabetes_MWAS$mz<-as.numeric(hilic_medcombat_diabetes_MWAS$mz)
c18_medcombat_diabetes_MWAS$mz<-as.numeric(c18_medcombat_diabetes_MWAS$mz)

hilic_medcombat_diabetes_MWAS$newp<--log(hilic_medcombat_diabetes_MWAS$FDR_pvalue)
c18_medcombat_diabetes_MWAS$newp<--log(c18_medcombat_diabetes_MWAS$FDR_pvalue)

p1<-ggplot(hilic_medcombat_diabetes_MWAS, aes(x=mz, y=newp))+ 
  geom_point(aes(colour =direction),show.legend = TRUE)+
  geom_hline(yintercept=3, linetype="dashed", color = "#D55E00")+ 
  scale_color_manual(values=c("#E69F00", "#0072B2", "#999999"))+
  labs(x ="m/z ratio", y="-log(FDR p value)", color = "Direction of Association")+ 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
p1

p2<-ggplot(c18_medcombat_diabetes_MWAS, aes(x=mz, y=newp))+ 
  geom_point(aes(colour =direction),show.legend = TRUE)+
  geom_hline(yintercept=3, linetype="dashed", color = "#D55E00")+ 
  scale_color_manual(values=c("#E69F00", "#0072B2", "#999999"))+
  labs(x ="m/z ratio", y="-log(FDR p value)", color = "Direction of Association")+ 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
p2

load("mwasResults/control_hilic_medcombat_diabetes_MWAS_new_sd.RData")
load("mwasResults/control_c18_medcombat_diabetes_MWAS_new_sd.RData")

hilic_medcombat_diabetes_MWAS$direction<-ifelse(hilic_medcombat_diabetes_MWAS$beta>0&hilic_medcombat_diabetes_MWAS$FDR_pvalue<0.05,"positive",
                                                ifelse(hilic_medcombat_diabetes_MWAS$beta<0&hilic_medcombat_diabetes_MWAS$FDR_pvalue<0.05,"negative","FDR-adjusted p-value<0.05"))
c18_medcombat_diabetes_MWAS$direction<-ifelse(c18_medcombat_diabetes_MWAS$beta>0&c18_medcombat_diabetes_MWAS$FDR_pvalue<0.05,"positive",
                                              ifelse(c18_medcombat_diabetes_MWAS$beta<0&c18_medcombat_diabetes_MWAS$FDR_pvalue<0.05,"negative","FDR-adjusted p-value<0.05"))
hilic_medcombat_diabetes_MWAS$direction <- factor(hilic_medcombat_diabetes_MWAS$direction, levels = c("positive", "negative", "FDR-adjusted p-value<0.05"))
c18_medcombat_diabetes_MWAS$direction <- factor(c18_medcombat_diabetes_MWAS$direction, levels = c("positive", "negative", "FDR-adjusted p-value<0.05"))


hilic_medcombat_diabetes_MWAS[c('drop1', 'drop2','mz','rt')] <- str_split_fixed(hilic_medcombat_diabetes_MWAS$metabo, '_', 4)
c18_medcombat_diabetes_MWAS[c('drop1', 'drop2','mz','rt')] <- str_split_fixed(c18_medcombat_diabetes_MWAS$metabo, '_', 4)

hilic_medcombat_diabetes_MWAS$mz<-as.numeric(hilic_medcombat_diabetes_MWAS$mz)
c18_medcombat_diabetes_MWAS$mz<-as.numeric(c18_medcombat_diabetes_MWAS$mz)

hilic_medcombat_diabetes_MWAS$newp<--log(hilic_medcombat_diabetes_MWAS$FDR_pvalue)
c18_medcombat_diabetes_MWAS$newp<--log(c18_medcombat_diabetes_MWAS$FDR_pvalue)

p3<-ggplot(hilic_medcombat_diabetes_MWAS, aes(x=mz, y=newp))+ 
  geom_point(aes(colour =direction),show.legend = TRUE)+
  geom_hline(yintercept=3, linetype="dashed", color = "#D55E00")+ 
  scale_color_manual(values=c("#E69F00", "#0072B2", "#999999"))+
  labs(x ="m/z ratio", y="-log(FDR p value)", color = "Direction of Association")+ 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
p3

p4<-ggplot(c18_medcombat_diabetes_MWAS, aes(x=mz, y=newp))+ 
  geom_point(aes(colour =direction),show.legend = TRUE)+
  geom_hline(yintercept=3, linetype="dashed", color = "#D55E00")+ 
  scale_color_manual(values=c("#E69F00", "#0072B2", "#999999"))+
  labs(x ="m/z ratio", y="-log(FDR p value)", color = "Direction of Association")+ 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
p4
#Combine plots
p5<-ggarrange(p1,p2,p3,p4, 
              labels = c("A", "B","C", "D"),
              heights = c(3, 3),
              ncol = 2, nrow = 2, common.legend = TRUE, legend="top")
p5
png(filename="Manuscript/figure/figure1.png",
    units="in", 
    width=12, 
    height=6, 
    res=300)
print(p5)
dev.off()



