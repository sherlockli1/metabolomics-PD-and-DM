library(readxl)
library(VennDiagram)
library(ggplot2)
library(dplyr)
library(data.table)
library(ggpubr)
library(gridExtra)
library("ggpmisc")
#Load mummichog output
path_hilic_pd<-read_excel("pathway_result/tables.xlsx", sheet = 1)
path_hilic_control<-read_excel("pathway_result/tables.xlsx", sheet = 2)
path_c18_pd<-read_excel("pathway_result/tables.xlsx", sheet = 3)
path_c18_control<-read_excel("pathway_result/tables.xlsx", sheet = 4)

#combine hilic and c18 and remove duplicated pathways
pd_all<-rbind(path_hilic_pd,path_c18_pd)
pd_unique<-pd_all[!duplicated(pd_all$Pathways),]

control_all<-rbind(path_hilic_control,path_c18_control)
control_unique<-control_all[!duplicated(control_all$Pathways),]

#combine pd and control
all_pd<-rbind(pd_unique,control_unique)

#count number of unique and shared pathways
counts<-as.data.frame(table(all_pd$Pathways, all_pd$PD))
counts_PD<-counts[which(counts$Var2=="PD"),]
colnames(counts_PD)[3]<-"freq_pd"
counts_control<-counts[which(counts$Var2=="non-PD"),]
colnames(counts_control)[3]<-"freq_control"
counts_final<-merge(counts_PD,counts_control,by="Var1")
counts_final$unique_pd<-ifelse(counts_final$freq_pd==1&counts_final$freq_control==0,1,0)
counts_final$unique_control<-ifelse(counts_final$freq_pd==0&counts_final$freq_control==1,1,0)
counts_final$shared<-ifelse(counts_final$freq_pd==1&counts_final$freq_control==1,1,0)

table(counts_final$unique_pd)
table(counts_final$unique_control)
table(counts_final$shared)

#create Venn diagram
# venn.plot <- draw.pairwise.venn(44, 20, 16, 
#                    category = c("PD Patient", "Control"), 
#                    lty = rep("blank",2), 
#                    fill = c("light blue", "pink"), 
#                    alpha = rep(0.5, 2), 
#                    cat.pos = c(0,0), 
#                    cat.dist = rep(0.025, 2))
# grid.newpage()
# p1 <- grobTree(venn.plot)
#

df <- data.frame(row.names=c("Shared Pathways", 
                             "Pathways only in PD patients", 
                             "Pathways only in non-PD participants"),
                 freq=c(16,28,4))
df <- as.data.frame(t(df))
tt <- ttheme_default(colhead=list(fg_params = list(parse=TRUE)),
                     base_size = 8,
                     padding = unit(c(2, 4), "mm"))
tbl <- tableGrob(df, rows=NULL, theme=tt)


# all<-rbind(path_hilic_pd,path_hilic_control,path_c18_pd, path_c18_control)
# write.csv(all,"Manuscript/figure/all_path_raw.csv",row.names=FALSE)

#Make Plots
table<-read.csv("Manuscript/figure/all_path.csv")
table$logpvalue<--log(table$p.value)

table$Column = factor(table$Column, levels=c('HILIC','C18'))
table$PD = factor(table$PD, levels=c('PD Patient','non-PD'))
table(table$Group)

table$Group <- factor(table$Group, levels = c("Simple Sugar",
                                              "Complex Sugar",
                                              "Glycosamino glycan",
                                              "N-Glycan",
                                              "Aminosugars",
                                              "Fatty Acid",
                                              "Glycerophospho lipid",
                                              "Glycosphingo lipid",
                                              "Amino Acid",
                                              "Nucleic Acid",
                                              "Vitamin",
                                              "Other"))
#Create figure
#plot by molecule type
p2<-ggplot(table, aes(x = reorder(Pathways, +Order), y = logpvalue, fill=Group)) + 
  geom_col(show.legend = FALSE)+
  facet_grid(Column+PD ~ Group, scales = "free", space='free',
             labeller = labeller(Group = label_wrap_gen(4)))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        axis.text.y = element_text(size = 5),
        plot.margin = margin(1, 1, 1, 2, unit = "cm"),
        strip.background.x =element_rect(fill=NA, colour = 'black', size = 0.5),
        strip.background.y =element_rect(fill=NA, colour = 'black', size = 0.5),
        panel.spacing = unit(.05, "lines"),
        strip.text.x = element_text(colour = 'black',size=5, angle = 90),
        strip.text.y = element_text(angle = 0, colour = 'black', size=5),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        legend.position="top")+
  labs(x ="Pathways", y = "-Log(p value)")
p2

#Plot by functional group
table_1<-read.csv("Manuscript/figure/all_path_function.csv")
colnames(table_1)
new_table<-table_1 %>%
  group_by(PD) %>%
  summarize(sum_Energy.Metabolism = sum(Energy.Metabolism== 1),
            sum_DNA.RNA = sum(DNA.RNA== 1),
            sum_Nitrogen.balance = sum(Nitrogen.balance== 1),
            sum_Connective.Tissue = sum(Connective.Tissue== 1),
            sum_Antioxidant = sum(Antioxidant== 1),
            sum_Immune.System = sum(Immune.System== 1),
            sum_Signaling.Pathway = sum(Signaling.Pathway== 1),
            sum_Glycolysis = sum(Glycolysis== 1),
            sum_Urea.Cycle = sum(Urea.Cycle== 1),
            sum_TCA.Cycle = sum(TCA.Cycle== 1)
            )

new_table<-melt(setDT(new_table),id.vars = c("PD"), variable.name = "Pathway")

new_table$Pathway_new<-
  c("Energy Metabolism",
    "Energy Metabolism",
    "DNA/RNA Damage and Repair",
    "DNA/RNA Damage and Repair",
    "Nitrogen Balance",
    "Nitrogen Balance",
    "Connective Tissue",
    "Connective Tissue",
    "Antioxidant/Oxidative Stress",
    "Antioxidant/Oxidative Stress",
    "Immune System",
    "Immune System",
    "Signaling Pathway",
    "Signaling Pathway",
    "Glycolysis",
    "Glycolysis",
    "Urea Cycle",
    "Urea Cycle",
    "Tricarboxylic Acid Cycle",
    "Tricarboxylic Acid Cycle")
new_table$PD = factor(new_table$PD, levels=c('PD Patient','non-PD'))

new_table$value<-ifelse(new_table$PD=="PD Patient",new_table$value/45,new_table$value/20)
new_table$value<-round(new_table$value,2)

# Create a bar chart using ggplot2
p3<-ggplot(data=new_table, aes(x=Pathway_new, y=value, fill=PD)) +
  geom_bar(stat="identity", position=position_dodge())+
  geom_text(aes(label=value), vjust=-0.3, color="black",
            position = position_dodge(0.9), size=1)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10), 
        plot.margin = margin(10, 10, 10, 100),
        strip.background.x =element_rect(fill=NA, colour = "black", size = 0.5),
        strip.background.y =element_rect(fill=NA, colour = "black", size = 0.5),
        panel.spacing = unit(.05, "lines"),
        strip.text.x = element_text(colour = 'black', size=5,face="bold"),
        strip.text.y = element_text(angle = 0, colour = 'black', size=8,face="bold"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.position="top")+
  labs(x ="Pathways' Functional Group", y = "Proportion of Pathways")
p3

# combined<-ggarrange(ggarrange(p2, p3,
#                     ncol = 2, nrow = 1,
#                     labels=c("A","B"),
#                     widths = c(2.5, 1)),
#                     tbl,
#                     nrow=2,
#                     heights=c(2,0.3),
#                     labels=c("","C"))
# combined


png(filename="Manuscript/figure/figure3A.png",
    units="in", 
    width=16, 
    height=6, 
    res=400)
print(p2)
dev.off()

png(filename="Manuscript/figure/figure3B.png",
    units="in", 
    width=16, 
    height=6, 
    res=400)
print(p3)
dev.off()











