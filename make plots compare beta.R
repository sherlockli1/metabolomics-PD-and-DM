library(ggplot2)
library(tibble)
library(ggpubr)
library(wesanderson)
library(gridExtra)
library(gtsummary)
library("ggpmisc")

#Load MWAS results
load("mwasResults/hilic_medcombat_diabetes_MWAS_new_sd.Rdata")
pd_hilic<-hilic_medcombat_diabetes_MWAS
load("mwasResults/c18_medcombat_diabetes_MWAS_new_sd.Rdata")
pd_c18<-c18_medcombat_diabetes_MWAS

load("mwasResults/control_hilic_medcombat_diabetes_MWAS_new_sd.Rdata")
control_hilic<-hilic_medcombat_diabetes_MWAS
load("mwasResults/control_c18_medcombat_diabetes_MWAS_new_sd.Rdata")
control_c18<-c18_medcombat_diabetes_MWAS

#Add indicator for significance
#HILIC
pd_hilic$indicator_crudep<-ifelse(pd_hilic$pvalue<0.05&control_hilic$pvalue<0.05,"Significant in both",
                             ifelse(pd_hilic$pvalue<0.05&control_hilic$pvalue>=0.05,"Significant only in PD patients",
                                    ifelse(pd_hilic$pvalue>=0.05&control_hilic$pvalue<0.05,"Significant only in non-PD participants","Not significant in either group")))
pd_hilic$indicator_crudep <- factor(pd_hilic$indicator_crudep, 
                                    levels = c("Significant in both", 
                                               "Significant only in PD patients", 
                                               "Significant only in non-PD participants",
                                               "Not significant in either group"))

pd_hilic$indicator_fdrp<-ifelse(pd_hilic$FDR_pvalue<0.05&control_hilic$FDR_pvalue<0.05,"Significant in both",
                                  ifelse(pd_hilic$FDR_pvalue<0.05&control_hilic$FDR_pvalue>=0.05,"Significant only in PD patients",
                                         ifelse(pd_hilic$FDR_pvalue>=0.05&control_hilic$FDR_pvalue<0.05,"Significant only in non-PD participants","Not significant in either group")))
pd_hilic$indicator_fdrp <- factor(pd_hilic$indicator_fdrp, 
                                    levels = c("Significant in both", 
                                               "Significant only in PD patients", 
                                               "Significant only in non-PD participants",
                                               "Not significant in either group"))

#C18
pd_c18$indicator_crudep<-ifelse(pd_c18$pvalue<0.05&control_c18$pvalue<0.05,"Significant in both",
                                  ifelse(pd_c18$pvalue<0.05&control_c18$pvalue>=0.05,"Significant only in PD patients",
                                         ifelse(pd_c18$pvalue>=0.05&control_c18$pvalue<0.05,"Significant only in non-PD participants","Not significant in either group")))
pd_c18$indicator_crudep <- factor(pd_c18$indicator_crudep, 
                                    levels = c("Significant in both", 
                                               "Significant only in PD patients", 
                                               "Significant only in non-PD participants",
                                               "Not significant in either group"))


pd_c18$indicator_fdrp<-ifelse(pd_c18$FDR_pvalue<0.05&control_c18$FDR_pvalue<0.05,"Significant in both",
                                ifelse(pd_c18$FDR_pvalue<0.05&control_c18$FDR_pvalue>=0.05,"Significant only in PD patients",
                                       ifelse(pd_c18$FDR_pvalue>=0.05&control_c18$FDR_pvalue<0.05,"Significant only in non-PD participants","Not significant in either group")))
pd_c18$indicator_fdrp <- factor(pd_c18$indicator_fdrp, 
                                  levels = c("Significant in both", 
                                             "Significant only in PD patients", 
                                             "Significant only in non-PD participants",
                                             "Not significant in either group"))


#Make Plots
#Colored by crude p value, beta from HILIC PD and Control
plot <- ggplot(pd_hilic, aes(x = beta, y = control_hilic$beta)) +
  geom_point(aes(color = indicator_crudep)) +
  scale_colour_manual(values=c("#E69F00","#56B4E9","#009E73","#999999")) +
  geom_smooth(method = "lm", color = "black") +
  stat_cor(method = "pearson", size = 3)+
  labs(title="Correlation between metabolic feature coefficients among \n PD patients and non-PD participants in HILIC column",
       x ="Metabolic feature coefficients among PD patients", 
       y = "Metabolic feature coefficients among non-PD participants",
       color= "Crude P value significance level")+ 
  geom_vline(xintercept = 0)+
  geom_hline(yintercept = 0)+
  theme(
         plot.title = element_text(color="Black", size=10, face="bold.italic"),
         axis.title.x = element_text(color="Black", size=8, face="bold"),
         axis.title.y = element_text(color="Black", size=8, face="bold"),
         legend.title = element_text(color="Black", size=8, face="bold")
       )+ 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                panel.background = element_blank(), axis.line = element_line(colour = "black"))
plot

#Make descriptive statistics
table(pd_hilic$indicator_crudep)
df <- data.frame(row.names=c("Significant in both", 
                        "Significant only in PD patients", 
                        "Significant only in non-PD participants",
                        "Not significant in either group"),
                 freq=c(139,409,247,2118))
df <- as.data.frame(t(df))
tt <- ttheme_default(colhead=list(fg_params = list(parse=TRUE)),
                     base_size = 6,
                     padding = unit(c(2, 4), "mm"))
tbl <- tableGrob(df, rows=NULL, theme=tt)
combined<-ggarrange(plot, tbl,
          ncol = 1, nrow = 2,
          heights = c(1, 0.5))
combined



#Colored by crude p value, beta from C18 PD and Control
plot <- ggplot(pd_c18, aes(x = beta, y = control_c18$beta)) +
  geom_point(aes(color = indicator_crudep)) +
  scale_colour_manual(values=c("#E69F00","#56B4E9","#009E73","#999999")) +
  geom_smooth(method = "lm", color = "black") +
  stat_cor(method = "pearson", size = 3)+
  labs(title="Correlation between metabolic feature coefficients among \n PD patients and non-PD participants in C18 column",
       x ="Metabolic feature coefficients among PD patients", 
       y = "Metabolic feature coefficients among non-PD participants",
       color= "Crude P value significance level")+
  geom_vline(xintercept = 0)+
  geom_hline(yintercept = 0)+
  theme(
    plot.title = element_text(color="Black", size=10, face="bold.italic"),
    axis.title.x = element_text(color="Black", size=8, face="bold"),
    axis.title.y = element_text(color="Black", size=8, face="bold"),
    legend.title = element_text(color="Black", size=8, face="bold")
  )+ 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
plot

#Make descriptive statistics
table(pd_c18$indicator_crudep)
df <- data.frame(row.names=c("Significant in both", 
                             "Significant only in PD patients", 
                             "Significant only in non-PD participants",
                             "Not significant in either group"),
                 freq=c(88,323,130,1681))
df <- as.data.frame(t(df))
tt <- ttheme_default(colhead=list(fg_params = list(parse=TRUE)),
                     base_size = 6,
                     padding = unit(c(2, 4), "mm"))
tbl <- tableGrob(df, rows=NULL, theme=tt)
combined2<-ggarrange(plot, tbl,
                    ncol = 1, nrow = 2,
                    heights = c(1, 0.5))
combined2

combined_final<-ggarrange(combined, combined2,
                     ncol = 2, nrow = 1,
                     labels = c("A", "B"),
                     heights = c(1, 1))
combined_final



png(filename="Manuscript/figure/cor_crude_p.png",
    units="in", 
    width=12, 
    height=6, 
    res=300)
print(combined_final)
dev.off()



################################################################################
#Colored by FDR p value, beta from HILIC PD and Control
plot <- ggplot(pd_hilic, aes(x = beta, y = control_hilic$beta)) +
  geom_point(aes(color = indicator_fdrp)) +
  scale_colour_manual(values=c("#E69F00","#56B4E9","#009E73","#999999")) +
  geom_smooth(method = "lm", color = "black") +
  stat_cor(method = "pearson", size = 3)+
  labs(title="Correlation between metabolic feature coefficients among \n PD patients and non-PD participants in HILIC column",
       x ="Metabolic feature coefficients among PD patients", 
       y = "Metabolic feature coefficients among non-PD participants",
       color= "FDR-adjusted P value significance level")+
  geom_vline(xintercept = 0)+
  geom_hline(yintercept = 0)+
  theme(
    plot.title = element_text(color="Black", size=10, face="bold.italic"),
    axis.title.x = element_text(color="Black", size=8, face="bold"),
    axis.title.y = element_text(color="Black", size=8, face="bold"),
    legend.title = element_text(color="Black", size=8, face="bold")
  )+ 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
plot

#Make descriptive statistics
table(pd_hilic$indicator_fdrp)
df <- data.frame(row.names=c("Significant in both", 
                             "Significant only in PD patients", 
                             "Significant only in non-PD participants",
                             "Not significant in either group"),
                 freq=c(30,159,30,2694))
df <- as.data.frame(t(df))
tt <- ttheme_default(colhead=list(fg_params = list(parse=TRUE)),
                     base_size = 6,
                     padding = unit(c(2, 4), "mm"))
tbl <- tableGrob(df, rows=NULL, theme=tt)
combined<-ggarrange(plot, tbl,
                    ncol = 1, nrow = 2,
                    heights = c(1, 0.5))
combined



#Colored by crude p value, beta from C18 PD and Control
plot <- ggplot(pd_c18, aes(x = beta, y = control_c18$beta)) +
  geom_point(aes(color = indicator_fdrp)) +
  scale_colour_manual(values=c("#E69F00","#56B4E9","#009E73","#999999")) +
  geom_smooth(method = "lm", color = "black") +
  stat_cor(method = "pearson", size = 3)+
  labs(title="Correlation between metabolic feature coefficients among \n PD patients and non-PD participants in C18 column",
       x ="Metabolic feature coefficients among PD patients", 
       y = "Metabolic feature coefficients among non-PD participants",
       color= "FDR-adjusted P value significance level")+ 
  geom_vline(xintercept = 0)+
  geom_hline(yintercept = 0)+
  theme(
    plot.title = element_text(color="Black", size=10, face="bold.italic"),
    axis.title.x = element_text(color="Black", size=8, face="bold"),
    axis.title.y = element_text(color="Black", size=8, face="bold"),
    legend.title = element_text(color="Black", size=8, face="bold")
  )+ 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
plot

#Make descriptive statistics
table(pd_c18$indicator_fdrp)
df <- data.frame(row.names=c("Significant in both", 
                             "Significant only in PD patients", 
                             "Significant only in non-PD participants",
                             "Not significant in either group"),
                 freq=c(25,113,7,2077))
df <- as.data.frame(t(df))
tt <- ttheme_default(colhead=list(fg_params = list(parse=TRUE)),
                     base_size = 6,
                     padding = unit(c(2, 4), "mm"))
tbl <- tableGrob(df, rows=NULL, theme=tt)
combined2<-ggarrange(plot, tbl,
                     ncol = 1, nrow = 2,
                     heights = c(1, 0.5))
combined2

combined_final<-ggarrange(combined, combined2,
                          ncol = 2, nrow = 1,
                          labels = c("C", "D"),
                          heights = c(1, 1))
combined_final



png(filename="Manuscript/figure/cor_fdr_p.png",
    units="in", 
    width=12, 
    height=6, 
    res=300)
print(combined_final)
dev.off()


# 
# #Multinomial Comparison
# #3: both, 2: PD, 1: diabetes
# load("mwasResults/hilic_medcombat_diabetes_MWAS_multi_1.Rdata")
# load("mwasResults/hilic_medcombat_diabetes_MWAS_multi_2.Rdata")
# load("mwasResults/hilic_medcombat_diabetes_MWAS_multi_3.Rdata")
# 
# load("mwasResults/c18_medcombat_diabetes_MWAS_multi_1.Rdata")
# load("mwasResults/c18_medcombat_diabetes_MWAS_multi_2.Rdata")
# load("mwasResults/c18_medcombat_diabetes_MWAS_multi_3.Rdata")
# 
# #Hilic
# #P value indicator
# final_3_1_hilic$indicator_crudep_3v2<-
#   ifelse(final_3_1_hilic$pvalue<0.05&final_2_1_hilic$pvalue<0.05,"Significant in both",
#          ifelse(final_3_1_hilic$pvalue<0.05&final_2_1_hilic$pvalue>=0.05,"Significant only in diabetic PD patients",
#                 ifelse(final_3_1_hilic$pvalue>=0.05&final_2_1_hilic$pvalue<0.05,"Significant only in non-diabetic PD patients",
#                        "Not significant in either group")))
# 
# final_3_1_hilic$indicator_crudep_3v1<-
#   ifelse(final_3_1_hilic$pvalue<0.05&final_1_1_hilic$pvalue<0.05,"Significant in both",
#          ifelse(final_3_1_hilic$pvalue<0.05&final_1_1_hilic$pvalue>=0.05,"Significant only in diabetic PD patients",
#                 ifelse(final_3_1_hilic$pvalue>=0.05&final_1_1_hilic$pvalue<0.05,"Significant only in Diabetic non-PD participants",
#                        "Not significant in either group")))
# 
# 
# 
# final_3_1_hilic$indicator_crudep_3v2 <- factor(final_3_1_hilic$indicator_crudep_3v2, 
#                                     levels = c("Significant in both", 
#                                                "Significant only in diabetic PD patients", 
#                                                "Significant only in non-diabetic PD patients",
#                                                "Not significant in either group"))
# final_3_1_hilic$indicator_crudep_3v1 <- factor(final_3_1_hilic$indicator_crudep_3v1, 
#                                                levels = c("Significant in both", 
#                                                           "Significant only in diabetic PD patients", 
#                                                           "Significant only in Diabetic non-PD participants",
#                                                           "Not significant in either group"))
# 
# #3v2
# plot <- ggplot(final_3_1_hilic, aes(x = beta, y = final_2_1_hilic$beta)) +
#   geom_point(aes(color = indicator_crudep_3v2)) +
#   scale_colour_manual(values=c("#E69F00","#56B4E9","#009E73","#999999")) +
#   geom_smooth(method = "lm", color = "black") +
#   stat_cor(method = "pearson", size=3)+
#   labs(x ="Diabetic PD patients vs Non-diabetic non-PD participants", 
#        y = "Non-diabetic PD patients vs Non-diabetic non-PD participants",
#        color= "Crude P value significance level")+ 
#   theme(
#     axis.title.x = element_text(color="Black", size=6, face="bold"),
#     axis.title.y = element_text(color="Black", size=6, face="bold"),
#     legend.title = element_text(color="Black", size=6, face="bold"),
#     legend.text  = element_text(color="Black", size=6))+ 
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#         panel.background = element_blank(), axis.line = element_line(colour = "black"))
# plot
# table(final_3_1_hilic$indicator_crudep_3v2)
# df <- data.frame(row.names=c("Significant in both", 
#                              "Significant only in diabetic PD patients", 
#                              "Significant only in non-diabetic PD patients",
#                              "Not significant in either group"),
#                  freq=c(212,1301,198,1202))
# df <- as.data.frame(t(df))
# tt <- ttheme_default(colhead=list(fg_params = list(parse=TRUE)),
#                      base_size = 6,
#                      padding = unit(c(2, 4), "mm"))
# tbl <- tableGrob(df, rows=NULL, theme=tt)
# combined<-ggarrange(plot, tbl,
#                      ncol = 1, nrow = 2,
#                      heights = c(1, 0.5))
# combined
# 
# #3vs1
# plot <- ggplot(final_3_1_hilic, aes(x = beta, y = final_1_1_hilic$beta)) +
#   geom_point(aes(color = indicator_crudep_3v1)) +
#   scale_colour_manual(values=c("#E69F00","#56B4E9","#009E73","#999999")) +
#   geom_smooth(method = "lm", color = "black") +
#   stat_cor(method = "pearson", size=3)+
#   labs(x ="PD patients vs Non-diabetic non-PD participants", 
#        y = "Diabetic Diabetic non-PD participants vs Non-diabetic non-PD participants",
#        color= "Crude P value significance level")+ 
#   theme(
#     axis.title.x = element_text(color="Black", size=6, face="bold"),
#     axis.title.y = element_text(color="Black", size=6, face="bold"),
#     legend.title = element_text(color="Black", size=6, face="bold"),
#     legend.text  = element_text(color="Black", size=6))+ 
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#         panel.background = element_blank(), axis.line = element_line(colour = "black"))
# plot
# table(final_3_1_hilic$indicator_crudep_3v1)
# df <- data.frame(row.names=c("Significant in both", 
#                              "Significant only in diabetic PD patients", 
#                              "Significant only in diabetic non-PD participants",
#                              "Not significant in either group"),
#                  freq=c(654,859,263,1137))
# df <- as.data.frame(t(df))
# tt <- ttheme_default(colhead=list(fg_params = list(parse=TRUE)),
#                      base_size = 6,
#                      padding = unit(c(2, 4), "mm"))
# tbl <- tableGrob(df, rows=NULL, theme=tt)
# combined1<-ggarrange(plot, tbl,
#                     ncol = 1, nrow = 2,
#                     heights = c(1, 0.5))
# combined1
# 
# 
# combined2<-ggarrange(combined, combined1,
#                      ncol = 2, nrow = 1,
#                      heights = c(1, 1))
# combined2
# 
# 
# #C18
# #P value indicator
# final_3_1_c18$indicator_crudep_3v2<-
#   ifelse(final_3_1_c18$pvalue<0.05&final_2_1_c18$pvalue<0.05,"Significant in both",
#          ifelse(final_3_1_c18$pvalue<0.05&final_2_1_c18$pvalue>=0.05,"Significant only in diabetic PD patients",
#                 ifelse(final_3_1_c18$pvalue>=0.05&final_2_1_c18$pvalue<0.05,"Significant only in non-diabetic PD patients",
#                        "Not significant in either group")))
# 
# final_3_1_c18$indicator_crudep_3v1<-
#   ifelse(final_3_1_c18$pvalue<0.05&final_1_1_c18$pvalue<0.05,"Significant in both",
#          ifelse(final_3_1_c18$pvalue<0.05&final_1_1_c18$pvalue>=0.05,"Significant only in diabetic PD patients",
#                 ifelse(final_3_1_c18$pvalue>=0.05&final_1_1_c18$pvalue<0.05,"Significant only in Diabetic non-PD participants",
#                        "Not significant in either group")))
# 
# 
# 
# final_3_1_c18$indicator_crudep_3v2 <- factor(final_3_1_c18$indicator_crudep_3v2, 
#                                                levels = c("Significant in both", 
#                                                           "Significant only in diabetic PD patients", 
#                                                           "Significant only in non-diabetic PD patients",
#                                                           "Not significant in either group"))
# final_3_1_c18$indicator_crudep_3v1 <- factor(final_3_1_c18$indicator_crudep_3v1, 
#                                                levels = c("Significant in both", 
#                                                           "Significant only in diabetic PD patients", 
#                                                           "Significant only in Diabetic non-PD participants",
#                                                           "Not significant in either group"))
# 
# #3v2
# plot <- ggplot(final_3_1_c18, aes(x = beta, y = final_2_1_c18$beta)) +
#   geom_point(aes(color = indicator_crudep_3v2)) +
#   scale_colour_manual(values=c("#E69F00","#56B4E9","#009E73","#999999")) +
#   geom_smooth(method = "lm", color = "black") +
#   stat_cor(method = "pearson", size=3)+
#   labs(x ="Diabetic PD patients vs Non-diabetic non-PD participants", 
#        y = "Non-diabetic PD patients vs Non-diabetic non-PD participants",
#        color= "Crude P value significance level")+ 
#   theme(
#     axis.title.x = element_text(color="Black", size=6, face="bold"),
#     axis.title.y = element_text(color="Black", size=6, face="bold"),
#     legend.title = element_text(color="Black", size=6, face="bold"),
#     legend.text  = element_text(color="Black", size=6))+ 
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#         panel.background = element_blank(), axis.line = element_line(colour = "black"))
# plot
# table(final_3_1_c18$indicator_crudep_3v2)
# df <- data.frame(row.names=c("Significant in both", 
#                              "Significant only in diabetic PD patients", 
#                              "Significant only in non-diabetic PD patients",
#                              "Not significant in either group"),
#                  freq=c(197,993,138,894))
# df <- as.data.frame(t(df))
# tt <- ttheme_default(colhead=list(fg_params = list(parse=TRUE)),
#                      base_size = 6,
#                      padding = unit(c(2, 4), "mm"))
# tbl <- tableGrob(df, rows=NULL, theme=tt)
# combined<-ggarrange(plot, tbl,
#                     ncol = 1, nrow = 2,
#                     heights = c(1, 0.5))
# combined
# 
# #3vs1
# plot <- ggplot(final_3_1_c18, aes(x = beta, y = final_1_1_c18$beta)) +
#   geom_point(aes(color = indicator_crudep_3v1)) +
#   scale_colour_manual(values=c("#E69F00","#56B4E9","#009E73","#999999")) +
#   geom_smooth(method = "lm", color = "black") +
#   stat_cor(method = "pearson", size=3)+
#   labs(x ="PD patients vs Non-diabetic non-PD participants", 
#        y = "Diabetic Diabetic non-PD participants vs Non-diabetic non-PD participants",
#        color= "Crude P value significance level")+ 
#   theme(
#     axis.title.x = element_text(color="Black", size=6, face="bold"),
#     axis.title.y = element_text(color="Black", size=6, face="bold"),
#     legend.title = element_text(color="Black", size=6, face="bold"),
#     legend.text  = element_text(color="Black", size=6))+ 
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#         panel.background = element_blank(), axis.line = element_line(colour = "black"))
# plot
# table(final_3_1_c18$indicator_crudep_3v1)
# df <- data.frame(row.names=c("Significant in both", 
#                              "Significant only in diabetic PD patients", 
#                              "Significant only in diabetic non-PD participants",
#                              "Not significant in either group"),
#                  freq=c(602,588,194,838))
# df <- as.data.frame(t(df))
# tt <- ttheme_default(colhead=list(fg_params = list(parse=TRUE)),
#                      base_size = 6,
#                      padding = unit(c(2, 4), "mm"))
# tbl <- tableGrob(df, rows=NULL, theme=tt)
# combined1<-ggarrange(plot, tbl,
#                      ncol = 1, nrow = 2,
#                      heights = c(1, 0.5))
# combined1
# 
# 
# combined3<-ggarrange(combined, combined1,
#                      ncol = 2, nrow = 1,
#                      heights = c(1, 1))
# combined3
# 
# combined_final<-ggarrange(combined2, NULL,combined3,
#                           ncol = 1, nrow = 3,vjust=-4,hjust=0.2,
#                           labels=c("HILIC column","", "C18 column"),
#                           font.label=list(color="black",size=10,face="bold"),
#                           heights = c(1,0.2, 1))+
#   theme(plot.margin = margin(2,1,1,1, "cm"))
# combined_final
# combined_final_new<-annotate_figure(combined_final,
#                 top = text_grob("Correlation between metabolic feature coefficients from multinomial models", color = "black", face = "bold", size = 14))
# combined_final_new
# png(filename="Manuscript/figure/cor_crude_p_multinom.png",
#     units="in", 
#     width=12, 
#     height=12, 
#     res=300)
# print(combined_final_new)
# dev.off()
# 
# 
# ########################################
# #FDR P value
# #Hilic
# #P value indicator
# final_3_1_hilic$indicator_fdrp_3v2<-
#   ifelse(final_3_1_hilic$FDR_pvalue<0.05&final_2_1_hilic$FDR_pvalue<0.05,"Significant in both",
#          ifelse(final_3_1_hilic$FDR_pvalue<0.05&final_2_1_hilic$FDR_pvalue>=0.05,"Significant only in diabetic PD patients",
#                 ifelse(final_3_1_hilic$FDR_pvalue>=0.05&final_2_1_hilic$FDR_pvalue<0.05,"Significant only in non-diabetic PD patients",
#                        "Not significant in either group")))
# 
# final_3_1_hilic$indicator_fdrp_3v1<-
#   ifelse(final_3_1_hilic$FDR_pvalue<0.05&final_1_1_hilic$FDR_pvalue<0.05,"Significant in both",
#          ifelse(final_3_1_hilic$FDR_pvalue<0.05&final_1_1_hilic$FDR_pvalue>=0.05,"Significant only in diabetic PD patients",
#                 ifelse(final_3_1_hilic$FDR_pvalue>=0.05&final_1_1_hilic$FDR_pvalue<0.05,"Significant only in Diabetic non-PD participants",
#                        "Not significant in either group")))
# 
# 
# 
# final_3_1_hilic$indicator_fdrp_3v2 <- factor(final_3_1_hilic$indicator_fdrp_3v2, 
#                                                levels = c("Significant in both", 
#                                                           "Significant only in diabetic PD patients", 
#                                                           "Significant only in non-diabetic PD patients",
#                                                           "Not significant in either group"))
# final_3_1_hilic$indicator_fdrp_3v1 <- factor(final_3_1_hilic$indicator_fdrp_3v1, 
#                                                levels = c("Significant in both", 
#                                                           "Significant only in diabetic PD patients", 
#                                                           "Significant only in Diabetic non-PD participants",
#                                                           "Not significant in either group"))
# 
# #3v2
# plot <- ggplot(final_3_1_hilic, aes(x = beta, y = final_2_1_hilic$beta)) +
#   geom_point(aes(color = indicator_fdrp_3v2)) +
#   scale_colour_manual(values=c("#E69F00","#56B4E9","#009E73","#999999")) +
#   geom_smooth(method = "lm", color = "black") +
#   stat_cor(method = "pearson", size=3)+
#   labs(x ="Diabetic PD patients vs Non-diabetic non-PD participants", 
#        y = "Non-diabetic PD patients vs Non-diabetic non-PD participants",
#        color= "FDR-adjusted P value significance level")+ 
#   theme(
#     axis.title.x = element_text(color="Black", size=6, face="bold"),
#     axis.title.y = element_text(color="Black", size=6, face="bold"),
#     legend.title = element_text(color="Black", size=6, face="bold"),
#     legend.text  = element_text(color="Black", size=6))+ 
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#         panel.background = element_blank(), axis.line = element_line(colour = "black"))
# plot
# table(final_3_1_hilic$indicator_fdrp_3v2)
# df <- data.frame(row.names=c("Significant in both", 
#                              "Significant only in diabetic PD patients", 
#                              "Significant only in non-diabetic PD patients",
#                              "Not significant in either group"),
#                  freq=c(64,1169,55,1625))
# df <- as.data.frame(t(df))
# tt <- ttheme_default(colhead=list(fg_params = list(parse=TRUE)),
#                      base_size = 6,
#                      padding = unit(c(2, 4), "mm"))
# tbl <- tableGrob(df, rows=NULL, theme=tt)
# combined<-ggarrange(plot, tbl,
#                     ncol = 1, nrow = 2,
#                     heights = c(1, 0.5))
# combined
# 
# #3vs1
# plot <- ggplot(final_3_1_hilic, aes(x = beta, y = final_1_1_hilic$beta)) +
#   geom_point(aes(color = indicator_fdrp_3v1)) +
#   scale_colour_manual(values=c("#E69F00","#56B4E9","#009E73","#999999")) +
#   geom_smooth(method = "lm", color = "black") +
#   stat_cor(method = "pearson", size=3)+
#   labs(x ="PD patients vs Non-diabetic non-PD participants", 
#        y = "Diabetic Diabetic non-PD participants vs Non-diabetic non-PD participants",
#        color= "FDR-adjusted P value significance level")+ 
#   theme(
#     axis.title.x = element_text(color="Black", size=6, face="bold"),
#     axis.title.y = element_text(color="Black", size=6, face="bold"),
#     legend.title = element_text(color="Black", size=6, face="bold"),
#     legend.text  = element_text(color="Black", size=6))+ 
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#         panel.background = element_blank(), axis.line = element_line(colour = "black"))
# plot
# table(final_3_1_hilic$indicator_fdrp_3v1)
# df <- data.frame(row.names=c("Significant in both", 
#                              "Significant only in diabetic PD patients", 
#                              "Significant only in diabetic non-PD participants",
#                              "Not significant in either group"),
#                  freq=c(232,1001,113,1567))
# df <- as.data.frame(t(df))
# tt <- ttheme_default(colhead=list(fg_params = list(parse=TRUE)),
#                      base_size = 6,
#                      padding = unit(c(2, 4), "mm"))
# tbl <- tableGrob(df, rows=NULL, theme=tt)
# combined1<-ggarrange(plot, tbl,
#                      ncol = 1, nrow = 2,
#                      heights = c(1, 0.5))
# combined1
# 
# 
# combined2<-ggarrange(combined, combined1,
#                      ncol = 2, nrow = 1,
#                      heights = c(1, 1))
# combined2
# 
# 
# #C18
# #P value indicator
# final_3_1_c18$indicator_fdrp_3v2<-
#   ifelse(final_3_1_c18$FDR_pvalue<0.05&final_2_1_c18$FDR_pvalue<0.05,"Significant in both",
#          ifelse(final_3_1_c18$FDR_pvalue<0.05&final_2_1_c18$FDR_pvalue>=0.05,"Significant only in diabetic PD patients",
#                 ifelse(final_3_1_c18$FDR_pvalue>=0.05&final_2_1_c18$FDR_pvalue<0.05,"Significant only in non-diabetic PD patients",
#                        "Not significant in either group")))
# 
# final_3_1_c18$indicator_fdrp_3v1<-
#   ifelse(final_3_1_c18$FDR_pvalue<0.05&final_1_1_c18$FDR_pvalue<0.05,"Significant in both",
#          ifelse(final_3_1_c18$FDR_pvalue<0.05&final_1_1_c18$FDR_pvalue>=0.05,"Significant only in diabetic PD patients",
#                 ifelse(final_3_1_c18$FDR_pvalue>=0.05&final_1_1_c18$FDR_pvalue<0.05,"Significant only in Diabetic non-PD participants",
#                        "Not significant in either group")))
# 
# 
# 
# final_3_1_c18$indicator_fdrp_3v2 <- factor(final_3_1_c18$indicator_fdrp_3v2, 
#                                              levels = c("Significant in both", 
#                                                         "Significant only in diabetic PD patients", 
#                                                         "Significant only in non-diabetic PD patients",
#                                                         "Not significant in either group"))
# final_3_1_c18$indicator_fdrp_3v1 <- factor(final_3_1_c18$indicator_fdrp_3v1, 
#                                              levels = c("Significant in both", 
#                                                         "Significant only in diabetic PD patients", 
#                                                         "Significant only in Diabetic non-PD participants",
#                                                         "Not significant in either group"))
# 
# #3v2
# plot <- ggplot(final_3_1_c18, aes(x = beta, y = final_2_1_c18$beta)) +
#   geom_point(aes(color = indicator_fdrp_3v2)) +
#   scale_colour_manual(values=c("#E69F00","#56B4E9","#009E73","#999999")) +
#   geom_smooth(method = "lm", color = "black") +
#   stat_cor(method = "pearson", size=3)+
#   labs(x ="Diabetic PD patients vs Non-diabetic non-PD participants", 
#        y = "Non-diabetic PD patients vs Non-diabetic non-PD participants",
#        color= "FDR-adjusted P value significance level")+ 
#   theme(
#     axis.title.x = element_text(color="Black", size=6, face="bold"),
#     axis.title.y = element_text(color="Black", size=6, face="bold"),
#     legend.title = element_text(color="Black", size=6, face="bold"),
#     legend.text  = element_text(color="Black", size=6))+ 
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#         panel.background = element_blank(), axis.line = element_line(colour = "black"))
# plot
# table(final_3_1_c18$indicator_fdrp_3v2)
# df <- data.frame(row.names=c("Significant in both", 
#                              "Significant only in diabetic PD patients", 
#                              "Significant only in non-diabetic PD patients",
#                              "Not significant in either group"),
#                  freq=c(47,932,13,1230))
# df <- as.data.frame(t(df))
# tt <- ttheme_default(colhead=list(fg_params = list(parse=TRUE)),
#                      base_size = 6,
#                      padding = unit(c(2, 4), "mm"))
# tbl <- tableGrob(df, rows=NULL, theme=tt)
# combined<-ggarrange(plot, tbl,
#                     ncol = 1, nrow = 2,
#                     heights = c(1, 0.5))
# combined
# 
# #3vs1
# plot <- ggplot(final_3_1_c18, aes(x = beta, y = final_1_1_c18$beta)) +
#   geom_point(aes(color = indicator_fdrp_3v1)) +
#   scale_colour_manual(values=c("#E69F00","#56B4E9","#009E73","#999999")) +
#   geom_smooth(method = "lm", color = "black") +
#   stat_cor(method = "pearson", size=3)+
#   labs(x ="PD patients vs Non-diabetic non-PD participants", 
#        y = "Diabetic Diabetic non-PD participants vs Non-diabetic non-PD participants",
#        color= "FDR-adjusted P value significance level")+ 
#   theme(
#     axis.title.x = element_text(color="Black", size=6, face="bold"),
#     axis.title.y = element_text(color="Black", size=6, face="bold"),
#     legend.title = element_text(color="Black", size=6, face="bold"),
#     legend.text  = element_text(color="Black", size=6))+ 
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#         panel.background = element_blank(), axis.line = element_line(colour = "black"))
# plot
# table(final_3_1_c18$indicator_fdrp_3v1)
# df <- data.frame(row.names=c("Significant in both", 
#                              "Significant only in diabetic PD patients", 
#                              "Significant only in diabetic non-PD participants",
#                              "Not significant in either group"),
#                  freq=c(176,803,66,1177))
# df <- as.data.frame(t(df))
# tt <- ttheme_default(colhead=list(fg_params = list(parse=TRUE)),
#                      base_size = 6,
#                      padding = unit(c(2, 4), "mm"))
# tbl <- tableGrob(df, rows=NULL, theme=tt)
# combined1<-ggarrange(plot, tbl,
#                      ncol = 1, nrow = 2,
#                      heights = c(1, 0.5))
# combined1
# 
# 
# combined3<-ggarrange(combined, combined1,
#                      ncol = 2, nrow = 1,
#                      heights = c(1, 1))
# combined3
# 
# combined_final<-ggarrange(combined2, NULL,combined3,
#                           ncol = 1, nrow = 3,vjust=-4,hjust=0.2,
#                           labels=c("HILIC column","", "C18 column"),
#                           font.label=list(color="black",size=10,face="bold"),
#                           heights = c(1,0.2, 1))+
#   theme(plot.margin = margin(2,1,1,1, "cm"))
# combined_final
# combined_final_new<-annotate_figure(combined_final,
#                                     top = text_grob("Correlation between metabolic feature coefficients from multinomial models", color = "black", face = "bold", size = 14))
# combined_final_new
# png(filename="Manuscript/figure/cor_fdr_p_multinom.png",
#     units="in", 
#     width=12, 
#     height=12, 
#     res=300)
# print(combined_final_new)
# dev.off()
