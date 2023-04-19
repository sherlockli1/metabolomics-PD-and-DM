##################################################################################################################
# Code used to conduct analyses for "Untargeted serum metabolic profiling of diabetes mellitus among             #
#Parkinson’s disease patients"                                                                                   #
# Analytic approach developed by: Sherlock Li and Yuyuan Lin                                                     #
# Code by: Sherlock Li and Yuyuan Lin                                                                            #               
# Last Update: 4/19/2023                                                                                         #
##################################################################################################################


require(haven)
library("readxl")
library(tidyverse)
library(lubridate)
library(gtsummary)
library(mice)

#Data Processing
#medical checklists
#import peg1
peg1.1 <- read_sas("Medical Checklist/peg1_medical_clean.sas7bdat")
peg1.2 <- read_sas("Medical Checklist/peg1_medical2_clean.sas7bdat")
peg1<-merge(peg1.1,peg1.2,by="PEGID")

#import peg2
peg2 <- read_sas("Medical Checklist/peg2_medical_1_2.sas7bdat")

#import cgep
cgep1.12<-read_sas("Medical Checklist/cgep1_medical12_clean.sas7bdat")
cgep1.3<-read_sas("Medical Checklist/cgep1_medical3_clean.sas7bdat")
cgep1<-merge(cgep1.12,cgep1.3,by="PEGID")

peg1_1<-peg1[,c("PEGID",
                "D7_AdultOnsetDiabetes",
                "D7_AdultOnsetDiabetesAge",
                "M_Diabetes1")]
colnames(peg1_1)[4]<-"M_Diabetes"

peg2_1<-peg2[,c("PEGID",
                "D7_AdultOnsetDiabetes",
                "D7_AdultOnsetDiabetesAge",
                "M_Diabetes")]
cgep1_1<-cgep1[,c("PEGID",
                  "D7_AdultOnsetDiabetes",
                  "D7_AdultOnsetDiabetesAge",
                  "M_Diabetes")]

all_med<-rbind(peg1_1,peg2_1,cgep1_1)

covar<-read.csv("KeyVariable/Pest_cov_more.csv")
load("KeyVariable/DP_metabo_keyVar.RData")
additional<-data_keyVar_dp[!data_keyVar_dp$pegid %in% covar$pegid, ]
additional<-additional[,colnames(additional) %in% colnames(covar)]
covar<-bind_rows(covar,additional)
covar$Study<-ifelse(grepl('^1', covar$pegid),1,2)

#import HH
all_key_med<-merge(covar,all_med,by.x="pegid",by.y="PEGID")

hh<-read.csv("KeyVariable/hh_keyvar_new.csv")
hh<-hh[,c(-8:-9)]

all_key_med<-bind_rows(all_key_med,hh)

#start with PD patients, start with people diagnosis
all_key_med$D7_AdultOnsetDiabetes<-ifelse(all_key_med$D7_AdultOnsetDiabetes==-6|all_key_med$D7_AdultOnsetDiabetes==77,NA,all_key_med$D7_AdultOnsetDiabetes)
all_key_med$Diabetes_med_or_diag<-ifelse(all_key_med$M_Diabetes==1|all_key_med$D7_AdultOnsetDiabetes==0,1,0)
all_key_med$Diabetes_med_and_diag<-ifelse(all_key_med$M_Diabetes==1&all_key_med$D7_AdultOnsetDiabetes==1,1,0)
all_key_med$D7_AdultOnsetDiabetes<-as.numeric(all_key_med$D7_AdultOnsetDiabetes)
all_key_med$D7_AdultOnsetDiabetesAge<-ifelse(all_key_med$D7_AdultOnsetDiabetes==0|all_key_med$D7_AdultOnsetDiabetesAge<0,NA,all_key_med$D7_AdultOnsetDiabetesAge)

all_key_med$diabetes_dur=all_key_med$agenew-all_key_med$D7_AdultOnsetDiabetesAge
all_key_med$D7_AdultOnsetDiabetes<-ifelse(all_key_med$diabetes_dur<=0&!is.na(all_key_med$diabetes_dur),0,
                                         ifelse(all_key_med$D7_AdultOnsetDiabetes==0,0,
                                                ifelse(!is.na(all_key_med$D7_AdultOnsetDiabetes),1,NA)))

all_key_med$Hisapnic<-ifelse(all_key_med$Race==3,1,
                             ifelse(is.na(all_key_med$Race),NA,0))
all_key_med$edu_cat<-ifelse(all_key_med$A1_SchYrs<12,0,
                            ifelse(is.na(all_key_med$A1_SchYrs),NA,1))
all_key_med$Smoker<-ifelse(all_key_med$Smoker==2,1,all_key_med$Smoker)

load("data_combat_for_date/matrix_PEG_hilic_med_combat_full_dateadj.RData")
load("Data/hilic/PEG_hilic_link_new2.Rdata")

load("data_combat_for_date/matrix_PEG_c18_med_combat_full_dateadj.RData")
load("Data/c18/PEG_c18_link_new2.Rdata")

#Imputation
#The patient refused to answer the question
all_key_med_pd<-all_key_med[which(all_key_med$PD==1),]
table(all_key_med_pd$D7_AdultOnsetDiabetes)
#Add the imputed value
all_key_med_pd$edu_cat<-ifelse(all_key_med_pd$pegid=="11454RJ30",1,all_key_med_pd$edu_cat)
all_key_med_pd$edu_cat<-ifelse(all_key_med_pd$pegid=="12426GC26",1,all_key_med_pd$edu_cat)
all_key_med_pd$edu_cat<-ifelse(all_key_med_pd$pegid=="12502JB39",1,all_key_med_pd$edu_cat)
all_key_med_pd$edu_cat<-ifelse(all_key_med_pd$pegid=="12532TL49",1,all_key_med_pd$edu_cat)
all_key_med_pd$edu_cat<-ifelse(all_key_med_pd$pegid=="86766NH47",1,all_key_med_pd$edu_cat)
all_key_med_pd$edu_cat<-ifelse(all_key_med_pd$pegid=="89832MW70",1,all_key_med_pd$edu_cat)
all_key_med_pd$edu_cat<-ifelse(all_key_med_pd$pegid=="89833DC42",1,all_key_med_pd$edu_cat)

# tidy data --------------------------------------------------------------------
#hilic
data_link <- PEG_hilic_link_new2 %>%
  # extract pegid and sample date
  filter(Emory_qc == 0) %>%
  mutate(pegid = str_extract(Sample.ID.hilic, 
                             pattern = "\\d{5}[:upper:]{2}\\d{2}"),
         date = str_extract(Sample.ID.hilic, pattern = "\\d{8}"),
         date = mdy(date)) %>%
  select(c("FullRunName", "pegid", "date")) %>%
  # drop quality control samples
  filter(!is.na(pegid)) %>%
  # get the baseline only (early sample date)
  arrange(date) %>%
  filter(!duplicated(pegid)) %>%
  select(-c(date))

data_metabo_hilic_medcombat_case<-matrix_PEG_hilic_med_combat_full_dateadj %>%
  as.data.frame() %>%
  rownames_to_column(var="ID") %>%
  dplyr::filter(grepl('_baseline',ID))%>%
  dplyr::mutate(ID = str_extract(ID, "[^_]+"))%>%
  dplyr::filter(ID %in% all_key_med_pd$pegid) %>%
  column_to_rownames("ID")

all_key_med_pd_1 <- all_key_med_pd %>%
  filter(pegid %in% rownames(data_metabo_hilic_medcombat_case))
table(rownames(data_metabo_hilic_medcombat_case) == all_key_med_pd_1$pegid)
X <- as.matrix(data_metabo_hilic_medcombat_case)

# save data
save(data_link, file = "tmpData/patient_link_hilic_new.rdata")
save(X, file = "tmpData/X_hilic_medcombat_new.rdata")
save(all_key_med_pd_1, file = "tmpData/all_key_med_pd_1_new.rdata")

#c18
data_link_2 <- PEG_c18_link_new2 %>%
  # extract pegid and sample date
  filter(Emory_qc == 0) %>%
  mutate(pegid = str_extract(Sample.ID, 
                             pattern = "\\d{5}[:upper:]{2}\\d{2}"),
         date = str_extract(Sample.ID, pattern = "\\d{8}"),
         date = mdy(date)) %>%
  select(c("FullRunName", "pegid", "date")) %>%
  # drop quality control samples
  filter(!is.na(pegid)) %>%
  # get the baseline only (early sample date)
  arrange(date) %>%
  filter(!duplicated(pegid)) %>%
  select(-c(date))


data_metabo_c18_medcombat_case<-matrix_PEG_c18_med_combat_full_dateadj %>%
  as.data.frame() %>%
  rownames_to_column(var="ID") %>%
  dplyr::filter(grepl('_baseline',ID))%>%
  dplyr::mutate(ID = str_extract(ID, "[^_]+"))%>%
  dplyr::filter(ID %in% all_key_med_pd$pegid) %>%
  column_to_rownames("ID")

all_key_med_pd_2 <- all_key_med_pd %>%
  filter(pegid %in% rownames(data_metabo_c18_medcombat_case))
table(rownames(data_metabo_c18_medcombat_case) == all_key_med_pd_2$pegid)
X_2 <- as.matrix(data_metabo_c18_medcombat_case)

# save data
save(data_link_2, file = "tmpData/patient_link_c18_new.rdata")
save(X_2, file = "tmpData/X_c18_medcombat_new.rdata")
save(all_key_med_pd_2, file = "tmpData/all_key_med_pd_2_new.rdata")


#Get Table 1
all_key_med_pd_2 %>% 
  select(D7_AdultOnsetDiabetes,agenew, Sex, edu_cat,Smoker,Hisapnic,Study) %>% 
  tbl_summary(by=D7_AdultOnsetDiabetes, 
              statistic = list(all_continuous() ~ "{mean} ({sd})",
                               all_categorical() ~ "{n} ({p}%)"),
              digits = all_continuous() ~ 2,
              percent= "column",
              missing_text = "(Missing)")

out <- all_key_med_pd %>%
  filter(!pegid %in% rownames(data_metabo_c18_medcombat_case))
out %>% 
  select(D7_AdultOnsetDiabetes,agenew, Sex, edu_cat,Smoker,Hisapnic,Study) %>% 
  tbl_summary(by=D7_AdultOnsetDiabetes, 
              statistic = list(all_continuous() ~ "{mean} ({sd})",
                               all_categorical() ~ "{n} ({p}%)"),
              digits = all_continuous() ~ 2,
              percent= "column",
              missing_text = "(Missing)")


# univariate MWAS --------------------------------------------------------------
dvList <- colnames(data_metabo_hilic_medcombat_case)

#hilic
MWAS_result <- lapply(dvList, function(x) {
  summary(glm(substitute(all_key_med_pd_1$D7_AdultOnsetDiabetes ~                     
                           i +   
                           all_key_med_pd_1$agenew + 
                           all_key_med_pd_1$Sex + 
                           all_key_med_pd_1$edu_cat +
                           all_key_med_pd_1$Smoker + 
                           all_key_med_pd_1$Hisapnic+
                           all_key_med_pd_1$Study,
                         list(i = as.name(x))),
              data = as.data.frame(X), 
              family = binomial(link = "logit")))$coefficients[c(2), 1:4]
})

hilic_medcombat_diabetes_MWAS <- MWAS_result %>%
  bind_rows() %>%
  mutate(metabo = dvList) %>%
  dplyr::rename(
    beta = `Estimate`,
    se   = `Std. Error`,
    zvalue = `z value`,
    pvalue = `Pr(>|z|)`
  ) %>%
  mutate(FDR_pvalue = p.adjust(pvalue, method = "fdr")) %>%
  relocate(metabo)

save(hilic_medcombat_diabetes_MWAS, 
     file = "mwasResults/hilic_medcombat_diabetes_MWAS_new.rdata")

#############################################################################
#c18
dvList <- colnames(data_metabo_c18_medcombat_case)

MWAS_result <- lapply(dvList, function(x) {
  summary(glm(substitute(all_key_med_pd_2$D7_AdultOnsetDiabetes ~ i + 
                           all_key_med_pd_2$agenew + 
                           all_key_med_pd_2$Sex + 
                           all_key_med_pd_2$edu_cat +
                           all_key_med_pd_2$Smoker + 
                           all_key_med_pd_2$Hisapnic+
                           all_key_med_pd_2$Study,
                         list(i = as.name(x))),
              data = as.data.frame(X_2), 
              family = binomial(link = "logit")))$coefficients[c(2), 1:4]
})


c18_medcombat_diabetes_MWAS <- MWAS_result %>%
  bind_rows() %>%
  mutate(metabo = dvList) %>%
  dplyr::rename(
    beta = `Estimate`,
    se   = `Std. Error`,
    zvalue = `z value`,
    pvalue = `Pr(>|z|)`
  ) %>%
  mutate(FDR_pvalue = p.adjust(pvalue, method = "fdr")) %>%
  relocate(metabo)

save(c18_medcombat_diabetes_MWAS, 
     file = "mwasResults/c18_medcombat_diabetes_MWAS_new.rdata")



#########
#Mummichog input
# save file for Mummichog pathway analysis -------------------------------------
# hilic
load("data/hilic/Comb_features_hilic.RData")

tbl_mum_input <- hilic_medcombat_diabetes_MWAS %>%
  left_join(Comb_features, by = c("metabo" = "match")) %>%
  dplyr::select(c(mz.adj, time.adj, pvalue, zvalue)) %>%
  dplyr::rename(`mz` = mz.adj, 
         `rtime` = time.adj, 
         `z-score` = zvalue, 
         `p-value` = pvalue)
write.table(tbl_mum_input, file = "mumInput/hilic_diabetes_new.txt",
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

test <- hilic_medcombat_diabetes_MWAS %>%
  left_join(Comb_features, by = c("metabo" = "match")) %>%
  dplyr::select(c(mz.adj, time.adj, pvalue, zvalue)) %>%
  dplyr::rename(`mz` = mz.adj, 
                `rtime` = time.adj, 
                `z-score` = zvalue, 
                `p-value` = pvalue)%>%
  dplyr::filter(`p-value`<0.05)
  
write.table(test, file = "mumInput/testhilic_diabetes_new.txt",
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

# tbl_metapone_input_hilic <- hilic_medcombat_diabetes_MWAS %>%
#   left_join(Comb_features, by = c("metabo" = "match")) %>%
#   dplyr::select(c(mz.adj, time.adj, beta, se, pvalue, zvalue)) %>%
#   dplyr::rename(`mz` = mz.adj, 
#                 `rtime` = time.adj, 
#                 `beta` = beta,
#                 `se` = se,
#                 `z-score` = zvalue, 
#                 `p-value` = pvalue)
# save(tbl_metapone_input_hilic, 
#      file = "mwasResults/tbl_metapone_input_hilic.rdata")


# c18
load("data/c18/Comb_features_c18.RData")

tbl_mum_input <- c18_medcombat_diabetes_MWAS %>%
  left_join(Comb_features, by = c("metabo" = "match")) %>%
  dplyr::select(c(mz.adj, time.adj, pvalue, zvalue)) %>%
  dplyr::rename(`mz` = mz.adj, 
         `rtime` = time.adj, 
         `z-score` = zvalue, 
         `p-value` = pvalue)
write.table(tbl_mum_input, file = "mumInput/c18_diabetes_new.txt",
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")


# tbl_metapone_input_c18 <- c18_medcombat_diabetes_MWAS %>%
#   left_join(Comb_features, by = c("metabo" = "match")) %>%
#   dplyr::select(c(mz.adj, time.adj, beta, se, pvalue, zvalue)) %>%
#   dplyr::rename(`mz` = mz.adj, 
#                 `rtime` = time.adj, 
#                 `beta` = beta,
#                 `se` = se,
#                 `z-score` = zvalue, 
#                 `p-value` = pvalue)
# save(tbl_metapone_input_c18, 
#      file = "mwasResults/tbl_metapone_input_c18.rdata")

#PLS Models
library(doParallel)
library(foreach)
library(mixOmics)
library(tibble)
no_cores <- 7
cl <- makeCluster(no_cores)

##############################################################################################################
Y_dia <- all_key_med_pd_1$D7_AdultOnsetDiabetes

# calculate residuals
registerDoParallel(cl)
X_resid <- foreach(i = 1:ncol(X), .combine = cbind) %dopar% {
  residuals(lm(X[, i] ~ Sex + edu_cat+Smoker+Hisapnic+Study, 
               na.action = na.exclude, data = all_key_med_pd_1))
}
stopImplicitCluster()
colnames(X_resid) <- colnames(X)
rownames(X_resid) <- rownames(X)

# original PLS -----------------------------------------------------------------
optimal.ncomp <- 3

pls_hilic_dia_una <- plsda(X_resid, Y_dia,
                             ncomp = optimal.ncomp)

####################################
Y_dia <- all_key_med_pd_2$D7_AdultOnsetDiabetes

# calculate residuals
registerDoParallel(cl)
X_resid_2 <- foreach(i = 1:ncol(X_2), .combine = cbind) %dopar% {
  residuals(lm(X_2[, i] ~ Sex + edu_cat+Smoker+Hisapnic+Study, 
               na.action = na.exclude, data = all_key_med_pd_2))
}
stopImplicitCluster()
colnames(X_resid_2) <- colnames(X_2)
rownames(X_resid_2) <- rownames(X_2)

pls_c18_dia_una <- plsda(X_resid_2, Y_dia,
                           ncomp = optimal.ncomp)

# save -------------------------------------------------------------------------
save(pls_hilic_dia_una, file = "plsModel/pls_hilic_dia_una_new.rdata")
save(pls_c18_dia_una, file = "plsModel/pls_c18_dia_una_new.rdata")

#Make Mummichog Input
hilic_pd_vip <- as.data.frame(vip(pls_hilic_dia_una))
c18_pd_vip <- as.data.frame(vip(pls_c18_dia_una))

hilic_pd_vip$fakepvalue<-ifelse(hilic_pd_vip$comp1>2,0.001,0.08)

