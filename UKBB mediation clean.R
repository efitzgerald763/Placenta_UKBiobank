###########################################################################
## exclude dropouts
## exclude related subjects
## exclude cases with inconsistent genetic and self-reported sex
## adjust for genotyping array + assessment center + 40PCs
###########################################################################

library(dplyr)
library(tidyverse)
library(data.table)
library(tidyr)


## keep a list of data-fields you need, otherwise will take too long to loaad
key <- fread("/share/projects/uk_biobank/UKBioBank_Key_all_14Aug2019.csv", header = T, stringsAsFactors = F) %>%
  filter(field.showcase %in% c(20544, 41270, 21003, 21022, 22027, 22021, 21001, 50, 12144, 20160, 22000, 22006, 54, 21000, 22001,
                               20022, 20116, 189, 22009, 23127))

XX  <- fread("/share/projects/uk_biobank/UKB_pheno_part1_part2_14Aug2019.csv", header = T, stringsAsFactors = F, 
             select = c("eid", "sex_f31_0_0", 
                        key$col.name))

###################################################### 

my_vars <- fread("/n0/uk_biobank/created.variables/f20544_f41270I20.25_f41270I25.9.csv")

####################################################
## 40 PCs are usually considered in all analyses
colnames(XX) <- gsub("genetic_principal_components_f22009_0_", "PC", colnames(XX))


## renaming to use a shorter name in analyses
setnames(XX, "age_at_recruitment_f21022_0_0", "age_f21022")


###### PRSs ######

##### standardize and rename PGS #####

## Cyan PGS
PRS_cyan <- fread("/n9/cortexjobs/UKbio_score/u_eprs/jobs/u163_eamon_placenta.rpkm.cyan22.ukbio.placenta_u163/Final_genotype/PRS_result/u163_prs.score.csv", header = T, stringsAsFactors = F) %>% 
  mutate(PRS_1.0 = as.numeric(scale(PRS_1.0, center = T, scale = T))) %>%
  dplyr::rename(eid=ID1)

colnames(PRS_cyan) <- c('eid', 'SNP_count_cyan', 'PRS_cyan')

## MDD and Cardio PRS
names.PRS <- c("Final_PRS_0.0001", "Final_PRS_0.001", "Final_PRS_0.01", "Final_PRS_0.05", "Final_PRS_0.1", "Final_PRS_0.2", "Final_PRS_0.5")

MDD_PRS <- fread("/n9/cortexjobs/UKbio_score/v_prs/jobs/v41_eamon_mdd.2019.ukbio_v41/Final_genotype/PRS_result/v41_prs_score.csv", header = T, stringsAsFactors = F) %>% 
  mutate_at(names.PRS, ~(scale(.) %>% as.vector))%>%
  dplyr::rename(eid=ID1)
colnames(MDD_PRS) <- c('eid', paste0(colnames(MDD_PRS)[-1], '_', 'MDD'))

Cardio_PRS <- fread("/n9/cortexjobs/UKbio_score/v_prs/jobs/v42_eamon_cardiovascular.disease.ukbio_v42/Final_genotype/PRS_result/v42_prs_score.csv", header = T, stringsAsFactors = F) %>% 
  mutate_at(names.PRS, ~(scale(.) %>% as.vector))%>%
  dplyr::rename(eid=ID1)
colnames(Cardio_PRS) <- c('eid', paste0(colnames(Cardio_PRS)[-1], '_', 'Cardio'))


# Merge
PRS <- merge(PRS_cyan, MDD_PRS, by = "eid", all = T) %>%
  merge(      Cardio_PRS, by = "eid", all = T) 


## Combine all data sets
XX <- full_join(XX,  my_vars,     by = "eid")%>%
  full_join(PRS,    by = "eid")



#######################################
########## exclude drop-outs ##########
#######################################
#dropouts <- fread("/share/projects/uk_biobank/w41975_20210809.csv")
zihan_used_dropouts <- fread("/n0/uk_biobank/w41975_20210201.csv")

XX <- XX %>% filter(!(eid %in% zihan_used_dropouts$V1))

###########################################
## remove cases with inconsistent sexes ##
###########################################
XX <- XX %>% filter(sex_f31_0_0 == genetic_sex_f22001_0_0)


########################################
### get only unrelated individuals  ####
## exclude outliers_for_heterozygosity #
########################################
# EXCLUDE subjects according to:
# Outliers for heterozygosity or missing rate	22027-0.0 - Indicates samples identified as outliers in heterozygosity and 
# missing rates, which implies that the genotypes for these samples are of poor quality.
#
# Genetic kinship to other participants	22021-0.0 > 0.044 - a further 131,790 related individuals based on a shared relatedness of up to the
# third degree using kinship coefficients (>0.044) calculated using the KING toolset

XX <- XX[-which(XX$outliers_for_heterozygosity_or_missing_rate_f22027_0_0 == "Yes"),]



## include unrelated subjects:

## run INSTEAD XX <- XX[which(XX$genetic_kinship_to_other_participants_f22021_0_0 == "No kinship found"),]

unrelated.id <- fread("/share/projects/uk_biobank/Unrelated.subset_of.related.subjects/UKB_66096_unrelated_subjects_in_related_sample_28July2020.csv")

XX$unrelated_66k<-0

XX$unrelated_66k[which(XX$eid %in% unrelated.id$q.id.V1)]<-1

XX$unrelated_66k[which(XX$genetic_kinship_to_other_participants_f22021_0_0=='Ten or more third-degree relatives identified')]<-0

XX <- XX[which(XX$genetic_kinship_to_other_participants_f22021_0_0 == "No kinship found" | XX$unrelated_66k == 1),]

# Assessment centre and genotyping arrray - usually included in all analyses
XX$uk_biobank_assessment_centre_f54_0_0   <- relevel(as.factor(XX$uk_biobank_assessment_centre_f54_0_0), ref = "11010")

 
# UKBiLEVEAX [-11,-1]
# Axiom [1,95]
XX <- XX %>% dplyr::mutate(genotype_array=ifelse(genotype_measurement_batch_f22000_0_0 %in% c(-11:-1), yes=1, no=0)) # yes="BiLEVEAX", no="Axiom"

XX$sex <- NA
XX$sex[which(XX$sex_f31_0_0 == "Male")] <- 1
XX$sex[which(XX$sex_f31_0_0 == "Female")] <- 2

XX_male <- XX[XX$sex == 1,]
XX_female <- XX[XX$sex == 2,]

######################## 
###### Regression ###### 
######################## 
XX = as.data.frame(XX)
rownames(XX_female) <- XX_female$eid
XX.cont <- XX_female

sink(paste0("~/UKBB_Med_Results_", Sys.Date() ,".txt",sep=""))

Final_PRS_0.00000001_Cardio <- glm(f41270_I25.9_Chronic_ischaemic_heart_disease ~ age_f21022 + PC1 + 
                     PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 +
                     PC11 + PC12 + PC13 + PC14 + PC15 + PC16 + PC17 + PC18 + PC19 + PC20 + PC21 + PC22 + PC23 +
                     PC24 + PC25 + PC26 + PC27 + PC28 + PC29 + PC30 + PC31 + PC32 + PC33 + PC34 + PC35 + PC36 + 
                     PC37 + PC38 + PC39 + PC40 + as.factor(genotype_array) + uk_biobank_assessment_centre_f54_0_0 + Final_PRS_0.00000001_Cardio +
                     PRS_cyan, family = binomial(link = "logit"), data = XX.cont, na.action=na.exclude)
summary(Final_PRS_0.00000001_Cardio)

Final_PRS_0.0000001_Cardio <- glm(f41270_I25.9_Chronic_ischaemic_heart_disease ~ age_f21022 + PC1 + 
                     PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 +
                     PC11 + PC12 + PC13 + PC14 + PC15 + PC16 + PC17 + PC18 + PC19 + PC20 + PC21 + PC22 + PC23 +
                     PC24 + PC25 + PC26 + PC27 + PC28 + PC29 + PC30 + PC31 + PC32 + PC33 + PC34 + PC35 + PC36 + 
                     PC37 + PC38 + PC39 + PC40 + as.factor(genotype_array) + uk_biobank_assessment_centre_f54_0_0 + Final_PRS_0.0000001_Cardio +
                     PRS_cyan, family = binomial(link = "logit"), data = XX.cont, na.action=na.exclude)
summary(Final_PRS_0.0000001_Cardio)

Final_PRS_0.000001_Cardio <- glm(f41270_I25.9_Chronic_ischaemic_heart_disease ~ age_f21022 + PC1 + 
                     PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 +
                     PC11 + PC12 + PC13 + PC14 + PC15 + PC16 + PC17 + PC18 + PC19 + PC20 + PC21 + PC22 + PC23 +
                     PC24 + PC25 + PC26 + PC27 + PC28 + PC29 + PC30 + PC31 + PC32 + PC33 + PC34 + PC35 + PC36 + 
                     PC37 + PC38 + PC39 + PC40 + as.factor(genotype_array) + uk_biobank_assessment_centre_f54_0_0 + Final_PRS_0.000001_Cardio +
                     PRS_cyan, family = binomial(link = "logit"), data = XX.cont, na.action=na.exclude)
summary(Final_PRS_0.000001_Cardio)

Final_PRS_0.00001_Cardio <- glm(f41270_I25.9_Chronic_ischaemic_heart_disease ~ age_f21022 + PC1 + 
                     PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 +
                     PC11 + PC12 + PC13 + PC14 + PC15 + PC16 + PC17 + PC18 + PC19 + PC20 + PC21 + PC22 + PC23 +
                     PC24 + PC25 + PC26 + PC27 + PC28 + PC29 + PC30 + PC31 + PC32 + PC33 + PC34 + PC35 + PC36 + 
                     PC37 + PC38 + PC39 + PC40 + as.factor(genotype_array) + uk_biobank_assessment_centre_f54_0_0 + Final_PRS_0.00001_Cardio +
                     PRS_cyan, family = binomial(link = "logit"), data = XX.cont, na.action=na.exclude)
summary(Final_PRS_0.00001_Cardio)

Final_PRS_0.0001_Cardio <- glm(f41270_I25.9_Chronic_ischaemic_heart_disease ~ age_f21022 + PC1 + 
                     PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 +
                     PC11 + PC12 + PC13 + PC14 + PC15 + PC16 + PC17 + PC18 + PC19 + PC20 + PC21 + PC22 + PC23 +
                     PC24 + PC25 + PC26 + PC27 + PC28 + PC29 + PC30 + PC31 + PC32 + PC33 + PC34 + PC35 + PC36 + 
                     PC37 + PC38 + PC39 + PC40 + as.factor(genotype_array) + uk_biobank_assessment_centre_f54_0_0 + Final_PRS_0.0001_Cardio +
                     PRS_cyan, family = binomial(link = "logit"), data = XX.cont, na.action=na.exclude)
summary(Final_PRS_0.0001_Cardio)


Final_PRS_0.001_Cardio <- glm(f41270_I25.9_Chronic_ischaemic_heart_disease ~ age_f21022 + PC1 + 
                     PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 +
                     PC11 + PC12 + PC13 + PC14 + PC15 + PC16 + PC17 + PC18 + PC19 + PC20 + PC21 + PC22 + PC23 +
                     PC24 + PC25 + PC26 + PC27 + PC28 + PC29 + PC30 + PC31 + PC32 + PC33 + PC34 + PC35 + PC36 + 
                     PC37 + PC38 + PC39 + PC40 + as.factor(genotype_array) + uk_biobank_assessment_centre_f54_0_0 + Final_PRS_0.001_Cardio +
                     PRS_cyan, family = binomial(link = "logit"), data = XX.cont, na.action=na.exclude)
summary(Final_PRS_0.001_Cardio)

Final_PRS_0.01_Cardio <- glm(f41270_I25.9_Chronic_ischaemic_heart_disease ~ age_f21022 + PC1 + 
                     PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 +
                     PC11 + PC12 + PC13 + PC14 + PC15 + PC16 + PC17 + PC18 + PC19 + PC20 + PC21 + PC22 + PC23 +
                     PC24 + PC25 + PC26 + PC27 + PC28 + PC29 + PC30 + PC31 + PC32 + PC33 + PC34 + PC35 + PC36 + 
                     PC37 + PC38 + PC39 + PC40 + as.factor(genotype_array) + uk_biobank_assessment_centre_f54_0_0 + Final_PRS_0.01_Cardio +
                     PRS_cyan, family = binomial(link = "logit"), data = XX.cont, na.action=na.exclude)
summary(Final_PRS_0.01_Cardio)

Final_PRS_0.1_Cardio <- glm(f41270_I25.9_Chronic_ischaemic_heart_disease ~ age_f21022 + PC1 + 
                     PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 +
                     PC11 + PC12 + PC13 + PC14 + PC15 + PC16 + PC17 + PC18 + PC19 + PC20 + PC21 + PC22 + PC23 +
                     PC24 + PC25 + PC26 + PC27 + PC28 + PC29 + PC30 + PC31 + PC32 + PC33 + PC34 + PC35 + PC36 + 
                     PC37 + PC38 + PC39 + PC40 + as.factor(genotype_array) + uk_biobank_assessment_centre_f54_0_0 + Final_PRS_0.1_Cardio +
                     PRS_cyan, family = binomial(link = "logit"), data = XX.cont, na.action=na.exclude)
summary(Final_PRS_0.1_Cardio)

Final_PRS_0.2_Cardio <- glm(f41270_I25.9_Chronic_ischaemic_heart_disease ~ age_f21022 + PC1 + 
                     PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 +
                     PC11 + PC12 + PC13 + PC14 + PC15 + PC16 + PC17 + PC18 + PC19 + PC20 + PC21 + PC22 + PC23 +
                     PC24 + PC25 + PC26 + PC27 + PC28 + PC29 + PC30 + PC31 + PC32 + PC33 + PC34 + PC35 + PC36 + 
                     PC37 + PC38 + PC39 + PC40 + as.factor(genotype_array) + uk_biobank_assessment_centre_f54_0_0 + Final_PRS_0.2_Cardio +
                     PRS_cyan, family = binomial(link = "logit"), data = XX.cont, na.action=na.exclude)
summary(Final_PRS_0.2_Cardio)

Final_PRS_0.5_Cardio <- glm(f41270_I25.9_Chronic_ischaemic_heart_disease ~ age_f21022 + PC1 + 
                     PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 +
                     PC11 + PC12 + PC13 + PC14 + PC15 + PC16 + PC17 + PC18 + PC19 + PC20 + PC21 + PC22 + PC23 +
                     PC24 + PC25 + PC26 + PC27 + PC28 + PC29 + PC30 + PC31 + PC32 + PC33 + PC34 + PC35 + PC36 + 
                     PC37 + PC38 + PC39 + PC40 + as.factor(genotype_array) + uk_biobank_assessment_centre_f54_0_0 + Final_PRS_0.5_Cardio +
                     PRS_cyan, family = binomial(link = "logit"), data = XX.cont, na.action=na.exclude)
summary(Final_PRS_0.5_Cardio)

MDD_in_CVD <- glm(f41270_I25.9_Chronic_ischaemic_heart_disease ~ age_f21022 + PC1 + 
                              PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 +
                              PC11 + PC12 + PC13 + PC14 + PC15 + PC16 + PC17 + PC18 + PC19 + PC20 + PC21 + PC22 + PC23 +
                              PC24 + PC25 + PC26 + PC27 + PC28 + PC29 + PC30 + PC31 + PC32 + PC33 + PC34 + PC35 + PC36 + 
                              PC37 + PC38 + PC39 + PC40 + as.factor(genotype_array) + uk_biobank_assessment_centre_f54_0_0 + f20544_Mental_health_problems_diagnosed_Depression +
                              PRS_cyan, family = binomial(link = "logit"), data = XX.cont, na.action=na.exclude)
summary(MDD_in_CVD)

#####################################
############### MDD #################
#####################################
Final_PRS_0.00000001_MDD <- glm(f20544_Mental_health_problems_diagnosed_Depression ~ age_f21022 + PC1 + 
                                     PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 +
                                     PC11 + PC12 + PC13 + PC14 + PC15 + PC16 + PC17 + PC18 + PC19 + PC20 + PC21 + PC22 + PC23 +
                                     PC24 + PC25 + PC26 + PC27 + PC28 + PC29 + PC30 + PC31 + PC32 + PC33 + PC34 + PC35 + PC36 + 
                                     PC37 + PC38 + PC39 + PC40 + as.factor(genotype_array) + uk_biobank_assessment_centre_f54_0_0 + Final_PRS_0.00000001_MDD +
                                     PRS_cyan, family = binomial(link = "logit"), data = XX.cont, na.action=na.exclude)
summary(Final_PRS_0.00000001_MDD)

Final_PRS_0.0000001_MDD <- glm(f20544_Mental_health_problems_diagnosed_Depression ~ age_f21022 + PC1 + 
                                    PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 +
                                    PC11 + PC12 + PC13 + PC14 + PC15 + PC16 + PC17 + PC18 + PC19 + PC20 + PC21 + PC22 + PC23 +
                                    PC24 + PC25 + PC26 + PC27 + PC28 + PC29 + PC30 + PC31 + PC32 + PC33 + PC34 + PC35 + PC36 + 
                                    PC37 + PC38 + PC39 + PC40 + as.factor(genotype_array) + uk_biobank_assessment_centre_f54_0_0 + Final_PRS_0.0000001_MDD +
                                    PRS_cyan, family = binomial(link = "logit"), data = XX.cont, na.action=na.exclude)
summary(Final_PRS_0.0000001_MDD)

Final_PRS_0.000001_MDD <- glm(f20544_Mental_health_problems_diagnosed_Depression ~ age_f21022 + PC1 + 
                                   PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 +
                                   PC11 + PC12 + PC13 + PC14 + PC15 + PC16 + PC17 + PC18 + PC19 + PC20 + PC21 + PC22 + PC23 +
                                   PC24 + PC25 + PC26 + PC27 + PC28 + PC29 + PC30 + PC31 + PC32 + PC33 + PC34 + PC35 + PC36 + 
                                   PC37 + PC38 + PC39 + PC40 + as.factor(genotype_array) + uk_biobank_assessment_centre_f54_0_0 + Final_PRS_0.000001_MDD +
                                   PRS_cyan, family = binomial(link = "logit"), data = XX.cont, na.action=na.exclude)
summary(Final_PRS_0.000001_MDD)

Final_PRS_0.00001_MDD <- glm(f20544_Mental_health_problems_diagnosed_Depression ~ age_f21022 + PC1 + 
                                  PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 +
                                  PC11 + PC12 + PC13 + PC14 + PC15 + PC16 + PC17 + PC18 + PC19 + PC20 + PC21 + PC22 + PC23 +
                                  PC24 + PC25 + PC26 + PC27 + PC28 + PC29 + PC30 + PC31 + PC32 + PC33 + PC34 + PC35 + PC36 + 
                                  PC37 + PC38 + PC39 + PC40 + as.factor(genotype_array) + uk_biobank_assessment_centre_f54_0_0 + Final_PRS_0.00001_MDD +
                                  PRS_cyan, family = binomial(link = "logit"), data = XX.cont, na.action=na.exclude)
summary(Final_PRS_0.00001_MDD)

Final_PRS_0.0001_MDD <- glm(f20544_Mental_health_problems_diagnosed_Depression ~ age_f21022 + PC1 + 
                                 PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 +
                                 PC11 + PC12 + PC13 + PC14 + PC15 + PC16 + PC17 + PC18 + PC19 + PC20 + PC21 + PC22 + PC23 +
                                 PC24 + PC25 + PC26 + PC27 + PC28 + PC29 + PC30 + PC31 + PC32 + PC33 + PC34 + PC35 + PC36 + 
                                 PC37 + PC38 + PC39 + PC40 + as.factor(genotype_array) + uk_biobank_assessment_centre_f54_0_0 + Final_PRS_0.0001_MDD +
                                 PRS_cyan, family = binomial(link = "logit"), data = XX.cont, na.action=na.exclude)
summary(Final_PRS_0.0001_MDD)


Final_PRS_0.001_MDD <- glm(f20544_Mental_health_problems_diagnosed_Depression ~ age_f21022 + PC1 + 
                                PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 +
                                PC11 + PC12 + PC13 + PC14 + PC15 + PC16 + PC17 + PC18 + PC19 + PC20 + PC21 + PC22 + PC23 +
                                PC24 + PC25 + PC26 + PC27 + PC28 + PC29 + PC30 + PC31 + PC32 + PC33 + PC34 + PC35 + PC36 + 
                                PC37 + PC38 + PC39 + PC40 + as.factor(genotype_array) + uk_biobank_assessment_centre_f54_0_0 + Final_PRS_0.001_MDD +
                                PRS_cyan, family = binomial(link = "logit"), data = XX.cont, na.action=na.exclude)
summary(Final_PRS_0.001_MDD)

Final_PRS_0.01_MDD <- glm(f20544_Mental_health_problems_diagnosed_Depression ~ age_f21022 + PC1 + 
                               PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 +
                               PC11 + PC12 + PC13 + PC14 + PC15 + PC16 + PC17 + PC18 + PC19 + PC20 + PC21 + PC22 + PC23 +
                               PC24 + PC25 + PC26 + PC27 + PC28 + PC29 + PC30 + PC31 + PC32 + PC33 + PC34 + PC35 + PC36 + 
                               PC37 + PC38 + PC39 + PC40 + as.factor(genotype_array) + uk_biobank_assessment_centre_f54_0_0 + Final_PRS_0.01_MDD +
                               PRS_cyan, family = binomial(link = "logit"), data = XX.cont, na.action=na.exclude)
summary(Final_PRS_0.01_MDD)

Final_PRS_0.1_MDD <- glm(f20544_Mental_health_problems_diagnosed_Depression ~ age_f21022 + PC1 + 
                              PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 +
                              PC11 + PC12 + PC13 + PC14 + PC15 + PC16 + PC17 + PC18 + PC19 + PC20 + PC21 + PC22 + PC23 +
                              PC24 + PC25 + PC26 + PC27 + PC28 + PC29 + PC30 + PC31 + PC32 + PC33 + PC34 + PC35 + PC36 + 
                              PC37 + PC38 + PC39 + PC40 + as.factor(genotype_array) + uk_biobank_assessment_centre_f54_0_0 + Final_PRS_0.1_MDD +
                              PRS_cyan, family = binomial(link = "logit"), data = XX.cont, na.action=na.exclude)
summary(Final_PRS_0.1_MDD)

Final_PRS_0.2_MDD <- glm(f20544_Mental_health_problems_diagnosed_Depression ~ age_f21022 + PC1 + 
                              PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 +
                              PC11 + PC12 + PC13 + PC14 + PC15 + PC16 + PC17 + PC18 + PC19 + PC20 + PC21 + PC22 + PC23 +
                              PC24 + PC25 + PC26 + PC27 + PC28 + PC29 + PC30 + PC31 + PC32 + PC33 + PC34 + PC35 + PC36 + 
                              PC37 + PC38 + PC39 + PC40 + as.factor(genotype_array) + uk_biobank_assessment_centre_f54_0_0 + Final_PRS_0.2_MDD +
                              PRS_cyan, family = binomial(link = "logit"), data = XX.cont, na.action=na.exclude)
summary(Final_PRS_0.2_MDD)

Final_PRS_0.5_MDD <- glm(f20544_Mental_health_problems_diagnosed_Depression ~ age_f21022 + PC1 + 
                              PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 +
                              PC11 + PC12 + PC13 + PC14 + PC15 + PC16 + PC17 + PC18 + PC19 + PC20 + PC21 + PC22 + PC23 +
                              PC24 + PC25 + PC26 + PC27 + PC28 + PC29 + PC30 + PC31 + PC32 + PC33 + PC34 + PC35 + PC36 + 
                              PC37 + PC38 + PC39 + PC40 + as.factor(genotype_array) + uk_biobank_assessment_centre_f54_0_0 + Final_PRS_0.5_MDD +
                              PRS_cyan, family = binomial(link = "logit"), data = XX.cont, na.action=na.exclude)
summary(Final_PRS_0.5_MDD)

CVD_in_MDD <- glm(f20544_Mental_health_problems_diagnosed_Depression ~ age_f21022 + PC1 + 
                           PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 +
                           PC11 + PC12 + PC13 + PC14 + PC15 + PC16 + PC17 + PC18 + PC19 + PC20 + PC21 + PC22 + PC23 +
                           PC24 + PC25 + PC26 + PC27 + PC28 + PC29 + PC30 + PC31 + PC32 + PC33 + PC34 + PC35 + PC36 + 
                           PC37 + PC38 + PC39 + PC40 + as.factor(genotype_array) + uk_biobank_assessment_centre_f54_0_0 + f41270_I25.9_Chronic_ischaemic_heart_disease +
                           PRS_cyan, family = binomial(link = "logit"), data = XX.cont, na.action=na.exclude)
summary(CVD_in_MDD)

sink()