library(haven)
library(dplyr)
library(interactions)
library(data.table)
library(jtools)
library(xlsx)
library(car)
library(ggplot2)

## full pheno
XX <- read_sav("/share/projects/alspac/Data/RawDataB2730/Parent_25Nov2019.sav") %>%
  mutate(UniqueID = paste0(cidB2730, qlet)) %>%    ## combine Unique pregnancy identifier (cidB2730) and Birth order (qlet) to get a unique ID for every kid
  as.data.frame()  

PCs <- fread("/share/projects/genodata/ALSPAC/DOWNSTREAM_ANALYSIS/PCA/PCA_2020Feb19/EIGENSOFT_2020Feb19/ALSPAC_unrelated_KIDS_2020Feb19/ALSPAC_unrelated_KIDS_ALL_only_genotyped_SNPs_2020Mar03/indep_100_5_1.01_by_alspac_without_MHC/alspac_kids_nodup_genotyped_without_MHC_maf05_pruned.pca_allSample.evec",
             header = F) %>% 
  mutate(UniqueID = substr(V1, 1, regexpr(":", V1)-1)) %>%
  dplyr::select(-c(V1, V12))
colnames(PCs)[1:10] <- paste0("PC", 1:10)


##### standardize and rename (e)PRSes #####
c782 <- fread("/n9/cortexjobs/c_eprs/jobs/c1285_eamon_acc.micro.alspac.kids.acc_c1285/result/c1285_prs.score.csv", header = T, stringsAsFactors = F) %>% 
  mutate(PRS_1.0 = as.numeric(scale(PRS_1.0, center = T, scale = T)))
c830 <- fread("/n9/cortexjobs/c_eprs/jobs/c830_eamon_outer.radial.glia.alspac.kids.fetaleQTLs_c830/scores_without_any_start_end_filter/result/c830_prs.score.csv", header = T, stringsAsFactors = F) %>% 
  mutate(PRS_1.0 = as.numeric(scale(PRS_1.0, center = T, scale = T)))
c829 <- fread("/n9/cortexjobs/c_eprs/jobs/c829_eamon_oligodendrocyte.precursor.cells.alspac.kids.fetaleQTLs_c829/scores_without_any_start_end_filter/result/c829_prs.score.csv", header = T, stringsAsFactors = F) %>% 
  mutate(PRS_1.0 = as.numeric(scale(PRS_1.0, center = T, scale = T)))
c828 <- fread("/n9/cortexjobs/c_eprs/jobs/c828_eamon_inner.layer.neurons.alspac.kids.fetaleQTLs_c828/scores_without_any_start_end_filter/result/c828_prs.score.csv", header = T, stringsAsFactors = F) %>% 
  mutate(PRS_1.0 = as.numeric(scale(PRS_1.0, center = T, scale = T)))

colnames(c782) <- c('UniqueID', 'SNP_count_microglia', 'ePRS_microglia')
colnames(c830) <- c('UniqueID', 'SNP_count_oRG', 'ePRS_oRG')
colnames(c829) <- c('UniqueID', 'SNP_count_OPCs', 'ePRS_OPCs')
colnames(c828) <- c('UniqueID', 'SNP_count_inner_neurons', 'ePRS_inner_neurons')

names.PRS <- c("PRS_0.0001", "PRS_0.001", "PRS_0.01", "PRS_0.05", "PRS_0.1", "PRS_0.2", "PRS_0.5")

b248 <- fread('/share/projects/cortexjobs/b_prs/jobs/b248_eamon_crossdisorderpgc.alspac.kids_b248/result/b248_prs.score.csv', header = T, stringsAsFactors = F) %>% 
  mutate_at(names.PRS, ~(scale(.) %>% as.vector))
b250 <- fread("/n9/cortexjobs/b_prs/jobs/b250_eamon_mddHoward19.alspac.kids_b250/result/b250_prs.score.csv", header = T, stringsAsFactors = F) %>% 
  mutate_at(names.PRS, ~(scale(.) %>% as.vector))

colnames(b248) <- c('UniqueID', paste0(colnames(b248)[-1], '_', 'crossdisorderpgc'))
colnames(b250) <- c('UniqueID', paste0(colnames(b250)[-1], '_', 'depression'))

ePRS <- merge(c782, c830, by = "UniqueID", all = T) %>%
  merge(      c829, by = "UniqueID", all = T) %>%
  merge(      c828, by = "UniqueID", all = T)

XX <- full_join(XX, ePRS,         by = "UniqueID") %>%
  full_join(    b250,          by = "UniqueID") %>%
  full_join(    PCs,          by = "UniqueID") 

#### get averages of prenatal maternal mood ####
#### CCEI and EPDS ####

XX <- XX %>% mutate(CCEI_Anxiety_mean_prenatal = rowMeans(cbind(b351, c573), na.rm = T)) %>%
  mutate(CCEI_Depression_mean_prenatal = rowMeans(cbind(b353, c579), na.rm = T)) %>%
  mutate(CCEI_Somatic_mean_prenatal = rowMeans(cbind(b355, c576), na.rm = T)) %>%
  mutate(CCEI_Total_mean_prenatal = rowMeans(cbind(b357, c582), na.rm = T)) %>%
  mutate(EPDS_mean_prenatal = rowMeans(cbind(b370, c600), na.rm = T))

#----- log CRP -----
XX$logCRP_TF3_15yrs <- log(XX$crp_TF3)
XX$logCRP_TF4_17.5yrs <- log(XX$CRP_TF4)


###########################
### exclusion criteria ####
###########################

XX <- XX %>% dplyr::filter ((bestgest >= 37) & (bestgest <= 42) &
                              (mz028b >= 18) & (kz030 >= 2000) & (mz010 == 1) &
                              (mz013 > 1) & (mz014 > 1) & (kz011b == 1))

# kz021,  ## gender: 1=male 2=female
XX$Gender <- NA
XX$Gender[which(XX$kz021 == 1)] <- "Male"
XX$Gender[which(XX$kz021 == 2)] <- "Female"

# Any dep variable
XX <- XX %>%
  rowwise() %>%
  mutate(Any_dep = as.integer(FKDQ1000 == 1 | FKDQ1010 == 1 | FKDQ1020 == 1)) %>%
  ungroup()

XX_female <- subset(XX, Gender =="Female")

g <- data.frame(mild = XX$FKDQ1000,
                mod = XX$FKDQ1020,
                sev = XX$FKDQ1020,
                toget = XX$Any_dep)
rownames(XX_female) <- XX_female$UniqueID
XX.cont <- XX_female


#----- Any Dep 24 years Females --------
glm.cont.old <- glm(Any_dep ~ b023 + bestgest + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + 
                      PC7 + PC8 + PC9 + PC10 + ePRS_microglia * EPDS_mean_prenatal, data = XX.cont, family = binomial)
summary(glm.cont.old)
nrow_coef <- nrow(summary(glm.cont.old)$coef)

im.cont <- data.frame(influence.measures(glm.cont.old)$is.inf)
case.exclude.cont <- rownames(im.cont[which(im.cont[,nrow_coef] == "TRUE"),])
if (length(case.exclude.cont)>0){XX.cont <- XX.cont[-which(XX.cont$UniqueID %in% case.exclude.cont),]}

glm.cont     <- glm(Any_dep ~ b023 + bestgest + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + 
                      PC7 + PC8 + PC9 + PC10 + ePRS_microglia * EPDS_mean_prenatal, data = XX.cont, family = binomial)
summary(glm.cont)



#Call:
#  glm(formula = Any_dep ~ b023 + bestgest + PC1 + PC2 + PC3 + PC4 + 
#        PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + ePRS_microglia * EPDS_mean_prenatal, 
#      family = binomial, data = XX.cont)

#Deviance Residuals: 
#  Min       1Q   Median       3Q      Max  
#-1.0006  -0.5469  -0.4682  -0.3950   2.4650  

#Coefficients:
#  Estimate Std. Error z value Pr(>|z|)   
#(Intercept)                         1.542651   2.622233   0.588  0.55633   
#b023                               -0.056700   0.017883  -3.171  0.00152 **
#  bestgest                           -0.060760   0.064205  -0.946  0.34397   
#PC1                                 1.907524   7.503753   0.254  0.79933   
#PC2                                 4.692534   7.196250   0.652  0.51435   
#PC3                                -4.904125   7.156773  -0.685  0.49319   
#PC4                                 0.624935   7.071126   0.088  0.92958   
#PC5                                 7.702378   7.157411   1.076  0.28186   
#PC6                                 4.213223   6.958846   0.605  0.54488   
#PC7                                -7.813324   7.447819  -1.049  0.29414   
#PC8                               -14.124721   7.150398  -1.975  0.04823 * 
#  PC9                                 7.030063   7.245760   0.970  0.33193   
#PC10                               -5.326249   7.290702  -0.731  0.46505   
#ePRS_microglia                      0.141740   0.151036   0.938  0.34801   
#EPDS_mean_prenatal                  0.047878   0.018568   2.579  0.00992 **
#ePRS_microglia:EPDS_mean_prenatal  -0.009657   0.018363  -0.526  0.59894   
