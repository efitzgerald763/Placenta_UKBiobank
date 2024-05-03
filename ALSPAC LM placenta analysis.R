setwd("***")
library(haven)
library(dplyr)
library(interactions)
library(data.table)
library(jtools)
library(xlsx)
library(car)
library(ggplot2)

getwd()
#####################################################################
########## PART I: LOAD, MERGE, AND PROCESS DATA SETS ###############
#####################################################################

## full pheno
XX <- read_spss("/share/projects/alspac/Data/RawDataB2730/Parent_25Nov2019.sav") %>%
  mutate(UniqueID = paste0(cidB2730, qlet)) %>%    ## combine Unique pregnancy identifier (cidB2730) and Birth order (qlet) to get a unique ID for every kid
  as.data.frame()  

B_weight <- read_spss("/n0/alspac/Analysis/ABC/Parent_ABC_08Oct2019.sav") %>%
  mutate(UniqueID = paste0(cidB2730, qlet)) %>% 
  dplyr::select(UniqueID, BWeightZScore, BWeightCentile) %>% 
  as.data.frame()  
  
##########################
## get PCs kids without MHC
## 13Jan2021
##########################
PCs <- fread("/share/projects/genodata/ALSPAC/DOWNSTREAM_ANALYSIS/PCA/PCA_2020Feb19/EIGENSOFT_2020Feb19/ALSPAC_unrelated_KIDS_2020Feb19/ALSPAC_unrelated_KIDS_ALL_only_genotyped_SNPs_2020Mar03/indep_100_5_1.01_by_alspac_without_MHC/alspac_kids_nodup_genotyped_without_MHC_maf05_pruned.pca_allSample.evec",
             header = F) %>% 
  mutate(UniqueID = substr(V1, 1, regexpr(":", V1)-1)) %>%
  dplyr::select(-c(V1, V12))
colnames(PCs)[1:10] <- paste0("PC", 1:10)



##### standardize and rename (e)PRSes #####
Cyan <- fread("/n9/cortexjobs/c_eprs/jobs/c841_eamon_placenta.rpkm.cyan22.alspac.kids.placenta_c841/result/c841_prs.score.csv", header = T, stringsAsFactors = F) %>% 
  mutate(PRS_1.0 = as.numeric(scale(PRS_1.0, center = T, scale = T)))
FDR05 <- fread("/n9/cortexjobs/c_eprs/jobs/c843_eamon_placenta.modules.fdr.0.05.alspac.kids.placenta_c843/result/c843_prs.score.csv", header = T, stringsAsFactors = F) %>% 
  mutate(PRS_1.0 = as.numeric(scale(PRS_1.0, center = T, scale = T)))
FDR10 <- fread("/n9/cortexjobs/c_eprs/jobs/c845_eamon_placenta.modules.fdr.0.1.alspac.kids.placenta_c845/result/c845_prs.score.csv", header = T, stringsAsFactors = F) %>% 
  mutate(PRS_1.0 = as.numeric(scale(PRS_1.0, center = T, scale = T)))


colnames(Cyan) <- c('UniqueID', 'SNP_count_Cyan', 'ePRS_Cyan')
colnames(FDR05) <- c('UniqueID', 'SNP_count_FDR05', 'ePRS_FDR05')
colnames(FDR10) <- c('UniqueID', 'SNP_count_FDR10', 'ePRS_FDR10')

ePRS <- merge(Cyan, FDR05, by = "UniqueID", all = T) %>%
  merge(      FDR10, by = "UniqueID", all = T) 


## ABC ## description available in alspac/Analysis/ABC
ABC_Oct2019 <- read_spss("/share/projects/alspac/Analysis/ABC/Parent_ABC_08Oct2019_neat.sav")[, c(1, 5:10)]



XX <- full_join(XX, ePRS,         by = "UniqueID") %>%
  full_join(    PCs,          by = "UniqueID") %>%
  full_join(    B_weight,          by = "UniqueID") %>%
  full_join(ABC_Oct2019,      by = "UniqueID")


# kz021,  ## gender: 1=male 2=female
XX$Gender <- NA
XX$Gender[which(XX$kz021 == 1)] <- "Male"
XX$Gender[which(XX$kz021 == 2)] <- "Female"



list.ePRS <- colnames(ePRS)[which(grepl("PRS", colnames(ePRS)))]
list.ePRS.names <- list.ePRS



###########################
### exclusion criteria ####
###########################
#               e111 any infection- data derived from mother completed questionaire at 8 weeks postnatal
#               c755 social class mother 
#b_seg_m   B: socio economic group (mother) 18 weeks
#c_seg_m   C: socio economic group (mother) 32 weeks
#f_seg_m    F: socio economic group (mother) 8 months
#g_seg_m   G: socio economic group (mother) 21 months
# 							bestgest, gestational length
# 							mz028b, ## Grouped age of mother at delivery
# 							kz030,  ## Preferred birthweight
# 							mz010,  ## Pregancy size (i.e. singletons, twins..)
# 							mz013,  ## 28-day survivors
# 							mz014,  ## 1-year survivors
# 							kz011b, ## Participant was alive at 1 year of age
#XX <- XX %>% dplyr::filter ((bestgest >= 37) & (bestgest <= 42) & (mz028b >= 18) & (kz030 >= 2000) & (mz010 == 1) & (mz013 > 1) & (mz014 > 1) & (kz011b == 1))
XX <- XX %>% dplyr::filter ((mz010 == 1) & (mz013 > 1) & (mz014 > 1) & (kz011b == 1))
XX <- XX %>% mutate(Perinatal_SES = rowMeans(cbind(b_seg_m, c_seg_m, f_seg_m, g_seg_m), na.rm = T))
XX <- XX %>% mutate(Perinatal_EPDS = rowMeans(cbind(b370, e390, f200, g290, c600), na.rm = T))

# Split males and females
XX_female <- subset(XX, Gender =="Female")
XX_male <- subset(XX, Gender =="Male")


XX <- XX_female

#####################################################################
########## PART II: EXAMPLE 1 ONE INTERACTION/REGRESSION ############
#####################################################################
#CRP_f9
#crp_FOM1
#crp_TF3
#CRP_TF4
#IL6_f9
#Insulin_BBS
#Insulin_cord
#insulin_F9
#insulin_TF3
#insulin_TF4
#insulingp_TF4
#Glucose_BBS
#glucose_TF3
#glucose_TF4
#Leptin_CIF61
#LEPTIN_f9







#ep
j557a
kq346c
n8365a
ku707b
kw6602b
ta7025a
tc4025a

#pp
j557d
kq346e
n8365d
ku709b
kw6604b
ta7025d
tc4025d

#ps
j557e
kq346a
n8365e
ku705b
kw6600b
ta7025e
tc4025e

rownames(XX) <- XX$UniqueID
XX.cont <- XX
lm.cont.old <- lm(tc4025e ~  bestgest + e111 + Perinatal_SES + bestgest + 
                    Perinatal_EPDS + mz028b + kz030 + PC1 + PC2 + PC3 + PC4 + 
                    PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + ePRS_Cyan, data = XX.cont)
summary(lm.cont.old)

nrow_coef <- nrow(summary(lm.cont.old)$coef)
im.cont <- data.frame(influence.measures(lm.cont.old)$is.inf)
case.exclude.cont <- rownames(im.cont[which(im.cont[,nrow_coef] == "TRUE"),])
if (length(case.exclude.cont)>0){XX.cont <- XX[-which(XX$UniqueID %in% case.exclude.cont),]}

lm.cont     <- lm(tc4025e ~ bestgest + e111 + Perinatal_SES + bestgest + 
                    Perinatal_EPDS + mz028b + kz030 + PC1 + PC2 + PC3 + PC4 + 
                    PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + ePRS_Cyan, data = XX.cont)
summary(lm.cont)
confint(lm.cont, 'ePRS_Cyan', level=0.95)

#sink("example_full_results.txt")  ## start - save full output to a text file
## diagnostics, check assumptions ##
par(mfrow = c(2, 2))
plot(lm.cont)
par(mfrow=c(1,1))
boxCox(lm.cont,lambda=c(-2,2,0.5))

lm.cont     <- lm(log (kw6602b +1) ~ Gender + bestgest + e111 + Perinatal_SES + bestgest + 
                    Perinatal_EPDS + mz028b + kz030 + PC1 + PC2 + PC3 + PC4 + 
                    PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + ePRS_Cyan, data = XX.cont)
summary(lm.cont)

######################################################################
####### PART III: EXAMPLE 2 LOOP ALL INTERACTIONS/REGRESSIONS ########
######################################################################

list.Y <- c('j557a',
            'j557b',
            'j557c',
            'j557d',
            'j557e',
            'j557f',
            'kq346a',
            'kq346b',
            'kq346c',
            'kq346d',
            'kq346e',
            'kq346f',
            'ku705b',
            'ku706b',
            'ku707b',
            'ku708b',
            'ku709b',
            'ku710b',
            'kw6600b',
            'kw6601b',
            'kw6602b',
            'kw6603b',
            'kw6604b',
            'kw6605b',
            'ta7025a',
            'ta7025b',
            'ta7025c',
            'ta7025d',
            'ta7025e',
            'ta7025f',
            'tc4025a',
            'tc4025b',
            'tc4025c',
            'tc4025d',
            'tc4025e',
            'tc4025f',
            'n8365a',
            'n8365b',
            'n8365c',
            'n8365d',
            'n8365e',
            'n8365f')

list.Y.names <- c('SDQ - emotional symptoms score (prorated) - 47m',
                  'SDQ - conduct problems score (prorated) - 47m',
                  'SDQ - hyperactivity score (prorated) - 47m',
                  'SDQ - peer problems score (prorated) - 47m',
                  'SDQ - prosocial score (prorated) - 47m',
                  'SDQ - total difficulties score (prorated) - 47m',
                  'SDQ - Prosocial score (prorated) 81m',
                  'SDQ - Hyperactivity score (prorated) 81m',
                  'SDQ - Emotional symptoms score (prorated) 81m',
                  'SDQ - Conduct problems score (prorated) 81m',
                  'SDQ - Peer problems score (prorated) 81m',
                  'SDQ - Total behavioural Difficulties score (prorated) 81m',
                  'SDQ - Prosocial score (prorated) - 9y7m',
                  'SDQ - Hyperactivity score (prorated) - 9y7m',
                  'SDQ - Emotional symptoms score (prorated) - 9y7m',
                  'SDQ - Conduct problems score (prorated) - 9y7m',
                  'SDQ - Peer problems score (prorated) - 9y7m',
                  'SDQ - Total difficulties score (prorated) - 9y7m',
                  'SDQ - Prosocial score (prorated) - 11y8m',
                  'SDQ - Hyperactivity score (prorated) - 11y8m',
                  'SDQ - Emotional symptoms score (prorated) - 11y8m',
                  'SDQ - Conduct problems score (prorated) - 11y8m',
                  'SDQ - Peer problems score (prorated) - 11y8m',
                  'SDQ - Total difficulties score (prorated) - 11y8m',
                  'SDQ - emotional symptoms score (prorated) - 13y1m',
                  'SDQ - conduct problems score (prorated) - 13y1m',
                  'SDQ - hyperactivity score (prorated) - 13y1m',
                  'SDQ - peer problems score (prorated) - 13y1m',
                  'SDQ - prosocial score (prorated) - 13y1m',
                  'SDQ - total difficulties score (prorated) - 13y1m',
                  'SDQ - emotional symptoms score (prorated) - 16y6m',
                  'SDQ - conduct problems score (prorated) - 16y6m',
                  'SDQ - hyperactivity score (prorated) - 16y6m',
                  'SDQ - peer problems score (prorated) - 16y6m',
                  'SDQ - prosocial score (prorated) - 16y6m',
                  'SDQ - total difficulties score (prorated) - 16y6m',
                  'SDQ - emotional symptoms score (prorated) - 8y1m',
                  'SDQ - conduct problems score (prorated) - 8y1m',
                  'SDQ - hyperactivity score (prorated) - 8y1m',
                  'SDQ - peer problems score (prorated) - 8y1m',
                  'SDQ - prosocial score (prorated) - 8y1m',
                  'SDQ - total difficulties score (prorated) - 8y1m')


#############################
#### Interactions   #########
#############################
options(digits=3)


warnings(2) ## to output warnings too with sink()

PP <- data.frame(matrix(NA, nrow = length(list.ePRS) * length(list.Y), ncol = 9)) 
colnames(PP) <- c('Main_effect_ID',
                  'Outcome',
                  'Predictor',
                  'Sample_size',
                  'Sig_before_inf_measure_Out',
                  'Cases_excluded',
                  'p_value',
                  'Beta',
                  'SE_Beta')

sink(paste0("~/Results_", Sys.Date() ,".txt",sep=""))

XX = as.data.frame(XX)
rownames(XX) <- XX$UniqueID
iID <- 0
########################### loop starts here ##############################
for (k in 1:length(list.ePRS)){  
  XX$ePRS.cur <-eval(parse(text=paste("XX$", list.ePRS[k], sep = "")))
  summary(XX$X)
  
  for (i in 1:length(list.Y)){  
    XX$Y <- eval(parse(text=paste0("XX$", list.Y[i])))
    
    
    XX.cont <- XX[complete.cases(XX[, c("Y", "e111", "Perinatal_SES", "Perinatal_EPDS" , "mz028b", "kz030", "bestgest", "PC1","ePRS.cur")]), ]
    
    
    lm.cont.old <- lm(Y ~ bestgest + e111 + Perinatal_SES + Perinatal_EPDS + mz028b + kz030 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + ePRS.cur, data = XX.cont)
    
    nrow_coef <- nrow(summary(lm.cont.old)$coef)
    
    im.cont <- data.frame(influence.measures(lm.cont.old)$is.inf)
    case.exclude <- rownames(im.cont[which(im.cont$dffit == "TRUE"),])
    if (length(case.exclude)>0){XX.cont <- XX.cont[-which(XX.cont$UniqueID %in% case.exclude),]}
    
    
    lm.cont     <- lm(Y ~ bestgest + e111 + Perinatal_SES + Perinatal_EPDS + mz028b + kz030 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + ePRS.cur, data = XX.cont)
    
    print(summary(lm.cont))
    
    nrow.cur <- nrow(summary(lm.cont)$coef)
    sum.cur <- summary(lm.cont)
    df.cur  <- ifelse(is.null(sum.cur$df.residual), sum.cur$df[2], sum.cur$df.residual)
    

    print(list.Y.names[i])
    print(list.ePRS.names[k])
    print(iID+1)
    print(summary(lm.cont))

    ## one line for each main effect
    iID <- iID + 1
    PP[iID, 1] <- iID
    PP$Predictor[iID] <- list.ePRS.names[k]
    PP$Outcome[iID]   <- list.Y.names[i]
    PP$p_value[iID] <- summary(lm.cont)$coefficients[nrow.cur,4]
    PP$Beta[iID]    <- summary(lm.cont)$coefficients[nrow.cur,1]
    PP$SE_Beta[iID] <- summary(lm.cont)$coefficients[nrow.cur,2]
    PP$Sample_size[iID] <- df.cur + nrow.cur
    PP$Cases_excluded[iID] <- length(case.exclude)
    PP$Sig_before_inf_measure_Out[iID] <- ifelse(summary(lm.cont.old)$coefficients[nrow.cur,4]<0.05,"Yes","No")
  }
}
sink()
########################### loop stops here ##############################

Cyan_res <- subset(PP, Predictor =="ePRS_Cyan")
Sub5 <- subset(PP, Predictor =="ePRS_FDR05")
Sub1 <- subset(PP, Predictor =="ePRS_FDR10")

#Adjust p-values for multiple comparisons
Cyan_res = Cyan_res[order(Cyan_res$p_value),]
Sub5 = Sub5[order(Sub5$p_value),]
Sub1 = Sub1[order(Sub1$p_value),]

Cyan_res$Bonferroni =
  p.adjust(Cyan_res$p_value,
           method = "bonferroni")
Sub5$Bonferroni =
  p.adjust(Sub5$p_value,
           method = "bonferroni")
Sub1$Bonferroni =
  p.adjust(Sub1$p_value,
           method = "bonferroni")
Cyan_res$BH =
  p.adjust(Cyan_res$p_value,
           method = "BH")
Sub5$BH =
  p.adjust(Sub5$p_value,
           method = "BH")
Sub1$BH =
  p.adjust(Sub1$p_value,
           method = "BH")

write_csv(Cyan_res, "ALSPAC_SDQ_Placenta_Cyan_Male_Results_Sept_23rd.csv")

full <-fread("/data1/meaneylab/eamon/ALSPAC_SDQ_Placenta_Cyan_MandF_Results_Sept_23rd.csv")
Male <-fread("/data1/meaneylab/eamon/ALSPAC_SDQ_Placenta_Cyan_Male_Results_Sept_23rd.csv")
Female <-fread("/data1/meaneylab/eamon/ALSPAC_SDQ_Placenta_Cyan_Female_Results_Sept_23rd.csv")

SD <- function(x) full$SE_Beta * sqrt(full$Sample_size)

full$SD <- SD(full)
# SE=SD/sqrt(n)

SD=SE*sqrt(n)



Left.CI=beta-1.96*SD

Right.CI=beta+1.96*SD