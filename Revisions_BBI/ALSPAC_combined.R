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
c782 <- fread("/n9/cortexjobs/c_eprs/jobs/c782_eamon_inflammatory_fetalMicroglia.alspackids.fetaleQTLs_c782/scores_without_any_start_end_filter/result/c782_prs.score.csv", header = T, stringsAsFactors = F) %>% 
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
#                                     Estimate Std. Error z value Pr(>|z|) 
# ePRS_microglia:EPDS_mean_prenatal   0.04004    0.01818   2.202  0.02767 * 
0.04004 - (0.01818*1.96) 
0.04004 + (0.01818*1.96) 

table(glm.cont$model$Any_dep)
#0    1 
#1273  179 

## simple slope analysis and plot (when there is a significant interaction)
sim_slopes(glm.cont, pred = "EPDS_mean_prenatal", modx = "ePRS_microglia", johnson_neyman = TRUE)$slopes

#Value of ePRS_microglia       Est.       S.E.        2.5%      97.5%    z val.            p
#1             -1.02099154 0.00715244 0.02643362 -0.04465650 0.05896138 0.2705812 0.7867131636
#2              0.02581637 0.04906592 0.01880765  0.01220361 0.08592823 2.6088285 0.0090852767
#3              1.07262427 0.09097940 0.02708079  0.03790204 0.14405677 3.3595554 0.0007806801

# Adjusted interact_plot code
plotme <- interact_plot(glm.cont, pred = "EPDS_mean_prenatal", modx = "ePRS_microglia", 
                      x.label = "Prenatal EPDS", 
                      y.label = "Depression (Probability)",
                      vary.lty = FALSE,
                      line.thickness = 0.75,
                      main.title = "CIS-R Depression (Females)",
                      legend.main = "Microglia PGS", colors = c("tomato3", "grey80", "#3399FF")) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 16),
        text = element_text(size = 16), 
        legend.position = 'right',
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) 


ggsave("/data1/meaneylab/eamon/Microglia_ePRS/Revisions_BBI/ALSPAC_Combined.pdf", 
       plot = plotme, width = 15, height = 10, units = "cm", bg="white")


#----- Any Dep 24 years Males --------
XX_male <- subset(XX, Gender =="Male")

rownames(XX_male) <- XX_male$UniqueID
XX.cont <- XX_male

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
# ePRS_microglia:EPDS_mean_prenatal  -0.027427   0.034306  -0.799    0.424

length(glm.cont$residuals)
#925

table(XX.cont$Any_dep)
#0    1 
#1174   92 


