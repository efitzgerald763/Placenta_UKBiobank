################################################ 
############ Not used until line 85 ############ 
################################################ 
library(tidyverse)

placenta <- fread("/data1/meaneylab/eamon/Placenta_ful_seq_geneID.csv")


placenta_eQTLs <- fread("/n9/cortexjobs/common/resources/expr_database/processed/expr.75/expr.75.txt")
placenta_snps_used <- fread("/n9/cortexjobs/UKbio_score/u_eprs/jobs/u163_eamon_placenta.rpkm.cyan22.ukbio.placenta_u163/Final_genotype/PRS_result/u163_prs.snplog.csv")%>%
  dplyr :: select(PRS_1.0)
placenta_genes <- fread("/n9/cortexjobs/UKbio_score/u_eprs/jobs/u163_eamon_placenta.rpkm.cyan22.ukbio.placenta_u163/logTABLE.csv")

# Filter first for genes,then eQTLs used in ePRS
placenta_trimmed <- filter(placenta_eQTLs, trait_id %in% placenta_genes$gtexENS_selectedForGWAS)
placenta_trimmed <- filter(placenta_trimmed, rsid %in% placenta_snps_used$PRS_1.0)

# Rename columns so they're recognized by TwosampleMR
colnames(placenta_trimmed) <- c('Study_ID', 'Feature', 'Variant_ID', 'Chr', 'Position', 'other_allele', 'effect_allele', 'SNP', 'direction', 
                                'beta', 'se', 'pval', 'sample_size', 'eaf', 'r_squared', 'platform_ID', 'Gene_name', 'cis_or_trans', 'q_value')

# Read in expression data and filter for the cyan module
Cyan <- fread("/data1/meaneylab/eamon/Placenta_ful_seq_geneID.csv") %>%
  dplyr::filter((moduleColors == "cyan"))

# Keep only genes in the PGS
Cyan <- filter(Cyan, V1 %in% placenta_trimmed$Feature)

# Calculate descriptive expression stats
mysummary2 <- function(x){
  if (is.numeric(x)){
    c(mean(x),median(x),sd(x), mad(x), IQR(x))
  }else{ print("pass")}
}
samp <- data.frame(colnames(Cyan))%>%
  as.data.frame()
stats <- data.frame(lapply(Cyan,mysummary2))

t_stats <- transpose(stats)%>%
  as.data.frame()
# get row and colnames in order
t_stats$PSCID <- samp$colnames.Cyan.
colnames(t_stats) <- c('mean', 'median', 'SD', 'Mad', 'IQR', 'PSCID')

# Remove first 2 and last 2 rows
t_stats <- slice(t_deets, 1:(n()-2))
t_stats <- t_deets[-(1:2), , drop = FALSE]

PRS_cyan <- fread("/n9/cortexjobs/c_eprs/jobs/c806_eamon_placenta.rpkm.cyan22.gusto.kids.placenta_c806/result/c806_prs.score.csv", header = T, stringsAsFactors = F) %>% 
  mutate(PRS_1.0 = as.numeric(scale(PRS_1.0, center = T, scale = T))) %>%
  mutate(ID1 = gsub("B", "", ID1))
colnames(PRS_cyan) <- c('PSCID', 'SNP_count_cyan', 'ePRS_cyan')

# Filter PGS dataframe for samples which have placenta seq data
samples <- colnames(placenta)
samples <- as.data.frame(samples)
PRS_cyan_Seq_samples <- filter(PRS_cyan, PSCID %in% samples$samples)
PRS_cyan_Seq_samples$ePRS_cyan_cut <- cut(PRS_cyan_Seq_samples$ePRS_cyan, 2, labels=FALSE)

PRS_cyan_Seq_samples <- PRS_cyan_Seq_samples %>%
  group_by(ePRS_cyan) %>% 
  mutate(median_split = c("above_median", "below_median")[1 + 
                                                            (ePRS_cyan <= median(ePRS_cyan))], 
         quartile_split = cut(ePRS_cyan, breaks = quantile(ePRS_cyan), 
                              labels = paste0(1:4, "_quartile")))

total <- full_join(t_stats, PRS_cyan_Seq_samples, by = "PSCID")
total$mean <- as.numeric(total$mean)
total$median_split <- as.character(total$median_split)
total$quartile_split <- as.character(total$quartile_split)
total$ePRS_cyan <- as.numeric(total$ePRS_cyan)



ggplot(total, aes(x=quartile_split, y=mean)) + 
  geom_boxplot(outlier.colour="red", outlier.shape=8,
               outlier.size=4)
library(ggpubr)
ggplot(total, aes(x=quartile_split, y=mean)) + 
  geom_boxplot(notch=TRUE)

library(ggpubr)
ggscatter(XX, x = 'ePRS_cyan', y = 'ssGSEA', 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "PGS", ylab = "ssGSEA")

#################################
########## GVEA #################
#################################
library(GSVA)
library(shinydashboard)
library(shinybusy)
library(tibble)
library(stringr) 
library(readxl)    


Modules <- fread("/data1/meaneylab/eamon/Placenta_study/Modules_GSVA.csv")

length(colnames(Modules))
nrow(Modules)


#Cyan_module <- fread("/data1/meaneylab/eamon/Placenta_study/placenta_rpkm_cyan22.txt")

Expr <- fread("/data1/meaneylab/eamon/Placenta_ful_seq_geneID.csv")

# Set rownames
Expr <- data.frame(column_to_rownames(Expr, var = "V1"))

# Remove "X" from sample names
colnames(Expr)<-gsub("X","",colnames(Expr))

# Remove surplus columns
Expr[1] <- NULL
Expr <- Expr[1:(length(Expr)-2)]

# Run analysis
ssGSEA_res <- gsva(as.matrix(Expr), Modules, method="ssgsea")%>%
  as.data.frame()
#rownames(ssGSEA_res) = "ssGSEA"
ssGSEA_res <- as.data.frame(t(ssGSEA_res))
ssGSEA_res <- tibble::rownames_to_column(ssGSEA_res, "PSCID")
ssGSEA_res$PSCID <- stringr::str_replace(ssGSEA_res$PSCID, '\\.', '-')


plage_res <- gsva(as.matrix(Expr), Modules, method="plage")%>%
  as.data.frame()
#rownames(plage_res) = "PLAGE"
plage_res <- as.data.frame(t(plage_res))
plage_res <- tibble::rownames_to_column(plage_res, "PSCID")
plage_res$PSCID <- stringr::str_replace(plage_res$PSCID, '\\.', '-')

zscore_res <- gsva(as.matrix(Expr), Modules, method="zscore")%>%
  as.data.frame()
#rownames(zscore_res) = "zscore"
zscore_res <- as.data.frame(t(zscore_res))
zscore_res <- tibble::rownames_to_column(zscore_res, "PSCID")
zscore_res$PSCID <- stringr::str_replace(zscore_res$PSCID, '\\.', '-')

GSVA_res <- gsva(as.matrix(Expr), Modules)%>%
  as.data.frame()
#rownames(GSVA_res) = "GSVA"
GSVA_res <- as.data.frame(t(GSVA_res))
GSVA_res <- tibble::rownames_to_column(GSVA_res, "PSCID")
GSVA_res$PSCID <- stringr::str_replace(GSVA_res$PSCID, '\\.', '-')

### Merge Results ###
#results <- cbind(ssGSEA_res, plage_res, zscore_res, GSVA_res)
#results <- tibble::rownames_to_column(results, "PSCID")

#results$PSCID <- stringr::str_replace(results$PSCID, '\\.', '-')

# Read in PRS
PRS_cyan <- fread("/n9/cortexjobs/c_eprs/jobs/c806_eamon_placenta.rpkm.cyan22.gusto.kids.placenta_c806/result/c806_prs.score.csv", header = T, stringsAsFactors = F) %>% 
  mutate(PRS_1.0 = as.numeric(scale(PRS_1.0, center = T, scale = T))) %>%
  mutate(ID1 = gsub("B", "", ID1))
colnames(PRS_cyan) <- c('PSCID', 'SNP_count_cyan', 'ePRS_cyan')

PRS_cyan_cortex <- fread("/n9/cortexjobs/c_eprs/jobs/c908_eamon_cyan.fetal.cortical.eqtl.gusto.kids.fetal_c908/scores_without_any_start_end_filter/result/c908_prs.score.csv", header = T, stringsAsFactors = F) %>% 
  mutate(PRS_1.0 = as.numeric(scale(PRS_1.0, center = T, scale = T))) %>%
  mutate(ID1 = gsub("B", "", ID1))
colnames(PRS_cyan_cortex) <- c('PSCID', 'SNP_count_cyan_cortex', 'ePRS_cyan_cortex')

PRS_random_2 <- fread("/n9/cortexjobs/c_eprs/jobs/c907_eamon_random2.placenta.eqtl.gusto.kids.placenta_c907/result/c907_prs.score.csv", header = T, stringsAsFactors = F) %>% 
  mutate(PRS_1.0 = as.numeric(scale(PRS_1.0, center = T, scale = T))) %>%
  mutate(ID1 = gsub("B", "", ID1))
colnames(PRS_random_2) <- c('PSCID', 'SNP_count_random_2', 'ePRS_random_2')

PRS_random_1 <- fread("/n9/cortexjobs/c_eprs/jobs/c906_eamon_random1.placenta.eqtl.gusto.kids.placenta_c906/result/c906_prs.score.csv", header = T, stringsAsFactors = F) %>% 
  mutate(PRS_1.0 = as.numeric(scale(PRS_1.0, center = T, scale = T))) %>%
  mutate(ID1 = gsub("B", "", ID1))
colnames(PRS_random_1) <- c('PSCID', 'SNP_count_random_1', 'ePRS_random_1')

PRS_light_cyan <- fread("/n9/cortexjobs/c_eprs/jobs/c905_eamon_lightcyan.placenta.eqtl.gusto.kids.placenta_c905/result/c905_prs.score.csv", header = T, stringsAsFactors = F) %>% 
  mutate(PRS_1.0 = as.numeric(scale(PRS_1.0, center = T, scale = T))) %>%
  mutate(ID1 = gsub("B", "", ID1))
colnames(PRS_light_cyan) <- c('PSCID', 'SNP_count_light_cyan', 'ePRS_light_cyan')

PRS_salmon <- fread("/n9/cortexjobs/c_eprs/jobs/c904_eamon_salmon.placenta.eqtl.gusto.kids.placenta_c904/result/c904_prs.score.csv", header = T, stringsAsFactors = F) %>% 
  mutate(PRS_1.0 = as.numeric(scale(PRS_1.0, center = T, scale = T))) %>%
  mutate(ID1 = gsub("B", "", ID1))
colnames(PRS_salmon) <- c('PSCID', 'SNP_count_salmon', 'ePRS_salmon')

PRS_midnight_blue <- fread("/n9/cortexjobs/c_eprs/jobs/c903_eamon_midnightblue.placenta.eqtl.gusto.kids.placenta_c903/result/c903_prs.score.csv", header = T, stringsAsFactors = F) %>% 
  mutate(PRS_1.0 = as.numeric(scale(PRS_1.0, center = T, scale = T))) %>%
  mutate(ID1 = gsub("B", "", ID1))
colnames(PRS_midnight_blue) <- c('PSCID', 'SNP_count_midnight_blue', 'ePRS_midnight_blue')

ePRS <- merge(PRS_cyan, PRS_cyan_cortex, by = "PSCID", all = T) %>%
  merge(      PRS_random_2, by = "PSCID", all = T) %>%
  merge(      PRS_random_1, by = "PSCID", all = T) %>%
  merge(      PRS_light_cyan, by = "PSCID", all = T) %>%
  merge(      PRS_salmon, by = "PSCID", all = T) %>%
  merge(      PRS_midnight_blue, by = "PSCID", all = T)
  

# Filter PRS dataframe for samples which have placenta seq data
#PSCID <- ssGSEA_res$PSCID
#PSCID <- as.data.frame(PSCID)

#colnames(Expr)<-gsub("X","",colnames(Expr))

PRS_cyan_Seq_samples <- full_join(PRS_cyan, ssGSEA_res, by = "PSCID")%>%
  na.omit
PRS_cyan_Seq_samples <- PRS_cyan_Seq_samples[,-1]

###### Make Corrplot ######
# Get p-values
corrp.mat <- cor_pmat(PRS_cyan_Seq_samples)
# Make correlation matrix
correlation_matrix <- round(cor(PRS_cyan_Seq_samples),2)

corrplot(correlation_matrix, type="full", tl.col="black", tl.srt=45,
         p.mat = corrp.mat, tl.cex = 0.7, sig.level = 0.05, insig = "blank",
         order = "hclust")

#Combine PRS_cyan_Seq_samples and results, then run regression
XX <- full_join(PRS_cyan_Seq_samples, ssGSEA_res, by = "PSCID")
XX <- XX[,-1]

###### Scatterplot ######  
library(ggpubr)
ggscatter(PRS_cyan_Seq_samples, x = 'ePRS_cyan', y = 'Grey60', 
          add = "none", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "Polygenic score", ylab = "ssGSEA")

###### Make Corrplot ######
# Get p-values
corrp.mat <- cor_pmat(XX)
# Make correlation matrix
correlation_matrix <- round(cor(XX),2)

corrplot(correlation_matrix, type="upper", tl.col="black", tl.srt=45,
         p.mat = corrp.mat, sig.level = 0.05, insig = "blank")

###########################################################################
### Birth metrics
BMI <- read_spss("/n0/gusto/DATA/BMI/GUSTO - ZBMI all 06Aug2019 - restructured - post_ageQC.sav")%>%
  dplyr::rename(PSCID = SubjectID)

### Demographics
demo <- read_spss("/n0/gusto/DATA/GUSTO for demog 23Aug2019.sav") %>%
  dplyr :: select(PSCID, mother_age_delivery, mother_highest_education_2cat, household_income_2cat) 

### Sex
sex <- openxlsx::read.xlsx("/share/projects/gusto/DATA/FORMS/FORMA_(340)_CS_29-SEPT-2018UPDATE/compiled data/FormA340_20181013.xlsx", sheet = 1, check.names = T) %>%
  dplyr :: select(SubjectID, sex) %>%
  dplyr :: rename(PSCID = SubjectID)


### Genetic PCs (kids)
PC_kid <- fread("/share/projects/genodata/GUSTO/toplot.child_allsamp_PCs.txt")

#### STAI and EPDS 26 weeks
STAI_EPDS <- readxl::read_excel ("/n0/gusto/DATA/FORMS/FORMA_(219)_CS_25-SEPT-2018UPDATE/clean dataset/FormA219_20180925 contains EPD, STAI.xlsx") %>%
  dplyr :: rename(PSCID = SubjectID)

### Make numeric 
STAI_EPDS$EPDS_imputed_pw26 <-as.numeric(STAI_EPDS$EPDS_imputed_pw26)

### Combine data
XX <- full_join(STAI_EPDS, PC_kid, by = "PSCID") %>%
  full_join(ePRS,         by = "PSCID") %>%
  full_join(ssGSEA_res,         by = "PSCID") %>%
  full_join(demo,         by = "PSCID") %>%
  full_join(sex,         by = "PSCID") %>%
  full_join(BMI,         by = "PSCID") %>%
  dplyr:: filter(substring(PSCID, 1, 3) %in% c("010", "020"))  ## only keep singletons


## Sex 1=male 2=female
XX$Sex <- NA
XX$Sex[which(XX$sex.x == "Male")] <- "1"
XX$Sex[which(XX$sex.x == "Female")] <- "2"

######### Multile Regression ######### 
rownames(XX) <- XX$PSCID
XX.cont <- XX
lm.cont.old <- lm(Cyan ~  PC1 + PC2 + PC3 + GA + ePRS_cyan, data = XX.cont, na.action=na.exclude)
summary(lm.cont.old)

nrow_coef <- nrow(summary(lm.cont.old)$coef)
im.cont <- data.frame(influence.measures(lm.cont.old)$is.inf)
case.exclude.cont <- rownames(im.cont[which(im.cont[,nrow_coef] == "TRUE"),])
if (length(case.exclude.cont)>0){XX.cont <- XX[-which(XX$UniqueID %in% case.exclude.cont),]}

lm.cont     <- lm((Cyan + 1) ~  PC1 + PC2 + PC3 + GA + ePRS_cyan, data = XX.cont, na.action=na.exclude)
summary(lm.cont)
confint(lm.cont, 'ePRS_cyan', level=0.95)

## diagnostics, check assumptions ##
par(mfrow = c(2, 2))
plot(lm.cont)
par(mfrow=c(1,1))
boxCox(lm.cont,lambda=c(-2,2,0.5))

########################################################

# List all column names except PSCID
list.Y <- colnames(ssGSEA_res)[colnames(ssGSEA_res) != "PSCID"] 


list.Y.names <- list.Y

list.ePRS <- colnames(ePRS)[which(grepl("PRS", colnames(ePRS)))]
list.ePRS.names <- list.ePRS
######################



library("questionr")
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

#sink(paste0("~/Results_", Sys.Date() ,".txt",sep=""))
XX = as.data.frame(XX)
rownames(XX) <- XX$PSCID
iID <- 0

for (k in 1:length(list.ePRS)){  
  XX$ePRS.cur <-eval(parse(text=paste("XX$", list.ePRS[k], sep = "")))
  summary(XX$X)
  
  for (i in 1:length(list.Y)){  
    XX$Y <- eval(parse(text=paste0("XX$", list.Y[i])))
    
    
    XX.cont <- XX[complete.cases(XX[, c("Y", "PC1", "ePRS.cur", "Sex")]), ]
    
    
    lm.cont.old <- lm(Y ~ Sex + PC1 + PC2 + PC3 + ePRS.cur, data = XX.cont,na.action=na.exclude)
    
    nrow_coef <- nrow(summary(lm.cont.old)$coef)
    
    im.cont <- data.frame(influence.measures(lm.cont.old)$is.inf)
    case.exclude <- rownames(im.cont[which(im.cont$dffit == "TRUE"),])
    if (length(case.exclude)>0){XX.cont <- XX.cont[-which(XX.cont$PSCID %in% case.exclude),]}
    
    
    lm.cont     <- lm(Y ~ Sex + PC1 + PC2 + PC3 + ePRS.cur, data = XX.cont, na.action=na.exclude)
    expr
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

###########

#Adjust p-values for multiple comparisons
PP = PP[order(PP$p_value),]

PP$Bonferroni =
  p.adjust(PP$p_value,
           method = "bonferroni")

PP$BH =
  p.adjust(PP$p_value,
           method = "BH")



write_csv(PP, "/data1/meaneylab/eamon/Placenta_study/NegControl_ssGSEA_GUSTO_reg_PCs_sex_26Jan2022.csv")

