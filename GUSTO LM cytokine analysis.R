## set working directory
setwd("")
getwd()
library(haven)
library(dplyr)
library(interactions)
library(data.table)
library(jtools)
library(xlsx)
library(readxl)
library(openxlsx)
library(car)
library(ggplot2)
library(readxl)
library(stargazer)
#################### where to find files ####################
## GUSTO pheno folder: /share/projects/gusto/DATA/

## GUSTO geno folder: /share/projects/genodata/GUSTO/

## (e)PRS scores folder: /share/projects/cortexjobs
## (e)PRS scores's info: /share/projects/cortexjobs/jobs_index.csv


#####################################################################
########## PART I: LOAD, MERGE, AND PROCESS DATA SETS ###############
#####################################################################

#Demographics
demo <- read_spss("/n0/gusto/DATA/GUSTO for demog 23Aug2019.sav") %>%
  dplyr :: select(PSCID, mother_age_delivery, mother_highest_education_2cat, household_income_2cat) 

#Birth metrics
BMI <- read_spss("/n0/gusto/DATA/BMI/GUSTO - ZBMI all 06Aug2019 - restructured - post_ageQC.sav")%>%
  dplyr::rename(PSCID = SubjectID) %>%
  dplyr :: select(PSCID, zbmi.day1, GA, weight.day1) 


#####################################
###### Get cord cytokine data #######
#####################################
# Remove adiponectin dilutions, leaving only the final data
ELISA <- openxlsx::read.xlsx("/data1/meaneylab/eamon/Placenta_study/Cord cytokines/ELISA/SIgN-ELISA_data.xlsx")%>%
  dplyr::select(-c(`Adipo_ngml.(2000X.dilution)`,`Adipo_ngml.(neat)`))

# Remove data with CV greater than 20% between plates
Luminex <- openxlsx::read.xlsx("/data1/meaneylab/eamon/Placenta_study/Cord cytokines/Luminex/SIgN-Luminex_data.xlsx")%>%
  dplyr::select(-c(`Eotaxin.(pg/ml)`,`Ghrelin.(pg/ml)`, `IL-2.(pg/ml)`, `IL-7.(pg/ml)`, `MIG.(pg/ml)`, `MIP-3alpha.(pg/ml)`,
                   `T4.(pg/ml)`,`Estradiol.(pg/ml)`, `Cortisol.(pg/ml)`, `IGFBP-4.(pg/ml)`, `PIGF-1.(pg/ml)`, `PAI-1.(pg/ml)`,
                   `TGF-beta1.(pg/ml)`, `hIgG4.(pg/ml)`, `hIgG3.(pg/ml)`, `hIgG2.(pg/ml)`, `hIgG1.(pg/ml)`, `IgM.(pg/ml)`,
                   `IgA.(pg/ml)`))

# Remove data with CV greater than 20% between plates
Quanterix <- openxlsx::read.xlsx("/data1/meaneylab/eamon/Placenta_study/Cord cytokines/Quanterix/Quanterix_cytokines_data.xlsx")%>%
  dplyr::rename(PSCID = SubjectID)%>%
  dplyr::select(-c("IFNg","IL-4"))

############# Use na.omit() if running corrplot ############# 
Cord_cyto <- full_join(ELISA, Luminex, by = "PSCID") %>%
  full_join(Quanterix,         by = "PSCID")
colnames(Cord_cyto) <- c('PSCID', 'Free_testosterone', 'Total_testosterone', 'Adiponectin', 'C_peptide',
                         'IL1RA', 'IP10', 'MCP1', 'MIP1a', 'MIP1b', 'VEGFA', 'GLP1', 'IL12p40', 'Leptin',
                         'Glucagon', 'Growth_hormone', 'IGFBP3', 'IGFBP7', 'Insulin', 'Prolactin', 'FSH',
                         'LH', 'TSH', 'IgE', 'CRP', 'IL10', 'IL6', 'TNFa')
Cord_cyto[, -c(1)] <- scale(Cord_cyto[, -c(1)])

###### Make Corrplot ######
# Remove the PSCID column
Cord_cyto <- Cord_cyto[,-1]

# Get p-values
corrp.mat <- cor_pmat(Cord_cyto)

# Make correlation matrix
correlation_matrix <- round(cor(Cord_cyto),2)

# Plot
corrplot(correlation_matrix, type="full", tl.col="black", tl.srt=45, tl.cex = 0.7,
         p.mat = corrp.mat, sig.level = 0.00006858, insig = "blank", order="hclust")
dev.off()
#################################################################
###### 7YO cytokines- scale so the betas are comparable ########
#################################################################

Cytokines <- read_spss("/n0/gusto/DATA/Samples/GUSTO cytokines CRP.sav")%>%
  dplyr::rename(PSCID = ID)%>%
  mutate(MeanpgmL.IL8 = as.numeric(scale(MeanpgmL.IL8, center = T, scale = T)))%>%
  mutate(MeanpgmL.IL1B = as.numeric(scale(MeanpgmL.IL1B, center = T, scale = T)))%>%
  mutate(MeanpgmL.IL6 = as.numeric(scale(MeanpgmL.IL6, center = T, scale = T)))%>%
  mutate(MeanpgmL.TNFa = as.numeric(scale(MeanpgmL.TNFa, center = T, scale = T)))%>%
  mutate(MeanpgmL.CRP = as.numeric(scale(MeanpgmL.CRP, center = T, scale = T)))

### Remove individuals that reported any sickness on day of sampling ###

Unwell_at_collection <- read_xlsx("/n0/gusto/DATA/Samples/7Y cytokine collection remarks_Medications_Unwell.xlsx")

TNF <- TNF[Unwell_at_collection$`Child Unwell (Y/N:what symptoms)` == 'N',]
IL1B <- IL1B[Unwell_at_collection$`Child Unwell (Y/N:what symptoms)` == 'N',]
IL8 <- IL8[Unwell_at_collection$`Child Unwell (Y/N:what symptoms)` == 'N',]
CRP <- CRP[Unwell_at_collection$`Child Unwell (Y/N:what symptoms)` == 'N',]
IL6 <- IL6[Unwell_at_collection$`Child Unwell (Y/N:what symptoms)` == 'N',]


#Delivery info
Delivery_info <- openxlsx::read.xlsx("/data1/meaneylab/eamon/FORMA_(332)_CS_28-DEC-2016_DeliveryGUSTOCaseReportForm.xlsx", sheet = 1)%>%
  dplyr :: select(1, 32, 47, 30, 117, 118) %>%
  dplyr::rename(PSCID = 1) %>%
  dplyr::rename(Chorio = "32") %>%
  dplyr::rename (Maternal_steroids = "47") %>%
  dplyr::rename (Ruptured_membranes = "30") %>%
  dplyr::rename (APGAR_1_min = "117") %>%
  dplyr::rename (APGAR_5_min = "118") 

Delivery_info = Delivery_info[-1,]

#Sex
sex <- openxlsx::read.xlsx("/share/projects/gusto/DATA/FORMS/FORMA_(340)_CS_29-SEPT-2018UPDATE/compiled data/FormA340_20181013.xlsx", sheet = 1, check.names = T) %>%
  dplyr :: select(SubjectID, sex) %>%
  dplyr :: rename(PSCID = SubjectID)


##########################
## get PCs (kids)
##########################
PC_kid <- fread("/share/projects/genodata/GUSTO/toplot.child_allsamp_PCs.txt")
# PC_mom <- fread("/share/projects/genodata/GUSTO/toplot.mother_allsamp_PCs.txt")



###### PRSs######

##### standardize and rename (e)PRSes #####

PRS_cyan <- fread("/n9/cortexjobs/c_eprs/jobs/c806_eamon_placenta.rpkm.cyan22.gusto.kids.placenta_c806/result/c806_prs.score.csv", header = T, stringsAsFactors = F) %>% 
  mutate(PRS_1.0 = as.numeric(scale(PRS_1.0, center = T, scale = T))) %>%
  mutate(ID1 = gsub("B", "", ID1))

colnames(PRS_cyan) <- c('PSCID', 'SNP_count_cyan', 'ePRS_cyan')



#################

ePRS <- PRS_cyan 


XX <- full_join(ePRS, PC_kid, by = "PSCID") %>%
  full_join(Cord_cyto,         by = "PSCID") %>%
  full_join(BMI,         by = "PSCID") %>%
  full_join(demo,         by = "PSCID") %>%
  full_join(Delivery_info,         by = "PSCID") %>%
  full_join(sex,          by = "PSCID") %>%
  dplyr:: filter(substring(PSCID, 1, 3) %in% c("010", "020"))  ## only keep singletons

#Remove samples who had chorioamnionitis or maternal steroids
XX <- XX %>% dplyr::filter ((Chorio == "0_no") & (Maternal_steroids == "1_No")  & (Ruptured_membranes == "0_no"))  

## Gender: 1=male 2=female
XX$Gender <- NA
XX$Gender[which(XX$sex == "Male")] <- "1"
XX$Gender[which(XX$sex == "Female")] <- "2"


# Split males and females
XX_female <- subset(XX, Gender =="2")
XX_male <- subset(XX, Gender =="1")


XX <- XX_female
#####################################################################
########## PART II: EXAMPLE ONE INTERACTION/REGRESSION   ############
#####################################################################
# Add to birth regression mother_age_delivery + mother_highest_education_2cat + household_income_2cat + zbmi.day1 + GA
# Add to 7 yo regression mother_age_delivery + mother_highest_education_2cat + household_income_2cat
# For cord cytokines run with PCs then introduce co-variates
list.Y

XX = as.data.frame(XX)
rownames(XX) <- XX$PSCID
XX.cont <- XX
lm.cont.old <- lm(IP10 ~  Gender + mother_age_delivery + household_income_2cat + 
                    GA + PC1 + PC2 + PC3 + ePRS_cyan, data = XX.cont, na.action=na.exclude)

summary(lm.cont.old)

nrow_coef <- nrow(summary(lm.cont.old)$coef)
im.cont <- data.frame(influence.measures(lm.cont.old)$is.inf)
case.exclude.cont <- rownames(im.cont[which(im.cont[,"dffit"] == "TRUE"),])
if (length(case.exclude.cont)>0){XX.cont <- XX[-which(XX$PSCID %in% case.exclude.cont),]}

lm.cont  <- lm(IP10 ~  Gender + mother_age_delivery + household_income_2cat + 
                 GA + PC1 + PC2 + PC3 + ePRS_cyan, data = XX.cont, na.action=na.exclude)

summary(lm.cont)

## diagnostics, check assumptions ##
par(mfrow = c(2, 2))
plot(lm.cont)
par(mfrow=c(1,1))
boxCox(lm.cont,lambda=c(-2,2,0.5))

lm.cont <- lm (log (IP10 + 1) ~  Gender + mother_age_delivery + household_income_2cat + 
                 GA + PC1 + PC2 + PC3 + ePRS_cyan, data = XX.cont, na.action=na.exclude)
summary(lm.cont)

######################################################################
####### PART III: EXAMPLE 2 LOOP ALL INTERACTIONS/REGRESSIONS ########
######################################################################

# List all column names except PSCID
list.Y <- colnames(Cord_cyto)[colnames(Cord_cyto) != "PSCID"] 

list.Y.names <- list.Y

list.ePRS <- colnames(PRS_cyan)[which(grepl("PRS", colnames(PRS_cyan)))]
list.ePRS.names <- list.ePRS
######################

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

sink(paste0("~/Cord_cytokine_results_", Sys.Date() ,".txt",sep=""))

XX = as.data.frame(XX)
rownames(XX) <- XX$PSCID
iID <- 0

for (k in 1:length(list.ePRS)){  
  XX$ePRS.cur <-eval(parse(text=paste("XX$", list.ePRS[k], sep = "")))
  summary(XX$X)
  
  for (i in 1:length(list.Y)){  
    XX$Y <- eval(parse(text=paste0("XX$", list.Y[i])))
    
      
      XX.cont <- XX[complete.cases(XX[, c("Y", "GA", "Gender", "PC1","ePRS.cur")]), ]
      
      
      lm.cont.old <- lm(Y ~ Gender + mother_age_delivery + household_income_2cat + 
                          GA + PC1 + PC2 + PC3 + ePRS.cur, data = XX.cont,na.action=na.exclude)
      
      nrow_coef <- nrow(summary(lm.cont.old)$coef)
      
      im.cont <- data.frame(influence.measures(lm.cont.old)$is.inf)
      case.exclude <- rownames(im.cont[which(im.cont$dffit == "TRUE"),])
      if (length(case.exclude)>0){XX.cont <- XX.cont[-which(XX.cont$PSCID %in% case.exclude),]}
      
      
      lm.cont     <- lm(log(Y + 1) ~ Gender + mother_age_delivery + household_income_2cat + 
                          GA + PC1 + PC2 + PC3 + ePRS.cur, data = XX.cont, na.action=na.exclude)
      
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

warnings()


write_csv(PP, "/data1/meaneylab/eamon/Placenta_study/Cord_cytokine_results_31Jan2022.csv")
