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
library(ggcorrplot)
library(corrplot)
library(RColorBrewer)


placenta_seq_sample_ids <- fread("/data1/meaneylab/eamon/Placenta_study/sample inflammatory markers.txt") %>%
  dplyr :: select(subjectid) %>%
  dplyr::rename(PSCID = subjectid)


#Birth metrics
BMI <- read_spss("/n0/gusto/DATA/BMI/GUSTO - ZBMI all 06Aug2019 - restructured - post_ageQC.sav")%>%
  dplyr::rename(PSCID = SubjectID)


#Cord cytokines- scale so the betas are comparable
Cytokines <- read_spss("/n0/gusto/DATA/Samples/GUSTO cytokines CRP.sav")%>%
  dplyr::rename(PSCID = ID)%>%
  mutate(MeanpgmL.IL8 = as.numeric(scale(MeanpgmL.IL8, center = T, scale = T)))%>%
  mutate(MeanpgmL.IL1B = as.numeric(scale(MeanpgmL.IL1B, center = T, scale = T)))%>%
  mutate(MeanpgmL.IL6 = as.numeric(scale(MeanpgmL.IL6, center = T, scale = T)))%>%
  mutate(MeanpgmL.TNFa = as.numeric(scale(MeanpgmL.TNFa, center = T, scale = T)))%>%
  mutate(MeanpgmL.CRP = as.numeric(scale(MeanpgmL.CRP, center = T, scale = T)))
  

#Delivery info
Delivery_info <- openxlsx::read.xlsx("/data1/meaneylab/eamon/FORMA_(332)_CS_28-DEC-2016_DeliveryGUSTOCaseReportForm.xlsx", sheet = 1)%>%
  dplyr :: select(1, 32, 47, 30, 117, 118, 62) %>%
  dplyr::rename(PSCID = 1) %>%
  dplyr::rename(Chorio = "32") %>%
  dplyr::rename (Maternal_steroids = "47") %>%
  dplyr::rename (Ruptured_membranes = "30") %>%
  dplyr::rename (APGAR_1_min = "117") %>%
  dplyr::rename (APGAR_5_min = "118") %>%
    dplyr::rename (Induced_Spon = "62")

Delivery_info = Delivery_info[-1,]

Delivery_info$Induced_Spon <- gsub("1_Spontaneous","1",Delivery_info$Induced_Spon)
Delivery_info$Induced_Spon <- gsub("2_Induced","2",Delivery_info$Induced_Spon)
Delivery_info$Induced_Spon <- gsub("3_na","0",Delivery_info$Induced_Spon)
Delivery_info$Induced_Spon <- gsub("not_answered","0",Delivery_info$Induced_Spon)
Delivery_info$Induced_Spon <- gsub("not_applicable","0",Delivery_info$Induced_Spon)

Delivery_info$Induced_Spon <- as.numeric(Delivery_info$Induced_Spon)
str(Delivery_info$Induced_Spon)
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
PRS_0.05 <- fread("/n9/cortexjobs/c_eprs/jobs/c832_eamon_placenta.modules.fdr.0.05.gusto.kids.placenta_c832/result/c832_prs.score.csv", header = T, stringsAsFactors = F) %>% 
 mutate(PRS_1.0 = as.numeric(scale(PRS_1.0, center = T, scale = T))) %>%
  mutate(ID1 = gsub("B", "", ID1))
PRS_0.1 <- fread("/n9/cortexjobs/c_eprs/jobs/c831_eamon_placenta.modules.fdr.0.1.gusto.kids.placenta_c831/result/c831_prs.score.csv", header = T, stringsAsFactors = F) %>% 
  mutate(PRS_1.0 = as.numeric(scale(PRS_1.0, center = T, scale = T))) %>%
  mutate(ID1 = gsub("B", "", ID1))
PRS_cyan <- fread("/n9/cortexjobs/c_eprs/jobs/c806_eamon_placenta.rpkm.cyan22.gusto.kids.placenta_c806/result/c806_prs.score.csv", header = T, stringsAsFactors = F) %>% 
  mutate(PRS_1.0 = as.numeric(scale(PRS_1.0, center = T, scale = T))) %>%
  mutate(ID1 = gsub("B", "", ID1))

colnames(PRS_0.05) <- c('PSCID', 'SNP_count_0.05', 'ePRS_0.05')
colnames(PRS_0.1) <- c('PSCID', 'SNP_count_0.1', 'ePRS_0.1')
colnames(PRS_cyan) <- c('PSCID', 'SNP_count_cyan', 'ePRS_cyan')



#################
ePRS <- merge(PRS_0.05, PRS_0.1, by = "PSCID", all = T) %>%
  merge(PRS_cyan, by = "PSCID", all = T)
#ePRS <- PRS_cyan 


####### STAI and EPDS 26 weeks#####
STAI_EPDS <- readxl::read_excel ("/n0/gusto/DATA/FORMS/FORMA_(219)_CS_25-SEPT-2018UPDATE/clean dataset/FormA219_20180925 contains EPD, STAI.xlsx") %>%
  dplyr :: rename(PSCID = SubjectID)

### Make the following numeric variables (as opposed to string), check using str()
STAI_EPDS$EPDS_imputed_pw26 <-as.numeric(STAI_EPDS$EPDS_imputed_pw26)

## ABC ## description available in /projects/gusto/Analysis/ABC
ABC <- read_spss("/share/projects/gusto/Analysis/ABC/GUSTO A,B,C 27Sep2018.sav")


XX <- full_join(BMI, PC_kid, by = "PSCID") %>%
  full_join(ePRS,         by = "PSCID") %>%
  full_join(Delivery_info,         by = "PSCID") %>%
  full_join(Cytokines,         by = "PSCID") %>%
  full_join(ABC,         by = "PSCID") %>%
  full_join(STAI_EPDS,         by = "PSCID") %>%
  full_join(sex,          by = "PSCID") %>%
  dplyr:: filter(substring(PSCID, 1, 3) %in% c("010", "020"))  ## only keep singletons

placenta_trimmed <- dplyr :: filter(XX, PSCID %in% placenta_seq_sample_ids$PSCID)

table(placenta_trimmed$Induced_Spon)

## Gender: 1=male 2=female
placenta_trimmed$Sex <- NA
placenta_trimmed$Sex[which(placenta_trimmed$sex.x == "Male")] <- "1"
placenta_trimmed$Sex[which(placenta_trimmed$sex.x == "Female")] <- "2"

x <- placenta_trimmed %>%
  dplyr::select(EPDS_imputed_pw26, mother_age_delivery, parity, 
                GA,weight.day1, length.day1, zbmi.day1, PC1, PC2, PC3, 
                SNP_count_cyan, ePRS_cyan, Sex, Induced_Spon) %>%
  na.omit()
colnames(x) <- c('Prenatal EPDS', 'Mother age at delivery', 'Parity',
                 'Gestational age at birth', 'Birth weight', 'Birth length',
                 'Birth BMI (z-score)', 'PC1', 'PC2', 'PC3', 'SNP_Count', 'Fetoplacent', 'Sex', 'Induced or Spontaneous labor')

x$Sex <-as.numeric(x$Sex)

###### Make Corrplot ######
# Get p-values
corrp.mat <- cor_pmat(x)
corrp.mat <- column_to_rownames(corrp.mat, var = 'rowname')%>%
  as.matrix()
# Make correlation matrix
correlation_matrix <- round(cor(x),2)

corrplot(correlation_matrix, type="upper", tl.col="black", tl.srt=45,
         p.mat = corrp.mat, sig.level = 0.05, insig = "blank")

