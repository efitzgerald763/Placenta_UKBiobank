
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

#################### where to find files ####################
## GUSTO pheno folder: /share/projects/gusto/DATA/

## GUSTO geno folder: /share/projects/genodata/GUSTO/

## (e)PRS scores folder: /share/projects/cortexjobs
## (e)PRS scores's info: /share/projects/cortexjobs/jobs_index.csv


#####################################################################
########## PART I: LOAD, MERGE, AND PROCESS DATA SETS ###############
#####################################################################


#Birth metrics
BMI <- read_spss("/n0/gusto/DATA/BMI/GUSTO - ZBMI all 06Aug2019 - restructured - post_ageQC.sav")%>%
  dplyr::rename(PSCID = SubjectID)

#Demographics
demo <- read_spss("/n0/gusto/DATA/GUSTO for demog 23Aug2019.sav") %>%
  dplyr :: select(PSCID, mother_highest_education_2cat, household_income_2cat) 


#Delivery info
Delivery_info <- openxlsx::read.xlsx("/data1/meaneylab/eamon/FORMA_(332)_CS_28-DEC-2016_DeliveryGUSTOCaseReportForm.xlsx", sheet = 1)%>%
  dplyr :: select(1, 32, 47, 30, 117, 118, 62,59) %>%
  dplyr::rename(PSCID = 1) %>%
  dplyr::rename(Chorio = "32") %>%
  dplyr::rename (Maternal_steroids = "47") %>%
  dplyr::rename (Ruptured_membranes = "30") %>%
  dplyr::rename (APGAR_1_min = "117") %>%
  dplyr::rename (APGAR_5_min = "118") %>%
  dplyr::rename (Delivery_mode = "59") %>%
  dplyr::rename (Induced_Spon = "62")

Delivery_info = Delivery_info[-1,]

Delivery_info$Induced_Spon <- gsub("1_Spontaneous","1",Delivery_info$Induced_Spon)
Delivery_info$Induced_Spon <- gsub("2_Induced","2",Delivery_info$Induced_Spon)
Delivery_info$Induced_Spon <- gsub("3_na","0",Delivery_info$Induced_Spon)
Delivery_info$Induced_Spon <- gsub("not_answered","0",Delivery_info$Induced_Spon)
Delivery_info$Induced_Spon <- gsub("not_applicable","0",Delivery_info$Induced_Spon)

Delivery_info$Induced_Spon <- as.numeric(Delivery_info$Induced_Spon)

table(Delivery_info$Delivery_mode)
Delivery_info$Delivery_mode <- gsub("1_Normal_vaginal_delivery","1",Delivery_info$Delivery_mode)
Delivery_info$Delivery_mode <- gsub("2_Emergency_LSCS","2",Delivery_info$Delivery_mode)
Delivery_info$Delivery_mode <- gsub("3_Elective_LSCS","2",Delivery_info$Delivery_mode)
Delivery_info$Delivery_mode <- gsub("5_Vacuum","NA",Delivery_info$Delivery_mode)
Delivery_info$Delivery_mode <- gsub("6_Forceps","NA",Delivery_info$Delivery_mode)
Delivery_info$Delivery_mode <- gsub("not_answered","NA",Delivery_info$Delivery_mode)
Delivery_info$Delivery_mode <- gsub("qn_13_mode_delivery","NA",Delivery_info$Delivery_mode)

Delivery_info$Delivery_mode <- as.numeric(Delivery_info$Delivery_mode)

#Sex
sex <- openxlsx::read.xlsx("/share/projects/gusto/DATA/FORMS/FORMA_(340)_CS_29-SEPT-2018UPDATE/compiled data/FormA340_20181013.xlsx", sheet = 1, check.names = T) %>%
  dplyr :: select(SubjectID, sex) %>%
  dplyr :: rename(PSCID = SubjectID)


##########################
## get PCs (kids)
##########################
PC_kid <- fread("/share/projects/genodata/GUSTO/toplot.child_allsamp_PCs.txt")
PC_mom <- fread("/share/projects/genodata/GUSTO/toplot.mother_allsamp_PCs.txt")



###### PRSs######

##### standardize and rename (e)PRSes #####
PRS_cyan <- fread("/n9/cortexjobs/c_eprs/jobs/c806_eamon_placenta.rpkm.cyan22.gusto.kids.placenta_c806/result/c806_prs.score.csv", header = T, stringsAsFactors = F) %>% 
  mutate(PRS_1.0 = as.numeric(scale(PRS_1.0, center = T, scale = T))) %>%
  mutate(ID1 = gsub("B", "", ID1))

colnames(PRS_cyan) <- c('PSCID', 'SNP_count_cyan', 'ePRS_cyan')


ePRS <- PRS_cyan 


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
  full_join(demo,         by = "PSCID") %>%
  full_join(ABC,         by = "PSCID") %>%
  full_join(STAI_EPDS,         by = "PSCID") %>%
  full_join(sex,          by = "PSCID") %>%
  dplyr:: filter(substring(PSCID, 1, 3) %in% c("010", "020"))  ## only keep singletons


## Gender: 1=male 2=female
XX$Sex <- NA
XX$Sex[which(XX$sex.x == "Male")] <- "1"
XX$Sex[which(XX$sex.x == "Female")] <- "2"

x <- XX %>%
  dplyr::select(ePRS_cyan, 
                EPDS_imputed_pw26,
                STAI_imputed_total_pw26,
                STAI_imputed_state_pw26,
                STAI_imputed_trait_pw26,
                parity,
                Induced_Spon,
                Delivery_mode,
                mother_age_delivery, 
                GA,
                weight.day1, 
                length.day1, 
                zbmi.day1,
                mother_highest_education_2cat,
                household_income_2cat,
                Sex) %>%
  na.omit()
colnames(x) <- c('Fetoplacental PGS',  
                 'EPDS PCW26', 
                 'STAI total PCW26',
                 'STAI state PCW26',
                 'STAI trait PCW26',
                 'Parity',
                 'Induced labor',
                 'Delivery mode',
                 'Mother age at delivery', 
                 'Gestational age at birth', 
                 'Birth weight', 
                 'Birth length',
                 'Birth BMI (z-score)', 
                 'Mother education', 
                 'Household income', 
                 'Sex')

x$Sex <-as.numeric(x$Sex)

###### Make Corrplot ######
# Get p-values
corrp.mat <- cor_pmat(x)
# Make correlation matrix
correlation_matrix <- round(cor(x),2)

pdf(file = "/data1/meaneylab/eamon/Placenta_study/Revisions/PGS_enviroment_corrplot.pdf")
corrplot(correlation_matrix, type="full", tl.col="black", tl.srt=45,
         p.mat = corrp.mat, sig.level = 0.05, insig = "blank", 
         number.cex=0.5)
dev.off()

