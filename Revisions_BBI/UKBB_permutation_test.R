library(dplyr)
library(tidyverse)
library(data.table)
library(tidyr)
library(interactions)

key <- fread("/share/projects/uk_biobank/UKBioBank_Key_all_02May2022.csv", header = T, stringsAsFactors = F) %>%
  filter(field.showcase %in% c(189, 21003, 21022, 20016, 20127, 22009, 21001, 22027, 22021, 22000, 22006, 54,34,52, 
                               20400, 21000, 22001, 20544, 20483, 25005, 25886, 25887,
                               21001, 20022, 20116, 26410, 26427, 26426, 21000, 21022, 30710, 884,
                               20525, 6141, 709, 6160))



XX  <- fread("/n0/uk_biobank/UKB_pheno_part1234_502506_subjects_20Feb2020.csv", header = T, stringsAsFactors = F,
             select = unique(c("eid", "sex_f31_0_0",
                               key$col.name,
                               "able_to_pay_rentmortgage_as_an_adult_f20525_0_0.x",
                               "sexual_interference_by_partner_or_expartner_without_consent_as_an_adult_f20524_0_0.x",
                               "physical_violence_by_partner_or_expartner_as_an_adult_f20523_0_0.x",
                               "been_in_a_confiding_relationship_as_an_adult_f20522_0_0.x",
                               "belittlement_by_partner_or_expartner_as_an_adult_f20521_0_0.x",
                               "been_in_serious_accident_believed_to_be_lifethreatening_f20526_0_0.x",
                               "been_involved_in_combat_or_exposed_to_warzone_f20527_0_0.x",
                               "diagnosed_with_lifethreatening_illness_f20528_0_0.x",
                               "victim_of_physically_violent_crime_f20529_0_0.x",
                               "witnessed_sudden_violent_death_f20530_0_0.x",
                               "victim_of_sexual_assault_f20531_0_0.x",
                               "age_at_first_episode_of_depression_f20433_0_0.x",
                               "recent_feelings_of_inadequacy_f20507_0_0.x",
                               "recent_trouble_concentrating_on_things_f20508_0_0.x",
                               "recent_feelings_of_depression_f20510_0_0.x",
                               "recent_poor_appetite_or_overeating_f20511_0_0.x",
                               "recent_thoughts_of_suicide_or_selfharm_f20513_0_0.x",
                               "recent_lack_of_interest_or_pleasure_in_doing_things_f20514_0_0.x",
                               "trouble_falling_or_staying_asleep_or_sleeping_too_much_f20517_0_0.x",
                               "recent_changes_in_speedamount_of_moving_or_speaking_f20518_0_0.x",
                               "recent_feelings_of_tiredness_or_low_energy_f20519_0_0.x")))


#----------- Format exposures ---------------

# Since you've been 16 have you..
# Reverse coding so high is worse
XX$rentmortgage_f20525_0_0 <- NA
XX$rentmortgage_f20525_0_0 <- dplyr::recode(XX$able_to_pay_rentmortgage_as_an_adult_f20525_0_0.x,
                                            "Prefer not to answer"="", "Never true"="4", "Rarely true"="3",
                                            "Sometimes true"="2", "Often"="1", "Very often true"="0")
table(XX$able_to_pay_rentmortgage_as_an_adult_f20525_0_0.x,XX$rentmortgage_f20525_0_0)
XX$rentmortgage_f20525_0_0<-as.numeric(XX$rentmortgage_f20525_0_0)


#---- PHQ9------


XX$depressed_mood_f20510_0_0 <- NA
XX$depressed_mood_f20510_0_0 <- dplyr::recode(XX$recent_feelings_of_depression_f20510_0_0.x,
                                              "Prefer not to answer"="", "Not at all"="0",
                                              "Several days"="1", "More than half the days"="2", "Nearly every day" = "3")
table(XX$recent_feelings_of_depression_f20510_0_0.x,XX$depressed_mood_f20510_0_0)
XX$depressed_mood_f20510_0_0<-as.numeric(XX$depressed_mood_f20510_0_0)


XX$anhedonia_f20514_0_0 <- NA
XX$anhedonia_f20514_0_0 <- dplyr::recode(XX$recent_lack_of_interest_or_pleasure_in_doing_things_f20514_0_0.x,
                                         "Prefer not to answer"="", "Not at all"="0",
                                         "Several days"="1", "More than half the days"="2", "Nearly every day" = "3")
table(XX$recent_lack_of_interest_or_pleasure_in_doing_things_f20514_0_0.x,XX$anhedonia_f20514_0_0)
XX$anhedonia_f20514_0_0<-as.numeric(XX$anhedonia_f20514_0_0)


XX$tiredness_f20519_0_0 <- NA
XX$tiredness_f20519_0_0 <- dplyr::recode(XX$recent_feelings_of_tiredness_or_low_energy_f20519_0_0.x,
                                         "Prefer not to answer"="", "Not at all"="0",
                                         "Several days"="1", "More than half the days"="2", "Nearly every day" = "3")
table(XX$recent_feelings_of_tiredness_or_low_energy_f20519_0_0.x,XX$tiredness_f20519_0_0)
XX$tiredness_f20519_0_0<-as.numeric(XX$tiredness_f20519_0_0)


#--------------- General formatting
## 40 PCs are usually considered in all analyses
colnames(XX) <- gsub("genetic_principal_components_f22009_0_", "PC", colnames(XX))

## renaming to use a shorter name in analyses
setnames(XX, "age_at_recruitment_f21022_0_0", "age_f21022")

###### PRSs ######
##### standardize and rename PGS #####


PRS_ACC_micro_GTEx <- fread("/n9/cortexjobs/UKbio_score/u_eprs/jobs/u206_eamon_acc.micro.ukbio.acc_u206/Final_genotype/PRS_result/u206_prs.score.csv")%>%
  mutate(PRS_1.0 = as.numeric(scale(PRS_1.0, center = T, scale = T))) %>%
  dplyr::rename(eid=ID1) %>%
  dplyr::rename(SNP_count_ACC_micro_GTEx =SNP_count_1.0) %>%
  dplyr::rename(PRS_ACC_micro_GTEx =PRS_1.0)

# Merge
PRS <- PRS_ACC_micro_GTEx



## Combine all data sets
XX <- full_join(XX, PRS,     by = "eid")




#######################################
########## exclude drop-outs ##########
#######################################
dropouts <- fread("/share/projects/uk_biobank/w41975_20210809.csv")

XX <- XX %>% filter(!(eid %in% dropouts$V1))

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


XX <- XX[XX$sex == 2,]

XX = as.data.frame(XX)
rownames(XX) <- XX$eid
XX$genotype_array <- as.factor(XX$genotype_array)

#------------------------------ Permutation fun ----------------------------------------------------
run_interaction_permutation_test <- function(data, formula, interaction_term, n_perm) {
  original_fit <- lm(formula, data = data)
  original_stat <- summary(original_fit)$coefficients[interaction_term, "t value"]
  
  perm_stats <- numeric(n_perm)
  
  for(i in 1:n_perm) {
    perm_data <- data
    
    # Everyday I'm shufflin
    perm_data[['PRS_ACC_micro_GTEx']] <- sample(perm_data[['PRS_ACC_micro_GTEx']])
    
    perm_fit <- lm(formula, data = perm_data)
    perm_stats[i] <- summary(perm_fit)$coefficients[interaction_term, "t value"]
  }
  
  p_value <- sum(abs(perm_stats) >= abs(original_stat)) / n_perm
  list(original_stat = original_stat, p_value = p_value, perm_stats = perm_stats)
}


#---------------------------- Rent*anhedonia --------------------------------------------------------
XX.cont <- XX
lm.cont <- lm(anhedonia_f20514_0_0 ~ age_f21022 +PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 +
                PC11 + PC12 + PC13 + PC14 + PC15 + PC16 + PC17 + PC18 + PC19 + PC20 + PC21 + PC22 + PC23 +
                PC24 + PC25 + PC26 + PC27 + PC28 + PC29 + PC30 + PC31 + PC32 + PC33 + PC34 + PC35 + PC36 + 
                PC37 + PC38 + PC39 + PC40 + genotype_array + uk_biobank_assessment_centre_f54_0_0 +  
                PRS_ACC_micro_GTEx * rentmortgage_f20525_0_0, data = XX.cont)

original_coefficient <- summary(lm.cont)$coefficients["PRS_ACC_micro_GTEx:rentmortgage_f20525_0_0", "Estimate"]
original_t_value <- summary(lm.cont)$coefficients["PRS_ACC_micro_GTEx:rentmortgage_f20525_0_0", "t value"]


perm_test_results <- run_interaction_permutation_test(data = XX.cont, formula = anhedonia_f20514_0_0 ~ age_f21022 +PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 +
                                                        PC11 + PC12 + PC13 + PC14 + PC15 + PC16 + PC17 + PC18 + PC19 + PC20 + PC21 + PC22 + PC23 +
                                                        PC24 + PC25 + PC26 + PC27 + PC28 + PC29 + PC30 + PC31 + PC32 + PC33 + PC34 + PC35 + PC36 + 
                                                        PC37 + PC38 + PC39 + PC40 + genotype_array + uk_biobank_assessment_centre_f54_0_0 +  
                                                        PRS_ACC_micro_GTEx * rentmortgage_f20525_0_0, 
                                                      interaction_term = "PRS_ACC_micro_GTEx:rentmortgage_f20525_0_0", n_perm = 1000)

# Prepare data for ggplot
perm_data_1 <- data.frame(T_values = perm_test_results$perm_stats)
original_stat_1 <- perm_test_results$original_stat
Perm_P_1 <- perm_test_results$p_value
Perm_P_1 <- ifelse(Perm_P_1 == 0, ">2.2e-16", as.character(Perm_P_1))
Perm_P_1 <- paste("P =",Perm_P_1)
Perm_P_1 <- "P < 2.2e-16"
# Plot
P1 <- ggplot(perm_data_1, aes(x = T_values)) +
  geom_histogram(aes(y = ..density..), binwidth = .2, fill = "skyblue1", color = "grey69") +
  geom_vline(xintercept = original_stat_1, color = "indianred4", linetype = "dashed", size = 1) +
  geom_density(alpha = .2, fill = "#FF6666") +
  labs(title = "PGS*Ability_to_pay_rent, Anhedonia (outcome)",
       x = "Permuted T-statistic",
       y = "Density") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        text = element_text(size = 16),
        plot.title = element_text(hjust = 0.5, size = 16)) +
  annotate("text", x = -3, y = 0.4, label = Perm_P_1)





#---------- Rent*Depressed-----------
XX.cont <- XX
lm.cont <- lm(depressed_mood_f20510_0_0 ~ age_f21022 +PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 +
                PC11 + PC12 + PC13 + PC14 + PC15 + PC16 + PC17 + PC18 + PC19 + PC20 + PC21 + PC22 + PC23 +
                PC24 + PC25 + PC26 + PC27 + PC28 + PC29 + PC30 + PC31 + PC32 + PC33 + PC34 + PC35 + PC36 + 
                PC37 + PC38 + PC39 + PC40 + genotype_array + uk_biobank_assessment_centre_f54_0_0 +  
                PRS_ACC_micro_GTEx * rentmortgage_f20525_0_0, data = XX.cont)

original_coefficient <- summary(lm.cont)$coefficients["PRS_ACC_micro_GTEx:rentmortgage_f20525_0_0", "Estimate"]
original_t_value <- summary(lm.cont)$coefficients["PRS_ACC_micro_GTEx:rentmortgage_f20525_0_0", "t value"]


perm_test_results <- run_interaction_permutation_test(data = XX.cont, formula = depressed_mood_f20510_0_0 ~ age_f21022 +PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 +
                                                        PC11 + PC12 + PC13 + PC14 + PC15 + PC16 + PC17 + PC18 + PC19 + PC20 + PC21 + PC22 + PC23 +
                                                        PC24 + PC25 + PC26 + PC27 + PC28 + PC29 + PC30 + PC31 + PC32 + PC33 + PC34 + PC35 + PC36 + 
                                                        PC37 + PC38 + PC39 + PC40 + genotype_array + uk_biobank_assessment_centre_f54_0_0 +  
                                                        PRS_ACC_micro_GTEx * rentmortgage_f20525_0_0, 
                                                      interaction_term = "PRS_ACC_micro_GTEx:rentmortgage_f20525_0_0", n_perm = 1000)

# Prepare data for ggplot
perm_data_2 <- data.frame(T_values = perm_test_results$perm_stats)
original_stat_2 <- perm_test_results$original_stat
Perm_P_2 <- perm_test_results$p_value
Perm_P_2 <- ifelse(Perm_P_2 == 0, ">2.2e-16", as.character(Perm_P_2))
Perm_P_2 <- paste("P =",Perm_P_2)

# Plot
P2 <- ggplot(perm_data_2, aes(x = T_values)) +
  geom_histogram(aes(y = ..density..), binwidth = .2, fill = "skyblue1", color = "grey69") +
  geom_vline(xintercept = original_stat_2, color = "indianred4", linetype = "dashed", size = 1) +
  geom_density(alpha = .2, fill = "#FF6666") +
  labs(title = "PGS*Ability_to_pay_rent, Depressed Mood (outcome)",
       x = "Permuted T-statistic",
       y = "Density") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        text = element_text(size = 16),
        plot.title = element_text(hjust = 0.5, size = 16)) +
  annotate("text", x = -3.7, y = 0.45, label = Perm_P_2)


#------------------------ BMI*fatigue------------------------
XX.cont <- XX
lm.cont <- lm(tiredness_f20519_0_0 ~ age_f21022 +PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 +
                PC11 + PC12 + PC13 + PC14 + PC15 + PC16 + PC17 + PC18 + PC19 + PC20 + PC21 + PC22 + PC23 +
                PC24 + PC25 + PC26 + PC27 + PC28 + PC29 + PC30 + PC31 + PC32 + PC33 + PC34 + PC35 + PC36 + 
                PC37 + PC38 + PC39 + PC40 + genotype_array + uk_biobank_assessment_centre_f54_0_0 +  
                PRS_ACC_micro_GTEx * body_mass_index_bmi_f21001_0_0, data = XX.cont)

original_coefficient <- summary(lm.cont)$coefficients["PRS_ACC_micro_GTEx:body_mass_index_bmi_f21001_0_0", "Estimate"]
original_t_value <- summary(lm.cont)$coefficients["PRS_ACC_micro_GTEx:body_mass_index_bmi_f21001_0_0", "t value"]


perm_test_results <- run_interaction_permutation_test(data = XX.cont, formula = tiredness_f20519_0_0 ~ age_f21022 +PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 +
                                                        PC11 + PC12 + PC13 + PC14 + PC15 + PC16 + PC17 + PC18 + PC19 + PC20 + PC21 + PC22 + PC23 +
                                                        PC24 + PC25 + PC26 + PC27 + PC28 + PC29 + PC30 + PC31 + PC32 + PC33 + PC34 + PC35 + PC36 + 
                                                        PC37 + PC38 + PC39 + PC40 + genotype_array + uk_biobank_assessment_centre_f54_0_0 +  
                                                        PRS_ACC_micro_GTEx * body_mass_index_bmi_f21001_0_0, 
                                                      interaction_term = "PRS_ACC_micro_GTEx:body_mass_index_bmi_f21001_0_0", n_perm = 1000)

# Prepare data for ggplot
perm_data_3 <- data.frame(T_values = perm_test_results$perm_stats)
original_stat_3 <- perm_test_results$original_stat
Perm_P_3 <- perm_test_results$p_value
Perm_P_3 <- ifelse(Perm_P_3 == 0, ">2.2e-16", as.character(Perm_P_3))
Perm_P_3 <- paste("P =",Perm_P_3)

# Plot
P3 <- ggplot(perm_data_3, aes(x = T_values)) +
  geom_histogram(aes(y = ..density..), binwidth = .2, fill = "skyblue1", color = "grey69") +
  geom_vline(xintercept = original_stat_3, color = "indianred4", linetype = "dashed", size = 1) +
  geom_density(alpha = .2, fill = "#FF6666") +
  labs(title = "PGS*BMI, Fatigue (Outcome)",
       x = "Permuted T-statistic",
       y = "Density") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        text = element_text(size = 16),
        plot.title = element_text(hjust = 0.5, size = 16)) +
  annotate("text", x = -4, y = 0.45, label = Perm_P_3)

#------------------------ BMI*anhedonia------------------------
XX.cont <- XX
lm.cont <- lm(anhedonia_f20514_0_0 ~ age_f21022 +PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 +
                PC11 + PC12 + PC13 + PC14 + PC15 + PC16 + PC17 + PC18 + PC19 + PC20 + PC21 + PC22 + PC23 +
                PC24 + PC25 + PC26 + PC27 + PC28 + PC29 + PC30 + PC31 + PC32 + PC33 + PC34 + PC35 + PC36 + 
                PC37 + PC38 + PC39 + PC40 + genotype_array + uk_biobank_assessment_centre_f54_0_0 +  
                PRS_ACC_micro_GTEx * body_mass_index_bmi_f21001_0_0, data = XX.cont)

original_coefficient <- summary(lm.cont)$coefficients["PRS_ACC_micro_GTEx:body_mass_index_bmi_f21001_0_0", "Estimate"]
original_t_value <- summary(lm.cont)$coefficients["PRS_ACC_micro_GTEx:body_mass_index_bmi_f21001_0_0", "t value"]


perm_test_results <- run_interaction_permutation_test(data = XX.cont, formula = anhedonia_f20514_0_0 ~ age_f21022 +PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 +
                                                        PC11 + PC12 + PC13 + PC14 + PC15 + PC16 + PC17 + PC18 + PC19 + PC20 + PC21 + PC22 + PC23 +
                                                        PC24 + PC25 + PC26 + PC27 + PC28 + PC29 + PC30 + PC31 + PC32 + PC33 + PC34 + PC35 + PC36 + 
                                                        PC37 + PC38 + PC39 + PC40 + genotype_array + uk_biobank_assessment_centre_f54_0_0 +  
                                                        PRS_ACC_micro_GTEx * body_mass_index_bmi_f21001_0_0, 
                                                      interaction_term = "PRS_ACC_micro_GTEx:body_mass_index_bmi_f21001_0_0", n_perm = 1000)

# Prepare data for ggplot
perm_data_4 <- data.frame(T_values = perm_test_results$perm_stats)
original_stat_4 <- perm_test_results$original_stat
Perm_P_4 <- perm_test_results$p_value
Perm_P_4 <- ifelse(Perm_P_4 == 0, ">2.2e-16", as.character(Perm_P_4))
Perm_P_4 <- paste("P =",Perm_P_4)

# Plot
P4 <- ggplot(perm_data_4, aes(x = T_values)) +
  geom_histogram(aes(y = ..density..), binwidth = .2, fill = "skyblue1", color = "grey69") +
  geom_vline(xintercept = original_stat_4, color = "indianred4", linetype = "dashed", size = 1) +
  geom_density(alpha = .2, fill = "#FF6666") +
  labs(title = "PGS*BMI, Anhedonia (Outcome)",
       x = "Permuted T-statistic",
       y = "Density") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        text = element_text(size = 16),
        plot.title = element_text(hjust = 0.5, size = 16)) +
  annotate("text", x = -3, y = 0.43, label = Perm_P_4)





save(perm_data_1,Perm_P_1, original_stat_1,
     perm_data_2,Perm_P_2, original_stat_2,
     perm_data_3,Perm_P_3, original_stat_3,
     perm_data_4,Perm_P_4, original_stat_4,
     file = "/data1/meaneylab/eamon/Microglia_ePRS/Revisions_BBI/UKBB_perm_res.RDS")

load("/data1/meaneylab/eamon/Microglia_ePRS/Revisions_BBI/UKBB_perm_res.RDS")

plotme <- cowplot::plot_grid(P1,P2,P3,P4, nrow = 2, labels = "AUTO")

ggsave("/data1/meaneylab/eamon/Microglia_ePRS/Revisions_BBI/UKBB_perm.pdf", 
       plot = plotme, width = 30, height = 25, units = "cm", bg="white")
