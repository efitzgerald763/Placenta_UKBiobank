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

# Assuming XX is your dataframe and BMI is the column of interest
# First, create intervals of 10 for BMI values
# Ensure that BMI contains only finite values
BMI_finite <- XX$body_mass_index_bmi_f21001_0_0[is.finite(XX$body_mass_index_bmi_f21001_0_0)]

# Check if there are any finite values to work with
if(length(BMI_finite) > 0) {
  # Create intervals of 10 for BMI values, considering only finite values
  BMI_intervals <- cut(BMI_finite,
                       breaks = seq(from = floor(min(BMI_finite)/10)*10,
                                    to = ceiling(max(BMI_finite)/10)*10,
                                    by = 10),
                       include.lowest = TRUE)
  
  # Create a table to count the number of values in each interval
  BMI_table <- table(BMI_intervals)
  
  # Print the table
  print(BMI_table)
} else {
  cat("No finite BMI values to process.")
}

write.csv(BMI_table_df,"/data1/meaneylab/eamon/Microglia_ePRS/Revisions_BBI/BMI_table.csv")

# Now, create a table to count the number of values in each interval
BMI_table_df <- as.data.frame(BMI_table)
max(BMI_finite)
min(BMI_finite)
mean(BMI_finite)
median(BMI_finite)

BMI_finite_df <- as.data.frame(BMI_finite)
plotme <- ggplot(BMI_finite_df, aes(x = BMI_finite)) +
  geom_histogram(binwidth = 1, fill = "lavender", color = "grey69") +
  labs(title = "BMI (females)",
       x = "BMI",
       y = "Count") +
  scale_y_log10()+  
  scale_x_continuous(breaks = seq(from = 10, to = 80, by = 10)) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        text = element_text(size = 16),
        plot.title = element_text(hjust = 0.5, size = 16)) 
ggsave("/data1/meaneylab/eamon/Microglia_ePRS/Revisions_BBI/BMI_histogram.png", 
       plot = plotme, width = 15, height = 12.5, units = "cm", bg="white")
