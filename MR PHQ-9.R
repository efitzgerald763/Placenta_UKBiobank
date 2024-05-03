library(R.utils)
library(mr.raps)
library(data.table)
library(TwoSampleMR)
library(MRPRESSO)
library(devtools)
library(readr)
library(tidyverse)



###########################
## Files were unzipped ####
gunzip("/data1/meaneylab/eamon/Placenta_study/PHQ9/tiredness.20519.gwas.imputed_v3.male.tsv.bgz", "/data1/meaneylab/eamon/Placenta_study/PHQ9/tiredness.20519.gwas.imputed_v3.male.tsv")
###########################

Anhedonia_total <- read_tsv("/data1/meaneylab/eamon/Placenta_study/PHQ9/anhedonia.20514.gwas.imputed_v3.both_sexes.tsv") %>%
  dplyr::select(variant, pval, se, beta, n_complete_samples)
Anhedonia_female <- read_tsv("/data1/meaneylab/eamon/Placenta_study/PHQ9/anhedonia.20514.gwas.imputed_v3.female.tsv")%>%
  dplyr::select(variant, pval, se, beta, n_complete_samples)
Anhedonia_male <- read_tsv("/data1/meaneylab/eamon/Placenta_study/PHQ9/anhedonia.20514.gwas.imputed_v3.male.tsv")%>%
  dplyr::select(variant, pval, se, beta, n_complete_samples)

Depressed_mood_total <- read_tsv("/data1/meaneylab/eamon/Placenta_study/PHQ9/depressedmood.20510.gwas.imputed_v3.both_sexes.tsv")%>%
  dplyr::select(variant, pval, se, beta, n_complete_samples)
Depressed_mood_female <- read_tsv("/data1/meaneylab/eamon/Placenta_study/PHQ9/depressedmood.20510.gwas.imputed_v3.female.tsv")%>%
  dplyr::select(variant, pval, se, beta, n_complete_samples)
Depressed_mood_male <- read_tsv("/data1/meaneylab/eamon/Placenta_study/PHQ9/depressedmood.20510.gwas.imputed_v3.male.tsv")%>%
  dplyr::select(variant, pval, se, beta, n_complete_samples)

Sleep_problems_total <- read_tsv("/data1/meaneylab/eamon/Placenta_study/PHQ9/sleepproblems.20517.gwas.imputed_v3.both_sexes.tsv")%>%
  dplyr::select(variant, pval, se, beta, n_complete_samples)
Sleep_problems_female <- read_tsv("/data1/meaneylab/eamon/Placenta_study/PHQ9/sleepproblems.20517.gwas.imputed_v3.female.tsv")%>%
  dplyr::select(variant, pval, se, beta, n_complete_samples)
Sleep_problems_male <- read_tsv("/data1/meaneylab/eamon/Placenta_study/PHQ9/sleepproblems.20517.gwas.imputed_v3.male.tsv")%>%
  dplyr::select(variant, pval, se, beta, n_complete_samples)

Tiredness_total <- read_tsv("/data1/meaneylab/eamon/Placenta_study/PHQ9/tiredness.20519.gwas.imputed_v3.both_sexes.tsv")%>%
  dplyr::select(variant, pval, se, beta, n_complete_samples)
Tiredness_female <- read_tsv("/data1/meaneylab/eamon/Placenta_study/PHQ9/tiredness.20519.gwas.imputed_v3.female.tsv")%>%
  dplyr::select(variant, pval, se, beta, n_complete_samples)
Tiredness_male <- read_tsv("/data1/meaneylab/eamon/Placenta_study/PHQ9/tiredness.20519.gwas.imputed_v3.male.tsv")%>%
  dplyr::select(variant, pval, se, beta, n_complete_samples)

Changed_appetite_total <- read_tsv("/data1/meaneylab/eamon/Placenta_study/PHQ9/appetite.20511.gwas.imputed_v3.both_sexes.tsv")%>%
  dplyr::select(variant, pval, se, beta, n_complete_samples)
Changed_appetite_female <- read_tsv("/data1/meaneylab/eamon/Placenta_study/PHQ9/appetite.20511.gwas.imputed_v3.female.tsv")%>%
  dplyr::select(variant, pval, se, beta, n_complete_samples)
Changed_appetite_male <- read_tsv("/data1/meaneylab/eamon/Placenta_study/PHQ9/appetite.20511.gwas.imputed_v3.male.tsv")%>%
  dplyr::select(variant, pval, se, beta, n_complete_samples)

Inadequacy_total <- read_tsv("/data1/meaneylab/eamon/Placenta_study/PHQ9/inadequacy.20507.gwas.imputed_v3.both_sexes.tsv")%>%
  dplyr::select(variant, pval, se, beta, n_complete_samples)
Inadequacy_female <- read_tsv("/data1/meaneylab/eamon/Placenta_study/PHQ9/inadequacy.20507.gwas.imputed_v3.female.tsv")%>%
  dplyr::select(variant, pval, se, beta, n_complete_samples)
Inadequacy_male <- read_tsv("/data1/meaneylab/eamon/Placenta_study/PHQ9/inadequacy.20507.gwas.imputed_v3.male.tsv")%>%
  dplyr::select(variant, pval, se, beta, n_complete_samples)

Concentration_problems_total <- read_tsv("/data1/meaneylab/eamon/Placenta_study/PHQ9/concentration.20508.gwas.imputed_v3.both_sexes.tsv")%>%
  dplyr::select(variant, pval, se, beta, n_complete_samples)
Concentration_problems_female <- read_tsv("/data1/meaneylab/eamon/Placenta_study/PHQ9/concentration.20508.gwas.imputed_v3.female.tsv")%>%
  dplyr::select(variant, pval, se, beta, n_complete_samples)
Concentration_problems_male <- read_tsv("/data1/meaneylab/eamon/Placenta_study/PHQ9/concentration.20508.gwas.imputed_v3.male.tsv")%>%
  dplyr::select(variant, pval, se, beta, n_complete_samples)

Psychomotor_changes_total <- read_tsv("/data1/meaneylab/eamon/Placenta_study/PHQ9/psychomotor.20518.gwas.imputed_v3.both_sexes.tsv")%>%
  dplyr::select(variant, pval, se, beta, n_complete_samples)
Psychomotor_changes_female <- read_tsv("/data1/meaneylab/eamon/Placenta_study/PHQ9/psychomotor.20518.gwas.imputed_v3.female.tsv")%>%
  dplyr::select(variant, pval, se, beta, n_complete_samples)
Psychomotor_changes_male <- read_tsv("/data1/meaneylab/eamon/Placenta_study/PHQ9/psychomotor.20518.gwas.imputed_v3.male.tsv")%>%
  dplyr::select(variant, pval, se, beta, n_complete_samples)

Suicidality_total <- read_tsv("/data1/meaneylab/eamon/Placenta_study/PHQ9/suicidality.20513.gwas.imputed_v3.both_sexes.tsv")%>%
  dplyr::select(variant, pval, se, beta, n_complete_samples)
Suicidality_female <- read_tsv("/data1/meaneylab/eamon/Placenta_study/PHQ9/suicidality.20513.gwas.imputed_v3.female.tsv")%>%
  dplyr::select(variant, pval, se, beta, n_complete_samples)
Suicidality_male <- read_tsv("/data1/meaneylab/eamon/Placenta_study/PHQ9/suicidality.20513.gwas.imputed_v3.male.tsv")%>%
  dplyr::select(variant, pval, se, beta, n_complete_samples)


variants <- read_tsv ("/data1/meaneylab/eamon/Placenta_study/PHQ9/variants.tsv") %>%
  dplyr::select(variant, rsid, alt, ref, AF, minor_allele, minor_AF)

#####################################################################
#### Merge GWAS summary stats with variants file and save as .csv ###
#####################################################################

#Anhedonia
Anhedonia_total_var <- merge(variants, Anhedonia_total, by = "variant")
Anhedonia_female_var <- merge(variants, Anhedonia_female, by = "variant")
Anhedonia_male_var <- merge(variants, Anhedonia_male, by = "variant")

write.table(Anhedonia_total_var, file='/data1/meaneylab/eamon/Placenta_study/PHQ9/MR input files/Anhedonia_total_var.csv', sep=',', row.names = FALSE)
write.table(Anhedonia_female_var, file='/data1/meaneylab/eamon/Placenta_study/PHQ9/MR input files/Anhedonia_female_var.csv', sep=',', row.names = FALSE)
write.table(Anhedonia_male_var, file='/data1/meaneylab/eamon/Placenta_study/PHQ9/MR input files/Anhedonia_male_var.csv', sep=',', row.names = FALSE)

#Depressed_mood
Depressed_mood_total_var <- merge(variants, Depressed_mood_total, by = "variant")
Depressed_mood_female_var <- merge(variants, Depressed_mood_female, by = "variant")
Depressed_mood_male_var <- merge(variants, Depressed_mood_male, by = "variant")


write.table(Depressed_mood_total_var, file='/data1/meaneylab/eamon/Placenta_study/PHQ9/MR input files/Depressed_mood_total_var.csv', sep=',', row.names = FALSE)
write.table(Depressed_mood_female_var, file='/data1/meaneylab/eamon/Placenta_study/PHQ9/MR input files/Depressed_mood_female_var.csv', sep=',', row.names = FALSE)
write.table(Depressed_mood_male_var, file='/data1/meaneylab/eamon/Placenta_study/PHQ9/MR input files/Depressed_mood_male_var.csv', sep=',', row.names = FALSE)

# Sleep_problems
Sleep_problems_total_var <- merge(variants, Sleep_problems_total, by = "variant")
Sleep_problems_female_var <- merge(variants, Sleep_problems_female, by = "variant")
Sleep_problems_male_var <- merge(variants, Sleep_problems_male, by = "variant")


write.table(Sleep_problems_total_var, file='/data1/meaneylab/eamon/Placenta_study/PHQ9/MR input files/Sleep_problems_total_var.csv', sep=',', row.names = FALSE)
write.table(Sleep_problems_female_var, file='/data1/meaneylab/eamon/Placenta_study/PHQ9/MR input files/Sleep_problems_female_var.csv', sep=',', row.names = FALSE)
write.table(Sleep_problems_male_var, file='/data1/meaneylab/eamon/Placenta_study/PHQ9/MR input files/Sleep_problems_male_var.csv', sep=',', row.names = FALSE)

# Tiredness
Tiredness_total_var <- merge(variants, Tiredness_total, by = "variant")
Tiredness_female_var <- merge(variants, Tiredness_female, by = "variant")
Tiredness_male_var <- merge(variants, Tiredness_male, by = "variant")


write.table(Tiredness_total_var, file='/data1/meaneylab/eamon/Placenta_study/PHQ9/MR input files/Tiredness_total_var.csv', sep=',', row.names = FALSE)
write.table(Tiredness_female_var, file='/data1/meaneylab/eamon/Placenta_study/PHQ9/MR input files/Tiredness_female_var.csv', sep=',', row.names = FALSE)
write.table(Tiredness_male_var, file='/data1/meaneylab/eamon/Placenta_study/PHQ9/MR input files/Tiredness_male_var.csv', sep=',', row.names = FALSE)

# Changed_appetite
Changed_appetite_total_var <- merge(variants, Changed_appetite_total, by = "variant")
Changed_appetite_female_var <- merge(variants, Changed_appetite_female, by = "variant")
Changed_appetite_male_var <- merge(variants, Changed_appetite_male, by = "variant")


write.table(Changed_appetite_total_var, file='/data1/meaneylab/eamon/Placenta_study/PHQ9/MR input files/Changed_appetite_total_var.csv', sep=',', row.names = FALSE)
write.table(Changed_appetite_female_var, file='/data1/meaneylab/eamon/Placenta_study/PHQ9/MR input files/Changed_appetite_female_var.csv', sep=',', row.names = FALSE)
write.table(Changed_appetite_male_var, file='/data1/meaneylab/eamon/Placenta_study/PHQ9/MR input files/Changed_appetite_male_var.csv', sep=',', row.names = FALSE)

# Inadequacy
Inadequacy_total_var <- merge(variants, Inadequacy_total, by = "variant")
Inadequacy_female_var <- merge(variants, Inadequacy_female, by = "variant")
Inadequacy_male_var <- merge(variants, Inadequacy_male, by = "variant")


write.table(Inadequacy_total_var, file='/data1/meaneylab/eamon/Placenta_study/PHQ9/MR input files/Inadequacy_total_var.csv', sep=',', row.names = FALSE)
write.table(Inadequacy_female_var, file='/data1/meaneylab/eamon/Placenta_study/PHQ9/MR input files/Inadequacy_female_var.csv', sep=',', row.names = FALSE)
write.table(Inadequacy_male_var, file='/data1/meaneylab/eamon/Placenta_study/PHQ9/MR input files/Inadequacy_male_var.csv', sep=',', row.names = FALSE)

# Concentration_problems
Concentration_problems_total_var <- merge(variants, Concentration_problems_total, by = "variant")
Concentration_problems_female_var <- merge(variants, Concentration_problems_female, by = "variant")
Concentration_problems_male_var <- merge(variants, Concentration_problems_male, by = "variant")


write.table(Concentration_problems_total_var, file='/data1/meaneylab/eamon/Placenta_study/PHQ9/MR input files/Concentration_problems_total_var.csv', sep=',', row.names = FALSE)
write.table(Concentration_problems_female_var, file='/data1/meaneylab/eamon/Placenta_study/PHQ9/MR input files/Concentration_problems_female_var.csv', sep=',', row.names = FALSE)
write.table(Concentration_problems_male_var, file='/data1/meaneylab/eamon/Placenta_study/PHQ9/MR input files/Concentration_problems_male_var.csv', sep=',', row.names = FALSE)

# Psychomotor_changes
Psychomotor_changes_total_var <- merge(variants, Psychomotor_changes_total, by = "variant")
Psychomotor_changes_female_var <- merge(variants, Psychomotor_changes_female, by = "variant")
Psychomotor_changes_male_var <- merge(variants, Psychomotor_changes_male, by = "variant")


write.table(Psychomotor_changes_total_var, file='/data1/meaneylab/eamon/Placenta_study/PHQ9/MR input files/Psychomotor_changes_total_var.csv', sep=',', row.names = FALSE)
write.table(Psychomotor_changes_female_var, file='/data1/meaneylab/eamon/Placenta_study/PHQ9/MR input files/Psychomotor_changes_female_var.csv', sep=',', row.names = FALSE)
write.table(Psychomotor_changes_male_var, file='/data1/meaneylab/eamon/Placenta_study/PHQ9/MR input files/Psychomotor_changes_male_var.csv', sep=',', row.names = FALSE)

# Suicidality
Suicidality_total_var <- merge(variants, Suicidality_total, by = "variant")
Suicidality_female_var <- merge(variants, Suicidality_female, by = "variant")
Suicidality_male_var <- merge(variants, Suicidality_male, by = "variant")


write.table(Suicidality_total_var, file='/data1/meaneylab/eamon/Placenta_study/PHQ9/MR input files/Suicidality_total_var.csv', sep=',', row.names = FALSE)
write.table(Suicidality_female_var, file='/data1/meaneylab/eamon/Placenta_study/PHQ9/MR input files/Suicidality_female_var.csv', sep=',', row.names = FALSE)
write.table(Suicidality_male_var, file='/data1/meaneylab/eamon/Placenta_study/PHQ9/MR input files/Suicidality_male_var.csv', sep=',', row.names = FALSE)


###########################################################
################# Extract exposure data ###################
###########################################################

# Get instruments
placenta_eQTLs <- fread("/n9/cortexjobs/common/resources/expr_database/processed/expr.75/expr.75.txt")
placenta_snps_used <- fread("/n9/cortexjobs/UKbio_score/u_eprs/jobs/u163_eamon_placenta.rpkm.cyan22.ukbio.placenta_u163/Final_genotype/PRS_result/u163_prs.snplog.csv")%>%
  dplyr :: select(PRS_1.0)
placenta_genes <- fread("/n9/cortexjobs/UKbio_score/u_eprs/jobs/u163_eamon_placenta.rpkm.cyan22.ukbio.placenta_u163/logTABLE.csv")

# Filter first for genes,then eQTLs used in ePRS
placenta_trimmed <- dplyr :: filter(placenta_eQTLs, trait_id %in% placenta_genes$gtexENS_selectedForGWAS)
placenta_trimmed <- dplyr :: filter(placenta_trimmed, rsid %in% placenta_snps_used$PRS_1.0)

# Rename columns so they're recognized by TwosampleMR
colnames(placenta_trimmed) <- c('Study_ID', 'Feature', 'Variant_ID', 'Chr', 'Position', 'other_allele', 'effect_allele', 'SNP', 'direction', 
                                'beta', 'se', 'pval', 'sample_size', 'eaf', 'r_squared', 'platform_ID', 'Gene_name', 'cis_or_trans', 'q_value')


exposure_dat <- format_data(placenta_trimmed, type="exposure")

###########################################################
################ Extract outcome data #####################
###########################################################

#Anhedonia

Ahendonia_total_outcome_dat<- read_outcome_data(
  snps = exposure_dat$SNP,
  filename = "/data1/meaneylab/eamon/Placenta_study/PHQ9/MR input files/Anhedonia_total_var.csv",
  sep = ",",
  snp_col = "rsid",
  beta_col = "beta",
  se_col = "se",
  effect_allele_col = "alt",
  other_allele_col = "ref",
  eaf_col = "AF",
  pval_col = "pval",
  units_col = FALSE,
  gene_col = FALSE,
  samplesize_col = "n_complete_samples"
)
Ahendonia_female_outcome_dat<- read_outcome_data(
  snps = exposure_dat$SNP,
  filename = "/data1/meaneylab/eamon/Placenta_study/PHQ9/MR input files/Anhedonia_female_var.csv",
  sep = ",",
  snp_col = "rsid",
  beta_col = "beta",
  se_col = "se",
  effect_allele_col = "alt",
  other_allele_col = "ref",
  eaf_col = "AF",
  pval_col = "pval",
  units_col = FALSE,
  gene_col = FALSE,
  samplesize_col = "n_complete_samples"
)
Ahendonia_male_outcome_dat<- read_outcome_data(
  snps = exposure_dat$SNP,
  filename = "/data1/meaneylab/eamon/Placenta_study/PHQ9/MR input files/Anhedonia_male_var.csv",
  sep = ",",
  snp_col = "rsid",
  beta_col = "beta",
  se_col = "se",
  effect_allele_col = "alt",
  other_allele_col = "ref",
  eaf_col = "AF",
  pval_col = "pval",
  units_col = FALSE,
  gene_col = FALSE,
  samplesize_col = "n_complete_samples"
)

#Depressed_mood

Depressed_mood_total_outcome_dat<- read_outcome_data(
  snps = exposure_dat$SNP,
  filename = "/data1/meaneylab/eamon/Placenta_study/PHQ9/MR input files/Depressed_mood_total_var.csv",
  sep = ",",
  snp_col = "rsid",
  beta_col = "beta",
  se_col = "se",
  effect_allele_col = "alt",
  other_allele_col = "ref",
  eaf_col = "AF",
  pval_col = "pval",
  units_col = FALSE,
  gene_col = FALSE,
  samplesize_col = "n_complete_samples"
)
Depressed_mood_female_outcome_dat<- read_outcome_data(
  snps = exposure_dat$SNP,
  filename = "/data1/meaneylab/eamon/Placenta_study/PHQ9/MR input files/Depressed_mood_female_var.csv",
  sep = ",",
  snp_col = "rsid",
  beta_col = "beta",
  se_col = "se",
  effect_allele_col = "alt",
  other_allele_col = "ref",
  eaf_col = "AF",
  pval_col = "pval",
  units_col = FALSE,
  gene_col = FALSE,
  samplesize_col = "n_complete_samples"
)
Depressed_mood_male_outcome_dat<- read_outcome_data(
  snps = exposure_dat$SNP,
  filename = "/data1/meaneylab/eamon/Placenta_study/PHQ9/MR input files/Depressed_mood_male_var.csv",
  sep = ",",
  snp_col = "rsid",
  beta_col = "beta",
  se_col = "se",
  effect_allele_col = "alt",
  other_allele_col = "ref",
  eaf_col = "AF",
  pval_col = "pval",
  units_col = FALSE,
  gene_col = FALSE,
  samplesize_col = "n_complete_samples"
)

# Sleep_problems

Sleep_problems_total_outcome_dat<- read_outcome_data(
  snps = exposure_dat$SNP,
  filename = "/data1/meaneylab/eamon/Placenta_study/PHQ9/MR input files/Sleep_problems_total_var.csv",
  sep = ",",
  snp_col = "rsid",
  beta_col = "beta",
  se_col = "se",
  effect_allele_col = "alt",
  other_allele_col = "ref",
  eaf_col = "AF",
  pval_col = "pval",
  units_col = FALSE,
  gene_col = FALSE,
  samplesize_col = "n_complete_samples"
)
Sleep_problems_female_outcome_dat<- read_outcome_data(
  snps = exposure_dat$SNP,
  filename = "/data1/meaneylab/eamon/Placenta_study/PHQ9/MR input files/Sleep_problems_female_var.csv",
  sep = ",",
  snp_col = "rsid",
  beta_col = "beta",
  se_col = "se",
  effect_allele_col = "alt",
  other_allele_col = "ref",
  eaf_col = "AF",
  pval_col = "pval",
  units_col = FALSE,
  gene_col = FALSE,
  samplesize_col = "n_complete_samples"
)
Sleep_problems_male_outcome_dat<- read_outcome_data(
  snps = exposure_dat$SNP,
  filename = "/data1/meaneylab/eamon/Placenta_study/PHQ9/MR input files/Sleep_problems_male_var.csv",
  sep = ",",
  snp_col = "rsid",
  beta_col = "beta",
  se_col = "se",
  effect_allele_col = "alt",
  other_allele_col = "ref",
  eaf_col = "AF",
  pval_col = "pval",
  units_col = FALSE,
  gene_col = FALSE,
  samplesize_col = "n_complete_samples"
)

# Tiredness
Tiredness_total_outcome_dat<- read_outcome_data(
  snps = exposure_dat$SNP,
  filename = "/data1/meaneylab/eamon/Placenta_study/PHQ9/MR input files/Tiredness_total_var.csv",
  sep = ",",
  snp_col = "rsid",
  beta_col = "beta",
  se_col = "se",
  effect_allele_col = "alt",
  other_allele_col = "ref",
  eaf_col = "AF",
  pval_col = "pval",
  units_col = FALSE,
  gene_col = FALSE,
  samplesize_col = "n_complete_samples"
)
Tiredness_female_outcome_dat<- read_outcome_data(
  snps = exposure_dat$SNP,
  filename = "/data1/meaneylab/eamon/Placenta_study/PHQ9/MR input files/Tiredness_female_var.csv",
  sep = ",",
  snp_col = "rsid",
  beta_col = "beta",
  se_col = "se",
  effect_allele_col = "alt",
  other_allele_col = "ref",
  eaf_col = "AF",
  pval_col = "pval",
  units_col = FALSE,
  gene_col = FALSE,
  samplesize_col = "n_complete_samples"
)
Tiredness_male_outcome_dat<- read_outcome_data(
  snps = exposure_dat$SNP,
  filename = "/data1/meaneylab/eamon/Placenta_study/PHQ9/MR input files/Tiredness_male_var.csv",
  sep = ",",
  snp_col = "rsid",
  beta_col = "beta",
  se_col = "se",
  effect_allele_col = "alt",
  other_allele_col = "ref",
  eaf_col = "AF",
  pval_col = "pval",
  units_col = FALSE,
  gene_col = FALSE,
  samplesize_col = "n_complete_samples"
)

# Changed_appetite
Changed_appetite_total_outcome_dat<- read_outcome_data(
  snps = exposure_dat$SNP,
  filename = "/data1/meaneylab/eamon/Placenta_study/PHQ9/MR input files/Changed_appetite_total_var.csv",
  sep = ",",
  snp_col = "rsid",
  beta_col = "beta",
  se_col = "se",
  effect_allele_col = "alt",
  other_allele_col = "ref",
  eaf_col = "AF",
  pval_col = "pval",
  units_col = FALSE,
  gene_col = FALSE,
  samplesize_col = "n_complete_samples"
)
Changed_appetite_female_outcome_dat<- read_outcome_data(
  snps = exposure_dat$SNP,
  filename = "/data1/meaneylab/eamon/Placenta_study/PHQ9/MR input files/Changed_appetite_female_var.csv",
  sep = ",",
  snp_col = "rsid",
  beta_col = "beta",
  se_col = "se",
  effect_allele_col = "alt",
  other_allele_col = "ref",
  eaf_col = "AF",
  pval_col = "pval",
  units_col = FALSE,
  gene_col = FALSE,
  samplesize_col = "n_complete_samples"
)
Changed_appetite_male_outcome_dat<- read_outcome_data(
  snps = exposure_dat$SNP,
  filename = "/data1/meaneylab/eamon/Placenta_study/PHQ9/MR input files/Changed_appetite_male_var.csv",
  sep = ",",
  snp_col = "rsid",
  beta_col = "beta",
  se_col = "se",
  effect_allele_col = "alt",
  other_allele_col = "ref",
  eaf_col = "AF",
  pval_col = "pval",
  units_col = FALSE,
  gene_col = FALSE,
  samplesize_col = "n_complete_samples"
)

# Inadequacy
Inadequacy_total_outcome_dat<- read_outcome_data(
  snps = exposure_dat$SNP,
  filename = "/data1/meaneylab/eamon/Placenta_study/PHQ9/MR input files/Inadequacy_total_var.csv",
  sep = ",",
  snp_col = "rsid",
  beta_col = "beta",
  se_col = "se",
  effect_allele_col = "alt",
  other_allele_col = "ref",
  eaf_col = "AF",
  pval_col = "pval",
  units_col = FALSE,
  gene_col = FALSE,
  samplesize_col = "n_complete_samples"
)
Inadequacy_female_outcome_dat<- read_outcome_data(
  snps = exposure_dat$SNP,
  filename = "/data1/meaneylab/eamon/Placenta_study/PHQ9/MR input files/Inadequacy_female_var.csv",
  sep = ",",
  snp_col = "rsid",
  beta_col = "beta",
  se_col = "se",
  effect_allele_col = "alt",
  other_allele_col = "ref",
  eaf_col = "AF",
  pval_col = "pval",
  units_col = FALSE,
  gene_col = FALSE,
  samplesize_col = "n_complete_samples"
)
Inadequacy_male_outcome_dat<- read_outcome_data(
  snps = exposure_dat$SNP,
  filename = "/data1/meaneylab/eamon/Placenta_study/PHQ9/MR input files/Inadequacy_male_var.csv",
  sep = ",",
  snp_col = "rsid",
  beta_col = "beta",
  se_col = "se",
  effect_allele_col = "alt",
  other_allele_col = "ref",
  eaf_col = "AF",
  pval_col = "pval",
  units_col = FALSE,
  gene_col = FALSE,
  samplesize_col = "n_complete_samples"
)

# Concentration_problems
Concentration_problems_total_outcome_dat<- read_outcome_data(
  snps = exposure_dat$SNP,
  filename = "/data1/meaneylab/eamon/Placenta_study/PHQ9/MR input files/Concentration_problems_total_var.csv",
  sep = ",",
  snp_col = "rsid",
  beta_col = "beta",
  se_col = "se",
  effect_allele_col = "alt",
  other_allele_col = "ref",
  eaf_col = "AF",
  pval_col = "pval",
  units_col = FALSE,
  gene_col = FALSE,
  samplesize_col = "n_complete_samples"
)
Concentration_problems_female_outcome_dat<- read_outcome_data(
  snps = exposure_dat$SNP,
  filename = "/data1/meaneylab/eamon/Placenta_study/PHQ9/MR input files/Concentration_problems_female_var.csv",
  sep = ",",
  snp_col = "rsid",
  beta_col = "beta",
  se_col = "se",
  effect_allele_col = "alt",
  other_allele_col = "ref",
  eaf_col = "AF",
  pval_col = "pval",
  units_col = FALSE,
  gene_col = FALSE,
  samplesize_col = "n_complete_samples"
)
Concentration_problems_male_outcome_dat<- read_outcome_data(
  snps = exposure_dat$SNP,
  filename = "/data1/meaneylab/eamon/Placenta_study/PHQ9/MR input files/Concentration_problems_male_var.csv",
  sep = ",",
  snp_col = "rsid",
  beta_col = "beta",
  se_col = "se",
  effect_allele_col = "alt",
  other_allele_col = "ref",
  eaf_col = "AF",
  pval_col = "pval",
  units_col = FALSE,
  gene_col = FALSE,
  samplesize_col = "n_complete_samples"
)

# Psychomotor_changes
Psychomotor_changes_total_outcome_dat<- read_outcome_data(
  snps = exposure_dat$SNP,
  filename = "/data1/meaneylab/eamon/Placenta_study/PHQ9/MR input files/Psychomotor_changes_total_var.csv",
  sep = ",",
  snp_col = "rsid",
  beta_col = "beta",
  se_col = "se",
  effect_allele_col = "alt",
  other_allele_col = "ref",
  eaf_col = "AF",
  pval_col = "pval",
  units_col = FALSE,
  gene_col = FALSE,
  samplesize_col = "n_complete_samples"
)
Psychomotor_changes_female_outcome_dat<- read_outcome_data(
  snps = exposure_dat$SNP,
  filename = "/data1/meaneylab/eamon/Placenta_study/PHQ9/MR input files/Psychomotor_changes_female_var.csv",
  sep = ",",
  snp_col = "rsid",
  beta_col = "beta",
  se_col = "se",
  effect_allele_col = "alt",
  other_allele_col = "ref",
  eaf_col = "AF",
  pval_col = "pval",
  units_col = FALSE,
  gene_col = FALSE,
  samplesize_col = "n_complete_samples"
)
Psychomotor_changes_male_outcome_dat<- read_outcome_data(
  snps = exposure_dat$SNP,
  filename = "/data1/meaneylab/eamon/Placenta_study/PHQ9/MR input files/Psychomotor_changes_male_var.csv",
  sep = ",",
  snp_col = "rsid",
  beta_col = "beta",
  se_col = "se",
  effect_allele_col = "alt",
  other_allele_col = "ref",
  eaf_col = "AF",
  pval_col = "pval",
  units_col = FALSE,
  gene_col = FALSE,
  samplesize_col = "n_complete_samples"
)

# Suicidality
Suicidality_total_outcome_dat<- read_outcome_data(
  snps = exposure_dat$SNP,
  filename = "/data1/meaneylab/eamon/Placenta_study/PHQ9/MR input files/Suicidality_total_var.csv",
  sep = ",",
  snp_col = "rsid",
  beta_col = "beta",
  se_col = "se",
  effect_allele_col = "alt",
  other_allele_col = "ref",
  eaf_col = "AF",
  pval_col = "pval",
  units_col = FALSE,
  gene_col = FALSE,
  samplesize_col = "n_complete_samples"
)
Suicidality_female_outcome_dat<- read_outcome_data(
  snps = exposure_dat$SNP,
  filename = "/data1/meaneylab/eamon/Placenta_study/PHQ9/MR input files/Suicidality_female_var.csv",
  sep = ",",
  snp_col = "rsid",
  beta_col = "beta",
  se_col = "se",
  effect_allele_col = "alt",
  other_allele_col = "ref",
  eaf_col = "AF",
  pval_col = "pval",
  units_col = FALSE,
  gene_col = FALSE,
  samplesize_col = "n_complete_samples"
)
Suicidality_male_outcome_dat<- read_outcome_data(
  snps = exposure_dat$SNP,
  filename = "/data1/meaneylab/eamon/Placenta_study/PHQ9/MR input files/Suicidality_male_var.csv",
  sep = ",",
  snp_col = "rsid",
  beta_col = "beta",
  se_col = "se",
  effect_allele_col = "alt",
  other_allele_col = "ref",
  eaf_col = "AF",
  pval_col = "pval",
  units_col = FALSE,
  gene_col = FALSE,
  samplesize_col = "n_complete_samples"
)

#######################################
########### Harmonize #################
#######################################

Anhedonia_total_dat <- harmonise_data(exposure_dat, Ahendonia_total_outcome_dat)
Anhedonia_female_dat <- harmonise_data(exposure_dat, Ahendonia_female_outcome_dat)
Anhedonia_male_dat <- harmonise_data(exposure_dat, Ahendonia_male_outcome_dat)

Depressed_mood_total_dat <- harmonise_data(exposure_dat, Depressed_mood_total_outcome_dat)
Depressed_mood_female_dat <- harmonise_data(exposure_dat, Depressed_mood_female_outcome_dat)
Depressed_mood_male_dat <- harmonise_data(exposure_dat, Depressed_mood_male_outcome_dat)

Sleep_problems_total_dat <- harmonise_data(exposure_dat, Sleep_problems_total_outcome_dat)
Sleep_problems_female_dat <- harmonise_data(exposure_dat, Sleep_problems_female_outcome_dat)
Sleep_problems_male_dat <- harmonise_data(exposure_dat, Sleep_problems_male_outcome_dat)

Tiredness_total_dat <- harmonise_data(exposure_dat, Tiredness_total_outcome_dat)
Tiredness_female_dat <- harmonise_data(exposure_dat, Tiredness_female_outcome_dat)
Tiredness_male_dat <- harmonise_data(exposure_dat, Tiredness_male_outcome_dat)

Changed_appetite_total_dat <- harmonise_data(exposure_dat, Changed_appetite_total_outcome_dat)
Changed_appetite_female_dat <- harmonise_data(exposure_dat, Changed_appetite_female_outcome_dat)
Changed_appetite_male_dat <- harmonise_data(exposure_dat, Changed_appetite_male_outcome_dat)

Inadequacy_total_dat <- harmonise_data(exposure_dat, Inadequacy_total_outcome_dat)
Inadequacy_female_dat <- harmonise_data(exposure_dat, Inadequacy_female_outcome_dat)
Inadequacy_male_dat <- harmonise_data(exposure_dat, Inadequacy_male_outcome_dat)

Concentration_problems_total_dat <- harmonise_data(exposure_dat, Concentration_problems_total_outcome_dat)
Concentration_problems_female_dat <- harmonise_data(exposure_dat, Concentration_problems_female_outcome_dat)
Concentration_problems_male_dat <- harmonise_data(exposure_dat, Concentration_problems_male_outcome_dat)

Psychomotor_total_dat <- harmonise_data(exposure_dat, Psychomotor_changes_total_outcome_dat)
Psychomotor_female_dat <- harmonise_data(exposure_dat, Psychomotor_changes_female_outcome_dat)
Psychomotor_male_dat <- harmonise_data(exposure_dat, Psychomotor_changes_male_outcome_dat)

Suicidality_total_dat <- harmonise_data(exposure_dat, Suicidality_total_outcome_dat)
Suicidality_female_dat <- harmonise_data(exposure_dat, Suicidality_female_outcome_dat)
Suicidality_male_dat <- harmonise_data(exposure_dat, Suicidality_male_outcome_dat)


############################################
########## Run MR ##########################
############################################

Anhedonia_total_dat_res <- mr(Anhedonia_total_dat, method_list=c("mr_weighted_median", "mr_ivw"))
Anhedonia_female_dat_res <- mr(Anhedonia_female_dat, method_list=c("mr_weighted_median", "mr_ivw"))
Anhedonia_male_dat_res <- mr(Anhedonia_male_dat, method_list=c("mr_weighted_median", "mr_ivw"))

Anhedonia_total_dat_res_ORs <- generate_odds_ratios(Anhedonia_total_dat_res)
Anhedonia_female_dat_res_ORs <- generate_odds_ratios(Anhedonia_female_dat_res)
Anhedonia_male_dat_res_ORs <- generate_odds_ratios(Anhedonia_male_dat_res)

Anhedonia_total_dat_res_ORs$Disease <- 'Anhedonia'
Anhedonia_female_dat_res_ORs$Disease <- 'Anhedonia'
Anhedonia_male_dat_res_ORs$Disease <- 'Anhedonia'

Anhedonia_total_dat_res_ORs$Sex <- 'Combined'
Anhedonia_female_dat_res_ORs$Sex <- 'Female'
Anhedonia_male_dat_res_ORs$Sex <- 'Male'



Depressed_mood_total_dat_res <- mr(Depressed_mood_total_dat, method_list=c("mr_weighted_median", "mr_ivw"))
Depressed_mood_female_dat_res <- mr(Depressed_mood_female_dat, method_list=c("mr_weighted_median", "mr_ivw"))
Depressed_mood_male_dat_res <- mr(Depressed_mood_male_dat, method_list=c("mr_weighted_median", "mr_ivw"))

Depressed_mood_total_dat_res_ORs <- generate_odds_ratios(Depressed_mood_total_dat_res)
Depressed_mood_female_dat_res_ORs <- generate_odds_ratios(Depressed_mood_female_dat_res)
Depressed_mood_male_dat_res_ORs <- generate_odds_ratios(Depressed_mood_male_dat_res)

Depressed_mood_total_dat_res_ORs$Disease <- 'Depressed mood'
Depressed_mood_female_dat_res_ORs$Disease <- 'Depressed mood'
Depressed_mood_male_dat_res_ORs$Disease <- 'Depressed mood'

Depressed_mood_total_dat_res_ORs$Sex <- 'Combined'
Depressed_mood_female_dat_res_ORs$Sex <- 'Female'
Depressed_mood_male_dat_res_ORs$Sex <- 'Male'



Sleep_problems_total_dat_res <- mr(Sleep_problems_total_dat, method_list=c("mr_weighted_median", "mr_ivw"))
Sleep_problems_female_dat_res <- mr(Sleep_problems_female_dat, method_list=c("mr_weighted_median", "mr_ivw"))
Sleep_problems_male_dat_res <- mr(Sleep_problems_male_dat, method_list=c("mr_weighted_median", "mr_ivw"))

Sleep_problems_total_dat_res_ORs <- generate_odds_ratios(Sleep_problems_total_dat_res)
Sleep_problems_female_dat_res_ORs <- generate_odds_ratios(Sleep_problems_female_dat_res)
Sleep_problems_male_dat_res_ORs <- generate_odds_ratios(Sleep_problems_male_dat_res)

Sleep_problems_total_dat_res_ORs$Disease <- 'Sleep problems'
Sleep_problems_female_dat_res_ORs$Disease <- 'Sleep problems'
Sleep_problems_male_dat_res_ORs$Disease <- 'Sleep problems'

Sleep_problems_total_dat_res_ORs$Sex <- 'Combined'
Sleep_problems_female_dat_res_ORs$Sex <- 'Female'
Sleep_problems_male_dat_res_ORs$Sex <- 'Male'


Tiredness_total_dat_res <- mr(Tiredness_total_dat, method_list=c("mr_weighted_median", "mr_ivw"))
Tiredness_female_dat_res <- mr(Tiredness_female_dat, method_list=c("mr_weighted_median", "mr_ivw"))
Tiredness_male_dat_res <- mr(Tiredness_male_dat, method_list=c("mr_weighted_median", "mr_ivw"))

Tiredness_total_dat_res_ORs <- generate_odds_ratios(Tiredness_total_dat_res)
Tiredness_female_dat_res_ORs <- generate_odds_ratios(Tiredness_female_dat_res)
Tiredness_male_dat_res_ORs <- generate_odds_ratios(Tiredness_male_dat_res)

Tiredness_total_dat_res_ORs$Disease <- 'Tiredness'
Tiredness_female_dat_res_ORs$Disease <- 'Tiredness'
Tiredness_male_dat_res_ORs$Disease <- 'Tiredness'

Tiredness_total_dat_res_ORs$Sex <- 'Combined'
Tiredness_female_dat_res_ORs$Sex <- 'Female'
Tiredness_male_dat_res_ORs$Sex <- 'Male'


Changed_appetite_total_dat_res <- mr(Changed_appetite_total_dat, method_list=c("mr_weighted_median", "mr_ivw"))
Changed_appetite_female_dat_res <- mr(Changed_appetite_female_dat, method_list=c("mr_weighted_median", "mr_ivw"))
Changed_appetite_male_dat_res <- mr(Changed_appetite_male_dat, method_list=c("mr_weighted_median", "mr_ivw"))

Changed_appetite_total_dat_res_ORs <- generate_odds_ratios(Changed_appetite_total_dat_res)
Changed_appetite_female_dat_res_ORs <- generate_odds_ratios(Changed_appetite_female_dat_res)
Changed_appetite_male_dat_res_ORs <- generate_odds_ratios(Changed_appetite_male_dat_res)

Changed_appetite_total_dat_res_ORs$Disease <- 'Change in appetite'
Changed_appetite_female_dat_res_ORs$Disease <- 'Change in appetite'
Changed_appetite_male_dat_res_ORs$Disease <- 'Change in appetite'

Changed_appetite_total_dat_res_ORs$Sex <- 'Combined'
Changed_appetite_female_dat_res_ORs$Sex <- 'Female'
Changed_appetite_male_dat_res_ORs$Sex <- 'Male'


Inadequacy_total_dat_res <- mr(Inadequacy_total_dat, method_list=c("mr_weighted_median", "mr_ivw"))
Inadequacy_female_dat_res <- mr(Inadequacy_female_dat, method_list=c("mr_weighted_median", "mr_ivw"))
Inadequacy_male_dat_res <- mr(Inadequacy_male_dat, method_list=c("mr_weighted_median", "mr_ivw"))

Inadequacy_total_dat_res_ORs <- generate_odds_ratios(Inadequacy_total_dat_res)
Inadequacy_female_dat_res_ORs <- generate_odds_ratios(Inadequacy_female_dat_res)
Inadequacy_male_dat_res_ORs <- generate_odds_ratios(Inadequacy_male_dat_res)

Inadequacy_total_dat_res_ORs$Disease <- 'Feelings of inadequacy'
Inadequacy_female_dat_res_ORs$Disease <- 'Feelings of inadequacy'
Inadequacy_male_dat_res_ORs$Disease <- 'Feelings of inadequacy'

Inadequacy_total_dat_res_ORs$Sex <- 'Combined'
Inadequacy_female_dat_res_ORs$Sex <- 'Female'
Inadequacy_male_dat_res_ORs$Sex <- 'Male'


Concentration_problems_total_dat_res <- mr(Concentration_problems_total_dat, method_list=c("mr_weighted_median", "mr_ivw"))
Concentration_problems_female_dat_res <- mr(Concentration_problems_female_dat, method_list=c("mr_weighted_median", "mr_ivw"))
Concentration_problems_male_dat_res <- mr(Concentration_problems_male_dat, method_list=c("mr_weighted_median", "mr_ivw"))

Concentration_problems_total_dat_res_ORs <- generate_odds_ratios(Concentration_problems_total_dat_res)
Concentration_problems_female_dat_res_ORs <- generate_odds_ratios(Concentration_problems_female_dat_res)
Concentration_problems_male_dat_res_ORs <- generate_odds_ratios(Concentration_problems_male_dat_res)

Concentration_problems_total_dat_res_ORs$Disease <- 'Concentration problems'
Concentration_problems_female_dat_res_ORs$Disease <- 'Concentration problems'
Concentration_problems_male_dat_res_ORs$Disease <- 'Concentration problems'

Concentration_problems_total_dat_res_ORs$Sex <- 'Combined'
Concentration_problems_female_dat_res_ORs$Sex <- 'Female'
Concentration_problems_male_dat_res_ORs$Sex <- 'Male'


Psychomotor_total_dat_res <- mr(Psychomotor_total_dat, method_list=c("mr_weighted_median", "mr_ivw"))
Psychomotor_female_dat_res <- mr(Psychomotor_female_dat, method_list=c("mr_weighted_median", "mr_ivw"))
Psychomotor_male_dat_res <- mr(Psychomotor_male_dat, method_list=c("mr_weighted_median", "mr_ivw"))

Psychomotor_total_dat_res_ORs <- generate_odds_ratios(Psychomotor_total_dat_res)
Psychomotor_female_dat_res_ORs <- generate_odds_ratios(Psychomotor_female_dat_res)
Psychomotor_male_dat_res_ORs <- generate_odds_ratios(Psychomotor_male_dat_res)

Psychomotor_total_dat_res_ORs$Disease <- 'Psychomotor problems'
Psychomotor_female_dat_res_ORs$Disease <- 'Psychomotor problems'
Psychomotor_male_dat_res_ORs$Disease <- 'Psychomotor problems'

Psychomotor_total_dat_res_ORs$Sex <- 'Combined'
Psychomotor_female_dat_res_ORs$Sex <- 'Female'
Psychomotor_male_dat_res_ORs$Sex <- 'Male'


Suicidality_total_dat_res <- mr(Suicidality_total_dat, method_list=c("mr_weighted_median", "mr_ivw"))
Suicidality_female_dat_res <- mr(Suicidality_female_dat, method_list=c("mr_weighted_median", "mr_ivw"))
Suicidality_male_dat_res <- mr(Suicidality_male_dat, method_list=c("mr_weighted_median", "mr_ivw"))

Suicidality_total_dat_res_ORs <- generate_odds_ratios(Suicidality_total_dat_res)
Suicidality_female_dat_res_ORs <- generate_odds_ratios(Suicidality_female_dat_res)
Suicidality_male_dat_res_ORs <- generate_odds_ratios(Suicidality_male_dat_res)

Suicidality_total_dat_res_ORs$Disease <- 'Suicidality'
Suicidality_female_dat_res_ORs$Disease <- 'Suicidality'
Suicidality_male_dat_res_ORs$Disease <- 'Suicidality'

Suicidality_total_dat_res_ORs$Sex <- 'Combined'
Suicidality_female_dat_res_ORs$Sex <- 'Female'
Suicidality_male_dat_res_ORs$Sex <- 'Male'

mr_report(
  Tiredness_total_dat,
    output_path = "/data1/meaneylab/eamon/Placenta_study/PHQ9/MR PHQ9 reports/Tiredness_total",
    output_type = "pdf",
    author = "Analyst",
    study = "Suicidality_female",
    path = system.file("reports", package = "TwoSampleMR"))

Total <- rbind(Anhedonia_total_dat_res_ORs, Anhedonia_female_dat_res_ORs, Anhedonia_male_dat_res_ORs, 
               Depressed_mood_total_dat_res_ORs,Depressed_mood_female_dat_res_ORs, Depressed_mood_male_dat_res_ORs, 
               Sleep_problems_total_dat_res_ORs, Sleep_problems_female_dat_res_ORs, Sleep_problems_male_dat_res_ORs, 
               Tiredness_total_dat_res_ORs, Tiredness_female_dat_res_ORs, Tiredness_male_dat_res_ORs,
               Changed_appetite_total_dat_res_ORs,Changed_appetite_female_dat_res_ORs,Changed_appetite_male_dat_res_ORs,
               Inadequacy_total_dat_res_ORs, Inadequacy_female_dat_res_ORs, Inadequacy_male_dat_res_ORs,
               Concentration_problems_total_dat_res_ORs, Concentration_problems_female_dat_res_ORs, Concentration_problems_male_dat_res_ORs,
               Psychomotor_total_dat_res_ORs, Psychomotor_female_dat_res_ORs, Psychomotor_male_dat_res_ORs,
               Suicidality_total_dat_res_ORs, Suicidality_female_dat_res_ORs, Suicidality_male_dat_res_ORs)
IVW_res <- Total[Total$method == 'Inverse variance weighted',]

write_csv(Total, "/data1/meaneylab/eamon/Placenta_study/PHQ9/Total_results.csv")
write_csv(IVW_res, "/data1/meaneylab/eamon/Placenta_study/PHQ9/IVW_results.csv")


IVW_res <- read_csv("/data1/meaneylab/eamon/Placenta_study/PHQ9/IVW_results.csv")
IVW_resM <- IVW_res[IVW_res$Sex == 'Male',]
IVW_resF <- IVW_res[IVW_res$Sex == 'Female',]
IVW_comb <- IVW_res[IVW_res$Sex == 'Combined',]

#Adjust p-values for multiple comparisons
IVW_resM = IVW_resM[order(IVW_resM$pval),]
IVW_resF = IVW_resF[order(IVW_resF$pval),]
IVW_comb = IVW_comb[order(IVW_comb$pval),]

IVW_resM$Bonferroni =
  p.adjust(IVW_resM$pval,
           method = "bonferroni")
IVW_resF$Bonferroni =
  p.adjust(IVW_resF$pval,
           method = "bonferroni")
IVW_comb$Bonferroni =
  p.adjust(IVW_comb$pval,
           method = "bonferroni")

IVW_resM$BH =
  p.adjust(IVW_resM$pval,
           method = "BH")
IVW_resF$BH =
  p.adjust(IVW_resF$pval,
           method = "BH")
IVW_comb$BH =
  p.adjust(IVW_comb$pval,
           method = "BH")

write_csv(IVW_resM, "/data1/meaneylab/eamon/Placenta_study/PHQ9/IVW_resultsMF_Male.csv")
write_csv(IVW_resF, "/data1/meaneylab/eamon/Placenta_study/PHQ9/IVW_results_Female.csv")
write_csv(IVW_comb, "/data1/meaneylab/eamon/Placenta_study/PHQ9/IVW_results_Combined.csv")

Graph <- read_csv("/data1/meaneylab/eamon/Placenta_study/PHQ9/IVW_resultsMF.csv")

#######################
###### MR PRESSO ######
#######################
data(SummaryStats)

# Run MR-PRESSO global method
trial <- mr_presso(BetaOutcome = "Y_effect", BetaExposure = "E1_effect", SdOutcome = "Y_se", SdExposure = "E1_se", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = SummaryStats, NbDistribution = 1000,  SignifThreshold = 0.05)

# Run MR-PRESSO on a multi-variable MR (MMR) model specifying several exposures
trial <- mr_presso(BetaOutcome = "Y_effect", BetaExposure = c("E1_effect", "E2_effect"), SdOutcome = "Y_se", SdExposure = c("E1_se", "E2_se"), OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = SummaryStats, NbDistribution = 1000,  SignifThreshold = 0.05)

mix <- run_mrmix(Suicidality_female_dat)

MR_PRESSO <- run_mr_presso(Suicidality_female_dat, NbDistribution = 1000, SignifThreshold = 0.05)

# Work in progress
MR_RAPS <- mr(Suicidality_female_dat, method_list = c("mr_raps"))
MR_RAPS <- mr(Suicidality_female_dat, method_list = c("mr_raps"), parameters = list(over.dispersion = FALSE))


## Forest Plot ##
IVW_res <- read_csv("/data1/meaneylab/eamon/Placenta_study/PHQ9/IVW_results.csv")

# To order as in file
IVW_res$Disease <- as.character(IVW_res$Disease)
IVW_res$Disease  <- factor(IVW_res$Disease, levels=unique(IVW_res$Disease))

ggplot(IVW_res, aes(x=Disease, y=b, ymin=lo_ci, ymax=up_ci, col = Sex)) +
  geom_pointrange(position=position_dodge(width = 0.5)) +
  geom_hline(yintercept=0, lty=2) +  # add a dotted line at x=0 after flip
  coord_flip() +  # flip coordinates (puts labels on y axis)
  ylab("Estimate (95% CI)")  +
  theme_bw() + # use a white background
  theme(axis.text.y =element_text(hjust=0.5,vjust=0.5)) +
  scale_color_paletteer_d(`"awtools::mpalette"`) +
  theme(legend.title=element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y=element_blank(),
        axis.ticks.y=element_blank(),
        panel.spacing = unit(2, "lines")) +
  theme(text = element_text(size=15)) 

### Facet plot ###
ggplot(IVW_res, aes(x=Disease, y=b, ymin=lo_ci, ymax=up_ci, col = Sex)) +
  geom_pointrange(position=position_dodge(width = 2)) +
  facet_wrap(~Disease, ncol = 3, scales = 'free_y') +
  geom_hline(yintercept=0, lty=2) +  # add a dotted line at x=0 after flip
  coord_flip() +  # flip coordinates (puts labels on y axis)
  ylab("Estimate (95% CI)")  +
  theme_bw() + # use a white background
  theme(axis.text.y =element_text(hjust=0.5,vjust=0.5)) +
  scale_color_paletteer_d(`"awtools::mpalette"`) +
  theme(legend.title=element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        panel.spacing = unit(2, "lines")) +
  theme(text = element_text(size=15)) 

#####################################################
res_single <- mr_singlesnp(Suicidality_female_dat)
res_single <- res_single[-c(30), ]

res_loo <- mr_leaveoneout(Suicidality_female_dat)


p1 <- mr_scatter_plot(Suicidality_female_dat_res, Suicidality_female_dat)
p1[[1]] + 
  theme(legend.title=element_blank()) +
  theme_bw() +
  theme(legend.position = "top") +
  theme(legend.title=element_blank()) +
  geom_hline(yintercept=0, color = "grey44", lty=2) +
  labs(y = "SNP effect on suicidality (females only)") +
  labs(x = "SNP effect on placenta cyan module eQTLs") 


p2 <- mr_forest_plot(res_single)
p2[[1]] + theme_bw() + theme(legend.position="none") +
  labs(x = "MR single SNP effect size for Placenta eQTLs on suicidality (females only)") 


p3 <- mr_leaveoneout_plot(res_loo)
p3[[1]] + theme_bw() + theme(legend.position="none") +
  labs(x = "MR leave one out; Outcome- suicidality (females only)") 
