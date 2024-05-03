
Anhedonia_total <- read_tsv("/data1/meaneylab/eamon/Placenta_study/PHQ9/anhedonia.20514.gwas.imputed_v3.both_sexes.tsv") %>%
  dplyr::select(n_complete_samples) 
Anhedonia_total <- mean(Anhedonia_total$n_complete_samples)

Anhedonia_female <- read_tsv("/data1/meaneylab/eamon/Placenta_study/PHQ9/anhedonia.20514.gwas.imputed_v3.female.tsv")%>%
  dplyr::select(n_complete_samples) 
Anhedonia_female <- mean(Anhedonia_female$n_complete_samples)

Anhedonia_male <- read_tsv("/data1/meaneylab/eamon/Placenta_study/PHQ9/anhedonia.20514.gwas.imputed_v3.male.tsv")%>%
  dplyr::select(n_complete_samples) 
Anhedonia_male <- mean(Anhedonia_male$n_complete_samples)

Depressed_mood_total <- read_tsv("/data1/meaneylab/eamon/Placenta_study/PHQ9/depressedmood.20510.gwas.imputed_v3.both_sexes.tsv")%>%
  dplyr::select(n_complete_samples) 
Depressed_mood_total <- mean(Depressed_mood_total$n_complete_samples)

Depressed_mood_female <- read_tsv("/data1/meaneylab/eamon/Placenta_study/PHQ9/depressedmood.20510.gwas.imputed_v3.female.tsv")%>%
  dplyr::select(n_complete_samples) 
Depressed_mood_female <- mean(Depressed_mood_female$n_complete_samples)

Depressed_mood_male <- read_tsv("/data1/meaneylab/eamon/Placenta_study/PHQ9/depressedmood.20510.gwas.imputed_v3.male.tsv")%>%
  dplyr::select(n_complete_samples) 
Depressed_mood_male <- mean(Depressed_mood_male$n_complete_samples)

Sleep_problems_total <- read_tsv("/data1/meaneylab/eamon/Placenta_study/PHQ9/sleepproblems.20517.gwas.imputed_v3.both_sexes.tsv")%>%
  dplyr::select(n_complete_samples) 
Sleep_problems_total <- mean(Sleep_problems_total$n_complete_samples)
Sleep_problems_female <- read_tsv("/data1/meaneylab/eamon/Placenta_study/PHQ9/sleepproblems.20517.gwas.imputed_v3.female.tsv")%>%
  dplyr::select(n_complete_samples) 
Sleep_problems_female <- mean(Sleep_problems_female$n_complete_samples)
Sleep_problems_male <- read_tsv("/data1/meaneylab/eamon/Placenta_study/PHQ9/sleepproblems.20517.gwas.imputed_v3.male.tsv")%>%
  dplyr::select(n_complete_samples) 
Sleep_problems_male <- mean(Sleep_problems_male$n_complete_samples)

Tiredness_total <- read_tsv("/data1/meaneylab/eamon/Placenta_study/PHQ9/tiredness.20519.gwas.imputed_v3.both_sexes.tsv")%>%
  dplyr::select(n_complete_samples) 
Tiredness_total <- mean(Tiredness_total$n_complete_samples)
Tiredness_female <- read_tsv("/data1/meaneylab/eamon/Placenta_study/PHQ9/tiredness.20519.gwas.imputed_v3.female.tsv")%>%
  dplyr::select(n_complete_samples) 
Tiredness_female <- mean(Tiredness_female$n_complete_samples)
Tiredness_male <- read_tsv("/data1/meaneylab/eamon/Placenta_study/PHQ9/tiredness.20519.gwas.imputed_v3.male.tsv")%>%
  dplyr::select(n_complete_samples) 
Tiredness_male <- mean(Tiredness_male$n_complete_samples)

Changed_appetite_total <- read_tsv("/data1/meaneylab/eamon/Placenta_study/PHQ9/appetite.20511.gwas.imputed_v3.both_sexes.tsv")%>%
  dplyr::select(n_complete_samples) 
Changed_appetite_total <- mean(Changed_appetite_total$n_complete_samples)
Changed_appetite_female <- read_tsv("/data1/meaneylab/eamon/Placenta_study/PHQ9/appetite.20511.gwas.imputed_v3.female.tsv")%>%
  dplyr::select(n_complete_samples) 
Changed_appetite_female <- mean(Changed_appetite_female$n_complete_samples)
Changed_appetite_male <- read_tsv("/data1/meaneylab/eamon/Placenta_study/PHQ9/appetite.20511.gwas.imputed_v3.male.tsv")%>%
  dplyr::select(n_complete_samples) 
Changed_appetite_male <- mean(Changed_appetite_male$n_complete_samples)

Inadequacy_total <- read_tsv("/data1/meaneylab/eamon/Placenta_study/PHQ9/inadequacy.20507.gwas.imputed_v3.both_sexes.tsv")%>%
  dplyr::select(n_complete_samples) 
Inadequacy_total <- mean(Inadequacy_total$n_complete_samples)
Inadequacy_female <- read_tsv("/data1/meaneylab/eamon/Placenta_study/PHQ9/inadequacy.20507.gwas.imputed_v3.female.tsv")%>%
  dplyr::select(n_complete_samples) 
Inadequacy_female <- mean(Inadequacy_female$n_complete_samples)
Inadequacy_male <- read_tsv("/data1/meaneylab/eamon/Placenta_study/PHQ9/inadequacy.20507.gwas.imputed_v3.male.tsv")%>%
  dplyr::select(n_complete_samples) 
Inadequacy_male <- mean(Inadequacy_male$n_complete_samples)

Concentration_problems_total <- read_tsv("/data1/meaneylab/eamon/Placenta_study/PHQ9/concentration.20508.gwas.imputed_v3.both_sexes.tsv")%>%
  dplyr::select(n_complete_samples) 
Concentration_problems_total <- mean(Concentration_problems_total$n_complete_samples)
Concentration_problems_female <- read_tsv("/data1/meaneylab/eamon/Placenta_study/PHQ9/concentration.20508.gwas.imputed_v3.female.tsv")%>%
  dplyr::select(n_complete_samples) 
Concentration_problems_female <- mean(Concentration_problems_female$n_complete_samples)
Concentration_problems_male <- read_tsv("/data1/meaneylab/eamon/Placenta_study/PHQ9/concentration.20508.gwas.imputed_v3.male.tsv")%>%
  dplyr::select(n_complete_samples) 
Concentration_problems_male <- mean(Concentration_problems_male$n_complete_samples)

Psychomotor_changes_total <- read_tsv("/data1/meaneylab/eamon/Placenta_study/PHQ9/psychomotor.20518.gwas.imputed_v3.both_sexes.tsv")%>%
  dplyr::select(n_complete_samples) 
Psychomotor_changes_total <- mean(Psychomotor_changes_total$n_complete_samples)
Psychomotor_changes_female <- read_tsv("/data1/meaneylab/eamon/Placenta_study/PHQ9/psychomotor.20518.gwas.imputed_v3.female.tsv")%>%
  dplyr::select(n_complete_samples) 
Psychomotor_changes_female <- mean(Psychomotor_changes_female$n_complete_samples)
Psychomotor_changes_male <- read_tsv("/data1/meaneylab/eamon/Placenta_study/PHQ9/psychomotor.20518.gwas.imputed_v3.male.tsv")%>%
  dplyr::select(n_complete_samples) 
Psychomotor_changes_male <- mean(Psychomotor_changes_male$n_complete_samples)

Suicidality_total <- read_tsv("/data1/meaneylab/eamon/Placenta_study/PHQ9/suicidality.20513.gwas.imputed_v3.both_sexes.tsv")%>%
  dplyr::select(n_complete_samples) 
Suicidality_total <- mean(Suicidality_total$n_complete_samples)
Suicidality_female <- read_tsv("/data1/meaneylab/eamon/Placenta_study/PHQ9/suicidality.20513.gwas.imputed_v3.female.tsv")%>%
  dplyr::select(n_complete_samples) 
Suicidality_female <- mean(Suicidality_female$n_complete_samples)
Suicidality_male <- read_tsv("/data1/meaneylab/eamon/Placenta_study/PHQ9/suicidality.20513.gwas.imputed_v3.male.tsv")%>%
  dplyr::select(n_complete_samples) 
Suicidality_male <- mean(Suicidality_male$n_complete_samples)

All2gether <- rbind(Anhedonia_total, Depressed_mood_total, Sleep_problems_total, Tiredness_total, Changed_appetite_total,
                    Inadequacy_total, Concentration_problems_total, Psychomotor_changes_total, Suicidality_total,
                    
                    Anhedonia_female, Depressed_mood_female, Sleep_problems_female, Tiredness_female, Changed_appetite_female,
                    Inadequacy_female, Concentration_problems_female, Psychomotor_changes_female, Suicidality_female,
                    
                    Anhedonia_male, Depressed_mood_male, Sleep_problems_male, Tiredness_male, Changed_appetite_male,
                    Inadequacy_male, Concentration_problems_male, Psychomotor_changes_male, Suicidality_male)
write.csv(All2gether, "/data1/meaneylab/eamon/Placenta_study/Revisons_2.0/PHQ9_N.csv")
