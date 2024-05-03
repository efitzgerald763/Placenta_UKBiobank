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

XX_female <- subset(XX, Gender =="Female")

rownames(XX_female) <- XX_female$UniqueID
XX.cont <- XX_female


#----- Mild Dep 24 years --------
glm.cont.old <- glm(FKDQ1000 ~ b023 + bestgest + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + 
                      PC7 + PC8 + PC9 + PC10 + ePRS_microglia * EPDS_mean_prenatal, data = XX.cont, family = binomial)
nrow_coef <- nrow(summary(glm.cont.old)$coef)

im.cont <- data.frame(influence.measures(glm.cont.old)$is.inf)
case.exclude.cont <- rownames(im.cont[which(im.cont[,nrow_coef] == "TRUE"),])
if (length(case.exclude.cont)>0){XX.cont <- XX.cont[-which(XX.cont$UniqueID %in% case.exclude.cont),]}

# Original model
glm.cont <- glm(FKDQ1000 ~ b023 + bestgest + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + 
                  PC7 + PC8 + PC9 + PC10 + ePRS_microglia * EPDS_mean_prenatal, 
                data = XX.cont, family = binomial)

# Extract the original coefficient of interest
original_coef <- coef(summary(glm.cont))["ePRS_microglia:EPDS_mean_prenatal", "Estimate"]

# Set up number of permutations
n_perm <- 1000 # You can increase this for more robust results

# Initialize a vector to store the permuted coefficients
perm_coefs <- numeric(n_perm)

# Permutation test
set.seed(123) # For reproducibility
for(i in 1:n_perm) {
  # Shuffle the outcome
  shuffled_data <- XX.cont
  shuffled_data$FKDQ1000 <- sample(XX.cont$FKDQ1000)
  
  # Fit the model with the shuffled data
  perm_model <- glm(FKDQ1000 ~ b023 + bestgest + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + 
                      PC7 + PC8 + PC9 + PC10 + ePRS_microglia * EPDS_mean_prenatal, 
                    data = shuffled_data, family = binomial)
  
  # Store the coefficient of interest
  perm_coefs[i] <- coef(summary(perm_model))["ePRS_microglia:EPDS_mean_prenatal", "Estimate"]
}

# Calculate p-value
p_value <- mean(abs(perm_coefs) >= abs(original_coef))

# Create a data frame for plotting
perm_data <- data.frame(T_values = perm_coefs)

# Original coefficient as a data frame for annotation
original_coef_df <- data.frame(original_coef = original_coef)

# Calculate p-value for annotation
p_value_text <- paste("P-value:", format(p_value, digits = 3))

# Plotting
ggplot(perm_data, aes(x = T_values)) +
  geom_histogram(aes(y = ..density..), binwidth = .001, fill = "skyblue1", color = "grey69") +
  geom_vline(data = original_coef_df, aes(xintercept = original_coef), color = "indianred4", linetype = "dashed", size = 1) +
  geom_density(alpha = .2, fill = "#FF6666") +
  labs(title = "PGS*EPDS_prenatal; Mild Dep (outcome)",
       x = "Permuted Coefficient",
       y = "Density") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        text = element_text(size = 16),
        plot.title = element_text(hjust = 0.5, size = 16)) +
  annotate("text", x = min(perm_coefs), y = max(density(perm_coefs)$y), label = p_value_text, hjust = 0, vjust = 1.5)







#----- Moderate Dep 24 years --------
XX.cont <- XX_female

glm.cont.old <- glm(FKDQ1010 ~  b023 + bestgest + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + 
                      PC7 + PC8 + PC9 + PC10 + ePRS_microglia * EPDS_mean_prenatal, data = XX.cont, family = binomial)
summary(glm.cont.old)
nrow_coef <- nrow(summary(glm.cont.old)$coef)

im.cont <- data.frame(influence.measures(glm.cont.old)$is.inf)
case.exclude.cont <- rownames(im.cont[which(im.cont[,nrow_coef] == "TRUE"),])
if (length(case.exclude.cont)>0){XX.cont <- XX.cont[-which(XX.cont$UniqueID %in% case.exclude.cont),]}

glm.cont     <- glm(FKDQ1010 ~ b023 + bestgest + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + 
                      PC7 + PC8 + PC9 + PC10 + ePRS_microglia * EPDS_mean_prenatal, data = XX.cont, family = binomial)

# Extract the original coefficient of interest
original_coef <- coef(summary(glm.cont))["ePRS_microglia:EPDS_mean_prenatal", "Estimate"]

# Set up number of permutations
n_perm <- 1000 # You can increase this for more robust results

# Initialize a vector to store the permuted coefficients
perm_coefs <- numeric(n_perm)

# Permutation test
set.seed(123) # For reproducibility
for(i in 1:n_perm) {
  # Shuffle the outcome
  shuffled_data <- XX.cont
  shuffled_data$FKDQ1010 <- sample(XX.cont$FKDQ1010)
  
  # Fit the model with the shuffled data
  perm_model <- glm(FKDQ1010 ~ b023 + bestgest + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + 
                      PC7 + PC8 + PC9 + PC10 + ePRS_microglia * EPDS_mean_prenatal, 
                    data = shuffled_data, family = binomial)
  
  # Store the coefficient of interest
  perm_coefs[i] <- coef(summary(perm_model))["ePRS_microglia:EPDS_mean_prenatal", "Estimate"]
}

# Calculate p-value
p_value <- mean(abs(perm_coefs) >= abs(original_coef))

# Create a data frame for plotting
perm_data <- data.frame(T_values = perm_coefs)

# Original coefficient as a data frame for annotation
original_coef_df <- data.frame(original_coef = original_coef)

# Calculate p-value for annotation
p_value_text <- paste("P-value:", format(p_value, digits = 3))

# Plotting
ggplot(perm_data, aes(x = T_values)) +
  geom_histogram(aes(y = ..density..), binwidth = .001, fill = "skyblue1", color = "grey69") +
  geom_vline(data = original_coef_df, aes(xintercept = original_coef), color = "indianred4", linetype = "dashed", size = 1) +
  geom_density(alpha = .2, fill = "#FF6666") +
  labs(title = "PGS*EPDS_prenatal; Moderate Dep (outcome)",
       x = "Permuted Coefficient",
       y = "Density") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        text = element_text(size = 16),
        plot.title = element_text(hjust = 0.5, size = 16)) +
  annotate("text", x = min(perm_coefs), y = max(density(perm_coefs)$y), label = p_value_text, hjust = 0, vjust = 1.5)



#----- Severe Dep 24 years --------
XX.cont <- XX_female

glm.cont.old <- glm(FKDQ1020 ~ b023 + bestgest + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + 
                      PC7 + PC8 + PC9 + PC10 + ePRS_microglia * EPDS_mean_prenatal, data = XX.cont, family = binomial)
summary(glm.cont.old)
nrow_coef <- nrow(summary(glm.cont.old)$coef)

im.cont <- data.frame(influence.measures(glm.cont.old)$is.inf)
case.exclude.cont <- rownames(im.cont[which(im.cont[,nrow_coef] == "TRUE"),])
if (length(case.exclude.cont)>0){XX.cont <- XX.cont[-which(XX.cont$UniqueID %in% case.exclude.cont),]}

glm.cont     <- glm(FKDQ1020 ~ bestgest + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + 
                      PC7 + PC8 + PC9 + PC10 + ePRS_microglia * EPDS_mean_prenatal, data = XX.cont, family = binomial)
# Extract the original coefficient of interest
original_coef <- coef(summary(glm.cont))["ePRS_microglia:EPDS_mean_prenatal", "Estimate"]

# Set up number of permutations
n_perm <- 1000 # You can increase this for more robust results

# Initialize a vector to store the permuted coefficients
perm_coefs <- numeric(n_perm)

# Permutation test
set.seed(123) # For reproducibility
for(i in 1:n_perm) {
  # Shuffle the outcome
  shuffled_data <- XX.cont
  shuffled_data$FKDQ1020 <- sample(XX.cont$FKDQ1020)
  
  # Fit the model with the shuffled data
  perm_model <- glm(FKDQ1020 ~ b023 + bestgest + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + 
                      PC7 + PC8 + PC9 + PC10 + ePRS_microglia * EPDS_mean_prenatal, 
                    data = shuffled_data, family = binomial)
  
  # Store the coefficient of interest
  perm_coefs[i] <- coef(summary(perm_model))["ePRS_microglia:EPDS_mean_prenatal", "Estimate"]
}

# Calculate p-value
p_value <- mean(abs(perm_coefs) >= abs(original_coef))

# Create a data frame for plotting
perm_data <- data.frame(T_values = perm_coefs)

# Original coefficient as a data frame for annotation
original_coef_df <- data.frame(original_coef = original_coef)

# Calculate p-value for annotation
p_value_text <- paste("P-value:", format(p_value, digits = 3))

# Plotting
ggplot(perm_data, aes(x = T_values)) +
  geom_histogram(aes(y = ..density..), binwidth = .001, fill = "skyblue1", color = "grey69") +
  geom_vline(data = original_coef_df, aes(xintercept = original_coef), color = "indianred4", linetype = "dashed", size = 1) +
  geom_density(alpha = .2, fill = "#FF6666") +
  labs(title = "PGS*EPDS_prenatal; Severe Dep (outcome)",
       x = "Permuted Coefficient",
       y = "Density") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        text = element_text(size = 16),
        plot.title = element_text(hjust = 0.5, size = 16)) +
  annotate("text", x = min(perm_coefs), y = max(density(perm_coefs)$y), label = p_value_text, hjust = 0, vjust = 1.5)













XX.cont <- XX_female


XX.cont <- XX.cont %>%
  rowwise() %>%
  mutate(Any_dep = as.integer(FKDQ1000 == 1 | FKDQ1010 == 1 | FKDQ1020 == 1)) %>%
  ungroup()


glm.cont.old <- glm(Any_dep ~ b023 + bestgest + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + 
                      PC7 + PC8 + PC9 + PC10 + ePRS_microglia * EPDS_mean_prenatal, data = XX.cont, family = binomial)
summary(glm.cont.old)
nrow_coef <- nrow(summary(glm.cont.old)$coef)

im.cont <- data.frame(influence.measures(glm.cont.old)$is.inf)
case.exclude.cont <- rownames(im.cont[which(im.cont[,nrow_coef] == "TRUE"),])
if (length(case.exclude.cont)>0){XX.cont <- XX.cont[-which(XX.cont$UniqueID %in% case.exclude.cont),]}

glm.cont     <- glm(Any_dep ~ bestgest + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + 
                      PC7 + PC8 + PC9 + PC10 + ePRS_microglia * EPDS_mean_prenatal, data = XX.cont, family = binomial)
# Extract the original coefficient of interest
original_coef <- coef(summary(glm.cont))["ePRS_microglia:EPDS_mean_prenatal", "Estimate"]

# Set up number of permutations
n_perm <- 1000 # You can increase this for more robust results

# Initialize a vector to store the permuted coefficients
perm_coefs <- numeric(n_perm)

# Permutation test
set.seed(123) # For reproducibility
for(i in 1:n_perm) {
  # Shuffle the outcome
  shuffled_data <- XX.cont
  shuffled_data$Any_dep <- sample(XX.cont$Any_dep)
  
  # Fit the model with the shuffled data
  perm_model <- glm(Any_dep ~ b023 + bestgest + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + 
                      PC7 + PC8 + PC9 + PC10 + ePRS_microglia * EPDS_mean_prenatal, 
                    data = shuffled_data, family = binomial)
  
  # Store the coefficient of interest
  perm_coefs[i] <- coef(summary(perm_model))["ePRS_microglia:EPDS_mean_prenatal", "Estimate"]
}

# Calculate p-value
p_value <- mean(abs(perm_coefs) >= abs(original_coef))

# Create a data frame for plotting
perm_data <- data.frame(T_values = perm_coefs)

# Original coefficient as a data frame for annotation
original_coef_df <- data.frame(original_coef = original_coef)

# Calculate p-value for annotation
p_value_text <- paste("P-value:", format(p_value, digits = 3))

# Plotting
plotme <- ggplot(perm_data, aes(x = T_values)) +
  geom_histogram(aes(y = ..density..), binwidth = .005, fill = "skyblue1", color = "grey69") +
  geom_vline(data = original_coef_df, aes(xintercept = original_coef), color = "indianred4", linetype = "dashed", size = 1) +
  geom_density(alpha = .2, fill = "#FF6666") +
  labs(title = "PGS*EPDS_prenatal; Depression (outcome)",
       x = "Permuted Coefficient",
       y = "Density") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        text = element_text(size = 16),
        plot.title = element_text(hjust = 0.5, size = 16)) +
  annotate("text", x = min(perm_coefs), y = max(density(perm_coefs)$y), label = p_value_text, hjust = 0, vjust = 1.5)


save(perm_data,p_value_text, original_coef_df,
     file = "/data1/meaneylab/eamon/Microglia_ePRS/Revisions_BBI/ALSPAC_perm_res.RDS")

load("/data1/meaneylab/eamon/Microglia_ePRS/Revisions_BBI/ALSPAC_perm_res.RDS")


ggsave("/data1/meaneylab/eamon/Microglia_ePRS/Revisions_BBI/ALSPAC_perm.pdf", 
       plot = plotme, width = 15, height = 12.5, units = "cm", bg="white")

