library(mr.raps)
library(data.table)
library(TwoSampleMR)
library(MRPRESSO)
library(devtools)

# List available GWASs in package
ao <- available_outcomes()

# Get instruments
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


exposure_dat <- format_data(placenta_trimmed, type="exposure")


# Get effects of instruments on outcome
# ieu-a-92 = Obesity
# ieu-a-1083 and ebi-a-GCST007557 = birth weight
# ieu-a-7 = coronoary heart disease
# ieu-a-798 = myocardial infarction
# ieu-a-2 = BMI
# ebi-a-GCST005902 = Broad dep
# ebi-a-GCST005903 = ICD10 MDD UKBB
# ieu-a-1188 = Wray et al 2018
# ieu-b-102 = Howard et al no 23 and me
# Howard et al use 23 and me data so only top 10,000 SNPs are available
# ieu-b-110 = LDL
# ieu-b-111 = Triglycerides
# ebi-a-GCST003837 = chronotype
# ukb-d-I50 = heart failure
# ieu-a-810 = schizophrenia
# Neurociticism || id:ebi-a-GCST006940
# Cerebrovascular disease || id:ukb-e-433_CSA

# ukb-b-11311 ukb-d-KRA_PSY_ANXIETY
# ukb-a-525 ieu-b-41    BP

CHD_outcome_dat <- extract_outcome_data(snps=exposure_dat$SNP, outcomes="ieu-a-7")
MI_outcome_dat <- extract_outcome_data(snps=exposure_dat$SNP, outcomes="ieu-a-798")
MDD_outcome_dat <- extract_outcome_data(snps=exposure_dat$SNP, outcomes="ieu-b-102")

# Harmonise the exposure and outcome data
CHD_dat <- harmonise_data(exposure_dat, CHD_outcome_dat)
MI_dat <- harmonise_data(exposure_dat, MI_outcome_dat)
MDD_dat <- harmonise_data(exposure_dat, MDD_outcome_dat)

# Perform MR
CHD_res <- mr(CHD_dat, method_list=c("mr_weighted_median", "mr_ivw"))
MI_res <- mr(MI_dat, method_list=c("mr_weighted_median", "mr_ivw"))
MDD_res <- mr(MDD_dat, method_list=c("mr_weighted_median", "mr_ivw"))


CHD_ORs <- generate_odds_ratios(CHD_res)
MDD_ORs <- generate_odds_ratios(MDD_res)
MI_ORs <- generate_odds_ratios(MI_res)

MDD$Study <- 'Nikpay et al 2015'
MDD_ORs$Disease <- 'Major Depression \n Howard et al 2019 \n (without 23 and me)'

Total <- rbind(CHD_ORs, MI_ORs, MDD_ORs)

TOT <- Total[Total$method %in% c("Inverse variance weighted", "Weighted median"), ]
write_csv(TOT, "/data1/meaneylab/eamon/MR_results_Sept_23rd.csv")
TOT <- read.csv("/data1/meaneylab/eamon/MR_results_Sept_23rd.csv")

### Forest plot
ggplot(TOT, aes(x=Disease, y=b, ymin=lo_ci, ymax=up_ci, col = method)) +
  geom_pointrange(position=position_dodge(width = 0.5)) +
  geom_hline(yintercept=0, lty=2) +  # add a dotted line at x=0 after flip
  coord_flip() +  # flip coordinates (puts labels on y axis)
  ylab("Estimate (95% CI)")  +
  theme_bw() + # use a white background
  theme(axis.text.y =element_text(hjust=0.5,vjust=0.5)) +
  scale_color_hue(l=40, c=40) +
  theme(legend.title=element_blank()) +
  theme(text = element_text(size=15)) +
  labs(x = "") 

# Pleiotropy and heterogeneity tests
CHD_egger <-mr_pleiotropy_test(CHD_dat)
MI_egger <- mr_pleiotropy_test(MI_dat)
MDD_egger <- mr_pleiotropy_test(MDD_dat)
Egger_intercept <- rbind(CHD_egger, MDD_egger, MI_egger)
write_csv(Egger_intercept, "MR_Egger_intercept_results.csv")

###############################################################

res_single <- mr_singlesnp(MDD_dat)
res_single <- res_single[-c(28), ]

res_loo <- mr_leaveoneout(MDD_dat)


p1 <- mr_scatter_plot(MDD_res, MDD_dat)
p1[[1]] + 
  theme(legend.title=element_blank()) +
  theme_bw() +
  theme(legend.position = "top") +
  theme(legend.title=element_blank()) +
  geom_hline(yintercept=0, color = "grey44", lty=2) +
  labs(y = "SNP effect on major depressive disorder") +
  labs(x = "SNP effect on placenta cyan module eQTLs") 
  

p2 <- mr_forest_plot(res_single)
p2[[1]] + theme_bw() + theme(legend.position="none") +
  labs(x = "MR single SNP effect size for Placenta eQTLs on major depressive disorder") 


p3 <- mr_leaveoneout_plot(res_loo)
p3[[1]] + theme_bw() + theme(legend.position="none") +
  labs(x = "MR leave one out; Outcome- major depressive disorder") 

#p4 <- mr_funnel_plot(res_single)
#p4[[1]]

# MR PRESSO
MR_PRESSO <- run_mr_presso(dat, NbDistribution = 1000, SignifThreshold = 0.05)

# Work in progress
MR_RAPS <- mr(dat, method_list = c("mr_raps"))
MR_RAPS <- mr(dat, method_list = c("mr_raps"), parameters = list(over.dispersion = FALSE, loss.function = "l2"))

######### Create Report ##############
getwd()
mr_report(
  MDD_dat,
  output_path = ".",
  output_type = "pdf",
  author = "Analyst",
  study = "Two Sample MR",
  path = system.file("reports", package = "TwoSampleMR"))

#######################
ao<-available_outcomes()
ao<-ao[ao$author=="Elliott LT",] #identify diseases
#ao<-ao[which(ao$ncase>100),]

dis_dat <- extract_outcome_data(
  snps = exposure_dat$SNP,
  outcomes = ao$id
)

dat3 <- harmonise_data(
  exposure_dat = exposure_dat, 
  outcome_dat = dis_dat
)


res<-mr(dat3,method_list=c("mr_wald_ratio","mr_ivw"))
write_csv(res, "placenta_inflam_full_MR_pheWAS_results.csv")

res<-split_outcome(res) # to keep the Y axis label clean we exclude the exposure ID labels from the exposure column 

res<-sort_1_to_many(res,b="b",sort_action=4) #this sorts results by decreasing effect size (largest effect at top of the plot)

res1<-res[1:52,]
res2<-res[53:138,]

plot1<-forest_plot_1_to_many(res1,b="b",se="se",
                             exponentiate=T,trans="log2",ao_slc=F,lo=0.004,
                             up=461,col1_width=2,TraitM="outcome",
                             col_text_size=3,xlab="")

plot2<-forest_plot_1_to_many(res2,b="b",se="se",
                             exponentiate=T,trans="log2",ao_slc=F,lo=0.004,
                             up=461,col1_width=2,TraitM="outcome",
                             subheading_size=11,col_text_size=3,xlab="")

plot1
plot2


