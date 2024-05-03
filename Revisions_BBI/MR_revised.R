library(R.utils)
library(mr.raps)
library(data.table)
library(TwoSampleMR)
library(MRPRESSO)
library(devtools)
library(readr)
library(tidyverse)
library(paletteer)
library(NatParksPalettes)

# Get instruments
Micro_eQTLs <- fread("/data1/meaneylab/eamon/Microglia_ePRS/New_stuff/MR/Micro_ACC_eQTLs.txt")%>%
  dplyr::rename(effect_allele=alt)%>%
  dplyr::rename(other_allele=ref)%>%
  dplyr::rename(eaf=maf)%>%
  dplyr::rename(SNP=rsID)%>%
  dplyr::rename(beta=slope)%>%
  dplyr::rename(se=slope_se)%>%
  dplyr::rename(pval=pval_nominal)%>%
  dplyr::rename(phenotype=gene_id)
exposure_dat <- format_data(Micro_eQTLs, type="exposure")

Female_dep <- read_outcome_data(
  snps = exposure_dat$SNP,
  filename = "/n9/cortexjobs/common/resources/gwas_database/original/mdd/Sex-specific_MDD_2019/GWAS.Broad_Females.SNPtest_FULL.txt",
  snp_col = "rsid",
  beta_col = "frequentist_add_beta_1:add/mdd_broad=1",
  se_col = "frequentist_add_se_1",
  effect_allele_col = "alleleA",
  other_allele_col = "alleleB",
  eaf_col = "all_maf",
  pval_col = "frequentist_add_wald_pvalue_1",
  units_col = FALSE,
  gene_col = FALSE,
  samplesize_col = FALSE
)
Female_dep <- harmonise_data(exposure_dat, Female_dep)
Female_dep_res <- mr(Female_dep, method_list=c("mr_weighted_median", "mr_ivw", "mr_weighted_mode", "mr_egger_regression", "mr_raps"))
Female_dep_res$GWAS <- 'Female_depression'

Female_dep_pleio <-mr_pleiotropy_test(Female_dep)
Female_dep_het <- mr_heterogeneity(Female_dep)
Female_dep_PRESSO <- run_mr_presso(Female_dep, NbDistribution = 10000, SignifThreshold = 0.05)
Female_dep_PRESSO$GWAS <- 'Female_depression'

Female_dep_PRESSO_res <- Female_dep_PRESSO[[1]][["Main MR results"]]
Female_dep_PRESSO_res$GWAS <- 'Female_depression'

Male_dep <- read_outcome_data(
  snps = exposure_dat$SNP,
  filename = "/n9/cortexjobs/common/resources/gwas_database/original/mdd/Sex-specific_MDD_2019/GWAS.Broad_Males.SNPtest_FULL.txt",
  sep = " ",
  snp_col = "rsid",
  beta_col = "frequentist_add_beta_1:add/mdd_broad=1",
  se_col = "frequentist_add_se_1",
  effect_allele_col = "alleleA",
  other_allele_col = "alleleB",
  eaf_col = "all_maf",
  pval_col = "frequentist_add_wald_pvalue_1",
  units_col = FALSE,
  gene_col = FALSE,
  samplesize_col = FALSE
)
Male_dep <- harmonise_data(exposure_dat, Male_dep)
Male_dep_res <- mr(Male_dep, method_list=c("mr_weighted_median", "mr_ivw", "mr_weighted_mode", "mr_egger_regression", "mr_raps"))
Male_dep_res$GWAS <- 'Male_depression'

Male_dep_pleio <-mr_pleiotropy_test(Male_dep)
Male_dep_het <- mr_heterogeneity(Male_dep)
Male_dep_PRESSO <- run_mr_presso(Male_dep, NbDistribution = 10000, SignifThreshold = 0.05)
Male_dep_PRESSO$GWAS <- 'Male_depression'

Male_dep_PRESSO_res <- Male_dep_PRESSO[[1]][["Main MR results"]]
Male_dep_PRESSO_res$GWAS <- 'Male_depression'


MDD_full <- extract_outcome_data(snps=exposure_dat$SNP, outcomes="ieu-b-102")
MDD_full <- harmonise_data(exposure_dat, MDD_full)
MDD_full_res <- mr(MDD_full, method_list=c("mr_weighted_median", "mr_ivw", "mr_weighted_mode", "mr_egger_regression", "mr_raps"))
MDD_full_res$GWAS <- 'Combined_depression'

MDD_full_pleio <-mr_pleiotropy_test(MDD_full)
MDD_full_het <- mr_heterogeneity(MDD_full)
MDD_full_PRESSO <- run_mr_presso(MDD_full, NbDistribution = 10000, SignifThreshold = 0.05)
MDD_full_PRESSO$GWAS <- 'Combined_depression'

MDD_full_PRESSO_res <- MDD_full_PRESSO[[1]][["Main MR results"]]
MDD_full_PRESSO_res$GWAS <- 'Combined_depression'

PRESSO_full <- rbind(MDD_full_PRESSO_res,
                     Female_dep_PRESSO_res,
                     Male_dep_PRESSO_res)

write.csv(PRESSO_full, "/data1/meaneylab/eamon/Microglia_ePRS/Revisions_BBI/PRESSO_results.csv")


MDD_full_PRESSO[[1]][["MR-PRESSO results"]]$`Distortion Test`
# P=0.84
Male_dep_PRESSO[[1]][["MR-PRESSO results"]]$`Distortion Test`
#No outliers, so not applicable
Female_dep_PRESSO[[1]][["MR-PRESSO results"]]$`Distortion Test`
# P=0.066

Total <- rbind(Female_dep_res, Male_dep_res, MDD_full_res)
write.csv(Total, "/data1/meaneylab/eamon/Microglia_ePRS/Revisions_BBI/MR_results_rev.csv")



Total <- read.csv("/data1/meaneylab/eamon/Microglia_ePRS/Revisions_BBI/MR_results_rev.csv")

Total$GWAS <- gsub('Combined_depression','Combined_GWAS',Total$GWAS)
Total$GWAS <- gsub('Female_depression','Female_GWAS',Total$GWAS)
Total$GWAS <- gsub('Male_depression','Male_GWAS',Total$GWAS)

MR <- ggplot(Total, aes(x=method, y= b, ymin= b-(se*1.96), ymax= b+(se*1.96), col = method)) +
  geom_pointrange(position=position_dodge(width = 0.75)) +
  geom_hline(yintercept=0, lty=2) + 
  facet_wrap(~GWAS)+
  coord_flip() +  
  scale_color_manual(values=natparks.pals("CapitolReef"))+
  ylab("Estimate (95% CI)")  +
  theme_bw() +
  xlab("")+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = 'none',
        strip.background = element_blank(),
        strip.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 16),
        axis.title.x = element_text(size = 20),
        axis.text.x = element_text(size=16))

ggsave("/data1/meaneylab/eamon/Microglia_ePRS/Revisions_BBI/MR_revised.pdf", 
       plot = MR, width = 40, height = 10, units = "cm", bg="white")


#----------------- LOO ------------------------------

Female_loo <- mr_leaveoneout(Female_dep)
p1 <- mr_leaveoneout_plot(Female_loo)
p1 <- p1$cISh3z.8ZnS9x + xlab('Leave One Out: effect of eQTLs on Female MDD') +
  theme_bw() +
  theme(legend.position = 'none',
        text = element_text(size = 14))


Male_loo <- mr_leaveoneout(Male_dep)
p2 <- mr_leaveoneout_plot(Male_loo)
p2 <- p2$cISh3z.gQVELn + xlab('Leave One Out: effect of eQTLs on Male MDD') +
  theme_bw() +
  theme(legend.position = 'none',
        text = element_text(size = 14))

MDD_loo <- mr_leaveoneout(MDD_full)
p3 <- mr_leaveoneout_plot(MDD_loo)
p3 <- p3$`cISh3z.ieu-b-102` + xlab('Leave One Out: effect of eQTLs on Combined MDD') +
  theme_bw() +
  theme(legend.position = 'none',
        text = element_text(size = 14))

plotme <- cowplot::plot_grid(p3,p1,p2,nrow = 1,align = 'h',labels = 'AUTO')
ggsave("/data1/meaneylab/eamon/Microglia_ePRS/Revisions_BBI/LOO.pdf", 
       plot = plotme, width = 45, height = 30, units = "cm", bg="white")

