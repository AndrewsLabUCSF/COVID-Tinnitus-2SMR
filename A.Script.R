#Packages
library(tidyverse)
library(TwoSampleMR)
library(RadialMR)
library(phenoscanner)
library(gt)
#Formatting Exposure Data
exposure_data <- read_tsv("hgi_A.tsv.gz", 
                          comment = "#",)
exposure_formatted <- format_data(
  exposure_data,
  type = "exposure",
  snps = NULL,
  header = TRUE,
  phenotype_col = "TRAIT",
  snp_col = "DBSNP_ID",
  beta_col = "BETA",
  se_col = "SE",
  eaf_col = "AF",
  effect_allele_col = "REF",
  other_allele_col = "ALT",
  pval_col = "P",  
  samplesize_col = "N",
  chr_col = "CHROM",
  pos_col = "POS")
#Clumping Exposure Data
exposure_clumped <- exposure_formatted %>% 
  filter(pval.exposure < 1e-6) %>% 
  clump_data(., )
#Formatting Outcome Data
outcome_data <- read_tsv("bhatt.tsv.gz")
outcome_formatted <- format_data(
  outcome_data,
  type = "outcome",
  snps = NULL,
  header = TRUE,
  snp_col = "SNP",
  beta_col = "BETA",
  se_col = "STDERR_BETA",
  eaf_col = "ALT_ALLELE_FREQ",
  effect_allele_col = "A1",
  other_allele_col = "A2",
  pval_col = "P",
  chr_col = "CHR",
  pos_col = "BP")
#Extracting/"Clumping" Outcome Data , I Think ?
outcome_clumped <- filter(outcome_formatted, SNP %in% exposure_clumped$SNP)

#Harmonization
harmony_data <- harmonise_data(
  exposure_dat = exposure_clumped,
  outcome_dat = outcome_clumped
) 
#This space will be used for filtering if I wanted to do that

#Performing MR
results <- mr(harmony_data, method_list = c("mr_egger_regression", 
                                            "mr_weighted_median", 
                                            "mr_ivw_fe", 
                                            "mr_weighted_mode"))
#Setting Plot Theme
theme_set(theme_bw())
#Scatter Plot
scatter <- mr_scatter_plot(results, harmony_data)
scatter[[1]]
#Forest Plot
trees <- mr_singlesnp(harmony_data, all_method = c("mr_egger_regression", 
                                                   "mr_weighted_median", 
                                                   "mr_ivw_fe", 
                                                   "mr_weighted_mode"))
forest <- mr_forest_plot(trees)
forest[[1]]
#Leave-One-Out Plot
one_out <- mr_leaveoneout(harmony_data)
leave <- mr_leaveoneout_plot(one_out)
leave[[1]] 
#Funnel Plot
funnel_cakes <- mr_singlesnp(harmony_data, all_method = c("mr_egger_regression", 
                                                          "mr_weighted_median", 
                                                          "mr_ivw_fe", 
                                                          "mr_weighted_mode"))
funnel <- mr_funnel_plot(funnel_cakes) 
funnel[[1]]
#Heterogeneity and Horizontal Pleiotropy Tests
mr_heterogeneity(harmony_data)
mr_pleiotropy_test(harmony_data)
#Radial Formatting
radial_data <- harmony_data %>% filter(mr_keep == T) %>% dat_to_RadialMR()
radial_mr_data <- radial_data$COVID_A2__EUR
#Radial IVW
bonff = 0.05/nrow(radial_mr_data) #Bonferonni Correction for reducing Type I Errors
radial_ivw_res <- ivw_radial(radial_mr_data, alpha = bonff)
#Radial Egger
radial_egger_res <- egger_radial(radial_mr_data, alpha = bonff)
#Radial Plots
ivw_radial_p <- plot_radial(radial_ivw_res, radial_scale = F, show_outliers = F)
egger_radial_p <- plot_radial(radial_egger_res, radial_scale = F, show_outliers = F)

cowplot::plot_grid(
  ivw_radial_p + coord_fixed(ratio=0.25) + theme(legend.position = 'bottom'), 
  egger_radial_p + theme(legend.position = 'bottom'), 
  align = 'h'
)
#The End .