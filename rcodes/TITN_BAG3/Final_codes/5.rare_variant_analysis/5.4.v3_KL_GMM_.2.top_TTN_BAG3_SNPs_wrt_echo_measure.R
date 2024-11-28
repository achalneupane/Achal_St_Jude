rm(list = ls())

library(dplyr)
library(haven)
library(lubridate)
library(stringr)
library(rlang)
library(lme4)
library(ggplot2)
library(ggpubr)
source("../utilities.R")
## Load echo data
load('data/EUR_common_variants.RData')
load('data/AFR_common_variants.RData')

diagDT <- readRDS("diagDT.rds") %>%
  rename(sjlid = IID)

# table(diagDT$era_variable)
# 
# (1968,1978] (1978,1988] (1988,1998] (1998,2008] [1958,1968] 
# 6635       12076       11101        3831         196

################################
## treatment era from SJLIFE
################################
sjlife_path <- "Z:/SJShare/SJCOMMON/ECC/SJLife/SJLIFE Data Freeze/2 Final Data SJLIFE/20200430/"

diagnosis <- read_sas(paste0(sjlife_path, "Clinical data/diagnosis.sas7bdat")) %>%
  rename_with(tolower) %>%
  select(sjlid, primdx, diagdt) %>%
  filter(primdx == 1)

sjlife_echo_data <- readRDS("results/sjlife_echo_analysis_dat.rds")

overlap_cohort <- sjlife_echo_data %>% filter(ccss %in% c("Expansion", "Original"))
sjonly_cohort <- sjlife_echo_data %>% filter(!sjlid %in% c(overlap_cohort$sjlid))

EUR_sjonly <- EUR_common_variants %>%
  rename_with(tolower) %>%
  filter(iid %in% sjonly_cohort$sjlid) %>%
  select(iid, pc1:pc10, rs3829747:rs3858340) %>%
  unique() %>%  ## N = 374
  mutate(sjlid = iid, ccssid = NA)


EUR_overlap_1 <- EUR_common_variants %>%
  rename_with(tolower) %>%
  filter(iid %in% overlap_cohort$sjlid) %>%
  select(iid, pc1:pc10, rs3829747:rs3858340) %>%
  unique() %>%  ## N = 1233
  mutate(sjlid = iid) %>%
  left_join(overlap_cohort %>% select(sjlid, ccssid), by = "sjlid")

EUR_overlap_2 <- EUR_common_variants %>%
  rename_with(tolower) %>%
  filter(iid %in% overlap_cohort$ccssid) %>%
  select(iid, pc1:pc10, rs3829747:rs3858340) %>%
  unique() %>%  ## N = 1233
  mutate(ccssid = iid) %>%
  left_join(overlap_cohort %>% select(sjlid, ccssid), by = "ccssid")

EUR_common_variant_data <- bind_rows(EUR_sjonly, EUR_overlap_1, EUR_overlap_2)  ## N = 3298
saveRDS(EUR_common_variant_data, "EUR_common_variant_data_clean.rds")


AFR_sjonly <- AFR_common_variants %>%
  rename_with(tolower) %>%
  filter(iid %in% sjonly_cohort$sjlid) %>%
  select(iid, pc1:pc10, rs3829747:rs3858340) %>%
  unique() %>%  ## N = 74
  mutate(sjlid = iid, ccssid = NA)


AFR_overlap_1 <- AFR_common_variants %>%
  rename_with(tolower) %>%
  filter(iid %in% overlap_cohort$sjlid) %>%
  select(iid, pc1:pc10, rs3829747:rs3858340) %>%
  unique() %>%  ## N = 401
  mutate(sjlid = iid) %>%
  left_join(overlap_cohort %>% select(sjlid, ccssid), by = "sjlid")

AFR_overlap_2 <- AFR_common_variants %>%
  rename_with(tolower) %>%
  filter(iid %in% overlap_cohort$ccssid) %>%
  select(iid, pc1:pc10, rs3829747:rs3858340) %>%
  unique() %>%  ## N = 0
  mutate(ccssid = iid) %>%
  left_join(overlap_cohort %>% select(sjlid, ccssid), by = "ccssid")

AFR_common_variant_data <- bind_rows(AFR_sjonly, AFR_overlap_1, AFR_overlap_2)  ## N = 475
saveRDS(AFR_common_variant_data, "AFR_common_variant_data_clean.rds")


EUR_common_variants2 <- EUR_common_variant_data %>%
  mutate(rs3829746_yn = as.numeric(rs3829746 > 0),
         rs2234962_yn = as.numeric(rs2234962 > 0))

AFR_common_variants2 <- AFR_common_variant_data %>%
  mutate(rs3829746_yn = as.numeric(rs3829746 > 0),
         rs2234962_yn = as.numeric(rs2234962 > 0))
## Need to merge echo data for SJLIFE EUR and SJLIFE AFR cohorts and we only want to analyze for:
# c("LV_Ejection_Fraction_3D", "LV_End_Diastolic_Volume_3D", "LV_End_Systolic_Volume_3D",
# "LV_Stroke_Volume_3D", "LVMassMM_Index", "LV_GLPS_AVG", "LV_Relative_Wall_Thickness")

# Clean the outliers
# Some EF values are in fractions; change them to percentages
sjlife_echo_data2 <- sjlife_echo_data %>%
  mutate(lv_ejection_fraction_3d = ifelse(lv_ejection_fraction_3d < 1, 
                                          lv_ejection_fraction_3d * 100,
                                          lv_ejection_fraction_3d),
         lv_end_diastolic_volume_3d = ifelse(lv_end_diastolic_volume_3d < 10,
                                             NA,
                                             lv_end_diastolic_volume_3d),
         lv_end_systolic_volume_3d = ifelse(lv_end_systolic_volume_3d < 10,
                                            10,
                                            lv_end_systolic_volume_3d),
         lv_stroke_volume_3d = ifelse(lv_stroke_volume_3d < 10,
                                      NA,
                                      lv_stroke_volume_3d),
         lvmassmm_index = ifelse(lvmassmm_index < 10,
                                 NA,
                                 lvmassmm_index)) %>%
  mutate(male = as.numeric(gender == "Male"),
         agevisit = yrdif(dob, studydatetime),
         agedx = yrdif(dob, diagdt))


EUR_data3 <- inner_join(sjlife_echo_data2, 
                        EUR_common_variants2 %>% select(-ccssid),
                        by = "sjlid")
AFR_data3 <- AFR_common_variants2 %>% inner_join(sjlife_echo_data2, by = "sjlid")


## Visualize the data
library(ggplot2)
library(rstatix)
EUR_data4 <- EUR_data3 %>%
  filter(!is.na(rs2234962) & !is.na(rs3829746)) %>%
  filter(visittype %in% paste0("SJLIFE Visit ", 1:7)) %>%
  filter(anthracyclines_5 == 1 | heartavg > 0) %>%
  left_join(diagDT %>% select(sjlid, era_variable), by = "sjlid") %>%
  mutate(era_variable = cut(year(diagdt), c(seq(1958, 1998, by = 10), 2009),
                            include.lowest = T, right = T))



AFR_data4 <- AFR_data3 %>%
  filter(!is.na(rs2234962) & !is.na(rs3829746)) %>%
  filter(visittype %in% paste0("SJLIFE Visit ", 1:7)) %>%
  filter(anthracyclines_5 == 1 | heartavg > 0) %>%
  mutate(era_variable = cut(year(diagdt), c(seq(1958, 1998, by = 10), 2009),
                            include.lowest = T, right = T))




param <- c("lv_ejection_fraction_3d", "lv_end_diastolic_volume_3d",    
           "lv_end_systolic_volume_3d", "lv_stroke_volume_3d",          
           "lvmassmm_index", "lv_glps_avg", "lv_relative_wall_thickness")
param_label <- c("LV ejection fraction", "LV end-diastolic volume",
                 "LV end-systolic volume", "LV stroke volume",
                 "LV mass index", "LV average global longitudinal strain",
                 "LV relative wall fitness")
snps <- c('rs3829746', 'rs2234962')
# dat = eur_baseline
res_EUR <- NULL


for (p in 1:length(param)){
  echo_p  <- param[p]
  for (s in (1:2)){
    snp_s <- snps[s]
    
    ## get the first visit where the measure of interest is non-missing
    dat <- EUR_data4 %>%
      filter(!is.na(!!sym(echo_p))) %>%
      arrange(sjlid, studydatetime) 
    
    fml <- as.formula(paste0(echo_p, " ~ ", snp_s, " + agedx + agevisit + male + anthracyclines_dose_5 + heartavg +
               pc1 + pc2 + pc3 + pc4 + pc5 + pc6 + pc7 + pc8 + pc9 + pc10 + era_variable + (1 | sjlid)"))
    fit <- lmer(formula = fml, data = dat)
    n <- as.numeric(summary(fit)$ngrps)
    beta <- summary(fit)$coef[2,1]
    se <- summary(fit)$coef[2,2]
    pval <- pchisq(beta ^ 2 / se ^ 2, 1, lower.tail = F)
    res_EUR <- rbind(res_EUR, data.frame(parameter = echo_p, snp = snp_s, n = n, beta, se, pval))
  }
}
print(res_EUR)
res_EUR$pval.adj <- p.adjust(res_EUR$pval, method = "BH")
res_EUR.ttn <- res_EUR %>% filter(snp == "rs3829746")
res_EUR.bag3 <- res_EUR %>% filter(snp == "rs2234962")

write.csv(res_EUR.ttn, "results/lmer_ECHO_param_common_variants_EUR_TTN.csv", row.names = F)
write.csv(res_EUR.bag3, "results/lmer_ECHO_param_common_variants_EUR_BAG3.csv", row.names = F)

## Fit regression models for each visit (AFR)
res_AFR <- NULL


for (p in 1:length(param)){
  echo_p  <- param[p]
  for (s in (1:2)){
    snp_s <- snps[s]
    
    ## get the first visit where the measure of interest is non-missing
    dat <- AFR_data4 %>%
      filter(!is.na(!!sym(echo_p))) %>%
      arrange(sjlid, studydatetime) 
    
    fml <- as.formula(paste0(echo_p, " ~ ", snp_s, " + agedx + agevisit + male + anthracyclines_dose_5 + heartavg +
               pc1 + pc2 + pc3 + pc4 + pc5 + pc6 + pc7 + pc8 + pc9 + pc10 + era_variable + (1 | sjlid)"))
    fit <- lmer(formula = fml, data=dat)
    n <- as.numeric(summary(fit)$ngrps)
    beta <- summary(fit)$coef[2,1]
    se <- summary(fit)$coef[2,2]
    pval <- pchisq(beta ^ 2 / se ^ 2, 1, lower.tail = F)
    res_AFR <- rbind(res_AFR, data.frame(parameter = echo_p, snp = snp_s, n = n, beta, se, pval))
  }
}
print(res_AFR)
res_AFR$pval.adj <- p.adjust(res_AFR$pval, method = "BH")
res_AFR.ttn <- res_AFR %>% filter(snp == "rs3829746")
res_AFR.bag3 <- res_AFR %>% filter(snp == "rs2234962")

write.csv(res_AFR.ttn, "results/lmer_ECHO_param_common_variants_AFR_TTN.csv", row.names = F)
write.csv(res_AFR.bag3, "results/lmer_ECHO_param_common_variants_AFR_BAG3.csv", row.names = F)

##########
## Male ##
##########
## Baseline
EUR_male_data <- EUR_data4 %>% filter(male == 1)
AFR_male_data <- AFR_data4 %>% filter(male == 1)

## Fit regression models for each visit (EUR)
res_EUR_male <- NULL


for (p in 1:length(param)){
  echo_p  <- param[p]
  for (s in (1:2)){
    snp_s <- snps[s]
    
    ## get the first visit where the measure of interest is non-missing
    dat <- EUR_male_data %>%
      filter(!is.na(!!sym(echo_p))) %>%
      arrange(sjlid, studydatetime) 
    
    fml <- as.formula(paste0(echo_p, " ~ ", snp_s, " + agedx + agevisit + anthracyclines_dose_5 + heartavg +
               pc1 + pc2 + pc3 + pc4 + pc5 + pc6 + pc7 + pc8 + pc9 + pc10 + era_variable + (1 | sjlid)"))
    fit <- lmer(formula = fml, data=dat)
    n <- as.numeric(summary(fit)$ngrps)
    beta <- summary(fit)$coef[2,1]
    se <- summary(fit)$coef[2,2]
    pval <- pchisq(beta ^ 2 / se ^ 2, 1, lower.tail = F)
    res_EUR_male <- rbind(res_EUR_male, data.frame(parameter = echo_p, snp = snp_s, n = n, beta, se, pval))
  }
}
print(res_EUR_male)
res_EUR_male$pval.adj <- p.adjust(res_EUR_male$pval, method = "BH")
res_EUR_male.ttn <- res_EUR_male %>% filter(snp == "rs3829746")
res_EUR_male.bag3 <- res_EUR_male %>% filter(snp == "rs2234962")

write.csv(res_EUR_male.ttn, "results/lmer_ECHO_param_common_variants_EUR_male_TTN.csv", row.names = F)
write.csv(res_EUR_male.bag3, "results/lmer_ECHO_param_common_variants_EUR_male_BAG3.csv", row.names = F)



## Fit regression models for each visit (AFR)
res_AFR_male <- NULL


for (p in 1:length(param)){
  echo_p  <- param[p]
  for (s in (1:2)){
    snp_s <- snps[s]
    
    ## get the first visit where the measure of interest is non-missing
    dat <- AFR_male_data %>%
      filter(!is.na(!!sym(echo_p))) %>%
      arrange(sjlid, studydatetime) 
    
    fml <- as.formula(paste0(echo_p, " ~ ", snp_s, " + agedx + agevisit + anthracyclines_dose_5 + heartavg +
               pc1 + pc2 + pc3 + pc4 + pc5 + pc6 + pc7 + pc8 + pc9 + pc10 + era_variable + (1 | sjlid)"))
    fit <- lmer(formula = fml, data = dat)
    n <- as.numeric(summary(fit)$ngrps)
    beta <- summary(fit)$coef[2,1]
    se <- summary(fit)$coef[2,2]
    pval <- pchisq(beta ^ 2 / se ^ 2, 1, lower.tail = F)
    res_AFR_male <- rbind(res_AFR_male, data.frame(parameter = echo_p, snp = snp_s, n = n, beta, se, pval))
  }
}
print(res_AFR_male)
res_AFR_male$pval.adj <- p.adjust(res_AFR_male$pval, method = "BH")
res_AFR_male.ttn <- res_AFR_male %>% filter(snp == "rs3829746")
res_AFR_male.bag3 <- res_AFR_male %>% filter(snp == "rs2234962")

write.csv(res_AFR_male.ttn, "results/lmer_ECHO_param_common_variants_AFR_male_TTN.csv", row.names = F)
write.csv(res_AFR_male.bag3, "results/lmer_ECHO_param_common_variants_AFR_male_BAG3.csv", row.names = F)




############
## Female ##
############
EUR_female_data <- EUR_data4 %>% filter(male == 0)
AFR_female_data <- AFR_data4 %>% filter(male == 0)

## Fit regression models for each visit (EUR)
res_EUR_female <- NULL


for (p in 1:length(param)){
  echo_p  <- param[p]
  for (s in (1:2)){
    snp_s <- snps[s]
    
    ## get the first visit where the measure of interest is non-missing
    dat <- EUR_female_data %>%
      filter(!is.na(!!sym(echo_p))) %>%
      arrange(sjlid, studydatetime) 
    
    fml <- as.formula(paste0(echo_p, " ~ ", snp_s, " + agedx + agevisit + anthracyclines_dose_5 + heartavg +
               pc1 + pc2 + pc3 + pc4 + pc5 + pc6 + pc7 + pc8 + pc9 + pc10 + era_variable + (1 | sjlid)"))
    fit <- lmer(formula = fml, data=dat)
    n <- as.numeric(summary(fit)$ngrps)
    beta <- summary(fit)$coef[2,1]
    se <- summary(fit)$coef[2,2]
    pval <- pchisq(beta ^ 2 / se ^ 2, 1, lower.tail = F)
    res_EUR_female <- rbind(res_EUR_female, data.frame(parameter = echo_p, snp = snp_s, n = n, beta, se, pval))
  }
}
print(res_EUR_female)
res_EUR_female$pval.adj <- p.adjust(res_EUR_female$pval, method = "BH")
res_EUR_female.ttn <- res_EUR_female %>% filter(snp == "rs3829746")
res_EUR_female.bag3 <- res_EUR_female %>% filter(snp == "rs2234962")

write.csv(res_EUR_female.ttn, "results/lmer_ECHO_param_common_variants_EUR_female_TTN.csv", row.names = F)
write.csv(res_EUR_female.bag3, "results/lmer_ECHO_param_common_variants_EUR_female_BAG3.csv", row.names = F)



## Fit regression models for each visit (AFR)
res_AFR_female <- NULL


for (p in 1:length(param)){
  echo_p  <- param[p]
  for (s in (1:2)){
    snp_s <- snps[s]
    
    ## get the first visit where the measure of interest is non-missing
    dat <- AFR_female_data %>%
      filter(!is.na(!!sym(echo_p))) %>%
      arrange(sjlid, studydatetime) 
    
    fml <- as.formula(paste0(echo_p, " ~ ", snp_s, " + agedx + agevisit + anthracyclines_dose_5 + heartavg +
               pc1 + pc2 + pc3 + pc4 + pc5 + pc6 + pc7 + pc8 + pc9 + pc10 + era_variable + (1 | sjlid)"))
    fit <- lmer(formula = fml, data=dat)
    n <- as.numeric(summary(fit)$ngrps)
    beta <- summary(fit)$coef[2,1]
    se <- summary(fit)$coef[2,2]
    pval <- pchisq(beta ^ 2 / se ^ 2, 1, lower.tail = F)
    res_AFR_female <- rbind(res_AFR_female, data.frame(parameter = echo_p, snp = snp_s, n = n, beta, se, pval))
  }
}
print(res_AFR_female)
res_AFR_female$pval.adj <- p.adjust(res_AFR_female$pval, method = "BH")
res_AFR_female.ttn <- res_AFR_female %>% filter(snp == "rs3829746")
res_AFR_female.bag3 <- res_AFR_female %>% filter(snp == "rs2234962")

write.csv(res_AFR_female.ttn, "results/lmer_ECHO_param_common_variants_AFR_female_TTN.csv", row.names = F)
write.csv(res_AFR_female.bag3, "results/lmer_ECHO_param_common_variants_AFR_female_BAG3.csv", row.names = F)


# # Regression for change in echo parameters
# library(dplyr)
# p = 'LV_Ejection_Fraction_3D'
# newdat = EUR_common_variants %>% group_by(sjlid) %>% mutate(diff=LV_Ejection_Fraction_3D-lag(LV_Ejection_Fraction_3D, default=first(LV_Ejection_Fraction_3D)))
