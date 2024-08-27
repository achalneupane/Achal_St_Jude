####################################
## Function to create demographic ##
####################################
# Function to create demographic table (Table 1)
get_demographic <- function(df, n){
## Ejection fraction
ejection_fraction_hrt <- paste0(round((median(df$ejection_fraction_hrt, na.rm = T)*100), 1), " (", round((quantile(df$ejection_fraction_hrt, prob=c(.25,.5,.75), type=1, na.rm = T)*100)[1], 1), "-" , round((quantile(df$ejection_fraction_hrt, prob=c(.25,.5,.75), type=1, na.rm = T)*100)[3],1), ")")

## Elapsed time
df$elapsedAGE <- df$agevent - df$agedx

# All elapsed negative age are converted to zero. ageevent probably used floor value so there are some negatives
df$elapsedAGE[df$elapsedAGE < 0] <- 0
df_elapsedAGE <- paste0(round(median(df$elapsedAGE,  na.rm = T), 1), " (", round(quantile(df$elapsedAGE, prob=c(.25,.5,.75), type=1, na.rm = T)[1], 1), "-" , round(quantile(df$elapsedAGE, prob=c(.25,.5,.75), type=1,  na.rm = T)[3],1), ")")



## Age at last contact: Median (IQR)
agelstcontact <- paste0(round(median(df$agelstcontact, na.rm = T), 1), " (", round(quantile(df$agelstcontact, prob=c(.25,.5,.75), type=1, na.rm = T)[1], 1), "-" , round(quantile(df$agelstcontact, prob=c(.25,.5,.75), type=1, na.rm = T)[3],1), ")")

## Age at diagnosis
agedx <- paste0(round(median(df$agedx, na.rm = T), 1), " (", round(quantile(df$agedx, prob=c(.25,.5,.75), type=1, na.rm = T)[1], 1), "-" , round(quantile(df$agedx, prob=c(.25,.5,.75), type=1, na.rm = T)[3],1), ")")

## Gender (0 = Male; 1 = Female)
male <- sum(df$gender == 0, na.rm = T)
male <- paste0(male, " (", round((male/n)*100,1), "%)")

female <- sum(df$gender == 1, na.rm = T)
female <- paste0(female, " (", round((female/n)*100,1), "%)")

## Anthracycline dose median (IQR)
jco_actual_dose <- df$anthra_jco_dose_any[df$anthra_jco_dose_any > 0] # keep only non-zeros
anthra_median_range <- paste0(round(median(jco_actual_dose, na.rm = T), 1), " (", round(quantile(jco_actual_dose, prob=c(.25,.5,.75), type=1, na.rm = T)[1], 1), "-" , round(quantile(jco_actual_dose, prob=c(.25,.5,.75), type=1, na.rm = T)[3],1), ")")

# ## Anthracycline dose median (range)
# jco_actual_dose <- df$anthra_jco_dose_any[df$anthra_jco_dose_any > 0] # keep only non-zeros
# anthra_median_range <- paste0(round(median(jco_actual_dose, na.rm = T), 1), " (", paste0(range(jco_actual_dose), collapse = "-"), ")")


## Anthracycline dose (%)
anthra_0 <- paste0(sum(df$anthra_jco_dose_any == 0, na.rm = T), 
                   " (", round(sum(df$anthra_jco_dose_any == 0, na.rm = T)/nrow(df)*100, 1), "%)")

anthra_1_250 <- paste0(sum(df$anthra_jco_dose_any >= 1 & df$anthra_jco_dose_any <= 250, na.rm = T),
                       " (", round(sum(df$anthra_jco_dose_any >= 1 & df$anthra_jco_dose_any <= 250, na.rm = T)/nrow(df)*100, 1), "%)")

anthra_gt_250 <- paste0(sum(df$anthra_jco_dose_any > 250, na.rm = T), 
                        " (", round(sum(df$anthra_jco_dose_any > 250, na.rm = T)/nrow(df)*100, 1), "%)")

## Average heart radiation dose median (IQR)
hrtavg_actual_dose <- df$hrtavg[df$hrtavg > 0]
hrtavg_median_range <- paste0(round(median(hrtavg_actual_dose, na.rm = T), 1), " (", round(quantile(hrtavg_actual_dose, prob=c(.25,.5,.75), type=1, na.rm = T)[1], 1), "-" , round(quantile(hrtavg_actual_dose, prob=c(.25,.5,.75), type=1, na.rm = T)[3],1), ")")

# ## Average heart radiation dose median (range)
# hrtavg_actual_dose <- df$hrtavg[df$hrtavg > 0]
# hrtavg_median_range <- paste0(round(median(hrtavg_actual_dose, na.rm = T), 1), " (", paste0(range(hrtavg_actual_dose), collapse = "-"), ")")



## Average heart radiation dose
heartrt_0 <- paste0(sum(df$hrtavg == 0, na.rm = T), 
                    " (", round(sum(df$hrtavg == 0, na.rm = T)/nrow(df)*100, 1), "%)")

heartrt_1_25 <- paste0(sum(df$hrtavg >= 1 & df$hrtavg <= 25, na.rm = T), 
                       " (", round(sum(df$hrtavg >= 1 & df$hrtavg <= 25, na.rm = T)/nrow(df)*100, 1), "%)")

heartrt_gt_25 <- paste0(sum(df$hrtavg > 25, na.rm = T), 
                        " (", round(sum(df$hrtavg > 25, na.rm = T)/nrow(df)*100, 1), "%)")

out <- setNames(rbind.data.frame(agedx, agelstcontact, df_elapsedAGE, ejection_fraction_hrt, male, female, anthra_median_range, anthra_0, anthra_1_250, anthra_gt_250, hrtavg_median_range, heartrt_0,
                                 heartrt_1_25, heartrt_gt_25), paste0("With CMP (n = ", n, ")"))

rownames(out) <- c("agedx", "agelstcontact", "df_elapsedAGE", "ejection_fraction_hrt", "male", "female", "anthra_median_range", "anthra_0", "anthra_1_250", "anthra_gt_250",
                   "hrtavg_median_range", "heartrt_0", "heartrt_1_25", "heartrt_gt_25")

out
}


###################################
## ccss_org, CCSS_exp and SJLIFE ##
###################################
setwd("Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3/pheno")
merged.dat <- read.table("sjlife_ccss_org_ccss_exp_ttn_bag3.pheno", header = T)
dim(merged.dat)

merged.dat_with_CMP <- merged.dat[merged.dat$CMP == 2,]
merged.dat_without_CMP <- merged.dat[merged.dat$CMP == 1,]

## With CMP
n=nrow(merged.dat_with_CMP)
get_demographic(merged.dat_with_CMP, n)

## Without CMP
n=nrow(merged.dat_without_CMP)
get_demographic(merged.dat_without_CMP, n)


# ALL
n = nrow(merged.dat)
get_demographic(merged.dat, n)

