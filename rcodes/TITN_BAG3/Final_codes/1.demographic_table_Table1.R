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
# merged.dat <- read.table("sjlife_ccss_org_ccss_exp_ttn_bag3.pheno", header = T)
merged.dat <- read.table("sjlife_ccss_org_ccss_exp_ttn_bag3_kendrick.pheno", header = T) # new data from Kendrick
dim(merged.dat)
# merged.dat$hrtavg <- merged.dat$hrtavg/100

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




############
## SJLIFE ##
############
setwd("Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3/pheno")
# merged.dat <- read.table("sjlife_ccss_org_ccss_exp_ttn_bag3.pheno", header = T)
merged.dat <- read.table("sjlife_ccss_org_ccss_exp_ttn_bag3_kendrick.pheno", header = T)
merged.dat <- merged.dat[merged.dat$cohort_two ==1,]
dim(merged.dat)
# merged.dat$hrtavg <- merged.dat$hrtavg/100

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

##########
## CCSS ##
##########
setwd("Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3/pheno")
# merged.dat <- read.table("sjlife_ccss_org_ccss_exp_ttn_bag3.pheno", header = T)
merged.dat <- read.table("sjlife_ccss_org_ccss_exp_ttn_bag3_kendrick.pheno", header = T)
merged.dat <- merged.dat[merged.dat$cohort_two ==2,]
dim(merged.dat)
# merged.dat$hrtavg <- merged.dat$hrtavg/100

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



#############
## African ##
#############

####################################
## Function to create demographic ##
####################################
# Function to create demographic table (Table 1) without elapsedAge
get_demographic <- function(df, n){
  ## Ejection fraction
  ejection_fraction_hrt <- paste0(round((median(df$ejection_fraction_hrt, na.rm = T)*100), 1), " (", round((quantile(df$ejection_fraction_hrt, prob=c(.25,.5,.75), type=1, na.rm = T)*100)[1], 1), "-" , round((quantile(df$ejection_fraction_hrt, prob=c(.25,.5,.75), type=1, na.rm = T)*100)[3],1), ")")
  
  ## Elapsed time
  # df$elapsedAGE <- df$agevent - df$agedx
  
  # All elapsed negative age are converted to zero. ageevent probably used floor value so there are some negatives
  # df$elapsedAGE[df$elapsedAGE < 0] <- 0
  # df_elapsedAGE <- paste0(round(median(df$elapsedAGE,  na.rm = T), 1), " (", round(quantile(df$elapsedAGE, prob=c(.25,.5,.75), type=1, na.rm = T)[1], 1), "-" , round(quantile(df$elapsedAGE, prob=c(.25,.5,.75), type=1,  na.rm = T)[3],1), ")")
  
  
  
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
  
  out <- setNames(rbind.data.frame(agedx, agelstcontact, ejection_fraction_hrt, male, female, anthra_median_range, anthra_0, anthra_1_250, anthra_gt_250, hrtavg_median_range, heartrt_0,
                                   heartrt_1_25, heartrt_gt_25), paste0("With CMP (n = ", n, ")"))
  
  rownames(out) <- c("agedx", "agelstcontact", "ejection_fraction_hrt", "male", "female", "anthra_median_range", "anthra_0", "anthra_1_250", "anthra_gt_250",
                     "hrtavg_median_range", "heartrt_0", "heartrt_1_25", "heartrt_gt_25")
  
  out
}

# > colnames(sjlife.afr.dat)
# [1] "FID"                       "IID"                       "ejection_fraction_hrt"     "gender"                   
# [5] "anthra_jco_dose_any"       "hrtavg"                    "PC1"                       "PC2"                      
# [9] "PC3"                       "PC4"                       "PC5"                       "PC6"                      
# [13] "PC7"                       "PC8"                       "PC9"                       "PC10"                     
# [17] "PC11"                      "PC12"                      "PC13"                      "PC14"                     
# [21] "PC15"                      "PC16"                      "PC17"                      "PC18"                     
# [25] "PC19"                      "PC20"                      "agedx_cat_5_10"            "agedx_cat_10_15"          
# [29] "agedx_cat_15_19.7"         "agelstcontact_cat_25_35"   "agelstcontact_cat_35_65.9" "CMP_EF_HEIRARCHY"         
# [33] "Final_Heirarchy_Grade"     "CMP"                       "agedx"                     "agelstcontact"            
# > colnames(merged.dat)
# [1] "FID"                   "IID"                   "CMP"                   "agedx"                 "agelstcontact"        
# [6] "gender"                "anthra_jco_dose_any"   "hrtavg"                "agevent"               "ejection_fraction_hrt"
# [11] "PC1"                   "PC2"                   "PC3"                   "PC4"                   "PC5"                  
# [16] "PC6"                   "PC7"                   "PC8"                   "PC9"                   "PC10"                 
# [21] "cohort"                "cohort_two"           

# cc <- merged.dat[match(pheno_gwas$IID, merged.dat$IID),]
# table(pheno_gwas$CMP == cc$CMP)
# cc <- merged.dat[(merged.dat$IID %in% pheno_gwas$IID),]
# cc.ca <-cc[cc$CMP == 2,]
# 
# gradefile<-  CTCAE[CTCAE$grade>0,]
# gradefile <- gradefile[!duplicated(gradefile$sjlid),]
# cc.ca$IID[(cc.ca$IID %in% gradefile$sjlid)] # grade 2
# cc.ca$IID[!(cc.ca$IID %in% gradefile$sjlid)] # grade 0

setwd("Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3/pheno")
# sjlife.afr.dat <- read.table("sjlife_ttn_bag3_afr.pheno", header = T)
sjlife.afr.dat <- read.table("sjlife_ttn_bag3_afr_kendrick.pheno", header = T)
dim(sjlife.afr.dat)
# sjlife.afr.dat$hrtavg <- sjlife.afr.dat$hrtavg/100

# ## Need to add ageevent to African Phenotype
# library(haven)
# CTCAE <- read_sas("Z:/SJShare/SJCOMMON/ECC/SJLife/SJLIFE Data Freeze/2 Final Data SJLIFE/20200430/Event Data/ctcaegrades.sas7bdat")
# CTCAE <- CTCAE[grepl("Cardiomyopathy", CTCAE$condition),]


sjlife.afr.dat_with_CMP <- sjlife.afr.dat[sjlife.afr.dat$CMP == 2,]
sjlife.afr.dat_without_CMP <- sjlife.afr.dat[sjlife.afr.dat$CMP == 1,]

## With CMP
n=nrow(sjlife.afr.dat_with_CMP)
get_demographic(sjlife.afr.dat_with_CMP, n)

## Without CMP
n=nrow(sjlife.afr.dat_without_CMP)
get_demographic(sjlife.afr.dat_without_CMP, n)


# ALL
n = nrow(sjlife.afr.dat)
get_demographic(sjlife.afr.dat, n)
