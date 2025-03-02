## Extract echo measures
rm(list=ls())
# setwd('/Volumes/Groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3//')
setwd("Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3/")

# Function to calculate odds ratio and confidence interval
calculate_odds_ratio <- function(coef, se) {
  odds_ratio <- exp(coef)
  ci_lower <- exp(coef - 1.96 * se)
  ci_upper <- exp(coef + 1.96 * se)
  return(list(odds_ratio = odds_ratio, ci_lower = ci_lower, ci_upper = ci_upper))
}

## Phenotype data for GWAS
pheno_gwas = read.table('pheno/sjlife_ccss_org_ccss_exp_ttn_bag3.pheno', header = TRUE)

#########
## ALL ##
#########
pheno_final <- pheno_gwas
## Genotypes of TTN and BAG3 top SNPs
bag3 = read.table('sjlife_ccss_org_ccss_exp_chr10_119670121_T_C.raw', header = TRUE)
bag3[c(2:6)] = NULL
ttn = read.table('sjlife_ccss_org_ccss_exp_chr2_178562809_T_C.raw', header = TRUE)
ttn[c(2:6)] = NULL
geno = merge(bag3, ttn, by="FID")

## Merge phenotype and genotype data
dat_final = merge(pheno_final, geno, by = 'FID')
table(dat_final$CMP)

# dat_final$chr10.119670121.T.C_C <- factor(ifelse(dat_final$chr10.119670121.T.C_C != 0, "1", "0"))
# dat_final$chr2.178562809.T.C_C <- factor(ifelse(dat_final$chr2.178562809.T.C_C != 0, "1", "0"))

## Get carrier status of TTN and BAG3
# dat_final$carrier <- factor(ifelse(rowSums(dat_final[c("chr2.178562809.T.C_C", "chr10.119670121.T.C_C")]) != 0, "1", "0"))
dat_final$carrier <- factor(ifelse(rowSums(dat_final[c("chr2.178562809.T.C_C", "chr10.119670121.T.C_C")]) != 0, "1", "0"))
dat_final$CMP <- factor(ifelse(dat_final$CMP == 2, 1, 0))

# fit <- glm(formula = CMP ~ chr10.119670121.T.C_C + chr2.178562809.T.C_C + agedx + agelstcontact + gender + anthra_jco_dose_any + hrtavg + PC1 + PC2 + PC3 +
#       PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, family = binomial,
#     data = dat_final)


#------------------1. Gender specific analysis



# Overall Analysis
all_data <- dat_final
overall_model <- glm(CMP ~ carrier + gender + agedx + agelstcontact + anthra_jco_dose_any + hrtavg + PC1 + PC2 + PC3 +
                       PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, family = binomial, data = all_data)
overall_summary <- summary(overall_model)

# Male-specific Analysis
male_pheno <- all_data[all_data$gender == 0, ]
male_model <- glm(CMP ~ carrier + agedx + agelstcontact + anthra_jco_dose_any + hrtavg + PC1 + PC2 + PC3 +
                    PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, family = binomial, data = male_pheno)
male_summary <- summary(male_model)

# Female-specific Analysis
female_pheno <- all_data[all_data$gender == 1, ]
female_model <- glm(CMP ~ carrier + agedx + agelstcontact + anthra_jco_dose_any + hrtavg + PC1 + PC2 + PC3 +
                      PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, family = binomial, data = female_pheno)
female_summary <- summary(female_model)


# Create the table
table_data <- data.frame(
  Gender = c("Both", "Male", "Female"),
  Estimate = c(overall_summary$coefficients[2, 1], male_summary$coefficients[2, 1], female_summary$coefficients[2, 1]),
  Std_Error = c(overall_summary$coefficients[2, 2], male_summary$coefficients[2, 2], female_summary$coefficients[2, 2]),
  Z_Value = c(overall_summary$coefficients[2, 3], male_summary$coefficients[2, 3], female_summary$coefficients[2, 3]),
  P_Value = c(overall_summary$coefficients[2, 4], male_summary$coefficients[2, 4], female_summary$coefficients[2, 4]),
  n = c(paste0(nrow(all_data), " (", sum(all_data$CMP==1), ")"), paste0(nrow(male_pheno), " (", sum(male_pheno$CMP==1), ")"), paste0(nrow(female_pheno), " (", sum(female_pheno$CMP==1), ")")),
  stringsAsFactors = FALSE
)

# Calculate odds ratio and confidence interval
odds_ratio_info <- calculate_odds_ratio(table_data$Estimate, table_data$Std_Error)
table_data$Odds_Ratio <- odds_ratio_info$odds_ratio
table_data$CI_Lower <- odds_ratio_info$ci_lower
table_data$CI_Upper <- odds_ratio_info$ci_upper

# Print the table
print(table_data)
# Gender    Estimate Std_Error    Z_Value    P_Value          n Odds_Ratio  CI_Lower  CI_Upper
# 1   Both -0.22010842 0.1021540 -2.1546733 0.03118741 6249 (461)  0.8024318 0.6568299 0.9803098
# 2   Male -0.36574462 0.1462691 -2.5004919 0.01240210 3109 (222)  0.6936799 0.5207777 0.9239870
# 3 Female -0.08349255 0.1439014 -0.5802067 0.56177525 3140 (239)  0.9198979 0.6938224 1.2196381



# Overall Analysis (carrier and gender in the same model):
#
# carrier: The coefficient estimate of carrier1 (genetic carrier status) is
# -0.4032, with a standard error of 0.1562. The z-value is -2.581, and the
# p-value is 0.00985. The p-value is less than the significance level of 0.05,
# indicating that the genetic carrier status (carrier) is statistically
# significant in predicting the odds of CMP. The negative coefficient suggests
# that being a carrier is associated with a decrease in the log odds of CMP.
# gender: The coefficient estimate of gender is -0.2944, with a standard error
# of 0.1574. The z-value is -1.870, and the p-value is 0.06148. The p-value is
# greater than the significance level of 0.05, indicating that gender is not
# statistically significant in predicting the odds of CMP. The negative
# coefficient suggests that being male is associated with a decrease in the log
# odds of CMP. 

# Analysis for Males (carrier in the male-specific model):
#
# carrier: The coefficient estimate of carrier1 (genetic carrier status) is
# -0.4469, with a standard error of 0.2059. The z-value is -2.171, and the
# p-value is 0.02994. The p-value is less than 0.05, indicating that genetic
# carrier status (carrier) is statistically significant in predicting the odds
# of CMP for males. The negative coefficient suggests that being a carrier is
# associated with a decrease in the log odds of CMP among males. 

# Analysis for Females (carrier in the female-specific model):
#
# carrier: The coefficient estimate of carrier1 (genetic carrier status) is
# -0.3245, with a standard error of 0.2438. The z-value is -1.331, and the
# p-value is 0.1831. The p-value is greater than 0.05, indicating that genetic
# carrier status (carrier) is not statistically significant in predicting the
# odds of CMP for females. The negative coefficient suggests that being a
# carrier is associated with a decrease in the log odds of CMP among females,
# but the effect is not statistically significant at the chosen significance
# level.
#
# In summary, the analysis suggests that genetic carrier status (carrier) is
# statistically significant in predicting the odds of CMP in the overall model
# and among males. However, the effect of genetic carrier status on CMP among
# females is not statistically significant. Gender (gender) itself does not
# appear to be a significant predictor of CMP in the overall model.



#------------------2. Anthracycline only, no heart radiation
all_data <- dat_final[dat_final$anthra_jco_dose_any > 0 & dat_final$hrtavg == 0,]
all_data$anthra_dose <- factor(ifelse(all_data$anthra_jco_dose_any < 250, "low", "high"), levels = c("low", "high"))


# Overall Analysis
overall_model <- glm(CMP ~ anthra_dose + carrier + gender + agedx + agelstcontact + hrtavg + PC1 + PC2 + PC3 +
                       PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, family = binomial, data = all_data)
overall_summary <- summary(overall_model)
overall_summary

# Male-specific Analysis
male_pheno <- all_data[all_data$gender == 0, ]
male_model <- glm(CMP ~ anthra_dose + carrier + agedx + agelstcontact + hrtavg + PC1 + PC2 + PC3 +
                    PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, family = binomial, data = male_pheno)
male_summary <- summary(male_model)

# Female-specific Analysis
female_pheno <- all_data[all_data$gender == 1, ]
female_model <- glm(CMP ~ anthra_dose + carrier + agedx + agelstcontact + hrtavg + PC1 + PC2 + PC3 +
                      PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, family = binomial, data = female_pheno)
female_summary <- summary(female_model)


# Create the table
table_data <- data.frame(
  Gender = c("Both", "Male", "Female"),
  Estimate = c(overall_summary$coefficients[2, 1], male_summary$coefficients[2, 1], female_summary$coefficients[2, 1]),
  Std_Error = c(overall_summary$coefficients[2, 2], male_summary$coefficients[2, 2], female_summary$coefficients[2, 2]),
  Z_Value = c(overall_summary$coefficients[2, 3], male_summary$coefficients[2, 3], female_summary$coefficients[2, 3]),
  P_Value = c(overall_summary$coefficients[2, 4], male_summary$coefficients[2, 4], female_summary$coefficients[2, 4]),
  n = c(paste0(nrow(all_data), " (", sum(all_data$CMP==1), ")"), paste0(nrow(male_pheno), " (", sum(male_pheno$CMP==1), ")"), paste0(nrow(female_pheno), " (", sum(female_pheno$CMP==1), ")")),
  stringsAsFactors = FALSE
)

# Calculate odds ratio and confidence interval
odds_ratio_info <- calculate_odds_ratio(table_data$Estimate, table_data$Std_Error)
table_data$Odds_Ratio <- odds_ratio_info$odds_ratio
table_data$CI_Lower <- odds_ratio_info$ci_lower
table_data$CI_Upper <- odds_ratio_info$ci_upper

# Print the table
print(table_data)
# Gender Estimate Std_Error  Z_Value      P_Value          n Odds_Ratio CI_Lower CI_Upper
# 1   Both 1.325270 0.2187225 6.059141 1.368503e-09 1804 (131)   3.763203 2.451192 5.777473
# 2   Male 1.393756 0.3044173 4.578438 4.684603e-06   906 (69)   4.029958 2.219092 7.318563
# 3 Female 1.346099 0.3238450 4.156615 3.229971e-05   898 (62)   3.842407 2.036765 7.248795


#------------------3. Heart radiation; no anthracycline
all_data <- dat_final[dat_final$anthra_jco_dose_any == 0 & dat_final$hrtavg > 0,]
all_data$heart_RT <- factor(ifelse(all_data$hrtavg < 15, "low", "high"), levels = c("low", "high"))

# Overall Analysis
overall_model <- glm(CMP ~ heart_RT + carrier + gender + agedx + agelstcontact + anthra_jco_dose_any + PC1 + PC2 + PC3 +
                       PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, family = binomial, data = all_data)
overall_summary <- summary(overall_model)
overall_summary

# Male-specific Analysis
male_pheno <- all_data[all_data$gender == 0, ]
male_model <- glm(CMP ~ heart_RT + carrier + agedx + agelstcontact + anthra_jco_dose_any + PC1 + PC2 + PC3 +
                    PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, family = binomial, data = male_pheno)
male_summary <- summary(male_model)

# Female-specific Analysis
female_pheno <- all_data[all_data$gender == 1, ]
female_model <- glm(CMP ~ heart_RT + carrier + agedx + agelstcontact + anthra_jco_dose_any + PC1 + PC2 + PC3 +
                      PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, family = binomial, data = female_pheno)
female_summary <- summary(female_model)


# Create the table
table_data <- data.frame(
  Gender = c("Both", "Male", "Female"),
  Estimate = c(overall_summary$coefficients[2, 1], male_summary$coefficients[2, 1], female_summary$coefficients[2, 1]),
  Std_Error = c(overall_summary$coefficients[2, 2], male_summary$coefficients[2, 2], female_summary$coefficients[2, 2]),
  Z_Value = c(overall_summary$coefficients[2, 3], male_summary$coefficients[2, 3], female_summary$coefficients[2, 3]),
  P_Value = c(overall_summary$coefficients[2, 4], male_summary$coefficients[2, 4], female_summary$coefficients[2, 4]),
  n = c(paste0(nrow(all_data), " (", sum(all_data$CMP==1), ")"), paste0(nrow(male_pheno), " (", sum(male_pheno$CMP==1), ")"), paste0(nrow(female_pheno), " (", sum(female_pheno$CMP==1), ")")),
  stringsAsFactors = FALSE
)

# Calculate odds ratio and confidence interval
odds_ratio_info <- calculate_odds_ratio(table_data$Estimate, table_data$Std_Error)
table_data$Odds_Ratio <- odds_ratio_info$odds_ratio
table_data$CI_Lower <- odds_ratio_info$ci_lower
table_data$CI_Upper <- odds_ratio_info$ci_upper

# Print the table
print(table_data)
# Gender Estimate Std_Error  Z_Value      P_Value          n Odds_Ratio CI_Lower CI_Upper
# 1   Both 2.207731 0.2924933 7.547973 4.420848e-14 2477 (118)   9.095061 5.126610 16.13544
# 2   Male 2.256660 0.4259181 5.298342 1.168589e-07  1205 (53)   9.551133 4.144824 22.00917
# 3 Female 2.177241 0.4064555 5.356654 8.477736e-08  1272 (65)   8.821935 3.977241 19.56797
