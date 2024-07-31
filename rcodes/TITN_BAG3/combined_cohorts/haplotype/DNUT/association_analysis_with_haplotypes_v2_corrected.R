setwd("Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3/") 

###############################
## Association of haplotypes ##
###############################
## After extracting haplotypes with extract_haplotypes.py, now check the association of the haplotypes
pheno <- read.table("pheno/sjlife_ccss_org_ccss_exp_ttn_bag3.pheno", header = T)
haplotypes <- read.table("haplotypes_ttn_r2_0.8_haplo.glm.txt", header = T, sep = "\t", stringsAsFactors = FALSE, colClasses = "character") # r2 > 0.8
# # haplotypes$haplo <- gsub('\\(|\\)|\\[|\\]',"", apply(haplotypes[2:11], 1, function(x) paste(x, collapse = "")))
# haplotypes$haplo <- gsub('\\(|\\)|\\[|\\]',"", apply(haplotypes[2:10], 1, function(x) paste(x, collapse = "")))
# # haplotypes$haplo <- apply(haplotypes[2:11], 1, function(x) paste(x, collapse = ""))
# length(table(haplotypes$haplo))
# # 28 


freq <- read.table(text="index  haplotype     freq         SE
         1   000000000    0.654829    0.000274
         2   000000001    0.009935    0.000147
         3   000000010    0.110620    0.000227
         4   000000011    0.000003    0.000017
         5   000000100    0.000107    0.000044
         6   000000101    0.000696    0.000174
         7   000000111    0.000173    0.000123
         8   000001100    0.000006    0.000030
         9   000001101    0.000075    0.000019
        10   000001110    0.000154    0.000030
        11   000010101    0.000148    0.000083
        12   000011101    0.000447    0.000054
        13   000100000    0.000022    0.000036
        14   000100010    0.000002    0.000014
        15   000111101    0.000689    0.000039
        16   010000000    0.000159    0.000042
        17   010000001    0.000012    0.000029
        18   010000100    0.001993    0.000123
        19   010000101    0.000112    0.000109
        20   010010100    0.000972    0.000075
        21   010010101    0.044120    0.000098
        22   010011101    0.000068    0.000029
        23   010100100    0.000075    0.000055
        24   010110100    0.000078    0.000016
        25   010110101    0.006944    0.000081
        26   010110111    0.000011    0.000066
        27   011110100    0.002223    0.000085
        28   011110101    0.018898    0.000088
        29   011110111    0.000002    0.000019
        30   011111100    0.000003    0.000016
        31   011111101    0.000157    0.000016
        32   100000000    0.000119    0.000040
        33   100000001    0.000018    0.000033
        34   110000000    0.000004    0.000017
        35   110000100    0.000013    0.000029
        36   110000101    0.000006    0.000022
        37   111000100    0.000038    0.000040
        38   111010101    0.000042    0.000040
        39   111011101    0.000081    0.000008
        40   111100000    0.000238    0.000011
        41   111100001    0.000002    0.000011
        42   111101101    0.000005    0.000019
        43   111111000    0.000002    0.000014
        44   111111001    0.000157    0.000016
        45   111111100    0.002293    0.000129
        46   111111101    0.142923    0.000117
        47   111111110    0.000284    0.000092
        48   111111111    0.000038    0.000076", header = T, colClasses = c("numeric", "character", "numeric", "numeric"))

dim(freq)


pheno$Best1 <- haplotypes$Best1[match(pheno$IID, haplotypes$sample)]
pheno$Best2 <- haplotypes$Best2[match(pheno$IID, haplotypes$sample)]

pheno.best1 <- cbind.data.frame(IID = pheno$IID, Best1 = pheno$Best1)
pheno.best2 <- cbind.data.frame(IID = pheno$IID, Best2 = pheno$Best2)

## Best 1
haplotypes <- unique(pheno$Best1)
for (i in 1:length(haplotypes)){
  pheno.best1[haplotypes[i]] <- ifelse(pheno.best1$Best1 %in% haplotypes[i], 1, 0)
}

## Best2
haplotypes <- unique(pheno$Best2)
for (i in 1:length(haplotypes)){
  pheno.best2[haplotypes[i]] <- ifelse(pheno.best2$Best2 %in% haplotypes[i], 1, 0)
}

all_columns <- unique(c(colnames(pheno.best1), colnames(pheno.best2)))
result_df <- data.frame(matrix(NA, nrow = nrow(pheno.best1), ncol = length(all_columns)))
colnames(result_df) <- all_columns
result_df$IID <- pheno.best1$IID
result_df$Best1 <- pheno.best1$Best1
result_df$Best2 <- pheno.best2$Best2

# Fill the new dataframe by summing or copying the values
for (col in all_columns) {
  if (col %in% c("IID", "Best1", "Best2")) next
  
  if (col %in% colnames(pheno.best1) && col %in% colnames(pheno.best2)) {
    result_df[[col]] <- pheno.best1[[col]] + pheno.best2[[col]]
  } else if (col %in% colnames(pheno.best1)) {
    result_df[[col]] <- pheno.best1[[col]]
  } else if (col %in% colnames(pheno.best2)) {
    result_df[[col]] <- pheno.best2[[col]]
  }
}

result_df <- result_df[!grepl("Best", colnames(result_df))]
colnames(result_df)[-1] <- paste0("haplo_", colnames(result_df))[-1]
merged_pheno <- merge(pheno, result_df, by = "IID")
save.pheno <- pheno
pheno <- merged_pheno

pheno$CMP <- as.factor(pheno$CMP)

haplos <- colnames(pheno)[grepl("haplo_", colnames(pheno))]

wanted.vars <- {}
## Check association of haplotypes
for (haplo in haplos) {
  cat("Doing:", haplo, "\n")
  formula <- as.formula(paste("CMP ~", haplo, "+ agedx + agelstcontact + gender + anthra_jco_dose_any + hrtavg + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10"))
  haplo.test.adj <- glm(formula, family = binomial, data = pheno)
  
  summary_table <- summary(haplo.test.adj)$coefficients
  row <- grep(haplo, rownames(summary_table))
  P <- summary_table[row, "Pr(>|z|)"]
  OR <- exp(summary_table[row, "Estimate"])
  CI1 <- exp(summary_table[row, "Estimate"] - 1.96 * summary_table[row, "Std. Error"])
  CI2 <- exp(summary_table[row, "Estimate"] + 1.96 * summary_table[row, "Std. Error"])
  
  # Round OR and CI values and handle values greater than 99
  OR_str <- sprintf("%.2f", OR)
  CI1_str <- ifelse(CI1 > 99, ">99", sprintf("%.2f", CI1))
  CI2_str <- ifelse(CI2 > 99, ">99", sprintf("%.2f", CI2))
  
  OR_CI <- paste0(OR_str, " (", CI1_str, " to ", CI2_str, ")")
  
  # Format p-value
  P_str <- ifelse(P <= 0.001, format(P, scientific = TRUE, digits = 2),
                  sprintf("%.2f", P))
  
  wanted.haplo <- gsub("haplo_", "", haplo)
  wanted.freq <- freq$freq[match(wanted.haplo, freq$haplotype)]
  wanted.freq.SE <- freq$SE[match(wanted.haplo, freq$haplotype)]
  wanted.freq_str <- paste0(wanted.freq, " (", wanted.freq.SE, ")")
  
  wanted.vars[[haplo]] <- c(haplotype = wanted.haplo, frequency = wanted.freq_str, OR_CI = OR_CI, P_value = P_str)
}

wanted.vars <- do.call(rbind, wanted.vars)
wanted.vars <- as.data.frame(wanted.vars, stringsAsFactors = FALSE)
colnames(wanted.vars) <- c("Haplotype", "Frequency (SE)", "OR (95% CI)", "P-value")

print(wanted.vars)

  
# haplo_000000010      2.898e-01  1.057e-01   2.743  0.00609 **
# haplo_111111101     -2.192e-01  1.066e-01  -2.055   0.0398 *  


haplo.test.adj <- glm(formula = CMP ~  haplo_000000010 + agedx + agelstcontact + gender + anthra_jco_dose_any + hrtavg + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, 
                      family = binomial,
                      data = pheno)
summary(haplo.test.adj)


haplo.test.adj <- glm(formula = CMP ~  haplo_111111101 + agedx + agelstcontact + gender + anthra_jco_dose_any + hrtavg + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, 
                      family = binomial,
                      data = pheno)
summary(haplo.test.adj)
##########################################################################
## Association of significant haplotypes with respect to echo in SJLIFE ##
##########################################################################
# haplo <- pheno[c("FID", "haplo_0111111101", "haplo_1000000010")]
haplo <- pheno[c("FID", "haplo_11110", "haplo_00001")]

## Extract echo measures
# rm(list=ls())
# setwd('/Volumes/Groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3//')
setwd("Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3/haplotype_analysis_v2/")
## Phenotype data for GWAS
pheno_gwas = read.table('../pheno/sjlife_ttn_bag3.pheno', header = TRUE)

## Echo measures with data from most recent visit
echo_data = read.table('../../echo_data/SJLIFE_echo_measures_most_recent_visit.txt', header = TRUE)

## Merge both data
pheno_final = merge(pheno_gwas, echo_data, by.x = 'FID', by.y = 'sjlid')

## Merge phenotype and significant haplotypes data
dat_final = merge(pheno_final, haplo, by = 'FID')

## Analyze both variants wrt each echo measures
# 0111101 haplotype
haplo_11110 = NULL
echo_phenotypes = colnames(echo_data)[-c(1:4)]
for (p in 1:length(echo_phenotypes)){
  dat_eff = subset(dat_final, !is.na(dat_final[echo_phenotypes[p]]))
  # Remove values with >3 sd away from the absolute value
  dat_eff = subset(dat_eff, abs(scale(dat_eff[,echo_phenotypes[p]])) <= 6)
  # Adjust for BSA; two (LVMass2D_Index and LA_Volume_Index) are already done
  if (!(echo_phenotypes[p] == "LVMass2D_Index" | echo_phenotypes[p] == "LA_Volume_Index"| 
        echo_phenotypes[p] == "LV_Ejection_Fraction_3D" | echo_phenotypes[p] == "LV_Cardiac_Output_3D" | echo_phenotypes[p] == "LV_GLPS_AVG")){
    fit = lm(scale(dat_eff[,echo_phenotypes[p]]/BodySurfaceArea)~haplo_11110+agedx+agelstcontact+gender+anthra_jco_dose_any+hrtavg+
               PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10,data=dat_eff)
  } else {
    fit = lm(scale(dat_eff[,echo_phenotypes[p]])~haplo_11110+agedx+agelstcontact+gender+anthra_jco_dose_any+hrtavg+
               PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10,data=dat_eff)
  }
  beta = summary(fit)$coef[2,1]
  se = summary(fit)$coef[2,2]
  pval = summary(fit)$coef[2,4]
  haplo_11110 = rbind(haplo_11110, data.frame(phenotype=echo_phenotypes[p], beta, se, pval, n=nrow(dat_eff)))
}
haplo_11110$p_BH = p.adjust(haplo_11110$pval, method = "BH")

# 1000010 haplotype
haplo_00001 = NULL
for (p in 1:length(echo_phenotypes)){
  dat_eff = subset(dat_final, !is.na(dat_final[echo_phenotypes[p]]))
  # Remove values with >3 sd away from the absolute value
  dat_eff = subset(dat_eff, abs(scale(dat_eff[,echo_phenotypes[p]])) <= 6)
  # Adjust for BSA; two (LVMass2D_Index and LA_Volume_Index) are already done
  if (!(echo_phenotypes[p] == "LVMass2D_Index" | echo_phenotypes[p] == "LA_Volume_Index"| 
        echo_phenotypes[p] == "LV_Ejection_Fraction_3D" | echo_phenotypes[p] == "LV_Cardiac_Output_3D" | echo_phenotypes[p] == "LV_GLPS_AVG")){
    fit = lm(scale(dat_eff[,echo_phenotypes[p]]/BodySurfaceArea)~haplo_00001+agedx+agelstcontact+gender+anthra_jco_dose_any+hrtavg+
               PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10,data=dat_eff)
  } else {
    fit = lm(scale(dat_eff[,echo_phenotypes[p]])~haplo_00001+agedx+agelstcontact+gender+anthra_jco_dose_any+hrtavg+
               PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10,data=dat_eff)
  }
  beta = summary(fit)$coef[2,1]
  se = summary(fit)$coef[2,2]
  pval = summary(fit)$coef[2,4]
  haplo_00001 = rbind(haplo_00001, data.frame(phenotype=echo_phenotypes[p], beta, se, pval, n=nrow(dat_eff)))
}
haplo_00001$p_BH = p.adjust(haplo_00001$pval, method = "BH")

print(haplo_11110); print(haplo_00001)
