setwd("Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3/") 

###############################
## Association of haplotypes ##
###############################
## After extracting haplotypes with extract_haplotypes.py, now check the association of the haplotypes
pheno <- read.table("pheno/sjlife_ccss_org_ccss_exp_ttn_bag3.pheno", header = T)
haplotypes <- read.table("haplotypes_ttn_r2_0.8.txt", header = F) # r2 > 0.8
# haplotypes$haplo <- gsub('\\(|\\)|\\[|\\]',"", apply(haplotypes[2:11], 1, function(x) paste(x, collapse = "")))
haplotypes$haplo <- gsub('\\(|\\)|\\[|\\]',"", apply(haplotypes[2:10], 1, function(x) paste(x, collapse = "")))
# haplotypes$haplo <- apply(haplotypes[2:11], 1, function(x) paste(x, collapse = ""))
length(table(haplotypes$haplo))
# 28 


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

freq$haplotype


pheno$haplotypes <- haplotypes$haplo[match(pheno$IID, haplotypes$V1)]

sum(unique(pheno$haplotypes) %in% freq$haplotype)
## 28 

freq <- freq[freq$haplotype %in% unique(haplotypes$haplo),]


haplotypes <- unique(pheno$haplotypes)
for (i in 1:length(haplotypes)){
  pheno[paste0("haplo_", haplotypes[i])] <- ifelse(pheno$haplotypes %in% haplotypes[i], 1, 0)
}

# freq.haplos <- c("0000000", "0000001", "0000010", "0000100", "0000101", "0001101", "0010000", "0010001", "0010100", "0011100", "0011101", "0100000", "0101101", "0110000", "0110101", "0111001", "0111100", "0111101", "0111110", "1000000", "1000010", "1000111", "1011101", "1111101")
# freq.haplos[!freq.haplos %in% haplotypes]

freq


# > colnames(pheno)
# [1] "FID"                   "IID"                   "CMP"                   "agedx"                 "agelstcontact"         "gender"                "anthra_jco_dose_any"  
# [8] "hrtavg"                "agevent"               "ejection_fraction_hrt" "PC1"                   "PC2"                   "PC3"                   "PC4"                  
# [15] "PC5"                   "PC6"                   "PC7"                   "PC8"                   "PC9"                   "PC10"                  "haplotypes"           
# [22] "haplo_0000000100"      "haplo_0000000100"      "haplo_1000000010"      "haplo_0000000000"      "haplo_0000000010"      "haplo_1000000000"      "haplo_0010110101"     
# [29] "haplo_0011110101"      "haplo_0010010101"      "haplo_0000000001"      "haplo_0010000100"      "haplo_1000000111"      "haplo_1010000000"      "haplo_0000000101"     
# [36] "haplo_0010010100"      "haplo_0010100100"      "haplo_0111011101"      "haplo_0000111101"      "haplo_0011110100"      "haplo_1010010101"      "haplo_1111111101"     
# [43] "haplo_0111111110"      "haplo_0010011101"      "haplo_0111111100"      "haplo_1000001110"      "haplo_0000011101"      "haplo_0111100000"      "haplo_0111000100"     
# [50] "haplo_0010000000"      "haplo_0010110100"      "haplo_0011111101"      "haplo_0100000000"      "haplo_0111111001"      "haplo_0000010101" 

pheno$CMP <- as.factor(pheno$CMP)

haplos <- colnames(pheno)[grepl("haplo_", colnames(pheno))]

wanted.vars <- {}
## Check association of haplotypes
for (i in 1:length(haplos)){
  print(paste0("Doing: ", haplos[i]))
  formulas <- paste0("CMP ~ ", haplos[i], "+ agedx + agelstcontact + gender + anthra_jco_dose_any + hrtavg + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10")
  haplo.test.adj <- glm(formula = formulas, family = binomial,
                        data = pheno)
  print(summary(haplo.test.adj))
  # Sys.sleep(5)
  wanted.haplo <- gsub("haplo_", "", haplos[i], )
  wanted.freq <- freq$freq[match(wanted.haplo, freq$haplotype)]
  wanted.freq.SD <- freq$SE[match(wanted.haplo, freq$haplotype)]
  wanted.freq <- paste0(wanted.freq, " (", wanted.freq.SD, ")")
  summary_table <- summary(haplo.test.adj)$coefficients[, c(1, 2, 4)]
  ROW <- match(haplos[i], rownames(summary_table))
  P <- round(summary_table[ROW,3],3)
  OR <- as.numeric(sprintf("%.2f", exp(summary_table[ROW,1])))
  CI1 <- paste0(" (",round(exp(summary_table[ROW,1]) -1.96*summary_table[ROW,2],2), 
                " to ", round(exp(summary_table[ROW,1]) +1.96*summary_table[ROW,2], 2),")")
  OR <- paste0(OR, CI1)
  
  wanted.vars.tmp <- c(wanted.haplo, wanted.freq, OR, P)
  wanted.vars <- rbind(wanted.vars, wanted.vars.tmp)
  Sys.sleep(10)
}

# ## Only these two haplotypes were found significant:
# haplo_11110; haplo_00001    ## r2 > 0.8
# haplo.test.adj <- glm(formula = CMP ~  haplo_0111111101 + agedx + agelstcontact + gender + anthra_jco_dose_any + hrtavg + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, 
#                       family = binomial,
#                       data = pheno)
# summary(haplo.test.adj)
# 
# 
# haplo.test.adj <- glm(formula = CMP ~  haplo_1000000010 + agedx + agelstcontact + gender + anthra_jco_dose_any + hrtavg + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, 
#                       family = binomial,
#                       data = pheno)
# summary(haplo.test.adj)

haplo.test.adj <- glm(formula = CMP ~  haplo_111101 + agedx + agelstcontact + gender + anthra_jco_dose_any + hrtavg + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, 
                      family = binomial,
                      data = pheno)
summary(haplo.test.adj)


haplo.test.adj <- glm(formula = CMP ~  haplo_000010 + agedx + agelstcontact + gender + anthra_jco_dose_any + hrtavg + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, 
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
