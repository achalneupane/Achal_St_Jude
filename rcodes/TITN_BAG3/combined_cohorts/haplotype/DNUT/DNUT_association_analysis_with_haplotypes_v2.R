# setwd("Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3") # (old)
setwd("Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3/haplotype_analysis_v2/") # (version_2)

###############################
## Association of haplotypes ##
###############################
## After extracting haplotypes with extract_haplotypes.py, now check the association of the haplotypes
pheno <- read.table("../pheno/sjlife_ccss_org_ccss_exp_ttn_bag3.pheno", header = T)
# haplotypes <- read.table("haplotypes_ttn.txt", header = F) # r2 < 0.8
haplotypes <- read.table("haplotypes_ttn_r2_0.2.txt", header = F) # r2 < 0.2
# haplotypes$haplo <- gsub('\\(|\\)|\\[|\\]',"", apply(haplotypes[2:11], 1, function(x) paste(x, collapse = "")))
haplotypes$haplo <- gsub('\\(|\\)|\\[|\\]',"", apply(haplotypes[2:7], 1, function(x) paste(x, collapse = "")))
# haplotypes$haplo <- apply(haplotypes[2:11], 1, function(x) paste(x, collapse = ""))
length(table(haplotypes$haplo))
# 34 # 0010000101 and 0000001101 haplotypes were not found

pheno$haplotypes <- haplotypes$haplo[match(pheno$IID, haplotypes$V1)]


freq <- read.table(text="index  haplotype     freq         SE
         1      000000    0.654865    0.000270
         2      000001    0.010008    0.000117
         3      000010    0.110590    0.000239
         4      000011    0.000001    0.000008
         5      000100    0.000168    0.000075
         6      000101    0.000485    0.000161
         7      000110    0.000069    0.000073
         8      000111    0.000279    0.000115
         9      001101    0.001490    0.000084
        10      010000    0.000158    0.000025
        11      010001    0.000005    0.000019
        12      010100    0.002179    0.000035
        13      010101    0.000001    0.000008
        14      011100    0.003316    0.000111
        15      011101    0.070140    0.000122
        16      011110    0.000001    0.000008
        17      011111    0.000011    0.000054
        18      100000    0.000090    0.000031
        19      110000    0.000250    0.000030
        20      110001    0.000001    0.000008
        21      110101    0.000071    0.000027
        22      111001    0.000159    0.000008
        23      111100    0.002257    0.000145
        24      111101    0.143082    0.000120
        25      111110    0.000309    0.000093
        26      111111    0.000014    0.000059", header = T, colClasses = c("numeric", "character", "numeric", "numeric"))

dim(freq)

freq$haplotype

# sum(unique(pheno$haplotypes) %in% freq$haplotype)

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
  # Sys.sleep(10)
}

# ## Only these two haplotypes were found significant:
# haplo_11110; haplo_00001    ## r2 < 0.2
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
