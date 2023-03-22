setwd("Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3")

###############################
## Association of haplotypes ##
###############################
## After extracting haplotypes with extract_haplotypes.py, now check the association of the haplotypes
pheno <- read.table("pheno/sjlife_ccss_org_ccss_exp_ttn_bag3.pheno", header = T)
haplotypes <- read.table("haplotypes_ttn.txt", header = F)
haplotypes$haplo <- gsub('\\(|\\)|\\[|\\]',"", apply(haplotypes[2:11], 1, function(x) paste(x, collapse = "")))
# haplotypes$haplo <- apply(haplotypes[2:11], 1, function(x) paste(x, collapse = ""))
length(table(haplotypes$haplo))
# 34 # 0010000101 and 0000001101 haplotypes were not found

pheno$haplotypes <- haplotypes$haplo[match(pheno$IID, haplotypes$V1)]

haplotypes <- unique(pheno$haplotypes)
for (i in 1:length(haplotypes)){
  pheno[paste0("haplo_", haplotypes[i])] <- ifelse(pheno$haplotypes %in% haplotypes[i], 1, 0)
}

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

## Check association of haplotypes
for (i in 1:length(haplos)){
  print(paste0("Doing: ", haplos[i]))
  formulas <- paste0("CMP ~ ", haplos[i], "+ agedx + agelstcontact + gender + anthra_jco_dose_any + hrtavg + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10")
  haplo.test.adj <- glm(formula = formulas, family = binomial,
                        data = pheno)
  print(summary(haplo.test.adj))
  Sys.sleep(3)
}

# ## Only these two haplotypes were found significant:
# haplo_0111111101, haplo_1000000010
haplo.test.adj <- glm(formula = CMP ~  haplo_0111111101 + agedx + agelstcontact + gender + anthra_jco_dose_any + hrtavg + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, 
                      family = binomial,
                      data = pheno)
summary(haplo.test.adj)


haplo.test.adj <- glm(formula = CMP ~  haplo_1000000010 + agedx + agelstcontact + gender + anthra_jco_dose_any + hrtavg + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, 
                      family = binomial,
                      data = pheno)
summary(haplo.test.adj)

##########################################################################
## Association of significant haplotypes with respect to echo in SJLIFE ##
##########################################################################
haplo <- pheno[c("FID", "haplo_0111111101", "haplo_1000000010")]

## Extract echo measures
# rm(list=ls())
# setwd('/Volumes/Groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3//')
setwd("Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3/")
## Phenotype data for GWAS
pheno_gwas = read.table('pheno/sjlife_ttn_bag3.pheno', header = TRUE)

## Echo measures with data from most recent visit
echo_data = read.table('../echo_data/SJLIFE_echo_measures_most_recent_visit.txt', header = TRUE)

## Merge both data
pheno_final = merge(pheno_gwas, echo_data, by.x = 'FID', by.y = 'sjlid')

## Merge phenotype and significant haplotypes data
dat_final = merge(pheno_final, haplo, by = 'FID')

## Analyze both variants wrt each echo measures
# 0111111101 haplotype
haplo.0111111101 = NULL
echo_phenotypes = colnames(echo_data)[-c(1:4)]
for (p in 1:length(echo_phenotypes)){
  dat_eff = subset(dat_final, !is.na(dat_final[echo_phenotypes[p]]))
  # Remove values with >3 sd away from the absolute value
  dat_eff = subset(dat_eff, abs(scale(dat_eff[,echo_phenotypes[p]])) <= 6)
  # Adjust for BSA; two (LVMass2D_Index and LA_Volume_Index) are already done
  if (!(echo_phenotypes[p] == "LVMass2D_Index" | echo_phenotypes[p] == "LA_Volume_Index"| 
        echo_phenotypes[p] == "LV_Ejection_Fraction_3D" | echo_phenotypes[p] == "LV_Cardiac_Output_3D" | echo_phenotypes[p] == "LV_GLPS_AVG")){
    fit = lm(scale(dat_eff[,echo_phenotypes[p]]/BodySurfaceArea)~haplo_0111111101+agedx+agelstcontact+gender+anthra_jco_dose_any+hrtavg+
               PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10,data=dat_eff)
  } else {
    fit = lm(scale(dat_eff[,echo_phenotypes[p]])~haplo_0111111101+agedx+agelstcontact+gender+anthra_jco_dose_any+hrtavg+
               PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10,data=dat_eff)
  }
  beta = summary(fit)$coef[2,1]
  se = summary(fit)$coef[2,2]
  pval = summary(fit)$coef[2,4]
  haplo.0111111101 = rbind(haplo.0111111101, data.frame(phenotype=echo_phenotypes[p], beta, se, pval, n=nrow(dat_eff)))
}
haplo.0111111101$p_BH = p.adjust(haplo.0111111101$pval, method = "BH")

# 1000000010 haplotype
haplo.1000000010 = NULL
for (p in 1:length(echo_phenotypes)){
  dat_eff = subset(dat_final, !is.na(dat_final[echo_phenotypes[p]]))
  # Remove values with >3 sd away from the absolute value
  dat_eff = subset(dat_eff, abs(scale(dat_eff[,echo_phenotypes[p]])) <= 6)
  # Adjust for BSA; two (LVMass2D_Index and LA_Volume_Index) are already done
  if (!(echo_phenotypes[p] == "LVMass2D_Index" | echo_phenotypes[p] == "LA_Volume_Index"| 
        echo_phenotypes[p] == "LV_Ejection_Fraction_3D" | echo_phenotypes[p] == "LV_Cardiac_Output_3D" | echo_phenotypes[p] == "LV_GLPS_AVG")){
    fit = lm(scale(dat_eff[,echo_phenotypes[p]]/BodySurfaceArea)~haplo_1000000010+agedx+agelstcontact+gender+anthra_jco_dose_any+hrtavg+
               PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10,data=dat_eff)
  } else {
    fit = lm(scale(dat_eff[,echo_phenotypes[p]])~haplo_1000000010+agedx+agelstcontact+gender+anthra_jco_dose_any+hrtavg+
               PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10,data=dat_eff)
  }
  beta = summary(fit)$coef[2,1]
  se = summary(fit)$coef[2,2]
  pval = summary(fit)$coef[2,4]
  haplo.1000000010 = rbind(haplo.1000000010, data.frame(phenotype=echo_phenotypes[p], beta, se, pval, n=nrow(dat_eff)))
}
haplo.1000000010$p_BH = p.adjust(haplo.1000000010$pval, method = "BH")

print(haplo.0111111101); print(haplo.1000000010)
