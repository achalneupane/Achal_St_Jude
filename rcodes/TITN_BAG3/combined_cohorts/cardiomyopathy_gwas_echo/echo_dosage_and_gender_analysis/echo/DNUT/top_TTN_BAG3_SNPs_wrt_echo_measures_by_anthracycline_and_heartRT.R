## Extract echo measures
rm(list=ls())
# setwd('/Volumes/Groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3//')
setwd("Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3/")
## Phenotype data for GWAS
pheno_gwas = read.table('pheno/sjlife_ttn_bag3.pheno', header = TRUE)

## Echo measures with data from most recent visit
echo_data = read.table('../echo_data/SJLIFE_echo_measures_most_recent_visit.txt', header = TRUE)

## Merge both data
pheno_final.all = merge(pheno_gwas, echo_data, by.x = 'FID', by.y = 'sjlid')


## 1. -----------------Anthracycline only
all_data <- pheno_final.all[pheno_final.all$anthra_jco_dose_any >0 & pheno_final.all$hrtavg == 0,]

pheno_final <- all_data
## Genotypes of TTN and BAG3 top SNPs
bag3 = read.table('sjlife_chr10_119670121_T_C.raw', header = TRUE)
bag3[c(2:6)] = NULL
ttn = read.table('sjlife_chr2_178562809_T_C.raw', header = TRUE)
ttn[c(2:6)] = NULL
geno = merge(bag3, ttn, by="FID")

## Merge phenotype and genotype data
dat_final = merge(pheno_final, geno, by = 'FID')

## Analyze both variants wrt each echo measures
# BAG3 SNP
res_bag3 = NULL
echo_phenotypes = colnames(echo_data)[-c(1:4)]
for (p in 1:length(echo_phenotypes)){
  dat_eff = subset(dat_final, !is.na(dat_final[echo_phenotypes[p]]))
  # Remove values with >3 sd away from the absolute value
  dat_eff = subset(dat_eff, abs(scale(dat_eff[,echo_phenotypes[p]])) <= 6)
  # Adjust for BSA; two (LVMass2D_Index and LA_Volume_Index) are already done
  if (!(echo_phenotypes[p] == "LVMass2D_Index" | echo_phenotypes[p] == "LA_Volume_Index"| 
        echo_phenotypes[p] == "LV_Ejection_Fraction_3D" | echo_phenotypes[p] == "LV_Cardiac_Output_3D" | echo_phenotypes[p] == "LV_GLPS_AVG")){
    fit = lm(scale(dat_eff[,echo_phenotypes[p]]/BodySurfaceArea)~chr10.119670121.T.C_C+agedx+agelstcontact+gender+anthra_jco_dose_any+hrtavg+
               PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10,data=dat_eff)
  } else {
    fit = lm(scale(dat_eff[,echo_phenotypes[p]])~chr10.119670121.T.C_C+agedx+agelstcontact+gender+anthra_jco_dose_any+hrtavg+
               PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10,data=dat_eff)
  }
  beta = summary(fit)$coef[2,1]
  se = summary(fit)$coef[2,2]
  pval = summary(fit)$coef[2,4]
  res_bag3 = rbind(res_bag3, data.frame(phenotype=echo_phenotypes[p], beta, se, pval, n=nrow(dat_eff)))
}
res_bag3$p_BH = p.adjust(res_bag3$pval, method = "BH")

# TTN SNP
res_ttn = NULL
for (p in 1:length(echo_phenotypes)){
  dat_eff = subset(dat_final, !is.na(dat_final[echo_phenotypes[p]]))
  # Remove values with >3 sd away from the absolute value
  dat_eff = subset(dat_eff, abs(scale(dat_eff[,echo_phenotypes[p]])) <= 6)
  # Adjust for BSA; two (LVMass2D_Index and LA_Volume_Index) are already done
  if (!(echo_phenotypes[p] == "LVMass2D_Index" | echo_phenotypes[p] == "LA_Volume_Index"| 
        echo_phenotypes[p] == "LV_Ejection_Fraction_3D" | echo_phenotypes[p] == "LV_Cardiac_Output_3D" | echo_phenotypes[p] == "LV_GLPS_AVG")){
    fit = lm(scale(dat_eff[,echo_phenotypes[p]]/BodySurfaceArea)~chr2.178562809.T.C_C+agedx+agelstcontact+gender+anthra_jco_dose_any+hrtavg+
               PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10,data=dat_eff)
  } else {
    fit = lm(scale(dat_eff[,echo_phenotypes[p]])~chr2.178562809.T.C_C+agedx+agelstcontact+gender+anthra_jco_dose_any+hrtavg+
               PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10,data=dat_eff)
  }
  beta = summary(fit)$coef[2,1]
  se = summary(fit)$coef[2,2]
  pval = summary(fit)$coef[2,4]
  res_ttn = rbind(res_ttn, data.frame(phenotype=echo_phenotypes[p], beta, se, pval, n=nrow(dat_eff)))
}
res_ttn$p_BH = p.adjust(res_ttn$pval, method = "BH")

print(res_bag3); print(res_ttn)


## 2. -----------------Heart RT only
all_data <- pheno_final.all[pheno_final.all$anthra_jco_dose_any ==0 & pheno_final.all$hrtavg > 0,]

pheno_final <- all_data
## Genotypes of TTN and BAG3 top SNPs
bag3 = read.table('sjlife_chr10_119670121_T_C.raw', header = TRUE)
bag3[c(2:6)] = NULL
ttn = read.table('sjlife_chr2_178562809_T_C.raw', header = TRUE)
ttn[c(2:6)] = NULL
geno = merge(bag3, ttn, by="FID")

## Merge phenotype and genotype data
dat_final = merge(pheno_final, geno, by = 'FID')

## Analyze both variants wrt each echo measures
# BAG3 SNP
res_bag3 = NULL
echo_phenotypes = colnames(echo_data)[-c(1:4)]
for (p in 1:length(echo_phenotypes)){
  dat_eff = subset(dat_final, !is.na(dat_final[echo_phenotypes[p]]))
  # Remove values with >3 sd away from the absolute value
  dat_eff = subset(dat_eff, abs(scale(dat_eff[,echo_phenotypes[p]])) <= 6)
  # Adjust for BSA; two (LVMass2D_Index and LA_Volume_Index) are already done
  if (!(echo_phenotypes[p] == "LVMass2D_Index" | echo_phenotypes[p] == "LA_Volume_Index"| 
        echo_phenotypes[p] == "LV_Ejection_Fraction_3D" | echo_phenotypes[p] == "LV_Cardiac_Output_3D" | echo_phenotypes[p] == "LV_GLPS_AVG")){
    fit = lm(scale(dat_eff[,echo_phenotypes[p]]/BodySurfaceArea)~chr10.119670121.T.C_C+agedx+agelstcontact+gender+anthra_jco_dose_any+hrtavg+
               PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10,data=dat_eff)
  } else {
    fit = lm(scale(dat_eff[,echo_phenotypes[p]])~chr10.119670121.T.C_C+agedx+agelstcontact+gender+anthra_jco_dose_any+hrtavg+
               PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10,data=dat_eff)
  }
  beta = summary(fit)$coef[2,1]
  se = summary(fit)$coef[2,2]
  pval = summary(fit)$coef[2,4]
  res_bag3 = rbind(res_bag3, data.frame(phenotype=echo_phenotypes[p], beta, se, pval, n=nrow(dat_eff)))
}
res_bag3$p_BH = p.adjust(res_bag3$pval, method = "BH")

# TTN SNP
res_ttn = NULL
for (p in 1:length(echo_phenotypes)){
  dat_eff = subset(dat_final, !is.na(dat_final[echo_phenotypes[p]]))
  # Remove values with >3 sd away from the absolute value
  dat_eff = subset(dat_eff, abs(scale(dat_eff[,echo_phenotypes[p]])) <= 6)
  # Adjust for BSA; two (LVMass2D_Index and LA_Volume_Index) are already done
  if (!(echo_phenotypes[p] == "LVMass2D_Index" | echo_phenotypes[p] == "LA_Volume_Index"| 
        echo_phenotypes[p] == "LV_Ejection_Fraction_3D" | echo_phenotypes[p] == "LV_Cardiac_Output_3D" | echo_phenotypes[p] == "LV_GLPS_AVG")){
    fit = lm(scale(dat_eff[,echo_phenotypes[p]]/BodySurfaceArea)~chr2.178562809.T.C_C+agedx+agelstcontact+gender+anthra_jco_dose_any+hrtavg+
               PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10,data=dat_eff)
  } else {
    fit = lm(scale(dat_eff[,echo_phenotypes[p]])~chr2.178562809.T.C_C+agedx+agelstcontact+gender+anthra_jco_dose_any+hrtavg+
               PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10,data=dat_eff)
  }
  beta = summary(fit)$coef[2,1]
  se = summary(fit)$coef[2,2]
  pval = summary(fit)$coef[2,4]
  res_ttn = rbind(res_ttn, data.frame(phenotype=echo_phenotypes[p], beta, se, pval, n=nrow(dat_eff)))
}
res_ttn$p_BH = p.adjust(res_ttn$pval, method = "BH")

print(res_bag3); print(res_ttn)

## 3. -----------------Anthracycline and Heart RT both
all_data <- pheno_final.all[pheno_final.all$anthra_jco_dose_any > 0 & pheno_final.all$hrtavg > 15,]

pheno_final <- all_data
## Genotypes of TTN and BAG3 top SNPs
bag3 = read.table('sjlife_chr10_119670121_T_C.raw', header = TRUE)
bag3[c(2:6)] = NULL
ttn = read.table('sjlife_chr2_178562809_T_C.raw', header = TRUE)
ttn[c(2:6)] = NULL
geno = merge(bag3, ttn, by="FID")

## Merge phenotype and genotype data
dat_final = merge(pheno_final, geno, by = 'FID')

## Analyze both variants wrt each echo measures
# BAG3 SNP
res_bag3 = NULL
echo_phenotypes = colnames(echo_data)[-c(1:4)]
for (p in 1:length(echo_phenotypes)){
  dat_eff = subset(dat_final, !is.na(dat_final[echo_phenotypes[p]]))
  # Remove values with >3 sd away from the absolute value
  dat_eff = subset(dat_eff, abs(scale(dat_eff[,echo_phenotypes[p]])) <= 6)
  # Adjust for BSA; two (LVMass2D_Index and LA_Volume_Index) are already done
  if (!(echo_phenotypes[p] == "LVMass2D_Index" | echo_phenotypes[p] == "LA_Volume_Index"| 
        echo_phenotypes[p] == "LV_Ejection_Fraction_3D" | echo_phenotypes[p] == "LV_Cardiac_Output_3D" | echo_phenotypes[p] == "LV_GLPS_AVG")){
    fit = lm(scale(dat_eff[,echo_phenotypes[p]]/BodySurfaceArea)~chr10.119670121.T.C_C+agedx+agelstcontact+gender+anthra_jco_dose_any+hrtavg+
               PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10,data=dat_eff)
  } else {
    fit = lm(scale(dat_eff[,echo_phenotypes[p]])~chr10.119670121.T.C_C+agedx+agelstcontact+gender+anthra_jco_dose_any+hrtavg+
               PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10,data=dat_eff)
  }
  beta = summary(fit)$coef[2,1]
  se = summary(fit)$coef[2,2]
  pval = summary(fit)$coef[2,4]
  res_bag3 = rbind(res_bag3, data.frame(phenotype=echo_phenotypes[p], beta, se, pval, n=nrow(dat_eff)))
}
res_bag3$p_BH = p.adjust(res_bag3$pval, method = "BH")

# TTN SNP
res_ttn = NULL
for (p in 1:length(echo_phenotypes)){
  dat_eff = subset(dat_final, !is.na(dat_final[echo_phenotypes[p]]))
  # Remove values with >3 sd away from the absolute value
  dat_eff = subset(dat_eff, abs(scale(dat_eff[,echo_phenotypes[p]])) <= 6)
  # Adjust for BSA; two (LVMass2D_Index and LA_Volume_Index) are already done
  if (!(echo_phenotypes[p] == "LVMass2D_Index" | echo_phenotypes[p] == "LA_Volume_Index"| 
        echo_phenotypes[p] == "LV_Ejection_Fraction_3D" | echo_phenotypes[p] == "LV_Cardiac_Output_3D" | echo_phenotypes[p] == "LV_GLPS_AVG")){
    fit = lm(scale(dat_eff[,echo_phenotypes[p]]/BodySurfaceArea)~chr2.178562809.T.C_C+agedx+agelstcontact+gender+anthra_jco_dose_any+hrtavg+
               PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10,data=dat_eff)
  } else {
    fit = lm(scale(dat_eff[,echo_phenotypes[p]])~chr2.178562809.T.C_C+agedx+agelstcontact+gender+anthra_jco_dose_any+hrtavg+
               PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10,data=dat_eff)
  }
  beta = summary(fit)$coef[2,1]
  se = summary(fit)$coef[2,2]
  pval = summary(fit)$coef[2,4]
  res_ttn = rbind(res_ttn, data.frame(phenotype=echo_phenotypes[p], beta, se, pval, n=nrow(dat_eff)))
}
res_ttn$p_BH = p.adjust(res_ttn$pval, method = "BH")

print(res_bag3); print(res_ttn)
