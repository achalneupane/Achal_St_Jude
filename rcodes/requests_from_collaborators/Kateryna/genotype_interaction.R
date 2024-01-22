setwd("Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/gwas")
df <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/gwas/pheno/sjlife_eur_dox_only.txt", header = T)
dim(df)
df <- df[df$IID %in% genotype$IID,]

genotype <- read.table("chr16_genotypes_eur.raw", header = T)
df <- cbind.data.frame(df,genotype[na.omit(match(df$IID, genotype$IID)), c("chr16.25608384.A.C_C", "chr16.25609966.C.T_T", "chr16.25610265.C.T_T", "chr16.25611595.C.T_T")])

df$CMP2plus <- factor(df$CMP2plus)

fit <- glm(CMP2plus ~ chr16.25611595.C.T_T + anthra_jco_dose_any + chr16.25611595.C.T_T*anthra_jco_dose_any, family = binomial, data = df)  
summary(fit)


# haplo.test.adj <- glm(formula = CMP ~  haplo_111101 + agedx + agelstcontact + gender + anthra_jco_dose_any + hrtavg + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, 
#                       family = binomial,
#                       data = pheno)