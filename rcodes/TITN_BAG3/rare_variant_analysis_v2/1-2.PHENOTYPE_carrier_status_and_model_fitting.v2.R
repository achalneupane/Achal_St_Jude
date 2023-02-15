## Read phenotype files
pheno.ccss_exp_eur <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3/pheno/ccss_exp_eur_cardiotoxic_exposed.pheno", header = T)
# pheno.ccss_org_eur <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3/pheno/ccss_org_eur_cardiotoxic_exposed.pheno", header = T)
pheno.sjlife_ttn_bag3 <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3/pheno/sjlife_ttn_bag3.pheno", header = T)
pheno.sjlife_ttn_bag3$CMP <- factor(pheno.sjlife_ttn_bag3$CMP)
pheno.ccss_exp_eur$CMP <- factor(pheno.ccss_exp_eur$CMP2plus)
#########################
## read carrier status ##
#########################
setwd("Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3/ALL_P_LP_combinations/rare_variant_analysis_v2/pablo_garcia_et_al_nine_genes")


load("common_p_LP_rare_variants_gnomad_all_gnomad_NFE_lt_0.01.RData")

setwd("Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3/ALL_P_LP_combinations/rare_variant_analysis_v2/pablo_garcia_et_al_nine_genes")
sjlife.raw <- read.table("sjlife_SNPS_maf_lt_0.01_gnomad_also_common_in_ccss_recodeA.raw", header = T)
# ccss.raw <- read.table("ccss_SNPS_maf_lt_0.01_gnomad_also_common_in_sjlife_recodeA.raw", header = T)

# exclude younger samples from ccss exp
ccss.raw <- read.table("ccss_SNPS_maf_lt_0.01_gnomad_also_common_in_sjlife_without_younger_samples_recodeA.raw", header = T)
sum(pheno.ccss_exp_eur$FID %in% ccss.raw$FID)
# 1457
pheno.ccss_exp_eur <- pheno.ccss_exp_eur[pheno.ccss_exp_eur$FID %in% ccss.raw$FID,]


genes <- unique(sjlife_vars_bim$GENE)





## SJLIFE analysis

for (i in 1:length(genes)){
  genes[i]
  rm(total.carriers, total.cases, carriers.cases, total.controls, carriers.controls, OR.CI, pvalue, OR.CI.adj, pvalue.adj)
  genes.carriers <-  sjlife.raw[c(2,grep(paste0(gsub(":" , "\\.",paste0("chr", sjlife_vars_bim$KEY[grepl(genes[i], sjlife_vars_bim$GENE)], ".")), collapse = "|"), colnames(sjlife.raw)))]
  genes.carriers$carrier <- factor(ifelse(rowSums(genes.carriers[-1]) != 0, "1", "0"))
  # table(genes.carriers$carrier)
  pheno.sjlife_ttn_bag3$carriers <- genes.carriers$carrier[match(pheno.sjlife_ttn_bag3$IID, genes.carriers$IID)]
  rm(gene.test)
  gene.test <- fisher.test(table(pheno.sjlife_ttn_bag3$CMP, pheno.sjlife_ttn_bag3$carriers)) # Got error; no carriers
  # gene.test <- chisq.test(table(pheno.sjlife_ttn_bag3$CMP, pheno.sjlife_ttn_bag3$carriers)) # Got error; no carriers
  
  pvalue <- gene.test$p.value
  OR.CI <- paste0(round(gene.test$estimate, 2), " (", paste0(round(gene.test$conf.int, 2), collapse = "-"), ")")
  
  rm(gene.test.adj)
  gene.test.adj <- glm(formula = CMP ~  pheno.sjlife_ttn_bag3$carriers + agedx + agelstcontact + gender + anthra_jco_dose_any + hrtavg + PC1 + PC2 + PC3 +
                         PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, family = binomial,
                       data = pheno.sjlife_ttn_bag3)
  
  
  CI1 <- paste0(" (",round(exp(c(coef(summary(gene.test.adj))[2,1]))-1.96*as.data.frame(summary(gene.test.adj)$coefficients)[c(1,2)][2,2],2), " to ", round(exp(c(coef(summary(gene.test.adj))[2,1]))+1.96*as.data.frame(summary(gene.test.adj)$coefficients)[c(1,2)][2,2], 2),")")
  
  total.carriers <- sum(pheno.sjlife_ttn_bag3$carriers == 1, na.rm = T)
  total.cases <- sum(pheno.sjlife_ttn_bag3$CMP == 2, na.rm = T)
  carriers.cases <- sum(pheno.sjlife_ttn_bag3$CMP == 2 & pheno.sjlife_ttn_bag3$carriers == 1, na.rm = T)
  total.controls <- sum(pheno.sjlife_ttn_bag3$CMP == 1, na.rm = T)
  carriers.controls <- sum(pheno.sjlife_ttn_bag3$CMP == 1 & pheno.sjlife_ttn_bag3$carriers == 1, na.rm = T)
  OR.CI.adj <- paste0(round(exp(c(coef(summary(gene.test.adj))[2,1])),2), CI1)
  pvalue.adj <- coef(summary(gene.test.adj))[2,4]
  # cbind.data.frame(genes[i],total.carriers, total.cases, carriers.cases, total.controls, carriers.controls, OR.CI, pvalue, OR.CI.adj, pvalue.adj)
  
  genes[i]
  total.carriers
  total.cases
  carriers.cases
  total.controls
  carriers.controls 
  OR.CI
  pvalue
  OR.CI.adj
  pvalue.adj
  
}

# ## Prevalence is defined as
# table.prevalence = table(pheno.sjlife_ttn_bag3$CMP, pheno.sjlife_ttn_bag3$carriers)
# #      0    1  # 1,2 is CMP and 0,1 is carrier status
# # 1 1402   28
# # 2  204    6
# prevalence = table.prevalence[1,2]/sum(table.prevalence[1,])
# # i.e., prevalence in CO is carriers in CO/ total CO; prevalence in CA is carriers in CA/ total CA

## CCSS analysis
for (i in 1:length(genes)){
  genes[i]
  rm(total.carriers, total.cases, carriers.cases, total.controls, carriers.controls, OR.CI, pvalue, OR.CI.adj, pvalue.adj)
  genes.carriers <-  ccss.raw[c(2,grep(paste0(gsub(":" , "\\.",paste0(sjlife_vars_bim$CCSS_equivalent[grepl(genes[i], sjlife_vars_bim$GENE)], ".")), collapse = "|"), colnames(ccss.raw)))]
  genes.carriers$carrier <- factor(ifelse(rowSums(genes.carriers[-1]) != 0, "1", "0"))
  # table(genes.carriers$carrier)
  pheno.ccss_exp_eur$carriers <- genes.carriers$carrier[match(pheno.ccss_exp_eur$IID, genes.carriers$IID)]
  rm(gene.test)
  gene.test <- fisher.test(table(pheno.ccss_exp_eur$CMP, pheno.ccss_exp_eur$carriers)) # Got error; no carriers
  # gene.test <- chisq.test(table(pheno.ccss_exp_eur$CMP, pheno.ccss_exp_eur$carriers)) # Got error; no carriers
  
  pvalue <- gene.test$p.value
  OR.CI <- paste0(round(gene.test$estimate, 2), " (", paste0(round(gene.test$conf.int, 2), collapse = "-"), ")")
  
  rm(gene.test.adj)
  gene.test.adj <- glm(formula = CMP ~  pheno.ccss_exp_eur$carriers + a_dx + a_end + SEX + anth_DED + HeartAvg + PC1 + PC2 + PC3 +
                       PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, family = binomial,
                       data = pheno.ccss_exp_eur)
  
  
  CI1 <- paste0(" (",round(exp(c(coef(summary(gene.test.adj))[2,1]))-1.96*as.data.frame(summary(gene.test.adj)$coefficients)[c(1,2)][2,2],2), " to ", round(exp(c(coef(summary(gene.test.adj))[2,1]))+1.96*as.data.frame(summary(gene.test.adj)$coefficients)[c(1,2)][2,2], 2),")")
  
  total.carriers <- sum(pheno.ccss_exp_eur$carriers == 1, na.rm = T)
  total.cases <- sum(pheno.ccss_exp_eur$CMP == 2, na.rm = T)
  carriers.cases <- sum(pheno.ccss_exp_eur$CMP == 2 & pheno.ccss_exp_eur$carriers == 1, na.rm = T)
  total.controls <- sum(pheno.ccss_exp_eur$CMP == 1, na.rm = T)
  carriers.controls <- sum(pheno.ccss_exp_eur$CMP == 1 & pheno.ccss_exp_eur$carriers == 1, na.rm = T)
  OR.CI.adj <- paste0(round(exp(c(coef(summary(gene.test.adj))[2,1])),2), CI1)
  pvalue.adj <- coef(summary(gene.test.adj))[2,4]
  # cbind.data.frame(genes[i],total.carriers, total.cases, carriers.cases, total.controls, carriers.controls, OR.CI, pvalue, OR.CI.adj, pvalue.adj)
  
  genes[i]
  total.carriers
  total.cases
  carriers.cases
  total.controls
  carriers.controls 
  OR.CI
  pvalue
  OR.CI.adj
  pvalue.adj
  
}

write.table(sjlife_vars_bim, "overlapping_rare_variants_in_sjlife_ccss_maf_0.01.txt", sep = "\t", col.names = T, row.names = F, quote = F)
