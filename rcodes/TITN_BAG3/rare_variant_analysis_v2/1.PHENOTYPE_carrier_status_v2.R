## Read phenotype files
pheno.ccss_exp_eur <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3/pheno/ccss_exp_eur_cardiotoxic_exposed.pheno", header = T)
# pheno.ccss_org_eur <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3/pheno/ccss_org_eur_cardiotoxic_exposed.pheno", header = T)
pheno.sjlife_ttn_bag3 <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3/pheno/sjlife_ttn_bag3.pheno", header = T)
pheno.sjlife_ttn_bag3$CMP <- factor(pheno.sjlife_ttn_bag3$CMP)
#########################
## read carrier status ##
#########################
setwd("Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3/ALL_P_LP_combinations/rare_variant_analysis_v2/pablo_garcia_et_al_nine_genes")


load("common_p_LP_rare_variants_gnomad_all_gnomad_NFE_lt_0.01.RData")

setwd("Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3/ALL_P_LP_combinations/rare_variant_analysis_v2/pablo_garcia_et_al_nine_genes")
sjlife.raw <- read.table("sjlife_SNPS_maf_lt_0.01_gnomad_also_common_in_ccss_recodeA.raw", header = T)

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


################
## CCSS::TITN ##
################
TITN_ccss_exp.Clinv.LoF <- read.table("TITN_ccss_exp.0.01maf_overlapping.clinvar.lof_recodeA.raw", header = T)
colnames(TITN_ccss_exp.Clinv.LoF)[-c(1:6)] <- str_split(gsub("\\.", ":", colnames(TITN_ccss_exp.Clinv.LoF)[-c(1:6)]), "_", simplify=T)[,1]
# Remove splice region that were included initially
TITN_ccss_exp.Clinv.LoF <- TITN_ccss_exp.Clinv.LoF[c(1:6,which(colnames(TITN_ccss_exp.Clinv.LoF) %in% splice.region.keep$KEY))]

TITN_ccss_exp.Clinv.LoF$Non.REF.Counts <- rowSums(TITN_ccss_exp.Clinv.LoF[!colnames(TITN_ccss_exp.Clinv.LoF) %in% c("FID", "IID", "PAT", "MAT", "SEX", "PHENOTYPE")])
TITN_ccss_exp.Clinv.LoF$carrierSTATUS <- factor(ifelse(TITN_ccss_exp.Clinv.LoF$Non.REF.Counts >= 1, "Y", "N"))
CCSS.carrier <- cbind.data.frame(IID = TITN_ccss_exp.Clinv.LoF$IID, TITN.Clinv.LoF.Non.REF.Counts = TITN_ccss_exp.Clinv.LoF$Non.REF.Counts, TITN.Clinv.LoF.carrierSTATUS = TITN_ccss_exp.Clinv.LoF$carrierSTATUS)

# TITN_ccss_exp.Clinv.LoF.REVEL <- read.table("TITN_ccss_exp.0.01maf_overlapping.clinvar.lof.REVEL_recodeA.raw", header = T)
# colnames(TITN_ccss_exp.Clinv.LoF.REVEL)[-c(1:6)] <- str_split(gsub("\\.", ":", colnames(TITN_ccss_exp.Clinv.LoF.REVEL)[-c(1:6)]), "_", simplify=T)[,1]
# # Remove splice region that were included initially
# TITN_ccss_exp.Clinv.LoF.REVEL <- TITN_ccss_exp.Clinv.LoF.REVEL[c(1:6,which(colnames(TITN_ccss_exp.Clinv.LoF.REVEL) %in% splice.region.keep$KEY))]
# 
# CCSS.carrier$TITN.Clinv.LoF.REVEL.Non.REF.Counts <- rowSums(TITN_ccss_exp.Clinv.LoF.REVEL[!colnames(TITN_ccss_exp.Clinv.LoF.REVEL) %in% c("FID", "IID", "PAT", "MAT", "SEX", "PHENOTYPE")])
# CCSS.carrier$TITN.Clinv.LoF.REVEL.carrierSTATUS <- factor(ifelse(CCSS.carrier$TITN.Clinv.LoF.REVEL.Non.REF.Counts >= 1, "Y", "N"))

## CCSS::BAG3
BAG3_ccss_exp.Clinv.LoF <- read.table("BAG3_ccss_exp.0.01maf_overlapping.clinvar.lof_recodeA.raw", header = T)
colnames(BAG3_ccss_exp.Clinv.LoF)[-c(1:6)] <- str_split(gsub("\\.", ":", colnames(BAG3_ccss_exp.Clinv.LoF)[-c(1:6)]), "_", simplify=T)[,1]
# Remove splice region that were included initially
BAG3_ccss_exp.Clinv.LoF <- BAG3_ccss_exp.Clinv.LoF[c(1:6,which(colnames(BAG3_ccss_exp.Clinv.LoF) %in% splice.region.keep$KEY))]


CCSS.carrier$BAG3.Clinv.LoF.Non.REF.Counts <- rowSums(BAG3_ccss_exp.Clinv.LoF[!colnames(BAG3_ccss_exp.Clinv.LoF) %in% c("FID", "IID", "PAT", "MAT", "SEX", "PHENOTYPE")])
CCSS.carrier$BAG3.Clinv.LoF.carrierSTATUS <- factor(ifelse(CCSS.carrier$BAG3.Clinv.LoF.Non.REF.Counts >= 1, "Y", "N"))


# BAG3_ccss_exp.Clinv.LoF.REVEL <- read.table("BAG3_ccss_exp.0.01maf_overlapping.clinvar.lof.REVEL_recodeA.raw", header = T)
# colnames(BAG3_ccss_exp.Clinv.LoF.REVEL)[-c(1:6)] <- str_split(gsub("\\.", ":", colnames(BAG3_ccss_exp.Clinv.LoF.REVEL)[-c(1:6)]), "_", simplify=T)[,1]
# # Remove splice region that were included initially
# BAG3_ccss_exp.Clinv.LoF.REVEL <- BAG3_ccss_exp.Clinv.LoF.REVEL[c(1:6,which(colnames(BAG3_ccss_exp.Clinv.LoF.REVEL) %in% splice.region.keep$KEY))]

 
# CCSS.carrier$BAG3.Clinv.LoF.REVEL.Non.REF.Counts <- rowSums(BAG3_ccss_exp.Clinv.LoF.REVEL[!colnames(BAG3_ccss_exp.Clinv.LoF.REVEL) %in% c("FID", "IID", "PAT", "MAT", "SEX", "PHENOTYPE")])
# CCSS.carrier$BAG3.Clinv.LoF.REVEL.carrierSTATUS <- factor(ifelse(CCSS.carrier$BAG3.Clinv.LoF.REVEL.Non.REF.Counts >= 1, "Y", "N"))

## Add carrier status to PHENO
pheno.ccss_exp_eur <- cbind.data.frame(pheno.ccss_exp_eur, CCSS.carrier[match(pheno.ccss_exp_eur$IID, CCSS.carrier$IID),-1])

## Add carrier status for TITN + BAG3
pheno.ccss_exp_eur$TITN_BAG3.clin.LoF.Counts <- pheno.ccss_exp_eur$TITN.Clinv.LoF.Non.REF.Counts + pheno.ccss_exp_eur$BAG3.Clinv.LoF.Non.REF.Counts
pheno.ccss_exp_eur$TITN_BAG3.clin.LoF.carrierSTATUS <- factor(ifelse(pheno.ccss_exp_eur$TITN_BAG3.clin.LoF.Counts >= 1, "Y", "N"))

# pheno.ccss_exp_eur$TITN_BAG3.clin.LoF.REVEL.Counts <- pheno.ccss_exp_eur$TITN.Clinv.LoF.REVEL.Non.REF.Counts + pheno.ccss_exp_eur$BAG3.Clinv.LoF.REVEL.Non.REF.Counts
# pheno.ccss_exp_eur$TITN_BAG3.clin.LoF.REVEL.carrierSTATUS <- factor(ifelse(pheno.ccss_exp_eur$TITN_BAG3.clin.LoF.REVEL.Counts >= 1, "Y", "N"))


##################
## SJLIFE::TITN ##
##################
TITN_sjlife.Clinv.LoF <- read.table("TITN_sjlife.0.01maf_overlapping.clinvar.lof_recodeA.raw", header = T)
colnames(TITN_sjlife.Clinv.LoF)[-c(1:6)] <- str_split(gsub("\\.", ":", colnames(TITN_sjlife.Clinv.LoF)[-c(1:6)]), "_", simplify=T)[,1]
# Remove splice region that were included initially
TITN_sjlife.Clinv.LoF <- TITN_sjlife.Clinv.LoF[c(1:6,which(colnames(TITN_sjlife.Clinv.LoF) %in% splice.region.keep$KEY))]


TITN_sjlife.Clinv.LoF$Non.REF.Counts <- rowSums(TITN_sjlife.Clinv.LoF[!colnames(TITN_sjlife.Clinv.LoF) %in% c("FID", "IID", "PAT", "MAT", "SEX", "PHENOTYPE")])
TITN_sjlife.Clinv.LoF$carrierSTATUS <- factor(ifelse(TITN_sjlife.Clinv.LoF$Non.REF.Counts >= 1, "Y", "N"))
SJLIFE.carrier <- cbind.data.frame(IID = TITN_sjlife.Clinv.LoF$IID, TITN.Clinv.LoF.Non.REF.Counts = TITN_sjlife.Clinv.LoF$Non.REF.Counts, TITN.Clinv.LoF.carrierSTATUS = TITN_sjlife.Clinv.LoF$carrierSTATUS)

# TITN_sjlife.Clinv.LoF.REVEL <- read.table("TITN_sjlife.0.01maf_overlapping.clinvar.lof.REVEL_recodeA.raw", header = T)
# colnames(TITN_sjlife.Clinv.LoF.REVEL)[-c(1:6)] <- str_split(gsub("\\.", ":", colnames(TITN_sjlife.Clinv.LoF.REVEL)[-c(1:6)]), "_", simplify=T)[,1]
# # Remove splice region that were included initially
# TITN_sjlife.Clinv.LoF.REVEL <- TITN_sjlife.Clinv.LoF.REVEL[c(1:6,which(colnames(TITN_sjlife.Clinv.LoF.REVEL) %in% splice.region.keep$KEY))]

# SJLIFE.carrier$TITN.Clinv.LoF.REVEL.Non.REF.Counts <- rowSums(TITN_sjlife.Clinv.LoF.REVEL[!colnames(TITN_sjlife.Clinv.LoF.REVEL) %in% c("FID", "IID", "PAT", "MAT", "SEX", "PHENOTYPE")])
# SJLIFE.carrier$TITN.Clinv.LoF.REVEL.carrierSTATUS <- factor(ifelse(SJLIFE.carrier$TITN.Clinv.LoF.REVEL.Non.REF.Counts >= 1, "Y", "N"))

## SJLIFE::BAG3
BAG3_sjlife.Clinv.LoF <- read.table("BAG3_sjlife.0.01maf_overlapping.clinvar.lof_recodeA.raw", header = T)
colnames(BAG3_sjlife.Clinv.LoF)[-c(1:6)] <- str_split(gsub("\\.", ":", colnames(BAG3_sjlife.Clinv.LoF)[-c(1:6)]), "_", simplify=T)[,1]
# Remove splice region that were included initially
BAG3_sjlife.Clinv.LoF <- BAG3_sjlife.Clinv.LoF[c(1:6,which(colnames(BAG3_sjlife.Clinv.LoF) %in% splice.region.keep$KEY))]


SJLIFE.carrier$BAG3.Clinv.LoF.Non.REF.Counts <- rowSums(BAG3_sjlife.Clinv.LoF[!colnames(BAG3_sjlife.Clinv.LoF) %in% c("FID", "IID", "PAT", "MAT", "SEX", "PHENOTYPE")])
SJLIFE.carrier$BAG3.Clinv.LoF.carrierSTATUS <- factor(ifelse(SJLIFE.carrier$BAG3.Clinv.LoF.Non.REF.Counts >= 1, "Y", "N"))


# # BAG3_sjlife.Clinv.LoF.REVEL <- read.table("BAG3_sjlife.0.01maf_overlapping.clinvar.lof.REVEL_recodeA.raw", header = T)
# # colnames(BAG3_sjlife.Clinv.LoF.REVEL)[-c(1:6)] <- str_split(gsub("\\.", ":", colnames(BAG3_sjlife.Clinv.LoF.REVEL)[-c(1:6)]), "_", simplify=T)[,1]
# # # Remove splice region that were included initially
# # BAG3_sjlife.Clinv.LoF.REVEL <- BAG3_sjlife.Clinv.LoF.REVEL[c(1:6,which(colnames(BAG3_sjlife.Clinv.LoF.REVEL) %in% splice.region.keep$KEY))]
# 
# 
# SJLIFE.carrier$BAG3.Clinv.LoF.REVEL.Non.REF.Counts <- rowSums(BAG3_sjlife.Clinv.LoF.REVEL[!colnames(BAG3_sjlife.Clinv.LoF.REVEL) %in% c("FID", "IID", "PAT", "MAT", "SEX", "PHENOTYPE")])
# SJLIFE.carrier$BAG3.Clinv.LoF.REVEL.carrierSTATUS <- factor(ifelse(SJLIFE.carrier$BAG3.Clinv.LoF.REVEL.Non.REF.Counts >= 1, "Y", "N"))

## Add carrier status to PHENO
pheno.sjlife_ttn_bag3 <- cbind.data.frame(pheno.sjlife_ttn_bag3, SJLIFE.carrier[match(pheno.sjlife_ttn_bag3$IID, SJLIFE.carrier$IID),-1])


## Add carrier status for TITN + BAG3
pheno.sjlife_ttn_bag3$TITN_BAG3.clin.LoF.Counts <- pheno.sjlife_ttn_bag3$TITN.Clinv.LoF.Non.REF.Counts + pheno.sjlife_ttn_bag3$BAG3.Clinv.LoF.Non.REF.Counts
pheno.sjlife_ttn_bag3$TITN_BAG3.clin.LoF.carrierSTATUS <- factor(ifelse(pheno.sjlife_ttn_bag3$TITN_BAG3.clin.LoF.Counts >= 1, "Y", "N"))

# pheno.sjlife_ttn_bag3$TITN_BAG3.clin.LoF.REVEL.Counts <- pheno.sjlife_ttn_bag3$TITN.Clinv.LoF.REVEL.Non.REF.Counts + pheno.sjlife_ttn_bag3$BAG3.Clinv.LoF.REVEL.Non.REF.Counts
# pheno.sjlife_ttn_bag3$TITN_BAG3.clin.LoF.REVEL.carrierSTATUS <- factor(ifelse(pheno.sjlife_ttn_bag3$TITN_BAG3.clin.LoF.REVEL.Counts >= 1, "Y", "N"))

# table(TITN_ccss_exp.Clinv.LoF$IID == TITN_ccss_exp.Clinv.LoF.REVEL$IID)
# table(BAG3_ccss_exp.Clinv.LoF$IID == BAG3_ccss_exp.Clinv.LoF.REVEL$IID)
# table(TITN_ccss_exp.Clinv.LoF$IID == BAG3_ccss_exp.Clinv.LoF.REVEL$IID)

# table(TITN_sjlife.Clinv.LoF$IID == TITN_sjlife.Clinv.LoF.REVEL$IID)
# table(BAG3_sjlife.Clinv.LoF$IID == BAG3_sjlife.Clinv.LoF.REVEL$IID)
# table(TITN_sjlife.Clinv.LoF$IID == BAG3_sjlife.Clinv.LoF$IID)

save.image("Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3/analysis_model_fitting/phenotype_carrier_status_v2.RDATA")

