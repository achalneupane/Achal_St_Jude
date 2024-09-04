## Read phenotype files
pheno.sjlife_ttn_bag3 <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3/pheno/sjlife_ccss_exp_ttn_bag3.pheno", header = T)
pheno.sjlife_ttn_bag3$CMP <- factor(pheno.sjlife_ttn_bag3$CMP)
#########################
## read carrier status ##
#########################
setwd("Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3/ALL_P_LP_combinations/rare_variant_analysis_v2/pablo_garcia_et_al_nine_genes")


load("p_LP_rare_variants_gnomad_all_gnomad_NFE_lt_0.01.RData")

setwd("Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3/rare_variant_sjlife_ccss_combined")
sjlife.raw <- read.table("sjlife_ccss_exp_SNPS_maf_lt_0.01_gnomad_recodeA.raw", header = T)


variants.raw <- colnames(sjlife.raw)[-c(1:6)]
HEADER = sub(pattern="_[T,A,G,C,*]+",replacement="",variants.raw)
HEADER = gsub(pattern=";rs\\d+",replacement="",HEADER)
HEADER = gsub("._...DEL", "",HEADER)
HEADER = gsub("....DEL.", ".<*:DEL>",HEADER)
HEADER = gsub("\\.", ":", HEADER)

sum(sjlife_vars_bim$SNP %in% HEADER)
# 30
sjlife_vars_bim <- sjlife_vars_bim[sjlife_vars_bim$SNP %in% HEADER,]

sum(HEADER %in% sjlife_vars_bim$SNP)
# HEADER [HEADER %in% sjlife_vars_bim$SNP]

df <- read.table(text = "chr3:38562422:C:A
chr1:156104763:A:C
chr1:156134954:G:A
chr1:156138821:C:T
chr1:201361214:C:T
chr1:201361215:G:A
chr1:201361215:G:C
chr1:201364335:CG:C
chr2:178542266:A:G
chr2:178577602:C:T
chr2:178615321:A:G
chr2:178630241:G:A
chr2:178646559:T:C
chr2:178651534:G:T
chr2:178664926:G:C
chr2:178665777:G:A
chr2:178667719:A:T
chr2:178669675:G:T
chr2:178692016:C:T
chr2:178693609:C:T
chr2:178701528:C:G
chr2:178702065:T:A
chr2:178706956:T:G
chr2:178728776:G:T
chr2:178746079:G:A
chr2:178756220:A:G
chr2:178777856:A:G
chr2:178778004:G:A
chr2:178782806:C:T
chr10:119672255:C:T
", header = F)

df$V1[!df$V1 %in% HEADER]


genes <- unique(sjlife_vars_bim$GENE)

i=1:5
## SJLIFE analysis
for (i in 1:length(genes)){
  genes[i]
  rm(total.carriers, total.cases, carriers.cases, total.controls, carriers.controls, OR.CI, pvalue, OR.CI.adj, pvalue.adj)
  genes.carriers <-  sjlife.raw[c(2,grep(paste0(gsub(":" , "\\.",paste0("chr", sjlife_vars_bim$KEY[grepl(genes[i], sjlife_vars_bim$GENE)], ".")), collapse = "|"), colnames(sjlife.raw)))]
  # If checking for All PLP
  # genes.carriers <-  sjlife.raw[c(2,grep(paste0(gsub(":" , "\\.",paste0("chr", sjlife_vars_bim$KEY[(sjlife_vars_bim$GENE %in% genes)], ".")), collapse = "|"), colnames(sjlife.raw)))]
  # See how many variants had carrier: Were there carriers in our data for all
  # these variants? For example, if no one had variant for chr10:119672255:C:T,
  # you would not have included in the generation of carrier status. This is what
  # we include in the total numbers of variants used to perform the analysis
  check.carriers <- genes.carriers[-1]
  check.carriers
  # df <- check.carriers[, complete.cases(check.carriers)]
  # df <- check.carriers[, colSums(is.na(check.carriers)) == 0]
  # "chr2.178542266.A.G_G" "chr2.178615321.A.G_G" "chr2.178630241.G.A_A" "chr2.178782806.C.T_T"
  cols_with_ones <- na.omit(apply(check.carriers == 1, 2, any))
  cols_with_ones
  # Extract the names of the columns with at least one 1
  n_vars_with_carriers <- length(names(cols_with_ones)[cols_with_ones])
  print(paste0(genes[i],"=",  n_vars_with_carriers))
# }
  ###############################################
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


# write.table(sjlife_vars_bim, "overlapping_rare_variants_in_sjlife_ccss_maf_0.01.txt", sep = "\t", col.names = T, row.names = F, quote = F)

#############################
## Run analysis in African ##
#############################
load("Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3/Rcodes/echo.PLP.afr.RData")


## Carriers in African
sjlife.LMNA.PLP <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3//ALL_P_LP_combinations/rare_variant_analysis_v2/pablo_garcia_et_al_nine_genes/LMNA_PLP_SJLIFE_recodeA.raw", header = T)
# sjlife.LMNA.PLP$carrier <- ifelse(rowSums(sjlife.LMNA.PLP[grepl("chr", colnames(sjlife.LMNA.PLP))]) > 0, 1, 0)
sjlife.LMNA.PLP <- sjlife.LMNA.PLP[sjlife.LMNA.PLP$IID %in% echo.PLP.afr$sjlid ,]
check.carriers <-  sjlife.LMNA.PLP[-c(1:6)]
cols_with_ones <- na.omit(apply(check.carriers == 1, 2, any))
cols_with_ones
# Extract the names of the columns with at least one 1
n_vars_with_carriers <- length(names(cols_with_ones)[cols_with_ones])
print(paste0("LMNA","=",  n_vars_with_carriers))
# "LMNA=0"


sjlife.TNTT2.PLP <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3//ALL_P_LP_combinations/rare_variant_analysis_v2/pablo_garcia_et_al_nine_genes/TNTT2_PLP_SJLIFE_recodeA.raw", header = T)
# sjlife.TNTT2.PLP$carrier <- ifelse(rowSums(sjlife.TNTT2.PLP[grepl("chr", colnames(sjlife.TNTT2.PLP))]) > 0, 1, 0)
sjlife.TNTT2.PLP <- sjlife.TNTT2.PLP[sjlife.TNTT2.PLP$IID %in% echo.PLP.afr$sjlid ,]
check.carriers <-  sjlife.TNTT2.PLP[-c(1:6)]
cols_with_ones <- na.omit(apply(check.carriers == 1, 2, any))
cols_with_ones
# Extract the names of the columns with at least one 1
n_vars_with_carriers <- length(names(cols_with_ones)[cols_with_ones])
print(paste0("TNTT2","=",  n_vars_with_carriers))
# "TNTT2=1"

sjlife.TTN.PLP <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3//ALL_P_LP_combinations/rare_variant_analysis_v2/pablo_garcia_et_al_nine_genes/TTN_PLP_SJLIFE_recodeA.raw", header = T)
# sjlife.TTN.PLP$carrier <- ifelse(rowSums(sjlife.TTN.PLP[grepl("chr", colnames(sjlife.TTN.PLP))]) > 0, 1, 0)
sjlife.TTN.PLP <- sjlife.TTN.PLP[sjlife.TTN.PLP$IID %in% echo.PLP.afr$sjlid ,]
check.carriers <-  sjlife.TTN.PLP[-c(1:6)]
cols_with_ones <- na.omit(apply(check.carriers == 1, 2, any))
cols_with_ones
# Extract the names of the columns with at least one 1
n_vars_with_carriers <- length(names(cols_with_ones)[cols_with_ones])
print(paste0("TTN","=",  n_vars_with_carriers))
# "TTN=2"

sjlife.BAG3.PLP <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3//ALL_P_LP_combinations/rare_variant_analysis_v2/pablo_garcia_et_al_nine_genes/BAG3_PLP_SJLIFE_recodeA.raw", header = T)
# sjlife.BAG3.PLP$carrier <- ifelse(rowSums(sjlife.BAG3.PLP[grepl("chr", colnames(sjlife.BAG3.PLP))]) > 0, 1, 0)
sjlife.BAG3.PLP <- sjlife.BAG3.PLP[sjlife.BAG3.PLP$IID %in% echo.PLP.afr$sjlid ,]
check.carriers <-  sjlife.BAG3.PLP[-c(1:6)]
cols_with_ones <- na.omit(apply(check.carriers == 1, 2, any))
cols_with_ones
# Extract the names of the columns with at least one 1
n_vars_with_carriers <- length(names(cols_with_ones)[cols_with_ones])
print(paste0("BAG3","=",  n_vars_with_carriers))
# "BAG3=0"

sjlife.ALL.PLP <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3//ALL_P_LP_combinations/rare_variant_analysis_v2/pablo_garcia_et_al_nine_genes/ALL_PLP_SJLIFE_recodeA.raw", header = T)
# sjlife.ALL.PLP$carrier <- ifelse(rowSums(sjlife.ALL.PLP[grepl("chr", colnames(sjlife.ALL.PLP))]) > 0, 1, 0)
sjlife.ALL.PLP <- sjlife.ALL.PLP[sjlife.ALL.PLP$IID %in% echo.PLP.afr$sjlid ,]
check.carriers <-  sjlife.ALL.PLP[-c(1:6)]
cols_with_ones <- na.omit(apply(check.carriers == 1, 2, any))
cols_with_ones
# Extract the names of the columns with at least one 1
n_vars_with_carriers <- length(names(cols_with_ones)[cols_with_ones])
print(paste0("ALL","=",  n_vars_with_carriers))


sjlife$LMNA.PLP_carrier <- sjlife.LMNA.PLP$carrier[match(sjlife$IID, sjlife.LMNA.PLP$IID)]
sjlife$TNTT2.PLP_carrier <- sjlife.TNTT2.PLP$carrier[match(sjlife$IID, sjlife.TNTT2.PLP$IID)]
sjlife$TTN.PLP_carrier <- sjlife.TTN.PLP$carrier[match(sjlife$IID, sjlife.TTN.PLP$IID)]
sjlife$BAG3.PLP_carrier <- sjlife.BAG3.PLP$carrier[match(sjlife$IID, sjlife.BAG3.PLP$IID)]
sjlife$ALL.PLP_carrier <- sjlife.ALL.PLP$carrier[match(sjlife$IID, sjlife.ALL.PLP$IID)]



pheno.sjlife_ttn_bag3 <- as.data.frame(echo.PLP.afr[echo.PLP.afr$visittype == "SJLIFE Visit 1",])
genes <- c("LMNA.PLP_carrier", "TNTT2.PLP_carrier", "TTN.PLP_carrier", "BAG3.PLP_carrier", "ALL.PLP_carrier")
pheno.sjlife_ttn_bag3$CMP[pheno.sjlife_ttn_bag3$CMP == 1 ] <- 2
pheno.sjlife_ttn_bag3$CMP[pheno.sjlife_ttn_bag3$CMP == 0 ] <- 1

i=1:5
## SJLIFE analysis
# for (i in 1:length(genes)){
  genes[i]
  rm(total.carriers, total.cases, carriers.cases, total.controls, carriers.controls, OR.CI, pvalue, OR.CI.adj, pvalue.adj)
  pheno.sjlife_ttn_bag3$carriers <-  pheno.sjlife_ttn_bag3[,genes[i]]
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