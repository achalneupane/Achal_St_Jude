library(data.table)
df <- fread("Z:/ResearchHome/Groups/sapkogrp/projects//Cardiotoxicity/common/ttn_bag3/ALL_P_LP_combinations/rare_variant_analysis_v2/pablo_garcia_et_al_nine_genes_SNPEFF/NINE_GENES_ANNOVAR", sep = "\t", header = T)
# df <- fread("Z:/ResearchHome/Groups/sapkogrp/projects//Cardiotoxicity/common/ttn_bag3/ALL_P_LP_combinations/rare_variant_analysis_v2/pablo_garcia_et_al_nine_genes/NINE_GENES_ANNOVAR", sep = "\t", header = T)
df$AF <- as.numeric(df$AF)
df$AF_nfe <- as.numeric(df$AF_nfe)
df$AF_afr <- as.numeric(df$AF_afr)

df$SNP <- sub(";.*", "", df$ID)
cc <- cbind.data.frame(df$SNP, df$`ANN[*].EFFECT`, df$AF, df$AF_nfe, df$AF_afr)
df$SNP_REF_ALT <- paste0(df$CHROM,":", df$POS, ":", df$REF, ":", df$ALT)


sjlife.bim <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3/rare_variant_sjlife_ccss_combined/merged_plink_PLP_sjlife.bim")
sjlife.bim$KEY1 <- paste0(sjlife.bim$V1, ":", sjlife.bim$V4)
sum(duplicated(sjlife.bim$KEY1))

ccss_exp.bim <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3/rare_variant_sjlife_ccss_combined/merged_plink_PLP_ccss_exp.bim")
ccss_exp.bim$KEY1 <- paste0(ccss_exp.bim$V1, ":", ccss_exp.bim$V4)
sum(duplicated(ccss_exp.bim$KEY1))

## Extract the overlapping variants
ccss_exp.bim$overlaps <- ccss_exp.bim$KEY1 %in% sjlife.bim$KEY1
overlaps <-  ccss_exp.bim[ccss_exp.bim$overlaps == TRUE,]
table(overlaps$V2 %in% df$SNP_REF_ALT)
table(overlaps$V2 %in% df$SNP)
overlaps <- overlaps[overlaps$V2 %in% df$SNP,]

df.overlap <- df[df$SNP %in% overlaps$V2,]


## NFE
AF_NFE.0.01 <- df.overlap[df.overlap$AF < 0.01 & df.overlap$AF_nfe < 0.01,]
table(AF_NFE.0.01$`ANN[*].EFFECT`)
AF_NFE.0.01 <- AF_NFE.0.01[grepl("stop|missense|frameshift", AF_NFE.0.01$`ANN[*].EFFECT`, ignore.case = T),]
AF_NFE.0.01 <- AF_NFE.0.01[!grepl("intron", AF_NFE.0.01$`ANN[*].EFFECT`, ignore.case = T),]
AF_NFE.0.01.patho <- AF_NFE.0.01[grepl("^Pathogenic", AF_NFE.0.01$CLNSIG),]
AF_NFE.0.01 <- AF_NFE.0.01[!AF_NFE.0.01$SNP %in% AF_NFE.0.01.patho$SNP ,]
AF_NFE.0.01 <- rbind.data.frame(AF_NFE.0.01, AF_NFE.0.01.patho)
AF_NFE.0.01 <- AF_NFE.0.01[!duplicated(AF_NFE.0.01$ID),]
table(AF_NFE.0.01$`ANN[*].EFFECT`, AF_NFE.0.01$`ANN[*].GENE`)
# BAG3 DSP LMNA MYH7 SCN5A TCAP TNNT2 TTN
# frameshift_variant                        0   0    0    0     0    0     1   0
# missense_variant                         19  41   12   10    32    7     6 717
# missense_variant&splice_region_variant    1   0    0    0     0    0     2  12
# stop_gained                               0   0    2    0     0    0     1   2
# stop_gained&splice_region_variant         0   0    0    0     0    0     1   0

table(AF_NFE.0.01$CLNSIG)

# write.table(as.data.frame(AF_NFE.0.01$SNP), "Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3/rare_variant_sjlife_ccss_combined/snpEFF_AF_nfe_0.01_vars.txt", col.names = F, row.names = F,  quote = F)

raw <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3/rare_variant_sjlife_ccss_combined/sjlife_ccss_exp_SNPS_maf_lt_0.01_gnomad_recodeA_EUR.raw", header = T)
rownames(raw) <- raw$IID
raw <- raw[-c(1:6)]
# HEADER <- colnames(raw)[-c(1:6)]
HEADER = sub(pattern="_[T,A,G,C,*]+",replacement="",HEADER)
HEADER = gsub(pattern=";rs\\d+",replacement="",HEADER)
HEADER = gsub("._...DEL", "",HEADER)
HEADER = gsub("....DEL.", ".<*:DEL>",HEADER)
HEADER = gsub("\\.", ":", HEADER)
colnames(raw) <- HEADER
table(AF_AFR.0.01$SNP %in% colnames(raw))
AF_NFE.0.01 <- AF_NFE.0.01[AF_NFE.0.01$SNP %in% colnames(raw),]
raw[raw == 2] <- 1
raw.eur <- raw
# dim(AF_NFE.0.01)
# [1] 854 246
keep.with.carriers <- names(colSums(raw.eur, na.rm = T))[colSums(raw.eur, na.rm = T) == 0]
raw.eur <- raw.eur[colnames(raw.eur) %in% keep.with.carriers]
dim(raw.eur)
# 3102  251
AF_NFE.0.01 <- AF_NFE.0.01[AF_NFE.0.01$SNP %in% colnames(raw.eur),]
dim(AF_NFE.0.01)
# 230  246
## Variants with at least one carrier
table(AF_NFE.0.01$`ANN[*].GENE`)
# BAG3   DSP  LMNA  MYH7 SCN5A  TCAP TNNT2   TTN 
# 7     9     1     1     8     1     3   200 

## AFR
AF_AFR.0.01 <- df.overlap[df.overlap$AF < 0.01 & df.overlap$AF_afr < 0.01,]
table(AF_AFR.0.01$`ANN[*].EFFECT`)
AF_AFR.0.01 <- AF_AFR.0.01[grepl("stop|missense|frameshift", AF_AFR.0.01$`ANN[*].EFFECT`),]
AF_AFR.0.01 <- AF_AFR.0.01[!grepl("intron", AF_AFR.0.01$`ANN[*].EFFECT`, ignore.case = T),]
AF_AFR.0.01.patho <- AF_AFR.0.01[grepl("^Pathogenic", AF_AFR.0.01$CLNSIG),]
AF_AFR.0.01 <- AF_AFR.0.01[!AF_AFR.0.01$SNP %in% AF_AFR.0.01.patho$SNP ,]
AF_AFR.0.01 <- rbind.data.frame(AF_AFR.0.01, AF_AFR.0.01.patho)
AF_AFR.0.01 <- AF_AFR.0.01[!duplicated(AF_AFR.0.01$ID),]
table(AF_AFR.0.01$`ANN[*].EFFECT`, AF_AFR.0.01$`ANN[*].GENE`)
# BAG3 DSP LMNA MYH7 SCN5A TCAP TNNT2 TTN
# frameshift_variant                        0   0    0    0     0    0     1   0
# missense_variant                         16  41   11   11    30    6     6 679
# missense_variant&splice_region_variant    1   0    0    0     0    0     0  12
# stop_gained                               0   0    2    0     0    0     1   2
# stop_gained&splice_region_variant         0   0    0    0     0    0     1   0
table(AF_AFR.0.01$CLNSIG)

# 14590+18500 = 33090
# write.table(as.data.frame(AF_AFR.0.01$SNP), "Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3/rare_variant_sjlife_ccss_combined/snpEFF_AF_afr_0.01_vars.txt", col.names = F, row.names = F,  quote = F)

raw <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3/rare_variant_sjlife_ccss_combined/sjlife_ccss_exp_SNPS_maf_lt_0.01_gnomad_recodeA_AFR.raw", header = T)
rownames(raw) <- raw$IID
raw <- raw[-c(1:6)]
HEADER <- colnames(raw)
HEADER = sub(pattern="_[T,A,G,C,*]+",replacement="",HEADER)
HEADER = gsub(pattern=";rs\\d+",replacement="",HEADER)
HEADER = gsub("._...DEL", "",HEADER)
HEADER = gsub("....DEL.", ".<*:DEL>",HEADER)
HEADER = gsub("\\.", ":", HEADER)
colnames(raw) <- HEADER
table(AF_AFR.0.01$SNP %in% colnames(raw))
AF_AFR.0.01 <- AF_AFR.0.01[AF_AFR.0.01$SNP %in% colnames(raw),]
raw[raw == 2] <- 1
raw.afr <- raw
dim(AF_AFR.0.01)
# [1] 802 246
keep.with.carriers <- names(colSums(raw.afr, na.rm = T))[colSums(raw.afr, na.rm = T) > 0]
raw.afr <- raw.afr[colnames(raw.afr) %in% keep.with.carriers]
dim(raw.afr)
# 3102  146
AF_AFR.0.01 <- AF_AFR.0.01[AF_AFR.0.01$SNP %in% colnames(raw.afr),]
dim(AF_AFR.0.01)
# 146  246
## Variants with at least one carrier
table(AF_AFR.0.01$`ANN[*].GENE`)
# BAG3   DSP  LMNA  MYH7 SCN5A TNNT2   TTN 
# 2     9     3     2     5     1   124 





