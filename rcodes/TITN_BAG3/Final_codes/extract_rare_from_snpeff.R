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
sjlife.bim <- sjlife.bim[!grepl(";", sjlife.bim$V2),]
sum(duplicated(sjlife.bim$KEY1))

ccss_exp.bim <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3/rare_variant_sjlife_ccss_combined/merged_plink_PLP_ccss_exp.bim")
ccss_exp.bim$KEY1 <- paste0(ccss_exp.bim$V1, ":", ccss_exp.bim$V4)
sum(duplicated(ccss_exp.bim$KEY1))

## Harmonize CCSS
sum(duplicated(ccss_exp.bim$KEY1))
## Harmonize
ccss_exp.bim$KEY <- paste0("chr", ccss_exp.bim$V1, ":", ccss_exp.bim$V4, ":", ccss_exp.bim$V6, ":", ccss_exp.bim$V5)
sum(duplicated(ccss_exp.bim$KEY))

keys <- ccss_exp.bim$KEY[duplicated(ccss_exp.bim$KEY)]
harmonize <- ccss_exp.bim[ccss_exp.bim$KEY %in% keys,]
harmonize <- harmonize[c("V2", "KEY")]



## Extract the overlapping variants
ccss_exp.bim$overlaps <- ccss_exp.bim$KEY1 %in% sjlife.bim$KEY1
ccss_exp.bim$overlapsV2 <- ccss_exp.bim$V2 %in% sjlife.bim$V2
ccss_exp.bim$sjlifeSNP <- sjlife.bim$V2[match(ccss_exp.bim$KEY1, sjlife.bim$KEY1)]
overlaps <-  ccss_exp.bim[ccss_exp.bim$overlaps == TRUE,]
table(overlaps$V2 %in% df$SNP_REF_ALT)
table(overlaps$V2 %in% df$SNP)
overlaps <- overlaps[overlaps$V2 %in% df$SNP,]

df.overlap <- df[df$SNP %in% overlaps$V2,]
dim(df.overlap)
# 128834    246

table(df.overlap$`ANN[*].GENE`)
df.overlap$`ANN[*].GENE`[df.overlap$`ANN[*].GENE` == "TTN-AS1"] <- "TTN"


## add TTN band and PSI
TTN <- df.overlap[df.overlap$`ANN[*].GENE` == "TTN",]
TTN <- cbind.data.frame(CHROM=TTN$CHROM, POS=TTN$POS, SNP=TTN$SNP, ID=TTN$ID)
ttn.pos <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects//Cardiotoxicity/common/ttn_bag3/TTN_band_and_PSI_roberts_et_al_2015.txt", sep = "\t", header = T)
split_data <- do.call(rbind, strsplit(as.character(ttn.pos$bed_GRCh38), "[:-]"))
ttn.pos <- data.frame(ttn.pos, chrom = split_data[,1], start = split_data[,2], end = split_data[,3])

TTN$Ave_PSI <- NA
TTN$Ave_PSI_GTEX <- NA
TTN$TTN_band <- NA

# Loop through each row of TTN
for(i in 1:nrow(TTN)) {
  print(paste0("Doing: ", i))
  # Find matching rows in ttn.pos where TTN$POS falls between ttn.pos$start and ttn.pos$end
  matching_row <- ttn.pos[TTN$POS[i] >= ttn.pos$start & TTN$POS[i] <= ttn.pos$end, ]
  
  # If there's a match, assign the corresponding Ave_PSI_GTEX value
  if(nrow(matching_row) > 0) {
    TTN$Ave_PSI[i] <- matching_row$Ave_PSI
    TTN$Ave_PSI_GTEX[i] <- matching_row$Ave_PSI_GTEX
    TTN$TTN_band[i] <- matching_row$Region
  }
}

# View the updated TTN dataframe
head(TTN)

df.overlap$Ave_PSI <- TTN$Ave_PSI[match(df.overlap$SNP, TTN$SNP)]
df.overlap$Ave_PSI_GTEX <- TTN$Ave_PSI_GTEX[match(df.overlap$SNP, TTN$SNP)]
df.overlap$TTN_band <- TTN$TTN_band[match(df.overlap$SNP, TTN$SNP)]

cc <- cbind.data.frame(df.overlap$SNP, df.overlap$`ANN[*].EFFECT`, df.overlap$AF, df.overlap$AF_nfe, df.overlap$AF_afr)

## NFE
AF_EUR.0.001 <- df.overlap[df.overlap$AF < 0.001 & df.overlap$AF_nfe < 0.001,]
table(AF_EUR.0.001$`ANN[*].EFFECT`)
AF_EUR.0.001 <- AF_EUR.0.001[grepl("stop|missense|frameshift|splice_acceptor", AF_EUR.0.001$`ANN[*].EFFECT`, ignore.case = T),]
# AF_EUR.0.001 <- AF_EUR.0.001[!grepl("intron", AF_EUR.0.001$`ANN[*].EFFECT`, ignore.case = T),]
AF_EUR.0.001.patho <- AF_EUR.0.001[grepl("^Pathogenic", AF_EUR.0.001$CLNSIG),]
AF_EUR.0.001 <- AF_EUR.0.001[!AF_EUR.0.001$SNP %in% AF_EUR.0.001.patho$SNP ,]
AF_EUR.0.001 <- rbind.data.frame(AF_EUR.0.001, AF_EUR.0.001.patho)
AF_EUR.0.001 <- AF_EUR.0.001[!duplicated(AF_EUR.0.001$ID),]
table(AF_EUR.0.001$`ANN[*].EFFECT`, AF_EUR.0.001$`ANN[*].GENE`)
#                                         BAG3 DSP LMNA MYH7 SCN5A TCAP TNNT2 TTN
# frameshift_variant                        0   0    0    0     0    0     1   0
# missense_variant                         13  28    8    9    23    5     4 484
# missense_variant&splice_region_variant    1   0    0    0     0    0     0  10
# splice_acceptor_variant&intron_variant    0   0    1    0     0    0     0   2
# stop_gained                               0   0    2    0     0    0     1   2
# stop_gained&splice_region_variant         0   0    0    0     0    0     1   0

table(AF_EUR.0.001$CLNSIG)
dim(AF_EUR.0.001)
# 595 249

# write.table(as.data.frame(AF_EUR.0.001$SNP), "Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3/rare_variant_sjlife_ccss_combined/snpEFF_AF_eur_0.001_vars.txt", col.names = F, row.names = F,  quote = F)

## AFR
AF_AFR.0.001 <- df.overlap[df.overlap$AF < 0.001 & df.overlap$AF_afr < 0.001,]
table(AF_AFR.0.001$`ANN[*].EFFECT`)
AF_AFR.0.001 <- AF_AFR.0.001[grepl("stop|missense|frameshift|splice_acceptor", AF_AFR.0.001$`ANN[*].EFFECT`),]
# AF_AFR.0.001 <- AF_AFR.0.001[!grepl("intron", AF_AFR.0.001$`ANN[*].EFFECT`, ignore.case = T),]
AF_AFR.0.001.patho <- AF_AFR.0.001[grepl("^Pathogenic", AF_AFR.0.001$CLNSIG),]
AF_AFR.0.001 <- AF_AFR.0.001[!AF_AFR.0.001$SNP %in% AF_AFR.0.001.patho$SNP ,]
AF_AFR.0.001 <- rbind.data.frame(AF_AFR.0.001, AF_AFR.0.001.patho)
AF_AFR.0.001 <- AF_AFR.0.001[!duplicated(AF_AFR.0.001$ID),]
table(AF_AFR.0.001$`ANN[*].EFFECT`, AF_AFR.0.001$`ANN[*].GENE`)
#                                         BAG3 DSP LMNA MYH7 SCN5A TCAP TNNT2 TTN
# frameshift_variant                        0   0    0    0     0    0     1   0
# missense_variant                         12  27    8   10    19    5     4 468
# missense_variant&splice_region_variant    1   0    0    0     0    0     0   9
# stop_gained                               0   0    1    0     0    0     1   2
# stop_gained&splice_region_variant         0   0    0    0     0    0     1   0
table(AF_AFR.0.001$CLNSIG)
dim(AF_AFR.0.001)
# 570 246

# 14590+18500 = 33090

# write.table(as.data.frame(AF_AFR.0.001$SNP), "Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3/rare_variant_sjlife_ccss_combined/snpEFF_AF_afr_0.001_vars.txt", col.names = F, row.names = F,  quote = F)


## EUR
raw <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3/rare_variant_sjlife_ccss_combined/sjlife_ccss_exp_SNPS_maf_lt_0.001_gnomad_recodeA_EUR.raw", header = T)
rownames(raw) <- raw$IID
raw <- raw[-c(1:6)]
# HEADER <- colnames(raw)[-c(1:6)]
HEADER=colnames(raw)
HEADER = sub(pattern="_[T,A,G,C,*]+",replacement="",HEADER)
HEADER = gsub(pattern=";rs\\d+",replacement="",HEADER)
HEADER = gsub("._...DEL", "",HEADER)
HEADER = gsub("....DEL.", ".<*:DEL>",HEADER)
HEADER = gsub("\\.", ":", HEADER)
colnames(raw) <- HEADER
table(AF_EUR.0.001$SNP %in% colnames(raw))
# FALSE  TRUE  ## removed due to geno 0.01
# 5   590

AF_EUR.0.001 <- AF_EUR.0.001[AF_EUR.0.001$SNP %in% colnames(raw),]
raw[raw == 2] <- 1
raw.eur <- raw
# dim(AF_EUR.0.001)
# [1] 854 246
keep.with.carriers <- names(colSums(raw.eur, na.rm = T))[colSums(raw.eur, na.rm = T) == 0]
raw.eur <- raw.eur[colnames(raw.eur) %in% keep.with.carriers]
dim(raw.eur)
# 3102  251
AF_EUR.0.001.carrier <- AF_EUR.0.001[AF_EUR.0.001$SNP %in% colnames(raw.eur),]
dim(AF_EUR.0.001)
# 169  249
## Variants with at least one carrier
table(AF_EUR.0.001.carrier$`ANN[*].GENE`)
# BAG3   DSP  LMNA  MYH7 SCN5A   TTN 
# 6     4     3     3     8   145 

AF_EUR.0.001$rsID <- sapply(strsplit(AF_EUR.0.001$ID, ";"), function(x) {
  rs_part <- x[grep("^rs", x)]  # Find the part that starts with 'rs'
  if(length(rs_part) > 0) return(rs_part) else return(NA)
})
AF_EUR.0.001.final <- AF_EUR.0.001[, c("CHROM", "POS", "rsID", "SNP", "REF", "ALT", "ANN[*].GENE", "ANN[*].EFFECT", 
                               "Ave_PSI", "Ave_PSI_GTEX", "TTN_band", "AF", "AF_nfe", "AF_afr")]


AF_EUR.0.001.final$carriers_in_EUR <- AF_EUR.0.001.final$SNP %in% AF_EUR.0.001.carrier$SNP

# same variants in AFR
write.table(as.data.frame(AF_EUR.0.001.final$SNP), "Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3/rare_variant_sjlife_ccss_combined/snpEFF_EUR_vars_final_0.001_vars.txt", col.names = F, row.names = F,  quote = F)
 
# ## AFR
# raw <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3/rare_variant_sjlife_ccss_combined/sjlife_SNPS_maf_lt_0.001_gnomad_recodeA_AFR.raw", header = T)
# rownames(raw) <- raw$IID
# raw <- raw[-c(1:6)]
# HEADER <- colnames(raw)
# HEADER = sub(pattern="_[T,A,G,C,*]+",replacement="",HEADER)
# HEADER = gsub(pattern=";rs\\d+",replacement="",HEADER)
# HEADER = gsub("._...DEL", "",HEADER)
# HEADER = gsub("....DEL.", ".<*:DEL>",HEADER)
# HEADER = gsub("\\.", ":", HEADER)
# colnames(raw) <- HEADER
# table(AF_AFR.0.001$SNP %in% colnames(raw))
# AF_AFR.0.001 <- AF_AFR.0.001[AF_AFR.0.001$SNP %in% colnames(raw),]
# raw[raw == 2] <- 1
# raw.afr <- raw
# dim(AF_AFR.0.001)
# # [1] 531 249
# keep.with.carriers <- names(colSums(raw.afr, na.rm = T))[colSums(raw.afr, na.rm = T) > 0]
# raw.afr <- raw.afr[colnames(raw.afr) %in% keep.with.carriers]
# dim(raw.afr)
# # 246  36
# AF_AFR.0.001.carrier <- AF_AFR.0.001[AF_AFR.0.001$SNP %in% colnames(raw.afr),]
# dim(AF_AFR.0.001.carrier)
# # 32  246
# ## Variants with at least one carrier
# table(AF_AFR.0.001.carrier$`ANN[*].GENE`)
# # BAG3   DSP  LMNA  MYH7 SCN5A TNNT2   TTN 
# # 2     9     3     2     5     1   124 

## AFR with the same variants
raw <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3/rare_variant_sjlife_ccss_combined/sjlife_SNPS_maf_lt_0.001_gnomad_recodeA_AFR.raw", header = T)
rownames(raw) <- raw$IID
raw <- raw[-c(1:6)]
HEADER <- colnames(raw)
HEADER = sub(pattern="_[T,A,G,C,*]+",replacement="",HEADER)
HEADER = gsub(pattern=";rs\\d+",replacement="",HEADER)
HEADER = gsub("._...DEL", "",HEADER)
HEADER = gsub("....DEL.", ".<*:DEL>",HEADER)
HEADER = gsub("\\.", ":", HEADER)
colnames(raw) <- HEADER
table(AF_EUR.0.001.final$SNP %in% colnames(raw))


raw[raw == 2] <- 1
raw.afr <- raw
dim(AF_EUR.0.001.final)
# [1] 531 249
keep.with.carriers <- names(colSums(raw.afr, na.rm = T))[colSums(raw.afr, na.rm = T) > 0]
raw.afr <- raw.afr[colnames(raw.afr) %in% keep.with.carriers]
dim(raw.afr)
# 246  36
AF_AFR.0.001.carrier <- AF_EUR.0.001.final[AF_EUR.0.001.final$SNP %in% colnames(raw.afr),]
dim(AF_AFR.0.001.carrier)
# 32  246
## Variants with at least one carrier
table(AF_AFR.0.001.carrier$`ANN[*].GENE`)
# BAG3   DSP  LMNA  MYH7 SCN5A  TTN
#   1     4     1     1     2    56

AF_EUR.0.001.final$carriers_in_AFR <- AF_EUR.0.001.final$SNP %in% AF_AFR.0.001.carrier$SNP

