# One technical concern is the MAF1% cutoff that the authors used for rare P/LP
# variants. MAF cutoff often used are 0.01%. Also, many of rare TTN variants
# reported in Table S4 are in I-band (in addition to having MAF>0.0001), which
# are not generally considered pathogenic/likely pathogenic. Mostly, variants in
# exons with high PSI, often located in A-band are considered P/LP. For example,
# chr2:178665777 A>G in Table S4 is in exon with PSI 1% and won’t be considered
# pathogenic. Therefore, I have concerns about their rare genetic variant
# analysis in the manuscript.
setwd("Z:/ResearchHome/Groups/sapkogrp/projects//Cardiotoxicity/common/ttn_bag3//rare_variant_sjlife_ccss_combined/")

library(data.table)
# df <- fread("Z:/ResearchHome/Groups/sapkogrp/projects//Cardiotoxicity/common/ttn_bag3/ALL_P_LP_combinations/rare_variant_analysis_v2/pablo_garcia_et_al_nine_genes_SNPEFF/NINE_GENES_ANNOVAR", sep = "\t", header = T)
## df <- fread("Z:/ResearchHome/Groups/sapkogrp/projects//Cardiotoxicity/common/ttn_bag3/ALL_P_LP_combinations/rare_variant_analysis_v2/pablo_garcia_et_al_nine_genes/NINE_GENES_ANNOVAR", sep = "\t", header = T)
df <- fread("Z:/ResearchHome/Groups/sapkogrp/projects//Cardiotoxicity/common/ttn_bag3//rare_variant_sjlife_ccss_combined/nine_genes_merged_snpeff_ann.txt", sep = "\t", header = T)
df$AF <- as.numeric(df$AF)
df$AF_nfe <- as.numeric(df$AF_nfe)
df$AF_afr <- as.numeric(df$AF_afr)

df$SNP <- sub(";.*", "", df$ID)
cc <- cbind.data.frame(df$SNP, df$`ANN[*].EFFECT`, df$AF, df$AF_nfe, df$AF_afr)
df$SNP_REF_ALT <- paste0(df$CHROM,":", df$POS, ":", df$REF, ":", df$ALT)

df$rsID <- sub(".*;(rs[0-9]+).*", "\\1", df$ID)


# sjlife.bim <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3/rare_variant_sjlife_ccss_combined/merged_plink_PLP_sjlife.bim")
# sjlife.bim$KEY1 <- paste0(sjlife.bim$V1, ":", sjlife.bim$V4)
# sjlife.bim <- sjlife.bim[!grepl(";", sjlife.bim$V2),]
# sum(duplicated(sjlife.bim$KEY1))
# 
# ccss_exp.bim <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3/rare_variant_sjlife_ccss_combined/merged_plink_PLP_ccss_exp.bim")
# ccss_exp.bim$KEY1 <- paste0(ccss_exp.bim$V1, ":", ccss_exp.bim$V4)
# sum(duplicated(ccss_exp.bim$KEY1))
# 
# ## Harmonize CCSS
# sum(duplicated(ccss_exp.bim$KEY1))
# ## Harmonize
# ccss_exp.bim$KEY <- paste0("chr", ccss_exp.bim$V1, ":", ccss_exp.bim$V4, ":", ccss_exp.bim$V6, ":", ccss_exp.bim$V5)
# sum(duplicated(ccss_exp.bim$KEY))
# 
# keys <- ccss_exp.bim$KEY[duplicated(ccss_exp.bim$KEY)]
# harmonize <- ccss_exp.bim[ccss_exp.bim$KEY %in% keys,]
# harmonize <- harmonize[c("V2", "KEY")]
# 
# 
# 
# ## Extract the overlapping variants
# ccss_exp.bim$overlaps <- ccss_exp.bim$KEY1 %in% sjlife.bim$KEY1
# ccss_exp.bim$overlapsV2 <- ccss_exp.bim$V2 %in% sjlife.bim$V2
# ccss_exp.bim$sjlifeSNP <- sjlife.bim$V2[match(ccss_exp.bim$KEY1, sjlife.bim$KEY1)]
# overlaps <-  ccss_exp.bim[ccss_exp.bim$overlaps == TRUE,]
# table(overlaps$V2 %in% df$SNP_REF_ALT)
# table(overlaps$V2 %in% df$SNP)
# overlaps <- overlaps[overlaps$V2 %in% df$SNP,]
# 
# df.overlap <- df[df$SNP %in% overlaps$V2,]
# dim(df.overlap)
# # 128834    246


df.overlap <- df
table(df.overlap$`ANN[*].GENE`)
df.overlap$`ANN[*].GENE`[df.overlap$`ANN[*].GENE` == "TTN-AS1"] <- "TTN"

# Remove those with either AF_afr >= 0.0001 0r AF_nfe >= 0.0001
df.overlap <- df.overlap[which(df.overlap$AF_nfe < 0.0001|df.overlap$AF_afr < 0.0001),]

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
TTN.saved <- TTN
head(TTN)
dim(TTN)
# 81220
# 185775 

df.overlap$Ave_PSI <- TTN$Ave_PSI[match(df.overlap$SNP, TTN$SNP)]
df.overlap$Ave_PSI_GTEX <- TTN$Ave_PSI_GTEX[match(df.overlap$SNP, TTN$SNP)]
df.overlap$TTN_band <- TTN$TTN_band[match(df.overlap$SNP, TTN$SNP)]

cc <- cbind.data.frame(df.overlap$SNP, df.overlap$`ANN[*].GENE`, df.overlap$`ANN[*].EFFECT`, df.overlap$AF, df.overlap$AF_nfe, df.overlap$AF_afr, df.overlap$rsID)




# 1. Protein-Altering Variants (BAG3, LMNA, TCAP, TNNC1, TNNT2)
# These include missense, nonsense, and frameshift variants:
# Missense variants: missense_variant
# Nonsense variants: stop_gained
# Frameshift variants: frameshift_variant
# In-frame insertions/deletions: inframe_insertion, inframe_deletion

# 2. Missense Variants and In-frame Insertion/Deletion (MYH7)
# Missense variants: missense_variant
# In-frame variants: inframe_insertion, inframe_deletion

# 3. Frameshift, Stop-Gained, Splice-Donor, and Splice-Acceptor Variants (DSP, SCN5A, TTN)
# Frameshift variants: frameshift_variant
# Stop-gained variants: stop_gained
# Splice-donor/acceptor variants: splice_donor_variant, splice_acceptor_variant

# 1. Protein-Altering Variants (BAG3, LMNA, TCAP, TNNC1, TNNT2)
df.overlap$POS <- as.numeric(df.overlap$POS) 
df.overlap.saved <- df.overlap
BAG3 <- df.overlap[(df.overlap$POS >=119651380 & df.overlap$POS <=119677819) & df.overlap$CHROM == "chr10",]
BAG3 <- BAG3[grepl("missense_variant|stop_gained|frameshift_variant|inframe_insertion|inframe_deletion", BAG3$`ANN[*].EFFECT`),]
BAG3 <- BAG3[!duplicated(BAG3$SNP),]

LMNA <- df.overlap[(df.overlap$POS >=156082573 & df.overlap$POS <=156140081) & df.overlap$CHROM == "chr1",]
LMNA <- LMNA[grepl("missense_variant|stop_gained|frameshift_variant|inframe_insertion|inframe_deletion", LMNA$`ANN[*].EFFECT`),]
LMNA <- LMNA[!duplicated(LMNA$SNP),]

TCAP <- df.overlap[(df.overlap$POS >=39665349 & df.overlap$POS <=39666554) & df.overlap$CHROM == "chr17",]
TCAP <- TCAP[grepl("missense_variant|stop_gained|frameshift_variant|inframe_insertion|inframe_deletion", TCAP$`ANN[*].EFFECT`),]
TCAP <- TCAP[!duplicated(TCAP$SNP),]

TNNC1 <- df.overlap[(df.overlap$POS >=52451100 & df.overlap$POS <=52454041) & df.overlap$CHROM == "chr3",]
TNNC1 <- TNNC1[grepl("missense_variant|stop_gained|frameshift_variant|inframe_insertion|inframe_deletion", TNNC1$`ANN[*].EFFECT`),]
TNNC1 <- TNNC1[!duplicated(TNNC1$SNP),]

TNNT2 <- df.overlap[(df.overlap$POS >=201359014 & df.overlap$POS <=201377680) & df.overlap$CHROM == "chr1",]
TNNT2 <- TNNT2[grepl("missense_variant|stop_gained|frameshift_variant|inframe_insertion|inframe_deletion", TNNT2$`ANN[*].EFFECT`),]
TNNT2 <- TNNT2[!duplicated(TNNT2$SNP),]

# 2. Missense Variants and In-frame Insertion/Deletion (MYH7)
MYH7 <- df.overlap[(df.overlap$POS >=23412740 & df.overlap$POS <=23435660) & df.overlap$CHROM == "chr14",]
MYH7 <- MYH7[grepl("missense_variant|stop_gained|inframe_insertion|inframe_deletion", MYH7$`ANN[*].EFFECT`),]
MYH7 <- MYH7[!duplicated(MYH7$SNP),]

# 3. Frameshift, Stop-Gained, Splice-Donor, and Splice-Acceptor Variants (DSP, SCN5A, TTN)
DSP <- df.overlap[(df.overlap$POS >=7541671 & df.overlap$POS <=7586714) & df.overlap$CHROM == "chr6",]
DSP <- DSP[grepl("frameshift_variant|stop_gained|splice_donor_variant|splice_acceptor_variant", DSP$`ANN[*].EFFECT`),]
DSP <- DSP[!duplicated(DSP$SNP),]

SCN5A <- df.overlap[(df.overlap$POS >=38548062 & df.overlap$POS <=38649687) & df.overlap$CHROM == "chr3",]
SCN5A <- SCN5A[grepl("frameshift_variant|stop_gained|splice_donor_variant|splice_acceptor_variant", SCN5A$`ANN[*].EFFECT`),]
SCN5A <- SCN5A[!duplicated(SCN5A$SNP),]

TTN <- df.overlap[(df.overlap$POS >=178525989 & df.overlap$POS <=178807423) & df.overlap$CHROM == "chr2",]
TTN <- TTN[grepl("frameshift_variant|stop_gained|splice_donor_variant|splice_acceptor_variant", TTN$`ANN[*].EFFECT`),]
TTN <- TTN[!duplicated(TTN$SNP),]
df.overlap <- rbind.data.frame(BAG3, LMNA, TCAP, TNNC1, TNNT2, MYH7, DSP, SCN5A, TTN)

cc <- cbind.data.frame(df.overlap$SNP, df.overlap$`ANN[*].GENE`, df.overlap$`ANN[*].EFFECT`, df.overlap$AF, df.overlap$AF_nfe, df.overlap$AF_afr, df.overlap$rsID, df.overlap$Ave_PSI_GTEX, df.overlap$TTN_band)

## NFE (apply only ancestry specific maf)
# AF_EUR.0.0001 <- df.overlap[df.overlap$AF < 0.0001 & df.overlap$AF_nfe < 0.0001,]
AF_EUR.0.0001 <- df.overlap[df.overlap$AF_nfe < 0.0001,]
table(AF_EUR.0.0001$`ANN[*].EFFECT`)
# AF_EUR.0.0001 <- AF_EUR.0.0001[grepl("stop|missense|frameshift|splice_acceptor", AF_EUR.0.0001$`ANN[*].EFFECT`, ignore.case = T),]
# AF_EUR.0.0001 <- AF_EUR.0.0001[!grepl("intron", AF_EUR.0.0001$`ANN[*].EFFECT`, ignore.case = T),]
# AF_EUR.0.0001.patho <- AF_EUR.0.0001[grepl("^Pathogenic", AF_EUR.0.0001$CLNSIG),]
# AF_EUR.0.0001 <- AF_EUR.0.0001[!AF_EUR.0.0001$SNP %in% AF_EUR.0.0001.patho$SNP ,]
# AF_EUR.0.0001 <- rbind.data.frame(AF_EUR.0.0001, AF_EUR.0.0001.patho)
# AF_EUR.0.0001 <- AF_EUR.0.0001[!duplicated(AF_EUR.0.0001$ID),]
table(AF_EUR.0.0001$`ANN[*].EFFECT`, AF_EUR.0.0001$`ANN[*].GENE`)
#                                         BAG3 DSP LMNA MYH7 SCN5A TCAP TNNC1 TNNT2 TTN
# disruptive_inframe_deletion                 1   0    0    0     0    0     0     0   0
# disruptive_inframe_insertion                1   0    0    0     0    0     0     0   0
# frameshift_variant                          0   0    1    0     0    1     0     0   4
# frameshift_variant&splice_region_variant    0   0    1    0     0    0     0     0   0
# missense_variant                           39   0   41   67     0   13     3    15   0
# missense_variant&splice_region_variant      1   0    0    4     0    1     0     0   0
# splice_acceptor_variant&intron_variant      0   0    0    0     0    0     0     0   2
# splice_donor_variant&intron_variant         0   1    0    0     0    0     0     0   6
# stop_gained                                 0   0    2    0     1    0     0     1  14

table(AF_EUR.0.0001$CLNSIG)
dim(AF_EUR.0.0001)
# 220 250
cc <- cbind.data.frame(AF_EUR.0.0001$SNP, AF_EUR.0.0001$`ANN[*].GENE`, AF_EUR.0.0001$`ANN[*].EFFECT`, AF_EUR.0.0001$AF, AF_EUR.0.0001$AF_nfe, AF_EUR.0.0001$AF_afr, AF_EUR.0.0001$rsID, AF_EUR.0.0001$Ave_PSI_GTEX, AF_EUR.0.0001$TTN_band)

write.table(as.data.frame(AF_EUR.0.0001$SNP), "Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3/rare_variant_sjlife_ccss_combined/snpEFF_AF_eur_0.0001_vars.txt", col.names = F, row.names = F,  quote = F)

## AFR
# AF_AFR.0.0001 <- df.overlap[df.overlap$AF < 0.0001 & df.overlap$AF_afr < 0.0001,]
AF_AFR.0.0001 <- df.overlap[df.overlap$AF_afr < 0.0001,]
table(AF_AFR.0.0001$`ANN[*].EFFECT`)
# AF_AFR.0.0001 <- AF_AFR.0.0001[grepl("stop|missense|frameshift|splice_acceptor", AF_AFR.0.0001$`ANN[*].EFFECT`),]
# AF_AFR.0.0001 <- AF_AFR.0.0001[!grepl("intron", AF_AFR.0.0001$`ANN[*].EFFECT`, ignore.case = T),]
# AF_AFR.0.0001.patho <- AF_AFR.0.0001[grepl("^Pathogenic", AF_AFR.0.0001$CLNSIG),]
# AF_AFR.0.0001 <- AF_AFR.0.0001[!AF_AFR.0.0001$SNP %in% AF_AFR.0.0001.patho$SNP ,]
# AF_AFR.0.0001 <- rbind.data.frame(AF_AFR.0.0001, AF_AFR.0.0001.patho)
# AF_AFR.0.0001 <- AF_AFR.0.0001[!duplicated(AF_AFR.0.0001$ID),]
table(AF_AFR.0.0001$`ANN[*].EFFECT`, AF_AFR.0.0001$`ANN[*].GENE`)
#                                                                               BAG3 DSP LMNA MYH7 SCN5A TCAP TNNT2 TTN
# frameshift_variant                                                              0   0    0    0     0    1     0     1   5
# frameshift_variant&splice_donor_variant&splice_region_variant&intron_variant    0   0    0    0     0    0     0     0   1
# frameshift_variant&splice_region_variant                                        0   0    1    0     0    0     0     0   0
# missense_variant                                                               38   0   41   71     0   12     3    17   0
# missense_variant&splice_region_variant                                          0   0    0    4     0    1     0     0   0
# splice_acceptor_variant&intron_variant                                          0   0    0    0     0    0     0     0   1
# splice_donor_variant&intron_variant                                             0   1    0    0     0    0     0     0   6
# stop_gained                                                                     0   0    1    0     1    0     0     1  16
# stop_gained&splice_region_variant                                               0   0    0    0     0    0     0     1   0

table(AF_AFR.0.0001$CLNSIG)
dim(AF_AFR.0.0001)
# 224 250

# 14590+18500 = 33090

write.table(as.data.frame(AF_AFR.0.0001$SNP), "Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3/rare_variant_sjlife_ccss_combined/snpEFF_AF_afr_0.0001_vars.txt", col.names = F, row.names = F,  quote = F)


## EUR
raw <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3/rare_variant_sjlife_ccss_combined/sjlife_ccss_exp_SNPS_maf_lt_0.0001_gnomad_recodeA_EUR.raw", header = T)
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
table(AF_EUR.0.0001$SNP %in% colnames(raw))
# TRUE 
# 220 

AF_EUR.0.0001 <- AF_EUR.0.0001[AF_EUR.0.0001$SNP %in% colnames(raw),]
raw[raw == 2] <- 1
raw.eur <- raw
dim(AF_EUR.0.0001)
# [1] 220 250
keep.with.carriers <- names(colSums(raw.eur, na.rm = T))[colSums(raw.eur, na.rm = T) > 0]
raw.eur <- raw.eur[colnames(raw.eur) %in% keep.with.carriers]
dim(raw.eur)
# 3102  50
AF_EUR.0.0001.carrier <- AF_EUR.0.0001[AF_EUR.0.0001$SNP %in% colnames(raw.eur),]
dim(AF_EUR.0.0001.carrier)
# 50  250
## Variants with at least one carrier
table(AF_EUR.0.0001.carrier$`ANN[*].GENE`)
# BAG3   DSP  LMNA  MYH7 SCN5A  TCAP TNNT2   TTN 
#    16     1    23    22     6     1     6    12

AF_EUR.0.0001.final <- AF_EUR.0.0001[, c("CHROM", "POS", "rsID", "SNP", "REF", "ALT", "ANN[*].GENE", "ANN[*].EFFECT", 
                               "Ave_PSI", "Ave_PSI_GTEX", "TTN_band", "AF", "AF_nfe", "AF_afr")]


AF_EUR.0.0001.final$var_in_EUR <- AF_EUR.0.0001.final$SNP %in% AF_EUR.0.0001.carrier$SNP
# N carriers
AF_EUR.0.0001.final$N_carriers_in_EUR <- NA
for (i in 1:nrow(AF_EUR.0.0001.final)){
AF_EUR.0.0001.final$N_carriers_in_EUR[i] <-  sum(raw[AF_EUR.0.0001.final$SNP[i]]==1, na.rm = T)
}
## same variants in AFR
# write.table(as.data.frame(AF_EUR.0.0001.final$SNP), "Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3/rare_variant_sjlife_ccss_combined/snpEFF_EUR_vars_final_0.0001_vars.txt", col.names = F, row.names = F,  quote = F)
 
# ## AFR
raw <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3/rare_variant_sjlife_ccss_combined/sjlife_SNPS_maf_lt_0.0001_gnomad_recodeA_AFR.raw", header = T)
rownames(raw) <- raw$IID
raw <- raw[-c(1:6)]
HEADER <- colnames(raw)
HEADER = sub(pattern="_[T,A,G,C,*]+",replacement="",HEADER)
HEADER = gsub(pattern=";rs\\d+",replacement="",HEADER)
HEADER = gsub("._...DEL", "",HEADER)
HEADER = gsub("....DEL.", ".<*:DEL>",HEADER)
HEADER = gsub("\\.", ":", HEADER)
colnames(raw) <- HEADER
table(AF_AFR.0.0001$SNP %in% colnames(raw))
# 143
AF_AFR.0.0001 <- AF_AFR.0.0001[AF_AFR.0.0001$SNP %in% colnames(raw),]
raw[raw == 2] <- 1
raw.afr <- raw
dim(AF_AFR.0.0001)
# [1] 143 250
keep.with.carriers <- names(colSums(raw.afr, na.rm = T))[colSums(raw.afr, na.rm = T) > 0]
raw.afr <- raw.afr[colnames(raw.afr) %in% keep.with.carriers]
dim(raw.afr)
# 246  11
AF_AFR.0.0001.carrier <- AF_AFR.0.0001[AF_AFR.0.0001$SNP %in% colnames(raw.afr),]
dim(AF_AFR.0.0001.carrier)
# 11  250
## Variants with at least one carrier
table(AF_AFR.0.0001.carrier$`ANN[*].GENE`)
# BAG3  LMNA  MYH7 TNNT2   TTN 
# 1     2     5     2     1 

AF_AFR.0.0001.final <- AF_AFR.0.0001[, c("CHROM", "POS", "rsID", "SNP", "REF", "ALT", "ANN[*].GENE", "ANN[*].EFFECT", 
                                         "Ave_PSI", "Ave_PSI_GTEX", "TTN_band", "AF", "AF_nfe", "AF_afr")]


AF_AFR.0.0001.final$var_in_AFR <- AF_AFR.0.0001.final$SNP %in% AF_AFR.0.0001.carrier$SNP
# N carriers
AF_AFR.0.0001.final$N_carriers_in_AFR <- NA
for (i in 1:nrow(AF_AFR.0.0001.final)){
  AF_AFR.0.0001.final$N_carriers_in_AFR[i] <-  sum(raw[AF_AFR.0.0001.final$SNP[i]]==1, na.rm = T)
}


# ## AFR with the same variants
# raw <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3/rare_variant_sjlife_ccss_combined/sjlife_SNPS_maf_lt_0.0001_gnomad_recodeA_AFR.raw", header = T)
# rownames(raw) <- raw$IID
# raw <- raw[-c(1:6)]
# HEADER <- colnames(raw)
# HEADER = sub(pattern="_[T,A,G,C,*]+",replacement="",HEADER)
# HEADER = gsub(pattern=";rs\\d+",replacement="",HEADER)
# HEADER = gsub("._...DEL", "",HEADER)
# HEADER = gsub("....DEL.", ".<*:DEL>",HEADER)
# HEADER = gsub("\\.", ":", HEADER)
# colnames(raw) <- HEADER
# table(AF_EUR.0.0001.final$SNP %in% colnames(raw))
# 
# 
# raw[raw == 2] <- 1
# raw.afr <- raw
# dim(raw.afr)
# # 246  181
# dim(AF_EUR.0.0001.final)
# # [1] 181 16
# keep.with.carriers <- names(colSums(raw.afr, na.rm = T))[colSums(raw.afr, na.rm = T) > 0]
# raw.afr <- raw.afr[colnames(raw.afr) %in% keep.with.carriers]
# dim(raw.afr)
# # 246  11
# AF_AFR.0.0001.carrier <- AF_EUR.0.0001.final[AF_EUR.0.0001.final$SNP %in% colnames(raw.afr),]
# dim(AF_AFR.0.0001.carrier)
# # 11  16
# ## Variants with at least one carrier
# table(AF_AFR.0.0001.carrier$`ANN[*].GENE`)
# # DSP TTN 
# #   1  10

# AF_EUR.0.0001.final$var_in_AFR <- AF_EUR.0.0001.final$SNP %in% AF_AFR.0.0001.carrier$SNP
# AF_EUR.0.0001.final$N_carriers_in_AFR <- NA
# for (i in 1:nrow(AF_EUR.0.0001.final)){
#   AF_EUR.0.0001.final$N_carriers_in_AFR[i] <-  sum(raw[AF_EUR.0.0001.final$SNP[i]]==1, na.rm = T)
# }
# 
# table(AF_EUR.0.0001.final$`ANN[*].EFFECT`)

AF_EUR.0.0001.final$`ANN[*].EFFECT`[AF_EUR.0.0001.final$`ANN[*].EFFECT` == "missense_variant&splice_region_variant"] <- "missense"
AF_EUR.0.0001.final$`ANN[*].EFFECT`[AF_EUR.0.0001.final$`ANN[*].EFFECT` == "splice_acceptor_variant&intron_variant"] <- "splice acceptor"
AF_EUR.0.0001.final$`ANN[*].EFFECT`[AF_EUR.0.0001.final$`ANN[*].EFFECT` == "stop_gained&splice_region_variant"] <- "stop-gained"
AF_EUR.0.0001.final$`ANN[*].EFFECT` <- gsub("disruptive-|&splice-region|&intron", "", AF_EUR.0.0001.final$`ANN[*].EFFECT`)
AF_EUR.0.0001.final$`ANN[*].EFFECT` <- gsub("_variant", "", AF_EUR.0.0001.final$`ANN[*].EFFECT`)
AF_EUR.0.0001.final$`ANN[*].EFFECT` <- gsub("_", "-", AF_EUR.0.0001.final$`ANN[*].EFFECT`)
table(AF_EUR.0.0001.final$`ANN[*].EFFECT`)

AF_AFR.0.0001.final$`ANN[*].EFFECT`[AF_AFR.0.0001.final$`ANN[*].EFFECT` == "missense_variant&splice_region_variant"] <- "missense"
AF_AFR.0.0001.final$`ANN[*].EFFECT`[AF_AFR.0.0001.final$`ANN[*].EFFECT` == "splice_acceptor_variant&intron_variant"] <- "splice acceptor"
AF_AFR.0.0001.final$`ANN[*].EFFECT`[AF_AFR.0.0001.final$`ANN[*].EFFECT` == "stop_gained&splice_region_variant"] <- "stop-gained"
AF_AFR.0.0001.final$`ANN[*].EFFECT` <- gsub("disruptive-|&splice-region|&intron|&splice-donor", "", AF_AFR.0.0001.final$`ANN[*].EFFECT`)
AF_AFR.0.0001.final$`ANN[*].EFFECT` <- gsub("_variant", "", AF_AFR.0.0001.final$`ANN[*].EFFECT`)
AF_AFR.0.0001.final$`ANN[*].EFFECT` <- gsub("_", "-", AF_AFR.0.0001.final$`ANN[*].EFFECT`)
table(AF_AFR.0.0001.final$`ANN[*].EFFECT`)

eur.freq <- fread("sjlife_ccss_exp_merged_EUR.freq.frq", header = T)
afr.freq <- fread("sjlife_AFR.freq.frq", header = T)


###################################################################
## EUR
AF_EUR.0.0001.final$cohort_EUR_MAF <- eur.freq$MAF[match(AF_EUR.0.0001.final$SNP, eur.freq$SNP)]
AF_EUR.0.0001.final$TTN_PSI <- ""
AF_EUR.0.0001.final$TTN_PSI [which(AF_EUR.0.0001.final$Ave_PSI_GTEX > 0.82)] <- "Yes"
AF_EUR.0.0001.final$TTN_PSI_A_Band <- ""
AF_EUR.0.0001.final$TTN_PSI_A_Band [which(AF_EUR.0.0001.final$Ave_PSI_GTEX > 0.82 &  AF_EUR.0.0001.final$TTN_band == "A-band")] <- "Yes"

# write.table(AF_EUR.0.0001.final, "Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3/rare_variant_sjlife_ccss_combined/snpEFF_EUR_vars_final_0.0001_annotation.txt", col.names = T, row.names = F,  quote = F, sep = "\t")


# with.carriers <- AF_EUR.0.0001.final[AF_EUR.0.0001.final$N_carriers_in_EUR >= 1,]
# table(with.carriers$`ANN[*].GENE`)
# # DSP  LMNA  MYH7 SCN5A  TCAP   TTN 
# # 4     4     1     2     1    42
# View(as.data.frame(with.carriers$SNP[with.carriers$`ANN[*].GENE` == "DSP"]))
# View(as.data.frame(with.carriers$SNP[with.carriers$`ANN[*].GENE` == "LMNA"]))
# View(as.data.frame(with.carriers$SNP[with.carriers$`ANN[*].GENE` == "MYH7"]))
# View(as.data.frame(with.carriers$SNP[with.carriers$`ANN[*].GENE` == "SCN5A"]))
# View(as.data.frame(with.carriers$SNP[with.carriers$`ANN[*].GENE` == "TCAP"]))
# View(as.data.frame(with.carriers$SNP[with.carriers$`ANN[*].GENE` == "TTN"]))
# View(as.data.frame(with.carriers$SNP[with.carriers$`ANN[*].GENE` == "TTN" & with.carriers$TTN_PSI=="Yes"]))


# extract variants by gene categories
table(AF_EUR.0.0001.final$`ANN[*].GENE`)
# BAG3   DSP  LMNA  MYH7 SCN5A  TCAP TNNT2   TTN 
# 42     1    45    71     1    15     3    16    26

genes <- names(table(AF_EUR.0.0001.final$`ANN[*].GENE`))
genes <- genes[!genes %in% "TTN"]
for (i in 1:length(genes)){
gene <- AF_EUR.0.0001.final$SNP[AF_EUR.0.0001.final$`ANN[*].GENE` == genes[i]]
write.table(as.data.frame(gene), paste0("Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3/rare_variant_sjlife_ccss_combined/", genes[i], "_EUR_vars.txt"), col.names = F, row.names = F,  quote = F, sep = "\t")
}

TTN_PSI <- AF_EUR.0.0001.final$SNP[AF_EUR.0.0001.final$TTN_PSI == "Yes"]
length(TTN_PSI)
# 11
write.table(as.data.frame(TTN_PSI), "Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3/rare_variant_sjlife_ccss_combined/TTN_PSI_EUR_vars.txt", col.names = F, row.names = F,  quote = F, sep = "\t")

TTN_PSI_A_Band <- AF_EUR.0.0001.final$SNP[AF_EUR.0.0001.final$TTN_PSI_A_Band == "Yes"]
length(TTN_PSI_A_Band)
# 1
write.table(as.data.frame(TTN_PSI_A_Band), "Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3/rare_variant_sjlife_ccss_combined/TTN_PSI_A_Band_EUR_vars.txt", col.names = F, row.names = F,  quote = F, sep = "\t")


BAG3.TTN_PSI <- AF_EUR.0.0001.final$SNP[AF_EUR.0.0001.final$`ANN[*].GENE` == "BAG3" | AF_EUR.0.0001.final$TTN_PSI == "Yes"]
length(BAG3.TTN_PSI)
# 53
write.table(as.data.frame(BAG3.TTN_PSI), "Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3/rare_variant_sjlife_ccss_combined/BAG3.TTN_PSI_EUR_vars.txt", col.names = F, row.names = F,  quote = F, sep = "\t")


BAG3.TTN_PSI_A_Band <- AF_EUR.0.0001.final$SNP[AF_EUR.0.0001.final$`ANN[*].GENE` == "BAG3" | AF_EUR.0.0001.final$TTN_PSI_A_Band == "Yes"]
length(BAG3.TTN_PSI_A_Band)
# 43
write.table(as.data.frame(BAG3.TTN_PSI_A_Band), "Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3/rare_variant_sjlife_ccss_combined/BAG3.TTN_PSI_A_Band_EUR_vars.txt", col.names = F, row.names = F,  quote = F, sep = "\t")


##################################################################
## AFR
AF_AFR.0.0001.final$cohort_AFR_MAF <- afr.freq$MAF[match(AF_AFR.0.0001.final$SNP, afr.freq$SNP)]

AF_AFR.0.0001.final$TTN_PSI <- ""
AF_AFR.0.0001.final$TTN_PSI [which(AF_AFR.0.0001.final$Ave_PSI_GTEX > 0.82)] <- "Yes"
AF_AFR.0.0001.final$TTN_PSI_A_Band <- ""
AF_AFR.0.0001.final$TTN_PSI_A_Band [which(AF_AFR.0.0001.final$Ave_PSI_GTEX > 0.82 &  AF_AFR.0.0001.final$TTN_band == "A-band")] <- "Yes"

# write.table(AF_AFR.0.0001.final, "Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3/rare_variant_sjlife_ccss_combined/snpEFF_AFR_vars_final_0.0001_annotation.txt", col.names = T, row.names = F,  quote = F, sep = "\t")

table(AF_AFR.0.0001.final$`ANN[*].GENE`)
# BAG3   DSP  LMNA  MYH7 SCN5A  TCAP TNNT2   TTN 
# 23     1    26    48     9     3    13    20 




# with.carriers <- AF_AFR.0.0001.final[AF_AFR.0.0001.final$N_carriers_in_AFR >= 1,]
# table(with.carriers$`ANN[*].GENE`)
# # DSP  LMNA  MYH7 SCN5A  TCAP   TTN 
# # 4     4     1     2     1    42
# View(as.data.frame(with.carriers$SNP[with.carriers$`ANN[*].GENE` == "DSP"]))
# View(as.data.frame(with.carriers$SNP[with.carriers$`ANN[*].GENE` == "LMNA"]))
# View(as.data.frame(with.carriers$SNP[with.carriers$`ANN[*].GENE` == "MYH7"]))
# View(as.data.frame(with.carriers$SNP[with.carriers$`ANN[*].GENE` == "SCN5A"]))
# View(as.data.frame(with.carriers$SNP[with.carriers$`ANN[*].GENE` == "TCAP"]))
# View(as.data.frame(with.carriers$SNP[with.carriers$`ANN[*].GENE` == "TTN"]))
# View(as.data.frame(with.carriers$SNP[with.carriers$`ANN[*].GENE` == "TTN" & with.carriers$TTN_PSI=="Yes"]))


# extract variants by gene categories
genes <- names(table(AF_AFR.0.0001.final$`ANN[*].GENE`))
genes <- genes[!genes %in% "TTN"]
for (i in 1:length(genes)){
  gene <- AF_AFR.0.0001.final$SNP[AF_AFR.0.0001.final$`ANN[*].GENE` == genes[i]]
  write.table(as.data.frame(gene), paste0("Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3/rare_variant_sjlife_ccss_combined/", genes[i], "_AFR_vars.txt"), col.names = F, row.names = F,  quote = F, sep = "\t")
}

TTN_PSI <- AF_AFR.0.0001.final$SNP[AF_AFR.0.0001.final$TTN_PSI == "Yes"]
length(TTN_PSI)
# 10
write.table(as.data.frame(TTN_PSI), "Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3/rare_variant_sjlife_ccss_combined/TTN_PSI_AFR_vars.txt", col.names = F, row.names = F,  quote = F, sep = "\t")

TTN_PSI_A_Band <- AF_AFR.0.0001.final$SNP[AF_AFR.0.0001.final$TTN_PSI_A_Band == "Yes"]
length(TTN_PSI_A_Band)
# 2
write.table(as.data.frame(TTN_PSI_A_Band), "Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3/rare_variant_sjlife_ccss_combined/TTN_PSI_A_Band_AFR_vars.txt", col.names = F, row.names = F,  quote = F, sep = "\t")


BAG3.TTN_PSI <- AF_AFR.0.0001.final$SNP[AF_AFR.0.0001.final$`ANN[*].GENE` == "BAG3" | AF_AFR.0.0001.final$TTN_PSI == "Yes"]
length(BAG3.TTN_PSI)
# 33
write.table(as.data.frame(BAG3.TTN_PSI), "Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3/rare_variant_sjlife_ccss_combined/BAG3.TTN_PSI_AFR_vars.txt", col.names = F, row.names = F,  quote = F, sep = "\t")

BAG3.TTN_PSI_A_Band <- AF_AFR.0.0001.final$SNP[AF_AFR.0.0001.final$`ANN[*].GENE` == "BAG3" | AF_AFR.0.0001.final$TTN_PSI_A_Band == "Yes"]
length(BAG3.TTN_PSI_A_Band)
# 25
write.table(as.data.frame(BAG3.TTN_PSI_A_Band), "Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3/rare_variant_sjlife_ccss_combined/BAG3.TTN_PSI_A_Band_AFR_vars.txt", col.names = F, row.names = F,  quote = F, sep = "\t")

