# save.image("SNPEFF_clinvar_metaSVM_LoF_from_R_filtering_process_PreQC_VCF.RData")
load("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/MERGED_sjlife1_2_PreQC/cleaned/annotation/snpEff/SNPEFF_clinvar_metaSVM_LoF_from_R_filtering_process_PreQC_VCF.RData")

setwd("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/MERGED_sjlife1_2_PreQC/cleaned/annotation/annovar/pablo_garcia_et_al_nine_genes")
## Extract variants regrions from nine genes: BAG3, DSP, LMNA, MYH7, SCN5A, TCAP, TNNC1, TNNT2, and TTN 
# BAG3
NINE_GENES.annovar <- read.table("NINE_GENES_ANNOVAR", sep = "\t", header = T)
dim(NINE_GENES.annovar)
# 21515

NINE_GENES.annovar$gnomAD_genome_ALL <- as.numeric(NINE_GENES.annovar$gnomAD_genome_ALL)
NINE_GENES.annovar$gnomAD_genome_NFE <- as.numeric(NINE_GENES.annovar$gnomAD_genome_NFE)




NINE_GENES.annovar <- NINE_GENES.annovar %>%
  filter(!(is.na(gnomAD_genome_ALL) & is.na(gnomAD_genome_NFE)))

NINE_GENES.annovar <- NINE_GENES.annovar[which(NINE_GENES.annovar$gnomAD_genome_ALL < 0.01 & NINE_GENES.annovar$gnomAD_genome_NFE < 0.01),]
dim(NINE_GENES.annovar)
# 11407

## Part 2

## Read bim files from ccss and sjlife maf less than 0.01 and get common variants bettween the cohorts to keep.
ccss_vars_bim <- read.table("all_ccss_exp_vars_lt_maf_0.01.txt")
sjlife_vars_bim <- read.table("all_sjlife_vars_lt_maf_0.01.txt")

colnames(ccss_vars_bim) <- c("CHROM", "SNP", "centi", "POS", "REF", "ALT")
colnames(sjlife_vars_bim) <- c("CHROM", "SNP", "centi", "POS", "REF", "ALT")


## Flip alleles
sjlife_vars_bim$REF_flipped <- chartr("acgtACGT", "tgcaTGCA", sjlife_vars_bim$REF)
sjlife_vars_bim$ALT_flipped <- chartr("acgtACGT", "tgcaTGCA", sjlife_vars_bim$ALT)

## Flip alleles
ccss_vars_bim$REF_flipped <- chartr("acgtACGT", "tgcaTGCA", ccss_vars_bim$REF)
ccss_vars_bim$ALT_flipped <- chartr("acgtACGT", "tgcaTGCA", ccss_vars_bim$ALT)

# Create positional keys
ccss_vars_bim$KEY <- paste(ccss_vars_bim$CHROM, ccss_vars_bim$POS, sep = ":")
sjlife_vars_bim$KEY <- paste(sjlife_vars_bim$CHROM, sjlife_vars_bim$POS, sep = ":")

## Keep gnomad rare variants
NINE_GENES.annovar$KEY <- gsub("chr", "", paste(NINE_GENES.annovar$Otherinfo4, NINE_GENES.annovar$Otherinfo5, sep = ":"))
# create SNPID for annovar
NINE_GENES.annovar$SNP <- paste0(NINE_GENES.annovar$Otherinfo4,":", NINE_GENES.annovar$Otherinfo5,":",
                                 NINE_GENES.annovar$Otherinfo7,":", NINE_GENES.annovar$Otherinfo8)
NINE_GENES.annovar$REF <- NINE_GENES.annovar$Otherinfo7
NINE_GENES.annovar$ALT <- NINE_GENES.annovar$Otherinfo8


sjlife_vars_bim <- sjlife_vars_bim[sjlife_vars_bim$KEY %in% NINE_GENES.annovar$KEY,]
ccss_vars_bim <- ccss_vars_bim[ccss_vars_bim$KEY %in% NINE_GENES.annovar$KEY,]

dim(sjlife_vars_bim)
# 10361     9
dim(ccss_vars_bim)
# [1] 5518    7


## gnomAD in SJLIFE
nine.genes.pablo.garcia <- cbind.data.frame(SNP=NINE_GENES.annovar$SNP, ref=NINE_GENES.annovar$Ref, alt=NINE_GENES.annovar$Alt, REF=NINE_GENES.annovar$REF, ALT=NINE_GENES.annovar$ALT, KEY=NINE_GENES.annovar$KEY)
# nine.genes.pablo.garcia$MATCH <- ifelse(nine.genes.pablo.garcia$ref == nine.genes.pablo.garcia$REF, "YES", "NO")

# i=470; chr10:119663435
# 119670135

for (i in 1:nrow(sjlife_vars_bim)){
  print(paste0("Doing iteration: ", i))
  if (sum(sjlife_vars_bim$KEY[i] %in% nine.genes.pablo.garcia$KEY) > 0){ # Only if position matches; do
    match.index <- grep(paste0("chr",sjlife_vars_bim$KEY[i],":"), nine.genes.pablo.garcia$SNP)
    for(j in 1:length(match.index)){ # direct match or match by swapping alleles
      if(sjlife_vars_bim$REF[i] == nine.genes.pablo.garcia$REF[match.index[j]] & sjlife_vars_bim$ALT[i] == nine.genes.pablo.garcia$ALT[match.index[j]]|
         sjlife_vars_bim$REF[i] == nine.genes.pablo.garcia$ALT[match.index[j]] & sjlife_vars_bim$ALT[i] == nine.genes.pablo.garcia$REF[match.index[j]]){
        sjlife_vars_bim$MATCH_gnomAD[i] <- "DIRECT_MATCH"
        sjlife_vars_bim$gnomAD_equivalent[i] <- nine.genes.pablo.garcia$SNP[match.index[j]]
      } else if # match by flipping alleles
      (sjlife_vars_bim$REF_flipped[i] == nine.genes.pablo.garcia$REF[match.index[j]] & sjlife_vars_bim$ALT_flipped[i] == nine.genes.pablo.garcia$ALT[match.index[j]]){
        sjlife_vars_bim$MATCH_gnomAD[i] <- "INDIRECT_MATCH1"
        sjlife_vars_bim$gnomAD_equivalent[i] <- nine.genes.pablo.garcia$SNP[match.index[j]]
      } else if # match by swapping flipped alleles
      (sjlife_vars_bim$REF_flipped[i] == nine.genes.pablo.garcia$ALT[match.index[j]] & sjlife_vars_bim$ALT_flipped[i] == nine.genes.pablo.garcia$REF[match.index[j]]){
        sjlife_vars_bim$MATCH_gnomAD[i] <- "INDIRECT_MATCH2"
        sjlife_vars_bim$gnomAD_equivalent[i] <- nine.genes.pablo.garcia$SNP[match.index[j]]
      } else if # match by flipping one allele or by swapping one of the flipped alleles
        ((sjlife_vars_bim$REF_flipped[i] == nine.genes.pablo.garcia$REF[match.index[j]] & sjlife_vars_bim$ALT[i] == nine.genes.pablo.garcia$ALT[match.index[j]])|
        (sjlife_vars_bim$ALT_flipped[i] == nine.genes.pablo.garcia$REF[match.index[j]] & sjlife_vars_bim$REF[i] == nine.genes.pablo.garcia$ALT[match.index[j]])){
        sjlife_vars_bim$MATCH_gnomAD[i] <- "INDIRECT_MATCH3"
        sjlife_vars_bim$gnomAD_equivalent[i] <- nine.genes.pablo.garcia$SNP[match.index[j]]
      } else{
        sjlife_vars_bim$MATCH_gnomAD[i] <- NA
        sjlife_vars_bim$gnomAD_equivalent[i] <- NA
      }
      
    }
  } else {
    sjlife_vars_bim$MATCH_gnomAD[i] <- NA
    sjlife_vars_bim$gnomAD_equivalent[i] <- NA
  }
  
}

# Keep gnomad rare variants only
sjlife_vars_bim <- sjlife_vars_bim[!is.na(sjlife_vars_bim$MATCH_gnomAD),]


## gnomAD in CCSS_Exp

for (i in 1:nrow(ccss_vars_bim)){
  print(paste0("Doing iteration: ", i))
  if (sum(ccss_vars_bim$KEY[i] %in% nine.genes.pablo.garcia$KEY) > 0){ # Only if position matches; do
    match.index <- grep(paste0("chr",ccss_vars_bim$KEY[i],":"), nine.genes.pablo.garcia$SNP)
    for(j in 1:length(match.index)){ # direct match or match by swapping alleles
      if(ccss_vars_bim$REF[i] == nine.genes.pablo.garcia$REF[match.index[j]] & ccss_vars_bim$ALT[i] == nine.genes.pablo.garcia$ALT[match.index[j]]|
         ccss_vars_bim$REF[i] == nine.genes.pablo.garcia$ALT[match.index[j]] & ccss_vars_bim$ALT[i] == nine.genes.pablo.garcia$REF[match.index[j]]){
        ccss_vars_bim$MATCH_gnomAD[i] <- "DIRECT_MATCH"
        ccss_vars_bim$gnomAD_equivalent[i] <- nine.genes.pablo.garcia$SNP[match.index[j]]
      } else if # match by flipping alleles
      (ccss_vars_bim$REF_flipped[i] == nine.genes.pablo.garcia$REF[match.index[j]] & ccss_vars_bim$ALT_flipped[i] == nine.genes.pablo.garcia$ALT[match.index[j]]){
        ccss_vars_bim$MATCH_gnomAD[i] <- "INDIRECT_MATCH1"
        ccss_vars_bim$gnomAD_equivalent[i] <- nine.genes.pablo.garcia$SNP[match.index[j]]
      } else if # match by swapping flipped alleles
      (ccss_vars_bim$REF_flipped[i] == nine.genes.pablo.garcia$ALT[match.index[j]] & ccss_vars_bim$ALT_flipped[i] == nine.genes.pablo.garcia$REF[match.index[j]]){
        ccss_vars_bim$MATCH_gnomAD[i] <- "INDIRECT_MATCH2"
        ccss_vars_bim$gnomAD_equivalent[i] <- nine.genes.pablo.garcia$SNP[match.index[j]]
      } else if # match by flipping one allele or by swapping one of the flipped alleles
      ((ccss_vars_bim$REF_flipped[i] == nine.genes.pablo.garcia$REF[match.index[j]] & ccss_vars_bim$ALT[i] == nine.genes.pablo.garcia$ALT[match.index[j]])|
       (ccss_vars_bim$ALT_flipped[i] == nine.genes.pablo.garcia$REF[match.index[j]] & ccss_vars_bim$REF[i] == nine.genes.pablo.garcia$ALT[match.index[j]])){
        ccss_vars_bim$MATCH_gnomAD[i] <- "INDIRECT_MATCH3"
        ccss_vars_bim$gnomAD_equivalent[i] <- nine.genes.pablo.garcia$SNP[match.index[j]]
      } else{
        ccss_vars_bim$MATCH_gnomAD[i] <- NA
        ccss_vars_bim$gnomAD_equivalent[i] <- NA
      }
      
    }
  } else {
    ccss_vars_bim$MATCH_gnomAD[i] <- NA
    ccss_vars_bim$gnomAD_equivalent[i] <- NA
  }
  
}


## Keep gnomad rare variants only
ccss_vars_bim <- ccss_vars_bim[!is.na(ccss_vars_bim$MATCH_gnomAD),]

## Add gnomAD maf
sjlife_vars_bim$gnomAD_genome_ALL <- NINE_GENES.annovar$gnomAD_genome_ALL[match(sjlife_vars_bim$gnomAD_equivalent , NINE_GENES.annovar$SNP)]
sjlife_vars_bim$gnomAD_genome_NFE <- NINE_GENES.annovar$gnomAD_genome_NFE[match(sjlife_vars_bim$gnomAD_equivalent , NINE_GENES.annovar$SNP)]

ccss_vars_bim$gnomAD_genome_ALL <- NINE_GENES.annovar$gnomAD_genome_ALL[match(ccss_vars_bim$gnomAD_equivalent , NINE_GENES.annovar$SNP)]
ccss_vars_bim$gnomAD_genome_NFE <- NINE_GENES.annovar$gnomAD_genome_NFE[match(ccss_vars_bim$gnomAD_equivalent , NINE_GENES.annovar$SNP)]




## Short list first by matching LOF and CLINVAR positions
CLINVAR.unique$KEY1 <- paste0(CLINVAR.unique$CHROM, ":", CLINVAR.unique$POS)
LoF.unique$KEY1 <- paste0(LoF.unique$CHROM, ":", LoF.unique$POS)

sjlife.CLINVAR.by.pos <- sjlife_vars_bim[sjlife_vars_bim$KEY %in% gsub("chr", "", CLINVAR.unique$KEY1),]
dim(sjlife.CLINVAR.by.pos)
# [1]  5 13

sjlife.LoF.by.pos <- sjlife_vars_bim[sjlife_vars_bim$KEY %in% gsub("chr", "", LoF.unique$KEY1),]
dim(sjlife.LoF.by.pos)
# [1] 204  17

sjlife.CLINVAR.by.pos$PRED_TYPE <- "CLINVAR"
sjlife.LoF.by.pos$PRED_TYPE <- "LOF"
sjlife_vars_bim <- rbind.data.frame(sjlife.CLINVAR.by.pos, sjlife.LoF.by.pos)
# remove LoF duplicated keeping Clinvar as is
sjlife_vars_bim <- sjlife_vars_bim[!(duplicated(sjlife_vars_bim$SNP) & sjlife_vars_bim$PRED_TYPE == "LOF"),]


CCSS.CLINVAR.by.pos <- ccss_vars_bim[ccss_vars_bim$KEY %in% gsub("chr", "", CLINVAR.unique$KEY1),]
dim(CCSS.CLINVAR.by.pos)
# [1]  1 13
CCSS.LoF.by.pos <- ccss_vars_bim[ccss_vars_bim$KEY %in% gsub("chr", "", LoF.unique$KEY1),]
dim(CCSS.LoF.by.pos)
# [1] 106  17

CCSS.CLINVAR.by.pos$PRED_TYPE <- "CLINVAR"
CCSS.LoF.by.pos$PRED_TYPE <- "LOF"
ccss_vars_bim <- rbind.data.frame(CCSS.CLINVAR.by.pos, CCSS.LoF.by.pos)
# remove LoF duplicated keeping Clinvar as is
ccss_vars_bim <- ccss_vars_bim[!(duplicated(ccss_vars_bim$SNP) & ccss_vars_bim$PRED_TYPE == "LOF"),]



## Match with P/LP by SNP ID
# SJLIFE
sjlife_vars_bim$SNP1 <- paste0("chr",sjlife_vars_bim$CHROM,":",sjlife_vars_bim$POS, ":", sjlife_vars_bim$REF, ":", sjlife_vars_bim$ALT)
sjlife_vars_bim$SNP2 <- paste0("chr",sjlife_vars_bim$CHROM,":",sjlife_vars_bim$POS, ":", sjlife_vars_bim$ALT, ":", sjlife_vars_bim$REF)
sum(sjlife_vars_bim$SNP1 %in% LoF.unique$KEY)
# 0
sum(sjlife_vars_bim$SNP2 %in% LoF.unique$KEY)
# 204
# sjlife_vars_bim$SNP2[!sjlife_vars_bim$SNP2 %in% LoF.unique$KEY] # not in LOF
sum(sjlife_vars_bim$SNP1 %in% CLINVAR.unique$KEY)
# 0
sum(sjlife_vars_bim$SNP2 %in% CLINVAR.unique$KEY)
# 5
sum(sjlife_vars_bim$SNP2 %in% CLINVAR.unique$KEY)
# sjlife_vars_bim$SNP2[sjlife_vars_bim$SNP2 %in% CLINVAR.unique$KEY] # Which is present in CLINVAR



# CCSS
ccss_vars_bim$SNP1 <- paste0("chr",ccss_vars_bim$CHROM,":",ccss_vars_bim$POS, ":", ccss_vars_bim$REF, ":", ccss_vars_bim$ALT)
ccss_vars_bim$SNP2 <- paste0("chr",ccss_vars_bim$CHROM,":",ccss_vars_bim$POS, ":", ccss_vars_bim$ALT, ":", ccss_vars_bim$REF)
sum(ccss_vars_bim$SNP1 %in% LoF.unique$KEY)
# 0
sum(ccss_vars_bim$SNP2 %in% LoF.unique$KEY)
# 106

sum(ccss_vars_bim$SNP1 %in% CLINVAR.unique$KEY)
# 0
sum(ccss_vars_bim$SNP2 %in% CLINVAR.unique$KEY)
# 1



## Add ANN_EFFECT
# SJLIFE
sjlife_vars_bim$ANN_EFFECT <- LoF.unique$`ANN[*].EFFECT`[match(sjlife_vars_bim$SNP2, LoF.unique$KEY)]
sum(sjlife_vars_bim$SNP %in% CLINVAR.unique$KEY)
# 5
sjlife_vars_bim$ANN_EFFECT[is.na(sjlife_vars_bim$ANN_EFFECT)] <- CLINVAR.unique$`ANN[*].EFFECT`[match(sjlife_vars_bim$SNP2, CLINVAR.unique$KEY)][is.na(sjlife_vars_bim$ANN_EFFECT)]

# CCSS
ccss_vars_bim$ANN_EFFECT <- LoF.unique$`ANN[*].EFFECT`[match(ccss_vars_bim$SNP2, LoF.unique$KEY)]
sum(ccss_vars_bim$SNP %in% CLINVAR.unique$KEY)
# 1
ccss_vars_bim$ANN_EFFECT[is.na(ccss_vars_bim$ANN_EFFECT)] <- CLINVAR.unique$`ANN[*].EFFECT`[match(ccss_vars_bim$SNP2, CLINVAR.unique$KEY)][is.na(ccss_vars_bim$ANN_EFFECT)]





sjlife_vars_bim <- sjlife_vars_bim[!sjlife_vars_bim$ANN_EFFECT %in% c("splice_region_variant&intron_variant",
                                                                      "splice_region_variant&non_coding_transcript_exon_variant",
                                                                      "splice_region_variant", 
                                                                      "splice_region_variant&synonymous_variant"),]

ccss_vars_bim <- ccss_vars_bim[!ccss_vars_bim$ANN_EFFECT %in% c("splice_region_variant&intron_variant",
                                                                      "splice_region_variant&non_coding_transcript_exon_variant",
                                                                      "splice_region_variant", 
                                                                      "splice_region_variant&synonymous_variant"),]



# Since all SNP ids match we can now use this dataframe for logistic regression analysis on nine genes


####################
## Annotate genes ##
####################

## SJLIFE 

# BAG3
sjlife_vars_bim$GENE[sjlife_vars_bim$CHROM == 10 &  (sjlife_vars_bim$POS >= 119651380 & sjlife_vars_bim$POS <= 119677819)] <- "BAG3"
# DSP
sjlife_vars_bim$GENE[sjlife_vars_bim$CHROM == 6 &  (sjlife_vars_bim$POS >= 7541671 & sjlife_vars_bim$POS <= 7586714)] <- "DSP"
# LMNA
sjlife_vars_bim$GENE[sjlife_vars_bim$CHROM == 1 &  (sjlife_vars_bim$POS >= 156082573 & sjlife_vars_bim$POS <= 156140081)] <- "LMNA"
# MYH7
sjlife_vars_bim$GENE[sjlife_vars_bim$CHROM == 14 &  (sjlife_vars_bim$POS >= 23412740 & sjlife_vars_bim$POS <= 23435660)] <- "MYH7"
# SCN5A
sjlife_vars_bim$GENE[sjlife_vars_bim$CHROM == 3 &  (sjlife_vars_bim$POS >= 38548062 & sjlife_vars_bim$POS <= 38649687)] <- "SCN5A"
# TCAP
sjlife_vars_bim$GENE[sjlife_vars_bim$CHROM == 17 &  (sjlife_vars_bim$POS >= 39665349 & sjlife_vars_bim$POS <= 39666554)] <- "TCAP"
# TNNC1
sjlife_vars_bim$GENE[sjlife_vars_bim$CHROM == 3 &  (sjlife_vars_bim$POS >= 52451100 & sjlife_vars_bim$POS <= 52454041)] <- "TNNC1"
# TNNT2
sjlife_vars_bim$GENE[sjlife_vars_bim$CHROM == 1 &  (sjlife_vars_bim$POS >= 201359014 & sjlife_vars_bim$POS <= 201377680)] <- "TNNT2"
# TTN
sjlife_vars_bim$GENE[sjlife_vars_bim$CHROM == 2 &  (sjlife_vars_bim$POS >= 178525989 & sjlife_vars_bim$POS <= 178807423)] <- "TTN"

table(sjlife_vars_bim$GENE)
# BAG3   DSP  LMNA  MYH7 SCN5A TNNT2   TTN 
# 2     2     4     1     2     3    46
## CCSS

ccss_vars_bim$GENE[ccss_vars_bim$CHROM == 10 &  (ccss_vars_bim$POS >= 119651380 & ccss_vars_bim$POS <= 119677819)] <- "BAG3"
# DSP
ccss_vars_bim$GENE[ccss_vars_bim$CHROM == 6 &  (ccss_vars_bim$POS >= 7541671 & ccss_vars_bim$POS <= 7586714)] <- "DSP"
# LMNA
ccss_vars_bim$GENE[ccss_vars_bim$CHROM == 1 &  (ccss_vars_bim$POS >= 156082573 & ccss_vars_bim$POS <= 156140081)] <- "LMNA"
# MYH7
ccss_vars_bim$GENE[ccss_vars_bim$CHROM == 14 &  (ccss_vars_bim$POS >= 23412740 & ccss_vars_bim$POS <= 23435660)] <- "MYH7"
# SCN5A
ccss_vars_bim$GENE[ccss_vars_bim$CHROM == 3 &  (ccss_vars_bim$POS >= 38548062 & ccss_vars_bim$POS <= 38649687)] <- "SCN5A"
# TCAP
ccss_vars_bim$GENE[ccss_vars_bim$CHROM == 17 &  (ccss_vars_bim$POS >= 39665349 & ccss_vars_bim$POS <= 39666554)] <- "TCAP"
# TNNC1
ccss_vars_bim$GENE[ccss_vars_bim$CHROM == 3 &  (ccss_vars_bim$POS >= 52451100 & ccss_vars_bim$POS <= 52454041)] <- "TNNC1"
# TNNT2
ccss_vars_bim$GENE[ccss_vars_bim$CHROM == 1 &  (ccss_vars_bim$POS >= 201359014 & ccss_vars_bim$POS <= 201377680)] <- "TNNT2"
# TTN
ccss_vars_bim$GENE[ccss_vars_bim$CHROM == 2 &  (ccss_vars_bim$POS >= 178525989 & ccss_vars_bim$POS <= 178807423)] <- "TTN"

table(ccss_vars_bim$GENE)
# BAG3  LMNA SCN5A TNNT2   TTN 
# 1     3     1     3    24 

save.image("allele_counts_of_rare_variants.RDATA")




