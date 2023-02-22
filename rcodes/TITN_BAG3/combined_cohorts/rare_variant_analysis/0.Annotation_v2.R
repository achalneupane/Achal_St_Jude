# save.image("SNPEFF_clinvar_metaSVM_LoF_from_R_filtering_process_PreQC_VCF.RData")
load("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/MERGED_sjlife1_2_PreQC/cleaned/annotation/snpEff/SNPEFF_clinvar_metaSVM_LoF_from_R_filtering_process_PreQC_VCF.RData")

setwd("Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3/ALL_P_LP_combinations/rare_variant_analysis_v2/pablo_garcia_et_al_nine_genes")
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

## Read bim files from ccss_exp and sjlife maf less than 0.01.
sjlife_ccss_exp_vars_bim <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3/rare_variant_sjlife_ccss_combined/sjlife_ccss_exp_merged_maxmaf_0.01.bim")

colnames(sjlife_ccss_exp_vars_bim) <- c("CHROM", "SNP", "centi", "POS", "REF", "ALT")
# colnames(sjlife_vars_bim) <- c("CHROM", "SNP", "centi", "POS", "REF", "ALT")


## Flip alleles
sjlife_ccss_exp_vars_bim$REF_flipped <- chartr("acgtACGT", "tgcaTGCA", sjlife_ccss_exp_vars_bim$REF)
sjlife_ccss_exp_vars_bim$ALT_flipped <- chartr("acgtACGT", "tgcaTGCA", sjlife_ccss_exp_vars_bim$ALT)

# Create positional keys
sjlife_ccss_exp_vars_bim$KEY <- paste(sjlife_ccss_exp_vars_bim$CHROM, sjlife_ccss_exp_vars_bim$POS, sep = ":")

## Keep gnomad rare variants
NINE_GENES.annovar$KEY <- gsub("chr", "", paste(NINE_GENES.annovar$Otherinfo4, NINE_GENES.annovar$Otherinfo5, sep = ":"))
# create SNPID for annovar
NINE_GENES.annovar$SNP <- paste0(NINE_GENES.annovar$Otherinfo4,":", NINE_GENES.annovar$Otherinfo5,":",
                                 NINE_GENES.annovar$Otherinfo7,":", NINE_GENES.annovar$Otherinfo8)
NINE_GENES.annovar$REF <- NINE_GENES.annovar$Otherinfo7
NINE_GENES.annovar$ALT <- NINE_GENES.annovar$Otherinfo8

sjlife_vars_bim <- sjlife_ccss_exp_vars_bim

sjlife_vars_bim <- sjlife_vars_bim[sjlife_vars_bim$KEY %in% NINE_GENES.annovar$KEY,]

dim(sjlife_vars_bim)

# cc <- NINE_GENES.annovar[NINE_GENES.annovar$Start >= 178525989 & NINE_GENES.annovar$Start <= 178807423,]
# sum(cc$gnomAD_genome_ALL < 0.01 & cc$gnomAD_genome_NFE < 0.01, na.rm = T)




## gnomAD
nine.genes.pablo.garcia <- cbind.data.frame(SNP=NINE_GENES.annovar$SNP, ref=NINE_GENES.annovar$Ref, alt=NINE_GENES.annovar$Alt, REF=NINE_GENES.annovar$REF, ALT=NINE_GENES.annovar$ALT, KEY=NINE_GENES.annovar$KEY)
# nine.genes.pablo.garcia$MATCH <- ifelse(nine.genes.pablo.garcia$ref == nine.genes.pablo.garcia$REF, "YES", "NO")

# i=470; chr10:119663435
# 119670135

for (i in 1:nrow(sjlife_vars_bim)){
  print(paste0("Doing iteration: ", i))
  if (sum(sjlife_vars_bim$KEY[i] %in% nine.genes.pablo.garcia$KEY) > 0){ # Only if position matches; do
    match.index <- grep(paste0("chr",sjlife_vars_bim$KEY[i],":"), nine.genes.pablo.garcia$SNP)
    for(j in 1:length(match.index)) { # direct match or match by swapping alleles
      if(sjlife_vars_bim$REF[i] == nine.genes.pablo.garcia$REF[match.index[j]] & sjlife_vars_bim$ALT[i] == nine.genes.pablo.garcia$ALT[match.index[j]]){
        sjlife_vars_bim$MATCH_gnomAD[i] <- "DIRECT_MATCH"
        sjlife_vars_bim$gnomAD_equivalent[i] <- nine.genes.pablo.garcia$SNP[match.index[j]]
      } else if
        (sjlife_vars_bim$REF[i] == nine.genes.pablo.garcia$ALT[match.index[j]] & sjlife_vars_bim$ALT[i] == nine.genes.pablo.garcia$REF[match.index[j]]){
        sjlife_vars_bim$MATCH_gnomAD[i] <- "INDIRECT_MATCH0"
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
        # sjlife_vars_bim$MATCH_gnomAD[i] <- NA
        # sjlife_vars_bim$gnomAD_equivalent[i] <- NA
        print(i,j)
      }
      
    }
  } else {
    sjlife_vars_bim$MATCH_gnomAD[i] <- NA
    sjlife_vars_bim$gnomAD_equivalent[i] <- NA
  }
  
}






sjlife_vars_bim <- sjlife_vars_bim[!is.na(sjlife_vars_bim$MATCH_CCSS),]
sjlife_vars_bim <- sjlife_vars_bim[!is.na(sjlife_vars_bim$MATCH_gnomAD),]

sjlife_vars_bim$gnomAD_genome_ALL <- NINE_GENES.annovar$gnomAD_genome_ALL[match(sjlife_vars_bim$gnomAD_equivalent , NINE_GENES.annovar$SNP)]
sjlife_vars_bim$gnomAD_genome_NFE <- NINE_GENES.annovar$gnomAD_genome_NFE[match(sjlife_vars_bim$gnomAD_equivalent , NINE_GENES.annovar$SNP)]



## Short list first by matching LOF and CLINVAR positions
CLINVAR.unique$KEY1 <- paste0(CLINVAR.unique$CHROM, ":", CLINVAR.unique$POS)
LoF.unique$KEY1 <- paste0(LoF.unique$CHROM, ":", LoF.unique$POS)

sjlife.CLINVAR.by.pos <- sjlife_vars_bim[sjlife_vars_bim$KEY %in% gsub("chr", "", CLINVAR.unique$KEY1),]
dim(sjlife.CLINVAR.by.pos)
# [1]  1 17
sjlife.LoF.by.pos <- sjlife_vars_bim[sjlife_vars_bim$KEY %in% gsub("chr", "", LoF.unique$KEY1),]
dim(sjlife.LoF.by.pos)
# [1] 100  17

sjlife_vars_bim <- rbind.data.frame(sjlife.CLINVAR.by.pos, sjlife.LoF.by.pos)


## Match with P/LP by SNP ID
sjlife_vars_bim$SNP1 <- paste0("chr",sjlife_vars_bim$CHROM,":",sjlife_vars_bim$POS, ":", sjlife_vars_bim$REF, ":", sjlife_vars_bim$ALT)
sjlife_vars_bim$SNP2 <- paste0("chr",sjlife_vars_bim$CHROM,":",sjlife_vars_bim$POS, ":", sjlife_vars_bim$ALT, ":", sjlife_vars_bim$REF)
sum(sjlife_vars_bim$SNP1 %in% LoF.unique$KEY)
# 0
sum(sjlife_vars_bim$SNP2 %in% LoF.unique$KEY)
# 100
sjlife_vars_bim$ANN_EFFECT <- LoF.unique$`ANN[*].EFFECT`[match(sjlife_vars_bim$SNP2, LoF.unique$KEY)]
sjlife_vars_bim$PRED_TYPE <- LoF.unique$PRED_TYPE[match(sjlife_vars_bim$SNP2, LoF.unique$KEY)]

sum(sjlife_vars_bim$SNP %in% CLINVAR.unique$KEY)
# 1
sjlife_vars_bim$ANN_EFFECT[is.na(sjlife_vars_bim$ANN_EFFECT)] <- CLINVAR.unique$`ANN[*].EFFECT`[match(sjlife_vars_bim$SNP2, CLINVAR.unique$KEY)][is.na(sjlife_vars_bim$ANN_EFFECT)]
sjlife_vars_bim$PRED_TYPE[is.na(sjlife_vars_bim$PRED_TYPE)] <- CLINVAR.unique$PRED_TYPE[match(sjlife_vars_bim$SNP2, CLINVAR.unique$KEY)][is.na(sjlife_vars_bim$PRED_TYPE)]

sjlife_vars_bim <- sjlife_vars_bim[!sjlife_vars_bim$ANN_EFFECT %in% c("splice_region_variant&intron_variant",
                                                                      "splice_region_variant&non_coding_transcript_exon_variant",
                                                                      "splice_region_variant", 
                                                                      "splice_region_variant&synonymous_variant"),]

# Since all SNP ids match we can now use this dataframe for logistic regression analysis on nine genes


####################
## Annotate genes ##
####################
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

# save.image("common_p_LP_rare_variants_gnomad_all_gnomad_NFE_lt_0.01.RData")



write.table(sjlife_vars_bim$SNP, "sjlife_SNPS_maf_lt_0.01_gnomad_also_common_in_ccss.txt", quote = FALSE, col.names = FALSE, row.names = F)
write.table(sjlife_vars_bim$CCSS_equivalent, "ccss_SNPS_maf_lt_0.01_gnomad_also_common_in_sjlife.txt", quote = FALSE, col.names = FALSE, row.names = F)


## Now run shell script rare_variant_extraction.sh to extract ccss and sjlife overlapping variants



