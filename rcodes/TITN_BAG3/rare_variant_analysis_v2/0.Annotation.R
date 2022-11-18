# save.image("SNPEFF_clinvar_metaSVM_LoF_from_R_filtering_process_PreQC_VCF.RData")
load("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/MERGED_sjlife1_2_PreQC/cleaned/annotation/snpEff/SNPEFF_clinvar_metaSVM_LoF_from_R_filtering_process_PreQC_VCF.RData")

setwd("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/MERGED_sjlife1_2_PreQC/cleaned/annotation/annovar/pablo_garcia_et_al_nine_genes")
## Extract variants regrions from nine genes: BAG3, DSP, LMNA, MYH7, SCN5A, TCAP, TNNC1, TNNT2, and TTN 
# BAG3
NINE_GENES.annovar <- read.table("NINE_GENES", sep = "\t", header = T)
dim(NINE_GENES.annovar)
NINE_GENES.annovar$gnomAD_genome_ALL <- as.numeric(NINE_GENES.annovar$gnomAD_genome_ALL)
NINE_GENES.annovar$gnomAD_genome_NFE <- as.numeric(NINE_GENES.annovar$gnomAD_genome_NFE)
NINE_GENES.annovar <- NINE_GENES.annovar[which(NINE_GENES.annovar$gnomAD_genome_ALL < 0.01 & NINE_GENES.annovar$gnomAD_genome_NFE < 0.01),]
dim(NINE_GENES.annovar)

## Part 2

## Read bim files from ccss and sjlife maf less than 0.01 and get common variants bettween the cohorts to keep.
ccss_vars_bim <- read.table("all_ccss_exp_vars_lt_maf_0.01.txt")
sjlife_vars_bim <- read.table("all_sjlife_vars_lt_maf_0.01.txt")

colnames(ccss_vars_bim) <- c("CHROM", "SNP", "centi", "POS", "REF", "ALT")
colnames(sjlife_vars_bim) <- c("CHROM", "SNP", "centi", "POS", "REF", "ALT")


## Flip alleles
sjlife_vars_bim$REF_flipped <- chartr("acgtACGT", "tgcaTGCA", sjlife_vars_bim$REF)
sjlife_vars_bim$ALT_flipped <- chartr("acgtACGT", "tgcaTGCA", sjlife_vars_bim$ALT)

# Create positional keys
ccss_vars_bim$KEY <- paste(ccss_vars_bim$CHROM, ccss_vars_bim$POS, sep = ":")
sjlife_vars_bim$KEY <- paste(sjlife_vars_bim$CHROM, sjlife_vars_bim$POS, sep = ":")

## Keep gnomad rare variants
NINE_GENES.annovar$KEY <- gsub("chr", "", paste(NINE_GENES.annovar$Otherinfo4, NINE_GENES.annovar$Otherinfo5, sep = ":"))
# create SNPID for annovar
NINE_GENES.annovar$SNP <- paste0(NINE_GENES.annovar$Otherinfo4,":", NINE_GENES.annovar$Otherinfo5,":",
                                 NINE_GENES.annovar$Otherinfo7,":", NINE_GENES.annovar$Otherinfo8)

sjlife_vars_bim <- sjlife_vars_bim[sjlife_vars_bim$KEY %in% NINE_GENES.annovar$KEY,]

dim(sjlife_vars_bim)
# 5263

for (i in 1:nrow(sjlife_vars_bim)){
  print(paste0("Doing iteration: ", i))
  if (sum(grepl(paste0("chr",sjlife_vars_bim$KEY[i],":"), ccss_vars_bim$SNP)) > 0){ # Only if position matches; do
    match.index <- grep(paste0("chr",sjlife_vars_bim$KEY[i],":"), ccss_vars_bim$SNP)
    for(j in 1:length(match.index)){
    if(sjlife_vars_bim$REF[i] == ccss_vars_bim$REF[match.index[j]] & sjlife_vars_bim$ALT[i] == ccss_vars_bim$ALT[match.index[j]]){
      sjlife_vars_bim$MATCH_CCSS[i] <- "DIRECT_MATCH"
      sjlife_vars_bim$CCSS_equivalent[i] <- ccss_vars_bim$SNP[match.index[j]]
    } else if
      ((sjlife_vars_bim$REF_flipped[i] == ccss_vars_bim$REF[match.index[j]] & sjlife_vars_bim$ALT_flipped[i] == ccss_vars_bim$ALT[match.index[j]])|
      (sjlife_vars_bim$REF_flipped[i] == ccss_vars_bim$REF[match.index[j]] & sjlife_vars_bim$ALT[i] == ccss_vars_bim$ALT[match.index[j]])|
      (sjlife_vars_bim$REF[i] == ccss_vars_bim$REF[j] & sjlife_vars_bim$ALT_flipped[i] == ccss_vars_bim$ALT[match.index[j]])){
      sjlife_vars_bim$MATCH_CCSS[i] <- "INDIRECT_MATCH"
      sjlife_vars_bim$CCSS_equivalent[i] <- ccss_vars_bim$SNP[match.index[j]]
    }else{
      sjlife_vars_bim$MATCH_CCSS[i] <- NA
      sjlife_vars_bim$CCSS_equivalent[i] <- NA
    }
      
    }
  } else {
    sjlife_vars_bim$MATCH_CCSS[i] <- NA
    sjlife_vars_bim$CCSS_equivalent[i] <- NA
  }
    
}


## Now get the match index from gnomAD variants with less than 0.01 maf
NINE_GENES.annovar_vars <- NINE_GENES.annovar$Otherinfo10
for (i in 1:nrow(sjlife_vars_bim)){
  print(paste0("Doing iteration: ", i))
  if (sum(grepl(paste0("chr",sjlife_vars_bim$KEY[i],":"), NINE_GENES.annovar_vars$SNP)) > 0){ # Only if position matches; do
    match.index <- grep(paste0("chr",sjlife_vars_bim$KEY[i],":"), NINE_GENES.annovar_vars$SNP)
    for(j in 1:length(match.index)){
      if(sjlife_vars_bim$REF[i] == NINE_GENES.annovar_vars$REF[match.index[j]] & sjlife_vars_bim$ALT[i] == NINE_GENES.annovar_vars$ALT[match.index[j]]){
        sjlife_vars_bim$MATCH_gnomad[i] <- "DIRECT_MATCH"
        sjlife_vars_bim$gnomad_equivalent[i] <- NINE_GENES.annovar_vars$SNP[match.index[j]]
      } else if
      ((sjlife_vars_bim$REF_flipped[i] == NINE_GENES.annovar_vars$REF[match.index[j]] & sjlife_vars_bim$ALT_flipped[i] == NINE_GENES.annovar_vars$ALT[match.index[j]])|
       (sjlife_vars_bim$REF_flipped[i] == NINE_GENES.annovar_vars$REF[match.index[j]] & sjlife_vars_bim$ALT[i] == NINE_GENES.annovar_vars$ALT[match.index[j]])|
       (sjlife_vars_bim$REF[i] == NINE_GENES.annovar_vars$REF[j] & sjlife_vars_bim$ALT_flipped[i] == NINE_GENES.annovar_vars$ALT[match.index[j]])){
        sjlife_vars_bim$MATCH_gnomad[i] <- "INDIRECT_MATCH"
        sjlife_vars_bim$gnomad_equivalent[i] <- NINE_GENES.annovar_vars$SNP[match.index[j]]
      }else{
        sjlife_vars_bim$MATCH_gnomad[i] <- NA
        sjlife_vars_bim$gnomad_equivalent[i] <- NA
      }
      
    }
  } else {
    sjlife_vars_bim$MATCH_gnomad[i] <- NA
    sjlife_vars_bim$gnomad_equivalent[i] <- NA
  }
  
}


############################
## Clinvar, LoF and Revel ##
############################

TITN_annovar <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/MERGED_sjlife1_2_PreQC/cleaned/annotation/annovar/TITN_BAG3_region_all/chr2_178_list_ALL.txt", sep = "\t", header = T)
dim(TITN_annovar)
colnames(TITN_annovar)
TITN_annovar <- TITN_annovar[TITN_annovar$Start >= 178525989 & TITN_annovar$Start <= 178807423,]
sum(TITN_annovar$gnomAD_genome_ALL < 0.01)
# 6398
TITN_annovar$gnomAD_genome_ALL <-  as.numeric(TITN_annovar$gnomAD_genome_ALL)
TITN_annovar$gnomAD_genome_NFE <-  as.numeric(TITN_annovar$gnomAD_genome_NFE)

saved.TITN_annovar <- TITN_annovar

sum(TITN_annovar$gnomAD_genome_ALL < 0.01 & TITN_annovar$gnomAD_genome_NFE < 0.01, na.rm = T)
# 5010

TITN_annovar <- TITN_annovar %>%
  filter(!(is.na(gnomAD_genome_ALL) & is.na(gnomAD_genome_NFE)))

TITN_annovar.lt.1.perc.maf.gnom.AD <- TITN_annovar[TITN_annovar$gnomAD_genome_ALL < 0.01 & TITN_annovar$gnomAD_genome_NFE < 0.01,]

TITN_annovar.lt.1.perc.maf.gnom.AD$KEY <- paste0(TITN_annovar.lt.1.perc.maf.gnom.AD$Otherinfo4,":", TITN_annovar.lt.1.perc.maf.gnom.AD$Otherinfo5,":",
                                                 TITN_annovar.lt.1.perc.maf.gnom.AD$Otherinfo7,":", TITN_annovar.lt.1.perc.maf.gnom.AD$Otherinfo8)

table(TITN_annovar.lt.1.perc.maf.gnom.AD$ExonicFunc.refGene)
# .     frameshift deletion    frameshift insertion  nonframeshift deletion nonframeshift insertion       nonsynonymous SNV 
# 3103                       2                       1                      18                       6                    1337 
# stopgain          synonymous SNV 
# 8                     535

TITN.1.per.maf.gnomad.ALL.and.NFE <- BAG3_TITN[BAG3_TITN$KEY %in% TITN_annovar.lt.1.perc.maf.gnom.AD$KEY,]
TITN.1.per.maf.gnomad.ALL.and.NFE$gnomAD_genome_ALL <- TITN_annovar.lt.1.perc.maf.gnom.AD$gnomAD_genome_ALL[match(TITN.1.per.maf.gnomad.ALL.and.NFE$KEY, TITN_annovar.lt.1.perc.maf.gnom.AD$KEY)]
TITN.1.per.maf.gnomad.ALL.and.NFE$gnomAD_genome.NFE <- TITN_annovar.lt.1.perc.maf.gnom.AD$gnomAD_genome_NFE[match(TITN.1.per.maf.gnomad.ALL.and.NFE$KEY, TITN_annovar.lt.1.perc.maf.gnom.AD$KEY)]
TITN.1.per.maf.gnomad.ALL.and.NFE$REVEL <- as.numeric(TITN_annovar.lt.1.perc.maf.gnom.AD$REVEL[match(TITN.1.per.maf.gnomad.ALL.and.NFE$KEY, TITN_annovar.lt.1.perc.maf.gnom.AD$KEY)])
TITN.1.per.maf.gnomad.ALL.and.NFE$REVEL[TITN.1.per.maf.gnomad.ALL.and.NFE$REVEL <= 0.5] <- NA

## Add Annovar annotation to this df
TITN.1.per.maf.gnomad.ALL.and.NFE$Annovar_ExonicFunc.refGene <- TITN_annovar.lt.1.perc.maf.gnom.AD$ExonicFunc.refGene[match(TITN.1.per.maf.gnomad.ALL.and.NFE$KEY, TITN_annovar.lt.1.perc.maf.gnom.AD$KEY)]

CLINVAR.Keys <- TITN.1.per.maf.gnomad.ALL.and.NFE$KEY[grepl("Clinvar", TITN.1.per.maf.gnomad.ALL.and.NFE$PRED_TYPE)]
TITN.1.per.maf.gnomad.ALL.and.NFE$CLINVAR <- ifelse(TITN.1.per.maf.gnomad.ALL.and.NFE$KEY %in% CLINVAR.Keys,"Y", "N")

LoF.Keys <- TITN.1.per.maf.gnomad.ALL.and.NFE$KEY[grepl("LoF", TITN.1.per.maf.gnomad.ALL.and.NFE$PRED_TYPE)]
TITN.1.per.maf.gnomad.ALL.and.NFE$LoF <- ifelse(TITN.1.per.maf.gnomad.ALL.and.NFE$KEY %in% LoF.Keys,"Y", "N")

MetaSVM.Keys <- TITN.1.per.maf.gnomad.ALL.and.NFE$KEY[grepl("MetaSVM", TITN.1.per.maf.gnomad.ALL.and.NFE$PRED_TYPE)]
TITN.1.per.maf.gnomad.ALL.and.NFE$MetaSVM <- ifelse(TITN.1.per.maf.gnomad.ALL.and.NFE$KEY %in% MetaSVM.Keys,"Y", "N")

TITN.1.per.maf.gnomad.ALL.and.NFE <- as.data.frame(TITN.1.per.maf.gnomad.ALL.and.NFE)
TITN.1.per.maf.gnomad.ALL.and.NFE <- TITN.1.per.maf.gnomad.ALL.and.NFE[!grepl("PRED_TYPE", colnames(TITN.1.per.maf.gnomad.ALL.and.NFE))]

# Remove duplicate rows
TITN.1.per.maf.gnomad.ALL.and.NFE <- TITN.1.per.maf.gnomad.ALL.and.NFE %>% 
  distinct(.keep_all = TRUE)

# TITN.1.per.maf.gnomad.ALL.and.NFE$KEY[duplicated(TITN.1.per.maf.gnomad.ALL.and.NFE$KEY)]
TITN.1.per.maf.gnomad.ALL.and.NFE$REVEL_TYPE <- ifelse(TITN.1.per.maf.gnomad.ALL.and.NFE$REVEL > 0.05, "REVEL_Pathogenic", "Not_pathogenic")

write.table(TITN.1.per.maf.gnomad.ALL.and.NFE, "Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/MERGED_sjlife1_2_PreQC/cleaned/annotation/snpEff/TTN_BAG3/TITN.1.per.maf.gnomad.ALL.and.NFE.txt", quote = F, col.names = T, sep = "\t", row.names = F)


TITN.1.per.maf.gnomad.ALL.and.NFE.clinvar.lof <- TITN.1.per.maf.gnomad.ALL.and.NFE$KEY[TITN.1.per.maf.gnomad.ALL.and.NFE$CLINVAR == "Y"| TITN.1.per.maf.gnomad.ALL.and.NFE$LoF =="Y" ]
write.table(TITN.1.per.maf.gnomad.ALL.and.NFE.clinvar.lof, "Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/MERGED_sjlife1_2_PreQC/cleaned/annotation/snpEff/TTN_BAG3/TITN.1.per.maf.gnomad.ALL.and.NFE.clinvar.lof.txt", quote = F, col.names = F, sep = "\t", row.names = F)

TITN.1.per.maf.gnomad.ALL.and.NFE.clinvar.lof.REVEL <- TITN.1.per.maf.gnomad.ALL.and.NFE$KEY[TITN.1.per.maf.gnomad.ALL.and.NFE$CLINVAR == "Y"| TITN.1.per.maf.gnomad.ALL.and.NFE$LoF =="Y"|TITN.1.per.maf.gnomad.ALL.and.NFE$REVEL_TYPE == "REVEL_Pathogenic" ]
write.table(TITN.1.per.maf.gnomad.ALL.and.NFE.clinvar.lof.REVEL, "Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/MERGED_sjlife1_2_PreQC/cleaned/annotation/snpEff/TTN_BAG3/TITN.1.per.maf.gnomad.ALL.and.NFE.clinvar.lof.REVEL.txt", quote = F, col.names = F, sep = "\t", row.names = F)


BAG3_annovar <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/MERGED_sjlife1_2_PreQC/cleaned/annotation/annovar/TITN_BAG3_region_all/chr10_119_list_ALL.txt", sep = "\t", header = T)
dim(BAG3_annovar)
colnames(BAG3_annovar)
BAG3_annovar <- BAG3_annovar[BAG3_annovar$Start >= 119651380 & BAG3_annovar$Start <= 119677819,]
sum(BAG3_annovar$gnomAD_genome_ALL < 0.01)
# 1063
BAG3_annovar$gnomAD_genome_ALL <-  as.numeric(BAG3_annovar$gnomAD_genome_ALL)
BAG3_annovar$gnomAD_genome_NFE <-  as.numeric(BAG3_annovar$gnomAD_genome_NFE)

saved.BAG3_annovar <- BAG3_annovar
sum(BAG3_annovar$gnomAD_genome_ALL < 0.01 & BAG3_annovar$gnomAD_genome_NFE < 0.01, na.rm = T)
# 817

BAG3_annovar <- BAG3_annovar %>%
  filter(!(is.na(gnomAD_genome_ALL) & is.na(gnomAD_genome_NFE)))

BAG3_annovar.lt.1.perc.maf.gnom.AD <- BAG3_annovar[BAG3_annovar$gnomAD_genome_ALL < 0.01 & BAG3_annovar$gnomAD_genome_NFE < 0.01,]

BAG3_annovar.lt.1.perc.maf.gnom.AD$KEY <- paste0(BAG3_annovar.lt.1.perc.maf.gnom.AD$Otherinfo4,":", BAG3_annovar.lt.1.perc.maf.gnom.AD$Otherinfo5,":",
                                                 BAG3_annovar.lt.1.perc.maf.gnom.AD$Otherinfo7,":", BAG3_annovar.lt.1.perc.maf.gnom.AD$Otherinfo8)

table(BAG3_annovar.lt.1.perc.maf.gnom.AD$ExonicFunc.refGene)
# .  nonframeshift deletion nonframeshift insertion       nonsynonymous SNV          synonymous SNV 
# 764                       1                       1                      32                      19 

BAG3.1.per.maf.gnomad.ALL.and.NFE <- BAG3_TITN[BAG3_TITN$KEY %in% BAG3_annovar.lt.1.perc.maf.gnom.AD$KEY,]
BAG3.1.per.maf.gnomad.ALL.and.NFE$gnomAD_genome_ALL <- BAG3_annovar.lt.1.perc.maf.gnom.AD$gnomAD_genome_ALL[match(BAG3.1.per.maf.gnomad.ALL.and.NFE$KEY, BAG3_annovar.lt.1.perc.maf.gnom.AD$KEY)]
BAG3.1.per.maf.gnomad.ALL.and.NFE$gnomAD_genome.NFE <- BAG3_annovar.lt.1.perc.maf.gnom.AD$gnomAD_genome_NFE[match(BAG3.1.per.maf.gnomad.ALL.and.NFE$KEY, BAG3_annovar.lt.1.perc.maf.gnom.AD$KEY)]
BAG3.1.per.maf.gnomad.ALL.and.NFE$REVEL <- as.numeric(BAG3_annovar.lt.1.perc.maf.gnom.AD$REVEL[match(BAG3.1.per.maf.gnomad.ALL.and.NFE$KEY, BAG3_annovar.lt.1.perc.maf.gnom.AD$KEY)])
BAG3.1.per.maf.gnomad.ALL.and.NFE$REVEL[BAG3.1.per.maf.gnomad.ALL.and.NFE$REVEL <= 0.5] <- NA

## Add Annovar annotation to this df
BAG3.1.per.maf.gnomad.ALL.and.NFE$Annovar_ExonicFunc.refGene <- BAG3_annovar.lt.1.perc.maf.gnom.AD$ExonicFunc.refGene[match(BAG3.1.per.maf.gnomad.ALL.and.NFE$KEY, BAG3_annovar.lt.1.perc.maf.gnom.AD$KEY)]


CLINVAR.Keys <- BAG3.1.per.maf.gnomad.ALL.and.NFE$KEY[grepl("Clinvar", BAG3.1.per.maf.gnomad.ALL.and.NFE$PRED_TYPE)]
BAG3.1.per.maf.gnomad.ALL.and.NFE$CLINVAR <- ifelse(BAG3.1.per.maf.gnomad.ALL.and.NFE$KEY %in% CLINVAR.Keys,"Y", "N")

LoF.Keys <- BAG3.1.per.maf.gnomad.ALL.and.NFE$KEY[grepl("LoF", BAG3.1.per.maf.gnomad.ALL.and.NFE$PRED_TYPE)]
BAG3.1.per.maf.gnomad.ALL.and.NFE$LoF <- ifelse(BAG3.1.per.maf.gnomad.ALL.and.NFE$KEY %in% LoF.Keys,"Y", "N")

MetaSVM.Keys <- BAG3.1.per.maf.gnomad.ALL.and.NFE$KEY[grepl("MetaSVM", BAG3.1.per.maf.gnomad.ALL.and.NFE$PRED_TYPE)]
BAG3.1.per.maf.gnomad.ALL.and.NFE$MetaSVM <- ifelse(BAG3.1.per.maf.gnomad.ALL.and.NFE$KEY %in% MetaSVM.Keys,"Y", "N")


BAG3.1.per.maf.gnomad.ALL.and.NFE <- as.data.frame(BAG3.1.per.maf.gnomad.ALL.and.NFE)
BAG3.1.per.maf.gnomad.ALL.and.NFE <- BAG3.1.per.maf.gnomad.ALL.and.NFE[!grepl("PRED_TYPE", colnames(BAG3.1.per.maf.gnomad.ALL.and.NFE))]


# Remove duplicate rows
BAG3.1.per.maf.gnomad.ALL.and.NFE <- BAG3.1.per.maf.gnomad.ALL.and.NFE %>% 
  distinct(.keep_all = TRUE)

BAG3.1.per.maf.gnomad.ALL.and.NFE$REVEL_TYPE <- ifelse(BAG3.1.per.maf.gnomad.ALL.and.NFE$REVEL > 0.05, "REVEL_Pathogenic", "Not_pathogenic")
write.table(BAG3.1.per.maf.gnomad.ALL.and.NFE, "Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/MERGED_sjlife1_2_PreQC/cleaned/annotation/snpEff/TTN_BAG3/BAG3.1.per.maf.gnomad.ALL.and.NFE.txt", quote = F, col.names = T, sep = "\t", row.names = F)


BAG3.1.per.maf.gnomad.ALL.and.NFE.clinvar.lof <- BAG3.1.per.maf.gnomad.ALL.and.NFE$KEY[BAG3.1.per.maf.gnomad.ALL.and.NFE$CLINVAR == "Y"| BAG3.1.per.maf.gnomad.ALL.and.NFE$LoF =="Y" ]
write.table(BAG3.1.per.maf.gnomad.ALL.and.NFE.clinvar.lof, "Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/MERGED_sjlife1_2_PreQC/cleaned/annotation/snpEff/TTN_BAG3/BAG3.1.per.maf.gnomad.ALL.and.NFE.clinvar.lof.txt", quote = F, col.names = F, sep = "\t", row.names = F)

BAG3.1.per.maf.gnomad.ALL.and.NFE.clinvar.lof.REVEL <- BAG3.1.per.maf.gnomad.ALL.and.NFE$KEY[BAG3.1.per.maf.gnomad.ALL.and.NFE$CLINVAR == "Y"| BAG3.1.per.maf.gnomad.ALL.and.NFE$LoF =="Y"|BAG3.1.per.maf.gnomad.ALL.and.NFE$REVEL_TYPE == "REVEL_Pathogenic" ]
write.table(BAG3.1.per.maf.gnomad.ALL.and.NFE.clinvar.lof.REVEL, "Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/MERGED_sjlife1_2_PreQC/cleaned/annotation/snpEff/TTN_BAG3/BAG3.1.per.maf.gnomad.ALL.and.NFE.clinvar.lof.REVEL.txt", quote = F, col.names = F, sep = "\t", row.names = F)


###########
## REVEL ##
###########
TITN_BAG3_annovar <- rbind.data.frame(saved.TITN_annovar, saved.BAG3_annovar)
TITN_BAG3_annovar$KEY <- paste0(TITN_BAG3_annovar$Otherinfo4,":", TITN_BAG3_annovar$Otherinfo5,":", 
                                TITN_BAG3_annovar$Otherinfo7,":", TITN_BAG3_annovar$Otherinfo8)
sum(BAG3_TITN$KEY %in% TITN_BAG3_annovar$KEY)

REVEL.gt.0.5.annovar <- TITN_BAG3_annovar[TITN_BAG3_annovar$REVEL > 0.5,]

# REVEL.gt.0.5 <- REVEL.gt.0.5.annovar$KEY

REVEL.gt.0.5.bed <- cbind.data.frame(REVEL.gt.0.5.annovar$Otherinfo4, REVEL.gt.0.5.annovar$Otherinfo5-1, REVEL.gt.0.5.annovar$Otherinfo5)



# write.table(REVEL.gt.0.5, "Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/MERGED_sjlife1_2_PreQC/cleaned/annotation/snpEff/TTN_BAG3/REVEL.gt.0.5.txt", quote = F, col.names = F, sep = "\t", row.names = F)

write.table(REVEL.gt.0.5.bed[grepl("chr2", REVEL.gt.0.5.bed$`REVEL.gt.0.5.annovar$Otherinfo4`),], "Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/MERGED_sjlife1_2_PreQC/cleaned/annotation/snpEff/TTN_BAG3/TTN_REVEL.gt.0.5.bed", quote = F, col.names = F, sep = "\t", row.names = F)
write.table(REVEL.gt.0.5.bed[grepl("chr10", REVEL.gt.0.5.bed$`REVEL.gt.0.5.annovar$Otherinfo4`),], "Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/MERGED_sjlife1_2_PreQC/cleaned/annotation/snpEff/TTN_BAG3/BAG3_REVEL.gt.0.5.bed", quote = F, col.names = F, sep = "\t", row.names = F)


################################################
## Extract common variants in SJLIFE and CCSS ##
################################################
# ttn_check <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/MERGED_sjlife1_2_PreQC/cleaned/annotation/snpEff/TTN_BAG3/check_ttn.txt")

## Check common vars between SJLIFE and CCSS_EXP
sjlife.maf.0.01 <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/MERGED_sjlife1_2_PreQC/cleaned/annotation/snpEff/TTN_BAG3/plink_out/all_var.sjlife.0.01.maf.txt")

ccss_exp <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3/ccss_exp.bim")
sum(sjlife.maf.0.01$V4 %in% ccss_exp$V4 )
sum(sjlife.maf.0.01$V2 %in% ccss_exp$V2 )

set.0 <- ccss_exp$V4[ccss_exp$V4 %in% sjlife.maf.0.01$V4]
sum(duplicated(set.0))
set.0[duplicated(set.0)]
set.1 <- sjlife.maf.0.01$V2[sjlife.maf.0.01$V4 %in% ccss_exp$V4] # match by position and extract variant ID
set.2 <- sjlife.maf.0.01$V2[sjlife.maf.0.01$V2 %in% ccss_exp$V2] # match by variant ID and extract matching variant ID

set.1[!set.1 %in% set.2] # print no direct match
# "chr2:178750594:G:A"          "chr2:178752047:A:G"          "chr2:178776816:C:A"          "chr2:178781235:C:T"          "chr10:119656670:GATGAAGGT:G"
# "chr2:178752047:A:G", "chr2:178781235:C:G", chr10:119656670:G:A, "chr10:119656670:GATGAAGGT:G # no match
direct.match <- set.1[set.1 %in% set.2]

replace.ccss_exp_id <- cbind.data.frame(old=c("chr2:178750594:G:T","chr2:178776816:C:T"), new=c("chr2:178750594:G:A", "chr2:178776816:C:A"))
# write.table(replace.ccss_exp_id, "Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3/replace.ccss_exp_id.txt", quote = F, col.names = F, sep = "\t", row.names = F)

ccss_bim <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3/ccss_exp.bim")

ccss_bim$V2[na.omit(match(replace.ccss_exp_id$old, ccss_bim$V2))] <- replace.ccss_exp_id$new
write.table(ccss_bim, "Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3/ccss_exp_edited.bim", quote = F, col.names = F, sep = "\t", row.names = F)


# Overlapping TITN P/LP and maf below 1% in gnomad_ALL and gnomad_NFE
ccss <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3/ALL_P_LP_combinations/TITN_ccss_exp.0.01maf.bim")
sjlife <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3/ALL_P_LP_combinations/TITN_sjlife.0.01maf.bim")

write.table(ccss$V2[ccss$V2 %in% sjlife$V2], "Z:/ResearchHome/Groups/sapkogrp/projects//Cardiotoxicity/common/ttn_bag3/ALL_P_LP_combinations/TITN_overlapping_vars_in_sjlife_ccss_exp_maf_0.01", quote = F, col.names = F, sep = "\t", row.names = F)

# Overlapping BAG3 P/LP
ccss <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3/ALL_P_LP_combinations/BAG3_ccss_exp.0.01maf.bim")
sjlife <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3/ALL_P_LP_combinations/BAG3_sjlife.0.01maf.bim")

write.table(ccss$V2[ccss$V2 %in% sjlife$V2], "Z:/ResearchHome/Groups/sapkogrp/projects//Cardiotoxicity/common/ttn_bag3/ALL_P_LP_combinations/BAG3_overlapping_vars_in_sjlife_ccss_exp_maf_0.01", quote = F, col.names = F, sep = "\t", row.names = F)
