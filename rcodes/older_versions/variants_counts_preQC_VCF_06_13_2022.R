###########################
## Achal Neupane         ##
## Date: 06/13/2022      ##
## VCF annotation PRE-QC ##
###########################
## First check how many variants in Zhaoming et al and Qin et al are also in our VCF dataset
# Read variant lists from Zhaoming et al 
library(data.table)
library(dplyr)
Sys.setlocale("LC_ALL", "C")

#############################
#############################
## SNPEFF ##### Annotation ##
#############################
#############################
# https://pcingola.github.io/SnpEff/adds/VCFannotationformat_v1.0.pdf
setwd("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/MERGED_sjlife1_2_PreQC/cleaned/annotation/snpEff/")
## read annotated SJLIFE annotated VCF 
# Loop over all chromosomes
chromosomes <- 1:22


## Now extract variants of ClinVar and MetaSVM significance
FINAL.VCF <- {}

capture.output (for( i in 1:length(chromosomes)){
print(paste0("Doing chromosome ", chromosomes[i]))
  
VCF <- fread(paste0("MERGED.SJLIFE.1.2.GATKv3.4.VQSR.chr", chromosomes[i], ".preQC_biallelic_renamed_ID_edited.vcf-annot-snpeff-dbnsfp-ExAC.0.3-clinvar.GRCh38.vcf.dbSNP155-FIELDS-simple.txt"))

#####################################
## Clinvar based P/LP VCF variants ##
#####################################
print(as.data.frame(table(VCF$CLNSIG)))
# Var1     Freq
# 1                                                          . 19388477
# 2                                                    Affects       51
# 3                                                     Benign   117561
# 4                                       Benign/Likely_benign    15102
# 5                                Benign/Likely_benign|_other        7
# 6                                              Benign|_other       14
# 7               Conflicting_interpretations_of_pathogenicity    19244
# 8  Conflicting_interpretations_of_pathogenicity|_association        6
# 9  Conflicting_interpretations_of_pathogenicity|_risk_factor        6
# 10                                             Likely_benign    62840
# 11                                         Likely_pathogenic      593
# 12                                                Pathogenic     2028
# 13                              Pathogenic/Likely_pathogenic      875
# 14                                   Pathogenic|_risk_factor        9
# 15                                    Uncertain_significance    47896
# 16                                               association        7
# 17                                              not_provided      535
# 18                                                protective        3
# 19                                   protective|_risk_factor        5
# 20                                               risk_factor       22

## Wanted Clinvar patterns
WANTED.types.clinvar <- c("^Pathogenic/Likely_pathogenic$|^Likely_pathogenic$|^Likely_pathogenic/Pathogenic$|^Pathogenic$|Pathogenic\\|_risk_factor|^Pathogenic\\|")
print(paste0("Total P or LP vars from clinvar: ", sum(grepl(WANTED.types.clinvar, VCF$CLNSIG, ignore.case = T))))

## Wanted MetaSVM patterns D (available patterns: D= Deleterious; T= Tolerated)
wanted.types.MetaSVM <- c("D")
print(paste0("Total Deleterious vars from MetaSVM: ", sum(grepl(wanted.types.MetaSVM, VCF$dbNSFP_MetaSVM_pred, ignore.case = T))))

VCF.clinvar <- VCF[grepl(WANTED.types.clinvar, VCF$CLNSIG, ignore.case = T),]
VCF.clinvar$PRED_TYPE <- "Clinvar"
VCF.MetaSVM <- VCF[grepl(wanted.types.MetaSVM, VCF$dbNSFP_MetaSVM_pred, ignore.case = T),]
VCF.MetaSVM$PRED_TYPE="MetaSVM"
VCF <- rbind.data.frame(VCF.clinvar, VCF.MetaSVM)
FINAL.VCF <- rbind.data.frame(FINAL.VCF, VCF)
}, file = "SNPEFF_clinvar_metaSVM_from_R_filtering_process.log")


# # First remove any leading and trailing spaces
FINAL.VCF$CHROM <- trimws(FINAL.VCF$CHROM, which = "both")
FINAL.VCF$POS <- trimws(FINAL.VCF$POS, which = "both")

FINAL.VCF$KEY <- paste(FINAL.VCF$CHROM, FINAL.VCF$POS, FINAL.VCF$REF, FINAL.VCF$ALT, sep =":")

# FINAL.VCF.unique <- FINAL.VCF[duplicated(FINAL.VCF$KEY),]

#############
## Clinvar ##
#############
CLINVAR <- FINAL.VCF[FINAL.VCF$PRED_TYPE == "Clinvar",]
table(CLINVAR$CLNSIG)
nrow(CLINVAR)
# 60892

#############
## MetaSVM ##
#############
MetaSVM <- FINAL.VCF[FINAL.VCF$PRED_TYPE == "MetaSVM",]

CLINVAR.unique <- CLINVAR[-which(duplicated(CLINVAR$KEY)),]
table(CLINVAR.unique$CLNSIG)

MetaSVM.unique <- MetaSVM[!duplicated(MetaSVM$KEY),]

nrow(CLINVAR.unique)
# 4705
nrow(MetaSVM.unique)
# 83479

######################
## Loss of Function ##
######################
LoF <- read.delim("LoF_variants.txt",  sep = "\t", header = T, check.names = F)
LoF$PRED_TYPE <- "LoF"
LoF$KEY <- paste(LoF$CHROM, LoF$POS, LoF$REF, LoF$ALT, sep = ":")
LoF.unique <- LoF[!duplicated(LoF$KEY),]
nrow(LoF.unique)
# 316386


# How many SNPs from Clinvar, MetaSVM and LoF
# Predicted.vars.inVCF <-  rbind.data.frame(CLINVAR, LoF, MetaSVM)
# table(Predicted.vars.inVCF$PRED_TYPE)
# # Clinvar     LoF MetaSVM 
# # 4705  316386   83479 
# 
# # remove duplicates
# Predicted.vars.inVCF.Unique <- Predicted.vars.inVCF[!duplicated(Predicted.vars.inVCF$KEY),]

Predicted.vars.inVCF.Unique <-  rbind.data.frame(CLINVAR.unique, LoF.unique, MetaSVM.unique)


#######################################################################
## Check which of the variants from the previous studies are present ##
#######################################################################
## Zhaoming et al ##
####################

zhaoming <- read.delim("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/MERGED_sjlife1_2_PreQC/cleaned/annotation/snpEff/zhaoming_et_al_SNV_INDEL_Search_list.txt", sep = "\t", header = T)
zhaoming.variants <- zhaoming[grepl("Y|N", zhaoming$Match_per_var),]

sum(zhaoming.variants$ID %in% Predicted.vars.inVCF.Unique$KEY)
# 291/299

###############
## Qin et al ##
###############
qin <- read.delim("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/MERGED_sjlife1_2_PreQC/cleaned/annotation/snpEff/qin_et_al_SNV_INDEL_Search_list.txt", sep = "\t", header = T)
qin.variants <- qin[grepl("Y|N", qin$Match_per_var),]

sum(qin.variants$ID %in% Predicted.vars.inVCF.Unique$KEY)
# 387/392

##################################
## How many variants per gene ? ##
##################################
df <- Predicted.vars.inVCF.Unique
group_and_concat <- df %>%
  dplyr::select(`ANN[*].GENE`, KEY) %>% 
  dplyr::group_by(`ANN[*].GENE`) %>%
  dplyr::summarise(all_variants_SJLIFE = paste(KEY, collapse = ","))

group_and_concat$counts_var_SJLIFE <- lengths(strsplit(group_and_concat$all_variants_SJLIFE, ","))

zhaoming.variants$Gene[!zhaoming.variants$Gene %in% group_and_concat$`ANN[*].GENE`]
# "MRE11A" == "MRE11"

qin.variants$Gene[!qin.variants$Gene %in% group_and_concat$`ANN[*].GENE`]
# "FAM175A" == "ABRAXAS1"
# "C1orf86" == "FAAP20"  
# "MRE11A" == "MRE11"
# "C19orf40" == "FAAP24"
# "SHFM1" == "SEM1"

## Recode above gene names
group_and_concat$`ANN[*].GENE`[group_and_concat$`ANN[*].GENE` == "MRE11"] <- "MRE11A"
group_and_concat$`ANN[*].GENE`[group_and_concat$`ANN[*].GENE` == "ABRAXAS1"] <- "FAM175A"
group_and_concat$`ANN[*].GENE`[group_and_concat$`ANN[*].GENE` == "FAAP20"] <- "C1orf86"
group_and_concat$`ANN[*].GENE`[group_and_concat$`ANN[*].GENE` == "FAAP24"] <- "C19orf40"
group_and_concat$`ANN[*].GENE`[group_and_concat$`ANN[*].GENE` == "SEM1"] <- "SHFM1"


sum(zhaoming.variants$Gene %in% group_and_concat$`ANN[*].GENE` )
# 299
sum(qin.variants$Gene %in% group_and_concat$`ANN[*].GENE`)
# 392

group_and_concat$Zhaoming_YN <- ifelse(group_and_concat$`ANN[*].GENE` %in% unique(zhaoming.variants$Gene), "Y", "N" )
group_and_concat$Qin_YN <- ifelse(group_and_concat$`ANN[*].GENE` %in% unique(qin.variants$Gene), "Y", "N" )


## Variants to be extracted from VCF file for genotypes
vars.to.extract.from.VCF <- Predicted.vars.inVCF.Unique[!duplicated(Predicted.vars.inVCF.Unique$KEY),]
vars.to.extract.from.VCF.bed <- cbind.data.frame(CHROM = vars.to.extract.from.VCF$CHROM, START = as.numeric(vars.to.extract.from.VCF$POS) -1, END = vars.to.extract.from.VCF$POS)
write.table(vars.to.extract.from.VCF.bed, "Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/annotated_variants/annotated_vars_from_PreQC_VCF_Clinvar_MetaSVM_LoF/Variants_from_annotation_Clinvar_MetaSVM_LoF_PreQC.bed", row.names = FALSE, quote = FALSE, col.names = F, sep = "\t")




## TTN and BAG 3 L/PL variants
# TITN: chr2:178525989-178807423
# BAG3: chr10:119651380-119677819

###################################
## With Clinvar and MetaSVM only ##
###################################

## TITN
clinvar.CHR2 <- CLINVAR.unique[grepl("^chr2$", CLINVAR.unique$CHROM),]
clinvar.CHR2$POS <- as.numeric(clinvar.CHR2$POS)
clinvar.CHR2 <- clinvar.CHR2[clinvar.CHR2$POS >= 178525989 & clinvar.CHR2$POS <= 178807423,]

MetaSVM.CHR2 <- MetaSVM.unique[grepl("^chr2$", MetaSVM.unique$CHROM),]
MetaSVM.CHR2$POS <- as.numeric(MetaSVM.CHR2$POS)
MetaSVM.CHR2 <- MetaSVM.CHR2[MetaSVM.CHR2$POS >= 178525989 & MetaSVM.CHR2$POS <= 178807423,]

TITN.df <- rbind.data.frame(clinvar.CHR2, MetaSVM.CHR2)

## BAG3
clinvar.CHR10 <- CLINVAR.unique[grepl("^chr10$", CLINVAR.unique$CHROM),]
clinvar.CHR10$POS <- as.numeric(clinvar.CHR10$POS)
clinvar.CHR10 <- clinvar.CHR10[clinvar.CHR10$POS >= 119651380 & clinvar.CHR10$POS <= 119677819,]

MetaSVM.CHR10 <- MetaSVM.unique[grepl("^chr10$", MetaSVM.unique$CHROM),]
MetaSVM.CHR10$POS <- as.numeric(MetaSVM.CHR10$POS)
MetaSVM.CHR10 <- MetaSVM.CHR10[MetaSVM.CHR10$POS >= 119651380 & MetaSVM.CHR10$POS <= 119677819,]

BAG3.df <- rbind.data.frame(clinvar.CHR10, MetaSVM.CHR10)


BAG3_TITN <- rbind.data.frame(BAG3.df,TITN.df)

BAG3_TITN_BED <- cbind.data.frame(BAG3_TITN$CHROM, BAG3_TITN$POS-1, BAG3_TITN$POS)

# write.table(BAG3_TITN, "Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/MERGED_sjlife1_2_PreQC/cleaned/annotation/snpEff/TTN_BAG3/TTN_BAG3_list_LP.txt", quote = F, col.names = T, sep = "\t")

TITN <- BAG3_TITN_BED[BAG3_TITN_BED$`BAG3_TITN$CHROM` == "chr2",]
TITN <- distinct(TITN)

BAG3 <- BAG3_TITN_BED[BAG3_TITN_BED$`BAG3_TITN$CHROM` == "chr10",]
BAG3 <- distinct(BAG3)

write.table(TITN, "Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/MERGED_sjlife1_2_PreQC/cleaned/annotation/snpEff/TTN_BAG3/TTN_LP_clinvar_metaSVM.bed", quote = F, col.names = F, sep = "\t", row.names = F)
write.table(BAG3, "Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/MERGED_sjlife1_2_PreQC/cleaned/annotation/snpEff/TTN_BAG3/BAG3_LP_clinvar_metaSVM.bed", quote = F, col.names = F, sep = "\t", row.names = F)


TITN_VCF.extracted <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/MERGED_sjlife1_2_PreQC/cleaned/annotation/snpEff/TTN_BAG3/TTN_VCF.txt", header = F)
TITN_VCF.extracted$V1[!TITN_VCF.extracted$V1 %in% TITN.df$KEY]
# [aneupane@noderome201 TTN_BAG3]$ grep chr2:178531620 TTN_VCF.txt
# chr2:178531620:G:A
# chr2:178531620:G:C
# [aneupane@noderome201 TTN_BAG3]$ grep chr2:178781235  TTN_VCF.txt
# chr2:178781235:C:T
# chr2:178781235:C:G


## Annotate file
TITN_BAG3.df <- read.delim("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/MERGED_sjlife1_2_PreQC/cleaned/annotation/snpEff/TTN_BAG3/TTN_BAG3_list_yadav.txt", header = T)
TITN_BAG3.df$KEY <- paste0("chr", TITN_BAG3.df$SNP)
TITN_BAG3.df$FULLKEY1 <- paste0("chr", TITN_BAG3.df$SNP, ":", TITN_BAG3.df$A2, ":", TITN_BAG3.df$A1)
TITN_BAG3.df$FULLKEY2 <- paste0("chr", TITN_BAG3.df$SNP, ":", TITN_BAG3.df$A1, ":", TITN_BAG3.df$A2)

# READ annotation VCF
chr2_178 <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/MERGED_sjlife1_2_PreQC/cleaned/annotation/snpEff/TTN_BAG3/chr2_178_list.txt", header = T)
chr10_119 <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/MERGED_sjlife1_2_PreQC/cleaned/annotation/snpEff/TTN_BAG3/chr10_119_list.txt", header = T)

chr2_10 <- rbind.data.frame(chr2_178, chr10_119)

chr2_10$KEY <- paste0(chr2_10$CHROM,":",  chr2_10$POS)
chr2_10$FULLKEY1 <- paste0(chr2_10$CHROM,":",  chr2_10$POS, ":", chr2_10$REF, ":", chr2_10$ALT)

sum(TITN_BAG3.df$KEY %in% chr2_10$KEY)
# 2601
sum(TITN_BAG3.df$FULLKEY1 %in% chr2_10$FULLKEY1)
# 2302
sum(TITN_BAG3.df$FULLKEY2 %in% chr2_10$FULLKEY1)
# 318
TITN_BAG3.df1 <- cbind(TITN_BAG3.df, chr2_10[match(TITN_BAG3.df$FULLKEY1, chr2_10$FULLKEY1),c("ANN....IMPACT", "ANN....FEATURE", "ANN....EFFECT")])
TITN_BAG3.df1 <-TITN_BAG3.df1[!is.na(TITN_BAG3.df1$ANN....IMPACT),]

TITN_BAG3.df2 <- cbind(TITN_BAG3.df, chr2_10[match(TITN_BAG3.df$FULLKEY2, chr2_10$FULLKEY1),c("ANN....IMPACT", "ANN....FEATURE", "ANN....EFFECT")])
TITN_BAG3.df2 <-TITN_BAG3.df2[!is.na(TITN_BAG3.df2$ANN....IMPACT),]

TITN_BAG3.df3 <- rbind.data.frame(TITN_BAG3.df1, TITN_BAG3.df2)
TITN_BAG3.df3 <- TITN_BAG3.df3[!duplicated(TITN_BAG3.df3$FULLKEY1),]

TITN_BAG3.df <- cbind.data.frame(TITN_BAG3.df, TITN_BAG3.df3[match(TITN_BAG3.df$FULLKEY1,TITN_BAG3.df3$FULLKEY1), c("ANN....IMPACT", "ANN....FEATURE", "ANN....EFFECT")])

write.table(TITN_BAG3.df, "Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/MERGED_sjlife1_2_PreQC/cleaned/annotation/snpEff/TTN_BAG3/Yadav_TITN_BAG3_list.txt", quote = F, col.names = T, sep = "\t", row.names = F)


###################################
## With Clinvar, MetaSVM and LoF ##
###################################
## TITN
clinvar.CHR2 <- CLINVAR.unique[grepl("^chr2$", CLINVAR.unique$CHROM),]
clinvar.CHR2$POS <- as.numeric(clinvar.CHR2$POS)
clinvar.CHR2 <- clinvar.CHR2[clinvar.CHR2$POS >= 178525989 & clinvar.CHR2$POS <= 178807423,]

MetaSVM.CHR2 <- MetaSVM.unique[grepl("^chr2$", MetaSVM.unique$CHROM),]
MetaSVM.CHR2$POS <- as.numeric(MetaSVM.CHR2$POS)
MetaSVM.CHR2 <- MetaSVM.CHR2[MetaSVM.CHR2$POS >= 178525989 & MetaSVM.CHR2$POS <= 178807423,]

LoF.CHR2 <- LoF.unique[grepl("^chr2$", LoF.unique$CHROM),]
LoF.CHR2$POS <- as.numeric(LoF.CHR2$POS)
LoF.CHR2 <- LoF.CHR2[LoF.CHR2$POS >= 178525989 & LoF.CHR2$POS <= 178807423,]

TITN.df <- rbind.data.frame(clinvar.CHR2, MetaSVM.CHR2, LoF.CHR2)

## BAG3
clinvar.CHR10 <- CLINVAR.unique[grepl("^chr10$", CLINVAR.unique$CHROM),]
clinvar.CHR10$POS <- as.numeric(clinvar.CHR10$POS)
clinvar.CHR10 <- clinvar.CHR10[clinvar.CHR10$POS >= 119651380 & clinvar.CHR10$POS <= 119677819,]

MetaSVM.CHR10 <- MetaSVM.unique[grepl("^chr10$", MetaSVM.unique$CHROM),]
MetaSVM.CHR10$POS <- as.numeric(MetaSVM.CHR10$POS)
MetaSVM.CHR10 <- MetaSVM.CHR10[MetaSVM.CHR10$POS >= 119651380 & MetaSVM.CHR10$POS <= 119677819,]

LoF.CHR10 <- LoF.unique[grepl("^chr10$", LoF.unique$CHROM),]
LoF.CHR10$POS <- as.numeric(LoF.CHR10$POS)
LoF.CHR10 <- LoF.CHR10[LoF.CHR10$POS >= 119651380 & LoF.CHR10$POS <= 119677819,]

BAG3.df <- rbind.data.frame(clinvar.CHR10, MetaSVM.CHR10, LoF.CHR10)


BAG3_TITN <- rbind.data.frame(BAG3.df,TITN.df)

# BAG3_TITN$SNP <- paste0(BAG3_TITN$CHROM, ":", BAG3_TITN$POS)

write.table(BAG3_TITN, "Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/MERGED_sjlife1_2_PreQC/cleaned/annotation/snpEff/TTN_BAG3/TTN_BAG3_list_LP.txt", quote = F, col.names = T, sep = "\t")

BAG3_TITN_BED <- cbind.data.frame(BAG3_TITN$CHROM, BAG3_TITN$POS-1, BAG3_TITN$POS)

TITN <- BAG3_TITN_BED[BAG3_TITN_BED$`BAG3_TITN$CHROM` == "chr2",]
TITN <- distinct(TITN)

BAG3 <- BAG3_TITN_BED[BAG3_TITN_BED$`BAG3_TITN$CHROM` == "chr10",]
BAG3 <- distinct(BAG3)

write.table(TITN, "Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/MERGED_sjlife1_2_PreQC/cleaned/annotation/snpEff/TTN_BAG3/TTN_LP_clinvar_metaSVM_LoF.bed", quote = F, col.names = F, sep = "\t", row.names = F)
write.table(BAG3, "Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/MERGED_sjlife1_2_PreQC/cleaned/annotation/snpEff/TTN_BAG3/BAG3_LP_clinvar_metaSVM_LoF.bed", quote = F, col.names = F, sep = "\t", row.names = F)


TITN_VCF.extracted <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/MERGED_sjlife1_2_PreQC/cleaned/annotation/snpEff/TTN_BAG3/TTN_LP_clinvar_metaSVM_LoF_v3.txt", header = F)
TITN_VCF.extracted$V1[!TITN_VCF.extracted$V1 %in% unique(TITN.df$KEY)]

unique(TITN.df$KEY)[unique(TITN.df$KEY) %in% TITN_VCF.extracted$V1]
TITN_VCF.extracted$V1[duplicated(TITN_VCF.extracted$V1)]



BAG3_VCF.extracted <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/MERGED_sjlife1_2_PreQC/cleaned/annotation/snpEff/TTN_BAG3/BAG3_VCF_clinvar_metaSVM_LoF_vars.txt", header = F)
BAG3_VCF.extracted$V1[!BAG3_VCF.extracted$V1 %in% unique(BAG3.df$KEY)]

# save.image("SNPEFF_clinvar_metaSVM_LoF_from_R_filtering_process_PreQC_VCF.RData")
load("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/MERGED_sjlife1_2_PreQC/cleaned/annotation/snpEff/SNPEFF_clinvar_metaSVM_LoF_from_R_filtering_process_PreQC_VCF.RData")


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

TITN.1.per.maf.gnomad.ALL.and.NFE <- BAG3_TITN[BAG3_TITN$KEY %in% TITN_annovar.lt.1.perc.maf.gnom.AD$KEY,]
TITN.1.per.maf.gnomad.ALL.and.NFE$gnomAD_genome_ALL <- TITN_annovar.lt.1.perc.maf.gnom.AD$gnomAD_genome_ALL[match(TITN.1.per.maf.gnomad.ALL.and.NFE$KEY, TITN_annovar.lt.1.perc.maf.gnom.AD$KEY)]
TITN.1.per.maf.gnomad.ALL.and.NFE$gnomAD_genome.NFE <- TITN_annovar.lt.1.perc.maf.gnom.AD$gnomAD_genome_NFE[match(TITN.1.per.maf.gnomad.ALL.and.NFE$KEY, TITN_annovar.lt.1.perc.maf.gnom.AD$KEY)]
TITN.1.per.maf.gnomad.ALL.and.NFE$REVEL <- as.numeric(TITN_annovar.lt.1.perc.maf.gnom.AD$REVEL[match(TITN.1.per.maf.gnomad.ALL.and.NFE$KEY, TITN_annovar.lt.1.perc.maf.gnom.AD$KEY)])
TITN.1.per.maf.gnomad.ALL.and.NFE$REVEL[TITN.1.per.maf.gnomad.ALL.and.NFE$REVEL <= 0.5] <- NA

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

BAG3.1.per.maf.gnomad.ALL.and.NFE <- BAG3_TITN[BAG3_TITN$KEY %in% BAG3_annovar.lt.1.perc.maf.gnom.AD$KEY,]
BAG3.1.per.maf.gnomad.ALL.and.NFE$gnomAD_genome_ALL <- BAG3_annovar.lt.1.perc.maf.gnom.AD$gnomAD_genome_ALL[match(BAG3.1.per.maf.gnomad.ALL.and.NFE$KEY, BAG3_annovar.lt.1.perc.maf.gnom.AD$KEY)]
BAG3.1.per.maf.gnomad.ALL.and.NFE$gnomAD_genome.NFE <- BAG3_annovar.lt.1.perc.maf.gnom.AD$gnomAD_genome_NFE[match(BAG3.1.per.maf.gnomad.ALL.and.NFE$KEY, BAG3_annovar.lt.1.perc.maf.gnom.AD$KEY)]
BAG3.1.per.maf.gnomad.ALL.and.NFE$REVEL <- as.numeric(BAG3_annovar.lt.1.perc.maf.gnom.AD$REVEL[match(BAG3.1.per.maf.gnomad.ALL.and.NFE$KEY, BAG3_annovar.lt.1.perc.maf.gnom.AD$KEY)])
BAG3.1.per.maf.gnomad.ALL.and.NFE$REVEL[BAG3.1.per.maf.gnomad.ALL.and.NFE$REVEL <= 0.5] <- NA


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
