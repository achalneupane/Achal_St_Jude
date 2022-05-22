######################
## Achal Neupane    ##
## Date: 05/12/2022 ##
######################
## Read variant lists from Zhaoming et al 
Sys.setlocale("LC_ALL", "C")
zhaoming.etal.vars <- read.delim("Z:/ResearchHome/ClusterHome/aneupane/St_Jude/Yadav_Sapkota/additional_papers/Zhaoming_Wang_Genetic_risk_for_subsequent_neoplasms_JCO.2018.77.8589/P-PL-sjlife-genetics-sn_SNV_INDELS.txt", sep = "\t", header =T, stringsAsFactors = F)
head(zhaoming.etal.vars)
dim(zhaoming.etal.vars)
var.classification <- table(zhaoming.etal.vars$Classification)
var.classification
# LP   P 
# 160 188

zhaoming.etal.vars$STUDY <- "Genetic_Risk_for_SN"
zhaoming.etal.vars$KEY.pos <- paste0("chr", zhaoming.etal.vars$Chr, ":", zhaoming.etal.vars$Pos_GRCh38)
zhaoming.etal.vars$KEY.varID <- paste0("chr", zhaoming.etal.vars$Chr, ":", zhaoming.etal.vars$Pos_GRCh38, ":", zhaoming.etal.vars$Reference_Allele, ":", zhaoming.etal.vars$Mutant_Allele)

## Reading variants from Qin et al (DNA-repair)
qin.etal.vars <- read.table("Z:/ResearchHome/ClusterHome/aneupane/St_Jude/Yadav_Sapkota/additional_papers/Na_Qin_Pathogenic Germline Mutations in DNA Repair Genes in Combination With Cancer Treatment Exposures/supplementary_DS_jco.19.02760.txt", sep = "\t", header =T, stringsAsFactors = F)
head(qin.etal.vars)

qin.etal.vars$STUDY <- "DNA_Repair"
qin.etal.vars$KEY.pos <- paste0("chr", qin.etal.vars$Chr, ":", qin.etal.vars$Pos_GRCh38)
qin.etal.vars$KEY.varID <- paste0("chr", qin.etal.vars$Chr, ":", qin.etal.vars$Pos_GRCh38, ":", qin.etal.vars$Reference_Allele, ":", qin.etal.vars$Mutant_Allele)

## Check if the variants in Zhaoming et al and Qin et al are also in our datasets
zhaoming.etal.vars$Chr <- trimws(zhaoming.etal.vars$Chr, which = "both")
zhaoming.etal.vars$Pos_GRCh38 <- trimws(zhaoming.etal.vars$Pos_GRCh38, which = "both")

qin.etal.vars$Chr <- trimws(qin.etal.vars$Chr, which = "both")
qin.etal.vars$Pos_GRCh38 <- trimws(qin.etal.vars$Pos_GRCh38, which = "both")


#############################
#############################
## SNPEFF ##### Annotation ##
#############################
#############################
setwd("Z:/ResearchHome/Groups/sapkogrp/projects/SJLIFE_WGS/common/sjlife/MERGED_SJLIFE_1_2/annotation/SNPEFF_ANNOTATION/")
## read annotated SJLIFE annotated VCF 
library(data.table)


# Loop over all chromosomes
chromosomes <- 1:22

# First check how many variants in Zhaoming et al and Qin et al are also in our VCF dataset
## Here 
FINAL.VCF <- {}

capture.output (for( i in 1:length(chromosomes)){
print(paste0("Doing chromosome ", chromosomes[i]))
  
VCF <- fread(paste0("MERGED.SJLIFE.1.2.GATKv3.4.VQSR.chr", chromosomes[i], ".PASS.decomposed.vcf-annot-snpeff-dbnsfp-ExAC.0.3-clinvar.GRCh38.vcf.dbSNP155-FIELDS-simple.txt"))


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

save.image("SNPEFF_clinvar_metaSVM_from_R_filtering_process.RData")

load("SNPEFF_clinvar_metaSVM_from_R_filtering_process.RData")

#######################################################################
## Check which of the variants from the previous studies are present ##
#######################################################################
# # First remove any leading and trailing spaces
FINAL.VCF$CHROM <- trimws(FINAL.VCF$CHROM, which = "both")
FINAL.VCF$POS <- trimws(FINAL.VCF$POS, which = "both")


as.data.frame(table(FINAL.VCF$`ANN[*].EFFECT`))


CLINVAR <- FINAL.VCF[grepl("Clinvar", FINAL.VCF$PRED_TYPE),]
MetaSVM <- FINAL.VCF[grepl("MetaSVM", FINAL.VCF$PRED_TYPE),]


## Clinvar and MetaSVM
FINAL.VCF$KEY.pos <- paste(FINAL.VCF$CHROM, FINAL.VCF$POS, sep = ":")
FINAL.VCF$KEY.varID <- paste(FINAL.VCF$CHROM, FINAL.VCF$POS, FINAL.VCF$REF, FINAL.VCF$ALT, sep = ":")


## Zhaoming et al
length(unique(zhaoming.etal.vars$KEY.pos)) # 295
sum(unique(zhaoming.etal.vars$KEY.pos) %in% FINAL.VCF$KEY.pos)
# 126

length(unique(zhaoming.etal.vars$KEY.varID)) # 299
sum(unique(zhaoming.etal.vars$KEY.varID) %in% FINAL.VCF$KEY.varID)
# 123

## Qin et al
length(unique(qin.etal.vars$KEY.pos)) # 389
sum(unique(qin.etal.vars$KEY.pos) %in% FINAL.VCF$KEY.pos)
# 94
length(unique(qin.etal.vars$KEY.varID)) # 392
sum(unique(qin.etal.vars$KEY.varID) %in% FINAL.VCF$KEY.varID)
# 95

############################
# check how many samples are from CLINVAR
## Zhaoming et al
length(unique(zhaoming.etal.vars$KEY.pos)) # 295
sum(unique(zhaoming.etal.vars$KEY.pos) %in% CLINVAR$KEY.pos)
# 115

length(unique(zhaoming.etal.vars$KEY.varID)) # 299
sum(unique(zhaoming.etal.vars$KEY.varID) %in% CLINVAR$KEY.varID)
# 114

## Qin et al
length(unique(qin.etal.vars$KEY.pos)) # 389
sum(unique(qin.etal.vars$KEY.pos) %in% CLINVAR$KEY.pos)
# 91
length(unique(qin.etal.vars$KEY.varID)) # 392
sum(unique(qin.etal.vars$KEY.varID) %in% CLINVAR$KEY.varID)
# 91

# check how many samples are from MetaSVM
length(unique(zhaoming.etal.vars$KEY.pos)) # 295
sum(unique(zhaoming.etal.vars$KEY.pos) %in% MetaSVM$KEY.pos)
# 30

length(unique(zhaoming.etal.vars$KEY.varID)) # 299
sum(unique(zhaoming.etal.vars$KEY.varID) %in% MetaSVM$KEY.varID)
# 28

## Qin et al
length(unique(qin.etal.vars$KEY.pos)) # 389
sum(unique(qin.etal.vars$KEY.pos) %in% MetaSVM$KEY.pos)
# 19
length(unique(qin.etal.vars$KEY.varID)) # 392
sum(unique(qin.etal.vars$KEY.varID) %in% MetaSVM$KEY.varID)
# 20

write.table(paste0(unique(zhaoming.etal.vars$KEY.varID),":"), "zhaoming_variants.txt", row.names = F, col.names = F, quote = F)
write.table(paste0(unique(zhaoming.etal.vars$KEY.pos),":"), "zhaoming_variant_sites.txt", row.names = F, col.names = F, quote = F)

