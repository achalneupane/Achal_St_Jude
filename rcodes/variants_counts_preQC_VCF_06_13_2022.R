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

CLINVAR.unique <- CLINVAR[!duplicated(CLINVAR$KEY),]
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

# save.image("SNPEFF_clinvar_metaSVM_from_R_filtering_process_PreQC_VCF.RData")
load("SNPEFF_clinvar_metaSVM_from_R_filtering_process_PreQC_VCF.RData")





