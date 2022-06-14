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

FINAL.VCF$KEY <- paste(FINAL.VCF$CHROM, FINAL.VCF$POS, FINAL.VCF$REF, FINAL.VCF$ALT, sep =":")

# FINAL.VCF.unique <- FINAL.VCF[duplicated(FINAL.VCF$KEY),]

CLINVAR <- FINAL.VCF[FINAL.VCF$PRED_TYPE == "Clinvar",]
table(CLINVAR$CLNSIG)
nrow(CLINVAR)
# 60892

MetaSVM <- FINAL.VCF[FINAL.VCF$PRED_TYPE == "MetaSVM",]

CLINVAR.unique <- CLINVAR[!duplicated(CLINVAR$KEY),]
table(CLINVAR.unique$CLNSIG)

MetaSVM.unique <- MetaSVM[!duplicated(MetaSVM$KEY),]

nrow(CLINVAR.unique)
# 4705
nrow(MetaSVM.unique)
# 83479

# How many from Clinvar, MetaSVM and LoF
table(FINAL.VCF.unique$PRED_TYPE)




## Read Qin and Zhaoming variants from the previous list
zhaoming <- read.delim("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/MERGED_sjlife1_2_PreQC/cleaned/annotation/snpEff/zhaoming_et_al_SNV_INDEL_Search_list.txt", sep = "\t", header = T)
zhaoming.variants <- zhaoming[grepl("Y|N", zhaoming$Match_per_var),]


#######################################################################
## Check which of the variants from the previous studies are present ##
#######################################################################
# # First remove any leading and trailing spaces
FINAL.VCF$CHROM <- trimws(FINAL.VCF$CHROM, which = "both")
FINAL.VCF$POS <- trimws(FINAL.VCF$POS, which = "both")


## Label indels and SNVs
qin.etal.vars$varTypes <- NULL
qin.etal.vars$varTypes[grepl("^A$|^T$|^G$|^C$", qin.etal.vars$Reference_Allele) & grepl("^A$|^T$|^G$|^C$", qin.etal.vars$Mutant_Allele)] <- "SNV"
qin.etal.vars$CHROM <- paste0("chr",qin.etal.vars$Chr)
qin.etal.vars$START <- as.numeric(qin.etal.vars$Pos_GRCh38)
qin.etal.vars$END <- as.numeric(qin.etal.vars$Pos_GRCh38)
qin.etal.vars$varTypes [!grepl("SNV", qin.etal.vars$varTypes)] <- "INDEL"
qin.etal.vars$END [grepl("INDEL", qin.etal.vars$varTypes)] <- NA
qin.etal.vars$END [grepl("bp", qin.etal.vars$Reference_Allele)] <- as.numeric(gsub("bp| ","", qin.etal.vars$Reference_Allele 
                                                                                   [grepl("bp", qin.etal.vars$Reference_Allele)])) +
  as.numeric(qin.etal.vars$START [grepl("bp", qin.etal.vars$Reference_Allele)]) + 1


qin.etal.vars$END [grepl("bp", qin.etal.vars$Mutant_Allele)] <- as.numeric(gsub("bp| ","", qin.etal.vars$Mutant_Allele 
                                                                                [grepl("bp", qin.etal.vars$Mutant_Allele)])) +
  as.numeric(qin.etal.vars$START [grepl("bp", qin.etal.vars$Mutant_Allele)]) + 1


qin.etal.vars$END[is.na(qin.etal.vars$END)&qin.etal.vars$Reference_Allele !="-"] <- nchar(gsub("-", "", qin.etal.vars$Reference_Allele))[is.na(qin.etal.vars$END)&qin.etal.vars$Reference_Allele !="-"] + 
  + as.numeric(qin.etal.vars$START[is.na(qin.etal.vars$END)&qin.etal.vars$Reference_Allele !="-"]) +1


qin.etal.vars$END[is.na(qin.etal.vars$END)&qin.etal.vars$Mutant_Allele !="-"] <- nchar(gsub("-", "", qin.etal.vars$Mutant_Allele))[is.na(qin.etal.vars$END)&qin.etal.vars$Mutant_Allele !="-"] + 
  + as.numeric(qin.etal.vars$START[is.na(qin.etal.vars$END)&qin.etal.vars$Mutant_Allele !="-"]) +1

qin.etal.vars$START[qin.etal.vars$varTypes == "INDEL"] <- as.numeric(qin.etal.vars$START) [qin.etal.vars$varTypes == "INDEL"]-1

## BED file
write.table(cbind.data.frame(qin.etal.vars$CHROM,qin.etal.vars$START, qin.etal.vars$END), "qin_et_al_variants.bed", row.names = F, col.names = F, quote = F, sep = "\t")

# Now check how many of the SNVs are in LOF and CLINVAR annotation in our dataset
qin.SNV <- qin.etal.vars[qin.etal.vars$varTypes == "SNV",]
qin.INDEL <- qin.etal.vars[qin.etal.vars$varTypes == "INDEL",]

length(unique(qin.SNV$KEY.varID))
# 211
sum(unique(qin.SNV$KEY.varID) %in% FINAL.VCF$KEY.varID[FINAL.VCF$PRED_TYPE == "Clinvar"])
# 91
sum(unique(qin.SNV$KEY.varID) %in% FINAL.VCF$KEY.varID[FINAL.VCF$PRED_TYPE == "MetaSVM"])
# 20
## CLINVAR+METASVM
sum(unique(qin.SNV$KEY.varID) %in% FINAL.VCF$KEY.varID)
# 95

# In LOF
LOF.df <- fread("LoF_variants_ID.txt", header = F)
LOF.df$Key.Pos <- sub('^([^:]+:[^:]+).*', '\\1', LOF.df$V1)
LOF <- as.character(LOF.df$V1)
sum(unique(qin.SNV$KEY.varID) %in% LOF)
# 188
# In everything
sum(unique(qin.SNV$KEY.varID) %in% c(LOF, FINAL.VCF$KEY.varID))
# 205




########################################
## Now repeat this for Zhaoming et al ##
########################################
## Label indels and SNVs
zhaoming.etal.vars$varTypes <- NULL
zhaoming.etal.vars$varTypes[grepl("^A$|^T$|^G$|^C$", zhaoming.etal.vars$Reference_Allele) & grepl("^A$|^T$|^G$|^C$", zhaoming.etal.vars$Mutant_Allele)] <- "SNV"
zhaoming.etal.vars$CHROM <- paste0("chr",zhaoming.etal.vars$Chr)
zhaoming.etal.vars$START <- as.numeric(zhaoming.etal.vars$Pos_GRCh38)
zhaoming.etal.vars$END <- as.numeric(zhaoming.etal.vars$Pos_GRCh38)
zhaoming.etal.vars$varTypes [!grepl("SNV", zhaoming.etal.vars$varTypes)] <- "INDEL"
zhaoming.etal.vars$END [grepl("INDEL", zhaoming.etal.vars$varTypes)] <- NA
zhaoming.etal.vars$END [grepl("bp", zhaoming.etal.vars$Reference_Allele)] <- as.numeric(gsub("bp| ","", zhaoming.etal.vars$Reference_Allele 
                                                                                   [grepl("bp", zhaoming.etal.vars$Reference_Allele)])) +
  as.numeric(zhaoming.etal.vars$START [grepl("bp", zhaoming.etal.vars$Reference_Allele)]) + 1


zhaoming.etal.vars$END [grepl("bp", zhaoming.etal.vars$Mutant_Allele)] <- as.numeric(gsub("bp| ","", zhaoming.etal.vars$Mutant_Allele 
                                                                                [grepl("bp", zhaoming.etal.vars$Mutant_Allele)])) +
  as.numeric(zhaoming.etal.vars$START [grepl("bp", zhaoming.etal.vars$Mutant_Allele)]) + 1


zhaoming.etal.vars$END[is.na(zhaoming.etal.vars$END)&zhaoming.etal.vars$Reference_Allele !="-"] <- nchar(gsub("-", "", zhaoming.etal.vars$Reference_Allele))[is.na(zhaoming.etal.vars$END)&zhaoming.etal.vars$Reference_Allele !="-"] + 
  + as.numeric(zhaoming.etal.vars$START[is.na(zhaoming.etal.vars$END)&zhaoming.etal.vars$Reference_Allele !="-"]) +1


zhaoming.etal.vars$END[is.na(zhaoming.etal.vars$END)&zhaoming.etal.vars$Mutant_Allele !="-"] <- nchar(gsub("-", "", zhaoming.etal.vars$Mutant_Allele))[is.na(zhaoming.etal.vars$END)&zhaoming.etal.vars$Mutant_Allele !="-"] + 
  + as.numeric(zhaoming.etal.vars$START[is.na(zhaoming.etal.vars$END)&zhaoming.etal.vars$Mutant_Allele !="-"]) +1

zhaoming.etal.vars$START[zhaoming.etal.vars$varTypes == "INDEL"] <- as.numeric(zhaoming.etal.vars$START) [zhaoming.etal.vars$varTypes == "INDEL"]-1

## BED file
write.table(cbind.data.frame(zhaoming.etal.vars$CHROM,zhaoming.etal.vars$START, zhaoming.etal.vars$END), "zhaoming_et_al_variants.bed", row.names = F, col.names = F, quote = F, sep = "\t")


# Now check how many of the SNVs are in LOF and CLINVAR annotation in our dataset

zhaoming.SNV <- zhaoming.etal.vars[zhaoming.etal.vars$varTypes == "SNV",]
zhaoming.INDEL <- zhaoming.etal.vars[zhaoming.etal.vars$varTypes == "INDEL",]

length(unique(zhaoming.SNV$KEY.varID))
# 166
sum(unique(zhaoming.SNV$KEY.varID) %in% FINAL.VCF$KEY.varID[FINAL.VCF$PRED_TYPE == "Clinvar"])
# 114
sum(unique(zhaoming.SNV$KEY.varID) %in% FINAL.VCF$KEY.varID[FINAL.VCF$PRED_TYPE == "MetaSVM"])
# 28
## CLINVAR+METASVM
sum(unique(zhaoming.SNV$KEY.varID) %in% FINAL.VCF$KEY.varID)
# 123

# In LOF
LOF <- fread("LoF_variants_ID.txt", header = F)
LOF <- as.character(LOF$V1)
sum(unique(zhaoming.SNV$KEY.varID) %in% LOF)
# 128
# In everything
sum(unique(zhaoming.SNV$KEY.varID) %in% c(LOF, FINAL.VCF$KEY.varID))
# 155

#############################
## NOw checking for INDELS ##
#############################
## Zhaoming et al ##
####################
unique(zhaoming.INDEL$KEY.varID)
zhaoming.INDEL <- zhaoming.INDEL %>% distinct(KEY.varID, .keep_all = TRUE)

zhaoming.INDEL$Pos_GRCh38 <- as.numeric(zhaoming.INDEL$Pos_GRCh38)
zhaoming.INDEL$Pos_GRCh38.minus1 <- zhaoming.INDEL$Pos_GRCh38-1
zhaoming.INDEL$KEY.pos.indels <- paste0("chr", zhaoming.INDEL$Chr, ":", zhaoming.INDEL$Pos_GRCh38.minus1)
zhaoming.INDEL.edited <- cbind.data.frame(KEY.varID = zhaoming.INDEL$KEY.varID, CHROM = zhaoming.INDEL$CHROM,
                                          START = zhaoming.INDEL$START, END = zhaoming.INDEL$END)
zhaoming.INDEL.edited$KEY.varID <- gsub(" ","",zhaoming.INDEL.edited$KEY.varID)
zhaoming.INDEL.edited$tabix_query <- paste0(zhaoming.INDEL.edited$CHROM, ":", zhaoming.INDEL.edited$START, "-", zhaoming.INDEL.edited$END)
zhaoming.INDEL.edited$tabix_query_PM_10bps <- paste0(zhaoming.INDEL.edited$CHROM, ":", zhaoming.INDEL.edited$START-9, "-", zhaoming.INDEL.edited$END+9)
write.table(zhaoming.INDEL.edited, "zhaoming_et_al_variants_INDEL.bed", row.names = F, col.names = F, quote = F, sep = "\t")
## Create another bed with INDELs +/- 10 and SNVs with START-1
BED.zhaoming.indel.snv <- cbind.data.frame(CHR = zhaoming.SNV$CHROM, START = zhaoming.SNV$START-1, END = zhaoming.SNV$END)
CHR <- sapply(strsplit(as.character(zhaoming.INDEL.edited$tabix_query_PM_10bps),':'), "[", 1)
START <- gsub("-.*","",sapply(strsplit(as.character(zhaoming.INDEL.edited$tabix_query_PM_10bps),':'), "[", 2))
END <- gsub(".*-","",sapply(strsplit(as.character(zhaoming.INDEL.edited$tabix_query_PM_10bps),':'), "[", 2))
BED.zhaoming.indel.snv <- rbind.data.frame(BED.zhaoming.indel.snv, cbind.data.frame(CHR,START, END))

BED.zhaoming.indel.snv$KEY <- paste(BED.zhaoming.indel.snv$CHR,BED.zhaoming.indel.snv$START, BED.zhaoming.indel.snv$END, sep= ":")
BED.zhaoming.indel.snv <- BED.zhaoming.indel.snv %>% distinct(KEY, .keep_all = TRUE)
write.table(BED.zhaoming.indel.snv[1:3], "zhaoming_et_al_SNV_start_minus_1_INDEL_plus_min_10.bed", row.names = F, col.names = F, quote = F, sep = "\t")


###############
## Qin et al ##
###############
unique(qin.INDEL$KEY.varID)
qin.INDEL <- qin.INDEL %>% distinct(KEY.varID, .keep_all = TRUE)

qin.INDEL$Pos_GRCh38 <- as.numeric(qin.INDEL$Pos_GRCh38)
qin.INDEL$Pos_GRCh38.minus1 <- qin.INDEL$Pos_GRCh38-1
qin.INDEL$KEY.pos.indels <- paste0("chr", qin.INDEL$Chr, ":", qin.INDEL$Pos_GRCh38.minus1)
qin.INDEL.edited <- cbind.data.frame(KEY.varID = qin.INDEL$KEY.varID, CHROM = qin.INDEL$CHROM,
                                          START = qin.INDEL$START, END = qin.INDEL$END)
qin.INDEL.edited$KEY.varID <- gsub(" ","",qin.INDEL.edited$KEY.varID)
qin.INDEL.edited$tabix_query <- paste0(qin.INDEL.edited$CHROM, ":", qin.INDEL.edited$START, "-", qin.INDEL.edited$END)
qin.INDEL.edited$tabix_query_PM_10bps <- paste0(qin.INDEL.edited$CHROM, ":", qin.INDEL.edited$START-9, "-", qin.INDEL.edited$END+9)
write.table(qin.INDEL.edited, "qin_et_al_variants_INDEL.bed", row.names = F, col.names = F, quote = F, sep = "\t")


## Create another bed with INDELs +/- 10 and SNVs with START-1
BED.qin.indel.snv <- cbind.data.frame(CHR = qin.SNV$CHROM, START = qin.SNV$START-1, END = qin.SNV$END)
CHR <- sapply(strsplit(as.character(qin.INDEL.edited$tabix_query_PM_10bps),':'), "[", 1)
START <- gsub("-.*","",sapply(strsplit(as.character(qin.INDEL.edited$tabix_query_PM_10bps),':'), "[", 2))
END <- gsub(".*-","",sapply(strsplit(as.character(qin.INDEL.edited$tabix_query_PM_10bps),':'), "[", 2))
BED.qin.indel.snv <- rbind.data.frame(BED.qin.indel.snv, cbind.data.frame(CHR,START, END))

BED.qin.indel.snv$KEY <- paste(BED.qin.indel.snv$CHR,BED.qin.indel.snv$START, BED.qin.indel.snv$END, sep= ":")
BED.qin.indel.snv <- BED.qin.indel.snv %>% distinct(KEY, .keep_all = TRUE)
write.table(BED.qin.indel.snv[1:3], "qin_et_al_SNV_start_minus_1_INDEL_plus_min_10.bed", row.names = F, col.names = F, quote = F, sep = "\t")


## Create second search list
# Zhaoming
zhaoming.SNV.edited <- zhaoming.SNV %>% distinct(KEY.varID, .keep_all = TRUE)
zhaoming.SNV.edited <- cbind.data.frame(KEY.varID = zhaoming.SNV.edited$KEY.varID, CHROM = zhaoming.SNV.edited$CHROM, START = zhaoming.SNV.edited$START, END = zhaoming.SNV.edited$END)
zhaoming.SNV.edited$tabix_query <- paste0(zhaoming.SNV.edited$CHROM, ":", zhaoming.SNV.edited$START, "-", zhaoming.SNV.edited$END)
zhaoming.SNV.edited$tabix_query_PM_10bps <- zhaoming.SNV.edited$tabix_query
zhaoming.SNV.edited$TYPE <- "SNV"
zhaoming.INDEL.edited$TYPE <- "INDEL" 
zhaoming.SNV.INDEL.search.list <- rbind.data.frame(zhaoming.INDEL.edited, zhaoming.SNV.edited)
write.table(zhaoming.SNV.INDEL.search.list, "zhaoming_et_al_SNV_INDEL_Search_list_V2.bed", row.names = F, col.names = F, quote = F, sep = "\t")


# Qin
qin.SNV.edited <- qin.SNV %>% distinct(KEY.varID, .keep_all = TRUE)
qin.SNV.edited <- cbind.data.frame(KEY.varID = qin.SNV.edited$KEY.varID, CHROM = qin.SNV.edited$CHROM, START = qin.SNV.edited$START, END = qin.SNV.edited$END)
qin.SNV.edited$tabix_query <- paste0(qin.SNV.edited$CHROM, ":", qin.SNV.edited$START, "-", qin.SNV.edited$END)
qin.SNV.edited$tabix_query_PM_10bps <- qin.SNV.edited$tabix_query
qin.SNV.edited$TYPE <- "SNV"
qin.INDEL.edited$TYPE <- "INDEL" 
qin.SNV.INDEL.search.list <- rbind.data.frame(qin.INDEL.edited, qin.SNV.edited)
write.table(qin.SNV.INDEL.search.list, "qin_et_al_SNV_INDEL_Search_list_V2.bed", row.names = F, col.names = F, quote = F, sep = "\t")


# KEY.varID CHROM     START       END              tabix_query     tabix_query_PM_10bps
# 1       chr2:47803501:C:-  chr2  47803500  47803503   chr2:47803500-47803503   chr2:47803491-47803512
# 2 chr2:68513118:AAACACC:-  chr2  68513117  68513126   chr2:68513117-68513126   chr2:68513108-68513135
# 3   chr14:75047853:-:16bp chr14  75047852  75047870  chr14:75047852-75047870  chr14:75047843-75047879
# 4   chr13:32339288:AAAG:- chr13  32339287  32339293  chr13:32339287-32339293  chr13:32339278-32339302
# 5   chr16:1775961:TCCCC:- chr16   1775960   1775967    chr16:1775960-1775967    chr16:1775951-1775976
# 6      chr2:127272935:T:-  chr2 127272934 127272937 chr2:127272934-127272937 chr2:127272925-127272946

# Now compare how many from zhaoming.INDEL.edited.bed were found in
# zhaoming_et_al_variants_INDEL.bed.out There is a excel file I have created for
# indel:
# //research_jude/rgs01_jude/groups/sapkogrp/projects/SJLIFE_WGS/common/sjlife/MERGED_SJLIFE_1_2/annotation/SNPEFF_ANNOTATION/annotated_indexed_vcf/zhaoming_and_qin_et_al_variant_INDEL_comparison_in_SJLIFE.xlsx


## Now that I have checked the concordance of these variants from Qin and Zhaoming in SJLIFE, I am going to group the variants by genes next
#####################
## Filter variants ##
#####################
CLINVAR <- FINAL.VCF[FINAL.VCF$PRED_TYPE == "Clinvar",]
CLINVAR <- CLINVAR[!duplicated(CLINVAR$KEY.varID),]

MetaSVM <- FINAL.VCF[FINAL.VCF$PRED_TYPE == "MetaSVM",]
MetaSVM <- MetaSVM[!duplicated(MetaSVM$KEY.varID),]

LOF.VCF <- fread("LoF_variants.txt", header = T, sep = "\t")
LOF.VCF$PRED_TYPE <- "LoF"
LOF.VCF$KEY.pos <- paste(LOF.VCF$CHROM, LOF.VCF$POS, sep = ":")        
LOF.VCF$KEY.varID <- paste(LOF.VCF$CHROM, LOF.VCF$POS, LOF.VCF$REF, LOF.VCF$ALT, sep = ":")     

non.intron.LoF <- LOF.VCF[grepl("splice_region_variant&intron_variant", LOF.VCF$`ANN[*].EFFECT`),]
sum(unique(zhaoming.SNV$KEY.varID) %in% non.intron.LoF$KEY.varID)
tt <- non.intron.LoF[non.intron.LoF$KEY.varID %in% unique(zhaoming.SNV$KEY.varID),]

# LOF.VCF <- LOF.VCF[!duplicated(LOF.VCF$KEY.varID),]



#####################################################################
#####################################################################
#####################################################################
#####################################################################

##################################################
## Now, check these in original VCFs from Yadav ##
##################################################
## Zhaoming ##
##############
zhaoming_in_vcf <- read.delim("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/Zhaoming_in_VCF_prior_to_any_QC.txt", sep = "\t", header = T, check.names = F)
zhaoming_in_vcf <- zhaoming_in_vcf[1:11]
head(zhaoming_in_vcf)
# zhaoming_in_vcf <- zhaoming_in_vcf[zhaoming_in_vcf$Match_per_var == "Y",]

zhaoming_in_vcf$MATCHED.IN.CLINVAR <- ""
zhaoming_in_vcf$MATCHED.IN.MetaSVM <- ""
zhaoming_in_vcf$MATCHED.IN.LOF <- ""
zhaoming_in_vcf$MATCHED.IN.ANY.ANNOTATION <- ""
zhaoming_in_vcf$MATCHED.IN.CLINVAR[zhaoming_in_vcf$ID %in% CLINVAR$KEY.varID] <- "Y"
zhaoming_in_vcf$MATCHED.IN.MetaSVM[zhaoming_in_vcf$ID %in% MetaSVM$KEY.varID] <- "Y"
zhaoming_in_vcf$MATCHED.IN.LOF[zhaoming_in_vcf$ID %in% LOF.VCF$KEY.varID] <- "Y"
zhaoming_in_vcf$MATCHED.IN.ANY.ANNOTATION[zhaoming_in_vcf$ID %in% unique(c(CLINVAR$KEY.varID, MetaSVM$KEY.varID, LOF.VCF$KEY.varID))] <- "Y"

# Now if zhaoming_in_vcf$Match_YN is not Y, then make all columns Not Y
zhaoming_in_vcf[!zhaoming_in_vcf$Match_YN == "Y",c("MATCHED.IN.CLINVAR", "MATCHED.IN.MetaSVM", "MATCHED.IN.LOF", "MATCHED.IN.ANY.ANNOTATION")] <- ""
zhaoming_in_vcf$notes <- ""
## Note: Match_per_var is where the variants in VCF match exactly (YES/NO). It is a unique Y/N value per variant whereas Match_YN could have duplicate YES/or NOs
library(dplyr)
library(tidyr)
# Return column names in a column where the value matches a given string
zhaoming_in_vcf <- zhaoming_in_vcf %>% 
  mutate(across(c(MATCHED.IN.CLINVAR, MATCHED.IN.MetaSVM, MATCHED.IN.LOF), ~case_when(. == "Y" ~ cur_column()), .names = 'new_{col}')) %>%
  unite(notes, starts_with('new'), na.rm = TRUE, sep = ';')
zhaoming_in_vcf$notes <- gsub("MATCHED.IN.", "", zhaoming_in_vcf$notes)

zhaoming_in_vcf$notes [zhaoming_in_vcf$QC_Dropped == "QC_Dropped"] <- "QC_Dropped"
# zhaoming_in_vcf$notes [zhaoming_in_vcf$Match_per_var =="Y" & !zhaoming_in_vcf$MATCHED.IN.ANY.ANNOTATION == "Y" ] <- "QC_DROPPED"
zhaoming_in_vcf$notes [is.na(zhaoming_in_vcf$ID)] <- "NOT_CALLED_BY_GATK"

zhaoming_in_vcf$notes [zhaoming_in_vcf$Match_per_var =="N" & !is.na(zhaoming_in_vcf$ID)] <- "Not_Matched"

write.table(zhaoming_in_vcf, "Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/Zhaoming_in_VCF_prior_to_any_QC_Final_list.txt", row.names = F, col.names = T, quote = F, sep = "\t")

#########
## Qin ##
#########
qin_in_vcf <- read.delim("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/qin_in_VCF_prior_to_any_QC.txt", sep = "\t", header = T, check.names = F)
qin_in_vcf <- qin_in_vcf[1:11]
head(qin_in_vcf)
# qin_in_vcf <- qin_in_vcf[qin_in_vcf$Match_per_var == "Y",]

qin_in_vcf$MATCHED.IN.CLINVAR <- ""
qin_in_vcf$MATCHED.IN.MetaSVM <- ""
qin_in_vcf$MATCHED.IN.LOF <- ""
qin_in_vcf$MATCHED.IN.ANY.ANNOTATION <- ""
qin_in_vcf$MATCHED.IN.CLINVAR[qin_in_vcf$ID %in% CLINVAR$KEY.varID] <- "Y"
qin_in_vcf$MATCHED.IN.MetaSVM[qin_in_vcf$ID %in% MetaSVM$KEY.varID] <- "Y"
qin_in_vcf$MATCHED.IN.LOF[qin_in_vcf$ID %in% LOF.VCF$KEY.varID] <- "Y"
qin_in_vcf$MATCHED.IN.ANY.ANNOTATION[qin_in_vcf$ID %in% unique(c(CLINVAR$KEY.varID, MetaSVM$KEY.varID, LOF.VCF$KEY.varID))] <- "Y"

# Now if qin_in_vcf$Match_YN is not Y, then make all columns Not Y
qin_in_vcf[!qin_in_vcf$Match_YN == "Y",c("MATCHED.IN.CLINVAR", "MATCHED.IN.MetaSVM", "MATCHED.IN.LOF", "MATCHED.IN.ANY.ANNOTATION")] <- ""
qin_in_vcf$notes <- ""
## Note: Match_per_var is where the variants in VCF match exactly (YES/NO). It is a unique Y/N value per variant whereas Match_YN could have duplicate YES/or NOs
library(dplyr)
library(tidyr)
# Return column names in a column where the value matches a given string
qin_in_vcf <- qin_in_vcf %>% 
  mutate(across(c(MATCHED.IN.CLINVAR, MATCHED.IN.MetaSVM, MATCHED.IN.LOF), ~case_when(. == "Y" ~ cur_column()), .names = 'new_{col}')) %>%
  unite(notes, starts_with('new'), na.rm = TRUE, sep = ';')
qin_in_vcf$notes <- gsub("MATCHED.IN.", "", qin_in_vcf$notes)

qin_in_vcf$notes [qin_in_vcf$QC_Dropped == "QC_Dropped"] <- "QC_Dropped"
# qin_in_vcf$notes [qin_in_vcf$Match_per_var =="Y" & !qin_in_vcf$MATCHED.IN.ANY.ANNOTATION == "Y" ] <- "QC_DROPPED"
qin_in_vcf$notes [is.na(qin_in_vcf$ID)] <- "NOT_CALLED_BY_GATK"

qin_in_vcf$notes [qin_in_vcf$Match_per_var =="N" & !is.na(qin_in_vcf$ID)] <- "Not_Matched"

write.table(qin_in_vcf, "Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/Qin_in_VCF_prior_to_any_QC_Final_list.txt", row.names = F, col.names = T, quote = F, sep = "\t")




##################################
## How many variants per gene ? ##
##################################
df <- CLINVAR
group_and_concat <- df %>%
  dplyr::select(`ANN[*].GENE`, KEY.varID) %>% 
  dplyr::group_by(`ANN[*].GENE`) %>%
  dplyr::summarise(all_variants_SJLIFE = paste(KEY.varID, collapse = ","))

group_and_concat$counts_var_SJLIFE <- lengths(strsplit(group_and_concat$all_variants_SJLIFE, ","))

qin.genes.vars <- qin.etal.vars
unique(qin.etal.vars$Gene) %in% group_and_concat$`ANN[*].GENE` 

save.image("SNPEFF_clinvar_metaSVM_from_R_filtering_process_PreQC_VCF.RData")
load("SNPEFF_clinvar_metaSVM_from_R_filtering_process_PreQC_VCF.RData")





