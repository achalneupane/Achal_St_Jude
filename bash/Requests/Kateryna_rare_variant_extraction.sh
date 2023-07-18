head -1 MERGED.SJLIFE.1.2.GATKv3.4.VQSR.chr10.preQC_biallelic_renamed_ID_edited.vcf-annot-snpeff-dbnsfp-ExAC.0.3-clinvar.GRCh38.vcf.dbSNP155-FIELDS-simple.txt > LoF_variants.txt
cat $(ls *.dbSNP155-FIELDS-simple.txt| sort -V)| egrep 'frameshift_variant|start_lost|stop_gained|splice|splice_region_variant' >> LoF_variants.txt
awk '{ print $3 }' LoF_variants.txt| cut -d";" -f1 > LoF_variants_ID.txt



## R code
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
# LoF.unique <- LoF[!duplicated(LoF$KEY),]
nrow(LoF.unique)
# 316386

splice.region <- LoF.unique[grepl("splice_region_variant", LoF.unique$`ANN[*].EFFECT`),]
dim(splice.region)
# 242541     25
non.splice.region <- LoF.unique[!grepl("splice_region_variant", LoF.unique$`ANN[*].EFFECT`),]
dim(non.splice.region)
# 73845    25

NO.splice.region <- splice.region[grepl("frameshift|stop|donor|acceptor", splice.region$`ANN[*].EFFECT`),]
dim(NO.splice.region)
# 4421   25

LOF.keep <- rbind.data.frame(non.splice.region, NO.splice.region)
dim(LOF.keep)

Predicted.vars.inVCF.Unique <-  rbind.data.frame(CLINVAR.unique, LOF.keep, MetaSVM.unique)

####################################
## Revel is in annovar annotation ##
####################################
REVEL > 0.5 is considered pathogenic