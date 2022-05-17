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
SNPEFF.VCF <- {}
for( i in 1:length(chromosomes)){

print(paste0("Doing chromosome ", chromosomes[i]))
    
VCF <- fread(paste0("MERGED.SJLIFE.1.2.GATKv3.4.VQSR.chr", chromosomes[i], ".PASS.decomposed.vcf-annot-snpeff-dbnsfp-ExAC.0.3-clinvar.GRCh38.vcf.dbSNP155-FIELDS-simple.txt"))

#####################################
## Clinvar based L/PL VCF variants ##
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


WANTED.types <- c("^Pathogenic/Likely_pathogenic$|^Likely_pathogenic$|^Pathogenic/Likely_pathogenic$|^Pathogenic$|Pathogenic\\|_risk_factor")
sum(grepl(WANTED.types, VCF$CLNSIG, ignore.case = T))
# 3505

VCF <- VCF[grepl(WANTED.types, VCF$CLNSIG, ignore.case = T),]
SNPEFF.VCF <- rbind.data.frame(SNPEFF.VCF, VCF)
}
as.data.frame(table(VCF$`ANN[*].EFFECT`))
                      