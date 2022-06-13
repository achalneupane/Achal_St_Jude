## Michigan overall; Take POS_START
OVERALL.37 <- fread("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/prs/breast/all_downloads/MichiganWeb/OVERALL_PRSWEB_PHECODE174.1_Onco-iCOGS-Overall-BRCA_PRS-CS_UKB_20200608_WEIGHTS_edited_1.txt", header = T)
OVERALL.38 <- fread("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/prs/breast/all_downloads/MichiganWeb/Michigan_OVERALL_GrCh38.bed", header = F)
colnames(OVERALL.38) <- c("CHROM_GRCh38", "POS_START_GRCh38", "POS_END_GRCh38", "KEY")

OVERALL.37$`ChrCHROM:POS`
FINAL.dat <- cbind.data.frame(OVERALL.38[match(OVERALL.37$`ChrCHROM:POS`, OVERALL.38$KEY),], OVERALL.37)

FINAL.dat <- cbind(CHROM = FINAL.dat$CHROM_GRCh38,POS_GRCh38 = FINAL.dat$POS_START_GRCh38, EA = FINAL.dat$EA, OA = FINAL.dat$OA, Weight = FINAL.dat$WEIGHT)
write.table(FINAL.dat, "Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/prs/breast/all_downloads/MichiganWeb/Michigan_OVERALL_GrCh38_edited.txt", row.names = F, col.names = T, quote = F, sep = "\t")



## Michigan ER_NEG; Take POS_START
ER_NEG.37 <- fread("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/prs/breast/all_downloads/MichiganWeb/ER_NEG_PRSWEB_PHECODE174.1_GWAS-Catalog-r2019-05-03-X174.1_PT_UKB_20200608_WEIGHTS_edited_1.txt", header = T)
ER_NEG.38 <- fread("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/prs/breast/all_downloads/MichiganWeb/Michigan_ER_NEG_GrCh38.bed", header = F)
colnames(ER_NEG.38) <- c("CHROM_GRCh38", "POS_START_GRCh38", "POS_END_GRCh38", "KEY")

ER_NEG.37$`ChrCHROM:POS`
FINAL.dat <- cbind.data.frame(ER_NEG.38[match(ER_NEG.37$`ChrCHROM:POS`, ER_NEG.38$KEY),], ER_NEG.37)

FINAL.dat <- cbind(CHROM = FINAL.dat$CHROM_GRCh38,POS_GRCh38 = FINAL.dat$POS_START_GRCh38, EA = FINAL.dat$EA, OA = FINAL.dat$OA, Weight = FINAL.dat$WEIGHT)
write.table(FINAL.dat, "Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/prs/breast/all_downloads/MichiganWeb/Michigan_ER_NEG_GrCh38_edited.txt", row.names = F, col.names = T, quote = F, sep = "\t")


## Michigan ER_POS; Take POS_START
ER_POS.37 <- fread("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/prs/breast/all_downloads/MichiganWeb/ER_POS_PRSWEB_PHECODE174.1_Onco-iCOGS-ER-positive-BRCA_PRS-CS_MGI_20200608_WEIGHTS_edited_1.txt", header = T)
ER_POS.38 <- fread("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/prs/breast/all_downloads/MichiganWeb/Michigan_ER_POS_GrCh38.bed", header = F)
colnames(ER_POS.38) <- c("CHROM_GRCh38", "POS_START_GRCh38", "POS_END_GRCh38", "KEY")

ER_POS.37$`ChrCHROM:POS`
FINAL.dat <- cbind.data.frame(ER_POS.38[match(ER_POS.37$`ChrCHROM:POS`, ER_POS.38$KEY),], ER_POS.37)

FINAL.dat <- cbind(CHROM = FINAL.dat$CHROM_GRCh38,POS_GRCh38 = FINAL.dat$POS_START_GRCh38, EA = FINAL.dat$EA, OA = FINAL.dat$OA, Weight = FINAL.dat$WEIGHT)
write.table(FINAL.dat, "Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/prs/breast/all_downloads/MichiganWeb/Michigan_ER_POS_GrCh38_edited.txt", row.names = F, col.names = T, quote = F, sep = "\t")

###################
## Breast Cancer ##
###################

## Mavaddat 2015
mavaddat_2015_ER_NEG <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/prs/breast/all_downloads/Mavaddat_2015/Mavaddat_2015_ER_NEG_edited.txt", header = T, sep = "\t")
mavaddat_2015_ER_NEG <- cbind.data.frame(CHROM = mavaddat_2015_ER_NEG$CHROM, POS_GRCh38 = mavaddat_2015_ER_NEG$POS_GRCh38, 
                 REF= mavaddat_2015_ER_NEG$other_allele, Effect_allele = mavaddat_2015_ER_NEG$effect_allele, 
                 Effect_size = mavaddat_2015_ER_NEG$effect_weight, TYPE = "Breast_Mavaddat_2015_ER_NEG", 
                 Cancer = "Breast", Significant_YN = "Y")

mavaddat_2015_ER_POS <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/prs/breast/all_downloads/Mavaddat_2015/Mavaddat_2015_ER_POS_edited.txt", header = T, sep = "\t")
mavaddat_2015_ER_POS <- cbind.data.frame(CHROM = mavaddat_2015_ER_POS$CHROM, POS_GRCh38 = mavaddat_2015_ER_POS$Pos_GRCh38, 
                                         REF= mavaddat_2015_ER_POS$other_allele, Effect_allele = mavaddat_2015_ER_POS$effect_allele, 
                                         Effect_size = mavaddat_2015_ER_POS$effect_weight, TYPE = "Breast_Mavaddat_2015_ER_POS", 
                                         Cancer = "Breast", Significant_YN = "Y")

mavaddat_2015_ER_OVERALL <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/prs/breast/all_downloads/Mavaddat_2015/Mavaddat_2015_ER_OVERALL_edited.txt", header = T, sep = "\t")
mavaddat_2015_ER_OVERALL <- cbind.data.frame(CHROM = mavaddat_2015_ER_OVERALL$CHROM, POS_GRCh38 = mavaddat_2015_ER_OVERALL$POS_GRCh38, 
                                         REF= mavaddat_2015_ER_OVERALL$other_allele, Effect_allele = mavaddat_2015_ER_OVERALL$effect_allele, 
                                         Effect_size = mavaddat_2015_ER_OVERALL$effect_weight, TYPE = "Breast_Mavaddat_2015_ER_OVERALL", 
                                         Cancer = "Breast", Significant_YN = "Y")

mavaddat_2015 <- rbind.data.frame(mavaddat_2015_ER_NEG, mavaddat_2015_ER_POS, mavaddat_2015_ER_OVERALL)

## Mavaddat 2019
mavaddat_2019_ER_NEG <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/prs/breast/all_downloads/Mavaddat_2019/Mavaddat_2019_ER_NEG_edited.txt", header = T, sep = "\t")
mavaddat_2019_ER_NEG <- cbind.data.frame(CHROM = mavaddat_2019_ER_NEG$CHROM, POS_GRCh38 = mavaddat_2019_ER_NEG$POS_GRCh38, 
                                         REF= mavaddat_2019_ER_NEG$REF, Effect_allele = mavaddat_2019_ER_NEG$effect_allele, 
                                         Effect_size = mavaddat_2019_ER_NEG$effect_weight, TYPE = "Breast_Mavaddat_2019_ER_NEG", 
                                         Cancer = "Breast", Significant_YN = "Y")

mavaddat_2019_ER_POS <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/prs/breast/all_downloads/Mavaddat_2019/Mavaddat_2019_ER_POS_edited.txt", header = T, sep = "\t")
mavaddat_2019_ER_POS <- cbind.data.frame(CHROM = mavaddat_2019_ER_POS$CHROM, POS_GRCh38 = mavaddat_2019_ER_POS$POS_GRCh38, 
                                         REF= mavaddat_2019_ER_POS$REF, Effect_allele = mavaddat_2019_ER_POS$effect_allele, 
                                         Effect_size = mavaddat_2019_ER_POS$effect_weight, TYPE = "Breast_Mavaddat_2019_ER_POS", 
                                         Cancer = "Breast", Significant_YN = "Y")

mavaddat_2019_ER_OVERALL <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/prs/breast/all_downloads/Mavaddat_2019/Mavaddat_2019_ER_OVERALL_edited.txt", header = T, sep = "\t")
mavaddat_2019_ER_OVERALL <- cbind.data.frame(CHROM = mavaddat_2019_ER_OVERALL$CHROM, POS_GRCh38 = mavaddat_2019_ER_OVERALL$POS_GRCh38, 
                                         REF= mavaddat_2019_ER_OVERALL$REF, Effect_allele = mavaddat_2019_ER_OVERALL$effect_allele, 
                                         Effect_size = mavaddat_2019_ER_OVERALL$effect_weight, TYPE = "Breast_Mavaddat_2019_ER_OVERALL", 
                                         Cancer = "Breast", Significant_YN = "Y")

mavaddat_2019 <- rbind.data.frame(mavaddat_2019_ER_NEG, mavaddat_2019_ER_POS, mavaddat_2019_ER_OVERALL)

## Khera_2018
Khera_2018 <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/prs/breast/all_downloads/Khera_2018/Khera_2018_edited.txt", header = T, sep = "\t")
Khera_2018 <- cbind.data.frame(CHROM = Khera_2018$CHROM, POS_GRCh38 = Khera_2018$POS_GRCh38, 
                                             REF= Khera_2018$other_allele, Effect_allele = Khera_2018$effect_allele, 
                                             Effect_size = Khera_2018$effect_weight, TYPE = "Khera_2018", 
                                             Cancer = "Breast", Significant_YN = "Y")



## Michigan Web
MichiganWeb_ER_NEG <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/prs/breast/all_downloads/MichiganWeb/Michigan_ER_NEG_GrCh38_edited.txt", header = T, sep = "\t")
MichiganWeb_ER_NEG <- cbind.data.frame(CHROM = MichiganWeb_ER_NEG$CHROM, POS_GRCh38 = MichiganWeb_ER_NEG$POS_GRCh38, 
                               REF= MichiganWeb_ER_NEG$OA, Effect_allele = MichiganWeb_ER_NEG$EA, 
                               Effect_size = MichiganWeb_ER_NEG$Weight, TYPE = "MichiganWeb_ER_NEG", 
                               Cancer = "Breast", Significant_YN = "Y")

MichiganWeb_ER_POS <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/prs/breast/all_downloads/MichiganWeb/Michigan_ER_POS_GrCh38_edited.txt", header = T, sep = "\t")
MichiganWeb_ER_POS <- cbind.data.frame(CHROM = MichiganWeb_ER_POS$CHROM, POS_GRCh38 = MichiganWeb_ER_POS$POS_GRCh38, 
                                       REF= MichiganWeb_ER_POS$OA, Effect_allele = MichiganWeb_ER_POS$EA, 
                                       Effect_size = MichiganWeb_ER_POS$Weight, TYPE = "MichiganWeb_ER_POS", 
                                       Cancer = "Breast", Significant_YN = "Y")

MichiganWeb_ER_OVERALL <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/prs/breast/all_downloads/MichiganWeb/Michigan_OVERALL_GrCh38_edited.txt", header = T, sep = "\t")
MichiganWeb_ER_OVERALL <- cbind.data.frame(CHROM = MichiganWeb_ER_OVERALL$CHROM, POS_GRCh38 = MichiganWeb_ER_OVERALL$POS_GRCh38, 
                                       REF= MichiganWeb_ER_OVERALL$OA, Effect_allele = MichiganWeb_ER_OVERALL$EA, 
                                       Effect_size = MichiganWeb_ER_OVERALL$Weight, TYPE = "MichiganWeb_ER_OVERALL", 
                                       Cancer = "Breast", Significant_YN = "Y")

MichiganWeb <- rbind.data.frame(MichiganWeb_ER_NEG, MichiganWeb_ER_POS, MichiganWeb_ER_OVERALL)

Breast <- rbind.data.frame(mavaddat_2015, mavaddat_2019, Khera_2018, MichiganWeb)


## Add rest of the cancer datasets
other_cancers <-  read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/prs/Final_list_1.txt", header = T, sep = "\t")

## Swap effect alleles with REF if Effect_size values are negative
all.cancers <- rbind.data.frame(Breast, other_cancers)
all.cancers[all.cancers$Effect_size < 0, c("REF", "Effect_allele")] <- all.cancers[all.cancers$Effect_size < 0, c("Effect_allele", "REF")]
all.cancers$Effect_size <- abs(all.cancers$Effect_size)

######################
## Create bed files ##
######################
all.cancers

all.cancers.INDELS <- all.cancers[nchar(all.cancers$REF) != 1 | nchar(all.cancers$Effect_allele) != 1 ,]
dim(all.cancers.INDELS)
library(stringr) # for str_remove function
fun <- function(a, b){
  a1 <- substr(a,1,1)
  b1 <- substr(b, 1, 1)
  d <- asplit(cbind(a, b), 1)
  ifelse(a1==b1, Recall(str_remove(a,a1), str_remove(b, b1)), d)
}

all.cancers.INDELS[c('REF', 'Effect_allele')] <-  do.call(rbind, fun(all.cancers.INDELS$REF, all.cancers.INDELS$Effect_allele))
all.cancers.INDELS

PRS.INDELS <- cbind.data.frame(CHROM = all.cancers.INDELS$CHROM, START = all.cancers.INDELS$POS_GRCh38, END = all.cancers.INDELS$POS_GRCh38 + nchar(all.cancers.INDELS$REF) + nchar(all.cancers.INDELS$Effect_allele))


all.cancers.SNVS <- all.cancers[!(nchar(all.cancers$REF) != 1 | nchar(all.cancers$Effect_allele) != 1),]
dim(all.cancers.SNVS)

PRS.SNVS <- cbind.data.frame(CHROM = all.cancers.SNVS$CHROM, START = all.cancers.SNVS$POS_GRCh38-1, END = all.cancers.SNVS$POS_GRCh38)


PRS.INDELS$KEY <- paste(PRS.INDELS$CHROM, PRS.INDELS$START, PRS.INDELS$END, sep = ":")
PRS.INDELS <- PRS.INDELS[!duplicated(PRS.INDELS$KEY),1:3]

PRS.SNVS$KEY <- paste(PRS.SNVS$CHROM, PRS.SNVS$START, PRS.SNVS$END, sep = ":")
PRS.SNVS <- PRS.SNVS[!duplicated(PRS.SNVS$KEY),1:3]

# all.cancers$KEY <- paste(all.cancers$CHROM, all.cancers$POS_GRCh38, sep = ":")
# all.cancers <- all.cancers[!duplicated(all.cancers$KEY),]
all.cancers$CHROM <- gsub("chr", "", all.cancers$CHROM)

all.cancers <- all.cancers[!is.na(all.cancers$CHROM),]

write.table(PRS.INDELS, "Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/prs/Indels_PRS.bed", row.names = F, col.names = F, quote = F, sep = "\t")
write.table(PRS.SNVS, "Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/prs/SNVs_PRS.bed", row.names = F, col.names = F, quote = F, sep = "\t")
write.table(all.cancers, "Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/prs/ALL_Cancers_PRS_data.txt", row.names = F, col.names = T, quote = F, sep = "\t")



setwd("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/prs/plink_data/")
## Combine Bim file and PRS file together
bim_file <- read.table("PRS_chr20.bim", header = T)

