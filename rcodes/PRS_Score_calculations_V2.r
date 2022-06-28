########################
## Intermediate edits ##
########################

## Michigan overall; Take POS_START
OVERALL.37 <- fread("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/prs/breast/all_downloads/MichiganWeb/OVERALL_PRSWEB_PHECODE174.1_Onco-iCOGS-Overall-BRCA_PRS-CS_UKB_20200608_WEIGHTS_edited_1.txt", header = T)
OVERALL.38 <- fread("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/prs/breast/all_downloads/MichiganWeb/Michigan_OVERALL_GrCh38_chr_only.bed", header = F)
colnames(OVERALL.38) <- c("CHROM_GRCh38", "POS_START_GRCh38", "POS_END_GRCh38", "KEY")

OVERALL.37$`ChrCHROM:POS`
FINAL.dat <- cbind.data.frame(OVERALL.38[match(OVERALL.37$`ChrCHROM:POS`, OVERALL.38$KEY),], OVERALL.37)

FINAL.dat_Michigan_overall_breast <- cbind.data.frame(CHROM = FINAL.dat$CHROM_GRCh38,POS_GRCh38 = FINAL.dat$POS_START_GRCh38, EA = FINAL.dat$EA, OA = FINAL.dat$OA, Weight = FINAL.dat$WEIGHT)
write.table(FINAL.dat_Michigan_overall_breast, "Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/prs/breast/all_downloads/MichiganWeb/Michigan_OVERALL_GrCh38_edited.txt", row.names = F, col.names = T, quote = F, sep = "\t")



## Michigan ER_NEG; Take POS_START
ER_NEG.37 <- fread("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/prs/breast/all_downloads/MichiganWeb/ER_NEG_PRSWEB_PHECODE174.1_GWAS-Catalog-r2019-05-03-X174.1_PT_UKB_20200608_WEIGHTS_edited_1.txt", header = T)
ER_NEG.38 <- fread("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/prs/breast/all_downloads/MichiganWeb/Michigan_ER_NEG_GrCh38.bed", header = F)
colnames(ER_NEG.38) <- c("CHROM_GRCh38", "POS_START_GRCh38", "POS_END_GRCh38", "KEY")

ER_NEG.37$`ChrCHROM:POS`
FINAL.dat <- cbind.data.frame(ER_NEG.38[match(ER_NEG.37$`ChrCHROM:POS`, ER_NEG.38$KEY),], ER_NEG.37)

FINAL.dat_Michigan_ER_NEG_breast <- cbind.data.frame(CHROM = FINAL.dat$CHROM_GRCh38,POS_GRCh38 = FINAL.dat$POS_START_GRCh38, EA = FINAL.dat$EA, OA = FINAL.dat$OA, Weight = FINAL.dat$WEIGHT)
FINAL.dat_Michigan_ER_NEG_breast[,2][grepl("145790095", FINAL.dat_Michigan_ER_NEG_breast[,2] )] <- "145790097"
write.table(FINAL.dat_Michigan_ER_NEG_breast, "Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/prs/breast/all_downloads/MichiganWeb/Michigan_ER_NEG_GrCh38_edited.txt", row.names = F, col.names = T, quote = F, sep = "\t")


## Michigan ER_POS; Take POS_START
ER_POS.37 <- fread("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/prs/breast/all_downloads/MichiganWeb/ER_POS_PRSWEB_PHECODE174.1_Onco-iCOGS-ER-positive-BRCA_PRS-CS_MGI_20200608_WEIGHTS_edited_1.txt", header = T)
ER_POS.38 <- fread("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/prs/breast/all_downloads/MichiganWeb/Michigan_ER_POS_GrCh38_chr_only.bed", header = F)
colnames(ER_POS.38) <- c("CHROM_GRCh38", "POS_START_GRCh38", "POS_END_GRCh38", "KEY")

ER_POS.37$`ChrCHROM:POS`
FINAL.dat <- cbind.data.frame(ER_POS.38[match(ER_POS.37$`ChrCHROM:POS`, ER_POS.38$KEY),], ER_POS.37)

FINAL.dat_Michigan_ER_POS_breast <- cbind.data.frame(CHROM = FINAL.dat$CHROM_GRCh38,POS_GRCh38 = FINAL.dat$POS_START_GRCh38, EA = FINAL.dat$EA, OA = FINAL.dat$OA, Weight = FINAL.dat$WEIGHT)
write.table(FINAL.dat_Michigan_ER_POS_breast, "Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/prs/breast/all_downloads/MichiganWeb/Michigan_ER_POS_GrCh38_edited.txt", row.names = F, col.names = T, quote = F, sep = "\t")


###################
## Breast Cancer ##
###################
options(scipen = 999)
## Mavaddat 2015
mavaddat_2015_ER_NEG <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/prs/breast/all_downloads/Mavaddat_2015/Mavaddat_2015_ER_NEG_edited.txt", header = T, sep = "\t")
mavaddat_2015_ER_NEG <- cbind.data.frame(CHROM = mavaddat_2015_ER_NEG$CHROM, POS_GRCh38 = mavaddat_2015_ER_NEG$POS_GRCh38, 
                 REF= mavaddat_2015_ER_NEG$other_allele, Effect_allele = mavaddat_2015_ER_NEG$effect_allele, 
                 Effect_size = mavaddat_2015_ER_NEG$effect_weight, TYPE = "Mavaddat_2015_ER_NEG_Breast", 
                 Cancer = "Breast", Significant_YN = "Y")

mavaddat_2015_ER_POS <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/prs/breast/all_downloads/Mavaddat_2015/Mavaddat_2015_ER_POS_edited.txt", header = T, sep = "\t")
mavaddat_2015_ER_POS <- cbind.data.frame(CHROM = mavaddat_2015_ER_POS$CHROM, POS_GRCh38 = mavaddat_2015_ER_POS$Pos_GRCh38, 
                                         REF= mavaddat_2015_ER_POS$other_allele, Effect_allele = mavaddat_2015_ER_POS$effect_allele, 
                                         Effect_size = mavaddat_2015_ER_POS$effect_weight, TYPE = "Mavaddat_2015_ER_POS_Breast", 
                                         Cancer = "Breast", Significant_YN = "Y")

mavaddat_2015_ER_OVERALL <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/prs/breast/all_downloads/Mavaddat_2015/Mavaddat_2015_ER_OVERALL_edited.txt", header = T, sep = "\t")
mavaddat_2015_ER_OVERALL <- cbind.data.frame(CHROM = mavaddat_2015_ER_OVERALL$CHROM, POS_GRCh38 = mavaddat_2015_ER_OVERALL$POS_GRCh38, 
                                         REF= mavaddat_2015_ER_OVERALL$other_allele, Effect_allele = mavaddat_2015_ER_OVERALL$effect_allele, 
                                         Effect_size = mavaddat_2015_ER_OVERALL$effect_weight, TYPE = "Mavaddat_2015_ER_OVERALL_Breast", 
                                         Cancer = "Breast", Significant_YN = "Y")

mavaddat_2015 <- rbind.data.frame(mavaddat_2015_ER_NEG, mavaddat_2015_ER_POS, mavaddat_2015_ER_OVERALL)

## Mavaddat 2019
mavaddat_2019_ER_NEG <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/prs/breast/all_downloads/Mavaddat_2019/Mavaddat_2019_ER_NEG_edited.txt", header = T, sep = "\t")
mavaddat_2019_ER_NEG <- cbind.data.frame(CHROM = mavaddat_2019_ER_NEG$CHROM, POS_GRCh38 = mavaddat_2019_ER_NEG$POS_GRCh38, 
                                         REF= mavaddat_2019_ER_NEG$REF, Effect_allele = mavaddat_2019_ER_NEG$effect_allele, 
                                         Effect_size = mavaddat_2019_ER_NEG$effect_weight, TYPE = "Mavaddat_2019_ER_NEG_Breast", 
                                         Cancer = "Breast", Significant_YN = "Y")

mavaddat_2019_ER_POS <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/prs/breast/all_downloads/Mavaddat_2019/Mavaddat_2019_ER_POS_edited.txt", header = T, sep = "\t")
mavaddat_2019_ER_POS <- cbind.data.frame(CHROM = mavaddat_2019_ER_POS$CHROM, POS_GRCh38 = mavaddat_2019_ER_POS$POS_GRCh38, 
                                         REF= mavaddat_2019_ER_POS$REF, Effect_allele = mavaddat_2019_ER_POS$effect_allele, 
                                         Effect_size = mavaddat_2019_ER_POS$effect_weight, TYPE = "Mavaddat_2019_ER_POS_Breast", 
                                         Cancer = "Breast", Significant_YN = "Y")

mavaddat_2019_ER_OVERALL <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/prs/breast/all_downloads/Mavaddat_2019/Mavaddat_2019_ER_OVERALL_edited.txt", header = T, sep = "\t")
mavaddat_2019_ER_OVERALL <- cbind.data.frame(CHROM = mavaddat_2019_ER_OVERALL$CHROM, POS_GRCh38 = mavaddat_2019_ER_OVERALL$POS_GRCh38, 
                                         REF= mavaddat_2019_ER_OVERALL$REF, Effect_allele = mavaddat_2019_ER_OVERALL$effect_allele, 
                                         Effect_size = mavaddat_2019_ER_OVERALL$effect_weight, TYPE = "Mavaddat_2019_ER_OVERALL_Breast", 
                                         Cancer = "Breast", Significant_YN = "Y")

mavaddat_2019 <- rbind.data.frame(mavaddat_2019_ER_NEG, mavaddat_2019_ER_POS, mavaddat_2019_ER_OVERALL)

## Khera_2018
Khera_2018 <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/prs/breast/all_downloads/Khera_2018/Khera_2018_edited.txt", header = T, sep = "\t")
Khera_2018 <- cbind.data.frame(CHROM = Khera_2018$CHROM, POS_GRCh38 = Khera_2018$POS_GRCh38, 
                                             REF= Khera_2018$other_allele, Effect_allele = Khera_2018$effect_allele, 
                                             Effect_size = Khera_2018$effect_weight, TYPE = "Khera_2018_Breast", 
                                             Cancer = "Breast", Significant_YN = "Y")



## Michigan Web
MichiganWeb_ER_NEG <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/prs/breast/all_downloads/MichiganWeb/Michigan_ER_NEG_GrCh38_edited.txt", header = T, sep = "\t")
MichiganWeb_ER_NEG <- cbind.data.frame(CHROM = MichiganWeb_ER_NEG$CHROM, POS_GRCh38 = MichiganWeb_ER_NEG$POS_GRCh38, 
                               REF= MichiganWeb_ER_NEG$OA, Effect_allele = MichiganWeb_ER_NEG$EA, 
                               Effect_size = MichiganWeb_ER_NEG$Weight, TYPE = "MichiganWeb_ER_NEG_Breast", 
                               Cancer = "Breast", Significant_YN = "Y")

MichiganWeb_ER_POS <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/prs/breast/all_downloads/MichiganWeb/Michigan_ER_POS_GrCh38_edited.txt", header = T, sep = "\t")
MichiganWeb_ER_POS <- cbind.data.frame(CHROM = MichiganWeb_ER_POS$CHROM, POS_GRCh38 = MichiganWeb_ER_POS$POS_GRCh38, 
                                       REF= MichiganWeb_ER_POS$OA, Effect_allele = MichiganWeb_ER_POS$EA, 
                                       Effect_size = MichiganWeb_ER_POS$Weight, TYPE = "MichiganWeb_ER_POS_Breast", 
                                       Cancer = "Breast", Significant_YN = "Y")

MichiganWeb_ER_OVERALL <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/prs/breast/all_downloads/MichiganWeb/Michigan_OVERALL_GrCh38_edited.txt", header = T, sep = "\t")
MichiganWeb_ER_OVERALL <- cbind.data.frame(CHROM = MichiganWeb_ER_OVERALL$CHROM, POS_GRCh38 = MichiganWeb_ER_OVERALL$POS_GRCh38, 
                                       REF= MichiganWeb_ER_OVERALL$OA, Effect_allele = MichiganWeb_ER_OVERALL$EA, 
                                       Effect_size = MichiganWeb_ER_OVERALL$Weight, TYPE = "MichiganWeb_ER_OVERALL_Breast", 
                                       Cancer = "Breast", Significant_YN = "Y")

MichiganWeb <- rbind.data.frame(MichiganWeb_ER_NEG, MichiganWeb_ER_POS, MichiganWeb_ER_OVERALL)


Wang_African <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/prs/breast/all_downloads/African_Hispanic/Wang_African_edited.txt", header = T, sep = "\t")
Wang_African <- cbind.data.frame(CHROM = Wang_African$CHROM, POS_GRCh38 = Wang_African$POS_GRCh38, 
                                           REF= Wang_African$REF, Effect_allele = Wang_African$Effect_allele, 
                                           Effect_size = Wang_African$Effect_size, TYPE = "Wang_African_Breast", 
                                           Cancer = "Breast", Significant_YN = "Y")

Allman_African <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/prs/breast/all_downloads/African_Hispanic/Allman_African_edited.txt", header = T, sep = "\t")
Allman_African <- cbind.data.frame(CHROM = Allman_African$CHROM, POS_GRCh38 = Allman_African$POS_GRCh38, 
                                 REF= Allman_African$REF, Effect_allele = Allman_African$Effect_allele, 
                                 Effect_size = Allman_African$Effect_size, TYPE = "Allman_African_Breast", 
                                 Cancer = "Breast", Significant_YN = "Y")


Allman_Hispanic <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/prs/breast/all_downloads/African_Hispanic/Allman_Hispanic_edited.txt", header = T, sep = "\t")
Allman_Hispanic <- cbind.data.frame(CHROM = Allman_Hispanic$CHROM, POS_GRCh38 = Allman_Hispanic$POS_GRCh38, 
                                   REF= Allman_Hispanic$REF, Effect_allele = Allman_Hispanic$Effect_allele, 
                                   Effect_size = Allman_Hispanic$Effect_size, TYPE = "Allman_Hispanic_Breast", 
                                   Cancer = "Breast", Significant_YN = "Y")


Breast <- rbind.data.frame(mavaddat_2015, mavaddat_2019, Khera_2018, MichiganWeb, Wang_African, Allman_African, Allman_Hispanic)


## Add rest of the cancer datasets
other_cancers <-  read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/prs/Final_list_1.txt", header = T, sep = "\t")

## Swap effect alleles with REF if Effect_size values are negative
all.cancers <- rbind.data.frame(Breast, other_cancers)

###############################
## fix effect size +/- signs ##
###############################
all.cancers.saved <- all.cancers
all.cancers[all.cancers$Effect_size < 0, c("REF", "Effect_allele")] <- all.cancers[all.cancers$Effect_size < 0, c("Effect_allele", "REF")]
all.cancers$Effect_size <- abs(all.cancers$Effect_size)

#############
## Fix Chr ##
#############
table(all.cancers$CHROM)
all.cancers$CHROM <- vapply(strsplit(all.cancers$CHROM, '_'), function(x) 
  paste(x[seq.int(1)], collapse='_'), character(1L))



CHR14 <- all.cancers[grepl("chr14", all.cancers$CHROM),]
######################
## Create bed files ##
######################
all.cancers

all.cancers <- all.cancers[all.cancers$CHROM != "NA",]


all.cancers.INDELS <- all.cancers[nchar(all.cancers$REF) != 1 | nchar(all.cancers$Effect_allele) != 1 ,]
nrow(all.cancers.INDELS)
# 185

# Calculating the base position for END in INDELS. First, I am removing the
# matching bases from REF and ALT and adding the nchar to the START to get the
# END position
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
nrow(all.cancers.SNVS)
# 2241532

PRS.SNVS <- cbind.data.frame(CHROM = all.cancers.SNVS$CHROM, START = all.cancers.SNVS$POS_GRCh38-1, END = all.cancers.SNVS$POS_GRCh38)


PRS.INDELS$KEY <- paste(PRS.INDELS$CHROM, PRS.INDELS$START, PRS.INDELS$END, sep = ":")
PRS.INDELS <- PRS.INDELS[!duplicated(PRS.INDELS$KEY),1:3]

PRS.SNVS$KEY <- paste(PRS.SNVS$CHROM, PRS.SNVS$START, PRS.SNVS$END, sep = ":")
PRS.SNVS <- PRS.SNVS[!duplicated(PRS.SNVS$KEY),1:3]

# all.cancers$KEY <- paste(all.cancers$CHROM, all.cancers$POS_GRCh38, sep = ":")
# all.cancers <- all.cancers[!duplicated(all.cancers$KEY),]
# all.cancers$CHROM <- gsub("chr", "", all.cancers$CHROM)






###################################
## Check Matches with plink file ##
###################################


# all.cancers <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/prs/ALL_Cancers_PRS_data.txt", header = T)
all.cancers$KEY <- paste(all.cancers$CHROM, all.cancers$POS_GRCh38, all.cancers$REF, all.cancers$Effect_allele, sep = ":")
all.cancers$KEY2 <- paste(all.cancers$CHROM, all.cancers$POS_GRCh38, sep = ":")
all.cancers$KEY3 <- paste(all.cancers$CHROM, all.cancers$POS_GRCh38, all.cancers$Effect_allele, all.cancers$REF, sep = ":")
dim(all.cancers) 
# 2241717      11

# Unique SNPIDs
length(unique(all.cancers$KEY))
# 1547024
# Unique Sites
length(unique(all.cancers$KEY2))
# 1121120

# Variants with missing REF
all.cancers[grepl("X",all.cancers$KEY),]

# extract plink subset
ALL.cancers.Keys <- as.data.frame(unique(c(all.cancers$KEY, all.cancers$KEY3)))
new.keys <- as.character(unique(c(all.cancers$KEY, all.cancers$KEY3)))

setwd("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/prs/plink_data/")
## Combine Bim file and PRS file together
bim_file_all <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/prs/plink_data/sjlife_all_PRS.bim", header = F)
bim_file_all$KEY <- paste0("chr", paste(bim_file_all$V1, bim_file_all$V4, bim_file_all$V5, bim_file_all$V6, sep = ":"))
bim_file_all$KEY2 <- paste0("chr", paste(bim_file_all$V1, bim_file_all$V4, sep = ":"))


# Unique SNPIDs
sum(unique(as.character(ALL.cancers.Keys$`unique(c(all.cancers$KEY, all.cancers$KEY3))`)) %in% as.character(bim_file_all$KEY)) 
# 1119974


## Extract the difference from the old PRS and run the bcftools on the difference again
OLD.all.cancers <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/prs/OLD_ALL_Cancers_PRS_data.txt", header = T)
OLD.all.cancers$KEY <- paste0("chr", paste(OLD.all.cancers$CHROM, OLD.all.cancers$POS_GRCh38, OLD.all.cancers$REF, OLD.all.cancers$Effect_allele, sep = ":"))
OLD.all.cancers$KEY3 <- paste0("chr", paste(OLD.all.cancers$CHROM, OLD.all.cancers$POS_GRCh38, OLD.all.cancers$Effect_allele, OLD.all.cancers$REF, sep = ":"))

old.keys <- as.character(c(OLD.all.cancers$KEY, OLD.all.cancers$KEY3))

sum(!new.keys %in% old.keys)
# 351
new.keys <- new.keys[!new.keys %in% old.keys]

new.set1 <- all.cancers[all.cancers$KEY %in% new.keys,]
new.set2 <- all.cancers[all.cancers$KEY3 %in% new.keys,]

new.set.all <- rbind.data.frame(new.set1, new.set2)
library(dplyr)
new.set.all <- distinct(new.set.all)

write.table(ALL.cancers.Keys, "PRS_all_cancers_vars.txt", col.names = F, row.names = F, quote = F)

write.table(PRS.INDELS, "Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/prs/Indels_PRS.bed", row.names = F, col.names = F, quote = F, sep = "\t")
write.table(PRS.SNVS, "Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/prs/SNVs_PRS.bed", row.names = F, col.names = F, quote = F, sep = "\t")

all.cancers1 <- all.cancers
all.cancers1$CHROM <- gsub("chr", "", all.cancers1$CHROM)

all.cancers2 <- all.cancers1

# Fix variants with X in reference
# i=1; chr2:145366886:X:A >> plink >> chr2:145366886:C:A
i=1
all.cancers2 [grepl("ALL_Vijayakrishnan", all.cancers2$TYPE),][i,] <- 
  apply(all.cancers2 [grepl("ALL_Vijayakrishnan", all.cancers2$TYPE),][i,], 2, function(x) gsub("X","C",x))  

i=2
all.cancers2 [grepl("ALL_Vijayakrishnan", all.cancers2$TYPE),][i,] <- 
  apply(all.cancers2 [grepl("ALL_Vijayakrishnan", all.cancers2$TYPE),][i,], 2, function(x) gsub("X","C",x))  

i=3
all.cancers2 [grepl("ALL_Vijayakrishnan", all.cancers2$TYPE),][i,] <- 
  apply(all.cancers2 [grepl("ALL_Vijayakrishnan", all.cancers2$TYPE),][i,], 2, function(x) gsub("X","T",x))  

i=4
all.cancers2 [grepl("ALL_Vijayakrishnan", all.cancers2$TYPE),][i,] <- 
  apply(all.cancers2 [grepl("ALL_Vijayakrishnan", all.cancers2$TYPE),][i,], 2, function(x) gsub("X","T",x))  

i=5
all.cancers2 [grepl("ALL_Vijayakrishnan", all.cancers2$TYPE),][i,] <- 
  apply(all.cancers2 [grepl("ALL_Vijayakrishnan", all.cancers2$TYPE),][i,], 2, function(x) gsub("X","A",x))  

i=6
all.cancers2 [grepl("ALL_Vijayakrishnan", all.cancers2$TYPE),][i,] <- 
  apply(all.cancers2 [grepl("ALL_Vijayakrishnan", all.cancers2$TYPE),][i,], 2, function(x) gsub("X","T",x))  

i=7
all.cancers2 [grepl("ALL_Vijayakrishnan", all.cancers2$TYPE),][i,] <- 
  apply(all.cancers2 [grepl("ALL_Vijayakrishnan", all.cancers2$TYPE),][i,], 2, function(x) gsub("X","G",x))  

i=8
all.cancers2 [grepl("ALL_Vijayakrishnan", all.cancers2$TYPE),][i,] <- 
  apply(all.cancers2 [grepl("ALL_Vijayakrishnan", all.cancers2$TYPE),][i,], 2, function(x) gsub("X","T",x))  

i=9
all.cancers2 [grepl("ALL_Vijayakrishnan", all.cancers2$TYPE),][i,] <- 
  apply(all.cancers2 [grepl("ALL_Vijayakrishnan", all.cancers2$TYPE),][i,], 2, function(x) gsub("X","C",x))  

i=10
all.cancers2 [grepl("ALL_Vijayakrishnan", all.cancers2$TYPE),][i,] <- 
  apply(all.cancers2 [grepl("ALL_Vijayakrishnan", all.cancers2$TYPE),][i,], 2, function(x) gsub("X","T",x))  

i=11
all.cancers2 [grepl("ALL_Vijayakrishnan", all.cancers2$TYPE),][i,] <- 
  apply(all.cancers2 [grepl("ALL_Vijayakrishnan", all.cancers2$TYPE),][i,], 2, function(x) gsub("X","T",x))  

i=12
all.cancers2 [grepl("ALL_Vijayakrishnan", all.cancers2$TYPE),][i,] <- 
  apply(all.cancers2 [grepl("ALL_Vijayakrishnan", all.cancers2$TYPE),][i,], 2, function(x) gsub("X","A",x))  

i=13
all.cancers2 [grepl("ALL_Vijayakrishnan", all.cancers2$TYPE),][i,] <- 
  apply(all.cancers2 [grepl("ALL_Vijayakrishnan", all.cancers2$TYPE),][i,], 2, function(x) gsub("X","G",x))  

i=14
all.cancers2 [grepl("ALL_Vijayakrishnan", all.cancers2$TYPE),][i,] <- 
  apply(all.cancers2 [grepl("ALL_Vijayakrishnan", all.cancers2$TYPE),][i,], 2, function(x) gsub("X","C",x))  

i=15
all.cancers2 [grepl("ALL_Vijayakrishnan", all.cancers2$TYPE),][i,] <- 
  apply(all.cancers2 [grepl("ALL_Vijayakrishnan", all.cancers2$TYPE),][i,], 2, function(x) gsub("X","A",x))  

i=1
all.cancers2 [grepl("Meningioma_Claus", all.cancers2$TYPE),][i,] <- 
  apply(all.cancers2 [grepl("Meningioma_Claus", all.cancers2$TYPE),][i,], 2, function(x) gsub("X","A",x)) 




# write.table(all.cancers1[1:8], "Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/prs/ALL_Cancers_PRS_data.txt", row.names = F, col.names = T, quote = F, sep = "\t")
write.table(all.cancers2[1:8], "Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/prs/ALL_Cancers_PRS_data.txt", row.names = F, col.names = T, quote = F, sep = "\t")



###################################################################
## create bed files for the new variants missed in first attempt ##
###################################################################
new.all.cancers.INDELS <- new.set.all[nchar(new.set.all$REF) != 1 | nchar(new.set.all$Effect_allele) != 1 ,]
nrow(new.all.cancers.INDELS)
# 0

new.all.cancers.SNVS <- new.set.all[!(nchar(new.set.all$REF) != 1 | nchar(new.set.all$Effect_allele) != 1),]
nrow(new.all.cancers.SNVS)
# 315

new.PRS.SNVS <- cbind.data.frame(CHROM = new.all.cancers.SNVS$CHROM, START = new.all.cancers.SNVS$POS_GRCh38-1, END = new.all.cancers.SNVS$POS_GRCh38)


new.PRS.SNVS$KEY <- paste(new.PRS.SNVS$CHROM, new.PRS.SNVS$START, new.PRS.SNVS$END, sep = ":")
new.PRS.SNVS <- new.PRS.SNVS[!duplicated(new.PRS.SNVS$KEY),1:3]

write.table(new.PRS.SNVS, "Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/prs/extract_batch2_SNVs_PRS.bed", row.names = F, col.names = F, quote = F, sep = "\t")


#################################################################################################
# Var1    Freq  Plink
# 1                    ALL_Vijayakrishnan      15 
# 2                 Allman_African_Breast      75 
# 3                Allman_Hispanic_Breast      71 
# 4                     Khera_2018_Breast      77 
# 5           Mavaddat_2015_ER_NEG_Breast      77 
# 6       Mavaddat_2015_ER_OVERALL_Breast      77 
# 7           Mavaddat_2015_ER_POS_Breast      77 
# 8           Mavaddat_2019_ER_NEG_Breast     313 
# 9       Mavaddat_2019_ER_OVERALL_Breast     313 
# 10          Mavaddat_2019_ER_POS_Breast     313 
# 11                     Meningioma_Claus       1 
# 12                   Meningioma_Dobbins       2 
# 13            MichiganWeb_ER_NEG_Breast      79 
# 14        MichiganWeb_ER_OVERALL_Breast 1120348 
# 15            MichiganWeb_ER_POS_Breast 1119079 
# 16            Pleiotropy_Bi_directional      21 
# 17             Pleiotropy_Meta_analysis      21 
# 18                Pleiotropy_One_cohort       9 
# 19           Pleiotropy_One_directional     137 
# 20                    Pleiotropy_PRSWEB     179 
# 21 Pleiotropy_Replication_prior_studies     308 
# 22                     Sarcoma_Machiela       6 
# 23                  Wang_African_Breast      98 

meningioma <- all.cancers [grepl("Meningioma", all.cancers$TYPE ),]
ALL <- all.cancers [grepl("Vijay", all.cancers$TYPE ),]
Allman_African <- all.cancers [grepl("Allman_African", all.cancers$TYPE ),]


# find that are missing
# first remove those that are present in final.set from all.cancers
all.cancers.not.matched <- all.cancers [!all.cancers$KEY %in% final.set$SNPId,]
all.cancers.not.matched <- all.cancers.not.matched [!all.cancers.not.matched$KEY3 %in% final.set$SNPId,]

head(plink.final.bim)
plink.final.bim$KEY2 <- paste0("chr", paste(plink.final.bim$V1, plink.final.bim$V4, sep = ":"))

sum(all.cancers.not.matched$KEY2 %in% plink.final.bim$KEY2)
# 183
plink.final.bim.matched2 <- plink.final.bim[plink.final.bim$KEY2 %in% all.cancers.not.matched$KEY2,]
plink.final.bim.matched2

TT <- merge(plink.final.bim.matched2, all.cancers.not.matched, by ="KEY2",  all = T)

plink.final.bim.not.found <- read.table("sjlife_found_v2.bim")

sum(all.cancers.not.matched$KEY %in% plink.final.bim.not.found$V2)
sum(all.cancers.not.matched$KEY3 %in% plink.final.bim.not.found$V2)

# #######################
# ## Harmonize alleles ##
# #######################
# 
# Process the variants to match their alleles

# # Read the input file including variants with no direct match of their alleles
# args = commandArgs(trailingOnly = TRUE)
# dat = read.table(args[1], header = FALSE, stringsAsFactors = FALSE)
# # dat = read.table("Z:/ResearchHome/ClusterHome/ysapkota/Work/CAD_PRS/y", header = FALSE, stringsAsFactors = FALSE)

dat <- TT[-c(1,8,9)]
dat <-  dat[!is.na(dat$V2),]
colnames(dat) <- paste0("V",1:14)
dat <- dat [! (grepl("DEL", dat$V5) | grepl("DEL", dat$V6)),]

## Function to flip alleles
flip_alleles = function(x){
  if(x=="A"){
    y="T"
  } else if (x=="T"){
    y="A"
  } else if (x=="C"){
    y="G"
  } else if (x=="G"){
    y="C"
  }
  return(y)
}

# Process each variant to check the alleles
dat.out = NULL
for (i in 1:nrow(dat)){
  print(paste0("Doing row ", i ))
  chr=dat$V1[i]
  pos = dat$V4[i]
  variant = dat$V2[i]
  khera_a1 = dat$V7[i]
  khera_a2 = dat$V8[i]
  khera_weight = dat$V9[i]
  wgs_a1 = dat$V5[i]
  wgs_a2 = dat$V6[i]
  # First find out if the wgs_a1 have more than one character
  if(nchar(wgs_a1)>1){
    wgs_a1_first = substr(wgs_a1, 1, 1)
    wgs_a1_last = substr(wgs_a1, nchar(wgs_a1), nchar(wgs_a1))
    wgs_a1_changed = ifelse(wgs_a1_first==wgs_a2, wgs_a1_last, wgs_a1_first)
  } else {
    wgs_a1_changed = wgs_a1
  }
  # Then do the same for wgs_a2 allele
  if (nchar(wgs_a2)>1){
    wgs_a2_first = substr(wgs_a2, 1, 1)
    wgs_a2_last = substr(wgs_a2, nchar(wgs_a2), nchar(wgs_a2))
    wgs_a2_changed = ifelse(wgs_a2_first==wgs_a1, wgs_a2_last, wgs_a2_first)
  } else {
    wgs_a2_changed = wgs_a2
  }
  # Now check if the changed wgs alleles match with those from Khera et al
  # No alleles flipped
  if ((wgs_a1_changed == khera_a1 & wgs_a2_changed == khera_a2) | (wgs_a1_changed == khera_a2 & wgs_a2_changed == khera_a1)){
    wgs_a1_new = wgs_a1_changed; wgs_a2_new = wgs_a2_changed; match=1
  } else if((flip_alleles(wgs_a1_changed) == khera_a1 & wgs_a2_changed == khera_a2) | (flip_alleles(wgs_a1_changed) == khera_a2 & wgs_a2_changed == khera_a1)) { # only a1 flipped
    wgs_a1_new = flip_alleles(wgs_a1_changed); wgs_a2_new = wgs_a2_changed; match=1
  } else if ((wgs_a1_changed == khera_a1 & flip_alleles(wgs_a2_changed) == khera_a2) | (wgs_a1_changed == khera_a2 & flip_alleles(wgs_a2_changed) == khera_a1)) { # only a2 flipped
    wgs_a1_new = wgs_a1_changed; wgs_a2_new = flip_alleles(wgs_a2_changed); match=1
  } else if ((flip_alleles(wgs_a1_changed) == khera_a1 & flip_alleles(wgs_a2_changed) == khera_a2) | (flip_alleles(wgs_a1_changed) == khera_a2 & flip_alleles(wgs_a2_changed) == khera_a1)){ # both a1 and a2 flipped
    wgs_a1_new = flip_alleles(wgs_a1_changed); wgs_a2_new = flip_alleles(wgs_a2_changed); match=1
  } else {
    wgs_a1_new = wgs_a1_changed; wgs_a2_new = wgs_a2_changed; match=0
  }
  dat.out = rbind(dat.out, data.frame(chr, pos, variant, khera_a1, khera_a2, wgs_a1, wgs_a2, wgs_a1_new, wgs_a2_new, match))
}

# Write data to disc
# write.table(dat.out, "PRS_SNVs_alleles_harmonized", row.names = FALSE, quote = FALSE)
save.image("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/prs/PRS_vars_search.RData")
# load("PRS_SNVs_alleles_harmonized.RData")
