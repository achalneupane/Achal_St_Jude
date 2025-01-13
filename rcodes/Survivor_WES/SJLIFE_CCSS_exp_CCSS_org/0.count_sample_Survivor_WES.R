rm(list=ls())
library(haven)
library(readxl)

#######################
## A. working on WES ##
#######################
# https://wiki.stjude.org/display/CAB/Genetic+Ancestry+Estimation+by+PCA
setwd("Z:/ResearchHome/Groups/sapkogrp/projects/CAB/common/SJLIFE_CCSS_WES_101724/GermlineQC/")


## 1. Check how many samples are available 
# SJLIFE
pop <- read_sas("Z:/SJShare/SJCOMMON/ECC/SJLife/SJLIFE Data Freeze/2 Final Data SJLIFE/20200430/Clinical Data/demographics.sas7bdat")
pop.survivor <- pop[pop$studypop == "Survivor",]
pop.survivor.control <- pop[grepl("Control", pop$studypop),]

## Read the specified sheet from the Excel file shared by Yadav. It contains the list of all samples that have been or will be sequenced.
sequencing.Record <- read_excel("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/Survivor_WES_QC/sample_mapping_files/SJLIFE_WGS_and_CCSS_SNP_WES_WGS_samples_and_overlap.xlsx", sheet = "SJLIFEWGS_4481", col_names = F)
table(sequencing.Record$...1 %in% pop.survivor.control$sjlid)
# FALSE  TRUE 
# 4375   106 

## List of 4402 SJLIFE samples
sjlife_4402 <- read_sas("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/PHENOTYPE/Data_from_Qi_Liu/toyadav.sas7bdat")
sjlife_4402 <- unique(sjlife_4402$sjlid)
table(sjlife_4402 %in% sequencing.Record$...1)
# FALSE  TRUE 
# 27  4375
missing.27.of.4402 <- sjlife_4402[!sjlife_4402 %in% sequencing.Record$...1]

## In preQC sjlife data
preQCsjlife <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/MERGED_sjlife1_2_PreQC/cleaned/MERGED.SJLIFE.1.2.GATKv3.4.VQSR.chr1.preQC_biallelic_renamed_ID_edited.vcf.gz.fam", header = F)
sjlife_4402[!(sjlife_4402 %in% preQCsjlife$V2)]
# "SJL5450006"

## This is TB ID conversion file from Yadav
SJLIFEwesTBID <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/Survivor_WES/SJLIFE_WESsamplelist_TBIDcheck_YSapkota_02Aug2022_FINAL.txt", header = T)

SJLIFEwesTBID$CompBioID_first_part <- sub("_.*", "", SJLIFEwesTBID$CompBioID)
#####################################################################
## This is to remove duplicate samples in TBID file: SJLIFEwesTBID ##
#####################################################################
## There are duplicate SJLIFE IDs, we can romove the ones with low call rate here:
# Check for duplicated IID values
imiss = read.table("/ResearchHome/Groups/sapkogrp/projects/Genomics/common/Survivor_WES_QC/biallelic2/plink_all/chrALL.SJLIFE_CCSS_WES_101724.GATK4180.hg38_biallelic.geno.0.1.hwe.1e-15_missing.imiss", header=TRUE)
# imiss <- imiss[!grepl("CCSS", imiss$IID),]
imiss$IID2 <- imiss$IID
imiss$IID2_first_part <- sub("_.*", "", imiss$IID2)
imiss$IID <- SJLIFEwesTBID$SJLID[match(imiss$IID2_first_part, SJLIFEwesTBID$CompBioID_first_part)]

imiss$CCSSID <- NA
imiss.CCSS.org <- imiss[grepl("SUBJECT", imiss$FID),]
imiss <- imiss[!grepl("SUBJECT", imiss$FID),]
# add CCSS exp
imiss$CCSSID <-  sub(".*CCSS-0", "", imiss$FID) 
imiss$CCSSID <- sub(".*CCSS-", "", imiss$CCSSID)
imiss$CCSSID[!grepl("CCSS",imiss$FID)] <- NA
imiss$IID[!is.na(imiss$CCSSID)] <- imiss$CCSSID[!is.na(imiss$CCSSID)]
# now merge CCSSorg 
imiss <- rbind.data.frame(imiss, imiss.CCSS.org)
# Keep remaining unmatched as is 
imiss$IID[is.na(imiss$IID)] <- imiss$IID2[is.na(imiss$IID)]

duplicate_samples <- imiss[duplicated(imiss$IID) | duplicated(imiss$IID, fromLast = TRUE), ]
dim(duplicate_samples)
# 86

# remove those labelled CCSS from duplicates
remove.dup.ccss <- duplicate_samples[grepl("CCSS", duplicate_samples$FID),]


# Sort and keep only the sample with the lowest F_MISS for each IID
library(dplyr)
samples_to_keep <- duplicate_samples %>%
  group_by(IID) %>%
  slice_min(F_MISS, with_ties = FALSE)  # with_ties = FALSE ensures only one sample is kept per group

# Identify the samples to remove by selecting all not in 'samples_to_keep'
samples_to_remove <- anti_join(duplicate_samples, samples_to_keep, by = c("FID", "IID"))

# View the samples to be removed
# View(samples_to_remove)
dim(samples_to_remove)
# 43  7

# keep SJLIFE and remove CCSS
samples_to_remove <- samples_to_remove[!samples_to_remove$IID %in% remove.dup.ccss$IID,]
samples_to_remove <- rbind.data.frame(samples_to_remove, remove.dup.ccss)
dim(samples_to_remove)
# 43  8
table(SJLIFEwesTBID$CompBioID %in% samples_to_remove$IID2)
# FALSE  TRUE 
# 4977    42 

## Remove these 42 duplicates
SJLIFEwesTBID <- SJLIFEwesTBID[!SJLIFEwesTBID$CompBioID %in% samples_to_remove$IID2,]
dim(SJLIFEwesTBID)
# 4977    4
#####################################################################

# table(sequencing.Record$...1 %in% unique(SJLIFEwesTBID$SJLID))
# # FALSE  TRUE 
# # 21  4460

## All WES samples
all.WES.samples <-  read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/Survivor_WES_QC/sample_mapping_files/WES_samples.txt", header = F)
all.WES.samples$CompBioID_first_part <- sub("_.*", "", all.WES.samples$V1)
dim(all.WES.samples)
# 13726     1
sum(duplicated(all.WES.samples$V1))
# 0
sum(duplicated(all.WES.samples$CompBioID_first_part))
# 57

## remove samples identified as duplicates
table(all.WES.samples$V1 %in% samples_to_remove$IID2)
# FALSE  TRUE 
# 13683    43

table(all.WES.samples$V1 %in% imiss$IID2)
# TRUE 
# 13726

all.WES.samples$renameVCF <- imiss$IID[match(all.WES.samples$V1, imiss$IID2)]
all.WES.samples$renameVCF [all.WES.samples$V1 %in% samples_to_remove$IID2] <- paste0(all.WES.samples$renameVCF [all.WES.samples$V1 %in% samples_to_remove$IID2], "_dups") 

## Survivor population
table(all.WES.samples$renameVCF %in% pop$sjlid)
# FALSE  TRUE 
# 9087  4639

table(all.WES.samples$renameVCF %in% pop.survivor.control$sjlid)
# FALSE  TRUE 
# 13620   106 

# all.WES.samples <- all.WES.samples[!all.WES.samples$V1 %in% samples_to_remove$IID2,]


# table(all.WES.samples$CompBioID_first_part %in% SJLIFEwesTBID$CompBioID_first_part)
# FALSE  TRUE 
# 9042  4630 

## AA sample list from Jenn
AA.90samples <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/RNAseq/common/RNAseq_sjlife/jenn_RNAseq/AA_samples.txt", header = F)
# AA.90samples <- read_excel("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/Survivor_WES_QC/sample_mapping_files/SJLIFE_WGS_and_CCSS_SNP_WES_WGS_samples_and_overlap.xlsx", sheet = "SJLIFEWGS_90AA", col_names = F)
# table(AA.90samples$...1 %in% SJLIFEwesTBID$SJLID)


## Classify population
all.WES.samples$pop <- NA
all.WES.samples$pop [grepl("CCSS_SUBJECT", all.WES.samples$V1)] <- "CCSS_org"
all.WES.samples$pop [all.WES.samples$renameVCF %in% pop.survivor.control$sjlid] <- "C.Control"
all.WES.samples$pop [all.WES.samples$renameVCF %in% pop.survivor$sjlid] <- "Survivor"
all.WES.samples$pop [all.WES.samples$renameVCF %in% AA.90samples$V1] <- "AA"
all.WES.samples$pop[grepl("-CCSS-", all.WES.samples$V1)] <-  "CCSS_exp"
table(all.WES.samples$pop)
# AA C.Control  CCSS_exp  CCSS_org  Survivor 
# 90       106      3034      5451      4443 

# ## combine rename variables
# all.WES.samples$all_rename <- all.WES.samples$renameSJLIFE
# all.WES.samples$all_rename [is.na(all.WES.samples$all_rename)]<- all.WES.samples$renameCCSS_exp[is.na(all.WES.samples$all_rename)]
# all.WES.samples$all_rename [is.na(all.WES.samples$all_rename)]<- all.WES.samples$renameCCSS_org[is.na(all.WES.samples$all_rename)]
# 
# cc <- all.WES.samples[is.na(all.WES.samples$all_rename),]

## 2. Check how many of SJLIFE are available in WES
table(sjlife_4402 %in% all.WES.samples$renameVCF)
# FALSE  TRUE 
# 21  4381 

View(as.data.frame(sjlife_4402[!sjlife_4402 %in% all.WES.samples$renameVCF]))
missing.sjlife.in.wes <- sjlife_4402[!sjlife_4402 %in% all.WES.samples$renameVCF]

sjlife.to.be.sequenced <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/Survivor_WES_QC/sample_mapping_files/27April2023_SJLIFE_WGSround4_SampleList.txt", header = T, sep = "\t")
table(sjlife.to.be.sequenced$sjlid %in% all.WES.samples$renameSJLIFE)
## Yet to be sequenced sjlife
table(sjlife.to.be.sequenced$sjlid %in% missing.sjlife.in.wes)
# FALSE 
# 569

## Missing SJLIFE of sequencing record from WES
table(sequencing.Record$...1 %in% all.WES.samples$renameVCF)
# FALSE  TRUE 
# 21  4460 
View(as.data.frame(sequencing.Record$...1[!sequencing.Record$...1 %in% all.WES.samples$renameVCF]))
missing.sjlife.fromseq.records.in.WES <- sequencing.Record$...1[!sequencing.Record$...1 %in% all.WES.samples$renameVCF]
table(missing.sjlife.fromseq.records.in.WES %in% missing.sjlife.in.wes)
# TRUE 
# 21 

## Additional SJLIFE from sequencing record
sequencing.Record.additional <- read_excel("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/Survivor_WES_QC/sample_mapping_files/SJLIFE_WGS_and_CCSS_SNP_WES_WGS_samples_and_overlap.xlsx", sheet = "SJLIFEWGS_additional817", col_names = F)
table(sequencing.Record.additional$...1 %in% sjlife.to.be.sequenced$sjlid)
# FALSE  TRUE 
# 248   569
table(all.WES.samples$renameVCF %in% sequencing.Record.additional$...1)
# FALSE  TRUE 
# 13649    32
table(missing.sjlife.in.wes %in% sequencing.Record.additional$...1)
# FALSE 
# 21

## 2. Check how many of CCSS_exp are available in WES
sequencing.Record.ccss_exp <- read_excel("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/Survivor_WES_QC/sample_mapping_files/SJLIFE_WGS_and_CCSS_SNP_WES_WGS_samples_and_overlap.xlsx", sheet = "CCSS_exp_WGS_2839", col_names = F)
table(sequencing.Record.ccss_exp$...1 %in% all.WES.samples$renameVCF)
# FALSE  TRUE 
# 2  2837
sequencing.Record.ccss_exp$...1[!sequencing.Record.ccss_exp$...1 %in% all.WES.samples$renameVCF]
# 3267622 2512001

## 3. Unidentified samples
unidentified.WES <- all.WES.samples[is.na(all.WES.samples$pop),]

unidentified.WES$TBID <- sub(".*(TB-[0-9-]+).*", "\\1", unidentified.WES$V1)

write.table(unidentified.WES, "Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/Survivor_WES_QC/sample_mapping_files/Unidentified_WES_samples.txt", col.names = T, row.names = F, sep = "\t", quote = F)
#######################
## B. working on WGS ##
#######################
## 1. check SJLIFE
vcf.samples.wgs <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects//Genomics/common/Survivor_WGS_QCed/Phenotypes/master_vcf_rename.txt", header = T)
vcf.samples.wgs$renameCCSS_exp <- NA

# Rename ccss_exp
vcf.samples.wgs$renameCCSS_exp [grepl("-CCSS-[0-9]+", vcf.samples.wgs$VCFID)] <- vcf.samples.wgs$VCFID [grepl("-CCSS-[0-9]+", vcf.samples.wgs$VCFID)]
vcf.samples.wgs$renameCCSS_exp <- sub(".*CCSS-0", "", vcf.samples.wgs$renameCCSS_exp) 
vcf.samples.wgs$renameCCSS_exp <- sub(".*CCSS-", "", vcf.samples.wgs$renameCCSS_exp) 


table(sequencing.Record$...1 %in% vcf.samples.wgs$VCFrename)
sequencing.Record$...1[!sequencing.Record$...1 %in% vcf.samples.wgs$VCFrename]
# "SJL5052915" "SJL5057815"
sjlife_4402[!sjlife_4402 %in% vcf.samples.wgs$VCFrename]

## 1. check CCSS_exp
table(sequencing.Record.ccss_exp$...1 %in% vcf.samples.wgs$VCFrename)
# FALSE  TRUE 
# 143  2696 


table(sequencing.Record.ccss_exp$...1 %in% vcf.samples.wgs$renameCCSS_exp)
# 2839


missing.ccss.exp <- sequencing.Record.ccss_exp$...1[!sequencing.Record.ccss_exp$...1 %in% vcf.samples.wgs$VCFrename]
missing.ccss.exp <-  vcf.samples.wgs[vcf.samples.wgs$renameCCSS_exp %in% missing.ccss.exp,]

table(missing.ccss.exp$VCFrename %in% sjlife_4402)
# 143
table(missing.ccss.exp$VCFrename %in% pop.survivor.control$sjlid)
# 143

vcf.samples.wgs$VCFrename[(vcf.samples.wgs$VCFrename %in% missing.ccss.exp$VCFrename)] <- vcf.samples.wgs$renameCCSS_exp[(vcf.samples.wgs$VCFrename %in% missing.ccss.exp$VCFrename)]
table(vcf.samples.wgs$VCFrename == vcf.samples.wgs$renameCCSS_exp)

# ## These 143 samples were supposed to be CCSS to need to correct the IDs
# write.table(vcf.samples.wgs[1:2], "Z:/ResearchHome/Groups/sapkogrp/projects//Genomics/common/Survivor_WGS_QCed/Phenotypes/master_vcf_rename_corrected.txt", col.names = F)

vcf.samples.wgs$pop <- NA
vcf.samples.wgs$pop [vcf.samples.wgs$VCFrename %in% pop.survivor$sjlid] <- "Survivor"
vcf.samples.wgs$pop [vcf.samples.wgs$VCFrename %in% pop.survivor.control$sjlid] <- "C.Control"
vcf.samples.wgs$pop [vcf.samples.wgs$VCFrename %in% AA.90samples$V1] <- "AA"
vcf.samples.wgs$pop [grepl("CCSS", vcf.samples.wgs$VCFID)] <- "CCSS_exp"
table((vcf.samples.wgs$pop))
# AA C.Control  CCSS_exp  Survivor 
# 90       450      2998      4425

vcf.samples.wgs[is.na(vcf.samples.wgs$pop),]

vcf.samples.wgs.renamed <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects//Genomics/common/Survivor_WGS_QCed/VCF_original/vcf_samples_renamed.txt", header = F)


