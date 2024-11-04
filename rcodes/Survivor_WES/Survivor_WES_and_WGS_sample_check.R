library(haven)
# https://wiki.stjude.org/display/CAB/Genetic+Ancestry+Estimation+by+PCA
setwd("Z:/ResearchHome/Groups/sapkogrp/projects/CAB/common/SJLIFE_CCSS_WES_101724/GermlineQC/")

load("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/PHENOTYPE/5_lifestyle_v11.RDATA")
overlaps <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/gwas/pheno/SJLIFE_WGS_samples_overlap_with_CCSS_org_SNP_samples.txt", header = T)
pop <- read_sas("Z:/SJShare/SJCOMMON/ECC/SJLife/SJLIFE Data Freeze/2 Final Data SJLIFE/20200430/Clinical Data/demographics.sas7bdat")
pop.survivor <- pop[pop$studypop == "Survivor",]
community.control <- pop[grepl("Control", pop$studypop),]
table(pop.survivor$sjlid %in% community.control$sjlid)
# FALSE 
# 9366 
######################
## Clean sample IDs ##
######################
SJLIFEwes.combio <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/Survivor_WES/SJLIFE_WESsamplelist_TBIDcheck_YSapkota_02Aug2022_FINAL.txt", header = T)
# SJLIFEwes.combio <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/Survivor_WES/biallelic2/plink_all/Survivors/chr.ALL.SURVIVORS_WES.GATK4180.hg38_biallelic.geno.0.1.hwe.1e-15.LCR.removed.MAC.ge.1_updated_unique.fam", header = F)
sjlife.4507 <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/MERGED_sjlife1_2_PreQC/cleaned/MERGED.SJLIFE.1.2.GATKv3.4.VQSR.chr9.preQC_biallelic_renamed_ID_edited.vcf.gz.fam", header = F)

all.WES <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/Survivor_WES/biallelic2/plink_all/chr.ALL.Survivor_WES.GATK4180.hg38_biallelic.geno.0.1.hwe.1e-15.LCR.removed.MAC.ge.1_ID_updated_unique.fam", header = T)

#########################
## Old WES VCF samples ##
#########################

#### preQC 

old.wes.samples = read.table("/ResearchHome/Groups/sapkogrp/projects/Genomics/common/Survivor_WES/biallelic2/plink_all/chr.ALL.Survivor_WES.GATK4180.hg38_biallelic.geno.0.1.hwe.1e-15.LCR.removed.MAC.ge.1_ID_updated.fam", header=F)
table(sjlife.4507$V2 %in% old.wes.samples$V2)
# FALSE  TRUE 
# 21  4486 

## 1. Check SJLIFE
ss1 <- sjlife.4507$V2[!sjlife.4507$V2 %in% old.wes.samples$V2]
ss1
# [1] "SJL1437501" "SJL1285201" "SJL1224901" "SJL1245501" "SJL1246701" "SJL5037906" "SJL1281212" "SJL5050105" "SJL1679008" "SJL5068708" "SJL1648308" "SJL5052915" "SJL5057815" "SJL5024018"
# [15] "SJL4749916" "SJL1750516" "SJL5058217" "SJL5271801" "SJL5385207" "SJL2521413" "SJL5321816"

table(PHENO.ANY_SN$sjlid %in% old.wes.samples$V2)
# FALSE  TRUE 
# 21  4380 
ss2 <- PHENO.ANY_SN$sjlid[!PHENO.ANY_SN$sjlid %in% old.wes.samples$V2]
ss2
# [1] "SJL1245501" "SJL1246701" "SJL1750516" "SJL5024018" "SJL1648308" "SJL1285201" "SJL1224901" "SJL1679008" "SJL5068708" "SJL1437501" "SJL1281212" "SJL5037906" "SJL5050105" "SJL5052915"
# [15] "SJL4749916" "SJL5057815" "SJL5058217" "SJL2521413" "SJL5271801" "SJL5321816" "SJL5385207"
table(ss1%in%ss2)
# TRUE 
# 21

table(unique(old.wes.samples$V2) %in% pop.survivor$sjlid)

## 2. Check AA samples
AA.90samples <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/RNAseq/common/RNAseq_sjlife/jenn_RNAseq/AA_samples.txt", header = F)
table(AA.90samples$V1 %in% old.wes.samples$V2)
# TRUE 
# 90 
table(AA.90samples$V1 %in% sjlife.4507$V2)
# FALSE 
# 90 

## 3. Check Community controls
table(old.wes.samples$V2 %in% community.control$sjlid)
# FALSE  TRUE 
# 7614   451 

table(sjlife.4507$V2 %in% community.control$sjlid)
# FALSE  TRUE 
# 4401   106 
table(PHENO.ANY_SN$sjlid %in% community.control$sjlid)
# FALSE 
# 4401 
table(AA.90samples$V1 %in% community.control$sjlid)
# FALSE 
# 90 

## 4. Check CCSS_exp
ccss_exp <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/ccss_exp_wgs/QCed_samples", header = F)
table(ccss_exp$V1 %in% old.wes.samples$V2)
# FALSE  TRUE 
# 2  2837 
ccss_exp$V1[!ccss_exp$V1 %in% old.wes.samples$V2]
# [1] 3267622 2512001

table(old.wes.samples$V2 %in% overlaps$V2)



#### postQC 

old.wes.samples.Survivor = read.table("/ResearchHome/Groups/sapkogrp/projects/Genomics/common/Survivor_WES/biallelic2/plink_all/Survivors/chr.ALL.SURVIVORS_WES.GATK4180.hg38_biallelic.geno.0.1.hwe.1e-15.LCR.removed.MAC.ge.1_updated_unique.fam", header=F)
table(sjlife.4507$V2 %in% old.wes.samples.Survivor$V2)
# FALSE  TRUE 
# 139  4368 

## 1. Check SJLIFE
ss3 <- sjlife.4507$V2[!sjlife.4507$V2 %in% old.wes.samples.Survivor$V2]
ss3
# [1] "SJL1437501" "SJL1285201" "SJL1224901" "SJL1245501" "SJL1246701" "SJL5037906" "SJL1281212" "SJL5050105" "SJL1679008" "SJL5068708" "SJL1648308" "SJL5052915" "SJL5057815" "SJL5024018"
# [15] "SJL4749916" "SJL1750516" "SJL5058217" "SJL5271801" "SJL5385207" "SJL2521413" "SJL5321816"

table(PHENO.ANY_SN$sjlid %in% old.wes.samples.Survivor$V2)
# FALSE  TRUE 
# 21  4380 
ss4 <- PHENO.ANY_SN$sjlid[!PHENO.ANY_SN$sjlid %in% old.wes.samples.Survivor$V2]
ss4
# [1] "SJL1245501" "SJL1246701" "SJL1750516" "SJL5024018" "SJL1648308" "SJL1285201" "SJL1224901" "SJL1679008" "SJL5068708" "SJL1437501" "SJL1281212" "SJL5037906" "SJL5050105" "SJL5052915"
# [15] "SJL4749916" "SJL5057815" "SJL5058217" "SJL2521413" "SJL5271801" "SJL5321816" "SJL5385207"
table(ss1%in%ss2)
# TRUE 
# 21

## 2. Check AA samples
AA.90samples <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/RNAseq/common/RNAseq_sjlife/jenn_RNAseq/AA_samples.txt", header = F)
table(AA.90samples$V1 %in% old.wes.samples$V2)
# TRUE 
# 90 
table(AA.90samples$V1 %in% sjlife.4507$V2)
# FALSE 
# 90 

## 3. Check Community controls
table(old.wes.samples$V2 %in% community.control$sjlid)
# FALSE  TRUE 
# 7614   451 

table(sjlife.4507$V2 %in% community.control$sjlid)
# FALSE  TRUE 
# 4401   106 
table(PHENO.ANY_SN$sjlid %in% community.control$sjlid)
# FALSE 
# 4401 
table(AA.90samples$V1 %in% community.control$sjlid)
# FALSE 
# 90 

## 4. Check CCSS_exp
ccss_exp <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/ccss_exp_wgs/QCed_samples", header = F)
table(ccss_exp$V1 %in% old.wes.samples$V2)
# FALSE  TRUE 
# 2  2837 
ccss_exp$V1[!ccss_exp$V1 %in% old.wes.samples$V2]
# [1] 3267622 2512001

table(old.wes.samples$V2 %in% overlaps$V2)


#############################
## WGS SJLIFE and CCSS exp ##
#############################
all.survivors.WGS <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/Survivor_WGS_QCed/QC/Survivor_WGS.GATK4180.hg38_renamed_chr1.PASS.decomposed.qced.fam", header = F)

## 1. Check SJLIFE
ss3 <- sjlife.4507$V2[!sjlife.4507$V2 %in% all.survivors.WGS$V2]
ss3
# [1] "SJL5207307" "SJL5188202" "SJL5128013" "SJL5094513" "SJL5052915" "SJL5057815"

table(PHENO.ANY_SN$sjlid %in% all.survivors.WGS$V2)
# FALSE  TRUE 
# 6  4395
ss2 <- PHENO.ANY_SN$sjlid[!PHENO.ANY_SN$sjlid %in% old.wes.samples$V2]
ss2
# [1] "SJL1245501" "SJL1246701" "SJL1750516" "SJL5024018" "SJL1648308" "SJL1285201" "SJL1224901" "SJL1679008" "SJL5068708" "SJL1437501" "SJL1281212" "SJL5037906" "SJL5050105" "SJL5052915"
# [15] "SJL4749916" "SJL5057815" "SJL5058217" "SJL2521413" "SJL5271801" "SJL5321816" "SJL5385207"
table(ss1%in%ss2)
# TRUE 
# 21

## 2. Check AA samples
AA.90samples <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/RNAseq/common/RNAseq_sjlife/jenn_RNAseq/AA_samples.txt", header = F)
table(AA.90samples$V1 %in% old.wes.samples$V2)
# TRUE 
# 90 
table(AA.90samples$V1 %in% sjlife.4507$V2)
# FALSE 
# 90 

## 3. Check Community controls
table(old.wes.samples$V2 %in% community.control$sjlid)
# FALSE  TRUE 
# 7614   451 

table(sjlife.4507$V2 %in% community.control$sjlid)
# FALSE  TRUE 
# 4401   106 
table(PHENO.ANY_SN$sjlid %in% community.control$sjlid)
# FALSE 
# 4401 
table(AA.90samples$V1 %in% community.control$sjlid)
# FALSE 
# 90 

## 4. Check CCSS_exp
ccss_exp <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/ccss_exp_wgs/QCed_samples", header = F)
table(ccss_exp$V1 %in% old.wes.samples$V2)
# FALSE  TRUE 
# 2  2837 
ccss_exp$V1[!ccss_exp$V1 %in% old.wes.samples$V2]
# [1] 3267622 2512001




## New WES VCF samples in SJLIFE, CCSS_exp and CCSS_org
new.wes.VCF.samples <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/Survivor_WES_QC/biallelic2/VCFsample_names.txt", header = F)
new.wes.VCF.samples$SJLID <- SJLIFEwes.combio$SJLID[match(new.wes.VCF.samples$V1, SJLIFEwes.combio$CompBioID )]

table(unique(new.wes.VCF.samples$SJLID) %in% sjlife.4507$V2)
# FALSE  TRUE 
# 152  4405 
table(unique(new.wes.VCF.samples$SJLID) %in% PHENO.ANY_SN$sjlid)
# FALSE  TRUE 
# 258  4299
table(PHENO.ANY_SN$sjlid %in% unique(new.wes.VCF.samples$SJLID))
# FALSE  TRUE 
# 102  4299 

AA.90samples <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/RNAseq/common/RNAseq_sjlife/jenn_RNAseq/AA_samples.txt", header = F)
table(AA.90samples$V1 %in% SJLIFEwes.combio$SJLID)
# TRUE 
# 90
table(AA.90samples$V1 %in% new.wes.VCF.samples$SJLID)
# TRUE 
# 90

table(unique(new.wes.VCF.samples$SJLID) %in% PHENO.ANY_SN$sjlid)
# FALSE  TRUE 
# 258  4299 
table(unique(SJLIFEwes.combio$SJLID) %in% PHENO.ANY_SN$sjlid)
# FALSE  TRUE 
# 604  4380 

table(SJLIFEwes.combio$CompBioID %in% new.wes.VCF.samples$V1)
# 4589
table(SJLIFEwes.combio$SJLID %in% unique(new.wes.VCF.samples$V2))

unique(SJLIFEwes.combio$V3)[!unique(SJLIFEwes.combio$V3) %in% PHENO.ANY_SN$sjlid]

PHENO.ANY_SN$sjlid[!PHENO.ANY_SN$sjlid %in% unique(SJLIFEwes.combio$V3)]
# [1] "SJL1245501" "SJL1246701" "SJL1750516" "SJL5024018" "SJL1648308" "SJL1285201" "SJL1224901" "SJL1679008" "SJL5068708" "SJL1437501" "SJL1281212" "SJL5037906" "SJL5050105" "SJL5052915"
# [15] "SJL4749916" "SJL5057815" "SJL5058217" "SJL2521413" "SJL5271801" "SJL5321816" "SJL5385207"

# all.samples[duplicated(all.samples$V2),]

# # Keep only those from Germline QC
# germline.QC <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/CAB/common/SJLIFE_CCSS_WES_101724/GermlineQC//Survivor_WES.fam")
# dim(germline.QC)
# # 8027
# sum(all.samples$V1 %in% germline.QC$V1)
# all.samples <- all.samples[all.samples$V1 %in% germline.QC$V1,]
dim(all.samples)
# 8055
# remove sex problem and het filtered
# all.samples.original <- all.samples
# all.samples <- all.samples[!all.samples$V1 %in% sex.problem.mind.and.3sd.outliers,]
# dim(all.samples)
sum(!all.samples$V1 %in% sex.problem.mind.and.3sd.outliers)
# 7907 # after sample level QC

## 8065 total samples in WES data
pop <- read_sas("Z:/SJShare/SJCOMMON/ECC/SJLife/SJLIFE Data Freeze/2 Final Data SJLIFE/20200430/Clinical Data/demographics.sas7bdat")
pop.survivor <- pop[pop$studypop == "Survivor",]
pop.survivor.control <- pop[grepl("Control", pop$studypop),]











#####################
## How many in WGS ##
#####################
WGS_SJLID <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/MERGED_sjlife1_2_PreQC/cleaned/MERGED.SJLIFE.1.2.GATKv3.4.VQSR.chr9.preQC_biallelic_renamed_ID_edited.vcf.gz.fam")
## Before QC
SJLID <- all.samples.before.QC$V2[all.samples.before.QC$cohort == "Survivor"]
sum(WGS_SJLID$V2 %in% SJLID)
## 4380 ## Survivors from WGS in Survivor WES

SJLID <- all.samples.before.QC$V2[all.samples.before.QC$cohort == "Community_control"]
sum(WGS_SJLID$V2 %in% SJLID)
# 106 ## Survivors from WGS in Community_control WES
## total 4486 in WGS if we include both Survivor and community control in SJLIFE before QC

## AFter QC
SJLID <- all.samples$V2[all.samples$cohort == "Survivor"]
sum(WGS_SJLID$V2 %in% SJLID)
## 4368 ## Survivors from WGS in Survivor WES

SJLID <- all.samples$V2[all.samples$cohort == "Community_control"]
sum(WGS_SJLID$V2 %in% SJLID)
# 106 ## Survivors from WGS in Community_control WES

## total 4474 in WGS if we include both Survivor and community control in SJLIFE

## Rename sample ID
setwd("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/Survivor_WES/biallelic/plink_all")
all.samples <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/Survivor_WES/biallelic/plink_all/chr19.Survivor_WES.GATK4180.hg38_biallelic.fam", header = F)

all.samples <- all.samples[1:2]

SJLIFE <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/Survivor_WES/SJLIFE_WESsamplelist_TBIDcheck_YSapkota_02Aug2022_FINAL.txt", header = T)
sum(duplicated(SJLIFE$SJLID))

all.samples$V3 <- all.samples$V1
all.samples$V4 <- SJLIFE$SJLID[match(all.samples$V2, SJLIFE$CompBioID)]
all.samples$V4[is.na(all.samples$V4)] <- all.samples$V2[is.na(all.samples$V4)]

all.samples$V4 <- sub(".*CCSS-0", "", all.samples$V4)
all.samples$V4 <- sub(".*CCSS-", "", all.samples$V4)

write.table(all.samples, "Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/Survivor_WES/rename_samples.txt", col.names = F, row.names = F, quote = F)


## There are duplicate SJLIFE IDs, we can romove the ones with low call rate here:
# Check for duplicated IID values
imiss = read.table("/ResearchHome/Groups/sapkogrp/projects/Genomics/common/Survivor_WES/biallelic2/plink_all/Survivors/chr.ALL.SURVIVORS_WES.GATK4180.hg38_biallelic.geno.0.1.hwe.1e-15.LCR.removed.MAC.ge.1_updated_missing.imiss", header=TRUE)
duplicate_samples <- imiss[duplicated(imiss$IID) | duplicated(imiss$IID, fromLast = TRUE), ]
dim(duplicate_samples)
# 69

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
## 35
write.table(samples_to_remove[, c("FID", "IID")], "Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/Survivor_WES/biallelic2/plink_all/Survivors/duplicate_samples_to_remove.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)

################################################################################################################
################################################################################################################
################################################################################################################
################################################################################################################
################################################################################################################
################################################################################################################
################################################################################################################
## heterozygosity
# # Load the data
# het_data <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/CAB/common/SJLIFE_CCSS_WES_101724/GermlineQC/heterozygosity.het", header=TRUE)
# het_data$F <- as.numeric(het_data$F) 
# het_data$IID <- factor(het_data$IID) 
# 
# # Calculate the mean and standard deviation of observed heterozygosity
# mean_het <- mean(het_data$F)
# sd_het <- sd(het_data$F)
# # Calculate the threshold
# threshold <- mean_het + 3 * sd_het
# # Identify outliers
# outliers <- het_data$IID[het_data$F > threshold]
# outliers
# # "SJRB056831_G1-TB-16-10618"
# 
# mean_het <- mean(het_data$F)
# sd_het <- sd(het_data$F)
# 
# 
# 
# # Load required libraries (if not already installed)
# library(ggplot2)
# 
# # Define the thresholds
# thresholds <- c(-8, -7, -6, -5, -4, -3, 3, 4, 5, 6, 7, 8)
# threshold_values <- mean_het + thresholds * sd_het
# 
# plot <- ggplot(het_data, aes(x = IID, y = F)) +
#   geom_point(aes(color = F < threshold_values[6] | F > threshold_values[7]), size = 3) +
#   geom_hline(yintercept = mean_het, linetype = "dashed", color = "white") +  # Add mean line
#   scale_color_manual(values = c("TRUE" = "red", "FALSE" = "blue")) +
#   geom_hline(yintercept = threshold_values, linetype = "dashed", color = "red") +
#   labs(title = "Scatter Plot of Samples by Cutoff Threshold",
#        x = "Samples",
#        y = "F Value",
#        color = "Outlier") +
#   
#   # Adjust plot appearance
#   theme_minimal()
# 
# # Label the threshold lines
# labels <- data.frame(
#   threshold = thresholds,
#   y = threshold_values,
#   x = rep(1, 2)  # Set x to 1 for both labels
# )
# 
# plot + geom_text(data = labels, aes(x = x, y = y, label = threshold), vjust = -0.5, hjust = 0, size = 4)
# 
# 
# write.table(outliers, file="Z:/ResearchHome/Groups/sapkogrp/projects/CAB/common/SJLIFE_CCSS_WES_101724/GermlineQC//heterozygosity_outside_8_std_outliers.txt", row.names=FALSE, col.names=FALSE, quote=FALSE)




