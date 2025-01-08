library(haven)
# https://wiki.stjude.org/display/CAB/Genetic+Ancestry+Estimation+by+PCA
setwd("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/Survivor_WES/GermlineQC/")

## Read kinship data
# kinship <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/Survivor_WES/GermlineQC/Survivor_WES.kinship", header = T)
sexcheck <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/Survivor_WES/GermlineQC/Survivor_WES.sexcheck", header = T)
sexcheck <- sexcheck[sexcheck$PEDSEX !=0 & sexcheck$SNPSEX !=0,]
sex.problem <- sexcheck[sexcheck$STATUS != "OK",]
write.table(sex.problem, "Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/Survivor_WES/GermlineQC/sex_mismatch.txt", col.names = T, row.names = F, quote = F)


#Read data
imiss = read.table("/ResearchHome/Groups/sapkogrp/projects/Genomics/common/Survivor_WES/biallelic/plink_all/chrALL.Survivor_WES.GATK4180.hg38_biallelic.geno.0.1.hwe.1e-15_missing.imiss", header=TRUE)
# imiss = read.table("/ResearchHome/Groups/sapkogrp/projects/Genomics/common/Survivor_WES/GermlineQC/Survivor_WES_missing.imiss", header=TRUE)

hetCalc = read.table("/ResearchHome/Groups/sapkogrp/projects/Genomics/common/Survivor_WES/biallelic/plink_all/chrALL.Survivor_WES.GATK4180.hg38_biallelic.geno.0.1.hwe.1e-15_maf0.05_indep_prune_heterozygosity.het", header=TRUE)
# hetCalc = read.table("/ResearchHome/Groups/sapkogrp/projects/Genomics/common/Survivor_WES/GermlineQC/Survivor_WES_heterozygosity.het", header=TRUE)
# hetCalc$INDV = paste(hetCalc$FID, hetCalc$IID, sep="_")
# hetCalc$INDV[1:4] = as.character(hetCalc$FID)

hetCalc$het<-(hetCalc$N.NM-hetCalc$O.HOM)/hetCalc$N.NM
mean_het = mean(hetCalc$het)
sd_het = sd(hetCalc$het)
lower_het = mean_het-3*sd_het
upper_het = mean_het+3*sd_het
outlier = hetCalc[hetCalc$het<lower_het | hetCalc$het>upper_het,]
dim(outlier)
# 33 at 5 sd
# 55 at 3 sd
germlineQC.fam <- read.table("/ResearchHome/Groups/sapkogrp/projects/Genomics/common/Survivor_WES/GermlineQC/Survivor_WES.fam")
sum(outlier$IID %in% germlineQC.fam$V2)

write.table(outlier, 'SJLIFE_WES_Per_sample_heterozygosity_outlier_check_3sd.txt', row.names=F, quote=F, sep="\t")

m = merge(imiss, hetCalc, by="IID")

#pdf("Per-sample_heterozygosity_rate_vs_proportion_of_missing_genotypes.pdf")
pdf("Per-sample_heterozygosity_rate_vs_proportion_of_missing_genotypes_3-5sd.pdf")
plot(m$F_MISS,m$het, pch=20, xlab="Proportion of missing genotypes", ylab="Heterozygosity rate")
abline(h=mean(hetCalc$het),col="red",lwd=2)
abline(h=mean(hetCalc$het)-(3*sd(hetCalc$het)),col="blue",lty=2)
abline(h=mean(hetCalc$het)+(3*sd(hetCalc$het)),col="blue",lty=2)
abline(h=mean(hetCalc$het)-(4*sd(hetCalc$het)),col="blue",lty=2)
abline(h=mean(hetCalc$het)+(4*sd(hetCalc$het)),col="blue",lty=2)
abline(h=mean(hetCalc$het)-(5*sd(hetCalc$het)),col="blue",lty=2)
abline(h=mean(hetCalc$het)+(5*sd(hetCalc$het)),col="blue",lty=2)
dev.off()


lower_het = mean_het-3*sd_het
upper_het = mean_het+3*sd_het
outlier.3sd = hetCalc[hetCalc$het<lower_het | hetCalc$het>upper_het,]

lower_het = mean_het-5*sd_het
upper_het = mean_het+5*sd_het
outlier.5sd = hetCalc[hetCalc$het<lower_het | hetCalc$het>upper_het,]

lower_het = mean_het-8*sd_het
upper_het = mean_het+8*sd_het
outlier.8sd = hetCalc[hetCalc$het<lower_het | hetCalc$het>upper_het,]


# Write the list of outliers to a file
sex.problem.mind.and.3sd.outliers <- unique(c(as.character(sex.problem$IID), outlier.3sd$IID, "SJALL019091_G1-TB-10-0710", "SJCNS017896_G1-TB-11-1618"))
write.table(sex.problem.mind.and.3sd.outliers, file="Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/Survivor_WES/GermlineQC/sex.problem.mind.and.3sd.outliers.txt", row.names=FALSE, col.names=FALSE, quote=FALSE)

all.samples <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/Survivor_WES/Survivor_WES.samplelist", header = F)
all.samples.sjlife <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/Survivor_WES/SJLIFE_WESsamplelist_TBIDcheck_YSapkota_02Aug2022_FINAL.txt", header = T)
all.samples.ccss <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/Survivor_WES/Survivor_WES.samplelist.ccss", header = F)
table(all.samples$V1 %in% all.samples.sjlife$CompBioID)
# FALSE  TRUE 
# 3046  5019

all.samples$SJLIFE <- all.samples.sjlife$SJLID[match(all.samples$V1, all.samples.sjlife$CompBioID)]
all.samples$CCSS_exp <- all.samples.ccss$V1[match(all.samples$V1, all.samples.ccss$V1)]
all.samples$CCSS_exp <- sub(".*CCSS-0", "", all.samples$CCSS_exp)
all.samples$CCSS_exp <- sub(".*CCSS-", "", all.samples$CCSS_exp)

all.samples.saved <- all.samples
all.samples$newID <- all.samples$SJLIFE
all.samples$newID[is.na(all.samples$newID)] <- all.samples$CCSS_exp[is.na(all.samples$newID)]
sum(is.na(all.samples$newID))
# 10 # these are not SJLIFE
all.samples <- all.samples[!is.na(all.samples$newID),]
dim(all.samples)
# 8055    4
all.samples <- all.samples[c("V1", "newID")]

WES.fam <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/Survivor_WES/biallelic2/plink_all/chr1.Survivor_WES.GATK4180.hg38_biallelic.fam", header = F)
WES.fam <- WES.fam[1:2]
WES.fam$V3 <- WES.fam$V1
WES.fam$V4 <- all.samples$newID[match(WES.fam$V3, all.samples$V1)]
all.samples <- WES.fam
all.samples <- all.samples[!is.na(all.samples$V4),]
write.table(all.samples, "Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/Survivor_WES/rename_samples.txt", col.names = F, row.names = F, quote = F)


## There are duplicate SJLIFE IDs, we can romove the ones with low call rate here:
# Check for duplicated IID values
imiss = read.table("/ResearchHome/Groups/sapkogrp/projects/Genomics/common/Survivor_WES/biallelic2/plink_all/chr.ALL.SURVIVORS_WES.GATK4180.hg38_biallelic.geno.0.1.hwe.1e-15.LCR.removed.MAC.ge.1_updated_missing.imiss", header=TRUE)
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

#########################
## Separate by cohorts ##
#########################

## 8065 total samples in WES data
pop <- read_sas("Z:/SJShare/SJCOMMON/ECC/SJLife/SJLIFE Data Freeze/2 Final Data SJLIFE/20200430/Clinical Data/demographics.sas7bdat")
pop.survivor <- pop[pop$studypop == "Survivor",]
pop.survivor.control <- pop[grepl("Control", pop$studypop),]



sum(all.samples$V4 %in% pop$sjlid)
# 5019
sum(unique(all.samples$V4) %in% pop$sjlid)
# 4984

all.samples$cohort <- NA
all.samples$cohort[grepl("CCSS|SJNPC018728_G1|SJNPC018729_G1", all.samples$V1)] <-  "CCSS"
all.samples$cohort[all.samples$V4 %in% pop.survivor.control$sjlid ] <- "Community_control"
all.samples$cohort[all.samples$V4 %in% pop.survivor$sjlid ] <- "Survivor"
table(all.samples$cohort)
# CCSS Community_control          Survivor 
# 3036               451              4568 
all.samples.before.QC <- all.samples

dim(all.samples)
# 8055

sum(!all.samples$V1 %in% sex.problem.mind.and.3sd.outliers)
# 7907 # after sample level QC

## OPTION 1: Separating cohorts after removing problematic samples from sample-level QC
# ## Remove --mind, heterozygosity, and discordant sex:
# all.samples <- all.samples[!all.samples$V1 %in% sex.problem.mind.and.3sd.outliers,]
# 
# 
# table(all.samples$cohort)
# # CCSS Community_control          Survivor 
# # 2911               445              4551 
# 
# CCSS <- as.data.frame(all.samples$V1[all.samples$cohort == "CCSS"])
# ## Not sure about these two samples in CCSS
# # SJNPC018728_G1  SJNPC018728_G1
# # SJNPC018729_G1  SJNPC018729_G1
# 
# dim(CCSS)
# # [1] 2911    1
# write.table(CCSS, "Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/Survivor_WES/biallelic/extract_CCSS.samples.txt", col.names = F, row.names = F, quote = F)
# 
# SJLIFE_control <- as.data.frame(all.samples$V1[all.samples$cohort == "Community_control"])
# dim(SJLIFE_control)
# # [1] 445   1
# 
# write.table(SJLIFE_control, "Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/Survivor_WES/biallelic/extract_SJLIFE_survivor_control.txt", col.names = F, row.names = F, quote = F)
# 
# SJLIFE_survivor <-  as.data.frame(all.samples$V1[all.samples$cohort == "Survivor"])
# dim(SJLIFE_survivor)
# # [1] 4551    1
# 
# write.table(SJLIFE_survivor, "Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/Survivor_WES/biallelic/extract_SJLIFE_survivor.txt", col.names = F, row.names = F, quote = F)



## OPTION 2: Separating cohorts without removing problematic samples from sample-level QC
table(all.samples$cohort)
# CCSS Community_control          Survivor 
# 3036               451              4568 

CCSS <- as.data.frame(all.samples[all.samples$cohort == "CCSS",c("V3", "V4")])
## Not sure about these two samples in CCSS
# SJNPC018728_G1  SJNPC018728_G1
# SJNPC018729_G1  SJNPC018729_G1

dim(CCSS)
# [1] 3036    1
write.table(CCSS, "Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/Survivor_WES/biallelic2/extract_CCSS.samples_iid_fid.no.sample.QC.txt", col.names = F, row.names = F, quote = F)
## remove SJLIFE overlaps from CCSS
sjlife_ccss_exp_overlaps <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/Survivor_WES_QC/sjlife_ccss_exp_overlaps.txt", header = TRUE, stringsAsFactors = FALSE, sep = "\t")
ccss.removed.sjlife <- CCSS[!(CCSS$V3 %in% sjlife_ccss_exp_overlaps$VCFID),]
write.table(ccss.removed.sjlife, "Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/Survivor_WES/biallelic2/extract_CCSS_without_sjlife_overlaps.samples_iid_fid.no.sample.QC.txt", col.names = F, row.names = F, quote = F)

SJLIFE_control <- as.data.frame(all.samples[all.samples$cohort == "Community_control",c("V3", "V4")])
dim(SJLIFE_control)
# [1] 451   2

write.table(SJLIFE_control, "Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/Survivor_WES/biallelic2/extract_SJLIFE_survivor_control_iid_fid.no.sample.QC.txt", col.names = F, row.names = F, quote = F)

SJLIFE_survivor <-  as.data.frame(all.samples[all.samples$cohort == "Survivor", c("V3", "V4")])
dim(SJLIFE_survivor)
# [1] 4568    2

write.table(SJLIFE_survivor, "Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/Survivor_WES/biallelic2/extract_SJLIFE_survivor_iid_fid.no.sample.QC.txt", col.names = F, row.names = F, quote = F)

write.table(all.samples, "Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/Survivor_WES/all_samples_with_cohort.txt", col.names = F, row.names = F, quote = F)

