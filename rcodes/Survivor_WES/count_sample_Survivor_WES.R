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




all.samples <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/Survivor_WES/biallelic/sample_mapping.txt", header = F)
# > dim(all.samples)
# [1] 8055    2

# all.samples[duplicated(all.samples$V2),]

# # Keep only those from Germline QC
# germline.QC <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/Survivor_WES/GermlineQC//Survivor_WES.fam")
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



sum(all.samples$V2 %in% pop$sjlid)
# 5019

all.samples$cohort <- NA
all.samples$cohort[grepl("CCSS|SJNPC018728_G1|SJNPC018729_G1", all.samples$V1)] <-  "CCSS"
all.samples$cohort[all.samples$V2 %in% pop.survivor.control$sjlid ] <- "Community_control"
all.samples$cohort[all.samples$V2 %in% pop.survivor$sjlid ] <- "Survivor"
table(all.samples$cohort)
# CCSS Community_control          Survivor 
# 3036               451              4568 
all.samples.before.QC <- all.samples

## Remove --mind, heterozygosity, and discordant sex:
all.samples <- all.samples[!all.samples$V1 %in% sex.problem.mind.and.3sd.outliers,]
table(all.samples$cohort)
# CCSS Community_control          Survivor 
# 2911               445              4551 

CCSS <- as.data.frame(all.samples$V1[all.samples$cohort == "CCSS"])
## Not sure about these two samples in CCSS
# SJNPC018728_G1  SJNPC018728_G1
# SJNPC018729_G1  SJNPC018729_G1

dim(CCSS)
# [1] 2911    1
write.table(CCSS, "Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/Survivor_WES/biallelic/extract_CCSS.samples.txt", col.names = F, row.names = F, quote = F)

SJLIFE_control <- as.data.frame(all.samples$V1[all.samples$cohort == "Community_control"])
dim(SJLIFE_control)
# [1] 445   1

write.table(SJLIFE_control, "Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/Survivor_WES/biallelic/extract_SJLIFE_survivor_control.txt", col.names = F, row.names = F, quote = F)

SJLIFE_survivor <-  as.data.frame(all.samples$V1[all.samples$cohort == "Survivor"])
dim(SJLIFE_survivor)
# [1] 4551    1

write.table(SJLIFE_survivor, "Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/Survivor_WES/biallelic/extract_SJLIFE_survivor.txt", col.names = F, row.names = F, quote = F)


## How many in WGS
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


################################################################################################################
################################################################################################################
################################################################################################################
################################################################################################################
################################################################################################################
################################################################################################################
################################################################################################################
## heterozygosity
# # Load the data
# het_data <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/Survivor_WES/GermlineQC/heterozygosity.het", header=TRUE)
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
# write.table(outliers, file="Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/Survivor_WES/GermlineQC//heterozygosity_outside_8_std_outliers.txt", row.names=FALSE, col.names=FALSE, quote=FALSE)




