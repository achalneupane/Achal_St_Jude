library(data.table)
setwd("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/WES_rare_variant/sjlife")

## 1. clinvar
clinvar <- read.delim("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/Survivor_WES_QC//annotation//snpEff_round3_preQC/all_new_clinvar_P_LP.txt", header = T, sep = "\t", stringsAsFactors = F)
dim(clinvar)
clinvar$SNP <- sub(";.*", "", clinvar$ID)
clinvar$AF <- as.numeric(clinvar$AF)
clinvar$AF_nfe <- as.numeric(clinvar$AF_nfe)
clinvar$AF_afr <- as.numeric(clinvar$AF_afr)


# make rare; .all is based on allee frequency from gnomAD global population; .eur is gnomAD EUR_nfe; .afr is gnomAD AFR
clinvar.all <- clinvar[which(clinvar$AF <0.01),]
clinvar.eur <- clinvar[which(clinvar$AF <0.01 & clinvar$AF_nfe < 0.01),]
clinvar.afr <- clinvar[which(clinvar$AF <0.01 & clinvar$AF_afr < 0.01),]


loftee <- read.delim("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/Survivor_WES_QC//annotation//snpEff_round3_preQC/loftee/loftee_HC_all_chr_with_gnomad.txt", header = T, sep = "\t", stringsAsFactors = F)
dim(loftee)
loftee$SNP <- sub(";.*", "", loftee$Uploaded_variation)
# LOF variants flagged by LOFTEE as dubious (e.g., affecting poorly conserved exons and splice variants affecting NAGNAG sites or non-canonical splice regions) will be excluded.
loftee <- loftee[!grepl("NAGNAG_SITE|NON_CAN_SPLICE", loftee$LoF_flags),]
table(loftee$LoF_flags)
# - 
# 70494 



loftee$AF <- as.numeric(loftee$AF.1)
loftee$AF_nfe <- as.numeric(loftee$AF_nfe)
loftee$AF_afr <- as.numeric(loftee$AF_afr)
loftee$CHROM  <- sub("([0-9XY]+):.+", "\\1", loftee$SNP)
loftee$POS <- sub("chr[0-9XY]+:(\\d+):.+", "\\1", loftee$SNP)

# make rare
loftee.all <- loftee[which(loftee$AF <0.01),]
loftee.eur <- loftee[which(loftee$AF <0.01 & loftee$AF_nfe < 0.01),]
loftee.afr <- loftee[which(loftee$AF <0.01 & loftee$AF_afr < 0.01),]


snpeff <- read.delim("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/Survivor_WES_QC//annotation//snpEff_round3_preQC/missense_variants_with_overlap_in_more_than_90_percent_of_prediction_tools_all_cols.txt", header = T, sep = "\t", stringsAsFactors = F)
dim(snpeff)
snpeff$SNP <- sub(";.*", "", snpeff$ID)

snpeff$AF <- as.numeric(snpeff$AF)
snpeff$AF_nfe <- as.numeric(snpeff$AF_nfe)
snpeff$AF_afr <- as.numeric(snpeff$AF_afr)


# make rare
snpeff.all <- snpeff[which(snpeff$AF <0.01),]
snpeff.eur <- snpeff[which(snpeff$AF <0.01 & snpeff$AF_nfe < 0.01),]
snpeff.afr <- snpeff[which(snpeff$AF <0.01 & snpeff$AF_afr < 0.01),]


length(unique(c(clinvar$SNP, loftee$SNP, snpeff$SNP)))
# 46076
cc <- as.data.frame(unique(c(clinvar$SNP, loftee$SNP, snpeff$SNP)))
dim(cc)
# 46076
write.table(cc, "Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/WES_rare_variant/sjlife/rare_variants_to_extract.txt", row.names = F, col.names = F, quote = F)

# ## Extract SJLIFE samples
# preQCsjlife <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/MERGED_sjlife1_2_PreQC/cleaned/MERGED.SJLIFE.1.2.GATKv3.4.VQSR.chr1.preQC_biallelic_renamed_ID_edited.vcf.gz.fam", header = F)
# all.WES.samples <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/Survivor_WES_QC/sample_mapping_files/WES_samples_after_Kubra.txt", header = T)
# table(all.WES.samples$V4 %in% preQCsjlife$V2)
# # FALSE  TRUE 
# # 9240  4486
# all.WES.samples.sjlife <- all.WES.samples[all.WES.samples$V4 %in% preQCsjlife$V2,3:4]
# write.table(all.WES.samples.sjlife, "Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common//WES_rare_variant/sjlife/sjlife_samples.txt", col.names = F, row.names = F, sep = "\t", quote = F)
all.WES.samples <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/Survivor_WES_QC/sample_mapping_files/WES_samples_after_Kubra.txt", header = T)
all.WES.samples.sjlife <- as.data.frame(all.WES.samples[which(all.WES.samples$pop == "Survivor"), c("V3", "V4")])
write.table(all.WES.samples.sjlife, "Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common//WES_rare_variant/sjlife/sjlife_samples.txt", col.names = F, row.names = F, sep = "\t", quote = F)

cc <- cc$`unique(c(clinvar$SNP, loftee$SNP, snpeff$SNP))`

## Final SJLIFE QCed data
bim.QC.sjlife <- fread("Z:/ResearchHome/Groups/sapkogrp/projects//Genomics/common/Survivor_WES_QC/biallelic2//plink_all/chrALL.SJLIFE_CCSS_WES_101724.GATK4180.hg38_biallelic.geno.0.1.hwe.1e-15.LCR.removed.MAC.ge.1_updated.bim")
table(cc %in% bim.QC.sjlife$V2)
# FALSE  TRUE 
# 49730 34697 


bim.QC.sjlife.PLP <- fread("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/WES_rare_variant//sjlife/all_rare_variants_maf0.01_all_sjlife.bim")
table(cc %in% bim.QC.sjlife.PLP$V2)
# FALSE  TRUE 
# 49888 34539
