library(readxl)

ccss <- read_excel("Z:/ResearchHome/ClusterHome/aneupane/data/request/cindy/colorectal/CCSS_colorectal_patients_v2.xlsx", trim_ws = TRUE)
ccss$ccssid <- trimws(ccss$ccssid)
sjlife <- read_excel("Z:/ResearchHome/ClusterHome/aneupane/data/request/cindy/colorectal/SJLIFE_colorectal_patients_v2.xlsx", trim_ws = TRUE)
overlaps <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/gwas/pheno/SJLIFE_WGS_samples_overlap_with_CCSS_org_SNP_samples.txt", header = T)

# When looking at the CCSS, please exclude any of the below SJLIFE survivors in the list of 90 or so.
to.remove <- c("SJL1246507", "SJL1295012", "SJL1417101", "SJL1645209", "SJL1685008", "SJL1727210", "SJL1730507", "SJL1740509", "SJL1741409", "SJL1742107", "SJL1742207", "SJL1757809", "SJL1760413", "SJL1774907", "SJL2505208", "SJL2519701", "SJL4182801", "SJL4195607", "SJL4743510", "SJL4793007", "SJL5003109", "SJL5005601", "SJL5007607", "SJL5011009", "SJL5039417", "SJL5045917", "SJL5048613", "SJL5066008", "SJL5112306", "SJL5117517", "SJL5343912", "SJL5404904")
overlaps$remove.ccss <- overlaps$V1 %in% to.remove
overlaps <- overlaps[overlaps$remove.ccss == TRUE,]

table(ccss$ccssid %in% overlaps$ccssid)
ccss <- ccss[!ccss$ccssid %in% overlaps$ccssid,]

ccss_org <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/ccss_org_hrc/plink/CCSS_org_GRCh37_chr22.fam", header = F)
ccss_exp <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/ccss_exp_wgs/preQC_VCF_per_chromosome/CCSS_exp_biallelic_chr22_ID_edited.fam", header = F)
WES <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/Survivor_WES/biallelic/plink_all/chr.ALL.Survivor_WES.GATK4180.hg38_biallelic_ID_updated.fam", header = F)
WES <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/Survivor_WES/biallelic/plink_all/chr.ALL.Survivor_WES.GATK4180.hg38_biallelic_ID_updated.fam", header = F)
library(stringr)

# Split the values in column v1 by "_" and keep the first item
ccss_org$ccssid <- str_split(ccss_org$V1, "_", simplify = TRUE)[, 1]

ccss$In_ccss_org <- ifelse(ccss$ccssid %in% ccss_org$ccssid, 1,0)
ccss$In_ccss_exp <- ifelse(ccss$ccssid %in% ccss_exp$V1, 1,0)
ccss$In_ccss_wes <- ifelse(ccss$ccssid %in% WES$V2, 1,0)

write.table(ccss, "Z:/ResearchHome/ClusterHome/aneupane/data/request/cindy/colorectal/ccss.genetic.data.check.txt", col.names = T, row.names = F, sep = "\t")

## On 6/5/2024
# Achal – can you please check both lists and see how many have genetic data in SJLIFE and CCSS? Exclude participants from CCSS if they are already in SJLIFE.
ccss <- read_excel("Z:/ResearchHome/ClusterHome/aneupane/data/request/cindy/colorectal/CCSS_colorectal_patients_2024-06-05.xlsx", trim_ws = TRUE)
ccss$CCSSID <- trimws(ccss$CCSSID)
sjlife <- read_excel("Z:/ResearchHome/ClusterHome/aneupane/data/request/cindy/colorectal/SJLIFE_colorectal_patients_2024-06-05.xlsx", trim_ws = TRUE)
overlaps <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/gwas/pheno/SJLIFE_WGS_samples_overlap_with_CCSS_org_SNP_samples.txt", header = T)

## Remove SJLIFE and CCSS overlaps from CCSS
ccss$SJLID_overlaps <- overlaps$V1[match(ccss$CCSSID,overlaps$ccssid)]
dim(ccss)
# 147

ccss$overlap_in_SJLIFE <- ccss$SJLID_overlaps %in% sjlife$studyid

# Check how many have genetic data

sjlife_fam <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/MERGED_sjlife1_2_PreQC/cleaned/MERGED.SJLIFE.1.2.GATKv3.4.VQSR.chr7.preQC_biallelic_renamed_ID_edited.vcf.gz.fam", header = F)
ccss_org <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/ccss_org_hrc/plink/CCSS_org_GRCh37_chr22.fam", header = F)
ccss_exp <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/ccss_exp_wgs/preQC_VCF_per_chromosome/CCSS_exp_biallelic_chr22_ID_edited.fam", header = F)
WES <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/Survivor_WES/biallelic/plink_all/chr.ALL.Survivor_WES.GATK4180.hg38_biallelic_ID_updated.fam", header = F)
library(stringr)

# Split the values in column v1 by "_" and keep the first item
ccss_org$ccssid <- str_split(ccss_org$V1, "_", simplify = TRUE)[, 1]

ccss$In_ccss_org <- ifelse(ccss$CCSSID %in% ccss_org$ccssid, 1,0)
ccss$In_ccss_exp <- ifelse(ccss$CCSSID %in% ccss_exp$V1, 1,0)
ccss$In_ccss_wes <- ifelse(ccss$CCSSID %in% WES$V2, 1,0)

write.table(ccss, "Z:/ResearchHome/ClusterHome/aneupane/data/request/cindy/colorectal/CCSS_colorectal_patients_2024-06-05_AN.txt", col.names = T, row.names = F, sep = "\t")

sjlife$in_wgs <- ifelse(sjlife$studyid %in% sjlife_fam$V2, 1,0)
sjlife$in_wes <- ifelse(sjlife$studyid %in% WES$V2,1,0)

write.table(sjlife, "Z:/ResearchHome/ClusterHome/aneupane/data/request/cindy/colorectal/SJLIFE_colorectal_patients_2024-06-05_AN.txt", col.names = T, row.names = F, sep = "\t")
