all_cancers <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/ccss_org_hrc/ccss_org_hrc_vcf_GRCh38/all_cancer_GrCh37.bed", sep='\t', header=F)
# Extract only those that are needed

all_cancers$KEY_GRCh37 <- paste0(all_cancers$V1, ":", all_cancers$V3)
all_cancers$GRCh38_POS <- sapply(strsplit(all_cancers$V4, "\\-"), `[`, 2)


all_cancers$KEY_GRCh38 <- paste0(all_cancers$V1, ":", all_cancers$GRCh38_POS)


bim_file <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/ccss_org_hrc/ccss_org_hrc_vcf_GRCh38/attr_fraction/prs/plink_data/merged.dat.bim", sep = '\t', header = F)
head(bim_file)
bim_file$KEY_GRCh37 <- paste0("chr", bim_file$V1, ":", bim_file$V4)


bim_file$GRCh38_POS <- all_cancers$GRCh38_POS[match(bim_file$KEY_GRCh37, all_cancers$KEY_GRCh37)]
sum(is.na(bim_file$GRCh38_POS))

bim_file$KEY_GRCh38 <- paste0('chr', bim_file$V1, ":", bim_file$GRCh38_POS)

bim_file <- bim_file[!is.na(bim_file$KEY_GRCh38),]

bim_file$ALLELES <- paste0(sapply(strsplit(bim_file$V2, "\\:"), `[`, 3), ":", sapply(strsplit(bim_file$V2, "\\:"), `[`, 4))

bim_file$GRCh38ID <- paste0(bim_file$KEY_GRCh38, ":", bim_file$ALLELES)


all_cancers_38 <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/ccss_exp_wgs/attr_fraction/prs/all_cancer.txt", sep='\t', header=T)
all_cancers_38$KEY_GRCh38 <- paste0('chr',all_cancers_38$CHROM, ":", all_cancers_38$POS_GRCh38)


# Find out which ones are missing
sum(all_cancers_38$KEY_GRCh38 %in% bim_file$KEY_GRCh38)
all_cancers_38.Not_matched <- all_cancers_38[!all_cancers_38$KEY_GRCh38 %in% bim_file$KEY_GRCh38,]

## Matched
all_cancers_38.matched <- all_cancers_38[all_cancers_38$KEY_GRCh38 %in% bim_file$KEY_GRCh38,]
all_cancers_38.matched$KEY_GRCh37 <- gsub("chr", "", all_cancers$KEY_GRCh37[match(all_cancers_38.matched$KEY_GRCh38, all_cancers$KEY_GRCh38)])

## Not matched or Not found in GWAS data
# Check all_cancers_38.Not_matched in all_cancers (GRCh37 version bed file)
all_cancers_38.Not_matched$GRCh37_POS <- all_cancers$V3[match(all_cancers_38.Not_matched$KEY_GRCh38, all_cancers$KEY_GRCh38)]
all_cancers_38.Not_matched$KEY_GRCh37 <- gsub("chr", "", all_cancers$KEY_GRCh37[match(all_cancers_38.Not_matched$KEY_GRCh38, all_cancers$KEY_GRCh38)])
all_cancers_38.Not_matched$BED_GRCh37 <- paste0(all_cancers_38.Not_matched$CHROM,":", all_cancers_38.Not_matched$GRCh37_POS-10,"-", all_cancers_38.Not_matched$GRCh37_POS+10) 

## Create bed file with 10 bp up/downstream
BED.NOT.FOUND <- cbind.data.frame(all_cancers_38.Not_matched$CHROM, all_cancers_38.Not_matched$GRCh37_POS-10, all_cancers_38.Not_matched$GRCh37_POS+10) 


write.table(all_cancers_38.Not_matched, "Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/ccss_org_hrc/ccss_org_hrc_vcf_GRCh38/attr_fraction/prs/all_cancers_38.Not_matched.txt", sep = "\t", col.names = T, quote = F, row.names = F)
write.table(BED.NOT.FOUND, "Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/ccss_org_hrc/ccss_org_hrc_vcf_GRCh38/attr_fraction/prs/all_cancers_38.Not_matched.bed", sep = "\t", col.names = T, quote = F, row.names = F)


## replace GRCh38 with GRCh37 in all_cancer.txt
all_cancers_38$POS_GRCh37 <- all_cancers$V3[match(all_cancers_38$KEY_GRCh38, all_cancers$KEY_GRCh38)]

colnames(all_cancers_38)

all_cancers_37 <- all_cancers_38[c("CHROM", "POS_GRCh37", "REF", "Effect_allele", "Effect_size", "TYPE", "Cancer", "Significant_YN")]

write.table(all_cancers_37, "Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/ccss_org_hrc/ccss_org_hrc_vcf_GRCh38/attr_fraction/prs/all_cancer.txt", sep = "\t", col.names = T, quote = F, row.names = F)

##########################################################################################
## Part 2: re-calculating PRS in SJLIFE with the variants present in CCSS_original data ##
##########################################################################################
## Extract matching variants from SJLIFE
SJLIFE_vars <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/prs/ALL_Cancers_PRS_data.txt", header = T)
SJLIFE_vars$KEY_GRCh38 <- paste0("chr", SJLIFE_vars$CHROM,":", SJLIFE_vars$POS_GRCh38)
all_cancers_38.matched$KEY_GRCh38 

SJLIFE_vars.in.CCSS_org <- SJLIFE_vars[SJLIFE_vars$KEY_GRCh38 %in% all_cancers_38.matched$KEY_GRCh38 ,]
## Add study identifier, so it's easy to identify the PRS files later on
SJLIFE_vars.in.CCSS_org$TYPE <- paste0(SJLIFE_vars.in.CCSS_org$TYPE, "_from_variants_also_in_CCSS_org")

write.table(SJLIFE_vars.in.CCSS_org[1:8], "Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/prs/ALL_Cancers_PRS_data_in_CCSS_org.txt", sep = "\t", col.names = T, quote = F, row.names = F)


##########################################################################################
## Part 3: re-calculating PRS in SJLIFE with the variants present in CCSS_original data ##
##########################################################################################
## Extract matching variants from SJLIFE
CCSS_exp_vars <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/ccss_exp_wgs/attr_fraction/prs/all_cancer.txt", header = T)
CCSS_exp_vars$KEY_GRCh38 <- paste0("chr", CCSS_exp_vars$CHROM,":", CCSS_exp_vars$POS_GRCh38)
all_cancers_38.matched$KEY_GRCh38 

CCSS_exp_vars.in.CCSS_org <- CCSS_exp_vars[CCSS_exp_vars$KEY_GRCh38 %in% all_cancers_38.matched$KEY_GRCh38 ,]
## Add study identifier, so it's easy to identify the PRS files later on
CCSS_exp_vars.in.CCSS_org$TYPE <- paste0(CCSS_exp_vars.in.CCSS_org$TYPE, "_from_variants_also_in_CCSS_org")

write.table(CCSS_exp_vars.in.CCSS_org[1:8], "Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/ccss_exp_wgs/attr_fraction/prs/ALL_Cancers_PRS_data_in_CCSS_org.txt", sep = "\t", col.names = T, quote = F, row.names = F)
