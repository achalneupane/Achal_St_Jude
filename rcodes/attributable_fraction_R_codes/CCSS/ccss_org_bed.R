all_cancers <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/ccss_org_hrc/ccss_org_hrc_vcf_GRCh38/all_cancer_GrCh37.bed", sep='\t', header=F)
# Extract only those that are needed

all_cancers$KEY_GRCh37 <- paste0(all_cancers$V1, ":", all_cancers$V3)
all_cancers$GRCh38_POS <- sapply(strsplit(all_cancers$V4, "\\-"), `[`, 2)


all_cancers$KEY_GRCh38 <- paste0(all_cancers$V1, ":", all_cancers$GRCh38_POS)


bim_file <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/ccss_org_hrc/ccss_org_hrc_vcf_GRCh38/attr_fraction/prs/plink_data/merged.dat.bim", sep = '\t', header = F)
head(bim_file)
bim_file$KEY_GRCh37 <- paste0("chr", bim_file$V1, ":", bim_file$V4)


bim_file$KEY_GRCh38 <- all_cancers$GRCh38_POS[match(bim_file$KEY_GRCh37, all_cancers$KEY_GRCh37)]
sum(is.na(bim_file$KEY_GRCh38))
