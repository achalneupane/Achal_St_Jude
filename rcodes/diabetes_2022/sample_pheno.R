## EUR
PCA <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics//common/sjlife/MERGED_SJLIFE_1_2/MERGED_SJLIFE_PLINK_PER_CHR/PCA//final_EUR.fam")
EUR <- read.delim("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics//common/diabetes/SJLIFE_T2D_GWAS_EUR.pheno", header = T, sep = " ")
head(EUR)

cc <- na.omit(EUR)

write.table(cc, "Z:/ResearchHome/Groups/sapkogrp/projects/Genomics//common/diabetes/final_data_submission/SJLIFE_T2D_GWAS_EUR_AN.pheno", col.names = T, row.names = F, sep = "\t", quote = F)
write.table(cc[1:2], "Z:/ResearchHome/Groups/sapkogrp/projects/Genomics//common/diabetes/final_data_submission/SJLIFE_T2D_GWAS_EUR_AN.samples.list", col.names = T, row.names = F, sep = "\t", quote = F)

## AFR
PCA <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics//common/sjlife/MERGED_SJLIFE_1_2/MERGED_SJLIFE_PLINK_PER_CHR/PCA//final_AFR.fam")
AFR <- read.delim("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics//common/diabetes/SJLIFE_T2D_GWAS_AFR.pheno", header = T, sep = " ")
head(EUR)

cc <- na.omit(AFR)
dim(cc)

write.table(cc, "Z:/ResearchHome/Groups/sapkogrp/projects/Genomics//common/diabetes/final_data_submission/SJLIFE_T2D_GWAS_AFR.pheno", col.names = T, row.names = F, sep = "\t", quote = F)
write.table(cc[1:2], "Z:/ResearchHome/Groups/sapkogrp/projects/Genomics//common/diabetes/final_data_submission/AFR_samples.list", col.names = T, row.names = F, sep = "\t", quote = F)
