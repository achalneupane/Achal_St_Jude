CCSS.samples <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/ccss_exp_wgs/CCSS_samples_Mingjuan.txt", header = F)
CCSS.samples$Sample_Num <- 1:nrow(CCSS.samples)
colnames(CCSS.samples) <- c("CCSSID", "SNUM")
# write.table(CCSS.samples, "Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/ccss_exp_wgs/CCSS_samples_Mingjuan_withSnum.txt", col.names = T, row.names = F, sep ="\t", quote = F)


ccss_org <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/ccss_org_hrc/CCSS_samples.txt", header = F)
ccss_exp <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/ccss_exp_wgs/QCed_samples", header = F)
merged.dat <- read.table("sjlife_ccss_org_ccss_exp_ttn_bag3_kendrick.pheno", header = T) # new data from Kendrick
table(merged.dat$IID %in% ccss_org$V1)
table(merged.dat$IID %in% ccss_exp$V1)

ccss_org_phenotype <- merged.dat[merged.dat$IID %in% ccss_org$V1,]
# write.table(ccss_org_phenotype[-1], "Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3/pheno/CCSS_org_ttn_Mingjuan.pheno", col.names = T, row.names = F, sep ="\t", quote = F)
