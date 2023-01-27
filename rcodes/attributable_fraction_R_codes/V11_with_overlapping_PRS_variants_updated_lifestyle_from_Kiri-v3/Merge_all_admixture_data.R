# Merge_all_admixture_data.R
SJLIFE1 <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/admixture/SJLIFE_and_1kGP_final_EUR_AFR_EAS.3.Q_sjlife1_samples", header = T)
head(SJLIFE1)
colnames(SJLIFE1)[1] <- "INDIVIDUAL"

SJLIFE2 <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/admixture/SJLIFE2_and_1kGP_final_EUR_AFR_EAS.3.Q_sjlife2_samples", header = T)
head(SJLIFE2)
colnames(SJLIFE2)[1] <- "INDIVIDUAL"

CCSS_exp <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/admixture/CCSS_and_1kGP_final_EUR_AFR_EAS.3.Q_ccssexp_samples", header = T)
head(CCSS_exp)
colnames(CCSS_exp)[1] <- "INDIVIDUAL"

CCSS_org <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/admixture/CCSS.SJLIFE.ancestry_CCSS.txt", header = T, sep = "\t")
head(CCSS_org)
CCSS_org <- CCSS_org[c("SAMPLE", "CEU", "YRI", "ASA")]
colnames(CCSS_org) <- c("INDIVIDUAL", "EUR", "AFR", "EAS")
CCSS_org$INDIVIDUAL <- paste(CCSS_org$INDIVIDUAL, CCSS_org$INDIVIDUAL, sep = "_")

ancestry <- rbind.data.frame(SJLIFE1, SJLIFE2, CCSS_exp, CCSS_org)
write.table(ancestry, "Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/admixture/merged.ancestry.file.txt", col.names = T, sep = "\t")
