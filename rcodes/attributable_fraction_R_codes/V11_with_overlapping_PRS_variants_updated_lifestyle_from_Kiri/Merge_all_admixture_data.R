# Merge_all_admixture_data.R
SJLIFE1 <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/admixture/SJLIFE_and_1kGP_final_EUR_AFR_EAS.3.Q_sjlife1_samples", header = T)
head(SJLIFE1)
colnames(SJLIFE1)[1] <- "INDIVIDUAL"

SJLIFE2 <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/admixture/SJLIFE2_and_1kGP_final_EUR_AFR_EAS.3.Q_sjlife2_samples", header = T)
head(SJLIFE2)
colnames(SJLIFE2)[1] <- "INDIVIDUAL"

CCSS_1 <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/admixture/CCSS_and_1kGP_final_EUR_AFR_EAS.3.Q_ccssexp_samples", header = T)
head(CCSS_1)
colnames(CCSS_1)[1] <- "INDIVIDUAL"

CCSS_2 <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/admixture/CCSS.SJLIFE.ancestry_CCSS.txt", header = T, sep = "\t")
head(CCSS_2)
CCSS_2 <- CCSS_2[c("SAMPLE", "CEU", "YRI", "ASA")]
colnames(CCSS_2) <- c("INDIVIDUAL", "EUR", "AFR", "EAS")

ancestry <- rbind.data.frame(SJLIFE1, SJLIFE2, CCSS_1, CCSS_2)
write.table(ancestry, "Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/admixture/merged.ancestry.file.txt", col.names = T, sep = "\t")
