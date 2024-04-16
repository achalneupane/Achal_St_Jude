## Admixture classification
admixture <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common//sjlife/MERGED_SJLIFE_1_2/MERGED_SJLIFE_PLINK_PER_CHR/PCA/SJLIFE_4481_Admixture_PCA_ethnicity.csv", sep = "\t", header = T)
AFR.admix <- admixture$INDIVIDUAL[admixture$AFR > 0.6]
EUR.admix <- admixture$INDIVIDUAL[admixture$EUR > 0.8]
