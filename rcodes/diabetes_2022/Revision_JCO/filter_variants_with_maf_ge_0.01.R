setwd("/research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/diabetes/Results_final_to_Cindy/Revision_task_JCO")

library(data.table)

#########
## EUR ##
#########
EUR <- fread("chr1_22_PCA_eur.assoc.logistic.clean.Psorted.withBETA_deposition", header = T)

EUR.frq <- fread("EUR_frq_chr_ALL.frq", header = T)
EUR.frq <- EUR.frq[EUR.frq$MAF >= 0.01,]

EUR$KEY <- paste0(EUR$SNP, ":", EUR$A1)
EUR.frq$KEY <- paste0(EUR.frq$SNP, ":", EUR.frq$A1)
EUR <- EUR[EUR$KEY %in% EUR.frq$KEY,]
EUR$effect_allele_frequency <- EUR.frq$MAF[match(EUR$KEY, EUR.frq$KEY)]

EUR$beta.lower_CI <- EUR$BETA - (1.96 * EUR$SE)
EUR$beta.upper_CI <- EUR$BETA + (1.96 * EUR$SE)


EUR$A1[grepl(";", EUR$SNP)] <- sub(";.*", "", EUR$A1[grepl(";", EUR$SNP)])
EUR$A2[grepl(";", EUR$SNP)] <- sub(";.*", "", EUR$A2[grepl(";", EUR$SNP)])

EUR$SNP_removed <- EUR$SNP
EUR$SNP[grepl(";", EUR$SNP)] <- paste(EUR$CHR, EUR$BP, EUR$A2, EUR$A1, sep =":")[grepl(";", EUR$SNP)]

# 13_40983974_A_G
EUR$variant_id <- paste(EUR$CHR, EUR$BP, EUR$A2, EUR$A1, sep = "_")

cc <- EUR[grepl(";", EUR$SNP_removed),]


EUR <- EUR[,c("CHR", "BP", "A1", "A2", "BETA", "SE", "effect_allele_frequency", "P", "variant_id", "beta.upper_CI", "beta.lower_CI", "OR")]
colnames(EUR) <- c("chromosome", "base_pair_location", "effect_allele", "other_allele", "beta", 
                   "standard_error", "effect_allele_frequency", "p_value", "variant_id", "ci_upper", "ci_lower", "odds_ratio")

write.table(EUR, "chr1_22_PCA_eur.assoc.logistic.clean.Psorted.withBETA_deposition.maf.ge.0.01.tsv", col.names = T, quote = F, row.names = F, sep = "\t")

#########
## AFR ##
#########
AFR <- fread("chr1_22_PCA_afr.assoc.logistic.clean.Psorted.withBETA_deposition", header = T)

AFR.frq <- fread("AFR_frq_chr_ALL.frq", header = T)
AFR.frq <- AFR.frq[AFR.frq$MAF >= 0.01,]

AFR$KEY <- paste0(AFR$SNP, ":", AFR$A1)
AFR.frq$KEY <- paste0(AFR.frq$SNP, ":", AFR.frq$A1)
AFR <- AFR[AFR$KEY %in% AFR.frq$KEY,]

AFR$effect_allele_frequency <- AFR.frq$MAF[match(AFR$KEY, AFR.frq$KEY)]

AFR$beta.lower_CI <- AFR$BETA - (1.96 * AFR$SE)
AFR$beta.upper_CI <- AFR$BETA + (1.96 * AFR$SE)

AFR$A1[grepl(";", AFR$SNP)] <- sub(";.*", "", AFR$A1[grepl(";", AFR$SNP)])
AFR$A2[grepl(";", AFR$SNP)] <- sub(";.*", "", AFR$A2[grepl(";", AFR$SNP)])

AFR$SNP_removed <- AFR$SNP
AFR$SNP[grepl(";", AFR$SNP)] <- paste(AFR$CHR, AFR$BP, AFR$A2, AFR$A1, sep =":")[grepl(";", AFR$SNP)]

# 13_40983974_A_G
AFR$variant_id <- paste(AFR$CHR, AFR$BP, AFR$A2, AFR$A1, sep = "_")

cc <- AFR[grepl(";", AFR$SNP_removed),]


AFR <- AFR[,c("CHR", "BP", "A1", "A2", "BETA", "SE", "effect_allele_frequency", "P", "variant_id", "beta.upper_CI", "beta.lower_CI", "OR")]
colnames(AFR) <- c("chromosome", "base_pair_location", "effect_allele", "other_allele", "beta", 
                   "standard_error", "effect_allele_frequency", "p_value", "variant_id", "ci_upper", "ci_lower", "odds_ratio")


write.table(AFR, "chr1_22_PCA_afr.assoc.logistic.clean.Psorted.withBETA_deposition.maf.ge.0.01.tsv", col.names = T, quote = F, row.names = F, sep = "\t")

