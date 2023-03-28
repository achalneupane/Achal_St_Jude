load("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/PHENOTYPE/00.CCSS_org_Genetic_data_P_LP_v14.Rdata")
ccss_org.PHENO.ANY_SN <- PHENO.ANY_SN
ccss_org.subneo <- subneo

load("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/PHENOTYPE/00.CCSS_exp_Genetic_data_P_LP_v14.Rdata")
ccss_exp.PHENO.ANY_SN <- PHENO.ANY_SN
ccss_exp.subneo <- subneo

sum(colnames(ccss_org.PHENO.ANY_SN) != colnames(ccss_exp.PHENO.ANY_SN))
sum(colnames(ccss_org.subneo) != colnames(ccss_exp.subneo))

PHENO.ANY_SN <- rbind.data.frame(ccss_org.PHENO.ANY_SN, ccss_exp.PHENO.ANY_SN)
subneo <- rbind.data.frame(ccss_org.subneo, ccss_exp.subneo)
