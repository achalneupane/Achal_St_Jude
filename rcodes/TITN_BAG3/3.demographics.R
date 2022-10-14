getwd()

setwd("Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3/pheno")

sjlife <- read.table("sjlife_ttn_bag3.pheno", header = T)
head(sjlife)




ccss_org <- read.table("ccss_org_eur_cardiotoxic_exposed.pheno", header = T)
head(ccss_org)




ccss_exp <- read.table("ccss_exp_eur_cardiotoxic_exposed.pheno", header = T)
head(ccss_exp)
