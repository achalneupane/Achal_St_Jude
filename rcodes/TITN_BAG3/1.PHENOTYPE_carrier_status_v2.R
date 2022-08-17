## Read phenotype files
pheno.ccss_exp_eur <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3/pheno/ccss_exp_eur_cardiotoxic_exposed.pheno", header = T)
# pheno.ccss_org_eur <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3/pheno/ccss_org_eur_cardiotoxic_exposed.pheno", header = T)
pheno.sjlife_ttn_bag3 <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3/pheno/sjlife_ttn_bag3.pheno", header = T)

#########################
## read carrier status ##
#########################
setwd("Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3/ALL_P_LP_combinations")

## On 08/15/2022, we decided to remove splice_region_variants. 
BAG3.splice_region.remove <-  read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/MERGED_sjlife1_2_PreQC/cleaned/annotation/snpEff/TTN_BAG3/BAG3.1.per.maf.gnomad.ALL.and.NFE.txt", header = T, sep = "\t")
## Remove Meta SVM and REVEL
BAG3.splice_region.remove <- BAG3.splice_region.remove[BAG3.splice_region.remove$CLINVAR == "Y" | BAG3.splice_region.remove$LoF == "Y",]


TITN.splice_region.remove <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/MERGED_sjlife1_2_PreQC/cleaned/annotation/snpEff/TTN_BAG3/TITN.1.per.maf.gnomad.ALL.and.NFE.txt", header = T, sep = "\t")
TITN.splice_region.remove <- TITN.splice_region.remove[TITN.splice_region.remove$CLINVAR == "Y" | TITN.splice_region.remove$LoF == "Y",]

splice_region.remove <- rbind.data.frame(BAG3.splice_region.remove, TITN.splice_region.remove)

NO.splice.region <- splice_region.remove[!grepl("splice_region_variant", splice_region.remove$ANN....EFFECT),]
splice.region <- splice_region.remove[grepl("splice_region_variant", splice_region.remove$ANN....EFFECT),]

splice.region.keep <- splice.region[grepl("frameshift|stop|donor|acceptor", splice.region$ANN....EFFECT),]
splice.region.keep <- rbind.data.frame(NO.splice.region, splice.region.keep)
# remove any missense that is not coming from clinvar
splice.region.keep <- splice.region.keep[!(splice.region.keep$CLINVAR == "N" & splice.region.keep$ANN....EFFECT == "missense_variant"),]

table(splice.region.keep$Annovar_ExonicFunc.refGene)
# .  frameshift deletion frameshift insertion    nonsynonymous SNV             stopgain       synonymous SNV 
# 9                    2                    1                    1                    8                    2 

cc <- cbind.data.frame(SNPId=splice.region.keep$KEY, SnpEff_annotation = splice.region.keep$ANN....EFFECT, Annovar_annotation = splice.region.keep$Annovar_ExonicFunc.refGene, splice.region.keep$CLINVAR)
write.table(cc, "Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3/ALL_P_LP_combinations/overlapping_sjlife_ccss_exp_maf_0.01.Clinvar.LoF.SNP_list.txt", sep = "\t", col.names = T, row.names = F, quote = F)

################
## CCSS::TITN ##
################
TITN_ccss_exp.Clinv.LoF <- read.table("TITN_ccss_exp.0.01maf_overlapping.clinvar.lof_recodeA.raw", header = T)
colnames(TITN_ccss_exp.Clinv.LoF)[-c(1:6)] <- str_split(gsub("\\.", ":", colnames(TITN_ccss_exp.Clinv.LoF)[-c(1:6)]), "_", simplify=T)[,1]
# Remove splice region that were included initially
TITN_ccss_exp.Clinv.LoF <- TITN_ccss_exp.Clinv.LoF[c(1:6,which(colnames(TITN_ccss_exp.Clinv.LoF) %in% splice.region.keep$KEY))]

TITN_ccss_exp.Clinv.LoF$Non.REF.Counts <- rowSums(TITN_ccss_exp.Clinv.LoF[!colnames(TITN_ccss_exp.Clinv.LoF) %in% c("FID", "IID", "PAT", "MAT", "SEX", "PHENOTYPE")])
TITN_ccss_exp.Clinv.LoF$carrierSTATUS <- factor(ifelse(TITN_ccss_exp.Clinv.LoF$Non.REF.Counts >= 1, "Y", "N"))
CCSS.carrier <- cbind.data.frame(IID = TITN_ccss_exp.Clinv.LoF$IID, TITN.Clinv.LoF.Non.REF.Counts = TITN_ccss_exp.Clinv.LoF$Non.REF.Counts, TITN.Clinv.LoF.carrierSTATUS = TITN_ccss_exp.Clinv.LoF$carrierSTATUS)

# TITN_ccss_exp.Clinv.LoF.REVEL <- read.table("TITN_ccss_exp.0.01maf_overlapping.clinvar.lof.REVEL_recodeA.raw", header = T)
# colnames(TITN_ccss_exp.Clinv.LoF.REVEL)[-c(1:6)] <- str_split(gsub("\\.", ":", colnames(TITN_ccss_exp.Clinv.LoF.REVEL)[-c(1:6)]), "_", simplify=T)[,1]
# # Remove splice region that were included initially
# TITN_ccss_exp.Clinv.LoF.REVEL <- TITN_ccss_exp.Clinv.LoF.REVEL[c(1:6,which(colnames(TITN_ccss_exp.Clinv.LoF.REVEL) %in% splice.region.keep$KEY))]
# 
# CCSS.carrier$TITN.Clinv.LoF.REVEL.Non.REF.Counts <- rowSums(TITN_ccss_exp.Clinv.LoF.REVEL[!colnames(TITN_ccss_exp.Clinv.LoF.REVEL) %in% c("FID", "IID", "PAT", "MAT", "SEX", "PHENOTYPE")])
# CCSS.carrier$TITN.Clinv.LoF.REVEL.carrierSTATUS <- factor(ifelse(CCSS.carrier$TITN.Clinv.LoF.REVEL.Non.REF.Counts >= 1, "Y", "N"))

## CCSS::BAG3
BAG3_ccss_exp.Clinv.LoF <- read.table("BAG3_ccss_exp.0.01maf_overlapping.clinvar.lof_recodeA.raw", header = T)
colnames(BAG3_ccss_exp.Clinv.LoF)[-c(1:6)] <- str_split(gsub("\\.", ":", colnames(BAG3_ccss_exp.Clinv.LoF)[-c(1:6)]), "_", simplify=T)[,1]
# Remove splice region that were included initially
BAG3_ccss_exp.Clinv.LoF <- BAG3_ccss_exp.Clinv.LoF[c(1:6,which(colnames(BAG3_ccss_exp.Clinv.LoF) %in% splice.region.keep$KEY))]


CCSS.carrier$BAG3.Clinv.LoF.Non.REF.Counts <- rowSums(BAG3_ccss_exp.Clinv.LoF[!colnames(BAG3_ccss_exp.Clinv.LoF) %in% c("FID", "IID", "PAT", "MAT", "SEX", "PHENOTYPE")])
CCSS.carrier$BAG3.Clinv.LoF.carrierSTATUS <- factor(ifelse(CCSS.carrier$BAG3.Clinv.LoF.Non.REF.Counts >= 1, "Y", "N"))


# BAG3_ccss_exp.Clinv.LoF.REVEL <- read.table("BAG3_ccss_exp.0.01maf_overlapping.clinvar.lof.REVEL_recodeA.raw", header = T)
# colnames(BAG3_ccss_exp.Clinv.LoF.REVEL)[-c(1:6)] <- str_split(gsub("\\.", ":", colnames(BAG3_ccss_exp.Clinv.LoF.REVEL)[-c(1:6)]), "_", simplify=T)[,1]
# # Remove splice region that were included initially
# BAG3_ccss_exp.Clinv.LoF.REVEL <- BAG3_ccss_exp.Clinv.LoF.REVEL[c(1:6,which(colnames(BAG3_ccss_exp.Clinv.LoF.REVEL) %in% splice.region.keep$KEY))]

 
# CCSS.carrier$BAG3.Clinv.LoF.REVEL.Non.REF.Counts <- rowSums(BAG3_ccss_exp.Clinv.LoF.REVEL[!colnames(BAG3_ccss_exp.Clinv.LoF.REVEL) %in% c("FID", "IID", "PAT", "MAT", "SEX", "PHENOTYPE")])
# CCSS.carrier$BAG3.Clinv.LoF.REVEL.carrierSTATUS <- factor(ifelse(CCSS.carrier$BAG3.Clinv.LoF.REVEL.Non.REF.Counts >= 1, "Y", "N"))

## Add carrier status to PHENO
pheno.ccss_exp_eur <- cbind.data.frame(pheno.ccss_exp_eur, CCSS.carrier[match(pheno.ccss_exp_eur$IID, CCSS.carrier$IID),-1])

## Add carrier status for TITN + BAG3
pheno.ccss_exp_eur$TITN_BAG3.clin.LoF.Counts <- pheno.ccss_exp_eur$TITN.Clinv.LoF.Non.REF.Counts + pheno.ccss_exp_eur$BAG3.Clinv.LoF.Non.REF.Counts
pheno.ccss_exp_eur$TITN_BAG3.clin.LoF.carrierSTATUS <- factor(ifelse(pheno.ccss_exp_eur$TITN_BAG3.clin.LoF.Counts >= 1, "Y", "N"))

# pheno.ccss_exp_eur$TITN_BAG3.clin.LoF.REVEL.Counts <- pheno.ccss_exp_eur$TITN.Clinv.LoF.REVEL.Non.REF.Counts + pheno.ccss_exp_eur$BAG3.Clinv.LoF.REVEL.Non.REF.Counts
# pheno.ccss_exp_eur$TITN_BAG3.clin.LoF.REVEL.carrierSTATUS <- factor(ifelse(pheno.ccss_exp_eur$TITN_BAG3.clin.LoF.REVEL.Counts >= 1, "Y", "N"))


##################
## SJLIFE::TITN ##
##################
TITN_sjlife.Clinv.LoF <- read.table("TITN_sjlife.0.01maf_overlapping.clinvar.lof_recodeA.raw", header = T)
colnames(TITN_sjlife.Clinv.LoF)[-c(1:6)] <- str_split(gsub("\\.", ":", colnames(TITN_sjlife.Clinv.LoF)[-c(1:6)]), "_", simplify=T)[,1]
# Remove splice region that were included initially
TITN_sjlife.Clinv.LoF <- TITN_sjlife.Clinv.LoF[c(1:6,which(colnames(TITN_sjlife.Clinv.LoF) %in% splice.region.keep$KEY))]


TITN_sjlife.Clinv.LoF$Non.REF.Counts <- rowSums(TITN_sjlife.Clinv.LoF[!colnames(TITN_sjlife.Clinv.LoF) %in% c("FID", "IID", "PAT", "MAT", "SEX", "PHENOTYPE")])
TITN_sjlife.Clinv.LoF$carrierSTATUS <- factor(ifelse(TITN_sjlife.Clinv.LoF$Non.REF.Counts >= 1, "Y", "N"))
SJLIFE.carrier <- cbind.data.frame(IID = TITN_sjlife.Clinv.LoF$IID, TITN.Clinv.LoF.Non.REF.Counts = TITN_sjlife.Clinv.LoF$Non.REF.Counts, TITN.Clinv.LoF.carrierSTATUS = TITN_sjlife.Clinv.LoF$carrierSTATUS)

# TITN_sjlife.Clinv.LoF.REVEL <- read.table("TITN_sjlife.0.01maf_overlapping.clinvar.lof.REVEL_recodeA.raw", header = T)
# colnames(TITN_sjlife.Clinv.LoF.REVEL)[-c(1:6)] <- str_split(gsub("\\.", ":", colnames(TITN_sjlife.Clinv.LoF.REVEL)[-c(1:6)]), "_", simplify=T)[,1]
# # Remove splice region that were included initially
# TITN_sjlife.Clinv.LoF.REVEL <- TITN_sjlife.Clinv.LoF.REVEL[c(1:6,which(colnames(TITN_sjlife.Clinv.LoF.REVEL) %in% splice.region.keep$KEY))]

# SJLIFE.carrier$TITN.Clinv.LoF.REVEL.Non.REF.Counts <- rowSums(TITN_sjlife.Clinv.LoF.REVEL[!colnames(TITN_sjlife.Clinv.LoF.REVEL) %in% c("FID", "IID", "PAT", "MAT", "SEX", "PHENOTYPE")])
# SJLIFE.carrier$TITN.Clinv.LoF.REVEL.carrierSTATUS <- factor(ifelse(SJLIFE.carrier$TITN.Clinv.LoF.REVEL.Non.REF.Counts >= 1, "Y", "N"))

## SJLIFE::BAG3
BAG3_sjlife.Clinv.LoF <- read.table("BAG3_sjlife.0.01maf_overlapping.clinvar.lof_recodeA.raw", header = T)
colnames(BAG3_sjlife.Clinv.LoF)[-c(1:6)] <- str_split(gsub("\\.", ":", colnames(BAG3_sjlife.Clinv.LoF)[-c(1:6)]), "_", simplify=T)[,1]
# Remove splice region that were included initially
BAG3_sjlife.Clinv.LoF <- BAG3_sjlife.Clinv.LoF[c(1:6,which(colnames(BAG3_sjlife.Clinv.LoF) %in% splice.region.keep$KEY))]


SJLIFE.carrier$BAG3.Clinv.LoF.Non.REF.Counts <- rowSums(BAG3_sjlife.Clinv.LoF[!colnames(BAG3_sjlife.Clinv.LoF) %in% c("FID", "IID", "PAT", "MAT", "SEX", "PHENOTYPE")])
SJLIFE.carrier$BAG3.Clinv.LoF.carrierSTATUS <- factor(ifelse(SJLIFE.carrier$BAG3.Clinv.LoF.Non.REF.Counts >= 1, "Y", "N"))


# # BAG3_sjlife.Clinv.LoF.REVEL <- read.table("BAG3_sjlife.0.01maf_overlapping.clinvar.lof.REVEL_recodeA.raw", header = T)
# # colnames(BAG3_sjlife.Clinv.LoF.REVEL)[-c(1:6)] <- str_split(gsub("\\.", ":", colnames(BAG3_sjlife.Clinv.LoF.REVEL)[-c(1:6)]), "_", simplify=T)[,1]
# # # Remove splice region that were included initially
# # BAG3_sjlife.Clinv.LoF.REVEL <- BAG3_sjlife.Clinv.LoF.REVEL[c(1:6,which(colnames(BAG3_sjlife.Clinv.LoF.REVEL) %in% splice.region.keep$KEY))]
# 
# 
# SJLIFE.carrier$BAG3.Clinv.LoF.REVEL.Non.REF.Counts <- rowSums(BAG3_sjlife.Clinv.LoF.REVEL[!colnames(BAG3_sjlife.Clinv.LoF.REVEL) %in% c("FID", "IID", "PAT", "MAT", "SEX", "PHENOTYPE")])
# SJLIFE.carrier$BAG3.Clinv.LoF.REVEL.carrierSTATUS <- factor(ifelse(SJLIFE.carrier$BAG3.Clinv.LoF.REVEL.Non.REF.Counts >= 1, "Y", "N"))

## Add carrier status to PHENO
pheno.sjlife_ttn_bag3 <- cbind.data.frame(pheno.sjlife_ttn_bag3, SJLIFE.carrier[match(pheno.sjlife_ttn_bag3$IID, SJLIFE.carrier$IID),-1])


## Add carrier status for TITN + BAG3
pheno.sjlife_ttn_bag3$TITN_BAG3.clin.LoF.Counts <- pheno.sjlife_ttn_bag3$TITN.Clinv.LoF.Non.REF.Counts + pheno.sjlife_ttn_bag3$BAG3.Clinv.LoF.Non.REF.Counts
pheno.sjlife_ttn_bag3$TITN_BAG3.clin.LoF.carrierSTATUS <- factor(ifelse(pheno.sjlife_ttn_bag3$TITN_BAG3.clin.LoF.Counts >= 1, "Y", "N"))

# pheno.sjlife_ttn_bag3$TITN_BAG3.clin.LoF.REVEL.Counts <- pheno.sjlife_ttn_bag3$TITN.Clinv.LoF.REVEL.Non.REF.Counts + pheno.sjlife_ttn_bag3$BAG3.Clinv.LoF.REVEL.Non.REF.Counts
# pheno.sjlife_ttn_bag3$TITN_BAG3.clin.LoF.REVEL.carrierSTATUS <- factor(ifelse(pheno.sjlife_ttn_bag3$TITN_BAG3.clin.LoF.REVEL.Counts >= 1, "Y", "N"))

# table(TITN_ccss_exp.Clinv.LoF$IID == TITN_ccss_exp.Clinv.LoF.REVEL$IID)
# table(BAG3_ccss_exp.Clinv.LoF$IID == BAG3_ccss_exp.Clinv.LoF.REVEL$IID)
# table(TITN_ccss_exp.Clinv.LoF$IID == BAG3_ccss_exp.Clinv.LoF.REVEL$IID)

# table(TITN_sjlife.Clinv.LoF$IID == TITN_sjlife.Clinv.LoF.REVEL$IID)
# table(BAG3_sjlife.Clinv.LoF$IID == BAG3_sjlife.Clinv.LoF.REVEL$IID)
# table(TITN_sjlife.Clinv.LoF$IID == BAG3_sjlife.Clinv.LoF$IID)

save.image("Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3/analysis_model_fitting/phenotype_carrier_status_v2.RDATA")

