library(data.table)

# load("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/MERGED_sjlife1_2_PreQC/cleaned/annotation/snpEff/SNPEFF_clinvar_metaSVM_LoF_from_R_filtering_process_PreQC_VCF.RData")
load("/research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/MERGED_sjlife1_2_PreQC/cleaned/annotation/snpEff/SNPEFF_clinvar_metaSVM_LoF_from_R_filtering_process_PreQC_VCF.RData")


splice.region <- LoF.unique[grepl("splice_region_variant", LoF.unique$`ANN[*].EFFECT`),]
dim(splice.region)
# 242541     25
non.splice.region <- LoF.unique[!grepl("splice_region_variant", LoF.unique$`ANN[*].EFFECT`),]
dim(non.splice.region)
# 73845    25



NO.splice.region <- splice.region[grepl("frameshift|stop|donor|acceptor", splice.region$`ANN[*].EFFECT`),]
dim(NO.splice.region)
# 4421   25



LOF.keep <- rbind.data.frame(non.splice.region, NO.splice.region)
dim(LOF.keep)
# 78266    25



Predicted.vars.inVCF.Unique <-  rbind.data.frame(CLINVAR.unique, LOF.keep)

# Predicted.vars.inVCF.Unique$KEY 


# gnomAD <- fread("Z:/ResearchHome/Groups/sapkogrp/projects//Genomics/common/diabetes/gnomAD_maf_lt_0.01_chrALL.txt")
gnomAD <- fread("/research_jude/rgs01_jude/groups/sapkogrp/projects//Genomics/common/diabetes/gnomAD_maf_lt_0.01_chrALL.txt")

head(gnomAD)

gnomAD$KEY <- paste(gnomAD$Otherinfo4, gnomAD$Otherinfo5, gnomAD$Otherinfo7, gnomAD$Otherinfo8, sep = ":")

Predicted.vars.inVCF.Unique.gnomad <- cbind.data.frame (gnomAD[match(Predicted.vars.inVCF.Unique$KEY, gnomAD$KEY), c("gnomAD_genome_ALL", "gnomAD_genome_AFR", "gnomAD_genome_NFE")], Predicted.vars.inVCF.Unique)


# write_rds(Predicted.vars.inVCF.Unique.gnomad, "Predicted.vars.inVCF.Unique.gnomad.rds")

setwd("Z:/ResearchHome/Groups/sapkogrp/projects//Genomics/common/diabetes")

df <- readRDS("Predicted.vars.inVCF.Unique.gnomad.rds")
