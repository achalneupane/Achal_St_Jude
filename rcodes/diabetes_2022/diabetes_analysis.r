# Main folder: /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/diabetes

load("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/MERGED_sjlife1_2_PreQC/cleaned/annotation/snpEff/SNPEFF_clinvar_metaSVM_LoF_from_R_filtering_process_PreQC_VCF.RData")
# load("/research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/MERGED_sjlife1_2_PreQC/cleaned/annotation/snpEff/SNPEFF_clinvar_metaSVM_LoF_from_R_filtering_process_PreQC_VCF.RData")




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

# remove any missense that is not coming from clinvar
Predicted.vars.inVCF.Unique <-  rbind.data.frame(CLINVAR.unique, LOF.keep)

FINAL <- {}

CHROM <- 1:22
for (i in 1:length(CHROM)){
print(paste0("DOING chr", CHROM[i]))
## Extract gnomAD MAF
annovar.chr <- read.delim(paste0("Z:/ResearchHome/Groups/sapkogrp/projects//Genomics/common/diabetes/gnomAD_maf_lt_0.01_chr",CHROM[i],".txt"), sep = "\t", header = T)
# annovar.chr <- read.delim(paste0("/research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/diabetes/gnomAD_maf_lt_0.01_chr",CHROM[i],".txt"), sep = "\t", stringsAsFactors = F, header = T)

annovar.chr$KEY <- paste(annovar.chr$Otherinfo4, annovar.chr$Otherinfo5, annovar.chr$Otherinfo7, annovar.chr$Otherinfo8, sep = ":")
# annovar.chr$KEY.chr.pos <- paste(annovar.chr$Otherinfo4, annovar.chr$Otherinfo5, sep = ":")
# Predicted.vars.inVCF.Unique$KEY.chr.pos <- paste(Predicted.vars.inVCF.Unique$CHROM, Predicted.vars.inVCF.Unique$POS, sep = ":")

sum(Predicted.vars.inVCF.Unique$KEY %in% annovar.chr$KEY)
# sum(Predicted.vars.inVCF.Unique$KEY.chr.pos %in% annovar.chr$KEY.chr.pos)

# gnomAD ALL
Predicted.vars.inVCF.Unique$gnomAD_genome_ALL <- annovar.chr$gnomAD_genome_ALL [match(Predicted.vars.inVCF.Unique$KEY, annovar.chr$KEY)]
Predicted.vars.inVCF.Unique$gnomAD_genome_NFE <- annovar.chr$gnomAD_genome_NFE [match(Predicted.vars.inVCF.Unique$KEY, annovar.chr$KEY)]
Predicted.vars.inVCF.Unique$gnomAD_genome_AFR <- annovar.chr$gnomAD_genome_AFR [match(Predicted.vars.inVCF.Unique$KEY, annovar.chr$KEY)]

wanted.vars.all <- Predicted.vars.inVCF.Unique[!is.na(Predicted.vars.inVCF.Unique$gnomAD_genome_ALL),]
wanted.vars.all <- wanted.vars.all[wanted.vars.all$gnomAD_genome_ALL < 0.01,]

FINAL <- rbind.data.frame(FINAL, wanted.vars.all)
print(paste0("NROWS of FINAL ", nrow(FINAL)))
}

# Keep the unique variants only
FINAL <- FINAL[!duplicated(FINAL$KEY),]

rm (annovar.chr)


write.table("Z:/ResearchHome/Groups/sapkogrp/projects//Genomics/common/diabetes/extract_variants_gnomAD_ALL_NFE_lt_0.01")

# gnomAD NFE
Predicted.vars.inVCF.Unique$gnomAD_genome_NFE <- annovar.chr$gnomAD_genome_NFE [match(Predicted.vars.inVCF.Unique$KEY, annovar.chr$KEY)]

# gnomAD AFR
Predicted.vars.inVCF.Unique$gnomAD_genome_AFR <- annovar.chr$gnomAD_genome_AFR [match(Predicted.vars.inVCF.Unique$KEY, annovar.chr$KEY)]



sum(Predicted.vars.inVCF.Unique$KEY.chr.pos %in% annovar.chr$KEY.chr.pos)
