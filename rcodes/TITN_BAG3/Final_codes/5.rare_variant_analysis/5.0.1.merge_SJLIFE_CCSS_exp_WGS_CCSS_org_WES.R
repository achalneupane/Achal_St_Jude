# One technical concern is the MAF1% cutoff that the authors used for rare P/LP
# variants. MAF cutoff often used are 0.01%. Also, many of rare TTN variants
# reported in Table S4 are in I-band (in addition to having MAF>0.0001), which
# are not generally considered pathogenic/likely pathogenic. Mostly, variants in
# exons with high PSI, often located in A-band are considered P/LP. For example,
# chr2:178665777 A>G in Table S4 is in exon with PSI 1% and won’t be considered
# pathogenic. Therefore, I have concerns about their rare genetic variant
# analysis in the manuscript.
library(data.table)
setwd("Z:/ResearchHome/Groups/sapkogrp/projects//Cardiotoxicity/common/ttn_bag3//rare_variant_sjlife_ccss_combined/rare_variants_from_ccss_org_WES")
dbgapID <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3/pheno/CCSS_org_ttn_Mingjuan_dbGAPID.txt", header = T, sep = "\t")
wes.samples <- fread("original_samples")
wes.samples$V2.cleaned <- sub(".*G1-", "", wes.samples$V2)
wes.samples$V2.cleaned <- gsub("_WES", "", wes.samples$V2.cleaned)

table(wes.samples$V2.cleaned %in% dbgapID$DBGAP_COMPBIO_ID)
sum(dbgapID$DBGAP_COMPBIO_ID %in% wes.samples$V2.cleaned)
# 3001
dbgapID$DBGAP_COMPBIO_ID[!dbgapID$DBGAP_COMPBIO_ID %in% wes.samples$V2.cleaned]
sum(dbgapID$DBGAP_COMPBIO_ID %in% wes.samples$V2.cleaned)

wes.samples <- wes.samples[wes.samples$V2.cleaned %in% dbgapID$DBGAP_COMPBIO_ID,]
update_IDs <- wes.samples[, c("V1", "V2", "V2.cleaned", "V2.cleaned")]

# write.table(update_IDs, "updated_Ids.txt", col.names = F, row.names = F, quote = F, sep = "\t")
# write.table(update_IDs[,3:4], "keep_updated_Ids.txt", col.names = F, row.names = F, quote = F, sep = "\t")

wgs.eur.PAV <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3/rare_variant_sjlife_ccss_combined/snpEFF_EUR_vars_final_0.0001_annotation.txt", header = T, sep = "\t")
wes.vars <- fread("chrALL.SJLIFE_CCSS_WES_101724.GATK4180.hg38_biallelic_updated_3001.bim")

table(wgs.eur.PAV$SNP %in% wes.vars$V2)
# FALSE  TRUE 
# 70   150 

wgs.eur.PAV$SNP[!wgs.eur.PAV$SNP %in% wes.vars$V2]
wgs.eur.PAV$SNP[!wgs.eur.PAV$SNP %in% wes.vars$V2]
# [1] "chr10:119651746:C:G" "chr10:119651808:C:T" "chr10:119669870:A:G" "chr10:119669881:C:T" "chr10:119669985:C:G" "chr10:119669986:C:T" "chr10:119670088:G:A" "chr10:119670130:G:C"
# [9] "chr10:119672297:T:G" "chr10:119676542:C:T" "chr10:119676720:G:A" "chr10:119677132:G:T" "chr1:156114941:G:A"  "chr1:156115268:A:G"  "chr1:156126893:C:T"  "chr1:156129890:C:T" 
# [17] "chr1:156134891:G:A"  "chr1:156134914:C:T"  "chr1:156135260:C:T"  "chr1:156136240:C:T"  "chr1:156137756:C:T"  "chr1:156138720:G:A"  "chr1:156138821:C:T"  "chr1:156138822:G:A" 
# [25] "chr1:156138830:C:G"  "chr1:156138846:C:G"  "chr1:156138860:C:T"  "chr1:156138894:G:T"  "chr1:156138902:C:G"  "chr17:39665792:C:T"  "chr17:39665814:G:A"  "chr17:39665831:C:T" 
# [33] "chr17:39665991:A:G"  "chr17:39666065:C:T"  "chr17:39666078:G:A"  "chr3:52452128:C:T"   "chr3:52452222:A:T"   "chr1:201359060:G:A"  "chr1:201359084:C:G"  "chr1:201359087:G:A" 
# [41] "chr1:201359177:G:A"  "chr1:201361274:C:T"  "chr1:201365656:T:C"  "chr14:23413770:T:A"  "chr14:23413860:G:A"  "chr14:23415083:T:C"  "chr14:23415084:T:C"  "chr14:23415832:C:A" 
# [49] "chr14:23416089:T:C"  "chr14:23416093:G:A"  "chr14:23416182:C:T"  "chr14:23417179:A:T"  "chr14:23418375:G:A"  "chr14:23419522:C:T"  "chr14:23419557:C:T"  "chr14:23419566:T:C" 
# [57] "chr14:23420164:C:T"  "chr14:23420234:C:T"  "chr14:23422315:G:T"  "chr14:23424148:T:C"  "chr14:23430591:A:G"  "chr6:7574797:T:C"    "chr2:178559309:A:T"  "chr2:178575970:G:A" 
# [65] "chr2:178586842:G:A"  "chr2:178675722:T:A"  "chr2:178706956:T:G"  "chr2:178731304:C:T"  "chr2:178750291:G:A"  "chr2:178756220:A:G" 

wgs.eur.PAV$KEY <- paste0(wgs.eur.PAV$CHROM, ":", wgs.eur.PAV$POS)
wes.vars$KEY <- paste0("chr", wes.vars$V1, ":", wes.vars$V4)
table(wgs.eur.PAV$KEY %in% wes.vars$KEY)
# 150 There are exactly 150 variants from WGS in WES. WE can keep these variants as is
wanted.snps <- as.data.frame(wgs.eur.PAV$SNP[wgs.eur.PAV$KEY %in% wes.vars$KEY])

write.table(wanted.snps, "wanted_snps.txt", col.names = F, row.names = F, quote = F, sep = "\t")

dbgapID.wanted <- dbgapID[dbgapID$DBGAP_COMPBIO_ID %in% wes.samples$V2.cleaned, ]
write.table(dbgapID.wanted, "Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3/pheno/ccss_org_mingjuan_pheno_found_in_WES.txt", col.names = T, row.names = F, quote = F, sep = "\t")
