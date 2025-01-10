sampleID <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/Survivor_WGS_QCed/Phenotypes/wgs_sampleIDs_SJLIDmatch_23March2023.txt", header = T, sep = "\t")
vcfID <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/Survivor_WGS_QCed/VCFsample_names.txt", header = F, sep = "\t")

sum(sampleID$CompBioID_WGS %in% vcfID$V1)
# 2
sampleID$CompBioID_WGS[sampleID$CompBioID_WGS %in% vcfID$V1]
# "SJNPC018728_G1" "SJNPC018729_G1"

vcfID$cleaned <- sub("-.*", "", vcfID$V1)
sum(sampleID$CompBioID_WGS %in% vcfID$cleaned)
sampleID$match <- sampleID$CompBioID_WGS %in% vcfID$cleaned

sampleID$VCFID <- vcfID$V1[match(sampleID$CompBioID_WGS, vcfID$cleaned)]

sum(duplicated(vcfID$cleaned))
vcfID$cleaned[duplicated(vcfID$cleaned)]

sampleID$CCSSID <- NA
sampleID$CCSSID <- sub(".*CCSS-0", "", sampleID$VCFID) 
sampleID$CCSSID <- sub(".*CCSS-", "", sampleID$CCSSID) 
sampleID$CCSSID[!grepl("CCSS",sampleID$VCFID)] <- NA

sampleID$VCFrename <- sampleID$SJLID
# sampleID$VCFrename[sampleID$VCFrename==""] <- sampleID$CCSSID[sampleID$VCFrename==""] 
sampleID$VCFrename[!is.na(sampleID$CCSSID)] <- sampleID$CCSSID[!is.na(sampleID$CCSSID)]
# cc <- sampleID
# pp <- sampleID[(sampleID$VCFrename!=cc$VCFrename),]
# gg <- pp[!pp$CCSSID %in% missing.ccss.exp$renameCCSS_exp,]
# # CompBioID_WGS      SJLID ERF_note match                        VCFID   CCSSID  VCFrename
# # 2577 SJCNS054399_G1 SJL2533606           TRUE SJCNS054399_G1-CCSS-15253362 15253362 SJL2533606
# # 4777  SJHL052240_G1 SJL5187207           TRUE  SJHL052240_G1-CCSS-15480903 15480903 SJL5187207
# # 4787  SJHL052295_G1 SJL2532507           TRUE  SJHL052295_G1-CCSS-15253253 15253253 SJL2532507
# # NA             <NA>       <NA>     <NA>    NA                         <NA>     <NA>       <NA>
# # NA.1           <NA>       <NA>     <NA>    NA                         <NA>     <NA>       <NA>

sampleID <- sampleID[c("CompBioID_WGS", "ERF_note", "SJLID", "CCSSID", "VCFID", "VCFrename")]

sampleID$dups <- duplicated(sampleID$VCFrename)
# sampleID$VCFrename2 <- ave(sampleID$VCFrename, sampleID$VCFrename, 
#                           FUN = function(x) ifelse(length(x) > 1, 
#                                                    ifelse(seq_along(x) == 1, x, paste0(x, "_", seq_along(x) - 1)), 
#                                                    x))
# 
# sampleID$VCFrename2 <- make.unique(as.character(sampleID$VCFrename), sep = "_")



write.table(sampleID, "Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/Survivor_WGS_QCed/Phenotypes/wgs_sampleIDs_cleaned_AN_corrected.txt", col.names = T, row.names = F, sep ="\t", quote = F)
# write.table(sampleID[c("VCFID", "VCFrename")], "Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/Survivor_WGS_QCed/Phenotypes/master_vcf_rename.txt", col.names = T, row.names = F, sep ="\t", quote = F)
cc <- sampleID
gg <- sampleID[match(cc$CompBioID_WGS, sampleID$CompBioID_WGS),]


## Recheck
sampleID <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/Survivor_WGS_QCed/Phenotypes/wgs_sampleIDs_cleaned_AN_recheck.txt", header = T, sep = "\t")
## these samples need to be renamed because I had labelled them as SJLIFE but they are CCSS.
dim(sampleID[grepl("SJL", sampleID$VCFrename) & grepl("CCSS", sampleID$VCFID),])
# 143
samples.relabel <- sampleID[grepl("SJL", sampleID$VCFrename) & grepl("CCSS", sampleID$VCFID),]

sampleID$VCFrename [grepl("CCSS", sampleID$VCFID)] <- sampleID$CCSSID [grepl("CCSS", sampleID$VCFID)] 
write.table(sampleID[c("VCFID", "VCFrename")], "Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/Survivor_WGS_QCed/Phenotypes/master_vcf_rename_corrected.txt", col.names = F, row.names = F, sep ="\t", quote = F)

write.table(samples.relabel[c("VCFrename", "CCSSID")], "Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/Survivor_WGS_QCed/Phenotypes/143_CCSS_rename_corrected.txt", col.names = F, row.names = F, sep ="\t", quote = F)

table(cc$VCFrename == sampleID$VCFrename)
cc[cc$VCFrename != sampleID$VCFrename,]

## Remove duplicate IDs and NOT-SJLIFE samples from qced Plink file
as.data.frame(sampleID[sampleID$dups == "TRUE"| sampleID$VCFrename == "REMOVE", "VCFrename"])
as.data.frame(sampleID[sampleID$dups == "TRUE"| sampleID$VCFrename == "REMOVE", "VCFrename"])
sampleID[sampleID$dups == "TRUE" | sampleID$VCFrename == "REMOVE", "VCFrename"]
## Remove these samples from Final VCF Survivor_WGS.GATK4180.hg38_renamed_chr10.PASS.decomposed.qced.fam
# remove_duplicate_IDs.txt (18 samples)
# SJL4790401 2
# 15253362 15253362
# SJL1656107 2
# 15480903 15480903
# 15253253 15480903
# REMOVE REMOVE
# SJNORM041129 G2-TB-15-6901
# SJNORM041130 G1-TB-15-6225
# SJNORM041130 G2-TB-15-6902
# SJNORM041131 G1-TB-15-6369
# SJNORM041131 G2-TB-15-6903
# SJNORM041132 G1-TB-15-6370
# SJNORM041132 G2-TB-15-6904
# SJNORM041133 G1-TB-15-6399
# SJNORM041133 G2-TB-15-6905
# SJL1656107 3
# SJNPC018729 G1
# SJL2053109 2
#############
## Post QC ##
#############
# Note: all.survivors.WES is preQC data
# all.survivors.WES <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/Survivor_WES/biallelic2/plink_all/chr.ALL.Survivor_WES.GATK4180.hg38_biallelic.geno.0.1.hwe.1e-15.LCR.removed.MAC.ge.1_ID_updated.fam", header = F)
all.survivors.WGS <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/Survivor_WGS_QCed/QC/Survivor_WGS.GATK4180.hg38_renamed_chr1.PASS.decomposed.qced.fam", header = F)
# c.control.WES <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/Survivor_WES/biallelic2/plink_all/Controls/chr.ALL.survivor.control_WES.GATK4180.hg38_biallelic.geno.0.1.hwe.1e-15.LCR.removed.MAC.ge.1_updated.fam", header = F)
# sample.not.found <- rownames(raw)[!rownames(raw) %in% all.survivors.WGS$V2]
# length(sample.not.found)
# ## [1] 219 samples were not found in WGS ..
# all.survivors.WGS$V2[!all.survivors.WGS$V2 %in% all.survivors.WES$V2]

## Check 4401 samples in SJLIFE
load("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/PHENOTYPE/5_lifestyle_v11.RDATA")
table(PHENO.ANY_SN$sjlid %in% all.survivors.WGS$V2)
# FALSE  TRUE 
# 6  4395
PHENO.ANY_SN$sjlid[!PHENO.ANY_SN$sjlid %in% all.survivors.WGS$V2]
# "SJL5052915" "SJL5057815" ## Not in the sample list
# SJOS042758_G1 = "SJL5094513"; SJAML041901_G1="SJL5188202"; SJOS042104_G1="SJL5128013"; SJHL042045_G1="SJL5207307"

AA.90samples <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/RNAseq/common/RNAseq_sjlife/jenn_RNAseq/AA_samples.txt", header = F)
table(AA.90samples$V1 %in% all.survivors.WGS$V2)
# TRUE 
# 90 

## Check all previous SJLIFE
sjlife.4507 <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/MERGED_sjlife1_2_PreQC/cleaned/MERGED.SJLIFE.1.2.GATKv3.4.VQSR.chr9.preQC_biallelic_renamed_ID_edited.vcf.gz.fam", header = F)
sjlife.4481 <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/sjlife/MERGED_SJLIFE_1_2/MERGED_SJLIFE_PLINK_PER_CHR/MERGED.SJLIFE.1.2.GATKv3.4.VQSR.chr2.PASS.decomposed.fam", header = F)
# table(sjlife.4481$V2 %in% AA.90samples$V1)
# table(sjlife.4507$V2 %in% AA.90samples$V1)
table(sjlife.4507$V2 %in% all.survivors.WGS$V2)
# FALSE  TRUE 
# 6  4501
sjlife.4507$V2[!sjlife.4507$V2 %in% all.survivors.WGS$V2]
# "SJL5207307" "SJL5188202" "SJL5128013" "SJL5094513" "SJL5052915" "SJL5057815"

table(rownames(raw) %in% all.survivors.WGS$V2)
# FALSE  TRUE 
# 219  7653

overlaps <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/gwas/pheno/SJLIFE_WGS_samples_overlap_with_CCSS_org_SNP_samples.txt", header = T)

## previously QCed CCSS_exp
ccss_exp <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/ccss_exp_wgs/QCed_samples", header = F)
table(ccss_exp$V1 %in% all.survivors.WGS$V2)
# FALSE  TRUE 
# 143  2696 
to.check <- ccss_exp$V1[!ccss_exp$V1 %in% all.survivors.WGS$V2]
to.check
# [1] 15479822 15417752 15478832 15251382 15476582 15476132 15481642 15478662 15283612 15482242 15250832 15253472 15476722 15480502 15251792 15419262 15283172 15473712 15483012 15417832 15485992
# [22] 15419092 15254562 15484582 15252622 15482792 15419492 15479242 15479772 15484882 15418202 15253362 15475772 15250628 15479965 15486305 15253425 15252735 15478945 15417675 15475785 15478325
# [43] 15479785 15284461 15477231 15479401 15482931 15253231 15473311 15477291 15254851 15472721 15283871 15474431 15417351 15478821 15480371 15250031 15472981 15253311 15477891 15482341 15475651
# [64] 15479021 15284211 15253641 15476541 15476561 15253971 15253291 15477371 15482821 15484781 15250941 15250411 15417511 15482301 15473601 15251551 15474541 15418751 15479131 15474711 15252431
# [85] 15483701 15254021 15283701 15486561 15479861 15252651 15476471 15284291 15250811 15474281 15283921 15481761 15484931 15250561 15480903 15482203 15419703 15478913 15251583 15251563 15283743
# [106] 15474933 15479203 15417493 15251663 15253753 15252753 15252893 15419593 15472713 15253683 15284093 15418073 15475343 15250083 15476263 15417613 15483333 15252926 15485966 15478684 15480654
# [127] 15283244 15480284 15486494 15478144 15475584 15479694 15418964 15476774 15252534 15480934 15475434 15283414 15251894 15250764 15283027 15485077 15251537 15419417 15283287