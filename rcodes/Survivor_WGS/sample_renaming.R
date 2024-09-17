sampleID <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/Survivor_WGS_QCed/Phenotypes/wgs_sampleIDs_SJLIDmatch_23March2023.txt", header = T, sep = "\t")
vcfID <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/Survivor_WGS_QCed/VCFsample_names.txt", header = F, sep = "\t")

sum(sampleID$CompBioID_WGS %in% vcfID$V1)

vcfID$cleaned <- sub("-.*", "", vcfID$V1)
sum(sampleID$CompBioID_WGS %in% vcfID$cleaned)
sampleID$match <- sampleID$CompBioID_WGS %in% vcfID$cleaned

sampleID$VCFID <- vcfID$V1[match(sampleID$CompBioID_WGS, vcfID$cleaned)]

sum(duplicated(vcfID$cleaned))
vcfID$cleaned[duplicated(vcfID$cleaned)]

sampleID$CCSSID <- sub(".*CCSS-", "", sampleID$VCFID) 
sampleID$CCSSID[!grepl("CCSS",sampleID$VCFID)] <- NA

sampleID$VCFrename <- sampleID$SJLID
sampleID$VCFrename[sampleID$VCFrename==""] <- sampleID$CCSSID[sampleID$VCFrename==""] 

sampleID <- sampleID[c("CompBioID_WGS", "ERF_note", "SJLID", "CCSSID", "VCFID", "VCFrename")]

sampleID$dups <- duplicated(sampleID$VCFrename)
# sampleID$VCFrename2 <- ave(sampleID$VCFrename, sampleID$VCFrename, 
#                           FUN = function(x) ifelse(length(x) > 1, 
#                                                    ifelse(seq_along(x) == 1, x, paste0(x, "_", seq_along(x) - 1)), 
#                                                    x))
# 
# sampleID$VCFrename2 <- make.unique(as.character(sampleID$VCFrename), sep = "_")



write.table(sampleID, "Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/Survivor_WGS_QCed/Phenotypes/wgs_sampleIDs_cleaned_AN.txt", col.names = T, row.names = F, sep ="\t", quote = F)
# write.table(sampleID[c("VCFID", "VCFrename")], "Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/Survivor_WGS_QCed/Phenotypes/master_vcf_rename.txt", col.names = T, row.names = F, sep ="\t", quote = F)


## Recheck
sampleID <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/Survivor_WGS_QCed/Phenotypes/wgs_sampleIDs_cleaned_AN_recheck.txt", header = T, sep = "\t")
write.table(sampleID[c("VCFID", "VCFrename")], "Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/Survivor_WGS_QCed/Phenotypes/master_vcf_rename.txt", col.names = T, row.names = F, sep ="\t", quote = F)
