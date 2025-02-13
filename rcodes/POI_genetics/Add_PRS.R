# Had to update the IDs after correcting the sample names in WGS
correct_ccss_exp <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/Survivor_WGS_QCed/Phenotypes/143_CCSS_rename_corrected.txt", header = F)
profile <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/POI_genetics/Ruth_et_al/prs/prs_out_wgs_survivor/POI_META_prs.profile", header = T)
# remove duplicates
profile <- profile[!((profile$IID == 3)|(profile$IID == 2)|(profile$IID == "REMOVE")),]

table(profile$FID %in% correct_ccss_exp$V1)
# FALSE  TRUE 
# 7834   143 
table(profile$IID %in% correct_ccss_exp$V1)
# 7834   143

profile$FID[profile$FID =="SJNPC018728"] <- "SJNPC018728_G1"
profile$FID[profile$FID =="SJNPC018729"] <- "SJNPC018729_G1"

profile$IID[profile$FID =="SJNPC018728_G1"] <- "SJNPC018728_G1"
profile$IID[profile$FID =="SJNPC018729_G1"] <- "SJNPC018729_G1"

# Identify matched indices in profile$IID
matched_indices <- match(profile$IID, correct_ccss_exp$V1)

# Replace FID with corresponding V2 values where matches exist
profile$FID[!is.na(matched_indices)] <- correct_ccss_exp$V2[matched_indices[!is.na(matched_indices)]]
profile$IID[!is.na(matched_indices)] <- correct_ccss_exp$V2[matched_indices[!is.na(matched_indices)]]

table(profile$IID %in% masterVCF$VCFrename)
profile[!profile$IID %in% masterVCF$VCFrename,]

replace_with <- paste0(profile$FID[grepl("SJNORM0", profile$FID)], profile$IID[grepl("SJNORM0", profile$FID)])
profile$FID[grepl("SJNORM0", profile$FID)] <- replace_with
profile$IID[grepl("SJNORM0", profile$FID)] <- replace_with
write.table(profile, "Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/POI_genetics/Ruth_et_al/prs/prs_out_wgs_survivor/POI_META_prs_updated.profile", col.names = T, row.names = F, sep = "\t", quote = F)

ccss.org <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/POI_genetics/Ruth_et_al/prs/prs_out_ccss_org/POI_META_prs.profile", header = T)
ccss.org$IID <- sub("_.+$", "", ccss.org$IID)
ccss.org$FID <- sub("_.+$", "", ccss.org$FID)
SurvivorWGS <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/POI_genetics/Ruth_et_al/prs/prs_out_wgs_survivor/POI_META_prs_updated.profile", header = T)

prs <- rbind(ccss.org, SurvivorWGS)
prs <- prs[c("FID", "IID", "SCORE")]
saveRDS(prs, file = "Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/POI_genetics/Ruth_et_al/prs/prs_wgs_survivor_and_ccss_org.rds")
