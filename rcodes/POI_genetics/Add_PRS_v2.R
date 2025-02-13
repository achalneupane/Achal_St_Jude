# Had to update the IDs after correcting the sample names in WGS
correct_ccss_exp <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/Survivor_WGS_QCed/Phenotypes/143_CCSS_rename_corrected.txt", header = F)
profile <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/POI_genetics/Ruth_et_al/prs/prs_out_wgs_survivor/POI_META_prs.profile", header = T)
# remove duplicates
profile <- profile[!((profile$IID == 3)|(profile$IID == 2)|(profile$IID == "REMOVE")),]
profile <- profile[!grepl("G2-|G1-", profile$IID),]

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

# replace_with <- paste0(profile$FID[grepl("SJNORM0", profile$FID)], profile$IID[grepl("SJNORM0", profile$FID)])
# profile$FID[grepl("SJNORM0", profile$FID)] <- replace_with
# profile$IID[grepl("SJNORM0", profile$FID)] <- replace_with
write.table(profile, "Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/POI_genetics/Ruth_et_al/prs/prs_out_wgs_survivor/POI_META_prs_updated.profile", col.names = T, row.names = F, sep = "\t", quote = F)

ccss.org <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/POI_genetics/Ruth_et_al/prs/prs_out_ccss_org/POI_META_prs.profile", header = T)
ccss.org$IID <- sub("_.+$", "", ccss.org$IID)
ccss.org$FID <- sub("_.+$", "", ccss.org$FID)
SurvivorWGS <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/POI_genetics/Ruth_et_al/prs/prs_out_wgs_survivor/POI_META_prs_updated.profile", header = T)

prs <- rbind(ccss.org, SurvivorWGS)
prs <- prs[c("FID", "IID", "SCORE")]
saveRDS(prs, file = "Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/POI_genetics/Ruth_et_al/prs/prs_wgs_survivor_and_ccss_org_v2.rds")

# summary(prs$SCORE)
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# -0.11594 -0.10530 -0.10037 -0.10123 -0.09735 -0.08736 


# Define a positive bin width
# Enhanced Histogram with Actual Frequencies
library(ggplot2)
library(scales)

ggplot(POI_meta_GRCh38_v2, aes(x = -1*V5)) + 
  geom_histogram(aes(y = after_stat(count)/sum(after_stat(count))), 
                 bins = 30, 
                 fill = "#4E79A7", 
                 color = "white") +
  geom_vline(xintercept = mean(-1*POI_meta_GRCh38_v2$V5), 
             color = "#E15759", 
             linetype = "dashed", 
             linewidth = 1) +
  scale_y_continuous(labels = percent_format(), 
                     name = "Frequency") +
  scale_x_continuous(name = "Effect Size (β)", 
                     breaks = seq(0, 1.5, 0.2)) +
  labs(title = "Size Distribution",
       subtitle = paste(nrow(POI_meta_GRCh38_v2), "SNPs | Mean β =", 
                        round(mean(-1*POI_meta_GRCh38_v2$V5), 3))) +
  theme_minimal(base_size = 14) +
  theme(panel.grid.minor = element_blank())


