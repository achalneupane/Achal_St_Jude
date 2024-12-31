## Check md5
MD5.downloaded <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/CAB/common/eQTL_RNAseq_NW_12_6_2024/usftp21.novogene.com/01.RawData/all_md5_downloaded.txt")
MD5.Achal <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/CAB/common/eQTL_RNAseq_NW_12_6_2024/usftp21.novogene.com/01.RawData/md5sums_Achal.txt")

table(MD5.Achal$V1 %in% MD5.downloaded$V1)


setwd("Z:/ResearchHome/Groups/sapkogrp/projects/CAB/common/eQTL_RNAseq_NW_12_6_2024/usftp21.novogene.com/01.RawData")
## Create sample info file 
DIR = "/research_jude/rgs01_jude/groups/sapkogrp/projects/CAB/common/eQTL_RNAseq_NW_12_6_2024/usftp21.novogene.com/01.RawData/"
all.dirs <- gsub("./","", list.dirs()[grepl("/", list.dirs())])

all.dirs <- as.data.frame(all.dirs)
colnames(all.dirs) <- "SubjectID"

duplicated_dirs <- all.dirs[rep(1:nrow(all.dirs), each = 2), ]
duplicated_dirs <- as.data.frame(duplicated_dirs)

colnames(duplicated_dirs) <- "SubjectID"
duplicated_dirs$Lab_sample_Ids <- duplicated_dirs$SubjectID

duplicated_dirs$SubjectID <- paste0(DIR, duplicated_dirs$SubjectID, "/", duplicated_dirs$SubjectID)

## R1
duplicated_dirs$SubjectID[seq(1, length(duplicated_dirs$SubjectID), by = 2)] <- paste0(duplicated_dirs$SubjectID[seq(1, length(duplicated_dirs$SubjectID), by = 2)], "_1.fq.gz")
## R2
duplicated_dirs$SubjectID[seq(2, length(duplicated_dirs$SubjectID), by = 2)] <- paste0(duplicated_dirs$SubjectID[seq(2, length(duplicated_dirs$SubjectID), by = 2)], "_2.fq.gz")

write.table(duplicated_dirs, "CAB_Transcriptomics_SampleInfo_2024.txt", col.names = T, row.names = F, sep = "\t", quote = F)
