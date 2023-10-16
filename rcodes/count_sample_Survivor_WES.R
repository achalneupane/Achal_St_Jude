library(haven)

## 8065 total samples in WES data
pop <- read_sas("Z:/SJShare/SJCOMMON/ECC/SJLife/SJLIFE Data Freeze/2 Final Data SJLIFE/20200430/Clinical Data/demographics.sas7bdat")
pop.survivor <- pop[pop$studypop == "Survivor",]
pop.survivor.control <- pop[grepl("Control", pop$studypop),]

all.samples <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/Survivor_WES/biallelic/sample_mapping.txt", header = F)
# > dim(all.samples)
# [1] 8055    2

all.samples[duplicated(all.samples$V2),]


sum(all.samples$V2 %in% pop$sjlid)
# 5019

CCSS <- as.data.frame(all.samples$V1[grepl("CCSS|SJNPC018728_G1|SJNPC018729_G1", all.samples$V1)])
## Not sure about these two samples in CCSS
# SJNPC018728_G1  SJNPC018728_G1
# SJNPC018729_G1  SJNPC018729_G1

dim(CCSS)
# [1] 3034    1
write.table(CCSS, "Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/Survivor_WES/biallelic/extract_CCSS.samples.txt", col.names = F, row.names = F, quote = F)

SJLIFE_survivor <- as.data.frame(all.samples$V1[all.samples$V2 %in% pop.survivor$sjlid ])
dim(SJLIFE_survivor)
# [1] 4568    1
write.table(SJLIFE_survivor, "Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/Survivor_WES/biallelic/extract_SJLIFE_survivor.txt", col.names = F, row.names = F, quote = F)

SJLIFE_control <- as.data.frame(all.samples$V1[all.samples$V2 %in% pop.survivor.control$sjlid ])
dim(SJLIFE_control)
# [1] 451   1
write.table(SJLIFE_control, "Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/Survivor_WES/biallelic/extract_SJLIFE_survivor_control.txt", col.names = F, row.names = F, quote = F)
