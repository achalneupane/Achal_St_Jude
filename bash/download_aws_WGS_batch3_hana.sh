cd /research_jude/rgs01_jude/groups/sapkogrp/projects/CAB/common/WGS_Northwestern/WGS_batch_3/F23A480000289_HOMvrpzR
## download files
aws s3 cp s3://homvrpzr-030971057479 . --recursive

## Check md5
find . -type f -exec md5sum {} \; | awk '{print $2, $1}' >>../../../F23A480000289_HOMvrpzR/achal_md5_sum.txt


# <R>
setwd("Z:/ResearchHome/Groups/sapkogrp/projects//CAB/common/WGS_Northwestern/WGS_batch_3/F23A480000289_HOMvrpzR")

achal.md5 <- read.table("achal_md5_sum.txt")
md5 <- read.table("md5sum_check.txt")

sum(md5$V1 %in% achal.md5$V2)
# 72

## Mismatch
md5$V2[!(md5$V1 %in% achal.md5$V2)]