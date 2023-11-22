cd /research_jude/rgs01_jude/groups/sapkogrp/projects/CAB/common/WGS_Northwestern/WGS_batch_3
## download files
module load aws-cli
aws configure
aws s3 cp s3://homvrpzr-030971057479 . --recursive
aws s3 cp s3://hombzsyr-598731762349 . --recursive
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

md5$Samples <- sub("^soapnuke/clean/", "", md5$V2)
md5$Samples  <- sub("^(.*?)/.*$", "\\1", md5$Samples)
table(md5$Samples)

found.file <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects//CAB/common/WGS_Northwestern/all.files.txt")
found.file$samples <- sub(".*/", "", found.file$V1)
found.file$sm <- sub("_.*", "", found.file$samples)
found.file$sm <- sub("-.*", "", found.file$sm)

## Check which ones are not found
eqtl.all <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects//CAB/common/WGS_Northwestern/eqtl_hiPSC.list")

eqtl.all$V1[!eqtl.all$V1 %in% found.file$sm]
## unique(md5$Samples)
as.data.frame(unique(md5$Samples))
KG39_C6_P20
KG23_C12_P19
SC29_C11_P18
JR66_C3_P24
JR64_C5_P26
SC37_C6_P34
JR23_C6_P18
CE01_C5_P26
JR53_C8_P24
JR20_C1_P24
JR57_C2_P33
KG33_C2_P26
SC36_C5_P21
KG13_C13_P27
SC41_C11_P41
SC16_C4_P20
SC32_C12_P36
JR68_C8_P42
JR18_C7_P20
SC39_C11_P44
JR45_C15_P20
JR59_C1_P20
SC31_C1_P15
KG10_C9_P23
KG42_C15_P24
KG24_C15_P22
JR56_C5_P12
SC20_C12_P24
JR55_C11_P9
JR65_C3_P16
SC34_C8_P24
KG19_C6_P26
KG34_C1_P20



for line in $(cat eqtl_hiPSC.txt); do
find . -type d | grep $line
done

while read -r line; do
  find . -type d | grep  "$line" 
done < eqtl_hiPSC.txt


while read -r line; do   find . -type d | grep -v eQTL_hiPSC| grep  "$line" ; done < eqtl_hiPSC.txt >> all.files.found.txt


"KG09/KG18" "SC06"      "SC28"      "SC30"      "SC35"      "JR01/JR04" "JR06/JR61" "JR54"      "JR63"      "SC38" 

Final missing to check with Hana: SC06, "SC28", "SC30", "SC35", "JR54", "JR63", "SC38"
Response from Hana:
SC06 is SC6
JR54 missing in BGI
SC28 missing in BGI


################################
# <R>
setwd("Z:/ResearchHome/Groups/sapkogrp/projects/CAB/common/WGS_Northwestern/")

df <- read.delim("WGS_Northwestern_eQTL.txt", header = T, sep = "\t", check.names = F)
head(df)
df$ShortName <- sapply(strsplit(df$FULLPATH, "/"), function(x) tail(x, 1))
table(df$ShortName == df$`Sample name`) ## all should be true

for(i in 1:nrow(df)){
directory_path <- gsub("/research_jude/rgs01_jude/groups", "/z/ResearchHome/Groups", df$FULLPATH[i])
files <- system(paste("ls", shQuote(directory_path)), intern = TRUE)
print(paste0("********DOing i : ", i, " has total files: " , length(files)))
Sys.sleep(1)
}


## Check if all samples are available for eQTL
df2 <- read.table("eqtl_hiPSC.txt")
dim(df2)
# 100

# Remove everything after the first occurrence of "_" or "-"
df$ShortName <- sub("[-_].*", "", df$ShortName)

sum(df2$V1 %in% df$ShortName)

df2$V1[!df2$V1 %in% df$ShortName]

df2$V1[df2$V1 == "KG09/KG18"] <- "KG09"
df2$V1[df2$V1 == "SC06" ] <- "SC6"
df2$V1[df2$V1 == "JR01/JR04" ] <- "JR01"
df2$V1[df2$V1 == "JR06/JR61" ] <- "JR06"

df2$V1[!df2$V1 %in% df$ShortName]

sum(df2$V1 %in% df$ShortName)
## 100
df2$match <- df$`Sample name`[match(df2$V1, df$ShortName)]

length(df$`Sample name`)
length(unique((df$`Sample name`)))
####################################

# for file in $(ls *L01_104_1.fq.gz*); do
# echo "Doing ${file}"
# file1=${file}
# file2=$(echo ${file1} | sed 's/_1\.fq\.gz/_2.fq.gz/')
# zcat ${file1}| wc -l; zcat ${file2}| wc -l
# done

for file in $(ls *L01_124_1.fq.gz); do
echo "Doing ${file}"
file1=${file}
file2=$(echo ${file1} | sed 's/_1\.fq\.gz/_2.fq.gz/')
md5sum ${file1}
md5sum ${file2}
done