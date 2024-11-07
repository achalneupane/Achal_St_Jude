library(haven)
# https://wiki.stjude.org/display/CAB/Genetic+Ancestry+Estimation+by+PCA
setwd("Z:/ResearchHome/Groups/sapkogrp/projects/CAB/common/SJLIFE_CCSS_WES_101724/GermlineQC/")

load("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/PHENOTYPE/5_lifestyle_v11.RDATA")
dim(PHENO.ANY_SN)
# [1] 4401  122 # we expect 4401 samples of old WGS in new WGS and WES data

overlaps <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/gwas/pheno/SJLIFE_WGS_samples_overlap_with_CCSS_org_SNP_samples.txt", header = T)
pop <- read_sas("Z:/SJShare/SJCOMMON/ECC/SJLife/SJLIFE Data Freeze/2 Final Data SJLIFE/20200430/Clinical Data/demographics.sas7bdat")
pop.survivor <- pop[pop$studypop == "Survivor",]
community.control <- pop[grepl("Control", pop$studypop),]
table(pop.survivor$sjlid %in% community.control$sjlid)
# FALSE 
# 9366 
######################
## Clean sample IDs ##
######################
SJLIFEwes.combio <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/Survivor_WES/SJLIFE_WESsamplelist_TBIDcheck_YSapkota_02Aug2022_FINAL.txt", header = T)
# SJLIFEwes.combio <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/Survivor_WES/biallelic2/plink_all/Survivors/chr.ALL.SURVIVORS_WES.GATK4180.hg38_biallelic.geno.0.1.hwe.1e-15.LCR.removed.MAC.ge.1_updated_unique.fam", header = F)
sjlife.4507 <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/MERGED_sjlife1_2_PreQC/cleaned/MERGED.SJLIFE.1.2.GATKv3.4.VQSR.chr9.preQC_biallelic_renamed_ID_edited.vcf.gz.fam", header = F)


#############################
## WGS SJLIFE and CCSS exp ##
#############################
all.survivors.WGS <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/Survivor_WGS_QCed/QC/Survivor_WGS.GATK4180.hg38_renamed_chr1.PASS.decomposed.qced.fam", header = F)
dim(all.survivors.WGS)
# [1] 7977    6

## 1. Check SJLIFE
ss1 <- sjlife.4507$V2[!sjlife.4507$V2 %in% all.survivors.WGS$V2]
ss1
# [1] "SJL5207307" "SJL5188202" "SJL5128013" "SJL5094513" "SJL5052915" "SJL5057815"

table(PHENO.ANY_SN$sjlid %in% all.survivors.WGS$V2)
# FALSE  TRUE 
# 6  4395
ss2 <- PHENO.ANY_SN$sjlid[!PHENO.ANY_SN$sjlid %in% all.survivors.WGS$V2]
ss2
# [1] "SJL5052915" "SJL5057815" "SJL5094513" "SJL5188202" "SJL5128013" "SJL5207307"
table(ss1%in%ss2)
# TRUE 
# 6

table(unique(all.survivors.WGS$V2) %in% pop.survivor$sjlid)
# FALSE  TRUE 
# 3320  4654
new.wgs.sjlife <- unique(all.survivors.WGS$V2)[unique(all.survivors.WGS$V2) %in% pop.survivor$sjlid]


## 2. Check AA samples
AA.90samples <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/RNAseq/common/RNAseq_sjlife/jenn_RNAseq/AA_samples.txt", header = F)
table(AA.90samples$V1 %in% all.survivors.WGS$V2)
# TRUE 
# 90 

## 3. Check Community controls
table(all.survivors.WGS$V2 %in% community.control$sjlid)
# FALSE  TRUE 
# 7527   450
wgs.Controls <- all.survivors.WGS$V2[all.survivors.WGS$V2 %in% community.control$sjlid]

## 4. Check CCSS_exp (This is the QCed CCSS_exp cohort, we would expect 2839 of these samples in both all.survivors.WGS and old.wes.samples datasets)
ccss_exp <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/ccss_exp_wgs/QCed_samples", header = F)
dim(ccss_exp)
# [1] 2839    1

table(ccss_exp$V1 %in% all.survivors.WGS$V2)
# FALSE  TRUE 
# 143  2696 
ccss.not.in.new.wgs <- ccss_exp$V1[!ccss_exp$V1 %in% all.survivors.WGS$V2]
ccss.not.in.new.wgs
# [1] 15479822 15417752 15478832 15251382 15476582 15476132 15481642 15478662 15283612 15482242 15250832 15253472 15476722 15480502 15251792 15419262 15283172 15473712 15483012 15417832 15485992
# [22] 15419092 15254562 15484582 15252622 15482792 15419492 15479242 15479772 15484882 15418202 15475772 15250628 15479965 15486305 15253425 15252735 15478945 15417675 15475785 15478325 15479785
# [43] 15284461 15477231 15479401 15482931 15253231 15473311 15477291 15254851 15472721 15283871 15474431 15417351 15478821 15480371 15250031 15472981 15253311 15477891 15482341 15475651 15479021
# [64] 15284211 15253641 15476541 15476561 15253971 15253291 15477371 15482821 15484781 15250941 15250411 15417511 15482301 15473601 15251551 15474541 15418751 15479131 15474711 15252431 15483701
# [85] 15254021 15283701 15486561 15479861 15252651 15476471 15284291 15250811 15474281 15283921 15481761 15484931 15250561 15482203 15419703 15478913 15251583 15251563 15283743 15474933 15479203
# [106] 15417493 15251663 15253753 15252753 15252893 15419593 15472713 15253683 15284093 15418073 15475343 15250083 15476263 15417613 15483333 15252926 15485966 15478684 15480654 15283244 15480284
# [127] 15486494 15478144 15475584 15479694 15418964 15476774 15252534 15480934 15475434 15283414 15251894 15250764 15283027 15485077 15251537 15419417 15283287
ccss.in.new.wgs <- ccss_exp$V1[ccss_exp$V1 %in% all.survivors.WGS$V2]

sampleID <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/Survivor_WGS_QCed/Phenotypes/wgs_sampleIDs_cleaned_AN_recheck.txt", header = T, sep = "\t")
sjlife_ccss_exp_overlaps <- sampleID[(sampleID$SJLID!="" & !is.na(sampleID$CCSSID)),]

table(ccss.not.in.new.wgs %in% sjlife_ccss_exp_overlaps$CCSSID)

table(sjlife_ccss_exp_overlaps$SJLID %in% PHENO.ANY_SN$sjlid)
write.table(sjlife_ccss_exp_overlaps, file="Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/Survivor_WES_QC/sjlife_ccss_exp_overlaps.txt", row.names=FALSE, col.names=TRUE, quote=FALSE, sep = "\t")

#########################
## Old WES VCF samples ##
#########################
#### preQC 
old.wes.samples = read.table("/ResearchHome/Groups/sapkogrp/projects/Genomics/common/Survivor_WES/biallelic2/plink_all/chr.ALL.Survivor_WES.GATK4180.hg38_biallelic.geno.0.1.hwe.1e-15.LCR.removed.MAC.ge.1_ID_updated_unique.fam", header=F)
dim(old.wes.samples) ## These are unique samples
# 8030    6
table(sjlife.4507$V2 %in% old.wes.samples$V2)
# FALSE  TRUE 
# 21  4486 

## 1. Check SJLIFE
ss3 <- sjlife.4507$V2[!sjlife.4507$V2 %in% old.wes.samples$V2]
ss3
# [1] "SJL1437501" "SJL1285201" "SJL1224901" "SJL1245501" "SJL1246701" "SJL5037906" "SJL1281212" "SJL5050105" "SJL1679008" "SJL5068708" "SJL1648308" "SJL5052915" "SJL5057815" "SJL5024018"
# [15] "SJL4749916" "SJL1750516" "SJL5058217" "SJL5271801" "SJL5385207" "SJL2521413" "SJL5321816"

table(PHENO.ANY_SN$sjlid %in% old.wes.samples$V2)
# FALSE  TRUE 
# 21  4380 
ss4 <- PHENO.ANY_SN$sjlid[!PHENO.ANY_SN$sjlid %in% old.wes.samples$V2]
ss4
# [1] "SJL1245501" "SJL1246701" "SJL1750516" "SJL5024018" "SJL1648308" "SJL1285201" "SJL1224901" "SJL1679008" "SJL5068708" "SJL1437501" "SJL1281212" "SJL5037906" "SJL5050105" "SJL5052915"
# [15] "SJL4749916" "SJL5057815" "SJL5058217" "SJL2521413" "SJL5271801" "SJL5321816" "SJL5385207"
table(ss3%in%ss4)
# TRUE 
# 21

# missing in both WGS and old WES
ss2[!ss2 %in% ss4]
# "SJL5094513" "SJL5188202" "SJL5128013" "SJL5207307"
ss2[ss2 %in% ss4]
# [1] "SJL5052915" "SJL5057815"

table(unique(old.wes.samples$V2) %in% pop.survivor$sjlid)
# FALSE  TRUE 
# 3497  4533 
old.wes.sjlife <- unique(old.wes.samples$V2)[unique(old.wes.samples$V2) %in% pop.survivor$sjlid]
table(old.wes.sjlife %in% new.wgs.sjlife)
# FALSE  TRUE 
# 30  4503 
table(new.wgs.sjlife %in% old.wes.sjlife)
# FALSE  TRUE 
# 151  4503 
## removing 90 AA samples from 4503 
# 4503-90 = 4413
new.wgs.sjlife.not.AA <- new.wgs.sjlife[!new.wgs.sjlife %in% AA.90samples$V1]
old.wes.sjlife.not.AA <- old.wes.sjlife[!old.wes.sjlife %in% AA.90samples$V1]
table(new.wgs.sjlife.not.AA %in% old.wes.sjlife.not.AA)
# FALSE  TRUE 
# 151  4413 
sjlife.overlaps.wgs.old.wes <- new.wgs.sjlife.not.AA[new.wgs.sjlife.not.AA %in% old.wes.sjlife.not.AA]
table(sjlife.overlaps.wgs.old.wes %in% PHENO.ANY_SN$sjlid)
# FALSE  TRUE 
# 37  4376 

table(sjlife_ccss_exp_overlaps$SJLID %in% sjlife.overlaps.wgs.old.wes)
# FALSE  TRUE 
# 132    14 

## 2. Check AA samples
AA.90samples <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/RNAseq/common/RNAseq_sjlife/jenn_RNAseq/AA_samples.txt", header = F)
table(AA.90samples$V1 %in% old.wes.samples$V2)
# TRUE 
# 90 
table(AA.90samples$V1 %in% sjlife.4507$V2)
# FALSE 
# 90 

## 3. Check Community controls
old.wes.samples.Controls = read.table("/ResearchHome/Groups/sapkogrp/projects/Genomics/common/Survivor_WES/biallelic2/plink_all/Controls/chr.ALL.survivor.control_WES.GATK4180.hg38_biallelic.geno.0.1.hwe.1e-15.LCR.removed.MAC.ge.1.fam", header=F)
table(old.wes.samples.Controls$V2 %in% community.control$sjlid)
# 451
table(old.wes.samples$V2 %in% community.control$sjlid)
# FALSE  TRUE 
# 7579   451 

table(sjlife.4507$V2 %in% community.control$sjlid)
# FALSE  TRUE 
# 4401   106 
table(PHENO.ANY_SN$sjlid %in% community.control$sjlid)
# FALSE 
# 4401 
table(AA.90samples$V1 %in% community.control$sjlid)
# FALSE 
# 90 

# In WGS vs old WES data
table(old.wes.samples.Controls$V2 %in% wgs.Controls)
# FALSE  TRUE 
# 2   449 
table(wgs.Controls %in% old.wes.samples.Controls$V2)
# FALSE  TRUE 
# 1   449
wgs.Controls[!wgs.Controls %in% old.wes.samples.Controls$V2]
# [1] SJL5144399

old.wes.samples.Controls$V2[!old.wes.samples.Controls$V2 %in% wgs.Controls]
# [1] "SJL5115599" "SJL5139899"


## 4. Check CCSS_exp
ccss_exp <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/ccss_exp_wgs/QCed_samples", header = F)
table(ccss_exp$V1 %in% old.wes.samples$V2)
# FALSE  TRUE 
# 2  2837 
ccss_exp$V1[!ccss_exp$V1 %in% old.wes.samples$V2]
# [1] 3267622 2512001

table(old.wes.samples$V2 %in% overlaps$V2)
# FALSE 
# 8030 

ccss_exp_missing_in_wgs <- ccss_exp$V1[!ccss_exp$V1 %in% all.survivors.WGS$V2]
ccss_exp_missing_in_wes <- ccss_exp$V1[!ccss_exp$V1 %in% old.wes.samples$V2]
table(ccss_exp_missing_in_wes %in% ccss_exp_missing_in_wgs)
# FALSE 
# 2 

ccss.in.old.wes <- ccss_exp$V1[ccss_exp$V1 %in% old.wes.samples$V2]
# overlap of ccss_exp in new WGS and old WES
table(ccss.in.old.wes %in% ccss.in.new.wgs)
# FALSE  TRUE 
# 143  2694 

table(sjlife_ccss_exp_overlaps$CCSSID %in% ccss.in.old.wes)
# TRUE 
# 146 


# Combine sample IDs with their sources
all.wes.samples <- data.frame(
  SampleID = c(PHENO.ANY_SN$sjlid, old.wes.samples.Controls$V2,  AA.90samples$V1, ccss_exp$V1),
  Source = c(
    rep("sjlife", length(PHENO.ANY_SN$sjlid)),
    rep("community.controls", length(old.wes.samples.Controls$V2)),
    rep("Survivor_African", length(AA.90samples$V1)),
    rep("ccss_exp", length(ccss_exp$V1))
  )
)

# Write to a tab-separated file
write.table(all.wes.samples, file = "Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/Survivor_WES/filtered_wes_samples_by_cohort.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)


#### postQC 
old.wes.samples.Survivor = read.table("/ResearchHome/Groups/sapkogrp/projects/Genomics/common/Survivor_WES/biallelic2/plink_all/Survivors/chr.ALL.SURVIVORS_WES.GATK4180.hg38_biallelic.geno.0.1.hwe.1e-15.LCR.removed.MAC.ge.1.fam", header=F)
table(sjlife.4507$V2 %in% old.wes.samples.Survivor$V2)
# FALSE  TRUE 
# 127  4380 

## 1. Check SJLIFE
ss5 <- sjlife.4507$V2[!sjlife.4507$V2 %in% old.wes.samples.Survivor$V2]

table(PHENO.ANY_SN$sjlid %in% old.wes.samples.Survivor$V2)
# FALSE  TRUE 
# 21  4380 
ss6 <- PHENO.ANY_SN$sjlid[!PHENO.ANY_SN$sjlid %in% old.wes.samples.Survivor$V2]
ss6
# [1] "SJL1245501" "SJL1246701" "SJL1750516" "SJL5024018" "SJL1648308" "SJL1285201" "SJL1224901" "SJL1679008" "SJL5068708" "SJL1437501" "SJL1281212" "SJL5037906" "SJL5050105" "SJL5052915"
# [15] "SJL4749916" "SJL5057815" "SJL5058217" "SJL2521413" "SJL5271801" "SJL5321816" "SJL5385207"
table(ss5%in%ss6)
# FALSE  TRUE 
# 106    21 

## 2. Check AA samples
AA.90samples <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/RNAseq/common/RNAseq_sjlife/jenn_RNAseq/AA_samples.txt", header = F)
table(AA.90samples$V1 %in% old.wes.samples.Survivor$V2)
# TRUE 
# 90 

## 3. Check Community controls
old.wes.samples.Controls = read.table("/ResearchHome/Groups/sapkogrp/projects/Genomics/common/Survivor_WES/biallelic2/plink_all/Controls/chr.ALL.survivor.control_WES.GATK4180.hg38_biallelic.geno.0.1.hwe.1e-15.LCR.removed.MAC.ge.1.fam", header=F)
table(old.wes.samples.Controls$V2 %in% community.control$sjlid)
# TRUE 
# 451 

table(old.wes.samples.Controls$V2 %in% community.control$sjlid)
# TRUE 
# 451
table(old.wes.samples.Controls$V2 %in% wgs.Controls)
# FALSE  TRUE 
# 2   449 
old.wes.samples.Controls$V2[!old.wes.samples.Controls$V2 %in% wgs.Controls]
# "SJL5115599" "SJL5139899"

## 4. Check CCSS_exp
old.wes.samples.CCSS_exp = read.table("/ResearchHome/Groups/sapkogrp/projects/Genomics/common/Survivor_WES/biallelic2/plink_all/CCSS_exp/chr.ALL.CCSS_exp_WES.GATK4180.hg38_biallelic.geno.0.1.hwe.1e-15.LCR.removed.MAC.ge.1.fam", header=F)
dim(old.wes.samples.CCSS_exp)
# 3036   6
table(ccss_exp$V1 %in% old.wes.samples.CCSS_exp$V2)
# FALSE  TRUE 
# 2  2837 
ccss_exp$V1[!ccss_exp$V1 %in% old.wes.samples$V2]
# [1] 3267622 2512001; These are probably SJNPC018728_G1 and SJNPC018729_G1?? Need to check with Kyla











##########################################################
## New WES VCF samples in SJLIFE, CCSS_exp and CCSS_org ##
##########################################################
new.wes.VCF.samples <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/Survivor_WES_QC/biallelic2/VCFsample_names.txt", header = F)
new.wes.VCF.samples$SJLID <- SJLIFEwes.combio$SJLID[match(new.wes.VCF.samples$V1, SJLIFEwes.combio$CompBioID )]

table(unique(new.wes.VCF.samples$SJLID) %in% sjlife.4507$V2)
# FALSE  TRUE 
# 152  4405 
table(unique(new.wes.VCF.samples$SJLID) %in% PHENO.ANY_SN$sjlid)
# FALSE  TRUE 
# 258  4299
table(PHENO.ANY_SN$sjlid %in% unique(new.wes.VCF.samples$SJLID))
# FALSE  TRUE 
# 102  4299 

AA.90samples <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/RNAseq/common/RNAseq_sjlife/jenn_RNAseq/AA_samples.txt", header = F)
table(AA.90samples$V1 %in% SJLIFEwes.combio$SJLID)
# TRUE 
# 90
table(AA.90samples$V1 %in% new.wes.VCF.samples$SJLID)
# TRUE 
# 90

table(unique(new.wes.VCF.samples$SJLID) %in% PHENO.ANY_SN$sjlid)
# FALSE  TRUE 
# 258  4299 
table(unique(SJLIFEwes.combio$SJLID) %in% PHENO.ANY_SN$sjlid)
# FALSE  TRUE 
# 604  4380 

table(SJLIFEwes.combio$CompBioID %in% new.wes.VCF.samples$V1)
# 4589
table(SJLIFEwes.combio$SJLID %in% unique(new.wes.VCF.samples$V2))

unique(SJLIFEwes.combio$V3)[!unique(SJLIFEwes.combio$V3) %in% PHENO.ANY_SN$sjlid]

PHENO.ANY_SN$sjlid[!PHENO.ANY_SN$sjlid %in% unique(SJLIFEwes.combio$V3)]
# [1] "SJL1245501" "SJL1246701" "SJL1750516" "SJL5024018" "SJL1648308" "SJL1285201" "SJL1224901" "SJL1679008" "SJL5068708" "SJL1437501" "SJL1281212" "SJL5037906" "SJL5050105" "SJL5052915"
# [15] "SJL4749916" "SJL5057815" "SJL5058217" "SJL2521413" "SJL5271801" "SJL5321816" "SJL5385207"

# all.samples[duplicated(all.samples$V2),]

# # Keep only those from Germline QC
# germline.QC <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/CAB/common/SJLIFE_CCSS_WES_101724/GermlineQC//Survivor_WES.fam")
# dim(germline.QC)
# # 8027
# sum(all.samples$V1 %in% germline.QC$V1)
# all.samples <- all.samples[all.samples$V1 %in% germline.QC$V1,]
dim(all.samples)
# 8055
# remove sex problem and het filtered
# all.samples.original <- all.samples
# all.samples <- all.samples[!all.samples$V1 %in% sex.problem.mind.and.3sd.outliers,]
# dim(all.samples)
sum(!all.samples$V1 %in% sex.problem.mind.and.3sd.outliers)
# 7907 # after sample level QC

## 8065 total samples in WES data
pop <- read_sas("Z:/SJShare/SJCOMMON/ECC/SJLife/SJLIFE Data Freeze/2 Final Data SJLIFE/20200430/Clinical Data/demographics.sas7bdat")
pop.survivor <- pop[pop$studypop == "Survivor",]
pop.survivor.control <- pop[grepl("Control", pop$studypop),]











#####################
## How many in WGS ##
#####################
WGS_SJLID <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/MERGED_sjlife1_2_PreQC/cleaned/MERGED.SJLIFE.1.2.GATKv3.4.VQSR.chr9.preQC_biallelic_renamed_ID_edited.vcf.gz.fam")
## Before QC
SJLID <- all.samples.before.QC$V2[all.samples.before.QC$cohort == "Survivor"]
sum(WGS_SJLID$V2 %in% SJLID)
## 4380 ## Survivors from WGS in Survivor WES

SJLID <- all.samples.before.QC$V2[all.samples.before.QC$cohort == "Community_control"]
sum(WGS_SJLID$V2 %in% SJLID)
# 106 ## Survivors from WGS in Community_control WES
## total 4486 in WGS if we include both Survivor and community control in SJLIFE before QC

## AFter QC
SJLID <- all.samples$V2[all.samples$cohort == "Survivor"]
sum(WGS_SJLID$V2 %in% SJLID)
## 4368 ## Survivors from WGS in Survivor WES

SJLID <- all.samples$V2[all.samples$cohort == "Community_control"]
sum(WGS_SJLID$V2 %in% SJLID)
# 106 ## Survivors from WGS in Community_control WES

## total 4474 in WGS if we include both Survivor and community control in SJLIFE

## Rename sample ID
setwd("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/Survivor_WES/biallelic/plink_all")
all.samples <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/Survivor_WES/biallelic/plink_all/chr19.Survivor_WES.GATK4180.hg38_biallelic.fam", header = F)

all.samples <- all.samples[1:2]

SJLIFE <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/Survivor_WES/SJLIFE_WESsamplelist_TBIDcheck_YSapkota_02Aug2022_FINAL.txt", header = T)
sum(duplicated(SJLIFE$SJLID))

all.samples$V3 <- all.samples$V1
all.samples$V4 <- SJLIFE$SJLID[match(all.samples$V2, SJLIFE$CompBioID)]
all.samples$V4[is.na(all.samples$V4)] <- all.samples$V2[is.na(all.samples$V4)]

all.samples$V4 <- sub(".*CCSS-0", "", all.samples$V4)
all.samples$V4 <- sub(".*CCSS-", "", all.samples$V4)

write.table(all.samples, "Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/Survivor_WES/rename_samples.txt", col.names = F, row.names = F, quote = F)


## There are duplicate SJLIFE IDs, we can romove the ones with low call rate here:
# Check for duplicated IID values
imiss = read.table("/ResearchHome/Groups/sapkogrp/projects/Genomics/common/Survivor_WES/biallelic2/plink_all/Survivors/chr.ALL.SURVIVORS_WES.GATK4180.hg38_biallelic.geno.0.1.hwe.1e-15.LCR.removed.MAC.ge.1_updated_missing.imiss", header=TRUE)
duplicate_samples <- imiss[duplicated(imiss$IID) | duplicated(imiss$IID, fromLast = TRUE), ]
dim(duplicate_samples)
# 69

# Sort and keep only the sample with the lowest F_MISS for each IID
library(dplyr)
samples_to_keep <- duplicate_samples %>%
  group_by(IID) %>%
  slice_min(F_MISS, with_ties = FALSE)  # with_ties = FALSE ensures only one sample is kept per group

# Identify the samples to remove by selecting all not in 'samples_to_keep'
samples_to_remove <- anti_join(duplicate_samples, samples_to_keep, by = c("FID", "IID"))

# View the samples to be removed
# View(samples_to_remove)
dim(samples_to_remove)
## 35
write.table(samples_to_remove[, c("FID", "IID")], "Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/Survivor_WES/biallelic2/plink_all/Survivors/duplicate_samples_to_remove.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)

################################################################################################################
################################################################################################################
################################################################################################################
################################################################################################################
################################################################################################################
################################################################################################################
################################################################################################################
## heterozygosity
# # Load the data
# het_data <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/CAB/common/SJLIFE_CCSS_WES_101724/GermlineQC/heterozygosity.het", header=TRUE)
# het_data$F <- as.numeric(het_data$F) 
# het_data$IID <- factor(het_data$IID) 
# 
# # Calculate the mean and standard deviation of observed heterozygosity
# mean_het <- mean(het_data$F)
# sd_het <- sd(het_data$F)
# # Calculate the threshold
# threshold <- mean_het + 3 * sd_het
# # Identify outliers
# outliers <- het_data$IID[het_data$F > threshold]
# outliers
# # "SJRB056831_G1-TB-16-10618"
# 
# mean_het <- mean(het_data$F)
# sd_het <- sd(het_data$F)
# 
# 
# 
# # Load required libraries (if not already installed)
# library(ggplot2)
# 
# # Define the thresholds
# thresholds <- c(-8, -7, -6, -5, -4, -3, 3, 4, 5, 6, 7, 8)
# threshold_values <- mean_het + thresholds * sd_het
# 
# plot <- ggplot(het_data, aes(x = IID, y = F)) +
#   geom_point(aes(color = F < threshold_values[6] | F > threshold_values[7]), size = 3) +
#   geom_hline(yintercept = mean_het, linetype = "dashed", color = "white") +  # Add mean line
#   scale_color_manual(values = c("TRUE" = "red", "FALSE" = "blue")) +
#   geom_hline(yintercept = threshold_values, linetype = "dashed", color = "red") +
#   labs(title = "Scatter Plot of Samples by Cutoff Threshold",
#        x = "Samples",
#        y = "F Value",
#        color = "Outlier") +
#   
#   # Adjust plot appearance
#   theme_minimal()
# 
# # Label the threshold lines
# labels <- data.frame(
#   threshold = thresholds,
#   y = threshold_values,
#   x = rep(1, 2)  # Set x to 1 for both labels
# )
# 
# plot + geom_text(data = labels, aes(x = x, y = y, label = threshold), vjust = -0.5, hjust = 0, size = 4)
# 
# 
# write.table(outliers, file="Z:/ResearchHome/Groups/sapkogrp/projects/CAB/common/SJLIFE_CCSS_WES_101724/GermlineQC//heterozygosity_outside_8_std_outliers.txt", row.names=FALSE, col.names=FALSE, quote=FALSE)




