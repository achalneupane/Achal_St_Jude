library(data.table)

## Read ACMG genes
ACMG <- read.delim("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/ACMG/ACMG_genes.txt", header = T, sep = "\t", stringsAsFactors = F)
ACMG$group <- ACMG$Disease_name_and_MIM_number
ACMG$group <- sub("\\(.*", "", ACMG$group)
ACMG$group <- sub(",.*", "", ACMG$group)
ACMG$group <- sub("[0-9].*", "", ACMG$group)
ACMG$group <- trimws(ACMG$group)
ACMG$group[ACMG$group == "RPE"] <- "RPE65-related retinopathy"
ACMG$group[ACMG$group == "Juvenile polyposis/hereditary hemorrhagic telangiectasia syndrome"] <- "Juvenile polyposis syndrome"
ACMG$group <- sub(" type*", "", ACMG$group)
as.data.frame(table(ACMG$group))
ACMG$GENE <- sub("<.*", "", ACMG$Gene_via_GTR)
ACMG$inheritence <- ""

## AD or AR
ACMG$inheritence[grepl("Adenomatous polyposis coli|Aortic aneurysm|Arrhythmogenic right ventricular cardiomyopathy|Breast-ovarian cancer|Brugada syndrome|Catecholaminergic|Dilated cardiomyopathy|Ehlers|Fabry|Familial hypercholesterolemia|Familial hypertrophic cardiomyopathy|Familial medullary thyroid carcinoma|Hereditary breast cancer|hemochromatosis|Hereditary hemorrhagic telangiectasia|Hereditary paraganglioma-pheochromocytoma|Hereditary transthyretin-related amyloidosis|Hypercholesterolemia|Juvenile polyposis|Li-Fraumeni syndrome|Loeys-Dietz syndrome|Long QT|Lynch|Malignant hyperthermia|Marfan|Maturity-Onset of Diabetes|Multiple endocrine neoplasia|Myofibrillar myopathy|Neurofibromatosis|Ornithine carbamoyltransferase deficiency|Paragangliomas|Peutz-Jeghers syndrome|Pheochromocytoma|PTEN hamartoma|Retinoblastoma|Tuberous sclerosis|Von Hippel-Lindau|Wilms", ACMG$group, ignore.case = T)] <- "AD"
ACMG$inheritence[grepl("Biotinidase deficiency|MYH-associated polyposis|Pompe disease|RPE65|Wilson disease", ACMG$group, ignore.case = T)] <- "AR"

setwd("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/ACMG/")

## 1. clinvar
clinvar <- read.delim("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/Survivor_WES/annotation/snpEff_round3/all_new_clinvar_P_LP.txt", header = T, sep = "\t", stringsAsFactors = F)
dim(clinvar)
clinvar$SNP <- sub(";.*", "", clinvar$ID)
clinvar$AF <- as.numeric(clinvar$AF)
clinvar$AF_nfe <- as.numeric(clinvar$AF_nfe)
clinvar$AF_afr <- as.numeric(clinvar$AF_afr)

tt <- cbind.data.frame(clinvar$SNP, clinvar$ID, clinvar$AF_nfe, clinvar$AF_afr, clinvar$AF)

# make rare; .all is based on allee frequency from gnomAD global population; .eur is gnomAD EUR_nfe; .afr is gnomAD AFR
clinvar$new_GENE.clinvar <- gsub("-.*", "", clinvar$ANN....GENE)
clinvar <- cbind.data.frame(CHROM=clinvar$CHROM, POS=clinvar$POS, REF=clinvar$REF, ALT=clinvar$ALT, SNP=clinvar$SNP, ID=clinvar$ID, GENE=clinvar$ANN....GENE, new_GENE.clinvar=clinvar$new_GENE.clinvar, AF_nfe=clinvar$AF_nfe, AF_afr=clinvar$AF_afr, AF=clinvar$AF, Effect=clinvar$ANN....EFFECT, CLNSIG=clinvar$CLNSIG)
clinvar.save <- clinvar

ACMG$GENE[!ACMG$GENE %in% clinvar$new_GENE.clinvar]
## Keep only those in ACMG
clinvar <- clinvar[clinvar$new_GENE.clinvar %in% ACMG$GENE,]

clinvar.all <- clinvar[which(clinvar$AF <0.01),]
clinvar.eur <- clinvar[which(clinvar$AF_nfe < 0.01),]
clinvar.afr <- clinvar[which(clinvar$AF_afr < 0.01),]




loftee <- read.delim("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/Survivor_WES/annotation/snpEff_round3/loftee/loftee_HC_all_chr_with_gnomad.txt", header = T, sep = "\t", stringsAsFactors = F)
dim(loftee)
loftee$SNP <- sub(";.*", "", loftee$Uploaded_variation)
# LOF variants flagged by LOFTEE as dubious (e.g., affecting poorly conserved exons and splice variants affecting NAGNAG sites or non-canonical splice regions) will be excluded.
loftee <- loftee[!grepl("NAGNAG_SITE|NON_CAN_SPLICE", loftee$LoF_flags),]
table(loftee$LoF_flags)
# - 
# 70494 


# make rare; .all is based on allee frequency from gnomAD global population; .eur is gnomAD EUR_nfe; .afr is gnomAD AFR
loftee$new_GENE.loftee <- gsub("-.*", "", loftee$SYMBOL)

ACMG$GENE[!ACMG$GENE %in% loftee$new_GENE.loftee]
## Keep only those in ACMG
loftee <- loftee[loftee$new_GENE.loftee %in% ACMG$GENE,]


loftee$AF <- as.numeric(loftee$AF.1)
loftee$AF_nfe <- as.numeric(loftee$AF_nfe)
loftee$AF_afr <- as.numeric(loftee$AF_afr)
loftee$CHROM  <- sub("([0-9XY]+):.+", "\\1", loftee$SNP)
loftee$POS <- sub("chr[0-9XY]+:(\\d+):.+", "\\1", loftee$SNP)

# make rare
loftee.all <- loftee[which(loftee$AF <0.01),]
loftee.eur <- loftee[which(loftee$AF_nfe < 0.01),]
loftee.afr <- loftee[which(loftee$AF_afr < 0.01),]


snpeff <- read.delim("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/Survivor_WES/annotation/snpEff_round3/missense_variants_with_overlap_in_more_than_90_percent_of_prediction_tools_all_cols.txt", header = T, sep = "\t", stringsAsFactors = F)
dim(snpeff)
snpeff$SNP <- sub(";.*", "", snpeff$ID)

snpeff$AF <- as.numeric(snpeff$AF)
snpeff$AF_nfe <- as.numeric(snpeff$AF_nfe)
snpeff$AF_afr <- as.numeric(snpeff$AF_afr)

snpeff$new_GENE.snpeff <- gsub("-.*", "", snpeff$ANN....GENE)
snpeff <- cbind.data.frame(CHROM=snpeff$CHROM, POS=snpeff$POS, REF=snpeff$REF, ALT=snpeff$ALT, SNP=snpeff$SNP, ID=snpeff$ID, GENE=snpeff$ANN....GENE, new_GENE.snpeff=snpeff$new_GENE.snpeff, AF_nfe=snpeff$AF_nfe, AF_afr=snpeff$AF_afr, AF=snpeff$AF, Effect=snpeff$ANN....EFFECT, CLNSIG=snpeff$CLNSIG)
snpeff.save <- snpeff

ACMG$GENE[!ACMG$GENE %in% snpeff$new_GENE.snpeff]
## Keep only those in ACMG
snpeff <- snpeff[snpeff$new_GENE.snpeff %in% ACMG$GENE,]

# make rare
snpeff.all <- snpeff[which(snpeff$AF <0.01),]
snpeff.eur <- snpeff[which(snpeff$AF_nfe < 0.01),]
snpeff.afr <- snpeff[which(snpeff$AF_afr < 0.01),]



length(unique(c(clinvar$SNP, loftee$SNP, snpeff$SNP)))
# 930
cc <- as.data.frame(unique(c(clinvar$SNP, loftee$SNP, snpeff$SNP)))
dim(cc)
# 930
colnames(cc) <- "SNP"

## QCed Bim
bim.QC <- fread("Z:/ResearchHome/Groups/sapkogrp/projects//Genomics/common/Survivor_WES/biallelic2/plink_all/chr.ALL.Survivor_WES.GATK4180.hg38_biallelic.geno.0.1.hwe.1e-15.LCR.removed.MAC.ge.1_ID_updated.bim")
cc$SNP[!cc$SNP %in% bim.QC$V2]
cc.final <- cc[cc$SNP %in% bim.QC$V2,]

write.table(cc.final, "Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/ACMG/ACMG_rare_variants_to_extract.txt", row.names = F, col.names = F, quote = F)

## preQC
bim.preQC <- fread("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/Survivor_WES/biallelic/plink_all/chr.ALL.Survivor_WES.GATK4180.hg38_biallelic_ID_updated.bim")
cc <- as.character(cc$SNP)
cc[!cc %in% bim.preQC$V2]
table(cc %in% bim.preQC$V2)
# TRUE 
# 930

## QCed for GQ, DP and VQSR
bim.QC <- fread("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/Survivor_WES/biallelic2/plink_all/chr.ALL.Survivor_WES.GATK4180.hg38_biallelic_ID_updated.bim")
table(cc %in% bim.QC$V2)
# FALSE  TRUE 
# 179   751

## Final SJLIFE QCed data
bim.QC.sjlife <- fread("Z:/ResearchHome/Groups/sapkogrp/projects//Genomics/common/Survivor_WES/biallelic2/plink_all/chr.ALL.Survivor_WES.GATK4180.hg38_biallelic.geno.0.1.hwe.1e-15.LCR.removed.MAC.ge.1_ID_updated.bim")
table(cc %in% bim.QC.sjlife$V2)
# FALSE  TRUE 
# 202   728 


## Create carrier status
raw <- fread("Z:/ResearchHome/Groups/sapkogrp/projects//Genomics/common/ACMG/ACMG_rare_variants_ALL_recodeA.raw")
raw <- as.data.frame(raw)
rownames(raw) <- raw$IID
raw <- raw[-c(1:6)]
# HEADER <- colnames(raw)[-c(1:6)]
HEADER=colnames(raw)
HEADER = sub(pattern="_[T,A,G,C,*]+",replacement="",HEADER)
HEADER = gsub(pattern=";rs\\d+",replacement="",HEADER)
HEADER = gsub("._...DEL", "",HEADER)
HEADER = gsub("....DEL.", ".<*:DEL>",HEADER)
HEADER = gsub("\\.", ":", HEADER)
colnames(raw) <- HEADER
table(colnames(raw) %in% cc)
# FALSE  TRUE 
# 6   728
## Looks good!


## Extract clinvar European
raw.clinvar.eur <- raw[which(colnames(raw) %in% clinvar.eur$SNP)]
clinvar.eur <- clinvar.eur[clinvar.eur$SNP %in% colnames(raw)]
genes <- unique(clinvar.eur$new_GENE.clinvar)
for (i in 1:length(genes))
clinvar.eur

