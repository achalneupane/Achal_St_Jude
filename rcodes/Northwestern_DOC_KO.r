setwd("Z:/ResearchHome/Groups/sapkogrp/projects/CAB/common/WGS_193/variantCalling/WGS_193/annotation")

## Part 1
## Missing rsIDs in ucsc browser:
# rs104642
# rs8187694

############
## Part 2 ##
#################################
## Read annotation of interest ##
#################################
# read original spreadsheet
dat <- read.table("Dox_KO_Supplementary_Table_1_v31_text.txt", sep = "\t", header = T)
dat$SNP.ID <- gsub(" ", "", dat$SNP.ID)
## read new bim file 
bim.extracted <- read.table("var_of_interest_plink.bim")  
bim.extracted$V1 <- paste0("chr", bim.extracted$V1)
bim.extracted$KEY <- paste(bim.extracted$V1, bim.extracted$V4, sep = ":")
bim.extracted$V2 [ bim.extracted$KEY == "chr7:99672916"] <- "rs776746"

sum(bim.extracted$KEY %in% bim.extracted$KEY)
# 144
## SNP eff
snpeff <- read.table("var_of_interest.vcf-annot-snpeff-dbnsfp-ExAC.0.3-clinvar.GRCh38.vcf.dbSNP155-FIELDS-simple.txt", sep ="\t", header = T)
# snpeff$ID <- sub(";.*", "", snpeff$ID)
snpeff$KEY <- paste(snpeff$CHROM, snpeff$POS, sep = ":")

sum(unique(dat$SNP.ID) %in% bim.extracted$V2)
unique(dat$SNP.ID)[!unique(dat$SNP.ID) %in% bim.extracted$V2]
length(c("rs2229109","rs104642","rs4148350","rs45511401","rs8053266","rs4148399","rs45451097","rs34381428","rs8187694","rs7910642","rs2231142","rs4340",
                "rs7801891","rs2229540","rs10836235","rs13397992","rs2069522","rs2069526","rs4646427","rs2294950","rs2020870","rs1736557",
                "rs1800562","rs7196087","rs2955159","rs34223702","rs78000710","rs564101364","rs148559410","rs7542939","rs11066861","rs723685","rs2305364",
                "rs17863783","rs4407290"))

# 35

snpeff$KEY1 <- paste(snpeff$CHROM, snpeff$POS, snpeff$REF, snpeff$ALT, sep = ":")
# cc <- snpeff[3:7]

snpeff$ID <- bim.extracted$V2[match(snpeff$KEY, bim.extracted$KEY)]

sum(unique(dat$SNP.ID) %in% snpeff$ID)
unique(dat$SNP.ID)[!unique(snpeff$ID) %in% unique(dat$SNP.ID)]
# rename chr7:99672916 to "rs776746" in bim file
snpeff$ID [ snpeff$KEY == "chr7:99672916"] <- "rs776746"

df <- snpeff[match(dat$SNP.ID, snpeff$ID),]


## Add gnomAD from annovar
annovar <- read.table("ANNOVAR_var_of_interest.vcf-annot-snpeff-dbnsfp-ExAC.0.3-clinvar.GRCh38.vcf.dbSNP155.vcf.hg38_multianno.txt", header = T, sep = "\t")
annovar$KEY <- paste0(annovar$Chr, ":", annovar$End)

sum(bim.extracted$KEY  %in% annovar$KEY)
# 142

df$MAF.gnomAD.genome.ALL <- annovar$gnomAD_genome_ALL[match(df$KEY, annovar$KEY)]
df$MAF.gnomAD.genome.NFE <- annovar$gnomAD_genome_NFE[match(df$KEY, annovar$KEY)]
df$MAF.1000genome <- annovar$X1000g2015aug_all[match(df$KEY, annovar$KEY)]

colnames(df)

colnames(df) <- c("CHROM", "POS", "ID", "REF", "ALT1", "ALT", "EFFECT", "IMPACT", "GENE", "GENEID", "FEATURE", "FEATUREID", "HGVS_C", "HGVS_P", "dbNSFP_CADD_phred", "dbNSFP_1000Gp3_AF", "dbNSFP_ExAC_AF", "dbNSFP_ExAC_Adj_AF",
"dbNSFP_MetaSVM_score", "dbNSFP_MetaSVM_rankscore", "dbNSFP_MetaSVM_pred", "dbNSFP_clinvar_clnsig", "CLNSIG", "KEY", "KEY1", "MAF.gnomAD.genome.ALL", "MAF.gnomAD.genome.NFE", "MAF.1000genome")

df <- df[, c("CHROM", "POS", "ID", "REF", "ALT", "EFFECT", "GENE", "dbNSFP_CADD_phred", "CLNSIG", "MAF.gnomAD.genome.ALL", "MAF.gnomAD.genome.NFE", "MAF.1000genome")]
dim(df)

## Part 3
# read recoded plink
raw <- read.delim("var_of_interest_plink_recodeA.raw", sep = " ", header = T, check.names = F)
colnames(raw)[grepl("^._T", colnames(raw))] <- "rs776746"
dim(raw)
head(raw)

HEADER = sub(pattern="_[T,A,G,C,*]+",replacement="",colnames(raw))
colnames(raw) = HEADER
raw <- raw[!grepl("FID|IID|PAT|MAT|SEX|PHENOTYPE", colnames(raw))]
colnames(raw) # check this with the .bim REF and NON-REFERENCE alleles
raw.t <- as.data.frame(t(raw))
raw.t$V1
raw.t$rsID <- colnames(raw)


df$WGS_193_genotype <- raw.t$V1[match(df$ID, raw.t$rsID)]
df$WGS_193_genotype[is.na(df$WGS_193_genotype)] <- ""

df$WGS_193_genotype[df$WGS_193_genotype == 0] <- paste0(df$REF[df$WGS_193_genotype == 0], "/", df$REF[df$WGS_193_genotype == 0])
df$WGS_193_genotype[df$WGS_193_genotype == 1] <- paste0(df$REF[df$WGS_193_genotype == 1], "/", df$ALT[df$WGS_193_genotype == 1])
df$WGS_193_genotype[df$WGS_193_genotype == 2] <- paste0(df$ALT[df$WGS_193_genotype == 2], "/",df$ALT[df$WGS_193_genotype == 2])

df$CADD <- sub(",.*", "", df$dbNSFP_CADD_phred)  
df <- df[!is.na(df$ID),]
df$EFFECT <- gsub("_", " ", df$EFFECT)

df$CLINVAR <- df$CLNSIG
df$CLINVAR[grepl("^Pathogenic", df$CLINVAR)] <- "Pathogenic/Likely-pathogenic"
df <- df[!duplicated(df$ID),]
dim(df)

write.table(df, "Final_list_of_variants_in_WGS_193.txt", row.names = F, col.names = T, quote = F, sep = "\t")


#############################################
## check missing 36 variants in Pre QC VCF ##
#############################################
setwd("Z:/ResearchHome/Groups/sapkogrp/projects/CAB/common/WGS_193/variantCalling/WGS_193/Merge_Intervals_ToChromVCF/annotation/")
## Part 1
## Missing rsIDs in ucsc browser:
# rs104642
# rs8187694

############
## Part 2 ##
#################################
## Read annotation of interest ##
#################################
# read original spreadsheet
dat <- read.table("Dox_KO_Supplementary_Table_1_v31_text.txt", sep = "\t", header = T)
dat$SNP.ID <- gsub(" ", "", dat$SNP.ID)
## read new bim file 
bim.extracted <- read.table("var_of_interest_plink.bim")  
bim.extracted$V1 <- paste0("chr", bim.extracted$V1)
bim.extracted$KEY <- paste(bim.extracted$V1, bim.extracted$V4, sep = ":")
bim.extracted$V2 [ bim.extracted$KEY == "chr7:99672916"] <- "rs776746"

sum(bim.extracted$KEY %in% bim.extracted$KEY)
# 144
## SNP eff
snpeff <- read.table("var_of_interest.vcf-annot-snpeff-dbnsfp-ExAC.0.3-clinvar.GRCh38.vcf.dbSNP155-FIELDS-simple.txt", sep ="\t", header = T)
# snpeff$ID <- sub(";.*", "", snpeff$ID)
snpeff$KEY <- paste(snpeff$CHROM, snpeff$POS, sep = ":")

sum(unique(dat$SNP.ID) %in% bim.extracted$V2)
unique(dat$SNP.ID)[!unique(dat$SNP.ID) %in% bim.extracted$V2]
length(c("rs2229109","rs104642","rs4148350","rs45511401","rs8053266","rs4148399","rs45451097","rs34381428","rs8187694","rs7910642","rs2231142","rs4340",
         "rs7801891","rs2229540","rs10836235","rs13397992","rs2069522","rs2069526","rs4646427","rs2294950","rs2020870","rs1736557",
         "rs1800562","rs7196087","rs2955159","rs34223702","rs78000710","rs564101364","rs148559410","rs7542939","rs11066861","rs723685","rs2305364",
         "rs17863783","rs4407290"))

# 35

snpeff$KEY1 <- paste(snpeff$CHROM, snpeff$POS, snpeff$REF, snpeff$ALT, sep = ":")
# cc <- snpeff[3:7]

snpeff$ID <- bim.extracted$V2[match(snpeff$KEY, bim.extracted$KEY)]

sum(unique(dat$SNP.ID) %in% snpeff$ID)
unique(dat$SNP.ID)[!unique(snpeff$ID) %in% unique(dat$SNP.ID)]
# rename chr7:99672916 to "rs776746" in bim file
snpeff$ID [ snpeff$KEY == "chr7:99672916"] <- "rs776746"

df <- snpeff[match(dat$SNP.ID, snpeff$ID),]


## Add gnomAD from annovar
annovar <- read.table("ANNOVAR_var_of_interest.vcf-annot-snpeff-dbnsfp-ExAC.0.3-clinvar.GRCh38.vcf.dbSNP155.vcf.hg38_multianno.txt", header = T, sep = "\t")
annovar$KEY <- paste0(annovar$Chr, ":", annovar$End)

sum(bim.extracted$KEY  %in% annovar$KEY)
# 142

df$MAF.gnomAD.genome.ALL <- annovar$gnomAD_genome_ALL[match(df$KEY, annovar$KEY)]
df$MAF.gnomAD.genome.NFE <- annovar$gnomAD_genome_NFE[match(df$KEY, annovar$KEY)]
df$MAF.1000genome <- annovar$X1000g2015aug_all[match(df$KEY, annovar$KEY)]

colnames(df)

colnames(df) <- c("CHROM", "POS", "ID", "REF", "ALT1", "ALT", "EFFECT", "IMPACT", "GENE", "GENEID", "FEATURE", "FEATUREID", "HGVS_C", "HGVS_P", "dbNSFP_CADD_phred", "dbNSFP_1000Gp3_AF", "dbNSFP_ExAC_AF", "dbNSFP_ExAC_Adj_AF",
                  "dbNSFP_MetaSVM_score", "dbNSFP_MetaSVM_rankscore", "dbNSFP_MetaSVM_pred", "dbNSFP_clinvar_clnsig", "CLNSIG", "KEY", "KEY1", "MAF.gnomAD.genome.ALL", "MAF.gnomAD.genome.NFE", "MAF.1000genome")

df <- df[, c("CHROM", "POS", "ID", "REF", "ALT", "EFFECT", "GENE", "dbNSFP_CADD_phred", "CLNSIG", "MAF.gnomAD.genome.ALL", "MAF.gnomAD.genome.NFE", "MAF.1000genome")]
dim(df)

## Part 3
# read recoded plink
raw <- read.delim("var_of_interest_plink_recodeA.raw", sep = " ", header = T, check.names = F)
colnames(raw)[grepl("^._T", colnames(raw))] <- "rs776746"
dim(raw)
head(raw)

HEADER = sub(pattern="_[T,A,G,C,*]+",replacement="",colnames(raw))
colnames(raw) = HEADER
raw <- raw[!grepl("FID|IID|PAT|MAT|SEX|PHENOTYPE", colnames(raw))]
colnames(raw) # check this with the .bim REF and NON-REFERENCE alleles
raw.t <- as.data.frame(t(raw))
raw.t$V1
raw.t$rsID <- colnames(raw)


df$WGS_193_genotype <- raw.t$V1[match(df$ID, raw.t$rsID)]
df$WGS_193_genotype[is.na(df$WGS_193_genotype)] <- ""

df$WGS_193_genotype[df$WGS_193_genotype == 0] <- paste0(df$REF[df$WGS_193_genotype == 0], df$REF[df$WGS_193_genotype == 0])
df$WGS_193_genotype[df$WGS_193_genotype == 1] <- paste0(df$REF[df$WGS_193_genotype == 1], df$ALT[df$WGS_193_genotype == 1])
df$WGS_193_genotype[df$WGS_193_genotype == 2] <- paste0(df$ALT[df$WGS_193_genotype == 2], df$ALT[df$WGS_193_genotype == 2])

df$CADD <- sub(",.*", "", df$dbNSFP_CADD_phred)  
df <- df[!is.na(df$ID),]
df$EFFECT <- gsub("_", " ", df$EFFECT)

df$CLINVAR <- df$CLNSIG
df$CLINVAR[grepl("^Pathogenic", df$CLINVAR)] <- "Pathogenic/Likely-pathogenic"
df <- df[!duplicated(df$ID),]
dim(df)


write.table(df, "Final_list_of_variants_in_WGS_193.txt", row.names = F, col.names = T, quote = F, sep = "\t")