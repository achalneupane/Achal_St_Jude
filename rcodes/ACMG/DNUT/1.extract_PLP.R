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

## Add start and end
library(biomaRt)
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# Define gene list
genes <- as.character(ACMG$GENE)
# Query Ensembl for start and end positions
results <- getBM(attributes = c('hgnc_symbol', 'chromosome_name', 'start_position', 'end_position'), 
                 filters = 'hgnc_symbol', 
                 values = genes, 
                 mart = ensembl)

# View the results
results
results <- results[!grepl("PATCH|HSC", results$chromosome_name),]
results$chromosome_name <- paste0("chr",results$chromosome_name)

ACMG$geneStart <- results$start_position[match(ACMG$GENE, results$hgnc_symbol)]
ACMG$geneEnd <- results$end_position[match(ACMG$GENE, results$hgnc_symbol)]


setwd("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/ACMG/")

## 1. clinvar
clinvar <- read.delim("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/Survivor_WES/annotation/snpEff_round3/all_new_clinvar_P_LP.txt", header = T, sep = "\t", stringsAsFactors = F)
dim(clinvar)
clinvar$SNP <- sub(";.*", "", clinvar$ID)
clinvar$AF <- as.numeric(clinvar$AF)
clinvar$AF_nfe <- as.numeric(clinvar$AF_nfe)
clinvar$AF_afr <- as.numeric(clinvar$AF_afr)

tt <- cbind.data.frame(clinvar$ID, clinvar$AF_nfe, clinvar$AF_afr, clinvar$AF)

# make rare; .all is based on allee frequency from gnomAD global population; .eur is gnomAD EUR_nfe; .afr is gnomAD AFR
clinvar.all <- clinvar[which(clinvar$AF <0.01),]
clinvar.eur <- clinvar[which(clinvar$AF_nfe < 0.01),]
clinvar.afr <- clinvar[which(clinvar$AF_afr < 0.01),]

matching_row_final <- {}
# Loop through each row of TTN
for(i in 1:nrow(ACMG)) {
  print(paste0("Doing: ", i))
  # Find matching rows in ttn.pos where TTN$POS falls between ttn.pos$start and ttn.pos$end
  matching_row <- clinvar.eur[clinvar.eur$POS >= as.numeric(ACMG$geneStart[i]) & 
                                clinvar.eur$POS <= as.numeric(ACMG$geneEnd[i]), ]
  matching_row_final <- rbind.data.frame(matching_row_final, matching_row)
}













loftee <- read.delim("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/Survivor_WES/annotation/snpEff_round3/loftee/loftee_HC_all_chr_with_gnomad.txt", header = T, sep = "\t", stringsAsFactors = F)
dim(loftee)
loftee$SNP <- sub(";.*", "", loftee$Uploaded_variation)
# LOF variants flagged by LOFTEE as dubious (e.g., affecting poorly conserved exons and splice variants affecting NAGNAG sites or non-canonical splice regions) will be excluded.
loftee <- loftee[!grepl("NAGNAG_SITE|NON_CAN_SPLICE", loftee$LoF_flags),]
table(loftee$LoF_flags)
# - 
# 70494 


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


# make rare
snpeff.all <- snpeff[which(snpeff$AF <0.01),]
snpeff.eur <- snpeff[which(snpeff$AF_nfe < 0.01),]
snpeff.afr <- snpeff[which(snpeff$AF_afr < 0.01),]




length(unique(c(clinvar$SNP, loftee$SNP, snpeff$SNP)))
# 46076
cc <- as.data.frame(unique(c(clinvar$SNP, loftee$SNP, snpeff$SNP)))
dim(cc)
# 46076
# write.table(cc, "Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/WES_rare_variant/sjlife/rare_variants_to_extract.txt", row.names = F, col.names = F, quote = F)

cc <- cc$`unique(c(clinvar$SNP, loftee$SNP, snpeff$SNP))`

## preQC
bim.preQC <- fread("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/Survivor_WES/biallelic/plink_all/chr.ALL.Survivor_WES.GATK4180.hg38_biallelic_ID_updated.bim")
table(cc %in% bim.preQC$V2)
# TRUE 
# 46076

## QCed for GQ, DP and VQSR
bim.QC <- fread("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/Survivor_WES/biallelic2/plink_all/chr.ALL.Survivor_WES.GATK4180.hg38_biallelic_ID_updated.bim")
table(cc %in% bim.QC$V2)
# FALSE  TRUE 
# 12848 33228 

## Final SJLIFE QCed data
bim.QC.sjlife <- fread("Z:/ResearchHome/Groups/sapkogrp/projects//Genomics/common/Survivor_WES/biallelic2/plink_all/chr.ALL.Survivor_WES.GATK4180.hg38_biallelic.geno.0.1.hwe.1e-15.LCR.removed.MAC.ge.1_ID_updated.bim")
table(cc %in% bim.QC.sjlife$V2)
# FALSE  TRUE 
# 16676 29400


bim.QC.sjlife.PLP <- fread("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/WES_rare_variant//sjlife/all_rare_variants_maf0.01_all_sjlife.bim")
table(cc %in% bim.QC.sjlife.PLP$V2)
# FALSE  TRUE 
# 16818 29254
