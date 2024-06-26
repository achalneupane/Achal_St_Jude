pathway.genes <- c("PROCR", "CCCA", "CCD41", "EPCR",
           "F5", "FVL", "PCCF", "RPRGL1", "THPH2",
           "SERPINC1", "AT3", "AT3D", "ATIII", "ATIII-R2", "ATIII-T1", "ATIII-T2", "THPH7",
           "SERPINA5", "PAI-3", "PAI3", "PCI", "PCI-B", "PLANH3", "PROCI",
           "PROC", "APC", "PC", "PROC1", "THPH3", "THPH4",
           "PROS1", "PROS", "PS21", "PS22", "PS23", "PS24", "PS25", "PSA", "THPH5", "THPH6")


tier1.genes <- c(
  "F12", "F10", "F11", "F12", "F13A1", "F13B", "F2", "F5", "F7", "F8", "F9",
  "FGA", "FGB", "FGG", "GGCX", "KNG1", "LMAN1", "MCFD2", "SERPINE1", "SERPINF2",
  "THBD", "VKORC1", "VWF", "ADAMTS13", "F2", "F5", "HRG", "PIGA", "PLG", "PROC",
  "PROS1", "SERPINC1", "SERPIND1", "THBD", "ABCG5", "ABCG8", "ACTB", "ACTN1",
  "ANKRD26", "ANO6", "AP3B1", "AP3D1", "ARPC1B", "BLOC1S3", "BLOC1S6", "CDC42",
  "CYCS", "DIAPH1", "DTNBP1", "ETV6", "FERMT3", "FLI1", "FLNA", "FYB1", "GATA1",
  "GFI1B", "GNE", "GP1BA", "GP1BA", "GP1BA", "GP1BB", "GP1BB", "GP6", "GP9",
  "HOXA11", "HPS1", "HPS3", "HPS4", "HPS5", "HPS6", "ITGA2B", "ITGA2B", "ITGB3",
  "ITGB3", "KDSR", "LYST", "MECOM", "MPIG6B", "MPL", "MYH9", "NBEA", "NBEAL2",
  "P2RY12", "PLA2G4A", "PLAU", "RASGRP2", "RBM8A", "RNU4ATAC", "RUNX1",
  "SLFN14", "SRC", "STIM1", "STXBP2", "TBXA2R", "TBXAS1", "THPO", "TUBB1",
  "VIPAS39", "VPS33B", "VWF", "WAS"
)

genes <- c(tier1.genes, pathway.genes)
genes <- unique(genes)
length(genes)
# 127

## Read PLP genes
clinvar <- read.delim("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/MERGED_sjlife1_2_PreQC/cleaned/annotation/snpEff/all_new_clinvar_P_LP.txt", header = T, sep = "\t", stringsAsFactors = F)
clinvar$ID <- sub(";.*", "", clinvar$ID)
clinvar <- clinvar[!duplicated(clinvar$ID),]
sum(clinvar$ANN....GENE %in% genes)
# 126

clinvar <- clinvar[clinvar$ANN....GENE %in% genes,]

extracted <- cbind.data.frame(CHROM=clinvar$CHROM, POS=clinvar$POS, SNPID=clinvar$ID, GENE=clinvar$ANN....GENE, CLINVAR=clinvar$CLNSIG, TYPE=clinvar$ANN....EFFECT, AF=clinvar$AF, AF_nfe=clinvar$AF_nfe, AF_afr=clinvar$AF_afr)

library(writexl)
write_xlsx(extracted, "Z:/ResearchHome/ClusterHome/aneupane/data/request/Rohith_jesudas/extracted_variants.xlsx", col_names = TRUE)


## Extract PLP variants in raw plink format
write.table(extracted$SNPID, "Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/MERGED_sjlife1_2_PreQC/cleaned/rohit_jesudas/rohith_jesudas_variants.txt", quote = F, col.names = F, row.names = F)

# UNIX command
# cd /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/MERGED_sjlife1_2_PreQC/cleaned/rohit_jesudas
# module load plink/1.90b
# for i in {1..22}; do
# plink --bfile /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/MERGED_sjlife1_2_PreQC/cleaned/MERGED.SJLIFE.1.2.GATKv3.4.VQSR.chr${i}.preQC_biallelic_renamed_ID_edited.vcf.gz --extract rohith_jesudas_variants.txt --keep-allele-order --make-bed --out tmp_clinvar_chr${i}.dat
# done
# 
# ## EUR
# ls -1 tmp_clinvar_chr*.dat*bed| sort -V | sed 's/\.bed//' > merge_list.txt
# plink --merge-list merge_list.txt --keep-allele-order --make-bed --out merged_data_plink
# plink --bfile merged_data_plink --recodeA --out merged_data_plink_recodeA

library(data.table)
raw <- as.data.frame(fread("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/MERGED_sjlife1_2_PreQC/cleaned/rohit_jesudas/merged_data_plink_recodeA.raw", header = T))
rownames(raw) <- raw$IID
raw <- raw[,-c(1:6)]
HEADER = sub(pattern="_[T,A,G,C,*]+",replacement="",colnames(raw))
colnames(raw) <- gsub("\\.", ":", HEADER)
raw  <- cbind.data.frame(samples= rownames(raw), raw)
write_xlsx(raw, "Z:/ResearchHome/ClusterHome/aneupane/data/request/Rohith_jesudas/carrier_status_extracted_variants.xlsx", col_names = TRUE)
