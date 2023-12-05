# Specify the file path
setwd("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/Survivor_WES/annotation/snpEff_round2/loftee")
file_path <- "chr22.Survivor_WES.GATK4180.hg38_biallelic.vcf-annot-snpeff-dbnsfp-ExAC.0.3-clinvar.GRCh38.with.gnomAD.revel.loftee.tsv"

rl <- readLines(file_path, n=1000)
header = read.table(rl[grep('^#Uploaded_variation.*', rl)], header = T, sep = "\t")
df <- read.table(file_path, skip = 108, header = FALSE, sep ='\t')
colnames(df) <- header

##################
## Extract P/LP ##
##################
table(df$CLIN_SIG)
clinvar <- (grepl("pathogenic|likely_pathogenic", df$CLIN_SIG)

#########
## NHW ##
#########




#############
## African ##
#############