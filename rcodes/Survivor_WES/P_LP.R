# Specify the file path
setwd("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/Survivor_WES/annotation/snpEff_round2/loftee")
file_path <- "chr22.Survivor_WES.GATK4180.hg38_biallelic.vcf-annot-snpeff-dbnsfp-ExAC.0.3-clinvar.GRCh38.with.gnomAD.revel.loftee.tsv"

rl <- readLines(file_path, n=1000)
header = read.table(rl[grep('^#Uploaded_variation.*', rl)], header = T, sep = "\t")

Final.DF <- {}
for(i in 1:22){
CHR=paste0("chr", i)
print(CHR)
file_path <- paste0(CHR, ".Survivor_WES.GATK4180.hg38_biallelic.vcf-annot-snpeff-dbnsfp-ExAC.0.3-clinvar.GRCh38.with.gnomAD.revel.loftee.tsv")
df <- read.table(file_path, skip = 108, header = FALSE, sep ='\t')
colnames(df) <- header

##################
## Extract P/LP ##
##################
patho <- as.data.frame(table(df$CLIN_SIG))
## Remove uncertain, benign
clinvar <- df[!grepl("uncertain|benign|conflicting|not_provided|drug_response|risk_factor|^\\-$", df$CLIN_SIG, ignore.case = T),]
# table(clinvar$CLIN_SIG)
clinvar$Prediction <- "Clinvar"
######################
## Loss of Function ##
######################
# LoF <- df[grepl("frameshift|start_lost|stop_gained|donor|acceptor", df$Consequence, ignore.case = T),]
LoF <- df[grepl("frameshift|start_lost|stop_gained|gain|donor|acceptor", df$Consequence, ignore.case = T),]
# dim(LoF)
LoF$Prediction <- "LoF"
tmp.df <- rbind.data.frame(clinvar, LoF)
Final.DF <- rbind.data.frame(Final.DF, tmp.df)
print(paste0("final df rows:: ", nrow(Final.DF)))
}

#########
## NHW ##
#########




#############
## African ##
#############