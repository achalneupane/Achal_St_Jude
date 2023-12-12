# Specify the file path
setwd("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/Survivor_WES/annotation/snpEff_round2/")
##-------------------------------------POI option 1
# Option 1: Variants in selected genes with the following 3 rare variant masks will be included :
# a.	Predicted deleterious missense variants: Annotation will be performed using SnpEff27. We will use dbNSFP28 (version 4.1a), which uses 30 in silico prediction tools for annotation . Missense variants will be classified as deleterious if >90% of collated annotations (across all tools) predict deleteriousness.
# b.	Predicted loss-of-function (LOF) variants: We will use Loss-of-Function Transcript Effect Estimator (LOFTEE; plug-in implemented in the Variant Effect Predictor or VEP29 (version 108), see https://github.com/konradjk/loftee). This tool uses VEP to annotate the most severe consequence of a given variant for each gene transcript, while LOFTEE annotates high-confidence loss-of-function (LOF) variants, which include frameshift indels, stop-gain variants and splice site disrupting variants. LOF variants flagged by LOFTEE as dubious (e.g., affecting poorly conserved exons and splice variants affecting NAGNAG sites or non-canonical splice regions) will be excluded.
# c.	Pathogenic or likely pathogenic (P/LP) variants: NCBI ClinVar30 will be accessed and search. The most recent ClinVar adjudications from clinical testing laboratories (2015 onwards) for variants without “conflicting interpretations” will be used, regardless of phenotype reported in ClinVar (given that these are sometimes vague or broad).

# version dbNSFP4.4a; downloaded on 10/25/2023
library(data.table)

Final.DF <- {}
for(i in 1:22){
  CHR=paste0("chr", i)
  print(CHR)
  file_path <- paste0("new_", CHR, ".Survivor_WES.GATK4180.hg38_biallelic.vcf-annot-snpeff-dbnsfp.vcf-FIELDS-simple.txt")
  df <- fread(file_path, header = TRUE, sep ='\t')
  df[df == "."] <- NA
  colnames(df)
  
  ## Keep unique
  # df <- df[!duplicated(df$ID),]
  ################################################
  ## a.	Predicted deleterious missense variants ##
  ################################################
  missense <- df[grepl("missense", df$`ANN[*].EFFECT`),]
  missense <- missense[!duplicated(missense$ID),]
  
  Aloft <- missense[is.na(missense$dbNSFP_Aloft_pred),]
  
  BayesDel <- missense[grepl("D", missense$dbNSFP_BayesDel_noAF_pred),] # this field may indicate the prediction of BayesDel without considering allele frequency as opposed to dbNSFP_BayesDel_addAF_pred. 
  
  missense$dbNSFP_CADD_phred <- as.numeric(missense$dbNSFP_CADD_phred)
  CADD <- missense[missense$dbNSFP_CADD_phred > 20,]
  
  ClinPred <- missense$dbNSFP_ClinPred_pred
  DANN
  DEOGEN2
  Eigen/Eigen-PC
  fathmm-MKL
  fathmm-XF
  FATHMM
  GenoCanyon
  LIST-S2
  LRT
  M-CAP
  MetaRNN
  MetaSVM/MetaLR
  Missense
  MPC
  MVP
  MutPred
  MutationAssessor
  MutationTaster
  Polyphen2
  PrimateAI
  PROVEAN
  REVEL
  SIFT
  VEST4
  
  
  
  ################
  ## c. Clinvar ##
  ################
  clinvar <- df[!grepl("uncertain|benign|conflicting|not_provided|drug_response|risk_factor|^\\-$|^\\.$", df$CLNSIG, ignore.case = T),]
  cc <- cbind.data.frame(clinvar$ID, clinvar$CLNSIG, clinvar$dbNSFP_clinvar_clnsig)
  clinvar <- clinvar[!duplicated(clinvar$ID),]


  

############
## Loftee ##
############
## file_path <- "chr22.Survivor_WES.GATK4180.hg38_biallelic.vcf-annot-snpeff-dbnsfp-ExAC.0.3-clinvar.GRCh38.with.gnomAD.revel.loftee.tsv"
#
# rl <- readLines(file_path, n=1000)
# matching_line <- grep('^#Uploaded_variation.*', rl)
# header = gsub("X.", "", names(read.table(text = rl[matching_line], header = TRUE, sep = "\t", comment.char = "")))
#
# Final.DF <- {}
# for(i in 1:22){
# CHR=paste0("chr", i)
# print(CHR)
# file_path <- paste0(CHR, ".Survivor_WES.GATK4180.hg38_biallelic.vcf-annot-snpeff-dbnsfp-ExAC.0.3-clinvar.GRCh38.with.gnomAD.revel.loftee.tsv")
# df <- read.table(file_path, skip = 108, header = FALSE, sep ='\t')
# colnames(df) <- header
# ##################
# ## Extract P/LP ##
# ##################
# patho <- as.data.frame(table(df$CLIN_SIG))
# ## Remove uncertain, benign
# clinvar <- df[!grepl("uncertain|benign|conflicting|not_provided|drug_response|risk_factor|^\\-$", df$CLIN_SIG, ignore.case = T),]
# # table(clinvar$CLIN_SIG)
# clinvar$Prediction <- "Clinvar"
# ######################
# ## Loss of Function ##
# ######################
# # LoF <- df[grepl("frameshift|start_lost|stop_gained|donor|acceptor", df$Consequence, ignore.case = T),]
# LoF <- df[grepl("frameshift|start_lost|stop_gained|gain|donor|acceptor", df$Consequence, ignore.case = T),] ## Use ANN[*].IMPACT column instead of Consequence; df$`ANN[*].IMPACT`
# # dim(LoF)
# LoF$Prediction <- "LoF"
# tmp.df <- rbind.data.frame(clinvar, LoF)
# Final.DF <- rbind.data.frame(Final.DF, tmp.df)
# print(paste0("final df rows:: ", nrow(Final.DF)))
# }
# 
# # save.image("All_p_lP_vars_before_maf_filter.Rdata")
# load("All_p_lP_vars_before_maf_filter.Rdata")
# 
# clinvar <- Final.DF[Final.DF$Prediction == "Clinvar",]
# # filter more clinvar
# clinvar <- clinvar[!grepl("^affects$|^affects,other$|^confers_sensitivity$|^other$|^protective$|^likely_risk_allele$|^association$|association,affects", clinvar$CLIN_SIG),]
# LoF <- Final.DF[Final.DF$Prediction == "LoF",]
# # Filter more for LoF
# # cc <- LoF[grepl("^non_coding_transcript_exon_variant$|^splice_region_variant,non_coding_transcript_exon_variant$", LoF$Consequence, ignore.case = T),]

###############
## Merge all ##
###############

Final.DF <- rbind.data.frame(clinvar, LoF)

####################
## All population ##
####################
Final.DF.rare.PLP <- Final.DF[Final.DF$AF < 0.01,]

saveRDS(Final.DF.rare.PLP, file = "Final.DF.rare.PLP.rds")
#########
## NHW ##
#########




#############
## African ##
#############