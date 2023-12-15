# Specify the file path
setwd("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/Survivor_WES/annotation/snpEff_round2/")
##-------------------------------------POI option 1
# Option 1: Variants in selected genes with the following 3 rare variant masks will be included :
# a.	Predicted deleterious missense variants: Annotation will be performed using SnpEff27. We will use dbNSFP28 (version 4.1a), which uses 30 in silico prediction tools for annotation . Missense variants will be classified as deleterious if >90% of collated annotations (across all tools) predict deleteriousness.
# b.	Predicted loss-of-function (LOF) variants: We will use Loss-of-Function Transcript Effect Estimator (LOFTEE; plug-in implemented in the Variant Effect Predictor or VEP29 (version 108), see https://github.com/konradjk/loftee). This tool uses VEP to annotate the most severe consequence of a given variant for each gene transcript, while LOFTEE annotates high-confidence loss-of-function (LOF) variants, which include frameshift indels, stop-gain variants and splice site disrupting variants. LOF variants flagged by LOFTEE as dubious (e.g., affecting poorly conserved exons and splice variants affecting NAGNAG sites or non-canonical splice regions) will be excluded.
# c.	Pathogenic or likely pathogenic (P/LP) variants: NCBI ClinVar30 will be accessed and search. The most recent ClinVar adjudications from clinical testing laboratories (2015 onwards) for variants without “conflicting interpretations” will be used, regardless of phenotype reported in ClinVar (given that these are sometimes vague or broad).

# version dbNSFP4.4a; downloaded on 10/25/2023
library(data.table)

Aloft.final <- {}
BayesDel.final <- {}
CADD.final <- {}
ClinPred.final <- {}
DANN.final <- {}
DEOGEN2.final <- {}
FATHMM.final <- {}
fathmm_MKL.final <- {}
fathmm_XF.final <- {}
GenoCanyon.final <- {}
GERP_RS.final <- {}
LIST_S2.final <- {}
LRT.final <- {}
MCAP.final <- {}
MetaRNN.final <- {}
MetaSVM.final <- {}
MetaLR.final <- {}
MPC.final <- {}
MVP.final <- {}
MutPred.final <- {}
MutationAssessor.final <- {}
MutationTaster.final <- {}
Polyphen2HDIV.final <- {}
Polyphen2HVAR.final <- {}
PrimateAI.final <- {}
PROVEAN.final <- {}
REVEL.final <- {}
SIFT.final <- {}
SIFT4G.final <- {}
VEST4.final <- {}

for(i in 1:22){
  CHR=paste0("chr", i)
  print(CHR)
  file_path <- paste0("new_", CHR, ".Survivor_WES.GATK4180.hg38_biallelic.vcf-annot-snpeff-dbnsfp_clinvar_12_10_2023_clinvar_12_10_2023.vcf-FIELDS-simple.txt")
  df <- fread(file_path, header = TRUE, sep ='\t')
  df[df == "."] <- NA
  colnames(df)
  
  ## Keep unique
  # df <- df[!duplicated(df$ID),]
  ################################################
  ## a.	Predicted deleterious missense variants ##
  ################################################
  missense <- df[grepl("missense", df$`ANN[*].EFFECT`),]
  # missense <- missense[!duplicated(missense$ID),]
  
  # annotation of loss-of-function transcripts (ALoFT) can predict the impact of premature stop variants and classify them as dominant disease-causing, recessive disease-causing and benign variants.
  # https://www.nature.com/articles/s41467-017-00443-5
  Aloft <- missense[is.na(missense$dbNSFP_Aloft_pred),]
  Aloft <- Aloft[!duplicated(Aloft$ID),]
  Aloft.final <- rbind.data.frame(Aloft.final, Aloft)
  
  # BayesDel is a deleteriousness meta-score. It works for coding and non-coding variants, single nucleotide variants and small insertion / deletions. The range of the score is from -1.29334 to 0.75731. The higher the score, the more likely the variant is pathogenic.
  # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5299048/
  BayesDel <- missense[grepl("D", missense$dbNSFP_BayesDel_noAF_pred),] # this field may indicate the prediction of BayesDel without considering allele frequency as opposed to dbNSFP_BayesDel_addAF_pred. 
  BayesDel <- BayesDel[!duplicated(BayesDel$ID),]
  BayesDel.final <- rbind.data.frame(BayesDel.final, BayesDel)
  
  # Combined Annotation Dependent Depletion
  missense$dbNSFP_CADD_phred <- as.numeric(missense$dbNSFP_CADD_phred)
  CADD <- missense[missense$dbNSFP_CADD_phred > 20,]
  CADD <- CADD[!duplicated(CADD$ID),]
  CADD.final <- rbind.data.frame(CADD.final, CADD)
  
  # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6174354/
  # Another major strength of our approach is the use of ClinVar—a rapidly growing database that allows selection of confidently annotated disease-causing variants—as a training set
  # Our classifier, ClinPred, combines random forest and gradient boosting models. As predictive features, we combine commonly used and recently developed individual prediction tool scores, as well as allele frequencies (AFs) of the variant in different populations from the gnomAD database. 
  ClinPred <- missense[grepl("D", missense$dbNSFP_ClinPred_pred),] 
  ClinPred <- ClinPred[!duplicated(ClinPred$ID),]
  ClinPred.final <- rbind.data.frame(ClinPred.final, ClinPred)
  
  # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4341060/
  # DANN trains a DNN consisting of an input layer, a sigmoid function output layer, and three 1000-node hidden layers with hyperbolic tangent activation function.
  DANN <- missense[grepl("D", missense$dbNSFP_ClinPred_pred),]
  DANN <- DANN[!duplicated(DANN$ID),]
  DANN.final <- rbind.data.frame(DANN.final, DANN)
  
  # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5570203/
  # DEOGEN2 incorporates heterogeneous information about the molecular effects of the variants, the domains involved, the relevance of the gene and the interactions in which it participates.
  DEOGEN2 <- missense[grepl("D", missense$dbNSFP_DEOGEN2_pred),] # DEOGEN2 uses the MCC-optimal deleteriousness threshold >0.45 (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5570203/)
  DEOGEN2 <- DEOGEN2[!duplicated(DEOGEN2$ID),]
  DEOGEN2.final <- rbind.data.frame(DEOGEN2.final, DEOGEN2)
  # Eigen/Eigen-PC #The Eigen scores themselves, including Eigen Phred Score (Coding) and Eigen PC Phred Score (Coding), are generally not direct indicators of deleteriousness of a genetic variant. 
  
  # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3558800/
  # Original version of the tool. Predicts whether missense variants are deleterious or tolerated based on their functional impact.
  # FATHMM: Frameshift Aware Translated Hidden Markov Models
  # FATHMM is a sequence similarity search tool that produces accurate translated alignments between protein profile hidden Markov models and DNA sequences containing frameshifts.
  FATHMM <- missense[grepl("D", missense$dbNSFP_FATHMM_pred),]
  # FATHMM <- FATHMM[!grepl("T", FATHMM$dbNSFP_FATHMM_pred),]# Remove any with tolerable
  FATHMM <- FATHMM[!duplicated(FATHMM$ID),]
  FATHMM.final <- rbind.data.frame(FATHMM.final, FATHMM)
  
  # An extension of FATHMM that incorporates multiple information sources (kernels) for improved prediction accuracy. It uses a combination of different types of genomic and functional data to make predictions.
  # https://academic.oup.com/bioinformatics/article/31/10/1536/177080?login=true
  fathmm_MKL <- missense[grepl("D", missense$dbNSFP_fathmm_MKL_coding_pred),]
  # fathmm_MKL <- fathmm_MKL[!grepl("T", fathmm_MKL$dbNSFP_FATHMM_pred),]# Remove any with tolerable
  fathmm_MKL <- fathmm_MKL[!duplicated(fathmm_MKL$ID),]
  fathmm_MKL.final <- rbind.data.frame(fathmm_MKL.final, fathmm_MKL)
  
  # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5860356/
  # A specialized version of FATHMM designed specifically for X-linked genes. It takes into account the unique characteristics of X-chromosomal inheritance.
  fathmm_XF <- missense[grepl("D", missense$dbNSFP_fathmm_XF_coding_pred),]
  fathmm_XF <- fathmm_XF[!duplicated(fathmm_XF$ID),]
  fathmm_XF.final <- rbind.data.frame(fathmm_XF.final, fathmm_XF)
  
  # used threshold of 0.5 "All of this evidence suggests that important functional segments could still be detected locally even in a generally lower-scored region". https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4444969/
  # citation: https://www.nature.com/articles/srep10576
  # Here we present GenoCanyon, a whole-genome annotation method that performs unsupervised statistical learning using 22 computational and experimental annotations thereby inferring the functional potential of each position in the human genome. With GenoCanyon, we are able to predict many of the known functional regions. The ability of predicting functional regions as well as its generalizable statistical framework makes GenoCanyon a unique and powerful tool for whole-genome annotation.
  missense$dbNSFP_GenoCanyon_score <- as.numeric(missense$dbNSFP_GenoCanyon_score)
  GenoCanyon <- missense[missense$dbNSFP_GenoCanyon_score> 0.5,] 
  GenoCanyon <- GenoCanyon[!duplicated(GenoCanyon$ID),]
  GenoCanyon.final <- rbind.data.frame(GenoCanyon.final, GenoCanyon)
  
  # http://mendel.stanford.edu/SidowLab/downloads/gerp
  # Genomic Evolutionary Rate Profiling (GERP++10) neutral rate (NR) – the number of substitutions expected under conditions of neutrality.
  # citation : https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1001025
  # source : https://www.frontiersin.org/articles/10.3389/fgene.2022.1010327/full#B7
  # GERP’s score ranges from −12.3 to 6.17, the higher the score, the more conserved that nucleotide/region and more likely to be deleterious (Davydov et al., 2010).
  # dbNSFP (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4109890/): Although GERP RS scores were typically used to measure the conservation of a nucleotide site in Mendelian disease studies (e.g. Cooper et al. 2010), an alternative measure might be a scaled RS score with the corresponding neutral rate (NR) of the site (i.e. RS/NR ratio). Therefore, both NR and RS scores were included in the database. 
  # missense$dbNSFP_GERP___NR <- as.numeric(missense$dbNSFP_GERP___NR)
  # GERP_NR <- missense[missense$dbNSFP_GERP___NR > 4,]

  # Genomic Evolutionary Rate Profiling (GERP++10) rejected substitutions (RS) – the number of substitutions expected under neutrality minus the number of substitutions ‘‘observed’’ at the position. Scores range from 1 to 6.18.The larger the score, the more conserved the site.
  missense$dbNSFP_GERP___RS <- as.numeric(missense$dbNSFP_GERP___RS)
  GERP_RS <-   missense[missense$dbNSFP_GERP___RS > 4,]
  GERP_RS <- GERP_RS[!duplicated(GERP_RS$ID),]
  GERP_RS.final <- rbind.data.frame(GERP_RS.final, GERP_RS)
  
  # https://academic.oup.com/nar/article/48/W1/W154/5827198
  # LIST-S2, a successor to our previously developed approach LIST, which aims to exploit local sequence identity and taxonomy distances in quantifying the conservation of human protein sequences. Unlike its predecessor, LIST-S2 is not limited to human sequences but can assess conservation and make predictions for sequences from any organism.
  LIST_S2 <- missense[grepl("D", missense$dbNSFP_LIST_S2_pred),]
  # LIST_S2 <- LIST_S2[!grepl("T", LIST_S2$dbNSFP_LIST_S2_pred),] # Remove any with tolerable
  LIST_S2 <- LIST_S2[!duplicated(LIST_S2$ID),]
  LIST_S2.final <- rbind.data.frame(LIST_S2.final, LIST_S2)
  
  # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2752137/
  # LRT (Likelihood Ratio Test) uses comparative genomics to identify variants that disrupt highly conserved amino acids. Variants are predicted to be deleterious ("D"), neutral ("N") or unknown ("U").
  LRT <- missense[grepl("D", missense$dbNSFP_LRT_pred),]
  LRT <- LRT[!duplicated(LRT$ID),]
  LRT.final <- rbind.data.frame(LRT.final, LRT)
  
  
  # https://www.nature.com/articles/ng.3703
  # The Mendelian Clinically Applicable Pathogenicity (M-CAP) score is a pathogenicity likelihood score that aims to misclassify no more than 5% of pathogenic variants while aggressively reducing the list of variants of uncertain significance. Much like allele frequency, M-CAP is readily interpreted; if it classifies a variant as benign, then that variant can be trusted to be benign with high confidence. M-CAP uses gradient boosting trees, a supervised learning classifier that excels at analyzing nonlinear interactions between features, and has state-of-the-art performance in a variety of classification tasks17,18,19. The features M-CAP uses for classification are based on both existing pathogenicity likelihood scores and direct measures of evolutionary conservation, the cross-species analog to frequency within the human population. 
  MCAP <- missense[grepl("D", missense$dbNSFP_M_CAP_pred),]
  MCAP <- MCAP[!duplicated(MCAP$ID),]
  MCAP.final <- rbind.data.frame(MCAP.final, MCAP)
  
  # https://genomemedicine.biomedcentral.com/articles/10.1186/s13073-022-01120-z
  # This study developed the MetaRNN and MetaRNN-indel models to overcome these limitations, enabling users to easily annotate and score both nsSNVs and nfINDELs. As predictive features, our classifiers combine recently developed independent prediction algorithms, conservation scores, and allele frequency information from the 1000 Genomes Project (1000GP) [19], ExAC [20], and gnomAD [21]
  MetaRNN <- missense[grepl("D", missense$dbNSFP_MetaRNN_pred),]
  MetaRNN <- MetaRNN[!duplicated(MetaRNN$ID),]
  MetaRNN.final <- rbind.data.frame(MetaRNN.final, MetaRNN)
  
  # https://academic.oup.com/hmg/article/24/8/2125/651446
  # we propose a meta-analytic support vector machine (Meta-SVM) that can accommodate multiple omics data, making it possible to detect consensus genes associated with diseases across studies
  MetaSVM <- missense[grepl("D", missense$dbNSFP_MetaSVM_pred),]
  MetaSVM <- MetaSVM[!duplicated(MetaSVM$ID),]
  MetaSVM.final <- rbind.data.frame(MetaSVM.final, MetaSVM)
  
  # https://academic.oup.com/hmg/article/24/8/2125/651446
  # MetaLR uses logistic regression to integrate nine independent variant deleteriousness scores and allele frequency information to predict the deleteriousness of missense variants. Variants are classified as 'tolerated' or 'damaging'; a score between 0 and 1 is also provided and variants with higher scores are more likely to be deleterious.
  # MetaLR scores are calculated by the dbNSFP project
  MetaLR <- missense[grepl("D", missense$dbNSFP_MetaLR_pred),]
  MetaLR <- MetaLR[!duplicated(MetaLR$ID),]
  MetaLR.final <- rbind.data.frame(MetaLR.final, MetaLR)
  
  # https://www.biorxiv.org/content/10.1101/148353v1.full.pdf
  # We combined information from orthogonal deleteriousness metrics into one score, called MPC (for Missense badness, PolyPhen-2, and Constraint).
  missense$dbNSFP_MPC_score <- as.numeric(missense$dbNSFP_MPC_score)
  MPC <- missense[missense$dbNSFP_MPC_score>0.5,] ## source
  MPC <- MPC[!duplicated(MPC$ID),]
  MPC.final <- rbind.data.frame(MPC.final, MPC)
  
  # Missense Variant Pathogenicity prediction: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7820281/
  # a rank score of 0.75 indicates the missense variant is more likely to be pathogenic than 75% of all possible missense variants.
  # select the max rankscore
  # Function to extract maximum number from a string
  extract_max <- function(x) {
    numbers <- as.numeric(unlist(strsplit(x, ",")))
    max_number <- max(numbers, na.rm = TRUE)
    return(max_number)
  }
  # Apply the function to each element in the vector
  missense$dbNSFP_MVP_rankscore <- sapply(missense$dbNSFP_MVP_rankscore, extract_max)
  MVP <- missense[missense$dbNSFP_MVP_rankscore>= 0.5 ,]
  MVP <- MVP[!duplicated(MVP$ID),]
  MVP.final <- rbind.data.frame(MVP.final, MVP)
  
  # MutPred2: https://www.nature.com/articles/s41467-020-19669-x
  # Statistically significant odds ratios exceeding 1.22 were observed starting at a score threshold of 0.45
  # Interpreting the results
  # http://mutpred.mutdb.org/help.html#:~:text=The%20output%20of%20MutPred2%20consists,of%200.50%20would%20suggest%20pathogenicity.
  # The output of MutPred2 consists of a general score (g), i.e., the probability that the amino acid substitution is pathogenic. This score is the average of the scores from all neural networks in MutPred2. If interpreted as a probability, a score threshold of 0.50 would suggest pathogenicity. However, in our evaluations, we have estimated that a threshold of 0.68 yields a false positive rate (fpr) of 10% and that of 0.80 yields an fpr of 5%.
  missense$dbNSFP_MutPred_score <- as.numeric(missense$dbNSFP_MutPred_score)
  MutPred <- missense[missense$dbNSFP_MutPred_score>= 0.5 ,]
  MutPred <- MutPred[!duplicated(MutPred$ID),]
  MutPred.final <- rbind.data.frame(MutPred.final, MutPred)
  
  # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3177186/
  # MutationAssessor [17] has a more elaborate conservation-based approach. It distinguishes between conservation patterns within aligned families (conservation score) and sub-families (specificity score) of homologs and so attempts to account for functional shifts between subfamilies of proteins. 
  # H, highly likely; M, moderately likely; L, less likely; N, Neutral
  MutationAssessor <- missense[grepl("H|M", missense$dbNSFP_MutationAssessor_pred),]
  MutationAssessor <- MutationAssessor[!duplicated(MutationAssessor$ID),]
  MutationAssessor.final <- rbind.data.frame(MutationAssessor.final, MutationAssessor)
  
  # https://www.nature.com/articles/nmeth0810-575
  # To meet the challenges of handling high-throughput sequencing data, we developed MutationTaster, a free, web-based application for rapid evaluation of the disease-causing potential of DNA sequence alterations. MutationTaster integrates information from different biomedical databases and uses established analysis tools
  MutationTaster <- missense[grepl("D", missense$dbNSFP_MutationTaster_pred),]
  MutationTaster <- MutationTaster[!duplicated(MutationTaster$ID),]
  MutationTaster.final <- rbind.data.frame(MutationTaster.final, MutationTaster)
  
  # Polyphen2: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2855889/
  # benign, possibly damaging, or probably damaging: http://genetics.bwh.harvard.edu/pph2/dokuwiki/overview
  # HDIV (HumDiv): This score is trained on a set of known human disease-causing mutations. It is specifically designed to predict the impact of missense mutations that contribute to Mendelian diseases.
  Polyphen2HDIV <- missense[grepl("D|P", missense$dbNSFP_Polyphen2_HDIV_pred),]
  Polyphen2HDIV <- Polyphen2HDIV[!duplicated(Polyphen2HDIV$ID),]
  Polyphen2HDIV.final <- rbind.data.frame(Polyphen2HDIV.final, Polyphen2HDIV)
  
  # benign, possibly damaging, or probably damaging: http://genetics.bwh.harvard.edu/pph2/dokuwiki/overview
  # HVAR (HumVar): This score is trained on a larger set of human genetic variations, including both disease-causing mutations and common polymorphisms. It is intended to predict the impact of missense mutations across a broader spectrum of genetic variation.
  Polyphen2HVAR <- missense[grepl("D|P", missense$dbNSFP_Polyphen2_HVAR_pred),]
  Polyphen2HVAR <- Polyphen2HVAR[!duplicated(Polyphen2HVAR$ID),]
  Polyphen2HVAR.final <- rbind.data.frame(Polyphen2HVAR.final, Polyphen2HVAR)
  
  # For practical application of PrimateAI scores, we recommend a threshold of > 0.8 for likely pathogenic classification, < 0.6 for likely benign, and 0.6–0.8 as intermediate, based on the enrichment of de novo variants in cases compared to controls (Fig. 3d).
  # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6237276/
  PrimateAI <- missense[grepl("D", missense$dbNSFP_PrimateAI_pred),]
  PrimateAI <- PrimateAI[!duplicated(PrimateAI$ID),]
  PrimateAI.final <- rbind.data.frame(PrimateAI.final, PrimateAI)
  
  # PROVEAN (Protein Variation Effect Analyzer)
  # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4528627/
  # If the PROVEAN score is equal to or below a predefined threshold (e.g. -2.5), the protein variant is predicted to have a "deleterious" effect.
  PROVEAN <- missense[grepl("D", missense$dbNSFP_PROVEAN_pred),]
  PROVEAN <- PROVEAN[!duplicated(PROVEAN$ID),]
  PROVEAN.final <- rbind.data.frame(PROVEAN.final, PROVEAN)

  # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5065685/
  # REVEL score above 0.5, corresponding to a sensitivity of 0.754 and specificity of 0.891. Selecting a more stringent REVEL score threshold of 0.75 would result in higher specificity but lower sensitivity, with 52.1% of disease mutations, 3.3% of neutral variants, and 4.1% of all ESVs being classified as pathogenic
  REVEL <-  missense[missense$dbNSFP_REVEL_score> 0.5,]
  REVEL <- REVEL[!duplicated(REVEL$ID),]
  REVEL.final <- rbind.data.frame(REVEL.final, REVEL)
  
  
  # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC168916/
  # SIFT (Sorting Intolerant From Tolerant) is a program that predicts whether an amino acid substitution affects protein function so that users can prioritize substitutions for further study. We have shown that SIFT can distinguish between functionally neutral and deleterious amino acid changes in mutagenesis studies and on human polymorphisms.
  SIFT <- missense[grepl("D", missense$dbNSFP_SIFT_pred),]
  SIFT <- SIFT[!duplicated(SIFT$ID),]
  SIFT.final <- rbind.data.frame(SIFT.final, SIFT)
  
  
  # SIFT 4G (SIFT for genomes), which is a faster version of SIFT that enables practical computations on reference genomes
  # https://www.nature.com/articles/nprot.2015.123
  SIFT4G <- missense[grepl("D", missense$dbNSFP_SIFT4G_pred),]
  SIFT4G <- SIFT4G[!duplicated(SIFT4G$ID),]
  SIFT4G.final <- rbind.data.frame(SIFT4G.final, SIFT4G)
  
  
  # We have developed the Variant Effect Scoring Tool (VEST), a supervised machine learning-based classifier, to prioritize rare missense variants with likely involvement in human disease.
  # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3665549/
  # For VEST, we classified the variants based on the VEST score with a cut-off 0.5, below which the variants were classified as benign and otherwise harmful.
  # Source: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6386394/#:~:text=The%20score%20for%20VEST%20indicates,VEST%20score%20cutoff%20of%200.5.
  VEST4 <- missense[missense$dbNSFP_VEST4_score > 0.5,]
  VEST4 <- VEST4[!duplicated(VEST4$ID),]
  VEST4.final <- rbind.data.frame(VEST4.final, VEST4)
  
}
  
  
  
cc <- cbind.data.frame(missense[, 1:12], missense$dbNSFP_MetaLR_pred, missense$dbNSFP_MetaLR_score)
  
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