########################
## Finemap with CARMA ##
########################
# export LD_PRELOAD=/research_jude/rgs01_jude/groups/sapkogrp/projects/RNAseq/common/scRNAseq_Paul_cardiomyopathy/codes/yes/pkgs/mkl-2023.1.0-h213fc3f_46344/lib/libmkl_rt.so R
# devtools::install_github("ZikunY/CARMA")
lib_paths <- c(
  "/research/rgs01/home/clusterHome/aneupane/R/x86_64-pc-linux-gnu-library/4.3",
  "/research/rgs01/applications/hpcf/authorized_apps/rhel8_apps/R/4.3.1/install_with_lapack/lib64/R/library"
)
# Set the library paths
.libPaths(lib_paths)

library(CARMA)
library(data.table)
library(magrittr)
library(dplyr)
library(R.utils)
##### setting up the working directory or the wd where the data are stored
setwd('/research_jude/rgs01_jude/groups/sapkogrp/projects/Cardiotoxicity/common/gwas/CARMA/analysis/')
###### load the GWAS summary statistics
sumstat<- fread(file = "/research_jude/rgs01_jude/groups/sapkogrp/projects/Cardiotoxicity/common/gwas/CARMA/analysis/summary_stat_carma_chr16.txt",
                sep = "\t", header = T, check.names = F, data.table = F,
                stringsAsFactors = F)
dim(sumstat)
# 859
###### load the pair-wise LD matrix (assuming the variants are sorted in the same order
###### as the variants in sumstat file)
ld = fread(file = "/research_jude/rgs01_jude/groups/sapkogrp/projects/Cardiotoxicity/common/gwas/finemap_kateryna/samplesnp_gt_MAF_1_perc.ld",
           sep = " ", header = F, check.names = F, data.table = F,
           stringsAsFactors = F)
length(ld)
# 859
print(head(sumstat))
# > print(head(sumstat))
# ID CHR       POS Ref Alt        SNP     N       MAF           Z      Pval
# 1 1:200938029:G:A   1 200938029   G   A rs10494829 10000 0.2593850 -0.47490054 1.3651421
# 2 1:200938474:G:T   1 200938474   G   T  rs4915210 10000 0.2901358  0.55036595 0.5820684

# Next, we run CARMA with the input summary statistics Z and the LD matrix with the hyperparameter
# η = 1 as the default setting. Notice that given that the LD matrix is in-sample, we do not turn on the outlier
# detection ‘outlier.switch=F’.
z.list<-list()
ld.list<-list()
lambda.list<-list()
z.list[[1]]<-sumstat$Z
ld.list[[1]]<-as.matrix(ld)
lambda.list[[1]]<-1
CARMA.results<-CARMA(z.list,ld.list,lambda.list=lambda.list,
                     outlier.switch=F)

###### Posterior inclusion probability (PIP) and credible set (CS)
sumstat.result = sumstat %>% mutate(PIP = CARMA.results[[1]]$PIPs, CS = 0)
if(length(CARMA.results[[1]]$`Credible set`[[2]])!=0){
  for(l in 1:length(CARMA.results[[1]]$`Credible set`[[2]])){
    sumstat.result$CS[CARMA.results[[1]]$`Credible set`[[2]][[l]]]=l
  }
}

head(sumstat.result)

###### write the GWAS summary statistics with PIP and CS
fwrite(x = sumstat.result,
       file = "sumstats_chr16_carma_results.txt",
       sep = "\t", quote = F, na = "NA", row.names = F, col.names = T)


CARMA.results[[1]]$`Credible set`[[2]]
## list()
CARMA.results[[1]]$`Credible model`[[3]]


# Running CARMA with annotations We can include functional annotations into CARMA:
###### load the functional annotations for the variants included in GWAS summary
###### statistics (assuming the variants are sorted in the same order as the
###### variants in sumstat file)
annot=fread(file = "Sample_data/sumstats_chr1_200937832_201937832_annotations.txt.gz",
            sep="\t", header = T, check.names = F, data.table = F,
            stringsAsFactors = F)
###### z.list and ld.list stay the same with the previous setting,
###### and we add annotations this time.
annot.list<-list()
annot.list[[1]]<-annot
CARMA.results<-CARMA(z.list,ld.list,lambda.list=lambda.list,w.list=annot.list,
                     outlier.switch=F)
###### Posterior inclusion probability (PIP) and credible set (CS)
sumstat.result = sumstat %>% mutate(PIP = CARMA.results[[1]]$PIPs, CS = 0)
if(length(CARMA.results[[1]]$`Credible set`[[2]])!=0){
  for(l in 1:length(CARMA.results[[1]]$`Credible set`[[2]])){
    sumstat.result$CS[CARMA.results[[1]]$`Credible set`[[2]][[l]]]=l
  }
}
###### write the GWAS summary statistics with PIP and CS
fwrite(x = sumstat.result,
       file = "Sample_data/sumstats_chr1_200937832_201937832_carma_annot.txt.gz",
       sep = "\t", quote = F, na = "NA", row.names = F, col.names = T,
       compress = "gzip")

CARMA.results[[1]]$`Credible set`[[2]]

CARMA.results[[1]]$`Credible model`[[3]]


# Next we run CARMA with two settings: 1. without annotations, and 2. with annotations as described above.
# We use the function “CARMA” to run CARMA with the external LD and η = 1 as the default setting. Given
# by that the LD matrix is extracted from reference panel (UKBB) instead of in-sample LD, we turn on the
# outlier detection setting ‘outlier.switch=TRUE’.

###### load the GWAS summary statistics (part of AD GWAS sumstats from Jansen et al., 2019)
sumstat.1 = fread(file = "Sample_data/ADAMTS4_sumstats.txt.gz",
                  sep = "\t", header = T, check.names = F, data.table = F,
                  stringsAsFactors = F)
sumstat.2 = fread(file = "Sample_data/CR1_sumstats.txt.gz",
                  sep = "\t", header = T, check.names = F, data.table = F,
                  stringsAsFactors = F)
###### load the functional annotations for the variants included in
###### GWAS summary statistics (assuming the variants are sorted in
###### the same order as the variants in sumstat file)
annot.1 = fread(file = "Sample_data/ADAMTS4_annotations.txt.gz",
                sep = "\t", header = T, check.names = F, data.table = F,
                stringsAsFactors = F)
annot.2 = fread(file = "Sample_data/CR1_annotations.txt.gz",
                sep = "\t", header = T, check.names = F, data.table = F,
                stringsAsFactors = F)
###### load the pair-wise LD matrix (assuming the variants are sorted in
###### the same order as the variants in sumstat file)
ld.1 = fread(file = "Sample_data/ADAMTS4_ld.txt.gz",
             sep = "\t", header = F, check.names = F, data.table = F,
             stringsAsFactors = F)
ld.2 = fread(file = "Sample_data/CR1_ld.txt.gz",
             sep = "\t", header = F, check.names = F, data.table = F,
             stringsAsFactors = F)
z.list<-list()
ld.list<-list()
lambda.list<-list()
z.list[[1]]<-sumstat.1$Z
z.list[[2]]<-sumstat.2$Z
ld.list[[1]]<-as.matrix(ld.1)
ld.list[[2]]<-as.matrix(ld.2)
lambda.list[[1]]<-1
lambda.list[[2]]<-1
###### Without annotations
CARMA.results_no_annot<-CARMA(z.list,ld.list,lambda.list = lambda.list,
                              outlier.switch=T)
###### With annotations
###### Exclude the variant information columns in annotation file
###### such as positions and REF/ALT alleles.
annot.list<-list()
annot.list[[1]]<-as.matrix(cbind(1, annot.1 %>% select(-(uniqID.a1a2:SNP))))
annot.list[[2]]<-as.matrix(cbind(1, annot.2 %>% select(-(uniqID.a1a2:SNP))))
CARMA.results_annot<-CARMA(z.list,ld.list,w.list=annot.list,
                           lambda.list = lambda.list,
                           input.alpha=0, outlier.switch=T)
###### Posterior inclusion probability (PIP) and credible set (CS)
sumstat.1 = sumstat.1 %>% mutate(PIP = CARMA.results_annot[[1]]$PIPs, CS = 0)


# In the figure above, the credible sets are highlighted by colored shapes. Next, we can examine the SNPs
# included in the credible sets. For simplicity we only show the credible sets when including functional
# annotations.
carma_annot[[1]]$`Credible set` #the first element of the list
#is the result of the first locus ADAMTS4
## [1] "This is the first credible set of the locus ADAMTS4"
## CHR BP A1 A2 SNP Z PIPs
## 2184 1 161155392 A G rs4575098 6.4 1
carma_annot[[2]]$`Credible set` #the second element of the list
#is the result of the second locus CR1


# We can also examine the credible models.
carma_annot[[1]]$`Credible model`#the credible model of locus ADAMTS4

carma_annot[[2]]$`Credible model`#the credible model of locus CR1




















## Example data

########################
## Finemap with CARMA ##
########################
# export LD_PRELOAD=/research_jude/rgs01_jude/groups/sapkogrp/projects/RNAseq/common/scRNAseq_Paul_cardiomyopathy/codes/yes/pkgs/mkl-2023.1.0-h213fc3f_46344/lib/libmkl_rt.so R
# devtools::install_github("ZikunY/CARMA")
lib_paths <- c(
  "/research/rgs01/home/clusterHome/aneupane/R/x86_64-pc-linux-gnu-library/4.3",
  "/research/rgs01/applications/hpcf/authorized_apps/rhel8_apps/R/4.3.1/install_with_lapack/lib64/R/library"
)
# Set the library paths
.libPaths(lib_paths)

library(CARMA)
library(data.table)
library(magrittr)
library(dplyr)
library(R.utils)
##### setting up the working directory or the wd where the data are stored
setwd('/research_jude/rgs01_jude/groups/sapkogrp/projects/Cardiotoxicity/common/gwas/CARMA/')
###### load the GWAS summary statistics
sumstat<- fread(file = "Sample_data/sumstats_chr1_200937832_201937832.txt.gz",
                sep = "\t", header = T, check.names = F, data.table = F,
                stringsAsFactors = F)
###### load the pair-wise LD matrix (assuming the variants are sorted in the same order
###### as the variants in sumstat file)
ld = fread(file = "Sample_data/sumstats_chr1_200937832_201937832_ld.txt.gz",
           sep = "\t", header = F, check.names = F, data.table = F,
           stringsAsFactors = F)
print(head(sumstat))


# Next, we run CARMA with the input summary statistics Z and the LD matrix with the hyperparameter
# η = 1 as the default setting. Notice that given that the LD matrix is in-sample, we do not turn on the outlier
# detection ‘outlier.switch=F’.
z.list<-list()
ld.list<-list()
lambda.list<-list()
z.list[[1]]<-sumstat$Z
ld.list[[1]]<-as.matrix(ld)
lambda.list[[1]]<-1
CARMA.results<-CARMA(z.list,ld.list,lambda.list=lambda.list,
                     outlier.switch=F)
###### Posterior inclusion probability (PIP) and credible set (CS)
sumstat.result = sumstat %>% mutate(PIP = CARMA.results[[1]]$PIPs, CS = 0)
if(length(CARMA.results[[1]]$`Credible set`[[2]])!=0){
  for(l in 1:length(CARMA.results[[1]]$`Credible set`[[2]])){
    sumstat.result$CS[CARMA.results[[1]]$`Credible set`[[2]][[l]]]=l
  }
}
###### write the GWAS summary statistics with PIP and CS
fwrite(x = sumstat.result,
       file = "Sample_data/sumstats_chr1_200937832_201937832_carma.txt.gz",
       sep = "\t", quote = F, na = "NA", row.names = F, col.names = T,
       compress = "gzip")


CARMA.results[[1]]$`Credible set`[[2]]
## list()
CARMA.results[[1]]$`Credible model`[[3]]


# Running CARMA with annotations We can include functional annotations into CARMA:
###### load the functional annotations for the variants included in GWAS summary
###### statistics (assuming the variants are sorted in the same order as the
###### variants in sumstat file)
annot=fread(file = "Sample_data/sumstats_chr1_200937832_201937832_annotations.txt.gz",
              sep="\t", header = T, check.names = F, data.table = F,
              stringsAsFactors = F)
###### z.list and ld.list stay the same with the previous setting,
###### and we add annotations this time.
annot.list<-list()
annot.list[[1]]<-annot
CARMA.results<-CARMA(z.list,ld.list,lambda.list=lambda.list,w.list=annot.list,
                     outlier.switch=F)
###### Posterior inclusion probability (PIP) and credible set (CS)
sumstat.result = sumstat %>% mutate(PIP = CARMA.results[[1]]$PIPs, CS = 0)
if(length(CARMA.results[[1]]$`Credible set`[[2]])!=0){
  for(l in 1:length(CARMA.results[[1]]$`Credible set`[[2]])){
    sumstat.result$CS[CARMA.results[[1]]$`Credible set`[[2]][[l]]]=l
  }
}
###### write the GWAS summary statistics with PIP and CS
fwrite(x = sumstat.result,
       file = "Sample_data/sumstats_chr1_200937832_201937832_carma_annot.txt.gz",
       sep = "\t", quote = F, na = "NA", row.names = F, col.names = T,
       compress = "gzip")

CARMA.results[[1]]$`Credible set`[[2]]

CARMA.results[[1]]$`Credible model`[[3]]


# Next we run CARMA with two settings: 1. without annotations, and 2. with annotations as described above.
# We use the function “CARMA” to run CARMA with the external LD and η = 1 as the default setting. Given
# by that the LD matrix is extracted from reference panel (UKBB) instead of in-sample LD, we turn on the
# outlier detection setting ‘outlier.switch=TRUE’.

###### load the GWAS summary statistics (part of AD GWAS sumstats from Jansen et al., 2019)
sumstat.1 = fread(file = "Sample_data/ADAMTS4_sumstats.txt.gz",
                  sep = "\t", header = T, check.names = F, data.table = F,
                  stringsAsFactors = F)
sumstat.2 = fread(file = "Sample_data/CR1_sumstats.txt.gz",
                  sep = "\t", header = T, check.names = F, data.table = F,
                  stringsAsFactors = F)
###### load the functional annotations for the variants included in
###### GWAS summary statistics (assuming the variants are sorted in
###### the same order as the variants in sumstat file)
annot.1 = fread(file = "Sample_data/ADAMTS4_annotations.txt.gz",
                sep = "\t", header = T, check.names = F, data.table = F,
                stringsAsFactors = F)
annot.2 = fread(file = "Sample_data/CR1_annotations.txt.gz",
                sep = "\t", header = T, check.names = F, data.table = F,
                stringsAsFactors = F)
###### load the pair-wise LD matrix (assuming the variants are sorted in
###### the same order as the variants in sumstat file)
ld.1 = fread(file = "Sample_data/ADAMTS4_ld.txt.gz",
             sep = "\t", header = F, check.names = F, data.table = F,
             stringsAsFactors = F)
ld.2 = fread(file = "Sample_data/CR1_ld.txt.gz",
             sep = "\t", header = F, check.names = F, data.table = F,
             stringsAsFactors = F)
z.list<-list()
ld.list<-list()
lambda.list<-list()
z.list[[1]]<-sumstat.1$Z
z.list[[2]]<-sumstat.2$Z
ld.list[[1]]<-as.matrix(ld.1)
ld.list[[2]]<-as.matrix(ld.2)
lambda.list[[1]]<-1
lambda.list[[2]]<-1
###### Without annotations
CARMA.results_no_annot<-CARMA(z.list,ld.list,lambda.list = lambda.list,
                              outlier.switch=T)
###### With annotations
###### Exclude the variant information columns in annotation file
###### such as positions and REF/ALT alleles.
annot.list<-list()
annot.list[[1]]<-as.matrix(cbind(1, annot.1 %>% select(-(uniqID.a1a2:SNP))))
annot.list[[2]]<-as.matrix(cbind(1, annot.2 %>% select(-(uniqID.a1a2:SNP))))
CARMA.results_annot<-CARMA(z.list,ld.list,w.list=annot.list,
                           lambda.list = lambda.list,
                           input.alpha=0, outlier.switch=T)
###### Posterior inclusion probability (PIP) and credible set (CS)
sumstat.1 = sumstat.1 %>% mutate(PIP = CARMA.results_annot[[1]]$PIPs, CS = 0)


# In the figure above, the credible sets are highlighted by colored shapes. Next, we can examine the SNPs
# included in the credible sets. For simplicity we only show the credible sets when including functional
# annotations.
carma_annot[[1]]$`Credible set` #the first element of the list
#is the result of the first locus ADAMTS4
## [1] "This is the first credible set of the locus ADAMTS4"
## CHR BP A1 A2 SNP Z PIPs
## 2184 1 161155392 A G rs4575098 6.4 1
carma_annot[[2]]$`Credible set` #the second element of the list
#is the result of the second locus CR1


# We can also examine the credible models.
carma_annot[[1]]$`Credible model`#the credible model of locus ADAMTS4

carma_annot[[2]]$`Credible model`#the credible model of locus CR1