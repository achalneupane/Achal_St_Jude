########################
## Finemap with CARMA ##
########################

# devtools::install_github("ZikunY/CARMA")
library(CARMA)
library(data.table)
library(magrittr)
library(dplyr)
library(devtools)
library(R.utils)
##### setting up the working directory or the wd where the data are stored
setwd('CARMA')
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

