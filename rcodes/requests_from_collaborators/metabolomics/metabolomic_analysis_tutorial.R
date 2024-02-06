# install BiocManager if not present
if (!requireNamespace("BiocManager", quietly = TRUE)){
  install.packages("BiocManager")
}

# install structToolbox and dependencies
BiocManager::install("structToolbox")

## install additional bioc packages for vignette if needed
#BiocManager::install(c('pmp', 'ropls', 'BiocFileCache'))

## install additional CRAN packages if needed
#install.packages(c('cowplot', 'openxlsx'))

suppressPackageStartupMessages({
  # Bioconductor packages
  library(structToolbox)
  library(pmp)
  library(ropls)
  library(BiocFileCache)
  
  # CRAN libraries
  library(ggplot2)
  library(gridExtra)
  library(cowplot)
  library(openxlsx)
})


# use the BiocFileCache
bfc <- BiocFileCache(ask = FALSE)

#############
## Dataset ##
#############
# 1. ----- Iris dataset (comment if using MTBLS79 benchmark data)
D = iris_DatasetExperiment()
D$sample_meta$class = D$sample_meta$Species

## MTBLS (comment if using Iris data)
# D = MTBLS79_DatasetExperiment(filtered=TRUE)
# M = pqn_norm(qc_label='QC',factor_name='sample_type') + 
#   knn_impute(neighbours=5) +
#   glog_transform(qc_label='QC',factor_name='sample_type') +
#   filter_smeta(mode='exclude',levels='QC',factor_name='sample_type')
# M = model_apply(M,D)
# D = predicted(M)

# show info
D

# 2. ----- DatasetExperiment objects
head(D$data[,1:4])


# 3. ----- Statistical models
P = PCA(number_components=15)
P$number_components=5
P$number_components
param_ids(P)
P
