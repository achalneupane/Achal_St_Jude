library("variancePartition")
library("edgeR")
library("BiocParallel")
data(varPartDEdata)

# filter genes by number of counts
isexpr <- rowSums(cpm(countMatrix) > 0.1) >= 5

# Standard usage of limma/voom
dge <- DGEList(countMatrix[isexpr, ])
dge <- calcNormFactors(dge)

# make this vignette faster by analyzing a subset of genes
dge <- dge[1:1000, ]

# ###
# Limma Analysis
# Limma has a built-in approach for analyzing repeated measures data using duplicateCorrelation(). The model can handle a single random effect, and forces the magnitude of the random effect to be the same across all genes.

# apply duplicateCorrelation is two rounds
design <- model.matrix(~Disease, metadata)
vobj_tmp <- voom(dge, design, plot = FALSE)
dupcor <- duplicateCorrelation(vobj_tmp, design, block = metadata$Individual)

# run voom considering the duplicateCorrelation results
# in order to compute more accurate precision weights
# Otherwise, use the results from the first voom run
vobj <- voom(dge, design, plot = FALSE, block = metadata$Individual, correlation = dupcor$consensus)

# Estimate linear mixed model with a single variance component
# Fit the model for each gene,
dupcor <- duplicateCorrelation(vobj, design, block = metadata$Individual)

# But this step uses only the genome-wide average for the random effect
fitDupCor <- lmFit(vobj, design, block = metadata$Individual, correlation = dupcor$consensus)

# Fit Empirical Bayes for moderated t-statistics
fitDupCor <- eBayes(fitDupCor)
####################
## Dream analysis ##
####################

# Specify parallel processing parameters
# this is used implicitly by dream() to run in parallel
param <- SnowParam(4, "SOCK", progressbar = TRUE)

# The variable to be tested must be a fixed effect
form <- ~ Disease + (1 | Individual)

# estimate weights using linear mixed model of dream
vobjDream <- voomWithDreamWeights(dge, form, metadata, BPPARAM = param)

# Fit the dream model on each gene
# For the hypothesis testing, by default,
# dream() uses the KR method for <= 20 samples,
# otherwise it uses the Satterthwaite approximation
fitmm <- dream(vobjDream, form, metadata)
fitmm <- eBayes(fitmm)

# Examine design matrix
head(fitmm$design, 3)

# Get results of hypothesis test on coefficients of interest
topTable(fitmm, coef = "Disease1", number = 3)
