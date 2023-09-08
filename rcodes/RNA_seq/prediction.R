# glmnet provides various arguments for users to customize the fit: we introduce some commonly used arguments here. (For more information, type ?glmnet.)
# 
# alpha is for the elastic net mixing parameter α
# , with range α∈[0,1]
# . α=1
# is lasso regression (default) and α=0
# is ridge regression.
# 
# weights is for the observation weights, default is 1 for each observation. (Note: glmnet rescales the weights internally to sum to N, the sample size.)
# 
# nlambda is the number of λ
# values in the sequence (default is 100).
# 
# lambda can be provided if the user wants to specify the lambda sequence, but typical usage is for the program to construct the lambda sequence on its own. When automatically generated, the λ
# sequence is determined by lambda.max and lambda.min.ratio. The latter is the ratio of smallest value of the generated λ
# sequence (say lambda.min) to lambda.max. The program generates nlambda values linear on the log scale from lambda.max down to lambda.min. lambda.max is not user-specified but is computed from the input x
# and y
# : it is the smallest value for lambda such that all the coefficients are zero. For alpha = 0 (ridge) lambda.max would be ∞
# : in this case we pick a value corresponding to a small value for alpha close to zero.)
# 
# standardize is a logical flag for x variable standardization prior to fitting the model sequence. The coefficients are always returned on the original scale. Default is standardize = TRUE.
# 
# As an example, we set α=0.2
# (more like a ridge regression), and give double weight to the latter half of the observations. We set nlambda to 20 so that the model fit is only compute for 20 values of λ
# . In practice, we recommend nlambda to be 100 (default) or more. In most cases, it does not come with extra cost because of the warm-starts used in the algorithm, and for nonlinear models leads to better convergence properties.


# Load necessary libraries
library(glmnet)
library(dplyr)

# Load your gene expression data and clinical data (replace with actual paths)
# Create a dummy gene expression dataset
gene_expression_data <- data.frame(
  SampleID = paste("Sample", 1:100, sep = ""),
  Gene1 = rnorm(100, mean = 2, sd = 0.5),
  Gene2 = rnorm(100, mean = 1.5, sd = 0.3),
  Gene3 = rnorm(100, mean = 1, sd = 0.2)
)
# Create a dummy clinical dataset
clinical_data <- data.frame(
  SampleID = paste("Sample", 1:100, sep = ""),
  ClinicalFactor1 = rnorm(100, mean = 50, sd = 10),
  ClinicalFactor2 = rnorm(100, mean = 25, sd = 5),
  CardiomyopathyStatus = sample(0:1, 100, replace = TRUE)
)


# Merge gene expression data and clinical data based on a common identifier (e.g., sample ID)
merged_data <- merge(gene_expression_data, clinical_data, by = "SampleID")

# Select relevant columns (gene expression and clinical risk factors)
features <- merged_data %>%
  select(starts_with("Gene"), starts_with("ClinicalFactor1"), starts_with("ClinicalFactor2"))

# Define your outcome variable (cardiomyopathy status) and covariates (clinical risk factors)
y <- merged_data$CardiomyopathyStatus
covariates <- features %>% select(starts_with("ClinicalFactor"))


# Install and load the matrixStats package if not already done
# install.packages("matrixStats")
library(matrixStats)

# Calculate the row-wise standard deviations of gene expression data
gene_sd <- matrixStats::rowSds(as.matrix(features[, grepl("^Gene", names(features))]))

# Get the indices of top 10% and 20% genes based on standard deviations
num_genes <- length(gene_sd)
top_10_percent_idx <- order(gene_sd, decreasing = TRUE)[1:round(num_genes * 0.1)]
top_20_percent_idx <- order(gene_sd, decreasing = TRUE)[1:round(num_genes * 0.2)]

# Get the names of selected genes
top_10_percent_genes <- names(features)[top_10_percent_idx]
top_20_percent_genes <- names(features)[top_20_percent_idx]



# Filter selected genes from the features
features <- features %>%
  select(all_of(c(top_10_percent_genes, top_20_percent_genes)))

# Combine covariates and selected genes
predictors <- cbind(covariates, features)

# Split data into training and testing sets
set.seed(123)
train_indices <- sample(1:nrow(predictors), 0.7 * nrow(predictors))
train_data <- predictors[train_indices, ]
test_data <- predictors[-train_indices, ]
train_labels <- y[train_indices]
test_labels <- y[-train_indices]

# Apply Lasso and elastic net models
lasso_model <- cv.glmnet(as.matrix(train_data), train_labels, alpha = 1)
elastic_net_model <- cv.glmnet(as.matrix(train_data), train_labels, alpha = 0.5)

# Get the most influential genes using the models
lasso_genes <- coef(lasso_model, s = "lambda.min")
elastic_net_genes <- coef(elastic_net_model, s = "lambda.min")

# Print the selected genes
cat("Lasso selected genes:\n")
print(names(lasso_genes)[which(lasso_genes != 0)])
cat("\nElastic Net selected genes:\n")
print(names(elastic_net_genes)[which(elastic_net_genes != 0)])





#########
# # Yadav wanted to:
# 1. Fit a clinical model including the variables you adjusted in DE analysis. This will include age at diagnosis, sex, etc. Let’s call it Clinical model.
# 2. Add top X% genes to the clinical model in 1). Let’s call it Clinical model + top X% genes.
# 3. Calculate area under the ROC curve (AUC) for Clinical model with and without top X% genes.
# 4. Compare AUC of the Clinical model with top X% genes with the AUC of the Clinical model, using DeLong’s test.
# 5. Validate Clinical model with and without the top X% genes in (a) AA survivors and (b) additional White survivors (n=116). Please note we should validate the exact same model from 1) and 2) in the validation datasets, without re-fitting the models.
# 6. Repeat (3) and (4) for both validation datasets as well.


# interesting read: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6989986/

# Load required libraries
library(pROC)
library(DescTools)

# Generate example data
set.seed(123)
n <- 500  # Number of samples
X_percent <- 10  # Percentage of top genes to add
mydata <- data.frame(
  Age = runif(n, min = 20, max = 80),
  Sex = sample(c("Male", "Female"), n, replace = TRUE),
  Anthracycline = runif(n, min = 0, max = 100),
  Radiation = sample(0:1, n, replace = TRUE),
  Cardiomyopathy = sample(0:1, n, replace = TRUE)
)

# Generate random gene expression data (replace with your actual gene data)
top_genes <- as.data.frame(matrix(runif(n * X_percent), ncol = X_percent))

# Step 1: Fit Clinical Model
clinical_model <- glm(Cardiomyopathy ~ Age + Sex + Anthracycline + Radiation, data = mydata, family = "binomial")

# Step 2: Add top X% genes to the Clinical Model
mydata_with_genes <- cbind(mydata, top_genes)
clinical_model_with_genes <- glm(Cardiomyopathy ~ Age + Sex + Anthracycline + Radiation + ., data = mydata_with_genes, family = "binomial")

# Step 3: Calculate AUC for Clinical Model
roc_clinical <- roc(mydata$Cardiomyopathy, fitted(clinical_model), levels = c(0, 1))

# Calculate AUC for Clinical Model + top X% genes
roc_clinical_with_genes <- roc(mydata$Cardiomyopathy, fitted(clinical_model_with_genes), levels = c(0, 1))

# Step 4: Compare AUCs using DeLong's test
delong_test <- roc.test(roc_clinical, roc_clinical_with_genes, method = "delong")

# Print the results
print(paste("AUC Clinical Model:", round(auc(roc_clinical), 2)))
print(paste("AUC Clinical Model + top X% genes:", round(auc(roc_clinical_with_genes), 2)))
print(paste("p-value (DeLong's Test):", round(delong_test$p.value, 4)))
