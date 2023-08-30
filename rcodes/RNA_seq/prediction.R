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
