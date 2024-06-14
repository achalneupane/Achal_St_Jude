########################################################################################################################################################
# 1. Prevalence of Carrying Rare P/LP Variants
# To evaluate the prevalence of carrying rare P/LP variants in BCC panel genes among BCC unaffected/affected survivors, you can use simple summary statistics and proportions.

# Load necessary libraries
library(dplyr)

# Example data
set.seed(123)
data <- data.frame(
  status = c(rep("affected", 50), rep("unaffected", 50)),
  rare_variant = c(sample(c(0, 1), 50, replace = TRUE), sample(c(0, 1), 50, replace = TRUE))
)

# Calculate prevalence
prevalence <- data %>%
  group_by(status) %>%
  summarize(prevalence = mean(rare_variant))

print(prevalence)

########################################################################################################################################################
# 2. Cox Proportional Hazards Regression
# To evaluate PRS associations with genetic factors using Cox proportional hazards regression, adjust for sex, age at childhood cancer diagnosis, genetic ancestry, and maximum RT dose.

# Load necessary libraries
library(survival)

# Example data
set.seed(123)
data <- data.frame(
  time = sample(1:100, 100, replace = TRUE),
  status = sample(0:1, 100, replace = TRUE),
  PRS = rnorm(100),
  sex = sample(c("male", "female"), 100, replace = TRUE),
  age_at_diagnosis = sample(1:18, 100, replace = TRUE),
  ancestry = sample(1:3, 100, replace = TRUE),
  RT_dose = sample(1:50, 100, replace = TRUE)
)

# Fit Cox model
cox_model <- coxph(Surv(time, status) ~ PRS + sex + age_at_diagnosis + ancestry + RT_dose, data = data)
summary(cox_model)

########################################################################################################################################################
# 3. Cross-Validated Prediction Performance Metrics Using XGBoost
# To include PRS to model with BCC predicted risk from novel age-specific models, you can use the XGBoost algorithm and evaluate using cross-validation.

# Load necessary libraries
library(xgboost)
library(caret)
library(pROC)

# Example data
set.seed(123)
data <- data.frame(
  status = sample(0:1, 100, replace = TRUE),
  PRS = rnorm(100),
  age = sample(18:70, 100, replace = TRUE),
  sex = sample(c("male", "female"), 100, replace = TRUE),
  ancestry = sample(1:3, 100, replace = TRUE),
  RT_dose = sample(1:50, 100, replace = TRUE)
)

# Prepare data for XGBoost
data_matrix <- model.matrix(status ~ . - 1, data = data)
label <- data$status

# Train-test split
set.seed(123)
train_index <- createDataPartition(label, p = 0.8, list = FALSE)
train_data <- data_matrix[train_index, ]
train_label <- label[train_index]
test_data <- data_matrix[-train_index, ]
test_label <- label[-train_index]

# Convert to DMatrix
dtrain <- xgb.DMatrix(data = train_data, label = train_label)
dtest <- xgb.DMatrix(data = test_data, label = test_label)

# Train model
params <- list(objective = "binary:logistic", eval_metric = "auc")
bst_model <- xgb.train(params, dtrain, nrounds = 100)

# Predict and evaluate
preds <- predict(bst_model, dtest)
roc_curve <- roc(test_label, preds)
auc <- auc(roc_curve)

print(auc)

########################################################################################################################################################
# 4. Comparing Models Using Bootstrapped 2-Sided P-Values
# To compare prediction performance using bootstrapped 2-sided p-values, you can use bootstrapping techniques to assess statistical significance.

# Load necessary libraries
library(boot)

# Define function to calculate AUC
auc_function <- function(data, indices) {
  dtrain <- xgb.DMatrix(data = data[indices, -1], label = data[indices, 1])
  bst_model <- xgb.train(params, dtrain, nrounds = 100)
  preds <- predict(bst_model, dtrain)
  roc_curve <- roc(data[indices, 1], preds)
  auc(roc_curve)
}

# Bootstrapping
data_combined <- cbind(label, data_matrix)
boot_results <- boot(data_combined, auc_function, R = 1000)

# Calculate 2-sided p-values
p_value <- 2 * min(mean(boot_results$t >= mean(boot_results$t0)), mean(boot_results$t <= mean(boot_results$t0)))

print(p_value)