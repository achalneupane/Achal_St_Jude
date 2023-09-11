library(glmnet)
set.seed(123)
#read in data
top_5<-read.table(file="Z:/ResearchHome/ClusterHome/aneupane/St_Jude/Achal_St_Jude/rcodes/RNA_seq/merged_94_samples_5.txt", header=TRUE)
#create training indeces from data
train_indices<-sample(1:94,0.7*94)

#create and format training and testing data
train_data<-top_5[train_indices,]
test_data<-top_5[-train_indices,]
train_labels<-train_data$Cardiomyopathy_Grade
test_labels<-test_data$Cardiomyopathy_Grade
train_data <- subset(train_data, select = -c(sjlid,gender,agedx,anthra_jco_dose_any,Cardiomyopathy_Grade,sample_age))
test_data <- subset(test_data, select = -c(sjlid,gender,agedx,anthra_jco_dose_any,Cardiomyopathy_Grade,sample_age))

#create lasso model and fit
lasso_model<-cv.glmnet(as.matrix(train_data),train_labels,alpha=1,family="binomial")
plot(lasso_model)
fit=glmnet(as.matrix(train_data),train_labels,alpha=1,lambda = lasso_model$lambda.1se,family="binomial")
lasso_genes<-coef(lasso_model,s="lambda.min")

alpha1_predicted<-predict(fit,s=lasso_model$lambda.1se,newx = as.matrix(test_data),type = "response")
mean((test_labels-alpha1_predicted)^2)
threshold <- 0.5
# Convert predicted probabilities to binary labels based on the threshold
predicted_labels <- ifelse(alpha1_predicted > threshold, 1, 0)
# Calculate accuracy
accuracy <- mean(predicted_labels == test_labels)
# Print the accuracy
cat("Accuracy:", accuracy, "\n")

#create ridge model and fit
ridge_model<-cv.glmnet(as.matrix(train_data),train_labels,alpha=0,family="binomial")
plot(ridge_model)
alpha0_fit=glmnet(as.matrix(train_data),train_labels,alpha=0,lambda = ridge_model$lambda.1se,family="binomial")
ridge_genes<-coef(ridge_model,s="lambda.min")

alpha0_predicted<-predict(alpha0_fit,s=ridge_model$lambda.1se,newx = as.matrix(test_data),type = "response")
mean((test_labels-alpha0_predicted)^2)
threshold <- 0.5
# Convert predicted probabilities to binary labels based on the threshold
predicted_labels <- ifelse(alpha0_predicted > threshold, 1, 0)
# Calculate accuracy
accuracy <- mean(predicted_labels == test_labels)
# Print the accuracy
cat("Accuracy:", accuracy, "\n")

#create elastic net model and fit
elastic_model<-cv.glmnet(as.matrix(train_data),train_labels,alpha=0.5,family="binomial")
plot(elastic_model)
alpha0.5_fit=glmnet(as.matrix(train_data),train_labels,alpha=0.5,lambda = elastic_model$lambda.1se,family="binomial")
elastic_genes<-coef(elastic_model,s="lambda.min")

alpha0.5_predicted<-predict(alpha0.5_fit,s=ridge_model$lambda.1se,newx = as.matrix(test_data),type = "response")
mean((test_labels-alpha0.5_predicted)^2)
threshold <- 0.5
# Convert predicted probabilities to binary labels based on the threshold
predicted_labels <- ifelse(alpha0.5_predicted > threshold, 1, 0)
# Calculate accuracy
accuracy <- mean(predicted_labels == test_labels)
# Print the accuracy
cat("Accuracy:", accuracy, "\n")


#Clinical Model
library(pROC)
library(DescTools)
clinical_df<-subset(top_5, select = c(sjlid,gender,agedx,anthra_jco_dose_any,Cardiomyopathy_Grade,sample_age))
clinical_model<-glm(Cardiomyopathy_Grade~sample_age+gender+anthra_jco_dose_any+agedx,data = clinical_df, family = "binomial")
clinical_model_with_genes<-glm(Cardiomyopathy_Grade~sample_age+gender+anthra_jco_dose_any+agedx+.,data = top_5, family = "binomial")
roc_clinical<-roc(top_5$Cardiomyopathy_Grade,fitted(clinical_model),levels=c(0,1))
roc_clinical_with_genes<-roc(top_5$Cardiomyopathy_Grade,fitted(clinical_model_with_genes),levels=c(0,1))
delong_test<-roc.test(roc_clinical,roc_clinical_with_genes, method="delong")
print(paste("AUC Clinical Model:", round(auc(roc_clinical), 2)))
print(paste("AUC Clinical Model + top X% genes:", round(auc(roc_clinical_with_genes), 2)))
print(paste("p-value (DeLong's Test):", round(delong_test$p.value, 4)))



#Clinical Model - dummy data
n <- 500  
X_percent <- 10  
mydata <- data.frame(
  Age = runif(n, min = 20, max = 80),
  Sex = sample(c("Male", "Female"), n, replace = TRUE),
  Anthracycline = runif(n, min = 0, max = 100),
  Radiation = sample(0:1, n, replace = TRUE),
  Cardiomyopathy = sample(0:1, n, replace = TRUE)
)


clinical_model <- glm(Cardiomyopathy ~ Age + Sex + Anthracycline + Radiation, data = mydata, family = "binomial")
mydata_with_genes <- cbind(mydata, top_genes)
clinical_model_with_genes <- glm(Cardiomyopathy ~ Age + Sex + Anthracycline + Radiation + ., data = mydata_with_genes, family = "binomial")
roc_clinical <- roc(mydata$Cardiomyopathy, fitted(clinical_model), levels = c(0, 1))
roc_clinical_with_genes <- roc(mydata$Cardiomyopathy, fitted(clinical_model_with_genes), levels = c(0, 1))
delong_test <- roc.test(roc_clinical, roc_clinical_with_genes, method = "delong")