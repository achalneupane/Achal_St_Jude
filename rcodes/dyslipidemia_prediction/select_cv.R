rm(list = ls())
library(glmnet)
library(dplyr)
load("results/cv_models.Rdata")

clinical_coef <- prs_coef <- vector("list", 5)
clinical_order <- order(sapply(clinical_model_list, function(x) min(x$cvm)))
for (i in 1:5) {
  clinical_coef[[i]] <- coef(clinical_model_list[[which(clinical_order == i)]])
}


prs_order <- order(sapply(prs_model_list, function(x) min(x$cvm)))
for (i in 1:5) {
  prs_coef[[i]] <- coef(prs_model_list[[which(prs_order == i)]])
}



clinical_model_cv <- clinical_model_list[[which(clinical_order == 2)]]
clinical_model_coef <- coef(clinical_model_cv)
clinical_model_covariates <- rownames(coef(clinical_model_cv))


prs_model_cv <- prs_model_list[[which(prs_order == 4)]]
prs_model_coef <- coef(prs_model_cv)
prs_model_covariates <- rownames(coef(prs_model_cv))

save(analysis_data,
     analysis_data_final,
     clinical_model_data, 
     clinical_model_x, clinical_model_y,
     clinical_model_cv, clinical_model_covariates, clinical_model_coef,
     prs_model_data, 
     prs_model_x, prs_model_y,
     prs_model_cv, prs_model_covariates, prs_model_coef, 
     file = "results/prediction_model_cv_results.Rdata")
