rm(list = ls())

library(risksetROC)
library(timeROC)
library(glmnet)
library(dplyr)
library(survival)
library(hdnom)
library(RColorBrewer)
cols <- brewer.pal(3, "Dark2")
# load("results/prediction_model_uni_select_results.Rdata")
seed <- 2024
load(sprintf("results/prediction_model_cv_results_seed_%s.Rdata", seed))

eval_data_df <- analysis_data_final %>% filter(complete.cases(.))
eval_data <- eval_data_df %>% as.matrix


eval_data_y <- with(eval_data_df,
                    Surv(time = time_to_dysl, event = dyslipidemia))



clinical_model_marker <- as.numeric(eval_data[, clinical_model_covariates] %*%
                                      clinical_model_coef)

roc_clinical_model <- risksetROC(Stime = eval_data[, "time_to_dysl"],
                                 status = eval_data[, "dyslipidemia"],
                                 marker = clinical_model_marker,
                                 predict.time = 10,
                                 method = "Cox",lty = 2, col = "red")

auc_clinical_model <- risksetAUC(Stime = eval_data[, "time_to_dysl"],
                                 status = eval_data[, "dyslipidemia"],
                                 marker = clinical_model_marker,
                                 tmax = 10,
                                 method = "Cox",lty = 2, col = "red")

auc_clinical_model$Cindex # 0.7020067

prs_model_marker <- as.numeric(eval_data[, prs_model_covariates] %*%
                                 prs_model_coef)

roc_prs_model <- risksetROC(Stime = eval_data[, "time_to_dysl"],
                            status = eval_data[, "dyslipidemia"],
                            marker = prs_model_marker,
                            predict.time = 10,
                            method = "Cox",lty = 2, col = "red")

auc_prs_model <- risksetAUC(Stime = eval_data[, "time_to_dysl"],
                            status = eval_data[, "dyslipidemia"],
                            marker = prs_model_marker,
                            tmax = 10,
                            method = "Cox",lty = 2, col = "red")

auc_prs_model$Cindex ## CV method: 0.7373103

## 5/10 year AUC
clinical_model_cd_ROC <- timeROC(T = eval_data[, "time_to_dysl"],
                                 delta = eval_data[, "dyslipidemia"],
                                 marker = clinical_model_marker,
                                 # other_markers = eval_data[, eval_model_covariates],
                                 # weighting = "cox",
                                 cause = 1,
                                 times = c(5, 10), iid = T) # 74.47

prs_model_cd_ROC <- timeROC(T = eval_data[, "time_to_dysl"],
                            delta = eval_data[, "dyslipidemia"],
                            marker = prs_model_marker,
                            #other_markers = eval_data[, selected_covariates],
                            cause = 1,
                            times = c(5, 10), 
                            weighting = "marginal",
                            iid = T) # 77.87

png("ROC_5yr_04.02.png", width = 480, height = 480)
par(mar = c(5, 5, 4, 1))
plot(x = 1, y = 1, type = "n", 
     xlab = "1-Specificity", ylab = "Sensitivity",
     cex.lab = 2, cex.axis = 2, xlim = c(0, 1), ylim = c(0, 1))
plot(clinical_model_cd_ROC, time = 5, col = 1, lty = 3, lwd = 3, 
     add = T)
plot(prs_model_cd_ROC, time = 5,add = T, col = 2, lty = 1, lwd = 3)

legend("bottomright", col = 1:2, lty = c(3, 1), lwd = 3,
       legend = c("Clinical Model (AUC = 0.74)",
                  "Genetic Model (AUC = 0.78)"), cex = 1.6,
       bty = "n")
dev.off()
compare(prs_model_cd_ROC, clinical_model_cd_ROC) # 0.0007078673
## TODO:
## C/D AUC;

###### 04.02.2024
# Export to PDF with higher resolution
pdf("ROC_5yr_high_res_04.02.pdf", width = 12, height = 12)  # Set dimensions in inches
par(mar = c(8, 8, 5, 5))  # Increase margin to give more space to the plot
plot(x = 1, y = 1, type = "n",
     xlab = "1-Specificity", ylab = "Sensitivity",
     cex.lab = 1.8, cex.axis = 1.8, xlim = c(0, 1), ylim = c(0, 1))
plot(clinical_model_cd_ROC, time = 5, col = "blue", lty = 1, lwd = 5,  # Change color to blue
     add = TRUE)
plot(prs_model_cd_ROC, time = 5, add = TRUE, col = 2, lty = 1, lwd = 5)
legend("bottomright", col = c("blue", 2), lty = c(1, 1), lwd = 3,
       legend = c("Clinical Model (AUC = 0.74)",
                  paste("Clinical Model + PRSs (AUC = 0.78)\n", "p = 0.0007")), cex = 2.3,  # Adjust legend font size
       bty = "n", x.intersp = 1.5, xjust = 1)
dev.off()

#######################################
## Risk Stratification ################
#######################################
# clinical_model_hr <- exp(clinical_model_marker)
# prs_model_hr <- exp(prs_model_marker)
# 
# summary(clinical_model_hr)
# # Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# # 0.8575  0.8575  0.8711  1.0575  1.0131  3.2482 
# sum(clinical_model_hr > 3)  # 6
# 
# 
# ## Youden index
# 
# summary(prs_model_hr)
# sum(prs_model_hr > 3) # 129
# # Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# # 0.2960  0.7088  0.9028  1.2172  1.2196 10.5393
# 
# sf_clinical <- survival::survfit(clinical_model_cv,
#                                  s = clinical_model_cv$lambda.min,
#                                  x = clinical_model_x,
#                                  y = clinical_model_y,
#                                  newx = eval_data[, clinical_model_covariates])
# smr_clinical_5yr <- summary(sf_clinical, times = 5)
# 
# risk_clinical_5yr <- 1 - c(smr_clinical_5yr$surv)
# summary(risk_clinical_5yr)
# # Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# # 0.01447 0.01729 0.02212 0.02806 0.02962 0.19857 
# 
# sf_prs <- survival::survfit(prs_model_cv,
#                             s = prs_model_cv$lambda.min,
#                             x = prs_model_x,
#                             y = prs_model_y,
#                             newx = eval_data[, prs_model_covariates])
# smr_prs_5yr <- summary(sf_prs, times = 5)
# 
# risk_prs_5yr <- 1 - c(smr_prs_5yr$surv)
# summary(risk_prs_5yr)
# # Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# # 0.001866 0.009687 0.016510 0.028195 0.029492 0.559022 

clinical_risk_group <- factor(cut(clinical_model_marker, 
                                  quantile(clinical_model_marker, c(0, 0.25, 0.90, 1)),
                                  include.lowest = T, right = F),
                              labels = c("Low risk", "Medium risk", "High risk"))

prs_risk_group <- factor(cut(prs_model_marker, 
                             quantile(prs_model_marker, c(0, 0.25, 0.90, 1)),
                             include.lowest = T, right = F),
                         labels = c("Low risk", "Medium risk", "High risk"))
table(prs_risk_group, clinical_risk_group, eval_data[, "dyslipidemia"])

km_clinical <- survfit(eval_data_y ~ clinical_risk_group)
km_prs <- survfit(eval_data_y ~ prs_risk_group)

png("km_risk_stratification_04.02_t.png", width = 2400, height = 2400, res = 300)
cols <- c("green", "blue", "red")
par(mar = c(5, 5, 4, 1))
plot(km_clinical, col = cols, fun = "event", 
     lty = 2, lwd = 2, ylim = c(0, 0.25),
     xaxs = "S", xlab = "Time (years)", ylab = "Cumulative risk",
     cex.axis = 2, cex.lab = 2)
lines(km_prs, col = cols, fun = "event",
      lwd = 3)
legend("topleft", 
       legend = c("Low risk", "Moderate risk", "High risk"),
       fill = cols,
       bty = "n", cex = 1.6)
legend(x = 4.3, y = 0.5,
       lwd = 2, lty = c(1, 2),
       bty = "n",
       legend = c("Clinical Model + PRSs",
                  "Clinical Model"),
       cex = 1.6)
dev.off()

pdf("Risk_stratification_high_res_04.02.pdf", width = 12, height = 12)  # Set dimensions in inches
par(mar = c(5, 5, 4, 1))
plot(km_clinical, col = cols, fun = "event", 
     lty = 2, lwd = 3, ylim = c(0, 0.25),
     xaxs = "S", xlab = "Time (years)", ylab = "Cumulative risk",
     cex.axis = 1.8, cex.lab = 1.8)
lines(km_prs, col = cols, fun = "event",
      lwd = 3)
lines(km_prs, col = cols, fun = "event",
      lwd = 3)
legend("topleft", 
       legend = c("Low risk", "Moderate risk", "High risk"),
       fill = cols,
       bty = "n", cex = 1.8)
legend(x = 4.3, y = 0.5,
       lwd = 3, lty = c(1, 5),
       bty = "n",
       legend = c("Clinical Model + PRSs",
                  "Clinical Model"),
       cex = 2.3)
dev.off()

clinical_marker_cutoff <- quantile(clinical_model_marker, c(0, 0.25, 0.90, 1))
prs_marker_cutoff <- quantile(prs_model_marker, c(0, 0.25, 0.90, 1))
save(clinical_marker_cutoff, prs_marker_cutoff,
     file = "results/marker_cutoffs.RData")

#############################
# Report coefficients
############################
clinical_log_hr <- as.data.frame(as.matrix(clinical_model_coef[as.numeric(clinical_model_coef) != 0, ]))
prs_log_hr <- as.data.frame(as.matrix(prs_model_coef[as.numeric(prs_model_coef) != 0, ]))

write.csv(clinical_log_hr, "results/clinical_model_coef.csv")
write.csv(prs_log_hr, "results/prs_model_coef.csv")





## comparing ROC in subgroups
clinical_model_cd_ROC_subgroup <- 
  prs_model_cd_ROC_subgroup <- 
  vector("list", 3)

auc_pval <- rep(NA, 3)
## 5/10 year AUC
for (gp in 1:3) {
  risk_subset <- clinical_risk_group == c("Low risk", "Medium risk", "High risk")[gp]
  clinical_model_cd_ROC_subgroup[[gp]] <- timeROC(T = eval_data[risk_subset, "time_to_dysl"],
                                                  delta = eval_data[risk_subset, "dyslipidemia"],
                                                  marker = clinical_model_marker[risk_subset],
                                                  # other_markers = eval_data[, eval_model_covariates],
                                                  # weighting = "cox",
                                                  cause = 1,
                                                  times = 5, iid = T) 
  ## 50, 57.2, 59.2
  prs_model_cd_ROC_subgroup[[gp]] <- timeROC(T = eval_data[risk_subset, "time_to_dysl"],
                                             delta = eval_data[risk_subset, "dyslipidemia"],
                                             marker = prs_model_marker[risk_subset],
                                             #other_markers = eval_data[, selected_covariates],
                                             cause = 1,
                                             times = 5, 
                                             weighting = "marginal",
                                             iid = T) # 75.53
  ## 41, 67.5, 63.3
  auc_pval[gp] <- compare(prs_model_cd_ROC_subgroup[[gp]], clinical_model_cd_ROC_subgroup[[gp]])$p_values_AUC[2]
}

