rm(list = ls())
library(glmnet)
library(risksetROC)
library(lubridate)
library(survival)
library(mgcv)
library(dplyr)
library(splines)
source("../utilities.R")

set.seed(2000)

dysl_data <- readRDS("results/dysl_analytic_data.rds")

## remove those who only had baseline visit
dysl_data <- dysl_data %>%
  filter(lstcondt > baselinevisitfinish + 90)

dysl_data$age_base <- with(dysl_data,
                           yrdif(dob, baselinevisitfinish))

dysl_data$time_to_dysl <-
  with(dysl_data, 
       ifelse(dyslipidemia == 1, 
              yrdif(baselinevisitfinish, dyslipidemia_dt),
              yrdif(baselinevisitfinish, lstcondt)))


analysis_data <- dysl_data %>%
  dplyr::transmute(mrn = mrn, 
                   male = as.numeric(gender == "Male"),
                   white = as.numeric(racegrp == "White"),
                   bmi2 = as.numeric(bmiadj >= 18.5 & bmiadj < 25),
                   bmi3 = as.numeric(bmiadj >= 25 & bmiadj < 30),
                   bmi4 = as.numeric(bmiadj >= 30),
                   agedx2 = as.numeric(agedx > 5 & agedx <= 10),
                   agedx3 = as.numeric(agedx > 10 & agedx <= 15),
                   agedx4 = as.numeric(agedx > 15),
                   age_base2 = as.numeric(age_base >= 15 & age_base < 25),
                   age_base3 = as.numeric(age_base >= 25 & age_base < 35),
                   age_base4 = as.numeric(age_base >= 35 & age_base < 45),
                   age_base5 = as.numeric(age_base >= 45),
                   smoker_nowformer = as.numeric(smoker %in% c(1, 2)),
                   smoker_never = as.numeric(smoker == 3),
                   smoker_unknown = as.numeric(smoker == 4),
                   mets = mets, diabetes = diabetes, ghd = ghd, 
                   hypogonadism = hypogonadism,
                   hyperthyroidism = hyperthyroidism, 
                   cardiomyopathy = cardiomyopathy, 
                   ckd = ckd, htn = htn,
                   anthracyclines1 = as.numeric(anthracyclines_dose_5 > 0),
                   anthracyclines2 = as.numeric(anthracyclines_dose_5 >= 150),
                   cisplatin1 = as.numeric(cisplatin_dose_5 > 0),
                   cisplatin2 = as.numeric(cisplatin_dose_5 >= 400),
                   carboplatin1 = as.numeric(carboplatin_dose_5 > 0),
                   carboplatin2 = as.numeric(carboplatin_dose_5 >= 2850),
                   corticosteroids_5 = corticosteroids_5,
                   maxabdrtdose1 = as.numeric(maxabdrtdose > 200),
                   maxabdrtdose2 = as.numeric(maxabdrtdose >= 2000),
                   maxsegrtdose1 = as.numeric(maxsegrtdose > 200),
                   maxsegrtdose2 = as.numeric(maxsegrtdose >= 2500),
                   score677, score686, score688, score699, score888,
                   list_a_carrier, list_ab_carrier, list_abc_carrier,
                   time_to_dysl = time_to_dysl, 
                   dyslipidemia = dyslipidemia) 

analysis_data_x <- analysis_data %>%
  dplyr::select(-mrn, -time_to_dysl, -dyslipidemia)

analysis_data_y <- analysis_data %>%
  dplyr::select(time_to_dysl, dyslipidemia)

analysis_data_x_mean <- colMeans(analysis_data_x, na.rm = T)
analysis_data_x_sd <- apply(analysis_data_x, 2, sd, na.rm = T)

analysis_data_x_standardized <- (analysis_data_x - 
                                   matrix(analysis_data_x_mean, nrow = nrow(analysis_data_x), ncol = ncol(analysis_data_x),
                                          byrow = T)) / 
  matrix(analysis_data_x_sd, nrow = nrow(analysis_data_x), ncol = ncol(analysis_data_x),
         byrow = T)

save(analysis_data_x_mean, analysis_data_x_sd, 
     file = "results/standardization_params.RData")

analysis_data_final <- cbind(mrn = analysis_data$mrn,
                             analysis_data_x_standardized,
                             analysis_data_y)
## N = 2524
clinical_model_data <- analysis_data_final %>%
  dplyr::select(-mrn,
                -score677, -score686, -score688, 
                -score699, -score688, -score888,
                -list_a_carrier, -list_ab_carrier, -list_abc_carrier) %>%
  filter(complete.cases(.))

clinical_model_x <- clinical_model_data %>%
  dplyr::select(-time_to_dysl, -dyslipidemia) %>%
  as.matrix

clinical_model_y <- with(clinical_model_data, 
                         Surv(time = time_to_dysl,
                              event = dyslipidemia))

aa <- seq(0, 1, 0.05)
model_list <- vector("list", length = length(aa))
for (i in seq_along(aa)) {
  model_list[[i]] <- cv.glmnet(x = clinical_model_x,
                               y = clinical_model_y,
                               family = "cox",
                               type.measure = "C",
                               alpha = aa[i])
  
}

model_C <- sapply(model_list, function(mm) max(mm$cvm))
clinical_model_cv <- model_list[[which.max(model_C)]]
clinical_model_coef <- coef(clinical_model_cv)
clinical_model_covariates <- rownames(coef(clinical_model_cv))
## N = 2490
# prs_names <- c("score677", "score686", "score688", "score699", "score888")
prs_model_data <- analysis_data_final %>%
  filter(complete.cases(.))


prs_model_x <- prs_model_data %>%
  dplyr::select(-mrn, -time_to_dysl, -dyslipidemia) %>%
  as.matrix

prs_model_y <- with(prs_model_data,
                    Surv(time = time_to_dysl, event = dyslipidemia))


aa <- seq(0, 1, 0.05)
prs_model_list <- vector("list", length = length(aa))
for (i in seq_along(aa)) {
  prs_model_list[[i]] <- cv.glmnet(x = prs_model_x,
                                   y = prs_model_y,
                                   family = "cox",
                                   type.measure = "C",
                                   alpha = aa[i])
  
}
prs_model_C <- sapply(prs_model_list, function(mm) max(mm$cvm))
prs_model_cv <- prs_model_list[[which.max(prs_model_C)]]

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

