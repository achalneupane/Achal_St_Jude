pheno <- read.table(text = "SJLID	Sample_id	EF	Cardtox	Genotypes	Cumulative anthracycline dose	Ancestry	Sex	Age_at_treatment	diaggrp	Age_at_last_contact
SJL1684113	TB-16-02775	0.449	Yes	TC	450.459	AFR	Male	11.4739726	Osteosarcoma	52.18356164
SJL1793013	TB-12-5170	0.574	No	TT	473.846	AFR	Male	14.29315068	Osteosarcoma	60.90410959
SJL2508708	TB-16-06955	0.425	Yes	TC	101.841	AFR	Female	15.13448611	Non-Hodgkin lymphoma	46.38858447
SJL4181813	TB-19-04682	0.445	Yes	TC	224.074	AFR	Female	11.075582	Osteosarcoma	34.68928063
SJL4187108	TB-14-7597	0.559	No	TT	461.765	AFR	Male	7.583561644	Non-Hodgkin lymphoma	40.32525638
SJL4726013	TB-13-6124	0.577	No	TT	291.489	AFR	Female	16.12054795	Osteosarcoma	39.86849315
SJL5049617	TB-14-0242	0.567	No	TT	365.839	AFR	Female	17.85205479	Soft tissue sarcoma	46.11780822
SJL5154909	TB-14-1537	0.472	Yes	TC	185.185	AFR	Male	2.071232877	Wilms tumor	28.93683659
SJL5188113	TB-17-04138	0.408	Yes	TC	391.228	AFR	Female	14.56164384	Osteosarcoma	35.71506849
SJL5237613	TB-15-9152	0.617	No	TT	384.211	AFR	Female	13.2739726	Osteosarcoma	27.20821918
SJL1691010	TB-14-5054	0.436	Yes	AA	181.482	EUR	Female	0.008219178	Neuroblastoma	37.14185942
SJL1741608	TB-13-4405	0.36	Yes	AC	393.329	EUR	Male	5.35890411	Non-Hodgkin lymphoma	44.47945205
SJL4726613	TB-15-10427	0.395	Yes	AC	388.75	EUR	Female	13.76712329	Osteosarcoma	42.61047983
SJL4824312	TB-16-06040	0.634	No	CC	246.342	EUR	Female	15.23444869	Ewing sarcoma family of tumors	44.03005464
SJL5113211	TB-12-4295	0.463	Yes	AA	341.177	EUR	Male	1.233692642	Liver malignancies	26.61917808
SJL5140702	TB-12-3198	0.636	No	CC	292.053	EUR	Female	17.31040497	Acute myeloid leukemia	36.19452055
SJL5152913	TB-17-00238	0.396	Yes	AC	348	EUR	Male	17.67555206	Osteosarcoma	36.45363425
SJL5195202	TB-13-2431	0.659	No	CC	306.084	EUR	Male	17.2524665	Acute myeloid leukemia	34.88534321
SJL5208908	TB-20-01253	0.704	No	CC	240.909	EUR	Male	18.6739726	Non-Hodgkin lymphoma	37.58334456
SJL5234413	TB-20-00725	0.674	No	CC	372.832	EUR	Female	20.51232877	Osteosarcoma	37.28481174", header = T, sep = "\t")


# pheno$new_agelstcont <- data.lastcondt$agelstcontact[match(pheno$SJLID, data.lastcondt$sjlid)]

pheno$ID <- sub("^[^-]+-[^-]+-", "", pheno$Sample_id)
pheno$ID <- gsub("^0+", "", pheno$ID)
pheno$ID[pheno$ID == "6955"] <- "60955"
pheno$ID[pheno$ID == "4682"] <- "4682"


intensities <- read.delim("Z:/ResearchHome/Groups/sapkogrp/projects/CAB/common/metabolomics_Northwestern/Report-Bur-Dis-20231206-CHMP_edited.txt", header = T, sep = "\t", stringsAsFactors = F)
colnames(intensities) <- gsub("^0","",gsub("Bur.Dis.Gra.20231206.", "", colnames(intensities)))
colnames(intensities) <-
c("Compounds", "KEGG.ID", "1537.0","1537.1","1537.3","7597.0","7597.1","7597.3","242.0","242.1","242.3","2775.0","2775.1","2775.3",
  "5170.0","5170.1","5170.3","9152.0","9152.1","9152.3","6124.0","6124.1","6124.3","4682.0","4682.1","4682.3","60955.0","60955.1","60955.3")

# norm.factor <- read.delim("Z:/ResearchHome/Groups/sapkogrp/projects/CAB/common/metabolomics_Northwestern/normalization_factor.txt", header = T, sep = "\t", stringsAsFactors = F)

## Log2 transformation, normalization and scaling

rownames(intensities) <- intensities$Compounds
intensities <- intensities[-c(1:2)]

intensities[] <- lapply(intensities, function(x) as.numeric(gsub(",", "", x)))


library(preprocessCore)
# intensities_log2 = log(intensities, base = 2)

# Log Transformation:
# Reason: The logarithmic transformation is often applied to data that has a
# wide range of values, helping to reduce the impact of extreme values and
# making patterns in the data more visible. Here, log2 is used to perform a
# base-2 logarithm on intensities. Note: Adding 1 before taking the logarithm is
# a common practice to handle cases where the data contains zeros. Logarithm of
# zero is undefined, so adding 1 avoids this issue.
intensities_log2 = log(intensities + 1)  # Base-2 logarithm with added 1 to handle log(0)

intensities_t <- t(intensities_log2)

# Z-Score Standardization: 
# Reason: Z-score standardization (also known as
# scaling or normalization) is applied to the transposed data. This
# transformation ensures that each variable (or sample, depending on the
# orientation after transposition) has a mean of 0 and a standard deviation of
# 1. This is useful in cases where variables have different scales, making them
# comparable in subsequent analyses.
intensities_t = as.data.frame(scale(intensities_t, center = TRUE, scale = TRUE))

## Add phenotype
intensities_t$samples <- row.names(intensities_t)
intensities_t$Cardtox <- pheno$Cardtox[match(sub("\\..*", "", intensities_t$samples), pheno$ID)]
intensities_t$Genotypes <- pheno$Genotypes[match(sub("\\..*", "", intensities_t$samples), pheno$ID)]


metabolites <- intensities_t
############
## t.test ##
############
# If there are NAs, we can remove them before running the t-test
# Replace NaN with NA in the entire dataframe
metabolites <- apply(metabolites, MARGIN = 2, FUN = function(x) ifelse(is.nan(x), NA, x))
metabolites <- as.data.frame(metabolites[, colSums(!is.na(metabolites)) > 0, drop = FALSE])
dim(metabolites)
metabolites$dose <- sub("^[^.]+\\.", "", row.names(metabolites))

## fix colnames
colnames(metabolites)[!grepl("samples|Cardtox|Genotypes|dose", colnames(metabolites))] <-  
  paste0("M_",gsub("_*_", "_",gsub(" |-|/|\\+|\\(|\\)|,|\\'", "_", colnames(metabolites)[!grepl("samples|Cardtox|Genotypes|dose", colnames(metabolites))])))

colnames(metabolites)[!grepl("samples|Cardtox|Genotypes|dose", colnames(metabolites))] <- gsub('\\"', "", colnames(metabolites)[!grepl("samples|Cardtox|Genotypes|dose", colnames(metabolites))])
colnames(metabolites) <- gsub("_$", "", colnames(metabolites))

metabolites_list <- colnames(metabolites)
metabolites_list <- metabolites_list[!grepl("samples|Cardtox|Genotypes|dose", metabolites_list)]

doses <- c(0, 1, 3)
all_p_values <- c()
for (dose in doses) {
  # Empty vector to store t-test results for the current dose
  p_values <- c()
  # Loop through each metabolite
  for (metabolite in metabolites_list) {
    metabolites.dose <- metabolites[metabolites$dose == dose,]
    
    # Perform the t-test only if there are enough observations
    t_test_result <- try(t.test(as.numeric(metabolites.dose[[metabolite]]) ~ Genotypes, data = metabolites.dose, na.rm = TRUE, paired = FALSE), silent = TRUE)
      
      # Check if t-test was successful
      if (!inherits(t_test_result, "try-error")) {
        p_values <- c(p_values, t_test_result$p.value)
      } else {
        # If an error occurs, set p-value to NA
        p_values <- c(p_values, NA)
      }
  }
  # Assign names to p_values for the current dose
  names(p_values) <- metabolites_list
  # Create a unique identifier for each dose
  dose_identifier <- paste0(dose,"dose", "_")
  # Append p_values to the all_p_values vector
  all_p_values <- c(all_p_values, setNames(p_values, paste0(dose_identifier, metabolites_list)))
}

# Display results or further analysis with the p_values vector
all_p_values <- data.frame(all_p_values)
all_p_values
all_p_values$metabolites <- sub("*.dose_", "", rownames(all_p_values))
all_p_values$dose <- sub("dose.*", "dose", rownames(all_p_values))
# View(all_p_values)
colnames(all_p_values) <- c("P", "metabolites", "dose")
all_p_values$metabolites <- sub("^M_", "", all_p_values$metabolites)

library(ggplot2)

# Convert 'dose' column to a factor with the desired order
all_p_values$dose <- factor(all_p_values$dose, levels = c("0dose", "1dose", "3dose"))

# Plot using ggplot2
ggplot(all_p_values, aes(x = dose, y = P, color = dose)) +
  geom_point(size = 3) +
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "red") +
  labs(title = "P-values for Each Dose", x = "Dose", y = "P-values") +
  scale_color_manual(values = c("0dose" = "green", "1dose" = "red", "3dose" = "darkred")) +
  scale_shape_manual(values = c("0dose" = 16, "1dose" = 17, "3dose" = 18)) +
  theme_minimal()+
  scale_y_reverse()

#################
## Mixed model ##
#################
cc <- metabolites
cc$samples <- sub("\\..*", "", cc$samples)
## Install and load the necessary package
# install.packages("lme4")
# library(lme4)
library(lmerTest)
# lmerTest is an extension of lme4 and is used for hypothesis testing of fixed effects in linear mixed-effects models.
# Specifically designed to provide p-values for fixed effects, which are not provided by default in lme4.

# ## check for collinearity. High correlation between predictors can lead to multicollinearity issues
# cor(cc[, c("xylitol", "glutamylcysteine")])

## fix colnames

metabolite_columns <- colnames(cc) [!grepl("samples|Cardtox|Genotypes|dose", colnames(cc))]
cc[metabolite_columns] <- sapply(cc[metabolite_columns], as.numeric)

# gg <- cc[c(248,250:253)]

estimates <- c()
p_values <- c()
result_df <- data.frame()

models_list <- list()
# Loop through metabolite columns and fit mixed-effects model
for (metabolite in metabolite_columns) {
  # Construct the formula
  gg <- cc[c("samples", "Cardtox", "Genotypes", "dose", metabolite)]
  gg$wanted_var <- gg[,c(metabolite)]
  
  gg <- gg %>%
    arrange(samples, dose) %>%
    group_by(samples) %>%
    mutate(baseline_dose = replace(wanted_var, dose != "0", NA)) %>%
    fill(baseline_dose, .direction = "downup")
  
  ## remove baseline samples
  gg <- gg[gg$dose != 0,]
  gg$dose <- factor(gg$dose, levels = c(1,3))
  gg$Genotypes <- factor(gg$Genotypes, levels = c("TT","TC"))
  formula <- as.formula(paste0(metabolite,"~ dose", " + baseline_dose + dose*Genotypes + (1 | Genotypes)"))
  # Fit the mixed-effects model
  model <- lmerTest::lmer(formula, data = gg, REML= T)
  
  # Store the model in the list
  models_list[[metabolite]] <- model
  summary_model <- summary(model)
  
  # Extract Estimate and Pr(>|t|) values
  estimates <- summary_model$coefficients[, "Estimate"]
  p_values <- summary_model$coefficients[, "Pr(>|t|)"]
  
  result_tmp <- cbind.data.frame(metabolite = metabolite, Estimate = estimates, P_Value = p_values)
  result_tmp$variables <- rownames(result_tmp)
  rownames(result_tmp) <- NULL
  result_df <- rbind.data.frame(result_df, result_tmp)
}

result_df$Adjusted_P_Value <- p.adjust(result_df$P_Value, method = "BH")

# Given that the baseline_dose is significant, it suggests that there is a
# significant difference in the response variable (M_1_methylnicotinamide)
# between the baseline dose (0 um) and the other doses (1 um and 3 um) after
# accounting for other variables in the model. Positive Coefficient: If the
# coefficient for baseline_dose is positive, it indicates that, on average, the
# response variable increases as the baseline dose (dose = 0 um) increases. In
# other words, individuals with a higher baseline dose tend to have higher
# values for the response variable.

# Negative Coefficient: If the coefficient for baseline_dose is negative, it
# suggests that, on average, the response variable decreases as the baseline
# dose increases. Individuals with a higher baseline dose tend to have lower
# values for the response variable.

mm.dose3 <- result_df[result_df$variables == "dose3",]
mm.baseline_dose <- result_df[result_df$variables == "baseline_dose",]
mm.Genotypes <- result_df[result_df$variables == "GenotypesTC",]
mm.dose3_Genotypes <- result_df[result_df$variables == "dose3:GenotypesTC",]


## with dose and genotype interaction
models_list <- list()
# Loop through metabolite columns and fit mixed-effects model
for (metabolite in metabolite_columns) {
  # Construct the formula
  formula <- as.formula(paste0(metabolite,"~ dose", "+" , baseline_dose, "(1 | Genotypes)"))
  # Fit the mixed-effects model
  model <- lmer(formula, data = cc)
  
  # Store the model in the list
  models_list[[metabolite]] <- model
}


## ggplot(cc,aes(x=Cardtox,y=valine_norvaline,col=dose)) + geom_jitter() + geom_boxplot(alpha=0.2) + facet_wrap(~dose)
# models_list <- list()
# # Loop through metabolite columns and fit mixed-effects model
# for (metabolite in metabolite_columns) {
#   # Construct the formula
#   formula <- as.formula(paste0(metabolite,"~ dose", "+ (1 | Genotypes)"))
#   # Fit the mixed-effects model
#   model <- lmer(formula, data = cc)
#   
#   # Store the model in the list
#   models_list[[metabolite]] <- model
# }