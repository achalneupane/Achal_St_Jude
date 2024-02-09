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
colnames(intensities) <- gsub("^0","",gsub("Bur.Dis.Gra.20231206.", "", colnames(metabol)))
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
metabolites_list <- colnames(metabolites)
metabolites_list <- metabolites_list[!grepl("samples|Cardtox|Genotypes", metabolites_list)]

# Empty vector to store t-test results
p_values <- c()

# Loop through each metabolite
for (metabolite in metabolites_list) {
  # t_test_result <- t.test(metabolites[[metabolite]] ~ Cardtox, data = metabolites)
  t_test_result <- t.test(as.numeric(metabolites[[metabolite]]) ~ Cardtox, data = metabolites, na.rm = TRUE)
  p_values <- c(p_values, t_test_result$p.value)
}

names(p_values) <- metabolites_list

# Display results or further analysis with the p_values vector
metabolite_df <- data.frame(Metabolite = names(p_values), p_value = unname(p_values))


#################
## Mixed model ##
#################
cc <- metabolites
cc$dose <- sub(".*\\.", "", cc$samples)

## Install and load the necessary package
# install.packages("lme4")
library(lme4)
# ## check for collinearity. High correlation between predictors can lead to multicollinearity issues
# cor(cc[, c("xylitol", "glutamylcysteine")])

## fix colnames
colnames(cc)[!grepl("samples|Cardtox|Genotypes|dose", colnames(cc))] <-  
  paste0("M_",gsub("_*_", "_",gsub(" |-|/|\\+|\\(|\\)|,|\\'", "_", colnames(cc)[!grepl("samples|Cardtox|Genotypes|dose", colnames(cc))])))

colnames(cc)[!grepl("samples|Cardtox|Genotypes|dose", colnames(cc))] <- gsub('\\"', "", colnames(cc)[!grepl("samples|Cardtox|Genotypes|dose", colnames(cc))])


metabolite_columns <- colnames(cc) [!grepl("samples|Cardtox|Genotypes|dose", colnames(cc))]
cc[metabolite_columns] <- sapply(cc[metabolite_columns], as.numeric)


models_list <- list()
# Loop through metabolite columns and fit mixed-effects model
for (metabolite in metabolite_columns) {
  # Construct the formula
  formula <- as.formula(paste0(metabolite,"~ dose", "+ (1 | Cardtox)"))
  # Fit the mixed-effects model
  model <- lmer(formula, data = cc)
  
  # Store the model in the list
  models_list[[metabolite]] <- model
}

# ggplot(cc,aes(x=Cardtox,y=valine_norvaline,col=dose)) + geom_jitter() + geom_boxplot(alpha=0.2) + facet_wrap(~dose)
