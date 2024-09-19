## Rare variant P/LP analysis with fisher exact test

#########################
## Load carrier status ##
#########################
load("Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3/Rcodes/EUR.dat.PLP.RData")
genes <- colnames(EUR.dat.PLP)[grepl("carrier", colnames(EUR.dat.PLP))]
## EUR analysis

# Empty dataframe
results <- data.frame(
  Gene = character(),
  TotalCarriers = numeric(),
  TotalCases = numeric(),
  CarriersCases = numeric(),
  TotalControls = numeric(),
  CarriersControls = numeric(),
  OR_CI = character(),
  PValue = numeric(),
  stringsAsFactors = FALSE
)

# Looping through each gene set/carrier variables
genes <- colnames(EUR.dat.PLP)[grepl("carrier", colnames(EUR.dat.PLP))]

for (i in 1:length(genes)){
  EUR.dat.PLP$carriers <- EUR.dat.PLP[, genes[i]]
  
  # Try to run Fisher's test, handle errors if no carriers
  tryCatch({
    gene.test <- fisher.test(table(EUR.dat.PLP$CMP, EUR.dat.PLP$carriers))
    pvalue <- gene.test$p.value
    OR.CI <- paste0(round(gene.test$estimate, 2), " (", paste0(round(gene.test$conf.int, 2), collapse = "-"), ")")
  }, error = function(e) {
    pvalue <- NA
    OR.CI <- NA
  })
  
  total.carriers <- sum(EUR.dat.PLP$carriers == 1, na.rm = TRUE)
  total.cases <- sum(EUR.dat.PLP$CMP == 2, na.rm = TRUE)
  carriers.cases <- sum(EUR.dat.PLP$CMP == 2 & EUR.dat.PLP$carriers == 1, na.rm = TRUE)
  total.controls <- sum(EUR.dat.PLP$CMP == 1, na.rm = TRUE)
  carriers.controls <- sum(EUR.dat.PLP$CMP == 1 & EUR.dat.PLP$carriers == 1, na.rm = TRUE)
  
  # Append the results 
  results <- rbind(results, data.frame(
    Gene = genes[i],
    TotalCarriers = total.carriers,
    TotalCases = total.cases,
    CarriersCases = carriers.cases,
    TotalControls = total.controls,
    CarriersControls = carriers.controls,
    OR_CI = OR.CI,
    PValue = pvalue
  ))
}

# View the results
print(results)

## AFR analysis
load("Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3/Rcodes/AFR.dat.PLP.RData")
genes <- colnames(AFR.dat.PLP)[grepl("carrier", colnames(AFR.dat.PLP))]
AFR.dat.PLP$CMP <- ifelse(AFR.dat.PLP$CMP==1, 2, 1)
## AFR analysis

# Empty dataframe
results <- data.frame(
  Gene = character(),
  TotalCarriers = numeric(),
  TotalCases = numeric(),
  CarriersCases = numeric(),
  TotalControls = numeric(),
  CarriersControls = numeric(),
  OR_CI = character(),
  PValue = numeric(),
  stringsAsFactors = FALSE
)

# Looping through each gene set/carrier variables
genes <- colnames(AFR.dat.PLP)[grepl("carrier", colnames(AFR.dat.PLP))]

for (i in 1:length(genes)){
  AFR.dat.PLP$carriers <- AFR.dat.PLP[, genes[i]]
  
  # Try to run Fisher's test, handle errors if no carriers
  tryCatch({
    gene.test <- fisher.test(table(AFR.dat.PLP$CMP, AFR.dat.PLP$carriers))
    pvalue <- gene.test$p.value
    OR.CI <- paste0(round(gene.test$estimate, 2), " (", paste0(round(gene.test$conf.int, 2), collapse = "-"), ")")
  }, error = function(e) {
    pvalue <- NA
    OR.CI <- NA
  })
  
  total.carriers <- sum(AFR.dat.PLP$carriers == 1, na.rm = TRUE)
  total.cases <- sum(AFR.dat.PLP$CMP == 2, na.rm = TRUE)
  carriers.cases <- sum(AFR.dat.PLP$CMP == 2 & AFR.dat.PLP$carriers == 1, na.rm = TRUE)
  total.controls <- sum(AFR.dat.PLP$CMP == 1, na.rm = TRUE)
  carriers.controls <- sum(AFR.dat.PLP$CMP == 1 & AFR.dat.PLP$carriers == 1, na.rm = TRUE)
  
  # Append the results 
  results <- rbind(results, data.frame(
    Gene = genes[i],
    TotalCarriers = total.carriers,
    TotalCases = total.cases,
    CarriersCases = carriers.cases,
    TotalControls = total.controls,
    CarriersControls = carriers.controls,
    OR_CI = OR.CI,
    PValue = pvalue
  ))
}

# View the results
print(results)
