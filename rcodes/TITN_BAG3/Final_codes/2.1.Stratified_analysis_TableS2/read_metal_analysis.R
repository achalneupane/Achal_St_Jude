## Load the data
## SJLIFE + CCSS
data <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3/meta_analysis/Meta_analysis_sjlife_ccss_fixed_stratified_analysis_1.tbl", header = TRUE, sep = "\t")

# Calculate OR_CI as effect size and standard error
data$OR_CI <- paste(round(exp(data$Effect), 3), " (", round(exp(data$Effect - 1.96 * data$StdErr), 3), "-", round(exp(data$Effect + 1.96 * data$StdErr), 3), ")", sep = "")

# Select and display relevant columns
result <- data[, c("MarkerName", "P.value", "OR_CI")]
metal.SJLIFE <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3/meta_analysis/stratified_analysis_SJLIFE.txt", header = T, sep = " ")

# sor the results in the same order
result <- result[match(metal.SJLIFE$MarkerName, result$MarkerName),]
View(result)




