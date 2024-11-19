# Load the data
data <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3/meta_analysis/Meta_analysis_sjlife_ccss_fixed_stratified_analysis_1.tbl", header = TRUE, sep = "\t")

# Calculate OR_CI as effect size and standard error
data$OR_CI <- paste(round(exp(data$Effect), 3), " (", round(exp(data$Effect - 1.96 * data$StdErr), 3), "-", round(exp(data$Effect + 1.96 * data$StdErr), 3), ")", sep = "")

# Select and display relevant columns
result <- data[, c("MarkerName", "P.value", "OR_CI")]
metal.SJLIFE <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3/meta_analysis/stratified_analysis_SJLIFE.txt", header = T, sep = " ")

# sor the results in the same order
result <- result[match(metal.SJLIFE$MarkerName, result$MarkerName),]
View(result)
# MarkerName  P.value              OR_CI
# 1                        rs2234962_Male 0.009596 0.697(0.531-0.916)
# 2               rs3829746_Moderate-risk 0.686000 0.905(0.559-1.466)
# 3                   rs3829746_High-risk 0.100900 0.815(0.638-1.041)
# 4  rs2234962_Anthracycline_with_heartRT 0.032650 0.796(0.646-0.981)
# 5   rs3829746_heartRT_with_Anthracyline 0.190500   0.872(0.71-1.07)
# 6                    rs3829746_Low-risk 0.028190    0.46(0.23-0.92)
# 7                rs3829746_heartRT_only 0.602500 0.911(0.641-1.294)
# 8                         rs3829746_All 0.025310 0.812(0.677-0.975)
# 9          rs3829746_Anthracycline_only 0.098100 0.737(0.513-1.058)
# 10                        rs2234962_All 0.011130 0.786(0.653-0.947)
# 11 rs3829746_Anthracycline_with_heartRT 0.032650 0.796(0.646-0.981)
# 12                  rs2234962_High-risk 0.027330 0.745(0.573-0.967)
# 13         rs2234962_Anthracycline_only 0.422500 0.865(0.608-1.232)
# 14                       rs3829746_Male 0.143300 0.822(0.632-1.069)
# 15                     rs2234962_Female 0.303700 0.873(0.674-1.131)
# 16                     rs3829746_Female 0.071230  0.789(0.61-1.021)
# 17               rs2234962_heartRT_only 0.602500 0.911(0.641-1.294)
# 18  rs2234962_heartRT_with_Anthracyline 0.190500   0.872(0.71-1.07)
# 19                   rs2234962_Low-risk 0.469600 0.819(0.477-1.406)