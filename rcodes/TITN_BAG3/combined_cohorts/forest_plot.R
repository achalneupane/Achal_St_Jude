
echo.sig <- read.table(text = "Phenotype	β	SE	P	n	P_BH	Gender
LV_ejection_fraction	0.121	0.047	0.010	1232	0.042	All
LV_end-diastolic	-0.118	0.046	0.010	1232	0.042	All
LV_end-systolic	-0.149	0.046	0.001	1232	0.020	All
LV_global_longitudinal_peak_strain_(GLPS)	-0.126	0.043	0.004	1359	0.029	All
LV_ejection_fraction	0.169	0.062	0.006	677	0.025	Male
LV_end-diastolic	-0.172	0.061	0.005	677	0.025	Male
LV_end-systolic	-0.202	0.061	0.001	676	0.017	Male
LV_global_longitudinal_peak_strain_(GLPS)	-0.165	0.059	0.005	750	0.025	Male", header = T)

echo.sig$Phenotype <- gsub("_|-", " ", echo.sig$Phenotype)

# echo.sig$Significant = ifelse(echo.sig$P_BH < 0.05, "*", "")
echo.sig$OR <- exp(echo.sig$β)
echo.sig$CI_low <- exp(echo.sig$β - 1.96 * echo.sig$SE)
echo.sig$CI_high <- exp(echo.sig$β + 1.96 * echo.sig$SE)


# Create the forest plot
ggplot(echo.sig, aes(x = β, y = reorder(Phenotype, β))) +
  # geom_point(aes(color = Significant)) +
  # geom_errorbarh(aes(xmin = β - 1.96 * SE, xmax = β + 1.96 * SE)) +
  geom_linerange(aes(xmin = β - 1.96 * SE, xmax = β + 1.96 * SE)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray") +
  labs(x = "Effect Size (β)", y = "Echo measures") +
  scale_color_manual(values = c("black", "red")) +
  geom_text(aes(label = sprintf("p = %.3f\nOR = %.3f", P_BH, OR)),
            hjust = -0.1, vjust = 0.5, size = 3) +
  theme_classic() +
  facet_wrap(~ Gender, nrow = 1)



ggplot(echo.sig, aes(x = Phenotype, y = β, color = Gender)) +
  geom_linerange(aes(ymin = β - 1.96 * SE, ymax = β + 1.96 * SE), position = position_dodge(width = 0.8)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
  labs(x = "Echo measures", y = "Effect Size (β)") +
  geom_text(aes(label = sprintf("p = %.3f\nOR = %.3f", P_BH, OR)),
            vjust = -0.5, size = 3) +
  theme_classic() +
  scale_color_manual(values = c("Male" = "blue", "All" = "black")) +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
  theme(legend.position = "top")


library(ggplot2)

library(ggplot2)

ggplot(echo.sig, aes(x = Phenotype, y = β, color = Gender)) +
  geom_linerange(aes(ymin = β - 1.96 * SE, ymax = β + 1.96 * SE), position = position_dodge(width = 0.8)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
  labs(x = NULL, y = "Effect Size (β)") +
  geom_text(aes(x = as.numeric(factor(Phenotype)) - 0.2, label = sprintf("P = %.3f", P_BH)),
            vjust = -0.5, size = 3, color = "black") +
  geom_text(aes(x = as.numeric(factor(Phenotype)) + 0.2, label = sprintf("P = %.3f", P_BH)),
            vjust = -0.5, size = 3, color = "blue") +
  geom_point(aes(x = as.numeric(factor(Phenotype)) - 0.2, y = β, color = "black"), size = 3) +
  geom_point(aes(x = as.numeric(factor(Phenotype)) + 0.2, y = β, color = "blue"), size = 3) +
  theme_classic() +
  scale_color_manual(values = c("Male" = "blue", "All" = "black")) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(legend.position = "top") +
  scale_x_discrete(labels = NULL) +
  scale_color_manual(values = c("black" = "black", "blue" = "blue"))


library(ggplot2)

ggplot(echo.sig, aes(x = Phenotype, y = β, color = Gender)) +
  geom_linerange(aes(ymin = β - 1.96 * SE, ymax = β + 1.96 * SE), position = position_dodge(width = 0.8)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
  labs(x = NULL, y = "Effect Size (β)") +
  geom_text(aes(x = as.numeric(factor(Phenotype)), label = sprintf("P = %.3f\nOR = %.3f", P_BH, OR)),
            vjust = -0.5, size = 3, color = ifelse(Gender == "Male", "blue", "black")) +
  geom_point(size = 3) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(legend.position = "top") +
  scale_x_discrete(labels = NULL) +
  scale_color_manual(values = c("Male" = "blue", "All" = "black"))



##############



















# library(ggplot2)
# data <- data.frame(
#   Phenotype = c("LV ejection fraction", "LV end-diastolic volume", "LV end-systolic volume", 
#                 "LV stroke volume", "LV mass index", "LV mass", "LA systolic volume", 
#                 "LA volume index", "TAPSE in mm", "LV global longitudinal peak strain", 
#                 "LV relative wall thickness", "3D LV cardiac output", 
#                 "LV posterior wall thickness", "LA systolic diameter", 
#                 "LV diastolic diameter", "Interventricular septal thickness"),
#   Beta = c(0.121, -0.118, -0.149, -0.057, -0.023, -0.017, 0.011, 0.033, 0.043, 
#            -0.126, 0.053, 0.018, 0.068, 0.066, 0.011, 0.037),
#   SE = c(0.047, 0.046, 0.046, 0.046, 0.039, 0.039, 0.044, 0.047, 0.050, 0.043, 
#          0.041, 0.048, 0.042, 0.042, 0.040, 0.041),
#   P_Value = c(0.042, 0.042, 0.020, 0.432, 0.753, 0.806, 0.806, 0.701, 0.617, 0.029,
#               0.432, 0.806, 0.328, 0.328, 0.806, 0.617)
# )
# 
# # Define a threshold for significance (e.g., p < 0.05)
# significance_threshold <- 0.05
# 
# # Create a column for significance
# data$Significant <- ifelse(data$P_Value < significance_threshold, "*", "")
# 
# # Create the forest plot
# ggplot(data, aes(x = Beta, y = reorder(Phenotype, Beta))) +
#   geom_point(aes(color = Significant)) +
#   geom_errorbarh(aes(xmin = Beta - 1.96 * SE, xmax = Beta + 1.96 * SE)) +
#   scale_color_manual(values = c("black", "red")) +
#   labs(x = "Effect Size (β)", y = "Phenotype") +
#   theme_classic()
# 




# data <- data.frame(
#   Phenotype = c("LV ejection fraction", "LV end-diastolic volume", "LV end-systolic volume", 
#                 "LV stroke volume", "LV mass index", "LV mass", "LA systolic volume", 
#                 "LA volume index", "TAPSE in mm", "LV global longitudinal peak strain", 
#                 "LV relative wall thickness", "3D LV cardiac output", 
#                 "LV posterior wall thickness", "LA systolic diameter", 
#                 "LV diastolic diameter", "Interventricular septal thickness"),
#   Beta = c(0.121, -0.118, -0.149, -0.057, -0.023, -0.017, 0.011, 0.033, 0.043, 
#            -0.126, 0.053, 0.018, 0.068, 0.066, 0.011, 0.037),
#   SE = c(0.047, 0.046, 0.046, 0.046, 0.039, 0.039, 0.044, 0.047, 0.050, 0.043, 
#          0.041, 0.048, 0.042, 0.042, 0.040, 0.041),
#   P_Value = c(0.042, 0.042, 0.020, 0.432, 0.753, 0.806, 0.806, 0.701, 0.617, 0.029,
#               0.432, 0.806, 0.328, 0.328, 0.806, 0.617),
#   Odds_Ratio = c(1.234, 0.876, 0.654, 0.987, 1.123, 0.987, 1.345, 0.987, 0.876, 0.765,
#                  1.234, 0.987, 1.123, 0.987, 1.234, 0.876),
#   Significant = ifelse(data$P_Value < 0.05, "*", "")
# )

# # Create the forest plot
# ggplot(data, aes(x = Beta, y = reorder(Phenotype, Beta))) +
#   geom_point(aes(color = Significant)) +
#   geom_errorbarh(aes(xmin = Beta - 1.96 * SE, xmax = Beta + 1.96 * SE)) +
#   geom_vline(xintercept = 0, linetype = "dashed", color = "gray") +
#   labs(x = "Effect Size (β)", y = "Echo measures") +
#   scale_color_manual(values = c("black", "red")) +
#   geom_text(aes(label = sprintf("p = %.3f\nOR = %.3f", P_Value, Odds_Ratio)),
#             hjust = -0.1, vjust = 0.5, size = 3) +
#   theme_classic()
