df <- read.table(text = "rsID	Gene	EA	NEA	OR CI_low CI_high	P
rs6723526	TTN	C	T	1.36 1.11 1.67	0.0037
rs2303838	TTN	T	C	0.78 0.64 0.95	0.0120
rs3731746	TTN	A	G	0.78 0.64 0.95	0.0138
rs2042996	TTN	A	G	0.81 0.67 0.96	0.0166
rs9808377	TTN	G	A	0.81 0.68 0.96	0.0176
rs1001238	TTN	C	T	0.81 0.68 0.96	0.0180
rs744426	TTN	A	G	0.78 0.63 0.96	0.0202
rs3829747	TTN	T	C	0.78 0.64 0.96	0.0216
rs3829746	TTN	C	T	0.81 0.68 0.97	0.0222
rs3731749	TTN	T	C	0.78 0.64 0.96	0.0222
rs16866406	TTN	A	G	0.78 0.64 0.97	0.0226
rs2234962	BAG3	C	T	0.81 0.68 0.97	0.0241
rs2288569	TTN	T	C	0.8 0.65 0.98	0.0312
rs2042995	TTN	C	T	0.83 0.69 0.98	0.0320
", header = T)


library(gridExtra)

custom_face <- ifelse(df$rsID %in% c("rs3829746", "rs2234962"), "bold", "plain")

p <- ggplot(data = df, aes(x = OR, y = reorder(rsID, -P))) +
  geom_point(size = 3) +
  geom_errorbarh(aes(xmin = CI_low, xmax = CI_high), height = 0) +
  geom_vline(xintercept = 1, linetype = "dashed") +
  scale_x_continuous(limits = c(0.6, 2), breaks = seq(0.6, 2, 0.2)) +
  xlab("OR (95% CI)") +
  ylab("") +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    # panel.border = element_blank(),
    panel.background = element_blank(),
    axis.text.x = element_text(color = "black"),
    axis.text.y = element_text(face = rev(custom_face), color = "black")
  )

p
ggsave("Z:/ResearchHome/ClusterHome/aneupane/St_Jude/Papers/TTN_BAG3/Final_v8/common_variant_analysis.tiff", P, dpi = 600, width = 5, height = 3, units = "in")





echo.sig <- read.table(text = "Phenotype	β	SE	P	n	P_BH	Gender
LV_ejection_fraction	0.121	0.047	0.010	1232	0.042	Both
LV_end-diastolic	-0.118	0.046	0.010	1232	0.042	Both
LV_end-systolic	-0.149	0.046	0.001	1232	0.020	Both
LV_global_longitudinal_peak_strain_(GLPS)	-0.126	0.043	0.004	1359	0.029	Both
LV_ejection_fraction	0.169	0.062	0.006	677	0.025	Male
LV_end-diastolic	-0.172	0.061	0.005	677	0.025	Male
LV_end-systolic	-0.202	0.061	0.001	676	0.017	Male
LV_global_longitudinal_peak_strain_(GLPS)	-0.165	0.059	0.005	750	0.025	Male", header = T)

echo.sig$Phenotype <- gsub("_|-", " ", echo.sig$Phenotype)

# Create the forest plot
library(ggplot2)

## With beta
p <- ggplot(echo.sig, aes(x = Phenotype, y = β, color = Gender)) +
  geom_linerange(aes(ymin = β - 1.96 * SE, ymax = β + 1.96 * SE), position = position_dodge(width = 0.8), size = 1) +
  geom_point(aes(x = as.numeric(factor(Phenotype)), y = β), position = position_dodge(width = 0.8), size = 2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
  labs(x = NULL, y = "Effect Size (β)") +
  geom_text(
    aes(x = as.numeric(factor(Phenotype)), 
        label = sprintf("P = %.3f\nβ = %.2f", P_BH, β), 
        hjust = ifelse(Gender == "Male", 1.2, -0.2),
        y = (β - 1.96 * SE + β + 1.96 * SE) / 2), # Center text vertically
    color = "black", # Set text color to black
    size = 4, show.legend = FALSE) + # Increase text size to 5
  theme_classic() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1, size = 12, color = "black"), # Increase x-axis text size to 14
        axis.text.y = element_text(size = 12, color = "black"), # Increase y-axis text size to 14
        axis.title.y = element_text(size = 12, color = "black"), # Increase Y-axis title size to 16
        legend.text = element_text(size = 12, color = "black"), # Increase legend text size to 14
        legend.title = element_text(size = 12, color = "black"), # Increase legend title size to 14
        plot.title = element_text(size = 12, color = "black")) + # Increase plot title size to 16
  theme(legend.position = "right") +
  scale_x_discrete(labels = echo.sig$Phenotype)


# ## With OR
# # echo.sig$Significant = ifelse(echo.sig$P_BH < 0.05, "*", "")
# echo.sig$OR <- exp(echo.sig$β)
# echo.sig$CI_low <- exp(echo.sig$β - 1.96 * echo.sig$SE)
# echo.sig$CI_high <- exp(echo.sig$β + 1.96 * echo.sig$SE)
# 
# 
# p <- ggplot(echo.sig, aes(x = Phenotype, y = OR, color = Gender)) +
#   geom_linerange(aes(ymin = CI_low, ymax = CI_high), position = position_dodge(width = 0.8), size = 1) +
#   geom_point(aes(x = as.numeric(factor(Phenotype)), y = OR), position = position_dodge(width = 0.8), size = 2) +
#   geom_hline(yintercept = 1, linetype = "dashed", color = "gray") +
#   labs(x = NULL, y = "OR (95% CI") +
#   geom_text(
#     aes(x = as.numeric(factor(Phenotype)), 
#         label = sprintf("P = %.3f\nOR = %.2f", P_BH, OR), 
#         hjust = ifelse(Gender == "Male", 1.2, -0.2),
#         y = (CI_low + CI_high) / 2), # Center text vertically
#     color = "black", # Set text color to black
#     size = 3, show.legend = FALSE) +
#   theme_classic() +
#   theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
#   theme(legend.position = "right") +
#   scale_x_discrete(labels = echo.sig$Phenotype)
#   # scale_color_manual(values = c("Male" = "blue", "Both" = "black"))

p
ggsave("Z:/ResearchHome/ClusterHome/aneupane/St_Jude/Papers/TTN_BAG3/Final_v8/echo_measure_sig.tiff", p, dpi = 600, width = 9, height = 7, units = "in")

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
