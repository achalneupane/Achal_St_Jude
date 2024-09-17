# Load necessary libraries
library(ggplot2)
library(dplyr)
library(gridExtra)
library(ggpubr)

data_table2 <- data.frame(
  Cohorts = rep(c("SJLIFE EUR", "CCSS EUR", "SJLIFE+CCSS EUR", "SJLIFE AFR"), each = 2),
  Variant = rep(c("TTN (rs3829746)", "BAG3 (rs2234962)"), times = 4),
  EAF = c(0.22, 0.21, 0.22, 0.20, 0.22, 0.21, 0.54, 0.03),
  OR = c(0.71, 0.75, 0.88, 0.84, 0.81, 0.81, 1.55, 0.24),
  CI_lower = c(0.54, 0.57, 0.69, 0.66, 0.68, 0.67, 0.82, 0.03),
  CI_upper = c(0.94, 0.98, 1.11, 1.08, 0.97, 0.96, 2.96, 1.74),
  P = c(0.017, 0.038, 0.27, 0.17, 0.020, 0.018, 0.18, 0.16),
  Population = rep(c("EUR", "EUR", "EUR", "AFR"), each = 2)  # New column for populations
)

# Create the first panel for Table 2 results with facets for TTN and BAG3
p0 <- ggplot(data_table2, aes(x = Cohorts, y = OR)) +
  geom_point(aes(color = Variant), size = 4) +
  geom_errorbar(aes(ymin = CI_lower, ymax = CI_upper, color = Variant), width = 0.1, size = 1, position = position_dodge(width = 0.4)) +
  scale_y_log10() +
  labs(title = "Association of TTN and BAG3 with CCM Risk in Long-term Survivors",
       x = "Cohorts",
       y = "Odds Ratio (95% CI)") +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  facet_wrap(~ Variant) + # Facet by genetic variant (TTN and BAG3) 
  scale_color_manual(values = c("TTN (rs3829746)" = "#F8766D", "BAG3 (rs2234962)" = "#00BFC4"))  # Distinct colors for variants
p0




# Data preparation for genetic associations including both TTN and BAG3
genetic_data <- data.frame(
  Variant = rep(c("TTN (rs3829746)", "BAG3 (rs2234962)"), each = 5),
  Stratab = rep(c("Male", "Female", "High risk", "Moderate risk", "Low-risk"), times = 2),
  Total_N = c(3109, 3140, 1651, 1163, 2065,
              3109, 3140, 1651, 1163, 2065),
  N_with_CCM = c(222, 239, 255, 60, 30,
                 222, 239, 255, 60, 30),
  OR = c(0.81, 0.81, 0.89, 0.72, 0.38,
         0.71, 0.89, 0.80, 0.94, 0.87),
  Lower_CI = c(0.63, 0.63, 0.70, 0.42, 0.16,
               0.54, 0.70, 0.63, 0.58, 0.44),
  Upper_CI = c(1.10, 1.03, 1.13, 1.23, 0.91,
               0.92, 1.14, 1.03, 1.51, 1.71),
  P_value = c(0.11, 0.086, 0.33, 0.23, 0.030,
              9.7E-03, 0.37, 0.081, 0.78, 0.68)
)

# Reorder factors for proper plot arrangement
genetic_data$Stratab <- factor(genetic_data$Stratab, levels = c("Male", "Female", "High risk", "Moderate risk", "Low-risk"))

# Create a combined plot for TTN and BAG3 with facets, using a distinct color scheme
p1 <- ggplot(genetic_data, aes(x = Stratab, y = OR, color = Variant)) +
  geom_point(size = 4) +  # Increased size for better visibility
  geom_errorbar(aes(ymin = Lower_CI, ymax = Upper_CI), width = 0.2, size = 1, position = position_dodge(width = 0.4)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  scale_y_log10() +  # Use a logarithmic scale for the y-axis
  labs(title = "Genetic Associations for TTN and BAG3",
       x = "Stratification",
       y = "Odds Ratio (OR)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  facet_wrap(~ Variant) +  # Facet by genetic variant (TTN and BAG3)
  scale_color_manual(values = c("TTN (rs3829746)" = "#F8766D", "BAG3 (rs2234962)" = "#00BFC4"))  # Distinct colors for variants

# Print the combined plot
print(p1)



# Creating the SJLIFE EUR data frame
sjlife_eur <- data.frame(
  Parameter = c("LV ejection fraction", "LV end-diastolic volume", "LV end-systolic volume", 
                "LV stroke volume", "LV mass index", "Global longitudinal peak strain", 
                "LV relative wall thickness"),
  BAG3 = c(0.17, -0.11, -0.17, -0.03, 0.02, -0.09, -0.02),
  P_BAG3 = c(2.9e-3, 0.035, 2.1e-3, 0.61, 0.83, 0.083, 0.69),
  TTN = c(0.05, -0.08, -0.10, -0.05, -0.03, -0.15, 0.02),
  P_TTN = c(0.39, 0.11, 0.074, 0.39, 0.65, 4.3e-3, 0.66),
  Group = "EUR"
)

# Creating the SJLIFE AFR data frame
sjlife_afr <- data.frame(
  Parameter = c("LV ejection fraction", "LV end-diastolic volume", "LV end-systolic volume", 
                "LV stroke volume", "LV mass index", "Global longitudinal peak strain", 
                "LV relative wall thickness"),
  BAG3 = c(0.18, 0.08, 0.01, 0.27, -0.43, -0.30, 0.09),
  P_BAG3 = c(0.58, 0.79, 0.98, 0.37, 0.50, 0.29, 0.74),
  TTN = c(0.17, -0.02, 0.00, 0.07, 0.002, 0.08, -0.12),
  P_TTN = c(0.15, 0.86, 0.99, 0.55, 0.99, 0.47, 0.22),
  Group = "AFR"
)

# Combine both datasets
combined_data <- bind_rows(sjlife_eur, sjlife_afr)

library(dplyr)
library(tidyr)
library(ggrepel)

# Convert to long format for plotting
long_data <- combined_data %>%
  pivot_longer(cols = c(BAG3, TTN), names_to = "Gene", values_to = "Effect_Size") %>%
  mutate(P_Value = ifelse(Gene == "BAG3", P_BAG3, P_TTN),
         Significance = paste0("p=", round(P_Value, 3)),  # Show p-value for all
         Parameter = factor(Parameter, levels = c("LV ejection fraction", "LV end-diastolic volume",
                                                  "LV end-systolic volume", "LV stroke volume",
                                                  "LV mass index", "Global longitudinal peak strain",
                                                  "LV relative wall thickness")),
         Label_Pos = ifelse(Effect_Size >= 0, Effect_Size + 0.02, Effect_Size - 0.02),  # Adjust based on sign of Effect_Size
         vjust = ifelse(Effect_Size >= 0, 0, 1))  # Calculate vertical adjustment

# Rename genes for clarity
long_data$Gene[long_data$Gene == "BAG3"] <- "BAG3 (rs2234962)"
long_data$Gene[long_data$Gene == "TTN"] <- "TTN (rs3829746)"

# long_data <- long_data %>%
#   filter(P_Value < 0.05)  # Adjust this threshold as needed

# Plotting the data
p2 <- ggplot(long_data, aes(x = Parameter, y = Effect_Size, fill = Gene)) +
  geom_bar(stat = "identity", position = "dodge", alpha = 0.7) +
  labs(title = "Echocardiographic Parameters by Gene",
       x = "Parameters",
       y = "Effect Size (β)",
       fill = "Gene") +  # Legend title
  scale_fill_manual(values = c("TTN (rs3829746)" = "#F8766D", "BAG3 (rs2234962)" = "#00BFC4")) +  # Custom colors
  theme_minimal() +
  geom_text(aes(label = Significance, y = Label_Pos, vjust = vjust), 
            position = position_dodge(width = 0.9),  # Correctly position the labels
            size = 2.5, 
            color = "black") + 
  facet_wrap(~ Group) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels

p2

# Arrange plots in a grid
combined_plot <- ggarrange(p0, p1, p2,
                           ncol = 1, nrow = 4,
                           common.legend = TRUE, legend = "top")

combined_plot
# Add annotation
annotate_figure(combined_plot,
                top = text_grob("Central Illustration for JACC Publication", face = "bold", size = 14))

# Save the plot
ggsave("Z:/ResearchHome/ClusterHome/aneupane/St_Jude/Papers/central_illustration_separate_panels.png", plot = combined_plot, width = 10, height = 15, dpi = 300)














# sjlife_eur <- data.frame(
#   Parameter = c("LV ejection fraction", "LV end-diastolic volume", "LV end-systolic volume", 
#                 "LV stroke volume", "LV mass index", "Global longitudinal peak strain", 
#                 "LV relative wall thickness"),
#   BAG3 = c(0.17, -0.11, -0.17, -0.03, 0.02, -0.09, -0.02),
#   P_BAG3 = c(2.9e-3, 0.035, 2.1e-3, 0.61, 0.83, 0.083, 0.69),
#   TTN = c(0.05, -0.08, -0.10, -0.05, -0.03, -0.15, 0.02),
#   P_TTN = c(0.39, 0.11, 0.074, 0.39, 0.65, 4.3e-3, 0.66),
#   Group = "EUR"
# )
# 
# # Creating the SJLIFE AFR data frame
# sjlife_afr <- data.frame(
#   Parameter = c("LV ejection fraction", "LV end-diastolic volume", "LV end-systolic volume", 
#                 "LV stroke volume", "LV mass index", "Global longitudinal peak strain", 
#                 "LV relative wall thickness"),
#   BAG3 = c(0.18, 0.08, 0.01, 0.27, -0.43, -0.30, 0.09),
#   P_BAG3 = c(0.58, 0.79, 0.98, 0.37, 0.50, 0.29, 0.74),
#   TTN = c(0.17, -0.02, 0.00, 0.07, 0.002, 0.08, -0.12),
#   P_TTN = c(0.15, 0.86, 0.99, 0.55, 0.99, 0.47, 0.22),
#   Group = "AFR"
# )
# 
# # Combine both datasets
# combined_data <- bind_rows(sjlife_eur, sjlife_afr)
# 
# # Convert to long format for plotting without the names_sep
# long_data <- combined_data %>%
#   pivot_longer(cols = c(BAG3, TTN), names_to = "Gene", values_to = "Effect_Size") %>%
#   mutate(P_Value = ifelse(Gene == "BAG3", P_BAG3, P_TTN),
#          Significance = ifelse(P_Value < 0.05, paste0("p=", round(P_Value, 3)), ""),
#          Parameter = factor(Parameter, levels = c("LV ejection fraction", "LV end-diastolic volume", 
#                                                   "LV end-systolic volume", "LV stroke volume", 
#                                                   "LV mass index", "Global longitudinal peak strain", 
#                                                   "LV relative wall thickness")))
# 
# # Plotting the data
# ggplot(long_data, aes(x = Parameter, y = Effect_Size, fill = Gene)) +
#   geom_bar(stat = "identity", position = "dodge", alpha = 0.7) +
#   labs(title = "Echocardiographic Parameters by Gene",
#        x = "Parameters",
#        y = "Effect Size (β)",
#        fill = "Gene") +  # Legend title
#   theme_minimal() +
#   geom_text(aes(label = Significance), position = position_dodge(width = 0.9), vjust = -0.5, size = 4, color = "black") +
#   facet_wrap(~ Group)
p2
