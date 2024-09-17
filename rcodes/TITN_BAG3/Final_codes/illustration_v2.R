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

data_table2$Cohorts <- factor(data_table2$Cohorts, levels = c("SJLIFE EUR", "CCSS EUR", "SJLIFE+CCSS EUR", "SJLIFE AFR"))
# "Association of TTN and BAG3 with CCM Risk in Long-term Survivors"
# Create the first panel for Table 2 results with facets for TTN and BAG3
p0 <- ggplot(data_table2, aes(x = Cohorts, y = OR)) +
  geom_point(aes(color = Variant), size = 4) +
  geom_errorbar(aes(ymin = CI_lower, ymax = CI_upper, color = Variant), width = 0.1, size = 1, position = position_dodge(width = 0.4)) +
  # scale_y_log10() +
  labs(title = "",
       x = "Cohorts",
       y = "Odds Ratio (95% CI)") +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  facet_wrap(~ Variant) + # Facet by genetic variant (TTN and BAG3) 
  scale_color_manual(values = c("TTN (rs3829746)" = "#F8766D", "BAG3 (rs2234962)" = "#00BFC4")) +
  theme_bw() +  # White background with box
  theme(axis.text.x = element_text(angle = 60, hjust = 1, color = "black", size = 10),  # Set x-axis text size
        axis.text.y = element_text(color = "black", size = 10),  # Set y-axis text size
        axis.title.x = element_text(color = "black", size = 12),  # Set x-axis title size
        axis.title.y = element_text(color = "black", size = 12),  # Set y-axis title size
        legend.text = element_text(size = 12),  # Set legend text size
        legend.title = element_text(size = 14),  # Set legend title size
        # strip.text = element_text(color = "black"),  # Set facet labels color to black
        strip.text = element_blank(),  # Remove facet labels
        plot.title = element_text(color = "black"),  # Set plot title color to black
        panel.grid.major = element_blank(),  # Remove major grid lines
        legend.position = "top",  # Adjust legend position
        panel.grid.minor = element_blank(),
        text = element_text(color = "black")) # Set all other text to black
p0

ggsave("Z:/ResearchHome/ClusterHome/aneupane/St_Jude/Papers/TTN_BAG3/Final_v9//cohort_associations.tiff", plot = p0, width = 5, height = 4, dpi = 400)


# Data preparation for genetic associations including both TTN and BAG3
genetic_data <- data.frame(
  Variant = rep(c("TTN (rs3829746)", "BAG3 (rs2234962)"), each = 5),
  Stratab = rep(c("Male", "Female", "High-risk", "Moderate-risk", "Low-risk"), times = 2),
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
genetic_data$Stratab <- factor(genetic_data$Stratab, levels = c("Male", "Female", "High-risk", "Moderate-risk", "Low-risk"))

# Create a combined plot for TTN and BAG3 with facets, using a distinct color scheme
# "Genetic Associations for TTN and BAG3"
p1 <- ggplot(genetic_data, aes(x = Stratab, y = OR, color = Variant)) +
  geom_point(size = 4) +  # Increased size for better visibility
  geom_errorbar(aes(ymin = Lower_CI, ymax = Upper_CI), width = 0.2, size = 1, position = position_dodge(width = 0.4)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  # scale_y_log10() +  # Use a logarithmic scale for the y-axis
  labs(title = "",
       x = "Stratification",
       y = "Odds Ratio (95% CI)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  facet_wrap(~ Variant) +  # Facet by genetic variant (TTN and BAG3)
  scale_color_manual(values = c("TTN (rs3829746)" = "#F8766D", "BAG3 (rs2234962)" = "#00BFC4")) + # Distinct colors for variants
  theme_bw() +  # White background with box
  theme(axis.text.x = element_text(angle = 60, hjust = 1, color = "black", size = 10),  # Set x-axis text size
        axis.text.y = element_text(color = "black", size = 10),  # Set y-axis text size
        axis.title.x = element_text(color = "black", size = 12),  # Set x-axis title size
        axis.title.y = element_text(color = "black", size = 12),  # Set y-axis title size
        legend.text = element_text(size = 12),  # Set legend text size
        legend.title = element_text(size = 14),  # Set legend title size
        # strip.text = element_text(color = "black"),  # Set facet labels color to black
        strip.text = element_blank(),  # Remove facet labels
        plot.title = element_text(color = "black"),  # Set plot title color to black
        panel.grid.major = element_blank(),  # Remove major grid lines
        legend.position = "top",  # Adjust legend position
        panel.grid.minor = element_blank(),
        text = element_text(color = "black")) # Set all other text to black
# Print the combined plot
print(p1)

ggsave("Z:/ResearchHome/ClusterHome/aneupane/St_Jude/Papers/TTN_BAG3/Final_v9//stratified_associations.tiff", plot = p1, width = 5, height = 4, dpi = 400)


# Creating the SJLIFE EUR data frame

# data <- read.table(text="Survivors	Echocardiographic_parameters	B_BAG3_EUR	SE_BAG3_EUR	P_BAG3_EUR	B_TTN_EUR	SE_TTN_EUR	P_TTN_EUR	B_BAG3_AFR	SE_BAG3_AFR	P_BAG3_AFR	B_TTN_AFR	SE_TTN_AFR	P_TTN_AFR
# All	LV ejection fraction	0.17	0.06	0.0029	0.05	0.06	0.39	0.18	0.31	0.58	0.17	0.12	0.15
# All	LV end-diastolic volume	-0.11	0.05	0.035	-0.08	0.05	0.11	0.08	0.30	0.79	-0.02	0.12	0.86
# All	LV end-systolic volume	-0.17	0.05	0.0021	-0.10	0.06	0.074	0.01	0.31	0.98	0.00	0.12	0.99
# All	LV stroke volume	-0.03	0.05	0.61	-0.05	0.05	0.39	0.27	0.29	0.37	0.07	0.11	0.55
# All	LV mass index	0.02	0.07	0.83	-0.03	0.07	0.65	-0.43	0.63	0.5	0.002	0.19	0.99
# All	Global longitudinal peak strain	-0.09	0.05	0.083	-0.15	0.05	0.0043	-0.30	0.28	0.29	0.08	0.10	0.47
# All	LV relative wall thickness	-0.02	0.04	0.69	0.02	0.05	0.66	0.09	0.27	0.74	-0.12	0.10	0.22
# Male	LV ejection fraction	0.12	0.08	0.11	0.05	0.08	0.53	0.32	0.43	0.45	0.23	0.18	0.22
# Male	LV end-diastolic volume	-0.08	0.07	0.27	-0.11	0.08	0.19	-0.47	0.43	0.28	-0.07	0.19	0.7
# Male	LV end-systolic volume	-0.12	0.08	0.11	-0.13	0.08	0.13	-0.49	0.46	0.29	0.01	0.20	0.96
# Male	LV stroke volume	0.00	0.07	0.99	-0.05	0.08	0.55	-0.32	0.43	0.46	0.08	0.18	0.66
# Male	LV mass index	-0.27	0.11	0.013	0.07	0.11	0.54	0.00	0.05	0.99	-0.11	0.41	0.79
# Male	Global longitudinal peak strain	-0.05	0.07	0.46	-0.11	0.07	0.13	-0.37	0.39	0.34	0.01	0.14	0.96
# Male	LV relative wall thickness	-0.01	0.06	0.86	0.03	0.06	0.6	0.74	0.38	0.051	-0.08	0.15	0.6
# Female	LV ejection fraction	0.23	0.09	0.01	0.06	0.09	0.53	-0.11	0.51	0.83	0.09	0.18	0.64
# Female	LV end-diastolic volume	-0.19	0.09	0.043	-0.05	0.09	0.56	1.04	0.49	0.038	-0.03	0.19	0.89
# Female	LV end-systolic volume	-0.28	0.09	0.0036	-0.07	0.09	0.45	0.78	0.50	0.13	-0.03	0.19	0.89
# Female	LV stroke volume	-0.08	0.09	0.37	-0.03	0.09	0.74	1.04	0.49	0.04	-0.02	0.19	0.92
# Female	LV mass index	0.10	0.10	0.29	-0.07	0.10	0.51	-0.37	0.75	0.63	-0.16	0.25	0.53
# Female	Global longitudinal peak strain	-0.11	0.07	0.13	-0.21	0.08	0.0076	0.14	0.45	0.76	0.07	0.16	0.65
# Female	LV relative wall thickness	-0.04	0.07	0.57	0.01	0.07	0.94	-0.85	0.43	0.051	-0.13	0.14	0.37", header = T, sep = "\t")

data <- read.table(text="Survivors	Echocardiographic_parameters	B_BAG3_EUR	SE_BAG3_EUR	P_BAG3_EUR	B_TTN_EUR	SE_TTN_EUR	P_TTN_EUR	B_BAG3_AFR	SE_BAG3_AFR	P_BAG3_AFR	B_TTN_AFR	SE_TTN_AFR	P_TTN_AFR
All	LV ejection fraction	0.17	0.06	0.0029	0.05	0.06	0.39	0.18	0.31	0.58	0.17	0.12	0.15
All	LV end-diastolic volume	-0.11	0.05	0.035	-0.08	0.05	0.11	0.08	0.30	0.79	-0.02	0.12	0.86
All	LV end-systolic volume	-0.17	0.05	0.0021	-0.10	0.06	0.074	0.01	0.31	0.98	0.00	0.12	0.99
All	LV stroke volume	-0.03	0.05	0.61	-0.05	0.05	0.39	0.27	0.29	0.37	0.07	0.11	0.55
All	LV mass index	0.02	0.07	0.83	-0.03	0.07	0.65	-0.43	0.63	0.5	0.002	0.19	0.99
All	Global longitudinal peak strain	-0.09	0.05	0.083	-0.15	0.05	0.0043	-0.30	0.28	0.29	0.08	0.10	0.47
All	LV relative wall thickness	-0.02	0.04	0.69	0.02	0.05	0.66	0.09	0.27	0.74	-0.12	0.10	0.22", header = T, sep = "\t")


data <- data %>%
  mutate(OddsRatio_BAG3_EUR = exp(B_BAG3_EUR),
         UpperCI_BAG3_EUR = exp(B_BAG3_EUR + 1.96 * SE_BAG3_EUR),
         LowerCI_BAG3_EUR = exp(B_BAG3_EUR - 1.96 * SE_BAG3_EUR),
         OddsRatio_TTN_EUR = exp(B_TTN_EUR),
         UpperCI_TTN_EUR = exp(B_TTN_EUR + 1.96 * SE_TTN_EUR),
         LowerCI_TTN_EUR = exp(B_TTN_EUR - 1.96 * SE_TTN_EUR),
         OddsRatio_BAG3_AFR = exp(B_BAG3_AFR),
         UpperCI_BAG3_AFR = exp(B_BAG3_AFR + 1.96 * SE_BAG3_AFR),
         LowerCI_BAG3_AFR = exp(B_BAG3_AFR - 1.96 * SE_BAG3_AFR),
         OddsRatio_TTN_AFR = exp(B_TTN_AFR),
         UpperCI_TTN_AFR = exp(B_TTN_AFR + 1.96 * SE_TTN_AFR),
         LowerCI_TTN_AFR = exp(B_TTN_AFR - 1.96 * SE_TTN_AFR)) 

data <- data [!grepl("_AFR", colnames(data))]

data_long <- data %>%
  pivot_longer(
    cols = c(starts_with("OddsRatio"), starts_with("UpperCI"), starts_with("LowerCI"), starts_with("P_")),
    names_to = c("Metric", "Variant", "Population"),
    names_sep = "_",
    values_to = "Value"
  ) %>%
  pivot_wider(
    names_from = Metric,
    values_from = Value
  ) %>%
  filter(Population == "EUR")  # Filter for EUR population


data_long$Variant[data_long$Variant == "TTN"] <- "TTN (rs3829746)"
data_long$Variant[data_long$Variant == "BAG3"] <- "BAG3 (rs2234962)"

# data_long <- data_long[data_long$P < 0.05,]

library("stringr")
# "Odds Ratios with 95% Confidence Intervals"
p2 <- ggplot(data_long, aes(x = Echocardiographic_parameters, y = OddsRatio, color = Variant)) +
  geom_point(size = 4, alpha = 0.7) +
  geom_errorbar(aes(ymin = LowerCI, ymax = UpperCI), width = 0.2) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  facet_grid(Variant ~ Survivors, scales = "free_y") +
  scale_color_manual(values = c("TTN (rs3829746)" = "#F8766D", "BAG3 (rs2234962)" = "#00BFC4")) + # Distinct colors for variants
  labs(title = "",
       y = "Odds Ratio (95% CI)",
       x = "Echocardiographic parameters") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1, color = "black", size = 10),  # Set x-axis text size
        axis.text.y = element_text(color = "black", size = 10),  # Set y-axis text size
        axis.title.x = element_text(color = "black", size = 12),  # Set x-axis title size
        axis.title.y = element_text(color = "black", size = 12),  # Set y-axis title size
        legend.text = element_text(size = 12),  # Set legend text size
        legend.title = element_text(size = 14),  # Set legend title size
        strip.text = element_text(color = "black"),
        plot.title = element_text(color = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "top",
        legend.box.margin = margin(0, 10, 0, 0),
        plot.margin = margin(t = 10, r = 30, b = 10, l = 10)) +
  scale_x_discrete(labels = function(x) sapply(x, function(y) str_wrap(y, width = 20))) +
  guides(color = guide_legend(override.aes = list(size = 4)))



p2

ggsave("Z:/ResearchHome/ClusterHome/aneupane/St_Jude/Papers/TTN_BAG3/Final_v9//echo_associations.tiff", plot = p2, width = 7, height = 5, dpi = 600)
# Arrange plots in a grid
combined_plot <- ggarrange(p0, p1, p2,
                           ncol = 1, nrow = 4,
                           common.legend = T, legend = "top")

combined_plot
# Add annotation
annotate_figure(combined_plot,
                top = text_grob("", face = "bold", size = 14))

# Save the plot
ggsave("Z:/ResearchHome/ClusterHome/aneupane/St_Jude/Papers/TTN_BAG3/Final_v9//central_illustration_separate_panels.tiff", plot = combined_plot, width = 10, height = 15, dpi = 600)






library(gt)
library(ggplot2)

# Create the data frame
strata_data <- data.frame(
  Stratification = c("Male", "Female", "High risk", "Moderate risk", "Low risk"),
  Total_N = c(3109, 3140, 1651, 1163, 2065),
  CCM_N = c(222, 239, 255, 60, 30)
)

strata_data$Stratification <- factor(strata_data$Stratification, 
                                     levels = c("Male", "Female", "High risk", "Moderate risk", "Low risk"))

# Reshape the data for plotting
library(reshape2)
strata_long <- melt(strata_data, id.vars = "Stratification", 
                    variable.name = "Category", value.name = "Count")

# Create a bar plot
p3 <- ggplot(strata_long, aes(x = Stratification, y = Count, fill = Category)) +
  geom_bar(stat = "identity", position = "dodge", alpha = 0.7) +
  geom_text(aes(label = Count), vjust = -0.5, size = 4, color = "black", 
            position = position_dodge(0.9)) +
  labs(title = "", 
       x = "", 
       y = "Counts") +
  scale_fill_manual(values = c("Total_N" = "#00BFC4", "CCM_N" = "#F8766D"),
                    labels = c("Total Participants", "CCM Cases")) +  # Custom legend labels
  scale_y_continuous(limits = c(0, 4000)) +  # Set y-axis limits
  theme_bw() +  # White background with box
  theme(axis.text.x = element_text(angle = 30, hjust = 1, color = "black", size = 10),  # Set x-axis text size
        axis.text.y = element_text(color = "black", size = 10),  # Set y-axis text size
        axis.title.x = element_text(color = "black", size = 12),  # Set x-axis title size
        axis.title.y = element_text(color = "black", size = 12),  # Set y-axis title size
        legend.text = element_text(size = 12),  # Set legend text size
        legend.title = element_text(size = 14),  # Set legend title size
        strip.text = element_blank(),  # Remove facet labels
        plot.title = element_text(color = "black"),  # Set plot title color to black
        panel.grid.major = element_blank(),  # Remove major grid lines
        legend.position = "top",  # Adjust legend position
        panel.grid.minor = element_blank(),
        text = element_text(color = "black"))  # Set all other text to black
p3

ggsave("Z:/ResearchHome/ClusterHome/aneupane/St_Jude/Papers/TTN_BAG3/Final_v9//stratified_counts.tiff", plot = p3, width = 5, height = 3, dpi = 600)



# Create the data frame
cohort_data <- data.frame(
  Cohorts = c("SJLIFE EUR", "CCSS EUR", "SJLIFE+CCSS EUR", "SJLIFE AFR"),
  Total_N = c(1645, 4604, 6249, 246),
  CCM_N = c(211, 250, 461, 40)
)

cohort_data$Cohorts <- factor(cohort_data$Cohorts, 
                              levels = c("SJLIFE EUR", "CCSS EUR", "SJLIFE+CCSS EUR", "SJLIFE AFR"))

# Reshape the data for plotting
cohort_long <- melt(cohort_data, id.vars = "Cohorts", 
                    variable.name = "Category", value.name = "Count")

# Create a bar plot
p4 <- ggplot(cohort_long, aes(x = Cohorts, y = Count, fill = Category)) +
  geom_bar(stat = "identity", position = "dodge", alpha = 0.7) +
  geom_text(aes(label = Count), vjust = -0.5, size = 4, color = "black", 
            position = position_dodge(0.9)) +
  labs(title = "", 
       x = "", 
       y = "Count") +
  scale_fill_manual(values = c("Total_N" = "#00BFC4", "CCM_N" = "#F8766D"),
                    labels = c("Total Participants", "CCM Cases")) +  # Custom legend labels
  scale_y_continuous(limits = c(0, 8000)) +  # Set y-axis limits
  theme_bw() +  # White background with box
  theme(axis.text.x = element_text(angle = 30, hjust = 1, color = "black", size = 10),  # Set x-axis text size
        axis.text.y = element_text(color = "black", size = 10),  # Set y-axis text size
        axis.title.x = element_text(color = "black", size = 12),  # Set x-axis title size
        axis.title.y = element_text(color = "black", size = 12),  # Set y-axis title size
        legend.text = element_text(size = 12),  # Set legend text size
        legend.title = element_text(size = 14),  # Set legend title size
        plot.title = element_text(color = "black", size = 18, hjust = 0.5),  # Set plot title size
        strip.text = element_blank(),  # Remove facet labels
        panel.grid.major = element_blank(),  # Remove major grid lines
        panel.grid.minor = element_blank(),
        legend.position = "top",  # Adjust legend position
        text = element_text(color = "black"))  # Set all other text to black
p4
ggsave("Z:/ResearchHome/ClusterHome/aneupane/St_Jude/Papers/TTN_BAG3/Final_v9//cohort_counts.tiff", plot = p4, width = 5, height = 3, dpi = 600)
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
