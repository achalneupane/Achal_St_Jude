library(ggplot2)

# Define colors for each SN type
sn_colors <- c("#1f78b4", "#33a02c", "#e31a1c", "#ff7f00", "#6a3d9a", "#a6cee3", "#b2df8a") 

# Define darker colors for SJLIFE and CCSS
sjlife_color <- "#1f78b4"
ccss_color <- "#33a02c"
darker_sjlife <- adjustcolor(sn_colors, alpha.f = 0.6)
darker_ccss <- adjustcolor(sn_colors, alpha.f = 0.2)

# Define data
data_melted <- data.frame(
  Cohort = c("SJLIFE", "SJLIFE", "SJLIFE", "SJLIFE", "SJLIFE", "SJLIFE", "SJLIFE", 
             "CCSS", "CCSS", "CCSS", "CCSS", "CCSS", "CCSS", "CCSS"),
  SN_types = c("Any SN", "SMN", "NMSC", "Breast cancer", "Thyroid cancer", "Meningioma", "Sarcoma", 
               "Any SN", "SMN", "NMSC", "Breast cancer", "Thyroid cancer", "Meningioma", "Sarcoma"),
  Variables = rep("Combined", 14),
  Overall = c(0.536, 0.467, 0.688, 0.704, 0.866, 0.423, 0.316, 
              0.435, 0.323, 0.577, 0.743, 0.664, 0.444, 0.362)
)





# Reorder data so SJLIFE bars come first for each SN type
data_melted <- data_melted[order(data_melted$Cohort, decreasing = TRUE), ]



# Define colors for each SN type
sn_colors <- c("#1f78b4", "#33a02c", "#e31a1c", "#ff7f00", "#6a3d9a", "#a6cee3", "#b2df8a") 

# Define darker colors for SJLIFE and CCSS
sjlife_color <- "#1f78b4"
ccss_color <- "#33a02c"
darker_sjlife <- adjustcolor(sn_colors, alpha.f = 0.6)
darker_ccss <- adjustcolor(sn_colors, alpha.f = 0.2)

angle_increment.AF <- 360 / 14
angle_increment.SN_type <- 360 / 7


# Reorder levels of SN_types based on Cohort
data_melted$SN_types <- factor(data_melted$SN_types, levels = unique(paste0(data_melted$SN_types)))


# Plot
ggplot(data_melted, aes(x = SN_types, y = Overall, fill = Cohort, label = Overall)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
  # geom_text(aes(y = Overall , angle = 0), position = position_dodge(width = 1), vjust = 0) +
  # geom_vline(xintercept = seq(1.5, 19.5, by = 3),color = "black", linetype = "dashed") + # If you want to include a vertical line
  geom_hline(yintercept = 1) +
  scale_fill_manual(values = c(darker_sjlife, darker_ccss), labels = c("SJLIFE", "CCSS")) +  # Specify colors for SJLIFE and CCSS bars
  coord_polar() +
  theme_minimal() +
  labs(title = "",
       x = NULL, y = "Attributable fraction",
       fill = "Cohort")

ggplot(data_melted, aes(x = SN_types, y = Overall, fill = Cohort, label = Overall)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
  geom_text(aes(y = Overall , angle = 0), position = position_dodge(width = 1), vjust = 0) +
  geom_vline(xintercept = seq(0.5, length(unique(data_melted$SN_types)) + 0.5, by = 1), color = "black", linetype = "solid") + 
  geom_hline(yintercept = 1) +
  scale_fill_manual(values = c(darker_sjlife, darker_ccss), labels = c("SJLIFE", "CCSS")) + 
  coord_polar() +
  theme_minimal() +
  labs(title = "Overall Score by SN Types and Cohort",
       x = NULL, y = "Overall Score",
       fill = "Cohort")

#####################################




# Plot
ggplot(data_melted, aes(x = SN_types, y = Overall, fill = Cohort, label = Overall)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
  geom_vline(xintercept = seq(1.5, 19.5, by = 3), color = "black", linetype = "dashed") + # If you want to include a vertical line
  geom_hline(yintercept = 1) +
  scale_fill_manual(values = c(darker_sjlife, darker_ccss), labels = c("SJLIFE", "CCSS")) +  # Specify colors for SJLIFE and CCSS bars
  scale_y_continuous(breaks = seq(0, 1, by = 0.25)) +  # Custom grid lines from 0 to 1 by 0.25
  coord_polar() +
  theme_minimal() +
  theme(panel.grid.major.y = element_line(color = "gray", linetype = "dotted"),  # Adjust major grid lines
        panel.grid.minor = element_blank(),  # Remove minor grid lines
        axis.line = element_blank(),   # Remove axis lines
        axis.text.x = element_text(hjust = 0.5, angle = angle_increment),  # Center SN type labels within grid
        axis.text = element_text(size = 12),  # Adjust axis text size
        axis.title = element_text(size = 14),  # Adjust axis title size
        legend.title = element_text(size = 14),  # Adjust legend title size
        legend.text = element_text(size = 12)) +  # Adjust legend text size
  labs(title = "",
       x = NULL, y = "Attributable fraction",
       fill = "Cohort")


############################################
# Plot
# Calculate angle increment
angle_increment.AF <- 360 / 14
angle_increment.SN_type <- 360 / 7

# Calculate angle increment
angle_increment <- 360 / length(unique(data_melted$SN_types))

# Plot
ggplot(data_melted, aes(x = SN_types, y = Overall, fill = Cohort, label = Overall)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
  geom_vline(xintercept = seq(1.5, 19.5, by = 3), color = "black", linetype = "dashed") + # If you want to include a vertical line
  geom_hline(yintercept = 1) +
  geom_text(aes(angle = 0 + (angle_increment.AF * (as.numeric(as.factor(SN_types)) - 1))),  # Adjust angle of text labels
            position = position_dodge(width = 1), vjust = 0, size = 2) +  # Adjust label position and size
  scale_fill_manual(values = c(darker_sjlife, darker_ccss), labels = c("SJLIFE", "CCSS")) +  # Specify colors for SJLIFE and CCSS bars
  scale_y_continuous(breaks = seq(0, 1, by = 0.25)) +  # Custom grid lines from 0 to 1 by 0.25
  coord_polar() +
  theme_minimal() +
  theme(panel.grid.major.y = element_line(color = "gray", linetype = "dotted"),  # Adjust major grid lines
        panel.grid.minor = element_blank(),  # Remove minor grid lines
        axis.line = element_blank(),   # Remove axis lines
        axis.text.x = element_text(hjust = 0.5, angle = angle_increment),  # Center SN type labels within grid
        axis.text = element_text(size = 12),  # Adjust axis text size
        axis.title = element_text(size = 14),  # Adjust axis title size
        legend.title = element_text(size = 14),  # Adjust legend title size
        legend.text = element_text(size = 12)) +  # Adjust legend text size
  labs(title = "",
       x = NULL, y = "Attributable fraction",
       fill = "Cohort")


#################

ggplot(data_melted, aes(x = SN_types, y = Overall, fill = Cohort, label = Overall)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
  geom_vline(xintercept = seq(1.5, 19.5, by = 3), color = "black", linetype = "dashed") + # If you want to include a vertical line
  geom_hline(yintercept = 1) +
  geom_text(aes(angle = 0 + (angle_increment.AF * (as.numeric(as.factor(SN_types)) - 1))),  # Adjust angle of text labels
            position = position_dodge(width = 1), vjust = 0, size = 2) +  # Adjust label position and size
  scale_fill_manual(values = c(darker_sjlife, darker_ccss), labels = c("SJLIFE", "CCSS")) +  # Specify colors for SJLIFE and CCSS bars
  scale_y_continuous(breaks = NULL) +  # Remove Y-axis values
  coord_polar() +
  theme_minimal() +
  theme(panel.grid.major.y = element_line(color = "gray", linetype = "dotted"),  # Adjust major grid lines
        panel.grid.minor = element_blank(),  # Remove minor grid lines
        axis.line = element_blank(),   # Remove axis lines
        axis.text.x = element_text(hjust = 0.5, angle = angle_increment),  # Center SN type labels within grid
        axis.text = element_text(size = 12),  # Adjust axis text size
        axis.title = element_text(size = 14),  # Adjust axis title size
        legend.title = element_text(size = 14),  # Adjust legend title size
        legend.text = element_text(size = 12)) +  # Adjust legend text size
  guides(fill = guide_legend(title = NULL)) +  # Remove left legend values for attributable fraction
  labs(title = "",
       x = NULL, y = "Attributable fraction",
       fill = "")  # Empty legend title for attributable fraction


##########################################
ggplot(data_melted, aes(x = SN_types, y = new_value, fill = variable, label = new_value)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
  geom_vline(xintercept = seq(1.5, 19.5, by = 3), color = "black", linetype = "dashed") + # If you want to include a vertical line
  geom_hline(yintercept = 100) +
  geom_text(aes(angle = 0 + (angle_increment.AF * (as.numeric(as.factor(SN_types)) - 1))),  # Adjust angle of text labels
            position = position_dodge(width = 1), vjust = 0.9, size = 5) +  # Adjust label position and size
  scale_fill_manual(values = c(darker_sjlife, darker_ccss), labels = c("SJLIFE", "CCSS")) +  # Specify colors for SJLIFE and CCSS bars
  scale_y_continuous(breaks = seq(0, 100, by = 10)) +  # Custom grid lines from 0 to 1 by 0.25
  coord_polar() +
  theme_minimal() +
  theme(panel.grid.major.y = element_line(color = "gray", linetype = "dotted"),  # Adjust major grid lines
        panel.grid.minor = element_blank(),  # Remove minor grid lines
        axis.line = element_blank(),   # Remove axis lines
        axis.text.x = element_text(hjust = 0.5, angle = angle_increment),  # Center SN type labels within grid
        axis.text = element_text(size = 12),  # Adjust axis text size
        axis.title = element_text(size = 14),  # Adjust axis title size
        legend.title = element_text(size = 14),  # Adjust legend title size
        legend.text = element_text(size = 12)) +  # Adjust legend text size
  labs(title = "",
       x = NULL, y = "Attributable fraction (%)",
       fill = "Cohort")

###############################################
library(geomtextpath)
ggplot(data_melted, aes(x = SN_types, y = new_value, 
                        fill = variable, label = new_value)) +
  geom_col(position = position_dodge(width = 0.9)) +
  geom_vline(xintercept = seq(1.5, 19.5, by = 3), color = "black", linetype = "dashed") + # If you want to include a vertical line
  geom_hline(yintercept = seq(0, 100, by = 10), color = "grey", linetype = "dotted") +  # Add dotted horizontal lines at y-axis labels
  geom_textpath(position = position_dodge(width = 1), vjust = 1, size = 5) + 
            position = position_dodge(width = 1), vjust = 0.9, size = 5, color = "black") +  # Adjust label position and size
  scale_fill_manual(values = c(darker_sjlife, darker_ccss), labels = c("SJLIFE", "CCSS")) +  # Specify colors for SJLIFE and CCSS bars
  scale_y_continuous(breaks = seq(0, 100, by = 10)) +  # Custom grid lines from 0 to 1 by 0.25
  coord_polar() +
  theme_minimal() +
  theme(panel.grid.major.y = element_blank(),  # Remove major grid lines
        panel.grid.minor = element_blank(),  # Remove minor grid lines
        axis.line = element_blank(),   # Remove axis lines
        axis.text.x = element_text(hjust = 0.5, angle = angle_increment, color = "black", size = 14),  # Center SN type labels within grid
        axis.text.y = element_text(size = 12, color = "black"),  # Adjust axis text size and color
        axis.title = element_text(size = 14, color = "black"),  # Adjust axis title size and color
        legend.title = element_text(size = 14, color = "black"),  # Adjust legend title size and color
        legend.text = element_text(size = 12, color = "black")) +  # Adjust legend text size and color
  labs(title = "",
       x = NULL, y = "Attributable fraction (%)",
       fill = "Cohort")


###############################################
ggplot(data_melted, aes(x = SN_types, y = new_value, fill = variable, label = new_value)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
  geom_vline(xintercept = seq(1.5, 19.5, by = 3), color = "black", linetype = "dashed") + # If you want to include a vertical line
  geom_hline(yintercept = 100) +
  geom_text(
    aes(
      # angle = 0 + (angle_increment.AF * (as.numeric(as.factor(SN_types)) - 1))
    ), # Adjust angle of text labels
    position = position_dodge(width = 1), vjust = 0, size = 5, color = "black"
  ) + # Adjust label position and size
  # scale_fill_manual(values = c(darker_sjlife, darker_ccss), labels = c("SJLIFE", "CCSS")) + # Specify colors for SJLIFE and CCSS bars
  scale_y_continuous(breaks = seq(0, 100, by = 10)) + # Custom grid lines from 0 to 1 by 0.25
  # coord_polar() +
  coord_radial(rotate_angle = TRUE, expand = FALSE) +
  guides(
    theta = guide_axis_theta(angle = 0)
  ) +
  theme_minimal() +
  theme(
    panel.grid.major.y = element_line(color = "gray", linetype = "dotted"), # Adjust major grid lines
    panel.grid.minor = element_blank(), # Remove minor grid lines
    axis.line = element_blank(), # Remove axis lines
    # axis.text.x = element_text(hjust = 0.5, angle = angle_increment, color = "black"), # Center SN type labels within grid
    axis.text = element_text(size = 12, color = "black"), # Adjust axis text size
    axis.title = element_text(size = 14, color = "black"), # Adjust axis title size
    legend.title = element_text(size = 14, color = "black"), # Adjust legend title size
    legend.text = element_text(size = 12, color = "black")
  ) + # Adjust legend text size
  labs(
    title = "",
    x = NULL, y = "(%)",
    fill = "Cohort"
  )


############################################
## Get radial plot
data_melted <- data[grepl("Combined", data$Variables), c("Cohort", "SN_types", "Variables", group)]
colnames(data_melted) <- c("variable", "SN_types", "AF_by", "value")
data_melted$value <- as.numeric(data_melted$value); data_melted <- data_melted[complete.cases(data_melted),]
data_melted$variable <- factor(data_melted$variable, levels = c("SJLIFE", "CCSS"))
data_melted$SN_types <- factor(data_melted$SN_types, levels = unique(data_melted$SN_types))
data_melted$new_value <- round(data_melted$value,2)*100

library(geomtextpath)
# Some of the numbers have to be upside down in a radial plot. The numbers on the
# bottom half of the plot are flipped to make them readable. If you want them
# all to be in the same orientation as the bars you can add upright = FALSE
# inside geom_textpath

p <- ggplot(data_melted, aes(x = SN_types, y = new_value, 
                             fill = variable, label = new_value)) +
  geom_col(position = position_dodge(width = 0.9)) +
  geom_vline(xintercept = seq(1.5, 19.5, by = 3), linetype = "dashed") +
  geom_textpath(position = position_dodge(width = 1), vjust = 1, size = 5, upright = TRUE) + 
  # scale_fill_manual(values = c('deepskyblue4', 'orangered'),
  scale_fill_manual(values = c('#1E90FF', '#FF6347'),
                    labels = c("SJLIFE", "CCSS")) + 
  scale_y_continuous(breaks = seq(0, 100, by = 10)) + 
  coord_curvedpolar() +
  theme_minimal() +
  theme(panel.grid.major.y = element_line(color = "gray", linetype = "dotted"), 
        panel.grid.minor = element_blank(),
        axis.line = element_blank(),  
        axis.text = element_text(size = 14, color = "black"), 
        axis.title = element_text(size = 14, color = "black"), 
        legend.title = element_text(size = 14, color = "black"), 
        legend.text = element_text(size = 12, color = "black")) + 
  labs(x = NULL, y = "Attributable fraction (%)", fill = "Cohort")

p

plot_name <- "Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/ANALYSIS/results/plots/v19b/figure_1_radial_plot.tiff"


ggsave(plot_name, p, width = 12, height = 8, dpi = 600, device = "tiff", compression = "lzw")