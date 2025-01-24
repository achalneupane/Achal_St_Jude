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
  scale_fill_manual(values = c('#1E90FF', '#FF6347'), labels = c("SJLIFE", "CCSS")) + 
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
    axis.text.y = element_blank(), # Remove y-axis labels
    axis.text = element_text(size = 12, color = "black"), # Adjust axis text size
    axis.title = element_text(size = 14, color = "black"), # Adjust axis title size
    legend.title = element_text(size = 14, color = "black"), # Adjust legend title size
    legend.text = element_text(size = 12, color = "black")
  ) + # Adjust legend text size
  labs(
    title = "",
    x = NULL, y = "AF (%)",
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

#################


p = ggplot(data_melted, aes(x = SN_types, y = new_value, fill = variable, label = new_value)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
  geom_vline(xintercept = seq(1.5, 19.5, by = 3), color = "black", linetype = "dashed") + # If you want to include a vertical line
  geom_hline(yintercept = 0, color = "black", linetype = "solid") + # Add horizontal line at y = 0
  geom_text(aes(), position = position_dodge(width = 1), vjust = 0, size = 5, color = "black") + # Adjust label position and size
  scale_fill_manual(values = c('#1E90FF', '#FF6347'), labels = c("SJLIFE", "CCSS")) + 
  coord_radial(rotate_angle = TRUE, expand = FALSE) + # Align y-axis labels with 0
  guides(
    theta = guide_axis_theta(angle = 0)
  ) +
  theme_minimal() +
  theme(
    panel.grid.major.y = element_line(color = "gray", linetype = "dotted"), # Adjust major grid lines
    panel.grid.minor = element_blank(), # Remove minor grid lines
    axis.line = element_blank(), # Remove axis lines
    axis.text = element_text(size = 18, color = "black"), # Adjust axis text size
    axis.title = element_text(size = 18, color = "black"), # Adjust axis title size
    legend.title = element_text(size = 14, color = "black"), # Adjust legend title size
    legend.text = element_text(size = 14, color = "black")
  ) + # Adjust legend text size
  labs(
    title = "",
    x = NULL, y = "Attributable fraction (%)",
    fill = "Cohort"
  )

p



p = ggplot(data_melted, aes(x = SN_types, y = new_value, fill = variable, label = paste0(new_value, "%"))) +
  geom_col(position = position_dodge(width = 0.9)) +
  geom_vline(xintercept = seq(1.5, 19.5, by = 3), linetype = "dashed") +
  geom_textpath(position = position_dodge(width = 1), vjust = 1, size = 5, upright = TRUE) + 
  scale_fill_manual(values = c('#1E90FF', '#FF6347'),
                    labels = c("SJLIFE", "CCSS")) + 
  # scale_y_continuous(breaks = seq(0, 100, by = 10)) + # Format y-axis labels as percentages
  coord_curvedpolar() +
  theme_minimal() +
  theme(panel.grid.major.y = element_line(color = "gray", linetype = "dotted"), 
        panel.grid.minor = element_blank(),
        axis.line = element_blank(),  
        axis.text = element_text(size = 14, color = "black"), 
        axis.title = element_text(size = 14, color = "black"), 
        legend.title = element_text(size = 14, color = "black"), 
        legend.text = element_text(size = 12, color = "black")) + 
  labs(x = NULL, y = "Attributable fraction", fill = "Cohort") # Removed "%" from y-axis label
p

plot_name <- "Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/ANALYSIS/results/plots/v19b/figure_1_radial_plot.tiff"


ggsave(plot_name, p, width = 12, height = 8, dpi = 600, device = "tiff", compression = "lzw")



###############################################

# library
library(tidyverse)

# Create dataset
data <- data.frame(
  individual=paste( "Mister ", seq(1,60), sep=""),
  group=as.factor(c( rep('A', 10), rep('B', 30), rep('C', 14), rep('D', 6))) ,
  value=sample( seq(10,100), 60, replace=T)
)


data <- data_melted$variable

# Set a number of 'empty bar' to add at the end of each group
empty_bar <- 3
to_add <- data.frame( matrix(NA, empty_bar*nlevels(data$group), ncol(data)) )
colnames(to_add) <- colnames(data)
to_add$group <- rep(levels(data$group), each=empty_bar)
data <- rbind(data, to_add)
data <- data %>% arrange(group)
data$id <- seq(1, nrow(data))

# Get the name and the y position of each label
label_data <- data
number_of_bar <- nrow(label_data)
angle <- 90 - 360 * (label_data$id-0.5) /number_of_bar     # I substract 0.5 because the letter must have the angle of the center of the bars. Not extreme right(1) or extreme left (0)
label_data$hjust <- ifelse( angle < -90, 1, 0)
label_data$angle <- ifelse(angle < -90, angle+180, angle)

# prepare a data frame for base lines
base_data <- data %>% 
  dplyr::group_by(group) %>% 
  dplyr::summarize(start=min(id), end=max(id) - empty_bar) %>% 
  dplyr::rowwise() %>% 
  dplyr::mutate(title=mean(c(start, end)))

# prepare a data frame for grid (scales)
grid_data <- base_data
grid_data$end <- grid_data$end[ c( nrow(grid_data), 1:nrow(grid_data)-1)] + 1
grid_data$start <- grid_data$start - 1
grid_data <- grid_data[-1,]

# Make the plot
ggplot(data, aes(x=as.factor(id), y=value, fill=group)) +       # Note that id is a factor. If x is numeric, there is some space between the first bar
  
  geom_bar(aes(x=as.factor(id), y=value, fill=group), stat="identity", alpha=0.5) +
  
  # Add a val=100/75/50/25 lines. I do it at the beginning to make sur barplots are OVER it.
  geom_segment(data=grid_data, aes(x = end, y = 80, xend = start, yend = 80), colour = "white", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  geom_segment(data=grid_data, aes(x = end, y = 60, xend = start, yend = 60), colour = "white", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  geom_segment(data=grid_data, aes(x = end, y = 40, xend = start, yend = 40), colour = "white", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  geom_segment(data=grid_data, aes(x = end, y = 20, xend = start, yend = 20), colour = "white", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  
  # Add text showing the value of each 100/75/50/25 lines
  annotate("text", x = rep(max(data$id),4), y = c(20, 40, 60, 80), label = c("20", "40", "60", "80") , color="grey", size=3 , angle=0, fontface="bold", hjust=1) +
  
  geom_bar(aes(x=as.factor(id), y=value, fill=group), stat="identity", alpha=0.5) +
  ylim(-100,120) +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    plot.margin = unit(rep(-1,4), "cm") 
  ) +
  coord_polar() + 
  geom_text(data=label_data, aes(x=id, y=value+10, label=individual, hjust=hjust), color="black", fontface="bold",alpha=0.6, size=2.5, angle= label_data$angle, inherit.aes = FALSE ) +
  
  # Add base line information
  geom_segment(data=base_data, aes(x = start, y = -5, xend = end, yend = -5), colour = "black", alpha=0.8, size=0.6 , inherit.aes = FALSE )  +
  geom_text(data=base_data, aes(x = title, y = -18, label=group), hjust=c(1,1,0,0), colour = "black", alpha=0.8, size=4, fontface="bold", inherit.aes = FALSE)

##############################################
##########################
## Radial plot by group ##
##########################
# library
library(tidyverse)

data_melted <- structure(list(variable = c("SJLIFE", "SJLIFE", "SJLIFE", "SJLIFE", 
                                           "SJLIFE", "SJLIFE", "SJLIFE", "SJLIFE", "SJLIFE", "SJLIFE", "SJLIFE", 
                                           "SJLIFE", "SJLIFE", "SJLIFE", "CCSS", "CCSS", "CCSS", "CCSS", 
                                           "CCSS", "CCSS", "CCSS", "CCSS", "CCSS", "CCSS", "CCSS", "CCSS", 
                                           "CCSS", "CCSS"), SN_types = c("SNs", "SNs", "SMNs", "SMNs", "NMSC", 
                                                                         "NMSC", "Breast cancer", "Breast cancer", "Thyroid cancer", "Thyroid cancer", 
                                                                         "Meningioma", "Meningioma", "Sarcoma", "Sarcoma", "SNs", "SNs", 
                                                                         "SMNs", "SMNs", "NMSC", "NMSC", "Breast cancer", "Breast cancer", 
                                                                         "Thyroid cancer", "Thyroid cancer", "Meningioma", "Meningioma", 
                                                                         "Sarcoma", "Sarcoma"), AF_by = c("All_treatments", "PRS", "All_treatments", 
                                                                                                          "PRS", "All_treatments", "PRS", "All_treatments", "PRS", "All_treatments", 
                                                                                                          "PRS", "All_treatments", "PRS", "All_treatments", "PRS", "All_treatments", 
                                                                                                          "PRS", "All_treatments", "PRS", "All_treatments", "PRS", "All_treatments", 
                                                                                                          "PRS", "All_treatments", "PRS", "All_treatments", "PRS", "All_treatments", 
                                                                                                          "PRS"), value = c(0.469, 0.125, 0.392, 0.122, 0.437, 0.276, 0.602, 
                                                                                                                            0.255, 0.726, 0.517, 0.489, NA, 0.324, NA, 0.407, 0.048, 0.29, 
                                                                                                                            0.047, 0.393, 0.31, 0.593, 0.365, 0.476, 0.358, 0.426, 0.032, 
                                                                                                                            0.342, 0.024), AF_by_new = c("All treatments", "PRS", "All treatments", 
                                                                                                                                                         "PRS", "All treatments", "PRS", "All treatments", "PRS", "All treatments", 
                                                                                                                                                         "PRS", "All treatments", "PRS", "All treatments", "PRS", "All treatments", 
                                                                                                                                                         "PRS", "All treatments", "PRS", "All treatments", "PRS", "All treatments", 
                                                                                                                                                         "PRS", "All treatments", "PRS", "All treatments", "PRS", "All treatments", 
                                                                                                                                                         "PRS"), legend_group = structure(c(3L, 4L, 3L, 4L, 3L, 4L, 3L, 
                                                                                                                                                                                            4L, 3L, 4L, 3L, 4L, 3L, 4L, 1L, 2L, 1L, 2L, 1L, 2L, 1L, 2L, 1L, 
                                                                                                                                                                                            2L, 1L, 2L, 1L, 2L), levels = c("CCSS All treatments", "CCSS PRS", 
                                                                                                                                                                                                                            "SJLIFE All treatments", "SJLIFE PRS"), class = "factor"), new_value = c(47, 
                                                                                                                                                                                                                                                                                                     12, 39, 12, 44, 28, 60, 26, 73, 52, 49, NA, 32, NA, 41, 5, 29, 
                                                                                                                                                                                                                                                                                                     5, 39, 31, 59, 36, 48, 36, 43, 3, 34, 2)), row.names = c(3L, 
                                                                                                                                                                                                                                                                                                                                                              4L, 9L, 10L, 15L, 16L, 21L, 22L, 27L, 28L, 33L, 34L, 39L, 40L, 
                                                                                                                                                                                                                                                                                                                                                              45L, 46L, 51L, 52L, 57L, 58L, 63L, 64L, 69L, 70L, 75L, 76L, 81L, 
                                                                                                                                                                                                                                                                                                                                                              82L), class = "data.frame")


# Define the desired order of SN_types within each variable
custom_order_within_variable <- c("SNs", "SMNs", "NMSC", "Breast cancer", "Thyroid cancer", "Meningioma", "Sarcoma")

# Create a custom factor variable based on the desired order
data_melted <- data_melted %>%
  dplyr::mutate(SN_types = factor(SN_types, levels = custom_order_within_variable))
data_melted <- data_melted %>%
  arrange(SN_types, AF_by)

data <- cbind.data.frame(individual= paste0(data_melted$new_value, " %"), group= data_melted$SN_types, group2= data_melted$legend_group, value=data_melted$new_value)
# data$value[is.na(data$value)]

# Set a number of 'empty bar' to add at the end of each group
empty_bar <- 3
to_add <- data.frame( matrix(NA, empty_bar*nlevels(data$group), ncol(data)) )
colnames(to_add) <- colnames(data)
to_add$group <- rep(levels(data$group), each=empty_bar)
data <- rbind(data, to_add)
data <- data %>% arrange(group)
data$id <- seq(1, nrow(data))

# Get the name and the y position of each label
label_data <- data
number_of_bar <- nrow(label_data)
angle <- 90 - 360 * (label_data$id-0.5) /number_of_bar     # I substract 0.5 because the letter must have the angle of the center of the bars. Not extreme right(1) or extreme left (0)
label_data$hjust <- ifelse( angle < -90, 1, 0)
label_data$angle <- ifelse(angle < -90, angle+180, angle)


# prepare a data frame for base lines
base_data <- data %>% 
  dplyr::group_by(group) %>% 
  dplyr::summarize(start=min(id), end=max(id) - empty_bar) %>% 
  dplyr::rowwise() %>% 
  dplyr::mutate(title=mean(c(start, end)))

# prepare a data frame for grid (scales)
grid_data <- base_data
grid_data$end <- grid_data$end[ c( nrow(grid_data), 1:nrow(grid_data)-1)] + 1
grid_data$start <- grid_data$start - 1
grid_data <- grid_data[-1,]


custom_colors <- c("SJLIFE PRS" = "#87CEFA", "CCSS PRS" = "#FFA07A", "SJLIFE All treatments" = "#1E90FF", "CCSS All treatments" = "#FF6347")
legend_order <- c("SJLIFE PRS", "CCSS PRS", "SJLIFE All treatments", "CCSS All treatments")

# Make the plot
ggplot(data, aes(x=as.factor(id), y=value, fill=group2)) +       # Note that id is a factor. If x is numeric, there is some space between the first bar
  
  geom_bar(aes(x=as.factor(id), y=value, fill=group2), stat="identity", alpha=0.5) +
  
  # Add a val=100/75/50/25 lines. I do it at the beginning to make sur barplots are OVER it.
  geom_segment(data=grid_data, aes(x = end, y = 80, xend = start, yend = 80), colour = "white", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  geom_segment(data=grid_data, aes(x = end, y = 60, xend = start, yend = 60), colour = "white", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  geom_segment(data=grid_data, aes(x = end, y = 40, xend = start, yend = 40), colour = "white", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  geom_segment(data=grid_data, aes(x = end, y = 20, xend = start, yend = 20), colour = "white", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  
  ## Add text showing the value of each 100/75/50/25 lines
  # annotate("text", x = rep(max(data$id),5), y = c(0, 25, 50, 75, 100), label = c("0%", "25%", "50%", "75%", "100%") , color="black", size=3 , angle=0, fontface= 1, hjust=1) +
  
  geom_bar(aes(x=as.factor(id), y=value, fill=group2), stat="identity", alpha=0.5) +
  ylim(-100,120) +
  theme_minimal() +
  theme(
    legend.position = "right",
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    plot.margin = unit(rep(-1,4), "cm") 
  ) +
  coord_polar() + 
  geom_text(data=label_data, aes(x=id, y=value+10, label=individual, hjust=hjust), color="black", fontface=1,alpha=1, size=4.5, angle= label_data$angle, inherit.aes = FALSE ) +
  
  # Add base line information
  # geom_segment(data=base_data, aes(x = start, y = -5, xend = end, yend = -5), colour = "black", alpha=0.8, size=1.5 , inherit.aes = FALSE )  +
  # geom_text(data=base_data, aes(x = title, y = -18, label=group), hjust=c(1,1,1,0,0,0,0), colour = "black", alpha=0.8, size=4.5, fontface=1, inherit.aes = FALSE) +
  geomtextpath::coord_curvedpolar() +
  geomtextpath::geom_textsegment(
    data = base_data,
    aes(
      x = start, y = -.1,
      xend = end, yend = -.1,
      label = group
    ),
    colour = "black",
    linewidth = 3,
    size = 16 / .pt,
    inherit.aes = FALSE,
    gap = FALSE,
    offset = unit(-24, "pt")
  ) +
  scale_fill_manual(values = custom_colors, breaks = legend_order) +
  guides(
    theta = guide_axis_theta(angle = 0)
  ) +
  # theme_minimal() +
  theme(
    # panel.grid.major.y = element_line(color = "gray", linetype = "dotted"), # Adjust major grid lines
    # panel.grid.minor = element_blank(), # Remove minor grid lines
    # axis.line = element_blank(), # Remove axis lines
    # axis.text.y = element_blank(), # Remove y-axis labels,
    legend.position = "right",
    # axis.text = element_text(size = 18, color = "black"), # Adjust axis text size
    # axis.title = element_text(size = 18, color = "black"), # Adjust axis title size
    legend.title = element_text(size = 14, color = "black", face=2), # Adjust legend title size
    legend.text = element_text(size = 14, color = "black")
  ) + # Adjust legend text size
  labs(
    title = "",
    x = NULL, y = "",
    fill = "Cohort and variables",
  )




##########################################################################


#########################
## Figure 1 : Option 2 ##
#########################
group <- "Overall"
variables <- unique(data$Variables)

custom_colors <- c("SJLIFE" = "darkred", "CCSS" = "darkblue")
legend_order <- c("SJLIFE", "CCSS")


for(i in 1:length(variables)){
  data_melted <- data[grepl(variables[i], data$Variables), c("Cohort", "SN_types", "Variables", group)]
  colnames(data_melted) <- c("variable", "SN_types", "AF_by", "value")
  data_melted$value <- as.numeric(data_melted$value); data_melted <- data_melted[complete.cases(data_melted),]
  data_melted$new_value <- round(data_melted$value,2)*100
  
  # Reshape the data for ggplot2
  library(reshape2)
  order <- unique(data_melted$SN_types)
  AF.type <- data_melted$AF_by[i] 
  lifestyle <- "without_lifestyle"
  if(lifestyle== "without_lifestyle" && variables[i] == "Lifestyle"){
    next
  }
  
  data_melted$legend_group <- factor(data_melted$variable, levels = c("SJLIFE", "CCSS"))
  
  # Create a factor with the desired order
  data_melted$SN_types <- factor(data_melted$SN_types, levels = order)
  
  
  
  ## Radial plot
  i=6
  data_melted <- data[grepl(variables[i], data$Variables), c("Cohort", "SN_types", "Variables", group)]
  colnames(data_melted) <- c("variable", "SN_types", "AF_by", "value")
  data_melted$value <- as.numeric(data_melted$value); data_melted <- data_melted[complete.cases(data_melted),]
  data_melted$new_value <- round(data_melted$value,2)*100
  
  # Reshape the data for ggplot2
  library(reshape2)
  order <- unique(data_melted$SN_types)
  AF.type <- data_melted$AF_by[i] 
  lifestyle <- "without_lifestyle"
  if(lifestyle== "without_lifestyle" && variables[i] == "Lifestyle"){
    next
  }
  
  data_melted$legend_group <- factor(data_melted$variable, levels = c("SJLIFE", "CCSS"))
  
  # Create a factor with the desired order
  data_melted$SN_types <- factor(data_melted$SN_types, levels = order)
  
  p= ggplot(data_melted, aes(x = SN_types, y = new_value, fill = variable, label = paste0(new_value, "%"))) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
    geom_vline(xintercept = seq(1.5, 19.5, by = 3), color = "black", linetype = "dashed") + # If you want to include a vertical line
    geom_hline(yintercept = 100) +
    geom_text(
      aes(
        # angle = 0 + (angle_increment.AF * (as.numeric(as.factor(SN_types)) - 1))
      ), # Adjust angle of text labels
      position = position_dodge(width = 1), vjust = -0.5, size = 5, color = "black"
    ) + # Adjust label position and size
    # scale_fill_manual(values = c(darker_sjlife, darker_ccss), labels = c("SJLIFE", "CCSS")) + # Specify colors for SJLIFE and CCSS bars
    scale_fill_manual(values = c('darkred', 'darkblue'), labels = c("SJLIFE", "CCSS")) + 
    # scale_y_continuous(breaks = seq(0, 100, by = 10)) + # Custom grid lines from 0 to 1 by 0.25
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
      axis.text.y = element_blank(), # Remove y-axis labels
      axis.text = element_text(size = 18, color = "black"), # Adjust axis text size
      axis.title = element_text(size = 18, color = "black"), # Adjust axis title size
      legend.title = element_text(size = 14, color = "black"), # Adjust legend title size
      legend.text = element_text(size = 14, color = "black")
    ) + # Adjust legend text size
    labs(
      title = "",
      x = NULL, y = "",
      fill = "Cohort"
    )
  p
  
  plot_name <- "Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/ANALYSIS/results/plots/v19b/overall/Overall_Combined_without_lifestyle_radial_plot.tiff"
  ggsave(plot_name, p, width = 12, height = 8, dpi = 600, device = "tiff", compression = "lzw")
  


##############
## Figure 2 ##
##############
## Radial plot
group <- "Overall"

variables <- "treatment|PRS"

new_data <- saved.data[grepl(variables, saved.data$Variables), c("Cohort", "SN_types", "Variables", group)]
colnames(new_data) <- c("variable", "SN_types", "AF_by", "value")
new_data$value <- as.numeric(new_data$value)

# Reshape the data for ggplot2
library(reshape2)

data_melted <-new_data
order <- unique(data_melted$SN_types)
AF.type <- "treatment_PRS"
lifestyle <- "without_lifestyle"

data_melted$AF_by_new <-  gsub("_", " ", data_melted$AF_by)
data_melted$legend_group <- factor(paste(data_melted$variable, data_melted$AF_by_new, sep = "-"))

data_melted$new_value <- round(data_melted$value,2)*100


# Define the desired order of SN_types within each variable
custom_order_within_variable <- c("SNs", "SMNs", "NMSC", "Breast cancer", "Thyroid cancer", "Meningioma", "Sarcoma")

# Create a custom factor variable based on the desired order
data_melted <- data_melted %>%
  dplyr::mutate(SN_types = factor(SN_types, levels = custom_order_within_variable))
data_melted <- data_melted %>%
  arrange(SN_types, AF_by)

data <- cbind.data.frame(individual= paste0(data_melted$new_value, " %"), group= data_melted$SN_types, group2= data_melted$legend_group, value=data_melted$new_value)
# data$value[is.na(data$value)]

# Set a number of 'empty bar' to add at the end of each group
empty_bar <- 1
to_add <- data.frame( matrix(NA, empty_bar*nlevels(data$group), ncol(data)) )
colnames(to_add) <- colnames(data)
to_add$group <- rep(levels(data$group), each=empty_bar)
data <- rbind(data, to_add)
data <- data %>% arrange(group)
data$id <- seq(1, nrow(data))

# Get the name and the y position of each label
label_data <- data
number_of_bar <- nrow(label_data)
angle <- 90 - 360 * (label_data$id-0.5) /number_of_bar     # I substract 0.5 because the letter must have the angle of the center of the bars. Not extreme right(1) or extreme left (0)
label_data$hjust <- ifelse( angle < -90, 1, 0)
label_data$angle <- ifelse(angle < -90, angle+180, angle)


# prepare a data frame for base lines
base_data <- data %>% 
  dplyr::group_by(group) %>% 
  dplyr::summarize(start=min(id), end=max(id) - empty_bar) %>% 
  dplyr::rowwise() %>% 
  dplyr::mutate(title=mean(c(start, end)))

# prepare a data frame for grid (scales)
grid_data <- base_data
grid_data$end <- grid_data$end[ c( nrow(grid_data), 1:nrow(grid_data)-1)] + 1
grid_data$start <- grid_data$start - 1
grid_data <- grid_data[-1,]


custom_colors <- c("SJLIFE-All treatments" = "darkred", "CCSS-All treatments" = "darkblue", "SJLIFE-PRS" = "red", "CCSS-PRS" = "blue")
legend_order <- c("SJLIFE-All treatments", "CCSS-All treatments", "SJLIFE-PRS", "CCSS-PRS")

# Make the plot
p = ggplot(data, aes(x=as.factor(id), y=value, fill=group2)) +       # Note that id is a factor. If x is numeric, there is some space between the first bar
  
  geom_bar(aes(x=as.factor(id), y=value, fill=group2), stat="identity", alpha=0.5) +
  
  # Add a val=100/75/50/25 lines. I do it at the beginning to make sur barplots are OVER it.
  geom_segment(data=grid_data, aes(x = end, y = 80, xend = start, yend = 80), colour = "white", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  geom_segment(data=grid_data, aes(x = end, y = 60, xend = start, yend = 60), colour = "white", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  geom_segment(data=grid_data, aes(x = end, y = 40, xend = start, yend = 40), colour = "white", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  geom_segment(data=grid_data, aes(x = end, y = 20, xend = start, yend = 20), colour = "white", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  
  ## Add text showing the value of each 100/75/50/25 lines
  # annotate("text", x = rep(max(data$id),5), y = c(0, 25, 50, 75, 100), label = c("0%", "25%", "50%", "75%", "100%") , color="black", size=3 , angle=0, fontface= 1, hjust=1) +
  
  geom_bar(aes(x=as.factor(id), y=value, fill=group2), stat="identity", alpha=0.5) +
  ylim(-100,120) +
  theme_minimal() +
  theme(
    legend.position = "right",
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    plot.margin = unit(rep(-1,4), "cm") 
  ) +
  coord_polar() + 
  geom_text(data=label_data, aes(x=id, y=value+10, label=individual, hjust=hjust), color="black", fontface=1,alpha=1, size=4.5, angle= label_data$angle, inherit.aes = FALSE ) +
  
  # Add base line information
  # geom_segment(data=base_data, aes(x = start, y = -5, xend = end, yend = -5), colour = "black", alpha=0.8, size=1.5 , inherit.aes = FALSE )  +
  # geom_text(data=base_data, aes(x = title, y = -18, label=group), hjust=c(1,1,1,0,0,0,0), colour = "black", alpha=0.8, size=4.5, fontface=1, inherit.aes = FALSE) +
  geomtextpath::coord_curvedpolar() +
  geomtextpath::geom_textsegment(
    data = base_data,
    aes(
      x = start, y = -.1,
      xend = end, yend = -.1,
      label = group
    ),
    colour = "black",
    linewidth = 1.5,
    size = 16 / .pt,
    inherit.aes = FALSE,
    gap = FALSE,
    offset = unit(-24, "pt")
  ) +
  scale_fill_manual(values = custom_colors, breaks = legend_order) +
  guides(
    theta = guide_axis_theta(angle = 0)
  ) +
  # theme_minimal() +
  theme(
    # panel.grid.major.y = element_line(color = "gray", linetype = "dotted"), # Adjust major grid lines
    # panel.grid.minor = element_blank(), # Remove minor grid lines
    # axis.line = element_blank(), # Remove axis lines
    # axis.text.y = element_blank(), # Remove y-axis labels,
    legend.position="right",
    legend.justification="right", 
    legend.box.spacing = unit(-120, "pt"),# The spacing between the plotting area and the legend box (unit)
    legend.margin=margin(0,0,0,0),
    # axis.text = element_text(size = 18, color = "black"), # Adjust axis text size
    # axis.title = element_text(size = 18, color = "black"), # Adjust axis title size
    legend.title = element_text(size = 16, color = "black", face=2), # Adjust legend title size
    plot.margin = margin(r = -3.5, l = -3.5, t = -3.5, b = -3.5, unit = "cm"),
    legend.text = element_text(size = 14, color = "black")
  ) + # Adjust legend text size
  labs(
    title = "",
    x = NULL, y = "",
    fill = "Cohorts and variables"
  )
p
plot_name <- "Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/ANALYSIS/results/plots/v19b/overall/Overall_treatment_and_PRS_without_lifestyle_radial_plot.tiff"
ggsave(plot_name, p, width = 12, height = 8, dpi = 600, device = "tiff", compression = "lzw")

