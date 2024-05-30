# Load required libraries
library(ggplot2)
library(ggthemes)# For additional themes
library(ggrepel)
library(dplyr)



data <- read.table(text="Cohort	SN_types	Variables	Overall	Female	Male	Age.lt.35	Age.ge.35
SJLIFE	SNs (605)	Radiation	0.425	0.417	0.435	0.4	0.447
SJLIFE	SNs (605)	Chemotherapy	0.08	0.079	0.08	0.109	0.054
SJLIFE	SNs (605)	All_treatments	0.469	0.462	0.478	0.463	0.474
SJLIFE	SNs (605)	PRS	0.125	0.123	0.129	0.126	0.125
SJLIFE	SNs (605)	Lifestyle	-	-	-	-	-
SJLIFE	SNs (605)	Combined	0.536	0.53	0.545	0.531	0.541
SJLIFE	SMNs (463)	Radiation	0.37	0.369	0.371	0.339	0.394
SJLIFE	SMNs (463)	Chemotherapy	0.032	0.032	0.032	0.045	0.021
SJLIFE	SMNs (463)	All_treatments	0.39	0.39	0.39	0.368	0.407
SJLIFE	SMNs (463)	PRS	0.117	0.115	0.12	0.117	0.117
SJLIFE	SMNs (463)	Lifestyle	-	-	-	-	-
SJLIFE	SMNs (463)	Combined	0.462	0.461	0.464	0.442	0.478
SJLIFE	NMSCs (251)	Radiation	0.436	0.427	0.448	0.408	0.452
SJLIFE	NMSCs (251)	Chemotherapy	-	-	-	-	-
SJLIFE	NMSCs (251)	All_treatments	0.436	0.427	0.448	0.408	0.452
SJLIFE	NMSCs (251)	PRS	0.28	0.266	0.296	0.278	0.281
SJLIFE	NMSCs (251)	Lifestyle	-	-	-	-	-
SJLIFE	NMSCs (251)	Combined	0.593	0.581	0.607	0.571	0.605
SJLIFE	Breast cancer (76)	Radiation	0.508	0.508	-	0.491	0.513
SJLIFE	Breast cancer (76)	Chemotherapy	0.189	0.189	-	0.222	0.18
SJLIFE	Breast cancer (76)	All_treatments	0.613	0.613	-	0.604	0.615
SJLIFE	Breast cancer (76)	PRS	0.192	0.192	-	0.197	0.191
SJLIFE	Breast cancer (76)	Lifestyle	-	-	-	-	-
SJLIFE	Breast cancer (76)	Combined	0.687	0.687	-	0.683	0.689
SJLIFE	Thyroid cancer (87)	Radiation	0.62	0.623	0.617	0.594	0.651
SJLIFE	Thyroid cancer (87)	Chemotherapy	0.235	0.236	0.233	0.294	0.165
SJLIFE	Thyroid cancer (87)	All_treatments	0.728	0.729	0.726	0.727	0.73
SJLIFE	Thyroid cancer (87)	PRS	0.517	0.519	0.514	0.519	0.515
SJLIFE	Thyroid cancer (87)	Lifestyle	-	-	-	-	-
SJLIFE	Thyroid cancer (87)	Combined	0.868	0.868	0.867	0.867	0.869
SJLIFE	Meningioma (149)	Radiation	0.913	0.91	0.917	0.904	0.92
SJLIFE	Meningioma (149)	Chemotherapy	0.241	0.244	0.237	0.287	0.205
SJLIFE	Meningioma (149)	All_treatments	0.93	0.928	0.932	0.929	0.93
SJLIFE	Meningioma (149)	PRS	-	-	-	-	-
SJLIFE	Meningioma (149)	Lifestyle	-	-	-	-	-
SJLIFE	Meningioma (149)	Combined	0.92	0.918	0.922	0.92	0.92
SJLIFE	Sarcoma (33)	Radiation	-	-	-	-	-
SJLIFE	Sarcoma (33)	Chemotherapy	0.346	0.349	0.342	0.329	0.371
SJLIFE	Sarcoma (33)	All_treatments	0.346	0.349	0.342	0.329	0.371
SJLIFE	Sarcoma (33)	PRS	-	-	-	-	-
SJLIFE	Sarcoma (33)	Lifestyle	-	-	-	-	-
SJLIFE	Sarcoma (33)	Combined	0.298	0.301	0.293	0.28	0.324
CCSS	SNs (1611)	Radiation	0.387	0.38	0.398	0.372	0.398
CCSS	SNs (1611)	Chemotherapy	0.031	0.03	0.033	0.047	0.019
CCSS	SNs (1611)	All_treatments	0.407	0.4	0.418	0.403	0.41
CCSS	SNs (1611)	PRS	0.048	0.047	0.049	0.047	0.048
CCSS	SNs (1611)	Lifestyle	-	-	-	-	-
CCSS	SNs (1611)	Combined	0.435	0.428	0.446	0.431	0.438
CCSS	SMNs (762)	Radiation	0.255	0.254	0.257	0.213	0.283
CCSS	SMNs (762)	Chemotherapy	0.042	0.042	0.042	0.07	0.024
CCSS	SMNs (762)	All_treatments	0.29	0.289	0.293	0.271	0.303
CCSS	SMNs (762)	PRS	0.048	0.047	0.048	0.047	0.048
CCSS	SMNs (762)	Lifestyle	-	-	-	-	-
CCSS	SMNs (762)	Combined	0.324	0.323	0.326	0.305	0.336
CCSS	NMSCs (774)	Radiation	0.391	0.387	0.397	0.386	0.394
CCSS	NMSCs (774)	Chemotherapy	-	-	-	-	-
CCSS	NMSCs (774)	All_treatments	0.391	0.387	0.397	0.386	0.394
CCSS	NMSCs (774)	PRS	0.306	0.304	0.309	0.308	0.305
CCSS	NMSCs (774)	Lifestyle	-	-	-	-	-
CCSS	NMSCs (774)	Combined	0.579	0.575	0.583	0.577	0.58
CCSS	Breast cancer (290)	Radiation	0.475	0.475	-	0.454	0.479
CCSS	Breast cancer (290)	Chemotherapy	0.186	0.186	-	0.223	0.179
CCSS	Breast cancer (290)	All_treatments	0.593	0.593	-	0.589	0.594
CCSS	Breast cancer (290)	PRS	0.365	0.365	-	0.368	0.365
CCSS	Breast cancer (290)	Lifestyle	-	-	-	-	-
CCSS	Breast cancer (290)	Combined	0.743	0.743	-	0.741	0.744
CCSS	Thyroid cancer (163)	Radiation	0.442	0.444	0.438	0.413	0.488
CCSS	Thyroid cancer (163)	Chemotherapy	0.056	0.054	0.059	0.073	0.03
CCSS	Thyroid cancer (163)	All_treatments	0.475	0.476	0.472	0.456	0.503
CCSS	Thyroid cancer (163)	PRS	0.358	0.362	0.351	0.355	0.362
CCSS	Thyroid cancer (163)	Lifestyle	-	-	-	-	-
CCSS	Thyroid cancer (163)	Combined	0.663	0.664	0.66	0.649	0.683
CCSS	Meningioma (256)	Radiation	0.864	0.855	0.877	0.87	0.858
CCSS	Meningioma (256)	Chemotherapy	0.082	0.073	0.095	0.11	0.052
CCSS	Meningioma (256)	All_treatments	0.877	0.868	0.889	0.887	0.865
CCSS	Meningioma (256)	PRS	0.013	0.012	0.014	0.016	0.009
CCSS	Meningioma (256)	Lifestyle	-	-	-	-	-
CCSS	Meningioma (256)	Combined	0.878	0.869	0.891	0.889	0.867
CCSS	Sarcoma (61)	Radiation	-	-	-	-	-
CCSS	Sarcoma (61)	Chemotherapy	0.353	0.338	0.369	0.349	0.358
CCSS	Sarcoma (61)	All_treatments	0.353	0.338	0.369	0.349	0.358
CCSS	Sarcoma (61)	PRS	0.041	0.042	0.039	0.04	0.042
CCSS	Sarcoma (61)	Lifestyle	-	-	-	-	-
CCSS	Sarcoma (61)	Combined	0.379	0.367	0.394	0.375	0.385", header = T, sep = "\t")

data$SN_types <- trimws(gsub("\\([0-9]+\\)", "", data$SN_types))
data[data == "-"] <- NA
saved.data <- data

## 1.___________________________________________________________________________ Combined (overall) analysis
## 1.___________________________________________________________________________ Combined (overall) analysis
## 1.___________________________________________________________________________ Combined (overall) analysis
## 1.___________________________________________________________________________ Combined (overall) analysis
## 1.___________________________________________________________________________ Combined (overall) analysis



group <- "Overall"
variables <- unique(data$Variables)

custom_colors <- c("SJLIFE" = "#EE6A50", "CCSS" = "#9FB6CD")
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



p <- ggplot(data_melted, aes(x = SN_types, y = new_value, fill = legend_group)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.4), width = 0.4) +
  # Customize the theme and appearance
  # # geom_hline(yintercept = seq(-0.1, 0.1, by = 0.1), color = "black", linetype = "solid", size = 0.7) +
  # # geom_hline(yintercept = 0, color = "black", linetype = "solid", size = 0.7) +
  # Add Y lines
  # # geom_vline(xintercept = seq(0.5, length(unique(data_melted$SN_types)) + 0.5), color = "black", linetype = "solid", size = 0.7) +
  # # geom_vline(xintercept = 0.5, color = "black", linetype = "solid", size = 0.7) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 16, color = "black"),
        axis.text.y = element_text(size = 16, color = "black"),
        axis.title.y = element_text(size = 20, color = "black"),
        axis.line.x = element_line(color="black"),
        axis.line.y = element_line(color="black"),
        legend.title = element_blank(),
        legend.text = element_text(size = 16, color = "black"),
        plot.title = element_text(size = 16, hjust = 0.5, vjust = 1.5, face = "bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  
  # Customize colors
  scale_fill_manual(values = custom_colors, breaks = legend_order) +
  
  # Add labels with geom_text
  geom_text(
    data = data_melted %>% filter(!is.na(new_value)),
    aes(label = new_value, y = new_value),  # Adjust y position
    position = position_dodge(width = 0.4),
    vjust = -0.20,  # Adjust vertical justification
    hjust = 0.5,  # Center text horizontally
    size = 5.5,
    color = "black"
  ) +
  # Adjust legend position
  theme(legend.position = "top", legend.box = "horizontal") +
  # Additional adjustments
  scale_x_discrete(limits = order) +
  coord_cartesian(clip = "off")
p <- p + scale_y_continuous(limits = c(0, 100), expand = c(0, 0)) + labs(title = "", y = "Attributable fraction (%)", x = NULL) 

p
# Save the plot as a high-resolution image
plot_name <- paste0("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/ANALYSIS/results/plots/v20b/overall/", group, "_", AF.type, "_", lifestyle, ".tiff")

ggsave(plot_name, p, width = 12, height = 8, dpi = 600, device = "tiff", compression = "lzw")
}



################################
## PRS and treatment together ##
################################
group <- "Overall"

variables <- "treatment|PRS"

new_data <- data[grepl(variables, data$Variables), c("Cohort", "SN_types", "Variables", group)]
colnames(new_data) <- c("variable", "SN_types", "AF_by", "value")
new_data$value <- as.numeric(new_data$value)

# Reshape the data for ggplot2
library(reshape2)
data_melted <-new_data
data_melted$new_value <- round(data_melted$value,2)*100

order <- unique(data_melted$SN_types)
AF.type <- "treatment_PRS"
lifestyle <- "without_lifestyle"

data_melted$legend_group <- gsub("_", " ", factor(paste(data_melted$variable, data_melted$AF_by, sep = "-")))


# Define custom colors and legend order
custom_colors <- c("SJLIFE-All treatments" = "#EE6A50", "CCSS-All treatments" = "darkblue", "SJLIFE-PRS" = "red", "CCSS-PRS" = "blue")
legend_order <- c("SJLIFE-All treatments", "CCSS-All treatments", "SJLIFE-PRS", "CCSS-PRS")
# Order the levels of the legend_group factor
data_melted$legend_group <- factor(data_melted$legend_group, levels = legend_order)

# Create a factor with the desired order
data_melted$SN_types <- factor(data_melted$SN_types, levels = order)

# Create the grouped bar chart
# Create the plot
p <- ggplot(data_melted, aes(x = SN_types, y = new_value, fill = legend_group)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.8) +
  # Customize the theme and appearance
  # geom_hline(yintercept = 0, color = "black", linetype = "solid", size = 0.7) +
  # Add Y lines
  # geom_vline(xintercept = 0.5, color = "black", linetype = "solid", size = 0.7) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 16, color = "black"),
        axis.text.y = element_text(size = 16, color = "black"),
        axis.title.y = element_text(size = 20, color = "black"), 
        axis.line.x = element_line(color="black"), 
        axis.line.y = element_line(color="black"),
        legend.title = element_blank(),
        legend.text = element_text(size = 16, color = "black"),
        plot.title = element_text(size = 16, hjust = 0.5, vjust = 1.5, face = "bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  
  # Customize colors
  scale_fill_manual(values = custom_colors, breaks = legend_order) +
  
  # Add labels with geom_text
  geom_text(
    data = data_melted,
    aes(label = new_value, y = new_value),  # Adjust y position
    position = position_dodge(width = 0.8),
    vjust = -0.20,  # Adjust vertical justification
    hjust = 0.5,  # Center text horizontally
    size = 5.5,
    color = "black"
  ) +
  # Adjust legend position
  theme(legend.position = "top", legend.box = "horizontal") +
  # Additional adjustments
  scale_x_discrete(limits = order) +
  coord_cartesian(clip = "off")
p
p <- p + scale_y_continuous(limits = c(0, 100), expand = c(0, 0)) + 
  labs(title = "", y = "Attributable fraction (%)", x =  NULL) +
  theme(axis.title.x = element_text(size = 18, color = "black")) +
  theme(axis.title.y = element_text(size = 18, color = "black"))

p

# Save the plot as a high-resolution image
plot_name <- paste0("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/ANALYSIS/results/plots/v20b/overall/", group, "_", AF.type,"_", lifestyle, ".tiff")

ggsave(plot_name, p, width = 12, height = 8, dpi = 600, device = "tiff", compression = "lzw")



##################################
## Chemo and Radiation together ##
##################################
group <- "Overall"

variables <- "Radiation|Chemotherapy"

new_data <- data[grepl(variables, data$Variables), c("Cohort", "SN_types", "Variables", group)]
colnames(new_data) <- c("variable", "SN_types", "AF_by", "value")
new_data$value <- as.numeric(new_data$value)

# Reshape the data for ggplot2
library(reshape2)
data_melted <-new_data
data_melted$new_value <- round(data_melted$value,2)*100

order <- unique(data_melted$SN_types)
AF.type <- "Chemo_Radiation"
lifestyle <- "without_lifestyle"

data_melted$legend_group <- factor(paste(data_melted$variable, data_melted$AF_by, sep = "-"))


# Define custom colors and legend order
custom_colors <- c("SJLIFE-Radiation" = "#EE6A50", "CCSS-Radiation" = "darkblue", "SJLIFE-Chemotherapy" = "red", "CCSS-Chemotherapy" = "blue")
legend_order <- c("SJLIFE-Radiation", "CCSS-Radiation", "SJLIFE-Chemotherapy", "CCSS-Chemotherapy")
# Order the levels of the legend_group factor
data_melted$legend_group <- factor(data_melted$legend_group, levels = legend_order)

# Create a factor with the desired order
data_melted$SN_types <- factor(data_melted$SN_types, levels = order)

# Create the grouped bar chart
# Create the plot
p <- ggplot(data_melted, aes(x = SN_types, y = new_value, fill = legend_group)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.8) +
  # Customize the theme and appearance
  # geom_hline(yintercept = 0, color = "black", linetype = "solid", size = 0.7) +
  # Add Y lines
  # geom_vline(xintercept = 0.5, color = "black", linetype = "solid", size = 0.7) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 16, color = "black"),
        axis.text.y = element_text(size = 16, color = "black"),
        axis.title.y = element_text(size = 20, color = "black"), 
        axis.line.x = element_line(color="black"), 
        axis.line.y = element_line(color="black"),
        legend.title = element_blank(),
        legend.text = element_text(size = 16, color = "black"),
        plot.title = element_text(size = 16, hjust = 0.5, vjust = 1.5, face = "bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  
  # Customize colors
  scale_fill_manual(values = custom_colors, breaks = legend_order) +
  
  # Add labels with geom_text
  geom_text(
    data = data_melted,
    aes(label = new_value, y = new_value),  # Adjust y position
    position = position_dodge(width = 0.8),
    vjust = -0.20,  # Adjust vertical justification
    hjust = 0.5,  # Center text horizontally
    size = 5.5,
    color = "black"
  ) +
  # Adjust legend position
  theme(legend.position = "top", legend.box = "horizontal") +
  # Additional adjustments
  scale_x_discrete(limits = order) +
  coord_cartesian(clip = "off")
p
p <- p + scale_y_continuous(limits = c(0, 100), expand = c(0, 0)) + 
  labs(title = "", y = "Attributable fraction (%)", x =  NULL) +
  theme(axis.title.x = element_text(size = 18, color = "black")) +
  theme(axis.title.y = element_text(size = 18, color = "black"))

p

# Save the plot as a high-resolution image
plot_name <- paste0("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/ANALYSIS/results/plots/v20b/overall/", group, "_", AF.type,"_", lifestyle, ".tiff")

ggsave(plot_name, p, width = 12, height = 8, dpi = 600, device = "tiff", compression = "lzw")


## 2.___________________________________________________________________________ Female analysis
## 2.___________________________________________________________________________ Female analysis
## 2.___________________________________________________________________________ Female analysis
## 2.___________________________________________________________________________ Female analysis
## 2.___________________________________________________________________________ Female analysis


data <- saved.data

group <- "Female"
variables <- unique(data$Variables)

custom_colors <- c("SJLIFE" = "#EE6A50", "CCSS" = "darkblue")
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
  
  
  p <- ggplot(data_melted, aes(x = SN_types, y = new_value, fill = legend_group)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.4), width = 0.4) +
    # Customize the theme and appearance
    # # geom_hline(yintercept = seq(-0.1, 0.1, by = 0.1), color = "black", linetype = "solid", size = 0.7) +
    # geom_hline(yintercept = 0, color = "black", linetype = "solid", size = 0.7) +
    # Add Y lines
    # # geom_vline(xintercept = seq(0.5, length(unique(data_melted$SN_types)) + 0.5), color = "black", linetype = "solid", size = 0.7) +
    # geom_vline(xintercept = 0.5, color = "black", linetype = "solid", size = 0.7) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 16, color = "black"),
          axis.text.y = element_text(size = 16, color = "black"),
          axis.title.y = element_text(size = 20, color = "black"), 
          axis.line.x = element_line(color="black"), 
          axis.line.y = element_line(color="black"),
          legend.title = element_blank(),
          legend.text = element_text(size = 16, color = "black"),
          plot.title = element_text(size = 16, hjust = 0.5, vjust = 1.5, face = "bold"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    
    # Customize colors
    scale_fill_manual(values = custom_colors, breaks = legend_order) +
    
    # Add labels with geom_text
    geom_text(
      data = data_melted,
      aes(label = new_value, y = new_value),  # Adjust y position
      position = position_dodge(width = 0.4),
      vjust = -0.20,  # Adjust vertical justification
      hjust = 0.5,  # Center text horizontally
      size = 5.5,
      color = "black"
    ) +
    # Adjust legend position
    theme(legend.position = "top", legend.box = "horizontal") +
    # Additional adjustments
    scale_x_discrete(limits = order) +
    coord_cartesian(clip = "off")
  p <- p + scale_y_continuous(limits = c(0, 100), expand = c(0, 0)) + labs(title = "", y = "Attributable fraction (%)", x = NULL) 
  
  p
  # Save the plot as a high-resolution image
  plot_name <- paste0("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/ANALYSIS/results/plots/v20b/female/", group, "_", AF.type, "_", lifestyle, ".tiff")
  
  ggsave(plot_name, p, width = 12, height = 8, dpi = 600, device = "tiff", compression = "lzw")
}

################################
## PRS and treatment together ##
################################
group <- "Female"

variables <- "treatment|PRS"

new_data <- data[grepl(variables, data$Variables), c("Cohort", "SN_types", "Variables", group)]
colnames(new_data) <- c("variable", "SN_types", "AF_by", "value")
new_data$value <- as.numeric(new_data$value)

# Reshape the data for ggplot2
library(reshape2)
data_melted <-new_data
data_melted$new_value <- round(data_melted$value,2)*100

order <- unique(data_melted$SN_types)
AF.type <- "treatment_PRS"
lifestyle <- "without_lifestyle"

data_melted$legend_group <- gsub("_", " ", factor(paste(data_melted$variable, data_melted$AF_by, sep = "-")))


# Define custom colors and legend order
custom_colors <- c("SJLIFE-All treatments" = "#EE6A50", "CCSS-All treatments" = "darkblue", "SJLIFE-PRS" = "red", "CCSS-PRS" = "blue")
legend_order <- c("SJLIFE-All treatments", "CCSS-All treatments", "SJLIFE-PRS", "CCSS-PRS")
# Order the levels of the legend_group factor
data_melted$legend_group <- factor(data_melted$legend_group, levels = legend_order)

# Create a factor with the desired order
data_melted$SN_types <- factor(data_melted$SN_types, levels = order)

# Create the grouped bar chart
# Create the plot
p <- ggplot(data_melted, aes(x = SN_types, y = new_value, fill = legend_group)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.8) +
  # Customize the theme and appearance
  # geom_hline(yintercept = 0, color = "black", linetype = "solid", size = 0.7) +
  # Add Y lines
  # geom_vline(xintercept = 0.5, color = "black", linetype = "solid", size = 0.7) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 16, color = "black"),
        axis.text.y = element_text(size = 16, color = "black"),
        axis.title.y = element_text(size = 20, color = "black"), 
        axis.line.x = element_line(color="black"), 
        axis.line.y = element_line(color="black"),
        legend.title = element_blank(),
        legend.text = element_text(size = 16, color = "black"),
        plot.title = element_text(size = 16, hjust = 0.5, vjust = 1.5, face = "bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  
  # Customize colors
  scale_fill_manual(values = custom_colors, breaks = legend_order) +
  
  # Add labels with geom_text
  geom_text(
    data = data_melted,
    aes(label = new_value, y = new_value),  # Adjust y position
    position = position_dodge(width = 0.8),
    vjust = -0.20,  # Adjust vertical justification
    hjust = 0.5,  # Center text horizontally
    size = 5.5,
    color = "black"
  ) +
  # Adjust legend position
  theme(legend.position = "top", legend.box = "horizontal") +
  # Additional adjustments
  scale_x_discrete(limits = order) +
  coord_cartesian(clip = "off")
p
p <- p + scale_y_continuous(limits = c(0, 100), expand = c(0, 0)) + 
  labs(title = "", y = "Attributable fraction (%)", x =  NULL) +
  theme(axis.title.x = element_text(size = 18, color = "black")) +
  theme(axis.title.y = element_text(size = 18, color = "black"))

p

# Save the plot as a high-resolution image
plot_name <- paste0("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/ANALYSIS/results/plots/v20b/female/", group, "_", AF.type,"_", lifestyle, ".tiff")

ggsave(plot_name, p, width = 12, height = 8, dpi = 600, device = "tiff", compression = "lzw")



##################################
## Chemo and Radiation together ##
##################################
group <- "Female"

variables <- "Radiation|Chemotherapy"

new_data <- data[grepl(variables, data$Variables), c("Cohort", "SN_types", "Variables", group)]
colnames(new_data) <- c("variable", "SN_types", "AF_by", "value")
new_data$value <- as.numeric(new_data$value)

# Reshape the data for ggplot2
library(reshape2)
data_melted <-new_data
data_melted$new_value <- round(data_melted$value,2)*100

order <- unique(data_melted$SN_types)
AF.type <- "Chemo_Radiation"
lifestyle <- "without_lifestyle"

data_melted$legend_group <- factor(paste(data_melted$variable, data_melted$AF_by, sep = "-"))


# Define custom colors and legend order
custom_colors <- c("SJLIFE-Radiation" = "#EE6A50", "CCSS-Radiation" = "darkblue", "SJLIFE-Chemotherapy" = "red", "CCSS-Chemotherapy" = "blue")
legend_order <- c("SJLIFE-Radiation", "CCSS-Radiation", "SJLIFE-Chemotherapy", "CCSS-Chemotherapy")
# Order the levels of the legend_group factor
data_melted$legend_group <- factor(data_melted$legend_group, levels = legend_order)

# Create a factor with the desired order
data_melted$SN_types <- factor(data_melted$SN_types, levels = order)

# Create the grouped bar chart
# Create the plot
p <- ggplot(data_melted, aes(x = SN_types, y = new_value, fill = legend_group)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.8) +
  # Customize the theme and appearance
  # geom_hline(yintercept = 0, color = "black", linetype = "solid", size = 0.7) +
  # Add Y lines
  # geom_vline(xintercept = 0.5, color = "black", linetype = "solid", size = 0.7) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 16, color = "black"),
        axis.text.y = element_text(size = 16, color = "black"),
        axis.title.y = element_text(size = 20, color = "black"), 
        axis.line.x = element_line(color="black"), 
        axis.line.y = element_line(color="black"),
        legend.title = element_blank(),
        legend.text = element_text(size = 16, color = "black"),
        plot.title = element_text(size = 16, hjust = 0.5, vjust = 1.5, face = "bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  
  # Customize colors
  scale_fill_manual(values = custom_colors, breaks = legend_order) +
  
  # Add labels with geom_text
  geom_text(
    data = data_melted,
    aes(label = new_value, y = new_value),  # Adjust y position
    position = position_dodge(width = 0.8),
    vjust = -0.20,  # Adjust vertical justification
    hjust = 0.5,  # Center text horizontally
    size = 5.5,
    color = "black"
  ) +
  # Adjust legend position
  theme(legend.position = "top", legend.box = "horizontal") +
  # Additional adjustments
  scale_x_discrete(limits = order) +
  coord_cartesian(clip = "off")
p
p <- p + scale_y_continuous(limits = c(0, 100), expand = c(0, 0)) + 
  labs(title = "", y = "Attributable fraction (%)", x =  NULL) +
  theme(axis.title.x = element_text(size = 18, color = "black")) +
  theme(axis.title.y = element_text(size = 18, color = "black"))

p

# Save the plot as a high-resolution image
plot_name <- paste0("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/ANALYSIS/results/plots/v20b/female/", group, "_", AF.type,"_", lifestyle, ".tiff")

ggsave(plot_name, p, width = 12, height = 8, dpi = 600, device = "tiff", compression = "lzw")


## 3.___________________________________________________________________________ Male analysis
## 3.___________________________________________________________________________ Male analysis
## 3.___________________________________________________________________________ Male analysis
## 3.___________________________________________________________________________ Male analysis
## 3.___________________________________________________________________________ Male analysis


data <- saved.data
data <- data[!data$SN_types %in% "Breast cancer",]

group <- "Male"
variables <- unique(data$Variables)

custom_colors <- c("SJLIFE" = "#EE6A50", "CCSS" = "darkblue")
legend_order <- c("SJLIFE", "CCSS")

for(i in 1:length(variables)){
  data_melted <- data[grepl(variables[i], data$Variables), c("Cohort", "SN_types", "Variables", group)]
  colnames(data_melted) <- c("variable", "SN_types", "AF_by", "value")
  data_melted$value <- as.numeric(data_melted$value); data_melted <- data_melted[complete.cases(data_melted),]
  data_melted$new_value <- round(data_melted$value,2)*100
  
  # Reshape the data for ggplot2
  library(reshape2)
  
  data_melted <- data_melted[complete.cases(data_melted),]
  
  order <- unique(data_melted$SN_types)
  AF.type <- data_melted$AF_by[i] 
  lifestyle <- "without_lifestyle"
  if(lifestyle== "without_lifestyle" && variables[i] == "Lifestyle"){
    next
  }
  
  data_melted$legend_group <- factor(data_melted$variable, levels = c("SJLIFE", "CCSS"))
  
  # Create a factor with the desired order
  data_melted$SN_types <- factor(data_melted$SN_types, levels = order)
  
  
  
  p <- ggplot(data_melted, aes(x = SN_types, y = new_value, fill = legend_group)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.4), width = 0.4) +
    # Customize the theme and appearance
    # # geom_hline(yintercept = seq(-0.1, 0.1, by = 0.1), color = "black", linetype = "solid", size = 0.7) +
    # geom_hline(yintercept = 0, color = "black", linetype = "solid", size = 0.7) +
    # Add Y lines
    # # geom_vline(xintercept = seq(0.5, length(unique(data_melted$SN_types)) + 0.5), color = "black", linetype = "solid", size = 0.7) +
    # geom_vline(xintercept = 0.5, color = "black", linetype = "solid", size = 0.7) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 16, color = "black"),
          axis.text.y = element_text(size = 16, color = "black"),
          axis.title.y = element_text(size = 20, color = "black"), 
          axis.line.x = element_line(color="black"), 
          axis.line.y = element_line(color="black"),
          legend.title = element_blank(),
          legend.text = element_text(size = 16, color = "black"),
          plot.title = element_text(size = 16, hjust = 0.5, vjust = 1.5, face = "bold"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    
    # Customize colors
    scale_fill_manual(values = custom_colors, breaks = legend_order) +
    
    # Add labels with geom_text
    geom_text(
      data = data_melted,
      aes(label = new_value, y = new_value),  # Adjust y position
      position = position_dodge(width = 0.4),
      vjust = -0.20,  # Adjust vertical justification
      hjust = 0.5,  # Center text horizontally
      size = 5.5,
      color = "black"
    ) +
    # Adjust legend position
    theme(legend.position = "top", legend.box = "horizontal") +
    # Additional adjustments
    scale_x_discrete(limits = order) +
    coord_cartesian(clip = "off")
  p <- p + scale_y_continuous(limits = c(0, 100), expand = c(0, 0)) + labs(title = "", y = "Attributable fraction (%)", x = NULL) 
  
  p
  # Save the plot as a high-resolution image
  plot_name <- paste0("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/ANALYSIS/results/plots/v20b/male/", group, "_", AF.type, "_", lifestyle, ".tiff")
  
  ggsave(plot_name, p, width = 12, height = 8, dpi = 600, device = "tiff", compression = "lzw")
}

################################
## PRS and treatment together ##
################################
group <- "Male"

variables <- "treatment|PRS"

new_data <- data[grepl(variables, data$Variables), c("Cohort", "SN_types", "Variables", group)]
colnames(new_data) <- c("variable", "SN_types", "AF_by", "value")
new_data$value <- as.numeric(new_data$value)

# Reshape the data for ggplot2
library(reshape2)
data_melted <-new_data
data_melted$new_value <- round(data_melted$value,2)*100

order <- unique(data_melted$SN_types)
AF.type <- "treatment_PRS"
lifestyle <- "without_lifestyle"

data_melted$legend_group <- gsub("_", " ", factor(paste(data_melted$variable, data_melted$AF_by, sep = "-")))


# Define custom colors and legend order
custom_colors <- c("SJLIFE-All treatments" = "#EE6A50", "CCSS-All treatments" = "darkblue", "SJLIFE-PRS" = "red", "CCSS-PRS" = "blue")
legend_order <- c("SJLIFE-All treatments", "CCSS-All treatments", "SJLIFE-PRS", "CCSS-PRS")
# Order the levels of the legend_group factor
data_melted$legend_group <- factor(data_melted$legend_group, levels = legend_order)

# Create a factor with the desired order
data_melted$SN_types <- factor(data_melted$SN_types, levels = order)

# Create the grouped bar chart
# Create the plot
p <- ggplot(data_melted, aes(x = SN_types, y = new_value, fill = legend_group)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.8) +
  # Customize the theme and appearance
  # geom_hline(yintercept = 0, color = "black", linetype = "solid", size = 0.7) +
  # Add Y lines
  # geom_vline(xintercept = 0.5, color = "black", linetype = "solid", size = 0.7) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 16, color = "black"),
        axis.text.y = element_text(size = 16, color = "black"),
        axis.title.y = element_text(size = 20, color = "black"), 
        axis.line.x = element_line(color="black"), 
        axis.line.y = element_line(color="black"),
        legend.title = element_blank(),
        legend.text = element_text(size = 16, color = "black"),
        plot.title = element_text(size = 16, hjust = 0.5, vjust = 1.5, face = "bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  
  # Customize colors
  scale_fill_manual(values = custom_colors, breaks = legend_order) +
  
  # Add labels with geom_text
  geom_text(
    data = data_melted,
    aes(label = new_value, y = new_value),  # Adjust y position
    position = position_dodge(width = 0.8),
    vjust = -0.20,  # Adjust vertical justification
    hjust = 0.5,  # Center text horizontally
    size = 5.5,
    color = "black"
  ) +
  # Adjust legend position
  theme(legend.position = "top", legend.box = "horizontal") +
  # Additional adjustments
  scale_x_discrete(limits = order) +
  coord_cartesian(clip = "off")
p
p <- p + scale_y_continuous(limits = c(0, 100), expand = c(0, 0)) + 
  labs(title = "", y = "Attributable fraction (%)", x =  NULL) +
  theme(axis.title.x = element_text(size = 18, color = "black")) +
  theme(axis.title.y = element_text(size = 18, color = "black"))

p

# Save the plot as a high-resolution image
plot_name <- paste0("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/ANALYSIS/results/plots/v20b/male/", group, "_", AF.type,"_", lifestyle, ".tiff")

ggsave(plot_name, p, width = 12, height = 8, dpi = 600, device = "tiff", compression = "lzw")

##################################
## Chemo and Radiation together ##
##################################
group <- "Male"

variables <- "Radiation|Chemotherapy"

new_data <- data[grepl(variables, data$Variables), c("Cohort", "SN_types", "Variables", group)]
colnames(new_data) <- c("variable", "SN_types", "AF_by", "value")
new_data$value <- as.numeric(new_data$value)

# Reshape the data for ggplot2
library(reshape2)
data_melted <-new_data
data_melted$new_value <- round(data_melted$value,2)*100

order <- unique(data_melted$SN_types)
AF.type <- "Chemo_Radiation"
lifestyle <- "without_lifestyle"

data_melted$legend_group <- factor(paste(data_melted$variable, data_melted$AF_by, sep = "-"))


# Define custom colors and legend order
custom_colors <- c("SJLIFE-Radiation" = "#EE6A50", "CCSS-Radiation" = "darkblue", "SJLIFE-Chemotherapy" = "red", "CCSS-Chemotherapy" = "blue")
legend_order <- c("SJLIFE-Radiation", "CCSS-Radiation", "SJLIFE-Chemotherapy", "CCSS-Chemotherapy")
# Order the levels of the legend_group factor
data_melted$legend_group <- factor(data_melted$legend_group, levels = legend_order)

# Create a factor with the desired order
data_melted$SN_types <- factor(data_melted$SN_types, levels = order)

# Create the grouped bar chart
# Create the plot
p <- ggplot(data_melted, aes(x = SN_types, y = new_value, fill = legend_group)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.8) +
  # Customize the theme and appearance
  # geom_hline(yintercept = 0, color = "black", linetype = "solid", size = 0.7) +
  # Add Y lines
  # geom_vline(xintercept = 0.5, color = "black", linetype = "solid", size = 0.7) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 16, color = "black"),
        axis.text.y = element_text(size = 16, color = "black"),
        axis.title.y = element_text(size = 20, color = "black"), 
        axis.line.x = element_line(color="black"), 
        axis.line.y = element_line(color="black"),
        legend.title = element_blank(),
        legend.text = element_text(size = 16, color = "black"),
        plot.title = element_text(size = 16, hjust = 0.5, vjust = 1.5, face = "bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  
  # Customize colors
  scale_fill_manual(values = custom_colors, breaks = legend_order) +
  
  # Add labels with geom_text
  geom_text(
    data = data_melted,
    aes(label = new_value, y = new_value),  # Adjust y position
    position = position_dodge(width = 0.8),
    vjust = -0.20,  # Adjust vertical justification
    hjust = 0.5,  # Center text horizontally
    size = 5.5,
    color = "black"
  ) +
  # Adjust legend position
  theme(legend.position = "top", legend.box = "horizontal") +
  # Additional adjustments
  scale_x_discrete(limits = order) +
  coord_cartesian(clip = "off")
p
p <- p + scale_y_continuous(limits = c(0, 100), expand = c(0, 0)) + 
  labs(title = "", y = "Attributable fraction (%)", x =  NULL) +
  theme(axis.title.x = element_text(size = 18, color = "black")) +
  theme(axis.title.y = element_text(size = 18, color = "black"))

p
# Save the plot as a high-resolution image
plot_name <- paste0("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/ANALYSIS/results/plots/v20b/male/", group, "_", AF.type,"_", lifestyle, ".tiff")

ggsave(plot_name, p, width = 12, height = 8, dpi = 600, device = "tiff", compression = "lzw")

## 4.___________________________________________________________________________ lt35 analysis
## 4.___________________________________________________________________________ lt35 analysis
## 4.___________________________________________________________________________ lt35 analysis
## 4.___________________________________________________________________________ lt35 analysis
## 4.___________________________________________________________________________ lt35 analysis


data <- saved.data

group <- "Age.lt.35"
variables <- unique(data$Variables)

custom_colors <- c("SJLIFE" = "#EE6A50", "CCSS" = "darkblue")
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
  
  
  p <- ggplot(data_melted, aes(x = SN_types, y = new_value, fill = legend_group)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.4), width = 0.4) +
    # Customize the theme and appearance
    # # geom_hline(yintercept = seq(-0.1, 0.1, by = 0.1), color = "black", linetype = "solid", size = 0.7) +
    # geom_hline(yintercept = 0, color = "black", linetype = "solid", size = 0.7) +
    # Add Y lines
    # # geom_vline(xintercept = seq(0.5, length(unique(data_melted$SN_types)) + 0.5), color = "black", linetype = "solid", size = 0.7) +
    # geom_vline(xintercept = 0.5, color = "black", linetype = "solid", size = 0.7) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 16, color = "black"),
          axis.text.y = element_text(size = 16, color = "black"),
          axis.title.y = element_text(size = 20, color = "black"), 
          axis.line.x = element_line(color="black"), 
          axis.line.y = element_line(color="black"),
          legend.title = element_blank(),
          legend.text = element_text(size = 16, color = "black"),
          plot.title = element_text(size = 16, hjust = 0.5, vjust = 1.5, face = "bold"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    
    # Customize colors
    scale_fill_manual(values = custom_colors, breaks = legend_order) +
    
    # Add labels with geom_text
    geom_text(
      data = data_melted,
      aes(label = new_value, y = new_value),  # Adjust y position
      position = position_dodge(width = 0.4),
      vjust = -0.20,  # Adjust vertical justification
      hjust = 0.5,  # Center text horizontally
      size = 5.5,
      color = "black"
    ) +
    # Adjust legend position
    theme(legend.position = "top", legend.box = "horizontal") +
    # Additional adjustments
    scale_x_discrete(limits = order) +
    coord_cartesian(clip = "off")
  p <- p + scale_y_continuous(limits = c(0, 100), expand = c(0, 0)) + labs(title = "", y = "Attributable fraction (%)", x = NULL) 
  
  p
  # Save the plot as a high-resolution image
  plot_name <- paste0("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/ANALYSIS/results/plots/v20b/lt35/", group, "_", AF.type, "_", lifestyle, ".tiff")
  
  ggsave(plot_name, p, width = 12, height = 8, dpi = 600, device = "tiff", compression = "lzw")
}

################################
## PRS and treatment together ##
################################
group <- "Age.lt.35"

variables <- "treatment|PRS"

new_data <- data[grepl(variables, data$Variables), c("Cohort", "SN_types", "Variables", group)]
colnames(new_data) <- c("variable", "SN_types", "AF_by", "value")
new_data$value <- as.numeric(new_data$value)

# Reshape the data for ggplot2
library(reshape2)
data_melted <-new_data
data_melted$new_value <- round(data_melted$value,2)*100

order <- unique(data_melted$SN_types)
AF.type <- "treatment_PRS"
lifestyle <- "without_lifestyle"

data_melted$legend_group <- gsub("_", " ", factor(paste(data_melted$variable, data_melted$AF_by, sep = "-")))


# Define custom colors and legend order
custom_colors <- c("SJLIFE-All treatments" = "#EE6A50", "CCSS-All treatments" = "darkblue", "SJLIFE-PRS" = "red", "CCSS-PRS" = "blue")
legend_order <- c("SJLIFE-All treatments", "CCSS-All treatments", "SJLIFE-PRS", "CCSS-PRS")
# Order the levels of the legend_group factor
data_melted$legend_group <- factor(data_melted$legend_group, levels = legend_order)

# Create a factor with the desired order
data_melted$SN_types <- factor(data_melted$SN_types, levels = order)

# Create the grouped bar chart
# Create the plot
p <- ggplot(data_melted, aes(x = SN_types, y = new_value, fill = legend_group)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.8) +
  # Customize the theme and appearance
  # geom_hline(yintercept = 0, color = "black", linetype = "solid", size = 0.7) +
  # Add Y lines
  # geom_vline(xintercept = 0.5, color = "black", linetype = "solid", size = 0.7) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 16, color = "black"),
        axis.text.y = element_text(size = 16, color = "black"),
        axis.title.y = element_text(size = 20, color = "black"), 
        axis.line.x = element_line(color="black"), 
        axis.line.y = element_line(color="black"),
        legend.title = element_blank(),
        legend.text = element_text(size = 16, color = "black"),
        plot.title = element_text(size = 16, hjust = 0.5, vjust = 1.5, face = "bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  
  # Customize colors
  scale_fill_manual(values = custom_colors, breaks = legend_order) +
  
  # Add labels with geom_text
  geom_text(
    data = data_melted,
    aes(label = new_value, y = new_value),  # Adjust y position
    position = position_dodge(width = 0.8),
    vjust = -0.20,  # Adjust vertical justification
    hjust = 0.5,  # Center text horizontally
    size = 5.5,
    color = "black"
  ) +
  # Adjust legend position
  theme(legend.position = "top", legend.box = "horizontal") +
  # Additional adjustments
  scale_x_discrete(limits = order) +
  coord_cartesian(clip = "off")
p
p <- p + scale_y_continuous(limits = c(0, 100), expand = c(0, 0)) + 
  labs(title = "", y = "Attributable fraction (%)", x =  NULL) +
  theme(axis.title.x = element_text(size = 18, color = "black")) +
  theme(axis.title.y = element_text(size = 18, color = "black"))

p
# Save the plot as a high-resolution image
plot_name <- paste0("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/ANALYSIS/results/plots/v20b/lt35/", group, "_", AF.type,"_", lifestyle, ".tiff")

ggsave(plot_name, p, width = 12, height = 8, dpi = 600, device = "tiff", compression = "lzw")

##################################
## Chemo and Radiation together ##
##################################
group <- "Age.lt.35"

variables <- "Radiation|Chemotherapy"

new_data <- data[grepl(variables, data$Variables), c("Cohort", "SN_types", "Variables", group)]
colnames(new_data) <- c("variable", "SN_types", "AF_by", "value")
new_data$value <- as.numeric(new_data$value)

# Reshape the data for ggplot2
library(reshape2)
data_melted <-new_data
data_melted$new_value <- round(data_melted$value,2)*100

order <- unique(data_melted$SN_types)
AF.type <- "Chemo_Radiation"
lifestyle <- "without_lifestyle"

data_melted$legend_group <- factor(paste(data_melted$variable, data_melted$AF_by, sep = "-"))


# Define custom colors and legend order
custom_colors <- c("SJLIFE-Radiation" = "#EE6A50", "CCSS-Radiation" = "darkblue", "SJLIFE-Chemotherapy" = "red", "CCSS-Chemotherapy" = "blue")
legend_order <- c("SJLIFE-Radiation", "CCSS-Radiation", "SJLIFE-Chemotherapy", "CCSS-Chemotherapy")
# Order the levels of the legend_group factor
data_melted$legend_group <- factor(data_melted$legend_group, levels = legend_order)

# Create a factor with the desired order
data_melted$SN_types <- factor(data_melted$SN_types, levels = order)

# Create the grouped bar chart
# Create the plot
p <- ggplot(data_melted, aes(x = SN_types, y = new_value, fill = legend_group)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.8) +
  # Customize the theme and appearance
  # geom_hline(yintercept = 0, color = "black", linetype = "solid", size = 0.7) +
  # Add Y lines
  # geom_vline(xintercept = 0.5, color = "black", linetype = "solid", size = 0.7) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 16, color = "black"),
        axis.text.y = element_text(size = 16, color = "black"),
        axis.title.y = element_text(size = 20, color = "black"), 
        axis.line.x = element_line(color="black"), 
        axis.line.y = element_line(color="black"),
        legend.title = element_blank(),
        legend.text = element_text(size = 16, color = "black"),
        plot.title = element_text(size = 16, hjust = 0.5, vjust = 1.5, face = "bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  
  # Customize colors
  scale_fill_manual(values = custom_colors, breaks = legend_order) +
  
  # Add labels with geom_text
  geom_text(
    data = data_melted,
    aes(label = new_value, y = new_value),  # Adjust y position
    position = position_dodge(width = 0.8),
    vjust = -0.20,  # Adjust vertical justification
    hjust = 0.5,  # Center text horizontally
    size = 5.5,
    color = "black"
  ) +
  # Adjust legend position
  theme(legend.position = "top", legend.box = "horizontal") +
  # Additional adjustments
  scale_x_discrete(limits = order) +
  coord_cartesian(clip = "off")
p
p <- p + scale_y_continuous(limits = c(0, 100), expand = c(0, 0)) + 
  labs(title = "", y = "Attributable fraction (%)", x =  NULL) +
  theme(axis.title.x = element_text(size = 18, color = "black")) +
  theme(axis.title.y = element_text(size = 18, color = "black"))

p
# Save the plot as a high-resolution image
plot_name <- paste0("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/ANALYSIS/results/plots/v20b/lt35/", group, "_", AF.type,"_", lifestyle, ".tiff")

ggsave(plot_name, p, width = 12, height = 8, dpi = 600, device = "tiff", compression = "lzw")

## 5.___________________________________________________________________________ ge35 analysis
## 5.___________________________________________________________________________ ge35 analysis
## 5.___________________________________________________________________________ ge35 analysis
## 5.___________________________________________________________________________ ge35 analysis
## 5.___________________________________________________________________________ ge35 analysis

data <- saved.data

group <- "Age.ge.35"
variables <- unique(data$Variables)

custom_colors <- c("SJLIFE" = "#EE6A50", "CCSS" = "darkblue")
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
  
  
  p <- ggplot(data_melted, aes(x = SN_types, y = new_value, fill = legend_group)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.4), width = 0.4) +
    # Customize the theme and appearance
    # # geom_hline(yintercept = seq(-0.1, 0.1, by = 0.1), color = "black", linetype = "solid", size = 0.7) +
    # geom_hline(yintercept = 0, color = "black", linetype = "solid", size = 0.7) +
    # Add Y lines
    # # geom_vline(xintercept = seq(0.5, length(unique(data_melted$SN_types)) + 0.5), color = "black", linetype = "solid", size = 0.7) +
    # geom_vline(xintercept = 0.5, color = "black", linetype = "solid", size = 0.7) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 16, color = "black"),
          axis.text.y = element_text(size = 16, color = "black"),
          axis.title.y = element_text(size = 20, color = "black"), 
          axis.line.x = element_line(color="black"), 
          axis.line.y = element_line(color="black"),
          legend.title = element_blank(),
          legend.text = element_text(size = 16, color = "black"),
          plot.title = element_text(size = 16, hjust = 0.5, vjust = 1.5, face = "bold"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    
    # Customize colors
    scale_fill_manual(values = custom_colors, breaks = legend_order) +
    
    # Add labels with geom_text
    geom_text(
      data = data_melted,
      aes(label = new_value, y = new_value),  # Adjust y position
      position = position_dodge(width = 0.4),
      vjust = -0.20,  # Adjust vertical justification
      hjust = 0.5,  # Center text horizontally
      size = 5.5,
      color = "black"
    ) +
    # Adjust legend position
    theme(legend.position = "top", legend.box = "horizontal") +
    # Additional adjustments
    scale_x_discrete(limits = order) +
    coord_cartesian(clip = "off")
  p <- p + scale_y_continuous(limits = c(0, 100), expand = c(0, 0)) + labs(title = "", y = "Attributable fraction (%)", x = NULL) 
  
  p
  # Save the plot as a high-resolution image
  plot_name <- paste0("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/ANALYSIS/results/plots/v20b/ge35/", group, "_", AF.type, "_", lifestyle, ".tiff")
  
  ggsave(plot_name, p, width = 12, height = 8, dpi = 600, device = "tiff", compression = "lzw")
}

################################
## PRS and treatment together ##
################################
group <- "Age.ge.35"

variables <- "treatment|PRS"

new_data <- data[grepl(variables, data$Variables), c("Cohort", "SN_types", "Variables", group)]
colnames(new_data) <- c("variable", "SN_types", "AF_by", "value")
new_data$value <- as.numeric(new_data$value)

# Reshape the data for ggplot2
library(reshape2)
data_melted <-new_data
data_melted$new_value <- round(data_melted$value,2)*100

order <- unique(data_melted$SN_types)
AF.type <- "treatment_PRS"
lifestyle <- "without_lifestyle"

data_melted$legend_group <- gsub("_", " ", factor(paste(data_melted$variable, data_melted$AF_by, sep = "-")))


# Define custom colors and legend order
custom_colors <- c("SJLIFE-All treatments" = "#EE6A50", "CCSS-All treatments" = "darkblue", "SJLIFE-PRS" = "red", "CCSS-PRS" = "blue")
legend_order <- c("SJLIFE-All treatments", "CCSS-All treatments", "SJLIFE-PRS", "CCSS-PRS")
# Order the levels of the legend_group factor
data_melted$legend_group <- factor(data_melted$legend_group, levels = legend_order)

# Create a factor with the desired order
data_melted$SN_types <- factor(data_melted$SN_types, levels = order)

# Create the grouped bar chart
# Create the plot
p <- ggplot(data_melted, aes(x = SN_types, y = new_value, fill = legend_group)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.8) +
  # Customize the theme and appearance
  # geom_hline(yintercept = 0, color = "black", linetype = "solid", size = 0.7) +
  # Add Y lines
  # geom_vline(xintercept = 0.5, color = "black", linetype = "solid", size = 0.7) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 16, color = "black"),
        axis.text.y = element_text(size = 16, color = "black"),
        axis.title.y = element_text(size = 20, color = "black"), 
        axis.line.x = element_line(color="black"), 
        axis.line.y = element_line(color="black"),
        legend.title = element_blank(),
        legend.text = element_text(size = 16, color = "black"),
        plot.title = element_text(size = 16, hjust = 0.5, vjust = 1.5, face = "bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  
  # Customize colors
  scale_fill_manual(values = custom_colors, breaks = legend_order) +
  
  # Add labels with geom_text
  geom_text(
    data = data_melted,
    aes(label = new_value, y = new_value),  # Adjust y position
    position = position_dodge(width = 0.8),
    vjust = -0.20,  # Adjust vertical justification
    hjust = 0.5,  # Center text horizontally
    size = 5.5,
    color = "black"
  ) +
  # Adjust legend position
  theme(legend.position = "top", legend.box = "horizontal") +
  # Additional adjustments
  scale_x_discrete(limits = order) +
  coord_cartesian(clip = "off")
p
p <- p + scale_y_continuous(limits = c(0, 100), expand = c(0, 0)) + 
  labs(title = "", y = "Attributable fraction (%)", x =  NULL) +
  theme(axis.title.x = element_text(size = 18, color = "black")) +
  theme(axis.title.y = element_text(size = 18, color = "black"))

p
# Save the plot as a high-resolution image
plot_name <- paste0("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/ANALYSIS/results/plots/v20b/ge35/", group, "_", AF.type,"_", lifestyle, ".tiff")

ggsave(plot_name, p, width = 12, height = 8, dpi = 600, device = "tiff", compression = "lzw")

##################################
## Chemo and Radiation together ##
##################################
group <- "Age.ge.35"

variables <- "Radiation|Chemotherapy"

new_data <- data[grepl(variables, data$Variables), c("Cohort", "SN_types", "Variables", group)]
colnames(new_data) <- c("variable", "SN_types", "AF_by", "value")
new_data$value <- as.numeric(new_data$value)

# Reshape the data for ggplot2
library(reshape2)
data_melted <-new_data
data_melted$new_value <- round(data_melted$value,2)*100

order <- unique(data_melted$SN_types)
AF.type <- "Chemo_Radiation"
lifestyle <- "without_lifestyle"

data_melted$legend_group <- factor(paste(data_melted$variable, data_melted$AF_by, sep = "-"))


# Define custom colors and legend order
custom_colors <- c("SJLIFE-Radiation" = "#EE6A50", "CCSS-Radiation" = "darkblue", "SJLIFE-Chemotherapy" = "red", "CCSS-Chemotherapy" = "blue")
legend_order <- c("SJLIFE-Radiation", "CCSS-Radiation", "SJLIFE-Chemotherapy", "CCSS-Chemotherapy")
# Order the levels of the legend_group factor
data_melted$legend_group <- factor(data_melted$legend_group, levels = legend_order)

# Create a factor with the desired order
data_melted$SN_types <- factor(data_melted$SN_types, levels = order)

# Create the grouped bar chart
# Create the plot
p <- ggplot(data_melted, aes(x = SN_types, y = new_value, fill = legend_group)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.8) +
  # Customize the theme and appearance
  # geom_hline(yintercept = 0, color = "black", linetype = "solid", size = 0.7) +
  # Add Y lines
  # geom_vline(xintercept = 0.5, color = "black", linetype = "solid", size = 0.7) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 16, color = "black"),
        axis.text.y = element_text(size = 16, color = "black"),
        axis.title.y = element_text(size = 20, color = "black"), 
        axis.line.x = element_line(color="black"), 
        axis.line.y = element_line(color="black"),
        legend.title = element_blank(),
        legend.text = element_text(size = 16, color = "black"),
        plot.title = element_text(size = 16, hjust = 0.5, vjust = 1.5, face = "bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  
  # Customize colors
  scale_fill_manual(values = custom_colors, breaks = legend_order) +
  
  # Add labels with geom_text
  geom_text(
    data = data_melted,
    aes(label = new_value, y = new_value),  # Adjust y position
    position = position_dodge(width = 0.8),
    vjust = -0.20,  # Adjust vertical justification
    hjust = 0.5,  # Center text horizontally
    size = 5.5,
    color = "black"
  ) +
  # Adjust legend position
  theme(legend.position = "top", legend.box = "horizontal") +
  # Additional adjustments
  scale_x_discrete(limits = order) +
  coord_cartesian(clip = "off")
p
p <- p + scale_y_continuous(limits = c(0, 100), expand = c(0, 0)) + 
  labs(title = "", y = "Attributable fraction (%)", x =  NULL) +
  theme(axis.title.x = element_text(size = 18, color = "black")) +
  theme(axis.title.y = element_text(size = 18, color = "black"))
p
# Save the plot as a high-resolution image
plot_name <- paste0("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/ANALYSIS/results/plots/v20b/ge35/", group, "_", AF.type,"_", lifestyle, ".tiff")

ggsave(plot_name, p, width = 12, height = 8, dpi = 600, device = "tiff", compression = "lzw")