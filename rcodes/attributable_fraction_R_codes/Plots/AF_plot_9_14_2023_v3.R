# Load required libraries
library(ggplot2)
library(ggthemes)# For additional themes
library(ggrepel)
library(dplyr)



## V18 b
data <- read.table(text="Cohort	SN_types	Variables	Overall	Female	Male	Age.lt.35	Age.ge.35
SJLIFE	Any SN (605)	Radiation	0.425	0.417	0.435	0.400	0.447
SJLIFE	Any SN (605)	Chemo	0.080	0.079	0.080	0.109	0.054
SJLIFE	Any SN (605)	All_treatments	0.469	0.462	0.478	0.463	0.474
SJLIFE	Any SN (605)	PRS	0.125	0.123	0.129	0.126	0.125
SJLIFE	Any SN (605)	Lifestyle	-	-	-	-	-
SJLIFE	Any SN (605)	Combined	0.536	0.530	0.545	0.531	0.541
SJLIFE	Any SMN (462)	Radiation	0.372	0.371	0.373	0.340	0.396
SJLIFE	Any SMN (462)	Chemo	0.033	0.033	0.032	0.047	0.022
SJLIFE	Any SMN (462)	All_treatments	0.392	0.392	0.392	0.370	0.409
SJLIFE	Any SMN (462)	PRS	0.122	0.119	0.125	0.122	0.122
SJLIFE	Any SMN (462)	Lifestyle	-	-	-	-	-
SJLIFE	Any SMN (462)	Combined	0.467	0.466	0.468	0.447	0.482
SJLIFE	NMSC (249)	Radiation	0.441	0.431	0.453	0.413	0.457
SJLIFE	NMSC (249)	Chemo	-	-	-	-	-
SJLIFE	NMSC (249)	All_treatments	0.441	0.431	0.453	0.413	0.457
SJLIFE	NMSC (249)	PRS	0.441	0.432	0.452	0.440	0.442
SJLIFE	NMSC (249)	Lifestyle	-	-	-	-	-
SJLIFE	NMSC (249)	Combined	0.688	0.680	0.698	0.671	0.698
SJLIFE	Breast cancer (74)	Radiation	0.493	0.493	-	0.475	0.498
SJLIFE	Breast cancer (74)	Chemo	0.193	0.193	-	0.228	0.184
SJLIFE	Breast cancer (74)	All_treatments	0.602	0.602	-	0.594	0.604
SJLIFE	Breast cancer (74)	PRS	0.255	0.255	-	0.261	0.254
SJLIFE	Breast cancer (74)	Lifestyle	-	-	-	-	-
SJLIFE	Breast cancer (74)	Combined	0.704	0.704	-	0.701	0.705
SJLIFE	Thyroid cancer (86)	Radiation	0.620	0.622	0.616	0.594	0.650
SJLIFE	Thyroid cancer (86)	Chemo	0.233	0.235	0.231	0.293	0.163
SJLIFE	Thyroid cancer (86)	All_treatments	0.726	0.727	0.723	0.725	0.727
SJLIFE	Thyroid cancer (86)	PRS	0.517	0.519	0.514	0.519	0.514
SJLIFE	Thyroid cancer (86)	Lifestyle	-	-	-	-	-
SJLIFE	Thyroid cancer (86)	Combined	0.866	0.867	0.866	0.866	0.868
SJLIFE	Meningioma (149)	Radiation	0.197	0.172	0.230	0.223	0.176
SJLIFE	Meningioma (149)	Chemo	0.361	0.371	0.349	0.472	0.271
SJLIFE	Meningioma (149)	All_treatments	0.489	0.478	0.504	0.598	0.400
SJLIFE	Meningioma (149)	PRS	-0.125	-0.117	-0.135	-0.121	-0.128
SJLIFE	Meningioma (149)	Lifestyle	-	-	-	-	-
SJLIFE	Meningioma (149)	Combined	0.425	0.416	0.437	0.549	0.323
SJLIFE	Sarcoma (32)	Radiation	-	-	-	-	-
SJLIFE	Sarcoma (32)	Chemo	0.324	0.329	0.320	0.308	0.351
SJLIFE	Sarcoma (32)	All_treatments	0.324	0.329	0.320	0.308	0.351
SJLIFE	Sarcoma (32)	PRS	-0.020	-0.021	-0.020	-0.020	-0.020
SJLIFE	Sarcoma (32)	Lifestyle	-	-	-	-	-
SJLIFE	Sarcoma (32)	Combined	0.311	0.315	0.306	0.294	0.337
CCSS	Any SN (1611)	Radiation	0.387	0.380	0.398	0.372	0.398
CCSS	Any SN (1611)	Chemo	0.031	0.030	0.033	0.047	0.019
CCSS	Any SN (1611)	All_treatments	0.407	0.400	0.418	0.403	0.410
CCSS	Any SN (1611)	PRS	0.048	0.047	0.049	0.047	0.048
CCSS	Any SN (1611)	Lifestyle	-	-	-	-	-
CCSS	Any SN (1611)	Combined	0.435	0.428	0.446	0.431	0.438
CCSS	Any SMN (762)	Radiation	0.255	0.254	0.257	0.213	0.283
CCSS	Any SMN (762)	Chemo	0.042	0.042	0.042	0.070	0.024
CCSS	Any SMN (762)	All_treatments	0.290	0.289	0.293	0.271	0.303
CCSS	Any SMN (762)	PRS	0.047	0.046	0.047	0.046	0.047
CCSS	Any SMN (762)	Lifestyle	-	-	-	-	-
CCSS	Any SMN (762)	Combined	0.323	0.322	0.326	0.305	0.336
CCSS	NMSC (769)	Radiation	0.391	0.386	0.396	0.385	0.394
CCSS	NMSC (769)	Chemo	-	-	-	-	-
CCSS	NMSC (769)	All_treatments	0.391	0.386	0.396	0.385	0.394
CCSS	NMSC (769)	PRS	0.304	0.306	0.302	0.306	0.303
CCSS	NMSC (769)	Lifestyle	-	-	-	-	-
CCSS	NMSC (769)	Combined	0.577	0.575	0.578	0.574	0.578
CCSS	Breast cancer (289)	Radiation	0.474	0.474	-	0.452	0.478
CCSS	Breast cancer (289)	Chemo	0.187	0.187	-	0.224	0.180
CCSS	Breast cancer (289)	All_treatments	0.593	0.593	-	0.589	0.593
CCSS	Breast cancer (289)	PRS	0.365	0.365	-	0.369	0.365
CCSS	Breast cancer (289)	Lifestyle	-	-	-	-	-
CCSS	Breast cancer (289)	Combined	0.743	0.743	-	0.741	0.743
CCSS	Thyroid cancer (163)	Radiation	0.443	0.445	0.439	0.413	0.489
CCSS	Thyroid cancer (163)	Chemo	0.055	0.055	0.056	0.073	0.029
CCSS	Thyroid cancer (163)	All_treatments	0.476	0.477	0.472	0.457	0.504
CCSS	Thyroid cancer (163)	PRS	0.358	0.362	0.352	0.355	0.363
CCSS	Thyroid cancer (163)	Lifestyle	-	-	-	-	-
CCSS	Thyroid cancer (163)	Combined	0.664	0.665	0.661	0.650	0.684
CCSS	Meningioma (255)	Radiation	0.370	0.333	0.422	0.394	0.345
CCSS	Meningioma (255)	Chemo	0.082	0.076	0.092	0.118	0.045
CCSS	Meningioma (255)	All_treatments	0.426	0.391	0.475	0.471	0.378
CCSS	Meningioma (255)	PRS	0.032	0.032	0.033	0.036	0.028
CCSS	Meningioma (255)	Lifestyle	-	-	-	-	-
CCSS	Meningioma (255)	Combined	0.444	0.410	0.494	0.490	0.397
CCSS	Sarcoma (60)	Radiation	-	-	-	-	-
CCSS	Sarcoma (60)	Chemo	0.342	0.327	0.358	0.337	0.348
CCSS	Sarcoma (60)	All_treatments	0.342	0.327	0.358	0.337	0.348
CCSS	Sarcoma (60)	PRS	0.024	0.026	0.022	0.023	0.025
CCSS	Sarcoma (60)	Lifestyle	-	-	-	-	-
CCSS	Sarcoma (60)	Combined	0.358	0.345	0.372	0.353	0.364", header = T, sep = "\t")

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

custom_colors <- c("SJLIFE" = "#1E90FF", "CCSS" = "#FF6347")
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
    position = position_dodge(width = 0.3),
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
plot_name <- paste0("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/ANALYSIS/results/plots/v18b/overall/", group, "_", AF.type, "_", lifestyle, ".tiff")

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
order <- unique(data_melted$SN_types)
AF.type <- "treatment_PRS"
lifestyle <- "without_lifestyle"

data_melted$legend_group <- factor(paste(data_melted$variable, data_melted$AF_by, sep = "_"))


data_melted$new_value <- round(data_melted$value,2)*100


# Define custom colors and legend order
custom_colors <- c("SJLIFE_PRS" = "#87CEFA", "CCSS_PRS" = "#FFA07A", "SJLIFE_All_treatments" = "#1E90FF", "CCSS_All_treatments" = "#FF6347")
legend_order <- c("SJLIFE_PRS", "CCSS_PRS", "SJLIFE_All_treatments", "CCSS_All_treatments")
# Order the levels of the legend_group factor
data_melted$legend_group <- factor(data_melted$legend_group, levels = legend_order)

# Create a factor with the desired order
data_melted$SN_types <- factor(order, levels = order)

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
    data = data_melted %>% filter(!is.na(new_value)),
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
p <- p + scale_y_continuous(limits = c(0, 100), expand = c(0, 0)) + labs(title = "", y = "Attributable fraction (%)", x = NULL) 

p

# Save the plot as a high-resolution image
plot_name <- paste0("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/ANALYSIS/results/plots/v18b/overall/", group, "_", AF.type,"_", lifestyle, ".tiff")

ggsave(plot_name, p, width = 12, height = 8, dpi = 600, device = "tiff", compression = "lzw")

##################################
## Chemo and Radiation together ##
##################################
group <- "Overall"

variables <- "Chemo|Radiation"

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

data_melted$legend_group <- factor(paste(data_melted$variable, data_melted$AF_by, sep = "_"))


# Define custom colors and legend order
custom_colors <- c("SJLIFE_Chemo" = "#87CEFA", "CCSS_Chemo" = "#FFA07A", "SJLIFE_Radiation" = "#1E90FF", "CCSS_Radiation" = "#FF6347")
legend_order <- c("SJLIFE_Chemo", "CCSS_Chemo", "SJLIFE_Radiation", "CCSS_Radiation")
# Order the levels of the legend_group factor
data_melted$legend_group <- factor(data_melted$legend_group, levels = legend_order)

# Define the order of x-axis categories
order <- c("Any SN", "Any SMN", "NMSC", "Breast cancer", "Thyroid cancer", "Meningioma", "Sarcoma")

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
p <- p + scale_y_continuous(limits = c(0, 100), expand = c(0, 0)) + labs(title = "", y = "Attributable fraction (%)", x = NULL) 

p

# Save the plot as a high-resolution image
plot_name <- paste0("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/ANALYSIS/results/plots/v18b/overall/", group, "_", AF.type,"_", lifestyle, ".tiff")

ggsave(plot_name, p, width = 12, height = 8, dpi = 600, device = "tiff", compression = "lzw")


## 2.___________________________________________________________________________ Female analysis
## 2.___________________________________________________________________________ Female analysis
## 2.___________________________________________________________________________ Female analysis
## 2.___________________________________________________________________________ Female analysis
## 2.___________________________________________________________________________ Female analysis


data <- saved.data

group <- "Female"
variables <- unique(data$Variables)

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
      position = position_dodge(width = 0.3),
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
  plot_name <- paste0("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/ANALYSIS/results/plots/v18b/female/", group, "_", AF.type, "_", lifestyle, ".tiff")
  
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


order <- unique(data_melted$SN_types)
AF.type <- "treatment_PRS"
lifestyle <- "without_lifestyle"

data_melted$legend_group <- factor(paste(data_melted$variable, data_melted$AF_by, sep = "_"))


data_melted$new_value <- round(data_melted$value,2)*100

# Define custom colors and legend order
custom_colors <- c("SJLIFE_PRS" = "#87CEFA", "CCSS_PRS" = "#FFA07A", "SJLIFE_All_treatments" = "#1E90FF", "CCSS_All_treatments" = "#FF6347")
legend_order <- c("SJLIFE_PRS", "CCSS_PRS", "SJLIFE_All_treatments", "CCSS_All_treatments")
# Order the levels of the legend_group factor
data_melted$legend_group <- factor(data_melted$legend_group, levels = legend_order)

# Define the order of x-axis categories
order <- c("Any SN", "Any SMN", "NMSC", "Breast cancer", "Thyroid cancer", "Meningioma", "Sarcoma")

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
p <- p + scale_y_continuous(limits = c(0, 100), expand = c(0, 0)) + labs(title = "", y = "Attributable fraction (%)", x = NULL) 
p

# Save the plot as a high-resolution image
plot_name <- paste0("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/ANALYSIS/results/plots/v18b/female/", group, "_", AF.type,"_", lifestyle, ".tiff")

ggsave(plot_name, p, width = 12, height = 8, dpi = 600, device = "tiff", compression = "lzw")

##################################
## Chemo and Radiation together ##
##################################
group <- "Female"

variables <- "Chemo|Radiation"

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

data_melted$legend_group <- factor(paste(data_melted$variable, data_melted$AF_by, sep = "_"))


# Define custom colors and legend order
custom_colors <- c("SJLIFE_Chemo" = "#87CEFA", "CCSS_Chemo" = "#FFA07A", "SJLIFE_Radiation" = "#1E90FF", "CCSS_Radiation" = "#FF6347")
legend_order <- c("SJLIFE_Chemo", "CCSS_Chemo", "SJLIFE_Radiation", "CCSS_Radiation")
# Order the levels of the legend_group factor
data_melted$legend_group <- factor(data_melted$legend_group, levels = legend_order)

# Define the order of x-axis categories
order <- c("Any SN", "Any SMN", "NMSC", "Breast cancer", "Thyroid cancer", "Meningioma", "Sarcoma")

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
p <- p + scale_y_continuous(limits = c(0, 100), expand = c(0, 0)) + labs(title = "", y = "Attributable fraction (%)", x = NULL) 

p
# Save the plot as a high-resolution image
plot_name <- paste0("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/ANALYSIS/results/plots/v18b/female/", group, "_", AF.type,"_", lifestyle, ".tiff")

ggsave(plot_name, p, width = 12, height = 8, dpi = 600, device = "tiff", compression = "lzw")

## 3.___________________________________________________________________________ Male analysis
## 3.___________________________________________________________________________ Male analysis
## 3.___________________________________________________________________________ Male analysis
## 3.___________________________________________________________________________ Male analysis
## 3.___________________________________________________________________________ Male analysis


data <- saved.data

group <- "Male"
variables <- unique(data$Variables)

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
      position = position_dodge(width = 0.3),
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
  plot_name <- paste0("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/ANALYSIS/results/plots/v18b/male/", group, "_", AF.type, "_", lifestyle, ".tiff")
  
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
order <- unique(data_melted$SN_types)
AF.type <- "treatment_PRS"
lifestyle <- "without_lifestyle"

data_melted$legend_group <- factor(paste(data_melted$variable, data_melted$AF_by, sep = "_"))


data_melted$new_value <- round(data_melted$value,2)*100

# Define custom colors and legend order
custom_colors <- c("SJLIFE_PRS" = "#87CEFA", "CCSS_PRS" = "#FFA07A", "SJLIFE_All_treatments" = "#1E90FF", "CCSS_All_treatments" = "#FF6347")
legend_order <- c("SJLIFE_PRS", "CCSS_PRS", "SJLIFE_All_treatments", "CCSS_All_treatments")
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
p <- p + scale_y_continuous(limits = c(0, 100), expand = c(0, 0)) + labs(title = "", y = "Attributable fraction (%)", x = NULL) 
p
# Save the plot as a high-resolution image
plot_name <- paste0("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/ANALYSIS/results/plots/v18b/male/", group, "_", AF.type,"_", lifestyle, ".tiff")

ggsave(plot_name, p, width = 12, height = 8, dpi = 600, device = "tiff", compression = "lzw")

##################################
## Chemo and Radiation together ##
##################################
group <- "Male"

variables <- "Chemo|Radiation"

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

data_melted$legend_group <- factor(paste(data_melted$variable, data_melted$AF_by, sep = "_"))


# Define custom colors and legend order
custom_colors <- c("SJLIFE_Chemo" = "#87CEFA", "CCSS_Chemo" = "#FFA07A", "SJLIFE_Radiation" = "#1E90FF", "CCSS_Radiation" = "#FF6347")
legend_order <- c("SJLIFE_Chemo", "CCSS_Chemo", "SJLIFE_Radiation", "CCSS_Radiation")
# Order the levels of the legend_group factor
data_melted$legend_group <- factor(data_melted$legend_group, levels = legend_order)

# Define the order of x-axis categories
order <- c("Any SN", "Any SMN", "NMSC", "Breast cancer", "Thyroid cancer", "Meningioma", "Sarcoma")

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
p <- p + scale_y_continuous(limits = c(0, 100), expand = c(0, 0)) + labs(title = "", y = "Attributable fraction (%)", x = NULL) 

p
# Save the plot as a high-resolution image
plot_name <- paste0("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/ANALYSIS/results/plots/v18b/male/", group, "_", AF.type,"_", lifestyle, ".tiff")

ggsave(plot_name, p, width = 12, height = 8, dpi = 600, device = "tiff", compression = "lzw")

## 4.___________________________________________________________________________ lt35 analysis
## 4.___________________________________________________________________________ lt35 analysis
## 4.___________________________________________________________________________ lt35 analysis
## 4.___________________________________________________________________________ lt35 analysis
## 4.___________________________________________________________________________ lt35 analysis


data <- saved.data

group <- "Age.lt.35"
variables <- unique(data$Variables)

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
      position = position_dodge(width = 0.3),
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
  plot_name <- paste0("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/ANALYSIS/results/plots/v18b/lt35/", group, "_", AF.type, "_", lifestyle, ".tiff")
  
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
# order <- c("Any SN", "Any SMN", "NMSC", "Breast cancer", "Thyroid cancer", "Meningioma", "Sarcoma")
order <- unique(data_melted$SN_types)
AF.type <- "treatment_PRS"
lifestyle <- "without_lifestyle"

data_melted$legend_group <- factor(paste(data_melted$variable, data_melted$AF_by, sep = "_"))


data_melted$new_value <- round(data_melted$value,2)*100

# Define custom colors and legend order
custom_colors <- c("SJLIFE_PRS" = "#87CEFA", "CCSS_PRS" = "#FFA07A", "SJLIFE_All_treatments" = "#1E90FF", "CCSS_All_treatments" = "#FF6347")
legend_order <- c("SJLIFE_PRS", "CCSS_PRS", "SJLIFE_All_treatments", "CCSS_All_treatments")
# Order the levels of the legend_group factor
data_melted$legend_group <- factor(data_melted$legend_group, levels = legend_order)

# Define the order of x-axis categories
order <- c("Any SN", "Any SMN", "NMSC", "Breast cancer", "Thyroid cancer", "Meningioma", "Sarcoma")

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
p <- p + scale_y_continuous(limits = c(0, 100), expand = c(0, 0)) + labs(title = "", y = "Attributable fraction (%)", x = NULL) 
p
# Save the plot as a high-resolution image
plot_name <- paste0("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/ANALYSIS/results/plots/v18b/lt35/", group, "_", AF.type,"_", lifestyle, ".tiff")

ggsave(plot_name, p, width = 12, height = 8, dpi = 600, device = "tiff", compression = "lzw")

##################################
## Chemo and Radiation together ##
##################################
group <- "Age.lt.35"

variables <- "Chemo|Radiation"

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

data_melted$legend_group <- factor(paste(data_melted$variable, data_melted$AF_by, sep = "_"))


# Define custom colors and legend order
custom_colors <- c("SJLIFE_Chemo" = "#87CEFA", "CCSS_Chemo" = "#FFA07A", "SJLIFE_Radiation" = "#1E90FF", "CCSS_Radiation" = "#FF6347")
legend_order <- c("SJLIFE_Chemo", "CCSS_Chemo", "SJLIFE_Radiation", "CCSS_Radiation")
# Order the levels of the legend_group factor
data_melted$legend_group <- factor(data_melted$legend_group, levels = legend_order)

# Define the order of x-axis categories
order <- c("Any SN", "Any SMN", "NMSC", "Breast cancer", "Thyroid cancer", "Meningioma", "Sarcoma")

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
p <- p + scale_y_continuous(limits = c(0, 100), expand = c(0, 0)) + labs(title = "", y = "Attributable fraction (%)", x = NULL) 

p
# Save the plot as a high-resolution image
plot_name <- paste0("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/ANALYSIS/results/plots/v18b/lt35/", group, "_", AF.type,"_", lifestyle, ".tiff")

ggsave(plot_name, p, width = 12, height = 8, dpi = 600, device = "tiff", compression = "lzw")

## 5.___________________________________________________________________________ ge35 analysis
## 5.___________________________________________________________________________ ge35 analysis
## 5.___________________________________________________________________________ ge35 analysis
## 5.___________________________________________________________________________ ge35 analysis
## 5.___________________________________________________________________________ ge35 analysis

data <- saved.data

group <- "Age.ge.35"
variables <- unique(data$Variables)

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
      position = position_dodge(width = 0.3),
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
  plot_name <- paste0("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/ANALYSIS/results/plots/v18b/ge35/", group, "_", AF.type, "_", lifestyle, ".tiff")
  
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
order <- unique(data_melted$SN_types)
AF.type <- "treatment_PRS"
lifestyle <- "without_lifestyle"

data_melted$legend_group <- factor(paste(data_melted$variable, data_melted$AF_by, sep = "_"))


data_melted$new_value <- round(data_melted$value,2)*100

# Define custom colors and legend order
custom_colors <- c("SJLIFE_PRS" = "#87CEFA", "CCSS_PRS" = "#FFA07A", "SJLIFE_All_treatments" = "#1E90FF", "CCSS_All_treatments" = "#FF6347")
legend_order <- c("SJLIFE_PRS", "CCSS_PRS", "SJLIFE_All_treatments", "CCSS_All_treatments")
# Order the levels of the legend_group factor
data_melted$legend_group <- factor(data_melted$legend_group, levels = legend_order)

# Define the order of x-axis categories
order <- c("Any SN", "Any SMN", "NMSC", "Breast cancer", "Thyroid cancer", "Meningioma", "Sarcoma")

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
p <- p + scale_y_continuous(limits = c(0, 100), expand = c(0, 0)) + labs(title = "", y = "Attributable fraction (%)", x = NULL) 
p
# Save the plot as a high-resolution image
plot_name <- paste0("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/ANALYSIS/results/plots/v18b/ge35/", group, "_", AF.type,"_", lifestyle, ".tiff")

ggsave(plot_name, p, width = 12, height = 8, dpi = 600, device = "tiff", compression = "lzw")

##################################
## Chemo and Radiation together ##
##################################
group <- "Age.ge.35"

variables <- "Chemo|Radiation"

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

data_melted$legend_group <- factor(paste(data_melted$variable, data_melted$AF_by, sep = "_"))


# Define custom colors and legend order
custom_colors <- c("SJLIFE_Chemo" = "#87CEFA", "CCSS_Chemo" = "#FFA07A", "SJLIFE_Radiation" = "#1E90FF", "CCSS_Radiation" = "#FF6347")
legend_order <- c("SJLIFE_Chemo", "CCSS_Chemo", "SJLIFE_Radiation", "CCSS_Radiation")
# Order the levels of the legend_group factor
data_melted$legend_group <- factor(data_melted$legend_group, levels = legend_order)

# Define the order of x-axis categories
order <- c("Any SN", "Any SMN", "NMSC", "Breast cancer", "Thyroid cancer", "Meningioma", "Sarcoma")

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
p <- p + scale_y_continuous(limits = c(0, 100), expand = c(0, 0)) + labs(title = "", y = "Attributable fraction (%)", x = NULL) 

p
# Save the plot as a high-resolution image
plot_name <- paste0("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/ANALYSIS/results/plots/v18b/ge35/", group, "_", AF.type,"_", lifestyle, ".tiff")

ggsave(plot_name, p, width = 12, height = 8, dpi = 600, device = "tiff", compression = "lzw")