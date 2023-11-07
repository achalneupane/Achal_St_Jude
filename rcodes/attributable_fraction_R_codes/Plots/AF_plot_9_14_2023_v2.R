# Load required libraries
library(ggplot2)
library(ggthemes)  # For additional themes
library(ggrepel)
library(dplyr)

## 1. Combined analysis

## Create a data frame with your data
# data <- data.frame(
#   SN_types = c("Any SN", "Any SMN", "NMSC", "Breast cancer", "Thyroid cancer", "Meningioma", "Sarcoma"),
#   SJLIFE = c(0.518, 0.455, 0.652, 0.672, 0.862, 0.438, 0.371),
#   CCSS = c(0.419, 0.326, 0.555, 0.709, 0.659, 0.432, 0.368)
# )


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

data$SN_types <- gsub("\\([0-9]+\\)", "", data$SN_types)
data[data == "-"] <- NA


##___________________________________________________________________________



group <- "Overall"
variables <- unique(data$Variables)

for(i in 1:length(variables)){
new_data <- data[grepl(variables[i], data$Variables), c("Cohort", "SN_types", "Variables", group)]
colnames(new_data) <- c("variable", "SN_types", "AF_by", "value")
new_data$value <- as.numeric(new_data$value)

# Reshape the data for ggplot2
library(reshape2)
# data_melted <- melt(data, id.vars = "SN_types")
# data_melted
# SN_types variable value
# 1          Any SN   SJLIFE 0.518
# 2             Any SMN   SJLIFE 0.455
# 3            NMSC   SJLIFE 0.652
data_melted <-new_data
## Exclude rows with NAs
data_melted <- data_melted[complete.cases(data_melted),]

# order <- c("Any SN", "Any SMN", "NMSC", "Breast cancer", "Thyroid cancer", "Meningioma", "Sarcoma")
order <- unique(data_melted$SN_types)
AF.type <- data_melted$AF_by[i] 
lifestyle <- "without_lifestyle"
if(lifestyle== "without_lifestyle" && variables[i] == "Lifestyle"){
  next
}

data_melted$legend_group <- factor(data_melted$variable, levels = c("SJLIFE", "CCSS"))

order <- unique(data_melted$SN_types)
AF.type <- data_melted$AF_by[i] 
lifestyle <- "without_lifestyle"



# Create the grouped bar chart
p <- ggplot(data_melted, aes(x = SN_types, y = value, fill = legend_group)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.6) +
  # Customize the theme and appearance
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12, color = "black"),  # Adjust font size and angle
        axis.text.y = element_text(size = 12, color = "black"),
        axis.title.y = element_text(size = 14, color = "black"),  # Customize Y-axis title
        legend.title = element_blank(),  # Remove legend title
        legend.text = element_text(size = 12, color = "black"),  # Adjust legend font size
        plot.title = element_text(size = 16, hjust = 0.5, vjust = 1.5, face = "bold"),  # Title formatting
        panel.grid.major = element_blank(),  # Remove major gridlines
        panel.grid.minor = element_blank()   # Remove minor gridlines
        
  ) +
  
  # Customize colors
  scale_fill_manual(values = c("SJLIFE" = "#1E90FF", "CCSS" = "#FF6347")) +  # Use color codes
  # Add labels
  labs(title = "", y = "Attributable Fraction", x = NULL) +
  # Remove extra space at the bottom of the Y-axis
  scale_y_continuous(expand = c(0, 0)) +
  # Add data labels
  geom_text(aes(label = sprintf("%.3f", value)), position = position_dodge(width = 0.8), vjust = -0.25, size = 3.5, color = "black") +
  # Adjust legend position
  theme(legend.position="top", legend.box = "horizontal") +
  ## Add caption
  # labs(caption = "Data")
  # scale_x_discrete(limits = c("Any SN ", "Any SMN ", "NMSC ", "Breast cancer ", "Thyroid cancer ", "Meningioma ", "Sarcoma ")) +
  scale_x_discrete(limits = order) +
  ## Adjust the x-axis limits to ensure that the secondary y-axis label is fully visible.
  coord_cartesian(clip = "off")

p <- p + scale_y_continuous(limits = c(0, 1))
# Save the plot as a high-resolution image
plot_name <- paste0("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/ANALYSIS/results/plots/v18b/", group, "_", AF.type, "_", lifestyle, ".tiff")

ggsave(plot_name, p, width = 10, height = 6, dpi = 300)
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
  # data_melted <- melt(data, id.vars = "SN_types")
  # data_melted
  # SN_types variable value
  # 1          Any SN   SJLIFE 0.518
  # 2             Any SMN   SJLIFE 0.455
  # 3            NMSC   SJLIFE 0.652
  data_melted <-new_data
  ## Exclude rows with NAs
  data_melted <- data_melted[complete.cases(data_melted),]
  
  # order <- c("Any SN", "Any SMN", "NMSC", "Breast cancer", "Thyroid cancer", "Meningioma", "Sarcoma")
  order <- unique(data_melted$SN_types)
  AF.type <- "treatment_PRS"
  lifestyle <- "without_lifestyle"
  
  data_melted$legend_group <- factor(paste(data_melted$variable, data_melted$AF_by, sep = "_"))
  
  
  data_melted$new_value <- data_melted$value
  # data_melted <- data_melted %>%
  #   group_by(SN_types, value) %>%
  #   mutate(new_value = ifelse(row_number() == which.min(new_value), new_value, NA)) %>%
  #   ungroup()

  
  
  
  # Define custom colors and legend order
  custom_colors <- c("SJLIFE_PRS" = "#87CEFA", "CCSS_PRS" = "#FFA07A", "SJLIFE_All_treatments" = "#1E90FF", "CCSS_All_treatments" = "#FF6347")
  legend_order <- c("SJLIFE_PRS", "CCSS_PRS", "SJLIFE_All_treatments", "CCSS_All_treatments")
  # Order the levels of the legend_group factor
  data_melted$legend_group <- factor(data_melted$legend_group, levels = legend_order)
  
  # Define the order of x-axis categories
  x_order <- c("Any SN ", "Any SMN ", "NMSC ", "Breast cancer ", "Thyroid cancer ", "Meningioma ", "Sarcoma ")
  
  # Create a factor with the desired order
  data_melted$SN_types <- factor(data_melted$SN_types, levels = x_order)
  
  # Create the grouped bar chart
  # Create the plot
  p <- ggplot(data_melted, aes(x = SN_types, y = value, fill = legend_group)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.6) +
    # Customize the theme and appearance
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12, color = "black"),
          axis.text.y = element_text(size = 12, color = "black"),
          axis.title.y = element_text(size = 14, color = "black"),
          legend.title = element_blank(),
          legend.text = element_text(size = 12, color = "black"),
          plot.title = element_text(size = 16, hjust = 0.5, vjust = 1.5, face = "bold"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    
    # Customize colors
    scale_fill_manual(values = custom_colors, breaks = legend_order) +
    
    # Add labels with repelling
    geom_text_repel(
      data = data_melted %>% filter(!is.na(new_value)),
      aes(label = sprintf("%.3f", new_value)),
      position = position_dodge(width = 0.8),
      vjust = 1.35,
      hjust = 0.15,
      size = 3.5,
      color = "black",
      box.padding = 0.15,
      force = 1,
      segment.color = "transparent",
      segment.size = 0,
      fill = "white"
    ) +
    
    # Adjust legend position
    theme(legend.position = "top", legend.box = "horizontal") +
    
    # Additional adjustments
    scale_x_discrete(limits = x_order) +
    coord_cartesian(clip = "off")
  
  p <- p + scale_y_continuous(limits = c(0, 1))
  p
  # Save the plot as a high-resolution image
  plot_name <- paste0("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/ANALYSIS/results/plots/v18b/", group, "_", AF.type,  "_", lifestyle, ".tiff")
  
  ggsave(plot_name, p, width = 10, height = 6, dpi = 300)

  
  
  
  
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
  # data_melted <- melt(data, id.vars = "SN_types")
  # data_melted
  # SN_types variable value
  # 1          Any SN   SJLIFE 0.518
  # 2             Any SMN   SJLIFE 0.455
  # 3            NMSC   SJLIFE 0.652
  data_melted <-new_data
  ## Exclude rows with NAs
  data_melted <- data_melted[complete.cases(data_melted),]
  
  # order <- c("Any SN", "Any SMN", "NMSC", "Breast cancer", "Thyroid cancer", "Meningioma", "Sarcoma")
  order <- unique(data_melted$SN_types)
  AF.type <- "Chemo_Radiation"
  lifestyle <- "without_lifestyle"
  
  data_melted$legend_group <- factor(paste(data_melted$variable, data_melted$AF_by, sep = "_"))
  
  
  data_melted$new_value <- data_melted$value
  # data_melted <- data_melted %>%
  #   group_by(SN_types, value) %>%
  #   mutate(new_value = ifelse(row_number() == which.min(new_value), new_value, NA)) %>%
  #   ungroup()
  
  
  
  
  # Define custom colors and legend order
  custom_colors <- c("SJLIFE_Chemo" = "#87CEFA", "CCSS_Chemo" = "#FFA07A", "SJLIFE_Radiation" = "#1E90FF", "CCSS_Radiation" = "#FF6347")
  legend_order <- c("SJLIFE_Chemo", "CCSS_Chemo", "SJLIFE_Radiation", "CCSS_Radiation")
  # Order the levels of the legend_group factor
  data_melted$legend_group <- factor(data_melted$legend_group, levels = legend_order)
  
  # Define the order of x-axis categories
  x_order <- c("Any SN ", "Any SMN ", "NMSC ", "Breast cancer ", "Thyroid cancer ", "Meningioma ", "Sarcoma ")
  
  # Create a factor with the desired order
  data_melted$SN_types <- factor(data_melted$SN_types, levels = x_order)
  
  # Create the grouped bar chart
  # Create the plot
  p <- ggplot(data_melted, aes(x = SN_types, y = value, fill = legend_group)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.6) +
    # Customize the theme and appearance
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12, color = "black"),
          axis.text.y = element_text(size = 12, color = "black"),
          axis.title.y = element_text(size = 14, color = "black"),
          legend.title = element_blank(),
          legend.text = element_text(size = 12, color = "black"),
          plot.title = element_text(size = 16, hjust = 0.5, vjust = 1.5, face = "bold"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    
    # Customize colors
    scale_fill_manual(values = custom_colors, breaks = legend_order) +
    
    # Add labels with repelling
    geom_text_repel(
      data = data_melted %>% filter(!is.na(new_value)),
      aes(label = sprintf("%.3f", new_value)),
      position = position_dodge(width = 0.8),
      vjust = 1.35,
      hjust = 0.15,
      size = 3.5,
      color = "black",
      box.padding = 0.15,
      force = 1,
      segment.color = "transparent",
      segment.size = 0,
      fill = "white"
    ) +
    
    # Adjust legend position
    theme(legend.position = "top", legend.box = "horizontal") +
    
    # Additional adjustments
    scale_x_discrete(limits = x_order) +
    coord_cartesian(clip = "off")
  
  p <- p + scale_y_continuous(limits = c(0, 1))
  p
  # Save the plot as a high-resolution image
  plot_name <- paste0("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/ANALYSIS/results/plots/v18b/", group, "_", AF.type,  "_", lifestyle, ".tiff")
  
  ggsave(plot_name, p, width = 10, height = 6, dpi = 300)

  
  
  
  
  
############################
## Lifestlye without diet ##
############################
  
  ############################
  ## Lifestyle without diet ##
  ############################
  data <- read.table(text = "Cohort	SN_types	Variables	Overall	Female	Male	Age.lt.35	Age.ge.35
SJLIFE	Any SN (303)	Radiation	0.442	0.432	0.455	0.368	0.473
SJLIFE	Any SN (303)	Chemo	0.017	0.017	0.016	0.031	0.011
SJLIFE	Any SN (303)	All_treatments	0.450	0.442	0.461	0.385	0.477
SJLIFE	Any SN (303)	PRS	0.179	0.177	0.183	0.182	0.178
SJLIFE	Any SN (303)	Lifestyle	-0.063	-0.058	-0.069	-0.058	-0.065
SJLIFE	Any SN (303)	Combined	0.520	0.514	0.527	0.462	0.543
SJLIFE	Any SMN (234)	Radiation	0.379	0.378	0.379	0.289	0.413
SJLIFE	Any SMN (234)	Chemo	0.017	0.018	0.017	0.030	0.013
SJLIFE	Any SMN (234)	All_treatments	0.389	0.390	0.388	0.308	0.421
SJLIFE	Any SMN (234)	PRS	0.143	0.140	0.146	0.142	0.143
SJLIFE	Any SMN (234)	Lifestyle	-0.083	-0.073	-0.096	-0.082	-0.084
SJLIFE	Any SMN (234)	Combined	0.431	0.437	0.423	0.354	0.460
SJLIFE	NMSC (118)	Radiation	0.294	0.282	0.310	0.239	0.306
SJLIFE	NMSC (118)	Chemo	-	-	-	-	-
SJLIFE	NMSC (118)	All_treatments	0.294	0.282	0.310	0.239	0.306
SJLIFE	NMSC (118)	PRS	0.432	0.418	0.448	0.431	0.432
SJLIFE	NMSC (118)	Lifestyle	-0.275	-0.251	-0.307	-0.242	-0.283
SJLIFE	NMSC (118)	Combined	0.490	0.482	0.500	0.459	0.497
SJLIFE	Breast cancer (51)	Radiation	0.466	0.466	-	0.391	0.478
SJLIFE	Breast cancer (51)	Chemo	0.106	0.106	-	0.102	0.107
SJLIFE	Breast cancer (51)	All_treatments	0.526	0.526	-	0.452	0.538
SJLIFE	Breast cancer (51)	PRS	0.222	0.222	-	0.250	0.218
SJLIFE	Breast cancer (51)	Lifestyle	-0.872	-0.872	-	-0.581	-0.919
SJLIFE	Breast cancer (51)	Combined	0.335	0.335	-	0.337	0.335
SJLIFE	Thyroid cancer (43)	Radiation	0.610	0.633	0.585	0.555	0.637
SJLIFE	Thyroid cancer (43)	Chemo	0.238	0.237	0.238	0.332	0.192
SJLIFE	Thyroid cancer (43)	All_treatments	0.724	0.740	0.706	0.708	0.732
SJLIFE	Thyroid cancer (43)	PRS	0.589	0.587	0.592	0.601	0.583
SJLIFE	Thyroid cancer (43)	Lifestyle	-0.178	-0.189	-0.166	-0.146	-0.194
SJLIFE	Thyroid cancer (43)	Combined	0.870	0.875	0.863	0.863	0.873
SJLIFE	Meningioma (81)	Radiation	0.184	0.150	0.225	0.245	0.164
SJLIFE	Meningioma (81)	Chemo	0.321	0.325	0.316	0.462	0.274
SJLIFE	Meningioma (81)	All_treatments	0.456	0.436	0.480	0.623	0.400
SJLIFE	Meningioma (81)	PRS	-0.117	-0.113	-0.122	-0.114	-0.118
SJLIFE	Meningioma (81)	Lifestyle	-0.135	-0.126	-0.146	-0.138	-0.134
SJLIFE	Meningioma (81)	Combined	0.308	0.293	0.325	0.509	0.240
SJLIFE	Sarcoma (NA)	Radiation	-	-	-	-	-
SJLIFE	Sarcoma (NA)	Chemo	-	-	-	-	-
SJLIFE	Sarcoma (NA)	All_treatments	-	-	-	-	-
SJLIFE	Sarcoma (NA)	PRS	-	-	-	-	-
SJLIFE	Sarcoma (NA)	Lifestyle	-	-	-	-	-
SJLIFE	Sarcoma (NA)	Combined	-	-	-	-	-
CCSS	Any SN (1286)	Radiation	0.358	0.351	0.369	0.338	0.368
CCSS	Any SN (1286)	Chemo	0.038	0.037	0.041	0.061	0.027
CCSS	Any SN (1286)	All_treatments	0.384	0.376	0.396	0.380	0.387
CCSS	Any SN (1286)	PRS	0.057	0.056	0.058	0.056	0.057
CCSS	Any SN (1286)	Lifestyle	-0.074	-0.071	-0.078	-0.068	-0.077
CCSS	Any SN (1286)	Combined	0.376	0.370	0.386	0.374	0.378
CCSS	Any SMN (619)	Radiation	0.249	0.246	0.255	0.186	0.275
CCSS	Any SMN (619)	Chemo	0.056	0.055	0.058	0.099	0.038
CCSS	Any SMN (619)	All_treatments	0.297	0.293	0.306	0.271	0.308
CCSS	Any SMN (619)	PRS	0.047	0.047	0.048	0.047	0.047
CCSS	Any SMN (619)	Lifestyle	-0.048	-0.046	-0.052	-0.046	-0.049
CCSS	Any SMN (619)	Combined	0.297	0.295	0.303	0.272	0.308
CCSS	NMSC (632)	Radiation	0.390	0.386	0.394	0.400	0.386
CCSS	NMSC (632)	Chemo	-	-	-	-	-
CCSS	NMSC (632)	All_treatments	0.390	0.386	0.394	0.400	0.386
CCSS	NMSC (632)	PRS	0.303	0.305	0.300	0.302	0.303
CCSS	NMSC (632)	Lifestyle	-0.149	-0.148	-0.152	-0.132	-0.156
CCSS	NMSC (632)	Combined	0.514	0.515	0.513	0.530	0.508
CCSS	Breast cancer (259)	Radiation	0.444	0.444	-	0.392	0.450
CCSS	Breast cancer (259)	Chemo	0.194	0.194	-	0.234	0.189
CCSS	Breast cancer (259)	All_treatments	0.572	0.572	-	0.541	0.575
CCSS	Breast cancer (259)	PRS	0.346	0.346	-	0.348	0.345
CCSS	Breast cancer (259)	Lifestyle	-0.088	-0.088	-	-0.078	-0.089
CCSS	Breast cancer (259)	Combined	0.696	0.696	-	0.678	0.698
CCSS	Thyroid cancer (123)	Radiation	0.416	0.416	0.418	0.360	0.486
CCSS	Thyroid cancer (123)	Chemo	0.095	0.093	0.097	0.122	0.061
CCSS	Thyroid cancer (123)	All_treatments	0.479	0.477	0.482	0.443	0.523
CCSS	Thyroid cancer (123)	PRS	0.386	0.391	0.376	0.379	0.394
CCSS	Thyroid cancer (123)	Lifestyle	-0.042	-0.036	-0.054	-0.025	-0.063
CCSS	Thyroid cancer (123)	Combined	0.666	0.667	0.666	0.646	0.692
CCSS	Meningioma (209)	Radiation	0.294	0.260	0.345	0.315	0.279
CCSS	Meningioma (209)	Chemo	0.087	0.080	0.098	0.134	0.055
CCSS	Meningioma (209)	All_treatments	0.368	0.334	0.420	0.421	0.331
CCSS	Meningioma (209)	PRS	0.028	0.027	0.029	0.033	0.024
CCSS	Meningioma (209)	Lifestyle	0.087	0.099	0.068	0.150	0.043
CCSS	Meningioma (209)	Combined	0.433	0.409	0.469	0.524	0.370
CCSS	Sarcoma (41)	Radiation	-	-	-	-	-
CCSS	Sarcoma (41)	Chemo	0.279	0.266	0.300	0.279	0.278
CCSS	Sarcoma (41)	All_treatments	0.279	0.266	0.300	0.279	0.278
CCSS	Sarcoma (41)	PRS	0.035	0.038	0.030	0.033	0.036
CCSS	Sarcoma (41)	Lifestyle	-0.261	-0.256	-0.271	-0.396	-0.180
CCSS	Sarcoma (41)	Combined	0.116	0.108	0.130	0.017	0.176", header = T, sep = "\t")
  
  data$SN_types_new <- gsub("\\([0-9]+\\)", "", data$SN_types)
  data[data == "-"] <- NA