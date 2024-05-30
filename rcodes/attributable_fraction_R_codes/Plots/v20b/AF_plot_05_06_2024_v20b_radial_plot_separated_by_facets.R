# Load required libraries
library(ggplot2)
library(ggthemes)# For additional themes
library(ggrepel)
library(dplyr)
library(ggpattern)
# https://r-charts.com/colors/

## V20b (Without lifestyle)
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

## Radial plot
group <- "Overall"

variables <- "Combined"

new_data <- saved.data[grepl(variables, saved.data$Variables), c("Cohort", "SN_types", "Variables", group)]
colnames(new_data) <- c("variable", "SN_types", "AF_by", "value")
new_data$value <- as.numeric(new_data$value)

# Reshape the data for ggplot2
library(reshape2)

data_melted <-new_data
order <- unique(data_melted$SN_types)
AF.type <- "Combined"
lifestyle <- "without_lifestyle"

data_melted$AF_by_new <-  gsub("_", " ", data_melted$AF_by)
data_melted$legend_group <- factor(data_melted$variable)

data_melted$new_value <- round(data_melted$value,2)*100


# Define the desired order of SN_types within each variable
custom_order_within_variable <- c("SNs", "SMNs", "NMSCs", "Breast cancer", "Thyroid cancer", "Meningioma", "Sarcoma")

# Create a custom factor variable based on the desired order
data_melted <- data_melted %>%
  dplyr::mutate(SN_types = factor(SN_types, levels = custom_order_within_variable))
data_melted <- data_melted %>%
  arrange(SN_types, AF_by)

# data <- cbind.data.frame(Cohort=data_melted$variable, Risks=gsub("All_t", "T", data_melted$AF_by), individual= paste0(data_melted$new_value, "%"), group= data_melted$SN_types, group2= data_melted$legend_group, value=data_melted$new_value)
# # data$value[is.na(data$value)]
# data$Cohort <- factor(data$Cohort, levels= c("SJLIFE", "CCSS"))
# data$Risks <- factor(data$Risks, levels= c("Treatments", "PRS"))
# 
# custom_colors <- c("SJLIFE" = "#EE6A50", "CCSS" = "#9FB6CD")
# legend_order <- c("SJLIFE", "CCSS")
# 
# 
# # Set a number of 'empty bar' to add at the end of each group
# empty_bar <- 1
# to_add <- data.frame( matrix(NA, empty_bar*nlevels(data$group), ncol(data)) )
# colnames(to_add) <- colnames(data)
# to_add$group <- rep(levels(data$group), each=empty_bar)
# data <- rbind(data, to_add)
# data <- data %>% arrange(group)
# data$id <- seq(1, nrow(data))
# 
# # Get the name and the y position of each label
# label_data <- data
# number_of_bar <- nrow(label_data)
# angle <- 90 - 360 * (label_data$id-0.5) /number_of_bar     # I substract 0.5 because the letter must have the angle of the center of the bars. Not extreme right(1) or extreme left (0)
# label_data$hjust <- ifelse( angle < -90, 1, 0)
# label_data$angle <- ifelse(angle < -90, angle+180, angle)
# 
# 
# # prepare a data frame for base lines
# base_data <- data %>% 
#   dplyr::group_by(group) %>% 
#   dplyr::summarize(start=min(id), end=max(id) - empty_bar) %>% 
#   dplyr::rowwise() %>% 
#   dplyr::mutate(title=mean(c(start, end)))
# 
# # prepare a data frame for grid (scales)
# grid_data <- base_data
# grid_data$end <- grid_data$end[ c( nrow(grid_data), 1:nrow(grid_data)-1)] + 1
# grid_data$start <- grid_data$start - 1
# grid_data <- grid_data[-1,]
# 
# # Make the plot
# p = ggplot(data, aes(x=as.factor(id), y=value, fill=Cohort)) +       # Note that id is a factor. If x is numeric, there is some space between the first bar
#   
#   # geom_bar(aes(x=as.factor(id), y=value, fill=group2), stat="identity", alpha=1, width=0.95) +
# 
#   # # Add a val=100/75/50/25 lines. I do it at the beginning to make sur barplots are OVER it.
#   # geom_segment(data=grid_data, aes(x = end, y = 80, xend = start, yend = 80), colour = "white", alpha=1, size=0.3 , inherit.aes = FALSE ) +
#   # geom_segment(data=grid_data, aes(x = end, y = 60, xend = start, yend = 60), colour = "white", alpha=1, size=0.3 , inherit.aes = FALSE ) +
#   # geom_segment(data=grid_data, aes(x = end, y = 40, xend = start, yend = 40), colour = "white", alpha=1, size=0.3 , inherit.aes = FALSE ) +
#   # geom_segment(data=grid_data, aes(x = end, y = 20, xend = start, yend = 20), colour = "white", alpha=1, size=0.3 , inherit.aes = FALSE ) +
#   
#   ## Add text showing the value of each 100/75/50/25 lines
#   # annotate("text", x = rep(max(data$id),5), y = c(0, 25, 50, 75, 100), label = c("0%", "25%", "50%", "75%", "100%") , color="black", size=3 , angle=0, fontface= 1, hjust=1) +
#   
#   geom_bar(aes(x=as.factor(id), y=value, fill=Cohort), stat="identity", alpha=1, width=0.95) +
#   theme_bw() +
#   ylim(-100,120) +
#   theme_minimal() +
#   theme(
#     legend.position = "right",
#     axis.text = element_blank(),
#     axis.title = element_blank(),
#     panel.grid = element_blank(),
#     plot.margin = unit(rep(-1,4), "cm") 
#   ) +
#   coord_polar() + 
#   geom_text(data=label_data, aes(x=id, y=value+2, label=individual, hjust=hjust), color="black", fontface=1,alpha=1, size=6.5, angle= label_data$angle, inherit.aes = FALSE ) +
#   
#   # Add base line information
#   # geom_segment(data=base_data, aes(x = start, y = -5, xend = end, yend = -5), colour = "black", alpha=0.8, size=1.5 , inherit.aes = FALSE )  +
#   # geom_text(data=base_data, aes(x = title, y = -18, label=group), hjust=c(1,1,1,0,0,0,0), colour = "black", alpha=0.8, size=4.5, fontface=1, inherit.aes = FALSE) +
#   geomtextpath::coord_curvedpolar() +
#   geomtextpath::geom_textsegment(
#     data = base_data,
#     aes(
#       x = start, y = -.1,
#       xend = end, yend = -.1,
#       label = group
#     ),
#     colour = "black",
#     linewidth = 2.5,
#     size = 16 / .pt,
#     inherit.aes = FALSE,
#     gap = FALSE,
#     offset = unit(-24, "pt")
#   ) +
#   scale_fill_manual(values = custom_colors, breaks = legend_order) +
#   guides(
#     theta = guide_axis_theta(angle = 0)
#   ) +
#   theme(
#     legend.position="right",
#     legend.justification="right", 
#     legend.box.spacing = unit(-120, "pt"),# The spacing between the plotting area and the legend box (unit)
#     legend.margin=margin(0,0,0,0),
#     # axis.text = element_text(size = 18, color = "black"), # Adjust axis text size
#     # axis.title = element_text(size = 18, color = "black"), # Adjust axis title size
#     legend.title = element_text(size = 16, color = "black", face=2), # Adjust legend title size
#     plot.margin = margin(r = .5, l = -3.5, t = -3.5, b = -3.5, unit = "cm"),
#     legend.text = element_text(size = 14, color = "black")
#   ) + # Adjust legend text size
#   labs(
#     title = "",
#     x = NULL, y = "",
#     fill = "Cohorts"
#   )
# p
# plot_name <- "Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/ANALYSIS/results/plots/v20b/overall/Overall_Combined_without_lifestyle_radial_plot.tiff"
# ggsave(plot_name, p, width = 12, height = 8, dpi = 600, device = "tiff", compression = "lzw")


################################
## PRS and treatment together ##
################################
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
custom_order_within_variable <- c("SNs", "SMNs", "NMSCs", "Breast cancer", "Thyroid cancer", "Meningioma", "Sarcoma")

# Create a custom factor variable based on the desired order
data_melted <- data_melted %>%
  dplyr::mutate(SN_types = factor(SN_types, levels = custom_order_within_variable))
data_melted <- data_melted %>%
  arrange(SN_types, factor(AF_by, levels = c("All_treatments", "PRS")))

data.all <- cbind.data.frame(Cohort=data_melted$variable, Risks=gsub("All_t", "T", data_melted$AF_by), individual= paste0(data_melted$new_value, "%"), group= data_melted$SN_types, group2= data_melted$legend_group, value=data_melted$new_value)
# data$value[is.na(data$value)]

## SJLIFE
data <- data.all[data.all$Cohort == "SJLIFE",]

data$Cohort <- factor(data$Cohort, levels= c("SJLIFE"))
data$Risks <- factor(data$Risks, levels= c("Treatments", "PRS"))

custom_colors <- c("Treatments" = "#EE6A50", PRS = "#EE6A50")
legend_order <- c("Treatments", "PRS")

# Set a number of 'empty bar' to add at the end of each group
empty_bar <- 1
to_add <- data.frame( matrix(0, empty_bar*nlevels(data$group), ncol(data)) )
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

## Fix issues with the gaps
data$value[is.na(data$Cohort)] <- NA
data$individual[is.na(data$Cohort)] <- NA
label_data$value[is.na(label_data$Cohort)] <- NA
label_data$individual[is.na(label_data$Cohort)] <- NA
data$Cohort[is.na(data$Cohort)] <- "SJLIFE"


# Make the plot
p_sjlife = ggplot(data, aes(x=as.factor(id), y=value, fill = Risks)) +  
  geom_col(
    position = "dodge",
    width = 0.95
  ) +
  geom_col_pattern(
    data = subset(data, Risks == "Treatments"),
    aes(fill = Risks),
    pattern = "stripe",
    pattern_fill = "black",
    pattern_angle = 45,
    pattern_density = 0.1,
    pattern_spacing = 0.01,
    position = position_dodge2(preserve = 'single'),
    na.rm = TRUE
  ) +
  scale_fill_discrete(na.translate = F) +
  # theme(legend.key.size = unit(1.5, 'cm')) +
  scale_pattern_manual(
    values = c("none", "stripe"),
    guide = guide_legend(override.aes = list(fill = "gray")) # <- make lighter
  ) +
  scale_fill_manual(
    values = custom_colors, 
    breaks = legend_order,
    labels = c(
      "Treatments" = "Higher exposure levels\nof cancer treatments",
      "PRS" = "Elevated genetic risk"
    ),
    # guide = guide_legend(override.aes = list(pattern = "none")) # <- hide pattern
  ) +
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
  geom_text(data=label_data, aes(x=id, y=value+2, label=individual, hjust=hjust), color="black", fontface=1,alpha=1, size=6.5, angle= label_data$angle, inherit.aes = FALSE ) +
  
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
    linewidth = 2.5,
    size = 16 / .pt,
    inherit.aes = FALSE,
    gap = FALSE,
    offset = unit(-24, "pt")
  ) +
  # guides(
  #   theta = guide_axis_theta(angle = 0)
  # ) +
  theme(
    legend.position="right",
    legend.justification="right", 
    legend.box.spacing = unit(-300, "pt"),# The spacing between the plotting area and the legend box (unit)
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
    fill = "Risk factors"
  )

p_sjlife


## CCSSS
data <- data.all[data.all$Cohort == "CCSS",]

data$Cohort <- factor(data$Cohort, levels= c("CCSS"))
data$Risks <- factor(data$Risks, levels= c("Treatments", "PRS"))

custom_colors <- c("Treatments" = "#9FB6CD", PRS = "#9FB6CD")
legend_order <- c("Treatments", "PRS")

# Set a number of 'empty bar' to add at the end of each group
empty_bar <- 1
to_add <- data.frame( matrix(0, empty_bar*nlevels(data$group), ncol(data)) )
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

## Fix issues with the gaps
data$value[is.na(data$Cohort)] <- NA
data$individual[is.na(data$Cohort)] <- NA
label_data$value[is.na(label_data$Cohort)] <- NA
label_data$individual[is.na(label_data$Cohort)] <- NA
data$Cohort[is.na(data$Cohort)] <- "CCSS"


# Make the plot
p_ccss = ggplot(data, aes(x=as.factor(id), y=value, fill = Risks)) +  
  geom_col(
    position = "dodge",
    width = 0.95
  ) +
  geom_col_pattern(
    data = subset(data, Risks == "Treatments"),
    aes(fill = Risks),
    pattern = "stripe",
    pattern_fill = "black",
    pattern_angle = 45,
    pattern_density = 0.1,
    pattern_spacing = 0.01,
    position = position_dodge2(preserve = 'single'),
    na.rm = TRUE
  ) +
  scale_fill_discrete(na.translate = F) +
  # theme(legend.key.size = unit(1.5, 'cm')) +
  scale_pattern_manual(
    values = c("none", "stripe"),
    guide = guide_legend(override.aes = list(fill = "gray")) # <- make lighter
  ) +
  scale_fill_manual(
    values = custom_colors, 
    breaks = legend_order,
    labels = c(
      "Treatments" = "Higher exposure levels\nof cancer treatments",
      "PRS" = "Elevated genetic risk"
    ),
    # guide = guide_legend(override.aes = list(pattern = "none")) # <- hide pattern
  ) +
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
  geom_text(data=label_data, aes(x=id, y=value+2, label=individual, hjust=hjust), color="black", fontface=1,alpha=1, size=6.5, angle= label_data$angle, inherit.aes = FALSE ) +
  
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
    linewidth = 2.5,
    size = 16 / .pt,
    inherit.aes = FALSE,
    gap = FALSE,
    offset = unit(-24, "pt")
  ) +
  # guides(
  #   theta = guide_axis_theta(angle = 0)
  # ) +
  theme(
    legend.position="right",
    legend.justification="right", 
    legend.box.spacing = unit(-300, "pt"),# The spacing between the plotting area and the legend box (unit)
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
    fill = "Risk factors"
  )

p_ccss

# library(cowplot)
# cowplot::plot_grid(p_sjlife, p_ccss, labels = c("SJLIFE", "CCSS"))
# p = cowplot::plot_grid(p_sjlife, p_ccss, labels = c("AUTO"))

# p_sjlife <- p_sjlife + theme(plot.margin = margin(r = 0, l = 0, t = 0, b = 0, "cm"))
# p_ccss <- p_ccss + theme(plot.margin = margin(r = 0, l = 0, t = 0, b = 0, "cm"))

p_sjlife <- p_sjlife + theme(plot.margin = margin(r = -15, l = -2, t = 0, b = 0, "cm"))
p_ccss <- p_ccss + theme(plot.margin = margin(r = -5, l = -10, t = 0, b = 0, "cm"))

library(ggplot2)
library(gridExtra)

p_sjlife <- p_sjlife + 
  theme(legend.position = "bottom") +
  guides(fill = guide_legend(override.aes = list(fill = "gray"))) # Set fill color to white
# Remove legend from the other plot
p_ccss <- p_ccss + theme(legend.position = "none")



p_sjlife_label <- p_sjlife +
  annotate("text", x = Inf, y = -Inf, label = "SJLIFE", vjust = 0, hjust = 1, size = 6)

# Add label "B" to p_ccss
p_ccss_label <- p_ccss +
  annotate("text", x = Inf, y = -Inf, label = "CCSS", vjust = 0, hjust = 1, size = 6)

# Combine the labeled plots
p <- ggarrange(
  p_sjlife_label, p_ccss_label, 
  ncol = 2, nrow = 1, 
  common.legend = TRUE, legend = "right",
  align = "hv"
)

p

p = p + theme(plot.margin = margin(r = 0, l = 0, t = 0, b = 0, "cm"))
p
plot_name <- "Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/ANALYSIS/results/plots/v20b/overall/Overall_treatment_and_PRS_without_lifestyle_radial_plot_faceted.tiff"
ggsave(plot_name, p, width = 18, height = 12, dpi = 400, device = "tiff", compression = "lzw")

# plot_name <- "Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/ANALYSIS/results/plots/v20b/overall/Overall_treatment_and_PRS_without_lifestyle_radial_plot_faceted_SJLIFE.tiff"
# ggsave(plot_name, p_sjlife, width = 12, height = 8, dpi = 600, device = "tiff", compression = "lzw")
# 
# plot_name <- "Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/ANALYSIS/results/plots/v20b/overall/Overall_treatment_and_PRS_without_lifestyle_radial_plot_faceted_CCSS.tiff"
# ggsave(plot_name, p_ccss, width = 12, height = 8, dpi = 600, device = "tiff", compression = "lzw")
# 

#########################################
## Chemotherapy and Radiation together ##
#########################################
## Radial plot
group <- "Overall"

variables <- "Chemotherapy|Radiation"

new_data <- saved.data[grepl(variables, saved.data$Variables), c("Cohort", "SN_types", "Variables", group)]
colnames(new_data) <- c("variable", "SN_types", "AF_by", "value")
new_data$value <- as.numeric(new_data$value)

# Reshape the data for ggplot2
library(reshape2)

#########################################
## Chemotherapy and Radiation together ##
#########################################
## Radial plot
group <- "Overall"

variables <- "Chemotherapy|Radiation"

new_data <- saved.data[grepl(variables, saved.data$Variables), c("Cohort", "SN_types", "Variables", group)]
colnames(new_data) <- c("variable", "SN_types", "AF_by", "value")
new_data$value <- as.numeric(new_data$value)

# Reshape the data for ggplot2
library(reshape2)

#########################################
## Radiation and Chemotherapy together ##
#########################################
## Radial plot
group <- "Overall"

variables <- "Chemotherapy|Radiation"

new_data <- saved.data[grepl(variables, saved.data$Variables), c("Cohort", "SN_types", "Variables", group)]
colnames(new_data) <- c("variable", "SN_types", "AF_by", "value")
new_data$value <- as.numeric(new_data$value)

# Reshape the data for ggplot2
library(reshape2)

data_melted <-new_data
order <- unique(data_melted$SN_types)
AF.type <- "Chemotherapy_Radiation"
lifestyle <- "without_lifestyle"

data_melted$AF_by_new <-  gsub("_", " ", data_melted$AF_by)
data_melted$legend_group <- factor(paste(data_melted$variable, data_melted$AF_by_new, sep = "-"))

data_melted$new_value <- round(data_melted$value,2)*100


# Define the desired order of SN_types within each variable
custom_order_within_variable <- c("SNs", "SMNs", "NMSCs", "Breast cancer", "Thyroid cancer", "Meningioma", "Sarcoma")

# Create a custom factor variable based on the desired order
data_melted <- data_melted %>%
  dplyr::mutate(SN_types = factor(SN_types, levels = custom_order_within_variable))
data_melted <- data_melted %>%
  arrange(SN_types, factor(AF_by, levels = c("Radiation", "Chemotherapy")))

data.all <- cbind.data.frame(Cohort=data_melted$variable, Risks=gsub("All_t", "T", data_melted$AF_by), individual= paste0(data_melted$new_value, "%"), group= data_melted$SN_types, group2= data_melted$legend_group, value=data_melted$new_value)
# data$value[is.na(data$value)]

## SJLIFE
data <- data.all[data.all$Cohort == "SJLIFE",]

data$Cohort <- factor(data$Cohort, levels= c("SJLIFE"))
data$Risks <- factor(data$Risks, levels= c("Radiation", "Chemotherapy"))

custom_colors <- c(Radiation = "#EE6A50", Chemotherapy = "#EE6A50")
legend_order <- c("Radiation", "Chemotherapy")

# Set a number of 'empty bar' to add at the end of each group
empty_bar <- 1
to_add <- data.frame( matrix(0, empty_bar*nlevels(data$group), ncol(data)) )
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

## Fix issues with the gaps
data$value[is.na(data$Cohort)] <- NA
data$individual[is.na(data$Cohort)] <- NA
label_data$value[is.na(label_data$Cohort)] <- NA
label_data$individual[is.na(label_data$Cohort)] <- NA
data$Cohort[is.na(data$Cohort)] <- "SJLIFE"


# Make the plot
p_sjlife = ggplot(data, aes(x=as.factor(id), y=value, fill = Risks)) +  
  geom_col(
    position = "dodge",
    width = 0.95
  ) +
  geom_col_pattern(
    data = subset(data, Risks == "Radiation"),
    aes(fill = Risks),
    pattern = "stripe",
    pattern_fill = "black",
    pattern_angle = 45,
    pattern_density = 0.1,
    pattern_spacing = 0.01,
    position = position_dodge2(preserve = 'single'),
    na.rm = TRUE
  ) +
  scale_fill_discrete(na.translate = F) +
  # theme(legend.key.size = unit(1.5, 'cm')) +
  scale_pattern_manual(
    values = c("none", "stripe"),
    guide = guide_legend(override.aes = list(fill = "gray")) # <- make lighter
  ) +
  scale_fill_manual(
    values = custom_colors, 
    breaks = legend_order,
    labels = c(
      "Radiation" = "Higher exposure levels\nof radiotherapy",
      "Chemotherapy" = "Higher exposure levels\nof chemotherapy"
    ),
    # guide = guide_legend(override.aes = list(pattern = "none")) # <- hide pattern
  ) +
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
  geom_text(data=label_data, aes(x=id, y=value+2, label=individual, hjust=hjust), color="black", fontface=1,alpha=1, size=6.5, angle= label_data$angle, inherit.aes = FALSE ) +
  
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
    linewidth = 2.5,
    size = 16 / .pt,
    inherit.aes = FALSE,
    gap = FALSE,
    offset = unit(-24, "pt")
  ) +
  # guides(
  #   theta = guide_axis_theta(angle = 0)
  # ) +
  theme(
    legend.position="right",
    legend.justification="right", 
    legend.box.spacing = unit(-300, "pt"),# The spacing between the plotting area and the legend box (unit)
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
    fill = "Risk factors"
  )

p_sjlife


## CCSSS
data <- data.all[data.all$Cohort == "CCSS",]

data$Cohort <- factor(data$Cohort, levels= c("CCSS"))
data$Risks <- factor(data$Risks, levels= c("Radiation", "Chemotherapy"))

custom_colors <- c(Radiation = "#9FB6CD", Chemotherapy = "#9FB6CD")
legend_order <- c("Radiation", "Chemotherapy")

# Set a number of 'empty bar' to add at the end of each group
empty_bar <- 1
to_add <- data.frame( matrix(0, empty_bar*nlevels(data$group), ncol(data)) )
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

## Fix issues with the gaps
data$value[is.na(data$Cohort)] <- NA
data$individual[is.na(data$Cohort)] <- NA
label_data$value[is.na(label_data$Cohort)] <- NA
label_data$individual[is.na(label_data$Cohort)] <- NA
data$Cohort[is.na(data$Cohort)] <- "CCSS"


# Make the plot
p_ccss = ggplot(data, aes(x=as.factor(id), y=value, fill = Risks)) +  
  geom_col(
    position = "dodge",
    width = 0.95
  ) +
  geom_col_pattern(
    data = subset(data, Risks == "Radiation"),
    aes(fill = Risks),
    pattern = "stripe",
    pattern_fill = "black",
    pattern_angle = 45,
    pattern_density = 0.1,
    pattern_spacing = 0.01,
    position = position_dodge2(preserve = 'single'),
    na.rm = TRUE
  ) +
  scale_fill_discrete(na.translate = F) +
  # theme(legend.key.size = unit(1.5, 'cm')) +
  scale_pattern_manual(
    values = c("none", "stripe"),
    guide = guide_legend(override.aes = list(fill = "gray")) # <- make lighter
  ) +
  scale_fill_manual(
    values = custom_colors, 
    breaks = legend_order,
    labels = c(
      "Radiation" = "Higher exposure levels\nof radiotherapy",
      "Chemotherapy" = "Higher exposure levels\nof chemotherapy"
    ),
    # guide = guide_legend(override.aes = list(pattern = "none")) # <- hide pattern
  ) +
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
  geom_text(data=label_data, aes(x=id, y=value+2, label=individual, hjust=hjust), color="black", fontface=1,alpha=1, size=6.5, angle= label_data$angle, inherit.aes = FALSE ) +
  
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
    linewidth = 2.5,
    size = 16 / .pt,
    inherit.aes = FALSE,
    gap = FALSE,
    offset = unit(-24, "pt")
  ) +
  # guides(
  #   theta = guide_axis_theta(angle = 0)
  # ) +
  theme(
    legend.position="right",
    legend.justification="right", 
    legend.box.spacing = unit(-300, "pt"),# The spacing between the plotting area and the legend box (unit)
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
    fill = "Risk factors"
  )

p_ccss

# library(cowplot)
# cowplot::plot_grid(p_sjlife, p_ccss, labels = c("SJLIFE", "CCSS"))
# p = cowplot::plot_grid(p_sjlife, p_ccss, labels = c("AUTO"))

# p_sjlife <- p_sjlife + theme(plot.margin = margin(r = 0, l = 0, t = 0, b = 0, "cm"))
# p_ccss <- p_ccss + theme(plot.margin = margin(r = 0, l = 0, t = 0, b = 0, "cm"))

p_sjlife <- p_sjlife + theme(plot.margin = margin(r = -15, l = -2, t = 0, b = 0, "cm"))
p_ccss <- p_ccss + theme(plot.margin = margin(r = -5, l = -10, t = 0, b = 0, "cm"))

library(ggplot2)
library(gridExtra)

p_sjlife <- p_sjlife + 
  theme(legend.position = "bottom") +
  guides(fill = guide_legend(override.aes = list(fill = "gray"))) # Set fill color to white
# Remove legend from the other plot
p_ccss <- p_ccss + theme(legend.position = "none")



p_sjlife_label <- p_sjlife +
  annotate("text", x = Inf, y = -Inf, label = "SJLIFE", vjust = 0, hjust = 1, size = 6)

# Add label "B" to p_ccss
p_ccss_label <- p_ccss +
  annotate("text", x = Inf, y = -Inf, label = "CCSS", vjust = 0, hjust = 1, size = 6)

# Combine the labeled plots
p <- ggarrange(
  p_sjlife_label, p_ccss_label, 
  ncol = 2, nrow = 1, 
  common.legend = TRUE, legend = "right",
  align = "hv"
)

p

p = p + theme(plot.margin = margin(r = 0, l = 0, t = 0, b = 0, "cm"))
p
plot_name <- "Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/ANALYSIS/results/plots/v20b/overall/Overall_Chemotherapy_and_Radiation_without_lifestyle_radial_plot_faceted.tiff"
ggsave(plot_name, p, width = 18, height = 12, dpi = 600, device = "tiff", compression = "lzw")
