# Load required libraries
library(ggplot2)
library(ggthemes)# For additional themes
library(ggrepel)
library(dplyr)
library(ggpattern)
# https://r-charts.com/colors/

## v21b (Without lifestyle)
data <- read.table(text="Cohorts	SN_types	Variables	Overall	Female	Male	Age.lt.35	Age.ge.35	EUR	AFR
SJLIFE	SNs (605)	Radiation	43	42	44	40	45	43	36
SJLIFE	SNs (605)	Chemotherapy	8	8	8	11	5	8	7
SJLIFE	SNs (605)	All_treatments	47	46	48	46	47	48	40
SJLIFE	SNs (605)	PRS	13	12	13	13	13	13	14
SJLIFE	SNs (605)	Lifestyle	NA	NA	NA	NA	NA	NA	NA
SJLIFE	SNs (605)	Combined	54	53	55	53	54	54	48
SJLIFE	SMNs (463)	Radiation	37	37	37	34	39	37	33
SJLIFE	SMNs (463)	Chemotherapy	3	3	3	5	2	3	3
SJLIFE	SMNs (463)	All_treatments	39	39	39	37	41	39	35
SJLIFE	SMNs (463)	PRS	12	12	12	12	12	12	13
SJLIFE	SMNs (463)	Lifestyle	NA	NA	NA	NA	NA	NA	NA
SJLIFE	SMNs (463)	Combined	46	46	46	44	48	47	43
SJLIFE	NMSCs (251)	Radiation	44	43	45	41	45	44	37
SJLIFE	NMSCs (251)	Chemotherapy	NA	NA	NA	NA	NA	NA	NA
SJLIFE	NMSCs (251)	All_treatments	44	43	45	41	45	44	37
SJLIFE	NMSCs (251)	PRS	28	27	30	28	28	28	2
SJLIFE	NMSCs (251)	Lifestyle	NA	NA	NA	NA	NA	NA	NA
SJLIFE	NMSCs (251)	Combined	59	58	61	57	61	60	38
SJLIFE	Breast cancer (76)	Radiation	51	51	NA	49	51	53	40
SJLIFE	Breast cancer (76)	Chemotherapy	19	19	NA	22	18	19	20
SJLIFE	Breast cancer (76)	All_treatments	61	61	NA	60	62	63	53
SJLIFE	Breast cancer (76)	PRS	19	19	NA	20	19	19	18
SJLIFE	Breast cancer (76)	Lifestyle	NA	NA	NA	NA	NA	NA	NA
SJLIFE	Breast cancer (76)	Combined	69	69	NA	68	69	70	61
SJLIFE	Thyroid cancer (87)	Radiation	62	62	62	59	65	62	59
SJLIFE	Thyroid cancer (87)	Chemotherapy	24	24	23	29	17	24	21
SJLIFE	Thyroid cancer (87)	All_treatments	73	73	73	73	73	73	69
SJLIFE	Thyroid cancer (87)	PRS	52	52	51	52	52	52	49
SJLIFE	Thyroid cancer (87)	Lifestyle	NA	NA	NA	NA	NA	NA	NA
SJLIFE	Thyroid cancer (87)	Combined	87	87	87	87	87	87	84
SJLIFE	Meningioma (149)	Radiation	91	91	92	90	92	92	88
SJLIFE	Meningioma (149)	Chemotherapy	24	24	24	29	21	24	24
SJLIFE	Meningioma (149)	All_treatments	93	93	93	93	93	93	90
SJLIFE	Meningioma (149)	PRS	NA	NA	NA	NA	NA	NA	NA
SJLIFE	Meningioma (149)	Lifestyle	NA	NA	NA	NA	NA	NA	NA
SJLIFE	Meningioma (149)	Combined	92	92	92	92	92	92	88
SJLIFE	Sarcoma (33)	Radiation	NA	NA	NA	NA	NA	NA	NA
SJLIFE	Sarcoma (33)	Chemotherapy	35	35	34	33	37	35	32
SJLIFE	Sarcoma (33)	All_treatments	35	35	34	33	37	35	32
SJLIFE	Sarcoma (33)	PRS	NA	NA	NA	NA	NA	NA	NA
SJLIFE	Sarcoma (33)	Lifestyle	NA	NA	NA	NA	NA	NA	NA
SJLIFE	Sarcoma (33)	Combined	30	30	29	28	32	30	28
CCSS	SNs (1611)	Radiation	39	38	40	37	40	39	40
CCSS	SNs (1611)	Chemotherapy	3	3	3	5	2	3	6
CCSS	SNs (1611)	All_treatments	41	40	42	40	41	41	43
CCSS	SNs (1611)	PRS	5	5	5	5	5	5	4
CCSS	SNs (1611)	Lifestyle	NA	NA	NA	NA	NA	NA	NA
CCSS	SNs (1611)	Combined	44	43	45	43	44	44	46
CCSS	SMNs (762)	Radiation	26	25	26	21	28	26	20
CCSS	SMNs (762)	Chemotherapy	4	4	4	7	2	4	9
CCSS	SMNs (762)	All_treatments	29	29	29	27	30	29	27
CCSS	SMNs (762)	PRS	5	5	5	5	5	5	6
CCSS	SMNs (762)	Lifestyle	NA	NA	NA	NA	NA	NA	NA
CCSS	SMNs (762)	Combined	32	32	33	31	34	33	31
CCSS	NMSCs (774)	Radiation	39	39	40	39	39	39	43
CCSS	NMSCs (774)	Chemotherapy	NA	NA	NA	NA	NA	NA	NA
CCSS	NMSCs (774)	All_treatments	39	39	40	39	39	39	43
CCSS	NMSCs (774)	PRS	31	30	31	31	31	31	8
CCSS	NMSCs (774)	Lifestyle	NA	NA	NA	NA	NA	NA	NA
CCSS	NMSCs (774)	Combined	58	58	58	58	58	58	46
CCSS	Breast cancer (290)	Radiation	48	48	NA	45	48	48	36
CCSS	Breast cancer (290)	Chemotherapy	19	19	NA	22	18	18	27
CCSS	Breast cancer (290)	All_treatments	59	59	NA	59	59	59	54
CCSS	Breast cancer (290)	PRS	37	37	NA	37	37	36	36
CCSS	Breast cancer (290)	Lifestyle	NA	NA	NA	NA	NA	NA	NA
CCSS	Breast cancer (290)	Combined	74	74	NA	74	74	74	70
CCSS	Thyroid cancer (163)	Radiation	44	44	44	41	49	45	39
CCSS	Thyroid cancer (163)	Chemotherapy	6	5	6	7	3	5	8
CCSS	Thyroid cancer (163)	All_treatments	48	48	47	46	50	48	44
CCSS	Thyroid cancer (163)	PRS	36	36	35	36	36	36	35
CCSS	Thyroid cancer (163)	Lifestyle	NA	NA	NA	NA	NA	NA	NA
CCSS	Thyroid cancer (163)	Combined	66	66	66	65	68	67	63
CCSS	Meningioma (256)	Radiation	86	86	88	87	86	86	90
CCSS	Meningioma (256)	Chemotherapy	8	7	10	11	5	8	16
CCSS	Meningioma (256)	All_treatments	88	87	89	89	87	88	91
CCSS	Meningioma (256)	PRS	1	1	1	2	1	1	NA
CCSS	Meningioma (256)	Lifestyle	NA	NA	NA	NA	NA	NA	NA
CCSS	Meningioma (256)	Combined	88	87	89	89	87	88	91
CCSS	Sarcoma (61)	Radiation	NA	NA	NA	NA	NA	NA	NA
CCSS	Sarcoma (61)	Chemotherapy	35	34	37	35	36	35	35
CCSS	Sarcoma (61)	All_treatments	35	34	37	35	36	35	35
CCSS	Sarcoma (61)	PRS	4	4	4	4	4	4	3
CCSS	Sarcoma (61)	Lifestyle	NA	NA	NA	NA	NA	NA	NA
CCSS	Sarcoma (61)	Combined	38	37	39	38	39	38	38", header = T, sep = "\t")


data$SN_types <- trimws(gsub("\\([0-9]+\\)", "", data$SN_types))
data[data == "-"] <- NA
data <- data %>%
  dplyr::mutate(across(where(is.numeric), ~ . / 100))
saved.data <- data

## 1.___________________________________________________________________________ Combined (overall) analysis
## 1.___________________________________________________________________________ Combined (overall) analysis
## 1.___________________________________________________________________________ Combined (overall) analysis
## 1.___________________________________________________________________________ Combined (overall) analysis
## 1.___________________________________________________________________________ Combined (overall) analysis

## Radial plot
group <- "Overall"

variables <- "Combined"

new_data <- saved.data[grepl(variables, saved.data$Variables), c("Cohorts", "SN_types", "Variables", group)]
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

data <- cbind.data.frame(Cohorts=data_melted$variable, Risks=gsub("All_t", "T", data_melted$AF_by), individual= paste0(data_melted$new_value, "%"), group= data_melted$SN_types, group2= data_melted$legend_group, value=data_melted$new_value)
# data$value[is.na(data$value)]
data$Cohorts <- factor(data$Cohorts, levels= c("SJLIFE", "CCSS"))
data$Risks <- factor(data$Risks, levels= c("Treatments", "PRS"))

custom_colors <- c("SJLIFE" = "#EE6A50", "CCSS" = "#9FB6CD")
legend_order <- c("SJLIFE", "CCSS")


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

# Make the plot
p = ggplot(data, aes(x=as.factor(id), y=value, fill=Cohorts)) +       # Note that id is a factor. If x is numeric, there is some space between the first bar
  
  # geom_bar(aes(x=as.factor(id), y=value, fill=group2), stat="identity", alpha=1, width=0.95) +

  # # Add a val=100/75/50/25 lines. I do it at the beginning to make sur barplots are OVER it.
  # geom_segment(data=grid_data, aes(x = end, y = 80, xend = start, yend = 80), colour = "white", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  # geom_segment(data=grid_data, aes(x = end, y = 60, xend = start, yend = 60), colour = "white", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  # geom_segment(data=grid_data, aes(x = end, y = 40, xend = start, yend = 40), colour = "white", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  # geom_segment(data=grid_data, aes(x = end, y = 20, xend = start, yend = 20), colour = "white", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  
  ## Add text showing the value of each 100/75/50/25 lines
  # annotate("text", x = rep(max(data$id),5), y = c(0, 25, 50, 75, 100), label = c("0%", "25%", "50%", "75%", "100%") , color="black", size=3 , angle=0, fontface= 1, hjust=1) +
  
  geom_bar(aes(x=as.factor(id), y=value, fill=Cohorts), stat="identity", alpha=1, width=0.95) +
  theme_bw() +
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
  scale_fill_manual(values = custom_colors, breaks = legend_order) +
  guides(
    theta = guide_axis_theta(angle = 0)
  ) +
  theme(
    legend.position="right",
    legend.justification="right", 
    legend.box.spacing = unit(-120, "pt"),# The spacing between the plotting area and the legend box (unit)
    legend.margin=margin(0,0,0,0),
    # axis.text = element_text(size = 18, color = "black"), # Adjust axis text size
    # axis.title = element_text(size = 18, color = "black"), # Adjust axis title size
    legend.title = element_text(size = 16, color = "black", face=2), # Adjust legend title size
    plot.margin = margin(r = .5, l = -3.5, t = -3.5, b = -3.5, unit = "cm"),
    legend.text = element_text(size = 14, color = "black")
  ) + # Adjust legend text size
  labs(
    title = "",
    x = NULL, y = "",
    fill = "Cohorts"
  )
p
plot_name <- "Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/ANALYSIS/results/plots/v21b/overall/Overall_Combined_without_lifestyle_radial_plot.tiff"
ggsave(plot_name, p, width = 12, height = 8, dpi = 600, device = "tiff", compression = "lzw")


################################
## PRS and treatment together ##
################################
## Radial plot
group <- "Overall"

variables <- "treatment|PRS"

new_data <- saved.data[grepl(variables, saved.data$Variables), c("Cohorts", "SN_types", "Variables", group)]
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

data <- cbind.data.frame(Cohorts=data_melted$variable, Risks=gsub("All_t", "T", data_melted$AF_by), individual= paste0(data_melted$new_value, "%"), group= data_melted$SN_types, group2= data_melted$legend_group, value=data_melted$new_value)
# data$value[is.na(data$value)]
data$Cohorts <- factor(data$Cohorts, levels= c("SJLIFE", "CCSS"))
data$Risks <- factor(data$Risks, levels= c("Treatments", "PRS"))

custom_colors <- c("Treatments" = "#FF7F50", "PRS" = "#A2B5CD")
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
data$value[is.na(data$Cohorts)] <- NA
data$individual[is.na(data$Cohorts)] <- NA
label_data$value[is.na(label_data$Cohorts)] <- NA
label_data$individual[is.na(label_data$Cohorts)] <- NA
data$Cohorts[is.na(data$Cohorts)] <- c("SJLIFE", "CCSS", "SJLIFE", "CCSS", "SJLIFE", "CCSS", "SJLIFE")


# Make the plot
p = ggplot(data, aes(x=as.factor(id), y=value, fill = Risks)) +  
  # # Add a val=100/75/50/25 lines. I do it at the beginning to make sur barplots are OVER it.
  # geom_segment(data=grid_data, aes(x = end, y = 80, xend = start, yend = 80), colour = "white", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  # geom_segment(data=grid_data, aes(x = end, y = 60, xend = start, yend = 60), colour = "white", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  # geom_segment(data=grid_data, aes(x = end, y = 40, xend = start, yend = 40), colour = "white", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  # geom_segment(data=grid_data, aes(x = end, y = 20, xend = start, yend = 20), colour = "white", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  
  ## Add text showing the value of each 100/75/50/25 lines
  # annotate("text", x = rep(max(data$id),5), y = c(0, 25, 50, 75, 100), label = c("0%", "25%", "50%", "75%", "100%") , color="black", size=3 , angle=0, fontface= 1, hjust=1) +
  
  geom_bar(aes(x=as.factor(id), y=value, fill=Risks), stat="identity", alpha=1, width=0.95) +
  geom_col_pattern(
    aes(pattern = Cohorts),
    colour = "black",
    pattern_fill = "black",
    pattern_angle = 45,
    pattern_density = 0.1,
    pattern_spacing = 0.01,
    position = position_dodge2(preserve = 'single'),
    na.rm = T
  ) +
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
    guide = guide_legend(override.aes = list(pattern = "none")) # <- hide pattern
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
    legend.box.spacing = unit(-120, "pt"),# The spacing between the plotting area and the legend box (unit)
    legend.margin=margin(0,0,0,0),
    # axis.text = element_text(size = 18, color = "black"), # Adjust axis text size
    # axis.title = element_text(size = 18, color = "black"), # Adjust axis title size
    legend.title = element_text(size = 16, color = "black", face=2), # Adjust legend title size
    plot.margin = margin(r = .5, l = -3.5, t = -3.5, b = -3.5, unit = "cm"),
    legend.text = element_text(size = 14, color = "black")
  ) + # Adjust legend text size
  labs(
    title = "",
    x = NULL, y = "",
    fill = "Risk factors"
  )

# p = p  + scale_color_brewer(palette = 2)
p
plot_name <- "Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/ANALYSIS/results/plots/v21b/overall/Overall_treatment_and_PRS_without_lifestyle_radial_plot.tiff"
ggsave(plot_name, p, width = 12, height = 8, dpi = 600, device = "tiff", compression = "lzw")

#########################################
## Chemotherapy and Radiation together ##
#########################################
## Radial plot
group <- "Overall"

variables <- "Chemotherapy|Radiation"

new_data <- saved.data[grepl(variables, saved.data$Variables), c("Cohorts", "SN_types", "Variables", group)]
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

data <- cbind.data.frame(Cohorts=data_melted$variable, Risks=gsub("All_t", "T", data_melted$AF_by), individual= paste0(data_melted$new_value, "%"), group= data_melted$SN_types, group2= data_melted$legend_group, value=data_melted$new_value)
# data$value[is.na(data$value)]
data$Cohorts <- factor(data$Cohorts, levels= c("SJLIFE", "CCSS"))
data$Risks[data$Risks=="Radiation"] <- "Radiotherapies"
data$Risks[data$Risks=="Chemotherapy"] <- "Chemotherapies"
data$Risks <- factor(data$Risks, levels= c("Radiotherapies", "Chemotherapies"))

custom_colors <- c("Radiotherapies" = "#FF7F50", "Chemotherapies" = "#A2B5CD")
legend_order <- c("Radiotherapies", "Chemotherapies")

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
data$value[is.na(data$Cohorts)] <- NA
data$individual[is.na(data$Cohorts)] <- NA
label_data$value[is.na(label_data$Cohorts)] <- NA
label_data$individual[is.na(label_data$Cohorts)] <- NA
data$Cohorts[is.na(data$Cohorts)] <- c("SJLIFE", "CCSS", "SJLIFE", "CCSS", "SJLIFE", "CCSS", "SJLIFE")

# Make the plot
p = ggplot(data, aes(x=as.factor(id), y=value, fill = Risks)) +  
  # # Add a val=100/75/50/25 lines. I do it at the beginning to make sur barplots are OVER it.
  # geom_segment(data=grid_data, aes(x = end, y = 80, xend = start, yend = 80), colour = "white", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  # geom_segment(data=grid_data, aes(x = end, y = 60, xend = start, yend = 60), colour = "white", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  # geom_segment(data=grid_data, aes(x = end, y = 40, xend = start, yend = 40), colour = "white", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  # geom_segment(data=grid_data, aes(x = end, y = 20, xend = start, yend = 20), colour = "white", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  
  ## Add text showing the value of each 100/75/50/25 lines
  # annotate("text", x = rep(max(data$id),5), y = c(0, 25, 50, 75, 100), label = c("0%", "25%", "50%", "75%", "100%") , color="black", size=3 , angle=0, fontface= 1, hjust=1) +
  
  geom_bar(aes(x=as.factor(id), y=value, fill=Risks), stat="identity", alpha=1, width=0.95) +
  geom_col_pattern(
    aes(pattern = Cohorts),
    colour = "black",
    pattern_fill = "black",
    pattern_angle = 45,
    pattern_density = 0.1,
    pattern_spacing = 0.01,
    position = position_dodge2(preserve = 'single'),
    na.rm = T
  ) +
  scale_pattern_manual(
    values = c("none", "stripe"),
    guide = guide_legend(override.aes = list(fill = "gray")) # <- make lighter
  ) +
  scale_fill_manual(
    values = custom_colors, 
    breaks = legend_order,
    labels = c(
      "Radiotherapies" = "Higher exposure levels\nof radiotherapy",
      "Chemotherapies" = "Higher exposure levels\nof chemotherapy"
    ),
    guide = guide_legend(override.aes = list(pattern = "none")) # <- hide pattern
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
    legend.box.spacing = unit(-120, "pt"),# The spacing between the plotting area and the legend box (unit)
    legend.margin=margin(0,0,0,0),
    # axis.text = element_text(size = 18, color = "black"), # Adjust axis text size
    # axis.title = element_text(size = 18, color = "black"), # Adjust axis title size
    legend.title = element_text(size = 16, color = "black", face=2), # Adjust legend title size
    plot.margin = margin(r = .5, l = -3.5, t = -3.5, b = -3.5, unit = "cm"),
    legend.text = element_text(size = 14, color = "black")
  ) + # Adjust legend text size
  labs(
    title = "",
    x = NULL, y = "",
    fill = "Risk factors"
  )

# p = p  + scale_color_brewer(palette = 2)
p

plot_name <- "Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/ANALYSIS/results/plots/v21b/overall/Overall_Chemotherapy_and_Radiation_without_lifestyle_radial_plot.tiff"
ggsave(plot_name, p, width = 12, height = 8, dpi = 600, device = "tiff", compression = "lzw")

