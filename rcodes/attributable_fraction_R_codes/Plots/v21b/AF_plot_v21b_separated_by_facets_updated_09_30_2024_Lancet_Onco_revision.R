# Load required libraries
library(ggplot2)
library(ggthemes)# For additional themes
library(ggrepel)
library(dplyr)
library(ggpattern)
library("ggpubr")
library(RColorBrewer)
# https://r-charts.com/colors/

## Please add CI to the Overall column!!!!!!!!!!!

## v21b (Without lifestyle)
data <- read.table(text="Cohort	SN_types	Variables	Overall	Female	Male	Age.lt.35	Age.ge.35	EUR	AFR
SJLIFE	SNs (605)	Radiotherapy	43	42	44	40	45	43	36
SJLIFE	SNs (605)	Chemotherapy	8	8	8	11	5	8	7
SJLIFE	SNs (605)	Treatments	47	46	48	46	47	48	40
SJLIFE	SNs (605)	PRS	13	12	13	13	13	13	14
SJLIFE	SNs (605)	Lifestyle	NA	NA	NA	NA	NA	NA	NA
SJLIFE	SNs (605)	Combined	54	53	55	53	54	54	48
SJLIFE	SMNs (463)	Radiotherapy	37	37	37	34	39	37	33
SJLIFE	SMNs (463)	Chemotherapy	3	3	3	5	2	3	3
SJLIFE	SMNs (463)	Treatments	39	39	39	37	41	39	35
SJLIFE	SMNs (463)	PRS	12	12	12	12	12	12	13
SJLIFE	SMNs (463)	Lifestyle	NA	NA	NA	NA	NA	NA	NA
SJLIFE	SMNs (463)	Combined	46	46	46	44	48	47	43
SJLIFE	NMSCs (251)	Radiotherapy	44	43	45	41	45	44	37
SJLIFE	NMSCs (251)	Chemotherapy	NA	NA	NA	NA	NA	NA	NA
SJLIFE	NMSCs (251)	Treatments	44	43	45	41	45	44	37
SJLIFE	NMSCs (251)	PRS	28	27	30	28	28	28	2
SJLIFE	NMSCs (251)	Lifestyle	NA	NA	NA	NA	NA	NA	NA
SJLIFE	NMSCs (251)	Combined	59	58	61	57	61	60	38
SJLIFE	Breast cancer (76)	Radiotherapy	51	51	NA	49	51	53	40
SJLIFE	Breast cancer (76)	Chemotherapy	19	19	NA	22	18	19	20
SJLIFE	Breast cancer (76)	Treatments	61	61	NA	60	62	63	53
SJLIFE	Breast cancer (76)	PRS	19	19	NA	20	19	19	18
SJLIFE	Breast cancer (76)	Lifestyle	NA	NA	NA	NA	NA	NA	NA
SJLIFE	Breast cancer (76)	Combined	69	69	NA	68	69	70	61
SJLIFE	Thyroid cancer (87)	Radiotherapy	62	62	62	59	65	62	59
SJLIFE	Thyroid cancer (87)	Chemotherapy	24	24	23	29	17	24	21
SJLIFE	Thyroid cancer (87)	Treatments	73	73	73	73	73	73	69
SJLIFE	Thyroid cancer (87)	PRS	52	52	51	52	52	52	49
SJLIFE	Thyroid cancer (87)	Lifestyle	NA	NA	NA	NA	NA	NA	NA
SJLIFE	Thyroid cancer (87)	Combined	87	87	87	87	87	87	84
SJLIFE	Meningioma (149)	Radiotherapy	91	91	92	90	92	92	88
SJLIFE	Meningioma (149)	Chemotherapy	24	24	24	29	21	24	24
SJLIFE	Meningioma (149)	Treatments	93	93	93	93	93	93	90
SJLIFE	Meningioma (149)	PRS	NA	NA	NA	NA	NA	NA	NA
SJLIFE	Meningioma (149)	Lifestyle	NA	NA	NA	NA	NA	NA	NA
SJLIFE	Meningioma (149)	Combined	92	92	92	92	92	92	88
SJLIFE	Sarcoma (33)	Radiotherapy	NA	NA	NA	NA	NA	NA	NA
SJLIFE	Sarcoma (33)	Chemotherapy	35	35	34	33	37	35	32
SJLIFE	Sarcoma (33)	Treatments	35	35	34	33	37	35	32
SJLIFE	Sarcoma (33)	PRS	NA	NA	NA	NA	NA	NA	NA
SJLIFE	Sarcoma (33)	Lifestyle	NA	NA	NA	NA	NA	NA	NA
SJLIFE	Sarcoma (33)	Combined	30	30	29	28	32	30	28
CCSS	SNs (1611)	Radiotherapy	39	38	40	37	40	39	40
CCSS	SNs (1611)	Chemotherapy	3	3	3	5	2	3	6
CCSS	SNs (1611)	Treatments	41	40	42	40	41	41	43
CCSS	SNs (1611)	PRS	5	5	5	5	5	5	4
CCSS	SNs (1611)	Lifestyle	NA	NA	NA	NA	NA	NA	NA
CCSS	SNs (1611)	Combined	44	43	45	43	44	44	46
CCSS	SMNs (762)	Radiotherapy	26	25	26	21	28	26	20
CCSS	SMNs (762)	Chemotherapy	4	4	4	7	2	4	9
CCSS	SMNs (762)	Treatments	29	29	29	27	30	29	27
CCSS	SMNs (762)	PRS	5	5	5	5	5	5	6
CCSS	SMNs (762)	Lifestyle	NA	NA	NA	NA	NA	NA	NA
CCSS	SMNs (762)	Combined	32	32	33	31	34	33	31
CCSS	NMSCs (728)	Radiotherapy	42	41	42	41	42	42	45
CCSS	NMSCs (728)	Chemotherapy	NA	NA	NA	NA	NA	NA	NA
CCSS	NMSCs (728)	Treatments	42	41	42	41	42	42	45
CCSS	NMSCs (728)	PRS	33	33	33	33	33	34	8
CCSS	NMSCs (728)	Lifestyle	NA	NA	NA	NA	NA	NA	NA
CCSS	NMSCs (728)	Combined	61	61	61	61	61	62	49
CCSS	Breast cancer (290)	Radiotherapy	48	48	NA	45	48	48	36
CCSS	Breast cancer (290)	Chemotherapy	19	19	NA	22	18	18	27
CCSS	Breast cancer (290)	Treatments	59	59	NA	59	59	59	54
CCSS	Breast cancer (290)	PRS	37	37	NA	37	37	36	36
CCSS	Breast cancer (290)	Lifestyle	NA	NA	NA	NA	NA	NA	NA
CCSS	Breast cancer (290)	Combined	74	74	NA	74	74	74	70
CCSS	Thyroid cancer (163)	Radiotherapy	44	44	44	41	49	45	39
CCSS	Thyroid cancer (163)	Chemotherapy	6	5	6	7	3	5	8
CCSS	Thyroid cancer (163)	Treatments	48	48	47	46	50	48	44
CCSS	Thyroid cancer (163)	PRS	36	36	35	36	36	36	35
CCSS	Thyroid cancer (163)	Lifestyle	NA	NA	NA	NA	NA	NA	NA
CCSS	Thyroid cancer (163)	Combined	66	66	66	65	68	67	63
CCSS	Meningioma (256)	Radiotherapy	86	86	88	87	86	86	90
CCSS	Meningioma (256)	Chemotherapy	8	7	10	11	5	8	16
CCSS	Meningioma (256)	Treatments	88	87	89	89	87	88	91
CCSS	Meningioma (256)	PRS	1	1	1	2	1	1	NA
CCSS	Meningioma (256)	Lifestyle	NA	NA	NA	NA	NA	NA	NA
CCSS	Meningioma (256)	Combined	88	87	89	89	87	88	91
CCSS	Sarcoma (61)	Radiotherapy	NA	NA	NA	NA	NA	NA	NA
CCSS	Sarcoma (61)	Chemotherapy	35	34	37	35	36	35	35
CCSS	Sarcoma (61)	Treatments	35	34	37	35	36	35	35
CCSS	Sarcoma (61)	PRS	4	4	4	4	4	4	3
CCSS	Sarcoma (61)	Lifestyle	NA	NA	NA	NA	NA	NA	NA
CCSS	Sarcoma (61)	Combined	38	37	39	38	39	38	38", header = T, sep = "\t")


data$SN_types <- trimws(gsub("\\([0-9]+\\)", "", data$SN_types))
data[data == "-"] <- NA
data <- data %>%
  dplyr::mutate(across(where(is.numeric), ~ . / 100))

data <- data[!grepl("Lifestyle", data$Variables),]
saved.data <- data

##############
## Figure 1 ##
##############
## 1.___________________________________________________________________________ Combined (overall) analysis
## 1.___________________________________________________________________________ Combined (overall) analysis
## 1.___________________________________________________________________________ Combined (overall) analysis
## 1.___________________________________________________________________________ Combined (overall) analysis
## 1.___________________________________________________________________________ Combined (overall) analysis


# Reshape the data for ggplot2
library(reshape2)

group.all <- c("Overall", "Female", "Male", "Age.lt.35", "Age.ge.35", "EUR", "AFR")

sn_types <- unique(data$SN_types)

j=1
group <- group.all[j]
variables <- unique(data$Variables)
order = c("Combined", "Treatments", "Radiotherapy", "Chemotherapy", "PRS")

custom_colors <- brewer.pal(5, "Dark2")
legend_order <- order


## combined all SNs
data_melted <- data[, c("Cohort", "SN_types", "Variables", group)]
data_melted <- data_melted[grepl("Combined", data_melted$Variables),]

data_melted$value <- as.numeric(data_melted[,group]*100)
data_melted$legend_group <- factor(data_melted$Cohort, levels = c("SJLIFE", "CCSS"))

# Create a factor with the desired order
data_melted$Variables <- factor(data_melted$Variables, levels = order)

Cohort <- c("SJLIFE", "CCSS")
data_melted$Cohort <- factor(data_melted$Cohort, levels = Cohort)
data_melted <- data_melted[complete.cases(data_melted),]
data_melted$new_value <- paste0(round(data_melted$value,2), "%")

data_melted$SN_types <- factor(data_melted$SN_types, levels = unique(data_melted$SN_types))

p <- ggplot(data_melted, aes(x = SN_types, y = value, fill = legend_group)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.6), width = 0.6) +
  geom_text(
    data = data_melted %>% filter(!is.na(new_value)),
    aes(label = new_value, y = value),  # Adjust y position
    position = position_dodge(width = 0.6),
    vjust = -0.30,  # Adjust vertical justification
    hjust = 0.5,  # Center text horizontally
    size = 4.5,
    color = "black"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1, size = 16, color = "black"),  # Tilt x-axis labels 45 degrees
    axis.text.y = element_text(size = 16, color = "black"),  # Make y-axis text black
    axis.title.y = element_text(size = 20, color = "black"),  # Make y-axis title black
    axis.title.x = element_blank(),
    axis.line.x = element_line(color = "black"),  # Make x-axis black
    axis.line.y = element_line(color = "black"),  # Make y-axis black
    legend.title = element_text(size = 22, color = "black", face = "bold"),  # Increase legend title size
    legend.text = element_text(size = 18, color = "black"),
    legend.position = "top",
    legend.box = "horizontal",
    legend.margin = margin(t = 10, b = 10, l = 10, r = 200),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.key = element_rect(color = "black", size = 1),  # Add border around legend
    # legend.background = element_rect(fill = "white", color = "black", size = 1)  # Add background to legend
    legend.background = element_blank()  # Remove background and borders from legend box
  ) +
  labs(fill = "Cohorts: ") +
  coord_cartesian(clip = "off") +
  scale_y_continuous(limits = c(0, 100), expand = c(0, 0)) +
  labs(title = "", y = "Attributable fraction (%)", x = NULL) +
  scale_fill_manual(
    values = c("SJLIFE" = "#D95F02", "CCSS" = "#7570B3"),  # Custom colors: burgundy for SJLIFE, dark red for CCSS
    breaks = levels(data_melted$legend_group),  # Ensure breaks match the levels of legend_group
)

p
plot_name <- paste0("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/ANALYSIS/results/plots/v21b/combined_all_SNs_without_lifestyle_plot.tiff")
ggsave(plot_name, p, width = 10, height = 9, dpi = 600, device = "tiff", compression = "lzw")

##############
## Figure 2 ##
##############
# SN and SMN
data_melted <- data[, c("Cohort", "SN_types", "Variables", group)]
data_melted <- data_melted[grepl(paste0(unique(data_melted$SN_types)[1:2], collapse="|"), data_melted$SN_types),]

data_melted$value <- as.numeric(data_melted[,group]*100)
data_melted$legend_group <- factor(data_melted$Variables, levels = order)

# Create a factor with the desired order
data_melted$Variables <- factor(data_melted$Variables, levels = order)

Cohort <- c("SJLIFE", "CCSS")
data_melted$Cohort <- factor(data_melted$Cohort, levels = Cohort)
data_melted <- data_melted[complete.cases(data_melted),]
data_melted$new_value <- paste0(round(data_melted$value,2), "%")

data_melted$SN_types <- factor(data_melted$SN_types, levels = unique(data_melted$SN_types)[1:2])

# Create the plot
p <- ggplot(data_melted, aes(x = Cohort, y = value, fill = legend_group)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9), width = 0.8) +
  geom_text(
    data = data_melted %>% filter(!is.na(new_value)),
    aes(label = new_value, y = value),  # Adjust y position
    position = position_dodge(width = 0.9),
    vjust = -0.20,  # Adjust vertical justification
    hjust = 0.5,  # Center text horizontally
    size = 6.5,
    color = "black"
  ) +
  geom_vline(xintercept = 1.5, linetype = "dashed", color = "black", size = 1) +  # Add vertical line
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 1, vjust = 1, size = 24, color = "black"),
    axis.text.y = element_text(size = 20, color = "black"),
    axis.title.y = element_text(size = 24, color = "black"),
    axis.title.x = element_blank(),
    axis.line.x = element_line(color = "black"),
    axis.line.y = element_line(color = "black"),
    legend.title = element_text(size = 22, color = "black", face = "bold"),  # Increase legend title size
    legend.text = element_text(size = 18, color = "black"),
    legend.position = "top",
    legend.box = "horizontal",
    legend.margin = margin(t = 10, b = 10, l = 20, r = 180),
    plot.margin = margin(t = 20, b = 80, l = 20, r = 20),  # Adjust plot margins
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.key = element_rect(color = "black", size = 1),  # Add border around legend
    legend.background = element_rect(fill = "white", color = "black", size = 1),  # Add background to legend
    strip.text = element_blank(),  # Remove facet labels (SN types)
    panel.spacing = unit(3, "lines")  # Increase space between facets
  ) +
  # Customize colors
  scale_fill_manual(
    values = custom_colors,
    breaks = legend_order,
    labels = c("Higher exposure levels of cancer treatments + Elevated genetic risks", 
               "Higher exposure levels of cancer treatments", 
               "Higher exposure levels of radiotherapy", 
               "Higher exposure levels of chemotherapy", 
               "Elevated genetic risks")
  ) +
  labs(fill = "Risk factors:") +
  scale_x_discrete(expand = c(0, 0), labels = c("SJLIFE", "CCSS")) +  # Set custom x-axis labels
  scale_y_continuous(limits = c(0, 100), breaks = c(0, 50, 100), expand = c(0, 0)) +
  labs(title = "", y = "Attributable fraction (%)", x = NULL) +
  guides(fill = guide_legend(
    title.position = "top",  # Position the legend title at the top
    title.hjust = 0,  # Center align the legend title
    label.hjust = 0,  # Left align the legend labels
    ncol = 1  # Adjust number of columns in legend
  )) +
  facet_wrap(~ SN_types, ncol = 1, scales='free') +  # Add faceting by SN_types
  # Adding the SN_types label at the top right of each facet
  geom_text(data = unique(data_melted %>% select(SN_types)),
            aes(x = Inf, y = Inf, label = SN_types),  # Place at the top right
            size = 8,
            fontface = "bold",
            hjust = 1,
            vjust = 2,
            color = "black",
            inherit.aes = FALSE) +
  coord_cartesian(clip = "off")  # Prevent clipping of text

p
plot_name <- paste0("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/ANALYSIS/results/plots/v21b/combined_SNs_and_SMNs_without_lifestyle_plot.tiff")
ggsave(plot_name, p, width = 11.6, height = 10, dpi = 600, device = "tiff", compression = "lzw")

##############
## Figure 3 ## 
##############
## "NMSCs",  "Breast cancer",  "Thyroid cancer", "Meningioma", "Sarcoma" 
data_melted <- data[, c("Cohort", "SN_types", "Variables", group)]
data_melted <- data_melted[grepl(paste0(unique(data_melted$SN_types)[-c(1:2)], collapse="|"), data_melted$SN_types),]

data_melted$value <- as.numeric(data_melted[,group]*100)
data_melted$legend_group <- factor(data_melted$Variables, levels = order)

# Create a factor with the desired order
data_melted$Variables <- factor(data_melted$Variables, levels = order)

Cohort <- c("SJLIFE", "CCSS")
data_melted$Cohort <- factor(data_melted$Cohort, levels = Cohort)
data_melted <- data_melted[complete.cases(data_melted),]
data_melted$new_value <- paste0(round(data_melted$value,2), "%")

data_melted$SN_types <- factor(data_melted$SN_types, levels = unique(data_melted$SN_types))

p <- ggplot(data_melted, aes(x = Cohort, y = value, fill = legend_group)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9), width = 0.8) +
  geom_text(
    data = data_melted %>% filter(!is.na(new_value)),
    aes(label = new_value, y = value),  # Adjust y position
    position = position_dodge(width = 0.9),
    vjust = -0.20,  # Adjust vertical justification
    hjust = 0.5,  # Center text horizontally
    size = 6.5,
    color = "black"
  ) +
  geom_vline(xintercept = 1.5, linetype = "dashed", color = "black", size = 1) +  # Add vertical line
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 1, vjust = 1, size = 24, color = "black"),
    axis.text.y = element_text(size = 20, color = "black"),
    axis.title.y = element_text(size = 24, color = "black"),
    axis.title.x = element_blank(),
    axis.line.x = element_line(color = "black"),
    axis.line.y = element_line(color = "black"),
    legend.title = element_text(size = 22, color = "black", face = "bold"),  # Increase legend title size
    legend.text = element_text(size = 18, color = "black"),
    legend.position = "top",
    legend.box = "horizontal",
    legend.margin = margin(t = 30, b = 10, l = 20, r = 180),
    plot.margin = margin(t = 40, b = 80, l = 20, r = 20),  # Adjust plot margins
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.key = element_rect(color = "black", size = 1),  # Add border around legend
    legend.background = element_rect(fill = "white", color = "black", size = 1),  # Add background to legend
    strip.text = element_blank(),  # Remove facet labels (SN types)
    panel.spacing = unit(3, "lines")  # Increase space between facets
  ) +
  # Customize colors
  scale_fill_manual(
    values = custom_colors,
    breaks = legend_order,
    labels = c("Higher exposure levels of cancer treatments + Elevated genetic risks", 
               "Higher exposure levels of cancer treatments", 
               "Higher exposure levels of radiotherapy", 
               "Higher exposure levels of chemotherapy", 
               "Elevated genetic risks")
  ) +
  labs(fill = "Risk factors:") +
  scale_x_discrete(expand = c(0, 0), labels = c("SJLIFE", "CCSS")) +  # Set custom x-axis labels
  scale_y_continuous(limits = c(0, 100), breaks = c(0, 50, 100), expand = c(0, 0)) +
  labs(title = "", y = "Attributable fraction (%)", x = NULL) +
  guides(fill = guide_legend(
    title.position = "top",  # Position the legend title at the top
    title.hjust = 0,  # Center align the legend title
    label.hjust = 0,  # Left align the legend labels
    ncol = 1  # Adjust number of columns in legend
  )) +
  facet_wrap(~ SN_types, ncol = 1, scales='free') +  # Add faceting by SN_types
  # Adding the SN_types label at the top right of each facet
  geom_text(data = unique(data_melted %>% select(SN_types)),
            aes(x = Inf, y = Inf, label = SN_types),  # Place at the top right
            size = 8,
            fontface = "bold",
            hjust = 1,
            vjust = 1,
            color = "black",
            inherit.aes = FALSE) +
  coord_cartesian(clip = "off")  # Prevent clipping of text

p
plot_name <- paste0("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/ANALYSIS/results/plots/v21b/combined_5SNs_without_lifestyle_plot.tiff")
ggsave(plot_name, p, width = 12, height = 15, dpi = 400, device = "tiff", compression = "lzw")



