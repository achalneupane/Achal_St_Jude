# Load required libraries
library(ggplot2)
library(ggthemes)# For additional themes
library(ggrepel)
library(dplyr)
library(ggpattern)


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
plot_name <- paste0("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/ANALYSIS/results/plots/v21b/overall/", group, "_", AF.type, "_", lifestyle, ".tiff")

ggsave(plot_name, p, width = 12, height = 8, dpi = 600, device = "tiff", compression = "lzw")
}



################################
## PRS and treatment together ##
################################
group <- "Overall"

variables <- "Treatment|PRS"

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

data_melted$legend_group <- gsub("All_t", "T", factor(data_melted$AF_by))


data_melted$Cohort <- data_melted$variable
data_melted.saved <- data_melted

# Define custom colors and legend order
custom_colors <- c("Treatments" = "#EE6A50", "PRS" = "#EE6A50")
legend_order <- c("Treatments", "PRS")
data_melted <- data_melted.saved[data_melted.saved$Cohort == "SJLIFE",]
data_melted$legend_group <- factor(data_melted$legend_group, levels = c("Treatments", "PRS"))
# Create a factor with the desired order
data_melted$SN_types <- factor(data_melted$SN_types, levels = order)


# Create the grouped bar chart
p_sjlife <- ggplot(data_melted, aes(x = SN_types, y = new_value, fill = legend_group)) +
  geom_bar_pattern(
    stat = "identity",
    position = position_dodge(width = 0.8),
    width = 0.8,
    pattern = ifelse(data_melted$legend_group == "Treatments", "stripe", "none"),
    pattern_density = 0.1,
    pattern_spacing = 0.01,
    pattern_fill = "black",
    pattern_angle = 45,
    colour = "black"
  ) +
  scale_pattern_manual(
    values = c("Treatments" = "stripe", "PRS" = "none"),
    guide = guide_legend(override.aes = list(fill = "gray"))
  ) +
  scale_fill_manual(
    values = custom_colors,
    breaks = legend_order,
    labels = c(
      "Treatments" = "Higher exposure levels\nof cancer treatments",
      "PRS" = "Elevated genetic risk"
    ),
    guide = guide_legend(override.aes = list(pattern = c("Treatments" = "stripe", "PRS" = "none")))
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 16, color = "black"),
        axis.text.y = element_text(size = 16, color = "black"),
        axis.title.y = element_text(size = 20, color = "black"), 
        axis.line.x = element_line(color = "black"), 
        axis.line.y = element_line(color = "black"),
        legend.title = element_blank(),
        legend.text = element_text(size = 16, color = "black"),
        plot.title = element_text(size = 16, hjust = 0.5, vjust = 1.5, face = "bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "top", 
        legend.box = "horizontal") +
  scale_x_discrete(limits = levels(data_melted$SN_types)) +
  coord_cartesian(clip = "off") +
  labs(title = "", y = "Attributable fraction (%)", x = NULL) +
  geom_text(
    aes(label = new_value, y = new_value),  # Adjust y position
    position = position_dodge(width = 0.8),
    vjust = -0.5,  # Adjust vertical justification
    size = 5.5,
    color = "black"
  ) +
  scale_y_continuous(limits = c(0, 100), expand = c(0, 0))
  
# p
p_sjlife <- p_sjlife + scale_y_continuous(limits = c(0, 100), expand = c(0, 0)) + 
  labs(title = "", y = "Attributable fraction (%)", x =  NULL) +
  theme(axis.title.x = element_text(size = 18, color = "black")) +
  theme(axis.title.y = element_text(size = 18, color = "black"))

p_sjlife


# Define custom colors and legend order
custom_colors <- c("Treatments" = "#9FB6CD", "PRS" = "#9FB6CD")
legend_order <- c("Treatments", "PRS")
data_melted <- data_melted.saved[data_melted.saved$Cohort == "CCSS",]
data_melted$legend_group <- factor(data_melted$legend_group, levels = c("Treatments", "PRS"))
# Create a factor with the desired order
data_melted$SN_types <- factor(data_melted$SN_types, levels = order)

# Create the grouped bar chart
p_ccss <- ggplot(data_melted, aes(x = SN_types, y = new_value, fill = legend_group)) +
  geom_bar_pattern(
    stat = "identity",
    position = position_dodge(width = 0.8),
    width = 0.8,
    pattern = ifelse(data_melted$legend_group == "Treatments", "stripe", "none"),
    pattern_density = 0.1,
    pattern_spacing = 0.01,
    pattern_fill = "black",
    pattern_angle = 45,
    colour = "black"
  ) +
  scale_pattern_manual(
    values = c("Treatments" = "stripe", "PRS" = "none"),
    guide = guide_legend(override.aes = list(fill = "gray"))
  ) +
  scale_fill_manual(
    values = custom_colors,
    breaks = legend_order,
    labels = c(
      "Treatments" = "Higher exposure levels\nof cancer treatments",
      "PRS" = "Elevated genetic risk"
    ),
    guide = guide_legend(override.aes = list(pattern = c("Treatments" = "stripe", "PRS" = "none")))
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 16, color = "black"),
        axis.text.y = element_text(size = 16, color = "black"),
        axis.title.y = element_text(size = 20, color = "black"), 
        axis.line.x = element_line(color = "black"), 
        axis.line.y = element_line(color = "black"),
        legend.title = element_blank(),
        legend.text = element_text(size = 16, color = "black"),
        plot.title = element_text(size = 16, hjust = 0.5, vjust = 1.5, face = "bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "top", 
        legend.box = "horizontal") +
  scale_x_discrete(limits = levels(data_melted$SN_types)) +
  coord_cartesian(clip = "off") +
  labs(title = "", y = "Attributable fraction (%)", x = NULL) +
  geom_text(
    aes(label = new_value, y = new_value),  # Adjust y position
    position = position_dodge(width = 0.8),
    vjust = -0.5,  # Adjust vertical justification
    size = 5.5,
    color = "black"
  ) +
  scale_y_continuous(limits = c(0, 100), expand = c(0, 0))
  
# p
p_ccss <- p_ccss + scale_y_continuous(limits = c(0, 100), expand = c(0, 0)) + 
  labs(title = "", y = "Attributable fraction (%)", x =  NULL) +
  theme(axis.title.x = element_text(size = 18, color = "black")) +
  theme(axis.title.y = element_text(size = 18, color = "black"))

p_ccss

p <- cowplot::plot_grid(p_sjlife, p_ccss, labels = c("AUTO"))


p

# Save the plot as a high-resolution image
plot_name <- paste0("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/ANALYSIS/results/plots/v21b/overall/", group, "_", AF.type,"_", lifestyle, ".tiff")

ggsave(plot_name, p, width = 12, height = 8, dpi = 600, device = "tiff", compression = "lzw")



##################################
## Chemo and Radiotherapy together ##
##################################
group <- "Overall"

variables <- "Radiotherapy|Chemotherapy"

new_data <- data[grepl(variables, data$Variables), c("Cohort", "SN_types", "Variables", group)]
colnames(new_data) <- c("variable", "SN_types", "AF_by", "value")
new_data$value <- as.numeric(new_data$value)

# Reshape the data for ggplot2
library(reshape2)
data_melted <-new_data
data_melted$new_value <- round(data_melted$value,2)*100

order <- unique(data_melted$SN_types)
AF.type <- "Chemo_Radiotherapy"
lifestyle <- "without_lifestyle"

data_melted$Cohort <- data_melted$variable
data_melted.saved <- data_melted

# Define custom colors and legend order
custom_colors <- c("Radiotherapy" = "#EE6A50", "Chemotherapy" = "#EE6A50")
legend_order <- c("Radiotherapy", "Chemotherapy")
data_melted <- data_melted.saved[data_melted.saved$Cohort == "SJLIFE",]
data_melted$legend_group <- factor(data_melted$AF_by, levels = c("Radiotherapy", "Chemotherapy"))
# Create a factor with the desired order
data_melted$SN_types <- factor(data_melted$SN_types, levels = order)


# Create the grouped bar chart
p_sjlife <- ggplot(data_melted, aes(x = SN_types, y = new_value, fill = legend_group)) +
  geom_bar_pattern(
    stat = "identity",
    position = position_dodge(width = 0.8),
    width = 0.8,
    pattern = ifelse(data_melted$legend_group == "Radiotherapy", "stripe", "none"),
    pattern_density = 0.1,
    pattern_spacing = 0.01,
    pattern_fill = "black",
    pattern_angle = 45,
    colour = "black"
  ) +
  scale_pattern_manual(
    values = c("Radiotherapy" = "stripe", "Chemotherapy" = "none"),
    guide = guide_legend(override.aes = list(fill = "gray"))
  ) +
  scale_fill_manual(
    values = custom_colors,
    breaks = legend_order,
    labels = c(
      "Radiotherapy" = "Higher exposure levels\nof radiotherapy",
      "Chemotherapy" = "Higher exposure levels\nof chemotherapy"
    ),
    guide = guide_legend(override.aes = list(pattern = c("Radiotherapy" = "stripe", "Chemotherapy" = "none")))
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 16, color = "black"),
        axis.text.y = element_text(size = 16, color = "black"),
        axis.title.y = element_text(size = 20, color = "black"), 
        axis.line.x = element_line(color = "black"), 
        axis.line.y = element_line(color = "black"),
        legend.title = element_blank(),
        legend.text = element_text(size = 16, color = "black"),
        plot.title = element_text(size = 16, hjust = 0.5, vjust = 1.5, face = "bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "top", 
        legend.box = "horizontal") +
  scale_x_discrete(limits = levels(data_melted$SN_types)) +
  coord_cartesian(clip = "off") +
  labs(title = "", y = "Attributable fraction (%)", x = NULL) +
  geom_text(
    aes(label = new_value, y = new_value),  # Adjust y position
    position = position_dodge(width = 0.8),
    vjust = -0.5,  # Adjust vertical justification
    size = 5.5,
    color = "black"
  ) +
  scale_y_continuous(limits = c(0, 100), expand = c(0, 0))

# p
p_sjlife <- p_sjlife + scale_y_continuous(limits = c(0, 100), expand = c(0, 0)) + 
  labs(title = "", y = "Attributable fraction (%)", x =  NULL) +
  theme(axis.title.x = element_text(size = 18, color = "black")) +
  theme(axis.title.y = element_text(size = 18, color = "black"))

# p_sjlife


# Define custom colors and legend order
custom_colors <- c("Radiotherapy" = "#9FB6CD", "Chemotherapy" = "#9FB6CD")
legend_order <- c("Radiotherapy", "Chemotherapy")
data_melted <- data_melted.saved[data_melted.saved$Cohort == "CCSS",]
data_melted$legend_group <- factor(data_melted$AF_by, levels = c("Radiotherapy", "Chemotherapy"))
# Create a factor with the desired order
data_melted$SN_types <- factor(data_melted$SN_types, levels = order)


# Create the grouped bar chart
p_ccss <- ggplot(data_melted, aes(x = SN_types, y = new_value, fill = legend_group)) +
  geom_bar_pattern(
    stat = "identity",
    position = position_dodge(width = 0.8),
    width = 0.8,
    pattern = ifelse(data_melted$legend_group == "Radiotherapy", "stripe", "none"),
    pattern_density = 0.1,
    pattern_spacing = 0.01,
    pattern_fill = "black",
    pattern_angle = 45,
    colour = "black"
  ) +
  scale_pattern_manual(
    values = c("Radiotherapy" = "stripe", "Chemotherapy" = "none"),
    guide = guide_legend(override.aes = list(fill = "gray"))
  ) +
  scale_fill_manual(
    values = custom_colors,
    breaks = legend_order,
    labels = c(
      "Radiotherapy" = "Higher exposure levels\nof radiotherapy",
      "Chemotherapy" = "Higher exposure levels\nof chemotherapy"
    ),
    guide = guide_legend(override.aes = list(pattern = c("Radiotherapy" = "stripe", "Chemotherapy" = "none")))
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 16, color = "black"),
        axis.text.y = element_text(size = 16, color = "black"),
        axis.title.y = element_text(size = 20, color = "black"), 
        axis.line.x = element_line(color = "black"), 
        axis.line.y = element_line(color = "black"),
        legend.title = element_blank(),
        legend.text = element_text(size = 16, color = "black"),
        plot.title = element_text(size = 16, hjust = 0.5, vjust = 1.5, face = "bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "top", 
        legend.box = "horizontal") +
  scale_x_discrete(limits = levels(data_melted$SN_types)) +
  coord_cartesian(clip = "off") +
  labs(title = "", y = "Attributable fraction (%)", x = NULL) +
  geom_text(
    aes(label = new_value, y = new_value),  # Adjust y position
    position = position_dodge(width = 0.8),
    vjust = -0.5,  # Adjust vertical justification
    size = 5.5,
    color = "black"
  ) +
  scale_y_continuous(limits = c(0, 100), expand = c(0, 0))

# p
p_ccss <- p_ccss + scale_y_continuous(limits = c(0, 100), expand = c(0, 0)) + 
  labs(title = "", y = "Attributable fraction (%)", x =  NULL) +
  theme(axis.title.x = element_text(size = 18, color = "black")) +
  theme(axis.title.y = element_text(size = 18, color = "black"))

# p_ccss



p <- cowplot::plot_grid(p_sjlife, p_ccss, labels = c("AUTO"))


# p

# Save the plot as a high-resolution image
plot_name <- paste0("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/ANALYSIS/results/plots/v21b/overall/", group, "_", AF.type,"_", lifestyle, ".tiff")

ggsave(plot_name, p, width = 12, height = 8, dpi = 600, device = "tiff", compression = "lzw")


## 2.___________________________________________________________________________ Female analysis
## 2.___________________________________________________________________________ Female analysis
## 2.___________________________________________________________________________ Female analysis
## 2.___________________________________________________________________________ Female analysis
## 2.___________________________________________________________________________ Female analysis


data <- saved.data

group <- "Female"
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
  plot_name <- paste0("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/ANALYSIS/results/plots/v21b/female/", group, "_", AF.type, "_", lifestyle, ".tiff")
  
  ggsave(plot_name, p, width = 12, height = 8, dpi = 600, device = "tiff", compression = "lzw")
}

################################
## PRS and treatment together ##
################################
group <- "Female"

variables <- "Treatment|PRS"

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

data_melted$legend_group <- gsub("All_t", "T", factor(data_melted$AF_by))


data_melted$Cohort <- data_melted$variable
data_melted.saved <- data_melted

# Define custom colors and legend order
custom_colors <- c("Treatments" = "#EE6A50", "PRS" = "#EE6A50")
legend_order <- c("Treatments", "PRS")
data_melted <- data_melted.saved[data_melted.saved$Cohort == "SJLIFE",]
data_melted$legend_group <- factor(data_melted$legend_group, levels = c("Treatments", "PRS"))
# Create a factor with the desired order
data_melted$SN_types <- factor(data_melted$SN_types, levels = order)


# Create the grouped bar chart
p_sjlife <- ggplot(data_melted, aes(x = SN_types, y = new_value, fill = legend_group)) +
  geom_bar_pattern(
    stat = "identity",
    position = position_dodge(width = 0.8),
    width = 0.8,
    pattern = ifelse(data_melted$legend_group == "Treatments", "stripe", "none"),
    pattern_density = 0.1,
    pattern_spacing = 0.01,
    pattern_fill = "black",
    pattern_angle = 45,
    colour = "black"
  ) +
  scale_pattern_manual(
    values = c("Treatments" = "stripe", "PRS" = "none"),
    guide = guide_legend(override.aes = list(fill = "gray"))
  ) +
  scale_fill_manual(
    values = custom_colors,
    breaks = legend_order,
    labels = c(
      "Treatments" = "Higher exposure levels\nof cancer treatments",
      "PRS" = "Elevated genetic risk"
    ),
    guide = guide_legend(override.aes = list(pattern = c("Treatments" = "stripe", "PRS" = "none")))
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 16, color = "black"),
        axis.text.y = element_text(size = 16, color = "black"),
        axis.title.y = element_text(size = 20, color = "black"), 
        axis.line.x = element_line(color = "black"), 
        axis.line.y = element_line(color = "black"),
        legend.title = element_blank(),
        legend.text = element_text(size = 16, color = "black"),
        plot.title = element_text(size = 16, hjust = 0.5, vjust = 1.5, face = "bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "top", 
        legend.box = "horizontal") +
  scale_x_discrete(limits = levels(data_melted$SN_types)) +
  coord_cartesian(clip = "off") +
  labs(title = "", y = "Attributable fraction (%)", x = NULL) +
  geom_text(
    aes(label = new_value, y = new_value),  # Adjust y position
    position = position_dodge(width = 0.8),
    vjust = -0.5,  # Adjust vertical justification
    size = 5.5,
    color = "black"
  ) +
  scale_y_continuous(limits = c(0, 100), expand = c(0, 0))

# p
p_sjlife <- p_sjlife + scale_y_continuous(limits = c(0, 100), expand = c(0, 0)) + 
  labs(title = "", y = "Attributable fraction (%)", x =  NULL) +
  theme(axis.title.x = element_text(size = 18, color = "black")) +
  theme(axis.title.y = element_text(size = 18, color = "black"))

p_sjlife


# Define custom colors and legend order
custom_colors <- c("Treatments" = "#9FB6CD", "PRS" = "#9FB6CD")
legend_order <- c("Treatments", "PRS")
data_melted <- data_melted.saved[data_melted.saved$Cohort == "CCSS",]
data_melted$legend_group <- factor(data_melted$legend_group, levels = c("Treatments", "PRS"))
# Create a factor with the desired order
data_melted$SN_types <- factor(data_melted$SN_types, levels = order)

# Create the grouped bar chart
p_ccss <- ggplot(data_melted, aes(x = SN_types, y = new_value, fill = legend_group)) +
  geom_bar_pattern(
    stat = "identity",
    position = position_dodge(width = 0.8),
    width = 0.8,
    pattern = ifelse(data_melted$legend_group == "Treatments", "stripe", "none"),
    pattern_density = 0.1,
    pattern_spacing = 0.01,
    pattern_fill = "black",
    pattern_angle = 45,
    colour = "black"
  ) +
  scale_pattern_manual(
    values = c("Treatments" = "stripe", "PRS" = "none"),
    guide = guide_legend(override.aes = list(fill = "gray"))
  ) +
  scale_fill_manual(
    values = custom_colors,
    breaks = legend_order,
    labels = c(
      "Treatments" = "Higher exposure levels\nof cancer treatments",
      "PRS" = "Elevated genetic risk"
    ),
    guide = guide_legend(override.aes = list(pattern = c("Treatments" = "stripe", "PRS" = "none")))
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 16, color = "black"),
        axis.text.y = element_text(size = 16, color = "black"),
        axis.title.y = element_text(size = 20, color = "black"), 
        axis.line.x = element_line(color = "black"), 
        axis.line.y = element_line(color = "black"),
        legend.title = element_blank(),
        legend.text = element_text(size = 16, color = "black"),
        plot.title = element_text(size = 16, hjust = 0.5, vjust = 1.5, face = "bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "top", 
        legend.box = "horizontal") +
  scale_x_discrete(limits = levels(data_melted$SN_types)) +
  coord_cartesian(clip = "off") +
  labs(title = "", y = "Attributable fraction (%)", x = NULL) +
  geom_text(
    aes(label = new_value, y = new_value),  # Adjust y position
    position = position_dodge(width = 0.8),
    vjust = -0.5,  # Adjust vertical justification
    size = 5.5,
    color = "black"
  ) +
  scale_y_continuous(limits = c(0, 100), expand = c(0, 0))

# p
p_ccss <- p_ccss + scale_y_continuous(limits = c(0, 100), expand = c(0, 0)) + 
  labs(title = "", y = "Attributable fraction (%)", x =  NULL) +
  theme(axis.title.x = element_text(size = 18, color = "black")) +
  theme(axis.title.y = element_text(size = 18, color = "black"))

p_ccss

p <- cowplot::plot_grid(p_sjlife, p_ccss, labels = c("AUTO"))


p

# Save the plot as a high-resolution image
plot_name <- paste0("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/ANALYSIS/results/plots/v21b/female/", group, "_", AF.type,"_", lifestyle, ".tiff")

ggsave(plot_name, p, width = 12, height = 8, dpi = 600, device = "tiff", compression = "lzw")



##################################
## Chemo and Radiotherapy together ##
##################################
group <- "Female"

variables <- "Radiotherapy|Chemotherapy"

new_data <- data[grepl(variables, data$Variables), c("Cohort", "SN_types", "Variables", group)]
colnames(new_data) <- c("variable", "SN_types", "AF_by", "value")
new_data$value <- as.numeric(new_data$value)

# Reshape the data for ggplot2
library(reshape2)
data_melted <-new_data
data_melted$new_value <- round(data_melted$value,2)*100

order <- unique(data_melted$SN_types)
AF.type <- "Chemo_Radiotherapy"
lifestyle <- "without_lifestyle"

data_melted$Cohort <- data_melted$variable
data_melted.saved <- data_melted

# Define custom colors and legend order
custom_colors <- c("Radiotherapy" = "#EE6A50", "Chemotherapy" = "#EE6A50")
legend_order <- c("Radiotherapy", "Chemotherapy")
data_melted <- data_melted.saved[data_melted.saved$Cohort == "SJLIFE",]
data_melted$legend_group <- factor(data_melted$AF_by, levels = c("Radiotherapy", "Chemotherapy"))
# Create a factor with the desired order
data_melted$SN_types <- factor(data_melted$SN_types, levels = order)


# Create the grouped bar chart
p_sjlife <- ggplot(data_melted, aes(x = SN_types, y = new_value, fill = legend_group)) +
  geom_bar_pattern(
    stat = "identity",
    position = position_dodge(width = 0.8),
    width = 0.8,
    pattern = ifelse(data_melted$legend_group == "Radiotherapy", "stripe", "none"),
    pattern_density = 0.1,
    pattern_spacing = 0.01,
    pattern_fill = "black",
    pattern_angle = 45,
    colour = "black"
  ) +
  scale_pattern_manual(
    values = c("Radiotherapy" = "stripe", "Chemotherapy" = "none"),
    guide = guide_legend(override.aes = list(fill = "gray"))
  ) +
  scale_fill_manual(
    values = custom_colors,
    breaks = legend_order,
    labels = c(
      "Radiotherapy" = "Higher exposure levels\nof radiotherapy",
      "Chemotherapy" = "Higher exposure levels\nof chemotherapy"
    ),
    guide = guide_legend(override.aes = list(pattern = c("Radiotherapy" = "stripe", "Chemotherapy" = "none")))
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 16, color = "black"),
        axis.text.y = element_text(size = 16, color = "black"),
        axis.title.y = element_text(size = 20, color = "black"), 
        axis.line.x = element_line(color = "black"), 
        axis.line.y = element_line(color = "black"),
        legend.title = element_blank(),
        legend.text = element_text(size = 16, color = "black"),
        plot.title = element_text(size = 16, hjust = 0.5, vjust = 1.5, face = "bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "top", 
        legend.box = "horizontal") +
  scale_x_discrete(limits = levels(data_melted$SN_types)) +
  coord_cartesian(clip = "off") +
  labs(title = "", y = "Attributable fraction (%)", x = NULL) +
  geom_text(
    aes(label = new_value, y = new_value),  # Adjust y position
    position = position_dodge(width = 0.8),
    vjust = -0.5,  # Adjust vertical justification
    size = 5.5,
    color = "black"
  ) +
  scale_y_continuous(limits = c(0, 100), expand = c(0, 0))

# p
p_sjlife <- p_sjlife + scale_y_continuous(limits = c(0, 100), expand = c(0, 0)) + 
  labs(title = "", y = "Attributable fraction (%)", x =  NULL) +
  theme(axis.title.x = element_text(size = 18, color = "black")) +
  theme(axis.title.y = element_text(size = 18, color = "black"))

# p_sjlife


# Define custom colors and legend order
custom_colors <- c("Radiotherapy" = "#9FB6CD", "Chemotherapy" = "#9FB6CD")
legend_order <- c("Radiotherapy", "Chemotherapy")
data_melted <- data_melted.saved[data_melted.saved$Cohort == "CCSS",]
data_melted$legend_group <- factor(data_melted$AF_by, levels = c("Radiotherapy", "Chemotherapy"))
# Create a factor with the desired order
data_melted$SN_types <- factor(data_melted$SN_types, levels = order)


# Create the grouped bar chart
p_ccss <- ggplot(data_melted, aes(x = SN_types, y = new_value, fill = legend_group)) +
  geom_bar_pattern(
    stat = "identity",
    position = position_dodge(width = 0.8),
    width = 0.8,
    pattern = ifelse(data_melted$legend_group == "Radiotherapy", "stripe", "none"),
    pattern_density = 0.1,
    pattern_spacing = 0.01,
    pattern_fill = "black",
    pattern_angle = 45,
    colour = "black"
  ) +
  scale_pattern_manual(
    values = c("Radiotherapy" = "stripe", "Chemotherapy" = "none"),
    guide = guide_legend(override.aes = list(fill = "gray"))
  ) +
  scale_fill_manual(
    values = custom_colors,
    breaks = legend_order,
    labels = c(
      "Radiotherapy" = "Higher exposure levels\nof radiotherapy",
      "Chemotherapy" = "Higher exposure levels\nof chemotherapy"
    ),
    guide = guide_legend(override.aes = list(pattern = c("Radiotherapy" = "stripe", "Chemotherapy" = "none")))
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 16, color = "black"),
        axis.text.y = element_text(size = 16, color = "black"),
        axis.title.y = element_text(size = 20, color = "black"), 
        axis.line.x = element_line(color = "black"), 
        axis.line.y = element_line(color = "black"),
        legend.title = element_blank(),
        legend.text = element_text(size = 16, color = "black"),
        plot.title = element_text(size = 16, hjust = 0.5, vjust = 1.5, face = "bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "top", 
        legend.box = "horizontal") +
  scale_x_discrete(limits = levels(data_melted$SN_types)) +
  coord_cartesian(clip = "off") +
  labs(title = "", y = "Attributable fraction (%)", x = NULL) +
  geom_text(
    aes(label = new_value, y = new_value),  # Adjust y position
    position = position_dodge(width = 0.8),
    vjust = -0.5,  # Adjust vertical justification
    size = 5.5,
    color = "black"
  ) +
  scale_y_continuous(limits = c(0, 100), expand = c(0, 0))

# p
p_ccss <- p_ccss + scale_y_continuous(limits = c(0, 100), expand = c(0, 0)) + 
  labs(title = "", y = "Attributable fraction (%)", x =  NULL) +
  theme(axis.title.x = element_text(size = 18, color = "black")) +
  theme(axis.title.y = element_text(size = 18, color = "black"))

# p_ccss



p <- cowplot::plot_grid(p_sjlife, p_ccss, labels = c("AUTO"))


# p

# Save the plot as a high-resolution image
plot_name <- paste0("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/ANALYSIS/results/plots/v21b/female/", group, "_", AF.type,"_", lifestyle, ".tiff")

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

custom_colors <- c("SJLIFE" = "#EE6A50", "CCSS" = "#9FB6CD")
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
  plot_name <- paste0("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/ANALYSIS/results/plots/v21b/male/", group, "_", AF.type, "_", lifestyle, ".tiff")
  
  ggsave(plot_name, p, width = 12, height = 8, dpi = 600, device = "tiff", compression = "lzw")
}

################################
## PRS and treatment together ##
################################
group <- "Male"

variables <- "Treatment|PRS"

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

data_melted$legend_group <- gsub("All_t", "T", factor(data_melted$AF_by))


data_melted$Cohort <- data_melted$variable
data_melted.saved <- data_melted

# Define custom colors and legend order
custom_colors <- c("Treatments" = "#EE6A50", "PRS" = "#EE6A50")
legend_order <- c("Treatments", "PRS")
data_melted <- data_melted.saved[data_melted.saved$Cohort == "SJLIFE",]
data_melted$legend_group <- factor(data_melted$legend_group, levels = c("Treatments", "PRS"))
# Create a factor with the desired order
data_melted$SN_types <- factor(data_melted$SN_types, levels = order)


# Create the grouped bar chart
p_sjlife <- ggplot(data_melted, aes(x = SN_types, y = new_value, fill = legend_group)) +
  geom_bar_pattern(
    stat = "identity",
    position = position_dodge(width = 0.8),
    width = 0.8,
    pattern = ifelse(data_melted$legend_group == "Treatments", "stripe", "none"),
    pattern_density = 0.1,
    pattern_spacing = 0.01,
    pattern_fill = "black",
    pattern_angle = 45,
    colour = "black"
  ) +
  scale_pattern_manual(
    values = c("Treatments" = "stripe", "PRS" = "none"),
    guide = guide_legend(override.aes = list(fill = "gray"))
  ) +
  scale_fill_manual(
    values = custom_colors,
    breaks = legend_order,
    labels = c(
      "Treatments" = "Higher exposure levels\nof cancer treatments",
      "PRS" = "Elevated genetic risk"
    ),
    guide = guide_legend(override.aes = list(pattern = c("Treatments" = "stripe", "PRS" = "none")))
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 16, color = "black"),
        axis.text.y = element_text(size = 16, color = "black"),
        axis.title.y = element_text(size = 20, color = "black"), 
        axis.line.x = element_line(color = "black"), 
        axis.line.y = element_line(color = "black"),
        legend.title = element_blank(),
        legend.text = element_text(size = 16, color = "black"),
        plot.title = element_text(size = 16, hjust = 0.5, vjust = 1.5, face = "bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "top", 
        legend.box = "horizontal") +
  scale_x_discrete(limits = levels(data_melted$SN_types)) +
  coord_cartesian(clip = "off") +
  labs(title = "", y = "Attributable fraction (%)", x = NULL) +
  geom_text(
    aes(label = new_value, y = new_value),  # Adjust y position
    position = position_dodge(width = 0.8),
    vjust = -0.5,  # Adjust vertical justification
    size = 5.5,
    color = "black"
  ) +
  scale_y_continuous(limits = c(0, 100), expand = c(0, 0))

# p
p_sjlife <- p_sjlife + scale_y_continuous(limits = c(0, 100), expand = c(0, 0)) + 
  labs(title = "", y = "Attributable fraction (%)", x =  NULL) +
  theme(axis.title.x = element_text(size = 18, color = "black")) +
  theme(axis.title.y = element_text(size = 18, color = "black"))

p_sjlife


# Define custom colors and legend order
custom_colors <- c("Treatments" = "#9FB6CD", "PRS" = "#9FB6CD")
legend_order <- c("Treatments", "PRS")
data_melted <- data_melted.saved[data_melted.saved$Cohort == "CCSS",]
data_melted$legend_group <- factor(data_melted$legend_group, levels = c("Treatments", "PRS"))
# Create a factor with the desired order
data_melted$SN_types <- factor(data_melted$SN_types, levels = order)

# Create the grouped bar chart
p_ccss <- ggplot(data_melted, aes(x = SN_types, y = new_value, fill = legend_group)) +
  geom_bar_pattern(
    stat = "identity",
    position = position_dodge(width = 0.8),
    width = 0.8,
    pattern = ifelse(data_melted$legend_group == "Treatments", "stripe", "none"),
    pattern_density = 0.1,
    pattern_spacing = 0.01,
    pattern_fill = "black",
    pattern_angle = 45,
    colour = "black"
  ) +
  scale_pattern_manual(
    values = c("Treatments" = "stripe", "PRS" = "none"),
    guide = guide_legend(override.aes = list(fill = "gray"))
  ) +
  scale_fill_manual(
    values = custom_colors,
    breaks = legend_order,
    labels = c(
      "Treatments" = "Higher exposure levels\nof cancer treatments",
      "PRS" = "Elevated genetic risk"
    ),
    guide = guide_legend(override.aes = list(pattern = c("Treatments" = "stripe", "PRS" = "none")))
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 16, color = "black"),
        axis.text.y = element_text(size = 16, color = "black"),
        axis.title.y = element_text(size = 20, color = "black"), 
        axis.line.x = element_line(color = "black"), 
        axis.line.y = element_line(color = "black"),
        legend.title = element_blank(),
        legend.text = element_text(size = 16, color = "black"),
        plot.title = element_text(size = 16, hjust = 0.5, vjust = 1.5, face = "bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "top", 
        legend.box = "horizontal") +
  scale_x_discrete(limits = levels(data_melted$SN_types)) +
  coord_cartesian(clip = "off") +
  labs(title = "", y = "Attributable fraction (%)", x = NULL) +
  geom_text(
    aes(label = new_value, y = new_value),  # Adjust y position
    position = position_dodge(width = 0.8),
    vjust = -0.5,  # Adjust vertical justification
    size = 5.5,
    color = "black"
  ) +
  scale_y_continuous(limits = c(0, 100), expand = c(0, 0))

# p
p_ccss <- p_ccss + scale_y_continuous(limits = c(0, 100), expand = c(0, 0)) + 
  labs(title = "", y = "Attributable fraction (%)", x =  NULL) +
  theme(axis.title.x = element_text(size = 18, color = "black")) +
  theme(axis.title.y = element_text(size = 18, color = "black"))

p_ccss

p <- cowplot::plot_grid(p_sjlife, p_ccss, labels = c("AUTO"))


p

# Save the plot as a high-resolution image
plot_name <- paste0("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/ANALYSIS/results/plots/v21b/male/", group, "_", AF.type,"_", lifestyle, ".tiff")

ggsave(plot_name, p, width = 12, height = 8, dpi = 600, device = "tiff", compression = "lzw")

##################################
## Chemo and Radiotherapy together ##
##################################
group <- "Male"

variables <- "Radiotherapy|Chemotherapy"

new_data <- data[grepl(variables, data$Variables), c("Cohort", "SN_types", "Variables", group)]
colnames(new_data) <- c("variable", "SN_types", "AF_by", "value")
new_data$value <- as.numeric(new_data$value)

# Reshape the data for ggplot2
library(reshape2)
data_melted <-new_data
data_melted$new_value <- round(data_melted$value,2)*100

order <- unique(data_melted$SN_types)
AF.type <- "Chemo_Radiotherapy"
lifestyle <- "without_lifestyle"

data_melted$Cohort <- data_melted$variable
data_melted.saved <- data_melted

# Define custom colors and legend order
custom_colors <- c("Radiotherapy" = "#EE6A50", "Chemotherapy" = "#EE6A50")
legend_order <- c("Radiotherapy", "Chemotherapy")
data_melted <- data_melted.saved[data_melted.saved$Cohort == "SJLIFE",]
data_melted$legend_group <- factor(data_melted$AF_by, levels = c("Radiotherapy", "Chemotherapy"))
# Create a factor with the desired order
data_melted$SN_types <- factor(data_melted$SN_types, levels = order)


# Create the grouped bar chart
p_sjlife <- ggplot(data_melted, aes(x = SN_types, y = new_value, fill = legend_group)) +
  geom_bar_pattern(
    stat = "identity",
    position = position_dodge(width = 0.8),
    width = 0.8,
    pattern = ifelse(data_melted$legend_group == "Radiotherapy", "stripe", "none"),
    pattern_density = 0.1,
    pattern_spacing = 0.01,
    pattern_fill = "black",
    pattern_angle = 45,
    colour = "black"
  ) +
  scale_pattern_manual(
    values = c("Radiotherapy" = "stripe", "Chemotherapy" = "none"),
    guide = guide_legend(override.aes = list(fill = "gray"))
  ) +
  scale_fill_manual(
    values = custom_colors,
    breaks = legend_order,
    labels = c(
      "Radiotherapy" = "Higher exposure levels\nof radiotherapy",
      "Chemotherapy" = "Higher exposure levels\nof chemotherapy"
    ),
    guide = guide_legend(override.aes = list(pattern = c("Radiotherapy" = "stripe", "Chemotherapy" = "none")))
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 16, color = "black"),
        axis.text.y = element_text(size = 16, color = "black"),
        axis.title.y = element_text(size = 20, color = "black"), 
        axis.line.x = element_line(color = "black"), 
        axis.line.y = element_line(color = "black"),
        legend.title = element_blank(),
        legend.text = element_text(size = 16, color = "black"),
        plot.title = element_text(size = 16, hjust = 0.5, vjust = 1.5, face = "bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "top", 
        legend.box = "horizontal") +
  scale_x_discrete(limits = levels(data_melted$SN_types)) +
  coord_cartesian(clip = "off") +
  labs(title = "", y = "Attributable fraction (%)", x = NULL) +
  geom_text(
    aes(label = new_value, y = new_value),  # Adjust y position
    position = position_dodge(width = 0.8),
    vjust = -0.5,  # Adjust vertical justification
    size = 5.5,
    color = "black"
  ) +
  scale_y_continuous(limits = c(0, 100), expand = c(0, 0))

# p
p_sjlife <- p_sjlife + scale_y_continuous(limits = c(0, 100), expand = c(0, 0)) + 
  labs(title = "", y = "Attributable fraction (%)", x =  NULL) +
  theme(axis.title.x = element_text(size = 18, color = "black")) +
  theme(axis.title.y = element_text(size = 18, color = "black"))

# p_sjlife


# Define custom colors and legend order
custom_colors <- c("Radiotherapy" = "#9FB6CD", "Chemotherapy" = "#9FB6CD")
legend_order <- c("Radiotherapy", "Chemotherapy")
data_melted <- data_melted.saved[data_melted.saved$Cohort == "CCSS",]
data_melted$legend_group <- factor(data_melted$AF_by, levels = c("Radiotherapy", "Chemotherapy"))
# Create a factor with the desired order
data_melted$SN_types <- factor(data_melted$SN_types, levels = order)


# Create the grouped bar chart
p_ccss <- ggplot(data_melted, aes(x = SN_types, y = new_value, fill = legend_group)) +
  geom_bar_pattern(
    stat = "identity",
    position = position_dodge(width = 0.8),
    width = 0.8,
    pattern = ifelse(data_melted$legend_group == "Radiotherapy", "stripe", "none"),
    pattern_density = 0.1,
    pattern_spacing = 0.01,
    pattern_fill = "black",
    pattern_angle = 45,
    colour = "black"
  ) +
  scale_pattern_manual(
    values = c("Radiotherapy" = "stripe", "Chemotherapy" = "none"),
    guide = guide_legend(override.aes = list(fill = "gray"))
  ) +
  scale_fill_manual(
    values = custom_colors,
    breaks = legend_order,
    labels = c(
      "Radiotherapy" = "Higher exposure levels\nof radiotherapy",
      "Chemotherapy" = "Higher exposure levels\nof chemotherapy"
    ),
    guide = guide_legend(override.aes = list(pattern = c("Radiotherapy" = "stripe", "Chemotherapy" = "none")))
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 16, color = "black"),
        axis.text.y = element_text(size = 16, color = "black"),
        axis.title.y = element_text(size = 20, color = "black"), 
        axis.line.x = element_line(color = "black"), 
        axis.line.y = element_line(color = "black"),
        legend.title = element_blank(),
        legend.text = element_text(size = 16, color = "black"),
        plot.title = element_text(size = 16, hjust = 0.5, vjust = 1.5, face = "bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "top", 
        legend.box = "horizontal") +
  scale_x_discrete(limits = levels(data_melted$SN_types)) +
  coord_cartesian(clip = "off") +
  labs(title = "", y = "Attributable fraction (%)", x = NULL) +
  geom_text(
    aes(label = new_value, y = new_value),  # Adjust y position
    position = position_dodge(width = 0.8),
    vjust = -0.5,  # Adjust vertical justification
    size = 5.5,
    color = "black"
  ) +
  scale_y_continuous(limits = c(0, 100), expand = c(0, 0))

# p
p_ccss <- p_ccss + scale_y_continuous(limits = c(0, 100), expand = c(0, 0)) + 
  labs(title = "", y = "Attributable fraction (%)", x =  NULL) +
  theme(axis.title.x = element_text(size = 18, color = "black")) +
  theme(axis.title.y = element_text(size = 18, color = "black"))

# p_ccss

p <- cowplot::plot_grid(p_sjlife, p_ccss, labels = c("AUTO"))


# p

# Save the plot as a high-resolution image
plot_name <- paste0("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/ANALYSIS/results/plots/v21b/male/", group, "_", AF.type,"_", lifestyle, ".tiff")

ggsave(plot_name, p, width = 12, height = 8, dpi = 600, device = "tiff", compression = "lzw")

## 4.___________________________________________________________________________ lt35 analysis
## 4.___________________________________________________________________________ lt35 analysis
## 4.___________________________________________________________________________ lt35 analysis
## 4.___________________________________________________________________________ lt35 analysis
## 4.___________________________________________________________________________ lt35 analysis


data <- saved.data

group <- "Age.lt.35"
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
  plot_name <- paste0("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/ANALYSIS/results/plots/v21b/lt35/", group, "_", AF.type, "_", lifestyle, ".tiff")
  
  ggsave(plot_name, p, width = 12, height = 8, dpi = 600, device = "tiff", compression = "lzw")
}

################################
## PRS and treatment together ##
################################
group <- "Age.lt.35"

variables <- "Treatment|PRS"

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

data_melted$legend_group <- gsub("All_t", "T", factor(data_melted$AF_by))


data_melted$Cohort <- data_melted$variable
data_melted.saved <- data_melted

# Define custom colors and legend order
custom_colors <- c("Treatments" = "#EE6A50", "PRS" = "#EE6A50")
legend_order <- c("Treatments", "PRS")
data_melted <- data_melted.saved[data_melted.saved$Cohort == "SJLIFE",]
data_melted$legend_group <- factor(data_melted$legend_group, levels = c("Treatments", "PRS"))
# Create a factor with the desired order
data_melted$SN_types <- factor(data_melted$SN_types, levels = order)


# Create the grouped bar chart
p_sjlife <- ggplot(data_melted, aes(x = SN_types, y = new_value, fill = legend_group)) +
  geom_bar_pattern(
    stat = "identity",
    position = position_dodge(width = 0.8),
    width = 0.8,
    pattern = ifelse(data_melted$legend_group == "Treatments", "stripe", "none"),
    pattern_density = 0.1,
    pattern_spacing = 0.01,
    pattern_fill = "black",
    pattern_angle = 45,
    colour = "black"
  ) +
  scale_pattern_manual(
    values = c("Treatments" = "stripe", "PRS" = "none"),
    guide = guide_legend(override.aes = list(fill = "gray"))
  ) +
  scale_fill_manual(
    values = custom_colors,
    breaks = legend_order,
    labels = c(
      "Treatments" = "Higher exposure levels\nof cancer treatments",
      "PRS" = "Elevated genetic risk"
    ),
    guide = guide_legend(override.aes = list(pattern = c("Treatments" = "stripe", "PRS" = "none")))
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 16, color = "black"),
        axis.text.y = element_text(size = 16, color = "black"),
        axis.title.y = element_text(size = 20, color = "black"), 
        axis.line.x = element_line(color = "black"), 
        axis.line.y = element_line(color = "black"),
        legend.title = element_blank(),
        legend.text = element_text(size = 16, color = "black"),
        plot.title = element_text(size = 16, hjust = 0.5, vjust = 1.5, face = "bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "top", 
        legend.box = "horizontal") +
  scale_x_discrete(limits = levels(data_melted$SN_types)) +
  coord_cartesian(clip = "off") +
  labs(title = "", y = "Attributable fraction (%)", x = NULL) +
  geom_text(
    aes(label = new_value, y = new_value),  # Adjust y position
    position = position_dodge(width = 0.8),
    vjust = -0.5,  # Adjust vertical justification
    size = 5.5,
    color = "black"
  ) +
  scale_y_continuous(limits = c(0, 100), expand = c(0, 0))

# p
p_sjlife <- p_sjlife + scale_y_continuous(limits = c(0, 100), expand = c(0, 0)) + 
  labs(title = "", y = "Attributable fraction (%)", x =  NULL) +
  theme(axis.title.x = element_text(size = 18, color = "black")) +
  theme(axis.title.y = element_text(size = 18, color = "black"))

p_sjlife


# Define custom colors and legend order
custom_colors <- c("Treatments" = "#9FB6CD", "PRS" = "#9FB6CD")
legend_order <- c("Treatments", "PRS")
data_melted <- data_melted.saved[data_melted.saved$Cohort == "CCSS",]
data_melted$legend_group <- factor(data_melted$legend_group, levels = c("Treatments", "PRS"))
# Create a factor with the desired order
data_melted$SN_types <- factor(data_melted$SN_types, levels = order)

# Create the grouped bar chart
p_ccss <- ggplot(data_melted, aes(x = SN_types, y = new_value, fill = legend_group)) +
  geom_bar_pattern(
    stat = "identity",
    position = position_dodge(width = 0.8),
    width = 0.8,
    pattern = ifelse(data_melted$legend_group == "Treatments", "stripe", "none"),
    pattern_density = 0.1,
    pattern_spacing = 0.01,
    pattern_fill = "black",
    pattern_angle = 45,
    colour = "black"
  ) +
  scale_pattern_manual(
    values = c("Treatments" = "stripe", "PRS" = "none"),
    guide = guide_legend(override.aes = list(fill = "gray"))
  ) +
  scale_fill_manual(
    values = custom_colors,
    breaks = legend_order,
    labels = c(
      "Treatments" = "Higher exposure levels\nof cancer treatments",
      "PRS" = "Elevated genetic risk"
    ),
    guide = guide_legend(override.aes = list(pattern = c("Treatments" = "stripe", "PRS" = "none")))
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 16, color = "black"),
        axis.text.y = element_text(size = 16, color = "black"),
        axis.title.y = element_text(size = 20, color = "black"), 
        axis.line.x = element_line(color = "black"), 
        axis.line.y = element_line(color = "black"),
        legend.title = element_blank(),
        legend.text = element_text(size = 16, color = "black"),
        plot.title = element_text(size = 16, hjust = 0.5, vjust = 1.5, face = "bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "top", 
        legend.box = "horizontal") +
  scale_x_discrete(limits = levels(data_melted$SN_types)) +
  coord_cartesian(clip = "off") +
  labs(title = "", y = "Attributable fraction (%)", x = NULL) +
  geom_text(
    aes(label = new_value, y = new_value),  # Adjust y position
    position = position_dodge(width = 0.8),
    vjust = -0.5,  # Adjust vertical justification
    size = 5.5,
    color = "black"
  ) +
  scale_y_continuous(limits = c(0, 100), expand = c(0, 0))

# p
p_ccss <- p_ccss + scale_y_continuous(limits = c(0, 100), expand = c(0, 0)) + 
  labs(title = "", y = "Attributable fraction (%)", x =  NULL) +
  theme(axis.title.x = element_text(size = 18, color = "black")) +
  theme(axis.title.y = element_text(size = 18, color = "black"))

p_ccss

p <- cowplot::plot_grid(p_sjlife, p_ccss, labels = c("AUTO"))


p
# Save the plot as a high-resolution image
plot_name <- paste0("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/ANALYSIS/results/plots/v21b/lt35/", group, "_", AF.type,"_", lifestyle, ".tiff")

ggsave(plot_name, p, width = 12, height = 8, dpi = 600, device = "tiff", compression = "lzw")

##################################
## Chemo and Radiotherapy together ##
##################################
group <- "Age.lt.35"

variables <- "Radiotherapy|Chemotherapy"

new_data <- data[grepl(variables, data$Variables), c("Cohort", "SN_types", "Variables", group)]
colnames(new_data) <- c("variable", "SN_types", "AF_by", "value")
new_data$value <- as.numeric(new_data$value)

# Reshape the data for ggplot2
library(reshape2)
data_melted <-new_data
data_melted$new_value <- round(data_melted$value,2)*100

order <- unique(data_melted$SN_types)
AF.type <- "Chemo_Radiotherapy"
lifestyle <- "without_lifestyle"

data_melted$Cohort <- data_melted$variable
data_melted.saved <- data_melted

# Define custom colors and legend order
custom_colors <- c("Radiotherapy" = "#EE6A50", "Chemotherapy" = "#EE6A50")
legend_order <- c("Radiotherapy", "Chemotherapy")
data_melted <- data_melted.saved[data_melted.saved$Cohort == "SJLIFE",]
data_melted$legend_group <- factor(data_melted$AF_by, levels = c("Radiotherapy", "Chemotherapy"))
# Create a factor with the desired order
data_melted$SN_types <- factor(data_melted$SN_types, levels = order)


# Create the grouped bar chart
p_sjlife <- ggplot(data_melted, aes(x = SN_types, y = new_value, fill = legend_group)) +
  geom_bar_pattern(
    stat = "identity",
    position = position_dodge(width = 0.8),
    width = 0.8,
    pattern = ifelse(data_melted$legend_group == "Radiotherapy", "stripe", "none"),
    pattern_density = 0.1,
    pattern_spacing = 0.01,
    pattern_fill = "black",
    pattern_angle = 45,
    colour = "black"
  ) +
  scale_pattern_manual(
    values = c("Radiotherapy" = "stripe", "Chemotherapy" = "none"),
    guide = guide_legend(override.aes = list(fill = "gray"))
  ) +
  scale_fill_manual(
    values = custom_colors,
    breaks = legend_order,
    labels = c(
      "Radiotherapy" = "Higher exposure levels\nof radiotherapy",
      "Chemotherapy" = "Higher exposure levels\nof chemotherapy"
    ),
    guide = guide_legend(override.aes = list(pattern = c("Radiotherapy" = "stripe", "Chemotherapy" = "none")))
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 16, color = "black"),
        axis.text.y = element_text(size = 16, color = "black"),
        axis.title.y = element_text(size = 20, color = "black"), 
        axis.line.x = element_line(color = "black"), 
        axis.line.y = element_line(color = "black"),
        legend.title = element_blank(),
        legend.text = element_text(size = 16, color = "black"),
        plot.title = element_text(size = 16, hjust = 0.5, vjust = 1.5, face = "bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "top", 
        legend.box = "horizontal") +
  scale_x_discrete(limits = levels(data_melted$SN_types)) +
  coord_cartesian(clip = "off") +
  labs(title = "", y = "Attributable fraction (%)", x = NULL) +
  geom_text(
    aes(label = new_value, y = new_value),  # Adjust y position
    position = position_dodge(width = 0.8),
    vjust = -0.5,  # Adjust vertical justification
    size = 5.5,
    color = "black"
  ) +
  scale_y_continuous(limits = c(0, 100), expand = c(0, 0))

# p
p_sjlife <- p_sjlife + scale_y_continuous(limits = c(0, 100), expand = c(0, 0)) + 
  labs(title = "", y = "Attributable fraction (%)", x =  NULL) +
  theme(axis.title.x = element_text(size = 18, color = "black")) +
  theme(axis.title.y = element_text(size = 18, color = "black"))

# p_sjlife


# Define custom colors and legend order
custom_colors <- c("Radiotherapy" = "#9FB6CD", "Chemotherapy" = "#9FB6CD")
legend_order <- c("Radiotherapy", "Chemotherapy")
data_melted <- data_melted.saved[data_melted.saved$Cohort == "CCSS",]
data_melted$legend_group <- factor(data_melted$AF_by, levels = c("Radiotherapy", "Chemotherapy"))
# Create a factor with the desired order
data_melted$SN_types <- factor(data_melted$SN_types, levels = order)


# Create the grouped bar chart
p_ccss <- ggplot(data_melted, aes(x = SN_types, y = new_value, fill = legend_group)) +
  geom_bar_pattern(
    stat = "identity",
    position = position_dodge(width = 0.8),
    width = 0.8,
    pattern = ifelse(data_melted$legend_group == "Radiotherapy", "stripe", "none"),
    pattern_density = 0.1,
    pattern_spacing = 0.01,
    pattern_fill = "black",
    pattern_angle = 45,
    colour = "black"
  ) +
  scale_pattern_manual(
    values = c("Radiotherapy" = "stripe", "Chemotherapy" = "none"),
    guide = guide_legend(override.aes = list(fill = "gray"))
  ) +
  scale_fill_manual(
    values = custom_colors,
    breaks = legend_order,
    labels = c(
      "Radiotherapy" = "Higher exposure levels\nof radiotherapy",
      "Chemotherapy" = "Higher exposure levels\nof chemotherapy"
    ),
    guide = guide_legend(override.aes = list(pattern = c("Radiotherapy" = "stripe", "Chemotherapy" = "none")))
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 16, color = "black"),
        axis.text.y = element_text(size = 16, color = "black"),
        axis.title.y = element_text(size = 20, color = "black"), 
        axis.line.x = element_line(color = "black"), 
        axis.line.y = element_line(color = "black"),
        legend.title = element_blank(),
        legend.text = element_text(size = 16, color = "black"),
        plot.title = element_text(size = 16, hjust = 0.5, vjust = 1.5, face = "bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "top", 
        legend.box = "horizontal") +
  scale_x_discrete(limits = levels(data_melted$SN_types)) +
  coord_cartesian(clip = "off") +
  labs(title = "", y = "Attributable fraction (%)", x = NULL) +
  geom_text(
    aes(label = new_value, y = new_value),  # Adjust y position
    position = position_dodge(width = 0.8),
    vjust = -0.5,  # Adjust vertical justification
    size = 5.5,
    color = "black"
  ) +
  scale_y_continuous(limits = c(0, 100), expand = c(0, 0))

# p
p_ccss <- p_ccss + scale_y_continuous(limits = c(0, 100), expand = c(0, 0)) + 
  labs(title = "", y = "Attributable fraction (%)", x =  NULL) +
  theme(axis.title.x = element_text(size = 18, color = "black")) +
  theme(axis.title.y = element_text(size = 18, color = "black"))

# p_ccss



p <- cowplot::plot_grid(p_sjlife, p_ccss, labels = c("AUTO"))


# p

# Save the plot as a high-resolution image
plot_name <- paste0("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/ANALYSIS/results/plots/v21b/lt35/", group, "_", AF.type,"_", lifestyle, ".tiff")

ggsave(plot_name, p, width = 12, height = 8, dpi = 600, device = "tiff", compression = "lzw")

## 5.___________________________________________________________________________ ge35 analysis
## 5.___________________________________________________________________________ ge35 analysis
## 5.___________________________________________________________________________ ge35 analysis
## 5.___________________________________________________________________________ ge35 analysis
## 5.___________________________________________________________________________ ge35 analysis

data <- saved.data

group <- "Age.ge.35"
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
  plot_name <- paste0("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/ANALYSIS/results/plots/v21b/ge35/", group, "_", AF.type, "_", lifestyle, ".tiff")
  
  ggsave(plot_name, p, width = 12, height = 8, dpi = 600, device = "tiff", compression = "lzw")
}

################################
## PRS and treatment together ##
################################
group <- "Age.ge.35"

variables <- "Treatment|PRS"

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

data_melted$legend_group <- gsub("All_t", "T", factor(data_melted$AF_by))


data_melted$Cohort <- data_melted$variable
data_melted.saved <- data_melted

# Define custom colors and legend order
custom_colors <- c("Treatments" = "#EE6A50", "PRS" = "#EE6A50")
legend_order <- c("Treatments", "PRS")
data_melted <- data_melted.saved[data_melted.saved$Cohort == "SJLIFE",]
data_melted$legend_group <- factor(data_melted$legend_group, levels = c("Treatments", "PRS"))
# Create a factor with the desired order
data_melted$SN_types <- factor(data_melted$SN_types, levels = order)


# Create the grouped bar chart
p_sjlife <- ggplot(data_melted, aes(x = SN_types, y = new_value, fill = legend_group)) +
  geom_bar_pattern(
    stat = "identity",
    position = position_dodge(width = 0.8),
    width = 0.8,
    pattern = ifelse(data_melted$legend_group == "Treatments", "stripe", "none"),
    pattern_density = 0.1,
    pattern_spacing = 0.01,
    pattern_fill = "black",
    pattern_angle = 45,
    colour = "black"
  ) +
  scale_pattern_manual(
    values = c("Treatments" = "stripe", "PRS" = "none"),
    guide = guide_legend(override.aes = list(fill = "gray"))
  ) +
  scale_fill_manual(
    values = custom_colors,
    breaks = legend_order,
    labels = c(
      "Treatments" = "Higher exposure levels\nof cancer treatments",
      "PRS" = "Elevated genetic risk"
    ),
    guide = guide_legend(override.aes = list(pattern = c("Treatments" = "stripe", "PRS" = "none")))
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 16, color = "black"),
        axis.text.y = element_text(size = 16, color = "black"),
        axis.title.y = element_text(size = 20, color = "black"), 
        axis.line.x = element_line(color = "black"), 
        axis.line.y = element_line(color = "black"),
        legend.title = element_blank(),
        legend.text = element_text(size = 16, color = "black"),
        plot.title = element_text(size = 16, hjust = 0.5, vjust = 1.5, face = "bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "top", 
        legend.box = "horizontal") +
  scale_x_discrete(limits = levels(data_melted$SN_types)) +
  coord_cartesian(clip = "off") +
  labs(title = "", y = "Attributable fraction (%)", x = NULL) +
  geom_text(
    aes(label = new_value, y = new_value),  # Adjust y position
    position = position_dodge(width = 0.8),
    vjust = -0.5,  # Adjust vertical justification
    size = 5.5,
    color = "black"
  ) +
  scale_y_continuous(limits = c(0, 100), expand = c(0, 0))

# p
p_sjlife <- p_sjlife + scale_y_continuous(limits = c(0, 100), expand = c(0, 0)) + 
  labs(title = "", y = "Attributable fraction (%)", x =  NULL) +
  theme(axis.title.x = element_text(size = 18, color = "black")) +
  theme(axis.title.y = element_text(size = 18, color = "black"))

p_sjlife


# Define custom colors and legend order
custom_colors <- c("Treatments" = "#9FB6CD", "PRS" = "#9FB6CD")
legend_order <- c("Treatments", "PRS")
data_melted <- data_melted.saved[data_melted.saved$Cohort == "CCSS",]
data_melted$legend_group <- factor(data_melted$legend_group, levels = c("Treatments", "PRS"))
# Create a factor with the desired order
data_melted$SN_types <- factor(data_melted$SN_types, levels = order)

# Create the grouped bar chart
p_ccss <- ggplot(data_melted, aes(x = SN_types, y = new_value, fill = legend_group)) +
  geom_bar_pattern(
    stat = "identity",
    position = position_dodge(width = 0.8),
    width = 0.8,
    pattern = ifelse(data_melted$legend_group == "Treatments", "stripe", "none"),
    pattern_density = 0.1,
    pattern_spacing = 0.01,
    pattern_fill = "black",
    pattern_angle = 45,
    colour = "black"
  ) +
  scale_pattern_manual(
    values = c("Treatments" = "stripe", "PRS" = "none"),
    guide = guide_legend(override.aes = list(fill = "gray"))
  ) +
  scale_fill_manual(
    values = custom_colors,
    breaks = legend_order,
    labels = c(
      "Treatments" = "Higher exposure levels\nof cancer treatments",
      "PRS" = "Elevated genetic risk"
    ),
    guide = guide_legend(override.aes = list(pattern = c("Treatments" = "stripe", "PRS" = "none")))
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 16, color = "black"),
        axis.text.y = element_text(size = 16, color = "black"),
        axis.title.y = element_text(size = 20, color = "black"), 
        axis.line.x = element_line(color = "black"), 
        axis.line.y = element_line(color = "black"),
        legend.title = element_blank(),
        legend.text = element_text(size = 16, color = "black"),
        plot.title = element_text(size = 16, hjust = 0.5, vjust = 1.5, face = "bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "top", 
        legend.box = "horizontal") +
  scale_x_discrete(limits = levels(data_melted$SN_types)) +
  coord_cartesian(clip = "off") +
  labs(title = "", y = "Attributable fraction (%)", x = NULL) +
  geom_text(
    aes(label = new_value, y = new_value),  # Adjust y position
    position = position_dodge(width = 0.8),
    vjust = -0.5,  # Adjust vertical justification
    size = 5.5,
    color = "black"
  ) +
  scale_y_continuous(limits = c(0, 100), expand = c(0, 0))

# p
p_ccss <- p_ccss + scale_y_continuous(limits = c(0, 100), expand = c(0, 0)) + 
  labs(title = "", y = "Attributable fraction (%)", x =  NULL) +
  theme(axis.title.x = element_text(size = 18, color = "black")) +
  theme(axis.title.y = element_text(size = 18, color = "black"))

p_ccss

p <- cowplot::plot_grid(p_sjlife, p_ccss, labels = c("AUTO"))


p

# Save the plot as a high-resolution image
plot_name <- paste0("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/ANALYSIS/results/plots/v21b/ge35/", group, "_", AF.type,"_", lifestyle, ".tiff")

ggsave(plot_name, p, width = 12, height = 8, dpi = 600, device = "tiff", compression = "lzw")

##################################
## Chemo and Radiotherapy together ##
##################################
group <- "Age.ge.35"

variables <- "Radiotherapy|Chemotherapy"

new_data <- data[grepl(variables, data$Variables), c("Cohort", "SN_types", "Variables", group)]
colnames(new_data) <- c("variable", "SN_types", "AF_by", "value")
new_data$value <- as.numeric(new_data$value)

# Reshape the data for ggplot2
library(reshape2)
data_melted <-new_data
data_melted$new_value <- round(data_melted$value,2)*100

order <- unique(data_melted$SN_types)
AF.type <- "Chemo_Radiotherapy"
lifestyle <- "without_lifestyle"

data_melted$Cohort <- data_melted$variable
data_melted.saved <- data_melted

# Define custom colors and legend order
custom_colors <- c("Radiotherapy" = "#EE6A50", "Chemotherapy" = "#EE6A50")
legend_order <- c("Radiotherapy", "Chemotherapy")
data_melted <- data_melted.saved[data_melted.saved$Cohort == "SJLIFE",]
data_melted$legend_group <- factor(data_melted$AF_by, levels = c("Radiotherapy", "Chemotherapy"))
# Create a factor with the desired order
data_melted$SN_types <- factor(data_melted$SN_types, levels = order)


# Create the grouped bar chart
p_sjlife <- ggplot(data_melted, aes(x = SN_types, y = new_value, fill = legend_group)) +
  geom_bar_pattern(
    stat = "identity",
    position = position_dodge(width = 0.8),
    width = 0.8,
    pattern = ifelse(data_melted$legend_group == "Radiotherapy", "stripe", "none"),
    pattern_density = 0.1,
    pattern_spacing = 0.01,
    pattern_fill = "black",
    pattern_angle = 45,
    colour = "black"
  ) +
  scale_pattern_manual(
    values = c("Radiotherapy" = "stripe", "Chemotherapy" = "none"),
    guide = guide_legend(override.aes = list(fill = "gray"))
  ) +
  scale_fill_manual(
    values = custom_colors,
    breaks = legend_order,
    labels = c(
      "Radiotherapy" = "Higher exposure levels\nof radiotherapy",
      "Chemotherapy" = "Higher exposure levels\nof chemotherapy"
    ),
    guide = guide_legend(override.aes = list(pattern = c("Radiotherapy" = "stripe", "Chemotherapy" = "none")))
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 16, color = "black"),
        axis.text.y = element_text(size = 16, color = "black"),
        axis.title.y = element_text(size = 20, color = "black"), 
        axis.line.x = element_line(color = "black"), 
        axis.line.y = element_line(color = "black"),
        legend.title = element_blank(),
        legend.text = element_text(size = 16, color = "black"),
        plot.title = element_text(size = 16, hjust = 0.5, vjust = 1.5, face = "bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "top", 
        legend.box = "horizontal") +
  scale_x_discrete(limits = levels(data_melted$SN_types)) +
  coord_cartesian(clip = "off") +
  labs(title = "", y = "Attributable fraction (%)", x = NULL) +
  geom_text(
    aes(label = new_value, y = new_value),  # Adjust y position
    position = position_dodge(width = 0.8),
    vjust = -0.5,  # Adjust vertical justification
    size = 5.5,
    color = "black"
  ) +
  scale_y_continuous(limits = c(0, 100), expand = c(0, 0))

# p
p_sjlife <- p_sjlife + scale_y_continuous(limits = c(0, 100), expand = c(0, 0)) + 
  labs(title = "", y = "Attributable fraction (%)", x =  NULL) +
  theme(axis.title.x = element_text(size = 18, color = "black")) +
  theme(axis.title.y = element_text(size = 18, color = "black"))

# p_sjlife


# Define custom colors and legend order
custom_colors <- c("Radiotherapy" = "#9FB6CD", "Chemotherapy" = "#9FB6CD")
legend_order <- c("Radiotherapy", "Chemotherapy")
data_melted <- data_melted.saved[data_melted.saved$Cohort == "CCSS",]
data_melted$legend_group <- factor(data_melted$AF_by, levels = c("Radiotherapy", "Chemotherapy"))
# Create a factor with the desired order
data_melted$SN_types <- factor(data_melted$SN_types, levels = order)


# Create the grouped bar chart
p_ccss <- ggplot(data_melted, aes(x = SN_types, y = new_value, fill = legend_group)) +
  geom_bar_pattern(
    stat = "identity",
    position = position_dodge(width = 0.8),
    width = 0.8,
    pattern = ifelse(data_melted$legend_group == "Radiotherapy", "stripe", "none"),
    pattern_density = 0.1,
    pattern_spacing = 0.01,
    pattern_fill = "black",
    pattern_angle = 45,
    colour = "black"
  ) +
  scale_pattern_manual(
    values = c("Radiotherapy" = "stripe", "Chemotherapy" = "none"),
    guide = guide_legend(override.aes = list(fill = "gray"))
  ) +
  scale_fill_manual(
    values = custom_colors,
    breaks = legend_order,
    labels = c(
      "Radiotherapy" = "Higher exposure levels\nof radiotherapy",
      "Chemotherapy" = "Higher exposure levels\nof chemotherapy"
    ),
    guide = guide_legend(override.aes = list(pattern = c("Radiotherapy" = "stripe", "Chemotherapy" = "none")))
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 16, color = "black"),
        axis.text.y = element_text(size = 16, color = "black"),
        axis.title.y = element_text(size = 20, color = "black"), 
        axis.line.x = element_line(color = "black"), 
        axis.line.y = element_line(color = "black"),
        legend.title = element_blank(),
        legend.text = element_text(size = 16, color = "black"),
        plot.title = element_text(size = 16, hjust = 0.5, vjust = 1.5, face = "bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "top", 
        legend.box = "horizontal") +
  scale_x_discrete(limits = levels(data_melted$SN_types)) +
  coord_cartesian(clip = "off") +
  labs(title = "", y = "Attributable fraction (%)", x = NULL) +
  geom_text(
    aes(label = new_value, y = new_value),  # Adjust y position
    position = position_dodge(width = 0.8),
    vjust = -0.5,  # Adjust vertical justification
    size = 5.5,
    color = "black"
  ) +
  scale_y_continuous(limits = c(0, 100), expand = c(0, 0))

# p
p_ccss <- p_ccss + scale_y_continuous(limits = c(0, 100), expand = c(0, 0)) + 
  labs(title = "", y = "Attributable fraction (%)", x =  NULL) +
  theme(axis.title.x = element_text(size = 18, color = "black")) +
  theme(axis.title.y = element_text(size = 18, color = "black"))

# p_ccss



p <- cowplot::plot_grid(p_sjlife, p_ccss, labels = c("AUTO"))


# p
# Save the plot as a high-resolution image
plot_name <- paste0("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/ANALYSIS/results/plots/v21b/ge35/", group, "_", AF.type,"_", lifestyle, ".tiff")

ggsave(plot_name, p, width = 12, height = 8, dpi = 600, device = "tiff", compression = "lzw")

## 5.___________________________________________________________________________ EUR analysis
## 5.___________________________________________________________________________ EUR analysis
## 5.___________________________________________________________________________ EUR analysis
## 5.___________________________________________________________________________ EUR analysis
## 5.___________________________________________________________________________ EUR analysis

data <- saved.data

group <- "EUR"
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
  plot_name <- paste0("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/ANALYSIS/results/plots/v21b/EUR/", group, "_", AF.type, "_", lifestyle, ".tiff")
  
  ggsave(plot_name, p, width = 12, height = 8, dpi = 600, device = "tiff", compression = "lzw")
}

################################
## PRS and treatment together ##
################################
group <- "EUR"

variables <- "Treatment|PRS"

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

data_melted$legend_group <- gsub("All_t", "T", factor(data_melted$AF_by))


data_melted$Cohort <- data_melted$variable
data_melted.saved <- data_melted

# Define custom colors and legend order
custom_colors <- c("Treatments" = "#EE6A50", "PRS" = "#EE6A50")
legend_order <- c("Treatments", "PRS")
data_melted <- data_melted.saved[data_melted.saved$Cohort == "SJLIFE",]
data_melted$legend_group <- factor(data_melted$legend_group, levels = c("Treatments", "PRS"))
# Create a factor with the desired order
data_melted$SN_types <- factor(data_melted$SN_types, levels = order)


# Create the grouped bar chart
p_sjlife <- ggplot(data_melted, aes(x = SN_types, y = new_value, fill = legend_group)) +
  geom_bar_pattern(
    stat = "identity",
    position = position_dodge(width = 0.8),
    width = 0.8,
    pattern = ifelse(data_melted$legend_group == "Treatments", "stripe", "none"),
    pattern_density = 0.1,
    pattern_spacing = 0.01,
    pattern_fill = "black",
    pattern_angle = 45,
    colour = "black"
  ) +
  scale_pattern_manual(
    values = c("Treatments" = "stripe", "PRS" = "none"),
    guide = guide_legend(override.aes = list(fill = "gray"))
  ) +
  scale_fill_manual(
    values = custom_colors,
    breaks = legend_order,
    labels = c(
      "Treatments" = "Higher exposure levels\nof cancer treatments",
      "PRS" = "Elevated genetic risk"
    ),
    guide = guide_legend(override.aes = list(pattern = c("Treatments" = "stripe", "PRS" = "none")))
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 16, color = "black"),
        axis.text.y = element_text(size = 16, color = "black"),
        axis.title.y = element_text(size = 20, color = "black"), 
        axis.line.x = element_line(color = "black"), 
        axis.line.y = element_line(color = "black"),
        legend.title = element_blank(),
        legend.text = element_text(size = 16, color = "black"),
        plot.title = element_text(size = 16, hjust = 0.5, vjust = 1.5, face = "bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "top", 
        legend.box = "horizontal") +
  scale_x_discrete(limits = levels(data_melted$SN_types)) +
  coord_cartesian(clip = "off") +
  labs(title = "", y = "Attributable fraction (%)", x = NULL) +
  geom_text(
    aes(label = new_value, y = new_value),  # Adjust y position
    position = position_dodge(width = 0.8),
    vjust = -0.5,  # Adjust vertical justification
    size = 5.5,
    color = "black"
  ) +
  scale_y_continuous(limits = c(0, 100), expand = c(0, 0))

# p
p_sjlife <- p_sjlife + scale_y_continuous(limits = c(0, 100), expand = c(0, 0)) + 
  labs(title = "", y = "Attributable fraction (%)", x =  NULL) +
  theme(axis.title.x = element_text(size = 18, color = "black")) +
  theme(axis.title.y = element_text(size = 18, color = "black"))

p_sjlife


# Define custom colors and legend order
custom_colors <- c("Treatments" = "#9FB6CD", "PRS" = "#9FB6CD")
legend_order <- c("Treatments", "PRS")
data_melted <- data_melted.saved[data_melted.saved$Cohort == "CCSS",]
data_melted$legend_group <- factor(data_melted$legend_group, levels = c("Treatments", "PRS"))
# Create a factor with the desired order
data_melted$SN_types <- factor(data_melted$SN_types, levels = order)

# Create the grouped bar chart
p_ccss <- ggplot(data_melted, aes(x = SN_types, y = new_value, fill = legend_group)) +
  geom_bar_pattern(
    stat = "identity",
    position = position_dodge(width = 0.8),
    width = 0.8,
    pattern = ifelse(data_melted$legend_group == "Treatments", "stripe", "none"),
    pattern_density = 0.1,
    pattern_spacing = 0.01,
    pattern_fill = "black",
    pattern_angle = 45,
    colour = "black"
  ) +
  scale_pattern_manual(
    values = c("Treatments" = "stripe", "PRS" = "none"),
    guide = guide_legend(override.aes = list(fill = "gray"))
  ) +
  scale_fill_manual(
    values = custom_colors,
    breaks = legend_order,
    labels = c(
      "Treatments" = "Higher exposure levels\nof cancer treatments",
      "PRS" = "Elevated genetic risk"
    ),
    guide = guide_legend(override.aes = list(pattern = c("Treatments" = "stripe", "PRS" = "none")))
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 16, color = "black"),
        axis.text.y = element_text(size = 16, color = "black"),
        axis.title.y = element_text(size = 20, color = "black"), 
        axis.line.x = element_line(color = "black"), 
        axis.line.y = element_line(color = "black"),
        legend.title = element_blank(),
        legend.text = element_text(size = 16, color = "black"),
        plot.title = element_text(size = 16, hjust = 0.5, vjust = 1.5, face = "bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "top", 
        legend.box = "horizontal") +
  scale_x_discrete(limits = levels(data_melted$SN_types)) +
  coord_cartesian(clip = "off") +
  labs(title = "", y = "Attributable fraction (%)", x = NULL) +
  geom_text(
    aes(label = new_value, y = new_value),  # Adjust y position
    position = position_dodge(width = 0.8),
    vjust = -0.5,  # Adjust vertical justification
    size = 5.5,
    color = "black"
  ) +
  scale_y_continuous(limits = c(0, 100), expand = c(0, 0))

# p
p_ccss <- p_ccss + scale_y_continuous(limits = c(0, 100), expand = c(0, 0)) + 
  labs(title = "", y = "Attributable fraction (%)", x =  NULL) +
  theme(axis.title.x = element_text(size = 18, color = "black")) +
  theme(axis.title.y = element_text(size = 18, color = "black"))

p_ccss

p <- cowplot::plot_grid(p_sjlife, p_ccss, labels = c("AUTO"))


p

# Save the plot as a high-resolution image
plot_name <- paste0("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/ANALYSIS/results/plots/v21b/EUR/", group, "_", AF.type,"_", lifestyle, ".tiff")

ggsave(plot_name, p, width = 12, height = 8, dpi = 600, device = "tiff", compression = "lzw")

##################################
## Chemo and Radiotherapy together ##
##################################
group <- "EUR"

variables <- "Radiotherapy|Chemotherapy"

new_data <- data[grepl(variables, data$Variables), c("Cohort", "SN_types", "Variables", group)]
colnames(new_data) <- c("variable", "SN_types", "AF_by", "value")
new_data$value <- as.numeric(new_data$value)

# Reshape the data for ggplot2
library(reshape2)
data_melted <-new_data
data_melted$new_value <- round(data_melted$value,2)*100

order <- unique(data_melted$SN_types)
AF.type <- "Chemo_Radiotherapy"
lifestyle <- "without_lifestyle"

data_melted$Cohort <- data_melted$variable
data_melted.saved <- data_melted

# Define custom colors and legend order
custom_colors <- c("Radiotherapy" = "#EE6A50", "Chemotherapy" = "#EE6A50")
legend_order <- c("Radiotherapy", "Chemotherapy")
data_melted <- data_melted.saved[data_melted.saved$Cohort == "SJLIFE",]
data_melted$legend_group <- factor(data_melted$AF_by, levels = c("Radiotherapy", "Chemotherapy"))
# Create a factor with the desired order
data_melted$SN_types <- factor(data_melted$SN_types, levels = order)


# Create the grouped bar chart
p_sjlife <- ggplot(data_melted, aes(x = SN_types, y = new_value, fill = legend_group)) +
  geom_bar_pattern(
    stat = "identity",
    position = position_dodge(width = 0.8),
    width = 0.8,
    pattern = ifelse(data_melted$legend_group == "Radiotherapy", "stripe", "none"),
    pattern_density = 0.1,
    pattern_spacing = 0.01,
    pattern_fill = "black",
    pattern_angle = 45,
    colour = "black"
  ) +
  scale_pattern_manual(
    values = c("Radiotherapy" = "stripe", "Chemotherapy" = "none"),
    guide = guide_legend(override.aes = list(fill = "gray"))
  ) +
  scale_fill_manual(
    values = custom_colors,
    breaks = legend_order,
    labels = c(
      "Radiotherapy" = "Higher exposure levels\nof radiotherapy",
      "Chemotherapy" = "Higher exposure levels\nof chemotherapy"
    ),
    guide = guide_legend(override.aes = list(pattern = c("Radiotherapy" = "stripe", "Chemotherapy" = "none")))
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 16, color = "black"),
        axis.text.y = element_text(size = 16, color = "black"),
        axis.title.y = element_text(size = 20, color = "black"), 
        axis.line.x = element_line(color = "black"), 
        axis.line.y = element_line(color = "black"),
        legend.title = element_blank(),
        legend.text = element_text(size = 16, color = "black"),
        plot.title = element_text(size = 16, hjust = 0.5, vjust = 1.5, face = "bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "top", 
        legend.box = "horizontal") +
  scale_x_discrete(limits = levels(data_melted$SN_types)) +
  coord_cartesian(clip = "off") +
  labs(title = "", y = "Attributable fraction (%)", x = NULL) +
  geom_text(
    aes(label = new_value, y = new_value),  # Adjust y position
    position = position_dodge(width = 0.8),
    vjust = -0.5,  # Adjust vertical justification
    size = 5.5,
    color = "black"
  ) +
  scale_y_continuous(limits = c(0, 100), expand = c(0, 0))

# p
p_sjlife <- p_sjlife + scale_y_continuous(limits = c(0, 100), expand = c(0, 0)) + 
  labs(title = "", y = "Attributable fraction (%)", x =  NULL) +
  theme(axis.title.x = element_text(size = 18, color = "black")) +
  theme(axis.title.y = element_text(size = 18, color = "black"))

# p_sjlife


# Define custom colors and legend order
custom_colors <- c("Radiotherapy" = "#9FB6CD", "Chemotherapy" = "#9FB6CD")
legend_order <- c("Radiotherapy", "Chemotherapy")
data_melted <- data_melted.saved[data_melted.saved$Cohort == "CCSS",]
data_melted$legend_group <- factor(data_melted$AF_by, levels = c("Radiotherapy", "Chemotherapy"))
# Create a factor with the desired order
data_melted$SN_types <- factor(data_melted$SN_types, levels = order)


# Create the grouped bar chart
p_ccss <- ggplot(data_melted, aes(x = SN_types, y = new_value, fill = legend_group)) +
  geom_bar_pattern(
    stat = "identity",
    position = position_dodge(width = 0.8),
    width = 0.8,
    pattern = ifelse(data_melted$legend_group == "Radiotherapy", "stripe", "none"),
    pattern_density = 0.1,
    pattern_spacing = 0.01,
    pattern_fill = "black",
    pattern_angle = 45,
    colour = "black"
  ) +
  scale_pattern_manual(
    values = c("Radiotherapy" = "stripe", "Chemotherapy" = "none"),
    guide = guide_legend(override.aes = list(fill = "gray"))
  ) +
  scale_fill_manual(
    values = custom_colors,
    breaks = legend_order,
    labels = c(
      "Radiotherapy" = "Higher exposure levels\nof radiotherapy",
      "Chemotherapy" = "Higher exposure levels\nof chemotherapy"
    ),
    guide = guide_legend(override.aes = list(pattern = c("Radiotherapy" = "stripe", "Chemotherapy" = "none")))
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 16, color = "black"),
        axis.text.y = element_text(size = 16, color = "black"),
        axis.title.y = element_text(size = 20, color = "black"), 
        axis.line.x = element_line(color = "black"), 
        axis.line.y = element_line(color = "black"),
        legend.title = element_blank(),
        legend.text = element_text(size = 16, color = "black"),
        plot.title = element_text(size = 16, hjust = 0.5, vjust = 1.5, face = "bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "top", 
        legend.box = "horizontal") +
  scale_x_discrete(limits = levels(data_melted$SN_types)) +
  coord_cartesian(clip = "off") +
  labs(title = "", y = "Attributable fraction (%)", x = NULL) +
  geom_text(
    aes(label = new_value, y = new_value),  # Adjust y position
    position = position_dodge(width = 0.8),
    vjust = -0.5,  # Adjust vertical justification
    size = 5.5,
    color = "black"
  ) +
  scale_y_continuous(limits = c(0, 100), expand = c(0, 0))

# p
p_ccss <- p_ccss + scale_y_continuous(limits = c(0, 100), expand = c(0, 0)) + 
  labs(title = "", y = "Attributable fraction (%)", x =  NULL) +
  theme(axis.title.x = element_text(size = 18, color = "black")) +
  theme(axis.title.y = element_text(size = 18, color = "black"))

# p_ccss



p <- cowplot::plot_grid(p_sjlife, p_ccss, labels = c("AUTO"))


# p
# Save the plot as a high-resolution image
plot_name <- paste0("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/ANALYSIS/results/plots/v21b/EUR/", group, "_", AF.type,"_", lifestyle, ".tiff")

ggsave(plot_name, p, width = 12, height = 8, dpi = 600, device = "tiff", compression = "lzw")

## 5.___________________________________________________________________________ AFR analysis
## 5.___________________________________________________________________________ AFR analysis
## 5.___________________________________________________________________________ AFR analysis
## 5.___________________________________________________________________________ AFR analysis
## 5.___________________________________________________________________________ AFR analysis

data <- saved.data

group <- "AFR"
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
  plot_name <- paste0("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/ANALYSIS/results/plots/v21b/AFR/", group, "_", AF.type, "_", lifestyle, ".tiff")
  
  ggsave(plot_name, p, width = 12, height = 8, dpi = 600, device = "tiff", compression = "lzw")
}

################################
## PRS and treatment together ##
################################
group <- "AFR"

variables <- "Treatment|PRS"

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

data_melted$legend_group <- gsub("All_t", "T", factor(data_melted$AF_by))


data_melted$Cohort <- data_melted$variable
data_melted.saved <- data_melted

# Define custom colors and legend order
custom_colors <- c("Treatments" = "#EE6A50", "PRS" = "#EE6A50")
legend_order <- c("Treatments", "PRS")
data_melted <- data_melted.saved[data_melted.saved$Cohort == "SJLIFE",]
data_melted$legend_group <- factor(data_melted$legend_group, levels = c("Treatments", "PRS"))
# Create a factor with the desired order
data_melted$SN_types <- factor(data_melted$SN_types, levels = order)


# Create the grouped bar chart
p_sjlife <- ggplot(data_melted, aes(x = SN_types, y = new_value, fill = legend_group)) +
  geom_bar_pattern(
    stat = "identity",
    position = position_dodge(width = 0.8),
    width = 0.8,
    pattern = ifelse(data_melted$legend_group == "Treatments", "stripe", "none"),
    pattern_density = 0.1,
    pattern_spacing = 0.01,
    pattern_fill = "black",
    pattern_angle = 45,
    colour = "black"
  ) +
  scale_pattern_manual(
    values = c("Treatments" = "stripe", "PRS" = "none"),
    guide = guide_legend(override.aes = list(fill = "gray"))
  ) +
  scale_fill_manual(
    values = custom_colors,
    breaks = legend_order,
    labels = c(
      "Treatments" = "Higher exposure levels\nof cancer treatments",
      "PRS" = "Elevated genetic risk"
    ),
    guide = guide_legend(override.aes = list(pattern = c("Treatments" = "stripe", "PRS" = "none")))
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 16, color = "black"),
        axis.text.y = element_text(size = 16, color = "black"),
        axis.title.y = element_text(size = 20, color = "black"), 
        axis.line.x = element_line(color = "black"), 
        axis.line.y = element_line(color = "black"),
        legend.title = element_blank(),
        legend.text = element_text(size = 16, color = "black"),
        plot.title = element_text(size = 16, hjust = 0.5, vjust = 1.5, face = "bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "top", 
        legend.box = "horizontal") +
  scale_x_discrete(limits = levels(data_melted$SN_types)) +
  coord_cartesian(clip = "off") +
  labs(title = "", y = "Attributable fraction (%)", x = NULL) +
  geom_text(
    aes(label = new_value, y = new_value),  # Adjust y position
    position = position_dodge(width = 0.8),
    vjust = -0.5,  # Adjust vertical justification
    size = 5.5,
    color = "black"
  ) +
  scale_y_continuous(limits = c(0, 100), expand = c(0, 0))

# p
p_sjlife <- p_sjlife + scale_y_continuous(limits = c(0, 100), expand = c(0, 0)) + 
  labs(title = "", y = "Attributable fraction (%)", x =  NULL) +
  theme(axis.title.x = element_text(size = 18, color = "black")) +
  theme(axis.title.y = element_text(size = 18, color = "black"))

p_sjlife


# Define custom colors and legend order
custom_colors <- c("Treatments" = "#9FB6CD", "PRS" = "#9FB6CD")
legend_order <- c("Treatments", "PRS")
data_melted <- data_melted.saved[data_melted.saved$Cohort == "CCSS",]
data_melted$legend_group <- factor(data_melted$legend_group, levels = c("Treatments", "PRS"))
# Create a factor with the desired order
data_melted$SN_types <- factor(data_melted$SN_types, levels = order)

# Create the grouped bar chart
p_ccss <- ggplot(data_melted, aes(x = SN_types, y = new_value, fill = legend_group)) +
  geom_bar_pattern(
    stat = "identity",
    position = position_dodge(width = 0.8),
    width = 0.8,
    pattern = ifelse(data_melted$legend_group == "Treatments", "stripe", "none"),
    pattern_density = 0.1,
    pattern_spacing = 0.01,
    pattern_fill = "black",
    pattern_angle = 45,
    colour = "black"
  ) +
  scale_pattern_manual(
    values = c("Treatments" = "stripe", "PRS" = "none"),
    guide = guide_legend(override.aes = list(fill = "gray"))
  ) +
  scale_fill_manual(
    values = custom_colors,
    breaks = legend_order,
    labels = c(
      "Treatments" = "Higher exposure levels\nof cancer treatments",
      "PRS" = "Elevated genetic risk"
    ),
    guide = guide_legend(override.aes = list(pattern = c("Treatments" = "stripe", "PRS" = "none")))
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 16, color = "black"),
        axis.text.y = element_text(size = 16, color = "black"),
        axis.title.y = element_text(size = 20, color = "black"), 
        axis.line.x = element_line(color = "black"), 
        axis.line.y = element_line(color = "black"),
        legend.title = element_blank(),
        legend.text = element_text(size = 16, color = "black"),
        plot.title = element_text(size = 16, hjust = 0.5, vjust = 1.5, face = "bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "top", 
        legend.box = "horizontal") +
  scale_x_discrete(limits = levels(data_melted$SN_types)) +
  coord_cartesian(clip = "off") +
  labs(title = "", y = "Attributable fraction (%)", x = NULL) +
  geom_text(
    aes(label = new_value, y = new_value),  # Adjust y position
    position = position_dodge(width = 0.8),
    vjust = -0.5,  # Adjust vertical justification
    size = 5.5,
    color = "black"
  ) +
  scale_y_continuous(limits = c(0, 100), expand = c(0, 0))

# p
p_ccss <- p_ccss + scale_y_continuous(limits = c(0, 100), expand = c(0, 0)) + 
  labs(title = "", y = "Attributable fraction (%)", x =  NULL) +
  theme(axis.title.x = element_text(size = 18, color = "black")) +
  theme(axis.title.y = element_text(size = 18, color = "black"))

p_ccss

p <- cowplot::plot_grid(p_sjlife, p_ccss, labels = c("AUTO"))


p

# Save the plot as a high-resolution image
plot_name <- paste0("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/ANALYSIS/results/plots/v21b/AFR/", group, "_", AF.type,"_", lifestyle, ".tiff")

ggsave(plot_name, p, width = 12, height = 8, dpi = 600, device = "tiff", compression = "lzw")

##################################
## Chemo and Radiotherapy together ##
##################################
group <- "AFR"

variables <- "Radiotherapy|Chemotherapy"

new_data <- data[grepl(variables, data$Variables), c("Cohort", "SN_types", "Variables", group)]
colnames(new_data) <- c("variable", "SN_types", "AF_by", "value")
new_data$value <- as.numeric(new_data$value)

# Reshape the data for ggplot2
library(reshape2)
data_melted <-new_data
data_melted$new_value <- round(data_melted$value,2)*100

order <- unique(data_melted$SN_types)
AF.type <- "Chemo_Radiotherapy"
lifestyle <- "without_lifestyle"

data_melted$Cohort <- data_melted$variable
data_melted.saved <- data_melted

# Define custom colors and legend order
custom_colors <- c("Radiotherapy" = "#EE6A50", "Chemotherapy" = "#EE6A50")
legend_order <- c("Radiotherapy", "Chemotherapy")
data_melted <- data_melted.saved[data_melted.saved$Cohort == "SJLIFE",]
data_melted$legend_group <- factor(data_melted$AF_by, levels = c("Radiotherapy", "Chemotherapy"))
# Create a factor with the desired order
data_melted$SN_types <- factor(data_melted$SN_types, levels = order)


# Create the grouped bar chart
p_sjlife <- ggplot(data_melted, aes(x = SN_types, y = new_value, fill = legend_group)) +
  geom_bar_pattern(
    stat = "identity",
    position = position_dodge(width = 0.8),
    width = 0.8,
    pattern = ifelse(data_melted$legend_group == "Radiotherapy", "stripe", "none"),
    pattern_density = 0.1,
    pattern_spacing = 0.01,
    pattern_fill = "black",
    pattern_angle = 45,
    colour = "black"
  ) +
  scale_pattern_manual(
    values = c("Radiotherapy" = "stripe", "Chemotherapy" = "none"),
    guide = guide_legend(override.aes = list(fill = "gray"))
  ) +
  scale_fill_manual(
    values = custom_colors,
    breaks = legend_order,
    labels = c(
      "Radiotherapy" = "Higher exposure levels\nof radiotherapy",
      "Chemotherapy" = "Higher exposure levels\nof chemotherapy"
    ),
    guide = guide_legend(override.aes = list(pattern = c("Radiotherapy" = "stripe", "Chemotherapy" = "none")))
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 16, color = "black"),
        axis.text.y = element_text(size = 16, color = "black"),
        axis.title.y = element_text(size = 20, color = "black"), 
        axis.line.x = element_line(color = "black"), 
        axis.line.y = element_line(color = "black"),
        legend.title = element_blank(),
        legend.text = element_text(size = 16, color = "black"),
        plot.title = element_text(size = 16, hjust = 0.5, vjust = 1.5, face = "bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "top", 
        legend.box = "horizontal") +
  scale_x_discrete(limits = levels(data_melted$SN_types)) +
  coord_cartesian(clip = "off") +
  labs(title = "", y = "Attributable fraction (%)", x = NULL) +
  geom_text(
    aes(label = new_value, y = new_value),  # Adjust y position
    position = position_dodge(width = 0.8),
    vjust = -0.5,  # Adjust vertical justification
    size = 5.5,
    color = "black"
  ) +
  scale_y_continuous(limits = c(0, 100), expand = c(0, 0))

# p
p_sjlife <- p_sjlife + scale_y_continuous(limits = c(0, 100), expand = c(0, 0)) + 
  labs(title = "", y = "Attributable fraction (%)", x =  NULL) +
  theme(axis.title.x = element_text(size = 18, color = "black")) +
  theme(axis.title.y = element_text(size = 18, color = "black"))

# p_sjlife


# Define custom colors and legend order
custom_colors <- c("Radiotherapy" = "#9FB6CD", "Chemotherapy" = "#9FB6CD")
legend_order <- c("Radiotherapy", "Chemotherapy")
data_melted <- data_melted.saved[data_melted.saved$Cohort == "CCSS",]
data_melted$legend_group <- factor(data_melted$AF_by, levels = c("Radiotherapy", "Chemotherapy"))
# Create a factor with the desired order
data_melted$SN_types <- factor(data_melted$SN_types, levels = order)


# Create the grouped bar chart
p_ccss <- ggplot(data_melted, aes(x = SN_types, y = new_value, fill = legend_group)) +
  geom_bar_pattern(
    stat = "identity",
    position = position_dodge(width = 0.8),
    width = 0.8,
    pattern = ifelse(data_melted$legend_group == "Radiotherapy", "stripe", "none"),
    pattern_density = 0.1,
    pattern_spacing = 0.01,
    pattern_fill = "black",
    pattern_angle = 45,
    colour = "black"
  ) +
  scale_pattern_manual(
    values = c("Radiotherapy" = "stripe", "Chemotherapy" = "none"),
    guide = guide_legend(override.aes = list(fill = "gray"))
  ) +
  scale_fill_manual(
    values = custom_colors,
    breaks = legend_order,
    labels = c(
      "Radiotherapy" = "Higher exposure levels\nof radiotherapy",
      "Chemotherapy" = "Higher exposure levels\nof chemotherapy"
    ),
    guide = guide_legend(override.aes = list(pattern = c("Radiotherapy" = "stripe", "Chemotherapy" = "none")))
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 16, color = "black"),
        axis.text.y = element_text(size = 16, color = "black"),
        axis.title.y = element_text(size = 20, color = "black"), 
        axis.line.x = element_line(color = "black"), 
        axis.line.y = element_line(color = "black"),
        legend.title = element_blank(),
        legend.text = element_text(size = 16, color = "black"),
        plot.title = element_text(size = 16, hjust = 0.5, vjust = 1.5, face = "bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "top", 
        legend.box = "horizontal") +
  scale_x_discrete(limits = levels(data_melted$SN_types)) +
  coord_cartesian(clip = "off") +
  labs(title = "", y = "Attributable fraction (%)", x = NULL) +
  geom_text(
    aes(label = new_value, y = new_value),  # Adjust y position
    position = position_dodge(width = 0.8),
    vjust = -0.5,  # Adjust vertical justification
    size = 5.5,
    color = "black"
  ) +
  scale_y_continuous(limits = c(0, 100), expand = c(0, 0))

# p
p_ccss <- p_ccss + scale_y_continuous(limits = c(0, 100), expand = c(0, 0)) + 
  labs(title = "", y = "Attributable fraction (%)", x =  NULL) +
  theme(axis.title.x = element_text(size = 18, color = "black")) +
  theme(axis.title.y = element_text(size = 18, color = "black"))

# p_ccss



p <- cowplot::plot_grid(p_sjlife, p_ccss, labels = c("AUTO"))


# p
# Save the plot as a high-resolution image
plot_name <- paste0("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/ANALYSIS/results/plots/v21b/AFR/", group, "_", AF.type,"_", lifestyle, ".tiff")

ggsave(plot_name, p, width = 12, height = 8, dpi = 600, device = "tiff", compression = "lzw")