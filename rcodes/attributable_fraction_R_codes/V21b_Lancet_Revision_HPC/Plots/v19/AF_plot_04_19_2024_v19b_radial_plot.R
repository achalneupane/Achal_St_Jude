# Load required libraries
library(ggplot2)
library(ggthemes)# For additional themes
library(ggrepel)
library(dplyr)



## V19 b (Without lifestyle)
data <- read.table(text="Cohort	SN_types	Variables	Overall	Female	Male	Age.lt.35	Age.ge.35
SJLIFE	SNs (605)	Radiation	0.425	0.417	0.435	0.400	0.447
SJLIFE	SNs (605)	Chemotherapy	0.080	0.079	0.080	0.109	0.054
SJLIFE	SNs (605)	All_treatments	0.469	0.462	0.478	0.463	0.474
SJLIFE	SNs (605)	PRS	0.125	0.123	0.129	0.126	0.125
SJLIFE	SNs (605)	Lifestyle	-	-	-	-	-
SJLIFE	SNs (605)	Combined	0.536	0.530	0.545	0.531	0.541
SJLIFE	SMNs (462)	Radiation	0.372	0.371	0.373	0.340	0.396
SJLIFE	SMNs (462)	Chemotherapy	0.033	0.033	0.032	0.047	0.022
SJLIFE	SMNs (462)	All_treatments	0.392	0.392	0.392	0.370	0.409
SJLIFE	SMNs (462)	PRS	0.122	0.119	0.125	0.122	0.122
SJLIFE	SMNs (462)	Lifestyle	-	-	-	-	-
SJLIFE	SMNs (462)	Combined	0.467	0.466	0.468	0.447	0.482
SJLIFE	NMSC (249)	Radiation	0.437	0.427	0.448	0.408	0.453
SJLIFE	NMSC (249)	Chemotherapy	-	-	-	-	-
SJLIFE	NMSC (249)	All_treatments	0.437	0.427	0.448	0.408	0.453
SJLIFE	NMSC (249)	PRS	0.276	0.263	0.292	0.274	0.278
SJLIFE	NMSC (249)	Lifestyle	-	-	-	-	-
SJLIFE	NMSC (249)	Combined	0.591	0.579	0.605	0.569	0.604
SJLIFE	Breast cancer (74)	Radiation	0.493	0.493	-	0.475	0.498
SJLIFE	Breast cancer (74)	Chemotherapy	0.193	0.193	-	0.228	0.184
SJLIFE	Breast cancer (74)	All_treatments	0.602	0.602	-	0.594	0.604
SJLIFE	Breast cancer (74)	PRS	0.255	0.255	-	0.261	0.254
SJLIFE	Breast cancer (74)	Lifestyle	-	-	-	-	-
SJLIFE	Breast cancer (74)	Combined	0.704	0.704	-	0.701	0.705
SJLIFE	Thyroid cancer (86)	Radiation	0.62	0.622	0.616	0.594	0.65
SJLIFE	Thyroid cancer (86)	Chemotherapy	0.233	0.235	0.231	0.293	0.163
SJLIFE	Thyroid cancer (86)	All_treatments	0.726	0.727	0.723	0.725	0.727
SJLIFE	Thyroid cancer (86)	PRS	0.517	0.519	0.514	0.519	0.514
SJLIFE	Thyroid cancer (86)	Lifestyle	-	-	-	-	-
SJLIFE	Thyroid cancer (86)	Combined	0.866	0.867	0.866	0.866	0.868
SJLIFE	Meningioma (149)	Radiation	0.197	0.172	0.23	0.223	0.176
SJLIFE	Meningioma (149)	Chemotherapy	0.361	0.371	0.349	0.472	0.271
SJLIFE	Meningioma (149)	All_treatments	0.489	0.478	0.504	0.598	0.4
SJLIFE	Meningioma (149)	PRS	-	-	-	-	-
SJLIFE	Meningioma (149)	Lifestyle	-	-	-	-	-
SJLIFE	Meningioma (149)	Combined	0.425	0.416	0.437	0.549	0.323
SJLIFE	Sarcoma (32)	Radiation	-	-	-	-	-
SJLIFE	Sarcoma (32)	Chemotherapy	0.324	0.329	0.32	0.308	0.351
SJLIFE	Sarcoma (32)	All_treatments	0.324	0.329	0.32	0.308	0.351
SJLIFE	Sarcoma (32)	PRS	-	-	-	-	-
SJLIFE	Sarcoma (32)	Lifestyle	-	-	-	-	-
SJLIFE	Sarcoma (32)	Combined	0.311	0.315	0.306	0.294	0.337
CCSS	SNs (1611)	Radiation	0.387	0.380	0.398	0.372	0.398
CCSS	SNs (1611)	Chemotherapy	0.031	0.030	0.033	0.047	0.019
CCSS	SNs (1611)	All_treatments	0.407	0.400	0.418	0.403	0.410
CCSS	SNs (1611)	PRS	0.048	0.047	0.049	0.047	0.048
CCSS	SNs (1611)	Lifestyle	-	-	-	-	-
CCSS	SNs (1611)	Combined	0.435	0.428	0.446	0.431	0.438
CCSS	SMNs (762)	Radiation	0.255	0.254	0.257	0.213	0.283
CCSS	SMNs (762)	Chemotherapy	0.042	0.042	0.042	0.070	0.024
CCSS	SMNs (762)	All_treatments	0.290	0.289	0.293	0.271	0.303
CCSS	SMNs (762)	PRS	0.047	0.046	0.047	0.046	0.047
CCSS	SMNs (762)	Lifestyle	-	-	-	-	-
CCSS	SMNs (762)	Combined	0.323	0.322	0.326	0.305	0.336
CCSS	NMSC (769)	Radiation	0.393	0.388	0.398	0.387	0.396
CCSS	NMSC (769)	Chemotherapy	-	-	-	-	-
CCSS	NMSC (769)	All_treatments	0.393	0.388	0.398	0.387	0.396
CCSS	NMSC (769)	PRS	0.31	0.308	0.313	0.312	0.309
CCSS	NMSC (769)	Lifestyle	-	-	-	-	-
CCSS	NMSC (769)	Combined	0.582	0.578	0.587	0.58	0.583
CCSS	Breast cancer (289)	Radiation	0.474	0.474	-	0.452	0.478
CCSS	Breast cancer (289)	Chemotherapy	0.187	0.187	-	0.224	0.180
CCSS	Breast cancer (289)	All_treatments	0.593	0.593	-	0.589	0.593
CCSS	Breast cancer (289)	PRS	0.365	0.365	-	0.369	0.365
CCSS	Breast cancer (289)	Lifestyle	-	-	-	-	-
CCSS	Breast cancer (289)	Combined	0.743	0.743	-	0.741	0.743
CCSS	Thyroid cancer (163)	Radiation	0.443	0.445	0.439	0.413	0.489
CCSS	Thyroid cancer (163)	Chemotherapy	0.055	0.055	0.056	0.073	0.029
CCSS	Thyroid cancer (163)	All_treatments	0.476	0.477	0.472	0.457	0.504
CCSS	Thyroid cancer (163)	PRS	0.358	0.362	0.352	0.355	0.363
CCSS	Thyroid cancer (163)	Lifestyle	-	-	-	-	-
CCSS	Thyroid cancer (163)	Combined	0.664	0.665	0.661	0.65	0.684
CCSS	Meningioma (255)	Radiation	0.37	0.333	0.422	0.394	0.345
CCSS	Meningioma (255)	Chemotherapy	0.082	0.076	0.092	0.118	0.045
CCSS	Meningioma (255)	All_treatments	0.426	0.391	0.475	0.471	0.378
CCSS	Meningioma (255)	PRS	0.032	0.032	0.033	0.036	0.028
CCSS	Meningioma (255)	Lifestyle	-	-	-	-	-
CCSS	Meningioma (255)	Combined	0.444	0.41	0.494	0.49	0.397
CCSS	Sarcoma (60)	Radiation	-	-	-	-	-
CCSS	Sarcoma (60)	Chemotherapy	0.342	0.327	0.358	0.337	0.348
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
custom_order_within_variable <- c("SNs", "SMNs", "NMSC", "Breast cancer", "Thyroid cancer", "Meningioma", "Sarcoma")

# Create a custom factor variable based on the desired order
data_melted <- data_melted %>%
  dplyr::mutate(SN_types = factor(SN_types, levels = custom_order_within_variable))
data_melted <- data_melted %>%
  arrange(SN_types, AF_by)

data <- cbind.data.frame(individual= paste0(data_melted$new_value, "%"), group= data_melted$SN_types, group2= data_melted$legend_group, value=data_melted$new_value)
# data$value[is.na(data$value)]

custom_colors <- c("SJLIFE" = "darkred", "CCSS" = "darkblue")
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
p = ggplot(data, aes(x=as.factor(id), y=value, fill=group2)) +       # Note that id is a factor. If x is numeric, there is some space between the first bar
  
  geom_bar(aes(x=as.factor(id), y=value, fill=group2), stat="identity", alpha=1, width=0.95) +
  
  # # Add a val=100/75/50/25 lines. I do it at the beginning to make sur barplots are OVER it.
  # geom_segment(data=grid_data, aes(x = end, y = 80, xend = start, yend = 80), colour = "white", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  # geom_segment(data=grid_data, aes(x = end, y = 60, xend = start, yend = 60), colour = "white", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  # geom_segment(data=grid_data, aes(x = end, y = 40, xend = start, yend = 40), colour = "white", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  # geom_segment(data=grid_data, aes(x = end, y = 20, xend = start, yend = 20), colour = "white", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  
  ## Add text showing the value of each 100/75/50/25 lines
  # annotate("text", x = rep(max(data$id),5), y = c(0, 25, 50, 75, 100), label = c("0%", "25%", "50%", "75%", "100%") , color="black", size=3 , angle=0, fontface= 1, hjust=1) +
  
  geom_bar(aes(x=as.factor(id), y=value, fill=group2), stat="identity", alpha=1, width=0.95) +
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
plot_name <- "Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/ANALYSIS/results/plots/v19b/overall/Overall_Combined_without_lifestyle_radial_plot.tiff"
ggsave(plot_name, p, width = 12, height = 8, dpi = 600, device = "tiff", compression = "lzw")


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
custom_order_within_variable <- c("SNs", "SMNs", "NMSC", "Breast cancer", "Thyroid cancer", "Meningioma", "Sarcoma")

# Create a custom factor variable based on the desired order
data_melted <- data_melted %>%
  dplyr::mutate(SN_types = factor(SN_types, levels = custom_order_within_variable))
data_melted <- data_melted %>%
  arrange(SN_types, factor(AF_by, levels = c("All_treatments", "PRS")))

data <- cbind.data.frame(individual= paste0(data_melted$new_value, "%"), group= data_melted$SN_types, group2= data_melted$legend_group, value=data_melted$new_value)
# data$value[is.na(data$value)]

custom_colors <- c("SJLIFE-All treatments" = "darkred", "CCSS-All treatments" = "darkblue", "SJLIFE-PRS" = "red", "CCSS-PRS" = "blue")
legend_order <- c("SJLIFE-All treatments", "CCSS-All treatments", "SJLIFE-PRS", "CCSS-PRS")

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
p = ggplot(data, aes(x=as.factor(id), y=value, fill=group2)) +       # Note that id is a factor. If x is numeric, there is some space between the first bar
  
  geom_bar(aes(x=as.factor(id), y=value, fill=group2), stat="identity", alpha=1, width=0.95) +
  
  # # Add a val=100/75/50/25 lines. I do it at the beginning to make sur barplots are OVER it.
  # geom_segment(data=grid_data, aes(x = end, y = 80, xend = start, yend = 80), colour = "white", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  # geom_segment(data=grid_data, aes(x = end, y = 60, xend = start, yend = 60), colour = "white", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  # geom_segment(data=grid_data, aes(x = end, y = 40, xend = start, yend = 40), colour = "white", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  # geom_segment(data=grid_data, aes(x = end, y = 20, xend = start, yend = 20), colour = "white", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  
  ## Add text showing the value of each 100/75/50/25 lines
  # annotate("text", x = rep(max(data$id),5), y = c(0, 25, 50, 75, 100), label = c("0%", "25%", "50%", "75%", "100%") , color="black", size=3 , angle=0, fontface= 1, hjust=1) +
  
  geom_bar(aes(x=as.factor(id), y=value, fill=group2), stat="identity", alpha=1, width=0.95) +
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
    fill = "Cohorts and variables"
  )
p
plot_name <- "Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/ANALYSIS/results/plots/v19b/overall/Overall_treatment_and_PRS_without_lifestyle_radial_plot.tiff"
ggsave(plot_name, p, width = 12, height = 8, dpi = 600, device = "tiff", compression = "lzw")

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

data_melted <-new_data
order <- unique(data_melted$SN_types)
AF.type <- "Chemotherapy_Radiation"
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
  arrange(SN_types, factor(AF_by, levels = c("Radiation", "Chemotherapy")))

data <- cbind.data.frame(individual= paste0(data_melted$new_value, "%"), group= data_melted$SN_types, group2= data_melted$legend_group, value=data_melted$new_value)
# data$value[is.na(data$value)]

custom_colors <- c("SJLIFE-Radiation" = "darkred", "CCSS-Radiation" = "darkblue", "SJLIFE-Chemotherapy" = "red", "CCSS-Chemotherapy" = "blue")
legend_order <- c("SJLIFE-Radiation", "CCSS-Radiation", "SJLIFE-Chemotherapy", "CCSS-Chemotherapy")


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
p = ggplot(data, aes(x=as.factor(id), y=value, fill=group2)) +       # Note that id is a factor. If x is numeric, there is some space between the first bar
  
  geom_bar(aes(x=as.factor(id), y=value, fill=group2), stat="identity", alpha=1, width=0.95) +
  
  # # Add a val=100/75/50/25 lines. I do it at the beginning to make sur barplots are OVER it.
  # geom_segment(data=grid_data, aes(x = end, y = 80, xend = start, yend = 80), colour = "white", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  # geom_segment(data=grid_data, aes(x = end, y = 60, xend = start, yend = 60), colour = "white", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  # geom_segment(data=grid_data, aes(x = end, y = 40, xend = start, yend = 40), colour = "white", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  # geom_segment(data=grid_data, aes(x = end, y = 20, xend = start, yend = 20), colour = "white", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  
  ## Add text showing the value of each 100/75/50/25 lines
  # annotate("text", x = rep(max(data$id),5), y = c(0, 25, 50, 75, 100), label = c("0%", "25%", "50%", "75%", "100%") , color="black", size=3 , angle=0, fontface= 1, hjust=1) +
  
  geom_bar(aes(x=as.factor(id), y=value, fill=group2), stat="identity", alpha=1, width=0.95) +
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
    plot.margin = margin(r = -3.5, l = -3.5, t = -3.5, b = -3.5, unit = "cm"),
    legend.text = element_text(size = 14, color = "black")
  ) + # Adjust legend text size
  labs(
    title = "",
    x = NULL, y = "",
    fill = "Cohorts and variables"
  )
p
plot_name <- "Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/ANALYSIS/results/plots/v19b/overall/Overall_Chemotherapy_and_Radiation_without_lifestyle_radial_plot.tiff"
ggsave(plot_name, p, width = 12, height = 8, dpi = 600, device = "tiff", compression = "lzw")

