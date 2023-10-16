# Load required libraries
library(ggplot2)
library(ggthemes)  # For additional themes

## Create a data frame with your data
# data <- data.frame(
#   SN_types = c("Any SN", "SMN", "NMSC", "Breast cancer", "Thyroid cancer", "Meningioma", "Sarcoma"),
#   SJLIFE = c(0.518, 0.455, 0.652, 0.672, 0.862, 0.438, 0.371),
#   CCSS = c(0.419, 0.326, 0.555, 0.709, 0.659, 0.432, 0.368)
# )


## V18 b
data <- read.table(text="Cohort	SN_types	Variables	Overall	Female	Male	Age.lt.35	Age.ge.35
SJLIFE	Any SN (605)	Radiation	0.425	0.417	0.435	0.400	0.447
SJLIFE	Any SN (605)	Chemo	0.080	0.079	0.080	0.109	0.054
SJLIFE	Any SN (605)	All treatments	0.469	0.462	0.478	0.463	0.474
SJLIFE	Any SN (605)	PRS	0.125	0.123	0.129	0.126	0.125
SJLIFE	Any SN (605)	Lifestyle	-	-	-	-	-
SJLIFE	Any SN (605)	Combined	0.536	0.530	0.545	0.531	0.541
SJLIFE	SMN (462)	Radiation	0.372	0.371	0.373	0.340	0.396
SJLIFE	SMN (462)	Chemo	0.033	0.033	0.032	0.047	0.022
SJLIFE	SMN (462)	All treatments	0.392	0.392	0.392	0.370	0.409
SJLIFE	SMN (462)	PRS	0.122	0.119	0.125	0.122	0.122
SJLIFE	SMN (462)	Lifestyle	-	-	-	-	-
SJLIFE	SMN (462)	Combined	0.467	0.466	0.468	0.447	0.482
SJLIFE	NMSC (249)	Radiation	0.441	0.431	0.453	0.413	0.457
SJLIFE	NMSC (249)	Chemo	-	-	-	-	-
SJLIFE	NMSC (249)	All treatments	0.441	0.431	0.453	0.413	0.457
SJLIFE	NMSC (249)	PRS	0.441	0.432	0.452	0.440	0.442
SJLIFE	NMSC (249)	Lifestyle	-	-	-	-	-
SJLIFE	NMSC (249)	Combined	0.688	0.680	0.698	0.671	0.698
SJLIFE	Breast cancer (74)	Radiation	0.493	0.493	-	0.475	0.498
SJLIFE	Breast cancer (74)	Chemo	0.193	0.193	-	0.228	0.184
SJLIFE	Breast cancer (74)	All treatments	0.602	0.602	-	0.594	0.604
SJLIFE	Breast cancer (74)	PRS	0.255	0.255	-	0.261	0.254
SJLIFE	Breast cancer (74)	Lifestyle	-	-	-	-	-
SJLIFE	Breast cancer (74)	Combined	0.704	0.704	-	0.701	0.705
SJLIFE	Thyroid cancer (86)	Radiation	0.620	0.622	0.616	0.594	0.650
SJLIFE	Thyroid cancer (86)	Chemo	0.233	0.235	0.231	0.293	0.163
SJLIFE	Thyroid cancer (86)	All treatments	0.726	0.727	0.723	0.725	0.727
SJLIFE	Thyroid cancer (86)	PRS	0.517	0.519	0.514	0.519	0.514
SJLIFE	Thyroid cancer (86)	Lifestyle	-	-	-	-	-
SJLIFE	Thyroid cancer (86)	Combined	0.866	0.867	0.866	0.866	0.868
SJLIFE	Meningioma (149)	Radiation	0.197	0.172	0.230	0.223	0.176
SJLIFE	Meningioma (149)	Chemo	0.361	0.371	0.349	0.472	0.271
SJLIFE	Meningioma (149)	All treatments	0.489	0.478	0.504	0.598	0.400
SJLIFE	Meningioma (149)	PRS	-0.125	-0.117	-0.135	-0.121	-0.128
SJLIFE	Meningioma (149)	Lifestyle	-	-	-	-	-
SJLIFE	Meningioma (149)	Combined	0.425	0.416	0.437	0.549	0.323
SJLIFE	Sarcoma (32)	Radiation	-	-	-	-	-
SJLIFE	Sarcoma (32)	Chemo	0.324	0.329	0.320	0.308	0.351
SJLIFE	Sarcoma (32)	All treatments	0.324	0.329	0.320	0.308	0.351
SJLIFE	Sarcoma (32)	PRS	-0.020	-0.021	-0.020	-0.020	-0.020
SJLIFE	Sarcoma (32)	Lifestyle	-	-	-	-	-
SJLIFE	Sarcoma (32)	Combined	0.311	0.315	0.306	0.294	0.337
CCSS	Any SN (1611)	Radiation	0.387	0.380	0.398	0.372	0.398
CCSS	Any SN (1611)	Chemo	0.031	0.030	0.033	0.047	0.019
CCSS	Any SN (1611)	All treatments	0.407	0.400	0.418	0.403	0.410
CCSS	Any SN (1611)	PRS	0.048	0.047	0.049	0.047	0.048
CCSS	Any SN (1611)	Lifestyle	-	-	-	-	-
CCSS	Any SN (1611)	Combined	0.435	0.428	0.446	0.431	0.438
CCSS	SMN (762)	Radiation	0.255	0.254	0.257	0.213	0.283
CCSS	SMN (762)	Chemo	0.042	0.042	0.042	0.070	0.024
CCSS	SMN (762)	All treatments	0.290	0.289	0.293	0.271	0.303
CCSS	SMN (762)	PRS	0.047	0.046	0.047	0.046	0.047
CCSS	SMN (762)	Lifestyle	-	-	-	-	-
CCSS	SMN (762)	Combined	0.323	0.322	0.326	0.305	0.336
CCSS	NMSC (769)	Radiation	0.391	0.386	0.396	0.385	0.394
CCSS	NMSC (769)	Chemo	-	-	-	-	-
CCSS	NMSC (769)	All treatments	0.391	0.386	0.396	0.385	0.394
CCSS	NMSC (769)	PRS	0.304	0.306	0.302	0.306	0.303
CCSS	NMSC (769)	Lifestyle	-	-	-	-	-
CCSS	NMSC (769)	Combined	0.577	0.575	0.578	0.574	0.578
CCSS	Breast cancer (289)	Radiation	0.474	0.474	-	0.452	0.478
CCSS	Breast cancer (289)	Chemo	0.187	0.187	-	0.224	0.180
CCSS	Breast cancer (289)	All treatments	0.593	0.593	-	0.589	0.593
CCSS	Breast cancer (289)	PRS	0.365	0.365	-	0.369	0.365
CCSS	Breast cancer (289)	Lifestyle	-	-	-	-	-
CCSS	Breast cancer (289)	Combined	0.743	0.743	-	0.741	0.743
CCSS	Thyroid cancer (163)	Radiation	0.443	0.445	0.439	0.413	0.489
CCSS	Thyroid cancer (163)	Chemo	0.055	0.055	0.056	0.073	0.029
CCSS	Thyroid cancer (163)	All treatments	0.476	0.477	0.472	0.457	0.504
CCSS	Thyroid cancer (163)	PRS	0.358	0.362	0.352	0.355	0.363
CCSS	Thyroid cancer (163)	Lifestyle	-	-	-	-	-
CCSS	Thyroid cancer (163)	Combined	0.664	0.665	0.661	0.650	0.684
CCSS	Meningioma (255)	Radiation	0.370	0.333	0.422	0.394	0.345
CCSS	Meningioma (255)	Chemo	0.082	0.076	0.092	0.118	0.045
CCSS	Meningioma (255)	All treatments	0.426	0.391	0.475	0.471	0.378
CCSS	Meningioma (255)	PRS	0.032	0.032	0.033	0.036	0.028
CCSS	Meningioma (255)	Lifestyle	-	-	-	-	-
CCSS	Meningioma (255)	Combined	0.444	0.410	0.494	0.490	0.397
CCSS	Sarcoma (60)	Radiation	-	-	-	-	-
CCSS	Sarcoma (60)	Chemo	0.342	0.327	0.358	0.337	0.348
CCSS	Sarcoma (60)	All treatments	0.342	0.327	0.358	0.337	0.348
CCSS	Sarcoma (60)	PRS	0.024	0.026	0.022	0.023	0.025
CCSS	Sarcoma (60)	Lifestyle	-	-	-	-	-
CCSS	Sarcoma (60)	Combined	0.358	0.345	0.372	0.353	0.364", header = T, sep = "\t")

data$SN_types <- gsub("\\([0-9]+\\)", "", data$SN_types)
data[data == "-"] <- NA

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
# 2             SMN   SJLIFE 0.455
# 3            NMSC   SJLIFE 0.652
data_melted <-new_data
## Exclude rows with NAs
data_melted <- data_melted[complete.cases(data_melted),]

# order <- c("Any SN", "SMN", "NMSC", "Breast cancer", "Thyroid cancer", "Meningioma", "Sarcoma")
order <- unique(data_melted$SN_types)
AF.type <- data_melted$AF_by[1] 
lifestyle <- "without_lifestyle"
if(lifestyle== "without_lifestyle" && variables[i] == "Lifestyle"){
  next
}

# Create the grouped bar chart
p <- ggplot(data_melted, aes(x = SN_types, y = value, fill = variable)) +
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
  scale_fill_manual(values = c("SJLIFE" = "#0072B2", "CCSS" = "#D55E00")) +  # Use color codes
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
  # scale_x_discrete(limits = c("Any SN ", "SMN ", "NMSC ", "Breast cancer ", "Thyroid cancer ", "Meningioma ", "Sarcoma ")) +
  scale_x_discrete(limits = order) +
  ## Adjust the x-axis limits to ensure that the secondary y-axis label is fully visible.
  coord_cartesian(clip = "off")

# Save the plot as a high-resolution image
plot_name <- paste0("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/ANALYSIS/results/plots/v18b/", group, "_", AF.type, "_", lifestyle, ".tiff")

ggsave(plot_name, p, width = 10, height = 6, dpi = 300)
}
