library(ggplot2)

# Your data
df <- read.table(text = "Cohort	SN_types	Variables	Overall	Female	Male	Age<35	Age≥35
SJLIFE	Any SN (605)	Radiation	0.432	0.426	0.441	0.322	0.456
SJLIFE	Any SN (605)	Chemo	0.04	0.041	0.04	0.072	0.034
SJLIFE	Any SN (605)	All treatments	0.454	0.448	0.461	0.368	0.472
SJLIFE	Any SN (605)	PRS	0.118	0.115	0.12	0.117	0.118
SJLIFE	Any SN (605)	Lifestyle	NA	NA	NA	NA	NA
SJLIFE	Any SN (605)	Combined	0.518	0.513	0.526	0.443	0.534
SJLIFE	SMN (462)	Radiation	0.383	0.382	0.383	0.268	0.406
SJLIFE	SMN (462)	Chemo	0.003	0.003	0.003	0.012	0.001
SJLIFE	SMN (462)	All treatments	0.384	0.384	0.383	0.275	0.406
SJLIFE	SMN (462)	PRS	0.115	0.113	0.117	0.113	0.115
SJLIFE	SMN (462)	Lifestyle	NA	NA	NA	NA	NA
SJLIFE	SMN (462)	Combined	0.455	0.454	0.455	0.357	0.474
SJLIFE	NMSC (249)	Radiation	0.392	0.384	0.4	0.31	0.399
SJLIFE	NMSC (249)	Chemo	NA	NA	NA	NA	NA
SJLIFE	NMSC (249)	All treatments	0.392	0.384	0.4	0.31	0.399
SJLIFE	NMSC (249)	PRS	0.428	0.417	0.44	0.421	0.428
SJLIFE	NMSC (249)	Lifestyle	NA	NA	NA	NA	NA
SJLIFE	NMSC (249)	Combined	0.652	0.644	0.662	0.601	0.657
SJLIFE	Breast cancer (76)	Radiation	0.509	0.532	0.486	0.446	0.513
SJLIFE	Breast cancer (76)	Chemo	0.075	0.076	0.075	0.06	0.076
SJLIFE	Breast cancer (76)	All treatments	0.55	0.572	0.53	0.486	0.554
SJLIFE	Breast cancer (76)	PRS	0.27	0.276	0.265	0.269	0.27
SJLIFE	Breast cancer (76)	Lifestyle	NA	NA	NA	NA	NA
SJLIFE	Breast cancer (76)	Combined	0.672	0.689	0.655	0.622	0.675
SJLIFE	Thyroid cancer (86)	Radiation	0.624	0.626	0.622	0.545	0.647
SJLIFE	Thyroid cancer (86)	Chemo	0.2	0.203	0.195	0.312	0.168
SJLIFE	Thyroid cancer (86)	All treatments	0.714	0.716	0.712	0.68	0.724
SJLIFE	Thyroid cancer (86)	PRS	0.52	0.523	0.517	0.525	0.519
SJLIFE	Thyroid cancer (86)	Lifestyle	NA	NA	NA	NA	NA
SJLIFE	Thyroid cancer (86)	Combined	0.862	0.863	0.861	0.844	0.867
SJLIFE	Meningioma (149)	Radiation	0.196	0.176	0.224	0.266	0.183
SJLIFE	Meningioma (149)	Chemo	0.379	0.392	0.361	0.563	0.342
SJLIFE	Meningioma (149)	All treatments	0.501	0.492	0.512	0.704	0.46
SJLIFE	Meningioma (149)	PRS	-0.125	-0.116	-0.137	-0.13	-0.124
SJLIFE	Meningioma (149)	Lifestyle	NA	NA	NA	NA	NA
SJLIFE	Meningioma (149)	Combined	0.438	0.433	0.446	0.666	0.393
SJLIFE	Sarcoma (32)	Radiation	NA	NA	NA	NA	NA
SJLIFE	Sarcoma (32)	Chemo	0.367	0.37	0.365	0.319	0.395
SJLIFE	Sarcoma (32)	All treatments	0.367	0.37	0.365	0.319	0.395
SJLIFE	Sarcoma (32)	PRS	0.006	0.006	0.007	0.007	0.006
SJLIFE	Sarcoma (32)	Lifestyle	NA	NA	NA	NA	NA
SJLIFE	Sarcoma (32)	Combined	0.371	0.374	0.369	0.324	0.398
CCSS	Any SN (1611)	Radiation	0.378	0.372	0.387	0.3	0.397
CCSS	Any SN (1611)	Chemo	0.017	0.016	0.018	0.04	0.011
CCSS	Any SN (1611)	All treatments	0.389	0.383	0.399	0.328	0.404
CCSS	Any SN (1611)	PRS	0.053	0.052	0.053	0.052	0.053
CCSS	Any SN (1611)	Lifestyle	NA	NA	NA	NA	NA
CCSS	Any SN (1611)	Combined	0.421	0.415	0.43	0.363	0.435
CCSS	SMN (762)	Radiation	0.268	0.269	0.265	0.171	0.293
CCSS	SMN (762)	Chemo	0.023	0.023	0.023	0.055	0.015
CCSS	SMN (762)	All treatments	0.287	0.288	0.286	0.217	0.305
CCSS	SMN (762)	PRS	0.056	0.056	0.057	0.055	0.056
CCSS	SMN (762)	Lifestyle	NA	NA	NA	NA	NA
CCSS	SMN (762)	Combined	0.327	0.328	0.326	0.26	0.344
CCSS	NMSC (769)	Radiation	0.375	0.371	0.381	0.299	0.387
CCSS	NMSC (769)	Chemo	NA	NA	NA	NA	NA
CCSS	NMSC (769)	All treatments	0.375	0.371	0.381	0.299	0.387
CCSS	NMSC (769)	PRS	0.289	0.291	0.287	0.29	0.289
CCSS	NMSC (769)	Lifestyle	NA	NA	NA	NA	NA
CCSS	NMSC (769)	Combined	0.557	0.555	0.558	0.504	0.565
CCSS	Breast cancer (294)	Radiation	0.469	0.474	0.464	0.374	0.476
CCSS	Breast cancer (294)	Chemo	0.132	0.125	0.14	0.174	0.129
CCSS	Breast cancer (294)	All treatments	0.555	0.556	0.554	0.484	0.56
CCSS	Breast cancer (294)	PRS	0.319	0.318	0.32	0.316	0.319
CCSS	Breast cancer (294)	Lifestyle	NA	NA	NA	NA	NA
CCSS	Breast cancer (294)	Combined	0.697	0.698	0.697	0.646	0.701
CCSS	Thyroid cancer (163)	Radiation	0.448	0.448	0.448	0.345	0.484
CCSS	Thyroid cancer (163)	Chemo	0.041	0.041	0.042	0.096	0.022
CCSS	Thyroid cancer (163)	All treatments	0.473	0.472	0.474	0.405	0.497
CCSS	Thyroid cancer (163)	PRS	0.352	0.356	0.346	0.347	0.354
CCSS	Thyroid cancer (163)	Lifestyle	NA	NA	NA	NA	NA
CCSS	Thyroid cancer (163)	Combined	0.659	0.659	0.658	0.613	0.675
CCSS	Meningioma (255)	Radiation	0.365	0.333	0.411	0.385	0.359
CCSS	Meningioma (255)	Chemo	0.076	0.071	0.084	0.16	0.048
CCSS	Meningioma (255)	All treatments	0.418	0.387	0.462	0.487	0.394
CCSS	Meningioma (255)	PRS	0.027	0.026	0.027	0.031	0.025
CCSS	Meningioma (255)	Lifestyle	NA	NA	NA	NA	NA
CCSS	Meningioma (255)	Combined	0.433	0.402	0.477	0.503	0.41
CCSS	Sarcoma (60)	Radiation	NA	NA	NA	NA	NA
CCSS	Sarcoma (60)	Chemo	0.328	0.313	0.345	0.316	0.333
CCSS	Sarcoma (60)	All treatments	0.328	0.313	0.345	0.316	0.333
CCSS	Sarcoma (60)	PRS	0.024	0.027	0.022	0.021	0.025
CCSS	Sarcoma (60)	Lifestyle	NA	NA	NA	NA	NA
CCSS	Sarcoma (60)	Combined	0.345	0.331	0.36	0.331	0.35
", header = TRUE, sep = "\t", check.names = FALSE)

df <- df[!grepl("Lifestyle", df$Variables),]

df$n <- gsub("\\D", "", df$SN_types)
df$SN_types <- sub("\\s*\\(.*", "", df$SN_types)

# Melt the data to long format
df_melted <- reshape2::melt(df, id.vars = c("SN_types", "Cohort", "Variables", "n"))
df_melted$value <- df_melted$value *100

## Fix the plot positions
df_melted$SN_types <- factor(df_melted$SN_types, levels = c("Any SN", "SMN", "NMSC", "Breast cancer", "Thyroid cancer", "Meningioma", "Sarcoma"))
df_melted$Cohort <- factor(df_melted$Cohort, levels = c("SJLIFE", "CCSS"))

library(RColorBrewer)

# Define the number of colors you need (corresponding to the number of variables)
num_colors <- length(unique(df_melted$variable))
color_palette <- brewer.pal(num_colors, "Set1")

#########################################################
## With n
n_df <- df %>%
  group_by(SN_types, Cohort) %>%
  summarise(n = first(n))

ggplot(df_melted, aes(x = Variables, y = value, fill = factor(variable))) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.6) +
  labs(x = "Variables", y = "Attributable Fraction (%)", fill = "Variables") +
  scale_fill_manual(values = color_palette) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "top",
        strip.text = element_text(size = 10),  # Adjust the size of facet text
        strip.background = element_blank(),    # Remove facet background color
        strip.placement = "outside",           # Place facet labels outside the plot
        panel.border = element_rect(color = "black", fill = NA, size = 1)  # Add border between facets
  ) +
  facet_grid(SN_types ~ Cohort, scales = "free_x", space = "free_x", switch = "x", ) +
  # facet_grid(.~ Cohort, scales = "free_x", space = "free_x", switch = "x", ) +
  scale_y_continuous(labels = scales::percent_format(scale = 1),    # Format Y-axis labels as percentages
                     limits = c(0, 100), breaks = seq(0, 100, 30))


##########################################################
# Create the plot with facets on the left and right
ggplot(df_melted, aes(x = Variables, y = value, fill = factor(variable))) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.6) +
  labs(x = "Variables", y = "Attributable Fraction (%)", fill = "Variables") +
  scale_fill_manual(values = color_palette) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "top",
        strip.text = element_text(size = 10),  # Adjust the size of facet text
        strip.background = element_blank(),    # Remove facet background color
        strip.placement = "outside",           # Place facet labels outside the plot
        panel.border = element_rect(color = "black", fill = NA, size = 1)  # Add border between facets
  ) +
  facet_grid(SN_types ~ Cohort, scales = "free_x", space = "free_x", switch = "x", ) +
  # facet_grid(.~ Cohort, scales = "free_x", space = "free_x", switch = "x", ) +
  scale_y_continuous(labels = scales::percent_format(scale = 1),    # Format Y-axis labels as percentages
                     limits = c(0, 100), breaks = seq(0, 100, 30)) +
  guides(fill = guide_legend(nrow = 1)) +  # Limits the legend to one row
  geom_text(aes(label = paste0(round(value, 2), "%")),   # Add percentage symbol to the label
            position = position_dodge(width = 0.6), 
            hjust = -0.1,                                 # Adjust vertical position of the labels above the bars
            size = 3,
            angle = 90) 

##########################################################
# # Create the plot with facets on the left and right
# ggplot(df_melted, aes(x = Variables, y = value, fill = factor(variable))) +
#   geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.6) +
#   labs(x = "Variables", y = "Attributable Fraction", fill = "Variables") +
#   scale_fill_manual(values = color_palette) +
#   theme_minimal() +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1),
#         legend.position = "top",
#         strip.text = element_text(size = 10),  # Adjust the size of facet text
#         strip.background = element_blank(),    # Remove facet background color
#         strip.placement = "outside",           # Place facet labels outside the plot
#         panel.border = element_rect(color = "black", fill = NA, size = 1)  # Add border between facets
#   ) +
#   facet_grid(SN_types ~ Cohort, scales = "free_x", space = "free_x", switch = "x", ) +
#   # facet_grid(.~ Cohort, scales = "free_x", space = "free_x", switch = "x", ) +
#   scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 20)) +
#   guides(fill = guide_legend(nrow = 1)) +  # Limits the legend to one row
#   geom_text(aes(label = round(value, 2)), 
#             position = position_dodge(width = 0.8), 
#             hjust = -0.1, 
#             size = 3,
#             angle = 90)  # Add numeric labels above each bar and rotate them vertically


# ##########################################################
# # Create the plot with facets
# ggplot(df_melted, aes(x = Variables, y = value, fill = factor(variable))) +
#   geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.6) +
#   labs(x = "Variables", y = "Attributable Fraction", fill = "Variables") +
#   scale_fill_manual(values = color_palette) +
#   theme_minimal() +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1),
#         legend.position = "top",
#         strip.text = element_text(size = 10),  # Adjust the size of facet text
#         strip.background = element_blank(),    # Remove facet background color
#         strip.placement = "outside",           # Place facet labels outside the plot
#         panel.border = element_rect(color = "black", fill = NA, size = 1)  # Add border between facets
#   ) +
#   facet_grid(Cohort ~ SN_types, scales = "free_y", space = "free_y", switch = "y") +
#   scale_y_continuous(limits = c(-0.2, 1), breaks = seq(-0.2, 1, 0.2)) +
#   guides(fill = guide_legend(nrow = 1)) +  # Limits the legend to one row
#   geom_text(aes(label = round(value, 2)), 
#             position = position_dodge(width = 0.8), 
#             hjust = -0.1, 
#             size = 3,
#             angle = 90) 
