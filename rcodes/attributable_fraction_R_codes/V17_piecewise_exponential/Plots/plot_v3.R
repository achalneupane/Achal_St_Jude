# Load the required libraries
library(ggplot2)
library(scales)
library(ggstance)
library(dplyr)

setwd("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/ANALYSIS/results/plots")

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

###########################################################
# Create the plot with facets on the left and right and with %
P <- ggplot(df_melted, aes(x = Variables, y = value, fill = factor(variable))) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.2) +
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
            position = position_dodge(width = 0.8), 
            hjust = -0.1,                                 # Adjust vertical position of the labels above the bars
            size = 3,
            angle = 90) 

plotname <- gsub(" ","_", paste0("V17_plots_without_lifestyle", ".tiff"))
ggsave(filename = plotname, plot = P, width = 10, height = 10, dpi = 600)

##########################################################################
## V17 without diet
mydf <- read.table(text = "Lifestyles	Cohort	SN_types	Variables	Overall	Female	Male	Age<35	Age≥35
All_4_lifestyles	SJLIFE	Any SN (303)	Radiation	0.436	0.427	0.448	0.346	0.456
All_4_lifestyles	SJLIFE	Any SN (303)	Chemo	0.034	0.034	0.034	0.061	0.028
All_4_lifestyles	SJLIFE	Any SN (303)	All treatments	0.454	0.446	0.464	0.382	0.47
All_4_lifestyles	SJLIFE	Any SN (303)	PRS	0.184	0.181	0.187	0.186	0.183
All_4_lifestyles	SJLIFE	Any SN (303)	Lifestyle	-0.051	-0.047	-0.057	-0.043	-0.053
All_4_lifestyles	SJLIFE	Any SN (303)	Combined	0.531	0.526	0.537	0.47	0.544
All_4_lifestyles	SJLIFE	SMN (234)	Radiation	0.378	0.378	0.377	0.266	0.404
All_4_lifestyles	SJLIFE	SMN (234)	Chemo	0.027	0.027	0.028	0.045	0.023
All_4_lifestyles	SJLIFE	SMN (234)	All treatments	0.395	0.397	0.393	0.299	0.418
All_4_lifestyles	SJLIFE	SMN (234)	PRS	0.152	0.149	0.156	0.15	0.153
All_4_lifestyles	SJLIFE	SMN (234)	Lifestyle	-0.073	-0.063	-0.086	-0.066	-0.075
All_4_lifestyles	SJLIFE	SMN (234)	Combined	0.449	0.455	0.44	0.359	0.469
All_4_lifestyles	SJLIFE	NMSC (118)	Radiation	0.276	0.272	0.28	0.24	0.28
All_4_lifestyles	SJLIFE	NMSC (118)	Chemo	NA	NA	NA	NA	NA
All_4_lifestyles	SJLIFE	NMSC (118)	All treatments	0.276	0.272	0.28	0.24	0.28
All_4_lifestyles	SJLIFE	NMSC (118)	PRS	0.439	0.424	0.458	0.444	0.438
All_4_lifestyles	SJLIFE	NMSC (118)	Lifestyle	-0.342	-0.294	-0.402	-0.243	-0.352
All_4_lifestyles	SJLIFE	NMSC (118)	Combined	0.459	0.465	0.451	0.472	0.457
All_4_lifestyles	SJLIFE	Breast cancer (53)	Radiation	0.486	0.508	0.466	0.43	0.489
All_4_lifestyles	SJLIFE	Breast cancer (53)	Chemo	0.075	0.078	0.073	0.062	0.076
All_4_lifestyles	SJLIFE	Breast cancer (53)	All treatments	0.527	0.548	0.508	0.471	0.531
All_4_lifestyles	SJLIFE	Breast cancer (53)	PRS	0.201	0.194	0.208	0.22	0.2
All_4_lifestyles	SJLIFE	Breast cancer (53)	Lifestyle	-0.524	-0.586	-0.471	-0.393	-0.533
All_4_lifestyles	SJLIFE	Breast cancer (53)	Combined	0.42	0.424	0.418	0.402	0.422
All_4_lifestyles	SJLIFE	Thyroid cancer (43)	Radiation	0.609	0.632	0.582	0.55	0.625
All_4_lifestyles	SJLIFE	Thyroid cancer (43)	Chemo	0.255	0.256	0.254	0.379	0.22
All_4_lifestyles	SJLIFE	Thyroid cancer (43)	All treatments	0.732	0.748	0.714	0.715	0.737
All_4_lifestyles	SJLIFE	Thyroid cancer (43)	PRS	0.583	0.581	0.586	0.595	0.58
All_4_lifestyles	SJLIFE	Thyroid cancer (43)	Lifestyle	-0.159	-0.17	-0.146	-0.131	-0.166
All_4_lifestyles	SJLIFE	Thyroid cancer (43)	Combined	0.873	0.878	0.867	0.865	0.875
All_4_lifestyles	SJLIFE	Meningioma (81)	Radiation	0.187	0.15	0.232	0.283	0.173
All_4_lifestyles	SJLIFE	Meningioma (81)	Chemo	0.341	0.347	0.334	0.501	0.318
All_4_lifestyles	SJLIFE	Meningioma (81)	All treatments	0.474	0.453	0.499	0.676	0.445
All_4_lifestyles	SJLIFE	Meningioma (81)	PRS	-0.119	-0.114	-0.125	-0.117	-0.119
All_4_lifestyles	SJLIFE	Meningioma (81)	Lifestyle	-0.129	-0.12	-0.138	-0.114	-0.131
All_4_lifestyles	SJLIFE	Meningioma (81)	Combined	0.334	0.317	0.355	0.59	0.297
All_4_lifestyles	CCSS	Any SN (1404)	Radiation	0.37	0.364	0.378	0.299	0.384
All_4_lifestyles	CCSS	Any SN (1404)	Chemo	0.022	0.021	0.025	0.044	0.018
All_4_lifestyles	CCSS	Any SN (1404)	All treatments	0.385	0.378	0.394	0.329	0.396
All_4_lifestyles	CCSS	Any SN (1404)	PRS	0.062	0.062	0.063	0.062	0.062
All_4_lifestyles	CCSS	Any SN (1404)	Lifestyle	0.43	0.43	0.429	0.431	0.43
All_4_lifestyles	CCSS	Any SN (1404)	Combined	0.67	0.667	0.676	0.641	0.676
All_4_lifestyles	CCSS	SMN (662)	Radiation	0.275	0.276	0.274	0.185	0.294
All_4_lifestyles	CCSS	SMN (662)	Chemo	0.031	0.031	0.033	0.068	0.024
All_4_lifestyles	CCSS	SMN (662)	All treatments	0.302	0.301	0.304	0.242	0.315
All_4_lifestyles	CCSS	SMN (662)	PRS	0.066	0.065	0.066	0.064	0.066
All_4_lifestyles	CCSS	SMN (662)	Lifestyle	-0.039	-0.039	-0.038	-0.054	-0.036
All_4_lifestyles	CCSS	SMN (662)	Combined	0.321	0.32	0.324	0.25	0.336
All_4_lifestyles	CCSS	NMSC (707)	Radiation	0.38	0.376	0.386	0.312	0.389
All_4_lifestyles	CCSS	NMSC (707)	Chemo	NA	NA	NA	NA	NA
All_4_lifestyles	CCSS	NMSC (707)	All treatments	0.38	0.376	0.386	0.312	0.389
All_4_lifestyles	CCSS	NMSC (707)	PRS	0.29	0.293	0.287	0.292	0.29
All_4_lifestyles	CCSS	NMSC (707)	Lifestyle	0.669	0.671	0.667	0.674	0.669
All_4_lifestyles	CCSS	NMSC (707)	Combined	0.855	0.855	0.854	0.842	0.856
All_4_lifestyles	CCSS	Breast cancer (278)	Radiation	0.449	0.451	0.447	0.36	0.455
All_4_lifestyles	CCSS	Breast cancer (278)	Chemo	0.127	0.121	0.134	0.156	0.125
All_4_lifestyles	CCSS	Breast cancer (278)	All treatments	0.534	0.533	0.536	0.464	0.539
All_4_lifestyles	CCSS	Breast cancer (278)	PRS	0.32	0.318	0.322	0.321	0.32
All_4_lifestyles	CCSS	Breast cancer (278)	Lifestyle	1	1	1	1	1
All_4_lifestyles	CCSS	Breast cancer (278)	Combined	1	1	1	1	1
All_4_lifestyles	CCSS	Thyroid cancer (135)	Radiation	0.438	0.436	0.442	0.33	0.476
All_4_lifestyles	CCSS	Thyroid cancer (135)	Chemo	0.059	0.057	0.063	0.117	0.039
All_4_lifestyles	CCSS	Thyroid cancer (135)	All treatments	0.474	0.471	0.479	0.403	0.498
All_4_lifestyles	CCSS	Thyroid cancer (135)	PRS	0.363	0.367	0.354	0.353	0.366
All_4_lifestyles	CCSS	Thyroid cancer (135)	Lifestyle	-1.046	-1.052	-1.036	-1.047	-1.046
All_4_lifestyles	CCSS	Thyroid cancer (135)	Combined	0.312	0.307	0.321	0.209	0.348
All_4_lifestyles	CCSS	Meningioma (215)	Radiation	0.303	0.272	0.348	0.326	0.297
All_4_lifestyles	CCSS	Meningioma (215)	Chemo	0.082	0.074	0.094	0.17	0.06
All_4_lifestyles	CCSS	Meningioma (215)	All treatments	0.371	0.339	0.418	0.452	0.351
All_4_lifestyles	CCSS	Meningioma (215)	PRS	0.017	0.016	0.018	0.022	0.016
All_4_lifestyles	CCSS	Meningioma (215)	Lifestyle	0.105	0.103	0.108	0.085	0.11
All_4_lifestyles	CCSS	Meningioma (215)	Combined	0.447	0.416	0.491	0.508	0.432
All_4_lifestyles	CCSS	Sarcoma (42)	Radiation	NA	NA	NA	NA	NA
All_4_lifestyles	CCSS	Sarcoma (42)	Chemo	0.24	0.228	0.26	0.267	0.234
All_4_lifestyles	CCSS	Sarcoma (42)	All treatments	0.24	0.228	0.26	0.267	0.234
All_4_lifestyles	CCSS	Sarcoma (42)	PRS	0.063	0.066	0.058	0.057	0.065
All_4_lifestyles	CCSS	Sarcoma (42)	Lifestyle	1	1	1	1	1
All_4_lifestyles	CCSS	Sarcoma (42)	Combined	1	1	1	1	1
Smoking_only	SJLIFE	Any SN (303)	Radiation	0.438	0.428	0.45	0.347	0.458
Smoking_only	SJLIFE	Any SN (303)	Chemo	0.035	0.035	0.035	0.062	0.029
Smoking_only	SJLIFE	Any SN (303)	All treatments	0.456	0.448	0.466	0.384	0.472
Smoking_only	SJLIFE	Any SN (303)	PRS	0.183	0.181	0.187	0.186	0.183
Smoking_only	SJLIFE	Any SN (303)	Lifestyle	-0.041	-0.043	-0.039	-0.03	-0.044
Smoking_only	SJLIFE	Any SN (303)	Combined	0.538	0.53	0.548	0.48	0.551
Smoking_only	SJLIFE	SMN (234)	Radiation	0.378	0.378	0.377	0.266	0.404
Smoking_only	SJLIFE	SMN (234)	Chemo	0.03	0.029	0.031	0.049	0.026
Smoking_only	SJLIFE	SMN (234)	All treatments	0.398	0.399	0.396	0.302	0.42
Smoking_only	SJLIFE	SMN (234)	PRS	0.15	0.148	0.154	0.149	0.151
Smoking_only	SJLIFE	SMN (234)	Lifestyle	-0.032	-0.032	-0.031	-0.025	-0.033
Smoking_only	SJLIFE	SMN (234)	Combined	0.472	0.473	0.471	0.389	0.492
Smoking_only	SJLIFE	NMSC (118)	Radiation	0.275	0.28	0.268	0.243	0.278
Smoking_only	SJLIFE	NMSC (118)	Chemo	NA	NA	NA	NA	NA
Smoking_only	SJLIFE	NMSC (118)	All treatments	0.275	0.28	0.268	0.243	0.278
Smoking_only	SJLIFE	NMSC (118)	PRS	0.416	0.401	0.435	0.419	0.416
Smoking_only	SJLIFE	NMSC (118)	Lifestyle	0.015	0.016	0.013	0.013	0.015
Smoking_only	SJLIFE	NMSC (118)	Combined	0.586	0.58	0.594	0.566	0.588
Smoking_only	SJLIFE	Breast cancer (53)	Radiation	0.491	0.515	0.468	0.432	0.495
Smoking_only	SJLIFE	Breast cancer (53)	Chemo	0.054	0.057	0.051	0.04	0.055
Smoking_only	SJLIFE	Breast cancer (53)	All treatments	0.521	0.544	0.499	0.461	0.524
Smoking_only	SJLIFE	Breast cancer (53)	PRS	0.197	0.19	0.203	0.213	0.196
Smoking_only	SJLIFE	Breast cancer (53)	Lifestyle	-0.242	-0.262	-0.223	-0.192	-0.245
Smoking_only	SJLIFE	Breast cancer (53)	Combined	0.523	0.537	0.511	0.48	0.526
Smoking_only	SJLIFE	Thyroid cancer (43)	Radiation	0.611	0.633	0.586	0.555	0.627
Smoking_only	SJLIFE	Thyroid cancer (43)	Chemo	0.25	0.251	0.249	0.368	0.216
Smoking_only	SJLIFE	Thyroid cancer (43)	All treatments	0.73	0.746	0.713	0.714	0.735
Smoking_only	SJLIFE	Thyroid cancer (43)	PRS	0.575	0.575	0.576	0.587	0.572
Smoking_only	SJLIFE	Thyroid cancer (43)	Lifestyle	0.124	0.128	0.119	0.119	0.125
Smoking_only	SJLIFE	Thyroid cancer (43)	Combined	0.901	0.906	0.895	0.892	0.904
Smoking_only	SJLIFE	Meningioma (81)	Radiation	0.189	0.154	0.23	0.28	0.175
Smoking_only	SJLIFE	Meningioma (81)	Chemo	0.341	0.345	0.336	0.506	0.317
Smoking_only	SJLIFE	Meningioma (81)	All treatments	0.474	0.454	0.498	0.68	0.444
Smoking_only	SJLIFE	Meningioma (81)	PRS	-0.114	-0.108	-0.121	-0.112	-0.114
Smoking_only	SJLIFE	Meningioma (81)	Lifestyle	-0.088	-0.08	-0.098	-0.068	-0.091
Smoking_only	SJLIFE	Meningioma (81)	Combined	0.362	0.345	0.383	0.612	0.325
Smoking_only	CCSS	Any SN (1404)	Radiation	0.37	0.364	0.378	0.298	0.384
Smoking_only	CCSS	Any SN (1404)	Chemo	0.022	0.021	0.025	0.044	0.018
Smoking_only	CCSS	Any SN (1404)	All treatments	0.385	0.378	0.395	0.329	0.396
Smoking_only	CCSS	Any SN (1404)	PRS	0.062	0.062	0.063	0.062	0.062
Smoking_only	CCSS	Any SN (1404)	Lifestyle	-0.012	-0.011	-0.014	-0.016	-0.011
Smoking_only	CCSS	Any SN (1404)	Combined	0.415	0.409	0.424	0.36	0.426
Smoking_only	CCSS	SMN (662)	Radiation	0.276	0.277	0.275	0.186	0.295
Smoking_only	CCSS	SMN (662)	Chemo	0.031	0.031	0.033	0.068	0.024
Smoking_only	CCSS	SMN (662)	All treatments	0.303	0.302	0.305	0.242	0.316
Smoking_only	CCSS	SMN (662)	PRS	0.065	0.065	0.066	0.064	0.066
Smoking_only	CCSS	SMN (662)	Lifestyle	-0.015	-0.014	-0.017	-0.02	-0.014
Smoking_only	CCSS	SMN (662)	Combined	0.338	0.337	0.338	0.276	0.351
Smoking_only	CCSS	NMSC (707)	Radiation	0.38	0.377	0.385	0.312	0.389
Smoking_only	CCSS	NMSC (707)	Chemo	NA	NA	NA	NA	NA
Smoking_only	CCSS	NMSC (707)	All treatments	0.38	0.377	0.385	0.312	0.389
Smoking_only	CCSS	NMSC (707)	PRS	0.29	0.293	0.287	0.292	0.29
Smoking_only	CCSS	NMSC (707)	Lifestyle	-0.021	-0.018	-0.026	-0.031	-0.02
Smoking_only	CCSS	NMSC (707)	Combined	0.551	0.552	0.551	0.499	0.558
Smoking_only	CCSS	Breast cancer (278)	Radiation	0.449	0.451	0.447	0.355	0.455
Smoking_only	CCSS	Breast cancer (278)	Chemo	0.129	0.122	0.136	0.159	0.127
Smoking_only	CCSS	Breast cancer (278)	All treatments	0.535	0.534	0.537	0.461	0.54
Smoking_only	CCSS	Breast cancer (278)	PRS	0.321	0.319	0.323	0.322	0.32
Smoking_only	CCSS	Breast cancer (278)	Lifestyle	-0.021	-0.02	-0.023	-0.03	-0.021
Smoking_only	CCSS	Breast cancer (278)	Combined	0.678	0.678	0.679	0.622	0.682
Smoking_only	CCSS	Thyroid cancer (135)	Radiation	0.439	0.436	0.443	0.328	0.477
Smoking_only	CCSS	Thyroid cancer (135)	Chemo	0.059	0.057	0.063	0.117	0.038
Smoking_only	CCSS	Thyroid cancer (135)	All treatments	0.474	0.471	0.48	0.402	0.499
Smoking_only	CCSS	Thyroid cancer (135)	PRS	0.361	0.366	0.353	0.351	0.365
Smoking_only	CCSS	Thyroid cancer (135)	Lifestyle	-0.002	-0.001	-0.002	-0.002	-0.001
Smoking_only	CCSS	Thyroid cancer (135)	Combined	0.663	0.661	0.666	0.613	0.68
Smoking_only	CCSS	Meningioma (215)	Radiation	0.303	0.272	0.349	0.325	0.298
Smoking_only	CCSS	Meningioma (215)	Chemo	0.082	0.074	0.093	0.17	0.06
Smoking_only	CCSS	Meningioma (215)	All treatments	0.371	0.338	0.419	0.451	0.351
Smoking_only	CCSS	Meningioma (215)	PRS	0.015	0.015	0.016	0.02	0.014
Smoking_only	CCSS	Meningioma (215)	Lifestyle	0.001	0.001	0.001	0.001	0.001
Smoking_only	CCSS	Meningioma (215)	Combined	0.381	0.349	0.429	0.462	0.361
Smoking_only	CCSS	Sarcoma (42)	Radiation	NA	NA	NA	NA	NA
Smoking_only	CCSS	Sarcoma (42)	Chemo	0.239	0.227	0.257	0.256	0.235
Smoking_only	CCSS	Sarcoma (42)	All treatments	0.239	0.227	0.257	0.256	0.235
Smoking_only	CCSS	Sarcoma (42)	PRS	0.052	0.054	0.049	0.047	0.054
Smoking_only	CCSS	Sarcoma (42)	Lifestyle	-0.009	-0.009	-0.011	-0.011	-0.009
Smoking_only	CCSS	Sarcoma (42)	Combined	0.271	0.262	0.285	0.282	0.268
Drinking_only	SJLIFE	Any SN (303)	Radiation	0.438	0.429	0.449	0.346	0.458
Drinking_only	SJLIFE	Any SN (303)	Chemo	0.036	0.036	0.036	0.064	0.029
Drinking_only	SJLIFE	Any SN (303)	All treatments	0.456	0.449	0.466	0.384	0.472
Drinking_only	SJLIFE	Any SN (303)	PRS	0.179	0.177	0.182	0.181	0.179
Drinking_only	SJLIFE	Any SN (303)	Lifestyle	-0.024	-0.02	-0.029	-0.023	-0.024
Drinking_only	SJLIFE	Any SN (303)	Combined	0.543	0.538	0.55	0.482	0.557
Drinking_only	SJLIFE	SMN (234)	Radiation	0.378	0.38	0.377	0.266	0.405
Drinking_only	SJLIFE	SMN (234)	Chemo	0.031	0.03	0.032	0.05	0.026
Drinking_only	SJLIFE	SMN (234)	All treatments	0.399	0.401	0.396	0.303	0.421
Drinking_only	SJLIFE	SMN (234)	PRS	0.147	0.145	0.15	0.145	0.148
Drinking_only	SJLIFE	SMN (234)	Lifestyle	-0.011	-0.009	-0.013	-0.01	-0.011
Drinking_only	SJLIFE	SMN (234)	Combined	0.482	0.484	0.479	0.397	0.502
Drinking_only	SJLIFE	NMSC (118)	Radiation	0.263	0.26	0.267	0.235	0.266
Drinking_only	SJLIFE	NMSC (118)	Chemo	NA	NA	NA	NA	NA
Drinking_only	SJLIFE	NMSC (118)	All treatments	0.263	0.26	0.267	0.235	0.266
Drinking_only	SJLIFE	NMSC (118)	PRS	0.44	0.426	0.458	0.446	0.439
Drinking_only	SJLIFE	NMSC (118)	Lifestyle	-0.151	-0.111	-0.201	-0.103	-0.156
Drinking_only	SJLIFE	NMSC (118)	Combined	0.528	0.533	0.522	0.527	0.528
Drinking_only	SJLIFE	Breast cancer (53)	Radiation	0.493	0.517	0.47	0.423	0.498
Drinking_only	SJLIFE	Breast cancer (53)	Chemo	0.057	0.058	0.055	0.048	0.057
Drinking_only	SJLIFE	Breast cancer (53)	All treatments	0.524	0.547	0.502	0.455	0.528
Drinking_only	SJLIFE	Breast cancer (53)	PRS	0.163	0.158	0.168	0.171	0.163
Drinking_only	SJLIFE	Breast cancer (53)	Lifestyle	-0.055	-0.049	-0.061	-0.066	-0.054
Drinking_only	SJLIFE	Breast cancer (53)	Combined	0.581	0.602	0.562	0.518	0.585
Drinking_only	SJLIFE	Thyroid cancer (43)	Radiation	0.608	0.629	0.584	0.548	0.625
Drinking_only	SJLIFE	Thyroid cancer (43)	Chemo	0.244	0.245	0.243	0.36	0.211
Drinking_only	SJLIFE	Thyroid cancer (43)	All treatments	0.725	0.74	0.708	0.703	0.731
Drinking_only	SJLIFE	Thyroid cancer (43)	PRS	0.568	0.567	0.57	0.576	0.566
Drinking_only	SJLIFE	Thyroid cancer (43)	Lifestyle	-0.063	-0.057	-0.071	-0.067	-0.062
Drinking_only	SJLIFE	Thyroid cancer (43)	Combined	0.873	0.88	0.866	0.861	0.877
Drinking_only	SJLIFE	Meningioma (81)	Radiation	0.187	0.152	0.228	0.278	0.174
Drinking_only	SJLIFE	Meningioma (81)	Chemo	0.339	0.342	0.335	0.497	0.316
Drinking_only	SJLIFE	Meningioma (81)	All treatments	0.469	0.449	0.494	0.671	0.44
Drinking_only	SJLIFE	Meningioma (81)	PRS	-0.118	-0.113	-0.125	-0.116	-0.119
Drinking_only	SJLIFE	Meningioma (81)	Lifestyle	-0.168	-0.134	-0.208	-0.139	-0.172
Drinking_only	SJLIFE	Meningioma (81)	Combined	0.306	0.303	0.309	0.576	0.266
Drinking_only	CCSS	Any SN (1404)	Radiation	0.372	0.366	0.38	0.3	0.386
Drinking_only	CCSS	Any SN (1404)	Chemo	0.022	0.02	0.025	0.044	0.018
Drinking_only	CCSS	Any SN (1404)	All treatments	0.387	0.38	0.397	0.33	0.398
Drinking_only	CCSS	Any SN (1404)	PRS	0.062	0.062	0.062	0.061	0.062
Drinking_only	CCSS	Any SN (1404)	Lifestyle	-0.018	-0.018	-0.017	-0.011	-0.019
Drinking_only	CCSS	Any SN (1404)	Combined	0.415	0.408	0.425	0.365	0.424
Drinking_only	CCSS	SMN (662)	Radiation	0.278	0.279	0.277	0.187	0.297
Drinking_only	CCSS	SMN (662)	Chemo	0.032	0.031	0.033	0.068	0.024
Drinking_only	CCSS	SMN (662)	All treatments	0.305	0.304	0.306	0.243	0.318
Drinking_only	CCSS	SMN (662)	PRS	0.065	0.065	0.065	0.064	0.065
Drinking_only	CCSS	SMN (662)	Lifestyle	0.001	0.001	0.001	0.001	0.001
Drinking_only	CCSS	SMN (662)	Combined	0.35	0.349	0.352	0.292	0.362
Drinking_only	CCSS	NMSC (707)	Radiation	0.382	0.378	0.387	0.313	0.39
Drinking_only	CCSS	NMSC (707)	Chemo	NA	NA	NA	NA	NA
Drinking_only	CCSS	NMSC (707)	All treatments	0.382	0.378	0.387	0.313	0.39
Drinking_only	CCSS	NMSC (707)	PRS	0.292	0.294	0.289	0.293	0.292
Drinking_only	CCSS	NMSC (707)	Lifestyle	-0.001	-0.001	-0.001	-0.001	-0.001
Drinking_only	CCSS	NMSC (707)	Combined	0.563	0.561	0.564	0.517	0.568
Drinking_only	CCSS	Breast cancer (278)	Radiation	0.45	0.452	0.447	0.355	0.456
Drinking_only	CCSS	Breast cancer (278)	Chemo	0.128	0.121	0.135	0.156	0.126
Drinking_only	CCSS	Breast cancer (278)	All treatments	0.535	0.534	0.537	0.459	0.54
Drinking_only	CCSS	Breast cancer (278)	PRS	0.321	0.32	0.323	0.322	0.321
Drinking_only	CCSS	Breast cancer (278)	Lifestyle	-0.02	-0.02	-0.019	-0.015	-0.02
Drinking_only	CCSS	Breast cancer (278)	Combined	0.68	0.678	0.681	0.628	0.683
Drinking_only	CCSS	Thyroid cancer (135)	Radiation	0.439	0.436	0.444	0.329	0.477
Drinking_only	CCSS	Thyroid cancer (135)	Chemo	0.059	0.057	0.063	0.117	0.039
Drinking_only	CCSS	Thyroid cancer (135)	All treatments	0.474	0.471	0.48	0.403	0.499
Drinking_only	CCSS	Thyroid cancer (135)	PRS	0.362	0.366	0.354	0.352	0.365
Drinking_only	CCSS	Thyroid cancer (135)	Lifestyle	-0.024	-0.025	-0.022	-0.015	-0.027
Drinking_only	CCSS	Thyroid cancer (135)	Combined	0.656	0.655	0.66	0.609	0.673
Drinking_only	CCSS	Meningioma (215)	Radiation	0.303	0.272	0.349	0.326	0.298
Drinking_only	CCSS	Meningioma (215)	Chemo	0.082	0.074	0.093	0.17	0.06
Drinking_only	CCSS	Meningioma (215)	All treatments	0.371	0.339	0.419	0.451	0.351
Drinking_only	CCSS	Meningioma (215)	PRS	0.016	0.015	0.016	0.021	0.014
Drinking_only	CCSS	Meningioma (215)	Lifestyle	-0.007	-0.008	-0.007	-0.005	-0.008
Drinking_only	CCSS	Meningioma (215)	Combined	0.377	0.344	0.425	0.46	0.356
Drinking_only	CCSS	Sarcoma (42)	Radiation	NA	NA	NA	NA	NA
Drinking_only	CCSS	Sarcoma (42)	Chemo	0.242	0.23	0.262	0.265	0.237
Drinking_only	CCSS	Sarcoma (42)	All treatments	0.242	0.23	0.262	0.265	0.237
Drinking_only	CCSS	Sarcoma (42)	PRS	0.06	0.063	0.055	0.053	0.061
Drinking_only	CCSS	Sarcoma (42)	Lifestyle	0.111	0.115	0.104	0.075	0.12
Drinking_only	CCSS	Sarcoma (42)	Combined	0.37	0.364	0.379	0.355	0.373
Physical_activity_only	SJLIFE	Any SN (303)	Radiation	0.439	0.43	0.45	0.347	0.459
Physical_activity_only	SJLIFE	Any SN (303)	Chemo	0.036	0.036	0.036	0.065	0.03
Physical_activity_only	SJLIFE	Any SN (303)	All treatments	0.457	0.45	0.467	0.386	0.473
Physical_activity_only	SJLIFE	Any SN (303)	PRS	0.18	0.177	0.183	0.182	0.179
Physical_activity_only	SJLIFE	Any SN (303)	Lifestyle	0.027	0.03	0.022	0.024	0.027
Physical_activity_only	SJLIFE	Any SN (303)	Combined	0.568	0.562	0.574	0.508	0.581
Physical_activity_only	SJLIFE	SMN (234)	Radiation	0.378	0.38	0.376	0.266	0.404
Physical_activity_only	SJLIFE	SMN (234)	Chemo	0.03	0.029	0.031	0.049	0.025
Physical_activity_only	SJLIFE	SMN (234)	All treatments	0.398	0.4	0.395	0.302	0.42
Physical_activity_only	SJLIFE	SMN (234)	PRS	0.148	0.146	0.151	0.146	0.149
Physical_activity_only	SJLIFE	SMN (234)	Lifestyle	0.063	0.072	0.052	0.055	0.065
Physical_activity_only	SJLIFE	SMN (234)	Combined	0.52	0.526	0.512	0.435	0.539
Physical_activity_only	SJLIFE	NMSC (118)	Radiation	0.277	0.28	0.273	0.243	0.28
Physical_activity_only	SJLIFE	NMSC (118)	Chemo	NA	NA	NA	NA	NA
Physical_activity_only	SJLIFE	NMSC (118)	All treatments	0.277	0.28	0.273	0.243	0.28
Physical_activity_only	SJLIFE	NMSC (118)	PRS	0.418	0.404	0.437	0.42	0.418
Physical_activity_only	SJLIFE	NMSC (118)	Lifestyle	-0.064	-0.067	-0.059	-0.051	-0.065
Physical_activity_only	SJLIFE	NMSC (118)	Combined	0.556	0.547	0.568	0.54	0.558
Physical_activity_only	SJLIFE	Breast cancer (53)	Radiation	0.493	0.516	0.473	0.427	0.498
Physical_activity_only	SJLIFE	Breast cancer (53)	Chemo	0.058	0.059	0.056	0.05	0.058
Physical_activity_only	SJLIFE	Breast cancer (53)	All treatments	0.525	0.547	0.505	0.459	0.529
Physical_activity_only	SJLIFE	Breast cancer (53)	PRS	0.163	0.156	0.17	0.174	0.162
Physical_activity_only	SJLIFE	Breast cancer (53)	Lifestyle	-0.114	-0.139	-0.092	-0.083	-0.116
Physical_activity_only	SJLIFE	Breast cancer (53)	Combined	0.559	0.569	0.551	0.512	0.563
Physical_activity_only	SJLIFE	Thyroid cancer (43)	Radiation	0.61	0.631	0.586	0.552	0.627
Physical_activity_only	SJLIFE	Thyroid cancer (43)	Chemo	0.254	0.256	0.251	0.378	0.218
Physical_activity_only	SJLIFE	Thyroid cancer (43)	All treatments	0.732	0.747	0.716	0.716	0.737
Physical_activity_only	SJLIFE	Thyroid cancer (43)	PRS	0.572	0.57	0.575	0.582	0.57
Physical_activity_only	SJLIFE	Thyroid cancer (43)	Lifestyle	-0.158	-0.191	-0.121	-0.123	-0.168
Physical_activity_only	SJLIFE	Thyroid cancer (43)	Combined	0.868	0.871	0.865	0.861	0.87
Physical_activity_only	SJLIFE	Meningioma (81)	Radiation	0.19	0.155	0.233	0.275	0.178
Physical_activity_only	SJLIFE	Meningioma (81)	Chemo	0.341	0.343	0.338	0.506	0.317
Physical_activity_only	SJLIFE	Meningioma (81)	All treatments	0.474	0.451	0.501	0.678	0.444
Physical_activity_only	SJLIFE	Meningioma (81)	PRS	-0.111	-0.104	-0.119	-0.108	-0.111
Physical_activity_only	SJLIFE	Meningioma (81)	Lifestyle	-0.106	-0.123	-0.087	-0.089	-0.109
Physical_activity_only	SJLIFE	Meningioma (81)	Combined	0.357	0.324	0.396	0.612	0.32
Physical_activity_only	CCSS	Any SN (1404)	Radiation	0.371	0.365	0.38	0.299	0.385
Physical_activity_only	CCSS	Any SN (1404)	Chemo	0.022	0.021	0.025	0.044	0.018
Physical_activity_only	CCSS	Any SN (1404)	All treatments	0.386	0.379	0.397	0.33	0.397
Physical_activity_only	CCSS	Any SN (1404)	PRS	0.062	0.061	0.062	0.061	0.062
Physical_activity_only	CCSS	Any SN (1404)	Lifestyle	-0.002	-0.002	-0.001	-0.007	-0.001
Physical_activity_only	CCSS	Any SN (1404)	Combined	0.423	0.416	0.433	0.367	0.434
Physical_activity_only	CCSS	SMN (662)	Radiation	0.278	0.279	0.276	0.187	0.297
Physical_activity_only	CCSS	SMN (662)	Chemo	0.031	0.031	0.033	0.068	0.024
Physical_activity_only	CCSS	SMN (662)	All treatments	0.304	0.304	0.306	0.243	0.317
Physical_activity_only	CCSS	SMN (662)	PRS	0.065	0.065	0.066	0.064	0.066
Physical_activity_only	CCSS	SMN (662)	Lifestyle	-0.006	-0.008	-0.004	-0.027	-0.002
Physical_activity_only	CCSS	SMN (662)	Combined	0.344	0.343	0.348	0.271	0.36
Physical_activity_only	CCSS	NMSC (707)	Radiation	0.382	0.378	0.388	0.313	0.391
Physical_activity_only	CCSS	NMSC (707)	Chemo	NA	NA	NA	NA	NA
Physical_activity_only	CCSS	NMSC (707)	All treatments	0.382	0.378	0.388	0.313	0.391
Physical_activity_only	CCSS	NMSC (707)	PRS	0.291	0.294	0.289	0.293	0.291
Physical_activity_only	CCSS	NMSC (707)	Lifestyle	0.005	0.007	0.004	0.031	0.002
Physical_activity_only	CCSS	NMSC (707)	Combined	0.566	0.565	0.566	0.532	0.57
Physical_activity_only	CCSS	Breast cancer (278)	Radiation	0.449	0.452	0.446	0.359	0.455
Physical_activity_only	CCSS	Breast cancer (278)	Chemo	0.13	0.123	0.137	0.162	0.128
Physical_activity_only	CCSS	Breast cancer (278)	All treatments	0.537	0.536	0.538	0.466	0.541
Physical_activity_only	CCSS	Breast cancer (278)	PRS	0.321	0.319	0.323	0.321	0.321
Physical_activity_only	CCSS	Breast cancer (278)	Lifestyle	-0.006	-0.007	-0.004	-0.055	-0.002
Physical_activity_only	CCSS	Breast cancer (278)	Combined	0.684	0.683	0.686	0.613	0.689
Physical_activity_only	CCSS	Thyroid cancer (135)	Radiation	0.438	0.436	0.442	0.33	0.476
Physical_activity_only	CCSS	Thyroid cancer (135)	Chemo	0.059	0.057	0.062	0.117	0.038
Physical_activity_only	CCSS	Thyroid cancer (135)	All treatments	0.474	0.471	0.478	0.403	0.498
Physical_activity_only	CCSS	Thyroid cancer (135)	PRS	0.361	0.366	0.353	0.351	0.365
Physical_activity_only	CCSS	Thyroid cancer (135)	Lifestyle	-0.01	-0.012	-0.006	-0.028	-0.004
Physical_activity_only	CCSS	Thyroid cancer (135)	Combined	0.658	0.656	0.663	0.601	0.678
Physical_activity_only	CCSS	Meningioma (215)	Radiation	0.303	0.271	0.349	0.325	0.298
Physical_activity_only	CCSS	Meningioma (215)	Chemo	0.082	0.074	0.093	0.17	0.06
Physical_activity_only	CCSS	Meningioma (215)	All treatments	0.371	0.338	0.419	0.451	0.352
Physical_activity_only	CCSS	Meningioma (215)	PRS	0.015	0.015	0.016	0.02	0.014
Physical_activity_only	CCSS	Meningioma (215)	Lifestyle	-0.011	-0.014	-0.007	-0.041	-0.004
Physical_activity_only	CCSS	Meningioma (215)	Combined	0.375	0.34	0.425	0.439	0.359
Physical_activity_only	CCSS	Sarcoma (42)	Radiation	NA	NA	NA	NA	NA
Physical_activity_only	CCSS	Sarcoma (42)	Chemo	0.237	0.226	0.255	0.257	0.232
Physical_activity_only	CCSS	Sarcoma (42)	All treatments	0.237	0.226	0.255	0.257	0.232
Physical_activity_only	CCSS	Sarcoma (42)	PRS	0.054	0.056	0.051	0.049	0.055
Physical_activity_only	CCSS	Sarcoma (42)	Lifestyle	-0.018	-0.022	-0.012	-0.067	-0.006
Physical_activity_only	CCSS	Sarcoma (42)	Combined	0.264	0.252	0.284	0.242	0.27
Obesity_only	SJLIFE	Any SN (303)	Radiation	0.439	0.43	0.451	0.347	0.46
Obesity_only	SJLIFE	Any SN (303)	Chemo	0.036	0.036	0.036	0.065	0.03
Obesity_only	SJLIFE	Any SN (303)	All treatments	0.458	0.45	0.468	0.386	0.474
Obesity_only	SJLIFE	Any SN (303)	PRS	0.179	0.177	0.182	0.181	0.179
Obesity_only	SJLIFE	Any SN (303)	Lifestyle	0.001	0.001	0.001	0.001	0.001
Obesity_only	SJLIFE	Any SN (303)	Combined	0.557	0.55	0.566	0.497	0.57
Obesity_only	SJLIFE	SMN (234)	Radiation	0.382	0.382	0.382	0.269	0.408
Obesity_only	SJLIFE	SMN (234)	Chemo	0.03	0.029	0.031	0.049	0.025
Obesity_only	SJLIFE	SMN (234)	All treatments	0.401	0.402	0.4	0.305	0.424
Obesity_only	SJLIFE	SMN (234)	PRS	0.147	0.144	0.15	0.144	0.148
Obesity_only	SJLIFE	SMN (234)	Lifestyle	-0.047	-0.047	-0.046	-0.033	-0.05
Obesity_only	SJLIFE	SMN (234)	Combined	0.466	0.467	0.464	0.385	0.484
Obesity_only	SJLIFE	NMSC (118)	Radiation	0.278	0.287	0.267	0.248	0.282
Obesity_only	SJLIFE	NMSC (118)	Chemo	NA	NA	NA	NA	NA
Obesity_only	SJLIFE	NMSC (118)	All treatments	0.278	0.287	0.267	0.248	0.282
Obesity_only	SJLIFE	NMSC (118)	PRS	0.412	0.397	0.432	0.414	0.412
Obesity_only	SJLIFE	NMSC (118)	Lifestyle	-0.076	-0.077	-0.074	-0.055	-0.078
Obesity_only	SJLIFE	NMSC (118)	Combined	0.549	0.542	0.557	0.535	0.55
Obesity_only	SJLIFE	Breast cancer (53)	Radiation	0.49	0.514	0.466	0.421	0.494
Obesity_only	SJLIFE	Breast cancer (53)	Chemo	0.057	0.058	0.056	0.05	0.058
Obesity_only	SJLIFE	Breast cancer (53)	All treatments	0.521	0.544	0.498	0.453	0.525
Obesity_only	SJLIFE	Breast cancer (53)	PRS	0.159	0.153	0.165	0.168	0.158
Obesity_only	SJLIFE	Breast cancer (53)	Lifestyle	-0.072	-0.066	-0.077	-0.058	-0.073
Obesity_only	SJLIFE	Breast cancer (53)	Combined	0.566	0.587	0.547	0.515	0.57
Obesity_only	SJLIFE	Thyroid cancer (43)	Radiation	0.609	0.63	0.584	0.55	0.626
Obesity_only	SJLIFE	Thyroid cancer (43)	Chemo	0.248	0.249	0.248	0.366	0.214
Obesity_only	SJLIFE	Thyroid cancer (43)	All treatments	0.728	0.743	0.711	0.709	0.733
Obesity_only	SJLIFE	Thyroid cancer (43)	PRS	0.571	0.569	0.573	0.58	0.568
Obesity_only	SJLIFE	Thyroid cancer (43)	Lifestyle	-0.054	-0.053	-0.055	-0.039	-0.058
Obesity_only	SJLIFE	Thyroid cancer (43)	Combined	0.877	0.883	0.87	0.867	0.88
Obesity_only	SJLIFE	Meningioma (81)	Radiation	0.188	0.151	0.233	0.274	0.176
Obesity_only	SJLIFE	Meningioma (81)	Chemo	0.337	0.342	0.331	0.496	0.314
Obesity_only	SJLIFE	Meningioma (81)	All treatments	0.469	0.448	0.495	0.67	0.44
Obesity_only	SJLIFE	Meningioma (81)	PRS	-0.12	-0.113	-0.128	-0.117	-0.12
Obesity_only	SJLIFE	Meningioma (81)	Lifestyle	0.138	0.132	0.145	0.101	0.144
Obesity_only	SJLIFE	Meningioma (81)	Combined	0.488	0.467	0.512	0.666	0.462
Obesity_only	CCSS	Any SN (1404)	Radiation	0.371	0.365	0.38	0.299	0.385
Obesity_only	CCSS	Any SN (1404)	Chemo	0.022	0.021	0.025	0.044	0.018
Obesity_only	CCSS	Any SN (1404)	All treatments	0.386	0.379	0.396	0.33	0.397
Obesity_only	CCSS	Any SN (1404)	PRS	0.062	0.061	0.062	0.061	0.062
Obesity_only	CCSS	Any SN (1404)	Lifestyle	0.004	0.004	0.004	0.008	0.003
Obesity_only	CCSS	Any SN (1404)	Combined	0.426	0.419	0.436	0.376	0.436
Obesity_only	CCSS	SMN (662)	Radiation	0.277	0.278	0.276	0.185	0.297
Obesity_only	CCSS	SMN (662)	Chemo	0.032	0.031	0.033	0.068	0.024
Obesity_only	CCSS	SMN (662)	All treatments	0.304	0.303	0.306	0.242	0.317
Obesity_only	CCSS	SMN (662)	PRS	0.065	0.065	0.065	0.064	0.065
Obesity_only	CCSS	SMN (662)	Lifestyle	0.007	0.007	0.007	0.015	0.005
Obesity_only	CCSS	SMN (662)	Combined	0.354	0.353	0.355	0.301	0.365
Obesity_only	CCSS	NMSC (707)	Radiation	0.382	0.378	0.388	0.314	0.391
Obesity_only	CCSS	NMSC (707)	Chemo	NA	NA	NA	NA	NA
Obesity_only	CCSS	NMSC (707)	All treatments	0.382	0.378	0.388	0.314	0.391
Obesity_only	CCSS	NMSC (707)	PRS	0.292	0.294	0.289	0.293	0.291
Obesity_only	CCSS	NMSC (707)	Lifestyle	-0.005	-0.005	-0.005	-0.011	-0.004
Obesity_only	CCSS	NMSC (707)	Combined	0.561	0.56	0.563	0.512	0.567
Obesity_only	CCSS	Breast cancer (278)	Radiation	0.449	0.452	0.447	0.355	0.455
Obesity_only	CCSS	Breast cancer (278)	Chemo	0.13	0.123	0.137	0.16	0.128
Obesity_only	CCSS	Breast cancer (278)	All treatments	0.537	0.535	0.538	0.461	0.541
Obesity_only	CCSS	Breast cancer (278)	PRS	0.321	0.319	0.323	0.322	0.321
Obesity_only	CCSS	Breast cancer (278)	Lifestyle	0.001	0.001	0.001	0.002	0.001
Obesity_only	CCSS	Breast cancer (278)	Combined	0.687	0.686	0.688	0.635	0.69
Obesity_only	CCSS	Thyroid cancer (135)	Radiation	0.439	0.436	0.444	0.328	0.477
Obesity_only	CCSS	Thyroid cancer (135)	Chemo	0.059	0.057	0.062	0.117	0.038
Obesity_only	CCSS	Thyroid cancer (135)	All treatments	0.474	0.471	0.48	0.402	0.499
Obesity_only	CCSS	Thyroid cancer (135)	PRS	0.362	0.367	0.354	0.352	0.365
Obesity_only	CCSS	Thyroid cancer (135)	Lifestyle	0.008	0.008	0.008	0.015	0.005
Obesity_only	CCSS	Thyroid cancer (135)	Combined	0.667	0.665	0.67	0.62	0.683
Obesity_only	CCSS	Meningioma (215)	Radiation	0.303	0.271	0.348	0.326	0.297
Obesity_only	CCSS	Meningioma (215)	Chemo	0.082	0.074	0.093	0.17	0.06
Obesity_only	CCSS	Meningioma (215)	All treatments	0.371	0.339	0.418	0.452	0.351
Obesity_only	CCSS	Meningioma (215)	PRS	0.016	0.016	0.017	0.021	0.015
Obesity_only	CCSS	Meningioma (215)	Lifestyle	0.005	0.006	0.005	0.01	0.004
Obesity_only	CCSS	Meningioma (215)	Combined	0.384	0.352	0.431	0.468	0.363
Obesity_only	CCSS	Sarcoma (42)	Radiation	NA	NA	NA	NA	NA
Obesity_only	CCSS	Sarcoma (42)	Chemo	0.239	0.228	0.258	0.255	0.235
Obesity_only	CCSS	Sarcoma (42)	All treatments	0.239	0.228	0.258	0.255	0.235
Obesity_only	CCSS	Sarcoma (42)	PRS	0.052	0.054	0.049	0.048	0.054
Obesity_only	CCSS	Sarcoma (42)	Lifestyle	0.002	0.002	0.003	0.005	0.002
Obesity_only	CCSS	Sarcoma (42)	Combined	0.28	0.271	0.296	0.294	0.277", header = T, sep = "\t", check.names = FALSE)

# df <- mydf[grepl("Any SN", mydf$SN_types),]





# Melt the dataframe for plotting
df_melted <- reshape::melt(mydf, id.vars = c("Lifestyles", "Cohort", "SN_types", "Variables"))

## Values to percentage
df_melted$value <- df_melted$value *100 
# df_melted$value [df_melted$value < 0] <- NA ## Non-significant

df_melted$bar_labels <- ifelse(df_melted$value > 0, paste0(round(df_melted$value, 2), "%"), "NS") 

df_melted$n <- as.numeric(gsub("\\D", "", df_melted$SN_types))
df_melted$SN_types <- sub("\\s*\\(.*", "", df_melted$SN_types)



SNs <- unique(df_melted$SN_types)
for (i in 1:length(SNs)){

print(paste0("Doing ", SNs[i]))
df_melted.tmp <-   df_melted[grepl(SNs[i], df_melted$SN_types),]
  

# Define the color palette
color_palette <- c("#1F77B4", "#FF7F0E", "#2CA02C", "#D62728", "#9467BD")

# Create the plot
P <- ggplot(df_melted.tmp, aes(x = Variables, y = value, fill = factor(variable))) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.2) +
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
  facet_grid(Lifestyles~ Cohort, scales = "free_x", space = "free_x", switch = "x") +
  scale_y_continuous(labels = percent_format(scale = 1),    # Format Y-axis labels as percentages
                     limits = c(0, 100), breaks = seq(0, 100, 30)) +
  guides(fill = guide_legend(nrow = 1)) +  # Limits the legend to one row
  geom_text(aes(label = bar_labels),   # Add percentage symbol or "NS" to the label
            position = position_dodge(width = 0.8), 
            hjust = -0.4,                                 # Adjust vertical position of the labels above the bars
            size = 3,
            angle = 90) 

plotname <- gsub(" ","_", paste0("V17_plots_without_diet_", SNs[i], ".tiff"))
ggsave(filename = plotname, plot = P, width = 10, height = 10, dpi = 600)

}






