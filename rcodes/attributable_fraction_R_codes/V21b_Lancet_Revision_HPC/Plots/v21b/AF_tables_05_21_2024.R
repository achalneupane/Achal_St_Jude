# Load required libraries
library(ggplot2)
library(ggthemes)  # For additional themes
library(tidyverse)

# ## V18 b (without lifestyle)
# data <- read.table(text="Cohort	SN_types	Variables	Overall	Female	Male	Age.lt.35	Age.ge.35
# SJLIFE	SNs (605)	Radiation	0.425	0.417	0.435	0.4	0.447
# SJLIFE	SNs (605)	Chemotherapy	0.08	0.079	0.08	0.109	0.054
# SJLIFE	SNs (605)	All_treatments	0.469	0.462	0.478	0.463	0.474
# SJLIFE	SNs (605)	PRS	0.125	0.123	0.129	0.126	0.125
# SJLIFE	SNs (605)	Lifestyle	-	-	-	-	-
# SJLIFE	SNs (605)	Combined	0.536	0.53	0.545	0.531	0.541
# SJLIFE	SMNs (463)	Radiation	0.37	0.369	0.371	0.339	0.394
# SJLIFE	SMNs (463)	Chemotherapy	0.032	0.032	0.032	0.045	0.021
# SJLIFE	SMNs (463)	All_treatments	0.39	0.39	0.39	0.368	0.407
# SJLIFE	SMNs (463)	PRS	0.117	0.115	0.12	0.117	0.117
# SJLIFE	SMNs (463)	Lifestyle	-	-	-	-	-
# SJLIFE	SMNs (463)	Combined	0.462	0.461	0.464	0.442	0.478
# SJLIFE	NMSC (251)	Radiation	0.436	0.427	0.448	0.408	0.452
# SJLIFE	NMSC (251)	Chemotherapy	-	-	-	-	-
# SJLIFE	NMSC (251)	All_treatments	0.436	0.427	0.448	0.408	0.452
# SJLIFE	NMSC (251)	PRS	0.28	0.266	0.296	0.278	0.281
# SJLIFE	NMSC (251)	Lifestyle	-	-	-	-	-
# SJLIFE	NMSC (251)	Combined	0.593	0.581	0.607	0.571	0.605
# SJLIFE	Breast cancer (76)	Radiation	0.508	0.508	-	0.491	0.513
# SJLIFE	Breast cancer (76)	Chemotherapy	0.189	0.189	-	0.222	0.18
# SJLIFE	Breast cancer (76)	All_treatments	0.613	0.613	-	0.604	0.615
# SJLIFE	Breast cancer (76)	PRS	0.192	0.192	-	0.197	0.191
# SJLIFE	Breast cancer (76)	Lifestyle	-	-	-	-	-
# SJLIFE	Breast cancer (76)	Combined	0.687	0.687	-	0.683	0.689
# SJLIFE	Thyroid cancer (87)	Radiation	0.62	0.623	0.617	0.594	0.651
# SJLIFE	Thyroid cancer (87)	Chemotherapy	0.235	0.236	0.233	0.294	0.165
# SJLIFE	Thyroid cancer (87)	All_treatments	0.728	0.729	0.726	0.727	0.73
# SJLIFE	Thyroid cancer (87)	PRS	0.517	0.519	0.514	0.519	0.515
# SJLIFE	Thyroid cancer (87)	Lifestyle	-	-	-	-	-
# SJLIFE	Thyroid cancer (87)	Combined	0.868	0.868	0.867	0.867	0.869
# SJLIFE	Meningioma (149)	Radiation	0.901	0.899	0.904	0.895	0.907
# SJLIFE	Meningioma (149)	Chemotherapy	0.167	0.182	0.149	0.219	0.128
# SJLIFE	Meningioma (149)	All_treatments	0.914	0.912	0.916	0.914	0.914
# SJLIFE	Meningioma (149)	PRS	-	-	-	-	-
# SJLIFE	Meningioma (149)	Lifestyle	-	-	-	-	-
# SJLIFE	Meningioma (149)	Combined	0.904	0.903	0.906	0.904	0.904
# SJLIFE	Sarcoma (33)	Radiation	-	-	-	-	-
# SJLIFE	Sarcoma (33)	Chemotherapy	0.346	0.349	0.342	0.329	0.371
# SJLIFE	Sarcoma (33)	All_treatments	0.346	0.349	0.342	0.329	0.371
# SJLIFE	Sarcoma (33)	PRS	-	-	-	-	-
# SJLIFE	Sarcoma (33)	Lifestyle	-	-	-	-	-
# SJLIFE	Sarcoma (33)	Combined	0.298	0.301	0.293	0.28	0.324
# CCSS	SNs (1611)	Radiation	0.387	0.38	0.398	0.372	0.398
# CCSS	SNs (1611)	Chemotherapy	0.031	0.03	0.033	0.047	0.019
# CCSS	SNs (1611)	All_treatments	0.407	0.4	0.418	0.403	0.41
# CCSS	SNs (1611)	PRS	0.048	0.047	0.049	0.047	0.048
# CCSS	SNs (1611)	Lifestyle	-	-	-	-	-
# CCSS	SNs (1611)	Combined	0.435	0.428	0.446	0.431	0.438
# CCSS	SMNs (762)	Radiation	0.255	0.254	0.257	0.213	0.283
# CCSS	SMNs (762)	Chemotherapy	0.042	0.042	0.042	0.07	0.024
# CCSS	SMNs (762)	All_treatments	0.29	0.289	0.293	0.271	0.303
# CCSS	SMNs (762)	PRS	0.048	0.047	0.048	0.047	0.048
# CCSS	SMNs (762)	Lifestyle	-	-	-	-	-
# CCSS	SMNs (762)	Combined	0.324	0.323	0.326	0.305	0.336
# CCSS	NMSC (774)	Radiation	0.391	0.387	0.397	0.386	0.394
# CCSS	NMSC (774)	Chemotherapy	-	-	-	-	-
# CCSS	NMSC (774)	All_treatments	0.391	0.387	0.397	0.386	0.394
# CCSS	NMSC (774)	PRS	0.306	0.304	0.309	0.308	0.305
# CCSS	NMSC (774)	Lifestyle	-	-	-	-	-
# CCSS	NMSC (774)	Combined	0.579	0.575	0.583	0.577	0.58
# CCSS	Breast cancer (290)	Radiation	0.475	0.475	-	0.454	0.479
# CCSS	Breast cancer (290)	Chemotherapy	0.186	0.186	-	0.223	0.179
# CCSS	Breast cancer (290)	All_treatments	0.593	0.593	-	0.589	0.594
# CCSS	Breast cancer (290)	PRS	0.365	0.365	-	0.368	0.365
# CCSS	Breast cancer (290)	Lifestyle	-	-	-	-	-
# CCSS	Breast cancer (290)	Combined	0.743	0.743	-	0.741	0.744
# CCSS	Thyroid cancer (163)	Radiation	0.442	0.444	0.438	0.413	0.488
# CCSS	Thyroid cancer (163)	Chemotherapy	0.056	0.054	0.059	0.073	0.03
# CCSS	Thyroid cancer (163)	All_treatments	0.475	0.476	0.472	0.456	0.503
# CCSS	Thyroid cancer (163)	PRS	0.358	0.362	0.351	0.355	0.362
# CCSS	Thyroid cancer (163)	Lifestyle	-	-	-	-	-
# CCSS	Thyroid cancer (163)	Combined	0.663	0.664	0.66	0.649	0.683
# CCSS	Meningioma (256)	Radiation	0.864	0.855	0.877	0.87	0.858
# CCSS	Meningioma (256)	Chemotherapy	0.082	0.073	0.095	0.11	0.052
# CCSS	Meningioma (256)	All_treatments	0.877	0.868	0.889	0.887	0.865
# CCSS	Meningioma (256)	PRS	0.013	0.012	0.014	0.016	0.009
# CCSS	Meningioma (256)	Lifestyle	-	-	-	-	-
# CCSS	Meningioma (256)	Combined	0.878	0.869	0.891	0.889	0.867
# CCSS	Sarcoma (61)	Radiation	-	-	-	-	-
# CCSS	Sarcoma (61)	Chemotherapy	0.353	0.338	0.369	0.349	0.358
# CCSS	Sarcoma (61)	All_treatments	0.353	0.338	0.369	0.349	0.358
# CCSS	Sarcoma (61)	PRS	0.041	0.042	0.039	0.04	0.042
# CCSS	Sarcoma (61)	Lifestyle	-	-	-	-	-
# CCSS	Sarcoma (61)	Combined	0.379	0.367	0.394	0.375	0.385", header = T, sep = "\t")
# 
# data$SN_types_new <- gsub("\\([0-9]+\\)", "", data$SN_types)
# data[data == "-"] <- NA
# 
# 
# lifestyle <- "without_lifestyle"
# all.group <- c("Overall", "Female", "Male", "Age.lt.35", "Age.ge.35")
# variables <-  c("Combined", "Radiation", "Chemotherapy", "treatments", "PRS", "Lifestyle")
# 
# 
# # Define the text file path
# data_name <- paste0("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/ANALYSIS/results/plots/v21b/", lifestyle, "_all_table_19b.txt")
# 
# file_conn <- file(data_name, "a")
# 
# for (j in 1:length(all.group)){
# group <-  all.group [j]
# 
# for(i in 1:length(variables)){
# # print(paste0("Doing ", group, "_", AF.type, "_", lifestyle))
#   
# if(lifestyle== "without_lifestyle" && variables[i] == "Lifestyle"){
#   next
# }
#     
# AF.type <- variables[i]
# new_data <- data[grepl(AF.type, data$Variables), c("Cohort", "SN_types", "SN_types_new", "Variables", group)]
# 
# # Filter the data for SJLIFE and CCSS cohorts
# sjlife_data <- new_data %>% filter(Cohort == "SJLIFE")
# ccss_data <- new_data %>% filter(Cohort == "CCSS")
# 
# # Create a new data frame with the desired format
# result_data <- tibble(
#   "SN types (CA in SJLIFE/CA in CCSS)" = paste0(sjlife_data$SN_types_new, " (", as.numeric(gsub("[^0-9.]", "", sjlife_data$SN_types)), "/", as.numeric(gsub("[^0-9.]", "", ccss_data$SN_types)), ")"),
#   "SJLIFE" = as.numeric(sjlife_data[,group]),
#   "CCSS" = as.numeric(ccss_data[,group]))
# 
# # Print the result_data data frame
# table_name <- paste0(group, "_", AF.type, "_", lifestyle,"\tNA\tNA")
# # write.table(table_name, data_name, col.names = T, row.names = F, quote = F, append = T, sep = "\t")
# # write.table(result_data, data_name, col.names = T, row.names = F, quote = F, append = T, sep = "\t")
# # }
# 
# cat(table_name, "\n", file = file_conn)
# write.table(result_data, file = file_conn, col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")
# }
# }
# # Close the text file
# close(file_conn)
# 
# ###########################
# ## Lifestyle without HEI ##
# ###########################
# data <- read.table(text = "Cohort	SN_types	Variables	Overall	Female	Male	Age.lt.35	Age.ge.35
# SJLIFE	SNs (303)	Radiation	0.442	0.432	0.455	0.367	0.472
# SJLIFE	SNs (303)	Chemotherapy	0.015	0.016	0.014	0.027	0.01
# SJLIFE	SNs (303)	All_treatments	0.449	0.441	0.46	0.382	0.476
# SJLIFE	SNs (303)	PRS	0.171	0.169	0.174	0.174	0.17
# SJLIFE	SNs (303)	Lifestyle	-	-	-	-	-
# SJLIFE	SNs (303)	Combined	0.525	0.521	0.53	0.461	0.55
# SJLIFE	SMNs (234)	Radiation	0.378	0.377	0.379	0.289	0.412
# SJLIFE	SMNs (234)	Chemotherapy	0.014	0.014	0.013	0.022	0.01
# SJLIFE	SMNs (234)	All_treatments	0.387	0.387	0.386	0.303	0.418
# SJLIFE	SMNs (234)	PRS	0.133	0.13	0.136	0.132	0.133
# SJLIFE	SMNs (234)	Lifestyle	-	-	-	-	-
# SJLIFE	SMNs (234)	Combined	0.436	0.443	0.427	0.351	0.468
# SJLIFE	NMSCs (119)	Radiation	0.298	0.286	0.312	0.239	0.31
# SJLIFE	NMSCs (119)	Chemotherapy	-	-	-	-	-
# SJLIFE	NMSCs (119)	All_treatments	0.298	0.286	0.312	0.239	0.31
# SJLIFE	NMSCs (119)	PRS	0.311	0.293	0.334	0.308	0.311
# SJLIFE	NMSCs (119)	Lifestyle	-	-	-	-	-
# SJLIFE	NMSCs (119)	Combined	0.416	0.403	0.434	0.363	0.428
# SJLIFE	Breast cancer (51)	Radiation	0.47	0.47	-	0.387	0.483
# SJLIFE	Breast cancer (51)	Chemotherapy	0.072	0.072	-	0.056	0.074
# SJLIFE	Breast cancer (51)	All_treatments	0.511	0.511	-	0.42	0.526
# SJLIFE	Breast cancer (51)	PRS	0.177	0.177	-	0.198	0.174
# SJLIFE	Breast cancer (51)	Lifestyle	-	-	-	-	-
# SJLIFE	Breast cancer (51)	Combined	0.341	0.341	-	0.314	0.345
# SJLIFE	Thyroid cancer (44)	Radiation	0.601	0.62	0.577	0.535	0.63
# SJLIFE	Thyroid cancer (44)	Chemotherapy	0.221	0.222	0.22	0.312	0.181
# SJLIFE	Thyroid cancer (44)	All_treatments	0.71	0.725	0.693	0.686	0.721
# SJLIFE	Thyroid cancer (44)	PRS	0.572	0.571	0.573	0.583	0.567
# SJLIFE	Thyroid cancer (44)	Lifestyle	-	-	-	-	-
# SJLIFE	Thyroid cancer (44)	Combined	0.849	0.854	0.842	0.833	0.856
# SJLIFE	Meningioma (81)	Radiation	0.185	0.151	0.225	0.251	0.162
# SJLIFE	Meningioma (81)	Chemotherapy	0.318	0.323	0.312	0.454	0.273
# SJLIFE	Meningioma (81)	All_treatments	0.456	0.436	0.479	0.623	0.4
# SJLIFE	Meningioma (81)	PRS	-	-	-	-	-
# SJLIFE	Meningioma (81)	Lifestyle	-	-	-	-	-
# SJLIFE	Meningioma (81)	Combined	0.217	0.21	0.226	0.461	0.136
# SJLIFE	Sarcoma (NA)	Radiation	-	-	-	-	-
# SJLIFE	Sarcoma (NA)	Chemotherapy	-	-	-	-	-
# SJLIFE	Sarcoma (NA)	All_treatments	-	-	-	-	-
# SJLIFE	Sarcoma (NA)	PRS	-	-	-	-	-
# SJLIFE	Sarcoma (NA)	Lifestyle	-	-	-	-	-
# SJLIFE	Sarcoma (NA)	Combined	-	-	-	-	-
# CCSS	SNs (1286)	Radiation	0.353	0.346	0.364	0.332	0.364
# CCSS	SNs (1286)	Chemotherapy	0.038	0.036	0.04	0.06	0.027
# CCSS	SNs (1286)	All_treatments	0.379	0.372	0.391	0.374	0.382
# CCSS	SNs (1286)	PRS	0.058	0.057	0.059	0.056	0.058
# CCSS	SNs (1286)	Lifestyle	-	-	-	-	-
# CCSS	SNs (1286)	Combined	0.352	0.345	0.363	0.352	0.352
# CCSS	SMNs (619)	Radiation	0.245	0.241	0.251	0.182	0.271
# CCSS	SMNs (619)	Chemotherapy	0.055	0.054	0.057	0.098	0.038
# CCSS	SMNs (619)	All_treatments	0.293	0.288	0.301	0.266	0.304
# CCSS	SMNs (619)	PRS	0.047	0.046	0.048	0.046	0.047
# CCSS	SMNs (619)	Lifestyle	-	-	-	-	-
# CCSS	SMNs (619)	Combined	0.273	0.27	0.28	0.252	0.282
# CCSS	NMSCs (637)	Radiation	0.385	0.383	0.387	0.393	0.382
# CCSS	NMSCs (637)	Chemotherapy	-	-	-	-	-
# CCSS	NMSCs (637)	All_treatments	0.385	0.383	0.387	0.393	0.382
# CCSS	NMSCs (637)	PRS	0.307	0.304	0.309	0.308	0.306
# CCSS	NMSCs (637)	Lifestyle	-	-	-	-	-
# CCSS	NMSCs (637)	Combined	0.485	0.48	0.49	0.504	0.477
# CCSS	Breast cancer (260)	Radiation	0.444	0.444	-	0.392	0.45
# CCSS	Breast cancer (260)	Chemotherapy	0.193	0.193	-	0.232	0.188
# CCSS	Breast cancer (260)	All_treatments	0.572	0.572	-	0.541	0.575
# CCSS	Breast cancer (260)	PRS	0.348	0.348	-	0.35	0.348
# CCSS	Breast cancer (260)	Lifestyle	-	-	-	-	-
# CCSS	Breast cancer (260)	Combined	0.69	0.69	-	0.673	0.692
# CCSS	Thyroid cancer (123)	Radiation	0.416	0.416	0.417	0.359	0.486
# CCSS	Thyroid cancer (123)	Chemotherapy	0.095	0.093	0.098	0.122	0.062
# CCSS	Thyroid cancer (123)	All_treatments	0.479	0.477	0.482	0.442	0.523
# CCSS	Thyroid cancer (123)	PRS	0.385	0.389	0.376	0.378	0.393
# CCSS	Thyroid cancer (123)	Lifestyle	-	-	-	-	-
# CCSS	Thyroid cancer (123)	Combined	0.667	0.666	0.669	0.646	0.693
# CCSS	Meningioma (210)	Radiation	0.297	0.264	0.348	0.319	0.282
# CCSS	Meningioma (210)	Chemotherapy	0.087	0.079	0.098	0.132	0.055
# CCSS	Meningioma (210)	All_treatments	0.371	0.337	0.422	0.423	0.334
# CCSS	Meningioma (210)	PRS	0.023	0.023	0.024	0.028	0.02
# CCSS	Meningioma (210)	Lifestyle	0.058	0.07	0.04	0.123	0.012
# CCSS	Meningioma (210)	Combined	0.413	0.387	0.451	0.506	0.347
# CCSS	Sarcoma (41)	Radiation	-	-	-	-	-
# CCSS	Sarcoma (41)	Chemotherapy	0.222	0.213	0.237	0.236	0.214
# CCSS	Sarcoma (41)	All_treatments	0.222	0.213	0.237	0.236	0.214
# CCSS	Sarcoma (41)	PRS	0.046	0.048	0.042	0.046	0.046
# CCSS	Sarcoma (41)	Lifestyle	-	-	-	-	-
# CCSS	Sarcoma (41)	Combined	-0.035	-0.036	-0.032	-0.108	0.01", header = T, sep = "\t")
# 
# data$SN_types_new <- gsub("\\([0-9]+\\)", "", data$SN_types)
# data[data == "-"] <- NA
# 
# 
# lifestyle <- "without_diet"
# all.group <- c("Overall", "Female", "Male", "Age.lt.35", "Age.ge.35")
# variables <-  c("Combined", "Radiation", "Chemotherapy", "treatments", "PRS", "Lifestyle")
# 
# 
# # Define the text file path
# data_name <- paste0("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/ANALYSIS/results/plots/v21b/", lifestyle, "_all_table_19b.txt")
# 
# file_conn <- file(data_name, "a")
# 
# for (j in 1:length(all.group)){
#   group <-  all.group [j]
#   for(i in 1:length(variables)){
#     if(lifestyle== "without_diet" && variables[i] == "Lifestyle"){
#     new_data <- data[grepl(variables[i], data$Variables), c("Cohort", "SN_types", "SN_types_new", "Variables", group)]
#     
#     AF.type <- variables[i]
#     
#     # Filter the data for SJLIFE and CCSS cohorts
#     sjlife_data <- new_data %>% filter(Cohort == "SJLIFE")
#     ccss_data <- new_data %>% filter(Cohort == "CCSS")
#     
#     # Create a new data frame with the desired format
#     result_data <- tibble(
#       "SN types (CA in SJLIFE/CA in CCSS)" = paste0(sjlife_data$SN_types_new, " (", as.numeric(gsub("[^0-9.]", "", sjlife_data$SN_types)), "/", as.numeric(gsub("[^0-9.]", "", ccss_data$SN_types)), ")"),
#       "SJLIFE" = as.numeric(sjlife_data[,group]),
#       "CCSS" = as.numeric(ccss_data[,group]))
#     
#     # Print the result_data data frame
#     table_name <- paste0(group, "_", AF.type, "_", lifestyle,"\tNA\tNA")
#     # write.table(table_name, data_name, col.names = T, row.names = F, quote = F, append = T, sep = "\t")
#     # write.table(result_data, data_name, col.names = T, row.names = F, quote = F, append = T, sep = "\t")
#     # }
#     
#     cat(table_name, "\n", file = file_conn)
#     write.table(result_data, file = file_conn, col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")
#     }
#   }
# }
# # Close the text file
# close(file_conn)



######################################
######################################
## create individual tables
## Lifestyle without diet
data <- read.table(text="Cohort	SNtypes	Variables	Overall_all4	Female_all4	Male_all4	Age_lt_35_all4	Age_ge_35_all4	Overall_smk	Female_smk	Male_smk	Age_lt_35_smk	Age_ge_35_smk	Overall_drk	Female_drk	Male_drk	Age_lt_35_drk	Age_ge_35_drk	Overall_PA	Female_PA	Male_PA	Age_lt_35_PA	Age_ge_35_PA	Overall_obese	Female_obese	Male_obese	Age_lt_35_obese	Age_ge_35_obese
SJLIFE	SNs	Radiation	0.442	0.432	0.455	0.367	0.472	0.442	0.432	0.455	0.367	0.472	0.442	0.432	0.455	0.367	0.472	0.442	0.432	0.455	0.367	0.472	0.442	0.432	0.455	0.367	0.472
SJLIFE	SNs	Chemotherapy	0.015	0.016	0.014	0.027	0.01	0.015	0.016	0.014	0.027	0.01	0.015	0.016	0.014	0.027	0.01	0.015	0.016	0.014	0.027	0.01	0.015	0.016	0.014	0.027	0.01
SJLIFE	SNs	All_treatments	0.449	0.441	0.46	0.382	0.476	0.449	0.441	0.46	0.382	0.476	0.449	0.441	0.46	0.382	0.476	0.449	0.441	0.46	0.382	0.476	0.449	0.441	0.46	0.382	0.476
SJLIFE	SNs	PRS	0.171	0.169	0.174	0.174	0.17	0.171	0.169	0.174	0.174	0.17	0.171	0.169	0.174	0.174	0.17	0.171	0.169	0.174	0.174	0.17	0.171	0.169	0.174	0.174	0.17
SJLIFE	SNs	Lifestyle	-0.038	-0.03	-0.05	-0.048	-0.034	-0.027	-0.026	-0.029	-0.036	-0.023	-0.04	-0.036	-0.045	-0.043	-0.039	0.009	0.011	0.005	0.005	0.01	-0.019	-0.019	-0.018	-0.02	-0.018
SJLIFE	SNs	Combined	0.525	0.521	0.53	0.461	0.55	0.531	0.524	0.541	0.469	0.556	0.525	0.518	0.533	0.465	0.549	0.548	0.541	0.557	0.49	0.571	0.536	0.528	0.546	0.478	0.559
SJLIFE	SMNs	Radiation	0.378	0.377	0.379	0.289	0.412	0.378	0.377	0.379	0.289	0.412	0.378	0.377	0.379	0.289	0.412	0.378	0.377	0.379	0.289	0.412	0.378	0.377	0.379	0.289	0.412
SJLIFE	SMNs	Chemotherapy	0.014	0.014	0.013	0.022	0.01	0.014	0.014	0.013	0.022	0.01	0.014	0.014	0.013	0.022	0.01	0.014	0.014	0.013	0.022	0.01	0.014	0.014	0.013	0.022	0.01
SJLIFE	SMNs	All_treatments	0.387	0.387	0.386	0.303	0.418	0.387	0.387	0.386	0.303	0.418	0.387	0.387	0.386	0.303	0.418	0.387	0.387	0.386	0.303	0.418	0.387	0.387	0.386	0.303	0.418
SJLIFE	SMNs	PRS	0.133	0.13	0.136	0.132	0.133	0.133	0.13	0.136	0.132	0.133	0.133	0.13	0.136	0.132	0.133	0.133	0.13	0.136	0.132	0.133	0.133	0.13	0.136	0.132	0.133
SJLIFE	SMNs	Lifestyle	-0.057	-0.044	-0.073	-0.068	-0.052	-0.049	-0.049	-0.049	-0.064	-0.043	-0.055	-0.054	-0.057	-0.063	-0.052	0.032	0.04	0.021	0.017	0.037	-0.102	-0.104	-0.099	-0.096	-0.104
SJLIFE	SMNs	Combined	0.436	0.443	0.427	0.351	0.468	0.442	0.442	0.441	0.355	0.474	0.438	0.438	0.438	0.355	0.47	0.485	0.488	0.48	0.403	0.515	0.413	0.412	0.414	0.335	0.443
SJLIFE	NMSCs	Radiation	0.298	0.286	0.312	0.239	0.31	0.298	0.286	0.312	0.239	0.31	0.298	0.286	0.312	0.239	0.31	0.298	0.286	0.312	0.239	0.31	0.298	0.286	0.312	0.239	0.31
SJLIFE	NMSCs	Chemotherapy	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-
SJLIFE	NMSCs	All_treatments	0.298	0.286	0.312	0.239	0.31	0.298	0.286	0.312	0.239	0.31	0.298	0.286	0.312	0.239	0.31	0.298	0.286	0.312	0.239	0.31	0.298	0.286	0.312	0.239	0.31
SJLIFE	NMSCs	PRS	0.311	0.293	0.334	0.308	0.311	0.311	0.293	0.334	0.308	0.311	0.311	0.293	0.334	0.308	0.311	0.311	0.293	0.334	0.308	0.311	0.311	0.293	0.334	0.308	0.311
SJLIFE	NMSCs	Lifestyle	-0.207	-0.188	-0.231	-0.201	-0.208	0.042	0.038	0.048	0.003	0.051	-0.156	-0.137	-0.18	-0.148	-0.157	-0.054	-0.058	-0.048	-0.064	-0.051	-0.146	-0.145	-0.146	-0.134	-0.148
SJLIFE	NMSCs	Combined	0.416	0.403	0.434	0.363	0.428	0.537	0.517	0.563	0.476	0.55	0.416	0.403	0.434	0.363	0.428	0.49	0.468	0.518	0.438	0.501	0.445	0.424	0.472	0.401	0.454
SJLIFE	Breast_cancer	Radiation	0.47	0.47	-	0.387	0.483	0.47	0.47	-	0.387	0.483	0.47	0.47	-	0.387	0.483	0.47	0.47	-	0.387	0.483	0.47	0.47	-	0.387	0.483
SJLIFE	Breast_cancer	Chemotherapy	0.072	0.072	-	0.056	0.074	0.072	0.072	-	0.056	0.074	0.072	0.072	-	0.056	0.074	0.072	0.072	-	0.056	0.074	0.072	0.072	-	0.056	0.074
SJLIFE	Breast_cancer	All_treatments	0.511	0.511	-	0.42	0.526	0.511	0.511	-	0.42	0.526	0.511	0.511	-	0.42	0.526	0.511	0.511	-	0.42	0.526	0.511	0.511	-	0.42	0.526
SJLIFE	Breast_cancer	PRS	0.177	0.177	-	0.198	0.174	0.177	0.177	-	0.198	0.174	0.177	0.177	-	0.198	0.174	0.177	0.177	-	0.198	0.174	0.177	0.177	-	0.198	0.174
SJLIFE	Breast_cancer	Lifestyle	-0.67	-0.67	-	-0.456	-0.704	-0.166	-0.166	-	-0.094	-0.178	0.055	0.055	-	0.057	0.055	-0.262	-0.262	-	-0.163	-0.278	0.051	0.051	-	0.054	0.05
SJLIFE	Breast_cancer	Combined	0.341	0.341	-	0.314	0.345	0.538	0.538	-	0.493	0.545	0.623	0.623	-	0.564	0.632	0.502	0.502	-	0.456	0.51	0.619	0.619	-	0.563	0.628
SJLIFE	Thyroid_cancer	Radiation	0.601	0.62	0.577	0.535	0.63	0.601	0.62	0.577	0.535	0.63	0.601	0.62	0.577	0.535	0.63	0.601	0.62	0.577	0.535	0.63	0.601	0.62	0.577	0.535	0.63
SJLIFE	Thyroid_cancer	Chemotherapy	0.221	0.222	0.22	0.312	0.181	0.221	0.222	0.22	0.312	0.181	0.221	0.222	0.22	0.312	0.181	0.221	0.222	0.22	0.312	0.181	0.221	0.222	0.22	0.312	0.181
SJLIFE	Thyroid_cancer	All_treatments	0.71	0.725	0.693	0.686	0.721	0.71	0.725	0.693	0.686	0.721	0.71	0.725	0.693	0.686	0.721	0.71	0.725	0.693	0.686	0.721	0.71	0.725	0.693	0.686	0.721
SJLIFE	Thyroid_cancer	PRS	0.572	0.571	0.573	0.583	0.567	0.572	0.571	0.573	0.583	0.567	0.572	0.571	0.573	0.583	0.567	0.572	0.571	0.573	0.583	0.567	0.572	0.571	0.573	0.583	0.567
SJLIFE	Thyroid_cancer	Lifestyle	-0.23	-0.245	-0.213	-0.235	-0.228	0.031	0.033	0.028	-0.005	0.046	-0.134	-0.121	-0.149	-0.136	-0.132	-0.181	-0.206	-0.152	-0.158	-0.192	-0.093	-0.094	-0.093	-0.085	-0.097
SJLIFE	Thyroid_cancer	Combined	0.849	0.854	0.842	0.833	0.856	0.881	0.886	0.874	0.864	0.888	0.86	0.867	0.851	0.845	0.866	0.855	0.858	0.851	0.842	0.86	0.865	0.871	0.858	0.853	0.871
SJLIFE	Meningioma	Radiation	0.185	0.151	0.225	0.251	0.162	0.185	0.151	0.225	0.251	0.162	0.185	0.151	0.225	0.251	0.162	0.185	0.151	0.225	0.251	0.162	0.185	0.151	0.225	0.251	0.162
SJLIFE	Meningioma	Chemotherapy	0.318	0.323	0.312	0.454	0.273	0.318	0.323	0.312	0.454	0.273	0.318	0.323	0.312	0.454	0.273	0.318	0.323	0.312	0.454	0.273	0.318	0.323	0.312	0.454	0.273
SJLIFE	Meningioma	All_treatments	0.456	0.436	0.479	0.623	0.4	0.456	0.436	0.479	0.623	0.4	0.456	0.436	0.479	0.623	0.4	0.456	0.436	0.479	0.623	0.4	0.456	0.436	0.479	0.623	0.4
SJLIFE	Meningioma	PRS	-0.123	-0.119	-0.127	-0.12	-0.124	-0.123	-0.119	-0.127	-0.12	-0.124	-0.123	-0.119	-0.127	-0.12	-0.124	-0.123	-0.119	-0.127	-0.12	-0.124	-0.123	-0.119	-0.127	-0.12	-0.124
SJLIFE	Meningioma	Lifestyle	-0.264	-0.237	-0.295	-0.229	-0.275	-0.172	-0.151	-0.196	-0.131	-0.185	-0.126	-0.101	-0.157	-0.117	-0.13	-0.086	-0.101	-0.069	-0.074	-0.091	0.151	0.145	0.158	0.127	0.159
SJLIFE	Meningioma	Combined	0.217	0.21	0.226	0.461	0.136	0.277	0.265	0.292	0.512	0.198	0.309	0.304	0.315	0.52	0.238	0.338	0.309	0.373	0.547	0.268	0.482	0.461	0.506	0.628	0.432
CCSS	SNs	Radiation	0.353	0.346	0.364	0.332	0.364	0.353	0.346	0.364	0.332	0.364	0.353	0.346	0.364	0.332	0.364	0.353	0.346	0.364	0.332	0.364	0.353	0.346	0.364	0.332	0.364
CCSS	SNs	Chemotherapy	0.038	0.036	0.04	0.06	0.027	0.038	0.036	0.04	0.06	0.027	0.038	0.036	0.04	0.06	0.027	0.038	0.036	0.04	0.06	0.027	0.038	0.036	0.04	0.06	0.027
CCSS	SNs	All_treatments	0.379	0.372	0.391	0.374	0.382	0.379	0.372	0.391	0.374	0.382	0.379	0.372	0.391	0.374	0.382	0.379	0.372	0.391	0.374	0.382	0.379	0.372	0.391	0.374	0.382
CCSS	SNs	PRS	0.058	0.057	0.059	0.056	0.058	0.058	0.057	0.059	0.056	0.058	0.058	0.057	0.059	0.056	0.058	0.058	0.057	0.059	0.056	0.058	0.058	0.057	0.059	0.056	0.058
CCSS	SNs	Lifestyle	-0.104	-0.103	-0.106	-0.093	-0.11	-0.073	-0.07	-0.076	-0.062	-0.078	-0.005	-0.005	-0.005	-0.005	-0.005	-0.007	-0.007	-0.006	-0.007	-0.007	-0.017	-0.019	-0.016	-0.018	-0.017
CCSS	SNs	Combined	0.352	0.345	0.363	0.352	0.352	0.37	0.363	0.38	0.37	0.37	0.412	0.405	0.423	0.406	0.415	0.411	0.403	0.422	0.405	0.414	0.405	0.397	0.417	0.399	0.408
CCSS	SMNs	Radiation	0.245	0.241	0.251	0.182	0.271	0.245	0.241	0.251	0.182	0.271	0.245	0.241	0.251	0.182	0.271	0.245	0.241	0.251	0.182	0.271	0.245	0.241	0.251	0.182	0.271
CCSS	SMNs	Chemotherapy	0.055	0.054	0.057	0.098	0.038	0.055	0.054	0.057	0.098	0.038	0.055	0.054	0.057	0.098	0.038	0.055	0.054	0.057	0.098	0.038	0.055	0.054	0.057	0.098	0.038
CCSS	SMNs	All_treatments	0.293	0.288	0.301	0.266	0.304	0.293	0.288	0.301	0.266	0.304	0.293	0.288	0.301	0.266	0.304	0.293	0.288	0.301	0.266	0.304	0.293	0.288	0.301	0.266	0.304
CCSS	SMNs	PRS	0.047	0.046	0.048	0.046	0.047	0.047	0.046	0.048	0.046	0.047	0.047	0.046	0.048	0.046	0.047	0.047	0.046	0.048	0.046	0.047	0.047	0.046	0.048	0.046	0.047
CCSS	SMNs	Lifestyle	-0.078	-0.077	-0.08	-0.068	-0.082	-0.06	-0.059	-0.062	-0.049	-0.065	-0.004	-0.004	-0.005	-0.004	-0.004	-0.008	-0.008	-0.008	-0.007	-0.008	-0.001	-0.001	-0.001	-0.001	-0.002
CCSS	SMNs	Combined	0.273	0.27	0.28	0.252	0.282	0.286	0.282	0.292	0.265	0.294	0.323	0.319	0.331	0.297	0.334	0.32	0.316	0.329	0.295	0.331	0.325	0.32	0.333	0.3	0.335
CCSS	NMSCs	Radiation	0.408	0.407	0.411	0.415	0.406	0.408	0.407	0.411	0.415	0.406	0.408	0.407	0.411	0.415	0.406	0.408	0.407	0.411	0.415	0.406	0.408	0.407	0.411	0.415	0.406
CCSS	NMSCs	Chemotherapy	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-
CCSS	NMSCs	All_treatments	0.408	0.407	0.411	0.415	0.406	0.408	0.407	0.411	0.415	0.406	0.408	0.407	0.411	0.415	0.406	0.408	0.407	0.411	0.415	0.406	0.408	0.407	0.411	0.415	0.406
CCSS	NMSCs	PRS	0.335	0.333	0.338	0.336	0.335	0.335	0.333	0.338	0.336	0.335	0.335	0.333	0.338	0.336	0.335	0.335	0.333	0.338	0.336	0.335	0.335	0.333	0.338	0.336	0.335
CCSS	NMSCs	Lifestyle	-0.225	-0.23	-0.22	-0.21	-0.231	-0.117	-0.111	-0.124	-0.108	-0.12	0.007	0.005	0.01	0.002	0.009	-0.074	-0.078	-0.069	-0.079	-0.072	-0.061	-0.067	-0.054	-0.068	-0.059
CCSS	NMSCs	Combined	0.52	0.515	0.526	0.532	0.515	0.56	0.559	0.561	0.568	0.557	0.612	0.608	0.615	0.615	0.61	0.579	0.575	0.584	0.583	0.578	0.586	0.583	0.591	0.589	0.585
CCSS	Breast_cancer	Radiation	0.444	0.444	NaN	0.392	0.45	0.444	0.444	NaN	0.392	0.45	0.444	0.444	NaN	0.392	0.45	0.444	0.444	NaN	0.392	0.45	0.444	0.444	NaN	0.392	0.45
CCSS	Breast_cancer	Chemotherapy	0.193	0.193	NaN	0.232	0.188	0.193	0.193	NaN	0.232	0.188	0.193	0.193	NaN	0.232	0.188	0.193	0.193	NaN	0.232	0.188	0.193	0.193	NaN	0.232	0.188
CCSS	Breast_cancer	All_treatments	0.572	0.572	NaN	0.541	0.575	0.572	0.572	NaN	0.541	0.575	0.572	0.572	NaN	0.541	0.575	0.572	0.572	NaN	0.541	0.575	0.572	0.572	NaN	0.541	0.575
CCSS	Breast_cancer	PRS	0.348	0.348	NaN	0.35	0.348	0.348	0.348	NaN	0.35	0.348	0.348	0.348	NaN	0.35	0.348	0.348	0.348	NaN	0.35	0.348	0.348	0.348	NaN	0.35	0.348
CCSS	Breast_cancer	Lifestyle	-0.111	-0.111	NaN	-0.1	-0.112	-0.048	-0.048	NaN	-0.04	-0.049	0.013	0.013	NaN	0.019	0.013	-0.033	-0.033	NaN	-0.029	-0.034	-0.027	-0.027	NaN	-0.025	-0.027
CCSS	Breast_cancer	Combined	0.69	0.69	NaN	0.673	0.692	0.708	0.708	NaN	0.691	0.71	0.726	0.726	NaN	0.709	0.728	0.713	0.713	NaN	0.695	0.715	0.714	0.714	NaN	0.695	0.716
CCSS	Thyroid_cancer	Radiation	0.416	0.416	0.417	0.359	0.486	0.416	0.416	0.417	0.359	0.486	0.416	0.416	0.417	0.359	0.486	0.416	0.416	0.417	0.359	0.486	0.416	0.416	0.417	0.359	0.486
CCSS	Thyroid_cancer	Chemotherapy	0.095	0.093	0.098	0.122	0.062	0.095	0.093	0.098	0.122	0.062	0.095	0.093	0.098	0.122	0.062	0.095	0.093	0.098	0.122	0.062	0.095	0.093	0.098	0.122	0.062
CCSS	Thyroid_cancer	All_treatments	0.479	0.477	0.482	0.442	0.523	0.479	0.477	0.482	0.442	0.523	0.479	0.477	0.482	0.442	0.523	0.479	0.477	0.482	0.442	0.523	0.479	0.477	0.482	0.442	0.523
CCSS	Thyroid_cancer	PRS	0.385	0.389	0.376	0.378	0.393	0.385	0.389	0.376	0.378	0.393	0.385	0.389	0.376	0.378	0.393	0.385	0.389	0.376	0.378	0.393	0.385	0.389	0.376	0.378	0.393
CCSS	Thyroid_cancer	Lifestyle	-0.039	-0.039	-0.04	-0.024	-0.059	-0.015	-0.013	-0.02	-0.002	-0.032	0.05	0.05	0.051	0.061	0.037	0.009	0.009	0.008	0.018	-0.003	0.003	0.003	0.002	0.011	-0.007
CCSS	Thyroid_cancer	Combined	0.667	0.666	0.669	0.646	0.693	0.675	0.675	0.674	0.654	0.7	0.696	0.695	0.697	0.676	0.72	0.683	0.682	0.683	0.661	0.708	0.681	0.68	0.682	0.659	0.707
CCSS	Meningioma	Radiation	0.297	0.264	0.348	0.319	0.282	0.297	0.264	0.348	0.319	0.282	0.297	0.264	0.348	0.319	0.282	0.297	0.264	0.348	0.319	0.282	0.297	0.264	0.348	0.319	0.282
CCSS	Meningioma	Chemotherapy	0.087	0.079	0.098	0.132	0.055	0.087	0.079	0.098	0.132	0.055	0.087	0.079	0.098	0.132	0.055	0.087	0.079	0.098	0.132	0.055	0.087	0.079	0.098	0.132	0.055
CCSS	Meningioma	All_treatments	0.371	0.337	0.422	0.423	0.334	0.371	0.337	0.422	0.423	0.334	0.371	0.337	0.422	0.423	0.334	0.371	0.337	0.422	0.423	0.334	0.371	0.337	0.422	0.423	0.334
CCSS	Meningioma	PRS	0.023	0.023	0.024	0.028	0.02	0.023	0.023	0.024	0.028	0.02	0.023	0.023	0.024	0.028	0.02	0.023	0.023	0.024	0.028	0.02	0.023	0.023	0.024	0.028	0.02
CCSS	Meningioma	Lifestyle	0.058	0.07	0.04	0.123	0.012	0.074	0.081	0.064	0.138	0.029	0.054	0.063	0.041	0.115	0.011	0.11	0.116	0.1	0.168	0.069	0.122	0.129	0.111	0.182	0.08
CCSS	Meningioma	Combined	0.413	0.387	0.451	0.506	0.347	0.428	0.401	0.47	0.52	0.364	0.415	0.388	0.455	0.506	0.35	0.453	0.426	0.494	0.539	0.392	0.459	0.432	0.499	0.544	0.399
CCSS	Sarcoma	Radiation	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-
CCSS	Sarcoma	Chemotherapy	0.222	0.213	0.237	0.236	0.214	0.222	0.213	0.237	0.236	0.214	0.222	0.213	0.237	0.236	0.214	0.222	0.213	0.237	0.236	0.214	0.222	0.213	0.237	0.236	0.214
CCSS	Sarcoma	All_treatments	0.222	0.213	0.237	0.236	0.214	0.222	0.213	0.237	0.236	0.214	0.222	0.213	0.237	0.236	0.214	0.222	0.213	0.237	0.236	0.214	0.222	0.213	0.237	0.236	0.214
CCSS	Sarcoma	PRS	0.046	0.048	0.042	0.046	0.046	0.046	0.048	0.042	0.046	0.046	0.046	0.048	0.042	0.046	0.046	0.046	0.048	0.042	0.046	0.046	0.046	0.048	0.042	0.046	0.046
CCSS	Sarcoma	Lifestyle	-0.392	-0.384	-0.406	-0.515	-0.317	-0.241	-0.243	-0.236	-0.325	-0.19	-0.217	-0.215	-0.22	-0.35	-0.136	-0.098	-0.104	-0.086	-0.199	-0.037	-0.137	-0.145	-0.123	-0.241	-0.074
CCSS	Sarcoma	Combined	-0.035	-0.036	-0.032	-0.108	0.01	0.076	0.067	0.091	0.031	0.103	0.098	0.093	0.107	0.015	0.148	-0.035	-0.036	-0.032	-0.108	0.01	0.156	0.143	0.178	0.093	0.193", header = T, sep = "\t")

## combined
data.sjlife <- data[data$Cohort=="SJLIFE",]
data.ccss <- data[data$Cohort=="CCSS",]

data.sjlife.wanted <- data.sjlife[grepl("Lifestyle", data.sjlife$Variables),]
data.ccss.wanted <- data.ccss[grepl("Lifestyle", data.ccss$Variables),]

variables <- c("drk", "smk", "PA", "obese", 'all4')


i=5
sjlife <- data.sjlife.wanted[grepl(paste0("SNtypes|",variables[i]), colnames(data.sjlife.wanted))]
sjlife <- rbind.data.frame(sjlife, c("Sarcoma",NA,NA,NA,NA,NA))
colnames(sjlife)<- paste0(colnames(sjlife), "_sjlife")
ccss <- data.ccss.wanted[grepl(paste0("SNtypes|",variables[i]), colnames(data.ccss.wanted))]
colnames(ccss)<- paste0(colnames(ccss), "_ccss")

wanted <- cbind(sjlife,ccss)[order(c(seq_along(sjlife), seq_along(ccss)))]
View(wanted)

###########################################################
###########################################################
## Lifestyle with diet HEI with HEI score less than 60 binary (SJLIFE only)
data <- read.table(text="SNtypes	Variables	Overall_all4	Female_all4	Male_all4	Age_lt_35_all4	Age_ge_35_all4	Overall_smk	Female_smk	Male_smk	Age_lt_35_smk	Age_ge_35_smk	Overall_drk	Female_drk	Male_drk	Age_lt_35_drk	Age_ge_35_drk	Overall_PA	Female_PA	Male_PA	Age_lt_35_PA	Age_ge_35_PA	Overall_obese	Female_obese	Male_obese	Age_lt_35_obese	Age_ge_35_obese	Overall_diet	Female_diet	Male_diet	Age_lt_35_diet	Age_ge_35_diet
SNs	Radiation	0.442	0.433	0.455	0.367	0.472	0.442	0.433	0.455	0.367	0.472	0.442	0.433	0.455	0.367	0.472	0.442	0.433	0.455	0.367	0.472	0.442	0.433	0.455	0.367	0.472	0.442	0.433	0.455	0.367	0.472
SNs	Chemotherapy	0.015	0.016	0.014	0.027	0.01	0.015	0.016	0.014	0.027	0.01	0.015	0.016	0.014	0.027	0.01	0.015	0.016	0.014	0.027	0.01	0.015	0.016	0.014	0.027	0.01	0.015	0.016	0.014	0.027	0.01
SNs	All_treatments	0.449	0.441	0.46	0.382	0.476	0.449	0.441	0.46	0.382	0.476	0.449	0.441	0.46	0.382	0.476	0.449	0.441	0.46	0.382	0.476	0.449	0.441	0.46	0.382	0.476	0.449	0.441	0.46	0.382	0.476
SNs	PRS	0.171	0.169	0.174	0.174	0.17	0.171	0.169	0.174	0.174	0.17	0.171	0.169	0.174	0.174	0.17	0.171	0.169	0.174	0.174	0.17	0.171	0.169	0.174	0.174	0.17	0.171	0.169	0.174	0.174	0.17
SNs	Lifestyle	-0.037	-0.027	-0.05	-0.047	-0.033	-0.024	-0.023	-0.026	-0.032	-0.021	-0.037	-0.033	-0.042	-0.039	-0.036	0.011	0.014	0.008	0.008	0.013	-0.014	-0.014	-0.013	-0.015	-0.014	-0.015	-0.014	-0.016	-0.017	-0.014
SNs	Combined	0.525	0.522	0.53	0.462	0.551	0.533	0.526	0.542	0.471	0.557	0.526	0.52	0.535	0.467	0.55	0.549	0.543	0.558	0.492	0.572	0.538	0.53	0.548	0.48	0.561	0.538	0.53	0.547	0.479	0.561
SMNs	Radiation	0.378	0.377	0.379	0.288	0.412	0.378	0.377	0.379	0.288	0.412	0.378	0.377	0.379	0.288	0.412	0.378	0.377	0.379	0.288	0.412	0.378	0.377	0.379	0.288	0.412	0.378	0.377	0.379	0.288	0.412
SMNs	Chemotherapy	0.013	0.013	0.013	0.021	0.01	0.013	0.013	0.013	0.021	0.01	0.013	0.013	0.013	0.021	0.01	0.013	0.013	0.013	0.021	0.01	0.013	0.013	0.013	0.021	0.01	0.013	0.013	0.013	0.021	0.01
SMNs	All_treatments	0.386	0.387	0.386	0.302	0.418	0.386	0.387	0.386	0.302	0.418	0.386	0.387	0.386	0.302	0.418	0.386	0.387	0.386	0.302	0.418	0.386	0.387	0.386	0.302	0.418	0.386	0.387	0.386	0.302	0.418
SMNs	PRS	0.132	0.13	0.134	0.131	0.132	0.132	0.13	0.134	0.131	0.132	0.132	0.13	0.134	0.131	0.132	0.132	0.13	0.134	0.131	0.132	0.132	0.13	0.134	0.131	0.132	0.132	0.13	0.134	0.131	0.132
SMNs	Lifestyle	-0.09	-0.075	-0.11	-0.106	-0.084	-0.066	-0.067	-0.064	-0.084	-0.059	-0.077	-0.077	-0.077	-0.087	-0.073	0.016	0.023	0.006	-0.002	0.022	-0.121	-0.125	-0.115	-0.118	-0.122	-0.069	-0.068	-0.07	-0.08	-0.064
SMNs	Combined	0.417	0.425	0.405	0.327	0.45	0.431	0.431	0.432	0.342	0.465	0.425	0.425	0.426	0.339	0.458	0.475	0.478	0.471	0.391	0.507	0.402	0.4	0.404	0.321	0.432	0.43	0.43	0.431	0.344	0.463
NMSCs	Radiation	0.3	0.289	0.314	0.241	0.313	0.3	0.289	0.314	0.241	0.313	0.3	0.289	0.314	0.241	0.313	0.3	0.289	0.314	0.241	0.313	0.3	0.289	0.314	0.241	0.313	0.3	0.289	0.314	0.241	0.313
NMSCs	Chemotherapy	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-
NMSCs	All_treatments	0.3	0.289	0.314	0.241	0.313	0.3	0.289	0.314	0.241	0.313	0.3	0.289	0.314	0.241	0.313	0.3	0.289	0.314	0.241	0.313	0.3	0.289	0.314	0.241	0.313	0.3	0.289	0.314	0.241	0.313
NMSCs	PRS	0.322	0.304	0.346	0.32	0.323	0.322	0.304	0.346	0.32	0.323	0.322	0.304	0.346	0.32	0.323	0.322	0.304	0.346	0.32	0.323	0.322	0.304	0.346	0.32	0.323	0.322	0.304	0.346	0.32	0.323
NMSCs	Lifestyle	-0.451	-0.355	-0.574	-0.486	-0.443	0.028	0.021	0.037	-0.028	0.04	-0.205	-0.185	-0.231	-0.209	-0.204	-0.055	-0.058	-0.051	-0.083	-0.049	-0.174	-0.176	-0.172	-0.18	-0.173	-0.288	-0.232	-0.359	-0.335	-0.277
NMSCs	Combined	0.316	0.333	0.294	0.229	0.334	0.539	0.518	0.566	0.472	0.553	0.427	0.414	0.443	0.367	0.44	0.498	0.477	0.525	0.439	0.511	0.441	0.419	0.47	0.389	0.453	0.391	0.393	0.389	0.31	0.409
Breast_cancer	Radiation	0.471	0.471	-	0.39	0.484	0.471	0.471	-	0.39	0.484	0.471	0.471	-	0.39	0.484	0.471	0.471	-	0.39	0.484	0.471	0.471	-	0.39	0.484	0.471	0.471	-	0.39	0.484
Breast_cancer	Chemotherapy	0.071	0.071	-	0.055	0.074	0.071	0.071	-	0.055	0.074	0.071	0.071	-	0.055	0.074	0.071	0.071	-	0.055	0.074	0.071	0.071	-	0.055	0.074	0.071	0.071	-	0.055	0.074
Breast_cancer	All_treatments	0.512	0.512	-	0.422	0.526	0.512	0.512	-	0.422	0.526	0.512	0.512	-	0.422	0.526	0.512	0.512	-	0.422	0.526	0.512	0.512	-	0.422	0.526	0.512	0.512	-	0.422	0.526
Breast_cancer	PRS	0.179	0.179	-	0.2	0.176	0.179	0.179	-	0.2	0.176	0.179	0.179	-	0.2	0.176	0.179	0.179	-	0.2	0.176	0.179	0.179	-	0.2	0.176	0.179	0.179	-	0.2	0.176
Breast_cancer	Lifestyle	-0.697	-0.697	-	-0.475	-0.733	-0.195	-0.195	-	-0.12	-0.207	0.031	0.031	-	0.034	0.03	-0.297	-0.297	-	-0.192	-0.313	0.024	0.024	-	0.029	0.023	0.07	0.07	-	0.074	0.07
Breast_cancer	Combined	0.334	0.334	-	0.31	0.338	0.53	0.53	-	0.484	0.538	0.616	0.616	-	0.556	0.625	0.492	0.492	-	0.444	0.5	0.611	0.611	-	0.554	0.62	0.631	0.631	-	0.575	0.64
Thyroid_cancer	Radiation	0.6	0.621	0.575	0.53	0.631	0.6	0.621	0.575	0.53	0.631	0.6	0.621	0.575	0.53	0.631	0.6	0.621	0.575	0.53	0.631	0.6	0.621	0.575	0.53	0.631	0.6	0.621	0.575	0.53	0.631
Thyroid_cancer	Chemotherapy	0.22	0.221	0.219	0.31	0.18	0.22	0.221	0.219	0.31	0.18	0.22	0.221	0.219	0.31	0.18	0.22	0.221	0.219	0.31	0.18	0.22	0.221	0.219	0.31	0.18	0.22	0.221	0.219	0.31	0.18
Thyroid_cancer	All_treatments	0.709	0.725	0.69	0.68	0.721	0.709	0.725	0.69	0.68	0.721	0.709	0.725	0.69	0.68	0.721	0.709	0.725	0.69	0.68	0.721	0.709	0.725	0.69	0.68	0.721	0.709	0.725	0.69	0.68	0.721
Thyroid_cancer	PRS	0.569	0.568	0.57	0.58	0.564	0.569	0.568	0.57	0.58	0.564	0.569	0.568	0.57	0.58	0.564	0.569	0.568	0.57	0.58	0.564	0.569	0.568	0.57	0.58	0.564	0.569	0.568	0.57	0.58	0.564
Thyroid_cancer	Lifestyle	-0.494	-0.446	-0.552	-0.552	-0.468	0.005	0.009	-0.001	-0.04	0.024	-0.202	-0.185	-0.222	-0.205	-0.2	-0.211	-0.231	-0.186	-0.195	-0.218	-0.137	-0.138	-0.136	-0.133	-0.139	-0.285	-0.233	-0.347	-0.325	-0.267
Thyroid_cancer	Combined	0.813	0.827	0.795	0.786	0.825	0.876	0.882	0.869	0.857	0.885	0.85	0.859	0.839	0.832	0.857	0.85	0.854	0.844	0.834	0.856	0.859	0.865	0.85	0.843	0.865	0.839	0.852	0.823	0.817	0.848
Meningioma	Radiation	0.182	0.151	0.219	0.25	0.159	0.184	0.151	0.224	0.251	0.162	0.184	0.151	0.224	0.251	0.162	0.184	0.151	0.224	0.251	0.162	0.184	0.151	0.224	0.251	0.162	0.184	0.151	0.224	0.251	0.162
Meningioma	Chemotherapy	0.316	0.323	0.306	0.453	0.27	0.318	0.323	0.311	0.453	0.272	0.318	0.323	0.311	0.453	0.272	0.318	0.323	0.311	0.453	0.272	0.318	0.323	0.311	0.453	0.272	0.318	0.323	0.311	0.453	0.272
Meningioma	All_treatments	0.448	0.431	0.467	0.614	0.392	0.454	0.435	0.477	0.621	0.399	0.454	0.435	0.477	0.621	0.399	0.454	0.435	0.477	0.621	0.399	0.454	0.435	0.477	0.621	0.399	0.454	0.435	0.477	0.621	0.399
Meningioma	PRS	-0.122	-0.119	-0.126	-0.12	-0.123	-0.122	-0.119	-0.126	-0.12	-0.123	-0.122	-0.119	-0.126	-0.12	-0.123	-0.122	-0.119	-0.126	-0.12	-0.123	-0.122	-0.119	-0.126	-0.12	-0.123	-0.122	-0.119	-0.126	-0.12	-0.123
Meningioma	Lifestyle	-0.218	-0.199	-0.241	-0.186	-0.229	-0.158	-0.136	-0.183	-0.112	-0.173	-0.109	-0.084	-0.14	-0.096	-0.114	-0.074	-0.088	-0.057	-0.057	-0.079	0.166	0.161	0.172	0.145	0.173	0.036	0.034	0.038	0.043	0.033
Meningioma	Combined	0.23	0.225	0.237	0.468	0.15	0.283	0.272	0.297	0.518	0.205	0.317	0.312	0.323	0.527	0.247	0.344	0.315	0.378	0.553	0.273	0.489	0.469	0.512	0.634	0.44	0.408	0.388	0.432	0.592	0.347", header = T, sep = "\t")

## combined
data.sjlife <- data

data.sjlife.wanted <- data.sjlife[grepl("Lifestyle", data.sjlife$Variables),]

variables <- c("diet", "drk", "smk", "PA", "obese", 'all4')


i=6
sjlife <- data.sjlife.wanted[grepl(paste0("SNtypes|",variables[i]), colnames(data.sjlife.wanted))]
View(sjlife)
