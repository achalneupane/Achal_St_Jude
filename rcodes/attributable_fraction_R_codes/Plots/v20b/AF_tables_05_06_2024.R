# Load required libraries
library(ggplot2)
library(ggthemes)  # For additional themes
library(tidyverse)

## V18 b (without lifestyle)
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
SJLIFE	NMSC (251)	Radiation	0.436	0.427	0.448	0.408	0.452
SJLIFE	NMSC (251)	Chemotherapy	-	-	-	-	-
SJLIFE	NMSC (251)	All_treatments	0.436	0.427	0.448	0.408	0.452
SJLIFE	NMSC (251)	PRS	0.28	0.266	0.296	0.278	0.281
SJLIFE	NMSC (251)	Lifestyle	-	-	-	-	-
SJLIFE	NMSC (251)	Combined	0.593	0.581	0.607	0.571	0.605
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
SJLIFE	Meningioma (149)	Radiation	0.901	0.899	0.904	0.895	0.907
SJLIFE	Meningioma (149)	Chemotherapy	0.167	0.182	0.149	0.219	0.128
SJLIFE	Meningioma (149)	All_treatments	0.914	0.912	0.916	0.914	0.914
SJLIFE	Meningioma (149)	PRS	-	-	-	-	-
SJLIFE	Meningioma (149)	Lifestyle	-	-	-	-	-
SJLIFE	Meningioma (149)	Combined	0.904	0.903	0.906	0.904	0.904
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
CCSS	NMSC (774)	Radiation	0.391	0.387	0.397	0.386	0.394
CCSS	NMSC (774)	Chemotherapy	-	-	-	-	-
CCSS	NMSC (774)	All_treatments	0.391	0.387	0.397	0.386	0.394
CCSS	NMSC (774)	PRS	0.306	0.304	0.309	0.308	0.305
CCSS	NMSC (774)	Lifestyle	-	-	-	-	-
CCSS	NMSC (774)	Combined	0.579	0.575	0.583	0.577	0.58
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

data$SN_types_new <- gsub("\\([0-9]+\\)", "", data$SN_types)
data[data == "-"] <- NA


lifestyle <- "without_lifestyle"
all.group <- c("Overall", "Female", "Male", "Age.lt.35", "Age.ge.35")
variables <-  c("Combined", "Radiation", "Chemotherapy", "treatments", "PRS", "Lifestyle")


# Define the text file path
data_name <- paste0("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/ANALYSIS/results/plots/v20b/", lifestyle, "_all_table_19b.txt")

file_conn <- file(data_name, "a")

for (j in 1:length(all.group)){
group <-  all.group [j]

for(i in 1:length(variables)){
# print(paste0("Doing ", group, "_", AF.type, "_", lifestyle))
  
if(lifestyle== "without_lifestyle" && variables[i] == "Lifestyle"){
  next
}
    
AF.type <- variables[i]
new_data <- data[grepl(AF.type, data$Variables), c("Cohort", "SN_types", "SN_types_new", "Variables", group)]

# Filter the data for SJLIFE and CCSS cohorts
sjlife_data <- new_data %>% filter(Cohort == "SJLIFE")
ccss_data <- new_data %>% filter(Cohort == "CCSS")

# Create a new data frame with the desired format
result_data <- tibble(
  "SN types (CA in SJLIFE/CA in CCSS)" = paste0(sjlife_data$SN_types_new, " (", as.numeric(gsub("[^0-9.]", "", sjlife_data$SN_types)), "/", as.numeric(gsub("[^0-9.]", "", ccss_data$SN_types)), ")"),
  "SJLIFE" = as.numeric(sjlife_data[,group]),
  "CCSS" = as.numeric(ccss_data[,group]))

# Print the result_data data frame
table_name <- paste0(group, "_", AF.type, "_", lifestyle,"\tNA\tNA")
# write.table(table_name, data_name, col.names = T, row.names = F, quote = F, append = T, sep = "\t")
# write.table(result_data, data_name, col.names = T, row.names = F, quote = F, append = T, sep = "\t")
# }

cat(table_name, "\n", file = file_conn)
write.table(result_data, file = file_conn, col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")
}
}
# Close the text file
close(file_conn)

###########################
## Lifestyle without HEI ##
###########################
data <- read.table(text = "Cohort	SN_types	Variables	Overall	Female	Male	Age.lt.35	Age.ge.35
SJLIFE	SNs (303)	Radiation	0.442	0.432	0.455	0.368	0.473
SJLIFE	SNs (303)	Chemotherapy	0.017	0.017	0.016	0.031	0.011
SJLIFE	SNs (303)	All_treatments	0.45	0.442	0.461	0.385	0.477
SJLIFE	SNs (303)	PRS	0.179	0.177	0.183	0.182	0.178
SJLIFE	SNs (303)	Lifestyle	-0.063	-0.058	-0.069	-0.058	-0.065
SJLIFE	SNs (303)	Combined	0.52	0.514	0.527	0.462	0.543
SJLIFE	SMNs (234)	Radiation	0.379	0.378	0.379	0.29	0.413
SJLIFE	SMNs (234)	Chemotherapy	0.017	0.017	0.016	0.029	0.012
SJLIFE	SMNs (234)	All_treatments	0.389	0.39	0.388	0.308	0.421
SJLIFE	SMNs (234)	PRS	0.142	0.139	0.145	0.141	0.142
SJLIFE	SMNs (234)	Lifestyle	-0.083	-0.073	-0.095	-0.081	-0.083
SJLIFE	SMNs (234)	Combined	0.431	0.437	0.423	0.353	0.46
SJLIFE	NMSC (119)	Radiation	0.291	0.28	0.305	0.236	0.303
SJLIFE	NMSC (119)	Chemotherapy	-	-	-	-	-
SJLIFE	NMSC (119)	All_treatments	0.291	0.28	0.305	0.236	0.303
SJLIFE	NMSC (119)	PRS	0.311	0.293	0.334	0.308	0.311
SJLIFE	NMSC (119)	Lifestyle	-0.275	-0.253	-0.304	-0.241	-0.282
SJLIFE	NMSC (119)	Combined	0.376	0.364	0.393	0.336	0.385
SJLIFE	Breast cancer (51)	Radiation	0.466	0.466	-	0.392	0.478
SJLIFE	Breast cancer (51)	Chemotherapy	0.107	0.107	-	0.103	0.107
SJLIFE	Breast cancer (51)	All_treatments	0.527	0.527	-	0.454	0.538
SJLIFE	Breast cancer (51)	PRS	0.22	0.22	-	0.247	0.215
SJLIFE	Breast cancer (51)	Lifestyle	-0.882	-0.882	-	-0.596	-0.927
SJLIFE	Breast cancer (51)	Combined	0.331	0.331	-	0.333	0.331
SJLIFE	Thyroid cancer (44)	Radiation	0.61	0.634	0.583	0.554	0.638
SJLIFE	Thyroid cancer (44)	Chemotherapy	0.239	0.238	0.241	0.334	0.193
SJLIFE	Thyroid cancer (44)	All_treatments	0.725	0.741	0.706	0.709	0.733
SJLIFE	Thyroid cancer (44)	PRS	0.589	0.587	0.592	0.6	0.584
SJLIFE	Thyroid cancer (44)	Lifestyle	-0.175	-0.184	-0.164	-0.141	-0.191
SJLIFE	Thyroid cancer (44)	Combined	0.871	0.876	0.864	0.863	0.874
SJLIFE	Meningioma (81)	Radiation	0.185	0.15	0.226	0.246	0.164
SJLIFE	Meningioma (81)	Chemotherapy	0.321	0.324	0.317	0.461	0.273
SJLIFE	Meningioma (81)	All_treatments	0.456	0.436	0.481	0.623	0.4
SJLIFE	Meningioma (81)	PRS	-0.121	-0.116	-0.126	-0.118	-0.122
SJLIFE	Meningioma (81)	Lifestyle	-0.134	-0.125	-0.145	-0.136	-0.133
SJLIFE	Meningioma (81)	Combined	0.306	0.291	0.324	0.509	0.238
SJLIFE	Sarcoma (NA)	Radiation	-	-	-	-	-
SJLIFE	Sarcoma (NA)	Chemotherapy	-	-	-	-	-
SJLIFE	Sarcoma (NA)	All_treatments	-	-	-	-	-
SJLIFE	Sarcoma (NA)	PRS	-	-	-	-	-
SJLIFE	Sarcoma (NA)	Lifestyle	-	-	-	-	-
SJLIFE	Sarcoma (NA)	Combined	-	-	-	-	-
CCSS	SNs (1286)	Radiation	0.358	0.351	0.369	0.338	0.368
CCSS	SNs (1286)	Chemotherapy	0.038	0.037	0.041	0.061	0.027
CCSS	SNs (1286)	All_treatments	0.384	0.376	0.396	0.38	0.387
CCSS	SNs (1286)	PRS	0.057	0.056	0.058	0.056	0.057
CCSS	SNs (1286)	Lifestyle	-0.074	-0.071	-0.078	-0.068	-0.077
CCSS	SNs (1286)	Combined	0.376	0.37	0.386	0.374	0.378
CCSS	SMNs (619)	Radiation	0.249	0.245	0.255	0.186	0.275
CCSS	SMNs (619)	Chemotherapy	0.056	0.055	0.058	0.099	0.038
CCSS	SMNs (619)	All_treatments	0.297	0.292	0.306	0.271	0.308
CCSS	SMNs (619)	PRS	0.048	0.047	0.049	0.047	0.048
CCSS	SMNs (619)	Lifestyle	-0.048	-0.047	-0.052	-0.046	-0.049
CCSS	SMNs (619)	Combined	0.298	0.295	0.303	0.273	0.308
CCSS	NMSC (637)	Radiation	0.391	0.388	0.395	0.401	0.387
CCSS	NMSC (637)	Chemotherapy	-	-	-	-	-
CCSS	NMSC (637)	All_treatments	0.391	0.388	0.395	0.401	0.387
CCSS	NMSC (637)	PRS	0.308	0.306	0.311	0.309	0.308
CCSS	NMSC (637)	Lifestyle	-0.155	-0.154	-0.156	-0.138	-0.162
CCSS	NMSC (637)	Combined	0.518	0.515	0.521	0.533	0.512
CCSS	Breast cancer (260)	Radiation	0.445	0.445	-	0.394	0.451
CCSS	Breast cancer (260)	Chemotherapy	0.192	0.192	-	0.232	0.188
CCSS	Breast cancer (260)	All_treatments	0.572	0.572	-	0.542	0.576
CCSS	Breast cancer (260)	PRS	0.346	0.346	-	0.348	0.346
CCSS	Breast cancer (260)	Lifestyle	-0.085	-0.085	-	-0.076	-0.086
CCSS	Breast cancer (260)	Combined	0.697	0.697	-	0.68	0.699
CCSS	Thyroid cancer (123)	Radiation	0.416	0.416	0.416	0.359	0.485
CCSS	Thyroid cancer (123)	Chemotherapy	0.095	0.094	0.097	0.122	0.061
CCSS	Thyroid cancer (123)	All_treatments	0.479	0.477	0.481	0.443	0.522
CCSS	Thyroid cancer (123)	PRS	0.385	0.39	0.375	0.378	0.394
CCSS	Thyroid cancer (123)	Lifestyle	-0.042	-0.035	-0.054	-0.025	-0.062
CCSS	Thyroid cancer (123)	Combined	0.666	0.667	0.664	0.645	0.691
CCSS	Meningioma (210)	Radiation	0.297	0.264	0.348	0.319	0.282
CCSS	Meningioma (210)	Chemotherapy	0.087	0.079	0.098	0.132	0.055
CCSS	Meningioma (210)	All_treatments	0.371	0.337	0.422	0.424	0.334
CCSS	Meningioma (210)	PRS	0.023	0.022	0.024	0.028	0.02
CCSS	Meningioma (210)	Lifestyle	0.072	0.083	0.054	0.133	0.029
CCSS	Meningioma (210)	Combined	0.423	0.398	0.46	0.513	0.359
CCSS	Sarcoma (41)	Radiation	-	-	-	-	-
CCSS	Sarcoma (41)	Chemotherapy	0.227	0.217	0.243	0.243	0.217
CCSS	Sarcoma (41)	All_treatments	0.227	0.217	0.243	0.243	0.217
CCSS	Sarcoma (41)	PRS	0.05	0.053	0.045	0.048	0.051
CCSS	Sarcoma (41)	Lifestyle	-0.278	-0.273	-0.286	-0.426	-0.189
CCSS	Sarcoma (41)	Combined	0.062	0.058	0.068	-0.03	0.117", header = T, sep = "\t")

data$SN_types_new <- gsub("\\([0-9]+\\)", "", data$SN_types)
data[data == "-"] <- NA


lifestyle <- "without_diet"
all.group <- c("Overall", "Female", "Male", "Age.lt.35", "Age.ge.35")
variables <-  c("Combined", "Radiation", "Chemotherapy", "treatments", "PRS", "Lifestyle")


# Define the text file path
data_name <- paste0("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/ANALYSIS/results/plots/v20b/", lifestyle, "_all_table_19b.txt")

file_conn <- file(data_name, "a")

for (j in 1:length(all.group)){
  group <-  all.group [j]
  for(i in 1:length(variables)){
    if(lifestyle== "without_diet" && variables[i] == "Lifestyle"){
    new_data <- data[grepl(variables[i], data$Variables), c("Cohort", "SN_types", "SN_types_new", "Variables", group)]
    
    AF.type <- variables[i]
    
    # Filter the data for SJLIFE and CCSS cohorts
    sjlife_data <- new_data %>% filter(Cohort == "SJLIFE")
    ccss_data <- new_data %>% filter(Cohort == "CCSS")
    
    # Create a new data frame with the desired format
    result_data <- tibble(
      "SN types (CA in SJLIFE/CA in CCSS)" = paste0(sjlife_data$SN_types_new, " (", as.numeric(gsub("[^0-9.]", "", sjlife_data$SN_types)), "/", as.numeric(gsub("[^0-9.]", "", ccss_data$SN_types)), ")"),
      "SJLIFE" = as.numeric(sjlife_data[,group]),
      "CCSS" = as.numeric(ccss_data[,group]))
    
    # Print the result_data data frame
    table_name <- paste0(group, "_", AF.type, "_", lifestyle,"\tNA\tNA")
    # write.table(table_name, data_name, col.names = T, row.names = F, quote = F, append = T, sep = "\t")
    # write.table(result_data, data_name, col.names = T, row.names = F, quote = F, append = T, sep = "\t")
    # }
    
    cat(table_name, "\n", file = file_conn)
    write.table(result_data, file = file_conn, col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")
    }
  }
}
# Close the text file
close(file_conn)



###############################################
## Lifestyle with diet HEI2015 (SJLIFE only) ##
###############################################
data <- read.table(text="Cohort	SN_types	Variables	Overall	Female	Male	Age.lt.35	Age.ge.35
SJLIFE	Any SN (303)	Radiation	0.443	0.433	0.457	0.369	0.473
SJLIFE	Any SN (303)	Chemotherapy	0.02	0.02	0.019	0.036	0.013
SJLIFE	Any SN (303)	All_treatments	0.453	0.443	0.465	0.389	0.479
SJLIFE	Any SN (303)	PRS	0.18	0.178	0.183	0.183	0.179
SJLIFE	Any SN (303)	Lifestyle	-0.057	-0.062	-0.049	-0.05	-0.06
SJLIFE	Any SN (303)	Combined	0.525	0.515	0.538	0.47	0.548
SJLIFE	SMNs (234)	Radiation	0.378	0.378	0.379	0.289	0.413
SJLIFE	SMNs (234)	Chemotherapy	0.018	0.018	0.018	0.031	0.013
SJLIFE	SMNs (234)	All_treatments	0.39	0.39	0.389	0.309	0.421
SJLIFE	SMNs (234)	PRS	0.14	0.137	0.143	0.139	0.14
SJLIFE	SMNs (234)	Lifestyle	-0.13	-0.122	-0.139	-0.133	-0.128
SJLIFE	SMNs (234)	Combined	0.404	0.409	0.398	0.321	0.436
SJLIFE	NMSC (119)	Radiation	0.285	0.274	0.299	0.23	0.297
SJLIFE	NMSC (119)	Chemotherapy	-	-	-	-	-
SJLIFE	NMSC (119)	All_treatments	0.285	0.274	0.299	0.23	0.297
SJLIFE	NMSC (119)	PRS	0.328	0.311	0.349	0.326	0.328
SJLIFE	NMSC (119)	Lifestyle	-0.549	-0.459	-0.665	-0.544	-0.55
SJLIFE	NMSC (119)	Combined	0.255	0.269	0.238	0.191	0.269
SJLIFE	Breast cancer (51)	Radiation	0.47	0.47	-	0.4	0.481
SJLIFE	Breast cancer (51)	Chemotherapy	0.109	0.109	-	0.106	0.11
SJLIFE	Breast cancer (51)	All_treatments	0.532	0.532	-	0.462	0.543
SJLIFE	Breast cancer (51)	PRS	0.226	0.226	-	0.259	0.221
SJLIFE	Breast cancer (51)	Lifestyle	-0.917	-0.917	-	-0.602	-0.968
SJLIFE	Breast cancer (51)	Combined	0.342	0.342	-	0.35	0.341
SJLIFE	Thyroid cancer (44)	Radiation	0.611	0.635	0.584	0.55	0.641
SJLIFE	Thyroid cancer (44)	Chemotherapy	0.24	0.238	0.243	0.329	0.197
SJLIFE	Thyroid cancer (44)	All_treatments	0.727	0.743	0.708	0.704	0.738
SJLIFE	Thyroid cancer (44)	PRS	0.582	0.579	0.585	0.593	0.576
SJLIFE	Thyroid cancer (44)	Lifestyle	-0.677	-0.561	-0.811	-0.709	-0.662
SJLIFE	Thyroid cancer (44)	Combined	0.811	0.834	0.785	0.792	0.82
SJLIFE	Meningioma (81)	Radiation	0.182	0.15	0.22	0.245	0.161
SJLIFE	Meningioma (81)	Chemotherapy	0.319	0.325	0.312	0.463	0.271
SJLIFE	Meningioma (81)	All_treatments	0.449	0.431	0.469	0.617	0.392
SJLIFE	Meningioma (81)	PRS	-0.119	-0.115	-0.124	-0.116	-0.12
SJLIFE	Meningioma (81)	Lifestyle	-0.13	-0.122	-0.14	-0.132	-0.129
SJLIFE	Meningioma (81)	Combined	0.298	0.288	0.31	0.503	0.229", header = T, sep = "\t")


data$SN_types_new <- gsub("\\([0-9]+\\)", "", data$SN_types)
data[data == "-"] <- NA


lifestyle <- "lifestyle_with_diet_HEI2015"
all.group <- c("Overall", "Female", "Male", "Age.lt.35", "Age.ge.35")
variables <-  c("Combined", "Radiation", "Chemotherapy", "treatments", "PRS", "Lifestyle")


# Define the text file path
data_name <- paste0("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/ANALYSIS/results/plots/v20b/", lifestyle, "_all_table_19b.txt")

file_conn <- file(data_name, "a")

for (j in 1:length(all.group)){
  group <-  all.group [j]
  for(i in 1:length(variables)){
    if(lifestyle== "lifestyle_with_diet_HEI2015" && variables[i] == "Lifestyle"){
      new_data <- data[grepl(variables[i], data$Variables), c("Cohort", "SN_types", "SN_types_new", "Variables", group)]
      
      AF.type <- variables[i]
      
      # Filter the data for SJLIFE and CCSS cohorts
      sjlife_data <- new_data %>% filter(Cohort == "SJLIFE")
      ccss_data <- new_data %>% filter(Cohort == "CCSS")
      
      # Create a new data frame with the desired format
      result_data <- tibble(
        "SN types (CA in SJLIFE/CA in CCSS)" = sjlife_data$SN_types,
        "SJLIFE" = as.numeric(sjlife_data[,group]))
      
      # Print the result_data data frame
      table_name <- paste0(group, "_", AF.type, "_", lifestyle,"\tNA\tNA")
      # write.table(table_name, data_name, col.names = T, row.names = F, quote = F, append = T, sep = "\t")
      # write.table(result_data, data_name, col.names = T, row.names = F, quote = F, append = T, sep = "\t")
      # }
      
      cat(table_name, "\n", file = file_conn)
      write.table(result_data, file = file_conn, col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")
    }
  }
}
# Close the text file
close(file_conn)
######################################
######################################
## create individual tables
## Lifestyle without diet
data <- read.table(text="Cohort	SNtypes	Variables	Overall_all4	Female_all4	Male_all4	Age_lt_35_all4	Age_ge_35_all4	Overall_smk	Female_smk	Male_smk	Age_lt_35_smk	Age_ge_35_smk	Overall_drk	Female_drk	Male_drk	Age_lt_35_drk	Age_ge_35_drk	Overall_PA	Female_PA	Male_PA	Age_lt_35_PA	Age_ge_35_PA	Overall_obese	Female_obese	Male_obese	Age_lt_35_obese	Age_ge_35_obese
SJLIFE	SNs	Radiation	0.442	0.432	0.455	0.368	0.473	0.442	0.432	0.455	0.368	0.473	0.442	0.432	0.455	0.368	0.473	0.442	0.432	0.455	0.368	0.473	0.442	0.432	0.455	0.368	0.473
SJLIFE	SNs	Chemotherapy	0.017	0.017	0.016	0.031	0.011	0.017	0.017	0.016	0.031	0.011	0.017	0.017	0.016	0.031	0.011	0.017	0.017	0.016	0.031	0.011	0.017	0.017	0.016	0.031	0.011
SJLIFE	SNs	All_treatments	0.45	0.442	0.461	0.385	0.477	0.45	0.442	0.461	0.385	0.477	0.45	0.442	0.461	0.385	0.477	0.45	0.442	0.461	0.385	0.477	0.45	0.442	0.461	0.385	0.477
SJLIFE	SNs	PRS	0.179	0.177	0.183	0.182	0.178	0.179	0.177	0.183	0.182	0.178	0.179	0.177	0.183	0.182	0.178	0.179	0.177	0.183	0.182	0.178	0.179	0.177	0.183	0.182	0.178
SJLIFE	SNs	Lifestyle	-0.063	-0.058	-0.069	-0.058	-0.065	-0.057	-0.059	-0.054	-0.052	-0.059	-0.036	-0.033	-0.039	-0.039	-0.034	0.01	0.014	0.006	0.005	0.012	-0.027	-0.028	-0.026	-0.028	-0.027
SJLIFE	SNs	Combined	0.52	0.514	0.527	0.462	0.543	0.523	0.515	0.534	0.466	0.546	0.533	0.526	0.542	0.475	0.556	0.554	0.548	0.562	0.498	0.577	0.538	0.53	0.548	0.481	0.56
SJLIFE	SMNs	Radiation	0.379	0.378	0.379	0.29	0.413	0.379	0.378	0.379	0.29	0.413	0.379	0.378	0.379	0.29	0.413	0.379	0.378	0.379	0.29	0.413	0.379	0.378	0.379	0.29	0.413
SJLIFE	SMNs	Chemotherapy	0.017	0.017	0.016	0.029	0.012	0.017	0.017	0.016	0.029	0.012	0.017	0.017	0.016	0.029	0.012	0.017	0.017	0.016	0.029	0.012	0.017	0.017	0.016	0.029	0.012
SJLIFE	SMNs	All_treatments	0.389	0.39	0.388	0.308	0.421	0.389	0.39	0.388	0.308	0.421	0.389	0.39	0.388	0.308	0.421	0.389	0.39	0.388	0.308	0.421	0.389	0.39	0.388	0.308	0.421
SJLIFE	SMNs	PRS	0.142	0.139	0.145	0.141	0.142	0.142	0.139	0.145	0.141	0.142	0.142	0.139	0.145	0.141	0.142	0.142	0.139	0.145	0.141	0.142	0.142	0.139	0.145	0.141	0.142
SJLIFE	SMNs	Lifestyle	-0.083	-0.073	-0.095	-0.081	-0.083	-0.079	-0.081	-0.076	-0.081	-0.077	-0.05	-0.05	-0.051	-0.058	-0.047	0.033	0.042	0.021	0.017	0.039	-0.11	-0.112	-0.107	-0.104	-0.112
SJLIFE	SMNs	Combined	0.431	0.437	0.423	0.353	0.46	0.434	0.433	0.435	0.354	0.465	0.449	0.449	0.45	0.369	0.48	0.493	0.497	0.488	0.414	0.523	0.418	0.417	0.418	0.342	0.446
SJLIFE	NMSC	Radiation	0.291	0.28	0.305	0.236	0.303	0.291	0.28	0.305	0.236	0.303	0.291	0.28	0.305	0.236	0.303	0.291	0.28	0.305	0.236	0.303	0.291	0.28	0.305	0.236	0.303
SJLIFE	NMSC	Chemotherapy	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-
SJLIFE	NMSC	All_treatments	0.291	0.28	0.305	0.236	0.303	0.291	0.28	0.305	0.236	0.303	0.291	0.28	0.305	0.236	0.303	0.291	0.28	0.305	0.236	0.303	0.291	0.28	0.305	0.236	0.303
SJLIFE	NMSC	PRS	0.311	0.293	0.334	0.308	0.311	0.311	0.293	0.334	0.308	0.311	0.311	0.293	0.334	0.308	0.311	0.311	0.293	0.334	0.308	0.311	0.311	0.293	0.334	0.308	0.311
SJLIFE	NMSC	Lifestyle	-0.275	-0.253	-0.304	-0.241	-0.282	-0.019	-0.021	-0.017	-0.036	-0.016	-0.152	-0.133	-0.175	-0.147	-0.152	-0.053	-0.057	-0.049	-0.067	-0.051	-0.157	-0.156	-0.157	-0.144	-0.159
SJLIFE	NMSC	Combined	0.376	0.364	0.393	0.336	0.385	0.502	0.482	0.527	0.451	0.513	0.376	0.364	0.393	0.336	0.385	0.485	0.463	0.512	0.434	0.496	0.434	0.413	0.461	0.392	0.443
SJLIFE	Breast_cancer	Radiation	0.466	0.466	-	0.392	0.478	0.466	0.466	-	0.392	0.478	0.466	0.466	-	0.392	0.478	0.466	0.466	-	0.392	0.478	0.466	0.466	-	0.392	0.478
SJLIFE	Breast_cancer	Chemotherapy	0.107	0.107	-	0.103	0.107	0.107	0.107	-	0.103	0.107	0.107	0.107	-	0.103	0.107	0.107	0.107	-	0.103	0.107	0.107	0.107	-	0.103	0.107
SJLIFE	Breast_cancer	All_treatments	0.527	0.527	-	0.454	0.538	0.527	0.527	-	0.454	0.538	0.527	0.527	-	0.454	0.538	0.527	0.527	-	0.454	0.538	0.527	0.527	-	0.454	0.538
SJLIFE	Breast_cancer	PRS	0.22	0.22	-	0.247	0.215	0.22	0.22	-	0.247	0.215	0.22	0.22	-	0.247	0.215	0.22	0.22	-	0.247	0.215	0.22	0.22	-	0.247	0.215
SJLIFE	Breast_cancer	Lifestyle	-0.882	-0.882	-	-0.596	-0.927	-0.26	-0.26	-	-0.166	-0.275	0.034	0.034	-	0.031	0.034	-0.379	-0.379	-	-0.228	-0.403	0.002	0.002	-	0.003	0.002
SJLIFE	Breast_cancer	Combined	0.331	0.331	-	0.333	0.331	0.551	0.551	-	0.523	0.556	0.648	0.648	-	0.607	0.655	0.51	0.51	-	0.492	0.513	0.634	0.634	-	0.596	0.64
SJLIFE	Thyroid_cancer	Radiation	0.61	0.634	0.583	0.554	0.638	0.61	0.634	0.583	0.554	0.638	0.61	0.634	0.583	0.554	0.638	0.61	0.634	0.583	0.554	0.638	0.61	0.634	0.583	0.554	0.638
SJLIFE	Thyroid_cancer	Chemotherapy	0.239	0.238	0.241	0.334	0.193	0.239	0.238	0.241	0.334	0.193	0.239	0.238	0.241	0.334	0.193	0.239	0.238	0.241	0.334	0.193	0.239	0.238	0.241	0.334	0.193
SJLIFE	Thyroid_cancer	All_treatments	0.725	0.741	0.706	0.709	0.733	0.725	0.741	0.706	0.709	0.733	0.725	0.741	0.706	0.709	0.733	0.725	0.741	0.706	0.709	0.733	0.725	0.741	0.706	0.709	0.733
SJLIFE	Thyroid_cancer	PRS	0.589	0.587	0.592	0.6	0.584	0.589	0.587	0.592	0.6	0.584	0.589	0.587	0.592	0.6	0.584	0.589	0.587	0.592	0.6	0.584	0.589	0.587	0.592	0.6	0.584
SJLIFE	Thyroid_cancer	Lifestyle	-0.175	-0.184	-0.164	-0.141	-0.191	0.095	0.101	0.089	0.09	0.098	-0.145	-0.13	-0.163	-0.148	-0.144	-0.182	-0.208	-0.152	-0.158	-0.194	-0.103	-0.104	-0.101	-0.093	-0.108
SJLIFE	Thyroid_cancer	Combined	0.871	0.876	0.864	0.863	0.874	0.9	0.905	0.894	0.891	0.904	0.871	0.879	0.862	0.861	0.876	0.868	0.872	0.864	0.86	0.872	0.877	0.883	0.87	0.869	0.881
SJLIFE	Meningioma	Radiation	0.185	0.15	0.226	0.246	0.164	0.185	0.15	0.226	0.246	0.164	0.185	0.15	0.226	0.246	0.164	0.185	0.15	0.226	0.246	0.164	0.185	0.15	0.226	0.246	0.164
SJLIFE	Meningioma	Chemotherapy	0.321	0.324	0.317	0.461	0.273	0.321	0.324	0.317	0.461	0.273	0.321	0.324	0.317	0.461	0.273	0.321	0.324	0.317	0.461	0.273	0.321	0.324	0.317	0.461	0.273
SJLIFE	Meningioma	All_treatments	0.456	0.436	0.481	0.623	0.4	0.456	0.436	0.481	0.623	0.4	0.456	0.436	0.481	0.623	0.4	0.456	0.436	0.481	0.623	0.4	0.456	0.436	0.481	0.623	0.4
SJLIFE	Meningioma	PRS	-0.121	-0.116	-0.126	-0.118	-0.122	-0.121	-0.116	-0.126	-0.118	-0.122	-0.121	-0.116	-0.126	-0.118	-0.122	-0.121	-0.116	-0.126	-0.118	-0.122	-0.121	-0.116	-0.126	-0.118	-0.122
SJLIFE	Meningioma	Lifestyle	-0.134	-0.125	-0.145	-0.136	-0.133	-0.051	-0.049	-0.054	-0.045	-0.054	-0.155	-0.126	-0.189	-0.144	-0.158	-0.088	-0.102	-0.072	-0.077	-0.092	0.151	0.145	0.159	0.126	0.159
SJLIFE	Meningioma	Combined	0.306	0.291	0.324	0.509	0.238	0.359	0.338	0.384	0.556	0.293	0.295	0.29	0.302	0.51	0.223	0.339	0.309	0.375	0.548	0.27	0.484	0.462	0.509	0.629	0.435
CCSS	SNs	Radiation	0.358	0.351	0.369	0.338	0.368	0.358	0.351	0.369	0.338	0.368	0.358	0.351	0.369	0.338	0.368	0.358	0.351	0.369	0.338	0.368	0.358	0.351	0.369	0.338	0.368
CCSS	SNs	Chemotherapy	0.038	0.037	0.041	0.061	0.027	0.038	0.037	0.041	0.061	0.027	0.038	0.037	0.041	0.061	0.027	0.038	0.037	0.041	0.061	0.027	0.038	0.037	0.041	0.061	0.027
CCSS	SNs	All_treatments	0.384	0.376	0.396	0.38	0.387	0.384	0.376	0.396	0.38	0.387	0.384	0.376	0.396	0.38	0.387	0.384	0.376	0.396	0.38	0.387	0.384	0.376	0.396	0.38	0.387
CCSS	SNs	PRS	0.057	0.056	0.058	0.056	0.057	0.057	0.056	0.058	0.056	0.057	0.057	0.056	0.058	0.056	0.057	0.057	0.056	0.058	0.056	0.057	0.057	0.056	0.058	0.056	0.057
CCSS	SNs	Lifestyle	-0.074	-0.071	-0.078	-0.068	-0.077	-0.052	-0.048	-0.057	-0.045	-0.055	-0.003	-0.002	-0.003	-0.002	-0.003	0.002	0.002	0.002	0.003	0.002	-0.017	-0.018	-0.015	-0.016	-0.017
CCSS	SNs	Combined	0.376	0.37	0.386	0.374	0.378	0.389	0.383	0.397	0.386	0.39	0.418	0.41	0.429	0.413	0.42	0.42	0.413	0.432	0.416	0.422	0.41	0.402	0.422	0.406	0.412
CCSS	SMNs	Radiation	0.249	0.245	0.255	0.186	0.275	0.249	0.245	0.255	0.186	0.275	0.249	0.245	0.255	0.186	0.275	0.249	0.245	0.255	0.186	0.275	0.249	0.245	0.255	0.186	0.275
CCSS	SMNs	Chemotherapy	0.056	0.055	0.058	0.099	0.038	0.056	0.055	0.058	0.099	0.038	0.056	0.055	0.058	0.099	0.038	0.056	0.055	0.058	0.099	0.038	0.056	0.055	0.058	0.099	0.038
CCSS	SMNs	All_treatments	0.297	0.292	0.306	0.271	0.308	0.297	0.292	0.306	0.271	0.308	0.297	0.292	0.306	0.271	0.308	0.297	0.292	0.306	0.271	0.308	0.297	0.292	0.306	0.271	0.308
CCSS	SMNs	PRS	0.048	0.047	0.049	0.047	0.048	0.048	0.047	0.049	0.047	0.048	0.048	0.047	0.049	0.047	0.048	0.048	0.047	0.049	0.047	0.048	0.048	0.047	0.049	0.047	0.048
CCSS	SMNs	Lifestyle	-0.048	-0.047	-0.052	-0.046	-0.049	-0.037	-0.035	-0.04	-0.034	-0.038	-0.004	-0.004	-0.005	-0.004	-0.004	-0.002	-0.002	-0.002	-0.002	-0.003	-0.002	-0.002	-0.002	-0.001	-0.002
CCSS	SMNs	Combined	0.298	0.295	0.303	0.273	0.308	0.306	0.303	0.311	0.281	0.316	0.328	0.323	0.336	0.303	0.338	0.329	0.324	0.338	0.304	0.339	0.329	0.325	0.338	0.305	0.34
CCSS	NMSC	Radiation	0.391	0.388	0.395	0.401	0.387	0.391	0.388	0.395	0.401	0.387	0.391	0.388	0.395	0.401	0.387	0.391	0.388	0.395	0.401	0.387	0.391	0.388	0.395	0.401	0.387
CCSS	NMSC	Chemotherapy	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-
CCSS	NMSC	All_treatments	0.391	0.388	0.395	0.401	0.387	0.391	0.388	0.395	0.401	0.387	0.391	0.388	0.395	0.401	0.387	0.391	0.388	0.395	0.401	0.387	0.391	0.388	0.395	0.401	0.387
CCSS	NMSC	PRS	0.308	0.306	0.311	0.309	0.308	0.308	0.306	0.311	0.309	0.308	0.308	0.306	0.311	0.309	0.308	0.308	0.306	0.311	0.309	0.308	0.308	0.306	0.311	0.309	0.308
CCSS	NMSC	Lifestyle	-0.155	-0.154	-0.156	-0.138	-0.162	-0.071	-0.061	-0.082	-0.059	-0.075	0.021	0.02	0.022	0.026	0.019	-0.04	-0.041	-0.038	-0.034	-0.042	-0.046	-0.05	-0.041	-0.041	-0.048
CCSS	NMSC	Combined	0.518	0.515	0.521	0.533	0.512	0.55	0.55	0.55	0.562	0.546	0.589	0.585	0.593	0.6	0.585	0.563	0.559	0.568	0.574	0.559	0.562	0.558	0.567	0.573	0.559
CCSS	Breast_cancer	Radiation	0.445	0.445	-	0.394	0.451	0.445	0.445	-	0.394	0.451	0.445	0.445	-	0.394	0.451	0.445	0.445	-	0.394	0.451	0.445	0.445	-	0.394	0.451
CCSS	Breast_cancer	Chemotherapy	0.192	0.192	-	0.232	0.188	0.192	0.192	-	0.232	0.188	0.192	0.192	-	0.232	0.188	0.192	0.192	-	0.232	0.188	0.192	0.192	-	0.232	0.188
CCSS	Breast_cancer	All_treatments	0.572	0.572	-	0.542	0.576	0.572	0.572	-	0.542	0.576	0.572	0.572	-	0.542	0.576	0.572	0.572	-	0.542	0.576	0.572	0.572	-	0.542	0.576
CCSS	Breast_cancer	PRS	0.346	0.346	-	0.348	0.346	0.346	0.346	-	0.348	0.346	0.346	0.346	-	0.348	0.346	0.346	0.346	-	0.348	0.346	0.346	0.346	-	0.348	0.346
CCSS	Breast_cancer	Lifestyle	-0.085	-0.085	-	-0.076	-0.086	-0.041	-0.041	-	-0.036	-0.041	0.018	0.018	-	0.025	0.017	-0.021	-0.021	-	-0.016	-0.021	-0.026	-0.026	-	-0.024	-0.026
CCSS	Breast_cancer	Combined	0.697	0.697	-	0.68	0.699	0.71	0.71	-	0.692	0.712	0.727	0.727	-	0.71	0.729	0.716	0.716	-	0.698	0.719	0.714	0.714	-	0.695	0.716
CCSS	Thyroid_cancer	Radiation	0.416	0.416	0.416	0.359	0.485	0.416	0.416	0.416	0.359	0.485	0.416	0.416	0.416	0.359	0.485	0.416	0.416	0.416	0.359	0.485	0.416	0.416	0.416	0.359	0.485
CCSS	Thyroid_cancer	Chemotherapy	0.095	0.094	0.097	0.122	0.061	0.095	0.094	0.097	0.122	0.061	0.095	0.094	0.097	0.122	0.061	0.095	0.094	0.097	0.122	0.061	0.095	0.094	0.097	0.122	0.061
CCSS	Thyroid_cancer	All_treatments	0.479	0.477	0.481	0.443	0.522	0.479	0.477	0.481	0.443	0.522	0.479	0.477	0.481	0.443	0.522	0.479	0.477	0.481	0.443	0.522	0.479	0.477	0.481	0.443	0.522
CCSS	Thyroid_cancer	PRS	0.385	0.39	0.375	0.378	0.394	0.385	0.39	0.375	0.378	0.394	0.385	0.39	0.375	0.378	0.394	0.385	0.39	0.375	0.378	0.394	0.385	0.39	0.375	0.378	0.394
CCSS	Thyroid_cancer	Lifestyle	-0.042	-0.035	-0.054	-0.025	-0.062	-0.036	-0.027	-0.053	-0.022	-0.053	0.061	0.06	0.062	0.073	0.045	0.02	0.021	0.019	0.03	0.008	0.005	0.005	0.004	0.014	-0.006
CCSS	Thyroid_cancer	Combined	0.666	0.667	0.664	0.645	0.691	0.667	0.67	0.663	0.646	0.693	0.7	0.699	0.7	0.681	0.723	0.686	0.686	0.686	0.666	0.711	0.682	0.681	0.682	0.661	0.707
CCSS	Meningioma	Radiation	0.297	0.264	0.348	0.319	0.282	0.297	0.264	0.348	0.319	0.282	0.297	0.264	0.348	0.319	0.282	0.297	0.264	0.348	0.319	0.282	0.297	0.264	0.348	0.319	0.282
CCSS	Meningioma	Chemotherapy	0.087	0.079	0.098	0.132	0.055	0.087	0.079	0.098	0.132	0.055	0.087	0.079	0.098	0.132	0.055	0.087	0.079	0.098	0.132	0.055	0.087	0.079	0.098	0.132	0.055
CCSS	Meningioma	All_treatments	0.371	0.337	0.422	0.424	0.334	0.371	0.337	0.422	0.424	0.334	0.371	0.337	0.422	0.424	0.334	0.371	0.337	0.422	0.424	0.334	0.371	0.337	0.422	0.424	0.334
CCSS	Meningioma	PRS	0.023	0.022	0.024	0.028	0.02	0.023	0.022	0.024	0.028	0.02	0.023	0.022	0.024	0.028	0.02	0.023	0.022	0.024	0.028	0.02	0.023	0.022	0.024	0.028	0.02
CCSS	Meningioma	Lifestyle	0.072	0.083	0.054	0.133	0.029	0.089	0.095	0.079	0.149	0.046	0.051	0.059	0.037	0.111	0.008	0.109	0.115	0.1	0.168	0.069	0.122	0.129	0.111	0.181	0.08
CCSS	Meningioma	Combined	0.423	0.398	0.46	0.513	0.359	0.439	0.412	0.481	0.528	0.377	0.413	0.385	0.453	0.504	0.348	0.453	0.425	0.494	0.539	0.392	0.459	0.432	0.499	0.544	0.399
CCSS	Sarcoma	Radiation	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-
CCSS	Sarcoma	Chemotherapy	0.227	0.217	0.243	0.243	0.217	0.227	0.217	0.243	0.243	0.217	0.227	0.217	0.243	0.243	0.217	0.227	0.217	0.243	0.243	0.217	0.227	0.217	0.243	0.243	0.217
CCSS	Sarcoma	All_treatments	0.227	0.217	0.243	0.243	0.217	0.227	0.217	0.243	0.243	0.217	0.227	0.217	0.243	0.243	0.217	0.227	0.217	0.243	0.243	0.217	0.227	0.217	0.243	0.243	0.217
CCSS	Sarcoma	PRS	0.05	0.053	0.045	0.048	0.051	0.05	0.053	0.045	0.048	0.051	0.05	0.053	0.045	0.048	0.051	0.05	0.053	0.045	0.048	0.051	0.05	0.053	0.045	0.048	0.051
CCSS	Sarcoma	Lifestyle	-0.278	-0.273	-0.286	-0.426	-0.189	-0.134	-0.139	-0.124	-0.24	-0.07	-0.24	-0.238	-0.245	-0.385	-0.153	-0.101	-0.107	-0.089	-0.205	-0.038	-0.142	-0.149	-0.128	-0.25	-0.076
CCSS	Sarcoma	Combined	0.062	0.058	0.068	-0.03	0.117	0.166	0.155	0.184	0.104	0.203	0.09	0.085	0.098	0.001	0.144	0.062	0.058	0.068	-0.03	0.117	0.161	0.148	0.182	0.097	0.199", header = T, sep = "\t")

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
## Lifestyle with diet HEI (SJLIFE only)
data <- read.table(text="SNtypes	Variables	Overall_all4	Female_all4	Male_all4	Age_lt_35_all4	Age_ge_35_all4	Overall_smk	Female_smk	Male_smk	Age_lt_35_smk	Age_ge_35_smk	Overall_drk	Female_drk	Male_drk	Age_lt_35_drk	Age_ge_35_drk	Overall_PA	Female_PA	Male_PA	Age_lt_35_PA	Age_ge_35_PA	Overall_obese	Female_obese	Male_obese	Age_lt_35_obese	Age_ge_35_obese	Overall_diet	Female_diet	Male_diet	Age_lt_35_diet	Age_ge_35_diet
SNs	Radiation	0.443	0.433	0.457	0.369	0.473	0.443	0.433	0.457	0.369	0.473	0.443	0.433	0.457	0.369	0.473	0.443	0.433	0.457	0.369	0.473	0.443	0.433	0.457	0.369	0.473	0.443	0.433	0.457	0.369	0.473
SNs	Chemotherapy	0.02	0.02	0.019	0.036	0.013	0.02	0.02	0.019	0.036	0.013	0.02	0.02	0.019	0.036	0.013	0.02	0.02	0.019	0.036	0.013	0.02	0.02	0.019	0.036	0.013	0.02	0.02	0.019	0.036	0.013
SNs	All_treatments	0.453	0.443	0.465	0.389	0.479	0.453	0.443	0.465	0.389	0.479	0.453	0.443	0.465	0.389	0.479	0.453	0.443	0.465	0.389	0.479	0.453	0.443	0.465	0.389	0.479	0.453	0.443	0.465	0.389	0.479
SNs	PRS	0.18	0.178	0.183	0.183	0.179	0.18	0.178	0.183	0.183	0.179	0.18	0.178	0.183	0.183	0.179	0.18	0.178	0.183	0.183	0.179	0.18	0.178	0.183	0.183	0.179	0.18	0.178	0.183	0.183	0.179
SNs	Lifestyle	-0.057	-0.062	-0.049	-0.05	-0.06	-0.059	-0.061	-0.055	-0.052	-0.062	-0.03	-0.028	-0.033	-0.032	-0.029	0.014	0.017	0.01	0.01	0.016	-0.02	-0.021	-0.02	-0.021	-0.02	-0.012	-0.022	0.001	-0.012	-0.012
SNs	Combined	0.525	0.515	0.538	0.47	0.548	0.525	0.516	0.537	0.471	0.547	0.538	0.531	0.548	0.483	0.561	0.558	0.552	0.567	0.504	0.581	0.543	0.535	0.554	0.489	0.565	0.547	0.535	0.563	0.493	0.569
SMNs	Radiation	0.378	0.378	0.379	0.289	0.413	0.378	0.378	0.379	0.289	0.413	0.378	0.378	0.379	0.289	0.413	0.378	0.378	0.379	0.289	0.413	0.378	0.378	0.379	0.289	0.413	0.378	0.378	0.379	0.289	0.413
SMNs	Chemotherapy	0.018	0.018	0.018	0.031	0.013	0.018	0.018	0.018	0.031	0.013	0.018	0.018	0.018	0.031	0.013	0.018	0.018	0.018	0.031	0.013	0.018	0.018	0.018	0.031	0.013	0.018	0.018	0.018	0.031	0.013
SMNs	All_treatments	0.39	0.39	0.389	0.309	0.421	0.39	0.39	0.389	0.309	0.421	0.39	0.39	0.389	0.309	0.421	0.39	0.39	0.389	0.309	0.421	0.39	0.39	0.389	0.309	0.421	0.39	0.39	0.389	0.309	0.421
SMNs	PRS	0.14	0.137	0.143	0.139	0.14	0.14	0.137	0.143	0.139	0.14	0.14	0.137	0.143	0.139	0.14	0.14	0.137	0.143	0.139	0.14	0.14	0.137	0.143	0.139	0.14	0.14	0.137	0.143	0.139	0.14
SMNs	Lifestyle	-0.13	-0.122	-0.139	-0.133	-0.128	-0.102	-0.106	-0.097	-0.107	-0.1	-0.074	-0.075	-0.073	-0.085	-0.07	0.014	0.022	0.005	-0.004	0.021	-0.13	-0.134	-0.125	-0.126	-0.131	-0.082	-0.087	-0.076	-0.096	-0.077
SMNs	Combined	0.404	0.409	0.398	0.321	0.436	0.421	0.419	0.423	0.339	0.452	0.435	0.434	0.437	0.353	0.467	0.482	0.485	0.478	0.401	0.513	0.406	0.405	0.407	0.329	0.436	0.431	0.428	0.435	0.347	0.464
NMSC	Radiation	0.285	0.274	0.299	0.23	0.297	0.285	0.274	0.299	0.23	0.297	0.285	0.274	0.299	0.23	0.297	0.285	0.274	0.299	0.23	0.297	0.285	0.274	0.299	0.23	0.297	0.285	0.274	0.299	0.23	0.297
NMSC	Chemotherapy	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-
NMSC	All_treatments	0.285	0.274	0.299	0.23	0.297	0.285	0.274	0.299	0.23	0.297	0.285	0.274	0.299	0.23	0.297	0.285	0.274	0.299	0.23	0.297	0.285	0.274	0.299	0.23	0.297	0.285	0.274	0.299	0.23	0.297
NMSC	PRS	0.328	0.311	0.349	0.326	0.328	0.328	0.311	0.349	0.326	0.328	0.328	0.311	0.349	0.326	0.328	0.328	0.311	0.349	0.326	0.328	0.328	0.311	0.349	0.326	0.328	0.328	0.311	0.349	0.326	0.328
NMSC	Lifestyle	-0.549	-0.459	-0.665	-0.544	-0.55	-0.04	-0.045	-0.035	-0.07	-0.034	-0.197	-0.179	-0.219	-0.204	-0.195	-0.065	-0.069	-0.059	-0.093	-0.059	-0.184	-0.186	-0.181	-0.188	-0.183	-0.293	-0.249	-0.35	-0.341	-0.283
NMSC	Combined	0.255	0.269	0.238	0.191	0.269	0.499	0.479	0.525	0.445	0.511	0.423	0.411	0.439	0.367	0.435	0.487	0.466	0.513	0.431	0.499	0.429	0.408	0.456	0.382	0.439	0.378	0.374	0.382	0.302	0.394
Breast_cancer	Radiation	0.47	0.47	-	0.4	0.481	0.47	0.47	-	0.4	0.481	0.47	0.47	-	0.4	0.481	0.47	0.47	-	0.4	0.481	0.47	0.47	-	0.4	0.481	0.47	0.47	-	0.4	0.481
Breast_cancer	Chemotherapy	0.109	0.109	-	0.106	0.11	0.109	0.109	-	0.106	0.11	0.109	0.109	-	0.106	0.11	0.109	0.109	-	0.106	0.11	0.109	0.109	-	0.106	0.11	0.109	0.109	-	0.106	0.11
Breast_cancer	All_treatments	0.532	0.532	-	0.462	0.543	0.532	0.532	-	0.462	0.543	0.532	0.532	-	0.462	0.543	0.532	0.532	-	0.462	0.543	0.532	0.532	-	0.462	0.543	0.532	0.532	-	0.462	0.543
Breast_cancer	PRS	0.226	0.226	-	0.259	0.221	0.226	0.226	-	0.259	0.221	0.226	0.226	-	0.259	0.221	0.226	0.226	-	0.259	0.221	0.226	0.226	-	0.259	0.221	0.226	0.226	-	0.259	0.221
Breast_cancer	Lifestyle	-0.917	-0.917	-	-0.602	-0.968	-0.294	-0.294	-	-0.194	-0.31	0.001	0.001	-	0.001	0.001	-0.417	-0.417	-	-0.259	-0.442	-0.044	-0.044	-	-0.04	-0.045	0.041	0.041	-	0.051	0.04
Breast_cancer	Combined	0.342	0.342	-	0.35	0.341	0.551	0.551	-	0.525	0.555	0.645	0.645	-	0.607	0.651	0.509	0.509	-	0.496	0.511	0.626	0.626	-	0.592	0.631	0.663	0.663	-	0.626	0.669
Thyroid_cancer	Radiation	0.611	0.635	0.584	0.55	0.641	0.611	0.635	0.584	0.55	0.641	0.611	0.635	0.584	0.55	0.641	0.611	0.635	0.584	0.55	0.641	0.611	0.635	0.584	0.55	0.641	0.611	0.635	0.584	0.55	0.641
Thyroid_cancer	Chemotherapy	0.24	0.238	0.243	0.329	0.197	0.24	0.238	0.243	0.329	0.197	0.24	0.238	0.243	0.329	0.197	0.24	0.238	0.243	0.329	0.197	0.24	0.238	0.243	0.329	0.197	0.24	0.238	0.243	0.329	0.197
Thyroid_cancer	All_treatments	0.727	0.743	0.708	0.704	0.738	0.727	0.743	0.708	0.704	0.738	0.727	0.743	0.708	0.704	0.738	0.727	0.743	0.708	0.704	0.738	0.727	0.743	0.708	0.704	0.738	0.727	0.743	0.708	0.704	0.738
Thyroid_cancer	PRS	0.582	0.579	0.585	0.593	0.576	0.582	0.579	0.585	0.593	0.576	0.582	0.579	0.585	0.593	0.576	0.582	0.579	0.585	0.593	0.576	0.582	0.579	0.585	0.593	0.576	0.582	0.579	0.585	0.593	0.576
Thyroid_cancer	Lifestyle	-0.677	-0.561	-0.811	-0.709	-0.662	0.066	0.07	0.061	0.055	0.071	-0.276	-0.245	-0.312	-0.282	-0.273	-0.244	-0.267	-0.218	-0.231	-0.251	-0.188	-0.189	-0.186	-0.186	-0.189	-0.515	-0.413	-0.633	-0.593	-0.477
Thyroid_cancer	Combined	0.811	0.834	0.785	0.792	0.82	0.896	0.901	0.89	0.884	0.902	0.855	0.865	0.843	0.84	0.862	0.86	0.865	0.854	0.848	0.866	0.866	0.873	0.858	0.854	0.873	0.826	0.847	0.803	0.805	0.837
Meningioma	Radiation	0.182	0.15	0.22	0.245	0.161	0.185	0.15	0.226	0.245	0.165	0.185	0.15	0.226	0.245	0.165	0.185	0.15	0.226	0.245	0.165	0.185	0.15	0.226	0.245	0.165	0.185	0.15	0.226	0.245	0.165
Meningioma	Chemotherapy	0.319	0.325	0.312	0.463	0.271	0.321	0.324	0.317	0.462	0.274	0.321	0.324	0.317	0.462	0.274	0.321	0.324	0.317	0.462	0.274	0.321	0.324	0.317	0.462	0.274	0.321	0.324	0.317	0.462	0.274
Meningioma	All_treatments	0.449	0.431	0.469	0.617	0.392	0.456	0.436	0.481	0.623	0.4	0.456	0.436	0.481	0.623	0.4	0.456	0.436	0.481	0.623	0.4	0.456	0.436	0.481	0.623	0.4	0.456	0.436	0.481	0.623	0.4
Meningioma	PRS	-0.119	-0.115	-0.124	-0.116	-0.12	-0.12	-0.116	-0.126	-0.117	-0.122	-0.12	-0.116	-0.126	-0.117	-0.122	-0.12	-0.116	-0.126	-0.117	-0.122	-0.12	-0.116	-0.126	-0.117	-0.122	-0.12	-0.116	-0.126	-0.117	-0.122
Meningioma	Lifestyle	-0.13	-0.122	-0.14	-0.132	-0.129	-0.038	-0.036	-0.041	-0.027	-0.042	-0.139	-0.111	-0.173	-0.124	-0.144	-0.073	-0.086	-0.057	-0.058	-0.078	0.166	0.16	0.173	0.144	0.173	0.006	0.007	0.004	0.011	0.004
Meningioma	Combined	0.298	0.288	0.31	0.503	0.229	0.367	0.346	0.391	0.564	0.3	0.304	0.299	0.311	0.519	0.232	0.348	0.319	0.383	0.556	0.278	0.492	0.471	0.517	0.637	0.444	0.395	0.375	0.419	0.583	0.332", header = T, sep = "\t")

## combined
data.sjlife <- data

data.sjlife.wanted <- data.sjlife[grepl("Lifestyle", data.sjlife$Variables),]

variables <- c("diet", "drk", "smk", "PA", "obese")


i=1
sjlife <- data.sjlife.wanted[grepl(paste0("SNtypes|",variables[i]), colnames(data.sjlife.wanted))]
View(sjlife)


###########################################################
###########################################################
## Lifestyle with diet HEI with HEI score less than 60 binary (SJLIFE only)
data <- read.table(text="SNtypes	Variables	Overall_all4	Female_all4	Male_all4	Age_lt_35_all4	Age_ge_35_all4	Overall_smk	Female_smk	Male_smk	Age_lt_35_smk	Age_ge_35_smk	Overall_drk	Female_drk	Male_drk	Age_lt_35_drk	Age_ge_35_drk	Overall_PA	Female_PA	Male_PA	Age_lt_35_PA	Age_ge_35_PA	Overall_obese	Female_obese	Male_obese	Age_lt_35_obese	Age_ge_35_obese	Overall_diet	Female_diet	Male_diet	Age_lt_35_diet	Age_ge_35_diet
SNs	Radiation	0.442	0.432	0.455	0.368	0.473	0.442	0.432	0.455	0.368	0.473	0.442	0.432	0.455	0.368	0.473	0.442	0.432	0.455	0.368	0.473	0.442	0.432	0.455	0.368	0.473	0.442	0.432	0.455	0.368	0.473
SNs	Chemotherapy	0.017	0.017	0.016	0.031	0.011	0.017	0.017	0.016	0.031	0.011	0.017	0.017	0.016	0.031	0.011	0.017	0.017	0.016	0.031	0.011	0.017	0.017	0.016	0.031	0.011	0.017	0.017	0.016	0.031	0.011
SNs	All_treatments	0.45	0.442	0.461	0.385	0.477	0.45	0.442	0.461	0.385	0.477	0.45	0.442	0.461	0.385	0.477	0.45	0.442	0.461	0.385	0.477	0.45	0.442	0.461	0.385	0.477	0.45	0.442	0.461	0.385	0.477
SNs	PRS	0.18	0.177	0.183	0.183	0.178	0.18	0.177	0.183	0.183	0.178	0.18	0.177	0.183	0.183	0.178	0.18	0.177	0.183	0.183	0.178	0.18	0.177	0.183	0.183	0.178	0.18	0.177	0.183	0.183	0.178
SNs	Lifestyle	-0.054	-0.05	-0.059	-0.049	-0.056	-0.053	-0.055	-0.05	-0.047	-0.056	-0.03	-0.028	-0.033	-0.032	-0.029	0.014	0.017	0.01	0.01	0.016	-0.021	-0.022	-0.02	-0.021	-0.021	-0.01	-0.011	-0.009	-0.012	-0.009
SNs	Combined	0.524	0.518	0.531	0.467	0.547	0.525	0.517	0.536	0.469	0.548	0.536	0.529	0.545	0.479	0.559	0.556	0.55	0.564	0.5	0.579	0.541	0.533	0.551	0.485	0.563	0.545	0.537	0.556	0.49	0.568
SMNs	Radiation	0.38	0.383	0.376	0.274	0.423	0.38	0.383	0.376	0.274	0.423	0.38	0.383	0.376	0.274	0.423	0.38	0.383	0.376	0.274	0.423	0.38	0.383	0.376	0.274	0.423	0.38	0.383	0.376	0.274	0.423
SMNs	Chemotherapy	0.062	0.061	0.063	0.106	0.044	0.062	0.061	0.063	0.106	0.044	0.062	0.061	0.063	0.106	0.044	0.062	0.061	0.063	0.106	0.044	0.062	0.061	0.063	0.106	0.044	0.062	0.061	0.063	0.106	0.044
SMNs	All_treatments	0.424	0.429	0.417	0.349	0.454	0.424	0.429	0.417	0.349	0.454	0.424	0.429	0.417	0.349	0.454	0.424	0.429	0.417	0.349	0.454	0.424	0.429	0.417	0.349	0.454	0.424	0.429	0.417	0.349	0.454
SMNs	PRS	0.218	0.213	0.223	0.215	0.219	0.218	0.213	0.223	0.215	0.219	0.218	0.213	0.223	0.215	0.219	0.218	0.213	0.223	0.215	0.219	0.218	0.213	0.223	0.215	0.219	0.218	0.213	0.223	0.215	0.219
SMNs	Lifestyle	0.028	0.007	0.053	0.032	0.026	-0.085	-0.088	-0.082	-0.089	-0.084	-0.004	-0.012	0.006	-0.009	-0.001	-0.027	-0.025	-0.028	-0.037	-0.022	-0.058	-0.06	-0.055	-0.063	-0.055	0.014	-0.003	0.036	0.015	0.014
SMNs	Combined	0.561	0.556	0.568	0.503	0.585	0.509	0.511	0.507	0.44	0.537	0.548	0.547	0.549	0.483	0.574	0.536	0.538	0.533	0.468	0.564	0.522	0.523	0.52	0.454	0.549	0.555	0.55	0.561	0.494	0.579
NMSCs	Radiation	0.293	0.282	0.306	0.236	0.305	0.293	0.282	0.306	0.236	0.305	0.293	0.282	0.306	0.236	0.305	0.293	0.282	0.306	0.236	0.305	0.293	0.282	0.306	0.236	0.305	0.293	0.282	0.306	0.236	0.305
NMSCs	Chemotherapy	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-
NMSCs	All_treatments	0.293	0.282	0.306	0.236	0.305	0.293	0.282	0.306	0.236	0.305	0.293	0.282	0.306	0.236	0.305	0.293	0.282	0.306	0.236	0.305	0.293	0.282	0.306	0.236	0.305	0.293	0.282	0.306	0.236	0.305
NMSCs	PRS	0.321	0.302	0.344	0.319	0.321	0.321	0.302	0.344	0.319	0.321	0.321	0.302	0.344	0.319	0.321	0.321	0.302	0.344	0.319	0.321	0.321	0.302	0.344	0.319	0.321	0.321	0.302	0.344	0.319	0.321
NMSCs	Lifestyle	-0.529	-0.425	-0.662	-0.527	-0.529	-0.034	-0.037	-0.029	-0.063	-0.027	-0.2	-0.181	-0.226	-0.207	-0.199	-0.056	-0.059	-0.052	-0.086	-0.049	-0.185	-0.186	-0.182	-0.19	-0.183	-0.289	-0.234	-0.359	-0.338	-0.278
NMSCs	Combined	0.268	0.288	0.242	0.198	0.283	0.503	0.482	0.529	0.447	0.515	0.422	0.409	0.438	0.363	0.434	0.491	0.47	0.518	0.433	0.504	0.429	0.407	0.457	0.379	0.44	0.382	0.384	0.38	0.303	0.4
Breast_cancer	Radiation	0.471	0.471	-	0.401	0.482	0.471	0.471	-	0.401	0.482	0.471	0.471	-	0.401	0.482	0.471	0.471	-	0.401	0.482	0.471	0.471	-	0.401	0.482	0.471	0.471	-	0.401	0.482
Breast_cancer	Chemotherapy	0.098	0.098	-	0.09	0.1	0.098	0.098	-	0.09	0.1	0.098	0.098	-	0.09	0.1	0.098	0.098	-	0.09	0.1	0.098	0.098	-	0.09	0.1	0.098	0.098	-	0.09	0.1
Breast_cancer	All_treatments	0.527	0.527	-	0.454	0.538	0.527	0.527	-	0.454	0.538	0.527	0.527	-	0.454	0.538	0.527	0.527	-	0.454	0.538	0.527	0.527	-	0.454	0.538	0.527	0.527	-	0.454	0.538
Breast_cancer	PRS	0.23	0.23	-	0.258	0.226	0.23	0.23	-	0.258	0.226	0.23	0.23	-	0.258	0.226	0.23	0.23	-	0.258	0.226	0.23	0.23	-	0.258	0.226	0.23	0.23	-	0.258	0.226
Breast_cancer	Lifestyle	-0.84	-0.84	-	-0.523	-0.891	-0.333	-0.333	-	-0.221	-0.352	0.004	0.004	-	0.003	0.004	-0.465	-0.465	-	-0.294	-0.493	-0.04	-0.04	-	-0.037	-0.041	0.091	0.091	-	0.103	0.089
Breast_cancer	Combined	0.362	0.362	-	0.377	0.36	0.537	0.537	-	0.51	0.541	0.645	0.645	-	0.602	0.651	0.491	0.491	-	0.472	0.494	0.626	0.626	-	0.587	0.632	0.676	0.676	-	0.644	0.681
Thyroid_cancer	Radiation	0.61	0.635	0.581	0.55	0.64	0.61	0.635	0.581	0.55	0.64	0.61	0.635	0.581	0.55	0.64	0.61	0.635	0.581	0.55	0.64	0.61	0.635	0.581	0.55	0.64	0.61	0.635	0.581	0.55	0.64
Thyroid_cancer	Chemotherapy	0.238	0.236	0.241	0.332	0.192	0.238	0.236	0.241	0.332	0.192	0.238	0.236	0.241	0.332	0.192	0.238	0.236	0.241	0.332	0.192	0.238	0.236	0.241	0.332	0.192	0.238	0.236	0.241	0.332	0.192
Thyroid_cancer	All_treatments	0.724	0.742	0.703	0.703	0.734	0.724	0.742	0.703	0.703	0.734	0.724	0.742	0.703	0.703	0.734	0.724	0.742	0.703	0.703	0.734	0.724	0.742	0.703	0.703	0.734	0.724	0.742	0.703	0.703	0.734
Thyroid_cancer	PRS	0.584	0.582	0.587	0.595	0.579	0.584	0.582	0.587	0.595	0.579	0.584	0.582	0.587	0.595	0.579	0.584	0.582	0.587	0.595	0.579	0.584	0.582	0.587	0.595	0.579	0.584	0.582	0.587	0.595	0.579
Thyroid_cancer	Lifestyle	-0.493	-0.421	-0.576	-0.508	-0.485	0.07	0.077	0.062	0.059	0.075	-0.228	-0.207	-0.252	-0.231	-0.226	-0.216	-0.237	-0.191	-0.203	-0.222	-0.156	-0.158	-0.153	-0.152	-0.158	-0.362	-0.296	-0.439	-0.415	-0.336
Thyroid_cancer	Combined	0.831	0.848	0.812	0.815	0.839	0.896	0.902	0.889	0.885	0.901	0.859	0.869	0.848	0.846	0.866	0.862	0.868	0.856	0.851	0.868	0.869	0.877	0.861	0.858	0.875	0.843	0.859	0.826	0.825	0.852
Meningioma	Radiation	0.182	0.15	0.22	0.245	0.161	0.185	0.15	0.226	0.245	0.164	0.185	0.15	0.226	0.245	0.164	0.185	0.15	0.226	0.245	0.164	0.185	0.15	0.226	0.245	0.164	0.185	0.15	0.226	0.245	0.164
Meningioma	Chemotherapy	0.318	0.325	0.31	0.461	0.27	0.32	0.324	0.316	0.461	0.273	0.32	0.324	0.316	0.461	0.273	0.32	0.324	0.316	0.461	0.273	0.32	0.324	0.316	0.461	0.273	0.32	0.324	0.316	0.461	0.273
Meningioma	All_treatments	0.448	0.431	0.468	0.616	0.391	0.456	0.435	0.48	0.622	0.4	0.456	0.435	0.48	0.622	0.4	0.456	0.435	0.48	0.622	0.4	0.456	0.435	0.48	0.622	0.4	0.456	0.435	0.48	0.622	0.4
Meningioma	PRS	-0.12	-0.115	-0.125	-0.116	-0.121	-0.121	-0.116	-0.126	-0.117	-0.122	-0.121	-0.116	-0.126	-0.117	-0.122	-0.121	-0.116	-0.126	-0.117	-0.122	-0.121	-0.116	-0.126	-0.117	-0.122	-0.121	-0.116	-0.126	-0.117	-0.122
Meningioma	Lifestyle	-0.109	-0.103	-0.116	-0.11	-0.108	-0.038	-0.036	-0.041	-0.027	-0.042	-0.138	-0.11	-0.172	-0.123	-0.144	-0.074	-0.087	-0.058	-0.059	-0.079	-0.074	-0.087	-0.058	-0.059	-0.079	-0.138	-0.11	-0.172	-0.123	-0.144
Meningioma	Combined	0.309	0.297	0.322	0.511	0.241	0.366	0.345	0.39	0.563	0.299	0.304	0.298	0.31	0.518	0.232	0.346	0.317	0.381	0.554	0.277	0.346	0.317	0.381	0.554	0.277	0.304	0.298	0.31	0.518	0.232", header = T, sep = "\t")

## combined
data.sjlife <- data

data.sjlife.wanted <- data.sjlife[grepl("Lifestyle", data.sjlife$Variables),]

variables <- c("diet", "drk", "smk", "PA", "obese")


i=1
sjlife <- data.sjlife.wanted[grepl(paste0("SNtypes|",variables[i]), colnames(data.sjlife.wanted))]
View(sjlife)
