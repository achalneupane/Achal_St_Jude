# Load required libraries
library(ggplot2)
library(ggthemes)  # For additional themes
library(tidyverse)

## V18 b
data <- read.table(text="Cohort	SN_types	Variables	Overall	Female	Male	Age.lt.35	Age.ge.35
SJLIFE	Any SN (605)	Radiation	0.425	0.417	0.435	0.400	0.447
SJLIFE	Any SN (605)	Chemo	0.080	0.079	0.080	0.109	0.054
SJLIFE	Any SN (605)	All_treatments	0.469	0.462	0.478	0.463	0.474
SJLIFE	Any SN (605)	PRS	0.125	0.123	0.129	0.126	0.125
SJLIFE	Any SN (605)	Lifestyle	-	-	-	-	-
SJLIFE	Any SN (605)	Combined	0.536	0.530	0.545	0.531	0.541
SJLIFE	Any SMN (462)	Radiation	0.372	0.371	0.373	0.340	0.396
SJLIFE	Any SMN (462)	Chemo	0.033	0.033	0.032	0.047	0.022
SJLIFE	Any SMN (462)	All_treatments	0.392	0.392	0.392	0.370	0.409
SJLIFE	Any SMN (462)	PRS	0.122	0.119	0.125	0.122	0.122
SJLIFE	Any SMN (462)	Lifestyle	-	-	-	-	-
SJLIFE	Any SMN (462)	Combined	0.467	0.466	0.468	0.447	0.482
SJLIFE	NMSC (249)	Radiation	0.441	0.431	0.453	0.413	0.457
SJLIFE	NMSC (249)	Chemo	-	-	-	-	-
SJLIFE	NMSC (249)	All_treatments	0.441	0.431	0.453	0.413	0.457
SJLIFE	NMSC (249)	PRS	0.441	0.432	0.452	0.440	0.442
SJLIFE	NMSC (249)	Lifestyle	-	-	-	-	-
SJLIFE	NMSC (249)	Combined	0.688	0.680	0.698	0.671	0.698
SJLIFE	Breast cancer (74)	Radiation	0.493	0.493	-	0.475	0.498
SJLIFE	Breast cancer (74)	Chemo	0.193	0.193	-	0.228	0.184
SJLIFE	Breast cancer (74)	All_treatments	0.602	0.602	-	0.594	0.604
SJLIFE	Breast cancer (74)	PRS	0.255	0.255	-	0.261	0.254
SJLIFE	Breast cancer (74)	Lifestyle	-	-	-	-	-
SJLIFE	Breast cancer (74)	Combined	0.704	0.704	-	0.701	0.705
SJLIFE	Thyroid cancer (86)	Radiation	0.620	0.622	0.616	0.594	0.650
SJLIFE	Thyroid cancer (86)	Chemo	0.233	0.235	0.231	0.293	0.163
SJLIFE	Thyroid cancer (86)	All_treatments	0.726	0.727	0.723	0.725	0.727
SJLIFE	Thyroid cancer (86)	PRS	0.517	0.519	0.514	0.519	0.514
SJLIFE	Thyroid cancer (86)	Lifestyle	-	-	-	-	-
SJLIFE	Thyroid cancer (86)	Combined	0.866	0.867	0.866	0.866	0.868
SJLIFE	Meningioma (149)	Radiation	0.197	0.172	0.230	0.223	0.176
SJLIFE	Meningioma (149)	Chemo	0.361	0.371	0.349	0.472	0.271
SJLIFE	Meningioma (149)	All_treatments	0.489	0.478	0.504	0.598	0.400
SJLIFE	Meningioma (149)	PRS	-0.125	-0.117	-0.135	-0.121	-0.128
SJLIFE	Meningioma (149)	Lifestyle	-	-	-	-	-
SJLIFE	Meningioma (149)	Combined	0.425	0.416	0.437	0.549	0.323
SJLIFE	Sarcoma (32)	Radiation	-	-	-	-	-
SJLIFE	Sarcoma (32)	Chemo	0.324	0.329	0.320	0.308	0.351
SJLIFE	Sarcoma (32)	All_treatments	0.324	0.329	0.320	0.308	0.351
SJLIFE	Sarcoma (32)	PRS	-0.020	-0.021	-0.020	-0.020	-0.020
SJLIFE	Sarcoma (32)	Lifestyle	-	-	-	-	-
SJLIFE	Sarcoma (32)	Combined	0.311	0.315	0.306	0.294	0.337
CCSS	Any SN (1611)	Radiation	0.387	0.380	0.398	0.372	0.398
CCSS	Any SN (1611)	Chemo	0.031	0.030	0.033	0.047	0.019
CCSS	Any SN (1611)	All_treatments	0.407	0.400	0.418	0.403	0.410
CCSS	Any SN (1611)	PRS	0.048	0.047	0.049	0.047	0.048
CCSS	Any SN (1611)	Lifestyle	-	-	-	-	-
CCSS	Any SN (1611)	Combined	0.435	0.428	0.446	0.431	0.438
CCSS	Any SMN (762)	Radiation	0.255	0.254	0.257	0.213	0.283
CCSS	Any SMN (762)	Chemo	0.042	0.042	0.042	0.070	0.024
CCSS	Any SMN (762)	All_treatments	0.290	0.289	0.293	0.271	0.303
CCSS	Any SMN (762)	PRS	0.047	0.046	0.047	0.046	0.047
CCSS	Any SMN (762)	Lifestyle	-	-	-	-	-
CCSS	Any SMN (762)	Combined	0.323	0.322	0.326	0.305	0.336
CCSS	NMSC (769)	Radiation	0.391	0.386	0.396	0.385	0.394
CCSS	NMSC (769)	Chemo	-	-	-	-	-
CCSS	NMSC (769)	All_treatments	0.391	0.386	0.396	0.385	0.394
CCSS	NMSC (769)	PRS	0.304	0.306	0.302	0.306	0.303
CCSS	NMSC (769)	Lifestyle	-	-	-	-	-
CCSS	NMSC (769)	Combined	0.577	0.575	0.578	0.574	0.578
CCSS	Breast cancer (289)	Radiation	0.474	0.474	-	0.452	0.478
CCSS	Breast cancer (289)	Chemo	0.187	0.187	-	0.224	0.180
CCSS	Breast cancer (289)	All_treatments	0.593	0.593	-	0.589	0.593
CCSS	Breast cancer (289)	PRS	0.365	0.365	-	0.369	0.365
CCSS	Breast cancer (289)	Lifestyle	-	-	-	-	-
CCSS	Breast cancer (289)	Combined	0.743	0.743	-	0.741	0.743
CCSS	Thyroid cancer (163)	Radiation	0.443	0.445	0.439	0.413	0.489
CCSS	Thyroid cancer (163)	Chemo	0.055	0.055	0.056	0.073	0.029
CCSS	Thyroid cancer (163)	All_treatments	0.476	0.477	0.472	0.457	0.504
CCSS	Thyroid cancer (163)	PRS	0.358	0.362	0.352	0.355	0.363
CCSS	Thyroid cancer (163)	Lifestyle	-	-	-	-	-
CCSS	Thyroid cancer (163)	Combined	0.664	0.665	0.661	0.650	0.684
CCSS	Meningioma (255)	Radiation	0.370	0.333	0.422	0.394	0.345
CCSS	Meningioma (255)	Chemo	0.082	0.076	0.092	0.118	0.045
CCSS	Meningioma (255)	All_treatments	0.426	0.391	0.475	0.471	0.378
CCSS	Meningioma (255)	PRS	0.032	0.032	0.033	0.036	0.028
CCSS	Meningioma (255)	Lifestyle	-	-	-	-	-
CCSS	Meningioma (255)	Combined	0.444	0.410	0.494	0.490	0.397
CCSS	Sarcoma (60)	Radiation	-	-	-	-	-
CCSS	Sarcoma (60)	Chemo	0.342	0.327	0.358	0.337	0.348
CCSS	Sarcoma (60)	All_treatments	0.342	0.327	0.358	0.337	0.348
CCSS	Sarcoma (60)	PRS	0.024	0.026	0.022	0.023	0.025
CCSS	Sarcoma (60)	Lifestyle	-	-	-	-	-
CCSS	Sarcoma (60)	Combined	0.358	0.345	0.372	0.353	0.364", header = T, sep = "\t")

data$SN_types_new <- gsub("\\([0-9]+\\)", "", data$SN_types)
data[data == "-"] <- NA


lifestyle <- "without_lifestyle"
all.group <- c("Overall", "Female", "Male", "Age.lt.35", "Age.ge.35")
variables <-  c("Combined", "Radiation", "Chemo", "treatments", "PRS", "Lifestyle")


# Define the text file path
data_name <- paste0("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/ANALYSIS/results/plots/v18b/", lifestyle, "_all_table_18b.txt")

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

############################
## Lifestyle without diet ##
############################
data <- read.table(text = "Cohort	SN_types	Variables	Overall	Female	Male	Age.lt.35	Age.ge.35
SJLIFE	Any SN (303)	Radiation	0.442	0.432	0.455	0.368	0.473
SJLIFE	Any SN (303)	Chemo	0.017	0.017	0.016	0.031	0.011
SJLIFE	Any SN (303)	All_treatments	0.450	0.442	0.461	0.385	0.477
SJLIFE	Any SN (303)	PRS	0.179	0.177	0.183	0.182	0.178
SJLIFE	Any SN (303)	Lifestyle	-0.063	-0.058	-0.069	-0.058	-0.065
SJLIFE	Any SN (303)	Combined	0.520	0.514	0.527	0.462	0.543
SJLIFE	Any SMN (234)	Radiation	0.379	0.378	0.379	0.289	0.413
SJLIFE	Any SMN (234)	Chemo	0.017	0.018	0.017	0.030	0.013
SJLIFE	Any SMN (234)	All_treatments	0.389	0.390	0.388	0.308	0.421
SJLIFE	Any SMN (234)	PRS	0.143	0.140	0.146	0.142	0.143
SJLIFE	Any SMN (234)	Lifestyle	-0.083	-0.073	-0.096	-0.082	-0.084
SJLIFE	Any SMN (234)	Combined	0.431	0.437	0.423	0.354	0.460
SJLIFE	NMSC (118)	Radiation	0.294	0.282	0.310	0.239	0.306
SJLIFE	NMSC (118)	Chemo	-	-	-	-	-
SJLIFE	NMSC (118)	All_treatments	0.294	0.282	0.310	0.239	0.306
SJLIFE	NMSC (118)	PRS	0.432	0.418	0.448	0.431	0.432
SJLIFE	NMSC (118)	Lifestyle	-0.275	-0.251	-0.307	-0.242	-0.283
SJLIFE	NMSC (118)	Combined	0.490	0.482	0.500	0.459	0.497
SJLIFE	Breast cancer (51)	Radiation	0.466	0.466	-	0.391	0.478
SJLIFE	Breast cancer (51)	Chemo	0.106	0.106	-	0.102	0.107
SJLIFE	Breast cancer (51)	All_treatments	0.526	0.526	-	0.452	0.538
SJLIFE	Breast cancer (51)	PRS	0.222	0.222	-	0.250	0.218
SJLIFE	Breast cancer (51)	Lifestyle	-0.872	-0.872	-	-0.581	-0.919
SJLIFE	Breast cancer (51)	Combined	0.335	0.335	-	0.337	0.335
SJLIFE	Thyroid cancer (43)	Radiation	0.610	0.633	0.585	0.555	0.637
SJLIFE	Thyroid cancer (43)	Chemo	0.238	0.237	0.238	0.332	0.192
SJLIFE	Thyroid cancer (43)	All_treatments	0.724	0.740	0.706	0.708	0.732
SJLIFE	Thyroid cancer (43)	PRS	0.589	0.587	0.592	0.601	0.583
SJLIFE	Thyroid cancer (43)	Lifestyle	-0.178	-0.189	-0.166	-0.146	-0.194
SJLIFE	Thyroid cancer (43)	Combined	0.870	0.875	0.863	0.863	0.873
SJLIFE	Meningioma (81)	Radiation	0.184	0.150	0.225	0.245	0.164
SJLIFE	Meningioma (81)	Chemo	0.321	0.325	0.316	0.462	0.274
SJLIFE	Meningioma (81)	All_treatments	0.456	0.436	0.480	0.623	0.400
SJLIFE	Meningioma (81)	PRS	-0.117	-0.113	-0.122	-0.114	-0.118
SJLIFE	Meningioma (81)	Lifestyle	-0.135	-0.126	-0.146	-0.138	-0.134
SJLIFE	Meningioma (81)	Combined	0.308	0.293	0.325	0.509	0.240
SJLIFE	Sarcoma (NA)	Radiation	-	-	-	-	-
SJLIFE	Sarcoma (NA)	Chemo	-	-	-	-	-
SJLIFE	Sarcoma (NA)	All_treatments	-	-	-	-	-
SJLIFE	Sarcoma (NA)	PRS	-	-	-	-	-
SJLIFE	Sarcoma (NA)	Lifestyle	-	-	-	-	-
SJLIFE	Sarcoma (NA)	Combined	-	-	-	-	-
CCSS	Any SN (1286)	Radiation	0.358	0.351	0.369	0.338	0.368
CCSS	Any SN (1286)	Chemo	0.038	0.037	0.041	0.061	0.027
CCSS	Any SN (1286)	All_treatments	0.384	0.376	0.396	0.380	0.387
CCSS	Any SN (1286)	PRS	0.057	0.056	0.058	0.056	0.057
CCSS	Any SN (1286)	Lifestyle	-0.074	-0.071	-0.078	-0.068	-0.077
CCSS	Any SN (1286)	Combined	0.376	0.370	0.386	0.374	0.378
CCSS	Any SMN (619)	Radiation	0.249	0.246	0.255	0.186	0.275
CCSS	Any SMN (619)	Chemo	0.056	0.055	0.058	0.099	0.038
CCSS	Any SMN (619)	All_treatments	0.297	0.293	0.306	0.271	0.308
CCSS	Any SMN (619)	PRS	0.047	0.047	0.048	0.047	0.047
CCSS	Any SMN (619)	Lifestyle	-0.048	-0.046	-0.052	-0.046	-0.049
CCSS	Any SMN (619)	Combined	0.297	0.295	0.303	0.272	0.308
CCSS	NMSC (632)	Radiation	0.390	0.386	0.394	0.400	0.386
CCSS	NMSC (632)	Chemo	-	-	-	-	-
CCSS	NMSC (632)	All_treatments	0.390	0.386	0.394	0.400	0.386
CCSS	NMSC (632)	PRS	0.303	0.305	0.300	0.302	0.303
CCSS	NMSC (632)	Lifestyle	-0.149	-0.148	-0.152	-0.132	-0.156
CCSS	NMSC (632)	Combined	0.514	0.515	0.513	0.530	0.508
CCSS	Breast cancer (259)	Radiation	0.444	0.444	-	0.392	0.450
CCSS	Breast cancer (259)	Chemo	0.194	0.194	-	0.234	0.189
CCSS	Breast cancer (259)	All_treatments	0.572	0.572	-	0.541	0.575
CCSS	Breast cancer (259)	PRS	0.346	0.346	-	0.348	0.345
CCSS	Breast cancer (259)	Lifestyle	-0.088	-0.088	-	-0.078	-0.089
CCSS	Breast cancer (259)	Combined	0.696	0.696	-	0.678	0.698
CCSS	Thyroid cancer (123)	Radiation	0.416	0.416	0.418	0.360	0.486
CCSS	Thyroid cancer (123)	Chemo	0.095	0.093	0.097	0.122	0.061
CCSS	Thyroid cancer (123)	All_treatments	0.479	0.477	0.482	0.443	0.523
CCSS	Thyroid cancer (123)	PRS	0.386	0.391	0.376	0.379	0.394
CCSS	Thyroid cancer (123)	Lifestyle	-0.042	-0.036	-0.054	-0.025	-0.063
CCSS	Thyroid cancer (123)	Combined	0.666	0.667	0.666	0.646	0.692
CCSS	Meningioma (209)	Radiation	0.294	0.260	0.345	0.315	0.279
CCSS	Meningioma (209)	Chemo	0.087	0.080	0.098	0.134	0.055
CCSS	Meningioma (209)	All_treatments	0.368	0.334	0.420	0.421	0.331
CCSS	Meningioma (209)	PRS	0.028	0.027	0.029	0.033	0.024
CCSS	Meningioma (209)	Lifestyle	0.087	0.099	0.068	0.150	0.043
CCSS	Meningioma (209)	Combined	0.433	0.409	0.469	0.524	0.370
CCSS	Sarcoma (41)	Radiation	-	-	-	-	-
CCSS	Sarcoma (41)	Chemo	0.279	0.266	0.300	0.279	0.278
CCSS	Sarcoma (41)	All_treatments	0.279	0.266	0.300	0.279	0.278
CCSS	Sarcoma (41)	PRS	0.035	0.038	0.030	0.033	0.036
CCSS	Sarcoma (41)	Lifestyle	-0.261	-0.256	-0.271	-0.396	-0.180
CCSS	Sarcoma (41)	Combined	0.116	0.108	0.130	0.017	0.176", header = T, sep = "\t")

data$SN_types_new <- gsub("\\([0-9]+\\)", "", data$SN_types)
data[data == "-"] <- NA


lifestyle <- "without_diet"
all.group <- c("Overall", "Female", "Male", "Age.lt.35", "Age.ge.35")
variables <-  c("Combined", "Radiation", "Chemo", "treatments", "PRS", "Lifestyle")


# Define the text file path
data_name <- paste0("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/ANALYSIS/results/plots/v18b/", lifestyle, "_all_table_18b.txt")

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




## Lifestyle with diet (SJLIFE only)
data <- read.table(text="Cohort	SN_types	Variables	Overall	Female	Male	Age.lt.35	Age.ge.35
SJLIFE	Any SN (303)	Radiation	0.442	0.430	0.456	0.368	0.472
SJLIFE	Any SN (303)	Chemo	0.016	0.017	0.015	0.030	0.011
SJLIFE	Any SN (303)	All_treatments	0.449	0.439	0.462	0.384	0.476
SJLIFE	Any SN (303)	PRS	0.179	0.177	0.183	0.182	0.178
SJLIFE	Any SN (303)	Lifestyle	0.199	0.193	0.206	0.203	0.197
SJLIFE	Any SN (303)	Combined	0.637	0.629	0.648	0.594	0.655
SJLIFE	Any SMN (234)	Radiation	0.378	0.375	0.383	0.289	0.412
SJLIFE	Any SMN (234)	Chemo	0.016	0.017	0.015	0.028	0.012
SJLIFE	Any SMN (234)	All_treatments	0.388	0.386	0.391	0.307	0.420
SJLIFE	Any SMN (234)	PRS	0.139	0.137	0.143	0.139	0.140
SJLIFE	Any SMN (234)	Lifestyle	0.266	0.260	0.273	0.269	0.265
SJLIFE	Any SMN (234)	Combined	0.613	0.610	0.615	0.560	0.633
SJLIFE	NMSC (118)	Radiation	0.295	0.282	0.310	0.239	0.307
SJLIFE	NMSC (118)	Chemo	-	-	-	-	-
SJLIFE	NMSC (118)	All_treatments	0.295	0.282	0.310	0.239	0.307
SJLIFE	NMSC (118)	PRS	0.431	0.418	0.447	0.430	0.431
SJLIFE	NMSC (118)	Lifestyle	-0.401	-0.375	-0.435	-0.375	-0.407
SJLIFE	NMSC (118)	Combined	0.439	0.430	0.450	0.400	0.448
SJLIFE	Breast cancer (51)	Radiation	0.473	0.473	-	0.408	0.484
SJLIFE	Breast cancer (51)	Chemo	0.118	0.118	-	0.112	0.119
SJLIFE	Breast cancer (51)	All_treatments	0.540	0.540	-	0.474	0.551
SJLIFE	Breast cancer (51)	PRS	0.216	0.216	-	0.239	0.213
SJLIFE	Breast cancer (51)	Lifestyle	0.199	0.199	-	0.349	0.174
SJLIFE	Breast cancer (51)	Combined	0.725	0.725	-	0.734	0.723
SJLIFE	Thyroid cancer (43)	Radiation	0.609	0.632	0.583	0.551	0.638
SJLIFE	Thyroid cancer (43)	Chemo	0.237	0.238	0.236	0.330	0.191
SJLIFE	Thyroid cancer (43)	All_treatments	0.723	0.739	0.705	0.705	0.732
SJLIFE	Thyroid cancer (43)	PRS	0.587	0.585	0.589	0.598	0.581
SJLIFE	Thyroid cancer (43)	Lifestyle	0.252	0.239	0.268	0.267	0.245
SJLIFE	Thyroid cancer (43)	Combined	0.916	0.919	0.913	0.911	0.919
SJLIFE	Meningioma (81)	Radiation	0.184	0.149	0.224	0.244	0.163
SJLIFE	Meningioma (81)	Chemo	0.319	0.323	0.315	0.461	0.271
SJLIFE	Meningioma (81)	All_treatments	0.453	0.432	0.478	0.620	0.397
SJLIFE	Meningioma (81)	PRS	-0.117	-0.113	-0.123	-0.116	-0.118
SJLIFE	Meningioma (81)	Lifestyle	0.269	0.267	0.272	0.270	0.269
SJLIFE	Meningioma (81)	Combined	0.550	0.536	0.567	0.681	0.506", header = T, sep = "\t")


data$SN_types_new <- gsub("\\([0-9]+\\)", "", data$SN_types)
data[data == "-"] <- NA


lifestyle <- "lifestyle_with_diet"
all.group <- c("Overall", "Female", "Male", "Age.lt.35", "Age.ge.35")
variables <-  c("Combined", "Radiation", "Chemo", "treatments", "PRS", "Lifestyle")


# Define the text file path
data_name <- paste0("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/ANALYSIS/results/plots/v18b/", lifestyle, "_all_table_18b.txt")

file_conn <- file(data_name, "a")

for (j in 1:length(all.group)){
  group <-  all.group [j]
  for(i in 1:length(variables)){
    if(lifestyle== "lifestyle_with_diet" && variables[i] == "Lifestyle"){
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



## create individual tables
## Lifestyle without diet
data <- read.table(text="Cohort	SNtypes	Variables	Overall_all4	Female_all4	Male_all4	Age_lt_35_all4	Age_ge_35_all4	Overall_smk	Female_smk	Male_smk	Age_lt_35_smk	Age_ge_35_smk	Overall_drk	Female_drk	Male_drk	Age_lt_35_drk	Age_ge_35_drk	Overall_PA	Female_PA	Male_PA	Age_lt_35_PA	Age_ge_35_PA	Overall_obese	Female_obese	Male_obese	Age_lt_35_obese	Age_ge_35_obese
SJLIFE	Any_SN	Radiation	0.442	0.432	0.455	0.368	0.473	0.442	0.432	0.455	0.368	0.473	0.442	0.432	0.455	0.368	0.473	0.442	0.432	0.455	0.368	0.473	0.442	0.432	0.455	0.368	0.473
SJLIFE	Any_SN	Chemo	0.017	0.017	0.016	0.031	0.011	0.017	0.017	0.016	0.031	0.011	0.017	0.017	0.016	0.031	0.011	0.017	0.017	0.016	0.031	0.011	0.017	0.017	0.016	0.031	0.011
SJLIFE	Any_SN	All_treatments	0.450	0.442	0.461	0.385	0.477	0.450	0.442	0.461	0.385	0.477	0.450	0.442	0.461	0.385	0.477	0.450	0.442	0.461	0.385	0.477	0.450	0.442	0.461	0.385	0.477
SJLIFE	Any_SN	PRS	0.179	0.177	0.183	0.182	0.178	0.179	0.177	0.183	0.182	0.178	0.179	0.177	0.183	0.182	0.178	0.179	0.177	0.183	0.182	0.178	0.179	0.177	0.183	0.182	0.178
SJLIFE	Any_SN	Lifestyle	-0.063	-0.058	-0.069	-0.058	-0.065	-0.057	-0.059	-0.054	-0.052	-0.059	-0.036	-0.033	-0.039	-0.039	-0.034	0.010	0.014	0.006	0.005	0.012	-0.027	-0.028	-0.026	-0.028	-0.027
SJLIFE	Any_SN	Combined	0.520	0.514	0.527	0.462	0.543	0.523	0.515	0.534	0.466	0.546	0.533	0.526	0.542	0.475	0.556	0.554	0.548	0.562	0.498	0.577	0.538	0.530	0.548	0.481	0.560
SJLIFE	Any_SMN	Radiation	0.379	0.378	0.379	0.289	0.413	0.379	0.378	0.379	0.289	0.413	0.379	0.378	0.379	0.289	0.413	0.379	0.378	0.379	0.289	0.413	0.379	0.378	0.379	0.289	0.413
SJLIFE	Any_SMN	Chemo	0.017	0.018	0.017	0.030	0.013	0.017	0.018	0.017	0.030	0.013	0.017	0.018	0.017	0.030	0.013	0.017	0.018	0.017	0.030	0.013	0.017	0.018	0.017	0.030	0.013
SJLIFE	Any_SMN	All_treatments	0.389	0.390	0.388	0.308	0.421	0.389	0.390	0.388	0.308	0.421	0.389	0.390	0.388	0.308	0.421	0.389	0.390	0.388	0.308	0.421	0.389	0.390	0.388	0.308	0.421
SJLIFE	Any_SMN	PRS	0.143	0.140	0.146	0.142	0.143	0.143	0.140	0.146	0.142	0.143	0.143	0.140	0.146	0.142	0.143	0.143	0.140	0.146	0.142	0.143	0.143	0.140	0.146	0.142	0.143
SJLIFE	Any_SMN	Lifestyle	-0.083	-0.073	-0.096	-0.082	-0.084	-0.079	-0.081	-0.077	-0.082	-0.078	-0.051	-0.051	-0.052	-0.059	-0.048	0.033	0.043	0.021	0.017	0.040	-0.111	-0.114	-0.109	-0.105	-0.114
SJLIFE	Any_SMN	Combined	0.431	0.437	0.423	0.354	0.460	0.434	0.433	0.435	0.355	0.465	0.449	0.449	0.450	0.370	0.480	0.494	0.498	0.488	0.415	0.524	0.417	0.417	0.418	0.343	0.446
SJLIFE	NMSC	Radiation	0.294	0.282	0.310	0.239	0.306	0.294	0.282	0.310	0.239	0.306	0.294	0.282	0.310	0.239	0.306	0.294	0.282	0.310	0.239	0.306	0.294	0.282	0.310	0.239	0.306
SJLIFE	NMSC	Chemo	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-
SJLIFE	NMSC	All_treatments	0.294	0.282	0.310	0.239	0.306	0.294	0.282	0.310	0.239	0.306	0.294	0.282	0.310	0.239	0.306	0.294	0.282	0.310	0.239	0.306	0.294	0.282	0.310	0.239	0.306
SJLIFE	NMSC	PRS	0.432	0.418	0.448	0.431	0.432	0.432	0.418	0.448	0.431	0.432	0.432	0.418	0.448	0.431	0.432	0.432	0.418	0.448	0.431	0.432	0.432	0.418	0.448	0.431	0.432
SJLIFE	NMSC	Lifestyle	-0.275	-0.251	-0.307	-0.242	-0.283	-0.024	-0.026	-0.022	-0.039	-0.021	-0.155	-0.134	-0.182	-0.149	-0.157	-0.042	-0.044	-0.039	-0.054	-0.039	-0.149	-0.149	-0.149	-0.136	-0.152
SJLIFE	NMSC	Combined	0.490	0.482	0.500	0.459	0.497	0.591	0.576	0.609	0.551	0.599	0.490	0.482	0.500	0.459	0.497	0.584	0.568	0.603	0.544	0.592	0.540	0.525	0.560	0.508	0.547
SJLIFE	Breast_cancer	Radiation	0.466	0.466	-	0.391	0.478	0.466	0.466	-	0.391	0.478	0.466	0.466	-	0.391	0.478	0.466	0.466	-	0.391	0.478	0.466	0.466	-	0.391	0.478
SJLIFE	Breast_cancer	Chemo	0.106	0.106	-	0.102	0.107	0.106	0.106	-	0.102	0.107	0.106	0.106	-	0.102	0.107	0.106	0.106	-	0.102	0.107	0.106	0.106	-	0.102	0.107
SJLIFE	Breast_cancer	All_treatments	0.526	0.526	-	0.452	0.538	0.526	0.526	-	0.452	0.538	0.526	0.526	-	0.452	0.538	0.526	0.526	-	0.452	0.538	0.526	0.526	-	0.452	0.538
SJLIFE	Breast_cancer	PRS	0.222	0.222	-	0.250	0.218	0.222	0.222	-	0.250	0.218	0.222	0.222	-	0.250	0.218	0.222	0.222	-	0.250	0.218	0.222	0.222	-	0.250	0.218
SJLIFE	Breast_cancer	Lifestyle	-0.872	-0.872	-	-0.581	-0.919	-0.256	-0.256	-	-0.159	-0.272	0.036	0.036	-	0.034	0.037	-0.375	-0.375	-	-0.224	-0.399	0.004	0.004	-	0.005	0.004
SJLIFE	Breast_cancer	Combined	0.335	0.335	-	0.337	0.335	0.553	0.553	-	0.525	0.557	0.649	0.649	-	0.608	0.656	0.512	0.512	-	0.494	0.515	0.635	0.635	-	0.597	0.641
SJLIFE	Thyroid_cancer	Radiation	0.610	0.633	0.585	0.555	0.637	0.610	0.633	0.585	0.555	0.637	0.610	0.633	0.585	0.555	0.637	0.610	0.633	0.585	0.555	0.637	0.610	0.633	0.585	0.555	0.637
SJLIFE	Thyroid_cancer	Chemo	0.238	0.237	0.238	0.332	0.192	0.238	0.237	0.238	0.332	0.192	0.238	0.237	0.238	0.332	0.192	0.238	0.237	0.238	0.332	0.192	0.238	0.237	0.238	0.332	0.192
SJLIFE	Thyroid_cancer	All_treatments	0.724	0.740	0.706	0.708	0.732	0.724	0.740	0.706	0.708	0.732	0.724	0.740	0.706	0.708	0.732	0.724	0.740	0.706	0.708	0.732	0.724	0.740	0.706	0.708	0.732
SJLIFE	Thyroid_cancer	PRS	0.589	0.587	0.592	0.601	0.583	0.589	0.587	0.592	0.601	0.583	0.589	0.587	0.592	0.601	0.583	0.589	0.587	0.592	0.601	0.583	0.589	0.587	0.592	0.601	0.583
SJLIFE	Thyroid_cancer	Lifestyle	-0.178	-0.189	-0.166	-0.146	-0.194	0.093	0.099	0.085	0.086	0.096	-0.149	-0.134	-0.166	-0.151	-0.148	-0.189	-0.216	-0.157	-0.163	-0.201	-0.101	-0.102	-0.100	-0.092	-0.105
SJLIFE	Thyroid_cancer	Combined	0.870	0.875	0.863	0.863	0.873	0.899	0.905	0.894	0.891	0.904	0.870	0.878	0.861	0.860	0.875	0.867	0.870	0.863	0.859	0.871	0.877	0.883	0.870	0.869	0.881
SJLIFE	Meningioma	Radiation	0.184	0.150	0.225	0.245	0.164	0.184	0.150	0.225	0.245	0.164	0.184	0.150	0.225	0.245	0.164	0.184	0.150	0.225	0.245	0.164	0.184	0.150	0.225	0.245	0.164
SJLIFE	Meningioma	Chemo	0.321	0.325	0.316	0.462	0.274	0.321	0.325	0.316	0.462	0.274	0.321	0.325	0.316	0.462	0.274	0.321	0.325	0.316	0.462	0.274	0.321	0.325	0.316	0.462	0.274
SJLIFE	Meningioma	All_treatments	0.456	0.436	0.480	0.623	0.400	0.456	0.436	0.480	0.623	0.400	0.456	0.436	0.480	0.623	0.400	0.456	0.436	0.480	0.623	0.400	0.456	0.436	0.480	0.623	0.400
SJLIFE	Meningioma	PRS	-0.117	-0.113	-0.122	-0.114	-0.118	-0.117	-0.113	-0.122	-0.114	-0.118	-0.117	-0.113	-0.122	-0.114	-0.118	-0.117	-0.113	-0.122	-0.114	-0.118	-0.117	-0.113	-0.122	-0.114	-0.118
SJLIFE	Meningioma	Lifestyle	-0.135	-0.126	-0.146	-0.138	-0.134	-0.053	-0.051	-0.055	-0.046	-0.055	-0.157	-0.129	-0.191	-0.147	-0.161	-0.086	-0.100	-0.070	-0.076	-0.089	0.149	0.143	0.157	0.124	0.158
SJLIFE	Meningioma	Combined	0.308	0.293	0.325	0.509	0.240	0.361	0.340	0.385	0.557	0.295	0.296	0.291	0.302	0.510	0.224	0.343	0.314	0.377	0.549	0.274	0.484	0.464	0.509	0.629	0.436
CCSS	Any_SN	Radiation	0.358	0.351	0.369	0.338	0.368	0.358	0.351	0.369	0.338	0.368	0.358	0.351	0.369	0.338	0.368	0.358	0.351	0.369	0.338	0.368	0.358	0.351	0.369	0.338	0.368
CCSS	Any_SN	Chemo	0.038	0.037	0.041	0.061	0.027	0.038	0.037	0.041	0.061	0.027	0.038	0.037	0.041	0.061	0.027	0.038	0.037	0.041	0.061	0.027	0.038	0.037	0.041	0.061	0.027
CCSS	Any_SN	All_treatments	0.384	0.376	0.396	0.380	0.387	0.384	0.376	0.396	0.380	0.387	0.384	0.376	0.396	0.380	0.387	0.384	0.376	0.396	0.380	0.387	0.384	0.376	0.396	0.380	0.387
CCSS	Any_SN	PRS	0.057	0.056	0.058	0.056	0.057	0.057	0.056	0.058	0.056	0.057	0.057	0.056	0.058	0.056	0.057	0.057	0.056	0.058	0.056	0.057	0.057	0.056	0.058	0.056	0.057
CCSS	Any_SN	Lifestyle	-0.074	-0.071	-0.078	-0.068	-0.077	-0.052	-0.048	-0.057	-0.045	-0.055	-0.003	-0.002	-0.003	-0.002	-0.003	0.002	0.002	0.002	0.003	0.002	-0.017	-0.018	-0.015	-0.016	-0.017
CCSS	Any_SN	Combined	0.376	0.370	0.386	0.374	0.378	0.389	0.383	0.397	0.386	0.390	0.418	0.410	0.429	0.413	0.420	0.420	0.413	0.432	0.416	0.422	0.410	0.402	0.422	0.406	0.412
CCSS	Any_SMN	Radiation	0.249	0.246	0.255	0.186	0.275	0.249	0.246	0.255	0.186	0.275	0.249	0.246	0.255	0.186	0.275	0.249	0.246	0.255	0.186	0.275	0.249	0.246	0.255	0.186	0.275
CCSS	Any_SMN	Chemo	0.056	0.055	0.058	0.099	0.038	0.056	0.055	0.058	0.099	0.038	0.056	0.055	0.058	0.099	0.038	0.056	0.055	0.058	0.099	0.038	0.056	0.055	0.058	0.099	0.038
CCSS	Any_SMN	All_treatments	0.297	0.293	0.306	0.271	0.308	0.297	0.293	0.306	0.271	0.308	0.297	0.293	0.306	0.271	0.308	0.297	0.293	0.306	0.271	0.308	0.297	0.293	0.306	0.271	0.308
CCSS	Any_SMN	PRS	0.047	0.047	0.048	0.047	0.047	0.047	0.047	0.048	0.047	0.047	0.047	0.047	0.048	0.047	0.047	0.047	0.047	0.048	0.047	0.047	0.047	0.047	0.048	0.047	0.047
CCSS	Any_SMN	Lifestyle	-0.048	-0.046	-0.052	-0.046	-0.049	-0.037	-0.035	-0.040	-0.034	-0.038	-0.004	-0.004	-0.005	-0.004	-0.004	-0.002	-0.002	-0.002	-0.001	-0.002	-0.002	-0.002	-0.002	-0.001	-0.002
CCSS	Any_SMN	Combined	0.297	0.295	0.303	0.272	0.308	0.305	0.302	0.311	0.281	0.315	0.327	0.323	0.335	0.302	0.338	0.329	0.324	0.337	0.304	0.339	0.329	0.324	0.337	0.304	0.339
CCSS	NMSC	Radiation	0.390	0.386	0.394	0.400	0.386	0.390	0.386	0.394	0.400	0.386	0.390	0.386	0.394	0.400	0.386	0.390	0.386	0.394	0.400	0.386	0.390	0.386	0.394	0.400	0.386
CCSS	NMSC	Chemo	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-
CCSS	NMSC	All_treatments	0.390	0.386	0.394	0.400	0.386	0.390	0.386	0.394	0.400	0.386	0.390	0.386	0.394	0.400	0.386	0.390	0.386	0.394	0.400	0.386	0.390	0.386	0.394	0.400	0.386
CCSS	NMSC	PRS	0.303	0.305	0.300	0.302	0.303	0.303	0.305	0.300	0.302	0.303	0.303	0.305	0.300	0.302	0.303	0.303	0.305	0.300	0.302	0.303	0.303	0.305	0.300	0.302	0.303
CCSS	NMSC	Lifestyle	-0.149	-0.148	-0.152	-0.132	-0.156	-0.067	-0.057	-0.079	-0.055	-0.072	0.021	0.020	0.022	0.028	0.018	-0.033	-0.034	-0.031	-0.025	-0.036	-0.043	-0.047	-0.039	-0.037	-0.046
CCSS	NMSC	Combined	0.514	0.515	0.513	0.530	0.508	0.546	0.550	0.542	0.559	0.541	0.584	0.583	0.585	0.595	0.580	0.561	0.560	0.562	0.572	0.556	0.558	0.557	0.560	0.569	0.554
CCSS	Breast_cancer	Radiation	0.444	0.444	-	0.392	0.450	0.444	0.444	-	0.392	0.450	0.444	0.444	-	0.392	0.450	0.444	0.444	-	0.392	0.450	0.444	0.444	-	0.392	0.450
CCSS	Breast_cancer	Chemo	0.194	0.194	-	0.234	0.189	0.194	0.194	-	0.234	0.189	0.194	0.194	-	0.234	0.189	0.194	0.194	-	0.234	0.189	0.194	0.194	-	0.234	0.189
CCSS	Breast_cancer	All_treatments	0.572	0.572	-	0.541	0.575	0.572	0.572	-	0.541	0.575	0.572	0.572	-	0.541	0.575	0.572	0.572	-	0.541	0.575	0.572	0.572	-	0.541	0.575
CCSS	Breast_cancer	PRS	0.346	0.346	-	0.348	0.345	0.346	0.346	-	0.348	0.345	0.346	0.346	-	0.348	0.345	0.346	0.346	-	0.348	0.345	0.346	0.346	-	0.348	0.345
CCSS	Breast_cancer	Lifestyle	-0.088	-0.088	-	-0.078	-0.089	-0.039	-0.039	-	-0.035	-0.039	0.021	0.021	-	0.029	0.020	-0.029	-0.029	-	-0.024	-0.029	-0.026	-0.026	-	-0.024	-0.026
CCSS	Breast_cancer	Combined	0.696	0.696	-	0.678	0.698	0.710	0.710	-	0.692	0.712	0.727	0.727	-	0.711	0.729	0.714	0.714	-	0.695	0.716	0.714	0.714	-	0.695	0.716
CCSS	Thyroid_cancer	Radiation	0.416	0.416	0.418	0.360	0.486	0.416	0.416	0.418	0.360	0.486	0.416	0.416	0.418	0.360	0.486	0.416	0.416	0.418	0.360	0.486	0.416	0.416	0.418	0.360	0.486
CCSS	Thyroid_cancer	Chemo	0.095	0.093	0.097	0.122	0.061	0.095	0.093	0.097	0.122	0.061	0.095	0.093	0.097	0.122	0.061	0.095	0.093	0.097	0.122	0.061	0.095	0.093	0.097	0.122	0.061
CCSS	Thyroid_cancer	All_treatments	0.479	0.477	0.482	0.443	0.523	0.479	0.477	0.482	0.443	0.523	0.479	0.477	0.482	0.443	0.523	0.479	0.477	0.482	0.443	0.523	0.479	0.477	0.482	0.443	0.523
CCSS	Thyroid_cancer	PRS	0.386	0.391	0.376	0.379	0.394	0.386	0.391	0.376	0.379	0.394	0.386	0.391	0.376	0.379	0.394	0.386	0.391	0.376	0.379	0.394	0.386	0.391	0.376	0.379	0.394
CCSS	Thyroid_cancer	Lifestyle	-0.042	-0.036	-0.054	-0.025	-0.063	-0.037	-0.028	-0.053	-0.023	-0.053	0.062	0.061	0.064	0.075	0.046	0.019	0.020	0.017	0.029	0.007	0.005	0.005	0.004	0.014	-0.006
CCSS	Thyroid_cancer	Combined	0.666	0.667	0.666	0.646	0.692	0.668	0.670	0.664	0.646	0.694	0.701	0.700	0.702	0.682	0.724	0.687	0.686	0.687	0.666	0.712	0.682	0.682	0.683	0.661	0.708
CCSS	Meningioma	Radiation	0.294	0.260	0.345	0.315	0.279	0.294	0.260	0.345	0.315	0.279	0.294	0.260	0.345	0.315	0.279	0.294	0.260	0.345	0.315	0.279	0.294	0.260	0.345	0.315	0.279
CCSS	Meningioma	Chemo	0.087	0.080	0.098	0.134	0.055	0.087	0.080	0.098	0.134	0.055	0.087	0.080	0.098	0.134	0.055	0.087	0.080	0.098	0.134	0.055	0.087	0.080	0.098	0.134	0.055
CCSS	Meningioma	All_treatments	0.368	0.334	0.420	0.421	0.331	0.368	0.334	0.420	0.421	0.331	0.368	0.334	0.420	0.421	0.331	0.368	0.334	0.420	0.421	0.331	0.368	0.334	0.420	0.421	0.331
CCSS	Meningioma	PRS	0.028	0.027	0.029	0.033	0.024	0.028	0.027	0.029	0.033	0.024	0.028	0.027	0.029	0.033	0.024	0.028	0.027	0.029	0.033	0.024	0.028	0.027	0.029	0.033	0.024
CCSS	Meningioma	Lifestyle	0.087	0.099	0.068	0.150	0.043	0.093	0.099	0.083	0.156	0.049	0.056	0.065	0.042	0.119	0.012	0.121	0.128	0.112	0.181	0.080	0.127	0.135	0.116	0.189	0.084
CCSS	Meningioma	Combined	0.433	0.409	0.469	0.524	0.370	0.442	0.415	0.484	0.533	0.380	0.417	0.390	0.457	0.510	0.352	0.461	0.434	0.501	0.547	0.400	0.463	0.437	0.503	0.550	0.402
CCSS	Sarcoma	Radiation	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-
CCSS	Sarcoma	Chemo	0.279	0.266	0.300	0.279	0.278	0.279	0.266	0.300	0.279	0.278	0.279	0.266	0.300	0.279	0.278	0.279	0.266	0.300	0.279	0.278	0.279	0.266	0.300	0.279	0.278
CCSS	Sarcoma	All_treatments	0.279	0.266	0.300	0.279	0.278	0.279	0.266	0.300	0.279	0.278	0.279	0.266	0.300	0.279	0.278	0.279	0.266	0.300	0.279	0.278	0.279	0.266	0.300	0.279	0.278
CCSS	Sarcoma	PRS	0.035	0.038	0.030	0.033	0.036	0.035	0.038	0.030	0.033	0.036	0.035	0.038	0.030	0.033	0.036	0.035	0.038	0.030	0.033	0.036	0.035	0.038	0.030	0.033	0.036
CCSS	Sarcoma	Lifestyle	-0.261	-0.256	-0.271	-0.396	-0.180	-0.121	-0.126	-0.112	-0.215	-0.064	-0.229	-0.225	-0.235	-0.360	-0.149	-0.090	-0.096	-0.080	-0.183	-0.033	-0.129	-0.136	-0.116	-0.226	-0.071
CCSS	Sarcoma	Combined	0.116	0.108	0.130	0.017	0.176	0.213	0.199	0.238	0.144	0.255	0.140	0.131	0.155	0.043	0.198	0.116	0.108	0.130	0.017	0.176	0.208	0.192	0.235	0.137	0.251", header = T, sep = "\t")

## combined
data.sjlife <- data[data$Cohort=="SJLIFE",]
data.ccss <- data[data$Cohort=="CCSS",]

data.sjlife.wanted <- data.sjlife[grepl("Lifestyle", data.sjlife$Variables),]
data.ccss.wanted <- data.ccss[grepl("Lifestyle", data.ccss$Variables),]

variables <- c("drk", "smk", "PA", "obese")


i=4
sjlife <- data.sjlife.wanted[grepl(paste0("SNtypes|",variables[i]), colnames(data.sjlife.wanted))]
sjlife <- rbind.data.frame(sjlife, c("Sarcoma",NA,NA,NA,NA,NA))
colnames(sjlife)<- paste0(colnames(sjlife), "_sjlife")
ccss <- data.ccss.wanted[grepl(paste0("SNtypes|",variables[i]), colnames(data.ccss.wanted))]
colnames(ccss)<- paste0(colnames(ccss), "_ccss")

wanted <- cbind(sjlife,ccss)[order(c(seq_along(sjlife), seq_along(ccss)))]
View(wanted)

###########################################################
## Lifestyle with diet (SJLIFE only)
data <- read.table(text="SNtypes	Variables	Overall_all4	Female_all4	Male_all4	Age_lt_35_all4	Age_ge_35_all4	Overall_smk	Female_smk	Male_smk	Age_lt_35_smk	Age_ge_35_smk	Overall_drk	Female_drk	Male_drk	Age_lt_35_drk	Age_ge_35_drk	Overall_PA	Female_PA	Male_PA	Age_lt_35_PA	Age_ge_35_PA	Overall_obese	Female_obese	Male_obese	Age_lt_35_obese	Age_ge_35_obese	Overall_diet	Female_diet	Male_diet	Age_lt_35_diet	Age_ge_35_diet
Any_SN	Radiation	0.442	0.430	0.456	0.368	0.472	0.442	0.430	0.456	0.368	0.472	0.442	0.430	0.456	0.368	0.472	0.442	0.430	0.456	0.368	0.472	0.442	0.430	0.456	0.368	0.472	0.442	0.430	0.456	0.368	0.472
Any_SN	Chemo	0.016	0.017	0.015	0.030	0.011	0.016	0.017	0.015	0.030	0.011	0.016	0.017	0.015	0.030	0.011	0.016	0.017	0.015	0.030	0.011	0.016	0.017	0.015	0.030	0.011	0.016	0.017	0.015	0.030	0.011
Any_SN	All_treatments	0.449	0.439	0.462	0.384	0.476	0.449	0.439	0.462	0.384	0.476	0.449	0.439	0.462	0.384	0.476	0.449	0.439	0.462	0.384	0.476	0.449	0.439	0.462	0.384	0.476	0.449	0.439	0.462	0.384	0.476
Any_SN	PRS	0.179	0.177	0.183	0.182	0.178	0.179	0.177	0.183	0.182	0.178	0.179	0.177	0.183	0.182	0.178	0.179	0.177	0.183	0.182	0.178	0.179	0.177	0.183	0.182	0.178	0.179	0.177	0.183	0.182	0.178
Any_SN	Lifestyle	0.199	0.193	0.206	0.203	0.197	-0.087	-0.092	-0.082	-0.087	-0.087	-0.064	-0.062	-0.066	-0.072	-0.060	-0.064	-0.062	-0.066	-0.072	-0.060	-0.058	-0.060	-0.055	-0.064	-0.055	0.250	0.243	0.260	0.247	0.251
Any_SN	Combined	0.637	0.629	0.648	0.594	0.655	0.509	0.497	0.523	0.448	0.533	0.519	0.510	0.531	0.458	0.544	0.519	0.510	0.531	0.458	0.544	0.523	0.513	0.536	0.463	0.547	0.662	0.653	0.674	0.619	0.679
Any_SMN	Radiation	0.378	0.375	0.383	0.289	0.412	0.378	0.375	0.383	0.289	0.412	0.378	0.375	0.383	0.289	0.412	0.378	0.375	0.383	0.289	0.412	0.378	0.375	0.383	0.289	0.412	0.378	0.375	0.383	0.289	0.412
Any_SMN	Chemo	0.016	0.017	0.015	0.028	0.012	0.016	0.017	0.015	0.028	0.012	0.016	0.017	0.015	0.028	0.012	0.016	0.017	0.015	0.028	0.012	0.016	0.017	0.015	0.028	0.012	0.016	0.017	0.015	0.028	0.012
Any_SMN	All_treatments	0.388	0.386	0.391	0.307	0.420	0.388	0.386	0.391	0.307	0.420	0.388	0.386	0.391	0.307	0.420	0.388	0.386	0.391	0.307	0.420	0.388	0.386	0.391	0.307	0.420	0.388	0.386	0.391	0.307	0.420
Any_SMN	PRS	0.139	0.137	0.143	0.139	0.140	0.139	0.137	0.143	0.139	0.140	0.139	0.137	0.143	0.139	0.140	0.139	0.137	0.143	0.139	0.140	0.139	0.137	0.143	0.139	0.140	0.139	0.137	0.143	0.139	0.140
Any_SMN	Lifestyle	0.266	0.260	0.273	0.269	0.265	-0.128	-0.135	-0.119	-0.139	-0.124	-0.096	-0.099	-0.093	-0.114	-0.090	-0.023	-0.018	-0.028	-0.046	-0.013	-0.163	-0.169	-0.155	-0.164	-0.162	0.315	0.304	0.329	0.309	0.318
Any_SMN	Combined	0.613	0.610	0.615	0.560	0.633	0.405	0.399	0.412	0.317	0.439	0.422	0.417	0.428	0.334	0.456	0.461	0.460	0.462	0.374	0.494	0.387	0.382	0.393	0.304	0.419	0.640	0.634	0.647	0.586	0.660
NMSC	Radiation	0.295	0.282	0.310	0.239	0.307	0.295	0.282	0.310	0.239	0.307	0.295	0.282	0.310	0.239	0.307	0.295	0.282	0.310	0.239	0.307	0.295	0.282	0.310	0.239	0.307	0.295	0.282	0.310	0.239	0.307
NMSC	Chemo	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-
NMSC	All_treatments	0.295	0.282	0.310	0.239	0.307	0.295	0.282	0.310	0.239	0.307	0.295	0.282	0.310	0.239	0.307	0.295	0.282	0.310	0.239	0.307	0.295	0.282	0.310	0.239	0.307	0.295	0.282	0.310	0.239	0.307
NMSC	PRS	0.431	0.418	0.447	0.430	0.431	0.431	0.418	0.447	0.430	0.431	0.431	0.418	0.447	0.430	0.431	0.431	0.418	0.447	0.430	0.431	0.431	0.418	0.447	0.430	0.431	0.431	0.418	0.447	0.430	0.431
NMSC	Lifestyle	-0.401	-0.375	-0.435	-0.375	-0.407	-0.052	-0.056	-0.046	-0.077	-0.046	-0.191	-0.172	-0.216	-0.195	-0.190	-0.068	-0.073	-0.061	-0.092	-0.063	-0.184	-0.187	-0.180	-0.181	-0.185	-0.128	-0.132	-0.124	-0.154	-0.123
NMSC	Combined	0.439	0.430	0.450	0.400	0.448	0.579	0.563	0.600	0.534	0.589	0.523	0.514	0.534	0.478	0.533	0.573	0.556	0.594	0.527	0.583	0.526	0.509	0.547	0.489	0.534	0.549	0.532	0.570	0.500	0.559
Breast_cancer	Radiation	0.473	0.473	-	0.408	0.484	0.473	0.473	-	0.408	0.484	0.473	0.473	-	0.408	0.484	0.473	0.473	-	0.408	0.484	0.473	0.473	-	0.408	0.484	0.473	0.473	-	0.408	0.484
Breast_cancer	Chemo	0.118	0.118	-	0.112	0.119	0.118	0.118	-	0.112	0.119	0.118	0.118	-	0.112	0.119	0.118	0.118	-	0.112	0.119	0.118	0.118	-	0.112	0.119	0.118	0.118	-	0.112	0.119
Breast_cancer	All_treatments	0.540	0.540	-	0.474	0.551	0.540	0.540	-	0.474	0.551	0.540	0.540	-	0.474	0.551	0.540	0.540	-	0.474	0.551	0.540	0.540	-	0.474	0.551	0.540	0.540	-	0.474	0.551
Breast_cancer	PRS	0.216	0.216	-	0.239	0.213	0.216	0.216	-	0.239	0.213	0.216	0.216	-	0.239	0.213	0.216	0.216	-	0.239	0.213	0.216	0.216	-	0.239	0.213	0.216	0.216	-	0.239	0.213
Breast_cancer	Lifestyle	0.199	0.199	-	0.349	0.174	-0.338	-0.338	-	-0.216	-0.357	-0.007	-0.007	-	-0.007	-0.007	-0.522	-0.522	-	-0.327	-0.553	-0.028	-0.028	-	-0.027	-0.028	0.602	0.602	-	0.610	0.600
Breast_cancer	Combined	0.725	0.725	-	0.734	0.723	0.535	0.535	-	0.512	0.539	0.642	0.642	-	0.602	0.648	0.476	0.476	-	0.462	0.478	0.632	0.632	-	0.594	0.639	0.859	0.859	-	0.846	0.861
Thyroid_cancer	Radiation	0.609	0.632	0.583	0.551	0.638	0.609	0.632	0.583	0.551	0.638	0.609	0.632	0.583	0.551	0.638	0.609	0.632	0.583	0.551	0.638	0.609	0.632	0.583	0.551	0.638	0.609	0.632	0.583	0.551	0.638
Thyroid_cancer	Chemo	0.237	0.238	0.236	0.330	0.191	0.237	0.238	0.236	0.330	0.191	0.237	0.238	0.236	0.330	0.191	0.237	0.238	0.236	0.330	0.191	0.237	0.238	0.236	0.330	0.191	0.237	0.238	0.236	0.330	0.191
Thyroid_cancer	All_treatments	0.723	0.739	0.705	0.705	0.732	0.723	0.739	0.705	0.705	0.732	0.723	0.739	0.705	0.705	0.732	0.723	0.739	0.705	0.705	0.732	0.723	0.739	0.705	0.705	0.732	0.723	0.739	0.705	0.705	0.732
Thyroid_cancer	PRS	0.587	0.585	0.589	0.598	0.581	0.587	0.585	0.589	0.598	0.581	0.587	0.585	0.589	0.598	0.581	0.587	0.585	0.589	0.598	0.581	0.587	0.585	0.589	0.598	0.581	0.587	0.585	0.589	0.598	0.581
Thyroid_cancer	Lifestyle	0.252	0.239	0.268	0.267	0.245	0.043	0.049	0.035	0.024	0.052	-0.215	-0.198	-0.235	-0.233	-0.207	-0.276	-0.311	-0.236	-0.260	-0.283	-0.162	-0.163	-0.160	-0.166	-0.159	0.351	0.347	0.356	0.343	0.356
Thyroid_cancer	Combined	0.916	0.919	0.913	0.911	0.919	0.892	0.898	0.886	0.881	0.898	0.860	0.869	0.851	0.847	0.867	0.855	0.858	0.852	0.845	0.860	0.868	0.875	0.861	0.857	0.874	0.926	0.929	0.923	0.919	0.930
Meningioma	Radiation	0.184	0.149	0.224	0.244	0.163	0.184	0.149	0.224	0.244	0.163	0.184	0.149	0.224	0.244	0.163	0.184	0.149	0.224	0.244	0.163	0.184	0.149	0.224	0.244	0.163	0.184	0.149	0.224	0.244	0.163
Meningioma	Chemo	0.319	0.323	0.315	0.461	0.271	0.319	0.323	0.315	0.461	0.271	0.319	0.323	0.315	0.461	0.271	0.319	0.323	0.315	0.461	0.271	0.319	0.323	0.315	0.461	0.271	0.319	0.323	0.315	0.461	0.271
Meningioma	All_treatments	0.453	0.432	0.478	0.620	0.397	0.453	0.432	0.478	0.620	0.397	0.453	0.432	0.478	0.620	0.397	0.453	0.432	0.478	0.620	0.397	0.453	0.432	0.478	0.620	0.397	0.453	0.432	0.478	0.620	0.397
Meningioma	PRS	-0.117	-0.113	-0.123	-0.116	-0.118	-0.117	-0.113	-0.123	-0.116	-0.118	-0.117	-0.113	-0.123	-0.116	-0.118	-0.117	-0.113	-0.123	-0.116	-0.118	-0.117	-0.113	-0.123	-0.116	-0.118	-0.117	-0.113	-0.123	-0.116	-0.118
Meningioma	Lifestyle	0.269	0.267	0.272	0.270	0.269	-0.040	-0.038	-0.042	-0.028	-0.044	-0.134	-0.106	-0.168	-0.122	-0.139	-0.080	-0.095	-0.062	-0.064	-0.086	0.168	0.161	0.176	0.147	0.175	0.356	0.350	0.363	0.360	0.355
Meningioma	Combined	0.550	0.536	0.567	0.681	0.506	0.363	0.342	0.388	0.560	0.297	0.305	0.299	0.311	0.517	0.234	0.341	0.311	0.378	0.550	0.271	0.492	0.470	0.517	0.635	0.443	0.606	0.590	0.626	0.727	0.565", header = T, sep = "\t")

## combined
data.sjlife <- data

data.sjlife.wanted <- data.sjlife[grepl("Lifestyle", data.sjlife$Variables),]

variables <- c("diet", "drk", "smk", "PA", "obese")


i=1
sjlife <- data.sjlife.wanted[grepl(paste0("SNtypes|",variables[i]), colnames(data.sjlife.wanted))]
View(sjlife)
