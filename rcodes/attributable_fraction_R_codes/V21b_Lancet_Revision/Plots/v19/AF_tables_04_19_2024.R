# Load required libraries
library(ggplot2)
library(ggthemes)  # For additional themes
library(tidyverse)

## V18 b (without lifestyle)
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

data$SN_types_new <- gsub("\\([0-9]+\\)", "", data$SN_types)
data[data == "-"] <- NA


lifestyle <- "without_lifestyle"
all.group <- c("Overall", "Female", "Male", "Age.lt.35", "Age.ge.35")
variables <-  c("Combined", "Radiation", "Chemotherapy", "treatments", "PRS", "Lifestyle")


# Define the text file path
data_name <- paste0("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/ANALYSIS/results/plots/v19b/", lifestyle, "_all_table_19b.txt")

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
SJLIFE	SMNs (234)	Radiation	0.379	0.378	0.379	0.289	0.413
SJLIFE	SMNs (234)	Chemotherapy	0.017	0.018	0.017	0.03	0.013
SJLIFE	SMNs (234)	All_treatments	0.389	0.39	0.388	0.308	0.421
SJLIFE	SMNs (234)	PRS	0.143	0.14	0.146	0.142	0.143
SJLIFE	SMNs (234)	Lifestyle	-0.083	-0.073	-0.096	-0.082	-0.084
SJLIFE	SMNs (234)	Combined	0.431	0.437	0.423	0.354	0.46
SJLIFE	NMSC (118)	Radiation	0.288	0.277	0.301	0.234	0.299
SJLIFE	NMSC (118)	Chemotherapy	-	-	-	-	-
SJLIFE	NMSC (118)	All_treatments	0.288	0.277	0.301	0.234	0.299
SJLIFE	NMSC (118)	PRS	0.307	0.289	0.33	0.304	0.308
SJLIFE	NMSC (118)	Lifestyle	-0.279	-0.255	-0.31	-0.246	-0.287
SJLIFE	NMSC (118)	Combined	0.368	0.356	0.383	0.328	0.377
SJLIFE	Breast cancer (51)	Radiation	0.466	0.466	-	0.391	0.478
SJLIFE	Breast cancer (51)	Chemotherapy	0.106	0.106	-	0.102	0.107
SJLIFE	Breast cancer (51)	All_treatments	0.526	0.526	-	0.452	0.538
SJLIFE	Breast cancer (51)	PRS	0.222	0.222	-	0.25	0.218
SJLIFE	Breast cancer (51)	Lifestyle	-0.872	-0.872	-	-0.581	-0.919
SJLIFE	Breast cancer (51)	Combined	0.335	0.335	-	0.337	0.335
SJLIFE	Thyroid cancer (43)	Radiation	0.61	0.633	0.585	0.555	0.637
SJLIFE	Thyroid cancer (43)	Chemotherapy	0.238	0.237	0.238	0.332	0.192
SJLIFE	Thyroid cancer (43)	All_treatments	0.724	0.74	0.706	0.708	0.732
SJLIFE	Thyroid cancer (43)	PRS	0.589	0.587	0.592	0.601	0.583
SJLIFE	Thyroid cancer (43)	Lifestyle	-0.178	-0.189	-0.166	-0.146	-0.194
SJLIFE	Thyroid cancer (43)	Combined	0.87	0.875	0.863	0.863	0.873
SJLIFE	Meningioma (81)	Radiation	0.184	0.15	0.225	0.245	0.164
SJLIFE	Meningioma (81)	Chemotherapy	0.321	0.325	0.316	0.462	0.274
SJLIFE	Meningioma (81)	All_treatments	0.456	0.436	0.48	0.623	0.4
SJLIFE	Meningioma (81)	PRS	-0.117	-0.113	-0.122	-0.114	-0.118
SJLIFE	Meningioma (81)	Lifestyle	-0.135	-0.126	-0.146	-0.138	-0.134
SJLIFE	Meningioma (81)	Combined	0.308	0.293	0.325	0.509	0.24
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
CCSS	SMNs (619)	Radiation	0.249	0.246	0.255	0.186	0.275
CCSS	SMNs (619)	Chemotherapy	0.056	0.055	0.058	0.099	0.038
CCSS	SMNs (619)	All_treatments	0.297	0.293	0.306	0.271	0.308
CCSS	SMNs (619)	PRS	0.047	0.047	0.048	0.047	0.047
CCSS	SMNs (619)	Lifestyle	-0.048	-0.046	-0.052	-0.046	-0.049
CCSS	SMNs (619)	Combined	0.297	0.295	0.303	0.272	0.308
CCSS	NMSC (632)	Radiation	0.393	0.39	0.396	0.403	0.389
CCSS	NMSC (632)	Chemotherapy	-	-	-	-	-
CCSS	NMSC (632)	All_treatments	0.393	0.39	0.396	0.403	0.389
CCSS	NMSC (632)	PRS	0.313	0.31	0.315	0.314	0.312
CCSS	NMSC (632)	Lifestyle	-0.151	-0.15	-0.153	-0.132	-0.159
CCSS	NMSC (632)	Combined	0.524	0.521	0.527	0.54	0.517
CCSS	Breast cancer (259)	Radiation	0.444	0.444	-	0.392	0.45
CCSS	Breast cancer (259)	Chemotherapy	0.194	0.194	-	0.234	0.189
CCSS	Breast cancer (259)	All_treatments	0.572	0.572	-	0.541	0.575
CCSS	Breast cancer (259)	PRS	0.346	0.346	-	0.348	0.345
CCSS	Breast cancer (259)	Lifestyle	-0.088	-0.088	-	-0.078	-0.089
CCSS	Breast cancer (259)	Combined	0.696	0.696	-	0.678	0.698
CCSS	Thyroid cancer (123)	Radiation	0.416	0.416	0.418	0.36	0.486
CCSS	Thyroid cancer (123)	Chemotherapy	0.095	0.093	0.097	0.122	0.061
CCSS	Thyroid cancer (123)	All_treatments	0.479	0.477	0.482	0.443	0.523
CCSS	Thyroid cancer (123)	PRS	0.386	0.391	0.376	0.379	0.394
CCSS	Thyroid cancer (123)	Lifestyle	-0.042	-0.036	-0.054	-0.025	-0.063
CCSS	Thyroid cancer (123)	Combined	0.666	0.667	0.666	0.646	0.692
CCSS	Meningioma (209)	Radiation	0.294	0.26	0.345	0.315	0.279
CCSS	Meningioma (209)	Chemotherapy	0.087	0.08	0.098	0.134	0.055
CCSS	Meningioma (209)	All_treatments	0.368	0.334	0.42	0.421	0.331
CCSS	Meningioma (209)	PRS	0.028	0.027	0.029	0.033	0.024
CCSS	Meningioma (209)	Lifestyle	0.087	0.099	0.068	0.15	0.043
CCSS	Meningioma (209)	Combined	0.433	0.409	0.469	0.524	0.37
CCSS	Sarcoma (41)	Radiation	-	-	-	-	-
CCSS	Sarcoma (41)	Chemotherapy	0.279	0.266	0.3	0.279	0.278
CCSS	Sarcoma (41)	All_treatments	0.279	0.266	0.3	0.279	0.278
CCSS	Sarcoma (41)	PRS	0.035	0.038	0.03	0.033	0.036
CCSS	Sarcoma (41)	Lifestyle	-0.261	-0.256	-0.271	-0.396	-0.18
CCSS	Sarcoma (41)	Combined	0.116	0.108	0.13	0.017	0.176", header = T, sep = "\t")

data$SN_types_new <- gsub("\\([0-9]+\\)", "", data$SN_types)
data[data == "-"] <- NA


lifestyle <- "without_diet"
all.group <- c("Overall", "Female", "Male", "Age.lt.35", "Age.ge.35")
variables <-  c("Combined", "Radiation", "Chemotherapy", "treatments", "PRS", "Lifestyle")


# Define the text file path
data_name <- paste0("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/ANALYSIS/results/plots/v19b/", lifestyle, "_all_table_19b.txt")

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
SJLIFE	SMNs (234)	Chemotherapy	0.019	0.019	0.018	0.032	0.013
SJLIFE	SMNs (234)	All_treatments	0.39	0.39	0.389	0.31	0.421
SJLIFE	SMNs (234)	PRS	0.141	0.138	0.144	0.14	0.141
SJLIFE	SMNs (234)	Lifestyle	-0.13	-0.122	-0.139	-0.133	-0.128
SJLIFE	SMNs (234)	Combined	0.405	0.41	0.398	0.322	0.437
SJLIFE	NMSC (118)	Radiation	0.282	0.271	0.296	0.229	0.294
SJLIFE	NMSC (118)	Chemotherapy	-	-	-	-	-
SJLIFE	NMSC (118)	All_treatments	0.282	0.271	0.296	0.229	0.294
SJLIFE	NMSC (118)	PRS	0.323	0.306	0.344	0.322	0.323
SJLIFE	NMSC (118)	Lifestyle	-0.541	-0.453	-0.653	-0.536	-0.542
SJLIFE	NMSC (118)	Combined	0.25	0.264	0.233	0.187	0.264
SJLIFE	Breast cancer (51)	Radiation	0.47	0.47	-	0.399	0.481
SJLIFE	Breast cancer (51)	Chemotherapy	0.109	0.109	-	0.104	0.109
SJLIFE	Breast cancer (51)	All_treatments	0.531	0.531	-	0.461	0.542
SJLIFE	Breast cancer (51)	PRS	0.229	0.229	-	0.262	0.224
SJLIFE	Breast cancer (51)	Lifestyle	-0.908	-0.908	-	-0.589	-0.959
SJLIFE	Breast cancer (51)	Combined	0.345	0.345	-	0.353	0.344
SJLIFE	Thyroid cancer (43)	Radiation	0.611	0.633	0.586	0.551	0.641
SJLIFE	Thyroid cancer (43)	Chemotherapy	0.238	0.237	0.24	0.326	0.196
SJLIFE	Thyroid cancer (43)	All_treatments	0.726	0.741	0.708	0.703	0.737
SJLIFE	Thyroid cancer (43)	PRS	0.582	0.58	0.585	0.595	0.576
SJLIFE	Thyroid cancer (43)	Lifestyle	-0.691	-0.576	-0.824	-0.722	-0.676
SJLIFE	Thyroid cancer (43)	Combined	0.809	0.832	0.783	0.79	0.818
SJLIFE	Meningioma (81)	Radiation	0.184	0.15	0.225	0.244	0.164
SJLIFE	Meningioma (81)	Chemotherapy	0.321	0.325	0.316	0.463	0.274
SJLIFE	Meningioma (81)	All_treatments	0.456	0.436	0.48	0.623	0.4
SJLIFE	Meningioma (81)	PRS	-0.117	-0.112	-0.122	-0.114	-0.117
SJLIFE	Meningioma (81)	Lifestyle	-0.121	-0.111	-0.132	-0.12	-0.121
SJLIFE	Meningioma (81)	Combined	0.316	0.302	0.332	0.517	0.249", header = T, sep = "\t")


data$SN_types_new <- gsub("\\([0-9]+\\)", "", data$SN_types)
data[data == "-"] <- NA


lifestyle <- "lifestyle_with_diet_HEI2015"
all.group <- c("Overall", "Female", "Male", "Age.lt.35", "Age.ge.35")
variables <-  c("Combined", "Radiation", "Chemotherapy", "treatments", "PRS", "Lifestyle")


# Define the text file path
data_name <- paste0("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/ANALYSIS/results/plots/v19b/", lifestyle, "_all_table_19b.txt")

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
SJLIFE	SMNs	Radiation	0.379	0.378	0.379	0.289	0.413	0.379	0.378	0.379	0.289	0.413	0.379	0.378	0.379	0.289	0.413	0.379	0.378	0.379	0.289	0.413	0.379	0.378	0.379	0.289	0.413
SJLIFE	SMNs	Chemotherapy	0.017	0.018	0.017	0.03	0.013	0.017	0.018	0.017	0.03	0.013	0.017	0.018	0.017	0.03	0.013	0.017	0.018	0.017	0.03	0.013	0.017	0.018	0.017	0.03	0.013
SJLIFE	SMNs	All_treatments	0.389	0.39	0.388	0.308	0.421	0.389	0.39	0.388	0.308	0.421	0.389	0.39	0.388	0.308	0.421	0.389	0.39	0.388	0.308	0.421	0.389	0.39	0.388	0.308	0.421
SJLIFE	SMNs	PRS	0.143	0.14	0.146	0.142	0.143	0.143	0.14	0.146	0.142	0.143	0.143	0.14	0.146	0.142	0.143	0.143	0.14	0.146	0.142	0.143	0.143	0.14	0.146	0.142	0.143
SJLIFE	SMNs	Lifestyle	-0.083	-0.073	-0.096	-0.082	-0.084	-0.079	-0.081	-0.077	-0.082	-0.078	-0.051	-0.051	-0.052	-0.059	-0.048	0.033	0.043	0.021	0.017	0.04	-0.111	-0.114	-0.109	-0.105	-0.114
SJLIFE	SMNs	Combined	0.431	0.437	0.423	0.354	0.46	0.434	0.433	0.435	0.355	0.465	0.449	0.449	0.45	0.37	0.48	0.494	0.498	0.488	0.415	0.524	0.417	0.417	0.418	0.343	0.446
SJLIFE	NMSC	Radiation	0.288	0.277	0.301	0.234	0.299	0.288	0.277	0.301	0.234	0.299	0.288	0.277	0.301	0.234	0.299	0.288	0.277	0.301	0.234	0.299	0.288	0.277	0.301	0.234	0.299
SJLIFE	NMSC	Chemotherapy	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-
SJLIFE	NMSC	All_treatments	0.288	0.277	0.301	0.234	0.299	0.288	0.277	0.301	0.234	0.299	0.288	0.277	0.301	0.234	0.299	0.288	0.277	0.301	0.234	0.299	0.288	0.277	0.301	0.234	0.299
SJLIFE	NMSC	PRS	0.307	0.289	0.33	0.304	0.308	0.307	0.289	0.33	0.304	0.308	0.307	0.289	0.33	0.304	0.308	0.307	0.289	0.33	0.304	0.308	0.307	0.289	0.33	0.304	0.308
SJLIFE	NMSC	Lifestyle	-0.279	-0.255	-0.31	-0.246	-0.287	-0.028	-0.03	-0.025	-0.044	-0.024	-0.156	-0.137	-0.181	-0.152	-0.157	-0.044	-0.047	-0.04	-0.059	-0.041	-0.16	-0.16	-0.161	-0.147	-0.163
SJLIFE	NMSC	Combined	0.368	0.356	0.383	0.328	0.377	0.493	0.473	0.518	0.443	0.504	0.368	0.356	0.383	0.328	0.377	0.485	0.464	0.511	0.434	0.496	0.427	0.407	0.454	0.387	0.436
SJLIFE	Breast_cancer	Radiation	0.466	0.466	-	0.391	0.478	0.466	0.466	-	0.391	0.478	0.466	0.466	-	0.391	0.478	0.466	0.466	-	0.391	0.478	0.466	0.466	-	0.391	0.478
SJLIFE	Breast_cancer	Chemotherapy	0.106	0.106	-	0.102	0.107	0.106	0.106	-	0.102	0.107	0.106	0.106	-	0.102	0.107	0.106	0.106	-	0.102	0.107	0.106	0.106	-	0.102	0.107
SJLIFE	Breast_cancer	All_treatments	0.526	0.526	-	0.452	0.538	0.526	0.526	-	0.452	0.538	0.526	0.526	-	0.452	0.538	0.526	0.526	-	0.452	0.538	0.526	0.526	-	0.452	0.538
SJLIFE	Breast_cancer	PRS	0.222	0.222	-	0.25	0.218	0.222	0.222	-	0.25	0.218	0.222	0.222	-	0.25	0.218	0.222	0.222	-	0.25	0.218	0.222	0.222	-	0.25	0.218
SJLIFE	Breast_cancer	Lifestyle	-0.872	-0.872	-	-0.581	-0.919	-0.256	-0.256	-	-0.159	-0.272	0.036	0.036	-	0.034	0.037	-0.375	-0.375	-	-0.224	-0.399	0.004	0.004	-	0.005	0.004
SJLIFE	Breast_cancer	Combined	0.335	0.335	-	0.337	0.335	0.553	0.553	-	0.525	0.557	0.649	0.649	-	0.608	0.656	0.512	0.512	-	0.494	0.515	0.635	0.635	-	0.597	0.641
SJLIFE	Thyroid_cancer	Radiation	0.61	0.633	0.585	0.555	0.637	0.61	0.633	0.585	0.555	0.637	0.61	0.633	0.585	0.555	0.637	0.61	0.633	0.585	0.555	0.637	0.61	0.633	0.585	0.555	0.637
SJLIFE	Thyroid_cancer	Chemotherapy	0.238	0.237	0.238	0.332	0.192	0.238	0.237	0.238	0.332	0.192	0.238	0.237	0.238	0.332	0.192	0.238	0.237	0.238	0.332	0.192	0.238	0.237	0.238	0.332	0.192
SJLIFE	Thyroid_cancer	All_treatments	0.724	0.74	0.706	0.708	0.732	0.724	0.74	0.706	0.708	0.732	0.724	0.74	0.706	0.708	0.732	0.724	0.74	0.706	0.708	0.732	0.724	0.74	0.706	0.708	0.732
SJLIFE	Thyroid_cancer	PRS	0.589	0.587	0.592	0.601	0.583	0.589	0.587	0.592	0.601	0.583	0.589	0.587	0.592	0.601	0.583	0.589	0.587	0.592	0.601	0.583	0.589	0.587	0.592	0.601	0.583
SJLIFE	Thyroid_cancer	Lifestyle	-0.178	-0.189	-0.166	-0.146	-0.194	0.093	0.099	0.085	0.086	0.096	-0.149	-0.134	-0.166	-0.151	-0.148	-0.189	-0.216	-0.157	-0.163	-0.201	-0.101	-0.102	-0.1	-0.092	-0.105
SJLIFE	Thyroid_cancer	Combined	0.87	0.875	0.863	0.863	0.873	0.899	0.905	0.894	0.891	0.904	0.87	0.878	0.861	0.86	0.875	0.867	0.87	0.863	0.859	0.871	0.877	0.883	0.87	0.869	0.881
SJLIFE	Meningioma	Radiation	0.184	0.15	0.225	0.245	0.164	0.184	0.15	0.225	0.245	0.164	0.184	0.15	0.225	0.245	0.164	0.184	0.15	0.225	0.245	0.164	0.184	0.15	0.225	0.245	0.164
SJLIFE	Meningioma	Chemotherapy	0.321	0.325	0.316	0.462	0.274	0.321	0.325	0.316	0.462	0.274	0.321	0.325	0.316	0.462	0.274	0.321	0.325	0.316	0.462	0.274	0.321	0.325	0.316	0.462	0.274
SJLIFE	Meningioma	All_treatments	0.456	0.436	0.48	0.623	0.4	0.456	0.436	0.48	0.623	0.4	0.456	0.436	0.48	0.623	0.4	0.456	0.436	0.48	0.623	0.4	0.456	0.436	0.48	0.623	0.4
SJLIFE	Meningioma	PRS	-0.117	-0.113	-0.122	-0.114	-0.118	-0.117	-0.113	-0.122	-0.114	-0.118	-0.117	-0.113	-0.122	-0.114	-0.118	-0.117	-0.113	-0.122	-0.114	-0.118	-0.117	-0.113	-0.122	-0.114	-0.118
SJLIFE	Meningioma	Lifestyle	-0.135	-0.126	-0.146	-0.138	-0.134	-0.053	-0.051	-0.055	-0.046	-0.055	-0.157	-0.129	-0.191	-0.147	-0.161	-0.086	-0.1	-0.07	-0.076	-0.089	0.149	0.143	0.157	0.124	0.158
SJLIFE	Meningioma	Combined	0.308	0.293	0.325	0.509	0.24	0.361	0.34	0.385	0.557	0.295	0.296	0.291	0.302	0.51	0.224	0.343	0.314	0.377	0.549	0.274	0.484	0.464	0.509	0.629	0.436
CCSS	SNs	Radiation	0.358	0.351	0.369	0.338	0.368	0.358	0.351	0.369	0.338	0.368	0.358	0.351	0.369	0.338	0.368	0.358	0.351	0.369	0.338	0.368	0.358	0.351	0.369	0.338	0.368
CCSS	SNs	Chemotherapy	0.038	0.037	0.041	0.061	0.027	0.038	0.037	0.041	0.061	0.027	0.038	0.037	0.041	0.061	0.027	0.038	0.037	0.041	0.061	0.027	0.038	0.037	0.041	0.061	0.027
CCSS	SNs	All_treatments	0.384	0.376	0.396	0.38	0.387	0.384	0.376	0.396	0.38	0.387	0.384	0.376	0.396	0.38	0.387	0.384	0.376	0.396	0.38	0.387	0.384	0.376	0.396	0.38	0.387
CCSS	SNs	PRS	0.057	0.056	0.058	0.056	0.057	0.057	0.056	0.058	0.056	0.057	0.057	0.056	0.058	0.056	0.057	0.057	0.056	0.058	0.056	0.057	0.057	0.056	0.058	0.056	0.057
CCSS	SNs	Lifestyle	-0.074	-0.071	-0.078	-0.068	-0.077	-0.052	-0.048	-0.057	-0.045	-0.055	-0.003	-0.002	-0.003	-0.002	-0.003	0.002	0.002	0.002	0.003	0.002	-0.017	-0.018	-0.015	-0.016	-0.017
CCSS	SNs	Combined	0.376	0.37	0.386	0.374	0.378	0.389	0.383	0.397	0.386	0.39	0.418	0.41	0.429	0.413	0.42	0.42	0.413	0.432	0.416	0.422	0.41	0.402	0.422	0.406	0.412
CCSS	SMNs	Radiation	0.249	0.246	0.255	0.186	0.275	0.249	0.246	0.255	0.186	0.275	0.249	0.246	0.255	0.186	0.275	0.249	0.246	0.255	0.186	0.275	0.249	0.246	0.255	0.186	0.275
CCSS	SMNs	Chemotherapy	0.056	0.055	0.058	0.099	0.038	0.056	0.055	0.058	0.099	0.038	0.056	0.055	0.058	0.099	0.038	0.056	0.055	0.058	0.099	0.038	0.056	0.055	0.058	0.099	0.038
CCSS	SMNs	All_treatments	0.297	0.293	0.306	0.271	0.308	0.297	0.293	0.306	0.271	0.308	0.297	0.293	0.306	0.271	0.308	0.297	0.293	0.306	0.271	0.308	0.297	0.293	0.306	0.271	0.308
CCSS	SMNs	PRS	0.047	0.047	0.048	0.047	0.047	0.047	0.047	0.048	0.047	0.047	0.047	0.047	0.048	0.047	0.047	0.047	0.047	0.048	0.047	0.047	0.047	0.047	0.048	0.047	0.047
CCSS	SMNs	Lifestyle	-0.048	-0.046	-0.052	-0.046	-0.049	-0.037	-0.035	-0.04	-0.034	-0.038	-0.004	-0.004	-0.005	-0.004	-0.004	-0.002	-0.002	-0.002	-0.001	-0.002	-0.002	-0.002	-0.002	-0.001	-0.002
CCSS	SMNs	Combined	0.297	0.295	0.303	0.272	0.308	0.305	0.302	0.311	0.281	0.315	0.327	0.323	0.335	0.302	0.338	0.329	0.324	0.337	0.304	0.339	0.329	0.324	0.337	0.304	0.339
CCSS	NMSC	Radiation	0.393	0.39	0.396	0.403	0.389	0.393	0.39	0.396	0.403	0.389	0.393	0.39	0.396	0.403	0.389	0.393	0.39	0.396	0.403	0.389	0.393	0.39	0.396	0.403	0.389
CCSS	NMSC	Chemotherapy	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-
CCSS	NMSC	All_treatments	0.393	0.39	0.396	0.403	0.389	0.393	0.39	0.396	0.403	0.389	0.393	0.39	0.396	0.403	0.389	0.393	0.39	0.396	0.403	0.389	0.393	0.39	0.396	0.403	0.389
CCSS	NMSC	PRS	0.313	0.31	0.315	0.314	0.312	0.313	0.31	0.315	0.314	0.312	0.313	0.31	0.315	0.314	0.312	0.313	0.31	0.315	0.314	0.312	0.313	0.31	0.315	0.314	0.312
CCSS	NMSC	Lifestyle	-0.151	-0.15	-0.153	-0.132	-0.159	-0.07	-0.06	-0.081	-0.056	-0.075	0.021	0.02	0.022	0.028	0.019	-0.034	-0.035	-0.032	-0.026	-0.036	-0.043	-0.047	-0.038	-0.037	-0.046
CCSS	NMSC	Combined	0.524	0.521	0.527	0.54	0.517	0.555	0.555	0.554	0.568	0.55	0.593	0.589	0.597	0.605	0.588	0.57	0.566	0.574	0.582	0.565	0.568	0.564	0.572	0.579	0.563
CCSS	Breast_cancer	Radiation	0.444	0.444	-	0.392	0.45	0.444	0.444	-	0.392	0.45	0.444	0.444	-	0.392	0.45	0.444	0.444	-	0.392	0.45	0.444	0.444	-	0.392	0.45
CCSS	Breast_cancer	Chemotherapy	0.194	0.194	-	0.234	0.189	0.194	0.194	-	0.234	0.189	0.194	0.194	-	0.234	0.189	0.194	0.194	-	0.234	0.189	0.194	0.194	-	0.234	0.189
CCSS	Breast_cancer	All_treatments	0.572	0.572	-	0.541	0.575	0.572	0.572	-	0.541	0.575	0.572	0.572	-	0.541	0.575	0.572	0.572	-	0.541	0.575	0.572	0.572	-	0.541	0.575
CCSS	Breast_cancer	PRS	0.346	0.346	-	0.348	0.345	0.346	0.346	-	0.348	0.345	0.346	0.346	-	0.348	0.345	0.346	0.346	-	0.348	0.345	0.346	0.346	-	0.348	0.345
CCSS	Breast_cancer	Lifestyle	-0.088	-0.088	-	-0.078	-0.089	-0.039	-0.039	-	-0.035	-0.039	0.021	0.021	-	0.029	0.02	-0.029	-0.029	-	-0.024	-0.029	-0.026	-0.026	-	-0.024	-0.026
CCSS	Breast_cancer	Combined	0.696	0.696	-	0.678	0.698	0.71	0.71	-	0.692	0.712	0.727	0.727	-	0.711	0.729	0.714	0.714	-	0.695	0.716	0.714	0.714	-	0.695	0.716
CCSS	Thyroid_cancer	Radiation	0.416	0.416	0.418	0.36	0.486	0.416	0.416	0.418	0.36	0.486	0.416	0.416	0.418	0.36	0.486	0.416	0.416	0.418	0.36	0.486	0.416	0.416	0.418	0.36	0.486
CCSS	Thyroid_cancer	Chemotherapy	0.095	0.093	0.097	0.122	0.061	0.095	0.093	0.097	0.122	0.061	0.095	0.093	0.097	0.122	0.061	0.095	0.093	0.097	0.122	0.061	0.095	0.093	0.097	0.122	0.061
CCSS	Thyroid_cancer	All_treatments	0.479	0.477	0.482	0.443	0.523	0.479	0.477	0.482	0.443	0.523	0.479	0.477	0.482	0.443	0.523	0.479	0.477	0.482	0.443	0.523	0.479	0.477	0.482	0.443	0.523
CCSS	Thyroid_cancer	PRS	0.386	0.391	0.376	0.379	0.394	0.386	0.391	0.376	0.379	0.394	0.386	0.391	0.376	0.379	0.394	0.386	0.391	0.376	0.379	0.394	0.386	0.391	0.376	0.379	0.394
CCSS	Thyroid_cancer	Lifestyle	-0.042	-0.036	-0.054	-0.025	-0.063	-0.037	-0.028	-0.053	-0.023	-0.053	0.062	0.061	0.064	0.075	0.046	0.019	0.02	0.017	0.029	0.007	0.005	0.005	0.004	0.014	-0.006
CCSS	Thyroid_cancer	Combined	0.666	0.667	0.666	0.646	0.692	0.668	0.67	0.664	0.646	0.694	0.701	0.7	0.702	0.682	0.724	0.687	0.686	0.687	0.666	0.712	0.682	0.682	0.683	0.661	0.708
CCSS	Meningioma	Radiation	0.294	0.26	0.345	0.315	0.279	0.294	0.26	0.345	0.315	0.279	0.294	0.26	0.345	0.315	0.279	0.294	0.26	0.345	0.315	0.279	0.294	0.26	0.345	0.315	0.279
CCSS	Meningioma	Chemotherapy	0.087	0.08	0.098	0.134	0.055	0.087	0.08	0.098	0.134	0.055	0.087	0.08	0.098	0.134	0.055	0.087	0.08	0.098	0.134	0.055	0.087	0.08	0.098	0.134	0.055
CCSS	Meningioma	All_treatments	0.368	0.334	0.42	0.421	0.331	0.368	0.334	0.42	0.421	0.331	0.368	0.334	0.42	0.421	0.331	0.368	0.334	0.42	0.421	0.331	0.368	0.334	0.42	0.421	0.331
CCSS	Meningioma	PRS	0.028	0.027	0.029	0.033	0.024	0.028	0.027	0.029	0.033	0.024	0.028	0.027	0.029	0.033	0.024	0.028	0.027	0.029	0.033	0.024	0.028	0.027	0.029	0.033	0.024
CCSS	Meningioma	Lifestyle	0.087	0.099	0.068	0.15	0.043	0.093	0.099	0.083	0.156	0.049	0.056	0.065	0.042	0.119	0.012	0.121	0.128	0.112	0.181	0.08	0.127	0.135	0.116	0.189	0.084
CCSS	Meningioma	Combined	0.433	0.409	0.469	0.524	0.37	0.442	0.415	0.484	0.533	0.38	0.417	0.39	0.457	0.51	0.352	0.461	0.434	0.501	0.547	0.4	0.463	0.437	0.503	0.55	0.402
CCSS	Sarcoma	Radiation	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-
CCSS	Sarcoma	Chemotherapy	0.279	0.266	0.3	0.279	0.278	0.279	0.266	0.3	0.279	0.278	0.279	0.266	0.3	0.279	0.278	0.279	0.266	0.3	0.279	0.278	0.279	0.266	0.3	0.279	0.278
CCSS	Sarcoma	All_treatments	0.279	0.266	0.3	0.279	0.278	0.279	0.266	0.3	0.279	0.278	0.279	0.266	0.3	0.279	0.278	0.279	0.266	0.3	0.279	0.278	0.279	0.266	0.3	0.279	0.278
CCSS	Sarcoma	PRS	0.035	0.038	0.03	0.033	0.036	0.035	0.038	0.03	0.033	0.036	0.035	0.038	0.03	0.033	0.036	0.035	0.038	0.03	0.033	0.036	0.035	0.038	0.03	0.033	0.036
CCSS	Sarcoma	Lifestyle	-0.261	-0.256	-0.271	-0.396	-0.18	-0.121	-0.126	-0.112	-0.215	-0.064	-0.229	-0.225	-0.235	-0.36	-0.149	-0.09	-0.096	-0.08	-0.183	-0.033	-0.129	-0.136	-0.116	-0.226	-0.071
CCSS	Sarcoma	Combined	0.116	0.108	0.13	0.017	0.176	0.213	0.199	0.238	0.144	0.255	0.14	0.131	0.155	0.043	0.198	0.116	0.108	0.13	0.017	0.176	0.208	0.192	0.235	0.137	0.251", header = T, sep = "\t")

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
###########################################################
## Lifestyle with diet HEI (SJLIFE only)
data <- read.table(text="SNtypes	Variables	Overall_all4	Female_all4	Male_all4	Age_lt_35_all4	Age_ge_35_all4	Overall_smk	Female_smk	Male_smk	Age_lt_35_smk	Age_ge_35_smk	Overall_drk	Female_drk	Male_drk	Age_lt_35_drk	Age_ge_35_drk	Overall_PA	Female_PA	Male_PA	Age_lt_35_PA	Age_ge_35_PA	Overall_obese	Female_obese	Male_obese	Age_lt_35_obese	Age_ge_35_obese	Overall_diet	Female_diet	Male_diet	Age_lt_35_diet	Age_ge_35_diet
Any_SN	Radiation	0.443	0.433	0.457	0.369	0.473	0.443	0.433	0.457	0.369	0.473	0.443	0.433	0.457	0.369	0.473	0.443	0.433	0.457	0.369	0.473	0.443	0.433	0.457	0.369	0.473	0.443	0.433	0.457	0.369	0.473
Any_SN	Chemotherapy	0.020	0.020	0.019	0.036	0.013	0.020	0.020	0.019	0.036	0.013	0.020	0.020	0.019	0.036	0.013	0.020	0.020	0.019	0.036	0.013	0.020	0.020	0.019	0.036	0.013	0.020	0.020	0.019	0.036	0.013
Any_SN	All_treatments	0.453	0.443	0.465	0.389	0.479	0.453	0.443	0.465	0.389	0.479	0.453	0.443	0.465	0.389	0.479	0.453	0.443	0.465	0.389	0.479	0.453	0.443	0.465	0.389	0.479	0.453	0.443	0.465	0.389	0.479
Any_SN	PRS	0.180	0.178	0.183	0.183	0.179	0.180	0.178	0.183	0.183	0.179	0.180	0.178	0.183	0.183	0.179	0.180	0.178	0.183	0.183	0.179	0.180	0.178	0.183	0.183	0.179	0.180	0.178	0.183	0.183	0.179
Any_SN	Lifestyle	-0.057	-0.062	-0.049	-0.050	-0.060	-0.059	-0.061	-0.055	-0.052	-0.062	-0.030	-0.028	-0.033	-0.032	-0.029	0.014	0.017	0.010	0.010	0.016	-0.020	-0.021	-0.020	-0.021	-0.020	-0.012	-0.022	0.001	-0.012	-0.012
Any_SN	Combined	0.525	0.515	0.538	0.470	0.548	0.525	0.516	0.537	0.471	0.547	0.538	0.531	0.548	0.483	0.561	0.558	0.552	0.567	0.504	0.581	0.543	0.535	0.554	0.489	0.565	0.547	0.535	0.563	0.493	0.569
Any_SMN	Radiation	0.378	0.378	0.379	0.289	0.413	0.378	0.378	0.379	0.289	0.413	0.378	0.378	0.379	0.289	0.413	0.378	0.378	0.379	0.289	0.413	0.378	0.378	0.379	0.289	0.413	0.378	0.378	0.379	0.289	0.413
Any_SMN	Chemotherapy	0.019	0.019	0.018	0.032	0.013	0.019	0.019	0.018	0.032	0.013	0.019	0.019	0.018	0.032	0.013	0.019	0.019	0.018	0.032	0.013	0.019	0.019	0.018	0.032	0.013	0.019	0.019	0.018	0.032	0.013
Any_SMN	All_treatments	0.390	0.390	0.389	0.310	0.421	0.390	0.390	0.389	0.310	0.421	0.390	0.390	0.389	0.310	0.421	0.390	0.390	0.389	0.310	0.421	0.390	0.390	0.389	0.310	0.421	0.390	0.390	0.389	0.310	0.421
Any_SMN	PRS	0.141	0.138	0.144	0.140	0.141	0.141	0.138	0.144	0.140	0.141	0.141	0.138	0.144	0.140	0.141	0.141	0.138	0.144	0.140	0.141	0.141	0.138	0.144	0.140	0.141	0.141	0.138	0.144	0.140	0.141
Any_SMN	Lifestyle	-0.130	-0.122	-0.139	-0.133	-0.128	-0.103	-0.107	-0.097	-0.108	-0.101	-0.075	-0.076	-0.074	-0.085	-0.071	0.015	0.023	0.005	-0.004	0.022	-0.131	-0.135	-0.126	-0.127	-0.132	-0.082	-0.086	-0.076	-0.095	-0.077
Any_SMN	Combined	0.405	0.410	0.398	0.322	0.437	0.421	0.419	0.423	0.339	0.452	0.436	0.435	0.437	0.354	0.467	0.483	0.486	0.479	0.402	0.514	0.406	0.405	0.407	0.329	0.436	0.432	0.429	0.436	0.348	0.464
NMSC	Radiation	0.290	0.278	0.306	0.235	0.302	0.290	0.278	0.306	0.235	0.302	0.290	0.278	0.306	0.235	0.302	0.290	0.278	0.306	0.235	0.302	0.290	0.278	0.306	0.235	0.302	0.290	0.278	0.306	0.235	0.302
NMSC	Chemotherapy	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-
NMSC	All_treatments	0.290	0.278	0.306	0.235	0.302	0.290	0.278	0.306	0.235	0.302	0.290	0.278	0.306	0.235	0.302	0.290	0.278	0.306	0.235	0.302	0.290	0.278	0.306	0.235	0.302	0.290	0.278	0.306	0.235	0.302
NMSC	PRS	0.441	0.429	0.457	0.441	0.441	0.441	0.429	0.457	0.441	0.441	0.441	0.429	0.457	0.441	0.441	0.441	0.429	0.457	0.441	0.441	0.441	0.429	0.457	0.441	0.441	0.441	0.429	0.457	0.441	0.441
NMSC	Lifestyle	-0.528	-0.440	-0.640	-0.525	-0.529	-0.044	-0.047	-0.039	-0.070	-0.038	-0.198	-0.177	-0.225	-0.202	-0.197	-0.049	-0.052	-0.047	-0.076	-0.043	-0.172	-0.174	-0.171	-0.176	-0.172	-0.274	-0.233	-0.326	-0.319	-0.264
NMSC	Combined	0.395	0.409	0.378	0.344	0.407	0.587	0.572	0.606	0.543	0.597	0.526	0.518	0.536	0.481	0.536	0.584	0.569	0.603	0.539	0.594	0.535	0.520	0.555	0.497	0.544	0.496	0.495	0.498	0.437	0.509
Breast_cancer	Radiation	0.470	0.470	-	0.399	0.481	0.470	0.470	-	0.399	0.481	0.470	0.470	-	0.399	0.481	0.470	0.470	-	0.399	0.481	0.470	0.470	-	0.399	0.481	0.470	0.470	-	0.399	0.481
Breast_cancer	Chemotherapy	0.109	0.109	-	0.104	0.109	0.109	0.109	-	0.104	0.109	0.109	0.109	-	0.104	0.109	0.109	0.109	-	0.104	0.109	0.109	0.109	-	0.104	0.109	0.109	0.109	-	0.104	0.109
Breast_cancer	All_treatments	0.531	0.531	-	0.461	0.542	0.531	0.531	-	0.461	0.542	0.531	0.531	-	0.461	0.542	0.531	0.531	-	0.461	0.542	0.531	0.531	-	0.461	0.542	0.531	0.531	-	0.461	0.542
Breast_cancer	PRS	0.229	0.229	-	0.262	0.224	0.229	0.229	-	0.262	0.224	0.229	0.229	-	0.262	0.224	0.229	0.229	-	0.262	0.224	0.229	0.229	-	0.262	0.224	0.229	0.229	-	0.262	0.224
Breast_cancer	Lifestyle	-0.908	-0.908	-	-0.589	-0.959	-0.289	-0.289	-	-0.185	-0.306	0.004	0.004	-	0.004	0.004	-0.412	-0.412	-	-0.255	-0.438	-0.042	-0.042	-	-0.038	-0.042	0.043	0.043	-	0.052	0.041
Breast_cancer	Combined	0.345	0.345	-	0.353	0.344	0.552	0.552	-	0.527	0.556	0.646	0.646	-	0.609	0.653	0.511	0.511	-	0.498	0.514	0.627	0.627	-	0.593	0.633	0.664	0.664	-	0.626	0.670
Thyroid_cancer	Radiation	0.614	0.625	0.600	0.547	0.646	0.614	0.625	0.600	0.547	0.646	0.614	0.625	0.600	0.547	0.646	0.614	0.625	0.600	0.547	0.646	0.614	0.625	0.600	0.547	0.646	0.614	0.625	0.600	0.547	0.646
Thyroid_cancer	Chemotherapy	0.246	0.250	0.242	0.336	0.202	0.246	0.250	0.242	0.336	0.202	0.246	0.250	0.242	0.336	0.202	0.246	0.250	0.242	0.336	0.202	0.246	0.250	0.242	0.336	0.202	0.246	0.250	0.242	0.336	0.202
Thyroid_cancer	All_treatments	0.732	0.741	0.722	0.707	0.744	0.732	0.741	0.722	0.707	0.744	0.732	0.741	0.722	0.707	0.744	0.732	0.741	0.722	0.707	0.744	0.732	0.741	0.722	0.707	0.744	0.732	0.741	0.722	0.707	0.744
Thyroid_cancer	PRS	0.592	0.587	0.599	0.605	0.586	0.592	0.587	0.599	0.605	0.586	0.592	0.587	0.599	0.605	0.586	0.592	0.587	0.599	0.605	0.586	0.592	0.587	0.599	0.605	0.586	0.592	0.587	0.599	0.605	0.586
Thyroid_cancer	Lifestyle	-0.732	-0.634	-0.844	-0.780	-0.708	0.065	0.068	0.063	0.050	0.073	-0.260	-0.238	-0.286	-0.268	-0.256	-0.286	-0.316	-0.253	-0.273	-0.293	-0.200	-0.203	-0.195	-0.204	-0.198	-0.604	-0.495	-0.730	-0.688	-0.563
Thyroid_cancer	Combined	0.815	0.830	0.797	0.793	0.826	0.901	0.903	0.899	0.889	0.907	0.863	0.869	0.857	0.848	0.870	0.862	0.863	0.861	0.849	0.869	0.871	0.874	0.868	0.857	0.878	0.825	0.841	0.806	0.800	0.837
Meningioma	Radiation	0.185	0.150	0.226	0.245	0.165	0.185	0.150	0.226	0.245	0.165	0.185	0.150	0.226	0.245	0.165	0.185	0.150	0.226	0.245	0.165	0.185	0.150	0.226	0.245	0.165	0.185	0.150	0.226	0.245	0.165
Meningioma	Chemotherapy	0.324	0.329	0.319	0.470	0.276	0.324	0.329	0.319	0.470	0.276	0.324	0.329	0.319	0.470	0.276	0.324	0.329	0.319	0.470	0.276	0.324	0.329	0.319	0.470	0.276	0.324	0.329	0.319	0.470	0.276
Meningioma	All_treatments	0.459	0.440	0.482	0.629	0.403	0.459	0.440	0.482	0.629	0.403	0.459	0.440	0.482	0.629	0.403	0.459	0.440	0.482	0.629	0.403	0.459	0.440	0.482	0.629	0.403	0.459	0.440	0.482	0.629	0.403
Meningioma	PRS	-0.118	-0.114	-0.123	-0.117	-0.119	-0.118	-0.114	-0.123	-0.117	-0.119	-0.118	-0.114	-0.123	-0.117	-0.119	-0.118	-0.114	-0.123	-0.117	-0.119	-0.118	-0.114	-0.123	-0.117	-0.119	-0.118	-0.114	-0.123	-0.117	-0.119
Meningioma	Lifestyle	-0.127	-0.116	-0.140	-0.125	-0.127	-0.040	-0.038	-0.042	-0.029	-0.043	-0.138	-0.110	-0.171	-0.124	-0.143	-0.084	-0.098	-0.067	-0.069	-0.089	0.163	0.158	0.169	0.141	0.170	0.001	0.002	-0.001	0.005	0.000
Meningioma	Combined	0.316	0.303	0.331	0.522	0.247	0.370	0.351	0.393	0.569	0.304	0.310	0.306	0.315	0.526	0.238	0.346	0.318	0.380	0.557	0.276	0.494	0.475	0.517	0.641	0.445	0.396	0.378	0.418	0.586	0.333", header = T, sep = "\t")

## combined
data.sjlife <- data

data.sjlife.wanted <- data.sjlife[grepl("Lifestyle", data.sjlife$Variables),]

variables <- c("diet", "drk", "smk", "PA", "obese")


i=1
sjlife <- data.sjlife.wanted[grepl(paste0("SNtypes|",variables[i]), colnames(data.sjlife.wanted))]
View(sjlife)
