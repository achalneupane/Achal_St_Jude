# V17_without_lifestyle_variable
mydf <- read.table(text = "Cohort	SN_types	Variables	Overall	Female	Male	Age<35	Age≥35
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
SJLIFE	Meningioma (149)	PRS	NA0.125	NA0.116	NA0.137	NA0.13	NA0.124
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
", header = T, sep = "\t", check.names = F)




# df <-  mydf
# df_subset <- df[df$SN_Types == "Any SN (303)",]

# df_subset <- df[, c("SN_Types", "Variables", "Overall_all_lifestyles", "Female_all_lifestyles", "Male_all_lifestyles")]
# df_subset <- df[, c("SN_Types", "Variables", "Overall_all_lifestyles")]

df_subset <- mydf[grepl("Any SN", mydf$SN_types),]
df_subset$Cohort <- paste(df_subset$Cohort, df_subset$Variables, sep = "_")

# Melt the data to long format
library(reshape2)
df_melted <- melt(df_subset, id.vars = c("SN_types", "Cohort", "Variables"))
df_melted$value <- as.numeric(df_melted$value)

library(RColorBrewer)

# Define the number of colors you need (corresponding to the number of variables)
num_colors <- length(unique(df_melted$variable))
color_palette <- brewer.pal(num_colors, "Set1")
# Create a grouped bar plot using ggplot2
library(ggplot2)
ggplot(df_melted, aes(x = Cohort, y = value, fill = variable)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(x = "Conditions", y = "Attributable Fraction") +
  # facet_wrap(~SN_types, scales = "free_y", ncol = 1) +
  theme_minimal() +
  theme(legend.position = "top", axis.text.x = element_text(angle = 45, hjust = 0.2)) +
  scale_y_continuous(limits = c(-0.2, 1), breaks = seq(-0.2, 1, 0.4)) +
scale_fill_manual(values = color_palette)


###########
ggplot(df_melted, aes(x = Variables, y = value, fill = variable)) +
  geom_bar(stat = "identity", position = "group") +
  geom_text(aes(label = ifelse(value >= 0, round(value, 2), paste0("(", round(value, 2), ")"))),
            position = position_stack(vjust = 0.5), color = "black", size = 3) +
  labs(x = "Conditions", y = "Attributable Fraction") +
  facet_wrap(~ SN_types, scales = "free_y", ncol = 1) +
  theme_minimal() +
  theme(legend.position = "top", axis.text.x = element_text(angle = 45, hjust = 0.2)) +
  scale_y_continuous(limits = c(-1, 1), breaks = seq(-1, 1, 0.2)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#########################
























mydf <- read.table(text="Variables	Overall	Female	Male	Age<35	Age≥35	Overall	Female	Male	Age<35	Age≥35	Overall	Female	Male	Age<35	Age≥35	Overall	Female	Male	Age<35	Age≥35	Overall	Female	Male	Age<35	Age≥35	Overall	Female	Male	Age<35	Age≥35
Radiation	0.435	0.424	0.449	0.345	0.455	0.438	0.428	0.45	0.347	0.458	0.438	0.429	0.449	0.346	0.458	0.439	0.43	0.45	0.347	0.459	0.439	0.43	0.451	0.347	0.46	0.439	0.428	0.452	0.347	0.459
Chemo	0.034	0.034	0.034	0.06	0.028	0.035	0.035	0.035	0.062	0.029	0.036	0.036	0.036	0.064	0.029	0.036	0.036	0.036	0.065	0.03	0.036	0.036	0.036	0.065	0.03	0.036	0.037	0.036	0.064	0.03
All treatments	0.453	0.443	0.465	0.381	0.469	0.456	0.448	0.466	0.384	0.472	0.456	0.449	0.466	0.384	0.472	0.457	0.45	0.467	0.386	0.473	0.458	0.45	0.468	0.386	0.474	0.458	0.449	0.469	0.386	0.474
PRS	0.184	0.182	0.187	0.187	0.184	0.184	0.181	0.187	0.186	0.183	0.18	0.177	0.182	0.181	0.179	0.18	0.177	0.183	0.182	0.179	0.179	0.177	0.182	0.181	0.179	0.18	0.178	0.183	0.182	0.18
Lifestyle	0.213	0.207	0.22	0.218	0.212	-0.041	-0.043	-0.039	-0.03	-0.044	-0.024	-0.02	-0.029	-0.023	-0.024	0.027	0.03	0.022	0.024	0.027	0.001	0.001	0.001	0.001	0.002	0.267	0.262	0.275	0.271	0.267
Combined	0.648	0.64	0.658	0.602	0.658	0.538	0.53	0.548	0.48	0.551	0.543	0.538	0.55	0.483	0.557	0.568	0.563	0.574	0.508	0.581	0.557	0.55	0.566	0.497	0.57	0.675	0.668	0.685	0.633	0.685
Radiation	0.376	0.374	0.379	0.266	0.402	0.378	0.378	0.377	0.266	0.404	0.378	0.38	0.377	0.266	0.405	0.378	0.38	0.376	0.266	0.404	0.382	0.382	0.382	0.269	0.408	0.378	0.377	0.38	0.267	0.404
Chemo	0.027	0.026	0.027	0.043	0.023	0.03	0.029	0.031	0.049	0.026	0.031	0.03	0.032	0.05	0.026	0.03	0.029	0.031	0.049	0.025	0.03	0.029	0.031	0.049	0.025	0.031	0.031	0.032	0.051	0.027
All treatments	0.394	0.392	0.395	0.298	0.416	0.398	0.399	0.396	0.302	0.42	0.399	0.401	0.396	0.303	0.421	0.398	0.4	0.395	0.302	0.42	0.401	0.402	0.4	0.305	0.424	0.399	0.398	0.4	0.304	0.421
PRS	0.15	0.147	0.154	0.148	0.151	0.151	0.148	0.154	0.149	0.151	0.147	0.145	0.15	0.145	0.148	0.148	0.146	0.151	0.146	0.149	0.147	0.145	0.15	0.145	0.148	0.147	0.145	0.15	0.145	0.148
Lifestyle	0.284	0.277	0.292	0.287	0.283	-0.032	-0.032	-0.031	-0.025	-0.033	-0.011	-0.009	-0.013	-0.01	-0.011	0.063	0.073	0.052	0.055	0.065	-0.046	-0.047	-0.045	-0.032	-0.049	0.354	0.347	0.364	0.36	0.353
Combined	0.63	0.628	0.633	0.569	0.645	0.472	0.473	0.471	0.389	0.492	0.482	0.484	0.479	0.397	0.502	0.52	0.526	0.512	0.435	0.54	0.466	0.467	0.465	0.385	0.485	0.67	0.668	0.674	0.618	0.682
Radiation	0.276	0.273	0.28	0.24	0.28	0.275	0.28	0.268	0.243	0.278	0.263	0.26	0.267	0.235	0.266	0.277	0.28	0.273	0.243	0.28	0.278	0.287	0.267	0.248	0.282	0.273	0.279	0.266	0.243	0.277
Chemo	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-
All treatments	0.276	0.273	0.28	0.24	0.28	0.275	0.28	0.268	0.243	0.278	0.263	0.26	0.267	0.235	0.266	0.277	0.28	0.273	0.243	0.28	0.278	0.287	0.267	0.248	0.282	0.273	0.279	0.266	0.243	0.277
PRS	0.437	0.423	0.456	0.442	0.437	0.416	0.401	0.435	0.419	0.416	0.44	0.426	0.458	0.446	0.439	0.418	0.404	0.437	0.42	0.418	0.413	0.397	0.432	0.414	0.412	0.417	0.403	0.436	0.419	0.417
Lifestyle	-0.49	-0.437	-0.557	-0.398	-0.5	0.015	0.016	0.013	0.013	0.015	-0.151	-0.111	-0.201	-0.103	-0.156	-0.064	-0.067	-0.059	-0.051	-0.065	-0.076	-0.077	-0.074	-0.055	-0.078	-0.071	-0.069	-0.074	-0.071	-0.072
Combined	0.398	0.405	0.388	0.405	0.397	0.586	0.58	0.594	0.566	0.588	0.528	0.533	0.522	0.527	0.528	0.556	0.547	0.568	0.54	0.558	0.549	0.542	0.557	0.535	0.55	0.55	0.544	0.558	0.529	0.552
Radiation	0.486	0.507	0.468	0.432	0.489	0.491	0.515	0.468	0.432	0.495	0.493	0.517	0.47	0.423	0.498	0.494	0.516	0.473	0.427	0.498	0.49	0.515	0.466	0.421	0.494	0.494	0.516	0.473	0.426	0.498
Chemo	0.074	0.076	0.072	0.058	0.075	0.054	0.057	0.051	0.04	0.055	0.057	0.058	0.055	0.048	0.057	0.058	0.059	0.056	0.05	0.058	0.057	0.058	0.056	0.05	0.058	0.048	0.049	0.047	0.038	0.049
All treatments	0.526	0.546	0.51	0.471	0.53	0.521	0.544	0.499	0.461	0.524	0.524	0.547	0.502	0.455	0.528	0.525	0.547	0.505	0.459	0.529	0.521	0.545	0.498	0.453	0.525	0.52	0.542	0.499	0.452	0.524
PRS	0.198	0.193	0.203	0.216	0.197	0.197	0.19	0.203	0.214	0.196	0.164	0.158	0.169	0.171	0.163	0.163	0.156	0.17	0.174	0.163	0.16	0.154	0.165	0.169	0.159	0.158	0.154	0.163	0.166	0.158
Lifestyle	0.112	0.058	0.157	0.194	0.107	-0.242	-0.262	-0.223	-0.192	-0.245	-0.055	-0.049	-0.061	-0.066	-0.054	-0.114	-0.139	-0.092	-0.083	-0.116	-0.071	-0.066	-0.076	-0.057	-0.072	0.374	0.366	0.381	0.377	0.374
Combined	0.661	0.657	0.664	0.652	0.662	0.524	0.537	0.511	0.48	0.526	0.581	0.602	0.562	0.519	0.585	0.56	0.569	0.551	0.512	0.563	0.567	0.588	0.548	0.516	0.571	0.748	0.756	0.741	0.714	0.75
Radiation	0.608	0.632	0.58	0.545	0.625	0.611	0.633	0.586	0.555	0.627	0.608	0.629	0.584	0.548	0.625	0.61	0.631	0.586	0.552	0.627	0.609	0.63	0.584	0.55	0.626	0.608	0.629	0.584	0.546	0.626
Chemo	0.254	0.256	0.252	0.377	0.22	0.25	0.251	0.249	0.368	0.216	0.244	0.245	0.243	0.36	0.211	0.254	0.256	0.251	0.378	0.218	0.248	0.249	0.248	0.366	0.214	0.245	0.248	0.243	0.361	0.213
All treatments	0.731	0.747	0.712	0.711	0.736	0.73	0.746	0.713	0.714	0.735	0.725	0.74	0.708	0.703	0.731	0.732	0.747	0.716	0.716	0.737	0.728	0.743	0.711	0.709	0.733	0.726	0.741	0.71	0.704	0.732
PRS	0.581	0.579	0.582	0.591	0.578	0.575	0.575	0.576	0.587	0.572	0.568	0.567	0.57	0.576	0.566	0.572	0.57	0.575	0.582	0.57	0.571	0.569	0.573	0.58	0.568	0.565	0.565	0.565	0.572	0.563
Lifestyle	0.295	0.282	0.31	0.298	0.294	0.124	0.128	0.119	0.119	0.125	-0.063	-0.057	-0.071	-0.067	-0.062	-0.158	-0.191	-0.121	-0.123	-0.168	-0.054	-0.053	-0.055	-0.039	-0.058	0.419	0.416	0.423	0.418	0.419
Combined	0.921	0.924	0.919	0.915	0.923	0.901	0.906	0.895	0.892	0.904	0.873	0.88	0.866	0.861	0.877	0.868	0.871	0.865	0.861	0.87	0.877	0.883	0.87	0.868	0.88	0.93	0.933	0.927	0.924	0.932
Radiation	0.186	0.149	0.23	0.281	0.172	0.189	0.154	0.23	0.28	0.175	0.187	0.152	0.228	0.278	0.174	0.19	0.155	0.233	0.275	0.178	0.188	0.151	0.233	0.274	0.176	0.188	0.154	0.229	0.272	0.176
Chemo	0.34	0.344	0.334	0.501	0.316	0.341	0.345	0.336	0.506	0.317	0.339	0.342	0.335	0.497	0.316	0.341	0.343	0.338	0.506	0.317	0.337	0.342	0.331	0.496	0.314	0.338	0.34	0.336	0.503	0.314
All treatments	0.471	0.449	0.498	0.673	0.442	0.474	0.454	0.498	0.68	0.444	0.469	0.449	0.494	0.671	0.44	0.474	0.451	0.501	0.678	0.444	0.469	0.448	0.495	0.67	0.44	0.468	0.446	0.494	0.671	0.438
PRS	-0.118	-0.113	-0.125	-0.118	-0.119	-0.114	-0.108	-0.121	-0.112	-0.114	-0.118	-0.113	-0.125	-0.116	-0.119	-0.111	-0.104	-0.119	-0.108	-0.111	-0.12	-0.113	-0.128	-0.117	-0.12	-0.116	-0.108	-0.124	-0.114	-0.116
Lifestyle	0.257	0.254	0.261	0.272	0.255	-0.088	-0.08	-0.098	-0.068	-0.091	-0.168	-0.134	-0.208	-0.139	-0.172	-0.106	-0.123	-0.087	-0.089	-0.109	0.139	0.133	0.145	0.101	0.144	0.333	0.328	0.34	0.335	0.333
Combined	0.559	0.543	0.578	0.729	0.534	0.362	0.345	0.383	0.612	0.325	0.306	0.303	0.309	0.576	0.266	0.357	0.324	0.396	0.612	0.32	0.488	0.467	0.512	0.666	0.462	0.605	0.588	0.625	0.755	0.583
", header = T, check.names = F)


# Get the current column names
current_names <- colnames(mydf)[-1]

# Append the variable names in a cyclical manner to the first five columns
variable_names <- c("all_lifestyles", "smoking", "drinking", "physical_activity", "obesity", "diet")
new_names <- rep(variable_names, each = 5)

# Set the new unique column names to your data frame
colnames(mydf)[-1] <- paste(current_names, new_names, sep = "_")

sn_vars <- c("Any SN (303)", "SMN (234)", "NMSC (118)", "Breast cancer (53)", "Thyroid cancer (43)", "Meningioma (81)")
mydf$SN_Types <- rep(sn_vars, each = 6)

mydf[mydf == "-"] <- NA

mydf[!grepl("Variables|SN_Types", colnames(mydf))] <- sapply(mydf[!grepl("Variables|SN_Types", colnames(mydf))], as.numeric)

df <-  mydf


df_subset <- df[df$SN_Types == "Any SN (303)",]

df_subset <- df[, c("SN_Types", "Variables", "Overall_all_lifestyles", "Female_all_lifestyles", "Male_all_lifestyles")]
df_subset <- df[, c("SN_Types", "Variables", "Overall_all_lifestyles")]



# Melt the data to long format
library(reshape2)
df_melted <- melt(df_subset, id.vars = c("SN_Types", "Variables"))

# Create a grouped bar plot using ggplot2
library(ggplot2)
ggplot(df_melted, aes(x = Variables, y = value, fill = variable)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(x = "Conditions", y = "Attributable Fraction") +
  facet_wrap(~SN_Types, scales = "free_y", ncol = 1) +
  theme_minimal() +
  theme(legend.position = "top", axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_y_continuous(limits = c(-1, 1), breaks = seq(-1, 1, 0.4))











# Plot the bar plot
ggplot(mydf, aes(x = Variables, y = Overall_all_lifestyles, fill = SN_Types)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(x = "Variables", y = "Attributable Fraction", title = "Attributable Fraction for Different Variables Across Subsequent Neoplasms") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_discrete(name = "Subsequent Neoplasms") +
theme(axis.text.x = element_text(angle = 45, hjust = 1))


ggplot(mydf, aes(x = Variables, y = Overall_all_lifestyles, fill = SN_Types)) +
  geom_bar(stat = "identity") +
  facet_wrap(~ SN_Types, scales = "free") +
  # scale_fill_discrete(name = "Subsequent Neoplasms") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# Reshape the dataframe to long format
mydf_long <- pivot_longer(mydf, cols = -c(Variables, SN_Types), names_to = "Lifestyle", values_to = "Attributable_Fraction")

# Plot the bar plot
ggplot(mydf_long, aes(x = Variables, y = Attributable_Fraction, fill = SN_Types)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~Lifestyle, scales = "free", ncol = 3) +
  labs(x = "Variables", y = "Attributable Fraction", title = "Attributable Fraction for Different Variables Across Subsequent Neoplasms") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_discrete(name = "Subsequent Neoplasms")






