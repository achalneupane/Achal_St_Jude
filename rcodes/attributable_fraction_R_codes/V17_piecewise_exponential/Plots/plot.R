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

df <-  mydf

mydf[mydf == "-"] <- NA

mydf[!grepl("Variables|SN_Types", colnames(mydf))] <- sapply(mydf[!grepl("Variables|SN_Types", colnames(mydf))], as.numeric)

# Plot the bar plot
ggplot(mydf, aes(x = Variables, y = Overall_all_lifestyles, fill = SN_Types)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(x = "Variables", y = "Attributable Fraction", title = "Attributable Fraction for Different Variables Across Subsequent Neoplasms") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_discrete(name = "Subsequent Neoplasms") +
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




ggplot(mydf_long, aes(x = Variables, y = Attributable_Fraction, group = SN_Types, color = Variables)) +
  geom_line() +
  geom_point() +
  labs(x = "Subsequent Neoplasms", y = "Attributable Fraction", title = "Trends of Attributable Fraction for Different Variables Across Subsequent Neoplasms") +
  scale_x_continuous(breaks = 1:length(unique(mydf$SN_Types)), labels = unique(mydf$SN_Types)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_color_discrete(name = "Variables")
# library(ggplot2)
# 
# # Create a new dataframe with the data in long format
# mydf_long <- tidyr::pivot_longer(mydf, cols = c("Overall", "Female", "Male", "Age_less_than_35", "Age_greater_than_or_equal_to_35"), 
#                                  names_to = "Gender_and_Age_Group", values_to = "Value")
# 
# # Create a stacked bar chart
# ggplot(mydf_long, aes(x = SN_types_N_cases, y = Value, fill = Gender_and_Age_Group)) +
#   geom_bar(stat = "identity") +
#   labs(title = "Comparison of Different Variables by SN Types",
#        x = "SN Types (N Cases)", y = "Proportion",
#        fill = "Gender and Age Group") +
#   theme_minimal() +
#   theme(legend.position = "bottom")


mydf_long <- gather(mydf, key = "AgeGroup", value = "AttributableFraction", Overall:Age_greater_than_or_equal_to_35)

# Plot
ggplot(mydf_long, aes(x = SN_types_N_cases, y = AttributableFraction, fill = Variables)) +
  geom_bar(position = "dodge", stat = "identity", width = 0.7) +
  labs(title = "Attributable Fraction for Different Subsequent Neoplasms",
       x = "Subsequent Neoplasms (N cases)",
       y = "Attributable Fraction",
       fill = "Variables") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_brewer(palette = "Set1")




