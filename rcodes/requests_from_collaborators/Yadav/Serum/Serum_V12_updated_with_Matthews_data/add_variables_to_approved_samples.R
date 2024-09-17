df <- read.table("Z:/ResearchHome/ClusterHome/aneupane/data/Yadav_serum/v12_output/Approved Samples(N=1200).txt", header = T, sep = "\t")
df.core <- read.table("Z:/ResearchHome/ClusterHome/aneupane/data/Yadav_serum/v12_output/plasma.df.1200.with_additional_variables_toYadav.txt", header = T, sep = "\t")

df$racegrp <-  df.core$racegrp[match(df$tb_number, df.core$tb_number)]
df$Sex <-  df.core$Sex[match(df$tb_number, df.core$tb_number)]
df$selection_group <-  df.core$selection_group[match(df$tb_number, df.core$tb_number)]
df$Sample_age <-  df.core$Sample_age[match(df$tb_number, df.core$tb_number)]

write.table(df, "Z:/ResearchHome/ClusterHome/aneupane/data/Yadav_serum/v12_output/Approved_Samples(N=1200)_to_core.txt", col.names = T, row.names = F, sep = "\t", quote = F)
