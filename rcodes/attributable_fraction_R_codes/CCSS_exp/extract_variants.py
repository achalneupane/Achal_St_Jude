import numpy as np
import pandas as pd


all_cancers = pd.read_csv("Z:/ResearchHome/Groups/sapkogrp/projects//Genomics/common/attr_fraction/prs/ALL_Cancers_PRS_data.txt", sep='\t', header = 0) # specify the first row as the header
thyroid_df = pd.read_csv("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/prs/thyroid/thyroid.score", sep = ' ', header = None)


# thyroid_df[["CHROM", "POS_GRCh38", "REF", "Effect_allele"]] = thyroid_df[0].str.split(':', 3, expand=True)
new_cols = ["CHROM", "POS_GRCh38", "REF", "Effect_allele"]
ex_thyroid_df= thyroid_df[0].str.split(':',expand=True).rename(columns={i:new_cols[i] for i in range(4)})

new_cols = ["true_effect_allele", "Effect_size"]
ex_thyroid_df_tmp = thyroid_df.iloc[:, 1:3].reset_index(drop=True)

# ex_thyroid_df_tmp = ex_thyroid_df_tmp.set_axis(new_cols, axis=1, inplace=False)
ex_thyroid_df_tmp.columns = new_cols

df_thyroid = pd.concat([ex_thyroid_df, ex_thyroid_df_tmp], axis=1)

# df_thyroid.columns.to_list()loc['Effect_allele']

df_thyroid.loc[df_thyroid['Effect_allele'] != df_thyroid['true_effect_allele'], "Effect_allele" ] = df_thyroid.loc[df_thyroid['Effect_allele'] != df_thyroid['true_effect_allele'], "true_effect_allele" ]


df[df$Effect_allele != df$V2, c("REF", "Effect_allele")] <- df[df$Effect_allele != df$V2, c("Effect_allele", "REF")]


pd.concat([all_cancers, df_thyroid], axis = 1)


