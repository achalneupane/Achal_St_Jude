import numpy as np
import pandas as pd
import re

all_cancers = pd.read_csv("Z:/ResearchHome/Groups/sapkogrp/projects//Genomics/common/attr_fraction/prs/ALL_Cancers_PRS_data.txt", sep='\t', header = 0) # specify the first row as the header
# Extract only those that are needed
all_cancers.columns.values
all_cancers.TYPE.value_counts()
all_cancers.Cancer.value_counts()

# View(pd.crosstab(index=all_cancers['TYPE'], columns=all_cancers['Cancer']))
all_cancers = all_cancers[all_cancers['TYPE'].str.contains("Mavaddat_2019|Pleiotropy_PRSWEB|Basal_cell_carcinoma|Squamous_cell_carcinoma|Sarcoma|Meningioma|ALL_Vijayakrishnan", regex=True)] # grepl equivalent
all_cancers.TYPE.value_counts()
all_cancers.Cancer.value_counts()

thyroid_df = pd.read_csv("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/prs/thyroid/thyroid.score", sep = ' ', header = None)


# thyroid_df[["CHROM", "POS_GRCh38", "REF", "Effect_allele"]] = thyroid_df[0].str.split(':', 3, expand=True)
new_cols = ["CHROM", "POS_GRCh38", "REF", "Effect_allele"]
ex_thyroid_df= thyroid_df[0].str.split(':',expand=True).rename(columns={i:new_cols[i] for i in range(4)})

new_cols = ["true_effect_allele", "Effect_size"]
ex_thyroid_df_tmp = thyroid_df.iloc[:, 1:3].reset_index(drop=True)

# ex_thyroid_df_tmp = ex_thyroid_df_tmp.set_axis(new_cols, axis=1, inplace=False)
ex_thyroid_df_tmp.columns = new_cols

df_thyroid = pd.concat([ex_thyroid_df, ex_thyroid_df_tmp], axis=1)
#     CHROM POS_GRCh38 REF Effect_allele true_effect_allele  Effect_size
# 0    chr1  233276815   A             G                  A     0.277632
# 1    chr2  217427435   C             G                  C     0.357674
# 2    chr3  169800667   T             G                  T     0.207014
# 3    chr5    1279675   C             T                  T     0.182322
# 4    chr5  112150207   A             T                  A     0.314811
# 5    chr8   32575278   G             T                  G     0.277632
# 6    chr9   97775520   A             C                  A     0.524729
# 7   chr10  103934543   C             T                  T     0.343590
# 8   chr14   36063370   G             C                  G     0.329304
# 9   chr14   36269155   C             T                  T     0.593327
# 10  chr15   67165147   G             C                  C     0.207014
# 11  chr15   67163292   C             T                  T     0.215111
# df_thyroid.columns.to_list()loc['Effect_allele']


# Re-order Reference and Effect alleles extracted from CHR:POS:REF:ALT
m = df_thyroid['Effect_allele'] != df_thyroid['true_effect_allele']
df_thyroid.loc[m, ['REF', 'Effect_allele']] = (df_thyroid.loc[m, ['Effect_allele', 'REF']].values)
df_thyroid = df_thyroid.loc[:, ['CHROM', 'POS_GRCh38', 'REF', 'Effect_allele', 'Effect_size']]

# Add columns: TYPE, Cancer, Significant_YN
df_thyroid['TYPE'] = "THYROID_PGS"
df_thyroid['Cancer'] = "THYROID"
df_thyroid['Significant_YN'] = "Y"

# remove chr from df_thyroid
df_thyroid['CHROM'] = df_thyroid['CHROM'].str.replace(r'\D+', '', regex=True).astype('int')
# df_thyroid['CHROM'].replace("chr", "", regex=True)

# Row bind two dataframes
df_thyroid.shape
all_cancers.shape
all_cancers = pd.concat([all_cancers, df_thyroid])


