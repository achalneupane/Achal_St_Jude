import numpy as np
import pandas as pd
from sklearn import svm

all_cancers = pd.read_csv("Z:/ResearchHome/Groups/sapkogrp/projects//Genomics/common/attr_fraction/prs/ALL_Cancers_PRS_data.txt", sep='\t', header = 0) # specify the first row as the header
thyroid_df = pd.read_csv("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/prs/thyroid/thyroid.score", sep = ' ', header = None)

thyroid_df

ex_thyroid_df = thyroid_df[0].str.split(':', 3, expand=True).rename(columns = ["CHROM", "POS_GRCh38", "REF", "Effect_allele"])
ex_thyroid_df_tmp = thyroid_df.iloc[:, 1:3].reset_index(drop=True).rename("true_effect_allele", "Effect_size")
df_c = pd.concat([ex_thyroid_df, ex_thyroid_df_tmp, axis=1)
