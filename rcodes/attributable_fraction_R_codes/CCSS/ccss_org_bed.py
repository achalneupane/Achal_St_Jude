import numpy as np
import pandas as pd
import re
import os
import natsort
from natsort import natsort_keygen


all_cancers = pd.read_csv("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/ccss_org_hrc/ccss_org_hrc_vcf_GRCh38/all_cancer_GrCh37.bed", sep='\t', header=None) # specify no header
# Extract only those that are needed

all_cancers.head
all_cancers.iloc[:,3]
all_cancers['KEY_GRCh37'] = all_cancers.iloc[:,0].astype(str) + ":" + all_cancers.iloc[:,2].astype(str)


all_cancers['GRCh38_POS'] = all_cancers.iloc[:,3].apply(lambda x: re.match(r'.*:(.*)-', x).group(1))
# all_cancers.iloc[:,3].apply(lambda s: s[s.index(":") + 1: s.rindex("-")])
# all_cancers.iloc[:,3].apply(lambda x: re.match(r'.*:(.*)-', x).group(1))
# all_cancers.iloc[:,3].str.extract(r'.*:(.*)-')
all_cancers['KEY_GRCh38'] = all_cancers.iloc[:,0].astype(str) + ":" + all_cancers.iloc[:,5].astype(str)

bim_file = pd.read_csv("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/ccss_org_hrc/ccss_org_hrc_vcf_GRCh38/attr_fraction/prs/plink_data/merged.dat.bim", sep = '\t', header = None)
bim_file.head

bim_file.shape
bim_file['KEY_GRCh37'] = "chr"+bim_file.iloc[:,0].astype(str) + ":" + bim_file.iloc[:,3].astype(str)

# c = [5,4,3,2,1,'c']
# d = [2,3,'c']
# # [ b.index(x) if x in b else None for x in a ] # match(a,b)
# match = lambda a, b: [ b.index(x) if x in b else None for x in a ]
# d[match(c,d).astype(int)]

# lambda a, b: [ b.index(x) if x in b else None for x in a ]
# [d[ind] for ind in [d.index(x) if x in d else None  for x in c if x in d] if ind is not None]
# [x if x in d else None  for x in c]

# [all_cancers['KEY_GRCh38'] if x in all_cancers['KEY_GRCh37'] else None  for x in bim_file['KEY_GRCh37']]

# bim_file['KEY_GRCh38'] = all_cancers['KEY_GRCh38'][bim_file['KEY_GRCh37'].isin(all_cancers['KEY_GRCh37'])]



# bim_file.merge(all_cancers, left_on='KEY_GRCh37', right_on='KEY_GRCh37')[['0_x', '1_x', '2_x', '3_x', 'KEY_GRCh37', 'KEY_GRCh38']]


# cc = pd.merge(bim_file, all_cancers, on=['KEY_GRCh37'], how = 'left')
# cc['KEY_GRCh38'].isna().sum()


# df = pd.DataFrame(dict(
#         AUTHOR_NAME=list('AAABBCCCCDEEFGG'),
#         title=      list('zyxwvutsrqponml')
#     ))
# 
# df2 = pd.DataFrame(dict(
#         AUTHOR_NAME=list('AABCCEGG'),
#         title      =list('zwvtrpml'),
#         CATEGORY   =list('11223344')
#     ))
# 
# df_concat = pd.concat([df2, df]).reset_index().drop_duplicates(['Index', 'AUTHOR_NAME'])
# df_concat.set_index('Index', inplace=True)
# df_concat[df_concat.index.isin(df.index)]
