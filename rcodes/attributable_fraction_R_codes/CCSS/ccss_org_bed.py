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


bim_file['KEY_GRCh38'] = all_cancers.KEY_GRCh38[bim_file.KEY_GRCh37.isin(all_cancers.KEY_GRCh37.unique())]

c = [5,4,3,2,1,'c']
d = [2,3,'c']
# [ b.index(x) if x in b else None for x in a ] # match(a,b)
match = lambda a, b: [ b.index(x) if x in b else None for x in a ]
d[match(c,d).astype(int)]


