ln -s /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/attr_fraction/prs/thyroid/thyroid.score .
ln -s /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/attr_fraction/prs/thyroid/PGS001354.txt_edited_hg38.bed_indat .

plink --bfile /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/attr_fraction/prs/plink_data/sjlife_all_PRS_all_final_v3 --extract thyroid.score --make-bed --out sjlife_all_PRS_all_final_v3_thyroid
plink --bfile sjlife_all_PRS_all_final_v3_thyroid --out sjlife_thyroid --score thyroid.score
