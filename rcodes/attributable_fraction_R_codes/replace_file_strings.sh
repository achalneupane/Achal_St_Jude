find ./ -type f| grep V11_with_overlapping| egrep -v "Pheno|Rhistory|Merge_all_admixture" | xargs -I {} sed -i 's/\bsmoker_former_or_never_yn\b/Current_smoker_yn/g' {}
find ./ -type f| grep V11_with_overlapping| egrep -v "Pheno|Rhistory|Merge_all_admixture" | xargs -I {} sed -i 's/\bNOT_RiskyHeavyDrink_yn\b/RiskyHeavyDrink_yn/g' {}
find ./ -type f| grep V11_with_overlapping| egrep -v "Pheno|Rhistory|Merge_all_admixture" | xargs -I {} sed -i 's/\bNot_obese_yn\b/Obese_yn/g' {}


find ./ -type f| grep V11_with_overlapping| egrep -v "Pheno|Rhistory|Merge_all_admixture" | xargs -I {} sed -i 's/\bMavaddat_2019_ER_POS_Breast_PRS.tertile.category +//g' {}
find ./ -type f| grep V11_with_overlapping| egrep -v "Pheno|Rhistory|Merge_all_admixture" | xargs -I {} sed -i 's/\bMavaddat_2019_ER_NEG_Breast_PRS.tertile.category +//g' {}
find ./ -type f| grep V11_with_overlapping| egrep -v "Pheno|Rhistory|Merge_all_admixture" | xargs -I {} sed -i 's/\bSQUAMOUScell_PRS.tertile.category +//g' {}