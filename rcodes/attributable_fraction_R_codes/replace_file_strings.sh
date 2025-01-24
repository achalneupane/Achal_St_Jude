find ./ -type f| grep V11_with_overlapping| egrep -v "Pheno|Rhistory|Merge_all_admixture" | xargs -I {} sed -i 's/\bsmoker_former_or_never_yn\b/Current_smoker_yn/g' {}
find ./ -type f| grep V11_with_overlapping| egrep -v "Pheno|Rhistory|Merge_all_admixture" | xargs -I {} sed -i 's/\bNOT_RiskyHeavyDrink_yn\b/RiskyHeavyDrink_yn/g' {}
find ./ -type f| grep V11_with_overlapping| egrep -v "Pheno|Rhistory|Merge_all_admixture" | xargs -I {} sed -i 's/\bNot_obese_yn\b/Obese_yn/g' {}


find ./ -type f| grep V11_with_overlapping| egrep -v "Pheno|Rhistory|Merge_all_admixture" | xargs -I {} sed -i 's/\bMavaddat_2019_ER_POS_Breast_PRS.tertile.category +//g' {}
find ./ -type f| grep V11_with_overlapping| egrep -v "Pheno|Rhistory|Merge_all_admixture" | xargs -I {} sed -i 's/\bMavaddat_2019_ER_NEG_Breast_PRS.tertile.category +//g' {}
find ./ -type f| grep V11_with_overlapping| egrep -v "Pheno|Rhistory|Merge_all_admixture" | xargs -I {} sed -i 's/\bSQUAMOUScell_PRS.tertile.category +//g' {}



find ./ -type f| grep V11_with_overlapping| egrep -v "Pheno|Rhistory|Merge_all_admixture" | xargs -I {} sed -i 's/\bsmoker_former_or_never_yn\b/Current_smoker_yn/g' {}




find ./ -type f| grep V11_with_overlapping| egrep -v "Pheno|Rhistory|Merge_all_admixture" | xargs -I {} sed -i '/^## Nullify Lifestyle/,$s/dat_lifestyle/dat_tx.plp.prs.lifestyle/g' {}



## After Lancet oncology revision
## comment out rounding for bootstrapping:
find ./ -type f -name "*model_fit*" -exec sed -i '/round(/ s/^/# /' {} +
find ./ -type f -name "*model_fit*" -exec sed -i '/View(/ s/^/# /' {} +



job=test4_$(date +%Y%m%d_%H%M%S)
bsub -q priority \
     -P bootstr-J \
     -J $job \
     -eo logs/$job.err \
     -oo logs/$job.out \
     -R "rusage[mem=20000]" \
     bash -c "module load R/4.2.2 && Rscript test4.R"


job=model_fit1_$(date +%Y%m%d_%H%M%S)
bsub -q priority \
     -P bootstr-J \
     -J $job \
     -eo logs/$job.err \
     -oo logs/$job.out \
     -R "rusage[mem=20000]" \
     bash -c "module load R/4.2.2 && Rscript 1.model_fit_Any_SN.R"


job=model_fit2_$(date +%Y%m%d_%H%M%S)
bsub -q priority \
     -P bootstr-J \
     -J $job \
     -eo logs/$job.err \
     -oo logs/$job.out \
     -R "rusage[mem=20000]" \
     bash -c "module load R/4.2.2 && Rscript 2.model_fit_SMNs.R"



job=model_fit3_$(date +%Y%m%d_%H%M%S)
bsub -q priority \
     -P bootstr-J \
     -J $job \
     -eo logs/$job.err \
     -oo logs/$job.out \
     -R "rusage[mem=20000]" \
     bash -c "module load R/4.2.2 && Rscript 3.model_fit_NMSCs.R"



job=model_fit4_$(date +%Y%m%d_%H%M%S)
bsub -q priority \
     -P bootstr-J \
     -J $job \
     -eo logs/$job.err \
     -oo logs/$job.out \
     -R "rusage[mem=20000]" \
     bash -c "module load R/4.2.2 && Rscript 4.model_fit_BREASTcancer.R"


job=model_fit5_$(date +%Y%m%d_%H%M%S)
bsub -q priority \
     -P bootstr-J \
     -J $job \
     -eo logs/$job.err \
     -oo logs/$job.out \
     -R "rusage[mem=20000]" \
     bash -c "module load R/4.2.2 && Rscript 5.model_fit_THYROIDcancer.R"


job=model_fit6_$(date +%Y%m%d_%H%M%S)
bsub -q priority \
     -P bootstr-J \
     -J $job \
     -eo logs/$job.err \
     -oo logs/$job.out \
     -R "rusage[mem=20000]" \
     bash -c "module load R/4.2.2 && Rscript 6.model_fit_MENINGIOMA.R"


job=model_fit7_$(date +%Y%m%d_%H%M%S)
bsub -q priority \
     -P bootstr-J \
     -J $job \
     -eo logs/$job.err \
     -oo logs/$job.out \
     -R "rusage[mem=20000]" \
     bash -c "module load R/4.2.2 && Rscript 7.model_fit_SARCOMA.R"