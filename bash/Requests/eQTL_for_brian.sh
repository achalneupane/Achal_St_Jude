##-------1. parallel submission

for i in {1..22}; do \
export CHR="chr${i}"; \
echo "Annotating $CHR"; \
bsub \
-P "${CHR}_eQTL" \
-J "${CHR}_eQTL" \
-q "rhel8_gpu" \
-o "/research_jude/rgs01_jude/groups/sapkogrp/projects/QTL_projects/common/eQTL_sjlife/logs/${CHR}_eQTL.%J" \
-gpu "num=1/host" \
-n "1" \
-R "rusage[mem=100GB]" \
"./entrypoint_eQTL_achal_with_status.sh"; \
done;

##------2. entrypoint_eQTL_achal.sh: 
cd /research_jude/rgs01_jude/groups/sapkogrp/projects/QTL_projects/common/eQTL_sjlife/
module load R/4.3.1
export CHR="${CHR}"

head -1 all_chr.Survivor_WGS.GATK4180.hg38.PASS.decomposed.vcf.gz_pruned_all_white3.vcf > ${CHR}_Survivor_WGS.GATK4180.hg38.PASS.decomposed.vcf.gz_pruned_all_white3.vcf
grep "^${CHR}:" all_chr.Survivor_WGS.GATK4180.hg38.PASS.decomposed.vcf.gz_pruned_all_white3.vcf >> ${CHR}_Survivor_WGS.GATK4180.hg38.PASS.decomposed.vcf.gz_pruned_all_white3.vcf

Rscript ./run_matrix_eqtl_achal.r



#-------3. run_matrix_eqtl_achal.r: 
CHR <- Sys.getenv("CHR")
library("MatrixEQTL")
useModel = modelLINEAR
SNP_file_name = paste0(CHR, "_Survivor_WGS.GATK4180.hg38.PASS.decomposed.vcf.gz_pruned_all_white3.vcf")
expression_file_name = ("gene_expression_all_white.txt")
covariates_file_name = ("all_covariates_all_white_no_peer3.txt")
output_file_name = tempfile()
pvOutputThreshold = 1e-6
errorCovariance = numeric()

snps = SlicedData$new()
snps$fileDelimiter = "\t"      # the TAB character
snps$fileOmitCharacters = "./." # denote missing values;
snps$fileSkipRows = 1          # one row of column labels
snps$fileSkipColumns = 1       # one column of row labels
snps$fileSliceSize = 2000     # read file in slices of 2,000 rows
snps$LoadFile(SNP_file_name)

gene = SlicedData$new()
gene$fileDelimiter = "\t"      # the TAB character
gene$fileOmitCharacters = "NA" # denote missing values;
gene$fileSkipRows = 1          # one row of column labels
gene$fileSkipColumns = 1       # one column of row labels
gene$fileSliceSize = 2000      # read file in slices of 2,000 rows
gene$LoadFile(expression_file_name)

cvrt = SlicedData$new()
cvrt$fileDelimiter = "\t"      # the TAB character
cvrt$fileOmitCharacters = "NA" # denote missing values;
cvrt$fileSkipRows = 1          # one row of column labels
cvrt$fileSkipColumns = 1       # one column of row labels
if(length(covariates_file_name)>0) {
cvrt$LoadFile(covariates_file_name)
}

# outcome specific difference??
me = Matrix_eQTL_engine(
    snps = snps,
    gene = gene,
    cvrt = cvrt,
    output_file_name = output_file_name,
    pvOutputThreshold = pvOutputThreshold,
    useModel = useModel,
    errorCovariance = errorCovariance,
    verbose = TRUE,
    pvalue.hist = TRUE,
    min.pv.by.genesnp = FALSE,
    noFDRsaveMemory = FALSE)

final_file <- paste0(CHR, "final_eQTL_out.txt") 
write.table(me$all$eqtls, file = final_file, sep = "\t")

