## compile reproducible peaks
#take 3 replicates in WT as example
cd /research_jude/rgs01_jude/groups/sapkogrp/projects/sapkogrp_auto/common/cab/ATACSEQ/SAPKO-322103-Northwestern-ATACSEQ/
cd /research_jude/rgs01_jude/groups/sapkogrp/projects/CAB/common/ATAC_seq_Northwestern/ATAC_seq_Northwestern_CAB_processed/Peaks
module load bedtools
for sample in $(ls *FDR50_free_macs2_peaks.narrowPeak | cut -d'.' -f1); do
echo "Doing ${sample}"
bedtools intersect -a ${sample}.free_macs2_peaks.narrowPeak -b ${sample}.FDR50_free_macs2_peaks.narrowPeak -u > ${sample}.kept.tmp.bed
done

cat $(ls *.kept.tmp.bed| sort -V) | sortBed -i | mergeBed -i stdin > All.kept.Repro.bed
rm  .tmp.bed


#Paired-End before get reads.
module load biobambam # for bamsort
module load samtools
for sample in $(ls *FDR50_free_macs2_peaks.narrowPeak | cut -d'.' -f1); do
samtools view -F 1804 -b -q 1 /research_jude/rgs01_jude/groups/sapkogrp/projects/sapkogrp_auto/common/cab/ATACSEQ/SAPKO-322103-Northwestern-ATACSEQ/${sample}/*/${sample}.bam | bamsort SO=queryname > ${sample}.q1.bam
bedtools bamtobed -bedpe -i ${sample}.q1.bam  2> /dev/null | awk '$NF != $(NF-1) && $1 == $4{ sub(/^/,"chr",$1); sub("chrchr","chr",$1); sub(/^/,"chr",$4); sub("chrchr","chr",$4); $7="r_"NR; OFS="\t"; print $0}' | grep -v chrM > ${sample}.bedpe
done

###############################
## reproducible peak regions ##
###############################
#ATACseq Paired-End, recommend only use nucleosome free fragments. See "ATACseq QC and peak calling" for details
for sample in $(ls *FDR50_free_macs2_peaks.narrowPeak | cut -d'.' -f1); do
echo "Doing ${sample}"
cat ${sample}.bedpe | awk '$6 - $2 < 109{OFS="\t"; print $1,$2,$6}' > ${sample}.bed
bedtools intersect -a All.kept.Repro.bed -b ${sample}.bed -c | awk '{print $NF}' > ${sample}.bed.count
done



#Combine into count table
awk '{print $1":"$2"-"$3}' All.kept.Repro.bed > All.kept.Repro.bed.pre

ls *.bed.count| sort -V

echo "Region $(ls *.bed.count | sort -V | cut -d'.' -f1 | tr '\n' ' ')" > All.counts.dat
paste All.kept.Repro.bed.pre $(ls *.bed.count | sort -V)  >> All.counts.dat


# #######################
# ## With wanted sites ##
# #######################

# #ATACseq Paired-End, recommend only use nucleosome free fragments. See "ATACseq QC and peak calling" for details
# for sample in $(ls *FDR50_free_macs2_peaks.narrowPeak | cut -d'.' -f1); do
# echo "Doing ${sample}"
# # cat ${sample}.bedpe | awk '$6 - $2 < 109{OFS="\t"; print $1,$2,$6}' > ${sample}.bed
# bedtools intersect -a wanted.sites.bed -b ${sample}.bed -c | awk '{print $NF}' > ${sample}.wanted.sites.bed.count
# done


# #Combine into count table
# awk '{print $1":"$2"-"$3}' wanted.sites.bed > wanted.sites.bed.pre

# ls *.wanted.sites.bed.count| sort -V

# echo "Region $(ls *.wanted.sites.bed.count | sort -V | cut -d'.' -f1 | tr '\n' ' ')" > All.wanted.sites.counts.dat
# paste wanted.sites.bed.pre $(ls *.wanted.sites.bed.count | sort -V)  >> All.wanted.sites.counts.dat





#############################################
# Can you examine the ATAC-seq peak data for the variants of interest, such as the TTN or BAG region and the new locus from Kateryna’s new GWAS, along with previously reported GWAS variants for cardiomyopathy in the Cancer Research (chr1 locus) and JNCI (chr6 locus) papers?
# We need to see if the peaks differ by cardiomyopathy status, with and without the variants of interest.

Sapkota et al 2022; JNCI (EUR)
rs2815063=chr6:39294759:C:A


Sapkota et al 2020; Cancer Research (AFR)
rs6689879=chr1:113559097:T:C

Among survivors of African ancestry; the rs6689879 is an eQTL for MAGI3, PHTF1 (putative 
homeodomain transcription factor 1), and two non-coding RNAs (RP5–1073O3.2 and RP4–
730K3.3) and a sQTL for PHTF1 and DCLRE1B (DNA cross-link repair 1B) in multiple 
human tissues (Supplementary Figures 4 and 5).

# PHTF1 promoter region was hypomethylated (chr1:113696831-113759486);
grep -E 'chr1:11369[0-9]{4,}|chr1:1137[0-9]{5,}'  Cardiotox_No_VS_Yes_afr_dose_0_diff.txt
grep -E 'chr1:11369[0-9]{4,}|chr1:1137[0-9]{5,}'  Cardiotox_No_VS_Yes_afr_dose_1_diff.txt
grep -E 'chr1:11369[0-9]{4,}|chr1:1137[0-9]{5,}'  Cardiotox_No_VS_Yes_afr_dose_3_diff.txt
# MAGI3: chr1:113390515-113685923; 
grep -E 'chr1:11339[0-9]{4,}|chr1:1134[0-9]{5,}|chr1:1135[0-9]{5,}|chr1:1136[0-9]{5,}' Cardiotox_No_VS_Yes_afr_dose_0_diff.txt
grep -E 'chr1:11339[0-9]{4,}|chr1:1134[0-9]{5,}|chr1:1135[0-9]{5,}|chr1:1136[0-9]{5,}' Cardiotox_No_VS_Yes_afr_dose_1_diff.txt
grep -E 'chr1:11339[0-9]{4,}|chr1:1134[0-9]{5,}|chr1:1135[0-9]{5,}|chr1:1136[0-9]{5,}' Cardiotox_No_VS_Yes_afr_dose_3_diff.txt

######################################
## Check for TTN and BAG3 variants: ##
######################################
## TTN: chr2:178525989-178807423
grep -E 'chr2:1785[0-9]{5,}|chr2:1786[0-9]{5,}|chr2:1787[0-9]{5,}|chr2:1788[0-9]{5,}'  Cardiotox_No_VS_Yes_afr_dose_0_diff.txt
logFC					AveExpr					t 						P.Value					adj.P.Value				B 						loci	
-0.942179162507886      1.99889877589193        -1.98046140548291       0.0609632455183909      0.999961858498848       -4.53461256140843	chr2:178804518-178805206
grep -E 'chr2:1785[0-9]{5,}|chr2:1786[0-9]{5,}|chr2:1787[0-9]{5,}|chr2:1788[0-9]{5,}'  Cardiotox_No_VS_Yes_afr_dose_1_diff.txt
grep -E 'chr2:1785[0-9]{5,}|chr2:1786[0-9]{5,}|chr2:1787[0-9]{5,}|chr2:1788[0-9]{5,}'  Cardiotox_No_VS_Yes_afr_dose_3_diff.txt

## BAG3:chr10:119651380-119677819
grep -E 'chr10:1196[0-9]{5,}'  Cardiotox_No_VS_Yes_afr_dose_0_diff.txt
grep -E 'chr10:1196[0-9]{5,}'  Cardiotox_No_VS_Yes_afr_dose_1_diff.txt

logFC					AveExpr					t 						P.Value					adj.P.Value				B 						loci	
-1.15415992375253       0.743022077383263       -2.42802930160121       0.0255413433767928      0.999981350798048       -4.49792974263232       chr10:119644485-119644615
# logFC indicates that the feature is downregulated in the Cardtox "Yes" condition compared to the "No".

# Positive B: A positive B-statistic indicates that the gene is more likely to be upregulated in the condition of interest compared to the reference condition.
# Negative B: A negative B-statistic indicates that the gene is more likely to be downregulated in the condition of interest compared to the reference condition.
# Magnitude of B: The magnitude of B reflects the strength of evidence for differential expression. Larger magnitudes (positive or negative) indicate stronger evidence.

grep -E 'chr10:1196[0-9]{5,}'  Cardiotox_No_VS_Yes_afr_dose_3_diff.txt
#################################################
## Check for Kateryna's dyslipidemia variants: ##
#################################################
# HS2ST1:  chr1:86914635-87098445
grep -E 'chr1:869[0-9]{5,}|chr1:870[0-9]{5,}'  Cardiotox_No_VS_Yes_afr_dose_0_diff.txt
grep -E 'chr1:869[0-9]{5,}|chr1:870[0-9]{5,}'  Cardiotox_No_VS_Yes_afr_dose_1_diff.txt
grep -E 'chr1:869[0-9]{5,}|chr1:870[0-9]{5,}'  Cardiotox_No_VS_Yes_afr_dose_3_diff.txt

logFC					AveExpr					t 						P.Value					adj.P.Value				B 						loci	
0.379035881449452       3.71791309721367        2.29170452954945        0.0265005671642148      0.99999531325849        -4.08174094334353       chr1:86913850-86915848

# HS6ST1: chr2:128265480-128318868
grep -E 'chr2:12823[0-9]{4,}|chr2:12824[0-9]{4,}|chr2:12825[0-9]{4,}|chr2:12826[0-9]{4,}'  Cardiotox_No_VS_Yes_afr_dose_0_diff.txt
grep -E 'chr2:12823[0-9]{4,}|chr2:12824[0-9]{4,}|chr2:12825[0-9]{4,}|chr2:12826[0-9]{4,}'  Cardiotox_No_VS_Yes_afr_dose_1_diff.txt
grep -E 'chr2:12823[0-9]{4,}|chr2:12824[0-9]{4,}|chr2:12825[0-9]{4,}|chr2:12826[0-9]{4,}'  Cardiotox_No_VS_Yes_afr_dose_3_diff.txt

# HS6ST3: chr13:96090107-96839562
grep -E 'chr13:96[0-9]{6,}' Cardiotox_No_VS_Yes_afr_dose_0_diff.txt
grep -E 'chr13:96[0-9]{6,}' Cardiotox_No_VS_Yes_afr_dose_1_diff.txt
grep -E 'chr13:96[0-9]{6,}' Cardiotox_No_VS_Yes_afr_dose_3_diff.txt

logFC					AveExpr					t 						P.Value					adj.P.Value				B 						loci	
0.698169639078131       0.83467008480492        2.05888570329865        0.0451337812817118      0.99999531325849        -4.47363848856412     chr13:96680056-96680372
-0.663719383781101      0.485882288776097       -1.98106585498666       0.0535179218248778      0.99999531325849        -4.48925518475091     chr13:96574462-96574639

# SULT1C3: chr2:108239968-108265351
grep -E 'chr2:10823[0-9]{4,}|chr2:10824[0-9]{4,}|chr2:10825[0-9]{4,}|chr2:10826[0-9]{4,}'  Cardiotox_No_VS_Yes_afr_dose_0_diff.txt
grep -E 'chr2:10823[0-9]{4,}|chr2:10824[0-9]{4,}|chr2:10825[0-9]{4,}|chr2:10826[0-9]{4,}'  Cardiotox_No_VS_Yes_afr_dose_1_diff.txt
grep -E 'chr2:10823[0-9]{4,}|chr2:10824[0-9]{4,}|chr2:10825[0-9]{4,}|chr2:10826[0-9]{4,}'  Cardiotox_No_VS_Yes_afr_dose_3_diff.txt

# GPC6: chr13:93226807-94408020
grep -E 'chr13:93[0-9]{6,}|chr13:94[0-9]{6,}' Cardiotox_No_VS_Yes_afr_dose_0_diff.txt
grep -E 'chr13:93[0-9]{6,}|chr13:94[0-9]{6,}' Cardiotox_No_VS_Yes_afr_dose_1_diff.txt
grep -E 'chr13:93[0-9]{6,}|chr13:94[0-9]{6,}' Cardiotox_No_VS_Yes_afr_dose_3_diff.txt
logFC				AveExpr					t 						P.Value					adj.P.Value				B 						loci	
0.830087363997164       0.402579823639871       2.11098320877508        0.0401816329577544      0.99999531325849        -4.51131239655434	chr13:93296946-93297094


