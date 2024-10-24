#!/bin/bash

myDir=/home/ysapkota/Work/WGS_SJLIFE/Analyses
cd $myDir

for chr in {1..22};do
#for chr in 22;do
 bsub -q priority -P SJLFIEGWAS -J annovar_chr$chr -eo ../log/annovar_chr$chr.err -oo ../log/annovar_chr$chr.out -R "rusage[mem=3000]" bash scripts/run_annovar.sh $chr
 for g in PTV NS NS1 REGMOTIF PTV_NS_NS1_REGMOTIF GENIC;do
  #bsub -P SJLFIEGWAS -J group_chr$chr -eo log/group_chr$chr.err -oo log/group_chr$chr.out -R "rusage[mem=2000]" python bed2groupfile.py annotations/chr$chr.$g
  bsub -q priority -P SJLFIEGWAS -w "done(annovar_chr$chr)" -J group_chr$chr -eo ../log/group_chr$chr.err -oo ../log/group_chr$chr.out -R "rusage[mem=2000]" python scripts/bed2groupfile.py annotations/chr$chr.$g
 done
done
