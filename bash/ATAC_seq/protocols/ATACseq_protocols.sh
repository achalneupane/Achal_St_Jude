1.2 Reference data and annotation
In our standardized pipelines, the reference genome sequence and annotation files are detailed in the following table. The files are accessible to everyone at St Jude. See Reference Library for more. 

ls /research/rgs01/applications/hpcf/authorized_apps/cab/Automation/REF
Homo_sapiens
Mus_musculus
Drosophila_melanogaster
Danio_rerio
Mustela_putorius_furo
Schizosaccharomyces_pombe

1.3 Preprocessing of ATACseq reads
Trim_Galore(v0.4.4) is used to quality trim 3’ adapters from raw reads with the cutadapt (DOI: 10.14806/ej.17.1.200) program and perform FASTQC (Andrews 2010) analysis.  A quality score cutoff of Q20 is used. The first 15bp of each reads were also clipped to reduce the GC bias.


Key example commands Collapse source
trim_galore --paired --fastqc --gzip --clip_R1 15 --clip_R2 15 ATAC_R1.fq.gz ATAC_R2.fq.gz

1.4 Mapping of ATACseq reads
Quality trimmed reads are mapped to the selected reference genome by bwa(v0.7.12-r1039, PMID: 19451168) and convert to bam file (a binary version of sam file) by samtools (v1.2, PMID: 19505943). biobambam2 (v2.0.87, DOI: 10.1186/1751-0473-9-13) were used to mark duplicated reads. 

Key example commands Collapse source
#old bam files from automapper(before 2020 March) might use bwa-aln, check bam header to make sure
#bwa aln ref.fa ATAC_R1_val_1.fq.gz > ATAC_R1.fq.gz.sai
#bwa aln ref.fa ATAC_R2_val_2.fq.gz > ATAC_R2.fq.gz.sai
#bwa sampe ref.fa ATAC_R1.fq.gz.sai ATAC_R2.fq.gz.sai ATAC_R1_val_1.fq.gz ATAC_R2_val_2.fq.gz | samtools view -b - > ATAC.unmark.bam
 
#bwa mem is quicker and even more quicker with our GPU-based(20x~45x) speedup so we switched to bwa mem in automapper. In our comparisons the QC pipeline didn't show strong differences. 
bwa mem ref.fa ATAC_R1_val_1.fq.gz ATAC_R2_val_2.fq.gz | samtools view -b - > ATAC.unmark.bam
 
cat ATAC.unmark.bam | bamsormadup > ATAC.bam


1.5 Fragment Size Distribution Analysis
Properly paired uniquely mapped reads were extracted by samtools (v1.2, PMID: 19505943). Fragments were extracted by bedtools (v2.24.0, PMID: 20110278) and fragment sizes distribution were plot by custom R script(in atac_count.sh).


Key example commands Collapse source
samtools view -F 1804 -b -q 1 ATAC.bam | bamsort SO=queryname > Bam/ATAC.q1.bam
bamToBed -bedpe -i Bam/ATAC.q1.bam  2> /dev/null | awk '$NF == "-" && $(NF-1) == "+" && $1 == $4{ sub(/^/,"chr",$1); sub("chrchr","chr",$1); sub("chrMT","chrM",$4); sub(/^/,"chr",$4); sub("chrchr","chr",$4); sub("chrMT","chrM",$4); $7="r_"NR; OFS="\t"; print $0}' | grep -v chrM > Bam/ATAC.bedpe
atac_count.sh Bam/ATAC.bedpe
# this would create Bam/ATAC.bedpe.dat each column is number of fragments as:
# filename all_fragment <2kb nucleosome_free mono-nucleosome di-nucleosome tri-nucleosome



scalefactor=$(head -n 1 Bam/ATAC.bedpe.dat | awk '{print 20000000/$4}')
#nucleosome free
cat Bam/ATAC.bedpe | awk 'BEGIN{OFS="\t";} $6 - $2 < 109 && $2+$6 > 80 && $6>$2 {mid=int(($2+$6)/2); print $1,mid-40,mid+41,$7,$8,"+"}' > Bam/QC/Coverage/ATAC.bed # 109 == 100 + Tn5 adjustment
sort-bed --max-mem 9G Bam/QC/Coverage/ATAC.bed > Bam/QC/Coverage/ATAC.free.bed #module load bedops
genomeCoverageBed -bg -i Bam/QC/Coverage/ATAC.free.bed -g genome.sizes -scale $scalefactor > Bam/QC/Coverage/ATAC.sort.bedGraph
bedGraphToBigWig Bam/QC/Coverage/ATAC.sort.bedGraph genome.sizes Bam/QC/Coverage/ATAC.free.bw
 
#mono/di/tri-nucleosomes
cat Bam/ATAC.bedpe | awk 'BEGIN{OFS="\t";} $NF == "-" && $(NF-1) == "+" && $1 == $4 && $6 - $2 < 624 && $6-$2 > 171{mid=int(($2+$6)/2); print $1,mid-40,mid+41}' > Bam/QC/Coverage/ATAC.bed
sort-bed --max-mem 9G Bam/QC/Coverage/ATAC.bed > Bam/QC/Coverage/ATAC.sort.bed
genomeCoverageBed -bg -i Bam/QC/Coverage/ATAC.sort.bed -g genome.sizes -scale $scalefacto > Bam/QC/Coverage/ATAC.sort.bedGraph
bedGraphToBigWig Bam/QC/Coverage/ATAC.sort.bedGraph genome.sizes Bam/QC/Coverage/ATAC.nuc.bw


1.7 Peak Calling
Nucleosome free fragment were extracted by bedtools (v2.24.0, PMID: 20110278). Narrow peaks were called by MACS2 (v2.1.1.20160309, PMID: 18798982) on nucleosome free reads (fragment size < 100bp, find out why @2.2) using False Positive Rate corrected p-value 0.05 (find out more @2.3). To review the peaks called, download IGV Genomics Browser



Key example commands
### Updated version (after 2021-07-16): # using center 80bp, peaks called more consistent with bigwig
macs2 callpeak -t Bam/QC/Coverage/ATAC.free.bed -f BEDPE -g hs --keep-dup all -n Bam/Peaks/ATAC.free -q 0.05
macs2 callpeak -t Bam/QC/Coverage/ATAC.free.bed -f BEDPE -g hs --keep-dup all -n Bam/Peaks/ATAC.FDR50_free -q 0.5
 
### Old version (before 2021-07-16):
#cat Bam/ATAC.bedpe | awk '$NF == "-" && $6 - $2 < 109{OFS="\t"; print $1,$2,$3 "\n" $4,$5,$6}' > Bam/ATAC.freepe.bed
#macs2 callpeak -t Bam/ATAC.freepe.bed -f BED --nomodel --keep-dup all -g hs --extsize 200 -n Bam/Peaks/ATAC.free -q 0.05
#macs2 callpeak -t Bam/ATAC.freepe.bed -f BED --nomodel --keep-dup all -g hs --extsize 200 -n Bam/Peaks/ATAC.FDR50_free -q 0.5
 
### Old version (before 2020-01)
#cat Bam/ATAC.bedpe | awk '$NF == "-" && $(NF-1) == "+" && $1 == $4 && $6 - $2 < 109{OFS="\t"; print $1,$2,$6,$7,$8,"+"}' > Bam/ATAC.free.bedpe
#macs2 callpeak -t Bam/ATAC.free.bedpe -f BEDPE --extsize 200 --nomodel --keep-dup all -g hs -n Bam/ATAC.free -q 0.05