#!/bin/bash
#########################################################################
# Copyright (c) 2019-~ Beisi Xu
# This source code is released for free distribution under the terms of the
# CreativeCommons BY-NC-SA 4.0 International License
#*Author:       Beisi Xu < xubeisi [at] gmail DOT com >
# File Name: atac_count.sh
#########################################################################

fppe=$1 # sorted bedpe file
if [ $# -gt 1 ]
then
fpe=$2
else
fpe=${fppe}.dat
fi

if [ $# -gt 2 ]
then
bedpe=$mode
else
bedpe=bedpe
fi

gettype(){
ff=$1
thetype=$2
bedpe=$3

case $thetype in
    centi0)
    min=0
    max=150
    ;;
    centi1)
    min=0
    max=180
    ;;
    free)
    min=0
    max=100
    ;;
    mono)
    min=180
    max=247
    ;;
    di)
    min=315
    max=473
    ;;
    tri)
    min=558
    max=615
    ;;
    all)
    min=0
    max=20000
    ;;
    2k)
    min=0
    max=2000
    ;;
    *)
    min=0
    max=2000
    ;;
esac

if [[ $bedpe =~ bedse ]]
then
awk "BEGIN{count=0;} \$3 - \$2 > $min - 9 && \$3 - \$2 < $max + 9 && \$1 !~ /chrM/{count+=1;} END{printf count\"\t\";}" $ff
else
awk "BEGIN{count=0;} \$6 - \$2 > $min -9 && \$6 - \$2 < $max + 9 && \$1 !~ /chrM/{count+=1;} END{printf count\"\t\";}" $ff
fi
}

foo(){

printf $(basename $fppe)"\t" | sed "s/\.\/Bam\/SICER\///;s/\.bedpe//;s/\.bed//" > ${fpe}

for i in all 2k free mono di tri
do
gettype $fppe $i $bedpe
done >> ${fpe}

if [[ $bedpe =~ bedse ]]
then
awk 'BEGIN{count=0;} $3-$2 > 0 && $1 ~ /chrM/{count+=1;} END{printf count"\t";}' $fppe >> ${fpe}
awk '$3-$2 > 0 && $1 !~ /chrM/{print $3-$2}' $fppe > ${fpe}.fsizes
else
awk 'BEGIN{count=0;} $6-$2 > 0 && $1 ~ /chrM/{count+=1;} END{printf count"\t";}' $fppe >> ${fpe}
awk '$6-$2 > 0 && $1 !~ /chrM/{print $6-$2}' $fppe > ${fpe}.fsizes
fi

echo "" >> ${fpe}
}

pppp(){
ff=$1
fff=$(echo $ff | sed "s/.*\///;s/\.bedpe//;s/\.bed//;s/\.fsizes//;s/\.dat//")
echo "
aa <- read.table('$ff')
aa\$FragmentSize <- aa\$V1
bb <- aa[aa\$FragmentSize < 1200 & aa\$FragmentSize > 0,]
if(!require(\"ggplot2\",character.only = TRUE)){
png('${ff}.dis.png')
plot(density(bb\$FragmentSize), main='$fff', xlab = 'Fragment Size', ylab = 'Density')
dev.off()
} else {
p <- ggplot(bb, aes(FragmentSize)) + geom_density() + ggtitle('$fff')
ggsave(filename='${ff}.dis.png',plot=p)
}
"
}

foo

pppp ${fpe}.fsizes > ${fpe}.r
Rscript ${fpe}.r
rm ${fpe}.r

