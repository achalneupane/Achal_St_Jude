#!/usr/local/bin/R
#this script plots observed heterozygosity per sample as a function of proportion of missingness

#Read data

# system('load module R/4.2.2-rhel8')
setwd("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/Survivor_WGS_QCed/QC/")

imiss = read.table("chrALL.Survivor_WGS.GATK4180.hg38_missingness.imiss", header=TRUE)
hetCalc = read.table("Survivor_WGS.GATK4180.hg38_renamed_chr1-22.PASS.decomposed_common_pruned_indep_biallelic_varnames_updated_het.het", header=TRUE)
# hetCalc$INDV = paste(hetCalc$FID, hetCalc$IID, sep="_")
# hetCalc$INDV[1:4] = as.character(hetCalc$FID)

hetCalc$het<-(hetCalc$N.NM-hetCalc$O.HOM)/hetCalc$N.NM
mean_het = mean(hetCalc$het)
sd_het = sd(hetCalc$het)
lower_het = mean_het-3*sd_het
upper_het = mean_het+3*sd_het
outlier = hetCalc[hetCalc$het<lower_het | hetCalc$het>upper_het,]
# dim(outlier)
# [1] 128   7

write.table(outlier, 'Survivor_WGS_Per_sample_heterozygosity_outlier_check.txt', row.names=F, quote=F, sep="\t")

m = merge(imiss, hetCalc, by="IID")

#pdf("Per-sample_heterozygosity_rate_vs_proportion_of_missing_genotypes.pdf")
pdf("Per-sample_heterozygosity_rate_vs_proportion_of_missing_genotypes.pdf")
plot(m$F_MISS,m$het, pch=20, xlab="Proportion of missing genotypes", ylab="Heterozygosity rate")
abline(h=mean(hetCalc$het),col="red",lwd=2)
abline(h=mean(hetCalc$het)-(3*sd(hetCalc$het)),col="blue",lty=2)
abline(h=mean(hetCalc$het)+(3*sd(hetCalc$het)),col="blue",lty=2)
dev.off()

