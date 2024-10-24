rm(list=ls())
setwd("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics//common/Survivor_WGS_QCed/QC/pca")
# Population labels for 1kGP
refPopLabels <- read.table('Z:/ResearchHome/Groups/sapkogrp/projects/Genomics//common/1kGP/integrated_call_samples_v3.20130502.ALL.panel', header=TRUE)
refPopLabels.subset <- subset(refPopLabels,select=c('sample','super_pop'))
colnames(refPopLabels.subset) <- c('sample', 'pop')

# Load the PCA data including both test and ref populations
pcs <- read.table('Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common//Survivor_WGS_QCed/QC/pca/Survivor_WGS_and_1kGP_final_top_20_PCs.eigenvec', header=FALSE)
colnames(pcs) <- c('FID','IID',paste('PC',1:20,sep=""))
pcs1_2 <- subset(pcs, select=c('FID', 'IID','PC1','PC2'))
pcs1_2.label <- merge(pcs1_2,refPopLabels.subset,by.x='FID',by.y='sample',all.x =TRUE)
pcs1_2.label$pop <- as.character(pcs1_2.label$pop)
pcs1_2.label$pop[is.na(pcs1_2.label$pop)] <- "Survivor_WGS"

pcs1_2.label.EUR <- subset(pcs1_2.label, pcs1_2.label$pop=="EUR")
mean_PC1.EUR <- mean(pcs1_2.label.EUR$PC1, na.rm=TRUE)
sd_PC1.EUR <- sd(pcs1_2.label.EUR$PC1, na.rm=TRUE)
mean_PC2.EUR <- mean(pcs1_2.label.EUR$PC2, na.rm=TRUE)
sd_PC2.EUR <- sd(pcs1_2.label.EUR$PC2, na.rm=TRUE)
lower_PC1.EUR <- mean_PC1.EUR-3*sd_PC1.EUR
upper_PC1.EUR <- mean_PC1.EUR+3*sd_PC1.EUR
lower_PC2.EUR <- mean_PC2.EUR-3*sd_PC2.EUR
upper_PC2.EUR <- mean_PC2.EUR+3*sd_PC2.EUR

pcs1_2.label.AMR <- subset(pcs1_2.label, pcs1_2.label$pop=="AMR" | pcs1_2.label$pop=="EUR")
mean_PC1.AMR <- mean(pcs1_2.label.AMR$PC1, na.rm=TRUE)
sd_PC1.AMR <- sd(pcs1_2.label.AMR$PC1, na.rm=TRUE)
mean_PC2.AMR <- mean(pcs1_2.label.AMR$PC2, na.rm=TRUE)
sd_PC2.AMR <- sd(pcs1_2.label.AMR$PC2, na.rm=TRUE)
lower_PC1.AMR <- mean_PC1.AMR-3*sd_PC1.AMR
upper_PC1.AMR <- mean_PC1.AMR+3*sd_PC1.AMR
lower_PC2.AMR <- mean_PC2.AMR-3*sd_PC2.AMR
upper_PC2.AMR <- mean_PC2.AMR+3*sd_PC2.AMR


pcs1_2.label.AFR <- subset(pcs1_2.label, pcs1_2.label$pop=="AFR")
mean_PC1.AFR <- mean(pcs1_2.label.AFR$PC1, na.rm=TRUE)
sd_PC1.AFR <- sd(pcs1_2.label.AFR$PC1, na.rm=TRUE)
mean_PC2.AFR <- mean(pcs1_2.label.AFR$PC2, na.rm=TRUE)
sd_PC2.AFR <- sd(pcs1_2.label.AFR$PC2, na.rm=TRUE)
lower_PC1.AFR <- mean_PC1.AFR-3*sd_PC1.AFR
upper_PC1.AFR <- mean_PC1.AFR+3*sd_PC1.AFR
lower_PC2.AFR <- mean_PC2.AFR-3*sd_PC2.AFR
upper_PC2.AFR <- mean_PC2.AFR+3*sd_PC2.AFR

# SJLFIE samples that are not EURs
test = subset(pcs1_2.label, pop == "Survivor_WGS")
# Samples that are not EURs
tmp1.nonEUR <- test[test$PC1 < lower_PC1.EUR,]
tmp2.nonEUR <- test[test$PC1 > upper_PC1.EUR,]
tmp3.nonEUR <- test[test$PC2 < lower_PC2.EUR,]
tmp4.nonEUR <- test[test$PC2 > upper_PC2.EUR,]
tmp.nonEUR <- rbind(tmp1.nonEUR, tmp2.nonEUR, tmp3.nonEUR, tmp4.nonEUR)
tmp.nonEUR <- unique(tmp.nonEUR)
# Samples that are EURs
tmp.EUR = subset(test, !(test$FID %in% tmp.nonEUR$FID))
write.table(tmp.EUR, 'Survivor_WGS_EUR_based_on_1kGP_Phase_3_data.txt',row.names=FALSE,quote=FALSE)

# Samples that are not AMRs
tmp1.nonAMR <- test[test$PC1 < lower_PC1.AMR,]
tmp2.nonAMR <- test[test$PC1 > upper_PC1.AMR,]
tmp3.nonAMR <- test[test$PC2 < lower_PC2.AMR,]
tmp4.nonAMR <- test[test$PC2 > upper_PC2.AMR,]
tmp.nonAMR <- rbind(tmp1.nonAMR, tmp2.nonAMR, tmp3.nonAMR, tmp4.nonAMR)
tmp.nonAMR <- unique(tmp.nonAMR)
# Samples that are AMRs
tmp.AMR = subset(test, !(test$FID %in% tmp.nonAMR$FID))
write.table(tmp.AMR, 'Survivor_WGS_EUR_AMR_based_on_1kGP_Phase_3_data.txt',row.names=FALSE,quote=FALSE)

# Samples that are not AFRs
tmp1.nonAFR <- test[test$PC1 < lower_PC1.AFR,]
tmp2.nonAFR <- test[test$PC1 > upper_PC1.AFR,]
tmp3.nonAFR <- test[test$PC2 < lower_PC2.AFR,]
tmp4.nonAFR <- test[test$PC2 > upper_PC2.AFR,]
tmp.nonAFR <- rbind(tmp1.nonAFR, tmp2.nonAFR, tmp3.nonAFR, tmp4.nonAFR)
tmp.nonAFR <- unique(tmp.nonAFR)
# Samples that are AFRs
tmp.AFR = subset(test, !(test$FID %in% tmp.nonAFR$FID))
write.table(tmp.AFR, 'Survivor_WGS_AFR_based_on_1kGP_Phase_3_data.txt',row.names=FALSE,quote=FALSE)

# Create a PCA plot
#pdf('~/Work/WGS_SJLIFE/QC/pca/PCA_plot_SJLIFEWGS_1kGP.pdf')
pdf('PCA_plot_SJLIFEWGS_1kGP_AFR.pdf')
palette(c("black","lightgrey","purple","lightgreen","yellow","orange","skyblue2","red","skyblue2","slateblue"))
# All 3006 samples
plot(pcs1_2.label$PC1, pcs1_2.label$PC2, main="SJLIFE WGS samples", col="black", bg=as.factor(pcs1_2.label$pop), pch=21, cex=1.2, xlab="Principal Component 1", ylab="Principal Component 2")
legend('topright', levels(as.factor((pcs1_2.label$pop))), fill = 1:10)
abline(v=lower_PC1.AFR, lty=2, col="red")
abline(v=upper_PC1.AFR, lty=2, col="red")
abline(h=lower_PC2.AFR, lty=2, col="blue")
abline(h=upper_PC2.AFR, lty=2, col="blue")
#abline(v=lower_PC1.AFR, lty=2, col="cyan")
#abline(v=upper_PC1.AFR, lty=2, col="cyan")
#abline(h=lower_PC2.AFR, lty=2, col="orange")
#abline(h=upper_PC2.AFR, lty=2, col="orange")
# EURs alone
pcs1_2.label_SJLIFE_EUR = subset(pcs1_2.label, !(pcs1_2.label$FID %in% tmp.nonEUR$FID))
plot(pcs1_2.label_SJLIFE_EUR$PC1, pcs1_2.label_SJLIFE_EUR$PC2, main="Survivor_WGS samples [Europeans alone]", col="black", bg=as.factor(pcs1_2.label_SJLIFE_EUR$pop), pch=21, cex=1.2, xlab="Principal Component 1", ylab="Principal Component 2")
legend('topright', levels(as.factor((pcs1_2.label_SJLIFE_EUR$pop))), fill = 1:10)
# AFRs alone
pcs1_2.label_SJLIFE_AFR = subset(pcs1_2.label, !(pcs1_2.label$FID %in% tmp.nonAFR$FID))
plot(pcs1_2.label_SJLIFE_AFR$PC1, pcs1_2.label_SJLIFE_AFR$PC2, main="Survivor_WGS samples [Africans alone]", col="black", bg=as.factor(pcs1_2.label_SJLIFE_AFR$pop), pch=21, cex=1.2, xlab="Principal Component 1", ylab="Principal Component 2")
legend('topright', levels(as.factor((pcs1_2.label_SJLIFE_AFR$pop))), fill = 1:10)
dev.off()
