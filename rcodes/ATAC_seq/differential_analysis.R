setwd("Z:/ResearchHome/Groups/sapkogrp/projects//CAB/common/ATAC_seq_Northwestern/ATAC_seq_Northwestern_CAB_processed/Peaks")
library(haven)
library(benchmarkme)
library(dplyr)
library(plyr)
library(data.table)
library (birk)

data.lastcondt <- read_sas("Z:/SJShare/SJCOMMON/ECC/SJLife/SJLIFE Data Freeze/2 Final Data SJLIFE/20200430/Tracking Data/lstcondt.sas7bdat")


pheno <- read.table(text = "SJLID	Sample_id	EF	Cardtox	Genotypes	Cumulative anthracycline dose	Ancestry	Sex	Age_at_treatment	diaggrp	Age_at_last_contact
SJL1684113	TB-16-02775	0.449	Yes	TC	450.459	AFR	Male	11.4739726	Osteosarcoma	52.18356164
SJL1793013	TB-12-5170	0.574	No	TT	473.846	AFR	Male	14.29315068	Osteosarcoma	60.90410959
SJL2508708	TB-16-06955	0.425	Yes	TC	101.841	AFR	Female	15.13448611	Non-Hodgkin lymphoma	46.38858447
SJL4181813	TB-19-04682	0.445	Yes	TC	224.074	AFR	Female	11.075582	Osteosarcoma	34.68928063
SJL4187108	TB-14-7597	0.559	No	TT	461.765	AFR	Male	7.583561644	Non-Hodgkin lymphoma	40.32525638
SJL4726013	TB-13-6124	0.577	No	TT	291.489	AFR	Female	16.12054795	Osteosarcoma	39.86849315
SJL5049617	TB-14-0242	0.567	No	TT	365.839	AFR	Female	17.85205479	Soft tissue sarcoma	46.11780822
SJL5154909	TB-14-1537	0.472	Yes	TC	185.185	AFR	Male	2.071232877	Wilms tumor	28.93683659
SJL5188113	TB-17-04138	0.408	Yes	TC	391.228	AFR	Female	14.56164384	Osteosarcoma	35.71506849
SJL5237613	TB-15-9152	0.617	No	TT	384.211	AFR	Female	13.2739726	Osteosarcoma	27.20821918
SJL1691010	TB-14-5054	0.436	Yes	AA	181.482	EUR	Female	0.008219178	Neuroblastoma	37.14185942
SJL1741608	TB-13-4405	0.36	Yes	AC	393.329	EUR	Male	5.35890411	Non-Hodgkin lymphoma	44.47945205
SJL4726613	TB-15-10427	0.395	Yes	AC	388.75	EUR	Female	13.76712329	Osteosarcoma	42.61047983
SJL4824312	TB-16-06040	0.634	No	CC	246.342	EUR	Female	15.23444869	Ewing sarcoma family of tumors	44.03005464
SJL5113211	TB-12-4295	0.463	Yes	AA	341.177	EUR	Male	1.233692642	Liver malignancies	26.61917808
SJL5140702	TB-12-3198	0.636	No	CC	292.053	EUR	Female	17.31040497	Acute myeloid leukemia	36.19452055
SJL5152913	TB-17-00238	0.396	Yes	AC	348	EUR	Male	17.67555206	Osteosarcoma	36.45363425
SJL5195202	TB-13-2431	0.659	No	CC	306.084	EUR	Male	17.2524665	Acute myeloid leukemia	34.88534321
SJL5208908	TB-20-01253	0.704	No	CC	240.909	EUR	Male	18.6739726	Non-Hodgkin lymphoma	37.58334456
SJL5234413	TB-20-00725	0.674	No	CC	372.832	EUR	Female	20.51232877	Osteosarcoma	37.28481174", header = T, sep = "\t")


# pheno$new_agelstcont <- data.lastcondt$agelstcontact[match(pheno$SJLID, data.lastcondt$sjlid)]

pheno$ID <- sub("^[^-]+-[^-]+-", "", pheno$Sample_id)
pheno$ID[pheno$ID == "06955"] <- "60955"
pheno$ID[pheno$ID == "04682"] <- "4682"

###########################
## differential analysis ##
###########################
# load package
library(edgeR)
library(limma)
library(dplyr)
library(tidyr)

## read count data
counts <- read.table("All.counts.dat", header = TRUE, row.names = 1, check.names = F)
samples <- unique(gsub("-.*", "", colnames(counts)))
samples[!samples %in% pheno$ID]
# 0

# Rearrange pheno
pheno <- pheno[pheno$ID %in% samples,]

# Create a new dataframe with expanded rows; Replicate each row two more times
replicated_df <- pheno[rep(seq_len(nrow(pheno)), each = 3), ]

# Add ID-0, ID-1, and ID-3 columns
replicated_df$new_ID <- rep(paste(rep(pheno$ID, each = 3), c("0", "1", "3"), sep = "-"), times = 1)
pheno <- replicated_df; rm(replicated_df)

pheno$new_ID
# [1] "02775-0" "02775-1" "02775-3" "5170-0"  "5170-1"  "5170-3"  "60955-0" "60955-1" "60955-3" "4682-0"  "4682-1"  "4682-3"  "7597-0"  "7597-1"  "7597-3" 
# [16] "6124-0"  "6124-1"  "6124-3"  "0242-0"  "0242-1"  "0242-3"  "1537-0"  "1537-1"  "1537-3"  "9152-0"  "9152-1"  "9152-3" 

colnames(counts)
# [1] "0242-0"  "0242-1"  "0242-3"  "1537-0"  "1537-1"  "1537-3"  "02775-0" "02775-1" "02775-3" "4682-0"  "4682-1"  "4682-3"  "5170-0"  "5170-1"  "5170-3" 
# [16] "6124-0"  "6124-1"  "6124-3"  "7597-0"  "7597-1"  "7597-3"  "9152-0"  "9152-1"  "9152-3"  "60955-0" "60955-1" "60955-3"

sum(!colnames(counts) %in% pheno$new_ID)
# 0

pheno <- pheno[order(match(pheno$new_ID, colnames(counts))), ]
pheno$new_ID
# [1] "0242-0"  "0242-1"  "0242-3"  "1537-0"  "1537-1"  "1537-3"  "02775-0" "02775-1" "02775-3" "4682-0"  "4682-1"  "4682-3"  "5170-0"  "5170-1"  "5170-3" 
# [16] "6124-0"  "6124-1"  "6124-3"  "7597-0"  "7597-1"  "7597-3"  "9152-0"  "9152-1"  "9152-3"  "60955-0" "60955-1" "60955-3"
sum(pheno$new_ID != colnames(counts))
# 0

day1 <- c("60955-0", "60955-1", "60955-3", "5170-0", "5170-1", "5170-3", 
          "4682-0", "4682-1", "4682-3", "6124-0", "6124-1", "6124-3", "9152-0", 
          "9152-1", "9152-3")

day2 <- c("1537-0", "1537-1", "1537-3", "0242-0", "0242-1", "0242-3", 
          "7597-0", "7597-1", "7597-3", "02775-0", "02775-1", "02775-3"
)


pheno$dose <- sub(".*-", "", pheno$new_ID)

pheno$batch [pheno$new_ID %in% day1] <- "day1"
pheno$batch [pheno$new_ID %in% day2] <- "day2"

pheno.saved <- pheno
counts.saved <- counts

#####################
## Annotated peaks ##
#####################
annotated.peaks <- read.delim("All.counts.dat_annotated",  header = T, sep = "\t")
annotated.peaks <- annotated.peaks[-1]
annotated.peaks$Peaks <- paste0(annotated.peaks$Chr, ":", annotated.peaks$Start, "-", annotated.peaks$End)

####################


# pheno <- pheno.saved 
# counts <- counts.saved



# PHTF1 <- counts[grepl("chr1:11369|chr1:1137", rownames(counts)),]

doses <- c(0, 1, 3)
# dose = 1
for (dose in doses){
print (dose)
## Dox 0
pheno <- pheno.saved[pheno.saved$dose == dose,]
counts <- counts.saved[colnames(counts.saved) %in% pheno$new_ID]
sum(pheno$new_ID != colnames(counts))
# 0

# define group
groupLabels <- pheno$Cardtox
TS <- factor(groupLabels, levels = c("No", "Yes"))

#batch correction: batch can be added as a covariate if needed
#this code does not do batch correction
design <- model.matrix(~ 0 + TS )
colnames(design) <- levels(TS)

# ## adding covariates
# design <- model.matrix(~ 0 +TS + Sex + Age_at_treatment + Cumulative.anthracycline.dose, data = pheno)
# colnames(design) <- c(levels(TS), "Sex", "Age_at_treatment", "Cumulative.anthracycline.dose")

dge <- DGEList(counts = counts, group = groupLabels)

# filter out low expressed regions
cpm_cutoff=10
cutoff <- as.vector(cpm(cpm_cutoff,mean(dge$samples$lib.size) ) )
keep <- (rowSums(cpm(dge) > cutoff) >= min(as.numeric(table(groupLabels))))
dge <- dge[keep, keep.lib.sizes = FALSE]
dge <- calcNormFactors(dge)

v <- voom(dge, design, plot = F)
fit <- lmFit(v, design)
cont.matrix <- makeContrasts(NovsYes = (No - Yes), levels = design)
fitcon <- contrasts.fit(fit, cont.matrix)
fitcon <- eBayes(fitcon)
results1 <- topTable(fitcon, n = Inf, sort.by="P", coef="NovsYes")
results1$loci <- rownames(results1)

sum(results1$loci %in% annotated.peaks$Peaks)
results1 <- merge(results1, annotated.peaks, by.x = "loci", by.y ="Peaks")

outFile1 <- paste0("annotated_Cardiotox_No_VS_Yes_afr_dose_", dose, "_diff.txt")
write.table(results1, outFile1, col.names = T, row.names = F, sep = "\t", quote = F)

# significant ones
results1.sig <- results1[results1$P.Value < 0.05,]
outFile2 <- paste0("annotated_significant_Cardiotox_No_VS_Yes_afr_dose_", dose, "_diff.txt")
write.table(results1.sig, outFile2, col.names = T, row.names = F, sep = "\t", quote = F)

# ## With P
# gg <- ggplot(results1, aes(x = logFC, y = -log10(P.Value))) +
#   geom_point(aes(color = ifelse(P.Value < 0.05, "Yes", "No")), size = 3) +
#   scale_color_manual(values = c("Yes" = "red", "No" = "black")) +
#   theme_minimal() +
#   labs(title = "Volcano Plot", x = "logFC", y = "-log10(P.Value)", color = "Siginificant")
# 
# ggsave(paste0("Cardiotox_No_VS_Yes_afr_dose_", dose, "_diff.tiff"), gg, width = 8, height = 6, dpi = 300)

}


## Get 250 KB up/downstream of each gene of interest
dose0.sig <- read.delim("annotated_significant_Cardiotox_No_VS_Yes_afr_dose_0_diff.txt", header = T, sep = "\t")
dose0.sig$dose <- 0
dose1.sig <- read.delim("annotated_significant_Cardiotox_No_VS_Yes_afr_dose_1_diff.txt", header = T, sep = "\t")
dose1.sig$dose <- 1
dose3.sig <- read.delim("annotated_significant_Cardiotox_No_VS_Yes_afr_dose_3_diff.txt", header = T, sep = "\t")
dose3.sig$dose <- 3

all.doses.significant <- rbind.data.frame(dose0.sig, dose1.sig, dose3.sig)

# get genes expaning 250kb +-

all.doses.significant$chr=CHR & all.doses.significant$Start=START 


library(GenomicRanges)
library(IRanges)
# Assuming all.doses.significant$gene contains gene names
# HS2ST1:  chr1:86914635-87098445

wanted.genes <- read.table(text="GENE	CHR 	START	END
PHTF1	chr1	113696831	113759486
MAGI3	chr1	113390515	113685923
TTN	chr2	178525989	178807423
BAG3	chr10	119651380	119677819
HS2ST1	chr1	86914635	87098445
HS6ST1	chr2	128265480	128318868
HS6ST3	chr13	96090107	96839562
SULT1C3	chr2	108239968	108265351
GPC6	chr13	93226807	94408020", header= T)

all.doses.significant$gene <- NA
all.doses.significant$gene.start.end <- NA
for (i in 1:nrow(wanted.genes)){
GENE = wanted.genes$GENE[i]
CHR = wanted.genes$CHR[i]
START = wanted.genes$START[i]
END = wanted.genes$END[i]

Upstream <- pmax(1, START - 250000)
Downstream <- END + 250000
# Create a GRanges object for gene regions
gene_regions <- GRanges(seqnames = CHR,
                        ranges = IRanges(start = Upstream,
                                         end = Downstream))
# Create a GRanges object for genes
genes <- GRanges(seqnames = all.doses.significant$Chr,
                 ranges = IRanges(start = all.doses.significant$Start,
                                  end = all.doses.significant$End, loci = all.doses.significant$loci))  # Add gene information
# Find overlaps between gene regions and loci
overlaps <- findOverlaps(gene_regions, genes)
# Extract the corresponding gene names
overlapping_genes <- mcols(genes)$loci[subjectHits(overlaps)]

# Add the labels to all.doses.significant
all.doses.significant$gene[all.doses.significant$loci %in% overlapping_genes] <- GENE
all.doses.significant$gene.start.end[all.doses.significant$loci %in% overlapping_genes] <- paste0(CHR, ":", START, "-", END)
}

all.doses.significant <- all.doses.significant[!is.na(all.doses.significant$gene),]
