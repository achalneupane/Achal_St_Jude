setwd('Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/Dyslipidemia')



# list.files()
PGS000888 <- read.table('PGS000888.txt.gz_cleaned.txt', header = TRUE, sep = "\t")
PGS000888$study <- "PGS000888"
PGS000677 = read.delim('new_PGS000677.txt.gz_cleaned.txt', header = TRUE, sep = "\t")
PGS000677$study <- "PGS000677"
PGS000686 = read.delim('new_PGS000686.txt.gz_cleaned.txt', header = TRUE, sep = "\t")
PGS000686$study <- "PGS000686"
PGS000699 = read.delim('new_PGS000699.txt.gz_cleaned.txt', header = TRUE, sep = "\t")
PGS000699$study <- "PGS000699"
PGS000688 = read.delim('new_PGS000688.txt.gz_cleaned.txt', header = TRUE, sep = "\t")
PGS000688$study <- "PGS000688"
df = rbind.data.frame(PGS000888, PGS000677, PGS000686, PGS000699, PGS000688)

df.saved <- df



df[df$effect_weight < 0, c("other_allele", "effect_allele")] <- df[df$effect_weight < 0, c("effect_allele", "other_allele")]
df$effect_weight <- abs(df$effect_weight)

head(df)

df$KEY <- paste0(df$chr_name, ":", df$chr_position)
df$chr <- paste0("chr", df$chr_name)
df$start <- df$chr_position -1
df$end <- df$chr_position

df$start <- as.character(df$start)
df$end <- as.character(df$end)

GRCh37 <- df[c("chr", "start", "end", "KEY", "study")]
# write.table(GRCh37, "GRCh37_PRS_score_dyslipidemia.txt", sep = "\t", quote = F, col.names = F, row.names = F)

GRCh38 <- read.table("GRCh38_PRS_score_dyslipidemia.bed")
df$start <- GRCh38$V2[match(df$KEY, GRCh38$V4)]
df$end <- GRCh38$V3[match(df$KEY, GRCh38$V4)]
df$chr_position <- df$end 

prs.df <- df[c("chr_name", "chr_position", "effect_allele", "other_allele", "effect_weight", "study")]
write.table(prs.df, "GRCh38_PRS_score_dyslipidemia.txt", sep = "\t", quote = F, col.names = F, row.names = F)

bim <- read.table("dyslipidemia_found_in_all.bim", header = F)
bim$KEY <- paste0(bim$V1, ":", bim$V4)
prs.df$KEY <- paste(prs.df$chr_name, prs.df$chr_position, sep = ":")
sum(prs.df$KEY %in% bim$KEY)
prs.df$KEY[!prs.df$KEY %in% bim$KEY]


