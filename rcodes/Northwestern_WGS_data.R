df <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/CAB/common/WGS_Northwestern/mapping_stats.txt", sep = "\t", header = T)
dim(df)

table(df$STATUS)

coverage <- df[grepl("COVERAGE|SAMPLE", colnames(df), ignore.case = T)]
colnames(coverage)
coverage <- coverage[c("SAMPLE", "COVERAGE_MEAN", "COVERAGE_EXON_10X_perc", 
           "COVERAGE_EXON_20X_perc", "COVERAGE_EXON_30X_perc", 
           "COVERAGE_OVERALL_10X_perc", "COVERAGE_OVERALL_20X_perc", 
           "COVERAGE_OVERALL_30X_perc")]

library(reshape2)
library(ggplot2)
d <- melt(coverage, id.vars="SAMPLE")

ggplot(d, aes(SAMPLE,value, col=variable)) + 
  stat_smooth(method="loess", span=0.1, se=TRUE, aes(fill=variable), alpha=0.3) +
  theme_bw()

d$level <- factor(d$variable)
d$SAMPLE <- factor(d$SAMPLE)
ggplot(d, aes(x=SAMPLE,y=value, colour=variable)) +
  geom_line(aes(linetype=level)) +
  theme_bw()
