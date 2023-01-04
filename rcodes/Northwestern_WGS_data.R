df <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/CAB/common/WGS_Northwestern/mapping_stats.txt", sep = "\t", header = T)
dim(df)

table(df$STATUS)
# PASS WARNING 
# 91      19 

coverage <- df[grepl("COVERAGE|SAMPLE|DUPLICATION_perc|MEAN_MAPPING_QUALITY|MAPPED_perc|READS_RAW", colnames(df), ignore.case = T)]
colnames(coverage)
# coverage <- coverage[c("SAMPLE", "COVERAGE_MEAN", "COVERAGE_EXON_10X_perc", 
#            "COVERAGE_EXON_20X_perc", "COVERAGE_EXON_30X_perc", 
#            "COVERAGE_OVERALL_10X_perc", "COVERAGE_OVERALL_20X_perc", 
#            "COVERAGE_OVERALL_30X_perc")]
coverage <- coverage[c("SAMPLE", "READS_RAW")]


library(reshape2)
library(ggplot2)
d <- melt(coverage, id.vars="SAMPLE")

ggplot(data=d,
       aes(x=SAMPLE, y=value, group = 1, colour=variable)) +
  geom_point() +
  scale_colour_manual(values=c("black")) +
  geom_line(size = .5) +
  scale_y_log10() +
  # scale_y_continuous(breaks = seq(0, 100, len = 10)) +
  # scale_y_continuous(breaks = c( 20, 30, 50, 80, 90, 100), limits = c(1,100)) +
  # facet_wrap(~variable, ncol=1)
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(size = 12)
        ) 

# panel.background = element_blank()
