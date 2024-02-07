intensities <- read.delim("Z:/ResearchHome/Groups/sapkogrp/projects/CAB/common/metabolomics_Northwestern/Report-Bur-Dis-20231206-CHMP_edited.txt", header = T, sep = "\t", stringsAsFactors = F)
colnames(intensities) <- gsub("^0","",gsub("Bur.Dis.Gra.20231206.", "", colnames(metabol)))
colnames(intensities) <-
c("Compounds", "KEGG.ID", "1537.0","1537.1","1537.3","7597.0","7597.1","7597.3","242.0","242.1","242.3","2775.0","2775.1","2775.3",
  "5170.0","5170.1","5170.3","9152.0","9152.1","9152.3","6124.0","6124.1","6124.3","4682.0","4682.1","4682.3","60955.0","60955.1","60955.3")

# norm.factor <- read.delim("Z:/ResearchHome/Groups/sapkogrp/projects/CAB/common/metabolomics_Northwestern/normalization_factor.txt", header = T, sep = "\t", stringsAsFactors = F)

## Log2 transformation, normalization and scaling

rownames(intensities) <- intensities$Compounds
intensities <- intensities[-c(1:2)]

intensities[] <- lapply(intensities, function(x) as.numeric(gsub(",", "", x)))


library(preprocessCore)
# intensities_log2 = log(intensities, base = 2)
intensities_log2 = log(intensities + 1)  # Base-2 logarithm with added 1 to handle log(0)
intensities_t <- t(intensities_log2)
intensities_t = scale(intensities_t, center = TRUE, scale = TRUE)
