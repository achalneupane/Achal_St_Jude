install.packages("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/tools/easyQC/EasyQC_23.8.tar.gz")
library(EasyQC)
setwd("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/tools/easyQC/EasyQC_23.8/EasyQC/inst/extdata/")

# decompress text files
system("gzip -d Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/tools/easyQC/EasyQC_23.8/EasyQC/inst/extdata/studyX_file1.txt.gz")
system("gzip -d Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/tools/easyQC/EasyQC_23.8/EasyQC/inst/extdata/studyX_file2.txt.gz")

EasyQC("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/tools/easyQC/EasyQC_23.8/EasyQC/inst/extdata/example_qc.ecf")
