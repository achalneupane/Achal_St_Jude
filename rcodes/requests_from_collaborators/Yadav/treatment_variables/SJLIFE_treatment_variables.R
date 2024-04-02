library(dplyr)
library(tidyr)
library(haven)

diag <- read_sas("Z:/SJShare/SJCOMMON/ECC/SJLife/SJLIFE Data Freeze/2 Final Data SJLIFE/20200430/Clinical Data/diagnosis.sas7bdat")
diag <- subset(diag, primdx==1 & Newstatus==3)

chemo <- read_sas("Z:/SJShare/SJCOMMON/ECC/SJLife/SJLIFE Data Freeze/2 Final Data SJLIFE/20200430/Clinical Data/chemosum_dose.sas7bdat")
radiation.2 <- read_sas("Z:/SJShare/SJCOMMON/ECC/SJLife/SJLIFE Data Freeze/2 Final Data SJLIFE/20200430/Clinical Data/radiationsum_yn.sas7bdat")

table(radiation.2$anyrt_5)
chemo$anyRT_5 <- radiation.2$anyrt_5[match(chemo$sjlid, radiation.2$sjlid)]

## those who completed campus visit
chemo <- chemo[chemo$Newstatus == 3,]
table(chemo$sjlid %in% diag$sjlid)
# TRUE 
# 5229 
###################
## Without AnyRT ##
###################



cyclophosphamide_dose_5.Yes <- length((which(chemo$cyclophosphamide_dose_5 > 0 & chemo$anyRT_5 == 0)))
cyclophosphamide_dose_5.No <- length((which(chemo$cyclophosphamide_dose_5 == 0 & chemo$anyRT_5 == 0)))

busulfan_dose_5.Yes <- length((which(chemo$busulfan_dose_5 > 0 & chemo$anyRT_5 == 0)))
busulfan_dose_5.No <- length((which(chemo$busulfan_dose_5 == 0 & chemo$anyRT_5 == 0)))

carmustine_dose_5.Yes <- length((which(chemo$carmustine_dose_5 > 0 & chemo$anyRT_5 == 0)))
carmustine_dose_5.No <- length((which(chemo$carmustine_dose_5 == 0 & chemo$anyRT_5 == 0)))

ifosfamide_dose_5.Yes <- length((which(chemo$ifosfamide_dose_5 > 0 & chemo$anyRT_5 == 0)))
ifosfamide_dose_5.No <- length((which(chemo$ifosfamide_dose_5 == 0 & chemo$anyRT_5 == 0)))

lomustine_dose_5.Yes <- length((which(chemo$lomustine_dose_5 > 0 & chemo$anyRT_5 == 0)))
lomustine_dose_5.No <- length((which(chemo$lomustine_dose_5 == 0 & chemo$anyRT_5 == 0)))

mechlorethamine_dose_5.Yes <- length((which(chemo$mechlorethamine_dose_5 > 0 & chemo$anyRT_5 == 0)))
mechlorethamine_dose_5.No <- length((which(chemo$mechlorethamine_dose_5 == 0 & chemo$anyRT_5 == 0)))

melphalan_dose_5.Yes <- length((which(chemo$melphalan_dose_5 > 0 & chemo$anyRT_5 == 0)))
melphalan_dose_5.No <- length((which(chemo$melphalan_dose_5 == 0 & chemo$anyRT_5 == 0)))

procarbazine_dose_5.Yes <- length((which(chemo$procarbazine_dose_5 > 0 & chemo$anyRT_5 == 0)))
procarbazine_dose_5.No <- length((which(chemo$procarbazine_dose_5 == 0 & chemo$anyRT_5 == 0)))

thiotepa_dose_5.Yes <- length((which(chemo$thiotepa_dose_5 > 0 & chemo$anyRT_5 == 0)))
thiotepa_dose_5.No <- length((which(chemo$thiotepa_dose_5 == 0 & chemo$anyRT_5 == 0)))

alkylating_dose_5.Yes <- length((which(chemo$alkylating_dose_5 > 0 & chemo$anyRT_5 == 0)))
alkylating_dose_5.No <- length((which(chemo$alkylating_dose_5 == 0 & chemo$anyRT_5 == 0)))

doxorubicin_dose_5.Yes <- length((which(chemo$doxorubicin_dose_5 > 0 & chemo$anyRT_5 == 0)))
doxorubicin_dose_5.No <- length((which(chemo$doxorubicin_dose_5 == 0 & chemo$anyRT_5 == 0)))

daunorubicin_dose_5.Yes <- length((which(chemo$daunorubicin_dose_5 > 0 & chemo$anyRT_5 == 0)))
daunorubicin_dose_5.No <- length((which(chemo$daunorubicin_dose_5 == 0 & chemo$anyRT_5 == 0)))

epirubicin_dose_5.Yes <- length((which(chemo$epirubicin_dose_5 > 0 & chemo$anyRT_5 == 0)))
epirubicin_dose_5.No <- length((which(chemo$epirubicin_dose_5 == 0 & chemo$anyRT_5 == 0)))

idarubicin_dose_5.Yes <- length((which(chemo$idarubicin_dose_5 > 0 & chemo$anyRT_5 == 0)))
idarubicin_dose_5.No <- length((which(chemo$idarubicin_dose_5 == 0 & chemo$anyRT_5 == 0)))

mitoxantrone_dose_5.Yes <- length((which(chemo$mitoxantrone_dose_5 > 0 & chemo$anyRT_5 == 0)))
mitoxantrone_dose_5.No <- length((which(chemo$mitoxantrone_dose_5 == 0 & chemo$anyRT_5 == 0)))

anthracyclines_dose_5.Yes <- length((which(chemo$anthracyclines_dose_5 > 0 & chemo$anyRT_5 == 0)))
anthracyclines_dose_5.No <- length((which(chemo$anthracyclines_dose_5 == 0 & chemo$anyRT_5 == 0)))

cytarabine_dose_5.Yes <- length((which(chemo$cytarabine_dose_5 > 0 & chemo$anyRT_5 == 0)))
cytarabine_dose_5.No <- length((which(chemo$cytarabine_dose_5 == 0 & chemo$anyRT_5 == 0)))

cytarabine_HD_dose_5.Yes <- length((which(chemo$cytarabine_HD_dose_5 > 0 & chemo$anyRT_5 == 0)))
cytarabine_HD_dose_5.No <- length((which(chemo$cytarabine_HD_dose_5 == 0 & chemo$anyRT_5 == 0)))

mercaptopurine_dose_5.Yes <- length((which(chemo$mercaptopurine_dose_5 > 0 & chemo$anyRT_5 == 0)))
mercaptopurine_dose_5.No <- length((which(chemo$mercaptopurine_dose_5 == 0 & chemo$anyRT_5 == 0)))

methotrexate_im_dose_5.Yes <- length((which(chemo$methotrexate_im_dose_5 > 0 & chemo$anyRT_5 == 0)))
methotrexate_im_dose_5.No <- length((which(chemo$methotrexate_im_dose_5 == 0 & chemo$anyRT_5 == 0)))

methotrexate_it_io_dose_5.Yes <- length((which(chemo$methotrexate_it_io_dose_5 > 0 & chemo$anyRT_5 == 0)))
methotrexate_it_io_dose_5.No <- length((which(chemo$methotrexate_it_io_dose_5 == 0 & chemo$anyRT_5 == 0)))

methotrexate_HD_dose_5.Yes <- length((which(chemo$methotrexate_HD_dose_5 > 0 & chemo$anyRT_5 == 0)))
methotrexate_HD_dose_5.No <- length((which(chemo$methotrexate_HD_dose_5 == 0 & chemo$anyRT_5 == 0)))

thioguanine_dose_5.Yes <- length((which(chemo$thioguanine_dose_5 > 0 & chemo$anyRT_5 == 0)))
thioguanine_dose_5.No <- length((which(chemo$thioguanine_dose_5 == 0 & chemo$anyRT_5 == 0)))

dactinomycin_dose_5.Yes <- length((which(chemo$dactinomycin_dose_5 > 0 & chemo$anyRT_5 == 0)))
dactinomycin_dose_5.No <- length((which(chemo$dactinomycin_dose_5 == 0 & chemo$anyRT_5 == 0)))

bleomycin_dose_5.Yes <- length((which(chemo$bleomycin_dose_5 > 0 & chemo$anyRT_5 == 0)))
bleomycin_dose_5.No <- length((which(chemo$bleomycin_dose_5 == 0 & chemo$anyRT_5 == 0)))

e_asparaginase_dose_5.Yes <- length((which(chemo$e_asparaginase_dose_5 > 0 & chemo$anyRT_5 == 0)))
e_asparaginase_dose_5.No <- length((which(chemo$e_asparaginase_dose_5 == 0 & chemo$anyRT_5 == 0)))

l_asparaginase_dose_5.Yes <- length((which(chemo$l_asparaginase_dose_5 > 0 & chemo$anyRT_5 == 0)))
l_asparaginase_dose_5.No <- length((which(chemo$l_asparaginase_dose_5 == 0 & chemo$anyRT_5 == 0)))

peg_asparaginase_dose_5.Yes <- length((which(chemo$peg_asparaginase_dose_5 > 0 & chemo$anyRT_5 == 0)))
peg_asparaginase_dose_5.No <- length((which(chemo$peg_asparaginase_dose_5 == 0 & chemo$anyRT_5 == 0)))

enzyme_dose_5.Yes <- length((which(chemo$enzyme_dose_5 > 0 & chemo$anyRT_5 == 0)))
enzyme_dose_5.No <- length((which(chemo$enzyme_dose_5 == 0 & chemo$anyRT_5 == 0)))

cortico_dose_5.Yes <- length((which(chemo$cortico_dose_5 > 0 & chemo$anyRT_5 == 0)))
cortico_dose_5.No <- length((which(chemo$cortico_dose_5 == 0 & chemo$anyRT_5 == 0)))

vincristine_dose_5.Yes <- length((which(chemo$vincristine_dose_5 > 0 & chemo$anyRT_5 == 0)))
vincristine_dose_5.No <- length((which(chemo$vincristine_dose_5 == 0 & chemo$anyRT_5 == 0)))

vinblastine_dose_5.Yes <- length((which(chemo$vinblastine_dose_5 > 0 & chemo$anyRT_5 == 0)))
vinblastine_dose_5.No <- length((which(chemo$vinblastine_dose_5 == 0 & chemo$anyRT_5 == 0)))

vinorelbine_dose_5.Yes <- length((which(chemo$vinorelbine_dose_5 > 0 & chemo$anyRT_5 == 0)))
vinorelbine_dose_5.No <- length((which(chemo$vinorelbine_dose_5 == 0 & chemo$anyRT_5 == 0)))

vinca_dose_5.Yes <- length((which(chemo$vinca_dose_5 > 0 & chemo$anyRT_5 == 0)))
vinca_dose_5.No <- length((which(chemo$vinca_dose_5 == 0 & chemo$anyRT_5 == 0)))

etoposide_dose_5.Yes <- length((which(chemo$etoposide_dose_5 > 0 & chemo$anyRT_5 == 0)))
etoposide_dose_5.No <- length((which(chemo$etoposide_dose_5 == 0 & chemo$anyRT_5 == 0)))

teniposide_dose_5.Yes <- length((which(chemo$teniposide_dose_5 > 0 & chemo$anyRT_5 == 0)))
teniposide_dose_5.No <- length((which(chemo$teniposide_dose_5 == 0 & chemo$anyRT_5 == 0)))

cisplatin_dose_5.Yes <- length((which(chemo$cisplatin_dose_5 > 0 & chemo$anyRT_5 == 0)))
cisplatin_dose_5.No <- length((which(chemo$cisplatin_dose_5 == 0 & chemo$anyRT_5 == 0)))

carboplatin_dose_5.Yes <- length((which(chemo$carboplatin_dose_5 > 0 & chemo$anyRT_5 == 0)))
carboplatin_dose_5.No <- length((which(chemo$carboplatin_dose_5 == 0 & chemo$anyRT_5 == 0)))

platinum_dose_5.Yes <- length((which(chemo$platinum_dose_5 > 0 & chemo$anyRT_5 == 0)))
platinum_dose_5.No <- length((which(chemo$platinum_dose_5 == 0 & chemo$anyRT_5 == 0)))

oxaliplatin_dose_5.Yes <- length((which(chemo$oxaliplatin_dose_5 > 0 & chemo$anyRT_5 == 0)))
oxaliplatin_dose_5.No <- length((which(chemo$oxaliplatin_dose_5 == 0 & chemo$anyRT_5 == 0)))

