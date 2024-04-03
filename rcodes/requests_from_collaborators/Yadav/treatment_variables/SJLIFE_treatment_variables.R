rm(list=ls())
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

cytarabine_OR_cytarabine_HD_dose_5.Yes <- length(which((chemo$cytarabine_dose_5 > 0 & chemo$anyRT_5 == 0)|
      (chemo$cytarabine_HD_dose_5 > 0 & chemo$anyRT_5 == 0)))
cytarabine_OR_cytarabine_HD_dose_5.No <- length(which((chemo$cytarabine_dose_5 == 0 &
      chemo$cytarabine_HD_dose_5 == 0) & chemo$anyRT_5 == 0))

mercaptopurine_dose_5.Yes <- length((which(chemo$mercaptopurine_dose_5 > 0 & chemo$anyRT_5 == 0)))
mercaptopurine_dose_5.No <- length((which(chemo$mercaptopurine_dose_5 == 0 & chemo$anyRT_5 == 0)))

methotrexate_im_dose_5.Yes <- length((which(chemo$methotrexate_im_dose_5 > 0 & chemo$anyRT_5 == 0)))
methotrexate_im_dose_5.No <- length((which(chemo$methotrexate_im_dose_5 == 0 & chemo$anyRT_5 == 0)))

methotrexate_it_io_dose_5.Yes <- length((which(chemo$methotrexate_it_io_dose_5 > 0 & chemo$anyRT_5 == 0)))
methotrexate_it_io_dose_5.No <- length((which(chemo$methotrexate_it_io_dose_5 == 0 & chemo$anyRT_5 == 0)))

methotrexate_HD_dose_5.Yes <- length((which(chemo$methotrexate_HD_dose_5 > 0 & chemo$anyRT_5 == 0)))
methotrexate_HD_dose_5.No <- length((which(chemo$methotrexate_HD_dose_5 == 0 & chemo$anyRT_5 == 0)))

methotrexate_IM_IT_5.Yes <- length((which((chemo$methotrexate_im_dose_5 > 0|chemo$methotrexate_it_io_dose_5 > 0) & chemo$anyRT_5 == 0)))
methotrexate_IM_IT_5.No <- length((which((chemo$methotrexate_im_dose_5 == 0 & chemo$methotrexate_it_io_dose_5 == 0) & chemo$anyRT_5 == 0)))

methotrexate_IM_HD_5.Yes <- length((which((chemo$methotrexate_im_dose_5 > 0|chemo$methotrexate_HD_dose_5 > 0) & chemo$anyRT_5 == 0)))
methotrexate_IM_HD_5.No <- length((which((chemo$methotrexate_im_dose_5 == 0 & chemo$methotrexate_HD_dose_5 == 0) & chemo$anyRT_5 == 0)))

methotrexate_IT_HD_5.Yes <- length((which((chemo$methotrexate_it_io_dose_5 > 0|chemo$methotrexate_HD_dose_5 > 0) & chemo$anyRT_5 == 0)))
methotrexate_IT_HD_5.No <- length((which((chemo$methotrexate_it_io_dose_5 == 0 & chemo$methotrexate_HD_dose_5 == 0) & chemo$anyRT_5 == 0)))

methotrexate_IM_IT_HD_5.Yes <- length((which((chemo$methotrexate_im_dose_5 > 0|chemo$methotrexate_it_io_dose_5 > 0|chemo$methotrexate_HD_dose_5 > 0) & chemo$anyRT_5 == 0)))
methotrexate_IM_IT_HD_5.No <- length((which((chemo$methotrexate_im_dose_5 == 0 & chemo$methotrexate_it_io_dose_5 == 0 & chemo$methotrexate_HD_dose_5 == 0) & chemo$anyRT_5 == 0)))

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

etoposide_teniposide_dose_5.Yes <- length((which((chemo$etoposide_dose_5 > 0|chemo$teniposide_dose_5 > 0) & chemo$anyRT_5 == 0)))
etoposide_teniposide_dose_5.No <- length((which((chemo$etoposide_dose_5 == 0|chemo$teniposide_dose_5 == 0) & chemo$anyRT_5 == 0)))

cisplatin_dose_5.Yes <- length((which(chemo$cisplatin_dose_5 > 0 & chemo$anyRT_5 == 0)))
cisplatin_dose_5.No <- length((which(chemo$cisplatin_dose_5 == 0 & chemo$anyRT_5 == 0)))

carboplatin_dose_5.Yes <- length((which(chemo$carboplatin_dose_5 > 0 & chemo$anyRT_5 == 0)))
carboplatin_dose_5.No <- length((which(chemo$carboplatin_dose_5 == 0 & chemo$anyRT_5 == 0)))

platinum_dose_5.Yes <- length((which(chemo$platinum_dose_5 > 0 & chemo$anyRT_5 == 0)))
platinum_dose_5.No <- length((which(chemo$platinum_dose_5 == 0 & chemo$anyRT_5 == 0)))

oxaliplatin_dose_5.Yes <- length((which(chemo$oxaliplatin_dose_5 > 0 & chemo$anyRT_5 == 0)))
oxaliplatin_dose_5.No <- length((which(chemo$oxaliplatin_dose_5 == 0 & chemo$anyRT_5 == 0)))


withoutRT.Yes <- as.data.frame(t(cbind.data.frame(cyclophosphamide_dose_5=cyclophosphamide_dose_5.Yes, 
                   busulfan_dose_5=busulfan_dose_5.Yes, 
                   carmustine_dose_5=carmustine_dose_5.Yes,
                   ifosfamide_dose_5=ifosfamide_dose_5.Yes,
                   lomustine_dose_5=lomustine_dose_5.Yes,
                   mechlorethamine_dose_5=mechlorethamine_dose_5.Yes,
                   melphalan_dose_5=melphalan_dose_5.Yes,
                   procarbazine_dose_5=procarbazine_dose_5.Yes,
                   thiotepa_dose_5=thiotepa_dose_5.Yes,
                   alkylating_dose_5=alkylating_dose_5.Yes,
                   doxorubicin_dose_5=doxorubicin_dose_5.Yes,
                   daunorubicin_dose_5=daunorubicin_dose_5.Yes,
                   epirubicin_dose_5=epirubicin_dose_5.Yes,
                   idarubicin_dose_5=idarubicin_dose_5.Yes,
                   mitoxantrone_dose_5=mitoxantrone_dose_5.Yes,
                   anthracyclines_dose_5=anthracyclines_dose_5.Yes,
                   cytarabine_dose_5=cytarabine_dose_5.Yes,
                   cytarabine_HD_dose_5=cytarabine_HD_dose_5.Yes,
                   cytarabine_OR_cytarabine_HD_dose_5=cytarabine_OR_cytarabine_HD_dose_5.Yes,
                   mercaptopurine_dose_5=mercaptopurine_dose_5.Yes,
                   methotrexate_im_dose_5=methotrexate_im_dose_5.Yes,
                   methotrexate_it_io_dose_5=methotrexate_it_io_dose_5.Yes,
                   methotrexate_HD_dose_5=methotrexate_HD_dose_5.Yes,
                   methotrexate_IM_IT_5=methotrexate_IM_IT_5.Yes,
                   methotrexate_IM_HD_5=methotrexate_IM_HD_5.Yes,
                   methotrexate_IT_HD_5=methotrexate_IT_HD_5.Yes,
                   methotrexate_IM_IT_HD_5=methotrexate_IM_IT_HD_5.Yes,
                   thioguanine_dose_5=thioguanine_dose_5.Yes,
                   dactinomycin_dose_5=dactinomycin_dose_5.Yes,
                   bleomycin_dose_5=bleomycin_dose_5.Yes,
                   e_asparaginase_dose_5=e_asparaginase_dose_5.Yes,
                   l_asparaginase_dose_5=l_asparaginase_dose_5.Yes,
                   peg_asparaginase_dose_5=peg_asparaginase_dose_5.Yes,
                   enzyme_dose_5=enzyme_dose_5.Yes,
                   cortico_dose_5=cortico_dose_5.Yes,
                   vincristine_dose_5=vincristine_dose_5.Yes,
                   vinblastine_dose_5=vinblastine_dose_5.Yes,
                   vinorelbine_dose_5=vinorelbine_dose_5.Yes,
                   vinca_dose_5=vinca_dose_5.Yes,
                   etoposide_dose_5=etoposide_dose_5.Yes,
                   teniposide_dose_5=teniposide_dose_5.Yes,
                   etoposide_teniposide_dose_5=etoposide_teniposide_dose_5.Yes,
                   cisplatin_dose_5=cisplatin_dose_5.Yes,
                   carboplatin_dose_5=carboplatin_dose_5.Yes,
                   platinum_dose_5=platinum_dose_5.Yes,
                   oxaliplatin_dose_5=oxaliplatin_dose_5.Yes)))


withoutRT.No <- as.data.frame(t(cbind.data.frame(cyclophosphamide_dose_5=cyclophosphamide_dose_5.No, 
                                                 busulfan_dose_5=busulfan_dose_5.No, 
                                                 carmustine_dose_5=carmustine_dose_5.No,
                                                 ifosfamide_dose_5=ifosfamide_dose_5.No,
                                                 lomustine_dose_5=lomustine_dose_5.No,
                                                 mechlorethamine_dose_5=mechlorethamine_dose_5.No,
                                                 melphalan_dose_5=melphalan_dose_5.No,
                                                 procarbazine_dose_5=procarbazine_dose_5.No,
                                                 thiotepa_dose_5=thiotepa_dose_5.No,
                                                 alkylating_dose_5=alkylating_dose_5.No,
                                                 doxorubicin_dose_5=doxorubicin_dose_5.No,
                                                 daunorubicin_dose_5=daunorubicin_dose_5.No,
                                                 epirubicin_dose_5=epirubicin_dose_5.No,
                                                 idarubicin_dose_5=idarubicin_dose_5.No,
                                                 mitoxantrone_dose_5=mitoxantrone_dose_5.No,
                                                 anthracyclines_dose_5=anthracyclines_dose_5.No,
                                                 cytarabine_dose_5=cytarabine_dose_5.No,
                                                 cytarabine_HD_dose_5=cytarabine_HD_dose_5.No,
                                                 cytarabine_OR_cytarabine_HD_dose_5=cytarabine_OR_cytarabine_HD_dose_5.No,
                                                 mercaptopurine_dose_5=mercaptopurine_dose_5.No,
                                                 methotrexate_im_dose_5=methotrexate_im_dose_5.No,
                                                 methotrexate_it_io_dose_5=methotrexate_it_io_dose_5.No,
                                                 methotrexate_HD_dose_5=methotrexate_HD_dose_5.No,
                                                 methotrexate_IM_IT_5=methotrexate_IM_IT_5.No,
                                                 methotrexate_IM_HD_5=methotrexate_IM_HD_5.No,
                                                 methotrexate_IT_HD_5=methotrexate_IT_HD_5.No,
                                                 methotrexate_IM_IT_HD_5=methotrexate_IM_IT_HD_5.No,
                                                 thioguanine_dose_5=thioguanine_dose_5.No,
                                                 dactinomycin_dose_5=dactinomycin_dose_5.No,
                                                 bleomycin_dose_5=bleomycin_dose_5.No,
                                                 e_asparaginase_dose_5=e_asparaginase_dose_5.No,
                                                 l_asparaginase_dose_5=l_asparaginase_dose_5.No,
                                                 peg_asparaginase_dose_5=peg_asparaginase_dose_5.No,
                                                 enzyme_dose_5=enzyme_dose_5.No,
                                                 cortico_dose_5=cortico_dose_5.No,
                                                 vincristine_dose_5=vincristine_dose_5.No,
                                                 vinblastine_dose_5=vinblastine_dose_5.No,
                                                 vinorelbine_dose_5=vinorelbine_dose_5.No,
                                                 vinca_dose_5=vinca_dose_5.No,
                                                 etoposide_dose_5=etoposide_dose_5.No,
                                                 teniposide_dose_5=teniposide_dose_5.No,
                                                 etoposide_teniposide_dose_5=etoposide_teniposide_dose_5.No,
                                                 cisplatin_dose_5=cisplatin_dose_5.No,
                                                 carboplatin_dose_5=carboplatin_dose_5.No,
                                                 platinum_dose_5=platinum_dose_5.No,
                                                 oxaliplatin_dose_5=oxaliplatin_dose_5.No)))

table(rownames(withoutRT.Yes) == rownames(withoutRT.No))
withoutRT <- cbind.data.frame(var=rownames(withoutRT.Yes), YesChemo=withoutRT.Yes$V1, NoChemo=withoutRT.No$V1)

chemotherapy_agents <- c("Cyclophosphamide", "Busulfan", "Carmustine", "Ifosfamide", "Lomustine", "Mechlorethamine", 
                         "Melphalan", "Procarbazine", "Thiotepa", "Cumulative Alkylating Agent", 
                         "Doxorubicin (IV)", "Daunorubicin (IV)", "Epirubicin (IV)", "Idarubicin (IV)", "Mitoxantrone (IV)", 
                         "Cumulative Anthracyline Dose (doxorubicin equivalent per JAMA)",
                         "Cytarabine cumulative dose", "High Dose Cytarabine", "Cytarabine cumulative dose OR High Dose Cytarabine", "6-Mercaptopurine (IV/PO)", 
                         "Methotrexate (IV/PO/IM)", "Methotrexate (IT/IO)", "High Dose Methotrexate (IV)", "Methotrexate (IV/PO/IM) OR Methotrexate (IT/IO)",
                         "Methotrexate (IV/PO/IM) OR High Dose Methotrexate (IV)", "Methotrexate (IT/IO) OR High Dose Methotrexate (IV)",
                         "Methotrexate (IV/PO/IM) OR Methotrexate (IT/IO) OR High Dose Methotrexate (IV)",
                         "Thioguanine (IV/PO)", "Dactinomycin (IV, mv/m2)", 
                         "Bleomycin (IV, IU/m2)", "Erwinia Asparaginase", 
                         "L-Asparaginase", "Peg-Asparaginase", "Cumulative Asparaginase Enzymes (Peg-Asparaginase Equivalent Dose)", 
                         "Cumulative Corticosteroids (Prednisone Equivalent Dose)", 
                         "Vincristine (IV)", "Vinblastine (IV)", "Vinorelbine (IV)", 
                         "Cumulative Vinca Alkaloids", "Etoposide (IV)", "Teniposide (IV)", "Etoposide (IV) OR Teniposide (IV)",
                         "Cisplatin (IV/IA)", "Carboplatin (IV/IO)", "Cumulative Platinum Agent (Cisplatin Equivalent Dose)", 
                         "Oxaliplatin (IV)")

rownames(withoutRT) <- chemotherapy_agents

################
## With AnyRT ##
################

cyclophosphamide_dose_5.Yes <- length((which(chemo$cyclophosphamide_dose_5 > 0 & chemo$anyRT_5 == 1)))
cyclophosphamide_dose_5.No <- length((which(chemo$cyclophosphamide_dose_5 == 0 & chemo$anyRT_5 == 1)))

busulfan_dose_5.Yes <- length((which(chemo$busulfan_dose_5 > 0 & chemo$anyRT_5 == 1)))
busulfan_dose_5.No <- length((which(chemo$busulfan_dose_5 == 0 & chemo$anyRT_5 == 1)))

carmustine_dose_5.Yes <- length((which(chemo$carmustine_dose_5 > 0 & chemo$anyRT_5 == 1)))
carmustine_dose_5.No <- length((which(chemo$carmustine_dose_5 == 0 & chemo$anyRT_5 == 1)))

ifosfamide_dose_5.Yes <- length((which(chemo$ifosfamide_dose_5 > 0 & chemo$anyRT_5 == 1)))
ifosfamide_dose_5.No <- length((which(chemo$ifosfamide_dose_5 == 0 & chemo$anyRT_5 == 1)))

lomustine_dose_5.Yes <- length((which(chemo$lomustine_dose_5 > 0 & chemo$anyRT_5 == 1)))
lomustine_dose_5.No <- length((which(chemo$lomustine_dose_5 == 0 & chemo$anyRT_5 == 1)))

mechlorethamine_dose_5.Yes <- length((which(chemo$mechlorethamine_dose_5 > 0 & chemo$anyRT_5 == 1)))
mechlorethamine_dose_5.No <- length((which(chemo$mechlorethamine_dose_5 == 0 & chemo$anyRT_5 == 1)))

melphalan_dose_5.Yes <- length((which(chemo$melphalan_dose_5 > 0 & chemo$anyRT_5 == 1)))
melphalan_dose_5.No <- length((which(chemo$melphalan_dose_5 == 0 & chemo$anyRT_5 == 1)))

procarbazine_dose_5.Yes <- length((which(chemo$procarbazine_dose_5 > 0 & chemo$anyRT_5 == 1)))
procarbazine_dose_5.No <- length((which(chemo$procarbazine_dose_5 == 0 & chemo$anyRT_5 == 1)))

thiotepa_dose_5.Yes <- length((which(chemo$thiotepa_dose_5 > 0 & chemo$anyRT_5 == 1)))
thiotepa_dose_5.No <- length((which(chemo$thiotepa_dose_5 == 0 & chemo$anyRT_5 == 1)))

alkylating_dose_5.Yes <- length((which(chemo$alkylating_dose_5 > 0 & chemo$anyRT_5 == 1)))
alkylating_dose_5.No <- length((which(chemo$alkylating_dose_5 == 0 & chemo$anyRT_5 == 1)))

doxorubicin_dose_5.Yes <- length((which(chemo$doxorubicin_dose_5 > 0 & chemo$anyRT_5 == 1)))
doxorubicin_dose_5.No <- length((which(chemo$doxorubicin_dose_5 == 0 & chemo$anyRT_5 == 1)))

daunorubicin_dose_5.Yes <- length((which(chemo$daunorubicin_dose_5 > 0 & chemo$anyRT_5 == 1)))
daunorubicin_dose_5.No <- length((which(chemo$daunorubicin_dose_5 == 0 & chemo$anyRT_5 == 1)))

epirubicin_dose_5.Yes <- length((which(chemo$epirubicin_dose_5 > 0 & chemo$anyRT_5 == 1)))
epirubicin_dose_5.No <- length((which(chemo$epirubicin_dose_5 == 0 & chemo$anyRT_5 == 1)))

idarubicin_dose_5.Yes <- length((which(chemo$idarubicin_dose_5 > 0 & chemo$anyRT_5 == 1)))
idarubicin_dose_5.No <- length((which(chemo$idarubicin_dose_5 == 0 & chemo$anyRT_5 == 1)))

mitoxantrone_dose_5.Yes <- length((which(chemo$mitoxantrone_dose_5 > 0 & chemo$anyRT_5 == 1)))
mitoxantrone_dose_5.No <- length((which(chemo$mitoxantrone_dose_5 == 0 & chemo$anyRT_5 == 1)))

anthracyclines_dose_5.Yes <- length((which(chemo$anthracyclines_dose_5 > 0 & chemo$anyRT_5 == 1)))
anthracyclines_dose_5.No <- length((which(chemo$anthracyclines_dose_5 == 0 & chemo$anyRT_5 == 1)))

cytarabine_dose_5.Yes <- length((which(chemo$cytarabine_dose_5 > 0 & chemo$anyRT_5 == 1)))
cytarabine_dose_5.No <- length((which(chemo$cytarabine_dose_5 == 0 & chemo$anyRT_5 == 1)))

cytarabine_HD_dose_5.Yes <- length((which(chemo$cytarabine_HD_dose_5 > 0 & chemo$anyRT_5 == 1)))
cytarabine_HD_dose_5.No <- length((which(chemo$cytarabine_HD_dose_5 == 0 & chemo$anyRT_5 == 1)))

cytarabine_OR_cytarabine_HD_dose_5.Yes <- length(which((chemo$cytarabine_dose_5 > 0 & chemo$anyRT_5 == 1)|
                                                         (chemo$cytarabine_HD_dose_5 > 0 & chemo$anyRT_5 == 1)))
cytarabine_OR_cytarabine_HD_dose_5.No <- length(which((chemo$cytarabine_dose_5 == 0 &
                                                         chemo$cytarabine_HD_dose_5 == 0) & chemo$anyRT_5 == 1))

mercaptopurine_dose_5.Yes <- length((which(chemo$mercaptopurine_dose_5 > 0 & chemo$anyRT_5 == 1)))
mercaptopurine_dose_5.No <- length((which(chemo$mercaptopurine_dose_5 == 0 & chemo$anyRT_5 == 1)))

methotrexate_im_dose_5.Yes <- length((which(chemo$methotrexate_im_dose_5 > 0 & chemo$anyRT_5 == 1)))
methotrexate_im_dose_5.No <- length((which(chemo$methotrexate_im_dose_5 == 0 & chemo$anyRT_5 == 1)))

methotrexate_it_io_dose_5.Yes <- length((which(chemo$methotrexate_it_io_dose_5 > 0 & chemo$anyRT_5 == 1)))
methotrexate_it_io_dose_5.No <- length((which(chemo$methotrexate_it_io_dose_5 == 0 & chemo$anyRT_5 == 1)))

methotrexate_HD_dose_5.Yes <- length((which(chemo$methotrexate_HD_dose_5 > 0 & chemo$anyRT_5 == 1)))
methotrexate_HD_dose_5.No <- length((which(chemo$methotrexate_HD_dose_5 == 0 & chemo$anyRT_5 == 1)))

methotrexate_IM_IT_5.Yes <- length((which((chemo$methotrexate_im_dose_5 > 0|chemo$methotrexate_it_io_dose_5 > 0) & chemo$anyRT_5 == 1)))
methotrexate_IM_IT_5.No <- length((which((chemo$methotrexate_im_dose_5 == 0 & chemo$methotrexate_it_io_dose_5 == 0) & chemo$anyRT_5 == 1)))

methotrexate_IM_HD_5.Yes <- length((which((chemo$methotrexate_im_dose_5 > 0|chemo$methotrexate_HD_dose_5 > 0) & chemo$anyRT_5 == 1)))
methotrexate_IM_HD_5.No <- length((which((chemo$methotrexate_im_dose_5 == 0 & chemo$methotrexate_HD_dose_5 == 0) & chemo$anyRT_5 == 1)))

methotrexate_IT_HD_5.Yes <- length((which((chemo$methotrexate_it_io_dose_5 > 0|chemo$methotrexate_HD_dose_5 > 0) & chemo$anyRT_5 == 1)))
methotrexate_IT_HD_5.No <- length((which((chemo$methotrexate_it_io_dose_5 == 0 & chemo$methotrexate_HD_dose_5 == 0) & chemo$anyRT_5 == 1)))

methotrexate_IM_IT_HD_5.Yes <- length((which((chemo$methotrexate_im_dose_5 > 0|chemo$methotrexate_it_io_dose_5 > 0|chemo$methotrexate_HD_dose_5 > 0) & chemo$anyRT_5 == 1)))
methotrexate_IM_IT_HD_5.No <- length((which((chemo$methotrexate_im_dose_5 == 0 & chemo$methotrexate_it_io_dose_5 == 0 & chemo$methotrexate_HD_dose_5 == 0) & chemo$anyRT_5 == 1)))

thioguanine_dose_5.Yes <- length((which(chemo$thioguanine_dose_5 > 0 & chemo$anyRT_5 == 1)))
thioguanine_dose_5.No <- length((which(chemo$thioguanine_dose_5 == 0 & chemo$anyRT_5 == 1)))

dactinomycin_dose_5.Yes <- length((which(chemo$dactinomycin_dose_5 > 0 & chemo$anyRT_5 == 1)))
dactinomycin_dose_5.No <- length((which(chemo$dactinomycin_dose_5 == 0 & chemo$anyRT_5 == 1)))

bleomycin_dose_5.Yes <- length((which(chemo$bleomycin_dose_5 > 0 & chemo$anyRT_5 == 1)))
bleomycin_dose_5.No <- length((which(chemo$bleomycin_dose_5 == 0 & chemo$anyRT_5 == 1)))

e_asparaginase_dose_5.Yes <- length((which(chemo$e_asparaginase_dose_5 > 0 & chemo$anyRT_5 == 1)))
e_asparaginase_dose_5.No <- length((which(chemo$e_asparaginase_dose_5 == 0 & chemo$anyRT_5 == 1)))

l_asparaginase_dose_5.Yes <- length((which(chemo$l_asparaginase_dose_5 > 0 & chemo$anyRT_5 == 1)))
l_asparaginase_dose_5.No <- length((which(chemo$l_asparaginase_dose_5 == 0 & chemo$anyRT_5 == 1)))

peg_asparaginase_dose_5.Yes <- length((which(chemo$peg_asparaginase_dose_5 > 0 & chemo$anyRT_5 == 1)))
peg_asparaginase_dose_5.No <- length((which(chemo$peg_asparaginase_dose_5 == 0 & chemo$anyRT_5 == 1)))

enzyme_dose_5.Yes <- length((which(chemo$enzyme_dose_5 > 0 & chemo$anyRT_5 == 1)))
enzyme_dose_5.No <- length((which(chemo$enzyme_dose_5 == 0 & chemo$anyRT_5 == 1)))

cortico_dose_5.Yes <- length((which(chemo$cortico_dose_5 > 0 & chemo$anyRT_5 == 1)))
cortico_dose_5.No <- length((which(chemo$cortico_dose_5 == 0 & chemo$anyRT_5 == 1)))

vincristine_dose_5.Yes <- length((which(chemo$vincristine_dose_5 > 0 & chemo$anyRT_5 == 1)))
vincristine_dose_5.No <- length((which(chemo$vincristine_dose_5 == 0 & chemo$anyRT_5 == 1)))

vinblastine_dose_5.Yes <- length((which(chemo$vinblastine_dose_5 > 0 & chemo$anyRT_5 == 1)))
vinblastine_dose_5.No <- length((which(chemo$vinblastine_dose_5 == 0 & chemo$anyRT_5 == 1)))

vinorelbine_dose_5.Yes <- length((which(chemo$vinorelbine_dose_5 > 0 & chemo$anyRT_5 == 1)))
vinorelbine_dose_5.No <- length((which(chemo$vinorelbine_dose_5 == 0 & chemo$anyRT_5 == 1)))

vinca_dose_5.Yes <- length((which(chemo$vinca_dose_5 > 0 & chemo$anyRT_5 == 1)))
vinca_dose_5.No <- length((which(chemo$vinca_dose_5 == 0 & chemo$anyRT_5 == 1)))

etoposide_dose_5.Yes <- length((which(chemo$etoposide_dose_5 > 0 & chemo$anyRT_5 == 1)))
etoposide_dose_5.No <- length((which(chemo$etoposide_dose_5 == 0 & chemo$anyRT_5 == 1)))

teniposide_dose_5.Yes <- length((which(chemo$teniposide_dose_5 > 0 & chemo$anyRT_5 == 1)))
teniposide_dose_5.No <- length((which(chemo$teniposide_dose_5 == 0 & chemo$anyRT_5 == 1)))

etoposide_teniposide_dose_5.Yes <- length((which((chemo$etoposide_dose_5 > 0|chemo$teniposide_dose_5 > 0) & chemo$anyRT_5 == 1)))
etoposide_teniposide_dose_5.No <- length((which((chemo$etoposide_dose_5 == 0|chemo$teniposide_dose_5 == 0) & chemo$anyRT_5 == 1)))

cisplatin_dose_5.Yes <- length((which(chemo$cisplatin_dose_5 > 0 & chemo$anyRT_5 == 1)))
cisplatin_dose_5.No <- length((which(chemo$cisplatin_dose_5 == 0 & chemo$anyRT_5 == 1)))

carboplatin_dose_5.Yes <- length((which(chemo$carboplatin_dose_5 > 0 & chemo$anyRT_5 == 1)))
carboplatin_dose_5.No <- length((which(chemo$carboplatin_dose_5 == 0 & chemo$anyRT_5 == 1)))

platinum_dose_5.Yes <- length((which(chemo$platinum_dose_5 > 0 & chemo$anyRT_5 == 1)))
platinum_dose_5.No <- length((which(chemo$platinum_dose_5 == 0 & chemo$anyRT_5 == 1)))

oxaliplatin_dose_5.Yes <- length((which(chemo$oxaliplatin_dose_5 > 0 & chemo$anyRT_5 == 1)))
oxaliplatin_dose_5.No <- length((which(chemo$oxaliplatin_dose_5 == 0 & chemo$anyRT_5 == 1)))


withRT.Yes <- as.data.frame(t(cbind.data.frame(cyclophosphamide_dose_5=cyclophosphamide_dose_5.Yes, 
                                               busulfan_dose_5=busulfan_dose_5.Yes, 
                                               carmustine_dose_5=carmustine_dose_5.Yes,
                                               ifosfamide_dose_5=ifosfamide_dose_5.Yes,
                                               lomustine_dose_5=lomustine_dose_5.Yes,
                                               mechlorethamine_dose_5=mechlorethamine_dose_5.Yes,
                                               melphalan_dose_5=melphalan_dose_5.Yes,
                                               procarbazine_dose_5=procarbazine_dose_5.Yes,
                                               thiotepa_dose_5=thiotepa_dose_5.Yes,
                                               alkylating_dose_5=alkylating_dose_5.Yes,
                                               doxorubicin_dose_5=doxorubicin_dose_5.Yes,
                                               daunorubicin_dose_5=daunorubicin_dose_5.Yes,
                                               epirubicin_dose_5=epirubicin_dose_5.Yes,
                                               idarubicin_dose_5=idarubicin_dose_5.Yes,
                                               mitoxantrone_dose_5=mitoxantrone_dose_5.Yes,
                                               anthracyclines_dose_5=anthracyclines_dose_5.Yes,
                                               cytarabine_dose_5=cytarabine_dose_5.Yes,
                                               cytarabine_HD_dose_5=cytarabine_HD_dose_5.Yes,
                                               cytarabine_OR_cytarabine_HD_dose_5=cytarabine_OR_cytarabine_HD_dose_5.Yes,
                                               mercaptopurine_dose_5=mercaptopurine_dose_5.Yes,
                                               methotrexate_im_dose_5=methotrexate_im_dose_5.Yes,
                                               methotrexate_it_io_dose_5=methotrexate_it_io_dose_5.Yes,
                                               methotrexate_HD_dose_5=methotrexate_HD_dose_5.Yes,
                                               methotrexate_IM_IT_5=methotrexate_IM_IT_5.Yes,
                                               methotrexate_IM_HD_5=methotrexate_IM_HD_5.Yes,
                                               methotrexate_IT_HD_5=methotrexate_IT_HD_5.Yes,
                                               methotrexate_IM_IT_HD_5=methotrexate_IM_IT_HD_5.Yes,
                                               thioguanine_dose_5=thioguanine_dose_5.Yes,
                                               dactinomycin_dose_5=dactinomycin_dose_5.Yes,
                                               bleomycin_dose_5=bleomycin_dose_5.Yes,
                                               e_asparaginase_dose_5=e_asparaginase_dose_5.Yes,
                                               l_asparaginase_dose_5=l_asparaginase_dose_5.Yes,
                                               peg_asparaginase_dose_5=peg_asparaginase_dose_5.Yes,
                                               enzyme_dose_5=enzyme_dose_5.Yes,
                                               cortico_dose_5=cortico_dose_5.Yes,
                                               vincristine_dose_5=vincristine_dose_5.Yes,
                                               vinblastine_dose_5=vinblastine_dose_5.Yes,
                                               vinorelbine_dose_5=vinorelbine_dose_5.Yes,
                                               vinca_dose_5=vinca_dose_5.Yes,
                                               etoposide_dose_5=etoposide_dose_5.Yes,
                                               teniposide_dose_5=teniposide_dose_5.Yes,
                                               etoposide_teniposide_dose_5=etoposide_teniposide_dose_5.Yes,
                                               cisplatin_dose_5=cisplatin_dose_5.Yes,
                                               carboplatin_dose_5=carboplatin_dose_5.Yes,
                                               platinum_dose_5=platinum_dose_5.Yes,
                                               oxaliplatin_dose_5=oxaliplatin_dose_5.Yes)))


withRT.No <- as.data.frame(t(cbind.data.frame(cyclophosphamide_dose_5=cyclophosphamide_dose_5.No, 
                                              busulfan_dose_5=busulfan_dose_5.No, 
                                              carmustine_dose_5=carmustine_dose_5.No,
                                              ifosfamide_dose_5=ifosfamide_dose_5.No,
                                              lomustine_dose_5=lomustine_dose_5.No,
                                              mechlorethamine_dose_5=mechlorethamine_dose_5.No,
                                              melphalan_dose_5=melphalan_dose_5.No,
                                              procarbazine_dose_5=procarbazine_dose_5.No,
                                              thiotepa_dose_5=thiotepa_dose_5.No,
                                              alkylating_dose_5=alkylating_dose_5.No,
                                              doxorubicin_dose_5=doxorubicin_dose_5.No,
                                              daunorubicin_dose_5=daunorubicin_dose_5.No,
                                              epirubicin_dose_5=epirubicin_dose_5.No,
                                              idarubicin_dose_5=idarubicin_dose_5.No,
                                              mitoxantrone_dose_5=mitoxantrone_dose_5.No,
                                              anthracyclines_dose_5=anthracyclines_dose_5.No,
                                              cytarabine_dose_5=cytarabine_dose_5.No,
                                              cytarabine_HD_dose_5=cytarabine_HD_dose_5.No,
                                              cytarabine_OR_cytarabine_HD_dose_5=cytarabine_OR_cytarabine_HD_dose_5.No,
                                              mercaptopurine_dose_5=mercaptopurine_dose_5.No,
                                              methotrexate_im_dose_5=methotrexate_im_dose_5.No,
                                              methotrexate_it_io_dose_5=methotrexate_it_io_dose_5.No,
                                              methotrexate_HD_dose_5=methotrexate_HD_dose_5.No,
                                              methotrexate_IM_IT_5=methotrexate_IM_IT_5.No,
                                              methotrexate_IM_HD_5=methotrexate_IM_HD_5.No,
                                              methotrexate_IT_HD_5=methotrexate_IT_HD_5.No,
                                              methotrexate_IM_IT_HD_5=methotrexate_IM_IT_HD_5.No,
                                              thioguanine_dose_5=thioguanine_dose_5.No,
                                              dactinomycin_dose_5=dactinomycin_dose_5.No,
                                              bleomycin_dose_5=bleomycin_dose_5.No,
                                              e_asparaginase_dose_5=e_asparaginase_dose_5.No,
                                              l_asparaginase_dose_5=l_asparaginase_dose_5.No,
                                              peg_asparaginase_dose_5=peg_asparaginase_dose_5.No,
                                              enzyme_dose_5=enzyme_dose_5.No,
                                              cortico_dose_5=cortico_dose_5.No,
                                              vincristine_dose_5=vincristine_dose_5.No,
                                              vinblastine_dose_5=vinblastine_dose_5.No,
                                              vinorelbine_dose_5=vinorelbine_dose_5.No,
                                              vinca_dose_5=vinca_dose_5.No,
                                              etoposide_dose_5=etoposide_dose_5.No,
                                              teniposide_dose_5=teniposide_dose_5.No,
                                              etoposide_teniposide_dose_5=etoposide_teniposide_dose_5.No,
                                              cisplatin_dose_5=cisplatin_dose_5.No,
                                              carboplatin_dose_5=carboplatin_dose_5.No,
                                              platinum_dose_5=platinum_dose_5.No,
                                              oxaliplatin_dose_5=oxaliplatin_dose_5.No)))

table(rownames(withRT.Yes) == rownames(withRT.No))
withRT <- cbind.data.frame(var=rownames(withRT.Yes), YesChemo=withRT.Yes$V1, NoChemo=withRT.No$V1)


rownames(withRT) <- chemotherapy_agents



