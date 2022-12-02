# setwd("Z:/SJShare/SJCOMMON/ECC/Ness Research Team/CCSS/CCSS 20200205/original/fu2007")

# setwd("Z:/SJShare/SJCOMMON/ECC/Ness Research Team/CCSS/CCSS 20200205/expansion/baseline")

# setwd("Z:/SJShare/SJCOMMON/ECC/Ness Research Team/CCSS/CCSS OCT2022/original/fu3")

# setwd("Z:/SJShare/SJCOMMON/ECC/Ness Research Team/CCSS/CCSS OCT2022/original/fu2007")

setwd("Z:/SJShare/SJCOMMON/ECC/Ness Research Team/CCSS/CCSS Dosimetry Data")

setwd("Z:/SJShare/SJCOMMON/ECC/Ness Research Team/CCSS/CCSS OCT2022/original/baseline")

allfiles <- list.files(pattern = "sas7bdat")



variables_required <- c("ccssid", "gender", "dob", "agedx", "diaggrp", "agelstcontact",
                        "heavydrink", "riskydrink", "BMIadj",
                        "maxchestrtdose", "maxneckrtdose", "maxpelvisrtdose", "maxabdrtdose", "maxsegrtdose", 
                        "anthra_jco_dose_5", "aa_class_dose_5", "epitxn_dose_5", "cisplateq_dose_5", "aa_hvymtl_dose_5",
                        "HEI2015_TOTAL_SCORE", 
                        "vpa10", "nopa", "pa20", "vpadays", "vpamin", "mpadays", "mpamin", "mpa10",
                        "evsm", "smnow", "smnvr", "cigmo", "cigd", "smyr",
                        "FRUITSRV", "NUTSFREQ", "VEGSRV", "WGRAINS", "NOTFRIEDFISHFREQ", "DAIRYSRV", "MIXEDBEEFPORKFREQ", "DT_TFAT", "SOFTDRINKSFREQ", "DT_SODI",
                        "nohat", "mercapto_dose_5")

found.variables <- {}
for( i in 1: length(allfiles)){
SASfile <- read_sas(allfiles[i])
found.variables.tmp <- variables_required[variables_required %in% colnames(SASfile)]
found.variables <- c(found.variables, found.variables.tmp)
found.variables <- unique(found.variables)
print(paste0("Found ", length(found.variables), " variables: ", found.variables, " in iteration ", i))
}


allcols <- {}
for( i in 1: length(allfiles)){
  SASfile <- read_sas(allfiles[i])
  allcols.tmp <- colnames(SASfile)
  allcols <- c(allcols, allcols.tmp)
  allcols <- unique(allcols)
  print(i)
}





