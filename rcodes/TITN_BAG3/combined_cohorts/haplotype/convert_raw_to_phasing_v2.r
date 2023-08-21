setwd("Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3")
## ped <- read.table("haplotype_input.ped")
## map <- read.table("haplotype_input.map")
raw <- read.table("haplotype_input_0.8.raw", header = T) # r2 0.8

rownames(raw) <- raw$IID
raw <- raw[-c(1:6)]
HEADER = sub(pattern="_[T,A,G,C,*]+",replacement="",colnames(raw))
colnames(raw) <- gsub("\\.", ":", HEADER)
bim <- read.table("ttn_significant_ld.bim")

bim$V2 <- paste0("chr",bim$V1, ":", bim$V4, ":", bim$V6, ":", bim$V5)
# sum(colnames(raw) %in% bim$V2)
# 7

# map$REF <- bim$V6[match(map$V2, bim$V2)]
# map$ALT <- bim$V5[match(map$V2, bim$V2)]

# > map$REF
# [1] "C" "C" "T" "G" "C" "G" "G" "T" "T" "T"

# > map$ALT
# [1] "G" "T" "C" "A" "T" "A" "A" "C" "C" "C"

first <- nrow(raw)
second <- ncol(raw)
third <- c("P", sapply(strsplit(colnames(raw), split=":"), "[", 2))
fourth <- paste0(rep("S", second), collapse = "")

save.raw <- raw
raw[raw == 1] <- "hetero"
raw[raw == 0] <- "homo_ref"
raw[raw == 2] <- "homo_alt"

raw[raw == "hetero"] <- "0 1"
raw[raw == "homo_ref"] <- "0 0"
raw[raw == "homo_alt"] <- "1 1"
raw[is.na(raw)] <- "? ?"

write(first, file="haplotype_input_edited.txt", append = T)
write(second, file="haplotype_input_edited.txt", append = T)
write(paste0(third, collapse = " "), file="haplotype_input_edited.txt", append = T)
write(fourth, file="haplotype_input_edited.txt", append = T)


for (i in 1:nrow(raw)){
  print(paste0("Doing row: ", i))  
  sample <- paste0("#", rownames(raw)[i])
  write(sample, file="haplotype_input_edited.txt", append = T)
  raw.get <- gsub(" ", "", paste0(as.character(raw[i,]), collapse = ""))
  line1 <- paste(strsplit(raw.get, "")[[1]][seq(1, nchar(raw.get), by = 2)], collapse = " ")
  line2 <- paste(strsplit(raw.get, "")[[1]][seq(2, nchar(raw.get), by = 2)], collapse = " ")
  write(line1, file="haplotype_input_edited.txt", append = T)
  write(line2, file="haplotype_input_edited.txt", append = T)
}
