setwd("Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3")
raw <- read.table("haplotype_input_0.8.raw", header = T)
rownames(raw) <- raw$IID
raw <- raw[-c(1:6)]
HEADER = sub(pattern="_[T,A,G,C,*]+",replacement="",colnames(raw))
colnames(raw) <- gsub("\\.", ":", HEADER)

first <- nrow(raw)
second <- ncol(raw)
third <- c("P", sapply(strsplit(colnames(raw), split=":"), "[", 2))
fourth <- paste0(rep("S", 10), collapse = "")

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

