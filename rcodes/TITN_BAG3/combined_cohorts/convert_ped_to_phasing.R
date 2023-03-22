setwd("Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3")
ped <- read.table("haplotype_input.ped")
map <- read.table("haplotype_input.map")
bim <- read.table("ttn_significant_ld.bim")

sum(map$V2 %in% bim$V2)
# 10

map$REF <- bim$V6[match(map$V2, bim$V2)]
map$ALT <- bim$V5[match(map$V2, bim$V2)]

first <- nrow(ped)
second <- nrow(map)
third <- c("P", map$V4)
fourth <- paste0(rep("S", 10), collapse = "")

write(first, file="haplotype_input_edited.txt", append = T)
write(second, file="haplotype_input_edited.txt", append = T)
write(paste0(third, collapse = " "), file="haplotype_input_edited.txt", append = T)
write(fourth, file="haplotype_input_edited.txt", append = T)

rownames(ped) <- ped$V1
ped <- ped[-c(1:6)]
# View(ped)

for (i in 1:nrow(ped)){
print(paste0("Doing row: ", i))  
sample <- paste0("#", rownames(ped)[i])
write(sample, file="haplotype_input_edited.txt", append = T)
ped.get.1 <- gsub("0", "?", ped[i,][seq_along(ped[i,]) %% 2 > 0])
ped.get.1[ped.get.1 == map$REF] <- 0
ped.get.1[ped.get.1 == map$ALT] <- 1
line1 <- paste0(ped.get.1, collapse = " ")

ped.get.2 <- gsub("0", "?", ped[i,][seq_along(ped[i,]) %% 2 < 1])
ped.get.2[ped.get.2 == map$REF] <- 0
ped.get.2[ped.get.2 == map$ALT] <- 1
line2 <- paste0(ped.get.2, collapse = " ")
write(line1, file="haplotype_input_edited.txt", append = T)
write(line2, file="haplotype_input_edited.txt", append = T)
}

