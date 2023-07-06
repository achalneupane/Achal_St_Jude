df <- read.table("z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3/sjlife_ccss_org_ccss_exp_samples.bim")
warnings <- read.table("z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3/warnings.txt")
dim(df)
df$KEY <- paste0("chr",df$V1, ":", df$V4)

df$warnings <- df$KEY %in% warnings$V1

library(dplyr)

# Assuming your data frame is named "df" and the column is named "column"
df_with_count <- df %>%
  group_by(KEY) %>%
  mutate(count = n())

df$KEY2 <- df$V2
problems <- df[df$warnings,]  

## work on allele swap
dat <- problems

dat <- dat [! (grepl("DEL", dat$V5) | grepl("DEL", dat$V6)),]


flip_alleles <- function(x) {
  flipped <- gsub("A", "1", x)
  flipped <- gsub("T", "2", flipped)
  flipped <- gsub("C", "3", flipped)
  flipped <- gsub("G", "4", flipped)
  flipped <- gsub("1", "T", flipped)
  flipped <- gsub("2", "A", flipped)
  flipped <- gsub("3", "G", flipped)
  flipped <- gsub("4", "C", flipped)
  return(flipped)
}


dat$a1 <- dat$V5
dat$a2 <- dat$V6
dat$a1.flipped <- sapply(dat$a1, flip_alleles)
dat$a2.flipped <- sapply(dat$a2, flip_alleles)

dat <- dat[1:1000,]


dat$out.indicator <- ""  # Initialize the indicator column
dat$out.new.a1 <- ifelse(dat$a1 == dat$a2, dat$a1, NA)  # Initialize the new.a1 column
dat$out.new.a2 <- ifelse(dat$a1 == dat$a2, dat$a2, NA)  # Initialize the new.a2 column

for (i in 1:nrow(dat)) {
  print(paste0("Doing row ", i))
  
  chr <- dat$V1[i]
  pos <- dat$V4[i]
  KEYpos <- dat$KEY[i]
  
  new.dat <- dat[grepl(KEYpos, dat$KEY), ]
  
  for (j in 1:nrow(new.dat)) {
    if (j == 1) {
      next  # Skip the first row
    }
    
    # Check if a1 and a2 already match
    if (new.dat$a1[1] == new.dat$a1[j] && new.dat$a2[1] == new.dat$a2[j]) {
      new.dat$out.indicator[j] <- "match"  # Update the indicator column for the first row
      new.dat$out.new.a1[j] <- new.dat$a1[j]  # Update the new.a1 column with original a1 allele for the first row
      new.dat$out.new.a2[j] <- new.dat$a2[j]  # Update the new.a2 column with original a2 allele for the first row
    }
    
    if (new.dat$out.indicator[j] == "match") {
      next  # Skip the row if it's already a match
    }
    
    if (new.dat$a1[1] == new.dat$a2[j] && new.dat$a2[1] == new.dat$a1[j]) {
      new.dat$out.indicator[j] <- "swapped"  # Update the indicator column
      new.dat$out.new.a1[j] <- new.dat$a1[j]  # Update the new.a1 column with original a1 allele
      new.dat$out.new.a2[j] <- new.dat$a2[j]  # Update the new.a2 column with original a2 allele
    } else if (new.dat$out.indicator[j] != "swapped" &&
               new.dat$a1[1] == flip_alleles(new.dat$a2[j]) && new.dat$a2[1] == flip_alleles(new.dat$a1[j])) {
      new.dat$out.indicator[j] <- "swapped_flipped"  # Update the indicator column
      new.dat$out.new.a1[j] <- flip_alleles(new.dat$a2[j])  # Update the new.a1 column with flipped allele
      new.dat$out.new.a2[j] <- flip_alleles(new.dat$a1[j])  # Update the new.a2 column with flipped allele
    }
  }
  
  
  
  # Update the indicator, new.a1, and new.a2 columns in the original dataframe
  dat[grepl(KEYpos, dat$KEY), "out.indicator"] <- new.dat$out.indicator
  dat[grepl(KEYpos, dat$KEY), "out.new.a1"] <- new.dat$out.new.a1
  dat[grepl(KEYpos, dat$KEY), "out.new.a2"] <- new.dat$out.new.a2
}
