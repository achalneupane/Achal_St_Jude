setwd("Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3")
df <- read.table("z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3/sjlife_ccss_org_ccss_exp_samples.bim")
warnings <- read.table("z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3/warnings.txt")  ## These are the list of variants from plink warnings when merging datasets
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

# dat <- dat[1:1000,]


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
    if (new.dat$a1[1] == new.dat$a1[j] && new.dat$a2[1] == new.dat$a2[j]) { ## direct match
      new.dat$out.indicator[1] <- "match_row1"  # Update the indicator column for the first row
      new.dat$out.new.a1[1] <- new.dat$a1[1]  # Update the new.a1 column with original a1 allele for the first row
      new.dat$out.new.a2[1] <- new.dat$a2[1]  # Update the new.a2 column with original a2 allele for the first row
      new.dat$out.indicator[j] <- "match"  # Update the indicator column for the j-th row
      new.dat$out.new.a1[j] <- new.dat$a1[j]  # Update the new.a1 column with original a1 allele for the j-th row
      new.dat$out.new.a2[j] <- new.dat$a2[j]  # Update the new.a2 column with original a2 allele for the j-th row
    } else if (new.dat$a1[1] == new.dat$a2[j] && new.dat$a2[1] == new.dat$a1[j]) { ## swapped
      new.dat$out.indicator[1] <- "swapped_row1"  # Update the indicator column
      new.dat$out.new.a1[1] <- new.dat$a1[1]  # Update the new.a1 column with original a1 allele
      new.dat$out.new.a2[1] <- new.dat$a2[1]  # Update the new.a2 column with original a2 allele
      new.dat$out.indicator[j] <- "swapped"  # Update the indicator column
      new.dat$out.new.a2[j] <- new.dat$a1[j]  # Update the new.a1 column with original a1 allele
      new.dat$out.new.a1[j] <- new.dat$a2[j]  # Update the new.a2 column with original a2 allele
    } else if (new.dat$a1[1] == flip_alleles(new.dat$a2[j]) && new.dat$a2[1] == flip_alleles(new.dat$a1[j])) { ## both swapped and flipped
      new.dat$out.indicator[1] <- "swapped_flipped_row1"  # Update the indicator column
      new.dat$out.new.a2[1] <- new.dat$a2[1]  # Update the new.a1 column with flipped allele
      new.dat$out.new.a1[1] <- new.dat$a1[1]  # Update the new.a2 column with flipped allele
      new.dat$out.indicator[j] <- "swapped_flipped"  # Update the indicator column
      new.dat$out.new.a1[j] <- flip_alleles(new.dat$a2[j])  # Update the new.a1 column with flipped allele
      new.dat$out.new.a2[j] <- flip_alleles(new.dat$a1[j])  # Update the new.a2 column with flipped allele
    } else if (new.dat$a1[1] == flip_alleles(new.dat$a1[j]) && new.dat$a2[1] == new.dat$a2[j]) { # a1 flipped
    new.dat$out.indicator[1] <- "a1_flipped_row1"  # Update the indicator column
    new.dat$out.new.a2[1] <- new.dat$a2[1]  # Update the new.a1 column with flipped allele
    new.dat$out.new.a1[1] <- new.dat$a1[1]  # Update the new.a2 column with flipped allele
    new.dat$out.indicator[j] <- "a1_flipped"  # Update the indicator column
    new.dat$out.new.a1[j] <- flip_alleles(new.dat$a1[j])  # Update the new.a1 column with flipped allele
    new.dat$out.new.a2[j] <- new.dat$a2[j]  # Update the new.a2 column with flipped allele
    } else if (new.dat$a1[1] == new.dat$a1[j] && new.dat$a2[1] == flip_alleles(new.dat$a2[j])) { # a2 flipped
    new.dat$out.indicator[1] <- "a2_flipped"  # Update the indicator column
    new.dat$out.new.a2[1] <- new.dat$a2[1]  # Update the new.a1 column with flipped allele
    new.dat$out.new.a1[1] <- new.dat$a1[1]  # Update the new.a2 column with flipped allele
    new.dat$out.indicator[j] <- "a2_flipped"  # Update the indicator column
    new.dat$out.new.a1[j] <- new.dat$a1[j]  # Update the new.a1 column with flipped allele
    new.dat$out.new.a2[j] <- flip_alleles(new.dat$a2[j])  # Update the new.a2 column with flipped allele
    } else if (new.dat$a1[1] == new.dat$a2[j] && new.dat$a2[1] == flip_alleles(new.dat$a1[j])) { # a1 swapped and flipped
    new.dat$out.indicator[1] <- "a1_flipped_swapped_row1"  # Update the indicator column
    new.dat$out.new.a2[1] <- new.dat$a2[1]  # Update the new.a1 column with flipped allele
    new.dat$out.new.a1[1] <- new.dat$a1[1]  # Update the new.a2 column with flipped allele
    new.dat$out.indicator[j] <- "a1_flipped_swapped"  # Update the indicator column
    new.dat$out.new.a1[j] <- new.dat$a2[j]  # Update the new.a1 column with flipped allele
    new.dat$out.new.a2[j] <- flip_alleles(new.dat$a1[j])  # Update the new.a2 column with flipped allele
    } else if (new.dat$a1[1] == flip_alleles(new.dat$a2[j]) && new.dat$a2[1] == new.dat$a1[j]) { # a2 swapped and flipped
    new.dat$out.indicator[1] <- "a2_flipped_swapped_row1"  # Update the indicator column
    new.dat$out.new.a2[1] <- new.dat$a2[1]  # Update the new.a1 column with flipped allele
    new.dat$out.new.a1[1] <- new.dat$a1[1]  # Update the new.a2 column with flipped allele
    new.dat$out.indicator[j] <- "a2_flipped_swapped"  # Update the indicator column
    new.dat$out.new.a2[j] <- new.dat$a1[j]  # Update the new.a1 column with flipped allele
    new.dat$out.new.a1[j] <- flip_alleles(new.dat$a2[j])  # Update the new.a2 column with flipped allele
  }
  }
  
  
  
  # Update the indicator, new.a1, and new.a2 columns in the original dataframe
  dat[grepl(KEYpos, dat$KEY), "out.indicator"] <- new.dat$out.indicator
  dat[grepl(KEYpos, dat$KEY), "out.new.a1"] <- new.dat$out.new.a1
  dat[grepl(KEYpos, dat$KEY), "out.new.a2"] <- new.dat$out.new.a2
}

dat.final <- dat[!dat$out.indicator=="",]

first.row <- dat.final[grepl("_row", dat.final$out.indicator),]
dat.final <- dat.final[!grepl("_row", dat.final$out.indicator),]


## extract variants
write.table(dat.final$V2, "extract_vars.txt", col.names = F, quote = F, row.names = F)
## update alleles
update_allele <- cbind.data.frame(dat.final$V2, dat.final$a2, dat.final$a1, dat.final$out.new.a2, dat.final$out.new.a1)
write.table(update_allele, "update_alleles.txt", col.names = F, quote = F, row.names = F)

## update names
update_names <- cbind.data.frame(OLD=dat.final$V2, NEW=paste0(dat.final$KEY, ":", dat.final$out.new.a2, ":", dat.final$out.new.a1))

## check if names are same
# update_names$MATCH <- update_names[1] == update_names[2]

# update_names$NEW[grepl("chr2:178316698:",update_names$OLD)] <- "chr2:178316698:C:T"
# update_names$NEW[grepl("chr2:178313779:",update_names$OLD)] <- "chr2:178313779:A:G"

write.table(update_names, "update_names.txt", col.names = F, quote = F, row.names = F)

# ## update names in the good ones
update_names_good_snps <- cbind.data.frame(OLD = first.row$V2, NEW=paste0(first.row$KEY, ":", first.row$out.new.a2, ":", first.row$out.new.a1))
write.table(update_names_good_snps, "update_names_good_snps.txt", col.names = F, quote = F, row.names = F)
