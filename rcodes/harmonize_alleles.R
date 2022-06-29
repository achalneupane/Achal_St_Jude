# Process the variants to match their alleles

# Read the input file including variants with no direct match of their alleles
args = commandArgs(trailingOnly = TRUE)
dat = read.table(args[1], header = FALSE, stringsAsFactors = FALSE, colClasses = c("character"))
# dat = read.table("Z:/ResearchHome/ClusterHome/ysapkota/Work/CAD_PRS/y", header = FALSE, stringsAsFactors = FALSE)

## Function to flip alleles
flip_alleles = function(x){
  if(x=="A"){
    y="T"
  } else if (x=="T"){
    y="A"
  } else if (x=="C"){
    y="G"
  } else if (x=="G"){
    y="C"
  }
  return(y)
}

# Process each variant to check the alleles
dat.out = NULL
for (i in 1:nrow(dat)){
  chr=dat$V1[i]
  pos = dat$V3[i]
  variant = dat$V2[i]
  khera_a1 = dat$V6[i]
  khera_a2 = dat$V7[i]
  # khera_weight = dat$V9[i]
  wgs_a1 = dat$V4[i]
  wgs_a2 = dat$V5[i]
  # First find out if the wgs_a1 have more than one character
  if(nchar(wgs_a1)>1){
    wgs_a1_first = substr(wgs_a1, 1, 1)
    wgs_a1_last = substr(wgs_a1, nchar(wgs_a1), nchar(wgs_a1))
    wgs_a1_changed = ifelse(wgs_a1_first==wgs_a2, wgs_a1_last, wgs_a1_first)
  } else {
    wgs_a1_changed = wgs_a1
  }
  # Then do the same for wgs_a2 allele
  if (nchar(wgs_a2)>1){
    wgs_a2_first = substr(wgs_a2, 1, 1)
    wgs_a2_last = substr(wgs_a2, nchar(wgs_a2), nchar(wgs_a2))
    wgs_a2_changed = ifelse(wgs_a2_first==wgs_a1, wgs_a2_last, wgs_a2_first)
  } else {
    wgs_a2_changed = wgs_a2
  }
  # Now check if the changed wgs alleles match with those from Khera et al
  # No alleles flipped
  if ((wgs_a1_changed == khera_a1 & wgs_a2_changed == khera_a2) | (wgs_a1_changed == khera_a2 & wgs_a2_changed == khera_a1)){
    wgs_a1_new = wgs_a1_changed; wgs_a2_new = wgs_a2_changed; match=1
  } else if((flip_alleles(wgs_a1_changed) == khera_a1 & wgs_a2_changed == khera_a2) | (flip_alleles(wgs_a1_changed) == khera_a2 & wgs_a2_changed == khera_a1)) { # only a1 flipped
    wgs_a1_new = flip_alleles(wgs_a1_changed); wgs_a2_new = wgs_a2_changed; match=1
  } else if ((wgs_a1_changed == khera_a1 & flip_alleles(wgs_a2_changed) == khera_a2) | (wgs_a1_changed == khera_a2 & flip_alleles(wgs_a2_changed) == khera_a1)) { # only a2 flipped
    wgs_a1_new = wgs_a1_changed; wgs_a2_new = flip_alleles(wgs_a2_changed); match=1
  } else if ((flip_alleles(wgs_a1_changed) == khera_a1 & flip_alleles(wgs_a2_changed) == khera_a2) | (flip_alleles(wgs_a1_changed) == khera_a2 & flip_alleles(wgs_a2_changed) == khera_a1)){ # both a1 and a2 flipped
    wgs_a1_new = flip_alleles(wgs_a1_changed); wgs_a2_new = flip_alleles(wgs_a2_changed); match=1
  } else {
    wgs_a1_new = wgs_a1_changed; wgs_a2_new = wgs_a2_changed; match=0
  }
  dat.out = rbind(dat.out, data.frame(chr, pos, variant, khera_a1, khera_a2, wgs_a1, wgs_a2, wgs_a1_new, wgs_a2_new, match))
}

# Write data to disc
write.table(dat.out, paste0(args[1], "_alleles_harmonized"), row.names = FALSE, quote = FALSE)
