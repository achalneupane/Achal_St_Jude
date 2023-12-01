# Specify the file path
setwd("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/Survivor_WES/annotation/snpEff_round2/loftee")
file_path <- "chr20.Survivor_WES.GATK4180.hg38_biallelic.vcf-annot-snpeff-dbnsfp-ExAC.0.3-clinvar.GRCh38.with.gnomAD.revel.loftee.tsv"

# Specify the chunk size
chunk_size <- 1000

# Read the first 1000 lines from the file
lines <- readLines(file_path, n = 1000)
# Filter lines starting with a single "#"
header_lines <- grep("^#[^#]", lines)
# Print or use the extracted header lines
header_vector <- gsub("#", "", lines[header_lines])
header <- read.delim(text = header_vector, header = T, sep = "\t")
header <- colnames(header)

# Open the file connection
con <- file(file_path, open = "r")

# Read and process the file in chunks
while (length(lines <- readLines(con, n = chunk_size)) > 0) {
  # Skip lines starting with #
  lines <- lines[!grepl("^#", lines)]
  df <- read.table(text = lines, header = F, sep = "\t")
  # Perform an operation on each chunk
  # For example, print the first few lines of each chunk
  print(head(lines))
}

# Close the file connection
close(con)
