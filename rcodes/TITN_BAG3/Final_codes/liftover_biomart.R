library(biomaRt)

## Select Ensembl dataset
# ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl", version = "GRCh38")
## Connect to the Ensembl database for human variants
# For GRCh37 (old genome build)
mart_grch37 <- useEnsembl(biomart = "ENSEMBL_MART_SNP", dataset = "hsapiens_snp", GRCh = 37)

# For GRCh38 (new genome build)
mart_grch38 <- useEnsembl(biomart = "ENSEMBL_MART_SNP", dataset = "hsapiens_snp")
# mart_grch38 <- useEnsembl(biomart = "ENSEMBL_MART_SNP", dataset = "hsapiens_snp", mirror = "useast")


# List available filters for GRCh38 mart
listFilters(mart_grch38)

# List available attributes for GRCh38 mart
listAttributes(mart_grch38)

# Define chromosome and position ranges based on previous results
chromosome <- "1"  # Example chromosome
start_position <- 156138821
end_position <- 156138821

# Create a data frame in R
df <- data.frame(
  Chromosome = rep("2", 16),
  start_position = c(179399704, 179400742, 179410112, 179414849, 179422284, 179425091, 179428124, 
            179432234, 179435679, 179441250, 179446855, 179453355, 179478777, 
            179571683, 179604819, 179631116),
  end_position = c(179399704, 179400742, 179410112, 179414849, 179422284, 179425091, 179428124, 
          179432234, 179435679, 179441250, 179446855, 179453355, 179478777, 
          179571683, 179604819, 179631116)
)

df$end_position <- df$end_position+1
# Print the data frame
# print(df)

chromosome <- df$Chromosome
start_position <- df$start_position
end_position <- df$end_position

# Query variants in GRCh37
variants_grch37 <- getBM(attributes = c('chr_name', 'chrom_start', 'snp'),
                         filters = c('chr_name', 'start'),
                         values = list(chromosome, start_position),
                         mart = mart_grch37)


variants_grch38 <- getBM(attributes = c('chr_name', 'chrom_start', 'snp'),
                         filters = c('chr_name', 'start'),
                         values = list(chromosome, start_position),
                         mart = mart_grch38)


# Merge data based on 'refsnp_id' to compare positions between builds
merged_variants <- merge(variants_grch37, variants_grch38, by = "refsnp_id", suffixes = c("_grch37", "_grch38"))

# View merged variants
print(merged_variants)
