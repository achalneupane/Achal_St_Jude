library(biomaRt)

## Select Ensembl dataset
# ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl", version = "GRCh38")
## Connect to the Ensembl database for human variants
# For GRCh37 (old genome build)
mart_grch37 <- useEnsembl(biomart = "ENSEMBL_MART_SNP", dataset = "hsapiens_snp", GRCh = 37)

# For GRCh38 (new genome build)
mart_grch38 <- useEnsembl(biomart = "ENSEMBL_MART_SNP", dataset = "hsapiens_snp")

# List available filters for GRCh38 mart
listFilters(mart_grch38)

# List available attributes for GRCh38 mart
listAttributes(mart_grch38)

# Define chromosome and position ranges based on previous results
chromosome <- "1"  # Example chromosome
start_position <- 156138821
end_position <- 156138821



# Query variants in GRCh37
variants_grch37 <- getBM(attributes = c('chr_name', 'chrom_start'),
                         filters = c('chr_name', 'start', 'end'),
                         values = list(chromosome, start_position, end_position),
                         mart = mart_grch37)


variants_grch38 <- getBM(attributes = c('chr_name', 'chrom_start'),
                         filters = c('chr_name', 'start', 'end'),
                         values = list(chromosome, start_position, end_position),
                         mart = mart_grch38)


# Merge data based on 'refsnp_id' to compare positions between builds
merged_variants <- merge(variants_grch37, variants_grch38, by = "refsnp_id", suffixes = c("_grch37", "_grch38"))

# View merged variants
print(merged_variants)
