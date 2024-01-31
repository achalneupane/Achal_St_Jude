



# Install and load the biomaRt package
if (!requireNamespace("biomaRt", quietly = TRUE)) {
  BiocManager::install("biomaRt")
}

library(biomaRt)

# Specify the Ensembl database and dataset
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# Get gene coordinates for ABCB11 and ACD
genes_of_interest <- c("ABCB11", "ACD")
gene_info <- getBM(
  attributes = c("external_gene_name", "chromosome_name", "start_position", "end_position", "strand"),
  filters = "external_gene_name",
  values = genes_of_interest,
  mart = ensembl
)

print(gene_info)
ensembl <- useMart("ensembl")
ensembl <- useDataset("hsapiens_gene_ensembl",mart=ensembl)
getBM(attributes=c('chromosome_name', 'start_position', 'end_position', 'strand'),
      filters=c('hgnc_symbol'),
      values=list('TTN'),
      mart=ensembl)
