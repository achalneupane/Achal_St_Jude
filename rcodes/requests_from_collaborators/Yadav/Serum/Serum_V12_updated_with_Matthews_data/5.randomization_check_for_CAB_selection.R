# Load the libraries
library(dplyr)
library(ggplot2)
library(knitr)

cabsamples <- read.table("Z:/ResearchHome/ClusterHome/aneupane/data/Yadav_serum/v12_output/Approved_Samples(N=1200)_to_core.txt", header = T, sep = "\t")

df <- read.table("Z:/ResearchHome/ClusterHome/aneupane/data/Yadav_serum/v12_output/CAB_randomized_batch1_1200_samples/Re_ SJLIFE samples for proteomics_metabolomics experiments/Output_for_Experiment_Design__v0.0.1.txt", header = T, sep = "\t")
df <- df[!grepl("IR", df$sampleID),]
# df.check <- df
# df.check$selection_group <- cabsamples$selection_group[match(df.check$sampleID, cabsamples$tb_number)]
# df.check$sampleAge <- cabsamples$Sample_age[match(df.check$sampleID, cabsamples$tb_number)]
# df.check$Race <- cabsamples$racegrp[match(df.check$sampleID, cabsamples$tb_number)]
# df.check$Sex <- cabsamples$Sex[match(df.check$sampleID, cabsamples$tb_number)]
# table(df.check$factor1 == df.check$selection_group)
# table(df.check$factor2 == df.check$Sex)
# table(df.check$factor4 == df.check$Race)

df$batchID <- factor(df$batchID, levels = unique(c(df$batchID)))

# 1. Check distribution of factor1 across batches
factor1_batch_table <- table(df$factor1, df$batchID)

# print(factor1_batch_table)
knitr::kable(t(factor1_batch_table), caption = "Distribution of Factor 1 Across Batches", format = "markdown")

library(ggplot2)

# Creating proportional table for each factor
## Factor 1
prop_table <- prop.table(factor1_batch_table, margin = 2)
# Convert to a data frame for ggplot2
prop_df <- as.data.frame(as.table(prop_table))
colnames(prop_df) <- c("Factor1", "Batch", "Proportion")

# Create the stacked bar plot
f1 <- ggplot(prop_df, aes(x = Batch, y = Proportion, fill = Factor1)) +
  geom_bar(stat = "identity") +
  labs(title = "Proportion of Selection-group Across Batches", x = "Batch", y = "Proportion") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
f1

# 2. Check distribution of factor2 across batches
factor2_batch_table <- table(df$factor2, df$batchID)
# print(factor2_batch_table)
knitr::kable(t(factor2_batch_table), caption = "Distribution of Factor 2 Across Batches", format = "markdown")

## Factor 2
prop_table <- prop.table(factor2_batch_table, margin = 2)
# Convert to a data frame for ggplot2
prop_df <- as.data.frame(as.table(prop_table))
colnames(prop_df) <- c("Factor2", "Batch", "Proportion")

# Create the stacked bar plot
f2 <- ggplot(prop_df, aes(x = Batch, y = Proportion, fill = Factor2)) +
  geom_bar(stat = "identity") +
  labs(title = "Proportion of Sex Across Batches", x = "Batch", y = "Proportion") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
f2

# 3. Check distribution of factor3 across batches
df$factor3 <- factor(df$factor3, levels = c("[18,24]", "(24,30]", "(30,36]", "(36,42]", "(42,48]", "(48,54]", "(54,60]", "(60,66]"))
factor3_batch_table <- table(df$factor3, df$batchID)
# print(factor3_batch_table)
knitr::kable(t(factor3_batch_table), caption = "Distribution of Factor 3 Across Batches", format = "markdown")

## Factor 3
prop_table <- prop.table(factor3_batch_table, margin = 2)
# Convert to a data frame for ggplot2
prop_df <- as.data.frame(as.table(prop_table))
colnames(prop_df) <- c("Factor3", "Batch", "Proportion")

# Create the stacked bar plot
f3 <- ggplot(prop_df, aes(x = Batch, y = Proportion, fill = Factor3)) +
  geom_bar(stat = "identity") +
  labs(title = "Proportion of Age Across Batches", x = "Batch", y = "Proportion") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
f3

# 4. Check distribution of factor4 across batches
factor4_batch_table <- table(df$factor4, df$batchID)
# print(factor4_batch_table)
knitr::kable(t(factor4_batch_table), caption = "Distribution of Factor 4 Across Batches", format = "markdown")

## Factor 4
prop_table <- prop.table(factor4_batch_table, margin = 2)
# Convert to a data frame for ggplot2
prop_df <- as.data.frame(as.table(prop_table))
colnames(prop_df) <- c("Factor4", "Batch", "Proportion")

# Create the stacked bar plot
f4 <- ggplot(prop_df, aes(x = Batch, y = Proportion, fill = Factor4)) +
  geom_bar(stat = "identity") +
  labs(title = "Proportion of Race Across Batches", x = "Batch", y = "Proportion") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
f4
