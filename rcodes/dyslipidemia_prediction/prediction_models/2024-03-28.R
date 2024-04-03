setwd("D:/Dropbox/R")
library(dplyr)
library(openxlsx)


###################################################################################################
###################################################################################################
#################agedx categories

df <- read.xlsx("dyslp_merged_data_2522.xlsx")
df0 <- df[df$dyslipidemia == 0, ]
df1 <- df[df$dyslipidemia == 1, ]

# Define breaks for the categories
breaks <- c(-Inf, 5, 10, 15, Inf)  # Breaks for the ranges: [0, 5], (5, 10], (10, 15], (15, Inf)
labels = c("<=5", ">5-10", ">10-15", ">15") # Labels for the breaks

# Coerce "ant_count" column to numeric
df0$agedx <- as.numeric(df0$agedx)
df1$agedx <- as.numeric(df1$agedx)

# Create categories based on the breaks
df0$categories <- cut(df0$agedx, breaks = breaks, labels = labels)
df1$categories <- cut(df1$agedx, breaks = breaks, labels = labels)

# Calculate the count of each category
count_by_category <- table(df0$categories)
count_by_category1 <- table(df1$categories)

# Calculate the percentage of each category
percentage_by_category <- round(prop.table(count_by_category) * 100, 1)
percentage_by_category1 <- round(prop.table(count_by_category1) * 100, 1)

# Create a dataframe with results
#results_df_agedx <- data.frame(Breaks = labels, Count_Percentage_agedx0 = paste0(as.numeric(count_by_category), " (", percentage_by_category, "%)"),Count_Percentage_agedx1 = paste0(as.numeric(count_by_category1), " (", percentage_by_category1, "%)"))

results_df_agedx <- data.frame(Breaks = labels, Count_Percentage_agedx1 = paste0(as.numeric(count_by_category1), " (", percentage_by_category1, "%)"),Count_Percentage_agedx0 = paste0(as.numeric(count_by_category), " (", percentage_by_category, "%)"))



###################################################################################################
###################################################################################################
#################age_base

df <- read.xlsx("dyslp_merged_data_2522.xlsx")
df0 <- df[df$dyslipidemia == 0, ]
df1 <- df[df$dyslipidemia == 1, ]

# Define breaks for the categories
breaks <- c(-Inf, 15, 25, 35, 45, Inf)  # Breaks for the ranges: [0, 15], (15, 25], (25, 35], (35, 45], (45, Inf)
labels = c("<=15", ">15-25", ">25-35", ">35-45", ">45") # Labels for the breaks

# Coerce "ant_count" column to numeric
df0$age_base <- as.numeric(df0$age_base)
df1$age_base <- as.numeric(df1$age_base)

# Create categories based on the breaks
df0$categories <- cut(df0$age_base, breaks = breaks, labels = labels)
df1$categories <- cut(df1$age_base, breaks = breaks, labels = labels)

# Calculate the count of each category
count_by_category <- table(df0$categories)
count_by_category1 <- table(df1$categories)

# Calculate the percentage of each category
percentage_by_category <- round(prop.table(count_by_category) * 100, 1)
percentage_by_category1 <- round(prop.table(count_by_category1) * 100, 1)

# Create a dataframe with results
#results_df_age_base <- data.frame(Breaks = labels, Count_Percentage_age_base0 = paste0(as.numeric(count_by_category), " (", percentage_by_category, "%)"),Count_Percentage_age_base1 = paste0(as.numeric(count_by_category1), " (", percentage_by_category1, "%)"))

results_df_age_base <- data.frame(Breaks = labels, Count_Percentage_age_base1 = paste0(as.numeric(count_by_category1), " (", percentage_by_category1, "%)"),Count_Percentage_age_base0 = paste0(as.numeric(count_by_category), " (", percentage_by_category, "%)"))




###################################################################################################
###################################################################################################
#################write all to one file


# Create a list of dataframes
library(openxlsx)

# Create a list of dataframes
list_of_dfs <- list(results_df_agedx,
                    results_df_age_base
                    ) 

# Specify the filename
file_name <- "multiple_dataframes.xlsx"

# Specify the directory where you want to save the file
# Replace "/path/to/your/directory/" with the actual directory path
dir_path <- "D:/Dropbox/R/"

# Create the directory if it doesn't exist
if (!dir.exists(dir_path)) {
  dir.create(dir_path, recursive = TRUE)
}

# Combine the directory path and the file name
file_path <- file.path(dir_path, file_name)

# Create a new workbook
wb <- createWorkbook()

# Create a new sheet
addWorksheet(wb, sheetName = "Combined")

# Track the starting row
current_row <- 1

# Write each data frame to the same sheet
for (i in seq_along(list_of_dfs)) {
  writeData(wb, sheet = "Combined", x = list_of_dfs[[i]], startRow = current_row)
  current_row <- current_row + nrow(list_of_dfs[[i]]) + 1  # Move to the next empty row
  if (i < length(list_of_dfs)) {
    current_row <- current_row + 1  # Add an empty row between dataframes
  }
}

# Save the workbook to the Excel file
saveWorkbook(wb, file_path)

