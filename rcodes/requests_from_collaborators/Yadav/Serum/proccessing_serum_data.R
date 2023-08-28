setwd("Z:/ResearchHome/ClusterHome/aneupane/St_Jude/Achal_St_Jude/rcodes/requests_from_collaborators/Yadav/Serum")
all.df <- read.table("All_serum_samples_for_R01_22Aug2023_for_Achal_tab.txt", sep = "\t", header = T, check.names = F, fileEncoding = "UTF-8")
dim(all.df)
head(all.df)

all.df$original_Sample_age <- all.df$Sample_age
all.df$original_ageevent <- all.df$ageevent

all.df$Sample_age <- as.numeric(all.df$Sample_age)

# Round age down to one decimal place, so easier to compare
all.df$Sample_age <- floor(all.df$Sample_age * 10) / 10
all.df$ageevent <- floor(all.df$ageevent * 10) / 10


##################################################
## 1. First process df without ageevent missing ##
##################################################
# Remove rows with missing ageevent and process them separately
missing.age.samples <- all.df$sjlid[is.na(all.df$ageevent)]
## Process them separately
df <- all.df[!all.df$sjlid %in% missing.age.samples,]

# Get unique SJLIDs
SJLID <- unique(df$sjlid)

FINAL <- {}

# grep("SJL1748413", SJLID)

## Loop through each SJLID
# i = 732
for (i in 1:length(SJLID)) { 
tmp.wanted <- df[SJLID[i] == df$sjlid,]

# Only keep rows with ageevent >= Sample_age
tmp.wanted <- tmp.wanted[!is.na(tmp.wanted$Sample_age) & tmp.wanted$ageevent >= tmp.wanted$Sample_age, ]

# If all grades start with grade 2 or higher, move to the next sample
if (sum(is.na(tmp.wanted$grade)) > 0 || sum(tmp.wanted$grade < 2) == 0) {
  next
  print(paste0("!!! Skipping sample: ", SJLID[i], "; Total rows: ", sum(grepl(SJLID[i], df$sjlid))))
}


# Sort by ageevent, tb_number, and grade
sorted_data <- tmp.wanted[order(tmp.wanted$ageevent, tmp.wanted$tb_number, tmp.wanted$grade), ]

## Keep rows with max num_vials, for each TB number; then keep ALIVE over DEAD
subset.1 <- sorted_data %>%
  group_by(tb_number) %>%
  filter(num_vials == max(num_vials)) %>%
  filter(Survival_Status == "ALIVE" | all(Survival_Status == "DEAD")) %>%
  ungroup()

# subset.1

# If all grades are 0; keep rows with smallest ageevent
if (all(subset.1$grade == 0)) {
  smallest_age_rows <- subset.1 %>%
    group_by(tb_number) %>%
    filter(ageevent == min(ageevent)) %>%
    ungroup()
} else {
  # Initialize variables
  prev_grade <- subset.1$grade[1]
  rows_to_keep <- c(1)
  
  # Loop through rows
  for (j in 2:nrow(subset.1)) {
    if (subset.1$grade[j] >= prev_grade || prev_grade == 0) {
      rows_to_keep <- c(rows_to_keep, j)
    } else {
      break
    }
    prev_grade <- subset.1$grade[j]
  }
  
  # Find unique tb_numbers and filter rows
  unique_tb_numbers <- subset.1$tb_number[rows_to_keep]
  smallest_age_rows <- subset.1 %>%
    filter(tb_number %in% unique_tb_numbers) %>%
    group_by(tb_number) %>%
    filter(ageevent == min(ageevent)) %>%
    ungroup()
}

## If all first events have grades 2 or higher, nullify the dataframe
if (!any(smallest_age_rows$grade == 0)) {
  smallest_age_rows <- data.frame()
}

print(paste0("Processing sample: ", SJLID[i], " in i ", i, "; processing total rows: ", sum(grepl(SJLID[i], df$sjlid)), "; final row counts: ", nrow(smallest_age_rows)))
## Display results for the sample
# print(sorted_data)
# print(smallest_age_rows)
FINAL <- rbind.data.frame(FINAL, smallest_age_rows)
# Sys.sleep(2)
}



############################################
## 2. Now work on age event missing ones, ##
############################################
# if ageevent is missing and grade is zero, keep (eg. SJL5237316, SJL5257302)
df <- all.df[all.df$sjlid %in% missing.age.samples,]
sum(is.na(df$ageevent)) # 323
sum(is.na(df$ageevent) & df$grade == 0) # 323


# Get unique SJLIDs
SJLID <- unique(df$sjlid)


FINAL.missing <- {}

## Loop through each SJLID
## i=22
for (i in 1:length(SJLID)) {
  tmp.wanted <- df[grepl(SJLID[i], df$sjlid),]
  
  ## Only keep rows with ageevent >= Sample_age
  # tmp.wanted <- tmp.wanted[!is.na(tmp.wanted$Sample_age) & tmp.wanted$ageevent >= tmp.wanted$Sample_age, ]
  
  # # If all grades start with grade 2 or higher, move to the next sample
  if (sum(is.na(tmp.wanted$grade)) > 0 || sum(tmp.wanted$grade < 2) == 0) {
    next
    print(paste0("!!! Skipping sample: ", SJLID[i], "; Total rows: ", sum(grepl(SJLID[i], df$sjlid))))
  }
  
  # Sort by ageevent, tb_number, and grade
  sorted_data <- tmp.wanted[order(tmp.wanted$ageevent, tmp.wanted$tb_number, tmp.wanted$grade), ]
  
  ## Keep rows with max num_vials, for each TB number; then keep ALIVE over DEAD
  subset.1 <- sorted_data %>%
    group_by(tb_number) %>%
    filter(num_vials == max(num_vials)) %>%
    filter(Survival_Status == "ALIVE" | all(Survival_Status == "DEAD")) %>%
    ungroup()
  
  subset.1
  
  # If all grades are 0; keep rows with smallest Sample_age; this is different from 1
  if (all(subset.1$grade == 0)) {
    smallest_age_rows <- subset.1 %>%
      group_by(tb_number) %>%
      filter(Sample_age == min(Sample_age)) %>%
      ungroup()
  } else {
    # Initialize variables
    prev_grade <- subset.1$grade[1]
    rows_to_keep <- c(1)
    
    # Loop through rows
    for (j in 2:nrow(subset.1)) {
      if (subset.1$grade[j] >= prev_grade || prev_grade == 0) {
        rows_to_keep <- c(rows_to_keep, j)
      } else {
        break
      }
      prev_grade <- subset.1$grade[j]
    }
    
    # Find unique tb_numbers and filter rows
    unique_tb_numbers <- subset.1$tb_number[rows_to_keep]
    smallest_age_rows <- subset.1 %>%
      filter(tb_number %in% unique_tb_numbers) %>%
      group_by(tb_number) %>%
      filter(Sample_age == min(Sample_age)) %>%
      ungroup()
  }
  
  ## If all first events have grades 2 or higher, nullify the dataframe
  if (!any(smallest_age_rows$grade == 0)) {
    smallest_age_rows <- data.frame()
  }
  
  print(paste0("Processing sample: ", SJLID[i], " in i ", i, "; processing total rows: ", sum(grepl(SJLID[i], df$sjlid)), "; final row counts: ", nrow(smallest_age_rows)))
  ## Display results for the sample
  # print(sorted_data)
  # print(smallest_age_rows)
  FINAL.missing <- rbind.data.frame(FINAL.missing, smallest_age_rows)
  # Sys.sleep(2)
}


## Merge 1 and 2.
FINAL.ALL <- rbind.data.frame(FINAL, FINAL.missing) 
write.table(FINAL.ALL, file = "serum_data_processed.txt", sep = "\t",  row.names = FALSE, col.names = TRUE, quote = F)
