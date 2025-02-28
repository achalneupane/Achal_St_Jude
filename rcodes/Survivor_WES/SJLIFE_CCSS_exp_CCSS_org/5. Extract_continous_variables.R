library(haven)
library(dplyr)
setwd("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/WES_rare_variant/sjlife/")
SJLIFE <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/WES_rare_variant/sjlife/all_rare_variants_all_sjlife.fam", header = F)
sjlife_4402 <- read_sas("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/PHENOTYPE/Data_from_Qi_Liu/toyadav.sas7bdat")
sjlife_4402 <- unique(sjlife_4402$sjlid)
sjlife_4381 <- sjlife_4402[sjlife_4402 %in% SJLIFE$V2]

####################
## Milli Measures ##
####################
# mili.dat <- read_sas("Z:/SJShare/SJCOMMON/ECC/SJLife/SJLIFE Data Freeze/2 Final Data SJLIFE/20200430/Clinical Data/milli_clinical.sas7bdat")
# BP
milli <- as.data.frame(read_sas("Z:/SJShare/SJCOMMON/ECC/SJLife/SJLIFE Data Freeze/2 Final Data SJLIFE/20200430/Clinical Data/milli_measures.sas7bdat"))
milli <- milli[milli$sjlid %in% SJLIFE$V2,]

wanted.cols <- colnames(milli)[-c(1:6)]


# Initialize lists to store results
base_bp_list <- list()
recent_bp_list <- list()

for (i in seq_along(wanted.cols)) {

print(paste0("Doing: ", wanted.cols[i]))
  
# Select necessary columns
milli2 <- milli %>%
  select(sjlid, labdt, !!sym(wanted.cols[i])) %>%
  mutate(!!sym(wanted.cols[i]) := as.numeric(na_if(as.character(!!sym(wanted.cols[i])), ""))) %>%  # Convert empty strings to NA, then to numeric
  filter(!is.na(!!sym(wanted.cols[i])) & !!sym(wanted.cols[i]) != "")  # filter out any missing or "" invalid values

if (nrow(milli2) == 0 || ncol(milli2) < 3 || length(milli2[,3]) == 0| length(unique(milli2[[3]])) <=3) {
  print(paste("Skipping", wanted.cols[i], "as it has no data"))
  next  # Skip this iteration and move to the next column
}

# Find the baseline (oldest) for each individual
base_bp <- milli2 %>%
  group_by(sjlid) %>%
  slice_min(labdt, with_ties = FALSE) %>%  # this picks the row with the minimum (oldest) labdt for each individual
  ungroup()  # ungroup after operation to make sure you don't accidentally retain groupings

# Store the result in the list
base_bp_list[[wanted.cols[i]]] <- base_bp


recent_bp <- milli2 %>%
  group_by(sjlid) %>%
  slice_max(labdt, with_ties = FALSE) %>%  # This picks the row with the maximum (most recent) labdt for each individual
  ungroup()

# Store the result in the list
recent_bp_list[[wanted.cols[i]]] <- recent_bp

}


# Combine baseline data frames
combined_base_bp <- as.data.frame(bind_rows(base_bp_list, .id = "variable"))
combined_base_bp.milli.measure <- combined_base_bp
# Combine recent data frames
combined_recent_bp <- as.data.frame(bind_rows(recent_bp_list, .id = "variable"))
combined_recent_bp.milli.measure <- combined_recent_bp

# cc <- combined_base_bp[combined_base_bp$variable =="DBP_MAX",]
# dim(cc)
# # 5051

## Plot
# Load necessary libraries
library(ggplot2)
library(scales)
library(gridExtra)

# Create a PDF to save the plots
pdf(paste0("Milli_measures_Distributions_of_", length(unique(milli$sjlid)), "_SJLIFE.pdf"), height = 10, width = 16)  # Adjust width for 4-column layout

# Initialize a list to hold the plots
plots <- list()
plot_count <- 0  # Counter for plots

for (i in seq_along(wanted.cols)) {
  print(paste0("Doing: ", wanted.cols[i]))
  # Process for baseline
  df.tmp <- combined_base_bp[combined_base_bp$variable == wanted.cols[i],]
  
  if (all(is.na(df.tmp[[wanted.cols[i]]]))) {
    print(paste("Skipping", wanted.cols[i], "as all values are NA"))
    next  # Skip this iteration
  }
  
  # Extract the column as numeric
  x_vals <- as.numeric(df.tmp[[wanted.cols[i]]])
  x_min <- min(x_vals, na.rm = TRUE)
  x_max <- max(x_vals, na.rm = TRUE)
  step_size <- (x_max - x_min) / 5  # Adjust for number of breaks
  breaks_seq <- seq(x_min, x_max, by = step_size)
  
  baselineP <- ggplot(df.tmp, aes(x = .data[[wanted.cols[i]]])) + 
    geom_histogram(aes(y = after_stat(count) / sum(after_stat(count))), 
                   bins = 30, 
                   fill = "#4E79A7", 
                   color = "white") +
    geom_vline(xintercept = median(x_vals, na.rm = TRUE),  
               color = "#E15759", 
               linetype = "dashed", 
               linewidth = 1) +
    scale_y_continuous(labels = scales::percent_format(), name = "Frequency") +
    scale_x_continuous(name = wanted.cols[i], breaks = breaks_seq) +
    labs(title = paste0("Baseline ", wanted.cols[i]),
         subtitle = paste(nrow(df.tmp), "Participants | Median =", 
                          round(median(x_vals, na.rm = TRUE), 3))) +
    theme_minimal(base_size = 14) +
    theme(panel.grid.minor = element_blank())
  
  plots[[length(plots) + 1]] <- baselineP
  plot_count <- plot_count + 1
  
  # Process for recent
  df.tmp <- combined_recent_bp[combined_recent_bp$variable == wanted.cols[i],]
  
  # Extract the column as numeric
  x_vals <- as.numeric(df.tmp[[wanted.cols[i]]])
  x_min <- min(x_vals, na.rm = TRUE)
  x_max <- max(x_vals, na.rm = TRUE)
  step_size <- (x_max - x_min) / 5  # Adjust for number of breaks
  breaks_seq <- seq(x_min, x_max, by = step_size)
  
  recentP <- ggplot(df.tmp, aes(x = .data[[wanted.cols[i]]])) + 
    geom_histogram(aes(y = after_stat(count) / sum(after_stat(count))), 
                   bins = 30, 
                   fill = "#4E79A7", 
                   color = "white") +
    geom_vline(xintercept = median(x_vals, na.rm = TRUE),  
               color = "#E15759", 
               linetype = "dashed", 
               linewidth = 1) +
    scale_y_continuous(labels = scales::percent_format(), name = "Frequency") +
    scale_x_continuous(name = wanted.cols[i], breaks = breaks_seq) +
    labs(title = paste0("Recent ", wanted.cols[i]),
         subtitle = paste(nrow(df.tmp), "Participants | Median =", 
                          round(median(x_vals, na.rm = TRUE), 3))) +
    theme_minimal(base_size = 14) +
    theme(panel.grid.minor = element_blank())
  
  plots[[length(plots) + 1]] <- recentP
  plot_count <- plot_count + 1
  
  # Arrange and print plots in a 2x2 grid every 4 plots
  if (plot_count %% 4 == 0) {
    grid.arrange(grobs = plots, ncol = 4, nrow = 1)  # 4 columns, 1 row
    plots <- list()  # Clear the list for the next set
    plot_count <- 0  # Reset the counter
  }
}

# Print remaining plots if any
if (length(plots) > 0) {
  grid.arrange(grobs = plots, ncol = 4, nrow = ceiling(length(plots) / 4))
}

# Close the PDF device
dev.off()

################
## Milli Dexa ##
################
milli <- read_sas("Z:/SJShare/SJCOMMON/ECC/SJLife/SJLIFE Data Freeze/2 Final Data SJLIFE/20200430/Clinical Data/milli_dexa.sas7bdat")
milli <- milli[milli$sjlid %in% SJLIFE$V2,]

wanted.cols <- colnames(milli)[-c(1:6)]


# Initialize lists to store results
base_bp_list <- list()
recent_bp_list <- list()

for (i in seq_along(wanted.cols)) {
  
  print(paste0("Doing: ", wanted.cols[i]))
  
  # Select necessary columns
  milli2 <- milli %>%
    select(sjlid, labdt, !!sym(wanted.cols[i])) %>%
    mutate(!!sym(wanted.cols[i]) := as.numeric(na_if(as.character(!!sym(wanted.cols[i])), ""))) %>%  # Convert empty strings to NA, then to numeric
    filter(!is.na(!!sym(wanted.cols[i])) & !!sym(wanted.cols[i]) != "")  # filter out any missing or "" invalid values

  if (nrow(milli2) == 0 || ncol(milli2) < 3 || length(milli2[,3]) == 0| length(unique(milli2[[3]])) <=3) {
    print(paste("Skipping", wanted.cols[i], "as it has no data"))
    next  # Skip this iteration and move to the next column
  }
  
  
  # Find the baseline (oldest) for each individual
  base_bp <- milli2 %>%
    group_by(sjlid) %>%
    slice_min(labdt, with_ties = FALSE) %>%  # this picks the row with the minimum (oldest) labdt for each individual
    ungroup()  # ungroup after operation to make sure you don't accidentally retain groupings
  
  # Store the result in the list
  base_bp_list[[wanted.cols[i]]] <- base_bp
  
  
  recent_bp <- milli2 %>%
    group_by(sjlid) %>%
    slice_max(labdt, with_ties = FALSE) %>%  # This picks the row with the maximum (most recent) labdt for each individual
    ungroup()
  
  # Store the result in the list
  recent_bp_list[[wanted.cols[i]]] <- recent_bp
  
}


# Combine baseline data frames
combined_base_bp <- as.data.frame(bind_rows(base_bp_list, .id = "variable"))
combined_base_bp.milli.dexa <- combined_base_bp
# Combine recent data frames
combined_recent_bp <- as.data.frame(bind_rows(recent_bp_list, .id = "variable"))
combined_recent_bp.milli.dexa <- combined_recent_bp



# Load necessary libraries
library(ggplot2)
library(scales)
library(gridExtra)

# Create a PDF to save the plots
pdf(paste0("Milli_dexa_Distributions_of_", length(unique(milli$sjlid)), "_SJLIFE.pdf"), height = 10, width = 16)  # Adjust width for 4-column layout

# Initialize a list to hold the plots
plots <- list()
plot_count <- 0  # Counter for plots

for (i in seq_along(wanted.cols)) {
  print(paste0("Doing: ", wanted.cols[i]))
  
  # Process for baseline
  df.tmp <- combined_base_bp[combined_base_bp$variable == wanted.cols[i],]
  
  if (all(is.na(df.tmp[[wanted.cols[i]]]))) {
    print(paste("Skipping", wanted.cols[i], "as all values are NA"))
    next  # Skip this iteration
  }
  
  # Extract the column as numeric
  x_vals <- as.numeric(df.tmp[[wanted.cols[i]]])
  x_min <- min(x_vals, na.rm = TRUE)
  x_max <- max(x_vals, na.rm = TRUE)
  step_size <- (x_max - x_min) / 5  # Adjust for number of breaks
  breaks_seq <- seq(x_min, x_max, by = step_size)
  
  baselineP <- ggplot(df.tmp, aes(x = .data[[wanted.cols[i]]])) + 
    geom_histogram(aes(y = after_stat(count) / sum(after_stat(count))), 
                   bins = 30, 
                   fill = "#4E79A7", 
                   color = "white") +
    geom_vline(xintercept = median(x_vals, na.rm = TRUE),  
               color = "#E15759", 
               linetype = "dashed", 
               linewidth = 1) +
    scale_y_continuous(labels = scales::percent_format(), name = "Frequency") +
    scale_x_continuous(name = wanted.cols[i], breaks = breaks_seq) +
    labs(title = paste0("Baseline ", wanted.cols[i]),
         subtitle = paste(nrow(df.tmp), "Participants | Median =", 
                          round(median(x_vals, na.rm = TRUE), 3))) +
    theme_minimal(base_size = 14) +
    theme(panel.grid.minor = element_blank())
  
  plots[[length(plots) + 1]] <- baselineP
  plot_count <- plot_count + 1
  
  # Process for recent
  df.tmp <- combined_recent_bp[combined_recent_bp$variable == wanted.cols[i],]
  
  # Extract the column as numeric
  x_vals <- as.numeric(df.tmp[[wanted.cols[i]]])
  x_min <- min(x_vals, na.rm = TRUE)
  x_max <- max(x_vals, na.rm = TRUE)
  step_size <- (x_max - x_min) / 5  # Adjust for number of breaks
  breaks_seq <- seq(x_min, x_max, by = step_size)
  
  recentP <- ggplot(df.tmp, aes(x = .data[[wanted.cols[i]]])) + 
    geom_histogram(aes(y = after_stat(count) / sum(after_stat(count))), 
                   bins = 30, 
                   fill = "#4E79A7", 
                   color = "white") +
    geom_vline(xintercept = median(x_vals, na.rm = TRUE),  
               color = "#E15759", 
               linetype = "dashed", 
               linewidth = 1) +
    scale_y_continuous(labels = scales::percent_format(), name = "Frequency") +
    scale_x_continuous(name = wanted.cols[i], breaks = breaks_seq) +
    labs(title = paste0("Recent ", wanted.cols[i]),
         subtitle = paste(nrow(df.tmp), "Participants | Median =", 
                          round(median(x_vals, na.rm = TRUE), 3))) +
    theme_minimal(base_size = 14) +
    theme(panel.grid.minor = element_blank())
  
  plots[[length(plots) + 1]] <- recentP
  plot_count <- plot_count + 1
  
  # Arrange and print plots in a 2x2 grid every 4 plots
  if (plot_count %% 4 == 0) {
    grid.arrange(grobs = plots, ncol = 4, nrow = 1)  # 4 columns, 1 row
    plots <- list()  # Clear the list for the next set
    plot_count <- 0  # Reset the counter
  }
}

# Print remaining plots if any
if (length(plots) > 0) {
  grid.arrange(grobs = plots, ncol = 4, nrow = ceiling(length(plots) / 4))
}

# Close the PDF device
dev.off()


################
## Milli Labs ##
################
milli <- read_sas("Z:/SJShare/SJCOMMON/ECC/SJLife/SJLIFE Data Freeze/2 Final Data SJLIFE/20200430/Clinical Data/milli_labs.sas7bdat")
milli <- milli[milli$sjlid %in% SJLIFE$V2,]

wanted.cols <- colnames(milli)[-c(1:6)]


# Initialize lists to store results
base_bp_list <- list()
recent_bp_list <- list()

for (i in seq_along(wanted.cols)) {
  
  print(paste0("Doing: ", wanted.cols[i]))
  
  # Select necessary columns
  milli2 <- milli %>%
    select(sjlid, labdt, !!sym(wanted.cols[i])) %>%
    mutate(!!sym(wanted.cols[i]) := as.numeric(na_if(as.character(!!sym(wanted.cols[i])), ""))) %>%  # Convert empty strings to NA, then to numeric
    filter(!is.na(!!sym(wanted.cols[i])) & !!sym(wanted.cols[i]) != "")  # filter out any missing or "" invalid values
  
  if (nrow(milli2) == 0 || ncol(milli2) < 3 || length(milli2[,3]) == 0| length(unique(milli2[[3]])) <=3) {
    print(paste("Skipping", wanted.cols[i], "as it has no data"))
    next  # Skip this iteration and move to the next column
  }
  
  
  # Find the baseline (oldest) for each individual
  base_bp <- milli2 %>%
    group_by(sjlid) %>%
    slice_min(labdt, with_ties = FALSE) %>%  # this picks the row with the minimum (oldest) labdt for each individual
    ungroup()  # ungroup after operation to make sure you don't accidentally retain groupings
  
  # Store the result in the list
  base_bp_list[[wanted.cols[i]]] <- base_bp
  
  
  recent_bp <- milli2 %>%
    group_by(sjlid) %>%
    slice_max(labdt, with_ties = FALSE) %>%  # This picks the row with the maximum (most recent) labdt for each individual
    ungroup()
  
  # Store the result in the list
  recent_bp_list[[wanted.cols[i]]] <- recent_bp
  
}


# Combine baseline data frames
combined_base_bp <- as.data.frame(bind_rows(base_bp_list, .id = "variable"))
combined_base_bp.milli.labs <- combined_base_bp
# Combine recent data frames
combined_recent_bp <- as.data.frame(bind_rows(recent_bp_list, .id = "variable"))
combined_recent_bp.milli.labs <- combined_recent_bp



# Load necessary libraries
library(ggplot2)
library(scales)
library(gridExtra)

# Create a PDF to save the plots
pdf(paste0("Milli_labs_Distributions_of_", length(unique(milli$sjlid)), "_SJLIFE.pdf"), height = 10, width = 16)  # Adjust width for 4-column layout

# Initialize a list to hold the plots
plots <- list()
plot_count <- 0  # Counter for plots

for (i in seq_along(wanted.cols)) {
  print(paste0("Doing: ", wanted.cols[i]))
  
  # Process for baseline
  df.tmp <- combined_base_bp[combined_base_bp$variable == wanted.cols[i],]
  
  if (all(is.na(df.tmp[[wanted.cols[i]]]))) {
    print(paste("Skipping", wanted.cols[i], "as all values are NA"))
    next  # Skip this iteration
  }
  
  # Extract the column as numeric
  x_vals <- as.numeric(df.tmp[[wanted.cols[i]]])
  x_min <- min(x_vals, na.rm = TRUE)
  x_max <- max(x_vals, na.rm = TRUE)
  step_size <- (x_max - x_min) / 5  # Adjust for number of breaks
  breaks_seq <- seq(x_min, x_max, by = step_size)
  
  baselineP <- ggplot(df.tmp, aes(x = .data[[wanted.cols[i]]])) + 
    geom_histogram(aes(y = after_stat(count) / sum(after_stat(count))), 
                   bins = 30, 
                   fill = "#4E79A7", 
                   color = "white") +
    geom_vline(xintercept = median(x_vals, na.rm = TRUE),  
               color = "#E15759", 
               linetype = "dashed", 
               linewidth = 1) +
    scale_y_continuous(labels = scales::percent_format(), name = "Frequency") +
    scale_x_continuous(name = wanted.cols[i], breaks = breaks_seq) +
    labs(title = paste0("Baseline ", wanted.cols[i]),
         subtitle = paste(nrow(df.tmp), "Participants | Median =", 
                          round(median(x_vals, na.rm = TRUE), 3))) +
    theme_minimal(base_size = 14) +
    theme(panel.grid.minor = element_blank())
  
  plots[[length(plots) + 1]] <- baselineP
  plot_count <- plot_count + 1
  
  # Process for recent
  df.tmp <- combined_recent_bp[combined_recent_bp$variable == wanted.cols[i],]
  
  # Extract the column as numeric
  x_vals <- as.numeric(df.tmp[[wanted.cols[i]]])
  x_min <- min(x_vals, na.rm = TRUE)
  x_max <- max(x_vals, na.rm = TRUE)
  step_size <- (x_max - x_min) / 5  # Adjust for number of breaks
  breaks_seq <- seq(x_min, x_max, by = step_size)
  
  recentP <- ggplot(df.tmp, aes(x = .data[[wanted.cols[i]]])) + 
    geom_histogram(aes(y = after_stat(count) / sum(after_stat(count))), 
                   bins = 30, 
                   fill = "#4E79A7", 
                   color = "white") +
    geom_vline(xintercept = median(x_vals, na.rm = TRUE),  
               color = "#E15759", 
               linetype = "dashed", 
               linewidth = 1) +
    scale_y_continuous(labels = scales::percent_format(), name = "Frequency") +
    scale_x_continuous(name = wanted.cols[i], breaks = breaks_seq) +
    labs(title = paste0("Recent ", wanted.cols[i]),
         subtitle = paste(nrow(df.tmp), "Participants | Median =", 
                          round(median(x_vals, na.rm = TRUE), 3))) +
    theme_minimal(base_size = 14) +
    theme(panel.grid.minor = element_blank())
  
  plots[[length(plots) + 1]] <- recentP
  plot_count <- plot_count + 1
  
  # Arrange and print plots in a 2x2 grid every 4 plots
  if (plot_count %% 4 == 0) {
    grid.arrange(grobs = plots, ncol = 4, nrow = 1)  # 4 columns, 1 row
    plots <- list()  # Clear the list for the next set
    plot_count <- 0  # Reset the counter
  }
}

# Print remaining plots if any
if (length(plots) > 0) {
  grid.arrange(grobs = plots, ncol = 4, nrow = ceiling(length(plots) / 4))
}

# Close the PDF device
dev.off()


#####################
## Milli pulmonary ##
#####################
milli <- read_sas("Z:/SJShare/SJCOMMON/ECC/SJLife/SJLIFE Data Freeze/2 Final Data SJLIFE/20200430/Clinical Data/milli_pulmonary.sas7bdat")
milli <- milli[milli$sjlid %in% SJLIFE$V2,]

wanted.cols <- colnames(milli)[-c(1:6)]


# Initialize lists to store results
base_bp_list <- list()
recent_bp_list <- list()

for (i in seq_along(wanted.cols)) {
  
  print(paste0("Doing: ", wanted.cols[i]))
  
  # Select necessary columns
  milli2 <- milli %>%
    select(sjlid, labdt, !!sym(wanted.cols[i])) %>%
    mutate(!!sym(wanted.cols[i]) := as.numeric(na_if(as.character(!!sym(wanted.cols[i])), ""))) %>%  # Convert empty strings to NA, then to numeric
    filter(!is.na(!!sym(wanted.cols[i])) & !!sym(wanted.cols[i]) != "")  # filter out any missing or "" invalid values
  
  if (nrow(milli2) == 0 || ncol(milli2) < 3 || length(milli2[,3]) == 0| length(unique(milli2[[3]])) <=3) {
    print(paste("Skipping", wanted.cols[i], "as it has no data"))
    next  # Skip this iteration and move to the next column
  }
  
  
  # Find the baseline (oldest) for each individual
  base_bp <- milli2 %>%
    group_by(sjlid) %>%
    slice_min(labdt, with_ties = FALSE) %>%  # this picks the row with the minimum (oldest) labdt for each individual
    ungroup()  # ungroup after operation to make sure you don't accidentally retain groupings
  
  # Store the result in the list
  base_bp_list[[wanted.cols[i]]] <- base_bp
  
  
  recent_bp <- milli2 %>%
    group_by(sjlid) %>%
    slice_max(labdt, with_ties = FALSE) %>%  # This picks the row with the maximum (most recent) labdt for each individual
    ungroup()
  
  # Store the result in the list
  recent_bp_list[[wanted.cols[i]]] <- recent_bp
  
}


# Combine baseline data frames
combined_base_bp <- as.data.frame(bind_rows(base_bp_list, .id = "variable"))
combined_base_bp.milli.pulmonary <- combined_base_bp
# Combine recent data frames
combined_recent_bp <- as.data.frame(bind_rows(recent_bp_list, .id = "variable"))
combined_recent_bp.milli.pulmonary <- combined_recent_bp



# Load necessary libraries
library(ggplot2)
library(scales)
library(gridExtra)

# Create a PDF to save the plots
pdf(paste0("Milli_pulmonary_Distributions_of_", length(unique(milli$sjlid)), "_SJLIFE.pdf"), height = 10, width = 16)  # Adjust width for 4-column layout

# Initialize a list to hold the plots
plots <- list()
plot_count <- 0  # Counter for plots

for (i in seq_along(wanted.cols)) {
  print(paste0("Doing: ", wanted.cols[i]))
  
  # Process for baseline
  df.tmp <- combined_base_bp[combined_base_bp$variable == wanted.cols[i],]
  
  if (all(is.na(df.tmp[[wanted.cols[i]]]))) {
    print(paste("Skipping", wanted.cols[i], "as all values are NA"))
    next  # Skip this iteration
  }
  
  # Extract the column as numeric
  x_vals <- as.numeric(df.tmp[[wanted.cols[i]]])
  x_min <- min(x_vals, na.rm = TRUE)
  x_max <- max(x_vals, na.rm = TRUE)
  step_size <- (x_max - x_min) / 5  # Adjust for number of breaks
  breaks_seq <- seq(x_min, x_max, by = step_size)
  
  baselineP <- ggplot(df.tmp, aes(x = .data[[wanted.cols[i]]])) + 
    geom_histogram(aes(y = after_stat(count) / sum(after_stat(count))), 
                   bins = 30, 
                   fill = "#4E79A7", 
                   color = "white") +
    geom_vline(xintercept = median(x_vals, na.rm = TRUE),  
               color = "#E15759", 
               linetype = "dashed", 
               linewidth = 1) +
    scale_y_continuous(labels = scales::percent_format(), name = "Frequency") +
    scale_x_continuous(name = wanted.cols[i], breaks = breaks_seq) +
    labs(title = paste0("Baseline ", wanted.cols[i]),
         subtitle = paste(nrow(df.tmp), "Participants | Median =", 
                          round(median(x_vals, na.rm = TRUE), 3))) +
    theme_minimal(base_size = 14) +
    theme(panel.grid.minor = element_blank())
  
  plots[[length(plots) + 1]] <- baselineP
  plot_count <- plot_count + 1
  
  # Process for recent
  df.tmp <- combined_recent_bp[combined_recent_bp$variable == wanted.cols[i],]
  
  # Extract the column as numeric
  x_vals <- as.numeric(df.tmp[[wanted.cols[i]]])
  x_min <- min(x_vals, na.rm = TRUE)
  x_max <- max(x_vals, na.rm = TRUE)
  step_size <- (x_max - x_min) / 5  # Adjust for number of breaks
  breaks_seq <- seq(x_min, x_max, by = step_size)
  
  recentP <- ggplot(df.tmp, aes(x = .data[[wanted.cols[i]]])) + 
    geom_histogram(aes(y = after_stat(count) / sum(after_stat(count))), 
                   bins = 30, 
                   fill = "#4E79A7", 
                   color = "white") +
    geom_vline(xintercept = median(x_vals, na.rm = TRUE),  
               color = "#E15759", 
               linetype = "dashed", 
               linewidth = 1) +
    scale_y_continuous(labels = scales::percent_format(), name = "Frequency") +
    scale_x_continuous(name = wanted.cols[i], breaks = breaks_seq) +
    labs(title = paste0("Recent ", wanted.cols[i]),
         subtitle = paste(nrow(df.tmp), "Participants | Median =", 
                          round(median(x_vals, na.rm = TRUE), 3))) +
    theme_minimal(base_size = 14) +
    theme(panel.grid.minor = element_blank())
  
  plots[[length(plots) + 1]] <- recentP
  plot_count <- plot_count + 1
  
  # Arrange and print plots in a 2x2 grid every 4 plots
  if (plot_count %% 4 == 0) {
    grid.arrange(grobs = plots, ncol = 4, nrow = 1)  # 4 columns, 1 row
    plots <- list()  # Clear the list for the next set
    plot_count <- 0  # Reset the counter
  }
}

# Print remaining plots if any
if (length(plots) > 0) {
  grid.arrange(grobs = plots, ncol = 4, nrow = ceiling(length(plots) / 4))
}

# Close the PDF device
dev.off()

###################
## echo research ##
###################
milli <- as.data.frame(read_sas("Z:/SJShare/SJCOMMON/ECC/SJLife/SJLIFE Data Freeze/2 Final Data SJLIFE/20200430/Clinical Data/echo_research.sas7bdat"))
milli <- milli[milli$sjlid %in% SJLIFE$V2,]
colnames(milli)[colnames(milli) =="studydatetime"] <- "labdt"
wanted.cols <- colnames(milli)[-c(1:4)]


# Initialize lists to store results
base_bp_list <- list()
recent_bp_list <- list()

for (i in seq_along(wanted.cols)) {
  
  print(paste0("Doing: ", wanted.cols[i]))
  
  # Select necessary columns
  milli2 <- milli %>%
    select(sjlid, labdt, !!sym(wanted.cols[i])) %>%
    mutate(!!sym(wanted.cols[i]) := as.numeric(na_if(as.character(!!sym(wanted.cols[i])), ""))) %>%  # Convert empty strings to NA, then to numeric
    filter(!is.na(!!sym(wanted.cols[i])) & !!sym(wanted.cols[i]) != "")  # filter out any missing or "" invalid values
  
  if (nrow(milli2) == 0 || ncol(milli2) < 3 || length(milli2[,3]) == 0| length(unique(milli2[[3]])) <=3) {
    print(paste("Skipping", wanted.cols[i], "as it has no data"))
    next  # Skip this iteration and move to the next column
  }
  
  # Find the baseline (oldest) for each individual
  base_bp <- milli2 %>%
    group_by(sjlid) %>%
    slice_min(labdt, with_ties = FALSE) %>%  # this picks the row with the minimum (oldest) labdt for each individual
    ungroup()  # ungroup after operation to make sure you don't accidentally retain groupings
  
  # Store the result in the list
  base_bp_list[[wanted.cols[i]]] <- base_bp
  
  
  recent_bp <- milli2 %>%
    group_by(sjlid) %>%
    slice_max(labdt, with_ties = FALSE) %>%  # This picks the row with the maximum (most recent) labdt for each individual
    ungroup()
  
  # Store the result in the list
  recent_bp_list[[wanted.cols[i]]] <- recent_bp
  
}


# Combine baseline data frames
combined_base_bp <- as.data.frame(bind_rows(base_bp_list, .id = "variable"))
combined_base_bp.echo.research <- combined_base_bp
# Combine recent data frames
combined_recent_bp <- as.data.frame(bind_rows(recent_bp_list, .id = "variable"))
combined_recent_bp.echo.research <- combined_recent_bp

## Plot
# Load necessary libraries
library(ggplot2)
library(scales)
library(gridExtra)

# Create a PDF to save the plots
pdf(paste0("Echo_research_Distributions_of_", length(unique(milli$sjlid)), "_SJLIFE.pdf"), height = 10, width = 16)  # Adjust width for 4-column layout

# Initialize a list to hold the plots
plots <- list()
plot_count <- 0  # Counter for plots

for (i in seq_along(wanted.cols)) {
  print(paste0("Doing: ", wanted.cols[i]))
  # Process for baseline
  df.tmp <- combined_base_bp[combined_base_bp$variable == wanted.cols[i],]
  
  if (all(is.na(df.tmp[[wanted.cols[i]]]))) {
    print(paste("Skipping", wanted.cols[i], "as all values are NA"))
    next  # Skip this iteration
  }
  
  # Extract the column as numeric
  x_vals <- as.numeric(df.tmp[[wanted.cols[i]]])
  x_min <- min(x_vals, na.rm = TRUE)
  x_max <- max(x_vals, na.rm = TRUE)
  step_size <- (x_max - x_min) / 5  # Adjust for number of breaks
  breaks_seq <- seq(x_min, x_max, by = step_size)
  
  baselineP <- ggplot(df.tmp, aes(x = .data[[wanted.cols[i]]])) + 
    geom_histogram(aes(y = after_stat(count) / sum(after_stat(count))), 
                   bins = 30, 
                   fill = "#4E79A7", 
                   color = "white") +
    geom_vline(xintercept = median(x_vals, na.rm = TRUE),  
               color = "#E15759", 
               linetype = "dashed", 
               linewidth = 1) +
    scale_y_continuous(labels = scales::percent_format(), name = "Frequency") +
    scale_x_continuous(name = wanted.cols[i], breaks = breaks_seq) +
    labs(title = paste0("Baseline ", wanted.cols[i]),
         subtitle = paste(nrow(df.tmp), "Participants | Median =", 
                          round(median(x_vals, na.rm = TRUE), 3))) +
    theme_minimal(base_size = 14) +
    theme(panel.grid.minor = element_blank())
  
  plots[[length(plots) + 1]] <- baselineP
  plot_count <- plot_count + 1
  
  # Process for recent
  df.tmp <- combined_recent_bp[combined_recent_bp$variable == wanted.cols[i],]
  
  # Extract the column as numeric
  x_vals <- as.numeric(df.tmp[[wanted.cols[i]]])
  x_min <- min(x_vals, na.rm = TRUE)
  x_max <- max(x_vals, na.rm = TRUE)
  step_size <- (x_max - x_min) / 5  # Adjust for number of breaks
  breaks_seq <- seq(x_min, x_max, by = step_size)
  
  recentP <- ggplot(df.tmp, aes(x = .data[[wanted.cols[i]]])) + 
    geom_histogram(aes(y = after_stat(count) / sum(after_stat(count))), 
                   bins = 30, 
                   fill = "#4E79A7", 
                   color = "white") +
    geom_vline(xintercept = median(x_vals, na.rm = TRUE),  
               color = "#E15759", 
               linetype = "dashed", 
               linewidth = 1) +
    scale_y_continuous(labels = scales::percent_format(), name = "Frequency") +
    scale_x_continuous(name = wanted.cols[i], breaks = breaks_seq) +
    labs(title = paste0("Recent ", wanted.cols[i]),
         subtitle = paste(nrow(df.tmp), "Participants | Median =", 
                          round(median(x_vals, na.rm = TRUE), 3))) +
    theme_minimal(base_size = 14) +
    theme(panel.grid.minor = element_blank())
  
  plots[[length(plots) + 1]] <- recentP
  plot_count <- plot_count + 1
  
  # Arrange and print plots in a 2x2 grid every 4 plots
  if (plot_count %% 4 == 0) {
    grid.arrange(grobs = plots, ncol = 4, nrow = 1)  # 4 columns, 1 row
    plots <- list()  # Clear the list for the next set
    plot_count <- 0  # Reset the counter
  }
}

# Print remaining plots if any
if (length(plots) > 0) {
  grid.arrange(grobs = plots, ncol = 4, nrow = ceiling(length(plots) / 4))
}

# Close the PDF device
dev.off()


##################
## echo machine ##
##################
milli <- as.data.frame(read_sas("Z:/SJShare/SJCOMMON/ECC/SJLife/SJLIFE Data Freeze/2 Final Data SJLIFE/20200430/Clinical Data/echo_machine.sas7bdat"))
milli <- milli[milli$sjlid %in% SJLIFE$V2,]
colnames(milli)[colnames(milli) =="studydatetime"] <- "labdt"
wanted.cols <- colnames(milli)[-c(1:7)]


# Initialize lists to store results
base_bp_list <- list()
recent_bp_list <- list()

for (i in seq_along(wanted.cols)) {
  
  print(paste0("Doing: ", wanted.cols[i]))
  
  # Select necessary columns
  milli2 <- milli %>%
    select(sjlid, labdt, !!sym(wanted.cols[i])) %>%
    mutate(!!sym(wanted.cols[i]) := as.numeric(na_if(as.character(!!sym(wanted.cols[i])), ""))) %>%  # Convert empty strings to NA, then to numeric
    filter(!is.na(!!sym(wanted.cols[i])) & !!sym(wanted.cols[i]) != "")  # filter out any missing or "" invalid values
  
  if (nrow(milli2) == 0 || ncol(milli2) < 3 || length(milli2[,3]) == 0| length(unique(milli2[[3]])) <=3) {
    print(paste("Skipping", wanted.cols[i], "as it has no data"))
    next  # Skip this iteration and move to the next column
  }
  
  # Find the baseline (oldest) for each individual
  base_bp <- milli2 %>%
    group_by(sjlid) %>%
    slice_min(labdt, with_ties = FALSE) %>%  # this picks the row with the minimum (oldest) labdt for each individual
    ungroup()  # ungroup after operation to make sure you don't accidentally retain groupings
  
  # Store the result in the list
  base_bp_list[[wanted.cols[i]]] <- base_bp
  
  
  recent_bp <- milli2 %>%
    group_by(sjlid) %>%
    slice_max(labdt, with_ties = FALSE) %>%  # This picks the row with the maximum (most recent) labdt for each individual
    ungroup()
  
  # Store the result in the list
  recent_bp_list[[wanted.cols[i]]] <- recent_bp
  
}


# Combine baseline data frames
combined_base_bp <- as.data.frame(bind_rows(base_bp_list, .id = "variable"))
combined_base_bp.echo.machine <- combined_base_bp
# Combine recent data frames
combined_recent_bp <- as.data.frame(bind_rows(recent_bp_list, .id = "variable"))
combined_recent_bp.echo.machine <- combined_recent_bp

## Plot
# Load necessary libraries
library(ggplot2)
library(scales)
library(gridExtra)

# Create a PDF to save the plots
pdf(paste0("Echo_machine_Distributions_of_", length(unique(milli$sjlid)), "_SJLIFE.pdf"), height = 10, width = 16)  # Adjust width for 4-column layout

# Initialize a list to hold the plots
plots <- list()
plot_count <- 0  # Counter for plots

for (i in seq_along(wanted.cols)) {
  print(paste0("Doing: ", wanted.cols[i]))
  # Process for baseline
  df.tmp <- combined_base_bp[combined_base_bp$variable == wanted.cols[i],]
  
  if (all(is.na(df.tmp[[wanted.cols[i]]]))) {
    print(paste("Skipping", wanted.cols[i], "as all values are NA"))
    next  # Skip this iteration
  }
  
  # Extract the column as numeric
  x_vals <- as.numeric(df.tmp[[wanted.cols[i]]])
  x_min <- min(x_vals, na.rm = TRUE)
  x_max <- max(x_vals, na.rm = TRUE)
  step_size <- (x_max - x_min) / 5  # Adjust for number of breaks
  breaks_seq <- seq(x_min, x_max, by = step_size)
  
  baselineP <- ggplot(df.tmp, aes(x = .data[[wanted.cols[i]]])) + 
    geom_histogram(aes(y = after_stat(count) / sum(after_stat(count))), 
                   bins = 30, 
                   fill = "#4E79A7", 
                   color = "white") +
    geom_vline(xintercept = median(x_vals, na.rm = TRUE),  
               color = "#E15759", 
               linetype = "dashed", 
               linewidth = 1) +
    scale_y_continuous(labels = scales::percent_format(), name = "Frequency") +
    scale_x_continuous(name = wanted.cols[i], breaks = breaks_seq) +
    labs(title = paste0("Baseline ", wanted.cols[i]),
         subtitle = paste(nrow(df.tmp), "Participants | Median =", 
                          round(median(x_vals, na.rm = TRUE), 3))) +
    theme_minimal(base_size = 14) +
    theme(panel.grid.minor = element_blank())
  
  plots[[length(plots) + 1]] <- baselineP
  plot_count <- plot_count + 1
  
  # Process for recent
  df.tmp <- combined_recent_bp[combined_recent_bp$variable == wanted.cols[i],]
  
  # Extract the column as numeric
  x_vals <- as.numeric(df.tmp[[wanted.cols[i]]])
  x_min <- min(x_vals, na.rm = TRUE)
  x_max <- max(x_vals, na.rm = TRUE)
  step_size <- (x_max - x_min) / 5  # Adjust for number of breaks
  breaks_seq <- seq(x_min, x_max, by = step_size)
  
  recentP <- ggplot(df.tmp, aes(x = .data[[wanted.cols[i]]])) + 
    geom_histogram(aes(y = after_stat(count) / sum(after_stat(count))), 
                   bins = 30, 
                   fill = "#4E79A7", 
                   color = "white") +
    geom_vline(xintercept = median(x_vals, na.rm = TRUE),  
               color = "#E15759", 
               linetype = "dashed", 
               linewidth = 1) +
    scale_y_continuous(labels = scales::percent_format(), name = "Frequency") +
    scale_x_continuous(name = wanted.cols[i], breaks = breaks_seq) +
    labs(title = paste0("Recent ", wanted.cols[i]),
         subtitle = paste(nrow(df.tmp), "Participants | Median =", 
                          round(median(x_vals, na.rm = TRUE), 3))) +
    theme_minimal(base_size = 14) +
    theme(panel.grid.minor = element_blank())
  
  plots[[length(plots) + 1]] <- recentP
  plot_count <- plot_count + 1
  
  # Arrange and print plots in a 2x2 grid every 4 plots
  if (plot_count %% 4 == 0) {
    grid.arrange(grobs = plots, ncol = 4, nrow = 1)  # 4 columns, 1 row
    plots <- list()  # Clear the list for the next set
    plot_count <- 0  # Reset the counter
  }
}

# Print remaining plots if any
if (length(plots) > 0) {
  grid.arrange(grobs = plots, ncol = 4, nrow = ceiling(length(plots) / 4))
}

# Close the PDF device
dev.off()

####################
## Function combo ##
####################
milli <- as.data.frame(read_sas("Z:/SJShare/SJCOMMON/ECC/SJLife/SJLIFE Data Freeze/2 Final Data SJLIFE/20200430/Clinical Data/function_combo_basic.sas7bdat"))
milli <- milli[milli$sjlid %in% SJLIFE$V2,]
colnames(milli)[colnames(milli) =="DateVisitStart"] <- "labdt"
wanted.cols <- colnames(milli)[-c(1:4)]


# Initialize lists to store results
base_bp_list <- list()
recent_bp_list <- list()

for (i in seq_along(wanted.cols)) {
  
  print(paste0("Doing: ", wanted.cols[i]))
  
  # Select necessary columns
  milli2 <- milli %>%
    select(sjlid, labdt, !!sym(wanted.cols[i])) %>%
    mutate(!!sym(wanted.cols[i]) := as.numeric(na_if(as.character(!!sym(wanted.cols[i])), ""))) %>%  # Convert empty strings to NA, then to numeric
    filter(!is.na(!!sym(wanted.cols[i])) & !!sym(wanted.cols[i]) != "")  # filter out any missing or "" invalid values
  
  if (nrow(milli2) == 0 || ncol(milli2) < 3 || length(milli2[,3]) == 0| length(unique(milli2[[3]])) <=3) {
    print(paste("Skipping", wanted.cols[i], "as it has no data"))
    next  # Skip this iteration and move to the next column
  }
  
  # Find the baseline (oldest) for each individual
  base_bp <- milli2 %>%
    group_by(sjlid) %>%
    slice_min(labdt, with_ties = FALSE) %>%  # this picks the row with the minimum (oldest) labdt for each individual
    ungroup()  # ungroup after operation to make sure you don't accidentally retain groupings
  
  # Store the result in the list
  base_bp_list[[wanted.cols[i]]] <- base_bp
  
  
  recent_bp <- milli2 %>%
    group_by(sjlid) %>%
    slice_max(labdt, with_ties = FALSE) %>%  # This picks the row with the maximum (most recent) labdt for each individual
    ungroup()
  
  # Store the result in the list
  recent_bp_list[[wanted.cols[i]]] <- recent_bp
  
}


# Combine baseline data frames
combined_base_bp <- as.data.frame(bind_rows(base_bp_list, .id = "variable"))
combined_base_bp.function.combo <- combined_base_bp
# Combine recent data frames
combined_recent_bp <- as.data.frame(bind_rows(recent_bp_list, .id = "variable"))
combined_recent_bp.function.combo <- combined_recent_bp

## Plot
# Load necessary libraries
library(ggplot2)
library(scales)
library(gridExtra)

# Create a PDF to save the plots
pdf(paste0("Function_combo_Distributions_of_", length(unique(milli$sjlid)), "_SJLIFE.pdf"), height = 10, width = 16)  # Adjust width for 4-column layout

# Initialize a list to hold the plots
plots <- list()
plot_count <- 0  # Counter for plots

for (i in seq_along(wanted.cols)) {
  print(paste0("Doing: ", wanted.cols[i]))
  # Process for baseline
  df.tmp <- combined_base_bp[combined_base_bp$variable == wanted.cols[i],]
  
  if (all(is.na(df.tmp[[wanted.cols[i]]]))) {
    print(paste("Skipping", wanted.cols[i], "as all values are NA"))
    next  # Skip this iteration
  }
  
  # Extract the column as numeric
  x_vals <- as.numeric(df.tmp[[wanted.cols[i]]])
  x_min <- min(x_vals, na.rm = TRUE)
  x_max <- max(x_vals, na.rm = TRUE)
  step_size <- (x_max - x_min) / 5  # Adjust for number of breaks
  breaks_seq <- seq(x_min, x_max, by = step_size)
  
  baselineP <- ggplot(df.tmp, aes(x = .data[[wanted.cols[i]]])) + 
    geom_histogram(aes(y = after_stat(count) / sum(after_stat(count))), 
                   bins = 30, 
                   fill = "#4E79A7", 
                   color = "white") +
    geom_vline(xintercept = median(x_vals, na.rm = TRUE),  
               color = "#E15759", 
               linetype = "dashed", 
               linewidth = 1) +
    scale_y_continuous(labels = scales::percent_format(), name = "Frequency") +
    scale_x_continuous(name = wanted.cols[i], breaks = breaks_seq) +
    labs(title = paste0("Baseline ", wanted.cols[i]),
         subtitle = paste(nrow(df.tmp), "Participants | Median =", 
                          round(median(x_vals, na.rm = TRUE), 3))) +
    theme_minimal(base_size = 14) +
    theme(panel.grid.minor = element_blank())
  
  plots[[length(plots) + 1]] <- baselineP
  plot_count <- plot_count + 1
  
  # Process for recent
  df.tmp <- combined_recent_bp[combined_recent_bp$variable == wanted.cols[i],]
  
  # Extract the column as numeric
  x_vals <- as.numeric(df.tmp[[wanted.cols[i]]])
  x_min <- min(x_vals, na.rm = TRUE)
  x_max <- max(x_vals, na.rm = TRUE)
  step_size <- (x_max - x_min) / 5  # Adjust for number of breaks
  breaks_seq <- seq(x_min, x_max, by = step_size)
  
  recentP <- ggplot(df.tmp, aes(x = .data[[wanted.cols[i]]])) + 
    geom_histogram(aes(y = after_stat(count) / sum(after_stat(count))), 
                   bins = 30, 
                   fill = "#4E79A7", 
                   color = "white") +
    geom_vline(xintercept = median(x_vals, na.rm = TRUE),  
               color = "#E15759", 
               linetype = "dashed", 
               linewidth = 1) +
    scale_y_continuous(labels = scales::percent_format(), name = "Frequency") +
    scale_x_continuous(name = wanted.cols[i], breaks = breaks_seq) +
    labs(title = paste0("Recent ", wanted.cols[i]),
         subtitle = paste(nrow(df.tmp), "Participants | Median =", 
                          round(median(x_vals, na.rm = TRUE), 3))) +
    theme_minimal(base_size = 14) +
    theme(panel.grid.minor = element_blank())
  
  plots[[length(plots) + 1]] <- recentP
  plot_count <- plot_count + 1
  
  # Arrange and print plots in a 2x2 grid every 4 plots
  if (plot_count %% 4 == 0) {
    grid.arrange(grobs = plots, ncol = 4, nrow = 1)  # 4 columns, 1 row
    plots <- list()  # Clear the list for the next set
    plot_count <- 0  # Reset the counter
  }
}

# Print remaining plots if any
if (length(plots) > 0) {
  grid.arrange(grobs = plots, ncol = 4, nrow = ceiling(length(plots) / 4))
}

# Close the PDF device
dev.off()

###############
## ECG coded ##
###############
milli <- as.data.frame(read_sas("Z:/SJShare/SJCOMMON/ECC/SJLife/SJLIFE Data Freeze/2 Final Data SJLIFE/20200430/Clinical Data/ecg_coded.sas7bdat"))
milli <- milli[milli$sjlid %in% SJLIFE$V2,]
colnames(milli)[colnames(milli) =="Date"] <- "labdt"
wanted.cols <- colnames(milli)[-c(1:5)]


# Initialize lists to store results
base_bp_list <- list()
recent_bp_list <- list()

for (i in seq_along(wanted.cols)) {
  
  print(paste0("Doing: ", wanted.cols[i]))
  
  # Select necessary columns
  milli2 <- milli %>%
    select(sjlid, labdt, !!sym(wanted.cols[i])) %>%
    mutate(!!sym(wanted.cols[i]) := as.numeric(na_if(as.character(!!sym(wanted.cols[i])), ""))) %>%  # Convert empty strings to NA, then to numeric
    filter(!is.na(!!sym(wanted.cols[i])) & !!sym(wanted.cols[i]) != "")  # filter out any missing or "" invalid values
  
  if (nrow(milli2) == 0 || ncol(milli2) < 3 || length(milli2[,3]) == 0| length(unique(milli2[[3]])) <=3) {
    print(paste("Skipping", wanted.cols[i], "as it has no data"))
    next  # Skip this iteration and move to the next column
  }
  
  # Find the baseline (oldest) for each individual
  base_bp <- milli2 %>%
    group_by(sjlid) %>%
    slice_min(labdt, with_ties = FALSE) %>%  # this picks the row with the minimum (oldest) labdt for each individual
    ungroup()  # ungroup after operation to make sure you don't accidentally retain groupings
  
  # Store the result in the list
  base_bp_list[[wanted.cols[i]]] <- base_bp
  
  
  recent_bp <- milli2 %>%
    group_by(sjlid) %>%
    slice_max(labdt, with_ties = FALSE) %>%  # This picks the row with the maximum (most recent) labdt for each individual
    ungroup()
  
  # Store the result in the list
  recent_bp_list[[wanted.cols[i]]] <- recent_bp
  
}


# Combine baseline data frames
combined_base_bp <- as.data.frame(bind_rows(base_bp_list, .id = "variable"))
combined_base_bp.ECG.coded <- combined_base_bp
# Combine recent data frames
combined_recent_bp <- as.data.frame(bind_rows(recent_bp_list, .id = "variable"))
combined_recent_bp.ECG.coded <- combined_recent_bp

## Plot
# Load necessary libraries
library(ggplot2)
library(scales)
library(gridExtra)

# Create a PDF to save the plots
pdf(paste0("ECG_coded_Distributions_of_", length(unique(milli$sjlid)), "_SJLIFE.pdf"), height = 10, width = 16)  # Adjust width for 4-column layout

# Initialize a list to hold the plots
plots <- list()
plot_count <- 0  # Counter for plots

for (i in seq_along(wanted.cols)) {
  print(paste0("Doing: ", wanted.cols[i]))
  # Process for baseline
  df.tmp <- combined_base_bp[combined_base_bp$variable == wanted.cols[i],]
  
  if (all(is.na(df.tmp[[wanted.cols[i]]]))) {
    print(paste("Skipping", wanted.cols[i], "as all values are NA"))
    next  # Skip this iteration
  }
  
  # Extract the column as numeric
  x_vals <- as.numeric(df.tmp[[wanted.cols[i]]])
  x_min <- min(x_vals, na.rm = TRUE)
  x_max <- max(x_vals, na.rm = TRUE)
  step_size <- (x_max - x_min) / 5  # Adjust for number of breaks
  breaks_seq <- seq(x_min, x_max, by = step_size)
  
  baselineP <- ggplot(df.tmp, aes(x = .data[[wanted.cols[i]]])) + 
    geom_histogram(aes(y = after_stat(count) / sum(after_stat(count))), 
                   bins = 30, 
                   fill = "#4E79A7", 
                   color = "white") +
    geom_vline(xintercept = median(x_vals, na.rm = TRUE),  
               color = "#E15759", 
               linetype = "dashed", 
               linewidth = 1) +
    scale_y_continuous(labels = scales::percent_format(), name = "Frequency") +
    scale_x_continuous(name = wanted.cols[i], breaks = breaks_seq) +
    labs(title = paste0("Baseline ", wanted.cols[i]),
         subtitle = paste(nrow(df.tmp), "Participants | Median =", 
                          round(median(x_vals, na.rm = TRUE), 3))) +
    theme_minimal(base_size = 14) +
    theme(panel.grid.minor = element_blank())
  
  plots[[length(plots) + 1]] <- baselineP
  plot_count <- plot_count + 1
  
  # Process for recent
  df.tmp <- combined_recent_bp[combined_recent_bp$variable == wanted.cols[i],]
  
  # Extract the column as numeric
  x_vals <- as.numeric(df.tmp[[wanted.cols[i]]])
  x_min <- min(x_vals, na.rm = TRUE)
  x_max <- max(x_vals, na.rm = TRUE)
  step_size <- (x_max - x_min) / 5  # Adjust for number of breaks
  breaks_seq <- seq(x_min, x_max, by = step_size)
  
  recentP <- ggplot(df.tmp, aes(x = .data[[wanted.cols[i]]])) + 
    geom_histogram(aes(y = after_stat(count) / sum(after_stat(count))), 
                   bins = 30, 
                   fill = "#4E79A7", 
                   color = "white") +
    geom_vline(xintercept = median(x_vals, na.rm = TRUE),  
               color = "#E15759", 
               linetype = "dashed", 
               linewidth = 1) +
    scale_y_continuous(labels = scales::percent_format(), name = "Frequency") +
    scale_x_continuous(name = wanted.cols[i], breaks = breaks_seq) +
    labs(title = paste0("Recent ", wanted.cols[i]),
         subtitle = paste(nrow(df.tmp), "Participants | Median =", 
                          round(median(x_vals, na.rm = TRUE), 3))) +
    theme_minimal(base_size = 14) +
    theme(panel.grid.minor = element_blank())
  
  plots[[length(plots) + 1]] <- recentP
  plot_count <- plot_count + 1
  
  # Arrange and print plots in a 2x2 grid every 4 plots
  if (plot_count %% 4 == 0) {
    grid.arrange(grobs = plots, ncol = 4, nrow = 1)  # 4 columns, 1 row
    plots <- list()  # Clear the list for the next set
    plot_count <- 0  # Reset the counter
  }
}

# Print remaining plots if any
if (length(plots) > 0) {
  grid.arrange(grobs = plots, ncol = 4, nrow = ceiling(length(plots) / 4))
}

# Close the PDF device
dev.off()


##################################
## Get Summary of all variables ##
##################################

all.data <- list(
  combined_base_bp.milli.measure = combined_base_bp.milli.measure,
  combined_base_bp.milli.dexa = combined_base_bp.milli.dexa,
  combined_base_bp.milli.labs = combined_base_bp.milli.labs,
  combined_base_bp.milli.pulmonary = combined_base_bp.milli.pulmonary,
  combined_base_bp.echo.machine = combined_base_bp.echo.machine,
  combined_base_bp.echo.research = combined_base_bp.echo.research,
  combined_base_bp.function.combo = combined_base_bp.function.combo,
  combined_base_bp.ECG.coded = combined_base_bp.ECG.coded
)

df.summary <- data.frame(
  this.dat.name <- character(),
  variable = character(),
  n_samples = numeric(),
  n_sjlife = numeric(),
  sjlife.percent = numeric(),
  n_sjlife.4381 = numeric(),
  sjlife.4381.percent = numeric(),
  stringsAsFactors = FALSE
)

for (j in 1:length(all.data)){

this.dat.name <- names(all.data)[j]  # Get the dataset name as a string
this.dat <- all.data[[j]]  # Retrieve the actual dataset
  
print(paste("Processing:", this.dat.name))

wanted.cols <- unique(this.dat$variable)
for (i in 1:length(wanted.cols)){
variable = wanted.cols[i]
df.tmp <- this.dat[this.dat$variable == wanted.cols[i],]
n_samples = nrow(df.tmp)
n_sjlife = sum(SJLIFE$V2 %in% df.tmp$sjlid)
sjlife.percent = round(n_sjlife/nrow(SJLIFE)*100, 2)
n_sjlife.4381 = sum(sjlife_4381 %in% df.tmp$sjlid)
sjlife.4381.percent = round(n_sjlife.4381/length(sjlife_4381)*100, 2)
df.summary.tmp <- data.frame(this.dat.name, variable, n_samples, n_sjlife, sjlife.percent, n_sjlife.4381, sjlife.4381.percent)
df.summary <- rbind.data.frame(df.summary, df.summary.tmp)
}
} 

df.summary.sorted <- df.summary %>%
  mutate(sjlife.4381.percent = as.numeric(sub("%", "", sjlife.4381.percent))) %>%  # Convert to numeric
  group_by(this.dat.name) %>%  # Keep groups intact
  arrange(this.dat.name, desc(sjlife.4381.percent), .by_group = TRUE) %>%  # Sort within groups
  ungroup()


write.table(df.summary.sorted, "All_continous_variables.txt", row.names = F, col.names = T, quote = F, sep = "\t")
