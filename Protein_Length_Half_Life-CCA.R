# Load necessary libraries
library(tidyverse)
library(CCA)

# Load protein Data
protein_data <- read_csv("/Users/chrisfluta/Downloads/data.csv")

# Check for missing data in the raw dataset
colSums(is.na(protein_data))

# Impute missing values with the mean of each column
protein_data_imputed <- protein_data %>%
  mutate(across(where(is.numeric), ~ ifelse(is.na(.), mean(., na.rm = TRUE), .)))

# Identify the column indices containing "half_life"
half_life_cols <- grep("half_life", names(protein_data_imputed))

# Set 1: Subset the data frame using the identified indices
set1 <- protein_data_imputed[, half_life_cols]

# Set 2: Protein length (using base R)
set2 <- protein_data_imputed[, "Length", drop = FALSE]

# Perform CCA on the raw data
cca_results <- cancor(set1, set2)

# Summarize the CCA results
summary(cca_results)

# Canonical correlations
cca_results$cor

# Canonical weights for Set 1 (half-life)
cca_results$xcoef

# Canonical weights for Set 2 (length)
cca_results$ycoef
