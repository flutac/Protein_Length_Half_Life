# Load necessary libraries
library(tidyverse)
library(wordspace)
library(ggplot2)

# Load protein Data
protein_data <- read_csv("/Users/chrisfluta/Downloads/data.csv")

# Clean Protein data from EDA (removing NAs)
colSums(is.na(protein_data))
protein_data <- na.omit(protein_data)
dim(protein_data)

# Calculate the mean for each cell type's half_life
cleaned_protein_data <- protein_data %>%
  mutate(
    Bcells_mean_half_life = rowMeans(select(., contains("Bcells") & contains("half_life")), na.rm = TRUE),
    NK_cells_mean_half_life = rowMeans(select(., contains("NK cells") & contains("half_life")), na.rm = TRUE),
    Hepatocytes_mean_half_life = rowMeans(select(., contains("Hepatocytes") & contains("half_life")), na.rm = TRUE),
    Monocytes_mean_half_life = rowMeans(select(., contains("Monocytes") & contains("half_life")), na.rm = TRUE),
    Mouse_Neurons_mean_half_life = rowMeans(select(., contains("Mouse Neurons") & contains("half_life")), na.rm = TRUE)
  )

# Drop the original replicate columns and R_sq columns if they are no longer needed
cleaned_protein_data <- cleaned_protein_data %>%
  select(-contains("replicate"), -contains("R_sq"))

# Select relevant columns for PCA (excluding R_sq)
pca_data <- select(cleaned_protein_data, contains("mean_half_life"), Length)

# Normalize the data by the maximum value over each column
pca_data_matrix <- data.matrix(pca_data)
pca_data_norm <- colNorms(pca_data_matrix, method = "maximum", p = 2)
normed_pca_data <- pca_data_matrix / pca_data_norm

# Compute the covariance matrix of the normalized data
pca_cov <- cov(normed_pca_data)

# Perform PCA using the covariance matrix
pca_results <- princomp(pca_cov)
summary(pca_results, loadings=TRUE)

# Plot the Scree Plot
plot(pca_results$sdev^2, xlab="Components", ylab="Captured variance",
     type="l", main="Scree plot for Cleaned Protein Data")

# Transform the original data using the first two principal components
# Keep the first two principal components
pca_keep <- pca_results$loadings[,1:2]

# Multiply the loading values by the data
pca_transformed <- pca_data_matrix %*% pca_keep

# Make sure the shape is n_obs by two components
dim(pca_transformed)

# Plot the transformed data
plot(pca_transformed, main="Cleaned Protein Data Principal Components")



#### PCA for uncleaned data (unmerged reps)
# Impute missing values with the mean of each column
protein_data_imputed <- protein_data %>%
  mutate(across(where(is.numeric), ~ ifelse(is.na(.), mean(., na.rm = TRUE), .)))

# Select relevant columns for PCA (excluding R_sq)
pca_data_uncleaned <- select(protein_data_imputed, contains("half_life"), Length)

# Normalize the data by the maximum value over each column
pca_data_matrix_uncleaned <- data.matrix(pca_data_uncleaned)
pca_data_norm_uncleaned <- colNorms(pca_data_matrix_uncleaned, method = "maximum", p = 2)
normed_pca_data_uncleaned <- pca_data_matrix_uncleaned / pca_data_norm_uncleaned

# Compute the covariance matrix of the normalized data
pca_cov_uncleaned <- cov(normed_pca_data_uncleaned)

# Perform PCA using the covariance matrix
pca_results_uncleaned <- princomp(pca_cov_uncleaned)
summary(pca_results_uncleaned, loadings=TRUE)

# Plot the Scree Plot
plot(pca_results_uncleaned$sdev^2, xlab="Components", ylab="Captured variance",
     type="l", main="Scree plot for Uncleaned Protein Data")

# Transform the original data using the first two principal components
# Keep the first two principal components
pca_keep_uncleaned <- pca_results_uncleaned$loadings[,1:2]

# Multiply the loading values by the data
pca_transformed_uncleaned <- pca_data_matrix_uncleaned %*% pca_keep_uncleaned

# Make sure the shape is n_obs by two components
dim(pca_transformed_uncleaned)

# Plot the transformed data
plot(pca_transformed_uncleaned, main="Uncleaned Protein Data Principal Components")

