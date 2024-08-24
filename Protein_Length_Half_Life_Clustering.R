# Load necessary libraries
library(tidyverse)
library(cluster)  # For hierarchical clustering
library(ggplot2)

# Load protein Data
protein_data <- read_csv("/Users/chrisfluta/Downloads/data.csv")

# Check for missing data in the raw dataset
colSums(is.na(protein_data))

# Impute missing values with the mean of each column
protein_data_imputed <- protein_data %>%
  mutate(across(where(is.numeric), ~ ifelse(is.na(.), mean(., na.rm = TRUE), .)))

# Identify the column indices containing "half_life"
half_life_cols <- grep("half_life", names(protein_data_imputed))

# Set 1: Subset the data frame using the identified indices for half-lives
set1 <- protein_data_imputed[, half_life_cols]

# Set 2: Protein length
set2 <- protein_data_imputed[, "Length", drop = FALSE]

# Combine the half-lives and length into a single dataset for clustering
clustering_data_raw <- cbind(set1, set2)

# Standardize the data before clustering (important for K-Means)
clustering_data_raw <- scale(clustering_data_raw)

# Step 1: Clustering on Raw Data
# Elbow Method for K-Means
wcss_raw <- vector()
for (i in 1:10) {
  kmeans_result <- kmeans(clustering_data_raw, centers = i, nstart = 25)
  wcss_raw[i] <- kmeans_result$tot.withinss  # Total within-cluster sum of squares
}

# Plot the WCSS against the number of clusters
plot(1:10, wcss_raw, type = "b", pch = 19, frame = FALSE, 
     xlab = "Number of Clusters", 
     ylab = "Total Within-Cluster Sum of Squares (WCSS)",
     main = "Elbow Method for Optimal Clusters (Raw Data)")

# K-Means Clustering on Raw Data
set.seed(42)  # For reproducibility
kmeans_results_raw <- kmeans(clustering_data_raw, centers = 3, nstart = 25)

# Hierarchical Clustering on Raw Data
dist_matrix_raw <- dist(clustering_data_raw)  # Compute the distance matrix
hclust_results_raw <- hclust(dist_matrix_raw, method = "ward.D2")  # Perform hierarchical clustering

# Plot the dendrogram for Raw Data
plot(hclust_results_raw, main = "Hierarchical Clustering Dendrogram (Raw Data)", xlab = "", sub = "", cex = 0.6)

# Cut the dendrogram to form clusters
hclust_clusters_raw <- cutree(hclust_results_raw, k = 3)  # Adjust 'k' to match the number of clusters

# Add cluster assignments to the original data
protein_data_imputed$kmeans_cluster_raw <- kmeans_results_raw$cluster
protein_data_imputed$hclust_cluster_raw <- hclust_clusters_raw

# Step 2: Clustering on PCA Data
# Perform PCA on Raw Data
pca_results_raw <- prcomp(clustering_data_raw, scale = TRUE)
pca_data_raw <- pca_results_raw$x[, 1:2]  # Use first two principal components

# Elbow Method for K-Means on PCA Data
wcss_pca <- vector()
for (i in 1:10) {
  kmeans_result_pca <- kmeans(pca_data_raw, centers = i, nstart = 25)
  wcss_pca[i] <- kmeans_result_pca$tot.withinss  # Total within-cluster sum of squares
}

# Plot the WCSS against the number of clusters
plot(1:10, wcss_pca, type = "b", pch = 19, frame = FALSE, 
     xlab = "Number of Clusters", 
     ylab = "Total Within-Cluster Sum of Squares (WCSS)",
     main = "Elbow Method for Optimal Clusters (PCA Data)")

# K-Means Clustering on PCA Data
set.seed(42)  # For reproducibility
kmeans_results_pca <- kmeans(pca_data_raw, centers = 3, nstart = 25)

# Hierarchical Clustering on PCA Data
dist_matrix_pca <- dist(pca_data_raw)  # Compute the distance matrix
hclust_results_pca <- hclust(dist_matrix_pca, method = "ward.D2")  # Perform hierarchical clustering

# Plot the dendrogram for PCA Data
plot(hclust_results_pca, main = "Hierarchical Clustering Dendrogram (PCA Data)", xlab = "", sub = "", cex = 0.6)

# Cut the dendrogram to form clusters
hclust_clusters_pca <- cutree(hclust_results_pca, k = 3)  # Adjust 'k' to match the number of clusters

# Add cluster assignments to the PCA data
protein_data_imputed$kmeans_cluster_pca <- kmeans_results_pca$cluster
protein_data_imputed$hclust_cluster_pca <- hclust_clusters_pca

# Step 3: Compare Clustering Results
# Visualize K-Means Clustering Results for Raw Data
ggplot(protein_data_imputed, aes(x = Length, y = `Bcells replicate 1 half_life`, color = factor(kmeans_cluster_raw))) +
  geom_point() +
  labs(title = "K-Means Clustering (3 Clusters) - Raw Data", x = "Protein Length", y = "Bcells replicate 1 half_life") +
  theme_minimal()

# Visualize K-Means Clustering Results for PCA Data
ggplot(protein_data_imputed, aes(x = pca_data_raw[,1], y = pca_data_raw[,2], color = factor(kmeans_cluster_pca))) +
  geom_point() +
  labs(title = "K-Means Clustering (3 Clusters) - PCA Data", x = "PC1", y = "PC2") +
  theme_minimal()

# Visualize Hierarchical Clustering Results for Raw Data
ggplot(protein_data_imputed, aes(x = Length, y = `Bcells replicate 1 half_life`, color = factor(hclust_cluster_raw))) +
  geom_point() +
  labs(title = "Hierarchical Clustering (3 Clusters) - Raw Data", x = "Protein Length", y = "Bcells replicate 1 half_life") +
  theme_minimal()

# Visualize Hierarchical Clustering Results for PCA Data
ggplot(protein_data_imputed, aes(x = pca_data_raw[,1], y = pca_data_raw[,2], color = factor(hclust_cluster_pca))) +
  geom_point() +
  labs(title = "Hierarchical Clustering (3 Clusters) - PCA Data", x = "PC1", y = "PC2") +
  theme_minimal()


# Compare the Cluster Assignments
table(protein_data_imputed$kmeans_cluster_raw)
table(protein_data_imputed$kmeans_cluster_pca)
table(protein_data_imputed$hclust_cluster_raw)
table(protein_data_imputed$hclust_cluster_pca)
