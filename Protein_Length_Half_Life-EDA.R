#install.packages("tidyverse")
library(tidyverse)
library(ggplot2)
library(tidyr)
protein_data <- read_csv("/Users/chrisfluta/Downloads/data.csv")

# Peek at data to understand general shape and variables
head(protein_data)
dim(protein_data)
str(protein_data)
summary(protein_data)
#colNames(protein_data)

# Omit NA data and see how many rows are left
colSums(is.na(protein_data))
protein_data <- na.omit(protein_data)
dim(protein_data)

# Merge replicates for each cell type by mean
# Assuming your data is in a data frame called 'protein_data'

# Calculate the mean for each cell type's half_life and R_sq
protein_data <- protein_data %>%
  mutate(
    Bcells_mean_half_life = rowMeans(select(., contains("Bcells") & contains("half_life")), na.rm = TRUE),
    Bcells_mean_R_sq = rowMeans(select(., contains("Bcells") & contains("R_sq")), na.rm = TRUE),
    
    NK_cells_mean_half_life = rowMeans(select(., contains("NK cells") & contains("half_life")), na.rm = TRUE),
    NK_cells_mean_R_sq = rowMeans(select(., contains("NK cells") & contains("R_sq")), na.rm = TRUE),
    
    Hepatocytes_mean_half_life = rowMeans(select(., contains("Hepatocytes") & contains("half_life")), na.rm = TRUE),
    Hepatocytes_mean_R_sq = rowMeans(select(., contains("Hepatocytes") & contains("R_sq")), na.rm = TRUE),
    
    Monocytes_mean_half_life = rowMeans(select(., contains("Monocytes") & contains("half_life")), na.rm = TRUE),
    Monocytes_mean_R_sq = rowMeans(select(., contains("Monocytes") & contains("R_sq")), na.rm = TRUE),
    
    Mouse_Neurons_mean_half_life = rowMeans(select(., contains("Mouse Neurons") & contains("half_life")), na.rm = TRUE),
    Mouse_Neurons_mean_R_sq = rowMeans(select(., contains("Mouse Neurons") & contains("R_sq")), na.rm = TRUE)
  )

# Drop the original replicate columns if they are no longer needed
protein_data <- protein_data %>%
  select(-contains("replicate"))

# View the updated data frame
head(protein_data)


# Histograms of Protein length
hist(protein_data$Length, main="Protein Length Distribution", xlab="Length (Amino Acids)")

# Plot half life for each cell type
ggplot(protein_data_long, aes(x = Mean_Half_Life, fill = Cell_Type)) + 
  geom_histogram(binwidth = 0.5, alpha = 0.6) +
  labs(title = "Distribution of Mean Half-Life Across Cell Types",
       x = "Mean Half-Life (Hours)", 
       y = "Count") +
  theme_minimal() +
  facet_wrap(~ Cell_Type, scales = "free_y") +
  scale_fill_brewer(palette = "Set1")


# Create boxplots for each of the key variables to identify potential outliers.
ggplot(protein_data, aes(y = Length)) + 
  geom_boxplot() + 
  labs(title = "Boxplot of Protein Length", y = "Length (Amino Acids)")

ggplot(protein_data_long, aes(y = Mean_Half_Life, fill = Cell_Type)) + 
  geom_boxplot() + 
  labs(title = "Boxplot of Mean Half-Life by Cell Type", y = "Mean Half-Life (Hours)") +
  facet_wrap(~ Cell_Type)

# Inspect variability in mean half-lives for each cell type
ggplot(protein_data_long, aes(x = Mean_Half_Life, fill = Cell_Type)) + 
  geom_density(alpha = 0.5) +
  labs(title = "Density Plot of Mean Half-Life Across Cell Types",
       x = "Mean Half-Life (Hours)", 
       y = "Density") +
  theme_minimal() +
  scale_fill_brewer(palette = "Set1")
