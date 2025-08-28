#清空数据
rm(list = ls())
# Load necessary libraries
library(ggplot2)
library(pheatmap)
library(dplyr)
library(RColorBrewer)

gene <- c("ADA", "ATXN7L1", "ATXN7L2", "ATXN7L3", "KAT2A", "SGF29", 
           "SUPT7L", "TAF5L", "TAF6L", "SUPT20H", "TAF10", "TAF12", 
           "SUPT3H", "TADA2B", "TADA3", "USP22")

# Load the data
load("DEMap_Data/CRISPR_DATA.Rdata")
mutData <- readRDS(file = "DEMap_Data/Mut_Data.rds")
# Data preprocessing
commoncell <- intersect(rownames(mutData), df_cell$ModelID)
rownames(df_cell) <- df_cell$ModelID
df_cell <- df_cell[rownames(mutData), ]
rownames(mutData) <- df_cell$StrippedCellLineName
Tissue <- unique(df_cell$OncotreeLineage)

# Create tissue groups
Tissue_groups <- list()
Tissue_groups <- split(df_cell$StrippedCellLineName, df_cell$OncotreeLineage)

# Assuming your mutation data (`mutData`) is a data frame with rows as cells and columns as genes

# Step 1: Preprocess data for the bubble plot
# Extract mutations for selected genes
mutGeneData <- mutData[, gene]  # Subset the mutData to include only the selected genes

# Create a binary mutation matrix (1 for mutation, 0 for no mutation)
mutGeneBinary <- mutGeneData
mutGeneBinary[mutGeneBinary != 0] <- 1  # Convert all non-zero values to 1 (mutation present)

# Remove the "Other" tissue group from Tissue_groups
Tissue_groups[["Other"]] <- NULL

# Data preprocessing for mutation frequencies
mutationSummary <- data.frame(
  Tissue = rep(names(Tissue_groups), each = length(gene)),
  Gene = rep(gene, times = length(Tissue_groups)),
  MutationFrequency = numeric(length(gene) * length(Tissue_groups)),
  MutatedCells = numeric(length(gene) * length(Tissue_groups)),
  TotalCells = numeric(length(gene) * length(Tissue_groups))
)

# Loop through each tissue group to calculate mutation frequencies
for (tissue in names(Tissue_groups)) {
  tissue_cells <- Tissue_groups[[tissue]]
  total_cells <- length(tissue_cells)
  
  for (g in gene) {
    # Get mutation data for the tissue and gene
    mutated_cells <- sum(mutGeneBinary[tissue_cells, g] == 1)
    
    # Calculate mutation frequency as the ratio of mutated cells to total cells
    mutation_summary_index <- which(mutationSummary$Tissue == tissue & mutationSummary$Gene == g)
    mutationSummary$MutatedCells[mutation_summary_index] <- mutated_cells
    mutationSummary$TotalCells[mutation_summary_index] <- total_cells
    mutationSummary$MutationFrequency[mutation_summary_index] <- mutated_cells / total_cells
  }
}

# Filter out rows where MutationFrequency is zero
mutationSummary <- mutationSummary[mutationSummary$MutationFrequency > 0, ]

# Create the Bubble Plot with filtered data
ggplot(mutationSummary, aes(x = Gene, y = Tissue, size = MutationFrequency, color = MutationFrequency)) +
  geom_point(alpha = 0.7) +
  scale_size_continuous(range = c(1, 10)) +  # Adjust bubble size range
  scale_color_gradient(low = "#7EB5D5", high = "#F57948") +  # Color gradient for frequency
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Gene", y = "Tissue", title = "Mutation Frequency by Gene and Tissue")


# # Bubble plot visualization using ggplot
# ggplot(mutationSummary, aes(x = Gene, y = Tissue, size = MutationFrequency, color = MutationFrequency)) +
#   geom_point(alpha = 0.6) +
#   scale_size_continuous(name = "Mutation Frequency", range = c(1, 10)) +  # Adjust size range
#   scale_color_viridis_c(name = "Mutation Frequency") +  # Color based on mutation frequency
#   theme_minimal() +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1)) +  # Rotate gene labels for better readability
#   labs(title = "Mutation Frequency of Selected Genes by Tissue Type",
#        x = "Gene",
#        y = "Tissue") +
#   theme(axis.title = element_text(size = 12), axis.text = element_text(size = 10))  # Adjust text size



















# # Step 3: Preprocess Data for Heatmap
# # Create a heatmap matrix (genes as rows, tissues as columns)
# heatmapData <- matrix(0, nrow = length(gene), ncol = length(Tissue_groups), 
#                       dimnames = list(gene, names(Tissue_groups)))
# 
# # Fill heatmap data with binary mutation presence (0/1)
# for (tissue in names(Tissue_groups)) {
#   for (g in gene) {
#     heatmapData[g, tissue] <- mean(mutGeneBinary[Tissue_groups[[tissue]], g])  # Mutation presence as 0/1
#   }
# }
# 
# # Step 4: Create Heatmap with Genes on x-axis and Tissues on y-axis
# pheatmap(heatmapData, 
#          cluster_rows = TRUE, cluster_cols = TRUE, 
#          color = colorRampPalette(brewer.pal(9, "RdYlBu"))(100), 
#          display_numbers = TRUE, 
#          main = "Heatmap of Gene Mutations Across Tissues (Excluding 'Other')",
#          fontsize_row = 8, fontsize_col = 8, 
#          angle_col = 45)