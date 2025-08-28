#清空数据
rm(list = ls())
# Load necessary libraries
library(ggplot2)
library(RColorBrewer)
library(dplyr)

# Load the data
load("DEMap_Data/CRISPR_DATA.Rdata")

# Data preprocessing
commoncell <- intersect(rownames(df_exp), df_cell$ModelID)
rownames(df_cell) <- df_cell$ModelID
df_cell <- df_cell[rownames(df_exp), ]
rownames(df_exp) <- df_cell$StrippedCellLineName
Tissue <- unique(df_cell$OncotreeLineage)

# Create tissue groups
Tissue_groups <- list()
Tissue_groups <- split(df_cell$StrippedCellLineName, df_cell$OncotreeLineage)

# Create crispr_groups by splitting df_exp based on tissue groups
crispr_groups <- list()

for (Tissue in names(Tissue_groups)) {
  crispr_groups[[Tissue]] <- df_exp[Tissue_groups[[Tissue]], ]
}

# Specify the gene for plotting
gene <- "PKMYT1"
gene_values <- list()

# Extract the gene expression values for the specified gene by tissue
for (Tissue in names(crispr_groups)) {
  gene_values[[Tissue]] <- crispr_groups[[Tissue]]["PKMYT1"]
}

# Create an empty list to store vectors for each tissue
crispr_vectors <- list()

# Loop through each tissue to convert data frames to vectors
for (Tissue in names(gene_values)) {
  # Convert the data frame to a matrix
  matrix_df <- as.matrix(gene_values[[Tissue]])
  
  # Flatten the matrix into a vector
  vector_df <- as.vector(matrix_df)
  
  # Store the vector in the list
  crispr_vectors[[Tissue]] <- vector_df
}

# Create an empty data frame to store the data in long format
plotdata <- data.frame(Value = numeric(0), Tissue = character(0), CellName = character(0))

# Combine all vectors into the long format data frame
for (Tissue in names(crispr_vectors)) {
  temp_df <- data.frame(
    Value = crispr_vectors[[Tissue]],
    Tissue = Tissue,
    CellName = rownames(gene_values[[Tissue]])  # Add CellName from row names
  )
  plotdata <- rbind(plotdata, temp_df)
}

# Remove tissues labeled as "other"
plotdata <- plotdata %>% filter(!grepl("other", Tissue, ignore.case = TRUE))

# Reorder Tissue factor levels by the median of Value
plotdata$Tissue <- reorder(plotdata$Tissue, plotdata$Value, FUN = median)

# Generate a color palette
colors <- rev(colorRampPalette(rev(brewer.pal(9, "RdYlBu")))(29))

# Create the boxplot
ggplot(plotdata, aes(x = Tissue, y = Value, fill = Tissue)) +
  geom_boxplot() +
  scale_fill_manual(values = colors) +  # Apply the custom color palette
  theme_minimal() +  # Use a clean, minimal theme
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),  # Rotate x-axis labels
    plot.title = element_text(hjust = 0.5, face = "bold", color = "black")  # Center-align title and set to black
  ) +
  labs(x = "Tissue", y = "Dependency", title = paste0(gene, " — Primary Disease"))


# Sort plotdata by Value from smallest to largest
plotdata <- plotdata[order(plotdata$Value), ]

# Ensure CellName is a factor with levels ordered by Value
plotdata$CellName <- factor(plotdata$CellName, levels = unique(plotdata$CellName[order(plotdata$Value)]))

# Generate the dot plot
ggplot(plotdata, aes(x = CellName, y = Value, color = Tissue)) +
  geom_jitter(width = 0.1, height = 0, alpha = 0.7, size = 1.5) +  # Reduce point size and jitter width
  scale_color_manual(values = colors) +  # Apply the custom color palette
  theme_minimal(base_size = 12) +  # Set a base font size
  theme(
    axis.text.x = element_blank(),  # Remove x-axis labels
    axis.ticks.x = element_blank(),  # Remove x-axis ticks
    panel.grid.major.x = element_blank(),  # Remove major grid lines for x-axis
    panel.grid.minor.x = element_blank(),  # Remove minor grid lines for x-axis
    plot.title = element_text(hjust = 0.5, face = "bold", color = "black")  # Center-align and set title to black
  ) +
  labs(y = "Value", color = "Tissue", title = paste0(gene, " — Primary Disease")) +
  coord_cartesian(clip = "off")  # Ensure all points are within the plotting area

