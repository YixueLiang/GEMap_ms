#清空数据
rm(list = ls())
# Load necessary libraries
library(ggplot2)
library(ggpubr)
library(RColorBrewer)
library(dplyr)

# Load the data
load("DEMap_Data/CRISPR_DATA.Rdata")

# Data preprocessing
COM <- intersect(rownames(df_crispr), rownames(df_exp))
df_crispr <- df_crispr[COM, ]
df_exp <- df_exp[COM, ]

commoncell <- intersect(rownames(df_crispr), df_cell$ModelID)
rownames(df_cell) <- df_cell$ModelID
df_cell <- df_cell[rownames(df_crispr), ]
rownames(df_crispr) <- df_cell$StrippedCellLineName
Tissue <- unique(df_cell$OncotreeLineage)
Tissue <- Tissue[Tissue != "Other"]

gene <- c("CCNE1", "PKMYT1")

# Create tissue groups
Tissue_groups <- list()
Tissue_groups <- split(df_cell$StrippedCellLineName, df_cell$OncotreeLineage)

# Remove the "Other" tissue group from Tissue_groups
Tissue_groups[["Other"]] <- NULL

# df_crispr <- df_crispr[rownames(df_crispr) != "SCH", ]

# Define the custom color palette for tissues
colors <- rev(colorRampPalette(rev(brewer.pal(9, "RdYlBu")))(27))

# Extract gene expression data for CCNE1 from df_crispr and PKMYT1 from df_exp
gene_data <- data.frame(
  Cellline = rownames(df_crispr),
  CCNE1 = df_crispr[, "CCNE1"],
  PKMYT1 = df_exp[, "PKMYT1"]
)
rownames(gene_data) <- gene_data$Cellline

# Calculate correlation for each tissue group
# Initialize an empty list to store the correlation results
cor_results <- list()

# Loop through each tissue group
for(tissue in names(Tissue_groups)) {
  tissue_cells <- Tissue_groups[[tissue]]
  tissue_data <- gene_data[tissue_cells, gene]
  cor_results[[tissue]] <- cor(tissue_data[,1], tissue_data[,2])
}

# Create a data frame with correlation results and tissue types
cor_df <- data.frame(
  Tissue = names(cor_results),
  Gene1 = gene[1],
  Gene2 = gene[2],
  Correlation = unlist(cor_results)
)

cor_df <- cor_df[order(cor_df$Correlation), ]
cor_df$color <- colors

gene_data$Tissue <- df_cell$OncotreeLineage[match(rownames(gene_data), df_cell$StrippedCellLineName)]

# Assign colors to tissue types
gene_data$Color <- cor_df$color[match(gene_data$Tissue, cor_df$Tissue)]

CCNE1_expr <- as.numeric(gene_data[, "CCNE1"])
PKMYT1_expr <- as.numeric(gene_data[, "PKMYT1"])

# Perform Spearman correlation test
gene_corr <- cor.test(CCNE1_expr, PKMYT1_expr)

cor <- gene_corr$estimate
pValue <- gene_corr$p.value

# Get the range of x and y to position the annotation dynamically
x_range <- range(gene_data$CCNE1, na.rm = TRUE)
y_range <- range(gene_data$PKMYT1, na.rm = TRUE)

# Set annotation position dynamically (e.g., place it at the top-left)
x_pos <- x_range[1] + 0.05 * diff(x_range)  # 5% into the x-axis range
y_pos <- y_range[2] - 0.05 * diff(y_range)  # 5% down from the top y-axis value

# Create scatter plot with tissue-based colors
ggplot(gene_data, aes(x = CCNE1, y = PKMYT1, color = Tissue)) +
  geom_point(aes(fill = Tissue), shape = 21, size = 3, color = "black", alpha = 0.7) +  # Colored points with gray border
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE, color = "#525252", fill = "grey80", linetype = "solid") +  # Add regression line with shaded confidence interval
  scale_color_manual(values = colors) +  # Use the tissue-based colors for borders
  scale_fill_manual(values = colors) +  # Use the tissue-based colors for fill
  theme_minimal() +
  labs(x = "CCNE1 Dependency", y = "PKMYT1 Expression", 
       title = paste("Scatter Plot of", gene[1], "vs", gene[2])) +
  theme(legend.position = "right") +
  annotate("text", x = x_pos, y = y_pos, label = paste("cor =", round(cor(gene_data$CCNE1, gene_data$PKMYT1, use = "complete.obs"), 2)), 
           hjust = -0.1, vjust = 1.1, color = "black", size = 4)  # Display correlation value in upper-left corner


# tissue ------------------------------------------------------------------


# Extract gene expression data for CCNE1 and PKMYT1 from df_crispr
gene_data_tissue <- gene_data[Tissue_groups[["Lung"]], , drop = FALSE]

# Extract numeric vectors for the gene expressions
CCNE1_expr <- as.numeric(gene_data_tissue[, "CCNE1"])
PKMYT1_expr <- as.numeric(gene_data_tissue[, "PKMYT1"])

# Perform Spearman correlation test
gene_corr <- cor.test(CCNE1_expr, PKMYT1_expr)

cor <- gene_corr$estimate
pValue <- gene_corr$p.value

# Get the range of x and y to position the annotation dynamically
x_range <- range(gene_data_tissue$CCNE1, na.rm = TRUE)
y_range <- range(gene_data_tissue$PKMYT1, na.rm = TRUE)

# Set annotation position dynamically (e.g., place it at the top-left)
x_pos <- x_range[1] + 0.05 * diff(x_range)  # 5% into the x-axis range
y_pos <- y_range[2] - 0.05 * diff(y_range)  # 5% down from the top y-axis value

ggplot(gene_data_tissue, aes(x = CCNE1, y = PKMYT1, color = Color)) +  # Use gene_data_tissue$Color for point colors
  geom_point(aes(fill = Color), shape = 21, size = 3, color = "black", alpha = 0.7) +  # Colored points with black border
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE, color = "#525252", fill = "grey80", linetype = "solid") +  # Add regression line with shaded confidence interval
  scale_color_identity() +  # Use the exact colors from gene_data_tissue$Color
  scale_fill_identity() +  # Use the exact colors from gene_data_tissue$Color
  theme_minimal() +
  labs(x = "CCNE1 Expression", y = "PKMYT1 Expression", 
       title = paste("Scatter Plot of", gene[1], "vs", gene[2])) +
  theme(legend.position = "right") +
  annotate("text", x = x_pos, y = y_pos, label = paste("cor =", round(cor(gene_data_tissue$CCNE1, gene_data_tissue$PKMYT1, use = "complete.obs"), 2)), 
           hjust = -0.1, vjust = 1.1, color = "black", size = 4)  # Display correlation value in upper-left corner
