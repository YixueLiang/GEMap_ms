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
df_crispr <- df_crispr[rownames(df_crispr) != "ACH-002198", ]

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

gene_data <- df_crispr[, c("CCNE1", "PKMYT1"), drop = FALSE]
gene_data$Cellline <- rownames(gene_data)

# Calculate correlation for each tissue group
cor_results <- list()
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
# ggplot(gene_data, aes(x = CCNE1, y = PKMYT1, color = Tissue)) +
#   geom_point(aes(fill = Tissue), shape = 21, size = 3, color = "black", alpha = 0.7) +  # Colored points with gray border
#   geom_smooth(method = "lm", formula = y ~ x, se = TRUE, color = "#525252", fill = "grey80", linetype = "solid") +  # Add regression line with shaded confidence interval
#   scale_color_manual(values = colors) +  # Use the tissue-based colors for borders
#   scale_fill_manual(values = colors) +  # Use the tissue-based colors for fill
#   theme_minimal() +
#   labs(x = "CCNE1 Dependency", y = "PKMYT1 Dependency", 
#        title = paste("Scatter Plot of", gene[1], "vs", gene[2])) +
#   theme(legend.position = "right") +
#   annotate("text", x = x_pos, y = y_pos, label = paste("cor =", round(cor(gene_data$CCNE1, gene_data$PKMYT1, use = "complete.obs"), 2)), 
#            hjust = -0.1, vjust = 1.1, color = "black", size = 4)  # Display correlation value in upper-left corner

ggplot(gene_data, aes(x = CCNE1, y = PKMYT1, color = Tissue)) +
  geom_point(aes(fill = Tissue), shape = 21, size = 3, color = "black", alpha = 0.7) +  # 彩色点，黑色边框
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE, color = "#525252", fill = "grey80", linetype = "solid") +  # 添加回归线和置信区间
  scale_color_manual(values = colors) +  # 使用基于组织的颜色（边框）
  scale_fill_manual(values = colors) +  # 使用基于组织的颜色（填充）
  theme(
    panel.background = element_blank(),  # 设置背景为空白
    panel.border = element_rect(colour = "black", fill = NA, size = 1),  # 添加黑方框
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),  # 标题居中加粗
    legend.position = "right"  # 图例放在右侧
  ) +
  labs(x = "CCNE1 Dependency", y = "PKMYT1 Dependency", 
       title = "Overview") +  # 修改标题为 "Overview"
  annotate("text", x = x_pos, y = y_pos, label = paste("cor =", round(cor(gene_data$CCNE1, gene_data$PKMYT1, use = "complete.obs"), 2)), 
           hjust = -0.1, vjust = 1.1, color = "black", size = 4)  # 在左上角显示相关系数
# ggplot(gene_data, aes(x = CCNE1, y = PKMYT1, color = Tissue)) +
#   geom_point(shape = 21, fill = "white", color = "gray50", size = 3, alpha = 0.7) +  # Gray border with white fill
#   geom_point(aes(fill = Tissue), shape = 21, size = 3, alpha = 0.7) +  # Colored points
#   geom_smooth(method = "lm", formula = y ~ x, se = TRUE, color = "#525252", fill = "grey80", linetype = "solid") +  # Add regression line with shaded confidence interval
#   scale_color_manual(values = tissue_colors) +  # Use the tissue-based colors
#   scale_fill_manual(values = tissue_colors) +  # Ensure the fill matches the tissue-based colors
#   theme_minimal() +
#   labs(x = "CCNE1 Expression", y = "PKMYT1 Expression", 
#        title = paste("Scatter Plot of", gene[1], "vs", gene[2])) +
#   theme(legend.position = "right") +
#   annotate("text", x = -2, y = 0, label = paste("cor =", round(cor(gene_data$CCNE1, gene_data$PKMYT1), 2)), 
#            hjust = -0.1, vjust = 1.1, color = "black", size = 4)  # Display correlation value in upper-left corner


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
