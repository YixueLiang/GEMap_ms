#清空数据
rm(list = ls())
# Load necessary libraries
library(ggplot2)
library(RColorBrewer)
library(dplyr)

# Load the data
load("DEMap_Data/GeCKO_DATA.Rdata")
load("DEMap_Data/CRISPR_sanger_avana_gecko.Rdata")

# Data preprocessing
# rownames(df_GeCKO_cell) <- df_GeCKO_cell$CellName
commoncell <- intersect(rownames(df_GeCKO), df_cell$StrippedCellLineName)
df_cell <- df_cell[commoncell, ]
Tissue <- unique(df_cell$OncotreeLineage)


# Create tissue groups
Tissue_groups <- list()
Tissue_groups <- split(df_cell$StrippedCellLineName, df_cell$OncotreeLineage)
a <- Tissue_groups[["Lung"]]
Tissue_groups[["Lung"]] <- c(a, "COLO699")

# Create GeCKO_groups by splitting df_GeCKO based on tissue groups
GeCKO_groups <- list()

for (Tissue in names(Tissue_groups)) {
  GeCKO_groups[[Tissue]] <- df_GeCKO[Tissue_groups[[Tissue]], ]
}

# Specify the gene for plotting
gene <- "CDK4"
gene_values <- list()

# Extract the gene expression values for the specified gene by tissue
for (Tissue in names(GeCKO_groups)) {
  gene_values[[Tissue]] <- GeCKO_groups[[Tissue]]["CDK4"]
}

# Create an empty list to store vectors for each tissue
GeCKO_vectors <- list()

# Loop through each tissue to convert data frames to vectors
for (Tissue in names(gene_values)) {
  # Convert the data frame to a matrix
  matrix_df <- as.matrix(gene_values[[Tissue]])
  
  # Flatten the matrix into a vector
  vector_df <- as.vector(matrix_df)
  
  # Store the vector in the list
  GeCKO_vectors[[Tissue]] <- vector_df
}

# Create an empty data frame to store the data in long format
plotdata <- data.frame(Value = numeric(0), Tissue = character(0), CellName = character(0))

# Combine all vectors into the long format data frame
for (Tissue in names(GeCKO_vectors)) {
  temp_df <- data.frame(
    Value = GeCKO_vectors[[Tissue]],
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
# colors <- rev(colorRampPalette(rev(brewer.pal(9, "RdYlBu")))(27))
colors <- colorRampPalette(rev(brewer.pal(9, "RdYlBu")))(11)

# Create the boxplot
# ggplot(plotdata, aes(x = Tissue, y = Value, fill = Tissue)) +
#   geom_boxplot() +
#   scale_fill_manual(values = colors) +  # 应用自定义颜色调色板
#   theme(
#     panel.background = element_blank(),  # 设置背景为空白
#     panel.border = element_rect(colour = "black", fill = NA, size = 1),  # 添加黑方框
#     axis.text.x = element_text(angle = 45, hjust = 1),  # 旋转 x 轴标签
#     axis.text.y = element_text(),  # 恢复 y 轴标签不加粗
#     plot.title = element_text(hjust = 0.5, face = "bold", color = "black")  # 居中标题并设置为黑色
#   ) +
#   labs(x = "Tissue", y = "Essentiality Score", title = paste0(gene, " — Primary Disease"))

ggplot(plotdata, aes(x = Tissue, y = Value, fill = Tissue)) +
  geom_boxplot() +
  scale_fill_manual(values = colors) +  # 应用自定义颜色调色板
  theme(
    panel.background = element_blank(),  # 设置背景为空白
    panel.border = element_rect(colour = "black", fill = NA, size = 1),  # 添加黑方框
    axis.text.x = element_text(angle = 45, hjust = 1),  # 旋转 x 轴标签
    axis.text.y = element_text(),  # y 轴标签不加粗
    axis.title.y = element_text(face = "bold"),  # 加粗 y 轴标题（Essentiality Score）
    plot.title = element_text(hjust = 0.5, face = "bold", color = "black"),  # 居中标题并设置为黑色
    legend.key = element_rect(colour = NA)  # 移除图例（颜色注释）的黑框
  ) +
  labs(x = "Tissue", y = "Essentiality Score", title = paste0(gene, " — Primary Disease"))


# single tissue -----------------------------------------------------------

single_tissue <- "BONE"

# 过滤数据
filtered_data <- plotdata %>% filter(Tissue == single_tissue)

# 绘制单个组织的箱线图
ggplot(filtered_data, aes(x = Tissue, y = Value, fill = Tissue)) +
  geom_boxplot() +
  scale_fill_manual(values = colors) +  # 应用自定义颜色调色板
  theme(
    panel.background = element_blank(),  # 设置背景为空白
    panel.border = element_rect(colour = "black", fill = NA, size = 1),  # 添加黑方框
    axis.text.x = element_text(angle = 45, hjust = 1),  # 旋转 x 轴标签
    axis.text.y = element_text(),  # y 轴标签不加粗
    axis.title.y = element_text(face = "bold"),  # 加粗 y 轴标题（Essentiality Score）
    plot.title = element_text(hjust = 0.5, face = "bold", color = "black"),  # 居中标题并设置为黑色
    legend.key = element_rect(colour = NA)  # 移除图例（颜色注释）的黑框
  ) +
  labs(x = "Tissue", y = "Essentiality Score", title = paste0(gene, " — ", single_tissue))


# Sort plotdata by Value from smallest to largest
plotdata <- plotdata[order(plotdata$Value), ]

# Ensure CellName is a factor with levels ordered by Value
plotdata$CellName <- factor(plotdata$CellName, levels = unique(plotdata$CellName[order(plotdata$Value)]))

# Generate the dot plot
# ggplot(plotdata, aes(x = CellName, y = Value, color = Tissue)) +
#   geom_jitter(width = 0.1, height = 0, alpha = 0.7, size = 1.5) +  # Reduce point size and jitter width
#   scale_color_manual(values = colors) +  # Apply the custom color palette
#   theme_minimal(base_size = 12) +  # Set a base font size
#   theme(
#     axis.text.x = element_blank(),  # Remove x-axis labels
#     axis.ticks.x = element_blank(),  # Remove x-axis ticks
#     panel.grid.major.x = element_blank(),  # Remove major grid lines for x-axis
#     panel.grid.minor.x = element_blank(),  # Remove minor grid lines for x-axis
#     plot.title = element_text(hjust = 0.5, face = "bold", color = "black")  # Center-align and set title to black
#   ) +
#   labs(y = "Essentiality Score", color = "Tissue", title = paste0(gene, " — Primary Disease")) +
#   coord_cartesian(clip = "off")  # Ensure all points are within the plotting area

ggplot(plotdata, aes(x = CellName, y = Value, color = Tissue)) +
  geom_jitter(width = 0.1, height = 0, alpha = 0.7, size = 1.5) +  # 调整点的抖动、透明度和大小
  scale_color_manual(values = colors) +  # 应用自定义颜色调色板
  theme(
    panel.background = element_blank(),  # 设置背景为空白
    panel.border = element_rect(colour = "black", fill = NA, size = 1),  # 添加黑方框
    axis.text.x = element_blank(),  # 移除 x 轴标签
    axis.ticks.x = element_blank(),  # 移除 x 轴刻度线
    axis.title.x = element_text(face = "bold", size = 12),  # 设置 x 轴标题为 "Cell lines" 并加粗
    panel.grid.major.x = element_blank(),  # 移除 x 轴主要网格线
    panel.grid.minor.x = element_blank(),  # 移除 x 轴次要网格线
    plot.title = element_text(hjust = 0.5, face = "bold", color = "black")  # 居中标题并设置为黑色
  ) +
  labs(x = "Cell lines", y = "Essentiality Score", color = "Tissue", title = paste0(gene, " — Primary Disease")) +
  coord_cartesian(clip = "off")  # 确保所有点都在绘图区域内

ggplot(filtered_data, aes(x = CellName, y = Value, color = Tissue)) +
  geom_jitter(width = 0.1, height = 0, alpha = 0.7, size = 1.5) +  # 调整点的抖动、透明度和大小
  scale_color_manual(values = colors) +  # 应用自定义颜色调色板
  theme(
    panel.background = element_blank(),  # 设置背景为空白
    panel.border = element_rect(colour = "black", fill = NA, size = 1),  # 添加黑方框
    axis.text.x = element_blank(),  # 移除 x 轴标签
    axis.ticks.x = element_blank(),  # 移除 x 轴刻度线
    axis.title.x = element_text(face = "bold", size = 12),  # 设置 x 轴标题为 "Cell lines" 并加粗
    panel.grid.major.x = element_blank(),  # 移除 x 轴主要网格线
    panel.grid.minor.x = element_blank(),  # 移除 x 轴次要网格线
    plot.title = element_text(hjust = 0.5, face = "bold", color = "black")  # 居中标题并设置为黑色
  ) +
  labs(x = "Cell lines", y = "Essentiality Score", color = "Tissue", title = paste0(gene, " — ", single_tissue)) +
  coord_cartesian(clip = "off")  # 确保所有点都在绘图区域内
