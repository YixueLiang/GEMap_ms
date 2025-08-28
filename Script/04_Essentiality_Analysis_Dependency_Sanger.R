#清空数据
rm(list = ls())
library(ggplot2)
library(RColorBrewer)
library(dplyr)

load("DEMap_Data/Sanger_DATA.Rdata")
df_cell = data.table::fread('inputdata/DepMap_Public_23Q2/Model.csv')
df_cell <- as.data.frame(df_cell)
rownames(df_cell) <- df_cell$StrippedCellLineName

df_cell$OncotreeLineage <- gsub("Ovary/Fallopian Tube", "Ovary", df_cell$OncotreeLineage)
df_cell$OncotreeLineage <- gsub("Bladder/Urinary Tract", "Bladder", df_cell$OncotreeLineage)
df_cell$OncotreeLineage <- gsub("CNS/Brain", "CNS", df_cell$OncotreeLineage)
df_cell$OncotreeLineage <- gsub("Vulva/Vagina", "Vulva", df_cell$OncotreeLineage)
df_cell$OncotreeLineage <- gsub("Esophagus/Stomach", "Esophagus and Stomach", df_cell$OncotreeLineage)

COM_cell <- intersect(rownames(df_sanger), df_cell$StrippedCellLineName)
df_cell <- df_cell[COM_cell, ]

Tissue <- unique(df_cell$OncotreeLineage)

# 创建分组织向量
Tissue_groups <- list()
Tissue_groups <- split(df_cell$StrippedCellLineName, df_cell$OncotreeLineage)

# 创建分组织后的 df_sanger
sanger_groups <- list()

# 根据 Tissue_groups 分组 df_sanger
for (Tissue in names(Tissue_groups)) {
  sanger_groups[[Tissue]] <- df_sanger[Tissue_groups[[Tissue]], ]
}



# 变量gene画plot -------------------------------------------------------------

gene <- "CDK4"
gene_values <- list()

for (Tissue in names(sanger_groups)) {
  gene_values[[Tissue]] <- sanger_groups[[Tissue]]["CDK4"]
}

# 创建一个空列表来存储每个组织的向量
sanger_vectors <- list()

# 遍历每个组织的数据框
for (Tissue in names(gene_values)) {
  # 将数据框转换为矩阵
  matrix_df <- as.matrix(gene_values[[Tissue]])
  
  # 将矩阵展平成一个向量
  vector_df <- as.vector(matrix_df)
  
  # 将向量存储在列表中
  sanger_vectors[[Tissue]] <- vector_df
}


# 创建一个空的数据框来存储长格式数据
plotdata <- data.frame(Value = numeric(0), Tissue = character(0), CellName = character(0))

# 将每个组织的向量添加到长格式数据框中
for (Tissue in names(sanger_vectors)) {
  temp_df <- data.frame(
    Value = sanger_vectors[[Tissue]],
    Tissue = Tissue,
    CellName = rownames(gene_values[[Tissue]])  # Add CellName from row names
  )
  plotdata <- rbind(plotdata, temp_df)
}

# Remove tissues labeled as "other"
plotdata <- plotdata %>% filter(!grepl("other", Tissue, ignore.case = TRUE))

# Reorder Tissue factor levels by the median of Value
plotdata$Tissue <- reorder(plotdata$Tissue, plotdata$Value, FUN = median)

# Generate a color palette with 28 distinct colors
colors <- rev(colorRampPalette(brewer.pal(9, "RdYlBu"))(20))

ggplot(plotdata, aes(x = Tissue, y = Value, fill = Tissue)) +
  geom_boxplot() +
  scale_fill_manual(values = colors) +  # Apply the custom color palette
  theme_minimal() +  # Use a clean, minimal theme
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
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
    panel.grid.minor.x = element_blank()   # Remove minor grid lines for x-axis
  ) +
  labs(y = "Value", color = "Tissue", title = paste0(gene, " — Primary Disease")) +
  coord_cartesian(clip = "off")  # Ensure all points are within the plotting area
