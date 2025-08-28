rm(list = ls())
library(ggplot2)
library(ggrepel)
library(dplyr)
library(ggpubr)
library(RColorBrewer)

# mut_vs_wt ---------------------------------------------------------------

# Load the data
load("DEMap_Data/CRISPR_DATA.Rdata")
mutData <- readRDS(file = "DEMap_Data/Mut_Data.rds")

df_crispr <- df_crispr[rownames(df_crispr) != "ACH-002198", ]
common <- intersect(rownames(mutData), rownames(df_crispr))

df_crispr <- df_crispr[common, ]
mutData <- mutData[common, ]

commoncell <- intersect(rownames(df_crispr), df_cell$ModelID)
rownames(df_cell) <- df_cell$ModelID
df_cell <- df_cell[rownames(df_crispr), ]
rownames(df_crispr) <- df_cell$StrippedCellLineName
rownames(mutData) <- df_cell$StrippedCellLineName
Tissue <- unique(df_cell$OncotreeLineage)

genes <- c("CCNE1", "PKMYT1")

# Create tissue groups
Tissue_groups <- list()
Tissue_groups <- split(df_cell$StrippedCellLineName, df_cell$OncotreeLineage)

# Define the custom color palette for tissues
colors <- rev(colorRampPalette(rev(brewer.pal(9, "RdYlBu")))(27))


# 保存突变和野生型的细胞列表
mutCells_list <- list()
wtCells_list <- list()

for (gene in genes) {  # 遍历基因列
  mutCells <- rownames(mutData)[mutData[[gene]] > 0]
  wtCells <- rownames(mutData)[mutData[[gene]] == 0]
  
  mutCells_list[[gene]] <- mutCells
  wtCells_list[[gene]] <- wtCells
}

gene <- "CCNE1"
mutCells = mutCells_list[[gene]]
wtCells = wtCells_list[[gene]]

mutCellsDepData = df_crispr[mutCells,]
wtCellsDepData = df_crispr[wtCells,]

# 合并数据并添加分类标签
combined_data <- rbind(
  data.frame(Value = mutCellsDepData$PKMYT1, Group = "MUT"),
  data.frame(Value = wtCellsDepData$PKMYT1, Group = "WT")
)

# Add row names to combined_data
rownames(combined_data) <- c(rownames(mutCellsDepData), rownames(wtCellsDepData))


ggplot(combined_data, aes(x = Group, y = Value, fill = Group)) +
  geom_violin(trim = FALSE, alpha = 0.6, color = "black") +  # 使小提琴图半透明
  geom_boxplot(data = subset(combined_data, Group == "MUT"), 
               aes(x = Group, y = Value), 
               outlier.colour = "#6096C5", width = 0.1, outlier.size = 1.5, color = "black") +  # MUT 箱线图
  geom_boxplot(data = subset(combined_data, Group == "WT"), 
               aes(x = Group, y = Value), 
               outlier.colour = "#F0663F", width = 0.1, outlier.size = 1.5, color = "black") +  # WT 箱线图
  labs(
    x = "Mutation Status of CCNE1", 
    y = "Essentiality Score of PKMYT1", 
    title = paste0(gene, " —— Mut vs WT")
  ) +  # 添加标题
  theme_minimal(base_size = 15) +  # 使用最小主题，并设置较大的基础字体大小
  theme(
    axis.text = element_text(size = 12, color = "black"),  # 调整轴文本大小和颜色
    axis.title = element_text(size = 14, face = "bold", color = "black"),  # 使轴标题加粗
    panel.grid.major = element_line(color = "gray80", size = 0.5),  # 自定义网格线
    panel.grid.minor = element_blank(),  # 移除次要网格线
    legend.position = "none",  # 隐藏图例
    plot.title = element_text(size = 16, face = "bold", color = "black", hjust = 0.5),  # 居中标题并加粗
    panel.background = element_rect(fill = "white")  # 将背景设置为白色
  ) +
  scale_fill_manual(values = c("MUT" = "#7EB5D5", "WT" = "#F57948")) +  # 自定义 MUT 和 WT 的颜色
  stat_compare_means(
    method = "wilcox.test", 
    label = "p.format", 
    method.args=list(exact=F), 
    label.x = 1.5, 
    label.y = max(combined_data$Value) * 1.1  # 动态定位标签在数据上方
  )  # 显示格式化后的 p 值


# tissue ------------------------------------------------------------------


# Extract gene expression data for CCNE1 and PKMYT1 from df_crispr
combined_data_tissue <- combined_data[Tissue_groups[["Ovary/Fallopian Tube"]], , drop = FALSE]


ggplot(combined_data_tissue, aes(x = Group, y = Value, fill = Group)) +
  geom_violin(trim = FALSE, alpha = 0.6, color = "black") +  
  geom_boxplot(data = subset(combined_data, Group == "MUT"), 
               aes(x = Group, y = Value), 
               outlier.colour = "#6096C5", width = 0.1, outlier.size = 1.5, color = "black") +  
  geom_boxplot(data = subset(combined_data, Group == "WT"), 
               aes(x = Group, y = Value), 
               outlier.colour = "#F0663F", width = 0.1, outlier.size = 1.5, color = "black") +  
  labs(
    x = "Mutation Status", 
    y = "POLR1C Value", 
    title = paste0(gene, " —— Mut vs WT")
  ) +  # Add a title
  theme_minimal(base_size = 15) +  # Minimal theme with larger base font size
  theme(
    axis.text = element_text(size = 12, color = "black"),  # Adjust axis text size and color
    axis.title = element_text(size = 14, face = "bold", color = "black"),  # Make axis titles bold
    panel.grid.major = element_line(color = "gray80", size = 0.5),  # Customize grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    legend.position = "none",  # Hide legend if not needed
    plot.title = element_text(size = 16, face = "bold", color = "black", hjust = 0.5)  
  ) +
  scale_fill_manual(values = c("MUT" = "#7EB5D5", "WT" = "#F57948")) +  
  stat_compare_means(
    method = "wilcox.test", 
    label = "p.format", 
    label.x = 1.5, 
    label.y = max(combined_data$Value) * 1.1  
  ) 


