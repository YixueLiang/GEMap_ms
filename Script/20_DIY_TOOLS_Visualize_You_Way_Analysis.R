#清空数据
rm(list = ls())
library(ggplot2)
library(RColorBrewer)

load("DEMap_Data/CRISPR_DATA.Rdata")

genes <- c("A1BG", "A1CF")
Tissues <- c("Breast", "Lymphoid", "Uterus")

commoncell <- intersect(rownames(df_crispr), df_cell$ModelID)
rownames(df_cell) <- df_cell$ModelID
df_cell <- df_cell[rownames(df_crispr), ]
rownames(df_crispr) <- df_cell$StrippedCellLineName
Tissue <- unique(df_cell$OncotreeLineage)

# 创建分组织向量
Tissue_groups <- list()
Tissue_groups <- split(df_cell$StrippedCellLineName, df_cell$OncotreeLineage)

# 创建分组织后的 df_crispr
crispr_groups <- list()

# 根据 Tissue_groups 分组 df_crispr
for (Tissue in Tissues) {
  crispr_groups[[Tissue]] <- df_crispr[Tissue_groups[[Tissue]], ]
}


gene_values <- list()

# 遍历基因和组织，提取基因值
for (gene in genes) {
  for (Tissue in names(crispr_groups)) {
    # 为每个基因创建一个子列表
    if (!is.list(gene_values[[gene]])) {
      gene_values[[gene]] <- list()
    }
    # 直接通过嵌套列表索引
    if (!is.null(crispr_groups[[Tissue]][[gene]])) {
      gene_values[[gene]][[Tissue]] <- crispr_groups[[Tissue]][[gene]]
    } else {
      cat("Gene", gene, "not found in Tissue", Tissue, "\n")
      gene_values[[gene]][[Tissue]] <- NA
    }
  }
}

# 创建一个空的数据框来存储长格式数据
plotdata <- data.frame(
  Gene = character(0),
  Value = numeric(0),
  Tissue = character(0),
  CellName = character(0)
)

# 遍历每个基因和组织，将向量数据添加到长格式数据框中
for (gene in names(gene_values)) {
  for (Tissue in names(gene_values[[gene]])) {
    temp_df <- data.frame(
      Gene = gene,
      Value = gene_values[[gene]][[Tissue]],
      Tissue = Tissue,
      CellName = rownames(crispr_groups[[Tissue]]) # 从行名中提取 CellName
    )
    plotdata <- rbind(plotdata, temp_df)
  }
}


library(ggplot2)
library(patchwork)

colors <- colorRampPalette(rev(brewer.pal(9, "RdYlBu")))(length(unique(plotdata$Tissue)))
# 绘制第一个基因的箱线图
plot_gene1 <- ggplot(plotdata, aes(x = Tissue, y = Value, fill = Tissue)) +
  geom_boxplot(outlier.size = 1, outlier.color = "gray40", lwd = 0.5) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, color = "black", size = 10),
    axis.text.y = element_text(color = "black", size = 10),
    axis.title = element_text(face = "bold", size = 12),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16), # 标题居中
    legend.position = "right", # 显示图例
    legend.title = element_text(face = "bold"), # 图例标题加粗
    legend.text = element_text(size = 10) # 图例文本大小
  ) +
  labs(title = "Boxplot for A1BG", y = "Expression Value", x = "Tissue") +
  scale_fill_manual(values = colors)

# 绘制第二个基因的箱线图
plot_gene2 <- ggplot(plotdata, aes(x = Tissue, y = Value, fill = Tissue)) +
  geom_boxplot(outlier.size = 1, outlier.color = "gray40", lwd = 0.5) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, color = "black", size = 10),
    axis.text.y = element_text(color = "black", size = 10),
    axis.title = element_text(face = "bold", size = 12),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16), # 标题居中
    legend.position = "right", # 显示图例
    legend.title = element_text(face = "bold"), # 图例标题加粗
    legend.text = element_text(size = 10) # 图例文本大小
  ) +
  labs(title = "Boxplot for A1CF", y = "Expression Value", x = "Tissue") +
  scale_fill_manual(values = colors)

# 拼接两个图
combined_plot <- plot_gene1 / plot_gene2 +
  plot_annotation(
    title = "Expression Profiles of A1BG and A1CF Across Tissues",
    theme = theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 18)
    )
  )

# 显示拼接图
print(combined_plot)
