rm(list = ls())

# Load necessary libraries
library(ggplot2)
library(ggpubr)
library(ggrepel)
library(RColorBrewer)
library(dplyr)

load("DEMap_Data/CRISPR_DATA.Rdata")
load("GEMap_Outputdata/09_targets_10.Rdata")

com_target <- intersect(depmap_23q2, dgidb)
df_crispr <- df_crispr[, com_target]

###  CN  ########################################################################
COM <- intersect(rownames(df_crispr), rownames(df_cn))
df_cn <- df_cn[COM, ]

### CCNE1为例
cnGene = "CCNE1"

cnGeneData =  data.frame(DepMap_ID=rownames(df_cn),cnGene = df_cn[,cnGene])

ampCells = cnGeneData$DepMap_ID[cnGeneData$cnGene > 1.58]
wtCells = cnGeneData$DepMap_ID[cnGeneData$cnGene <= 1.58]

### 批量处理数据
ampCellsDepData = df_crispr[ampCells,]
wtCellsDepData = df_crispr[wtCells,]

#### 批量执行
results <- data.frame()
for (i in 1:ncol(df_crispr)) {
  ## A标记
  print(i)
  
  ## B提取数据
  # ampD = ampCellsDepData[,i]
  # wtD = wtCellsDepData[,i]
  ampD <- na.omit(ampCellsDepData[, i])
  wtD <- na.omit(wtCellsDepData[, i])
  ## C数据计算以及储存
  ## 1.GeneSymbol
  results[i,1] = colnames(df_crispr)[i]
  ## 2.med.dep amp
  results[i,2] = median(ampD,na.rm = T)
  ## 3.med.dep wt
  results[i,3] = median(wtD,na.rm = T)
  ## 4.dep FC
  results[i,4] = median(ampD,na.rm = T) - median(wtD,na.rm = T)
  
  ## 5.mean.dep amp
  results[i,5] = mean(ampD,na.rm = T)
  ## 6.mean.dep wt
  results[i,6] = mean(wtD,na.rm = T)
  ## 7.dep FC mean
  results[i,7] = mean(ampD,na.rm = T) - mean(wtD,na.rm = T)
  
  # 样本数量检查并进行统计测试
  if (length(ampD) > 1 & length(wtD) > 1) {
    results[i, 8] <- wilcox.test(ampD, wtD, alternative = "less")$p.value
    results[i, 9] <- t.test(ampD, wtD, alternative = "less")$p.value
  } else {
    results[i, 8] <- NA
    results[i, 9] <- NA
  }
}

colnames(results) = c("GeneSymbol",
                      "amp_median_dep",
                      "wt_median_dep",
                      "dep.FC",
                      "amp_mean_dep",
                      "wt_mean_dep",
                      "dep.FC.mean",
                      "Wilcox_p_value",
                      "Ttest_p_value"
)

### 增加矫正的p值
results$Wilcox_p_value_fdr = p.adjust(results$Wilcox_p_value,method ="fdr")
results$Ttest_p_value_fdr = p.adjust(results$Ttest_p_value,method ="fdr")

data <- results
data$rank <- -log10(data$Wilcox_p_value)*data$dep.FC

genes_to_show <- data %>% 
  arrange(rank) %>%
  dplyr::slice(1:10) %>%
  pull(GeneSymbol)

ggplot(data, aes(dep.FC, -log10(Wilcox_p_value))) +
  geom_point(shape = 21, colour = "gray30", aes(fill = -log10(Wilcox_p_value)), size = 3) + 
  scale_fill_gradient(low = "#9ECFE3", high = "#FB9F5A") +  # 从浅蓝到橙色的渐变填充
  geom_point(data = subset(data, GeneSymbol %in% genes_to_show), 
             shape = 21, colour = "gray30", fill = "#F0663F", size = 3, stroke = 1.2) +  # 高亮点的颜色为橙红色
  geom_text_repel(data = subset(data, GeneSymbol %in% genes_to_show), 
                  aes(label = GeneSymbol), size = 5, fontface = "bold", color = "black", 
                  box.padding = 0.35, point.padding = 0.5, max.overlaps = 10) +  # 文本位置调整
  theme(
    panel.background = element_blank(),  # 设置背景为空白
    panel.border = element_rect(colour = "black", fill = NA, size = 1),  # 添加黑方框
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),  # 标题居中加粗
    plot.subtitle = element_text(hjust = 0.5, size = 14, margin = margin(b = 10)),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    legend.position = "right", 
    legend.title = element_blank(),  # 移除图例标题
    legend.key.size = unit(0.8, "cm"),  # 调整图例键大小
    legend.text = element_text(size = 12)
  ) +
  xlab(expression(paste(Delta,"FC (CCNE1 amplification vs WT)"))) +
  ylab(expression(paste("-log"[10], "(Wilcoxon ", italic("P"), " value)"))) +
  labs(title = "Overview")

###  CN vs WT  ########################################################################
# CDK2、NAPG、EMC7、ORC6、TFIP11、UFM1、CWC22、ST13、ADAR和PKMYT1
gene <- "CDK2"

# 合并数据并添加分类标签
combined_data <- rbind(
  data.frame(Value = ampCellsDepData$CDK2, Group = "AMP"),
  data.frame(Value = wtCellsDepData$CDK2, Group = "WT")
)

ggplot(combined_data, aes(x = Group, y = Value, fill = Group)) +
  geom_violin(trim = FALSE, alpha = 0.6, color = "black") +  # 小提琴图，半透明，黑色边框
  geom_boxplot(data = subset(combined_data, Group == "AMP"), 
               aes(x = Group, y = Value), 
               outlier.colour = "#6096C5", width = 0.1, outlier.size = 1.5, color = "black") +  # AMP 箱线图
  geom_boxplot(data = subset(combined_data, Group == "WT"), 
               aes(x = Group, y = Value), 
               outlier.colour = "#F0663F", width = 0.1, outlier.size = 1.5, color = "black") +  # WT 箱线图
  labs(
    x = "Amplification Status based on CCNE1", 
    y = "Essentiality Score", 
    title = "Essentiality Partner of CCNE1",
    fill = "Group"  # 添加颜色标签的名称
  ) +
  theme(
    panel.background = element_blank(),  # 设置背景为完全透明
    plot.background = element_blank(),   # 设置整个图的背景为完全透明
    panel.border = element_rect(colour = "black", fill = NA, size = 1),  # 添加黑方框
    axis.text = element_text(size = 12, color = "black"),  # 调整坐标轴文字大小和颜色
    axis.title = element_text(size = 14, face = "bold", color = "black"),  # 坐标轴标题加粗
    panel.grid.major = element_blank(),  # 移除主要网格线
    panel.grid.minor = element_blank(),  # 移除次要网格线
    legend.position = "right",  # 显示图例，位置在右侧
    legend.title = element_text(size = 12, face = "bold"),  # 图例标题加粗
    legend.text = element_text(size = 12),  # 图例文字大小
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5)  # 居中标题，加粗
  ) +
  scale_fill_manual(
    values = c("AMP" = "#7EB5D5", "WT" = "#F57948"),  # 自定义 AMP 和 WT 的颜色
    labels = c("AMP" = "Amp", "WT" = "WT")  # 设置颜色标签的名称
  ) +
  stat_compare_means(
    method = "wilcox.test", 
    label = "p.format", 
    label.x = 1.5, 
    label.y = 0.8  # 动态定位 p 值标签
  )  # 显示格式化后的 p 值


###  MUT  ########################################################################

mutData <- readRDS(file = "DEMap_Data/Mut_Data.rds")

### CCNE1为例
mutGene = "CCNE1"
mutGeneData =  data.frame(DepMap_ID=rownames(mutData),mutGene = mutData[,mutGene])

mutCells = mutGeneData$DepMap_ID[mutGeneData$mutGene > 0]
wtCells = mutGeneData$DepMap_ID[mutGeneData$mutGene == 0]

### 批量处理数据
mutCellsDepData = df_crispr[mutCells,]
wtCellsDepData = df_crispr[wtCells,]

#### 批量执行
results <- data.frame()
for (i in 1:ncol(df_crispr)) {
  ## A标记
  print(i)
  
  ## B提取数据
  mutD = mutCellsDepData[,i]
  wtD = wtCellsDepData[,i]
  
  ## C数据计算以及储存
  ## 1.GeneSymbol
  results[i,1] = colnames(df_crispr)[i]
  ## 2.med.dep mut
  results[i,2] = median(mutD,na.rm = T)
  ## 3.med.dep wt
  results[i,3] = median(wtD,na.rm = T)
  ## 4.dep FC
  results[i,4] = median(mutD,na.rm = T) - median(wtD,na.rm = T)
  
  ## 5.mean.dep mut
  results[i,5] = mean(mutD,na.rm = T)
  ## 6.mean.dep wt
  results[i,6] = mean(wtD,na.rm = T)
  ## 7.dep FC mean
  results[i,7] = mean(mutD,na.rm = T) - mean(wtD,na.rm = T)
  
  ## 8.
  results[i,8] = wilcox.test(mutD,wtD,alternative = "less")$p.value
  ## 9.
  results[i,9] = t.test(mutD,wtD,alternative = "less")$p.value
}

colnames(results) = c("GeneSymbol",
                      "mut_median_dep",
                      "wt_median_dep",
                      "dep.FC",
                      "mut_mean_dep",
                      "wt_mean_dep",
                      "dep.FC.mean",
                      "Wilcox_p_value",
                      "Ttest_p_value"
)

### 增加矫正的p值
results$Wilcox_p_value_fdr = p.adjust(results$Wilcox_p_value,method ="fdr")
results$Ttest_p_value_fdr = p.adjust(results$Ttest_p_value,method ="fdr")

data <- results
data$rank <- -log10(data$Wilcox_p_value)*data$dep.FC

genes_to_show <- data %>% 
  arrange(rank) %>%
  dplyr::slice(1:10) %>%
  pull(GeneSymbol)

ggplot(data, aes(dep.FC, -log10(Wilcox_p_value))) +
  geom_point(shape = 21, colour = "gray30", aes(fill = -log10(Wilcox_p_value)), size = 3) + 
  scale_fill_gradient(low = "#9ECFE3", high = "#FB9F5A") + # Gradient from light blue to dark blue
  geom_point(data = subset(data, GeneSymbol %in% genes_to_show), 
             shape = 21, colour = "gray30", fill = "#F0663F", size = 3, stroke = 1.2) + # Orange-red for highlighted points
  geom_text_repel(data = subset(data, GeneSymbol %in% genes_to_show), 
                  aes(label = GeneSymbol), size = 5, fontface = "bold", color = "black", 
                  box.padding = 0.35, point.padding = 0.5, max.overlaps = 10) + # Text positioning adjustments
  theme_minimal(base_size = 14) + # Cleaner, minimalist theme with larger base font size
  xlab(expression(paste(Delta,"FC (CCNE1 MUT vs WT)"))) +
  ylab(expression(paste("-log"[10], "(Wilcoxon ", italic("P"), " value)"))) +
  labs(title = "DepMap Dependencies", subtitle = "CCNE1-MUT Tumor Cell Lines") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"), # Larger, bold title
    plot.subtitle = element_text(hjust = 0.5, size = 14, margin = margin(b = 10)),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    legend.position = "right", 
    legend.title = element_blank(), # Remove legend title for a cleaner look
    legend.key.size = unit(0.8, "cm"), # Smaller legend key size
    legend.text = element_text(size = 12)
  )

###  MUT vs WT  ########################################################################
gene <- "CCNE1"

combined_data <- rbind(
  data.frame(Value = mutCellsDepData$POLR3A, Group = "MUT"),
  data.frame(Value = wtCellsDepData$POLR3A, Group = "WT")
)

ggplot(combined_data, aes(x = Group, y = Value, fill = Group)) +
  geom_violin(trim = FALSE, alpha = 0.6, color = "black") +  # Make the violin plot semi-transparent
  geom_boxplot(data = subset(combined_data, Group == "MUT"), 
               aes(x = Group, y = Value), 
               outlier.colour = "#6096C5", width = 0.1, outlier.size = 1.5, color = "black") +  # mut boxplot
  geom_boxplot(data = subset(combined_data, Group == "WT"), 
               aes(x = Group, y = Value), 
               outlier.colour = "#F0663F", width = 0.1, outlier.size = 1.5, color = "black") +  # WT boxplot
  labs(x = "MUT Status", y = "POLR1C Value", 
       title = paste0(gene, " —— MUT vs WT")) +
  theme_minimal(base_size = 15) +  # Minimal theme with larger base font size
  theme(
    axis.text = element_text(size = 12, color = "black"),  # Adjust axis text size and color
    axis.title = element_text(size = 14, face = "bold", color = "black"),  # Make axis titles bold
    panel.grid.major = element_line(color = "gray80", size = 0.5),  # Customize grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    legend.position = "none",  # Hide legend if not needed
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5)  # Centered plot title with bold font
  ) +
  scale_fill_manual(values = c("MUT" = "#7EB5D5", "WT" = "#F57948")) +  # Use a pastel color palette for the violin plot
  stat_compare_means(method = "wilcox.test", 
                     method.args = list(alternative =  "less"), 
                     label = "p.format", 
                     label.x = 1.5, label.y = 0.8)  # Show formatted p-value



