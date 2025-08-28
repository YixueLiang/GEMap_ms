rm(list = ls())

# Load necessary libraries
library(ggplot2)
library(ggpubr)
library(ggrepel)
library(stringr)
library(ggfun)
library(RColorBrewer)
library(dplyr)
library(tidyr)
library(tibble)


# cn
###########################################################################
load("DEMap_Data/CRISPR_DATA.Rdata")
cellinfor <- readRDS(file = "DEMap_Data/cellinfor.rds")

### 细胞系取交集，1072
coID <- intersect(rownames(df_crispr),rownames(cellinfor)) %>% 
  intersect(rownames(df_cn))

Tissue = "lineage"
mutGene <- "CCNE1"
targetGene <- "PKMYT1"

### 优化流程
## 数据筛选，Mut >=3,WT>=5
TissueStatus <- data.frame(
  Tissue = cellinfor[coID,Tissue],
  CN = ifelse(df_cn[coID,mutGene]<=1.58,"WT","Mut")) %>% 
  select(Tissue,CN) %>%
  group_by(Tissue,CN) %>% 
  summarise(n =n(),.groups = "drop") %>% 
  ungroup() %>% 
  pivot_wider(names_from = "CN",
              values_from = n,
              values_fill = 0) %>% 
  filter(Mut >=1,WT >=1) %>%
  mutate(Sum=Mut+WT) %>% 
  mutate(Mut_Number = sum(Mut),WT_Number=sum(WT)) 


data <- data.frame(
  DepmapID = coID,
  Tissue = cellinfor[coID,Tissue],
  CN = ifelse(df_cn[coID,mutGene]<=1.58,"WT","Mut"),
  Dependency = df_crispr[coID,targetGene]
)

mydata <- data %>% 
  filter(Tissue %in% TissueStatus$Tissue)

## 作图数据整理
Mean_Dep <- mydata %>% 
  filter(CN=="Mut") %>% 
  group_by(Tissue) %>% 
  summarise(Sensitivity = mean(Dependency))

Diff_Dep <- mydata %>%
  group_by(Tissue) %>%
  summarize(Difference = mean(Dependency[CN == "Mut"]) - mean(Dependency[CN == "WT"]))

## 数据合并
PlotData <- merge(Mean_Dep,Diff_Dep,by="Tissue") %>% 
  merge(TissueStatus,by="Tissue")

# plot --------------------------------------------------------------------

ggplot(PlotData, aes(x = Sensitivity, y = Difference)) +
  geom_point(aes(size = Sum, fill = Difference > 0), shape = 21, color = "black", alpha = 0.7, stroke = 1.5) + 
  geom_text_repel(aes(label = Tissue), vjust = 1.5, hjust = 1.5, size = 3.5, color = "black") + 
  scale_fill_manual(values = c("TRUE" = "#E31A1C", "FALSE" = "#4575B4")) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0, ymax = Inf,
           fill = "lightyellow", alpha = 0.3, colour = "#984EA3", linetype = "dashed", linewidth = 1) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = 0,
           fill = "lightblue", alpha = 0.3, colour = "#66C2A5", linetype = "dashed", linewidth = 1) +
  labs(x = "Sensitivity Score for Cell Lines with CN", 
       y = "AMP vs. WT Sensitivity", 
       title = paste0("AMP vs. WT: ", mutGene, " & ", targetGene)) +
  theme_minimal(base_size = 15) + 
  guides(size = guide_legend(title = "Number of Cells"), fill = guide_legend(title = "Difference > 0")) +
  theme(
    panel.grid.major = element_line(linetype = "dashed", color = "grey80"),
    panel.grid.minor = element_line(linetype = "dotted", linewidth = 0.5, color = "grey85"),
    panel.background = element_blank(),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    legend.position = "right",  # Place legend to the right
    legend.justification = "center",
    legend.box = "vertical",  # Arrange legend items vertically
    legend.box.margin = margin(0, 10, 0, 0),  # Add spacing between plot and legend
    legend.margin = margin(10, 10, 10, 10),
    legend.background = element_rect(fill = "white", color = "#808080", linetype = 1),
    legend.title = element_text(face = "bold"),
    legend.text = element_text(size = 12)
  ) +
  theme(plot.margin = margin(10, 50, 10, 10))  # Add space for the legend outside the plot area


# mut
###########################################################################
mutData <- readRDS(file = "DEMap_Data/Mut_Data.rds")

### 细胞系取交集，1072
coID <- intersect(rownames(df_crispr),rownames(cellinfor)) %>% 
  intersect(rownames(mutData))
cellinfor$OncotreeLineage
Tissue = "OncotreeLineage"
mutGene <- "CCNE1"
targetGene <- "PKMYT1"

### 优化流程
## 数据筛选，Mut >=3,WT>=5
TissueStatus <- data.frame(
  Tissue = cellinfor[coID,Tissue],
  Mutation = ifelse(mutData[coID,mutGene]==0,"WT","Mut")
) %>% 
  select(Tissue,Mutation) %>%
  group_by(Tissue,Mutation) %>% 
  summarise(n =n(),.groups = "drop") %>% 
  ungroup() %>% 
  pivot_wider(names_from = "Mutation",
              values_from = n,
              values_fill = 0) %>% 
  filter(Mut >=0,WT>=0) %>%
  mutate(Sum=Mut+WT) %>% 
  mutate(Mut_Number = sum(Mut),WT_Number=sum(WT)) 


data <- data.frame(
  DepmapID = coID,
  Tissue = cellinfor[coID,Tissue],
  Mutation = ifelse(mutData[coID,mutGene]==0,"WT","Mut"),
  Dependency = df_crispr[coID,targetGene]
)

mydata <- data %>% 
  filter(Tissue %in% TissueStatus$Tissue)

## 作图数据整理
Mean_Dep <- mydata %>% 
  filter(Mutation=="Mut") %>% 
  group_by(Tissue) %>% 
  summarise(Sensitivity = mean(Dependency))

Diff_Dep <- mydata %>%
  group_by(Tissue) %>%
  summarize(Difference = mean(Dependency[Mutation == "Mut"]) - mean(Dependency[Mutation == "WT"]))

## 数据合并
PlotData <- merge(Mean_Dep,Diff_Dep,by="Tissue") %>% 
  merge(TissueStatus,by="Tissue")

# plot --------------------------------------------------------------------

ggplot(PlotData, aes(x = Sensitivity, y = Difference)) +
  geom_point(aes(size = Sum, fill = Difference > 0), shape = 21, color = "black", alpha = 0.7, stroke = 1.5) + 
  geom_text_repel(aes(label = Tissue), vjust = 1.5, hjust = 1.5, size = 3.5, color = "black") + 
  scale_fill_manual(values = c("TRUE" = "#E31A1C", "FALSE" = "#4575B4")) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0, ymax = Inf,
           fill = "lightyellow", alpha = 0.3, colour = "#984EA3", linetype = "dashed", linewidth = 1) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = 0,
           fill = "lightblue", alpha = 0.3, colour = "#66C2A5", linetype = "dashed", linewidth = 1) +
  labs(x = "Sensitivity Score for Cell Lines with Mutations", 
       y = "Mutant vs. Wildtype Sensitivity", 
       title = paste0("Mutant vs. Wildtype: ", mutGene, " & ", targetGene)) +
  theme_minimal(base_size = 15) + 
  guides(size = guide_legend(title = "Number of Cells"), fill = guide_legend(title = "Difference > 0")) +
  theme(
    panel.grid.major = element_line(linetype = "dashed", color = "grey80"),
    panel.grid.minor = element_line(linetype = "dotted", linewidth = 0.5, color = "grey85"),
    panel.background = element_blank(),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    legend.position = "right",  # Place legend to the right
    legend.justification = "center",
    legend.box = "vertical",  # Arrange legend items vertically
    legend.box.margin = margin(0, 10, 0, 0),  # Add spacing between plot and legend
    legend.margin = margin(10, 10, 10, 10),
    legend.background = element_rect(fill = "white", color = "#808080", linetype = 1),
    legend.title = element_text(face = "bold"),
    legend.text = element_text(size = 12)
  ) +
  theme(plot.margin = margin(10, 50, 10, 10))  # Add space for the legend outside the plot area


