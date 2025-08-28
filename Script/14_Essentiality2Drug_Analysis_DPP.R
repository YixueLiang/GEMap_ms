#清空数据
rm(list = ls())
library(ggplot2)
library(RColorBrewer)
library(dplyr)

Gene <- "DOCK5"
Tissue <- "Breast"

load("DEMap_Data/drug_data_list_20_new.RData")

Tissue_data <- drug_data_list[[Tissue]]
Tissue_data <- Tissue_data[,c("pert_iname", "cell_id", "OncotreeLineage", Gene)]

length(unique(Tissue_data$pert_iname))
Tissue_data <- Tissue_data %>%
  group_by(pert_iname) %>%            # 按pert_iname进行分组
  summarise(DOCK5 = mean(DOCK5, na.rm = TRUE), # 计算PKMYT1的平均值，na.rm = TRUE 处理NA值
            .groups = "drop")            # 取消分组


# Sort the data by the values of the specified column in ascending order
Tissue_data <- Tissue_data[order(Tissue_data[[Gene]]), ]

# Extract the first 15 rows (minimum values)
top_min <- Tissue_data[1:15, ]

# Extract the last 15 rows (maximum values)
top_max <- Tissue_data[(nrow(Tissue_data) - 14):nrow(Tissue_data), ]

# Combine the top and bottom rows
result <- rbind(top_min, top_max)



library(RColorBrewer)

# Create color palette
colors <- colorRampPalette(rev(brewer.pal(9, "RdYlBu")))(30)

# Assign colors based on expression level (not group)
result <- result %>%
  arrange(!!sym(Gene)) %>%
  mutate(color_index = cut(!!sym(Gene), breaks = 30, labels = FALSE))


## plot

# 计算数据的范围，确保0在中心
data_range <- range(result[[Gene]], na.rm = TRUE)
max_abs_value <- max(abs(data_range[1]))
min_abs_value <- max(abs(data_range[2]))
limits <- c(-max_abs_value, min_abs_value)


ggplot(result) +
  # Add a box around the entire plot (using annotate instead of geom_rect to avoid warning)
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, 
           color = "black", fill = NA, linewidth = 0.5) +
  
  # Make sure pert_iname is ordered by the Gene value
  geom_segment(
    aes(x = reorder(pert_iname, !!sym(Gene)), 
        xend = reorder(pert_iname, !!sym(Gene)), 
        y = 0, yend = !!sym(Gene)),
    color = "gray70",
    linewidth = 0.8
  ) +
  geom_point(
    aes(x = reorder(pert_iname, !!sym(Gene)), 
        y = !!sym(Gene), 
        fill = !!sym(Gene)),
    size = 6,
    shape = 21,
    color = "black",
    stroke = 0.5
  ) +
  scale_fill_gradientn(
    colors = colors,
    limits = limits,
    name = "Level",
    guide = guide_colorbar(
      title.position = "top",
      title.hjust = 0.5,
      barwidth = unit(0.3, "cm"),
      barheight = unit(5, "cm"),
      frame.colour = "black",
      ticks.colour = "black",
      direction = "vertical"
    )
  ) +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.5, linetype = "dashed") +
  
  labs(
    x = NULL, 
    y = paste(Gene, "Expression Value"),
    title = paste("Top 15 and Bottom 15", Gene, "Values in", Tissue, "Tissue")
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5, margin = margin(b = 15)),
    axis.text.y = element_text(size = 11, color = "black"),
    axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
    axis.title = element_text(size = 14),
    axis.title.x = element_text(margin = margin(t = 10)),
    legend.position = "right",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.margin = unit(c(1, 1, 1, 1), "cm"),
    plot.background = element_rect(fill = "white", color = NA),
    panel.border = element_blank()
  ) +
  geom_vline(xintercept = 15.5, linetype = "dashed", color = "gray30")





# LINCS APPROVED ----------------------------------------------------------

mydata <- read.delim("inputdata/lincs/repurposing_drugs_20200324.txt",
                     skip = 9,  
                     stringsAsFactors = FALSE,
                     quote = "\"",  
                     comment.char = "!",  
                     fill = TRUE)  


Lincs_Approved_drugs <- intersect(mydata$pert_iname, unique(drug_data$pert_iname))

drug_data2 <- drug_data[drug_data$pert_iname %in% Lincs_Approved_drugs,]
drug_data_list2 <- split(drug_data2, drug_data2$OncotreeLineage)

Tissue_data <- drug_data_list2[[Tissue]]
Tissue_data <- Tissue_data[,c("pert_iname", "cell_id", "OncotreeLineage", Gene)]

length(unique(Tissue_data$pert_iname))
Tissue_data <- Tissue_data %>%
  group_by(pert_iname) %>%            # 按pert_iname进行分组
  summarise(DOCK5 = mean(DOCK5, na.rm = TRUE), # 计算PKMYT1的平均值，na.rm = TRUE 处理NA值
            .groups = "drop")            # 取消分组


# Sort the data by the values of the specified column in ascending order
Tissue_data <- Tissue_data[order(Tissue_data[[Gene]]), ]

# Extract the first 15 rows (minimum values)
top_min <- Tissue_data[1:15, ]

# Extract the last 15 rows (maximum values)
top_max <- Tissue_data[(nrow(Tissue_data) - 14):nrow(Tissue_data), ]

# Combine the top and bottom rows
result <- rbind(top_min, top_max)



library(RColorBrewer)

# Create color palette
colors <- colorRampPalette(rev(brewer.pal(9, "RdYlBu")))(30)

# Assign colors based on expression level (not group)
result <- result %>%
  arrange(!!sym(Gene)) %>%
  mutate(color_index = cut(!!sym(Gene), breaks = 30, labels = FALSE))


## plot

# 计算数据的范围，确保0在中心
data_range <- range(result[[Gene]], na.rm = TRUE)
max_abs_value <- max(abs(data_range[1]))
min_abs_value <- max(abs(data_range[2]))
limits <- c(-max_abs_value, min_abs_value)


ggplot(result) +
  # Add a box around the entire plot (using annotate instead of geom_rect to avoid warning)
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, 
           color = "black", fill = NA, linewidth = 0.5) +
  
  # Make sure pert_iname is ordered by the Gene value
  geom_segment(
    aes(x = reorder(pert_iname, !!sym(Gene)), 
        xend = reorder(pert_iname, !!sym(Gene)), 
        y = 0, yend = !!sym(Gene)),
    color = "gray70",
    linewidth = 0.8
  ) +
  geom_point(
    aes(x = reorder(pert_iname, !!sym(Gene)), 
        y = !!sym(Gene), 
        fill = !!sym(Gene)),
    size = 6,
    shape = 21,
    color = "black",
    stroke = 0.5
  ) +
  scale_fill_gradientn(
    colors = colors,
    limits = limits,
    name = "Level",
    guide = guide_colorbar(
      title.position = "top",
      title.hjust = 0.5,
      barwidth = unit(0.3, "cm"),
      barheight = unit(5, "cm"),
      frame.colour = "black",
      ticks.colour = "black",
      direction = "vertical"
    )
  ) +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.5, linetype = "dashed") +
  
  labs(
    x = NULL, 
    y = paste(Gene, "Expression Value"),
    title = paste("Top 15 and Bottom 15", Gene, "Values in", Tissue, "Tissue")
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5, margin = margin(b = 15)),
    axis.text.y = element_text(size = 11, color = "black"),
    axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
    axis.title = element_text(size = 14),
    axis.title.x = element_text(margin = margin(t = 10)),
    legend.position = "right",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.margin = unit(c(1, 1, 1, 1), "cm"),
    plot.background = element_rect(fill = "white", color = NA),
    panel.border = element_blank()
  ) +
  geom_vline(xintercept = 15.5, linetype = "dashed", color = "gray30")





# ggplot(result) +
#   # Add a box around the entire plot
#   geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf), 
#             color = "black", fill = NA, linewidth = 0.5) +
#   
#   geom_segment(
#     aes(x = pert_iname, xend = pert_iname, y = 0, yend = !!sym(Gene)),
#     color = "gray70",
#     linewidth = 0.8
#   ) +
#   geom_point(
#     aes(x = pert_iname, y = !!sym(Gene), fill = !!sym(Gene)),
#     size = 6,
#     shape = 21,
#     color = "black",
#     stroke = 0.5
#   ) +
#   scale_fill_gradientn(
#     colors = colors,
#     limits = limits,  # 设置对称的limits
#     name = "Level",
#     guide = guide_colorbar(
#       title.position = "top",
#       title.hjust = 0.5,
#       barwidth = unit(0.3, "cm"),
#       barheight = unit(5, "cm"),
#       frame.colour = "black",
#       ticks.colour = "black",
#       direction = "vertical"
#     )
#   ) +
#   # Modified horizontal line at y=0 (now dashed)
#   geom_hline(yintercept = 0, color = "black", linewidth = 0.5, linetype = "dashed") +
#   
#   labs(
#     x = NULL, 
#     y = paste(Gene, "Expression Value"),
#     title = paste("Top 15 and Bottom 15", Gene, "Values in", Tissue, "Tissue")
#   ) +
#   theme_minimal() +
#   theme(
#     plot.title = element_text(size = 16, face = "bold", hjust = 0.5, margin = margin(b = 15)),
#     axis.text.y = element_text(size = 11, color = "black"),
#     axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
#     axis.title = element_text(size = 14),
#     axis.title.x = element_text(margin = margin(t = 10)),
#     legend.position = "right",
#     legend.title = element_text(size = 12),
#     legend.text = element_text(size = 10),
#     panel.grid.major = element_blank(),
#     panel.grid.minor = element_blank(),
#     plot.margin = unit(c(1, 1, 1, 1), "cm"),
#     plot.background = element_rect(fill = "white", color = NA),
#     panel.border = element_blank()
#   ) +
#   geom_vline(xintercept = 15.5, linetype = "dashed", color = "gray30") +
#   annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, 
#            fill = NA, color = "black", linewidth = 0.5)
# 
# 
# 
# 












# ggplot(result) +
#   # Add a box around the entire plot
#   geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf),
#             color = "black", fill = NA, linewidth = 0.5) +
# 
#   geom_segment(
#     aes(x = pert_iname, xend = pert_iname, y = 0, yend = !!sym(Gene)),
#     color = "gray70",
#     linewidth = 0.8
#   ) +
#   geom_point(
#     aes(x = pert_iname, y = !!sym(Gene), fill = !!sym(Gene)),
#     size = 6,
#     shape = 21,
#     color = "black",
#     stroke = 0.5
#   ) +
#   scale_fill_gradientn(
#     colors = colors,
#     name = "Level",
#     guide = guide_colorbar(
#       title.position = "top",
#       title.hjust = 0.5,
#       barwidth = unit(0.3, "cm"),
#       barheight = unit(5, "cm"),
#       frame.colour = "black",
#       ticks.colour = "black",
#       direction = "vertical"
#     )
#   ) +
#   # Modified horizontal line at y=0 (now dashed)
#   geom_hline(yintercept = 0, color = "black", linewidth = 0.5, linetype = "dashed") +
# 
#   labs(
#     x = NULL,
#     y = paste(Gene, "Expression Value"),
#     title = paste("Top 15 and Bottom 15", Gene, "Values in", Tissue, "Tissue")
#   ) +
#   theme_minimal() +
#   theme(
#     plot.title = element_text(size = 16, face = "bold", hjust = 0.5, margin = margin(b = 15)),
#     axis.text.y = element_text(size = 11, color = "black"),
#     axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
#     axis.title = element_text(size = 14),
#     axis.title.x = element_text(margin = margin(t = 10)),
#     legend.position = "right",
#     legend.title = element_text(size = 12),
#     legend.text = element_text(size = 10),
#     panel.grid.major = element_blank(),
#     panel.grid.minor = element_blank(),
#     plot.margin = unit(c(1, 1, 1, 1), "cm"),
#     plot.background = element_rect(fill = "white", color = NA),
#     panel.border = element_blank()
#   ) +
#   geom_vline(xintercept = 15.5, linetype = "dashed", color = "gray30") +
#   annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf,
#            fill = NA, color = "black", linewidth = 0.5)






# # RColorBrewer ------------------------------------------------------------
# 
# 
# 
# library(RColorBrewer)
# library(ggplot2)
# 
# # Create the plot with Spectral palette
# ggplot(result) +
#   geom_segment(
#     aes(x = pert_iname, xend = pert_iname, y = 0, yend = !!sym(Gene)),
#     color = "gray70",
#     linewidth = 0.8
#   ) +
#   geom_point(
#     aes(x = pert_iname, y = !!sym(Gene), color = !!sym(Gene)),
#     size = 6
#   ) +
#   scale_color_distiller(
#     palette = "Spectral",
#     name = "Expression Level",
#     guide = guide_colorbar(
#       title.position = "top",
#       barwidth = unit(5, "cm"),
#       barheight = unit(0.3, "cm")
#     )
#   ) +
#   labs(
#     x = NULL, 
#     y = paste(Gene, "Expression Value"),
#     title = paste("Drugs with Extreme", Gene, "Expression in", Tissue, "Tissue")
#   ) +
#   coord_flip() +
#   theme_minimal() +
#   theme(
#     plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
#     axis.text.y = element_text(size = 11, color = "black"),
#     axis.text.x = element_text(size = 12),
#     axis.title = element_text(size = 14),
#     legend.position = "top",
#     legend.title = element_text(size = 12),
#     legend.text = element_text(size = 10),
#     panel.grid.major.y = element_blank()
#   ) +
#   geom_vline(xintercept = 15.5, linetype = "dashed", color = "gray30") +
#   # Add group labels
#   annotate("text", 
#            x = c(8, 23), y = min(result[[Gene]]) * 0.9, 
#            label = c("Low PKMYT1", "High PKMYT1"), 
#            size = 5, fontface = "bold", color = "black")
