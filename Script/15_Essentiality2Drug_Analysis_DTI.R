
#清空数据
rm(list = ls())


# Load the data
load("GEMap_Outputdata/00_interactions2.Rdata")

interactions_approved <- interactions2[interactions2$approved %in% "TRUE", ]

DGIDB_Approved_drugs <- unique(interactions_approved$drug_name)


# 加载必要的包
library(igraph)

# 假设 interactions2 是包含药物和靶点相互作用的数据框
# 筛选出与靶点 CDK2 相关的相互作用
target <- "CDK2"
interactions_CDK2 <- interactions2[interactions2$gene_claim_name %in% target, ]
interactions_CDK2_APPROVED <- interactions_CDK2[interactions_CDK2$gene_claim_name %in% "TRUE", ]
unique(interactions_CDK2$drug_name)
# 检查筛选后的数据中是否有 NA 值
if (any(is.na(interactions_CDK2$drug_claim_name)) | any(is.na(interactions_CDK2$interaction_type))) {
  stop("数据中存在 NA 值，请检查并清理数据。")
}

# 删除 interactions_CDK2 中包含 NA 的行
interactions_CDK2 <- na.omit(interactions_CDK2)

# 筛选出药物和靶点的相互作用数据
edges <- data.frame(
  from = interactions_CDK2$drug_claim_name,
  to = target,
  interaction_type = interactions_CDK2$interaction_type,
  interaction_score = interactions_CDK2$interaction_score
)

# 根据 interaction_type 设置边的颜色
edges$color <- ifelse(
  edges$interaction_type %in% c("agonist", "activator", "positive modulator", "vaccine"),
  "#D73027",  # 激动作用为红色
  "#4575B4"  # 其他为蓝色
)

# 根据 interaction_score 设置边的权重
# 假设 interaction_score 值越大，相互作用越强，边的长度越短
edges$interaction_score <- as.numeric(edges$interaction_score)
edges$weight <- 1 / edges$interaction_score  # 权重与 interaction_score 成反比
# edges$weight <- scale(edges$interaction_score)

# 创建网络图
# 节点包括药物和靶点
nodes <- data.frame(
  name = unique(c(edges$from, edges$to)),
  type = ifelse(unique(c(edges$from, edges$to)) == target, "target", "drug")
)

# 创建图对象
g <- graph_from_data_frame(edges, directed = TRUE, vertices = nodes)

# 设置节点颜色
V(g)$color <- ifelse(V(g)$type == "target", "orange", "lightblue")  # 靶点为橙色，药物为浅蓝色

# 设置节点大小
V(g)$size <- ifelse(V(g)$type == "target", 20, 15)  # 靶点较大，药物较小


# 计算布局并保存
layout_fixed <- layout_with_fr(g, weights = E(g)$weight)  # 使用 Fruchterman-Reingold 布局

# 绘制网络图（使用固定的布局）
plot(g, 
     edge.arrow.size = 0.5,  # 箭头大小
     vertex.label.cex = 0.8,  # 节点标签大小
     vertex.label.color = "black",  # 节点标签颜色
     edge.color = E(g)$color,  # 边的颜色
     layout = layout_fixed,  # 使用固定的布局
     main = paste("Drug-Target Network for", target)  # 图标题
)

