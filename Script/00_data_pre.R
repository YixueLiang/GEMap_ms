#清空数据
rm(list = ls())


# depmap_23q2 -------------------------------------------------------------

# CRISPRGeneEffect.csv
# Model.csv
# OmicsExpressionProteinCodingGenesTPMLogp1.csv
# OmicsSomaticMutations.csv
# OmicsCNGene.csv

### CRISPR knockout screen across 1000+ cell lines
df_crispr = data.table::fread('inputdata/DepMap_Public_23Q2/CRISPRGeneEffect.csv')
df_cell = data.table::fread('inputdata/DepMap_Public_23Q2/Model.csv')
df_exp = data.table::fread('inputdata/DepMap_Public_23Q2/OmicsExpressionProteinCodingGenesTPMLogp1.csv')
df_mut = data.table::fread('inputdata/DepMap_Public_23Q2/OmicsSomaticMutations.csv')
df_cn = data.table::fread('inputdata/DepMap_Public_23Q2/OmicsCNGene.csv')


colnames(df_crispr) <- gsub("\\s+\\(\\d+\\)","",colnames(df_crispr))
colnames(df_exp) <- gsub("\\s+\\(\\d+\\)","",colnames(df_exp))
df_crispr <- as.data.frame(df_crispr)
df_exp <- as.data.frame(df_exp)
rownames(df_crispr) <- df_crispr$ModelID
rownames(df_exp) <- df_exp$V1
df_crispr <- df_crispr[,-1]
df_exp <- df_exp[,-1]

df_cell <- as.data.frame(df_cell)
rownames(df_cell) <- df_cell[,4]


colnames(df_cn) <- gsub("\\s+\\(\\d+\\)","",colnames(df_cn))
df_cn <- as.data.frame(df_cn)
rownames(df_cn) <- df_cn$V1
df_cn <- df_cn[,-1]
colnames(df_cn) <- gsub("\\.", "-", colnames(df_cn))

# 使用 names() 函数将新的列名应用到数据框中
df_mut$MutationType <- paste(df_mut$HugoSymbol, df_mut$ProteinChange, sep = "_")
colnames(df_mut)[55] <- "DepMap_ID"


save(df_cell, df_cn, df_crispr, df_exp, df_mut,
     file = "outputdata/df_5.Rdata")


# load("outputdata/df_5.Rdata")

unique(df_cell$OncotreeLineage)
# df_cell = data.table::fread('inputdata/DepMap_Public_23Q2/Model.csv')
# df_cell <- as.data.frame(df_cell)
# rownames(df_cell) <- df_cell$StrippedCellLineName

df_cell$OncotreeLineage <- gsub("Ovary/Fallopian Tube", "Ovary", df_cell$OncotreeLineage)
df_cell$OncotreeLineage <- gsub("Bladder/Urinary Tract", "Bladder", df_cell$OncotreeLineage)
df_cell$OncotreeLineage <- gsub("CNS/Brain", "CNS", df_cell$OncotreeLineage)
df_cell$OncotreeLineage <- gsub("Vulva/Vagina", "Vulva", df_cell$OncotreeLineage)
df_cell$OncotreeLineage <- gsub("Esophagus/Stomach", "Esophagus and Stomach", df_cell$OncotreeLineage)


save(df_cell,file = "DEMap_Data/CRISPR_sanger_avana_gecko.Rdata")

save(df_cell, df_cn, df_crispr, df_exp, df_mut,
     file = "DEMap_Data/CRISPR_DATA.Rdata")


load("outputdata_lincs/03_drug_data_list_20.RData")
merged_df <- do.call(rbind, drug_data_list)

save(merged_df, file = "outputdata/17_merged_df.Rdata")


# sanger ------------------------------------------------------------------


df_sanger = data.table::fread('inputdata/Sanger_CRISPR/Sanger CRISPR (Project Score, CERES) Files/gene_effect.csv')

df_sanger <- as.data.frame(df_sanger)
colnames(df_sanger) <- gsub("\\s+\\(\\d+\\)","",colnames(df_sanger))
rownames(df_sanger) <- df_sanger[,1]
df_sanger <- df_sanger[,-1]

common <- intersect(rownames(df_sanger), rownames(df_cell))
df_cell2 <- df_cell[common, ]
df_sanger <- df_sanger[common, ]
rownames(df_sanger) <- df_cell2$StrippedCellLineName

save(df_dep, df_sanger,
     file = "outputdata_Sanger/df_2_cers.Rdata")


# pubmed ------------------------------------------------------------------

load("outputdata_pubmed/merged_df.RData")
screen_info <- read.delim("inputdata/PubMed/screen_information.tsv", header = TRUE, sep = "\t")
unique(screen_info$PHENOTYPE)
unique(screen_info$LIBRARY_METHODOLOGY)
unique(screen_info$SCREEN_TYPE)
unique(screen_info$EXPERIMENTAL_SETUP)
# 1375-1260
screen_info <- screen_info[screen_info$PHENOTYPE %in% "Viability" & 
                             screen_info$LIBRARY_METHODOLOGY %in% "Knockout" &
                             screen_info$SCREEN_TYPE %in% "Negative selection" &
                             screen_info$EXPERIMENTAL_SETUP  %in% "Timecourse", ]

merged_df <- as.data.frame(merged_df)
colnames(merged_df) <- gsub("\\.BCSCORE", "", colnames(merged_df))
rownames(merged_df) <- merged_df$HGNC
merged_df <- merged_df[, -c(1, 2)]
merged_df <- as.data.frame(t(merged_df))

merged_df <- cbind(Cell_id = rownames(merged_df), merged_df)
COM <- intersect(rownames(merged_df), screen_info$EWHA_ID)
merged_df <- merged_df[COM, ]
merged_df$Cell_id <- screen_info$CELL_LINE[match(merged_df$Cell_id, screen_info$EWHA_ID)]
length(unique(merged_df$Cell_id))

# save(merged_df, file = "outputdata_pubmed/merged_df_cell.RData")
# save(screen_info, file = "outputdata_pubmed/screen_info_1260.RData")
# 

save(merged_df, screen_info, PubMed_info, file = "DEMap_Data/PubMed_Data.RData")


# AVANA -------------------------------------------------------------------

df_Avana = data.table::fread('inputdata/Achilles Avana /ceresgeneeffects.csv')
df_Avana_cell <- data.table::fread('inputdata/Achilles Avana /celllineinfo.csv')

df_Avana <- as.data.frame(df_Avana)
df_Avana_cell <- as.data.frame(df_Avana_cell)
rownames(df_Avana) <- df_Avana$V1
df_Avana <- df_Avana[, df_Avana_cell$Name]

df_Avana_cell$CellName <- sub("_.*", "", df_Avana_cell$Name)
# length(intersect(df_Avana_cell$CellName, df_cell$StrippedCellLineName))

colnames(df_Avana) <- df_Avana_cell$CellName
df_Avana <- as.data.frame(t(df_Avana))

rownames(df_Avana_cell) <- df_Avana_cell$CellName

save(df_Avana, df_Avana_cell, file = "DEMap_Data/Avana_DATA.Rdata")



# GeCKO -------------------------------------------------------------------

df_GeCKO = data.table::fread('inputdata/Achilles GeCKO/gene_effect.csv')
df_GeCKO <- as.data.frame(df_GeCKO)

rownames(df_GeCKO) <- df_GeCKO$V1
df_GeCKO <- df_GeCKO[,-1]

CellName <- sub("_.*", "", colnames(df_GeCKO))
Tissue <- sub("^[^_]*_", "", colnames(df_GeCKO))

colnames(df_GeCKO) <- sub("_.*", "", colnames(df_GeCKO))
df_GeCKO <- as.data.frame(t(df_GeCKO))


df_GeCKO_cell <- data.frame(CellName = CellName, 
                            Lineage = Tissue)

save(df_GeCKO, df_GeCKO_cell, file = "DEMap_Data/GeCKO_DATA.Rdata")






