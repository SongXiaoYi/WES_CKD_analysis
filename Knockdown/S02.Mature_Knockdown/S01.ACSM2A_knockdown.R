####################################
library(dplyr)
library(Seurat)
library(scTenifoldKnk)
library(ggplot2)
library(ggrepel)  # 用于防止标签重叠
library(data.table)
##########################################
setwd('H:\\CKDwork\\knokdown\\Data')
load('./Mature_Pseudo_Bulk.Rdata')

pseudo_Bulk__mature_normalize <- FindVariableFeatures(object = pseudo_Bulk__mature_normalize, selection.method = "vst", nfeatures = 10000)   #这里设置为前2500

Select_Features <- VariableFeatures(pseudo_Bulk__mature_normalize)

setwd('H:\\CKDwork\\knokdown\\Experiment')

Gene_list <- fread('./Gene_list.txt')

Gene_list <- as.character(Gene_list$Gene_id)

Final_gene_select <- unique(c(Select_Features, Gene_list))

pseudo_Bulk__mature_normalize_copy <- pseudo_Bulk__mature_normalize

pseudo_Bulk__mature_normalize_copy <- pseudo_Bulk__mature_normalize_copy[rownames(pseudo_Bulk__mature_normalize_copy) %in% Final_gene_select]

pseudo_Bulk__mature_normalize_copy <- NormalizeData(object = pseudo_Bulk__mature_normalize_copy, normalization.method = "LogNormalize", scale.factor = 10000)

countMatrix <- pseudo_Bulk__mature_normalize_copy@assays$RNA$data

#countMatrix <- pseudo_Bulk__mature_normalize@assays$RNA$data

result <- scTenifoldKnk(countMatrix = countMatrix, 
                        gKO = 'ACSM2A', #需要敲除的基因
                        qc = FALSE,#是否进行QC
                        qc_mtThreshold = 0.2,#mt阈值
                        qc_minLSize = 1000,#文库阈值(细胞测到的基因总数)
                        nc_nNet = 10, #子网络数量
                        nc_nCells = 500, #每个网络中随机抽取的细胞数
                        nc_nComp = 3 #PCA 的主成分数量
                        )


