####################################
library(dplyr)
library(Seurat)
library(scTenifoldKnk)
library(ggplot2)
library(ggrepel)  # 用于防止标签重叠
library(data.table)
##########################################
setwd('H:\\CKDwork\\knokdown\\Data')
load('./Fetal_Pseudo_Bulk.Rdata')

pseudo_Bulk__fetal_normalize <- FindVariableFeatures(object = pseudo_Bulk__fetal_normalize, selection.method = "vst", nfeatures = 10000)   #这里设置为前2500

Select_Features <- VariableFeatures(pseudo_Bulk__fetal_normalize)

setwd('H:\\CKDwork\\knokdown\\Experiment')

Gene_list <- fread('./Gene_list.txt')

Gene_list <- as.character(Gene_list$Gene_id)

Final_gene_select <- unique(c(Select_Features, Gene_list))

pseudo_Bulk__fetal_normalize_copy <- pseudo_Bulk__fetal_normalize

pseudo_Bulk__fetal_normalize_copy <- pseudo_Bulk__fetal_normalize_copy[rownames(pseudo_Bulk__fetal_normalize_copy) %in% Final_gene_select]

pseudo_Bulk__fetal_normalize_copy <- NormalizeData(object = pseudo_Bulk__fetal_normalize_copy, normalization.method = "LogNormalize", scale.factor = 10000)

countMatrix <- pseudo_Bulk__fetal_normalize_copy@assays$RNA$data

#countMatrix <- pseudo_Bulk__fetal_normalize@assays$RNA$data

########################################################### ACSM2A
result <- scTenifoldKnk(countMatrix = countMatrix, 
                        gKO = 'ACSM2A', #需要敲除的基因
                        qc = FALSE,#是否进行QC
                        qc_mtThreshold = 0.2,#mt阈值
                        qc_minLSize = 1000,#文库阈值(细胞测到的基因总数)
                        nc_nNet = 10, #子网络数量
                        nc_nCells = 500, #每个网络中随机抽取的细胞数
                        nc_nComp = 3 #PCA 的主成分数量
                        )

setwd('H:\\CKDwork\\knokdown\\Experiment\\ACSM2A_Fetal')
save(result, file = 'ACSM2A_result_G10000_C2000.Rdata')
