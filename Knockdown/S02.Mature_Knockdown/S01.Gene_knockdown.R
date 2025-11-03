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

setwd('H:\\CKDwork\\knokdown\\Experiment\\ACSM2A')
save(result, file = 'ACSM2A_result_G10000_C2000.Rdata')

########################################################### TPTE2
result <- scTenifoldKnk(countMatrix = countMatrix, 
                        gKO = 'TPTE2', #需要敲除的基因
                        qc = FALSE,#是否进行QC
                        qc_mtThreshold = 0.2,#mt阈值
                        qc_minLSize = 1000,#文库阈值(细胞测到的基因总数)
                        nc_nNet = 10, #子网络数量
                        nc_nCells = 500, #每个网络中随机抽取的细胞数
                        nc_nComp = 3 #PCA 的主成分数量
                        )

setwd('H:\\CKDwork\\knokdown\\Experiment\\TPTE2')
save(result, file = 'TPTE2_result_G10000_C2000.Rdata')

########################################################### OR2L8
result <- scTenifoldKnk(countMatrix = countMatrix, 
                        gKO = 'OR2L8', #需要敲除的基因
                        qc = FALSE,#是否进行QC
                        qc_mtThreshold = 0.2,#mt阈值
                        qc_minLSize = 1000,#文库阈值(细胞测到的基因总数)
                        nc_nNet = 10, #子网络数量
                        nc_nCells = 500, #每个网络中随机抽取的细胞数
                        nc_nComp = 3 #PCA 的主成分数量
                        )

setwd('H:\\CKDwork\\knokdown\\Experiment\\OR2L8')
save(result, file = 'OR2L8_result_G10000_C2000.Rdata')
gc()
########################################################### ZNF782
result <- scTenifoldKnk(countMatrix = countMatrix, 
                        gKO = 'ZNF782', #需要敲除的基因
                        qc = FALSE,#是否进行QC
                        qc_mtThreshold = 0.2,#mt阈值
                        qc_minLSize = 1000,#文库阈值(细胞测到的基因总数)
                        nc_nNet = 10, #子网络数量
                        nc_nCells = 500, #每个网络中随机抽取的细胞数
                        nc_nComp = 3 #PCA 的主成分数量
                        )

setwd('H:\\CKDwork\\knokdown\\Experiment\\ZNF782')
save(result, file = 'ZNF782_result_G10000_C2000.Rdata')
gc()
########################################################### ZNF510
result <- scTenifoldKnk(countMatrix = countMatrix, 
                        gKO = 'ZNF510', #需要敲除的基因
                        qc = FALSE,#是否进行QC
                        qc_mtThreshold = 0.2,#mt阈值
                        qc_minLSize = 1000,#文库阈值(细胞测到的基因总数)
                        nc_nNet = 10, #子网络数量
                        nc_nCells = 500, #每个网络中随机抽取的细胞数
                        nc_nComp = 3 #PCA 的主成分数量
                        )

setwd('H:\\CKDwork\\knokdown\\Experiment\\ZNF510')
save(result, file = 'ZNF782_result_G10000_C2000.Rdata')
gc()
########################################################### CYFIP1
result <- scTenifoldKnk(countMatrix = countMatrix, 
                        gKO = 'CYFIP1', #需要敲除的基因
                        qc = FALSE,#是否进行QC
                        qc_mtThreshold = 0.2,#mt阈值
                        qc_minLSize = 1000,#文库阈值(细胞测到的基因总数)
                        nc_nNet = 10, #子网络数量
                        nc_nCells = 500, #每个网络中随机抽取的细胞数
                        nc_nComp = 3 #PCA 的主成分数量
                        )

setwd('H:\\CKDwork\\knokdown\\Experiment\\CYFIP1')
save(result, file = 'CYFIP1_result_G10000_C2000.Rdata')
gc()
########################################################### SLC35G6
result <- scTenifoldKnk(countMatrix = countMatrix, 
                        gKO = 'SLC35G6', #需要敲除的基因
                        qc = FALSE,#是否进行QC
                        qc_mtThreshold = 0.2,#mt阈值
                        qc_minLSize = 1000,#文库阈值(细胞测到的基因总数)
                        nc_nNet = 10, #子网络数量
                        nc_nCells = 500, #每个网络中随机抽取的细胞数
                        nc_nComp = 3 #PCA 的主成分数量
                        )

setwd('H:\\CKDwork\\knokdown\\Experiment\\SLC35G6')
save(result, file = 'SLC35G6_result_G10000_C2000.Rdata')
gc()
########################################################### MUC22
result <- scTenifoldKnk(countMatrix = countMatrix, 
                        gKO = 'MUC22', #需要敲除的基因
                        qc = FALSE,#是否进行QC
                        qc_mtThreshold = 0.2,#mt阈值
                        qc_minLSize = 1000,#文库阈值(细胞测到的基因总数)
                        nc_nNet = 10, #子网络数量
                        nc_nCells = 500, #每个网络中随机抽取的细胞数
                        nc_nComp = 3 #PCA 的主成分数量
                        )

setwd('H:\\CKDwork\\knokdown\\Experiment\\MUC22')
save(result, file = 'MUC22_result_G10000_C2000.Rdata')
gc()



