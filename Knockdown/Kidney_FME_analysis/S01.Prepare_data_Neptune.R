####### Dietary Restriction GSVA Score
library(dplyr)
library(stringr)
library(patchwork)
library(ggplot2)
library(knitr)
library(msigdbr) 
library(fgsea)
library(ggplot2)
library(GSVA)
library(GSEABase)
library(limma)
library(data.table)
#######################
library(Hmisc)
library(tidyr)
library(tibble)
#######################
library(GSVA)
library(GSEABase)
library(ggpubr)
library(reshape2)
library(sva)
######################################
setwd('G:\\Mendilian\\RNA_experiment\\Data\\Bulkseq\\海王星信息整理\\Datasets\\GSE182380')
pbmc <- fread('GSE182380-normalize_data.csv') %>% as.data.frame()
Names <- pbmc$V1
pbmc <- as.matrix(pbmc[,2:200])
rownames(pbmc) <- Names
pbmc = avereps(pbmc)
pbmc = pbmc[rowMeans(pbmc)>0,]
pbmc = normalizeBetweenArrays(pbmc)
##################################
setwd('G:\\Mendilian\\RNA_experiment\\Kidney fibrosis')
gmt <- getGmt('./FME_signature.gmt')
##################################
gs.exp <- gsva(pbmc,gmt,method = "ssgsea",kcdf = "Gaussian", min.sz = 10)
Total_score <- t(gs.exp)
Total_score <- as.data.frame(Total_score)
Total_score$ID <- rownames(Total_score)
Total_score$ACSM2A <- pbmc['ACSM2A',]
Total_score <- Total_score[,c(2,3,1)]
######################################################
setwd('G:\\Mendilian\\RNA_experiment\\Data\\Bulkseq\\海王星信息整理\\Datasets\\Integrate_cohort')
meta <- fread('./Integrate_meta.csv')
Total_score$EGFR <- meta$EGFR[match(Total_score$ID,meta$ID_match)]
setwd('G:\\CKDwork\\knokdown\\Experiment\\ACSM2A_Mature')
write.csv(Total_score,file = 'Neptune_total_FME_score.csv',row.names = FALSE)







