####### Dietary Restriction GSVA Score
library(dplyr)
library(stringr)
library(patchwork)
library(ggplot2)
library(knitr)
#library(msigdbr) 
#library(fgsea)
library(ggplot2)
#library(GSVA)
library(GSEABase)
library(limma)
library(data.table)
#######################
library(Hmisc)
library(tidyr)
library(tibble)
#######################
#library(GSVA)
library(GSEABase)
library(ggpubr)
library(reshape2)
library(sva)
######################################
setwd('H:\\Mendilian\\RNA_experiment\\Kidney fibrosis')
#gmt <- getGmt('./FME_signature.gmt')
################################## Neptune meta 
###################################### Calculate GSE197307
setwd('H:\\Mendilian\\RNA_experiment\\Signature_characteristics\\Neptune')
GSE197307_matrix <- read.csv('./Neptune_GSE197307_filter.csv')
rownames(GSE197307_matrix) <- GSE197307_matrix$X
GSE197307_matrix <- GSE197307_matrix[,-1]
GSE197307_matrix <- as.matrix(GSE197307_matrix)

# log2 transform
qx <- as.numeric(quantile(GSE197307_matrix, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
          (qx[6]-qx[1] > 50 && qx[2] > 0)
if (LogC) { GSE197307_matrix[which(GSE197307_matrix <= 0)] <- NaN
  GSE197307_matrix <- log2(GSE197307_matrix) }

GSE197307_matrix[which(GSE197307_matrix == 'NaN')] <- 0
GSE197307_matrix <- avereps(GSE197307_matrix)
GSE197307_matrix <- GSE197307_matrix[rowMeans(GSE197307_matrix)>0,]

###################################### Calculate GSE133288
setwd('H:\\Mendilian\\RNA_experiment\\Signature_characteristics\\Neptune')
GSE133288_matrix <- read.csv('./Neptune_GSE133288_filter.csv')
rownames(GSE133288_matrix) <- GSE133288_matrix$X
GSE133288_matrix <- GSE133288_matrix[,-1]
GSE133288_matrix <- as.matrix(GSE133288_matrix)
###################################### Combine
Come_genes <- intersect(rownames(GSE133288_matrix),rownames(GSE197307_matrix))
GSE197307_matrix <- GSE197307_matrix[Come_genes,]
GSE133288_matrix <- GSE133288_matrix[Come_genes,]
############################################################## Combine
Combine_matrix <- cbind(GSE197307_matrix,GSE133288_matrix) %>% as.data.frame()
a <- t(Combine_matrix['ACSM2A',]) %>% as.data.frame()
a$Sample <- rownames(a)
################################### FME score
setwd('H:\\Mendilian\\RNA_experiment\\Kidney fibrosis')
Neptune_score <- read.csv('./Neptune_total_FME_score.csv')
a$FME <- Neptune_score$FME_signature[match(Neptune_score$ID, a$Sample)]
setwd('H:\\CKDwork\\knokdown\\Experiment')
write.csv(a, file = 'FME_score_analysis.csv', row.names = FALSE)


