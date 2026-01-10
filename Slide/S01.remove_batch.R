##############
library(limma)
##############

rm(list = ls())  ## 魔幻操作，一键清空~
options(stringsAsFactors = F)

library("FactoMineR")
library("factoextra")
###################

pca_plot = function(dddd,ggggg){
  library("FactoMineR")
  library("factoextra")
  df.pca <- PCA(t(dddd), graph = FALSE)
  fviz_pca_ind(df.pca,
               #axes = c(2,3),
               geom.ind = "point",
               col.ind = ggggg ,
               addEllipses = TRUE,
               legend.title = "Groups"
  )
}

##### Load data
setwd('H:\\Liman_ref\\imageJ处理过')
dat <- read.csv('./ACSM2A_compare.csv')
dat <- dat[,-3]
data <- dat[,2]
data_t <- t(data)
#dim(data_t)

##### Arrange phenotype
group_list <- dat$Group
# table(group_list)
g=factor(group_list)
design=model.matrix(~g)

###### Original PCA
pca_plot(data_t,g)

################################################# Eliminate the batch effect
batch = dat$Batch
ex_b_limma <- removeBatchEffect(data_t, batch = batch, design = design)

###### Eliminate batch effect PCA
#pca_plot(ex_b_limma,g)
#scale_data <- scale(ex_b_limma[1,])
#ex_b_limma <- as.data.frame(scale_data)
data_final <- dat
data_final$ACSM2A_mean <- ex_b_limma[1,]

write.csv(data_final,file = './Remove_batch_effect_ACSM2A.csv', row.names = FALSE)

