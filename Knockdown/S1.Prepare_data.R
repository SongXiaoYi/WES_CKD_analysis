##################### 转换
library(reticulate)
library(Seurat)
library(patchwork)
library(stringr)
library(Matrix)
library(sceasy)
library(dior)
library(schard)
use_python("G:/miniconda3/python.exe")

os <- import("os")
ad <- import("anndata")
loompy <- import('loompy')
#########################
setwd('H:\\CKDwork\\knokdown\\Data')
h5ad_file <- './Mature_Full_v3.h5ad'

sc <- schard::h5ad2seurat("test.h5ad",use.raw=F)

setwd('H:\\CKDwork\\knokdown\\Data')
saveRDS(sc, file = 'Mature_Full_v3.rds')

setwd('H:\\CKDwork\\knokdown\\Data')
h5ad_file <- './Fetal_full_v3.h5ad'

sc <- schard::h5ad2seurat(h5ad_file,use.raw=F)

setwd('H:\\CKDwork\\knokdown\\Data')
saveRDS(sc, file = 'Fetal_full_v3.rds')

