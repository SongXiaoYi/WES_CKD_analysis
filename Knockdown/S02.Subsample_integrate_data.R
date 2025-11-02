#####################
library(Seurat)
library(dplyr)
library(stringr)
library(harmony)
library(patchwork)
library(ggplot2)
library(ggbump)
library(knitr)
library(purrr)
library(scater)
library(tidyr)
###########################
############# PseudoCell function
GatherData <- function(object, ...) {

  UseMethod("GatherData")

}

GatherData.Seurat <- function(object,
                              assay,
                              slot_use,
                              ...) {

  assay <- assay %||% "RNA"
  slot_use <- slot_use %||% "data"
  obj_data <- GetAssayData(
    object = object,
    assay = assay,
    slot = slot_use
  ) %>%
    as.matrix()
  return(obj_data)
}

PseudoCell  <- function(object,
                        assay_use = NULL,
                        slot_use = NULL,
                        cluster_use =NULL,
                        pseudocell.size  =NULL){
  message("tips: 
  Cluster_use : one col in metadata
  pseudocell.size : how many cell will be pseudo")

  Inter<- GatherData(object = object,
                     assay = assay_use,
                     slot_use = slot_use) 
  Inter[Inter<0]=0
  idd<-object@meta.data
  Inter.id<-cbind(rownames(idd),as.vector(idd[,cluster_use]))

  rownames(Inter.id)<-rownames(idd)
  colnames(Inter.id)<-c("CellID","Celltype")

  Inter.id<-as.data.frame(Inter.id)
  Inter1<-Inter[,Inter.id$CellID]
  Inter<-as.matrix(Inter1)
  pseudocell.size = pseudocell.size ## 10 test
  new_ids_list = list()
  Inter.id$Celltype <- as.factor(Inter.id$Celltype)
  for (i in 1:length(levels(Inter.id$Celltype))) {
    cluster_id = levels(Inter.id$Celltype)[i]
    cluster_cells <- rownames(Inter.id[Inter.id$Celltype == cluster_id,])
    cluster_size <- length(cluster_cells)       
    pseudo_ids <- floor(seq_along(cluster_cells)/pseudocell.size)
    pseudo_ids <- paste0(cluster_id, "_Cell", pseudo_ids)
    names(pseudo_ids) <- sample(cluster_cells)  
    new_ids_list[[i]] <- pseudo_ids     
  }

  new_ids <- unlist(new_ids_list)
  new_ids <- as.data.frame(new_ids)
  new_ids_length <- table(new_ids)

  new_colnames <- rownames(new_ids)  ###add
  all.data<-Inter[,as.character(new_colnames)] ###add
  all.data <- t(all.data)###add

  new.data<-aggregate(list(all.data[,1:length(all.data[1,])]),
                      list(name=new_ids[,1]),FUN=mean)
  rownames(new.data)<-new.data$name
  new.data<-new.data[,-1]

  new_ids_length<-as.matrix(new_ids_length)##
  short<-which(new_ids_length< pseudocell.size -1 )##
  new_good_ids<-as.matrix(new_ids_length[-short,])##
  result<-t(new.data)[,rownames(new_good_ids)]
  rownames(result)<-rownames(Inter)

  newobject <- CreateSeuratObject(result)
  newobject@misc$idtrans <- new_ids
  newobject@commands$PseudoCell <- LogSeuratCommand(newobject, return.command = TRUE)
  return(newobject)
}

Provide_tag <- function(x){
  cell_list <- colnames(x)
  cell_ann <- as.data.frame(cell_list)
  cell_ann[,2:4] <- str_split_fixed(cell_ann$cell_list, "_", 3)
  colnames(cell_ann)[2:4] <- c('celltype','age','id')
  cell_ann$age <- as.numeric(cell_ann$age)
  cell_ann$group <- cell_ann$age
  cell_ann[which(cell_ann$age <= 6),'group'] <- 0
  cell_ann[which(cell_ann$age > 6),'group'] <- 1
  cell_ann$group
}

Seurat_preprocessing <- function(seurat, project = "Scissor_Single_Cell", 
                                 normalization.method = "LogNormalize", scale.factor = 10000,
                                 selection.method = "vst", resolution = 0.6,
                                 dims_Neighbors = 1:10, dims_TSNE = 1:10, dims_UMAP = 1:10,
                                 verbose = TRUE){
    data <- NormalizeData(object = seurat, normalization.method = normalization.method, scale.factor = scale.factor, verbose = verbose)
    data <- FindVariableFeatures(object = data, selection.method = selection.method, verbose = verbose)
    data <- ScaleData(object = data, verbose = verbose)
    data <- RunPCA(object = data, features = VariableFeatures(data), verbose = verbose)
    data <- FindNeighbors(object = data, dims = dims_Neighbors, verbose = verbose)
    data <- FindClusters( object = data, resolution = resolution, verbose = verbose)
    data <- RunTSNE(object = data, dims = dims_TSNE)
    data <- RunUMAP(object = data, dims = dims_UMAP, verbose = verbose)

    return(data)
}
################################################### Heree begin
########################################
setwd('H:\\CKDwork\\knokdown\\Data')
pbmc <- readRDS('./Mature_Full_v3.rds')
pseudo_Bulk <- PseudoCell(pbmc, "RNA", "counts", "celltype", 20)











