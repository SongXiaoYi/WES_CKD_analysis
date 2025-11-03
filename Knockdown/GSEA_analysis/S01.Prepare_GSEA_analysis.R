####################
library(Seurat)
library(dplyr)
library(stringr)
library(harmony)
library(patchwork)
library(ggplot2)
library(ggbump)
library(knitr)
library(clusterProfiler)
library(org.Mm.eg.db)
library(GseaVis)
library(dplyr)
library(stringr)
library(patchwork)
library(knitr)
library(msigdbr) 
library(fgsea)
library(ggplot2)
library(GSVA)
library(GSEABase)
library(limma)
library(clusterProfiler)
library(DESeq2)
library(ggplotify)
library(cowplot)
library(edgeR)
library(tinyarray)
library(gridExtra)
library(Seurat)
library(pheatmap)
library(ggplot2)
library(reshape2)
library(enrichplot)
###########################
setwd('G:\\WuXiang\\Experiment\\WGCNA\\Data')
############################################### Matrix
pbmc <- read.csv('./Raw_count.csv')
rownames(pbmc) <- pbmc[,1]
Rownames <- pbmc[,1]
pbmc <- pbmc[,-1]
rownames(pbmc) <- Rownames
pbmc <- avereps(pbmc)
pbmc <- pbmc[rowMeans(pbmc)>0,]
################################################ read meta
pbmc_meta <- read.csv('./WGCNA_meta.csv')
pbmc <- pbmc[,pbmc_meta$Sample]
data <- pbmc
################################################ DEG
pbmc_meta$Compare <- pbmc_meta$subtype
pbmc_meta[which(pbmc_meta$subtype != 'DandG'),'Compare'] <- 'Control'
description <- factor(pbmc_meta$Compare, levels = c('Control','DandG'))
matrix <- as.data.frame(pbmc)
design <- model.matrix(~description + 0, matrix)
colnames(design) <- c('Control','DandG')

dge <- DGEList(counts = pbmc)
v <- voom(dge, design, normalize = 'quantile')
fit <- lmFit(v, design)

cont.matrix <- makeContrasts(DandG-Control, levels = design)
fit2 = contrasts.fit(fit, cont.matrix)
fit2 = eBayes(fit2)

DEG = topTable(fit2, adjust.method = "fdr", sort.by = "B", number = nrow(matrix))
DEG = na.omit(DEG)

DEG <- data.frame(Gene = row.names(DEG),DEG)
#diff <- subset(DEG,DEG$P.Val < 0.05 & abs(DEG$logFC) > 2)
#diff$Trend <- "up"
#diff$Trend[diff$logFC < 0] <- "down" 
#diff <- diff[diff$Trend == 'up',]


h.gmt <- msigdbr(species = 'human', category = "H") %>%
            dplyr::select(gs_name, gene_symbol)

############## inflammatory Response
GeneSet <- h.gmt[which(h.gmt$gs_name == 'HALLMARK_INFLAMMATORY_RESPONSE'),'gene_symbol']
GeneSet <- GeneSet$gene_symbol

#################### Heatmap
GeneSet <- GeneSet[GeneSet %in% rownames(pbmc)]
gene_expression <- pbmc[GeneSet,]

dat <- gene_expression
n=t(scale(t(log(dat+1)))) # 'scale'可以对log-ratio数值进行归一化
n[n>2]=2
n[n< -2]= -2

n1 <- n
group_list <- pbmc_meta$subtype
ac=data.frame(g=group_list)
rownames(ac)=colnames(n1) 
#n1 <- melt(n1)
#p1 <- ggplot(n1,aes(x = Var1, y = Var2, fill = value)) + geom_tile(color = 'black') + coord_fixed()
p1 <- pheatmap(n1,show_colnames =F,show_rownames = F,scale = "none",cluster_cols = FALSE,cluster_rows = FALSE,
        annotation_col=ac,cellwidth = 13.5, cellheight = 0.5)

setwd('G:\\WuXiang\\Experiment\\GSEA')
pdf(file="./INFLAMMATORY_RESPONSE_H.pdf",width= 4,height= 2.6)
p1
dev.off()

#####################################################
############## INTERFERON_GAMMA_RESPONSE
GeneSet <- h.gmt[which(h.gmt$gs_name == 'HALLMARK_INTERFERON_GAMMA_RESPONSE'),'gene_symbol']
GeneSet <- GeneSet$gene_symbol

#################### Heatmap
GeneSet <- GeneSet[GeneSet %in% rownames(pbmc)]
gene_expression <- pbmc[GeneSet,]

dat <- gene_expression
n=t(scale(t(log(dat+1)))) # 'scale'可以对log-ratio数值进行归一化
n[n>2]=2
n[n< -2]= -2

n1 <- n
group_list <- pbmc_meta$subtype
ac=data.frame(g=group_list)
rownames(ac)=colnames(n1) 
#n1 <- melt(n1)
#p1 <- ggplot(n1,aes(x = Var1, y = Var2, fill = value)) + geom_tile(color = 'black') + coord_fixed()
p1 <- pheatmap(n1,show_colnames =F,show_rownames = F,scale = "none",cluster_cols = FALSE,cluster_rows = FALSE,
        annotation_col=ac,cellwidth = 13.5, cellheight = 0.5)

setwd('G:\\WuXiang\\Experiment\\GSEA')
pdf(file="./INTERFERON_GAMMA_H.pdf",width= 4,height= 2.6)
p1
dev.off()

#####################################################
############## HALLMARK_INTERFERON_ALPHA_RESPONSE
GeneSet <- h.gmt[which(h.gmt$gs_name == 'HALLMARK_INTERFERON_ALPHA_RESPONSE'),'gene_symbol']
GeneSet <- GeneSet$gene_symbol

#################### Heatmap
GeneSet <- GeneSet[GeneSet %in% rownames(pbmc)]
gene_expression <- pbmc[GeneSet,]

dat <- gene_expression
n=t(scale(t(log(dat+1)))) # 'scale'可以对log-ratio数值进行归一化
n[n>2]=2
n[n< -2]= -2

n1 <- n
group_list <- pbmc_meta$subtype
ac=data.frame(g=group_list)
rownames(ac)=colnames(n1) 
#n1 <- melt(n1)
#p1 <- ggplot(n1,aes(x = Var1, y = Var2, fill = value)) + geom_tile(color = 'black') + coord_fixed()
p1 <- pheatmap(n1,show_colnames =F,show_rownames = F,scale = "none",cluster_cols = FALSE,cluster_rows = FALSE,
        annotation_col=ac,cellwidth = 13.5, cellheight = 0.5)

setwd('G:\\WuXiang\\Experiment\\GSEA')
pdf(file="./INTERFERON_ALPHA_H.pdf",width= 4,height= 2.6)
p1
dev.off()

#####################################################
############## HALLMARK_INTERFERON_ALPHA_RESPONSE
GeneSet <- h.gmt[which(h.gmt$gs_name == 'HALLMARK_IL6_JAK_STAT3_SIGNALING'),'gene_symbol']
GeneSet <- GeneSet$gene_symbol

#################### Heatmap
GeneSet <- GeneSet[GeneSet %in% rownames(pbmc)]
gene_expression <- pbmc[GeneSet,]

dat <- gene_expression
n=t(scale(t(log(dat+1)))) # 'scale'可以对log-ratio数值进行归一化
n[n>2]=2
n[n< -2]= -2

n1 <- n
group_list <- pbmc_meta$subtype
ac=data.frame(g=group_list)
rownames(ac)=colnames(n1) 
#n1 <- melt(n1)
#p1 <- ggplot(n1,aes(x = Var1, y = Var2, fill = value)) + geom_tile(color = 'black') + coord_fixed()
p1 <- pheatmap(n1,show_colnames =F,show_rownames = F,scale = "none",cluster_cols = FALSE,cluster_rows = FALSE,
        annotation_col=ac,cellwidth = 13.5, cellheight = 0.5)

setwd('G:\\WuXiang\\Experiment\\GSEA')
pdf(file="./IL6_JAK_H.pdf",width= 4,height= 2.6)
p1
dev.off()


##################################### GSEA
IR <- DEG[!duplicated(DEG$Gene),c('logFC','adj.P.Val')]
#########################
library(org.Hs.eg.db)
gene <- bitr(rownames(IR),     #转换的列是nrDEG的列名
             fromType = "SYMBOL",     #需要转换ID类型
             toType =  "ENTREZID",    #转换成的ID类型
             OrgDb = org.Hs.eg.db)    #对应的物种，小鼠的是org.Mm.eg.db
#让基因名、ENTREZID、logFC对应起来
gene$logFC <- IR$logFC[match(gene$SYMBOL,rownames(IR))]
head(gene)
geneList=gene$logFC
names(geneList)=gene$SYMBOL 
geneList <- geneList[!duplicated(names(geneList))]
#按照logFC的值来排序geneList
geneList=sort(geneList,decreasing = T)
head(geneList)

Plot_GeneSet <- h.gmt
colnames(Plot_GeneSet) <- c('term','gene')

set.seed(1234)
egmt <- GSEA(geneList, TERM2GENE = Plot_GeneSet, 
             minGSSize = 10,maxGSSize = 10000,
             pvalueCutoff = 1,pAdjustMethod="fdr", 
                    seed=TRUE, by="fgsea")

gesa_res <- as.data.frame(egmt@result)
###############################################
#gseaplot2(egmt, "HALLMARK_INFLAMMATORY_RESPONSE", color = "firebrick",
#            rel_heights=c(1, .2, .6), subplots=1:2)
############################################### Inflammatory Res
p1 <- gseaplot2(egmt,#数据
            "HALLMARK_INFLAMMATORY_RESPONSE",#画哪一列的信号通路
            #title = "Prion disease",#标题
            base_size = 15,#字体大小
            color = "green",#线条的颜色
            #pvalue_table = TRUE,#加不加p值
            ES_geom="line",subplots=1:2)#是用线，还是用d点

setwd('G:\\WuXiang\\Experiment\\GSEA')
pdf(file="./InflammatoryRes.pdf",width= 3.25,height= 2.5)
p1
dev.off()

############################################### IL6_JAK_STAT3_SIGNALING
p1 <- gseaplot2(egmt,#数据
            "HALLMARK_IL6_JAK_STAT3_SIGNALING",#画哪一列的信号通路
            #title = "Prion disease",#标题
            base_size = 15,#字体大小
            color = "green",#线条的颜色
            #pvalue_table = TRUE,#加不加p值
            ES_geom="line",subplots=1:2)#是用线，还是用d点

setwd('G:\\WuXiang\\Experiment\\GSEA')
pdf(file="./IL6_JAK.pdf",width= 3.25,height= 2.5)
p1
dev.off()

############################################### INTERFERON_ALPHA_RESPONSE
p1 <- gseaplot2(egmt,#数据
            "HALLMARK_INTERFERON_ALPHA_RESPONSE",#画哪一列的信号通路
            #title = "Prion disease",#标题
            base_size = 15,#字体大小
            color = "green",#线条的颜色
            #pvalue_table = TRUE,#加不加p值
            ES_geom="line",subplots=1:2)#是用线，还是用d点

setwd('G:\\WuXiang\\Experiment\\GSEA')
pdf(file="./INTERFERON_ALPHA.pdf",width= 3.25,height= 2.5)
p1
dev.off()

############################################### HALLMARK_INTERFERON_GAMMA_RESPONSE
p1 <- gseaplot2(egmt,#数据
            "HALLMARK_INTERFERON_GAMMA_RESPONSE",#画哪一列的信号通路
            #title = "Prion disease",#标题
            base_size = 15,#字体大小
            color = "green",#线条的颜色
            #pvalue_table = TRUE,#加不加p值
            ES_geom="line",subplots=1:2)#是用线，还是用d点

setwd('G:\\WuXiang\\Experiment\\GSEA')
pdf(file="./INTERFERON_GAMMA.pdf",width= 3.25,height= 2.5)
p1
dev.off()

############################################### Viral mimicry
gene <- c('TP53','SETDB1‌','DNMT1‌','KDM1A','IFIH1','DDX58','DNMT1','IRF3','IRF7','OAS1','OAS2','OAS3','MX1','ISG15‌')
term <- rep('Viral mimicry', length(gene))
Gmt <- cbind(term,gene) %>% as.data.frame()
‌
egmt_viral <- GSEA(geneList, TERM2GENE = Gmt, 
             minGSSize = 10,maxGSSize = 10000,
             pvalueCutoff = 10,pAdjustMethod="fdr", 
                    seed=TRUE, by="fgsea")

p1 <- gseaplot2(egmt_viral,#数据
            "Viral mimicry",#画哪一列的信号通路
            #title = "Prion disease",#标题
            base_size = 15,#字体大小
            color = "green",#线条的颜色
            #pvalue_table = TRUE,#加不加p值
            ES_geom="line",subplots=1:2)#是用线，还是用d点

setwd('G:\\WuXiang\\Experiment\\GSEA')
pdf(file="./Viral mimicry.pdf",width= 3.25,height= 2.5)
p1
dev.off()

############################################### PRC2
gene <- c('EZH1','EZH2','SUZ12','EED','RBBP4','RBBP7','JARID2','AEBP2','PHF1','MTF2','PHF19','EPOP','DNMT1','UHRF1')
term <- rep('PRC2', length(gene))
Gmt <- cbind(term,gene) %>% as.data.frame()
‌
egmt_viral <- GSEA(geneList, TERM2GENE = Gmt, 
             minGSSize = 9,maxGSSize = 10000,
             pvalueCutoff = 10,pAdjustMethod="fdr", 
                    seed=TRUE, by="fgsea")

p1 <- gseaplot2(egmt_viral,#数据
            "PRC2",#画哪一列的信号通路
            #title = "Prion disease",#标题
            base_size = 15,#字体大小
            color = "green",#线条的颜色
            #pvalue_table = TRUE,#加不加p值
            ES_geom="line",subplots=1:2)#是用线，还是用d点

setwd('G:\\WuXiang\\Experiment\\GSEA')
pdf(file="./PRC2.pdf",width= 3.25,height= 2.5)
p1
dev.off()
library(enrichplot)
