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
library(enrichR)
##########################################
setwd('G:\\CKDwork\\knokdown\\Experiment\\ACSM2A_Mature')
load('ACSM2A_result_G10000_C2000.Rdata')
abc <- result$diffRegulation
abc <- abc[abc$distance > 1e-10 & abc$p.adj < 0.05,]
abc <- abc[,c(1,4,6)]
abc$FC <- log2(abc$FC)
abc <- abc[-1,]
#####################################
h.gmt <- getGmt('./gene sets for Proximal tubule.gmt')
##################################### GSEA
IR <- abc
#########################
geneList <- IR$FC
names(geneList) <- IR$gene 
#按照logFC的值来排序geneList

Plot_GeneSet <- cbind(h.gmt[[1]]@setName,h.gmt[[1]]@geneIds) %>% as.data.frame()
colnames(Plot_GeneSet) <- c('term','gene')

set.seed(1234)
egmt <- GSEA(geneList, TERM2GENE = Plot_GeneSet, 
             minGSSize = 1,maxGSSize = 10000,
             pvalueCutoff = 1,pAdjustMethod="fdr", 
                    seed=TRUE, by="fgsea")

gesa_res <- as.data.frame(egmt@result)
###############################################
#gseaplot2(egmt, "HALLMARK_INFLAMMATORY_RESPONSE", color = "firebrick",
#            rel_heights=c(1, .2, .6), subplots=1:2)
############################################### Inflammatory Res
p1 <- gseaplot2(egmt,#数据
            "Proximal tubule bicarbonate reclamation",#画哪一列的信号通路
            #title = "Prion disease",#标题
            base_size = 15,#字体大小
            color = "green",#线条的颜色
            #pvalue_table = TRUE,#加不加p值
            ES_geom="line",subplots=1:2)#是用线，还是用d点

setwd('G:\\WuXiang\\Experiment\\GSEA')
pdf(file="./InflammatoryRes.pdf",width= 3.25,height= 2.5)
p1
dev.off()

