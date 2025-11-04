####################################
library(dplyr)
library(Seurat)
library(ggplot2)
library(ggrepel)  # 用于防止标签重叠
library(data.table)
library(enrichR)
library(igraph)
library(ggbreak)
library(patchwork)
##########################################
##########################################
setwd('G:\\CKDwork\\knokdown\\Experiment\\ACSM2A_Mature')
load('ACSM2A_result_G10000_C2000.Rdata')
################################################################## enrich function
enrichFunction <- function(X, fdrThreshold = fdrThreshold){
    E <- enrichr(X, c('KEGG_2019_Human', 'GO_Biological_Process_2018','GO_Cellular_Component_2018', 'GO_Molecular_Function_2018','BioPlanet_2019', 'WikiPathways_2019_Human', 'Reactome_2016'))
    E <- do.call(rbind.data.frame, E)
    E <- E[E$Adjusted.P.value < fdrThreshold,]
    E <- E[order(E$Adjusted.P.value),]
    E$Term <- unlist(lapply(strsplit(E$Term,''), function(X){
        X[1] <- toupper(X[1])
        X <- paste0(X,collapse = '')
        X <- gsub('\\([[:print:]]+\\)|Homo[[:print:]]+|WP[[:digit:]]+','',X)
        X <- gsub("'s",'',X)
        X <- unlist(strsplit(X,','))[1]
        X <- gsub('[[:blank:]]$','',X)
        return(X)
    }))
    selectedSet <- rep(FALSE, nrow(E))
    for(i in seq_len(nrow(E))){
        if(i == 1){
            selectedSet[i] <- TRUE
        } else {
            A <- unique(unlist(strsplit(E[which(selectedSet[seq_len(i)]),'Genes'], ';')))
            B <- unlist(strsplit(E[i,'Genes'], ';'))
            selectedSet[i] <- !all(B %in% A)
        }
    }
    gSets <- table(toupper(E$Term))
    gSets <- names(gSets[gSets > 1])
    for(i in gSets){
        selectedSet[which(toupper(E$Term) %in% i)[-1]] <- FALSE
    }
    E <- E[selectedSet,]
    if(nrow(E) > nCategories){
        E <- E[seq_len(nCategories),]  
    }
    return(E)
}

#################################################
gKO <- 'ACSM2A'
q = 0.99
annotate = TRUE
nCategories = 100
fdrThreshold = 0.05
X <- result
################################### plotKO modify
gList <- unique(c(gKO, X$diffRegulation$gene[X$diffRegulation$distance > 1e-10 & X$diffRegulation$p.adj < 0.05]))
sCluster <- as.matrix(X$tensorNetworks$WT[gList,gList])
koInfo <- sCluster[gKO,]
gList <- gList[!grepl('^mt-|^Rpl|^Rps',gList, ignore.case = TRUE)]
sCluster[abs(sCluster) <= quantile(abs(sCluster), q)] <- 0
sCluster[gKO,] <- koInfo
diag(sCluster) <- 0
sCluster <-  reshape2::melt(as.matrix(sCluster))
colnames(sCluster) <- c('from', 'to', 'W')
sCluster <- sCluster[sCluster$W != 0,]
netPlot <- graph_from_data_frame(sCluster, directed = TRUE)
dPlot <- centr_degree(netPlot)$res
W <- rep(1,nrow(sCluster))
sG   <- (names(V(netPlot))[dPlot > 1])[-1]
W[sCluster$from %in% sG] <- 0.2
W[sCluster$to %in% sG] <- 0.2
W[sCluster$from %in% gKO] <- 1
W[sCluster$from %in% gKO & sCluster$to %in% sG] <- 0.8
set.seed(1)
layPlot <- layout_with_fr(netPlot, weights = W)
dPlot <- (dPlot/max(dPlot))*20
######################################################
gKO <- 'ACSM2A'
nCategories = 100
gList <- unique(c(gKO, result$diffRegulation$gene[result$diffRegulation$distance > 1e-10 & result$diffRegulation$p.adj < 0.05]))
E <- enrichFunction(gList, 1)
aaa <- c(1:4,6,7,9:11,14,15,19,21,29)
E <- E[aaa,]
#######################################################
tPlot <- strsplit(E$Genes, ';')
pPlot <- matrix(0,nrow = length(V(netPlot)), ncol = nrow(E))
rownames(pPlot) <- toupper(names(V(netPlot)))
for(i in seq_along(tPlot)){
    pPlot[unlist(tPlot[i]),i] <- 1
}
pPlot <- lapply(seq_len(nrow(pPlot)), function(X){as.vector(pPlot[X,])})
names(pPlot) <- names(V(netPlot))
tPlot <- unique(unlist(tPlot))
eGenes <- toupper(names(V(netPlot))) %in% tPlot
vColor <- rgb(195/255, 199/255, 198/255 ,0.3)
###########################################################
pieColors <- list(hcl.colors(nrow(E), palette = 'Zissou 1', alpha = 0.7))
################################################## Plot 
pdf(file = 'Gene_network.pdf', width = 20, height = 10)
par(mar=c(4,0,0,0), xpd = TRUE)
suppressWarnings(plot(netPlot,
                      layout = layPlot, 
                      edge.arrow.size=.2,
                      vertex.label.color="black", 
                      vertex.shape = ifelse(eGenes,'pie','circle'),
                      vertex.pie = pPlot,
                      vertex.size = 10+dPlot, 
                      vertex.pie.color=pieColors,
                      vertex.label.family="Times", 
                      vertex.label.font=ifelse(eGenes,2,1),
                      edge.color = ifelse(E(netPlot)$W > 0, 'red', 'blue'),
                      edge.curved = ifelse(W == 0.2, 0, 0.1),
                      vertex.color = vColor, 
                      vertex.frame.color = NA))

sigLevel <- formatC(E$Adjusted.P.value, digits = 2, format = 'g', width = 0, drop0trailing = TRUE)
gSetNames <- lengths(strsplit(E$Genes, ';'))
gSetNames <- paste0('(', gSetNames,') ', E$Term, ' FDR = ', sigLevel)
legend(x = -1.05, y = -1.05, legend = gSetNames, bty = 'n', ncol = 2, cex = 1, col = unlist(pieColors), pch = 16)
dev.off()

###################################################################
library(gg.gap)
top_genes <- head(result$diffRegulation[order(-result$diffRegulation$FC), ], 21)
top_genes <- top_genes[-1,]

p1 <- ggplot(top_genes, aes(x=reorder(gene, log2(FC)), y=log2(FC))) +
  geom_bar(stat='identity', fill='steelblue') + scale_y_continuous(expand = expansion(c(0, 0))) +
  coord_flip() +
  labs(title="Top 20 Differentially Regulated Genes",
       x="", y="Log2(Fold change)") + 
  theme_minimal()

pdf(file = 'Top20_diff_gene.pdf', width = 5.1, height = 3.5)
p1
dev.off()


########################################################
df <- result$diffRegulation
df <- df[-1,]
df <- df[df$Z > 0,]
df$color <- ifelse(df$FC > 2 & df$p.adj < 0.01, 'red','black')
df$log_pval <- -log10(df$p.adj)
label_genes <- subset(df, abs(Z) > 2 & p.adj < 0.01)
ggplot(df, aes(x=Z, y=log_pval)) +
  geom_point(alpha=0.5,color = df$color) +
  geom_hline(yintercept=-log10(0.05), linetype="dashed", color="red") +
  geom_vline(xintercept=c(2), linetype="dashed", color="blue") +
  geom_text_repel(data=label_genes, aes(label=gene),
                  size=3, max.overlaps=50) +
  labs(title="",
       x="Z-score", y="-log10(p-value)") + 
  theme_classic() 

pdf(file = 'Vocano_gene_plot.pdf', width = 4.15, height = 4.1)
ggplot(df, aes(x=Z, y=log_pval)) +
  geom_point(alpha=0.5,color = df$color) +
  geom_hline(yintercept=-log10(0.05), linetype="dashed", color="red") +
  geom_vline(xintercept=c(2), linetype="dashed", color="blue") +
  geom_text_repel(data=label_genes, aes(label=gene),
                  size=3, max.overlaps=50) +
  labs(title="",
       x="Z-score", y="-log10(p-value)") + 
  theme_classic() 
dev.off()

