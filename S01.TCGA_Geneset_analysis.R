############## Pathway
library(SummarizedExperiment)
library(TCGAbiolinks)
library(data.table)
library(dplyr)
##############
setwd('G:\\CKDwork\\Other_analysis\\GeneSet')
pbmc <- read.csv('./QTL_analysis.csv')

ALL <- pbmc$id

ansEA <- TCGAanalyze_EAcomplete(
    TFname = "CKD-related Genes",
    RegulonList = ALL
)  

####### Plot the pathway
TCGAvisualize_EAbarplot(
    tf = rownames(ansEA$ResBP),
    GOBPTab = ansEA$ResBP,
    PathTab = ansEA$ResPat,
    nRGTab = ALL,
    nBar = 10,
    text.size = 2,
    fig.width = 20,
    fig.height = 15
)
