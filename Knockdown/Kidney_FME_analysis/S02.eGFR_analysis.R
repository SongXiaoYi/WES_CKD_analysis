##### Prepare GTEx and TCGA row count
###############################
library(dplyr)
library(ComplexHeatmap)
library(data.table)
#######################
library(dplyr)
library(Hmisc)
library(tidyr)
library(tibble)
#######################
library(GSVA)
library(GSEABase)
library(limma)
library(ggpubr)
library(reshape2)
##################################
setwd('G:\\CKDwork\\knokdown\\Experiment\\ACSM2A_Mature')
Score <- fread('./Neptune_total_FME_score.csv') %>% as.data.frame()
Score <- Score[Score$EGFR != 'None',]
Score <- na.omit(Score)
Score$EGFR <- as.numeric(Score$EGFR)
Score$Group <- Score$EGFR > 45
#Score$ACSM2A <- log(Score$ACSM2A)

################################# Color assign
Type_colors <- c('#D9D9D9', '#CA93B5')   ####PT VS RT
ROI_colors <- c('#80C1C0', '#D695FE','#F0BEA6')      ######## ROI
                 
################################################ Plot
################################ Total compare
newdata <- Score
newdata$Group <- factor(newdata$Group, level = c('TRUE','FALSE'))
my_comparisons <- list(c("TRUE", "FALSE"))
a = ggboxplot(newdata, x="Group", y="ACSM2A", color = "Group", palette = Type_colors, add = "jitter", shape = 16,
            ylab="ACSM2A normalizaed experission", xlab="",outlier.shape = NA,bxp.errorbar = TRUE)
a = a+rotate_x_text(51) + theme(plot.title = element_text(size=12,hjust=0.5))
a + stat_compare_means() 

setwd('G:\\CKDwork\\knokdown\\Experiment\\ACSM2A_Mature')
pdf(file = 'eGFR_analysis.pdf', width = 3.1, height = 4)
a + stat_compare_means()
dev.off()

