library(ggplot2)
library(dplyr)
library(cowplot)
library(reshape2)
library(gghalves)
library(ggpubr)

# read data
immuneHCC <- read.csv2("/data/TCGA.Kallisto.fullIDs.cibersort.relative.tsv", header = TRUE, sep = '\t') 
immuneHCC <- type.convert(immuneHCC, as.is = TRUE)
immuneHCC <- immuneHCC %>% subset(CancerType == 'LIHC' & str_sub(SampleID, 14, 15) == '01')
rownames(immuneHCC) <- gsub('\\.', '-', str_sub(immuneHCC$SampleID, 1, 12))
immuneHCC <- immuneHCC %>% select(-SampleID, -CancerType)

evoSigScore <- readRDS(file = '/data/tcga_risk_score.rds')
evoSigScore <- evoSigScore %>% column_to_rownames(var = 'patient_id')

# Immune Cellular Fraction Estimates
immune22Score <- merge(evoSigScore, immuneHCC, by = 'row.names')[, 19:41]
immune22Score <- melt(data = immune22Score, id.vars = 'risk.categ', variable.name = 'cellType', value.name = 'fraction')


imm22p.values <- sapply(unique(immune22Score$cellType), function(ct){
  
  p.value <- wilcox.test(fraction~risk.categ, data = subset(immune22Score, cellType == ct), alternative = "two.sided")$p.value
  p.value <- setNames(p.value, ct)
  return(p.value)
})

# > imm22p.values[imm22p.values < 0.05]
# B.cells.naive   T.cells.CD4.memory.resting T.cells.CD4.memory.activated             NK.cells.resting 
# 1.263308e-02                 3.005270e-02                 2.344712e-02                 1.937435e-03 
# Monocytes               Macrophages.M0           Mast.cells.resting                  Neutrophils 
# 3.165378e-04                 6.894344e-09                 6.494022e-03                 2.184109e-08 


# plot
mytheme <- theme_bw() +
  theme(axis.text.x = element_text(angle = 15, hjust = 1, family="Times", size = 7), 
        axis.text.y = element_text(family = "Times",size = 7,face = "plain"), 
        axis.title.y = element_text(family = "Times",size = 7,face = "plain"), 
        axis.line = element_line(colour = "black",size = 1), 
        legend.position = "top",
        legend.text = element_text(face = "italic", family = "Times",size = 8),
        legend.title = element_blank())


showcellTypes <- names(imm22p.values)[imm22p.values < 0.05]
immune22Score <- subset(immune22Score, cellType %in% showcellTypes) %>% mutate(cellType = factor(cellType))

p1 <- ggplot() + 
  ###左半边 Group1
  geom_half_violin(data = subset(immune22Score, risk.categ == 'high risk'),
                   aes(x = cellType, y = fraction),colour = "white",fill = "#db6968",
                   side = "l",nudge = 0.01, scale = "width") +
  ###右半边 Group2
  geom_half_violin(data = subset(immune22Score, risk.categ == 'low risk'),
                   aes(x = cellType,y = fraction),colour = "white",fill = "#4d97cd",
                   side = "r",nudge = 0.01, scale = "width") + 
  ###添加均值点
  geom_point(data = immune22Score, aes(x = cellType,y = fraction, fill = risk.categ),
             stat = 'summary', fun = mean, size = 1.5,
             show.legend = FALSE,
             position = position_dodge(width = 0.3))+
  ###添加errorbar
  stat_summary(data = immune22Score, aes(x = cellType,y = fraction, fill = risk.categ),
               fun.min = function(x){quantile(x)[2]},
               fun.max = function(x){quantile(x)[4]},
               geom = 'errorbar', color = 'black',
               width = 0.1, size = 1,
               position = position_dodge(width = 0.3)) + 
  ylab("Cell fraction") + xlab(NULL) + mytheme



# Aggregate 
cellTypes <- list(
     T.cells.CD8 = 'T.cells.CD8', 
     T.cells.CD4 = c('T.cells.CD4.naive', 'T.cells.CD4.memory.resting', 'T.cells.CD4.memory.activated'),
     T.cells.follicular.helper = 'T.cells.follicular.helper',
     T.cells.regulatory..Tregs. = 'T.cells.regulatory..Tregs.',
     T.cells.gamma.delta = 'T.cells.gamma.delta',
     B.cells = c('B.cells.naive', 'B.cells.memory'),
     NK.cells = c('NK.cells.resting', 'NK.cells.activated'),
     Plasma.cells = 'Plasma.cells',
     Monocytes = 'Monocytes',
     Macrophage = c('Macrophages.M0', 'Macrophages.M1', 'Macrophages.M2'),
     Dendritic.cells = c('Dendritic.cells.resting', 'Dendritic.cells.activated'),
     Mast.cells = c('Mast.cells.resting', 'Mast.cells.activated'),
     Neutrophils = 'Neutrophils',
     Eosinophils = 'Eosinophils')


immune14Score <- lapply(cellTypes, function(celltype){
  
  cellFrac <- rowSums(immuneHCC[, celltype, FALSE])
  return(cellFrac)
})

immune14Score <- cbind.data.frame(immune14Score)



immune14Score <- merge(evoSigScore, immune14Score, by = 'row.names')[, 19:33]
immune14Score <- melt(data = immune14Score, id.vars = 'risk.categ', variable.name = 'cellType', value.name = 'fraction')

imm14p.values <- sapply(unique(immune14Score$cellType), function(ct){
  
  p.value <- wilcox.test(fraction~risk.categ, data = subset(immune14Score, cellType == ct), alternative = "two.sided")$p.value
  p.value <- setNames(p.value, ct)
  return(p.value)
})

# imm14p.values[imm14p.values < 0.05]
# T.cells.CD4     NK.cells    Monocytes   Macrophage   Mast.cells  Neutrophils 
# 3.067156e-02 6.117837e-03 3.165378e-04 1.550218e-04 1.610473e-02 2.184109e-08 

showcellTypes <- names(imm14p.values)[imm14p.values < 0.05]
immune14Score <- subset(immune14Score, cellType %in% showcellTypes) %>% mutate(cellType = factor(cellType))

p2 <- ggplot() + 
  ###左半边 Group1
  geom_half_violin(data = subset(immune14Score, risk.categ == 'high risk'),
                   aes(x = cellType, y = fraction),colour = "white",fill = "#db6968",
                   side = "l",nudge = 0.01, scale = "width") +
  ###右半边 Group2
  geom_half_violin(data = subset(immune14Score, risk.categ == 'low risk'),
                   aes(x = cellType,y = fraction),colour = "white",fill = "#4d97cd",
                   side = "r",nudge = 0.01, scale = "width") + 
  ###添加均值点
  geom_point(data = immune14Score, aes(x = cellType,y = fraction, fill = risk.categ),
             stat = 'summary', fun = mean, size = 1.5,
             show.legend = FALSE,
             position = position_dodge(width = 0.3))+
  ###添加errorbar
  stat_summary(data = immune14Score, aes(x = cellType,y = fraction, fill = risk.categ),
               fun.min = function(x){quantile(x)[2]},
               fun.max = function(x){quantile(x)[4]},
               geom = 'errorbar', color = 'black',
               width = 0.1, size = 1,
               position = position_dodge(width = 0.3)) + 
  ylab("Cell fraction") + xlab(NULL) + mytheme


combinedPlots <- plot_grid(p1, p2, labels = "AUTO", ncol = 2, nrow = 3)
ggsave("/result/Section6/CIBERSORT_cellfraction_tcga.pdf", combinedPlots)



# Aggregate 
Lymphocytes <- list(Lymphocytes = 
                    c('B.cells.naive', 'B.cells.memory', 'T.cells.CD4.naive', 'T.cells.CD4.memory.resting', 
                      'T.cells.CD4.memory.activated', 'T.cells.follicular.helper', 'T.cells.regulatory..Tregs.', 'T.cells.gamma.delta', 
                      'T.cells.CD8', 'NK.cells.resting', 'NK.cells.activated', 'Plasma.cells'))
                  

lymphocyteScore <- lapply(Lymphocytes, function(celltype){
  
  cellFrac <- rowSums(immuneHCC[, celltype, FALSE])
  return(cellFrac)
})

lymphocyteScore <- cbind.data.frame(lymphocyteScore)



lymphocyteScore <- merge(evoSigScore, lymphocyteScore, by = 'row.names')[, 19:20]
lymphocyteScore <- melt(data = lymphocyteScore, id.vars = 'risk.categ', variable.name = 'cellType', value.name = 'fraction')

# wilcox.test(fraction~risk.categ, data = lymphocyteScore, alternative = "two.sided")$p.value
# 0.06604536

p3 <- ggplot(lymphocyteScore, aes(x = risk.categ, y = fraction, fill = risk.categ)) +
  geom_violin(trim =TRUE) +
  geom_boxplot(width = 0.1, fill = "white", alpha = 0.7) + labs(x = NULL, y = 'Leukocyte fraction') + 
  geom_signif(comparisons = list(c("high risk", "low risk")),
              test = "wilcox.test", textsize = 4, step_increase = 0.1) + 
  scale_fill_brewer(palette="Blues") + 
  theme(
    panel.background = element_rect(fill = "white"),  # 白色背景
    panel.grid.major = element_blank(),               # 移除主要网格线
    panel.grid.minor = element_blank(),               # 移除次要网格线
    axis.line = element_line(color = "black"),        # 黑色坐标轴线
    plot.background = element_rect(fill = "white"),   # 整个绘图区域白色背景
    panel.border = element_blank(), # 移除面板边框
    legend.position = "top")



leukocyteHCC <- read.csv2("/data/TCGA_all_leuk_estimate.masked.20170107.tsv", header = FALSE, sep = '\t') 
colnames(leukocyteHCC) <- c('CancerType', 'sampleID', 'leukocyte')
leukocyteHCC <- type.convert(leukocyteHCC, as.is = TRUE)
leukocyteHCC <- leukocyteHCC %>% subset.data.frame(CancerType == 'LIHC' & str_sub(sampleID, 14, 15) == '01')
rownames(leukocyteHCC) <- str_sub(leukocyteHCC$sampleID, 1, 12)
leukocyteHCC <- leukocyteHCC %>% select(-CancerType, -sampleID)

leukocyteScore <- merge(evoSigScore, leukocyteHCC, by = 'row.names')[, 19:20]
leukocyteScore <- melt(data = leukocyteScore, id.vars = 'risk.categ', variable.name = 'cellType', value.name = 'fraction')

# wilcox.test(fraction~risk.categ, data = leukocyteScore, alternative = "two.sided")$p.value
# 0.37306

p4 <- ggplot(leukocyteScore, aes(x = risk.categ, y = fraction, fill = risk.categ)) +
  geom_violin(trim =TRUE) +
  geom_boxplot(width = 0.1, fill = "white", alpha = 0.7) + labs(x = NULL, y = 'Leukocyte fraction') + 
  geom_signif(comparisons = list(c("high risk", "low risk")),
              test = "wilcox.test", textsize = 4, step_increase = 0.1) + 
  scale_fill_brewer(palette="Blues") + 
  theme(
    panel.background = element_rect(fill = "white"),  # 白色背景
    panel.grid.major = element_blank(),               # 移除主要网格线
    panel.grid.minor = element_blank(),               # 移除次要网格线
    axis.line = element_line(color = "black"),        # 黑色坐标轴线
    plot.background = element_rect(fill = "white"),   # 整个绘图区域白色背景
    panel.border = element_blank(), # 移除面板边框
    legend.position = "top")



combinedPlots <- plot_grid(p3, p4, labels = "AUTO", ncol = 2, nrow = 3)
ggsave("/result/Section6/CIBERSORT_leukocytefraction.pdf", combinedPlots)




