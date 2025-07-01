
coInM1 <- c('PDCD1', 'PDCD1', 'HAVCR2', 'LAIR1', 'LAIR1', 'CTLA4', 'CTLA4', 
'TIGIT', 'TIGIT', 'TIGIT', 'TIGIT', 'TIGIT', 'CD160', 'BTLA')

coInM2 <- c('CD274', 'PDCD1LG2', 'LGALS9', 'PTPN6', 'PTPN11', 'CD80',
'CD86', 'PVR', 'NECTIN1', 'NECTIN2', 'NECTIN3', 'NECTIN4', 'TNFRSF14', 'TNFRSF14')


coStiM1 <- c('ICOS', 'TNFRSF4', 'TNFRSF8', 'TNFRSF9', 'TNFRSF18',
'CD40LG', 'CD27', 'CD28', 'CD28', 'CD226', 'CD226', 'CD226', 'CD226', 'CD226', 'TNFSF14')


coStiM2 <- c('ICOSLG', 'TNFSF4', 'TNFSF8', 'TNFSF9', 'TNFSF18',
'CD40', 'CD70', 'CD80', 'CD86', 'PVR', 'NECTIN1', 'NECTIN2', 'NECTIN3', 'NECTIN4', 'TNFRSF14')



# PDCD1, CD274, CTLA4, HAVCR2, LAG3, TIGIT, BTLA, SIRPA, SIGLEC7, CD276, TNFRSF14
# VSIR, TIGIT, PDCD1, LAG3, KLRD1, IAPP, HAVCR2, CTLA4, CD96, CD274, CD226
# IDO1, CD274, HAVCR2, PDCD1, CTLA4, LAG3, CD8A, CXCL10, CXCL9, GZMA, GZMB, PRF1, IFNG, TBX2, TNF, CD80, CD86

tcga.lihc.cli.data <- readRDS(file='/data/tcga_risk_score.rds')
tcga.lihc.cli.data <- tcga.lihc.cli.data %>% select(patient_id, risk.categ) %>% column_to_rownames(var = 'patient_id')


tcga.lihc.vst.exp <- readRDS(file='/data/tcga_vst_norm_geneExp.rds')
load(file='/data/exp_gene_anno.RData') 
rownames(tcga.lihc.vst.exp) <- exp.gene.anno$Symbol[match(rownames(tcga.lihc.vst.exp), exp.gene.anno$GeneID)]
colnames(tcga.lihc.vst.exp) <- str_sub(colnames(tcga.lihc.vst.exp), 1, 12)


ExpClusteringHeatmap <- function(exp.data, gene.sig, sample.anno, file.name){
  library(ComplexHeatmap)
  library(dplyr)
  library(tibble)
  library(circlize)
  library(tidyr)
  library(viridis)
  library(RColorBrewer)
  
  
  # mat1
  mat1 <- exp.data[gene.sig, ] %>% t() %>% as.data.frame() 
  
  # Z-score
  tmp <- sapply(X=mat1, FUN=scale)
  tmp <- as.data.frame(tmp)
  rownames(tmp) <- rownames(mat1)
  mat1 <- tmp
  
  
  # Annotation
  sample.anno <- sample.anno[rownames(mat1), , FALSE]
  row.col = structure(names = c('high risk', 'low risk'), c("black", "blue"))
  row.ha <- HeatmapAnnotation(df=sample.anno, name = 'risk.categ', which='row', col=list(risk.categ=row.col))
  
  
  # ht1: expression heatmap
  # make heatmap
  ht1 = Heatmap(
    matrix=as.matrix(mat1), 
    name = "ht1", 
    col = colorRamp2(seq(-3,3), viridis(7, option = "magma")),
    column_dend_height = unit(15, "mm"),
    row_dend_width = unit(15, "mm"),
    show_row_names = FALSE,
    width = 1,
    heatmap_legend_param = list(title = NULL, color_bar = "continuous"),
    column_names_gp = gpar(fontsize = 6), 
    cluster_columns = FALSE,
    cluster_rows = FALSE,
    left_annotation = row.ha
  )
  
  
  # plot ComplexHeatmaps 
  pdf(paste0(file.name,'.pdf'))
  print(ht1)
  
  # ht1: add border
  decorate_heatmap_body("ht1", {grid.rect(gp = gpar(fill = "transparent", col = "black", lwd = 0.75))})
  
  dev.off()
}


imScores <- lapply(list(coInM1, coInM2, coStiM1, coStiM2), function(geneSet){
  
  score <- colMeans(tcga.lihc.vst.exp[geneSet,])
  return(score)
})

imScores <- cbind.data.frame(imScores)
colnames(imScores) <- c('coInM1', 'coInM2', 'coStiM1', 'coStiM2')
imScores <- merge(tcga.lihc.cli.data, imScores, by = 'row.names') %>% column_to_rownames(var = 'Row.names')


setwd('/result/Section6')

ExpClusteringHeatmap(tcga.lihc.vst.exp[, imScores %>% select(risk.categ, coInM1) %>% arrange(desc(risk.categ), coInM1) %>% rownames()], 
                     coInM1, tcga.lihc.cli.data, 'tcga_coInM1_exp_heatmap')

ExpClusteringHeatmap(tcga.lihc.vst.exp[, imScores %>% select(risk.categ, coInM2) %>% arrange(desc(risk.categ), coInM2) %>% rownames()], 
                     coInM2, tcga.lihc.cli.data, 'tcga_coInM2_exp_heatmap')



ExpClusteringHeatmap(tcga.lihc.vst.exp[, imScores %>% select(risk.categ, coStiM1) %>% arrange(desc(risk.categ), coStiM1) %>% rownames()], 
                     coStiM1, tcga.lihc.cli.data, 'tcga_coStiM1_exp_heatmap')

ExpClusteringHeatmap(tcga.lihc.vst.exp[, imScores %>% select(risk.categ, coStiM2) %>% arrange(desc(risk.categ), coStiM2) %>% rownames()], 
                     coStiM2, tcga.lihc.cli.data, 'tcga_coStiM2_exp_heatmap')


###########
others <- c('LAG3', 'LAG3', 'SIRPA', 'SIGLEC7', 'CD276', 'VSIR', 'KLRD1', 'IAPP', 'CD96', 
            'IDO1', 'CD8A', 'CXCL10', 'CXCL9', 'GZMA', 'GZMB', 'PRF1', 'IFNG', 'TBX2', 'TNF')

imExp <- tcga.lihc.vst.exp[unique(c(coInM1, coInM2, coStiM1, coStiM2)), ] %>% t() %>% as.data.frame() %>% 
  merge(tcga.lihc.cli.data, by = 'row.names') %>% column_to_rownames(var = 'Row.names') %>% 
  melt(id.vars = 'risk.categ', variable.name = 'name', value.name = 'expr')


diffPvalues <- sapply(unique(imExp$name), function(gene){
  
  p.value <- wilcox.test(expr~risk.categ, data = subset(imExp, name == gene), alternative = "two.sided")$p.value
  p.value <- setNames(p.value, gene)
  return(p.value)
})

# diffPvalues <- p.adjust(diffPvalues, method = 'fdr')
# diffPvalues[diffPvalues < 0.05]
# HAVCR2        LAIR1        CTLA4     PDCD1LG2       LGALS9         CD80         CD86      NECTIN1      NECTIN4      TNFRSF4      TNFRSF9 
# 4.189157e-04 8.172453e-03 3.068738e-03 3.953541e-02 1.806498e-03 1.088787e-04 1.430822e-02 6.238815e-06 8.488892e-03 1.218883e-02 2.390105e-02 
# TNFRSF18        CD226       ICOSLG       TNFSF4       TNFSF9 
# 5.170723e-03 3.845761e-02 2.916431e-06 1.153168e-03 5.718508e-03 


imExp <- subset(imExp, name %in% names(diffPvalues)[diffPvalues < 0.05])

mytheme <- theme_bw() +
  theme(axis.text.x = element_text(angle = 15, hjust = 1, family="Times", size = 7), 
        axis.text.y = element_text(family = "Times",size = 7,face = "plain"), 
        axis.title.y = element_text(family = "Times",size = 7,face = "plain"), 
        axis.line = element_line(colour = "black",size = 1), 
        legend.position = "top",
        legend.text = element_text(face = "italic", family = "Times",size = 8),
        legend.title = element_blank())

p1 <- ggplot() + 
  ###左半边 Group1
  geom_half_violin(data = subset(imExp, risk.categ == 'high risk'),
                   aes(x = name, y = expr),colour = "white",fill = "#db6968",
                   side = "l",nudge = 0.01, scale = "width") +
  ###右半边 Group2
  geom_half_violin(data = subset(imExp, risk.categ == 'low risk'),
                   aes(x = name,y = expr),colour = "white",fill = "#4d97cd",
                   side = "r",nudge = 0.01, scale = "width") + 
  ###添加均值点
  geom_point(data = imExp, aes(x = name,y = expr, fill = risk.categ),
             stat = 'summary', fun = mean, size = 1.5,
             show.legend = FALSE,
             position = position_dodge(width = 0.3))+
  ###添加errorbar
  stat_summary(data = imExp, aes(x = name,y = expr, fill = risk.categ),
               fun.min = function(x){quantile(x)[2]},
               fun.max = function(x){quantile(x)[4]},
               geom = 'errorbar', color = 'black',
               width = 0.1, size = 1,
               position = position_dodge(width = 0.3)) + 
  ylab("Expression") + xlab(NULL) + mytheme


combinedPlots <- plot_grid(p1, labels = "AUTO", ncol = 1, nrow = 3)
ggsave("/result/Section6/tcga_immuneCheckPoint_expr.pdf", combinedPlots)



