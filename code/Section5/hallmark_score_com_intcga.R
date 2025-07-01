
library('msigdbr')
library("GSVA")

# load data
tcga.lihc.cli.data <- readRDS(file='/data/tcga_risk_score.rds')
load(file='/data/tcga_lihc_rnaseq.RData')

tcga.lihc.cli.data$risk.categ <- factor(tcga.lihc.cli.data$risk.categ, levels = c('low risk', 'high risk'))

tcga.lihc.cli.data$patient_id <- paste0(tcga.lihc.cli.data$patient_id, '-01')
tcga.lihc.exp <- tcga.lihc.tpm[, tcga.lihc.cli.data$patient_id]
tcga.lihc.exp <- log2(tcga.lihc.exp + 1)
rownames(tcga.lihc.exp) <- stringr::str_split_i(rownames(tcga.lihc.exp), pattern = '\\|', i = 2)


# gene set variation analysis
h_gene_sets = msigdbr(species = "human", category = "H")
h_gene_sets = split(x = h_gene_sets$entrez_gene, f = h_gene_sets$gs_name)


## build GSVA parameter object
ssgseapar <- ssgseaParam(
  exprData = as.matrix(tcga.lihc.exp),
  geneSets = h_gene_sets,
  minSize = 1,
  maxSize = Inf,
  alpha = 0.25,
  normalize = TRUE)

## estimate GSVA enrichment scores
ssgsea_es <- gsva(ssgseapar)

# saveRDS(ssgsea_es, file = '/data/tcga_lihc_hallmark_score.rds')


#####################compare high risk vs. low risk 
# load data
tcga.lihc.cli.data <- readRDS(file='/data/tcga_risk_score.rds')
tcga.lihc.cli.data$risk.categ <- factor(tcga.lihc.cli.data$risk.categ, levels = c('low risk', 'high risk'))
tcga.lihc.cli.data$patient_id <- paste0(tcga.lihc.cli.data$patient_id, '-01')

ssgsea_es <- readRDS(file = '/data/tcga_lihc_hallmark_score.rds')

tcga.lihc.cli.data <- cbind(tcga.lihc.cli.data, t(ssgsea_es[, tcga.lihc.cli.data$patient_id]))

p.values <- sapply(rownames(ssgsea_es), function(hallmark){
  
  gseascore <- tcga.lihc.cli.data[, c(hallmark, 'risk.categ')]
  colnames(gseascore) <- c('hallmark', 'risk.categ')
  
  meanScore <- gseascore %>% dplyr::select(hallmark, risk.categ) %>% group_by(risk.categ) %>% 
    summarise(meanscore = mean(hallmark))
  
  
  p.value <- wilcox.test(hallmark ~ risk.categ, alternative = "two.sided", data = gseascore)$p.value
  return(setNames(c(meanScore$meanscore, p.value), c('Lsocre', 'Hscore', 'pvalue')))
}, simplify = F)

p.values <- do.call(rbind, p.values) %>% as.data.frame() %>% 
  mutate(qvalue = p.adjust(pvalue, method = 'fdr')) %>% arrange(qvalue) %>% mutate(log2FC = log2(Lsocre/Hscore))


# write.table(p.values, file = '/result/Section5/hallmark_score_com_intcga.txt',
#             row.names = T, col.names = T, sep ='\t', quote = F)
# 

##################### plot

ExpClusteringHeatmap <- function(exp.data, gene.sig, sample.anno, file.name){
  library(ComplexHeatmap)
  library(dplyr)
  library(tibble)
  library(circlize)
  library(tidyr)
  library(viridis)
  library(RColorBrewer)
  
  exp.data <- exp.data[gene.sig, ] %>% t() %>% as.data.frame()
  exp.data <- exp.data[rownames(sample.anno), ]
  
  # Annotation
  row.col = structure(names = c('high risk', 'low risk'), c("black", "blue"))
  row.ha <- HeatmapAnnotation(df=sample.anno, name = 'risk.categ', which='row', col=list(risk.categ=row.col))
  
  
  # ht1: expression heatmap
  # make heatmap
  ht1 = Heatmap(
    matrix=as.matrix(exp.data), 
    name = "ht1", 
    col = colorRamp2(seq(-0.5, 1.0, by = 0.01), viridis(151, option = "magma")),
    column_dend_height = unit(15, "mm"),
    row_dend_width = unit(15, "mm"),
    show_row_names = FALSE,
    width = 1,
    heatmap_legend_param = list(title = NULL, color_bar = "continuous"),
    column_names_gp = gpar(fontsize = 6), 
    cluster_columns = T,
    cluster_rows = F,
    left_annotation = row.ha
  )
  
  
  # plot ComplexHeatmaps 
  pdf(paste0(file.name,'.pdf'))
  print(ht1)
  
  # ht1: add border
  decorate_heatmap_body("ht1", {grid.rect(gp = gpar(fill = "transparent", col = "black", lwd = 0.75))})
  
  dev.off()
}


setwd('/result/Section5/')
sam.anno <- tcga.lihc.cli.data[, 'risk.categ', F] %>% arrange(risk.categ)
hallmark <- rownames(p.values)[p.values$qvalue < 0.05]

ExpClusteringHeatmap(ssgsea_es, hallmark, sam.anno, 'tcga_hallmark_score')


hallmarks <- c('HALLMARK_G2M_CHECKPOINT', 'HALLMARK_E2F_TARGETS', 'HALLMARK_MITOTIC_SPINDLE',
'HALLMARK_FATTY_ACID_METABOLISM', 'HALLMARK_BILE_ACID_METABOLISM', 'HALLMARK_GLYCOLYSIS')



plots <- lapply(hallmarks, function(hallmark){
  
  p <- ggboxplot(tcga.lihc.cli.data, x = "risk.categ", y = hallmark, 
            color = "risk.categ", palette = "jco", add.params = list(size = 0.6), 
            add = "jitter", xlab = F, ylab = hallmark) +
    stat_compare_means(aes(label = ..p.format..), method = "wilcox.test", label.x = 1.5) + 
    scale_x_discrete(labels = c('high risk' = "High risk", 'low risk' = "Low risk")) + theme(legend.position = "none")
  
  return(p)
})
 

ggsave(cowplot::plot_grid(plotlist = plots, ncol=4, nrow=4), 
       file='/result/Section5/select_hallmark_score_com.pdf')


