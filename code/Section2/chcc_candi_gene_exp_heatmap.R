
# load data
# setwd('/data/OriginalData/Gao_Cell_2019')
# identifier <- readxl::read_xlsx(path='About the RNA and protein Identifier match.xlsx', sheet = 1, col_names = TRUE, skip=1)
# colnames(identifier) <- c('protein_T', 'protein_N', 'rna_T', 'rna_N')
# 
# 
# # expression data
# exp.data <- read.csv(file='HCC_UQ_RSEM.tsv', header=TRUE, sep='\t', stringsAsFactors=FALSE)
# exp.data <- tibble::column_to_rownames(exp.data, var = "protein")
# exp.data <- dplyr::mutate(exp.data, X=NULL)
# colnames(exp.data) <- gsub('X', '', colnames(exp.data))
# 
# exp.data <- exp.data[, as.character(c(identifier$rna_T, identifier$rna_N))]
# colnames(exp.data) <- c(paste0('T', identifier$rna_T), paste0('N', identifier$rna_N))
# 
# 
# exp.data <- round(exp.data)
# exp.data <- exp.data[rowSums(exp.data >= 1) >= 64, ]
# 
# 
# 
# # sample information
# sams.info <- data.frame(labels=c(rep('Tumor', nrow(identifier)), rep('Normal', nrow(identifier))))
# rownames(sams.info) <- c(paste0('T', identifier$rna_T), paste0('N', identifier$rna_N))
# 
# 
# # Differential expression analysis
# library(DESeq2)
# dds <- DESeqDataSetFromMatrix(countData = exp.data, colData = sams.info, design = ~ labels)
# 
# vsd <- vst(dds, blind = FALSE)
# gao.etal.vst.exp <- assay(vsd)
# 
# dds <- DESeq(dds)
# gao.etal.diff.res <- results(dds)
# 
# save(sams.info, gao.etal.vst.exp, gao.etal.diff.res, file='/data/Gao_Cell_2019_diff_exp.RData')


#####################diff gene
load(file='/data/Gao_Cell_2019_diff_exp.RData')
load(file='/data/exp_gene_anno.RData') 
candi.prog.gene <- readRDS(file='/data/candi_prog_sig_gene.rds')


rownames(gao.etal.vst.exp) <- stringr::str_split_i(rownames(gao.etal.vst.exp), pattern = '\\.', i = 1)
rownames(gao.etal.diff.res) <- stringr::str_split_i(rownames(gao.etal.diff.res), pattern = '\\.', i = 1)

gao.etal.vst.exp <- gao.etal.vst.exp[intersect(rownames(gao.etal.vst.exp), exp.gene.anno$EnsemblID), ]
gao.etal.diff.res <- gao.etal.diff.res[intersect(rownames(gao.etal.diff.res), exp.gene.anno$EnsemblID), ]

rownames(gao.etal.vst.exp) <- exp.gene.anno$Symbol[match(rownames(gao.etal.vst.exp), exp.gene.anno$EnsemblID)]
rownames(gao.etal.diff.res) <- exp.gene.anno$Symbol[match(rownames(gao.etal.diff.res), exp.gene.anno$EnsemblID)]


gao.etal.diff.res[intersect(candi.prog.gene, rownames(gao.etal.diff.res)), ] %>% as.data.frame()


# heatmap
setwd('/result/Section2')
colnames(sams.info) <- 'Type'
ExpClusteringHeatmap(gao.etal.vst.exp, candi.prog.gene, sams.info, 'chcc_candi_gene_exp_heatmap')

