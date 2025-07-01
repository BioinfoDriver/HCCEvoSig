
# setwd('/data/OriginalData/GSE144269')
# cli.data <- read.csv(file='patient_sample_metadata.txt', header=TRUE, sep='\t', stringsAsFactors=FALSE)
# cli.data <- subset(cli.data, !is.na(RNASeq_T))
# 
# 
# exp.count <- read.csv(file='GSE144269_RSEM_GeneCounts.txt', header=TRUE, sep='\t', stringsAsFactors=FALSE)
# 
# exp.count$entrez_id <- do.call(rbind, strsplit(exp.count$entrez_id, split='\\|'))[, 1]
# exp.count <- tibble::column_to_rownames(exp.count, var = "entrez_id")
# 
# exp.count <- exp.count[rowSums(exp.count >= 1) >= 28, ] 
# 
# 
# # sample information
# sams.info <- data.frame(labels=c(rep('Tumor', nrow(cli.data)), rep('Normal', nrow(cli.data))))
# rownames(sams.info) <- c(cli.data$RNASeq_T, cli.data$RNASeq_NT)
# 
# exp.count <- exp.count[, rownames(sams.info)]
# 
# # Differential expression analysis
# library(DESeq2)
# dds <- DESeqDataSetFromMatrix(countData = exp.count, colData = sams.info, design = ~ labels)
# 
# vsd <- vst(dds, blind = FALSE)
# GSE144269.vst.exp <- assay(vsd)
# 
# dds <- DESeq(dds)
# GSE144269.diff.res <- results(dds)
# 
# save(cli.data, sams.info, GSE144269.vst.exp, GSE144269.diff.res, file='/data/GSE144269_diff_exp.RData')


#####################diff gene
load(file='/data/GSE144269_diff_exp.RData')
load(file='/data/exp_gene_anno.RData') 
candi.prog.gene <- readRDS(file='/data/candi_prog_sig_gene.rds')

rownames(GSE144269.vst.exp) <- stringr::str_split_i(rownames(GSE144269.vst.exp), pattern = '\\.', i = 1)
rownames(GSE144269.diff.res) <- stringr::str_split_i(rownames(GSE144269.diff.res), pattern = '\\.', i = 1)

GSE144269.vst.exp <- GSE144269.vst.exp[intersect(rownames(GSE144269.vst.exp), exp.gene.anno$EnsemblID), ]
GSE144269.diff.res <- GSE144269.diff.res[intersect(rownames(GSE144269.diff.res), exp.gene.anno$EnsemblID), ]

rownames(GSE144269.vst.exp) <- exp.gene.anno$Symbol[match(rownames(GSE144269.vst.exp), exp.gene.anno$EnsemblID)]
rownames(GSE144269.diff.res) <- exp.gene.anno$Symbol[match(rownames(GSE144269.diff.res), exp.gene.anno$EnsemblID)]


GSE144269.diff.res[candi.prog.gene, ] %>% as.data.frame()


# heatmap
setwd('/result/Section2')
colnames(sams.info) <- 'Type'
ExpClusteringHeatmap(GSE144269.vst.exp, candi.prog.gene, sams.info, 'mongolian_candi_gene_exp_heatmap')


