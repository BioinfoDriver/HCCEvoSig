
load(file='/data/exp_gene_anno.RData')

gene.exp <- read.csv(file = 'data/GSE136711/GSE136711_readcount.txt', 
                     sep='\t', row.names = 1, header=TRUE)

gene.exp <- gene.exp[!duplicated(substr(rownames(gene.exp), 1, 15)), ]
rownames(gene.exp) <- substr(rownames(gene.exp), 1, 15)

gene.exp <- gene.exp[intersect(na.omit(exp.gene.anno$EnsemblID), rownames(gene.exp)), ]
rownames(gene.exp) <- exp.gene.anno$Symbol[match(rownames(gene.exp), exp.gene.anno$EnsemblID)]

gene.exp <- gene.exp[rowSums(gene.exp >= 1) >= 29, ]

# The variance stabilizing transformation
library(DESeq2)

# sample information
sam.info <- data.frame(patientID=substr(colnames(gene.exp), 1, 5), 
                       sampleID = colnames(gene.exp), regionID = substr(colnames(gene.exp), 7, 10))
rownames(sam.info) <- colnames(gene.exp)

dds <- DESeqDataSetFromMatrix(countData = gene.exp, colData=sam.info, design = ~ 1)
vsd <- vst(dds, blind = FALSE)
exp.vst <- assay(vsd)

saveRDS(exp.vst, file='/data/GSE136711_vst_exp.RData')

