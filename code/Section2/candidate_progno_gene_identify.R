

# load data
# differential gene expression selecction
load('/data/exp_gene_anno.RData')
load(file='/data/icgc_lihc_diff_exp.RData')

diff.exp.gene <- subset(icgc.paired.diff.res, abs(log2FoldChange) > 1 & padj < 0.05) # 2024(up), 1221(down)
diff.exp.gene.list <- rownames(diff.exp.gene) # 3245

diff.exp.gene.list <- exp.gene.anno$Symbol[match(diff.exp.gene.list, exp.gene.anno$GeneID)]



# clonal gene expression selection
clonal.exp.gene.list <- readRDS(file='/data/evo_genes.rds') # 580


# prognostic gene selecction
load('/data/univ_cox_pvalue.RData')
univ.cox.pvalue <- tcga.univ.cox.pvalue

prog.gene.list <- names(univ.cox.pvalue[univ.cox.pvalue < 0.05]) # 3495
# prog.gene.list <- names(univ.cox.pvalue[univ.cox.pvalue < 0.01]) # 1953

# univ.cox.pvalue <- p.adjust(univ.cox.pvalue, method='fdr')
# prog.gene.list <- names(univ.cox.pvalue[univ.cox.pvalue < 0.05]) # 1652
# prog.gene.list <- names(univ.cox.pvalue[univ.cox.pvalue < 0.01]) # 642




# candidate prognostic signature
candi.prog.sig.gene <- Reduce(intersect, list(diff.exp.gene.list, clonal.exp.gene.list, prog.gene.list))

saveRDS(candi.prog.sig.gene, file='/data/candi_prog_sig_gene.rds')


######Plot
library("VennDiagram")
library("RColorBrewer")

myCol <- brewer.pal(3, "Pastel2")
p <- venn.diagram(x = list(diff.exp.gene.list, clonal.exp.gene.list, prog.gene.list),
                  category.names = c("DE Set" , "Q1 Set" , "Surv Set"), filename=NULL, output=TRUE, 
                  fill=myCol, fontface = "bold", fontfamily = "sans")

pdf("/result/Section2/candidate_gene_venn_diagramm.pdf")
grid.draw(p)
dev.off()

