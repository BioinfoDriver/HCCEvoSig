
library('ggpubr')
library('ggplot2')
library('limma')
library('clusterProfiler')
library('org.Hs.eg.db')
library('enrichplot')

# load data
tcga.lihc.cli.data <- readRDS(file='/data/tcga_risk_score.rds')
load(file='/data/tcga_lihc_rnaseq.RData')

tcga.lihc.cli.data$patient_id <- paste0(tcga.lihc.cli.data$patient_id, '-01')
tcga.lihc.exp <- tcga.lihc.tpm[, tcga.lihc.cli.data$patient_id]
tcga.lihc.exp <- log2(tcga.lihc.exp + 1)

rownames(tcga.lihc.exp) <- stringr::str_split_i(rownames(tcga.lihc.exp), pattern = '\\|', i = 2)

# diff gene
tcga.lihc.cli.data$risk.categ <- factor(tcga.lihc.cli.data$risk.categ, levels = c('low risk', 'high risk'))
design <- model.matrix(~risk.categ, data = tcga.lihc.cli.data)


fit <- lmFit(tcga.lihc.exp, design)
fit <- eBayes(fit)
difRes <- topTable(fit, coef="risk.categhigh risk", adjust.method = "BH", n = Inf)



# GO,KEGG pathway over-representation analysis
diff.genes <- rownames(subset(difRes, abs(logFC) > 1.0 & adj.P.Val < 0.05))

ego <- enrichGO(gene = diff.genes, keyType = 'ENTREZID', OrgDb = org.Hs.eg.db, ont = "BP", pAdjustMethod = "BH", pvalueCutoff  = 0.05, qvalueCutoff  = 0.05, readable = TRUE)

ekegg <- enrichKEGG(gene = diff.genes, keyType = "ncbi-geneid", pAdjustMethod = "BH", organism = 'hsa', pvalueCutoff = 0.05, qvalueCutoff = 0.05)
ekegg <- setReadable(ekegg, OrgDb = org.Hs.eg.db, keyType="ENTREZID")


# KEGG Gene Set Enrichment Analysis
gene.list <- setNames(difRes$logFC, rownames(difRes))
gene.list = sort(gene.list, decreasing = TRUE)

gseakegg <- gseKEGG(geneList = gene.list, organism = 'hsa',
                    keyType = "ncbi-geneid", minGSSize = 10, maxGSSize = 500, pvalueCutoff = 0.05,
                    pAdjustMethod = "BH", verbose = FALSE)

gseakegg <- setReadable(gseakegg, OrgDb = org.Hs.eg.db, keyType="ENTREZID")

# save(difRes, ego, ekegg, gseakegg, file = '/data/tcga_lihc_gokegg_enrich.RData')


################ plot
# write.table(ego, file = '/result/Section5/goenrich.txt', sep = '\t', col.names = T, row.names = F, quote = F)
# write.table(ekegg, file = '/result/Section5/keggenrich.txt', sep = '\t', col.names = T, row.names = F, quote = F)
# write.table(gseakegg, file = '/result/Section5/kegggsea.txt', sep = '\t', col.names = T, row.names = F, quote = F)


categorys <- c('Retinol metabolism',
               'Metabolism of xenobiotics by cytochrome P450',
               'Fatty acid degradation',
               'Steroid hormone biosynthesis',
               'Tryptophan metabolism',
               'Bile secretion',
               'Valine, leucine and isoleucine degradation',
               'Tyrosine metabolism',
               'Primary bile acid biosynthesis',
               'Glycolysis / Gluconeogenesis',
               'Biosynthesis of amino acids',
               'Central carbon metabolism in cancer',
               'Cholesterol metabolism',
               'Fatty acid metabolism',
               'Cell cycle',
               'DNA replication', 
               'p53 signaling pathway',
               'Tight junction',
               'TNF signaling pathway',
               'IL-17 signaling pathway',
               'Hippo signaling pathway',
               'Regulation of actin cytoskeleton')


pdf('/result/Section5/gsea_kegg.pdf')
dotplot(gseakegg, showCategory=categorys, split=".sign") + facet_grid(.~.sign)

gseaplot2(gseakegg, geneSetID = c(1, 32, 49), pvalue_table = TRUE,
          color = c("#E495A5", "#86B875", "#7DB0DD"), ES_geom = "line")

gseaplot2(gseakegg, geneSetID = c(7, 19, 68), pvalue_table = TRUE,
          color = c("#E495A5", "#86B875", "#7DB0DD"), ES_geom = "line")

gseaplot2(gseakegg, geneSetID = c(13, 25, 71), pvalue_table = TRUE,
          color = c("#E495A5", "#86B875", "#7DB0DD"), ES_geom = "line")

gseaplot2(gseakegg, geneSetID = c(74, 93, 103), pvalue_table = TRUE,
          color = c("#E495A5", "#86B875", "#7DB0DD"), ES_geom = "line")

dev.off()
