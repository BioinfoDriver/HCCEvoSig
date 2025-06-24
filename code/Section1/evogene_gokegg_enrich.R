
library(clusterProfiler)
library(org.Hs.eg.db)
# load datahttp://127.0.0.1:10633/graphics/plot_zoom_png?width=1920&height=1017
evo.genes <- readRDS(file='/data/evo_genes.rds')
evo.gene.ids <- bitr(evo.genes, fromType="SYMBOL", toType=c("ENTREZID"), OrgDb="org.Hs.eg.db")

# GO,KEGG pathway over-representation analysis
ego <- enrichGO(gene = evo.genes, keyType = 'SYMBOL', OrgDb = org.Hs.eg.db, ont = "BP", pAdjustMethod = "BH", pvalueCutoff  = 0.05, qvalueCutoff  = 0.05, readable = TRUE)

ekegg <- enrichKEGG(gene = evo.gene.ids$ENTREZID, keyType = "ncbi-geneid", pAdjustMethod = "BH", organism = 'hsa', pvalueCutoff = 0.05, qvalueCutoff = 0.05)
ekegg <- setReadable(ekegg, OrgDb = org.Hs.eg.db, keyType="ENTREZID")


save(ego, ekegg, file = '/data/evo_genes_enrich.RData')


# write.table(ego, file = '/result/section1/evo_gene_goenrich.txt', sep = '\t', col.names = T, row.names = F, quote = F)
# write.table(ekegg, file = '/result/section1/evo_gene_keggenrich.txt', sep = '\t', col.names = T, row.names = F, quote = F)


# Enrichment Map
goEnrichMap <- pairwise_termsim(ego, showCategory = 50)
ggsave(plot = emapplot(goEnrichMap, cex_category=1.5, showCategory = 50), 
       filename = '/result/section1/goEnrichMap.pdf', width = 14, height = 14)


keggEnrichMap <- pairwise_termsim(ekegg, showCategory = 50)
ggsave(plot = emapplot(keggEnrichMap, cex_category=1.5, showCategory = 50), 
       filename = '/result/section1/keggEnrichMap.pdf', width = 14, height = 14)


# Dot plot

showKEGG <- c("Retinol metabolism"                    ,                          
"Bile secretion"                                      , "Steroid hormone biosynthesis"              ,                
"Biosynthesis of amino acids"                         ,
"Tyrosine metabolism"                                 , "ECM-receptor interaction"                  ,                
"Arginine biosynthesis"                               , "Complement and coagulation cascades"       ,       
"Primary bile acid biosynthesis"                      , "Focal adhesion"                            ,                
"Fructose and mannose metabolism"                     , "Cholesterol metabolism"                    ,                
"Cysteine and methionine metabolism"                  , "Glycolysis / Gluconeogenesis"              ,                
"TNF signaling pathway"                               , "Glycine, serine and threonine metabolism"  ,                
"Protein digestion and absorption"                    , "Glutathione metabolism"                    ,                
"PI3K-Akt signaling pathway"                          ,
"Cytokine-cytokine receptor interaction"              , "Pyruvate metabolism") 


ggsave(plot = dotplot(ekegg, showCategory = showKEGG), 
       filename = '/result/section1/keggDotplot.pdf')


showGO <- c("steroid metabolic process"                    ,  
"olefinic compound metabolic process"            ,"carboxylic acid biosynthetic process"       ,  
"organic acid biosynthetic process"              ,"fatty acid metabolic process"               ,"hormone metabolic process",  
"primary alcohol metabolic process"              ,  
"terpenoid metabolic process"                    ,  
"alcohol metabolic process"                      ,"cell chemotaxis"                            ,  
"icosanoid metabolic process"                    ,  
"small molecule catabolic process"               ,"amino acid metabolic process"               ,  
"positive regulation of chemotaxis"              ,  
"leukocyte chemotaxis"                           ,"positive regulation of leukocyte migration" ,  
"isoprenoid metabolic process"                   ,"monosaccharide biosynthetic process"        ,  
"gluconeogenesis"                                ,  
"positive regulation of leukocyte chemotaxis"    ,  
"steroid biosynthetic process")

ggsave(plot = dotplot(ego, showCategory = showGO), 
       filename = '/result/section1/goDotplot.pdf')

