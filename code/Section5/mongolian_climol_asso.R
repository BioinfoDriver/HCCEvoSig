
# load data
lihc.risk.score <- readRDS(file='/data/mongolian_risk_score.rds')

# data prepare
lihc.risk.score$risk.categ <- factor(lihc.risk.score$risk.categ, levels=c('low risk', 'high risk'))

cli.sig.char <- lihc.risk.score[, c('LHC_ID', 'WES_T', 'age_bin', 'sex', 'smoker', 'alcohol', 'obesity', 'cirrhosis', 
                                    'FHx_LC', 'hbv', 'hcv', 'hdv', 'stage', 'tumor_size', 'multinodular', 'alb', 'bil', 'alt', 'afp', 'risk.categ')]


# Fisher's Exact Test
p.values <- lapply(colnames(cli.sig.char)[3:19], function(char){
  print(rowSums(table(cli.sig.char[, c(char, 'risk.categ')])))
  print(table(cli.sig.char[, c(char, 'risk.categ')]))
  
  or <- fisher.test(table(cli.sig.char[, c(char, 'risk.categ')]))$estimate
  p.value <- fisher.test(table(cli.sig.char[, c(char, 'risk.categ')]))$p.value
  return(c(or, p.value))
})

p.values <- as.data.frame(do.call(rbind, p.values))
colnames(p.values) <- c('OR', 'Pvalue')
p.values$Qvalue <- p.adjust(p.values$Pvalue, 'fdr')
rownames(p.values) <- colnames(cli.sig.char)[3:19]


# OR      Pvalue     Qvalue
# tumor_size   3.6032119 0.048318263 0.41070524
# afp          6.7834992 0.001136466 0.01931992

####Plot
library('ggpubr')

tumorSize.p <- ggbarplot(as.data.frame(table(cli.sig.char[, c('tumor_size', 'risk.categ')])), 
                         "risk.categ", "Freq", fill = "tumor_size", color="tumor_size", ylab='Number of patients', xlab=FALSE, label=TRUE, 
                         lab.col="white", lab.pos="in", lab.size=2, position=position_dodge(0.9),
                         title='tumorSize', lab.vjust=2, palette="Paired") + theme(axis.text.y = element_text(size=7), axis.text.x = element_text(size=7))

afp.p <- ggbarplot(as.data.frame(table(cli.sig.char[, c('afp', 'risk.categ')])), 
                   "risk.categ", "Freq", fill = "afp", color="afp", ylab='Number of patients', xlab=FALSE, label=TRUE, 
                   lab.col="white", lab.pos="in", lab.size=2, position=position_dodge(0.9),
                   title='AFP', lab.vjust=2, palette="Paired") + theme(axis.text.y = element_text(size=7), axis.text.x = element_text(size=7))


ggsave(ggarrange(tumorSize.p, afp.p, ncol=4, nrow=3, legend = "none"), 
       file='/result/Section5/mongolian_cli_risk_com.pdf')




##########################################Mutation data
library('maftools')

colnames(cli.sig.char)[2] <- 'Tumor_Sample_Barcode'
mut.maf = read.maf(maf = '/data/OriginalData/GSE144269/mutect2_merged.maf', 
                   clinicalData = cli.sig.char, verbose = TRUE)


# Changing colors for variant classifications
vc_cols = RColorBrewer::brewer.pal(n = 8, name = 'Paired')
names(vc_cols) = c('Frame_Shift_Del', 'Missense_Mutation', 'Nonsense_Mutation', 'Multi_Hit',
                   'Frame_Shift_Ins', 'In_Frame_Ins', 'Splice_Site', 'In_Frame_Del')

show_genes <- c('TP53', 'NFE2L2','ACVR2A','KEAP1','ARID2','BAP1','RB1','AXIN1','ARID1A','APOB','ALB','CTNNB1',
                'BRD7', 'CDKN1A', 'CDKN2A', 'CREB3L3', 'DHX9', 'EEF1A1', 'GNAS','HIST1H1E', 'IDH1', 'IL6ST', 'KRAS', 'LZTR1', 'NRAS',
                'NUP113','PIK3CA', 'RPS6KA3', 'SMARCA4', 'TSC2', 'WHSC1', 'XPO1')

cols = RColorBrewer::brewer.pal(n = 8,name = 'Spectral')
risk.cols <- cols[1:2]
names(risk.cols) <- c('high risk', 'low risk')
stage.cols <- cols[3:4]
names(stage.cols) <- c('0', '1')
ts.cols <- cols[5:6]
names(ts.cols) <- c('0', '1')
afp.cols <- cols[7:8]
names(afp.cols) <- c('0', '1')


annocolors = list(risk.categ = risk.cols, stage = stage.cols, tumor_size = ts.cols, afp = afp.cols)

pdf('/result/Section5/mongolian_mut_oncoplot.pdf')
oncoplot(maf = mut.maf, genes = show_genes, colors = vc_cols, # minMut=9,
         clinicalFeatures = c('risk.categ','stage', 'tumor_size','afp'),
         annotationColor = annocolors, sortByAnnotation = TRUE)
dev.off()


# Clinical enrichment analysis
fab.ce = clinicalEnrichment(maf = mut.maf, clinicalFeature = 'risk.categ', minMut = 6)
# fab.ce$pairwise_comparision[fdr < 0.05]

pdf('/result/Section5/mongolian_mut_enrich.pdf')
plotEnrichmentResults(enrich_res = fab.ce, pVal = 0.03, geneFontSize = 0.5, annoFontSize = 0.6)
dev.off()



# Fisher's Exact Test
lihc.mut.maf <- read.csv(file='/data/OriginalData/GSE144269/mutect2_merged.maf', 
                         header=TRUE, sep='\t', stringsAsFactors=FALSE)

lihc.mut.maf <- subset(lihc.mut.maf, !Variant_Classification %in% 
                         c("3'Flank", "3'UTR", "5'Flank","5'UTR", "IGR", "Intron", "RNA", "Silent"))

mut.genes <- unique(lihc.mut.maf$Hugo_Symbol)
gene.mut.mat <- sapply(unique(lihc.mut.maf$Tumor_Sample_Barcode), function(patient.id){
  gene.sym <- subset(lihc.mut.maf, Tumor_Sample_Barcode == patient.id)$Hugo_Symbol
  index <- as.numeric(mut.genes %in% gene.sym)
  return(index)
})

rownames(gene.mut.mat) <- mut.genes
high.fre.mut <- gene.mut.mat[rowSums(gene.mut.mat)/ncol(gene.mut.mat) >= 0.03, ]
high.fre.mut <- as.data.frame(t(high.fre.mut))

cli.sig.char <- subset(cli.sig.char, !is.na(Tumor_Sample_Barcode))
high.fre.mut <- high.fre.mut[cli.sig.char$Tumor_Sample_Barcode, ]

# driver genes with 3% mutation frequence
high.fre.mut <- high.fre.mut[, intersect(show_genes, colnames(high.fre.mut))]
high.fre.mut <- cbind(high.fre.mut, cli.sig.char[, 'risk.categ', FALSE])


p.values <- sapply(setdiff(colnames(high.fre.mut), 'risk.categ'), function(driver.gene){
  
  num <- table(high.fre.mut[, driver.gene])
  char <- names(num)[num >= 3]
  dat <- subset(high.fre.mut, get(driver.gene) %in% char)
  
  if(dplyr::n_distinct(dat[, driver.gene])<=1){
    p <- NA
    
  }else{
    stat <- table(dat[, c(driver.gene, 'risk.categ')])
    p <- fisher.test(stat, alternative = "two.sided")$p.value
    
  }
  
  return(p)
  
})

# TP53     NFE2L2     ACVR2A      KEAP1      ARID2       BAP1        RB1      AXIN1     ARID1A       APOB        ALB     CTNNB1 
# 0.02124325 0.41327388 1.00000000 0.46309160 1.00000000 0.02519629 1.00000000 0.24426479 0.20372044 0.71089548 0.22188651 0.79478129 
# CDKN2A       DHX9       GNAS   HIST1H1E      IL6ST      LZTR1     PIK3CA    RPS6KA3    SMARCA4       TSC2       XPO1 
# 0.60185440 1.00000000 0.73589324 0.18377695 0.72611150 0.02519629 0.43011770 0.71089548 1.00000000 1.00000000 0.76806837

p.adjust(p.values, method='fdr')
# TP53    NFE2L2    ACVR2A     KEAP1     ARID2      BAP1       RB1     AXIN1    ARID1A      APOB       ALB    CTNNB1    CDKN2A 
# 0.1931715 1.0000000 1.0000000 1.0000000 1.0000000 0.1931715 1.0000000 0.8025843 0.8025843 1.0000000 0.8025843 1.0000000 1.0000000 
# DHX9      GNAS  HIST1H1E     IL6ST     LZTR1    PIK3CA   RPS6KA3   SMARCA4      TSC2      XPO1 
# 1.0000000 1.0000000 0.8025843 1.0000000 0.1931715 1.0000000 1.0000000 1.0000000 1.0000000 1.0000000

