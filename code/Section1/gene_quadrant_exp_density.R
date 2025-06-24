
load(file='/data/intra.inter.ith.quadrant.RData')

##################Losic
load(file='/data/losic.refseq.gene.exp.RData')

losic.gene.exp <- tran.exp$abundance
rownames(losic.gene.exp) <- stringr::str_split_i(rownames(losic.gene.exp), pattern = '_', i = 2)
losic.gene.exp <- losic.gene.exp[, grep(pattern = 'H', colnames(losic.gene.exp), value = TRUE)]


losic.gene.exp <- losic.ith.quadrant %>% rownames_to_column(var = 'Var1') %>% merge.data.frame(reshape2::melt(losic.gene.exp), by = 'Var1') 
losic.gene.exp <- losic.gene.exp %>% mutate(value = log2(value + 1))


ggplot(losic.gene.exp, aes(x = value, fill = quadrant)) + geom_density(alpha = 0.5)

ggplot(subset(losic.gene.exp, Var2 %in% c("H12.a", "H12.b", "H12.c", "H12.d", "H12.e")), aes(x = value, fill = quadrant)) +
  geom_density(alpha = 0.5) 


##################Shen
shen.gene.exp <- readRDS(file='/data/GSE136711_vst_exp.RData')
load(file='/data/exp_gene_anno.RData')

rownames(shen.gene.exp) <- exp.gene.anno$Symbol[match(rownames(shen.gene.exp), exp.gene.anno$GeneID)]

shen.gene.exp <- shen.ith.quadrant %>% rownames_to_column(var = 'Var1') %>% merge.data.frame(reshape2::melt(shen.gene.exp), by = 'Var1') 


ggplot(shen.gene.exp, aes(x = value, fill = quadrant)) + geom_density(alpha = 0.5) 

ggplot(subset(shen.gene.exp, Var2 %in% c("HCC04.IT2",  "HCC04.IT3",  "HCC04.IT4",  "HCC04.IT5")), aes(x = value, fill = quadrant)) +
  geom_density(alpha = 0.5) 


# Gene expression from Shi et al.
shi.gene.exp <- readRDS(file='/data/shi.lihc.gene.exp.rds')
rownames(shi.gene.exp) <- shi.gene.exp$GENE_SYMBOL
shi.gene.exp <- as.matrix(shi.gene.exp[, -c(1:3)])

shi.gene.exp <- shi.ith.quadrant %>% rownames_to_column(var = 'Var1') %>% merge.data.frame(reshape2::melt(shi.gene.exp), by = 'Var1') 


ggplot(shi.gene.exp, aes(x = value, fill = quadrant)) + geom_density(alpha = 0.5) 
ggplot(subset(shi.gene.exp, Var2 %in% c("H3.1", "H3.2", "H3.3", "H3.4", "H3.5")), aes(x = value, fill = quadrant)) +
  geom_density(alpha = 0.5) 


##################Renji
renji.gene.exp <- read.csv(file = '/data/Renji/Renji_cohort_MR_exp.matrix', 
                           sep='\t', row.names = 1, header=TRUE)

renji.gene.exp <- renji.gene.exp[, grep(pattern = 'T', colnames(renji.gene.exp), value = T)]
renji.gene.exp <- as.matrix(renji.gene.exp)

renji.gene.exp <- renji.ith.quadrant %>% rownames_to_column(var = 'Var1') %>% merge.data.frame(reshape2::melt(renji.gene.exp), by = 'Var1') 

ggplot(renji.gene.exp, aes(x = value, fill = quadrant)) + geom_density(alpha = 0.5) 
ggplot(subset(renji.gene.exp, Var2 %in% c("T02_01", "T02_02", "T02_03", "T02_04", "T02_05", "T02_06")), aes(x = value, fill = quadrant)) +
  geom_density(alpha = 0.5) 



