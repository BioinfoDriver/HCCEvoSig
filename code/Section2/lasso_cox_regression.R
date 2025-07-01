
# load data
load('/data/exp_gene_anno.RData')
tcga.lihc.vst.exp <- readRDS(file='/data/tcga_vst_norm_geneExp.rds')
tcga.lihc.cli.data <- readRDS(file='/data/tcga.lihc.cli.data.rds')
candi.prog.sig.gene <- readRDS(file='/data/candi_prog_sig_gene.rds')


# data integration
rownames(tcga.lihc.vst.exp) <- exp.gene.anno$Symbol[match(rownames(tcga.lihc.vst.exp), exp.gene.anno$GeneID)]


colnames(tcga.lihc.vst.exp) <- substr(colnames(tcga.lihc.vst.exp), 1, 12)
tcga.lihc.vst.exp <- t(tcga.lihc.vst.exp)

tcga.lihc.cli.data <- tcga.lihc.cli.data[, c('os_time', 'os')]
tcga.lihc.cli.data <- as.matrix(tcga.lihc.cli.data)



# Regularized Cox Regression
LassoCoxFunction <- function(candi.gene, exp.data, cli.data, n.cv, t.measure=c('deviance', 'C'), file.name){
  library(glmnet)
  library(survival)
  
  t.measure = match.arg(t.measure)
  
  colnames(cli.data) <- c('time', 'status')
  exp.data <- exp.data[, candi.gene]
  
  inters.sams <- intersect(rownames(exp.data), rownames(cli.data))
  exp.data <- exp.data[inters.sams, ]
  cli.data <- cli.data[inters.sams, ]
  
  
  # Does k-fold cross-validation for glmnet,
  # set.seed(123)
  cv.fit = cv.glmnet(x = exp.data, y = cli.data, type.measure = t.measure, nfolds=n.cv, family = "cox")
  
  
  # plot the cross-validation curve
  # pdf(file.name)
  # plot(cv.fit)
  # dev.off()
  
  return(cv.fit)
}


random.lasso.cox.res <- lapply(seq(1000), function(i){
  
  exp.data <- tcga.lihc.vst.exp[sample(seq(nrow(tcga.lihc.vst.exp)), round(nrow(tcga.lihc.vst.exp) * 0.75)), ]
  lasso.cox.res <- LassoCoxFunction(candi.prog.sig.gene, exp.data, 
                                    tcga.lihc.cli.data, n.cv=10, t.measure='deviance', 
                                    file.name='/result/Section2/cv_curve_deviance.pdf')
  
  return(lasso.cox.res)
})



# extract non-zero coefficients
gene.appear <-lapply(random.lasso.cox.res, function(lasso.cox.res){
  est.coef = coef(lasso.cox.res, s = lasso.cox.res$lambda.min)
  active.k.vals = est.coef[which(est.coef != 0), ]
  
  return(names(active.k.vals))
})

# gene.appear.1se <-lapply(random.lasso.cox.res, function(lasso.cox.res){
#   est.coef = coef(lasso.cox.res, s = lasso.cox.res$lambda.1se)
#   active.k.vals.1se = est.coef[which(est.coef != 0), ]
#   
#   return(names(active.k.vals.1se))
# })



gene.freq <- as.data.frame(table(unlist(gene.appear))) %>% arrange(desc(Freq)) %>% mutate(Var1 = factor(Var1, levels = Var1))
my_color <- colorRampPalette(c("#FF0000", "#FFFFFF"), bias = 5)(21)
names(my_color) <- gene.freq$Var1

p1 <- ggplot(data = gene.freq, aes(x = Var1, y = Freq, fill = Var1)) + geom_bar(stat = "identity", show.legend = FALSE) + 
  theme_minimal() + scale_fill_manual(values = my_color) + 
  geom_text(aes(label = Freq), vjust = -0.5, size = 4) + 
  theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 30, hjust = 1)) + labs(y = "Frequency")
  
ggsave(plot = p1, filename = '/result/Section2/geneFreq.pdf')



# Multivariate Cox regression
gene.appear <- table(unlist(gene.appear))
gene.appear <- names(gene.appear[gene.appear > 850])

cox.model <- coxph(Surv(os_time, os) ~ ADH4 + CDC20 + CFHR3 + CYP2C9 + RAMP3 + RDH16 + SERPINE1 + SLC16A11 + SPP1 + SPP2 + TRNP1, 
                   data = merge(tcga.lihc.vst.exp, tcga.lihc.cli.data, by = 'row.names'))
active.k.vals <- coef(cox.model)

save(random.lasso.cox.res, active.k.vals, file='/data/lasso_cox_res.RData')


