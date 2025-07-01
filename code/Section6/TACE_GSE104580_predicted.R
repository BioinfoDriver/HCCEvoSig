
library('GEOquery')
library('dplyr')
library('stringr')
library('pROC')
library('cowplot')

# load data
gse <- getGEO(filename = '/data/Kim_2022_GSE104580_HepatolCommun/GSE104580_series_matrix.txt.gz')
load(file='/data/lasso_cox_res.RData')


exprData <- exprs(gse)
phenoData <- pData(gse)
featureData <- fData(gse)

featureData <- featureData[, c("ID", "Gene Title", "Gene Symbol", "ENTREZ_GENE_ID")]
colnames(featureData) <- c("ID", "Gene_Title", "Gene_Symbol", "ENTREZ_GENE_ID")
featureData <- featureData %>% filter(across(everything(), ~!str_detect(., "///")))

# expr data
featureData <- subset(featureData, Gene_Symbol %in% names(active.k.vals))
exprData <- merge.data.frame(featureData[, 'Gene_Symbol', FALSE], exprData, by = 'row.names')
exprData <- lapply(split.data.frame(exprData, f = ~Gene_Symbol), function(geneExprData){
  
  geneExprData <- geneExprData[which.max(rowMeans(geneExprData[, -c(1:2)])), , FALSE]
  return(geneExprData)
})

exprData <- do.call(rbind.data.frame, exprData)
exprData <- exprData[, -c(1, 2)]


# risk evaluation
RiskEsti <- function(exp.dat, gene.set, risk.coef, cut.off=NULL){
  
  exp.dat <- as.matrix(exp.dat[gene.set, ])
  
  risk.score <- crossprod(exp.dat, matrix(risk.coef, nrow=length(risk.coef)))[, 1]
  if(!is.null(cut.off)){
    risk.categ <- ifelse(risk.score >= cut.off, 'high risk', 'low risk')
    
  }else{
    risk.categ <- ifelse(risk.score >= median(risk.score), 'high risk', 'low risk')
    
  }
  return(data.frame(risk.score, risk.categ))
}

riskScore <- RiskEsti(exp.dat=exprData, gene.set=names(active.k.vals), risk.coef=active.k.vals, cut.off=NULL)
cliData <- merge(phenoData[, 'subject subgroup:ch1',FALSE], riskScore, by='row.names')
colnames(cliData) <- c('patient_id', 'TACE', 'risk.score', 'risk.categ')




# plot
# wilcox.test(risk.score~TACE, data = cliData)$p.value
# [1] 3.144701e-10

p1 <- ggplot(cliData, aes(x = TACE, y = risk.score, fill = TACE)) +
  geom_violin(trim =TRUE) +
  geom_boxplot(width = 0.1, fill = "white", alpha = 0.7) + labs(x = NULL, y = 'EvoSig score') + 
  geom_signif(comparisons = list(c("TACE non-responders", "TACE responders")),
              test = "wilcox.test", textsize = 4, step_increase = 0.1) + 
  scale_fill_brewer(palette="Blues") + 
  theme(
    panel.background = element_rect(fill = "white"),  # 白色背景
    panel.grid.major = element_blank(),               # 移除主要网格线
    panel.grid.minor = element_blank(),               # 移除次要网格线
    axis.line = element_line(color = "black"),        # 黑色坐标轴线
    plot.background = element_rect(fill = "white"),   # 整个绘图区域白色背景
    panel.border = element_blank(), # 移除面板边框
    legend.position = "top")



# fisher.test(table(cliData$risk.categ, cliData$TACE))$p.value
# [1] 3.717338e-10
# fisher.test(table(cliData$risk.categ, cliData$TACE))$estimate
# odds ratio 
# 9.772727 

p2 <- ggplot(data = cliData %>% group_by(TACE, risk.categ) %>% summarise(n = n()), aes(x = risk.categ, y = n, fill = TACE)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_y_continuous(labels = scales::percent) +
  labs(x = NULL, y = "Percentage") + 
  scale_fill_brewer(palette="Blues") + 
  theme(
    panel.background = element_rect(fill = "white"),  # 白色背景
    panel.grid.major = element_blank(),               # 移除主要网格线
    panel.grid.minor = element_blank(),               # 移除次要网格线
    axis.line = element_line(color = "black"),        # 黑色坐标轴线
    plot.background = element_rect(fill = "white"),   # 整个绘图区域白色背景
    panel.border = element_blank(), # 移除面板边框
    legend.position = "top")




# Fitting Generalized Linear Models
cliData <- cliData %>% mutate(resp = ifelse(TACE == 'TACE responders', 1, 0))
model <- glm(resp ~ risk.score, data = cliData, family = binomial)
predProbs <- predict(model, type = "response")

rocObj <- roc(cliData$resp, predProbs)

# cat("AUC:", auc(rocObj), "\n")
# AUC: 0.8022821
# cat("95% CI:", ci(rocObj), "\n")
# 95% CI: 0.7279499 0.8022821 0.8766142 


rocDf <- data.frame(Sensitivity = rocObj$sensitivities, Specificity = rocObj$specificities)
rocDf$OneMinusSpec <- 1 - rocDf$Specificity

p3 <- ggplot(rocDf, aes(x = OneMinusSpec, y = Sensitivity)) +
  geom_line(color = "black", size = 1) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +
  geom_polygon(aes(x = OneMinusSpec, y = Sensitivity), 
               fill = "white", alpha = 0.2) +
  annotate("text", x = 0.6, y = 0.4, 
           label = paste("AUC =", round(auc(rocObj), 3)), 
           size = 5) +
  labs(x = "1 - Specificity (False Positive Rate)",
       y = "Sensitivity (True Positive Rate)") +
  theme_classic2() +
  coord_equal()


combinedPlots <- plot_grid(p1, p2, p3, labels = "AUTO", ncol = 2)
ggsave("/result/Section6/TACE_GSE104580_predicted.pdf", combinedPlots)



