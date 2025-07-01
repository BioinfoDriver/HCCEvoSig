library(ggplot2)
library(dplyr)
library(cowplot)
library(reshape2)
library(gghalves)
library(ggpubr)


evoSigScore <- readRDS(file = '/data/tcga_risk_score.rds')
evoSigScore <- evoSigScore %>% column_to_rownames(var = 'patient_id')


tcia <- read.csv2("/data/TCIA-ClinicalData.tsv", header = TRUE, sep = '\t') 
tcia <- tcia %>% column_to_rownames(var = 'barcode') %>% select(ips_ctla4_neg_pd1_neg, 
                                                                ips_ctla4_neg_pd1_pos, ips_ctla4_pos_pd1_neg, ips_ctla4_pos_pd1_pos)
tcia <- type.convert(tcia, as.is = TRUE)



tcia <- merge(evoSigScore, tcia, by = 'row.names')[, 19:23]
tcia <- melt(data = tcia, id.vars = 'risk.categ', variable.name = 'IPS', value.name = 'score')

sapply(c('ips_ctla4_neg_pd1_neg', 'ips_ctla4_neg_pd1_pos', 'ips_ctla4_pos_pd1_neg', 'ips_ctla4_pos_pd1_pos'), function(ct){
  
  p.value <- wilcox.test(score~risk.categ, data = subset(tcia, IPS == ct), alternative = "two.sided")$p.value
  p.value <- setNames(p.value, ct)
  return(p.value)
})
# ips_ctla4_neg_pd1_neg.ips_ctla4_neg_pd1_neg ips_ctla4_neg_pd1_pos.ips_ctla4_neg_pd1_pos ips_ctla4_pos_pd1_neg.ips_ctla4_pos_pd1_neg 
# 8.908442e-06                                7.041216e-04                                2.383746e-04 
# ips_ctla4_pos_pd1_pos.ips_ctla4_pos_pd1_pos 
# 6.184125e-03

p1 <- ggplot(tcia, aes(x = risk.categ, y = score, fill = risk.categ)) +
  geom_violin(trim =TRUE) +
  geom_boxplot(width = 0.1, fill = "white", alpha = 0.7) + 
  labs(x = NULL, y = 'Score') + 
  geom_signif(comparisons = list(c("high risk", "low risk")),
              test = "wilcox.test", textsize = 4, step_increase = 0.1) + 
  scale_fill_brewer(palette="Blues") + facet_grid(cols = vars(IPS))+
  theme(
    panel.background = element_rect(fill = "white"),  # 白色背景
    panel.grid.major = element_blank(),               # 移除主要网格线
    panel.grid.minor = element_blank(),               # 移除次要网格线
    axis.line = element_line(color = "black"),        # 黑色坐标轴线
    plot.background = element_rect(fill = "white"),   # 整个绘图区域白色背景
    panel.border = element_blank(), # 移除面板边框
    legend.position = "none")



combinedPlots <- plot_grid(p1, labels = "AUTO", ncol = 1, nrow = 3)
ggsave("/result/Section6/TCIA_IPS_score.pdf", combinedPlots)



