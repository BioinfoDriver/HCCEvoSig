
library(ggplot2)
library(dplyr)
library(cowplot)

# read data
tideHCC <- read.csv2("/data/tide_result.csv", header = TRUE, sep = ',', row.names = 1) 
rownames(tideHCC) <- str_sub(rownames(tideHCC), 1, 12)
tideHCC <- type.convert(tideHCC, as.is = TRUE)


evoSigScore <- readRDS(file = '/data/tcga_risk_score.rds')
evoSigScore <- evoSigScore %>% column_to_rownames(var = 'patient_id')

evoSigScore <- merge(evoSigScore, tideHCC, by = 'row.names')


# wilcox.test(Dysfunction ~ risk.categ, data = evoSigScore, alternative = "two.sided")$p.value
# [1] 0.01724197

# wilcox.test(Exclusion ~ risk.categ, data = evoSigScore, alternative = "two.sided")$p.value
# [1] 9.260288e-05

# wilcox.test(TIDE ~ risk.categ, data = evoSigScore, alternative = "two.sided")$p.value
# [1] 0.0008908428

p1 <- ggplot(evoSigScore, aes(x = risk.categ, y = TIDE, fill = risk.categ)) +
  geom_violin(trim =TRUE) +
  geom_boxplot(width = 0.1, fill = "white", alpha = 0.7) + labs(x = NULL, y = 'TIDE score') + 
  geom_signif(comparisons = list(c("high risk", "low risk")),
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


p2 <- ggplot(evoSigScore, aes(x = risk.categ, y = Dysfunction, fill = risk.categ)) +
  geom_violin(trim =TRUE) +
  geom_boxplot(width = 0.1, fill = "white", alpha = 0.7) + labs(x = NULL, y = 'Dysfunction score') + 
  geom_signif(comparisons = list(c("high risk", "low risk")),
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


p3 <- ggplot(evoSigScore, aes(x = risk.categ, y = Exclusion, fill = risk.categ)) +
  geom_violin(trim =TRUE) +
  geom_boxplot(width = 0.1, fill = "white", alpha = 0.7) + labs(x = NULL, y = 'Exclusion score') + 
  geom_signif(comparisons = list(c("high risk", "low risk")),
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



# cor.test(evoSigScore$TIDE, evoSigScore$risk.score, method = 'spearman')$estimate
# 0.2712361
# cor.test(evoSigScore$TIDE, evoSigScore$risk.score, method = 'spearman')$p.value
# [1] 8.24972e-07


# cor.test(evoSigScore$Dysfunction, evoSigScore$risk.score, method = 'spearman')$estimate
# -0.1074803
# cor.test(evoSigScore$Dysfunction, evoSigScore$risk.score, method = 'spearman')$p.value
# [1] 0.05367028


# cor.test(evoSigScore$Exclusion, evoSigScore$risk.score, method = 'spearman')$estimate
# 0.29589
# cor.test(evoSigScore$Exclusion, evoSigScore$risk.score, method = 'spearman')$p.value
# [1] 6.931023e-08



p4 <- ggplot(evoSigScore, aes(x = risk.score, y = TIDE)) +
  geom_point(
    size = 1,            # 点的大小
    shape = 21,          # 点的形状(21是带边框的圆)
    color = "black",     # 边框颜色
    fill = "#1c61b6",    # 填充颜色
    alpha = 0.7          # 透明度(0-1)
  ) +
  geom_smooth(
    method = "lm",       # 线性模型
    se = TRUE,           # 显示置信区间
    color = "#ff7f00",   # 线条颜色
    fill = "#ff7f00",    # 置信区间填充色
    alpha = 0.2          # 置信区间透明度
  ) + stat_cor(method = "spearman", 
               label.x = -4, label.y = 2,
               color = "red",
               label.sep = "\n") + 
  labs(x = 'Evosig score', y = 'TIDE score') + 
  theme_classic() 


p5 <- ggplot(evoSigScore, aes(x = risk.score, y = Dysfunction)) +
  geom_point(
    size = 1,            # 点的大小
    shape = 21,          # 点的形状(21是带边框的圆)
    color = "black",     # 边框颜色
    fill = "#1c61b6",    # 填充颜色
    alpha = 0.7          # 透明度(0-1)
  ) +
  geom_smooth(
    method = "lm",       # 线性模型
    se = TRUE,           # 显示置信区间
    color = "#ff7f00",   # 线条颜色
    fill = "#ff7f00",    # 置信区间填充色
    alpha = 0.2          # 置信区间透明度
  ) + stat_cor(method = "spearman", 
               label.x = -4, label.y = 2,
               color = "red",
               label.sep = "\n") + 
  labs(x = 'Evosig score', y = 'Dysfunction score') + 
  theme_classic() 


p6 <- ggplot(evoSigScore, aes(x = risk.score, y = Exclusion)) +
  geom_point(
    size = 1,            # 点的大小
    shape = 21,          # 点的形状(21是带边框的圆)
    color = "black",     # 边框颜色
    fill = "#1c61b6",    # 填充颜色
    alpha = 0.7          # 透明度(0-1)
  ) +
  geom_smooth(
    method = "lm",       # 线性模型
    se = TRUE,           # 显示置信区间
    color = "#ff7f00",   # 线条颜色
    fill = "#ff7f00",    # 置信区间填充色
    alpha = 0.2          # 置信区间透明度
  ) + stat_cor(method = "spearman", 
               label.x = -4, label.y = 2,
               color = "red",
               label.sep = "\n") + 
  labs(x = 'Evosig score', y = 'Exclusion score') + 
  theme_classic() 


combinedPlots <- plot_grid(p1, p4, p2, p5, p3, p6, labels = "AUTO", ncol = 2, nrow = 3)
ggsave("/result/Section6/TIDE_tcga.pdf", combinedPlots)


