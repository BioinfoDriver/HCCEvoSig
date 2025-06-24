
library("VennDiagram")
library("RColorBrewer")
library(dplyr)
library(ggplot2)

### load
load(file='/data/intra.inter.ith.score.RData')
load(file='/data/intra.inter.ith.quadrant.RData')


# overlap
{
pdf("/result/section1/top_right_overlap.pdf")

myCol <- brewer.pal(4, "Pastel2")
p <- venn.diagram(x = list(rownames(subset(losic.ith.quadrant, quadrant == 'top_right')),
                           rownames(subset(renji.ith.quadrant, quadrant == 'top_right')),
                           rownames(subset(shen.ith.quadrant, quadrant == 'top_right')),
                           rownames(subset(shi.ith.quadrant, quadrant == 'top_right'))),
                  category.names = c("Q1 Set(Losic et al.)", "Q1 Set(Yang et al.)", "Q1 Set(Shen et al.)", "Q1 Set(Shi et al.)"), filename=NULL, output=TRUE, 
                  fill=myCol, fontface = "bold", fontfamily = "sans")
grid.draw(p)
grid.newpage()

p <- venn.diagram(x = list(rownames(subset(losic.ith.quadrant, quadrant == 'top_right')),
                           rownames(subset(renji.ith.quadrant, quadrant == 'top_right'))),
                  category.names = c("Q1 Set(Losic et al.)", "Q1 Set(Yang et al.)"), filename=NULL, output=TRUE, 
                  fill=myCol[c(1, 2)], fontface = "bold", fontfamily = "sans")
grid.draw(p)
grid.newpage()


p <- venn.diagram(x = list(rownames(subset(losic.ith.quadrant, quadrant == 'top_right')),
                           rownames(subset(shen.ith.quadrant, quadrant == 'top_right'))),
                  category.names = c("Q1 Set(Losic et al.)", "Q1 Set(Shen et al.)"), filename=NULL, output=TRUE, 
                  fill=myCol[c(1, 3)], fontface = "bold", fontfamily = "sans")
grid.draw(p)
grid.newpage()

p <- venn.diagram(x = list(rownames(subset(losic.ith.quadrant, quadrant == 'top_right')),
                           rownames(subset(shi.ith.quadrant, quadrant == 'top_right'))),
                  category.names = c("Q1 Set(Losic et al.)", "Q1 Set(Shi et al.)"), filename=NULL, output=TRUE, 
                  fill=myCol[c(1, 4)], fontface = "bold", fontfamily = "sans")
grid.draw(p)
grid.newpage()

p <- venn.diagram(x = list(rownames(subset(renji.ith.quadrant, quadrant == 'top_right')),
                           rownames(subset(shen.ith.quadrant, quadrant == 'top_right'))),
                  category.names = c("Q1 Set(Yang et al.)", "Q1 Set(Shen et al.)"), filename=NULL, output=TRUE, 
                  fill=myCol[c(2, 3)], fontface = "bold", fontfamily = "sans")
grid.draw(p)
grid.newpage()

p <- venn.diagram(x = list(rownames(subset(renji.ith.quadrant, quadrant == 'top_right')),
                           rownames(subset(shi.ith.quadrant, quadrant == 'top_right'))),
                  category.names = c("Q1 Set(Yang et al.)", "Q1 Set(Shi et al.)"), filename=NULL, output=TRUE, 
                  fill=myCol[c(2, 4)], fontface = "bold", fontfamily = "sans")
grid.draw(p)
grid.newpage()

p <- venn.diagram(x = list(rownames(subset(shen.ith.quadrant, quadrant == 'top_right')),
                           rownames(subset(shi.ith.quadrant, quadrant == 'top_right'))),
                  category.names = c("Q1 Set(Shen et al.)", "Q1 Set(Shi et al.)"), filename=NULL, output=TRUE, 
                  fill=myCol[c(3, 4)], fontface = "bold", fontfamily = "sans")
grid.draw(p)

dev.off()
}


# 1 - phyper(q = length(intersect(rownames(subset(losic.ith.quadrant, quadrant == 'top_right')), rownames(subset(renji.ith.quadrant, quadrant == 'top_right')))), 
#            m = nrow(subset(renji.ith.quadrant, quadrant == 'top_right')), 
#            n = nrow(subset(renji.ith.quadrant, quadrant != 'top_right')), 
#            k = nrow(subset(losic.ith.quadrant, quadrant == 'top_right')), lower.tail = T, log.p = FALSE)
# 
# 
# 1 - phyper(q = length(intersect(rownames(subset(losic.ith.quadrant, quadrant == 'top_right')), rownames(subset(shen.ith.quadrant, quadrant == 'top_right')))), 
#            m = nrow(subset(shen.ith.quadrant, quadrant == 'top_right')), 
#            n = nrow(subset(shen.ith.quadrant, quadrant != 'top_right')), 
#            k = nrow(subset(losic.ith.quadrant, quadrant == 'top_right')), lower.tail = T, log.p = FALSE)
# 
# 
# 1 - phyper(q = length(intersect(rownames(subset(losic.ith.quadrant, quadrant == 'top_right')), rownames(subset(shi.ith.quadrant, quadrant == 'top_right')))), 
#            m = nrow(subset(shi.ith.quadrant, quadrant == 'top_right')), 
#            n = nrow(subset(shi.ith.quadrant, quadrant != 'top_right')), 
#            k = nrow(subset(losic.ith.quadrant, quadrant == 'top_right')), lower.tail = T, log.p = FALSE)
# 
# 
# 1 - phyper(q = length(intersect(rownames(subset(renji.ith.quadrant, quadrant == 'top_right')), rownames(subset(shen.ith.quadrant, quadrant == 'top_right')))), 
#            m = nrow(subset(shen.ith.quadrant, quadrant == 'top_right')), 
#            n = nrow(subset(shen.ith.quadrant, quadrant != 'top_right')), 
#            k = nrow(subset(renji.ith.quadrant, quadrant == 'top_right')), lower.tail = T, log.p = FALSE)
# 
# 
# 1 - phyper(q = length(intersect(rownames(subset(renji.ith.quadrant, quadrant == 'top_right')), rownames(subset(shi.ith.quadrant, quadrant == 'top_right')))), 
#            m = nrow(subset(shi.ith.quadrant, quadrant == 'top_right')), 
#            n = nrow(subset(shi.ith.quadrant, quadrant != 'top_right')), 
#            k = nrow(subset(renji.ith.quadrant, quadrant == 'top_right')), lower.tail = T, log.p = FALSE)
# 
# 
# 1 - phyper(q = length(intersect(rownames(subset(shen.ith.quadrant, quadrant == 'top_right')), rownames(subset(shi.ith.quadrant, quadrant == 'top_right')))), 
#            m = nrow(subset(shi.ith.quadrant, quadrant == 'top_right')), 
#            n = nrow(subset(shi.ith.quadrant, quadrant != 'top_right')), 
#            k = nrow(subset(shen.ith.quadrant, quadrant == 'top_right')), lower.tail = T, log.p = FALSE)
# 
# 




# scatter diagram
evo.genes <- Reduce(intersect, list(rownames(subset(losic.ith.quadrant, quadrant == 'top_right')),
                       rownames(subset(renji.ith.quadrant, quadrant == 'top_right')),
                       rownames(subset(shen.ith.quadrant, quadrant == 'top_right')),
                       rownames(subset(shi.ith.quadrant, quadrant == 'top_right'))))

# saveRDS(evo.genes, file='D:/LIHCSurvival/data/evo_genes.rds')


losic.ith.score$group <- factor(ifelse(rownames(losic.ith.score) %in% evo.genes, 'Q1 Gene', 'Other'), levels = c('Q1 Gene', 'Other'))
renji.ith.score$group <- factor(ifelse(rownames(renji.ith.score) %in% evo.genes, 'Q1 Gene', 'Other'), levels = c('Q1 Gene', 'Other'))
shen.ith.score$group <- factor(ifelse(rownames(shen.ith.score) %in% evo.genes, 'Q1 Gene', 'Other'), levels = c('Q1 Gene', 'Other'))
shi.ith.score$group <- factor(ifelse(rownames(shi.ith.score) %in% evo.genes, 'Q1 Gene', 'Other'), levels = c('Q1 Gene', 'Other'))


# MAD
{
pdf("/result/section1/top_right_mad_score.pdf")
losic.mad <- ggplot(losic.ith.score, aes(x = losic.intra.mad, y = losic.inter.mad, color = group)) + 
  geom_point(size = 1) + 
  scale_color_manual(values = c("Q1 Gene" = "orange", "Other" = "black"), name = "Group") + 
  labs(title = "Losic et al.", x = "Intra (MAD)", y = "Inter (MAD)") + theme_minimal() + 
  theme(legend.position = 'bottom')

ggExtra::ggMarginal(losic.mad, type = "boxplot", groupColour = TRUE, groupFill = TRUE)
grid.newpage()


renji.mad <- ggplot(renji.ith.score, aes(x = renji.intra.mad, y = renji.inter.mad, color = group)) + 
  geom_point(size = 1) + 
  scale_color_manual(values = c("Q1 Gene" = "orange", "Other" = "black"), name = "Group") + 
  labs(title = "Yang et al.", x = "Intra (MAD)", y = "Inter (MAD)") + theme_minimal() + 
  theme(legend.position = 'bottom')

ggExtra::ggMarginal(renji.mad, type = "boxplot", groupColour = TRUE, groupFill = TRUE)
grid.newpage()


shen.mad <- ggplot(shen.ith.score, aes(x = shen.intra.mad, y = shen.inter.mad, color = group)) + 
  geom_point(size = 1) + 
  scale_color_manual(values = c("Q1 Gene" = "orange", "Other" = "black"), name = "Group") + 
  labs(title = "Shen et al.", x = "Intra (MAD)", y = "Inter (MAD)") + theme_minimal() + 
  theme(legend.position = 'bottom')

ggExtra::ggMarginal(shen.mad, type = "boxplot", groupColour = TRUE, groupFill = TRUE)
grid.newpage()


shi.mad <- ggplot(shi.ith.score, aes(x = shi.intra.mad, y = shi.inter.mad, color = group)) + 
  geom_point(size = 1) + 
  scale_color_manual(values = c("Q1 Gene" = "orange", "Other" = "black"), name = "Group") + 
  labs(title = "Shi et al.", x = "Intra (MAD)", y = "Inter (MAD)") + theme_minimal() + 
  theme(legend.position = 'bottom')

ggExtra::ggMarginal(shi.mad, type = "boxplot", groupColour = TRUE, groupFill = TRUE)

dev.off()
}

# CV
{
pdf("/result/section1/top_right_cv_score.pdf")
losic.cv <- ggplot(losic.ith.score, aes(x = losic.intra.cv, y = losic.inter.cv, color = group)) + 
  geom_point(size = 1) + 
  scale_color_manual(values = c("Q1 Gene" = "orange", "Other" = "black"), name = "Group") + 
  labs(title = "Losic et al.", x = "Intra (CV)", y = "Inter (CV)") + theme_minimal() + 
  theme(legend.position = 'bottom')

ggExtra::ggMarginal(losic.cv, type = "boxplot", groupColour = TRUE, groupFill = TRUE)
grid.newpage()

renji.cv <- ggplot(renji.ith.score, aes(x = renji.intra.cv, y = renji.inter.cv, color = group)) + 
  geom_point(size = 1) + 
  scale_color_manual(values = c("Q1 Gene" = "orange", "Other" = "black"), name = "Group") + 
  labs(title = "Yang et al.", x = "Intra (CV)", y = "Inter (CV)") + theme_minimal() + 
  theme(legend.position = 'bottom')

ggExtra::ggMarginal(renji.cv, type = "boxplot", groupColour = TRUE, groupFill = TRUE)
grid.newpage()


shen.cv <- ggplot(shen.ith.score, aes(x = shen.intra.cv, y = shen.inter.cv, color = group)) + 
  geom_point(size = 1) + 
  scale_color_manual(values = c("Q1 Gene" = "orange", "Other" = "black"), name = "Group") + 
  labs(title = "Shen et al.", x = "Intra (CV)", y = "Inter (CV)") + theme_minimal() + 
  theme(legend.position = 'bottom')

ggExtra::ggMarginal(shen.cv, type = "boxplot", groupColour = TRUE, groupFill = TRUE)
grid.newpage()


shi.cv <- ggplot(shi.ith.score, aes(x = shi.intra.cv, y = shi.inter.cv, color = group)) + 
  geom_point(size = 1) + 
  scale_color_manual(values = c("Q1 Gene" = "orange", "Other" = "black"), name = "Group") + 
  labs(title = "Shi et al.", x = "Intra (CV)", y = "Inter (CV)") + theme_minimal() + 
  theme(legend.position = 'bottom')

ggExtra::ggMarginal(shi.cv, type = "boxplot", groupColour = TRUE, groupFill = TRUE)

dev.off()
}


# > wilcox.test(losic.intra.mad ~ group, data = losic.ith.score)$p.value
# [1] 2.418536e-185
# > wilcox.test(losic.inter.mad ~ group, data = losic.ith.score)$p.value
# [1] 8.946868e-207
# > wilcox.test(losic.intra.cv ~ group, data = losic.ith.score)$p.value
# [1] 2.085398e-79
# > wilcox.test(losic.inter.cv ~ group, data = losic.ith.score)$p.value
# [1] 1.782257e-108
# > wilcox.test(renji.intra.mad ~ group, data = renji.ith.score)$p.value
# [1] 1.872449e-206
# > wilcox.test(renji.inter.mad ~ group, data = renji.ith.score)$p.value
# [1] 1.205884e-212
# > wilcox.test(renji.intra.cv ~ group, data = renji.ith.score)$p.value
# [1] 3.994068e-83
# > wilcox.test(renji.inter.cv ~ group, data = renji.ith.score)$p.value
# [1] 1.041102e-117
# > wilcox.test(shen.intra.mad ~ group, data = shen.ith.score)$p.value
# [1] 6.867077e-179
# > wilcox.test(shen.inter.mad ~ group, data = shen.ith.score)$p.value
# [1] 7.050443e-216
# > wilcox.test(shen.intra.cv ~ group, data = shen.ith.score)$p.value
# [1] 4.462389e-163
# > wilcox.test(shen.inter.cv ~ group, data = shen.ith.score)$p.value
# [1] 3.735043e-183
# > wilcox.test(shi.intra.mad ~ group, data = shi.ith.score)$p.value
# [1] 7.951029e-190
# > wilcox.test(shi.inter.mad ~ group, data = shi.ith.score)$p.value
# [1] 2.406024e-187
# > wilcox.test(shi.intra.cv ~ group, data = shi.ith.score)$p.value
# [1] 4.802988e-172
# > wilcox.test(shi.inter.cv ~ group, data = shi.ith.score)$p.value
# [1] 2.724363e-202

