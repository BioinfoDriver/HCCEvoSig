
QuadrantPlot <- function(ith.score, file.name){
  
  ##  input: rna-ith scores
  colnames(ith.score) <- c('intra.score', 'inter.score')
  
  # within.ith <- ifelse(ith.score$intra.score > mean(ith.score$intra.score), 'top', 'bottom')
  # between.ith <- ifelse(ith.score$inter.score > mean(ith.score$inter.score), 'right', 'left')
  within.ith <- ifelse(ith.score$intra.score > quantile(ith.score$intra.score, 0.75), 'top', 'bottom')
  between.ith <- ifelse(ith.score$inter.score > quantile(ith.score$inter.score, 0.75), 'right', 'left')
  ith.score$quadrant <- paste(within.ith, between.ith, sep='_')
  
  ## Make new quadrant type figure - 2 density plots, one with four to scale circles
  a <- table(factor(ith.score[,'quadrant'],levels=c("top_left", "bottom_left", "top_right", "bottom_right")))
  b <- a/min(a)
  
  pdf(paste0(file.name, '.pdf'))
  # scale circles
  plot(c(100,100,400,400),c(400,100,400,100), pch=16, cex=b*2, ylim=c(0,500),xlim=c(0,500), 
       col=c("firebrick1", "darkorchid2", "gold1", "turquoise2"),axes=F,ylab='',xlab='')
  
  abline(h=250,v=250,col=1,lty=2,lwd=4)
  text(c(100,100,400,400),c(450,200,450,200), paste(a,'genes'))
  box()
  
  # density plots
  hist(ith.score[,'inter.score'], breaks=100, prob=TRUE,col='#6baed699',axes=F,ylab='', xlab='',main='Between tumour, LIHC')
  lines(density(ith.score[,'inter.score']), lwd = 3,col = '#08519c')
  abline(v=quantile(ith.score[,'inter.score'], 0.75), col=1,lwd=4, lty=2)
  
  
  hist(ith.score[,'intra.score'], breaks=100, prob=TRUE,col='#fd8d3c99',axes=F,ylab='', xlab='',main='Within tumour, LIHC')
  lines(density(ith.score[,'intra.score']), lwd = 3,col = '#a63603')
  abline(v=quantile(ith.score[,'intra.score'], 0.75), col=1,lwd=4, lty=2)
  
  dev.off()
}


# load data
load(file='/data/intra.inter.ith.score.RData')


setwd('/result/section1/')
QuadrantPlot(losic.ith.score[, c('losic.intra.sd', 'losic.inter.sd')], 'losic_quadrant_plot_sd')
QuadrantPlot(renji.ith.score[, c('renji.intra.sd', 'renji.inter.sd')], 'renji_quadrant_plot_sd')
QuadrantPlot(shen.ith.score[, c('shen.intra.sd', 'shen.inter.sd')], 'shen_quadrant_plot_sd')
QuadrantPlot(shi.ith.score[, c('shi.intra.sd', 'shi.inter.sd')], 'shi_quadrant_plot_sd')




ITHQuadrant <- function(ith.score){
  
  ##  input: rna-ith scores
  colnames(ith.score) <- c('intra.score', 'inter.score')
  
  # within.ith <- ifelse(ith.score$intra.score > mean(ith.score$intra.score), 'top', 'bottom')
  # between.ith <- ifelse(ith.score$inter.score > mean(ith.score$inter.score), 'right', 'left')
  within.ith <- ifelse(ith.score$intra.score > quantile(ith.score$intra.score, 0.75), 'top', 'bottom')
  between.ith <- ifelse(ith.score$inter.score > quantile(ith.score$inter.score, 0.75), 'right', 'left')
  ith.score$quadrant <- paste(within.ith, between.ith, sep='_')
  
  return(ith.score)
}


losic.ith.quadrant <- ITHQuadrant(losic.ith.score[, c('losic.intra.sd', 'losic.inter.sd')])
renji.ith.quadrant <- ITHQuadrant(renji.ith.score[, c('renji.intra.sd', 'renji.inter.sd')])
shen.ith.quadrant <- ITHQuadrant(shen.ith.score[, c('shen.intra.sd', 'shen.inter.sd')])
shi.ith.quadrant <- ITHQuadrant(shi.ith.score[, c('shi.intra.sd', 'shi.inter.sd')])


# intersect(rownames(subset(losic.ith.quadrant, quadrant == 'top_right')), rownames(subset(renji.ith.quadrant, quadrant == 'top_right')))
# intersect(rownames(subset(losic.ith.quadrant, quadrant == 'top_right')), rownames(subset(shen.ith.quadrant, quadrant == 'top_right')))
# intersect(rownames(subset(losic.ith.quadrant, quadrant == 'top_right')), rownames(subset(shi.ith.quadrant, quadrant == 'top_right')))
# 
# 
# Reduce(intersect, list(rownames(subset(losic.ith.quadrant, quadrant == 'top_right')),
#                        rownames(subset(renji.ith.quadrant, quadrant == 'top_right')),
#                        rownames(subset(shen.ith.quadrant, quadrant == 'top_right')),
#                        rownames(subset(shi.ith.quadrant, quadrant == 'top_right'))))


save(losic.ith.quadrant, renji.ith.quadrant, shen.ith.quadrant, shi.ith.quadrant, file='/data/intra.inter.ith.quadrant.RData')


