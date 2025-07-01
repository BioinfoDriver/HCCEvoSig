

# Function to fit meta-analytic equal-, fixed-, and random-effects models 
# https://stats.stackexchange.com/questions/343316/hazard-ratio-meta-analysis
# https://www.bmj.com/content/343/bmj.d2304
MetaforHrMetaAnalysis <- function(dat, method=c('EE', 'FE', 'REML'), ifplot, fname){
  method <- match.arg(method)
  
  library('metafor')
  library('dplyr')
  
  colnames(dat) <- c('study', 'hr', 'ci.lb', 'ci.ub', 'pval')
  dat <- dat %>% mutate(yi = log(hr), sei = (log(ci.ub)-log(ci.lb))/(2*1.96))
  
  # whether an equal-, a fixed or random-effects model should be fitted.
  res <- rma(yi=yi, sei=sei, data=dat, method=method, slab=study)
  
  meta.sum <- exp(c(res$b[1, 1], res$ci.lb, res$ci.ub))
  meta.sum <- c(meta.sum, res$pval)
  names(meta.sum) <- c('hr', 'ci.lb', 'ci.ub', 'pval')
  
  # plot
  if(ifplot){
    pdf(fname)
    forest(x=res, annotate=TRUE, header=c('Author(s) and Year', 'HR [95% CI]'), refline=1, 
           xlab='Hazard ratio', mlab='Overall', ilab=dat$pval, ilab.xpos=7, ilab.pos=2,
           colout='#3182bd', col='#a50f15', transf=exp)
    dev.off()
  }
  
  return(meta.sum)
}


RmetaHrMetaAnalysis <- function(dat, method=c("fixed", "random"), ifplot, fname){
  method <- match.arg(method)
  
  library('rmeta')
  library('dplyr')
  
  colnames(dat) <- c('study', 'hr', 'ci.lb', 'ci.ub', 'pval')
  dat <- dat %>% mutate(yi = log(hr), sei = (log(ci.ub)-log(ci.lb))/(2*1.96))
  
  # whether a fixed- or random-effects model should be fitted.
  res <- meta.summaries(d = yi, se = sei, method = method, logscale=TRUE, names=study, data=dat)
  
  meta.sum <- exp(c(res$summary, res$summary - 1.96 * res$se.summary, res$summary + 1.96 * res$se.summary))
  meta.sum <- c(meta.sum, res$test[2])
  names(meta.sum) <- c('hr', 'ci.lb', 'ci.ub', 'pval')
  
  # plot
  if(ifplot){
    pdf(fname)
    metaplot(mn=dat$yi, se=dat$sei, labels=dat$study, xlab='Hazard Ratio', 
             summn = res$summary, sumse = res$se.summary, sumnn= 1/res$se.summary^2, xlim=c(-1, 3), summlabel="Overall",
             zero=0, colors=meta.colors(box="#3182bd",lines="#a50f15", zero="red", summary="black",text="black"), xaxt='n')
    axis(1, at=log(c(0.5,1,2,4,8,16)), labels=c(0.5,1,2,4,8,16))
    dev.off()
  }
  
  return(meta.sum)
}

exam.data <- data.frame(dataset=c('TCGA-LIHC', 'ICGC-LIRI-JP', 'CHCC-HBV', 'Mongolian-HCC', 'FULCI-HCC', 'NCI-HCC'), 
                        HR=c(3.40, 5.06, 4.20, 2.88, 2.84, 1.76), 
                        lower=c(2.27, 2.20, 2.26, 1.18, 1.80, 1.08), upper=c(5.09, 11.64, 7.81, 7.02, 4.49, 2.87), 
                        pvalue=c('<0.0001', '0.0001', '<0.0001', '0.0197', '<0.0001', '0.0241'))

RmetaHrMetaAnalysis(exam.data, 'fixed', TRUE, '/result/Section4/metaHRmetafor-1.pdf')
MetaforHrMetaAnalysis(exam.data, 'FE', TRUE, '/result/Section4/metaHRrmeta-1.pdf')

# hr        ci.lb        ci.ub         pval 
# 2.977960e+00 2.386653e+00 3.715766e+00 4.345931e-22