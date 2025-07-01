# load data
clinical.data <- readRDS(file='/data/tcga_sig_risk_score.rds')

# Time-dependent ROC curve analysis
TimeDependentROC <- function(clinical.data, sur.time, sur.status, models, t.points){
  
  library(timeROC)
  library(survival)
  
  time.roc.auc <- lapply(models, function(x){
    
    time.roc <- timeROC(T = sur.time, delta = sur.status, marker = clinical.data[, x], cause=1,
                        times = t.points, ROC = TRUE, iid = TRUE)
    
    # Area under the curve
    auc <- time.roc$AUC
    # Confidence intervals
    auc.ci95 <- confint(time.roc, level = 0.95)$CI_AUC/100
    
    roc.auc <- data.frame(cbind(auc, auc.ci95))
    roc.auc$auc_ci <- paste0(sprintf('%0.2f', roc.auc$auc), 
                             '(', sprintf('%0.2f', roc.auc$X2.5.), '-', sprintf('%0.2f', roc.auc$X97.5.), ')')
    
    roc.auc$model <- x
    return(roc.auc)
  })
  
  return(time.roc.auc)
}

risk.sig <- c('risk.score', 'PMID_33251144', 'PMID_32198063', 'PMID_35123387', 'PMID_31335995', 'PMID_33033585',
              'PMID_33828988', 'PMID_34676211', 'PMID_32903581', 'PMID_34900672', 'PMID_33089373', 'PMID_35311113', 
              'PMID_34975331', 'PMID_25666192', 'PMID_22105560', 'PMID_23800896')


# risk.sig <- c('risk.score', grep('PMID_', colnames(clinical.data), value = T))
  
sig.td.roc <- TimeDependentROC(clinical.data, clinical.data$os_time, clinical.data$os, 
                               risk.sig, t.points = c(1*365, 3*365, 5*365))

sig.td.roc <- do.call(rbind, sig.td.roc)
sig.td.roc[seq(from = 1, to = 16*3, by = 3), ]
sig.td.roc[seq(from = 2, to = 16*3, by = 3), ]
sig.td.roc[seq(from = 3, to = 16*3, by = 3), ]