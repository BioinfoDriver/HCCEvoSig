# load data
tcga.lihc.cli.char <- readRDS(file='/data/tcga.lihc.cli.char.rds')
tcga.lihc.risk.score <- readRDS(file='/data/tcga_risk_score.rds')

# data prepare
sig.score <- tcga.lihc.risk.score[, c('patient_id', 'os', 'os_time', 'risk.score', 'risk.categ')]
sig.score$risk.categ <- factor(sig.score$risk.categ, levels=c('low risk', 'high risk'))

cli.char <- tcga.lihc.cli.char[, c('bcr_patient_barcode', 'age_category', 'gender_category', 'race_category', 
                                   'family_cancer_history_category', 'BMI_category',
                                   'pathologic_stage_category', 'histologic_grade_category', 'residual_tumor_category', 'ishak_score_category', 
                                   'child_pugh_grade_category', 'alpha_fetoprotein',
                                   'vascular_tumor_invasion_category','inflammation_extent_category')]

cli.sig.char <- merge(sig.score, cli.char, by.x='patient_id', by.y='bcr_patient_barcode')

# 单因素cox分析
UnivariateCox <- function(cli.data, covariates)
{
  library('survival')
  #STEP1:构建单因素分析的对象
  univ_formulas <- sapply(covariates,
                          function(x) as.formula(paste('Surv(os_time, os)~', x)));
  
  #STEP2:单因素Cox分析
  univ_models <- lapply(univ_formulas, function(x){coxph(x, data = cli.data)});
  
  #STEP3:提取有用信息
  univ_results <- lapply(univ_models, function(x)
  {                             
    tmp <-summary(x);
    
    #提取p值，保留两位有效数字
    p.value <- round(tmp$coefficients[ ,5], digits = 4);
    p.value[which(p.value < 0.0001)] <- "<0.0001";
    
    #提取beta值，这里的coefficients为矩阵，但是只有一行，所以可以这样取值
    #beta <- round(tmp$coefficients[ ,1], digits = 4);
    
    #提取风险比
    HR <- round(tmp$coefficients[ ,2], digits = 4);
    
    #提取95%置信区间上下界
    HR.confint.lower <- round(tmp$conf.int[,"lower .95"], 4);
    HR.confint.upper <- round(tmp$conf.int[,"upper .95"], 4);    
    
    #合并风险比HR和置信区间为一个内容
    HR <- paste0(HR, " (", HR.confint.lower, "-", HR.confint.upper, ")");
    
    variate <- rownames(tmp$coefficients);
    
    #将所有值合并在一个矩阵中
    all.data <- as.data.frame(cbind(variate, HR, p.value));
  }
  )
  univ_results <- do.call(rbind, univ_results)
  return(univ_results)
}

uni.cox.res <- UnivariateCox(cli.data = cli.sig.char, covariates=colnames(cli.sig.char)[5:18])
# sapply(cli.sig.char[, 5:18], function(x) sum(!is.na(x)))

# 多因素cox分析
cli.sig.char <- cli.sig.char[, c('patient_id', 'os', 'os_time', 'risk.categ', 'age_category', 'gender_category', 
                                 'pathologic_stage_category', 'histologic_grade_category', 'ishak_score_category', 'alpha_fetoprotein')]
cli.sig.char <- cli.sig.char[apply(cli.sig.char, 1, function(x) !any(is.na(x))), ]


source('/code/Rscript/Cox.function.R')
uni.mul.cox.res <- Cox.function(time=cli.sig.char$os_time, event=cli.sig.char$os, 
                                clinical.data=cli.sig.char, clinical.variate = c(4:10))

