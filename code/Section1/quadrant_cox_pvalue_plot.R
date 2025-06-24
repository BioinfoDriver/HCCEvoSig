library(dplyr)
library(tibble)

##  core function for quadrant box-plots
boxplot_DB <- function(data, x, y, fill_var, fill_scale=c("firebrick1", "darkorchid1", "gold1", "turquoise1"), title="", subtitle="", x_axis="x", y_axis="y", test="wilcox", aspect_ratio=1, add_violin=T, exclude_NA=T, log_axis=T, zero_axis=T) {
  
  # tidy input data
  data <- data %>% as.data.frame()
  if(exclude_NA) {data <- data %>% tidyr::drop_na()}
  
  # extract y-variable summary stats
  if(add_violin){ y_max <- ceiling(max(data[,y, drop=T], na.rm = T)) + 1 }
  if(!add_violin){ y_max <- ceiling(fivenum(data[,y, drop=T], na.rm = T)[4]) + 1 }
  
  # perform statistical test
  if(test=="wilcox") {
    
    ##  pairwise wilcox test
    # calculate p-values for all comparisons
    pv <- pairwise.wilcox.test(x = data[, y, drop=T], g = data[, x, drop=T], p.adjust.method="none", paired=F)$p.value
    # tidy: create "quad_1" and "quad_2" cols to show comparison for each test
    
    
    pv <- pv %>% as.data.frame() %>% tibble::rownames_to_column("quad_1")
    pv <- reshape2::melt(pv, id.vars="quad_1")
    
    colnames(pv)[2:3] <- c("quad_2", "p_val")
    # QC: output to console
    print(pv)
    # filter significant values, and label with asterisks
    pv <- pv[which(pv$p_val < 0.05), ]
    pv$symbol <- ifelse(pv$p_val > 0.01, "*", ifelse(pv$p_val > 0.001, "**", "***"))
    # fix class of "quad_2" col
    pv$quad_2 <- pv$quad_2 %>% as.character()
    
  }
  
  # plot
  gg <- ggplot(data, aes_string(x, y, fill=fill_var)) 
  gg <- gg + scale_fill_manual(values=fill_scale)
  # add violin
  if(add_violin) {gg <- gg + geom_violin(alpha=0.25, scale="width", width=0.8) }
  # add boxplot
  gg <- gg + geom_boxplot(width=0.15, outlier.shape = NA)
  # add significance bars
  if(test=="wilcox") {
    if(nrow(pv)>0) {
      tmp <- as.data.frame(t(pv[, 1:2]))
      tmp[] <- sapply(tmp, as.character)
      gg <- gg + geom_signif(comparisons = as.list(tmp), annotation=pv$symbol)
    }
  }
  # add themes: black border, origin at 0, inward y tickmarks
  gg <- gg + theme_classic() + theme(legend.position = "none") + theme(axis.line = element_blank(), panel.border = element_rect(colour = "black", fill=NA, size=1), aspect.ratio = aspect_ratio, axis.text.y = element_text(size = 8, margin = unit(c(t = 0, r = 2.5, b = 0, l = 0), "mm")), axis.text.x = element_text(size = 8, margin = unit(c(t = 2.5, r = 0, b = 0, l = 0), "mm")),  axis.ticks.length = unit(-1.4, "mm")) 
  
  if(log_axis) {
    sig_values <- c(1, 0.05, 0.01, 10^(-seq(3, y_max)))
    gg <- gg + scale_y_continuous(expand = c(0,0), limits = c(0, y_max), breaks = sapply(X=sig_values, FUN=function(x) {-log10(x)}), labels=sig_values)
  }
  if(zero_axis) {
    gg <- gg + scale_y_continuous(expand = c(0,0), limits = c(0, y_max))
  }
  
  # add plot labels
  gg <- gg + ggtitle(label=title, subtitle=subtitle)
  gg <- gg + xlab(x_axis) + ylab(y_axis) + theme(plot.title=element_text(hjust=0.5, face="bold"))
  # output
  gg %>% return()
  
}

## wrapper function for gene-wise survival association ~ RNA heterogeneity quadrant as box-plots
quadrant_hr_boxplot_DB <- function(prog, gene_hazard_ratios, histology="LIHC", cohort_name="") {
  ##  UPDATE: 17 SEPTEMBER 2018, gene_hazard_ratios = output from rapid Cox UVA, unique gene names
  
  df <- prog
  
  ##  PLOT: p-value ~ signature, coloured by heterogeneity quadrant
  
  ##  prep HR info for joining
  tmp <- dplyr::select(gene_hazard_ratios, Gene, P)
  colnames(tmp)[which(colnames(tmp)=="P")] <- "pval"
  
  ##  prep prog: filter histology & unique gene names
  df <- prog %>% dplyr::filter(Histology==histology) %>% dplyr::select(Gene, quadrant) %>% dplyr::distinct()
  
  ##  inner join: drop NA values (quadrant, pval)
  df <- dplyr::inner_join(x=df, y=tmp, by="Gene")
  
  ##  no. signatures
  n_sigs <- nrow(unique(prog[prog$Histology==histology, "signature"]))
  
  ##prep data for plotting
  # log-transform pvals
  df$log_p <- -log10(df$pval)
  # rename quadrants
  df$q <- c("top_right"="Q1", "top_left"="Q2", "bottom_left"="Q3", "bottom_right"="Q4")[df$quadrant]
  
  
  ##  box-plot
  # boxplot
  gg <- boxplot_DB(
    data = df, 
    x = "q", 
    y = "log_p", 
    fill_var = "q", 
    title = "", 
    subtitle = paste(n_sigs, histology, "prognostic signatures"), 
    x_axis = "", 
    y_axis = paste0("Gene-Wise Prognostic Value in\n", cohort_name, " cohort (", gene_hazard_ratios$n_patients, " patients)\n-log10(Cox UVA P-value)"), zero_axis = F
  ) 
  gg + geom_hline(yintercept = -log10(0.05), lty="dashed") + geom_hline(yintercept = -log10(0.01), lty="dotted")
  
}


quadrant_cox_pvalue_plot <- function(quadrant.gene, cox.pvalue, plot.file.name, cohort_name = ''){
  # input
  prog_histology <- quadrant.gene
  prog_histology$Histology <- "LIHC"
  
  #  run analysis
  library(ggplot2)
  gg <- quadrant_hr_boxplot_DB(prog_histology, cox.pvalue, "LIHC", cohort_name)
  
  ##  re-plot with greater y-axis upper limit
  y_max <- 12
  sig_values <- c(1, 0.05, 0.01, 10^(-seq(3, y_max)))
  gg <- gg + scale_y_continuous(expand = c(0,0), limits = c(0, y_max), 
                                breaks = sapply(X=sig_values, FUN=function(x) {-log10(x)}), labels=sig_values)
  # gg %>% print()
  
  ggsave(filename=plot.file.name, plot=gg, path='D:/LIHCSurvival/result/section1/')
}


quadrant_precog_zscore_plot <- function(quadrant.gene, pancancer=TRUE, plot.file.name){
  # prep pan-cancer z-score values from PRECOG 
  if(pancancer){
    data <- dplyr::inner_join(x=precog_z[, c("Gene", "Unweighted_meta.Z_of_all_cancers")], y=quadrant.gene[, c("Gene", "quadrant")], by="Gene")
  }else{
    data <- dplyr::inner_join(x=precog_z[, c("Gene", "Liver_cancer")], y=quadrant.gene[, c("Gene", "quadrant")], by="Gene")
  }
  
  colnames(data)[2] <- "pancan_z"
  # join NSCLC quadrants
  data$quadrant <- data$quadrant %>% gsub(pattern="top_right", replacement="Q1") %>% gsub(pattern="top_left", replacement="Q2") %>% gsub(pattern="bottom_left", replacement="Q3") %>% gsub(pattern="bottom_right", replacement="Q4")
  data$quadrant <- factor(data$quadrant, levels = c("Q1", "Q2", "Q3", "Q4"))
  
  # box-plot
  gg <- boxplot_DB(data = data, x = "quadrant", y = "pancan_z", fill_var = "quadrant", x_axis = "Quadrant", y_axis = "Pan-cancer prognostic value\n(z score)")
  
  ggsave(filename=plot.file.name, plot=gg, path='D:/LIHCSurvival/result/section1/')
}



#################################
load('/data/univ_cox_pvalue.RData')
evo_genes <- readRDS('/data/evo_genes.rds')

tcga.cox.pvalue <- data.frame(Gene=names(tcga.univ.cox.pvalue), P=tcga.univ.cox.pvalue, n_patients=323)
icgc.cox.pvalue <- data.frame(Gene=names(icgc.univ.cox.pvalue), P=icgc.univ.cox.pvalue, n_patients=203)
chcc.cox.pvalue <- data.frame(Gene=names(chcc.univ.cox.pvalue), P=chcc.univ.cox.pvalue, n_patients=159)
fulci.cox.pvalue <- data.frame(Gene=names(fulci.univ.cox.pvalue), P=fulci.univ.cox.pvalue, n_patients=221)
nci.cox.pvalue <- data.frame(Gene=names(nci.univ.cox.pvalue), P=nci.univ.cox.pvalue, n_patients=112)
mongolian.cox.pvalue <- data.frame(Gene=names(mongolian.univ.cox.pvalue), P=mongolian.univ.cox.pvalue, n_patients=70)


tcga.quadrant <- data.frame(Gene = names(tcga.univ.cox.pvalue), 
                            quadrant = ifelse(names(tcga.univ.cox.pvalue) %in% evo_genes, 'top_right', 'top_left'), signature = 'COX UVA')
icgc.quadrant <- data.frame(Gene = names(icgc.univ.cox.pvalue), 
                            quadrant = ifelse(names(icgc.univ.cox.pvalue) %in% evo_genes, 'top_right', 'top_left'), signature = 'COX UVA')
chcc.quadrant <- data.frame(Gene = names(chcc.univ.cox.pvalue), 
                            quadrant = ifelse(names(chcc.univ.cox.pvalue) %in% evo_genes, 'top_right', 'top_left'), signature = 'COX UVA')
fulci.quadrant <- data.frame(Gene = names(fulci.univ.cox.pvalue), 
                            quadrant = ifelse(names(fulci.univ.cox.pvalue) %in% evo_genes, 'top_right', 'top_left'), signature = 'COX UVA')
nci.quadrant <- data.frame(Gene = names(nci.univ.cox.pvalue), 
                            quadrant = ifelse(names(nci.univ.cox.pvalue) %in% evo_genes, 'top_right', 'top_left'), signature = 'COX UVA')
mongolian.quadrant <- data.frame(Gene = names(mongolian.univ.cox.pvalue), 
                            quadrant = ifelse(names(mongolian.univ.cox.pvalue) %in% evo_genes, 'top_right', 'top_left'), signature = 'COX UVA')


quadrant_cox_pvalue_plot(tcga.quadrant, tcga.cox.pvalue, 'evo_genes_tcga_coxp.pdf', cohort_name = 'TCGA') # 1.157704e-07
quadrant_cox_pvalue_plot(icgc.quadrant, icgc.cox.pvalue, 'evo_genes_icgc_coxp.pdf', cohort_name = 'ICGC') # 4.448762e-05
quadrant_cox_pvalue_plot(chcc.quadrant, chcc.cox.pvalue, 'evo_genes_chcc_coxp.pdf', cohort_name = 'CHCC') # 4.937998e-24
quadrant_cox_pvalue_plot(fulci.quadrant, fulci.cox.pvalue, 'evo_genes_fulci_coxp.pdf', cohort_name = 'FULCI') # 1.857932e-17
quadrant_cox_pvalue_plot(nci.quadrant, nci.cox.pvalue, 'evo_genes_nci_coxp.pdf', cohort_name = 'NCI') # 8.11698e-08
quadrant_cox_pvalue_plot(mongolian.quadrant, mongolian.cox.pvalue, 'evo_genes_mongolian_coxp.pdf', cohort_name = 'Mongolian') # 2.790114e-11


# Fisher test
fisher.test(table(merge(tcga.quadrant, tcga.cox.pvalue, by = 'Gene') %>% 
                    mutate(sig = ifelse(P < 0.05, 'sig', 'unsig'), quadrant = factor(quadrant, level = c('top_right', 'top_left'))) %>% 
                    dplyr::select(sig, quadrant))) # OR = 1.727904, p-value = 7.213e-08

fisher.test(table(merge(icgc.quadrant, icgc.cox.pvalue, by = 'Gene') %>% 
                    mutate(sig = ifelse(P < 0.05, 'sig', 'unsig'), quadrant = factor(quadrant, level = c('top_right', 'top_left'))) %>% 
                    dplyr::select(sig, quadrant))) # OR = 1.314747, p-value = 0.007065

fisher.test(table(merge(chcc.quadrant, chcc.cox.pvalue, by = 'Gene') %>% 
                    mutate(sig = ifelse(P < 0.05, 'sig', 'unsig'), quadrant = factor(quadrant, level = c('top_right', 'top_left'))) %>% 
                    dplyr::select(sig, quadrant))) # OR = 2.613234, p-value = 1.382206e-21

fisher.test(table(merge(fulci.quadrant, fulci.cox.pvalue, by = 'Gene') %>% 
                    mutate(sig = ifelse(P < 0.05, 'sig', 'unsig'), quadrant = factor(quadrant, level = c('top_right', 'top_left'))) %>% 
                    dplyr::select(sig, quadrant))) # OR = 2.759802, p-value = 3.778794e-18

fisher.test(table(merge(nci.quadrant, nci.cox.pvalue, by = 'Gene') %>% 
                    mutate(sig = ifelse(P < 0.05, 'sig', 'unsig'), quadrant = factor(quadrant, level = c('top_right', 'top_left'))) %>% 
                    dplyr::select(sig, quadrant))) # OR = 2.111864, p-value = 6.411e-11

fisher.test(table(merge(mongolian.quadrant, mongolian.cox.pvalue, by = 'Gene') %>% 
                    mutate(sig = ifelse(P < 0.05, 'sig', 'unsig'), quadrant = factor(quadrant, level = c('top_right', 'top_left'))) %>% 
                    dplyr::select(sig, quadrant))) # OR = 1.642333, p-value = 2.782e-06



p1 <- ggplot(data = merge(tcga.quadrant, tcga.cox.pvalue, by = 'Gene') %>% 
  mutate(sig = ifelse(P < 0.05, 'sig', 'unsig'), quadrant = factor(quadrant, level = c('top_right', 'top_left'))),
  aes(x = quadrant, fill = sig)) + geom_bar(position = "fill", stat = "count") +
  scale_fill_manual(values = c("unsig" = "blue", "sig" = "red"), labels = c("P < 0.05", "P > 0.05")) + 
  scale_x_discrete(labels = c("top_right" = "Q1", "top_left" = "Other")) + 
  theme_minimal() + labs(x = NULL, y = "Percentage", title = "TCGA", fill = 'Cox Pvalue') + 
  theme(legend.position = "bottom")


p2 <- ggplot(data = merge(icgc.quadrant, icgc.cox.pvalue, by = 'Gene') %>% 
         mutate(sig = ifelse(P < 0.05, 'sig', 'unsig'), quadrant = factor(quadrant, level = c('top_right', 'top_left'))),
       aes(x = quadrant, fill = sig)) + geom_bar(position = "fill", stat = "count") +
  scale_fill_manual(values = c("unsig" = "blue", "sig" = "red"), labels = c("P < 0.05", "P > 0.05")) + 
  scale_x_discrete(labels = c("top_right" = "Q1", "top_left" = "Other")) + 
  theme_minimal() + labs(x = NULL, y = "Percentage", title = "ICGC", fill = 'Cox Pvalue')+
  theme(legend.position = "bottom")


p3 <- ggplot(data = merge(chcc.quadrant, chcc.cox.pvalue, by = 'Gene') %>% 
         mutate(sig = ifelse(P < 0.05, 'sig', 'unsig'), quadrant = factor(quadrant, level = c('top_right', 'top_left'))),
       aes(x = quadrant, fill = sig)) + geom_bar(position = "fill", stat = "count") +
  scale_fill_manual(values = c("unsig" = "blue", "sig" = "red"), labels = c("P < 0.05", "P > 0.05")) + 
  scale_x_discrete(labels = c("top_right" = "Q1", "top_left" = "Other")) + 
  theme_minimal() + labs(x = NULL, y = "Percentage", title = "CHCC", fill = 'Cox Pvalue') + 
  theme(legend.position = "bottom")


p4 <- ggplot(data = merge(nci.quadrant, nci.cox.pvalue, by = 'Gene') %>% 
         mutate(sig = ifelse(P < 0.05, 'sig', 'unsig'), quadrant = factor(quadrant, level = c('top_right', 'top_left'))),
       aes(x = quadrant, fill = sig)) + geom_bar(position = "fill", stat = "count") +
  scale_fill_manual(values = c("unsig" = "blue", "sig" = "red"), labels = c("P < 0.05", "P > 0.05")) + 
  scale_x_discrete(labels = c("top_right" = "Q1", "top_left" = "Other")) + 
  theme_minimal() + labs(x = NULL, y = "Percentage", title = "NCI", fill = 'Cox Pvalue') + 
  theme(legend.position = "bottom")


p5 <- ggplot(data = merge(fulci.quadrant, fulci.cox.pvalue, by = 'Gene') %>% 
         mutate(sig = ifelse(P < 0.05, 'sig', 'unsig'), quadrant = factor(quadrant, level = c('top_right', 'top_left'))),
       aes(x = quadrant, fill = sig)) + geom_bar(position = "fill", stat = "count") +
  scale_fill_manual(values = c("unsig" = "blue", "sig" = "red"), labels = c("P < 0.05", "P > 0.05")) + 
  scale_x_discrete(labels = c("top_right" = "Q1", "top_left" = "Other")) + 
  theme_minimal() + labs(x = NULL, y = "Percentage", title = "FULCI", fill = 'Cox Pvalue') + 
  theme(legend.position = "bottom")


p6 <- ggplot(data = merge(mongolian.quadrant, mongolian.cox.pvalue, by = 'Gene') %>% 
         mutate(sig = ifelse(P < 0.05, 'sig', 'unsig'), quadrant = factor(quadrant, level = c('top_right', 'top_left'))),
       aes(x = quadrant, fill = sig)) + geom_bar(position = "fill", stat = "count") +
  scale_fill_manual(values = c("unsig" = "blue", "sig" = "red"), labels = c("P < 0.05", "P > 0.05")) + 
  scale_x_discrete(labels = c("top_right" = "Q1", "top_left" = "Other")) + 
  theme_minimal() + labs(x = NULL, y = "Percentage", title = "Mongolian", fill = 'Cox Pvalue') + 
  theme(legend.position = "bottom")


com.plot <- ggarrange(p1, p2, p3, p4, p5, p6, labels = c("A", "B", "C", "D", "E", "F"), ncol = 3, nrow = 2)

ggsave(plot = com.plot, filename = '/result/section1/evo_genes_cox_sig_enrich.pdf')




#################################
load('D:/LIHCSurvival/data/univ_cox_pvalue.RData')
load(file='D:/LIHCSurvival/data/intra.inter.ith.quadrant.RData')
ith.genes <- Reduce(intersect, list(rownames(subset(losic.ith.quadrant, quadrant == 'bottom_right')),
                                    rownames(subset(renji.ith.quadrant, quadrant == 'bottom_right')),
                                    rownames(subset(shen.ith.quadrant, quadrant == 'bottom_right')),
                                    rownames(subset(shi.ith.quadrant, quadrant == 'bottom_right'))))



tcga.cox.pvalue <- data.frame(Gene=names(tcga.univ.cox.pvalue), P=tcga.univ.cox.pvalue, n_patients=323)
icgc.cox.pvalue <- data.frame(Gene=names(icgc.univ.cox.pvalue), P=icgc.univ.cox.pvalue, n_patients=203)
chcc.cox.pvalue <- data.frame(Gene=names(chcc.univ.cox.pvalue), P=chcc.univ.cox.pvalue, n_patients=159)
fulci.cox.pvalue <- data.frame(Gene=names(fulci.univ.cox.pvalue), P=fulci.univ.cox.pvalue, n_patients=221)
nci.cox.pvalue <- data.frame(Gene=names(nci.univ.cox.pvalue), P=nci.univ.cox.pvalue, n_patients=112)
mongolian.cox.pvalue <- data.frame(Gene=names(mongolian.univ.cox.pvalue), P=mongolian.univ.cox.pvalue, n_patients=70)


tcga.quadrant <- data.frame(Gene = names(tcga.univ.cox.pvalue), 
                            quadrant = ifelse(names(tcga.univ.cox.pvalue) %in% ith.genes, 'bottom_right', 'top_left'), signature = 'COX UVA')
icgc.quadrant <- data.frame(Gene = names(icgc.univ.cox.pvalue), 
                            quadrant = ifelse(names(icgc.univ.cox.pvalue) %in% ith.genes, 'bottom_right', 'top_left'), signature = 'COX UVA')
chcc.quadrant <- data.frame(Gene = names(chcc.univ.cox.pvalue), 
                            quadrant = ifelse(names(chcc.univ.cox.pvalue) %in% ith.genes, 'bottom_right', 'top_left'), signature = 'COX UVA')
fulci.quadrant <- data.frame(Gene = names(fulci.univ.cox.pvalue), 
                             quadrant = ifelse(names(fulci.univ.cox.pvalue) %in% ith.genes, 'bottom_right', 'top_left'), signature = 'COX UVA')
nci.quadrant <- data.frame(Gene = names(nci.univ.cox.pvalue), 
                           quadrant = ifelse(names(nci.univ.cox.pvalue) %in% ith.genes, 'bottom_right', 'top_left'), signature = 'COX UVA')
mongolian.quadrant <- data.frame(Gene = names(mongolian.univ.cox.pvalue), 
                                 quadrant = ifelse(names(mongolian.univ.cox.pvalue) %in% ith.genes, 'bottom_right', 'top_left'), signature = 'COX UVA')


quadrant_cox_pvalue_plot(tcga.quadrant, tcga.cox.pvalue, 'ith_genes_tcga_coxp.pdf', cohort_name = 'TCGA')
quadrant_cox_pvalue_plot(icgc.quadrant, icgc.cox.pvalue, 'ith_genes_icgc_coxp.pdf', cohort_name = 'ICGC')
quadrant_cox_pvalue_plot(chcc.quadrant, chcc.cox.pvalue, 'ith_genes_chcc_coxp.pdf', cohort_name = 'CHCC')
quadrant_cox_pvalue_plot(fulci.quadrant, fulci.cox.pvalue, 'ith_genes_fulci_coxp.pdf', cohort_name = 'FULCI')
quadrant_cox_pvalue_plot(nci.quadrant, nci.cox.pvalue, 'ith_genes_nci_coxp.pdf', cohort_name = 'NCI')
quadrant_cox_pvalue_plot(mongolian.quadrant, mongolian.cox.pvalue, 'ith_genes_mongolian_coxp.pdf', cohort_name = 'Mongolian')


# Fisher test
fisher.test(table(merge(tcga.quadrant, tcga.cox.pvalue, by = 'Gene') %>% 
                    mutate(sig = ifelse(P < 0.05, 'sig', 'unsig'), quadrant = factor(quadrant, level = c('bottom_right', 'top_left'))) %>% 
                    dplyr::select(sig, quadrant)))

fisher.test(table(merge(icgc.quadrant, icgc.cox.pvalue, by = 'Gene') %>% 
                    mutate(sig = ifelse(P < 0.05, 'sig', 'unsig'), quadrant = factor(quadrant, level = c('bottom_right', 'top_left'))) %>% 
                    dplyr::select(sig, quadrant)))

fisher.test(table(merge(chcc.quadrant, chcc.cox.pvalue, by = 'Gene') %>% 
                    mutate(sig = ifelse(P < 0.05, 'sig', 'unsig'), quadrant = factor(quadrant, level = c('bottom_right', 'top_left'))) %>% 
                    dplyr::select(sig, quadrant)))

fisher.test(table(merge(fulci.quadrant, fulci.cox.pvalue, by = 'Gene') %>% 
                    mutate(sig = ifelse(P < 0.05, 'sig', 'unsig'), quadrant = factor(quadrant, level = c('bottom_right', 'top_left'))) %>% 
                    dplyr::select(sig, quadrant)))

fisher.test(table(merge(nci.quadrant, nci.cox.pvalue, by = 'Gene') %>% 
                    mutate(sig = ifelse(P < 0.05, 'sig', 'unsig'), quadrant = factor(quadrant, level = c('bottom_right', 'top_left'))) %>% 
                    dplyr::select(sig, quadrant)))

fisher.test(table(merge(mongolian.quadrant, mongolian.cox.pvalue, by = 'Gene') %>% 
                    mutate(sig = ifelse(P < 0.05, 'sig', 'unsig'), quadrant = factor(quadrant, level = c('bottom_right', 'top_left'))) %>% 
                    dplyr::select(sig, quadrant)))



p1 <- ggplot(data = merge(tcga.quadrant, tcga.cox.pvalue, by = 'Gene') %>% 
               mutate(sig = ifelse(P < 0.05, 'sig', 'unsig'), quadrant = factor(quadrant, level = c('bottom_right', 'top_left'))),
             aes(x = quadrant, fill = sig)) + geom_bar(position = "fill", stat = "count") +
  scale_fill_manual(values = c("unsig" = "blue", "sig" = "red"), labels = c("P < 0.05", "P > 0.05")) + 
  scale_x_discrete(labels = c("bottom_right" = "Q1", "top_left" = "Other")) + 
  theme_minimal() + labs(x = NULL, y = "Percentage", title = "TCGA", fill = 'Cox Pvalue') + 
  theme(legend.position = "bottom")


p2 <- ggplot(data = merge(icgc.quadrant, icgc.cox.pvalue, by = 'Gene') %>% 
               mutate(sig = ifelse(P < 0.05, 'sig', 'unsig'), quadrant = factor(quadrant, level = c('bottom_right', 'top_left'))),
             aes(x = quadrant, fill = sig)) + geom_bar(position = "fill", stat = "count") +
  scale_fill_manual(values = c("unsig" = "blue", "sig" = "red"), labels = c("P < 0.05", "P > 0.05")) + 
  scale_x_discrete(labels = c("bottom_right" = "Q1", "top_left" = "Other")) + 
  theme_minimal() + labs(x = NULL, y = "Percentage", title = "ICGC", fill = 'Cox Pvalue')+
  theme(legend.position = "bottom")


p3 <- ggplot(data = merge(chcc.quadrant, chcc.cox.pvalue, by = 'Gene') %>% 
               mutate(sig = ifelse(P < 0.05, 'sig', 'unsig'), quadrant = factor(quadrant, level = c('bottom_right', 'top_left'))),
             aes(x = quadrant, fill = sig)) + geom_bar(position = "fill", stat = "count") +
  scale_fill_manual(values = c("unsig" = "blue", "sig" = "red"), labels = c("P < 0.05", "P > 0.05")) + 
  scale_x_discrete(labels = c("bottom_right" = "Q1", "top_left" = "Other")) + 
  theme_minimal() + labs(x = NULL, y = "Percentage", title = "CHCC", fill = 'Cox Pvalue') + 
  theme(legend.position = "bottom")


p4 <- ggplot(data = merge(nci.quadrant, nci.cox.pvalue, by = 'Gene') %>% 
               mutate(sig = ifelse(P < 0.05, 'sig', 'unsig'), quadrant = factor(quadrant, level = c('bottom_right', 'top_left'))),
             aes(x = quadrant, fill = sig)) + geom_bar(position = "fill", stat = "count") +
  scale_fill_manual(values = c("unsig" = "blue", "sig" = "red"), labels = c("P < 0.05", "P > 0.05")) + 
  scale_x_discrete(labels = c("bottom_right" = "Q1", "top_left" = "Other")) + 
  theme_minimal() + labs(x = NULL, y = "Percentage", title = "NCI", fill = 'Cox Pvalue') + 
  theme(legend.position = "bottom")


p5 <- ggplot(data = merge(fulci.quadrant, fulci.cox.pvalue, by = 'Gene') %>% 
               mutate(sig = ifelse(P < 0.05, 'sig', 'unsig'), quadrant = factor(quadrant, level = c('bottom_right', 'top_left'))),
             aes(x = quadrant, fill = sig)) + geom_bar(position = "fill", stat = "count") +
  scale_fill_manual(values = c("unsig" = "blue", "sig" = "red"), labels = c("P < 0.05", "P > 0.05")) + 
  scale_x_discrete(labels = c("bottom_right" = "Q1", "top_left" = "Other")) + 
  theme_minimal() + labs(x = NULL, y = "Percentage", title = "FULCI", fill = 'Cox Pvalue') + 
  theme(legend.position = "bottom")


p6 <- ggplot(data = merge(fulci.quadrant, fulci.cox.pvalue, by = 'Gene') %>% 
               mutate(sig = ifelse(P < 0.05, 'sig', 'unsig'), quadrant = factor(quadrant, level = c('bottom_right', 'top_left'))),
             aes(x = quadrant, fill = sig)) + geom_bar(position = "fill", stat = "count") +
  scale_fill_manual(values = c("unsig" = "blue", "sig" = "red"), labels = c("P < 0.05", "P > 0.05")) + 
  scale_x_discrete(labels = c("bottom_right" = "Q1", "top_left" = "Other")) + 
  theme_minimal() + labs(x = NULL, y = "Percentage", title = "FULCI", fill = 'Cox Pvalue') + 
  theme(legend.position = "bottom")


com.plot <- ggarrange(p1, p2, p3, p4, p5, p6, labels = c("A", "B", "C", "D", "E", "F"), ncol = 3, nrow = 2)

ggsave(plot = com.plot, filename = '/result/section1/ith_genes_cox_sig_enrich.pdf')




#################################
# load(file='/data/intra.inter.ith.quadrant.RData')
# losic.ith.quadrant, renji.ith.quadrant, shen.ith.quadrant, shi.ith.quadrant

# losic.ith.quadrant <- losic.ith.quadrant %>% rownames_to_column(var = 'Gene') %>% dplyr::select(Gene, quadrant) %>% mutate(signature = 'COX UVA')
# renji.ith.quadrant <- renji.ith.quadrant %>% rownames_to_column(var = 'Gene') %>% dplyr::select(Gene, quadrant) %>% mutate(signature = 'COX UVA') 
# shen.ith.quadrant <- shen.ith.quadrant %>% rownames_to_column(var = 'Gene') %>% dplyr::select(Gene, quadrant) %>% mutate(signature = 'COX UVA')
# shi.ith.quadrant <- shi.ith.quadrant %>% rownames_to_column(var = 'Gene') %>% dplyr::select(Gene, quadrant) %>% mutate(signature = 'COX UVA')

# 
# #################################TCGA
# load('/data/univ_cox_pvalue.RData')
# tcga.cox.pvalue <- data.frame(Gene=names(tcga.univ.cox.pvalue), P=tcga.univ.cox.pvalue, n_patients=323)
# 
# quadrant_cox_pvalue_plot(losic.ith.quadrant, tcga.cox.pvalue, 'losic_sd_quadrant_tcga_coxp.pdf', cohort_name = 'TCGA')
# quadrant_cox_pvalue_plot(renji.ith.quadrant, tcga.cox.pvalue, 'renji_sd_quadrant_tcga_coxp.pdf', cohort_name = 'TCGA')
# quadrant_cox_pvalue_plot(shen.ith.quadrant, tcga.cox.pvalue, 'shen_sd_quadrant_tcga_coxp.pdf', cohort_name = 'TCGA')
# quadrant_cox_pvalue_plot(shi.ith.quadrant, tcga.cox.pvalue, 'shi_sd_quadrant_tcga_coxp.pdf', cohort_name = 'TCGA')
# 
# 
# 
# #################################ICGC
# icgc.cox.pvalue <- data.frame(Gene=names(icgc.univ.cox.pvalue), P=icgc.univ.cox.pvalue, n_patients=203)
# 
# quadrant_cox_pvalue_plot(losic.ith.quadrant, icgc.cox.pvalue, 'losic_sd_quadrant_icgc_coxp.pdf', cohort_name = 'ICGC')
# quadrant_cox_pvalue_plot(renji.ith.quadrant, icgc.cox.pvalue, 'renji_sd_quadrant_icgc_coxp.pdf', cohort_name = 'TCGA')
# quadrant_cox_pvalue_plot(shen.ith.quadrant, icgc.cox.pvalue, 'shen_sd_quadrant_icgc_coxp.pdf', cohort_name = 'ICGC')
# quadrant_cox_pvalue_plot(shi.ith.quadrant, icgc.cox.pvalue, 'shi_sd_quadrant_icgc_coxp.pdf', cohort_name = 'ICGC')
# 
# 
# #################################CHCC
# chcc.cox.pvalue <- data.frame(Gene=names(chcc.univ.cox.pvalue), P=chcc.univ.cox.pvalue, n_patients=159)
# 
# quadrant_cox_pvalue_plot(losic.ith.quadrant, chcc.cox.pvalue, 'losic_sd_quadrant_chcc_coxp.pdf', cohort_name = 'CHCC')
# quadrant_cox_pvalue_plot(renji.ith.quadrant, chcc.cox.pvalue, 'renji_sd_quadrant_chcc_coxp.pdf', cohort_name = 'CHCC')
# quadrant_cox_pvalue_plot(shen.ith.quadrant, chcc.cox.pvalue, 'shen_sd_quadrant_chcc_coxp.pdf', cohort_name = 'CHCC')
# quadrant_cox_pvalue_plot(shi.ith.quadrant, chcc.cox.pvalue, 'shi_sd_quadrant_chcc_coxp.pdf', cohort_name = 'CHCC')
# 
# 
# #################################FULCI
# fulci.cox.pvalue <- data.frame(Gene=names(fulci.univ.cox.pvalue), P=fulci.univ.cox.pvalue, n_patients=221)
# 
# quadrant_cox_pvalue_plot(losic.ith.quadrant, fulci.cox.pvalue, 'losic_sd_quadrant_fulci_coxp.pdf', cohort_name = 'FULCI')
# quadrant_cox_pvalue_plot(renji.ith.quadrant, fulci.cox.pvalue, 'renji_sd_quadrant_fulci_coxp.pdf', cohort_name = 'FULCI')
# quadrant_cox_pvalue_plot(shen.ith.quadrant, fulci.cox.pvalue, 'shen_sd_quadrant_fulci_coxp.pdf', cohort_name = 'FULCI')
# quadrant_cox_pvalue_plot(shi.ith.quadrant, fulci.cox.pvalue, 'shi_sd_quadrant_fulci_coxp.pdf', cohort_name = 'FULCI')
# 
# 
# #################################NCI
# nci.cox.pvalue <- data.frame(Gene=names(nci.univ.cox.pvalue), P=nci.univ.cox.pvalue, n_patients=112)
# 
# quadrant_cox_pvalue_plot(losic.ith.quadrant, nci.cox.pvalue, 'losic_sd_quadrant_nci_coxp.pdf', cohort_name = 'NCI')
# quadrant_cox_pvalue_plot(renji.ith.quadrant, nci.cox.pvalue, 'renji_sd_quadrant_nci_coxp.pdf', cohort_name = 'NCI')
# quadrant_cox_pvalue_plot(shen.ith.quadrant, nci.cox.pvalue, 'shen_sd_quadrant_nci_coxp.pdf', cohort_name = 'NCI')
# quadrant_cox_pvalue_plot(shi.ith.quadrant, nci.cox.pvalue, 'shi_sd_quadrant_fnci_coxp.pdf', cohort_name = 'NCI')
# 
# 
# 
# #################################Mongolian
# mongolian.cox.pvalue <- data.frame(Gene=names(mongolian.univ.cox.pvalue), P=mongolian.univ.cox.pvalue, n_patients=70)
# 
# quadrant_cox_pvalue_plot(losic.ith.quadrant, mongolian.cox.pvalue, 'losic_sd_quadrant_mongolian_coxp.pdf', cohort_name = 'Mongolian')
# quadrant_cox_pvalue_plot(renji.ith.quadrant, mongolian.cox.pvalue, 'renji_sd_quadrant_mongolian_coxp.pdf', cohort_name = 'Mongolian')
# quadrant_cox_pvalue_plot(shen.ith.quadrant, mongolian.cox.pvalue, 'shen_sd_quadrant_mongolian_coxp.pdf', cohort_name = 'Mongolian')
# quadrant_cox_pvalue_plot(shi.ith.quadrant, mongolian.cox.pvalue, 'shi_sd_quadrant_mongolian_coxp.pdf', cohort_name = 'Mongolian')
# 






evo.genes <- readRDS(file='/data/evo_genes.rds')
ith.quadrant <- data.frame(Gene = names(tcga.univ.cox.pvalue), quadrant = ifelse(names(tcga.univ.cox.pvalue) %in% evo.genes, 'top_right', 'top_left'), signature = 'COX UVA')

quadrant_cox_pvalue_plot(ith.quadrant, tcga.cox.pvalue, 'quadrant_tcga_coxp.pdf', cohort_name = 'TCGA')


fisher.test(table(merge(ith.quadrant, tcga.cox.pvalue, by = 'Gene') %>% 
  mutate(sig = ifelse(P < 0.05, 'sig', 'unsig'), quadrant = factor(quadrant, level = c('top_right', 'top_left'))) %>% dplyr::select(sig, quadrant)))




###############################
precog_z <- read.delim(file="/data/PRECOG-metaZ.pcl")

quadrant.gene$Gene <- gene.info$Symbol[match(quadrant.gene$Gene, gene.info$GeneID)]
quadrant_precog_zscore_plot(quadrant.gene, pancancer=FALSE, 'shi_cv_quadrant_precogz_liver.pdf')


