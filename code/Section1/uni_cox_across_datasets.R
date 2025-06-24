
######################################TCGA
# load data
tcga.lihc.vst.exp <- readRDS(file='/data/tcga_vst_norm_geneExp.rds')
tcga.lihc.cli.data <- readRDS(file='/data/tcga.lihc.cli.data.rds')
load('/data/exp_gene_anno.RData')

tcga.lihc.vst.exp <- tcga.lihc.vst.exp[as.character(na.omit(exp.gene.anno$TCGAExpGeneFilter)), ]

# data preparation
rownames(tcga.lihc.vst.exp) <- paste0('ID_', rownames(tcga.lihc.vst.exp))
colnames(tcga.lihc.vst.exp) <- substr(colnames(tcga.lihc.vst.exp), 1, 12)
tcga.lihc.vst.exp <- as.data.frame(t(tcga.lihc.vst.exp))


tcga.lihc.cli.data <- tcga.lihc.cli.data[, c('os_time', 'os')]
cli.exp.data <- merge(tcga.lihc.cli.data, tcga.lihc.vst.exp, by='row.names')


# Univariate Cox regression analysis
library(survival)
gene.list <- colnames(cli.exp.data)[-c(1:3)]
univ.formulas <- sapply(gene.list, function(x) as.formula(paste('Surv(os_time, os)~', x)))

univ.models <- lapply(univ.formulas, function(x){coxph(x, data = cli.exp.data)})
tcga.univ.cox.pvalue <- sapply(univ.models, function(x){summary(x)$coefficients[, 5]})

names(tcga.univ.cox.pvalue) <- gsub('ID_', '', names(tcga.univ.cox.pvalue))
names(tcga.univ.cox.pvalue) <- exp.gene.anno$Symbol[match(names(tcga.univ.cox.pvalue), exp.gene.anno$GeneID)]

# saveRDS(tcga.univ.cox.pvalue, file='/data/tcga.lihc.univ.cox.pvalue.rds')


######################################ICGC
# load data
icgc.linc.cli.data <- readRDS(file='/data/icgc_linc_cli_data.rds')
icgc.lihc.vst.exp <- readRDS(file='/data/icgc_vst_norm_geneExp.rds')

icgc.lihc.vst.exp <- icgc.lihc.vst.exp[as.character(subset(exp.gene.anno, !is.na(ICGCExpGeneFilter))$GeneID), ]


# data preparation
icgc.linc.cli.data <- icgc.linc.cli.data %>% dplyr::select(icgc_sample_id, donor_survival_time, donor_vital_status) %>% 
  subset(icgc_sample_id %in% colnames(icgc.lihc.vst.exp)) %>% remove_rownames %>% column_to_rownames(var = 'icgc_sample_id') %>% 
dplyr::rename(os_time = donor_survival_time, os = donor_vital_status) %>% mutate(os = ifelse(os == 'deceased', 1, 0))


rownames(icgc.lihc.vst.exp) <- paste0('ID_', rownames(icgc.lihc.vst.exp))
icgc.lihc.vst.exp <- as.data.frame(t(icgc.lihc.vst.exp))

cli.exp.data <- merge(icgc.linc.cli.data, icgc.lihc.vst.exp, by='row.names')


# Univariate Cox regression analysis
library(survival)
gene.list <- colnames(cli.exp.data)[-c(1:3)]
univ.formulas <- sapply(gene.list, function(x) as.formula(paste('Surv(os_time, os)~', x)))

univ.models <- lapply(univ.formulas, function(x){coxph(x, data = cli.exp.data)})
icgc.univ.cox.pvalue <- sapply(univ.models, function(x){summary(x)$coefficients[, 5]})

names(icgc.univ.cox.pvalue) <- gsub('ID_', '', names(icgc.univ.cox.pvalue))
names(icgc.univ.cox.pvalue) <- exp.gene.anno$Symbol[match(names(icgc.univ.cox.pvalue), exp.gene.anno$GeneID)]

# saveRDS(icgc.univ.cox.pvalue, file='/data/icgc.lihc.univ.cox.pvalue.rds')


######################################CHCC
# load data
chcc.lihc.vst.exp <- readRDS(file='/data/Gao_Cell_2019_expr_data.rds')
chcc.linc.cli.data <- readRDS(file='/data/Gao_Cell_2019_cli_data.rds')

rownames(chcc.lihc.vst.exp) <- stringr::str_split_i(rownames(chcc.lihc.vst.exp ), '\\.', i = 1)
chcc.lihc.vst.exp  <- chcc.lihc.vst.exp [intersect(rownames(chcc.lihc.vst.exp ), exp.gene.anno$EnsemblID), ]


# data preparation
chcc.linc.cli.data <- chcc.linc.cli.data %>% dplyr::select(`Overall survial (month)`, `Survial  (1, dead; 0, alive)`) %>% 
 dplyr::rename(os_time = `Overall survial (month)`, os = `Survial  (1, dead; 0, alive)`)

chcc.lihc.vst.exp <- as.data.frame(t(chcc.lihc.vst.exp))

cli.exp.data <- merge(chcc.linc.cli.data, chcc.lihc.vst.exp, by='row.names')


# Univariate Cox regression analysis
library(survival)
gene.list <- colnames(cli.exp.data)[-c(1:3)]
univ.formulas <- sapply(gene.list, function(x) as.formula(paste('Surv(os_time, os)~', x)))

univ.models <- lapply(univ.formulas, function(x){coxph(x, data = cli.exp.data)})
chcc.univ.cox.pvalue <- sapply(univ.models, function(x){summary(x)$coefficients[, 5]})

names(chcc.univ.cox.pvalue) <- exp.gene.anno$Symbol[match(names(chcc.univ.cox.pvalue), exp.gene.anno$EnsemblID)]

# saveRDS(chcc.univ.cox.pvalue, file='/data/chcc.lihc.univ.cox.pvalue.rds')




###################################### Mongolian
# load data
mongolian.lihc.vst.exp <- readRDS(file='/data/GSE144269_expr_data.rds')
mongolian.linc.cli.data <- readRDS(file='/data/GSE144269_cli_data.rds')

rownames(mongolian.lihc.vst.exp) <- stringr::str_split_i(rownames(mongolian.lihc.vst.exp), '\\.', i = 1)
mongolian.lihc.vst.exp <- mongolian.lihc.vst.exp[intersect(rownames(mongolian.lihc.vst.exp), exp.gene.anno$EnsemblID), ]

# data preparation
mongolian.linc.cli.data <- mongolian.linc.cli.data %>% dplyr::select(RNASeq_T, survival.time, survival.status) %>% 
  remove_rownames() %>% column_to_rownames(var = 'RNASeq_T') %>% 
  dplyr::rename(os_time = survival.time, os = survival.status)

mongolian.lihc.vst.exp <- as.data.frame(t(mongolian.lihc.vst.exp))

cli.exp.data <- merge(mongolian.linc.cli.data, mongolian.lihc.vst.exp, by='row.names')


# Univariate Cox regression analysis
library(survival)
gene.list <- colnames(cli.exp.data)[-c(1:3)]
univ.formulas <- sapply(gene.list, function(x) as.formula(paste('Surv(os_time, os)~', x)))

univ.models <- lapply(univ.formulas, function(x){coxph(x, data = cli.exp.data)})
mongolian.univ.cox.pvalue <- sapply(univ.models, function(x){summary(x)$coefficients[, 5]})


names(mongolian.univ.cox.pvalue) <- exp.gene.anno$Symbol[match(names(mongolian.univ.cox.pvalue), exp.gene.anno$EnsemblID)]

# saveRDS(mongolian.univ.cox.pvalue, file='/data/mongolian.lihc.univ.cox.pvalue.rds')



###################################### Fulci
# load data
load(file='/data/lci_xinweiwang_expr_data.RData')
fulci.linc.cli.data <- readRDS(file='/data/lci_xinweiwang_cli_data.rds')

fulci.lihc.vst.exp <- gene.max.exp.profile
fulci.lihc.vst.exp <- fulci.lihc.vst.exp %>% select(-ID, -ENTREZ_GENE_ID, -`Gene Symbol`)
fulci.lihc.vst.exp <- fulci.lihc.vst.exp[intersect(rownames(fulci.lihc.vst.exp), exp.gene.anno$GeneID), ]

# data preparation
fulci.linc.cli.data <- fulci.linc.cli.data %>% subset(Tissue.Type == 'Tumor') %>% dplyr::select(Affy_GSM, Survival.months, Survival.status) %>% 
  remove_rownames() %>% column_to_rownames(var = 'Affy_GSM') %>% 
  dplyr::rename(os_time = Survival.months, os = Survival.status)

fulci.lihc.vst.exp <- fulci.lihc.vst.exp[, rownames(fulci.linc.cli.data)]
rownames(fulci.lihc.vst.exp) <- paste0('ID_', rownames(fulci.lihc.vst.exp))
fulci.lihc.vst.exp <- as.data.frame(t(fulci.lihc.vst.exp))

cli.exp.data <- merge(fulci.linc.cli.data, fulci.lihc.vst.exp, by='row.names')


# Univariate Cox regression analysis
library(survival)
gene.list <- colnames(cli.exp.data)[-c(1:3)]
univ.formulas <- sapply(gene.list, function(x) as.formula(paste('Surv(os_time, os)~', x)))

univ.models <- lapply(univ.formulas, function(x){coxph(x, data = cli.exp.data)})
fulci.univ.cox.pvalue <- sapply(univ.models, function(x){summary(x)$coefficients[, 5]})


names(fulci.univ.cox.pvalue) <- gsub('ID_', '', names(fulci.univ.cox.pvalue))
names(fulci.univ.cox.pvalue) <- exp.gene.anno$Symbol[match(names(fulci.univ.cox.pvalue), exp.gene.anno$GeneID)]

# saveRDS(fulci.univ.cox.pvalue, file='/data/fulci.lihc.univ.cox.pvalue.rds')



###################################### NCI
# load data
load(file='/data/lec_snorri_s_thorgeirsson_expr_data.RData')
nci.linc.cli.data <- readRDS(file='/data/lec_snorri_s_thorgeirsson_cli_data.rds')

nci.lihc.vst.exp <- batch.expr.max.data
nci.lihc.vst.exp <- nci.lihc.vst.exp %>% select(-GeneID, -Symbol, -ID)
nci.lihc.vst.exp <- nci.lihc.vst.exp[intersect(rownames(nci.lihc.vst.exp), exp.gene.anno$GeneID), ]


# data preparation
nci.linc.cli.data <- nci.linc.cli.data %>% dplyr::select(Array, OS_Time, OS_Status) %>% 
  remove_rownames() %>% column_to_rownames(var = 'Array') %>% subset(!is.na(OS_Time) & !is.na(OS_Status)) %>% 
  dplyr::rename(os_time = OS_Time, os = OS_Status)

rownames(nci.lihc.vst.exp) <- paste0('ID_', rownames(nci.lihc.vst.exp))
nci.lihc.vst.exp <- as.data.frame(t(nci.lihc.vst.exp))

cli.exp.data <- merge(nci.linc.cli.data, nci.lihc.vst.exp, by='row.names')


# Univariate Cox regression analysis
library(survival)
gene.list <- colnames(cli.exp.data)[-c(1:3)]
univ.formulas <- sapply(gene.list, function(x) as.formula(paste('Surv(os_time, os)~', x)))

univ.models <- lapply(univ.formulas, function(x){coxph(x, data = cli.exp.data)})
nci.univ.cox.pvalue <- sapply(univ.models, function(x){summary(x)$coefficients[, 5]})


names(nci.univ.cox.pvalue) <- gsub('ID_', '', names(nci.univ.cox.pvalue))
names(nci.univ.cox.pvalue) <- exp.gene.anno$Symbol[match(names(nci.univ.cox.pvalue), exp.gene.anno$GeneID)]

# saveRDS(nci.univ.cox.pvalue, file='/data/nci.lihc.univ.cox.pvalue.rds')


###################################### INSERM
# load data
load(file='/data/E_TABM_36_expr_data.RData')
inserm.linc.cli.data <- readRDS(file='/data/E_TABM_36_cli_data.rds')

inserm.lihc.vst.exp <- gene.max.exp.profile
inserm.lihc.vst.exp <- inserm.lihc.vst.exp %>% select(-GeneID, -GeneSymbol, -ID)
inserm.lihc.vst.exp <- inserm.lihc.vst.exp[intersect(rownames(inserm.lihc.vst.exp), exp.gene.anno$GeneID), ]


# data preparation
inserm.linc.cli.data <- inserm.linc.cli.data %>% subset(DiseaseState == "HCC tumor") %>% dplyr::select(os_time, os_status) %>% 
  subset(!is.na(os_time) & !is.na(os_status)) %>% dplyr::rename(os_time = os_time, os = os_status)

rownames(inserm.lihc.vst.exp) <- paste0('ID_', rownames(inserm.lihc.vst.exp))
inserm.lihc.vst.exp <- as.data.frame(t(inserm.lihc.vst.exp))

cli.exp.data <- merge(inserm.linc.cli.data, inserm.lihc.vst.exp, by='row.names')


# Univariate Cox regression analysis
library(survival)
gene.list <- colnames(cli.exp.data)[-c(1:3)]
univ.formulas <- sapply(gene.list, function(x) as.formula(paste('Surv(os_time, os)~', x)))

univ.models <- lapply(univ.formulas, function(x){coxph(x, data = cli.exp.data)})
inserm.univ.cox.pvalue <- sapply(univ.models, function(x){summary(x)$coefficients[, 5]})


names(inserm.univ.cox.pvalue) <- gsub('ID_', '', names(inserm.univ.cox.pvalue))
names(inserm.univ.cox.pvalue) <- exp.gene.anno$Symbol[match(names(inserm.univ.cox.pvalue), exp.gene.anno$GeneID)]

# saveRDS(inserm.univ.cox.pvalue, file='/data/inserm.lihc.univ.cox.pvalue.rds')


save(chcc.univ.cox.pvalue, fulci.univ.cox.pvalue, icgc.univ.cox.pvalue, mongolian.univ.cox.pvalue, 
     nci.univ.cox.pvalue, tcga.univ.cox.pvalue, inserm.univ.cox.pvalue, file = '/data/univ_cox_pvalue.RData')



