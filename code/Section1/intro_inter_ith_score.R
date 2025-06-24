
# Function
IntraHetScore <- function(exp.dat, sample.list, method){
  
  ith.scores <- lapply(sample.list, function(sample.set){
    pat.exp.dat <- exp.dat[, sample.set]
    
    if(method == 'sd'){
      ith.score <- apply(pat.exp.dat, 1, sd)
    }else if(method == 'mad'){
      ith.score <- apply(pat.exp.dat, 1, mad)
    }else{
      ith.score <- apply(pat.exp.dat, 1, function(x) sd(x)/mean(x))
    }
    
    return(ith.score)
  })
  ith.scores <- do.call(cbind, ith.scores)
  return(rowMeans(ith.scores))
}

InterHetScore <- function(exp.dat, sample.list, method){
  
  set.seed(41)
  sample.list <- lapply(1:10, function(i){
    sample.set <- sapply(sample.list, function(sample.set) sample(sample.set, 1))
    return(sample.set)
  })
  
  ith.scores <- lapply(sample.list, function(sample.set){
    pat.exp.dat <- exp.dat[, sample.set]
    
    if(method == 'sd'){
      ith.score <- apply(pat.exp.dat, 1, sd)
    }else if(method == 'mad'){
      ith.score <- apply(pat.exp.dat, 1, mad)
    }else{
      ith.score <- apply(pat.exp.dat, 1, function(x) sd(x)/mean(x))
    }
    
    return(ith.score)
  })
  ith.scores <- do.call(cbind, ith.scores)
  return(rowMeans(ith.scores))
}


# Gene expression from Losic et al.
load(file='/data/losic.refseq.gene.exp.RData')

losic.gene.exp <- tran.exp$abundance
rownames(losic.gene.exp) <- stringr::str_split_i(rownames(losic.gene.exp), pattern = '_', i = 2)
losic.gene.exp <- losic.gene.exp[, grep(pattern = 'H', colnames(losic.gene.exp), value = TRUE)]

# keeping genes with an expression value of at least 1 TPM in at least 70% (9/43) of tumor samples
losic.gene.exp <- losic.gene.exp[rowSums(losic.gene.exp >= 1) >= 30, ]
losic.gene.exp <- as.data.frame(log2(losic.gene.exp + 1))
# 12429

losic.pats <- list(
  # H1 = c("H1.a", "H1.b"),
  H2 = c("H2.a1", "H2.b", "H2.c", "H2.d", "H2.e"), # "H2.a2", 
  # H3 = c("H3.a", "H3.b"),
  H4 = c("H4.a", "H4.b", "H4.c", "H4.d", "H4.e"), 
  # H6 = c("H6.a", "H6.b"),
  H7 = c("H7.a", "H7.b", "H7.c", "H7.d", "H7.e"), 
  H8 = c("H8.a", "H8.b", "H8.c"),
  H9 = c("H9.a", "H9.b", "H9.d", "H9.e", "H9.f"), # "H9.c", 
  H10 = c("H10.a", "H10.b", "H10.c", "H10.d", "H10.e"),
  # H11 = c("H11.a", "H11.b"),
  H12 = c("H12.a", "H12.b", "H12.c", "H12.d", "H12.e"))


losic.ith.score <- data.frame(losic.intra.sd = IntraHetScore(losic.gene.exp, losic.pats, 'sd'),
                              losic.intra.mad = IntraHetScore(losic.gene.exp, losic.pats, 'mad'),
                              losic.intra.cv = IntraHetScore(losic.gene.exp, losic.pats, 'cv'),
                              losic.inter.sd = InterHetScore(losic.gene.exp, losic.pats, 'sd'),
                              losic.inter.mad = InterHetScore(losic.gene.exp, losic.pats, 'mad'),
                              losic.inter.cv = InterHetScore(losic.gene.exp, losic.pats, 'cv'))



# Gene expression from Renji.
renji.gene.exp <- read.csv(file = '/data/Renji/Renji_cohort_MR_exp.matrix', 
                     sep='\t', row.names = 1, header=TRUE)

renji.gene.exp <- renji.gene.exp[, grep(pattern = 'T', colnames(renji.gene.exp), value = T)]

# keeping genes with an expression value of at least 1 TPM in at least 70% (53/75) of tumor samples
renji.gene.exp <- renji.gene.exp[rowSums(2^renji.gene.exp - 1 >= 1) >= 53, ]
# 10898

renji.pats <- list(T2 = c("T02_01", "T02_02", "T02_03", "T02_04", "T02_05", "T02_06"),
                 T4 = c("T04_01", "T04_02", "T04_03", "T04_04"),
                 T5 = c("T05_01", "T05_02", "T05_03", "T05_04"),
                 T6 = c("T06_01", "T06_02", "T06_04", "T06_05", "T06_06"),
                 T7 = c("T07_01", "T07_02", "T07_03", "T07_04"),
                 T8 = c("T08_01", "T08_02", "T08_03", "T08_04", "T08_05"), 
                 T9 = c("T09_01", "T09_02", "T09_03", "T09_04", "T09_05"), 
                 T10 = c("T10_01", "T10_02", "T10_03", "T10_04", "T10_06", "T10_07"),
                 T11 = c("T11_01", "T11_02", "T11_03", "T11_04"),
                 T13 = c("T13_01", "T13_02", "T13_03", "T13_04", "T13_05", "T13_06", "T13_07", "T13_08", "T13_09", "T13_10"),
                 # T14 = c("T14_01", "T14_02", "T14_03"),
                 T16 = c("T16_01", "T16_02", "T16_03", "T16_04"),
                 T17 = c("T17_01", "T17_02", "T17_03", "T17_04", "T17_05"), 
                 T18 = c("T18_01", "T18_02", "T18_03", "T18_04", "T18_05", "T18_06", "T18_07", "T18_08", "T18_09", "T18_10"))


renji.ith.score <- data.frame(renji.intra.sd = IntraHetScore(renji.gene.exp, renji.pats, 'sd'),
                              renji.intra.mad = IntraHetScore(renji.gene.exp, renji.pats, 'mad'),
                              renji.intra.cv = IntraHetScore(renji.gene.exp, renji.pats, 'cv'),
                              renji.inter.sd = InterHetScore(renji.gene.exp, renji.pats, 'sd'),
                              renji.inter.mad = InterHetScore(renji.gene.exp, renji.pats, 'mad'),
                              renji.inter.cv = InterHetScore(renji.gene.exp, renji.pats, 'cv'))


# Gene expression from Shi et al.
shi.gene.exp <- readRDS(file='/data/shi.lihc.gene.exp.rds')
load(file='/data/exp_gene_anno.RData')

rownames(shi.gene.exp) <- exp.gene.anno$Symbol[match(shi.gene.exp$GENE, exp.gene.anno$GeneID)]
shi.gene.exp <- shi.gene.exp[, -c(1:3)]
# 12789

shi.pats <- list(
  H1 = c("H1.1", "H1.2", "H1.3", "H1.4", "H1.5"),
  H2 = c("H2.1", "H2.2", "H2.3", "H2.4", "H2.5"), 
  H3 = c("H3.1", "H3.2", "H3.3", "H3.4", "H3.5"), 
  H4 = c("H4.1", "H4.2", "H4.3", "H4.4", "H4.5"), 
  H5 = c("H5.1", "H5.2", "H5.3", "H5.4", "H5.5"))

shi.ith.score <- data.frame(shi.intra.sd = IntraHetScore(shi.gene.exp, shi.pats, 'sd'),
                            shi.intra.mad = IntraHetScore(shi.gene.exp, shi.pats, 'mad'),
                            shi.intra.cv = IntraHetScore(shi.gene.exp, shi.pats, 'cv'),
                            shi.inter.sd = InterHetScore(shi.gene.exp, shi.pats, 'sd'),
                            shi.inter.mad = InterHetScore(shi.gene.exp, shi.pats, 'mad'),
                            shi.inter.cv = InterHetScore(shi.gene.exp, shi.pats, 'cv'))


# Gene expression from Shen et al.
shen.gene.exp <- readRDS(file='/data/GSE136711_vst_exp.RData')


shen.pats <- list(H1 = c("HCC01.IT2",  "HCC01.IT3",  "HCC01.IT4",  "HCC01.IT5"),  
                  H2 = c("HCC02.IT1",  "HCC02.IT3",  "HCC02.IT6",  "HCC02.IT9",  "HCC02.IT12"),
                  H3 = c("HCC03.IT4",  "HCC03.IT5",  "HCC03.IT8",  "HCC03.IT9"), 
                  H4 = c("HCC04.IT2",  "HCC04.IT3",  "HCC04.IT4",  "HCC04.IT5"), 
                  H5 = c("HCC05.IT1",  "HCC05.IT2",  "HCC05.IT4",  "HCC05.IT6"),  
                  H7 = c("HCC07.IT2",  "HCC07.IT4",  "HCC07.IT5"),  
                  H8 = c("HCC08.IT1",  "HCC08.IT4",  "HCC08.IT5"), 
                  H9 = c("HCC09.IT1",  "HCC09.IT3",  "HCC09.IT4"),  
                  H11 = c("HCC11.IT1",  "HCC11.IT4",  "HCC11.IT5"),  
                  H12 = c("HCC12.IT1",  "HCC12.IT2",  "HCC12.IT6"),  
                  # H13 = c("HCC13.IT1",  "HCC13.IT2"),
                  H14 = c("HCC14.IT1",  "HCC14.IT3",  "HCC14.IT5"))

shen.ith.score <- data.frame(shen.intra.sd = IntraHetScore(shen.gene.exp, shen.pats, 'sd'),
                             shen.intra.mad = IntraHetScore(shen.gene.exp, shen.pats, 'mad'),
                             shen.intra.cv = IntraHetScore(shen.gene.exp, shen.pats, 'cv'),
                             shen.inter.sd = InterHetScore(shen.gene.exp, shen.pats, 'sd'),
                             shen.inter.mad = InterHetScore(shen.gene.exp, shen.pats, 'mad'),
                             shen.inter.cv = InterHetScore(shen.gene.exp, shen.pats, 'cv'))


save(losic.ith.score, renji.ith.score, shi.ith.score, shen.ith.score, file='/data/intra.inter.ith.score.RData')


