
# Published signature
library(readxl)
lihc.signature <- read_excel(path = '/data/Signature.xlsx', sheet =1, col_names = TRUE)
evo_genes <- readRDS('/data/evo_genes.rds')

lihc.signature <- lihc.signature %>% mutate(Quadrant = ifelse(GeneName %in% evo_genes, 'Q1', 'Other'))

sig.q1.percentage <- lihc.signature %>% subset(Quadrant == 'Q1') %>% group_by(PMID) %>% count(Quadrant, name = 'Q1num') %>% 
  select(PMID, Q1num) %>% right_join(lihc.signature, by = 'PMID') %>% mutate(Q1num = ifelse(is.na(Q1num), 0, Q1num)) %>% 
  mutate(Q1perc = Q1num/Size) %>% select(PMID, Q1num, Q1perc, Signature, IF, Size, Name, Cited) %>% distinct()


# saveRDS(sig.q1.percentage, file = '/data/sig_q1_percentage.rds')

#####################################plot
load('/data/intra.inter.ith.quadrant.RData')

risk.sig <- c('PMID_33251144', 'PMID_32198063', 'PMID_35123387', 'PMID_31335995', 'PMID_33033585',
              'PMID_33828988', 'PMID_34676211', 'PMID_32903581', 'PMID_34900672', 'PMID_33089373', 'PMID_35311113', 
              'PMID_34975331', 'PMID_25666192', 'PMID_22105560', 'PMID_23800896')

lihc.signature <- lihc.signature %>% mutate(PMID = paste0('PMID_', PMID))
sig15.stat <- lihc.signature %>% subset(PMID %in% risk.sig) %>% distinct(GeneName, .keep_all = T) %>% count(Quadrant) %>% mutate(type = 'Sig15')
sig65.stat <- lihc.signature %>% distinct(GeneName, .keep_all = T) %>% count(Quadrant) %>% mutate(type = 'Sig65')

losic.quadrant.stat <- data.frame(Quadrant = c('Q1', 'Other'), 
                                  n = c(length(evo_genes), length(setdiff(rownames(losic.ith.quadrant), evo_genes))), type = 'Losic')

renji.quadrant.stat <- data.frame(Quadrant = c('Q1', 'Other'), 
                                  n = c(length(evo_genes), length(setdiff(rownames(renji.ith.quadrant), evo_genes))), type = 'Renji')

shen.quadrant.stat <- data.frame(Quadrant = c('Q1', 'Other'), 
                                  n = c(length(evo_genes), length(setdiff(rownames(shen.ith.quadrant), evo_genes))), type = 'Shen')

shi.quadrant.stat <- data.frame(Quadrant = c('Q1', 'Other'), 
                                  n = c(length(evo_genes), length(setdiff(rownames(shi.ith.quadrant), evo_genes))), type = 'Shi')


quadrant.stat <- rbind.data.frame(sig15.stat, sig65.stat, losic.quadrant.stat, renji.quadrant.stat, shen.quadrant.stat, shi.quadrant.stat)




p1 <- ggplot(data = quadrant.stat %>% mutate(Quadrant = factor(Quadrant, level = c('Q1', 'Other'))),
       aes(x = type, y = n, fill = Quadrant)) + geom_bar(position = "fill", stat = 'identity') +
  scale_fill_manual(values = c("Other" = "blue", "Q1" = "red"), labels = c("EvoGene", "Other")) + 
  scale_x_discrete(labels = c("Sig15" = "Observed (15)", "Sig65" = "Observed (65)", 
                              "Losic" = "Expected (Losic)", "Renji" = "Expected (Renji)",
                              "Shen" = "Observed (Shen)", "Shi" = "Observed (Shi)")) + 
  theme_minimal() + labs(x = NULL, y = "Percentage") + theme(legend.position = "bottom")

ggsave(plot = p1, filename = '/result/section1/sig_enrich_in_evogenes.pdf')


# fisher.test

quadrant.stat <- quadrant.stat %>% tidyr::pivot_wider(names_from = type, values_from = n, ) %>% 
  arrange(desc(Quadrant))


# > fisher.test(quadrant.stat[, c('Sig15', 'Losic')])$p.value
# [1] 7.556616e-10
# > fisher.test(quadrant.stat[, c('Sig15', 'Losic')])$estimate
# odds ratio 
# 4.553982 
# > fisher.test(quadrant.stat[, c('Sig15', 'Losic')])$conf.int
# [1] 2.903435 6.916121


# > fisher.test(quadrant.stat[, c('Sig15', 'Renji')])$p.value
# [1] 1.255977e-08
# > fisher.test(quadrant.stat[, c('Sig15', 'Renji')])$estimate
# odds ratio 
# 3.972447 
# > fisher.test(quadrant.stat[, c('Sig15', 'Renji')])$conf.int
# [1] 2.532076 6.033561


# > fisher.test(quadrant.stat[, c('Sig15', 'Shen')])$p.value
# [1] 8.74066e-12
# > fisher.test(quadrant.stat[, c('Sig15', 'Shen')])$estimate
# odds ratio 
# 5.595125 
# > fisher.test(quadrant.stat[, c('Sig15', 'Shen')])$conf.int
# [1] 3.567422 8.498287



# > fisher.test(quadrant.stat[, c('Sig15', 'Shi')])$p.value
# [1] 4.043636e-10
# > fisher.test(quadrant.stat[, c('Sig15', 'Shi')])$estimate
# odds ratio 
# 4.691021 
# > fisher.test(quadrant.stat[, c('Sig15', 'Shi')])$conf.int
# [1] 2.990481 7.125645


# 
# > fisher.test(quadrant.stat[, c('Sig65', 'Losic')])$p.value
# [1] 1.430474e-05
# > fisher.test(quadrant.stat[, c('Sig65', 'Losic')])$estimate
# odds ratio 
# 2.493982 
# > fisher.test(quadrant.stat[, c('Sig65', 'Losic')])$conf.int
# [1] 1.669435 3.619756
# 
# 
# > fisher.test(quadrant.stat[, c('Sig65', 'Renji')])$p.value
# [1] 0.0001525831
# > fisher.test(quadrant.stat[, c('Sig65', 'Renji')])$estimate
# odds ratio 
# 2.175265 
# > fisher.test(quadrant.stat[, c('Sig65', 'Renji')])$conf.int
# [1] 1.455990 3.157432
# 
# 
# > fisher.test(quadrant.stat[, c('Sig65', 'Shen')])$p.value
# [1] 1.427113e-07
# > fisher.test(quadrant.stat[, c('Sig65', 'Shen')])$estimate
# odds ratio 
# 3.064191 
# > fisher.test(quadrant.stat[, c('Sig65', 'Shen')])$conf.int
# [1] 2.051282 4.446841
# 
# > fisher.test(quadrant.stat[, c('Sig65', 'Shi')])$p.value
# [1] 5.756427e-06
# > fisher.test(quadrant.stat[, c('Sig65', 'Shi')])$estimate
# odds ratio 
# 2.568789 
# > fisher.test(quadrant.stat[, c('Sig65', 'Shi')])$conf.int
# [1] 1.719615 3.728452

