
adjuvantTACE <- c('LCS_007A','LCS_009A','LCS_019A','LCS_020A','LCS_025A','LCS_027A','LCS_028A','LCS_029A','LCS_031A','LCS_033A',
'LCS_034A','LCS_038A','LCS_043A','LCS_047A','LCS_049A','LCS_050A','LCS_062A','LCS_065A','LCS_068A','LCS_071A',
'LCS_075A','LCS_079A','LCS_085A','LCS_086A','LCS_092A','LCS_097A','LCS_100A','LCS_104A','LCS_110A','LCS_116A',
'LCS_117A','LCS_118A','LCS_121A','LCS_127A','LCS_134A','LCS_136A','LCS_140A','LCS_142A','LCS_146A','LCS_154A',
'LCS_158A','LCS_159A','LCS_161A','LCS_166A','LCS_167A','LCS_170A','LCS_171A','LCS_177A','LCS_185A','LCS_191A',
'LCS_192A','LCS_196A','LCS_197A','LCS_208A','LCS_209A','LCS_212A','LCS_213A','LCS_223A','LCS_228A','LCS_231A',
'LCS_240A','LCS_241A','LCS_245A','LCS_251A','LCS_259A','LCS_260A','LCS_263A','LCS_264A','LCS_265A','LCS_266A',
'LCS_270A','LCS_272A','LCS_284A','LCS_289A','LCS_393A')

postTACE <- c('LCS_008A','LCS_012A','LCS_023A','LCS_024A','LCS_032A','LCS_035A','LCS_067A','LCS_072A','LCS_088A','LCS_096A',
'LCS_120A','LCS_138A','LCS_139A','LCS_145A','LCS_178A','LCS_190A','LCS_194A','LCS_198A','LCS_200A','LCS_207A',
'LCS_224A','LCS_227A','LCS_234A','LCS_238A','LCS_267A','LCS_273A','LCS_274A','LCS_281A','LCS_333A','LCS_403A')

resection <- c('LCS_010A','LCS_014A','LCS_015A','LCS_016A','LCS_018A','LCS_022A','LCS_040A','LCS_041A','LCS_044A','LCS_045A',
'LCS_046A','LCS_048A','LCS_051A','LCS_056A','LCS_057A','LCS_061A','LCS_063A','LCS_064A','LCS_069A','LCS_073A',
'LCS_076A','LCS_078A','LCS_084A','LCS_090A','LCS_091A','LCS_094A','LCS_099A','LCS_101A','LCS_102A','LCS_105A',
'LCS_106A','LCS_108A','LCS_109A','LCS_119A','LCS_122A','LCS_130A','LCS_131A','LCS_132A','LCS_137A','LCS_144A',
'LCS_147A','LCS_150A','LCS_151A','LCS_156A','LCS_160A','LCS_163A','LCS_165A','LCS_169A','LCS_172A','LCS_174A',
'LCS_179A','LCS_180A','LCS_184A','LCS_189A','LCS_205A','LCS_210A','LCS_211A','LCS_215A','LCS_216A','LCS_222A',
'LCS_236A','LCS_237A','LCS_243A','LCS_247A','LCS_249A','LCS_253A','LCS_254A','LCS_261A','LCS_262A','LCS_268A',
'LCS_269A','LCS_275A','LCS_278A','LCS_279A','LCS_282A','LCS_285A','LCS_286A','LCS_291A','LCS_343A','LCS_344A',
'LCS_346A','LCS_400A','LCS_406A','LCS_415A','LCS_424A','LCS_426A')

fulci.linc.cli.data <- readRDS(file='/data/fulci_risk_score.rds')


TACEsams <- fulci.linc.cli.data %>% subset.data.frame(LCS.ID %in% c(adjuvantTACE, postTACE)) %>% 
  mutate(risk.categ = ifelse(risk.score >= median(risk.score), 'high risk', 'low risk'))

adjuvantTACEsams <- fulci.linc.cli.data %>% subset.data.frame(LCS.ID %in% adjuvantTACE) %>% 
  mutate(risk.categ = ifelse(risk.score >= median(risk.score), 'high risk', 'low risk'))

postTACEsams <- fulci.linc.cli.data %>% subset.data.frame(LCS.ID %in% postTACE) %>% 
  mutate(risk.categ = ifelse(risk.score >= median(risk.score), 'high risk', 'low risk'))

resectionSams <- fulci.linc.cli.data %>% subset.data.frame(LCS.ID %in% resection) %>% 
  mutate(risk.categ = ifelse(risk.score >= median(risk.score), 'high risk', 'low risk'))


# survival plot 
source('/code/Rscript/survival_plot.R')

SurvivalPlot(survival.data=TACEsams[, c('LCS.ID', 'Survival.months', 'Survival.status')], 
             sample.class=TACEsams[, c('LCS.ID', 'risk.categ')], filename='fulci_lihc_TACE_os.pdf', 
             out.file.path='/result/Section6/')

SurvivalPlot(survival.data=TACEsams[, c('LCS.ID', 'Recurr.months', 'Recurr.status')], 
             sample.class=TACEsams[, c('LCS.ID', 'risk.categ')], filename='fulci_lihc_TACE_rfs.pdf', 
             out.file.path='/result/Section6/')


SurvivalPlot(survival.data=adjuvantTACEsams[, c('LCS.ID', 'Survival.months', 'Survival.status')], 
             sample.class=adjuvantTACEsams[, c('LCS.ID', 'risk.categ')], filename='fulci_lihc_adjuvantTACE_os.pdf', 
             out.file.path='/result/Section6/')

SurvivalPlot(survival.data=adjuvantTACEsams[, c('LCS.ID', 'Recurr.months', 'Recurr.status')], 
             sample.class=adjuvantTACEsams[, c('LCS.ID', 'risk.categ')], filename='fulci_lihc_adjuvantTACE_rfs.pdf', 
             out.file.path='/result/Section6/')


SurvivalPlot(survival.data=postTACEsams[, c('LCS.ID', 'Survival.months', 'Survival.status')], 
             sample.class=postTACEsams[, c('LCS.ID', 'risk.categ')], filename='fulci_lihc_postTACE_os.pdf', 
             out.file.path='/result/Section6/')

SurvivalPlot(survival.data=postTACEsams[, c('LCS.ID', 'Recurr.months', 'Recurr.status')], 
             sample.class=postTACEsams[, c('LCS.ID', 'risk.categ')], filename='fulci_lihc_postTACE_rfs.pdf', 
             out.file.path='/result/Section6/')



SurvivalPlot(survival.data=resectionSams[, c('LCS.ID', 'Survival.months', 'Survival.status')], 
             sample.class=resectionSams[, c('LCS.ID', 'risk.categ')], filename='fulci_lihc_resection_os.pdf', 
             out.file.path='/result/Section6/')

SurvivalPlot(survival.data=resectionSams[, c('LCS.ID', 'Recurr.months', 'Recurr.status')], 
             sample.class=resectionSams[, c('LCS.ID', 'risk.categ')], filename='fulci_lihc_resection_rfs.pdf', 
             out.file.path='/result/Section6/')




