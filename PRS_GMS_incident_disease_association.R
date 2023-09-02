


#########################################################################
# 
# PRS association
# 
#########################################################################



# CHD
prs <- 'CAD_PGS000018'
ds <- 'CHD'; ds.incident <- paste0('INCIDENT_',ds); ds.agediff <- paste0(ds, '_AGEDIFF')
riskf <- c('MEN','BL_AGE','BMI', 'SYSTM','KOL','HDL','CURR_SMOKE', 'exercise','PREVAL_DIAB_GDM','MI_FAMILYHIST')
num <-  c('BL_AGE','BMI','SYSTM', 'KOL', 'HDL')
fac <-  c('MEN','CURR_SMOKE','exercise','PREVAL_DIAB_GDM','MI_FAMILYHIST')
covariates <- setdiff(riskf, c('MEN'))
covariates
strat_sex <- 'strata(MEN)'



# T2D
prs <- 'T2D_PGS000036'
ds <- 'DIAB_T2'; ds.incident <- paste0('INCIDENT_',ds); ds.agediff <- paste0(ds, '_AGEDIFF')
riskf <- c('MEN','BL_AGE','BMI','SYSTM','KOL','HDL','TRIG','CURR_SMOKE','exercise','DIAB_FAMILYHIST')
num <-  c('BL_AGE','BMI','SYSTM','KOL','HDL','TRIG')
fac <-  c('MEN','CURR_SMOKE','exercise','DIAB_FAMILYHIST')
covariates <- setdiff(riskf, c('MEN'))
covariates
strat_sex <- 'strata(MEN)'



# T2D include glucose
prs <- 'T2D_PGS000036'
ds <- 'DIAB_T2'; ds.incident <- paste0('INCIDENT_',ds); ds.agediff <- paste0(ds, '_AGEDIFF')
riskf <- c('MEN','BL_AGE','BMI','SYSTM','KOL','HDL','TRIG','CURR_SMOKE','exercise','DIAB_FAMILYHIST','Glc')
num <-  c('BL_AGE','BMI','SYSTM','KOL','HDL','TRIG','Glc')
fac <-  c('MEN','CURR_SMOKE','exercise','DIAB_FAMILYHIST')
covariates <- setdiff(riskf, c('MEN'))
covariates
strat_sex <- 'strata(MEN)'




# AD
prs <- 'AD_PGS000334'
ds <- 'AD'; ds.incident <- paste0('INCIDENT_',ds); ds.agediff <- paste0(ds, '_AGEDIFF')
riskf <- c('MEN','BL_AGE','BMI','SYSTM','DIASM','KOL','HDL','ALKI2_FR02','CURR_SMOKE','exercise',
           'PREVAL_DIAB_T2','PREVAL_STR_SAH_TIA','PREVAL_MENTAL')
num <-  c('BL_AGE','BMI','SYSTM','DIASM','KOL','HDL','ALKI2_FR02')
fac <-  c('MEN','CURR_SMOKE','exercise','PREVAL_DIAB_T2','PREVAL_STR_SAH_TIA','PREVAL_MENTAL')
covariates <- setdiff(riskf, c('MEN'))
covariates
strat_sex <- 'strata(MEN)'



# PC
prs <- 'prostate_PGS000662'
ds <- 'CR_PROSTCANC'; ds.incident <- paste0('INCIDENT_',ds); ds.agediff <- paste0(ds, '_AGEDIFF')
riskf <- c('BL_AGE','BMI','ALKI2_FR02','CURR_SMOKE','exercise','CANC_FAMILYHIST')
num <-  c('BL_AGE','BMI','ALKI2_FR02')
fac <-  c('CURR_SMOKE','exercise','CANC_FAMILYHIST')
covariates <- setdiff(riskf, c('MEN'))
covariates
strat_sex <- NULL
adjust_sex <- NULL



# --------------------------------------------------------



unlink('.RData')
rm(list = ls())
pcks <- c('data.table','survival','ggplot2','gridExtra','grid','cowplot')
sapply(pcks, require, character.only = TRUE)
sapply(pcks, packageVersion)


dim(pheno <- readRDS(get(paste0('fn.phenodata_', ds)))) 
# pheno[,(fac):=lapply(.SD, as.factor), .SDcols=fac]


lapply(pheno[,..riskf],summary)
summary(pheno[,c(ds.incident,ds.agediff),with=F])
covariates


dim(dt <- Reduce(function(x,y) merge(x = x, y = y, by = "ID_variable", sort=F), 
                 list(pheno, dt.pcs, dt.prs)))
dt[,(ds.incident):=lapply(.SD, as.numeric), .SDcol=ds.incident]


fml.cox1 <- list(as.formula(paste0(paste0('Surv(', ds.agediff,',',ds.incident,')~'),
                                   paste0(c( paste0('scale(',prs,')'), strat_sex), collapse = "+")  )) )
names(fml.cox1) <- prs


fml.cox2 <- c(
  sapply(num, function(x){
    as.formula(paste0(paste0('Surv(', ds.agediff,',',ds.incident,')~'), paste0(c(paste0('scale(',x,')'),strat_sex), collapse = "+")  ) )   
  }) ,
  sapply(setdiff(fac, c('MEN')), function(x){
    as.formula(paste0(paste0('Surv(', ds.agediff,',',ds.incident,')~'), paste0(c(x,strat_sex), collapse = "+")  ) )   
  })
)


fml.cox3 <- list(as.formula(paste0('Surv(', ds.agediff,',',ds.incident,')~', 
                                   paste0(c( paste0('scale(',num,')', collapse = ' + '), 
                                             paste0(setdiff(fac, 'MEN'), collapse = ' + '),
                                             strat_sex), collapse = "+")  )) )  
names(fml.cox3) <- 'all_crf'


fml.cox4 <- list(as.formula(paste0('Surv(', ds.agediff,',',ds.incident,')~', 
                                   paste0(c( paste0('scale(',prs,')'), 
                                             paste0('scale(',num,')', collapse = ' + '), 
                                             paste0(setdiff(fac, 'MEN'), collapse = ' + '),
                                             strat_sex), collapse = "+")  )) )  
names(fml.cox4) <- 'PRS_and_all_crf'


fml.cox5 <- list(as.formula(paste0('Surv(', ds.agediff,',',ds.incident,')~', 
                                   paste0(c( paste0('scale(',prs,')'), 
                                             paste0('scale(',num,')', collapse = ' + '), 
                                             paste0(setdiff(fac, 'MEN'), collapse = ' + '),
                                             strat_sex,
                                             'BATCH',
                                             paste0(paste0('PC',1:10), collapse = ' + ')), 
                                          collapse = "+")  )) )  
names(fml.cox5) <- 'PRS_and_all_crf_and_10pc_BATCH'
fml.cox5

fml.prs <-  c(fml.cox1,fml.cox2,fml.cox3, fml.cox4, fml.cox5)
writeLines(c('','','','fml.prs',''), sep = "\n")
print(fml.prs)


md.prs <- lapply(fml.prs, function(x){
  output.fit <- coxph(x, data = dt)
  output.fit$call$formula <- x
  return(output.fit)
})

# lapply(md.prs, function(x){
#   print(x$call$formula)
# })


res.prs <- lapply(md.prs, get_multivar_cox_stats)

numericcol <- c("beta", "pval", 
                "HR", "HR.confint.lower.95","HR.confint.upper.95",  
                "concordance", "concordance.se", 
                "p.waldTest", "p.likelihoodTest")

writeLines(c('','','','res.prs',''), sep = "\n")

invisible(lapply(res.prs, function(x){
  x[, (numericcol) := lapply(.SD, as.numeric), .SDcols = numericcol]
}))



# --------------------------------------------------------

dt.newvarnm <- data.table(varname = c("BL_AGE","BMI","SYSTM","DIASM","KOL","HDL","TRIG","ALKI2_FR02","Glc",
                                      "fit_preval","CAD_PGS000018","T2D_PGS000036","AD_PGS000334","prostate_PGS000662",
                                      "CURR_SMOKE","exercise","PREVAL_DIAB_GDM","PREVAL_DIAB_T2","PREVAL_MENTAL","PREVAL_STR_SAH_TIA","HOMAIR_fac","HOMAIR_fac",
                                      "MI_FAMILYHIST","DIAB_FAMILYHIST","CANC_FAMILYHIST",
                                      "all_crf","PRS_and_all_crf","PRS_and_all_crf_and_10pc"),
                          labelname = c("Baseline age", "BMI", "Systolic BP", "Diastolic BP","Total cholesterol","HDL","Triglyceride","Alcohol consumption", "Glucose",
                                        "Gut microbiome score","PRS","PRS","PRS","PRS",
                                        "Smoking","Exercise","Prevalent diabetes", "Prevalent T2D","Prevalent psychiatric disorders","Prevalent stroke","HOMA-IR (1.9,2.9]","HOMA-IR > 2.9",
                                        "Family history","Family history","Family history",
                                        "Conventional risk factors","PRS + Conventional risk factors","PRS + Conventional risk factors + 10PCs"))


# dt.cad <- readRDS('')
numericcol <- c('concordance', 'CI_low', 'CI_up')
dt.cad[, (numericcol) := lapply(.SD, as.numeric), .SDcols = numericcol]
dt.cad <- merge(dt.cad, dt.newvarnm, by='varname', all.x = T, sort = F)
dt.cad <- dt.cad[order(concordance,decreasing=F),]
dt.cad
cad.var.order <- setdiff(dt.cad$varname, 'PRS_and_all_crf_and_10pc_BATCH')
# [1] "CURR_SMOKE"      "PREVAL_DIAB_GDM" "exercise"        "HDL"            
# [5] "MI_FAMILYHIST"   "BMI"             "KOL"             "CAD_PGS000018"  
# [9] "SYSTM"           "BL_AGE"          "all_crf"         "PRS_and_all_crf"



# plot c-statistics
cad.p <- 
  ggplot(subset(dt.cad, varname %in% cad.var.order ), 
         aes(x=factor(varname, levels = cad.var.order), 
             y=concordance,
             color=factor(model,levels = c('individual_model','all_crf','PRS_and_all_crf')))) +         
  geom_point(size = 2.5, position=position_dodge(0.3)) +
  geom_errorbar(aes(ymin = CI_low, ymax = CI_up), width=0.3, size=0.8, position=position_dodge(0.3)) +
  ylab("C-statistic (95% CI)") +
  ggtitle ('CAD') +
  scale_y_continuous(limits = c(0.47, 0.85), breaks = seq(0.5, 0.9, 0.05)) +
  scale_x_discrete(labels = dt.cad$labelname, guide = guide_axis(angle = 45)) +
  theme_classic() +
  theme(aspect.ratio=1/1, 
        axis.text.x = element_text(size = 13),
        axis.title.x = element_blank(), 
        axis.text.y = element_text(size = 13),
        axis.title.y = element_text(size = 16),  
        plot.title= element_text(size = 16, face = 'bold.italic'),
        legend.position="none",
        legend.text = element_text(size=14),
        legend.title = element_text(size=14))+
  scale_color_manual(name='Model', 
                     labels=c('Individual risk factor',
                              'All conventional risk factorss',
                              'PRS + all conventional risk factors'),
                     values = c('black','darkorange2','blue')) +
  guides(color=guide_legend(title="",nrow=4, byrow=TRUE))



# dt.t2d <- readRDS('')
dt.t2d[, (numericcol) := lapply(.SD, as.numeric), .SDcols = numericcol]
dt.t2d <- merge(dt.t2d, dt.newvarnm, by='varname', all.x = T, sort = F)
dt.t2d <- dt.t2d[order(concordance,decreasing=F),]
t2d.var.order <- setdiff(dt.t2d$varname, 'PRS_and_all_crf_and_10pc_BATCH')
t2d.var.order
# [1] "CURR_SMOKE"      "exercise"        "KOL"             "DIAB_FAMILYHIST"
# [5] "T2D_PGS000036"   "SYSTM"           "BL_AGE"          "HDL"            
# [9] "TRIG"            "BMI"             "all_crf"         "PRS_and_all_crf"


# # subanalysis
# t2d.var.order <- setdiff(dt.t2d$varname, 'PRS_and_all_crf_and_10pc_BATCH')
# t2d.var.order
# # [1] "CURR_SMOKE"      "KOL"             "exercise"        "DIAB_FAMILYHIST"
# # [5] "T2D_PGS000036"   "SYSTM"           "BL_AGE"          "HDL"            
# # [9] "Glc"             "TRIG"            "BMI"             "all_crf"        
# # [13] "PRS_and_all_crf"



t2d.p <- 
  ggplot(subset(dt.t2d, varname %in% t2d.var.order ), 
         aes(x=factor(varname, levels = t2d.var.order), 
             y=concordance,
             color=factor(model,levels = c('individual_model','all_crf','PRS_and_all_crf')))) +         
  geom_point(size = 2.5, position=position_dodge(0.3)) +
  geom_errorbar(aes(ymin = CI_low, ymax = CI_up), width=0.3, size=0.8, position=position_dodge(0.3)) +
  ylab("C-statistic (95% CI)") +
  ggtitle ('T2D') +
  scale_y_continuous(limits = c(0.47, 0.85), breaks = seq(0.5, 0.9, 0.05)) +
  scale_x_discrete(labels = dt.t2d$labelname, guide = guide_axis(angle = 45)) +
  theme_classic() +
  theme(aspect.ratio=1/1, 
        axis.text.x = element_text(size = 13),
        axis.title.x = element_blank(), 
        axis.text.y = element_text(size = 13),
        axis.title.y = element_text(size = 16),  
        plot.title= element_text(size = 16, face = 'bold.italic'),
        legend.position="none",
        legend.text = element_text(size=14),
        legend.title = element_text(size=14))+
  scale_color_manual(name='Model', 
                     labels=c('Individual risk factor',
                              'All conventional risk factorss',
                              'PRS + all conventional risk factors'),
                     values = c('black','darkorange2','blue')) +
  guides(color=guide_legend(title="",nrow=4, byrow=TRUE))



# dt.ad <- readRDS('')
numericcol <- c('concordance', 'CI_low', 'CI_up')
dt.ad[, (numericcol) := lapply(.SD, as.numeric), .SDcols = numericcol]
dt.ad <- merge(dt.ad, dt.newvarnm, by='varname', all.x = T, sort = F)
dt.ad <- dt.ad[order(concordance,decreasing=F),]
ad.var.order <- setdiff(dt.ad$varname, 'PRS_and_all_crf_and_10pc_BATCH')
ad.var.order
# [1] "PREVAL_MENTAL"      "PREVAL_STR_SAH_TIA" "HDL"               
# [4] "DIASM"              "PREVAL_DIAB_T2"     "CURR_SMOKE"        
# [7] "exercise"           "KOL"                "BMI"               
# [10] "ALKI2_FR02"         "AD_PGS000334"       "SYSTM"             
# [13] "BL_AGE"             "all_crf"            "PRS_and_all_crf" 


ad.p <- 
  ggplot(subset(dt.ad, varname %in% ad.var.order ), 
         aes(x=factor(varname, levels = ad.var.order), 
             y=concordance,
             color=factor(model,levels = c('individual_model','all_crf','PRS_and_all_crf')))) +         
  geom_point(size = 2.5, position=position_dodge(0.3)) +
  geom_errorbar(aes(ymin = CI_low, ymax = CI_up), width=0.3, size=0.8, position=position_dodge(0.3)) +
  ylab("C-statistic (95% CI)") +
  ggtitle ("Alzheimer's disease") +
  scale_y_continuous(limits = c(0.47, 0.92), breaks = seq(0.5, 0.9, 0.05)) +
  scale_x_discrete(labels = dt.ad$labelname, guide = guide_axis(angle = 45)) +
  theme_classic() +
  theme(aspect.ratio=1/1, 
        axis.text.x = element_text(size = 12),
        axis.title.x = element_blank(), 
        axis.text.y = element_text(size = 13),
        axis.title.y = element_text(size = 16),  
        plot.title= element_text(size = 16, face = 'bold.italic'),
        legend.position="none",
        legend.text = element_text(size=14),
        legend.title = element_text(size=14))+
  scale_color_manual(name='Model', 
                     labels=c('Individual risk factor',
                              'All conventional risk factorss',
                              'PRS + all conventional risk factors'),
                     values = c('black','darkorange2','blue')) +
  guides(color=guide_legend(title="",nrow=4, byrow=TRUE))




# dt.pc <- readRDS('')
numericcol <- c('concordance', 'CI_low', 'CI_up')
dt.pc[, (numericcol) := lapply(.SD, as.numeric), .SDcols = numericcol]
dt.pc <- merge(dt.pc, dt.newvarnm, by='varname', all.x = T, sort = F)
dt.pc <- dt.pc[order(concordance,decreasing=F),]
pc.var.order <- setdiff(dt.pc$varname, 'PRS_and_all_crf_and_10pc_BATCH')
pc.var.order
# [1] "ALKI2_FR02"         "exercise"           "BMI"               
# [4] "CANC_FAMILYHIST"    "CURR_SMOKE"         "prostate_PGS000662"
# [7] "BL_AGE"             "all_crf"            "PRS_and_all_crf"


pc.p <- 
  ggplot(subset(dt.pc, varname %in% pc.var.order ), 
         aes(x=factor(varname, levels = pc.var.order), 
             y=concordance,
             color=factor(model,levels = c('individual_model','all_crf','PRS_and_all_crf')))) +  
  geom_point(size = 2.5, position=position_dodge(0.3)) +
  geom_errorbar(aes(ymin = CI_low, ymax = CI_up), width=0.3, size=0.8, position=position_dodge(0.3)) +
  ylab("C-statistic (95% CI)") +
  ggtitle ('Prostate cancer') +
  scale_y_continuous(limits = c(0.46, 0.85), breaks = seq(0.5, 0.9, 0.05)) +
  scale_x_discrete(labels = dt.pc$labelname, guide = guide_axis(angle = 45)) +
  theme_classic() +
  theme(aspect.ratio=1/1, 
        axis.text.x = element_text(size = 13),
        axis.title.x = element_blank(), 
        axis.text.y = element_text(size = 13),
        axis.title.y = element_text(size = 16),  
        plot.title= element_text(size = 16, face = 'bold.italic'),
        legend.position="bottom",
        # legend.position="",
        legend.text = element_text(size=12),
        legend.title = element_text(size=12))   +
  scale_color_manual(name='Model', 
                     labels=c('Individual risk factor',
                              'All conventional risk factorss',
                              'PRS + all conventional risk factors'),
                     values = c('black','darkorange2','blue')) +
  guides(color=guide_legend(title="Model",nrow=4, byrow=TRUE))


p <- plot_grid(cad.p, 
               t2d.p, 
               ad.p, 
               pc.p, 
               ncol = 2, align = "h",axis = "bt",
               labels = c("A", "B","C","D"), label_size = 20,
               label_x = 0.01,  hjust = -2, vjust = 1)




# plot HR


dt.newvarnm <- data.table(var = c(paste0(paste0('scale(',c("BL_AGE","BMI","SYSTM","DIASM","KOL","HDL","TRIG","ALKI2_FR02",
                                                           "fit_preval","CAD_PGS000018","T2D_PGS000036","AD_PGS000334","prostate_PGS000662","Glc")),')'),
                                  paste0(c("CURR_SMOKE","exercise","PREVAL_DIAB_GDM","PREVAL_DIAB_T2","PREVAL_MENTAL","PREVAL_STR_SAH_TIA",
                                           "MI_FAMILYHIST","DIAB_FAMILYHIST","CANC_FAMILYHIST"),1)), 
                          labelname = c(paste0(c("Baseline age", "BMI", "Systolic BP", "Diastolic BP","Total cholesterol","HDL","Triglyceride","Alcohol consumption",
                                                 "Gut microbiome score", "PRS","PRS","PRS","PRS","Glucose"),' per s.d.'),
                                        paste0(c("Smoking","Exercise","Prevalent diabetes", "Prevalent T2D","Prevalent psychiatric disorders","Prevalent stroke",
                                                 "Family history","Family history","Family history"), " (yes/no)")))



# dt.cad <- readRDS('')
dt.cad <- dt.cad[order(HR,decreasing=F),]
cad.var.order <- dt.cad$var
cad.p <-
  ggplot(dt.cad, 
         aes(x=factor(var, levels = cad.var.order), 
             y=HR) ) +
  geom_point(size = 2.5, position=position_dodge(0.3)) +
  geom_errorbar(aes(ymin = HR.confint.lower.95, ymax = HR.confint.upper.95), width=0.3, size=0.8,  position=position_dodge(0.3)) +
  geom_hline(yintercept = 1, linetype='dashed', col = 'red', size=0.8) +
  ylab("HR (95% CI)") +
  ggtitle ('CAD') +
  scale_y_continuous(limits = c(0, ceiling(dt.cad[, max(HR.confint.upper.95)])), 
                     breaks = seq(0, ceiling(dt.cad[, max(HR.confint.upper.95)]), 1)) +
  scale_x_discrete(labels = dt.cad$labelname, guide = guide_axis(angle = 45)) +
  theme_classic() +
  theme(aspect.ratio=1/1, 
        plot.title = element_text(size=16, face="bold.italic"),
        legend.position="none",
        legend.title = element_text(size=14), 
        legend.text=element_text(size=14), 
        axis.text.x = element_text(size = 14), 
        axis.title.x = element_blank(), 
        axis.text.y = element_text(size = 14),
        axis.title.y = element_text(size = 16)) 



# dt.t2d<- readRDS('')
dt.t2d <- dt.t2d[order(HR,decreasing=F),]
t2d.var.order <- dt.t2d$var
t2d.p <-
  ggplot(dt.t2d, 
         aes(x=factor(var, levels = t2d.var.order), 
             y=HR) ) +
  geom_point(size = 2.5, position=position_dodge(0.3)) +
  geom_errorbar(aes(ymin = HR.confint.lower.95, ymax = HR.confint.upper.95), width=0.3, size=0.8, position=position_dodge(0.3)) +
  geom_hline(yintercept = 1, linetype='dashed', col = 'red', size=0.8) +
  ylab("HR (95% CI)") +
  ggtitle ('T2D') +
  scale_y_continuous(limits = c(0, 2.1  ), 
                     breaks = seq(0, 2, 1)) +
  scale_x_discrete(labels = dt.t2d$labelname, guide = guide_axis(angle = 45)) +
  theme_classic() +
  theme(aspect.ratio=1/1, 
        plot.title = element_text(size=16, face="bold.italic"),
        legend.position="none",
        legend.title = element_text(size=14), 
        legend.text=element_text(size=14), 
        axis.text.x = element_text(size = 14), 
        axis.title.x = element_blank(), 
        axis.text.y = element_text(size = 14),
        axis.title.y = element_text(size = 16))


# dt.ad<-readRDS('')
dt.ad <- dt.ad[order(HR,decreasing=F),]
ad.var.order <- dt.ad$var
ad.p <-
  ggplot(dt.ad,
         aes(x=factor(var, levels = ad.var.order),
             y=HR) ) +
  geom_point(size = 2.5, position=position_dodge(0.3)) +
  geom_errorbar(aes(ymin = HR.confint.lower.95, ymax = HR.confint.upper.95), width=0.3, size=0.8, position=position_dodge(0.3)) +
  geom_hline(yintercept = 1, linetype='dashed', col = 'red', size=0.8) +
  ylab("HR (95% CI)") +
  ggtitle ("Alzheimer's disease") +
  scale_y_continuous(limits = c(0, ceiling(dt.ad[, max(HR.confint.upper.95)]))  ,
                     breaks = seq(1, ceiling(dt.ad[, max(HR.confint.upper.95)]), 2)) +
  scale_x_discrete(labels = dt.ad$labelname, guide = guide_axis(angle = 45)) +
  theme_classic() +
  theme(aspect.ratio=1/1,
        plot.title = element_text(size=16, face="bold.italic"),
        legend.position="none",
        legend.title = element_text(size=14),
        legend.text=element_text(size=14),
        axis.text.x = element_text(size = 12),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 14),
        axis.title.y = element_text(size = 16))


# dt.pc <-readRDS('')
dt.pc <- dt.pc[order(HR,decreasing=F),]
pc.var.order <- dt.pc$var

pc.p <-
  ggplot(dt.pc, 
         aes(x=factor(var, levels = pc.var.order), 
             y=HR) ) +
  geom_point(size = 2.5, position=position_dodge(0.3)) +
  geom_errorbar(aes(ymin = HR.confint.lower.95, ymax = HR.confint.upper.95), width=0.3, size=0.8,  position=position_dodge(0.3)) +
  geom_hline(yintercept = 1, linetype='dashed', col = 'red', size=0.8) +
  ylab("HR (95% CI)") +
  ggtitle ('Prostate cancer') +
  scale_y_continuous(limits = c(0, ceiling(dt.pc[, max(HR.confint.upper.95)])), 
                     breaks = seq(0, ceiling(dt.pc[, max(HR.confint.upper.95)]), 1)) +
  scale_x_discrete(labels = dt.pc$labelname, guide = guide_axis(angle = 45)) +
  theme_classic() +
  theme(aspect.ratio=1/1, 
        plot.title = element_text(size=16, face="bold.italic"),
        legend.position="none",
        legend.title = element_text(size=14), 
        legend.text=element_text(size=14), 
        axis.text.x = element_text(size = 14), 
        axis.title.x = element_blank(), 
        axis.text.y = element_text(size = 14),
        axis.title.y = element_text(size = 16)) 






#########################################################################
# 
# PRS and gut microbiome score association
# 
#########################################################################


# CHD
ds <- 'CHD'
ds.incident <- paste0('INCIDENT_',ds)
ds.agediff <- paste0(ds, '_AGEDIFF')
ds;ds.incident;ds.agediff
prs <- 'CAD_PGS000018'
riskf <- c('MEN','BL_AGE','BMI', 'SYSTM','KOL','HDL','CURR_SMOKE', 'exercise','PREVAL_DIAB_GDM','MI_FAMILYHIST')
num <-  c('BL_AGE','BMI','SYSTM', 'KOL', 'HDL')
fac <-  c('MEN','CURR_SMOKE','exercise','PREVAL_DIAB_GDM','MI_FAMILYHIST')
covariates <- setdiff(riskf, c('MEN'))
covariates
strat_sex <- 'strata(MEN)'




# T2D
prs <- 'T2D_PGS000036'
ds <- 'DIAB_T2'
ds.incident <- paste0('INCIDENT_',ds)
ds.agediff <- paste0(ds, '_AGEDIFF')
riskf <- c('MEN','BL_AGE','BMI','SYSTM','KOL','HDL','TRIG','CURR_SMOKE','exercise', 'DIAB_FAMILYHIST')
num <-  c('BL_AGE','BMI','SYSTM','KOL','HDL','TRIG')
fac <-  c('MEN','CURR_SMOKE','exercise','DIAB_FAMILYHIST')
covariates <- setdiff(riskf, c('MEN'))
covariates
strat_sex <- 'strata(MEN)'




# AD
prs <- 'AD_PGS000334'
ds <- 'AD'; ds.incident <- paste0('INCIDENT_',ds); ds.agediff <- paste0(ds, '_AGEDIFF')
riskf <- c('MEN','BL_AGE','BMI','SYSTM','DIASM','KOL','HDL','ALKI2_FR02','CURR_SMOKE','exercise',
           'PREVAL_DIAB_T2','PREVAL_STR_SAH_TIA','PREVAL_MENTAL')
num <-  c('BL_AGE','BMI','SYSTM','DIASM','KOL','HDL','ALKI2_FR02')
fac <-  c('MEN','CURR_SMOKE','exercise','PREVAL_DIAB_T2','PREVAL_STR_SAH_TIA','PREVAL_MENTAL')
covariates <- setdiff(riskf, c('MEN'))
covariates
strat_sex <- 'strata(MEN)'



# PC
prs <- 'prostate_PGS000662'
ds <- 'CR_PROSTCANC'; ds.incident <- paste0('INCIDENT_',ds); ds.agediff <- paste0(ds, '_AGEDIFF')
riskf <- c('BL_AGE','BMI','ALKI2_FR02','CURR_SMOKE','exercise','CANC_FAMILYHIST')
num <-  c('BL_AGE','BMI','ALKI2_FR02')
fac <-  c('CURR_SMOKE','exercise','CANC_FAMILYHIST')
covariates <- setdiff(riskf, c('MEN'))
covariates
strat_sex <- NULL
adjust_sex <- NULL



unlink('.RData')
rm(list = ls())
pcks <- c('data.table','survival','caret','glmnet','ggplot2','gridExtra','grid','cowplot')
sapply(pcks, require, character.only = TRUE)
sapply(pcks, packageVersion)




get_multivar_cox_stats <- function(cmd){
  
  x <- summary(cmd)
  
  res <- data.table(x$coef, keep.rownames = 'var')
  res[,pval:=formatC(`Pr(>|z|)`, digits=2, format = "e")]
  setnames(res, old = c('coef','exp(coef)'), new = c('beta','HR'))
  conf.int <- data.table(x$conf.int, keep.rownames = 'var')
  res <- merge(res, conf.int, sort=F)
  setnames(res, old = c('lower .95' , 'upper .95'), new = c('HR.confint.lower.95' , 'HR.confint.upper.95'))
  res[,`HR (95% CI)`:=paste0(round(HR, digits=2), " (", round(HR.confint.lower.95, digits=2), "-", round(HR.confint.upper.95, digits=2), ")")]
  rmcol <- c('se(coef)','Pr(>|z|)','z','exp(coef)','exp(-coef)')
  res[,(rmcol):=NULL]
  
  
  p.wald <- formatC(x$wald["pvalue"], digits=2, format = "e") 
  p.likelihood <- formatC(x$logtest["pvalue"], digits=2, format = "e")  
  concordance <- x$concordance[1]
  concordance.se <- x$concordance[2]
  res[,p.waldTest:=p.wald]
  res[,p.likelihoodTest:=p.likelihood]
  res[,concordance:=concordance]
  res[,concordance.se:=concordance.se]
  
  return(res)
}


dt.prs <- fread(fn.PRSice_all_score)
dim(clrdata <- fread('...')) 
dim(pheno <- readRDS('...')) 
pheno[,(fac):=lapply(.SD, as.factor), .SDcols=fac]
length(foldids <- readRDS('...'))  



dim(ds.data <- merge(pheno[,c('Barcode', ds.incident),with=F], 
                     clrdata, 
                     by='Barcode',
                     sort=F)) 
system.time(
  ridge_cvmod <- 
    lapply( names(foldids) , function(i){
      ridge_mod = cv.glmnet(x = data.matrix(ds.data[,!c('Barcode',ds.incident), with=F]), 
                            y = ds.data[,get(ds.incident)], 
                            alpha = 0,        
                            lambda = 10^seq(-4, 2, length = 200), 
                            family = "binomial", 
                            type.measure = "auc", 
                            intercept = T, 
                            foldid = foldids[[i]],
                            keep=TRUE, 
                            # parallel = T, 
                            standardize = FALSE) 
      return(ridge_mod)
    })
)
names(ridge_cvmod) <- names(foldids); names(ridge_cvmod)


# ridge_cvm_lambda_min <-
#   sapply(names(ridge_cvmod), USE.NAMES = TRUE, function(x){
# 
#     ridge_mod <- ridge_cvmod[[x]]
#     cvm.res <- data.table(cvm = ridge_mod$cvm,
#                           lambda = ridge_mod$lambda)
#     n.min <- which(ridge_mod$lambda==ridge_mod$lambda.min)
# 
#     return(data.table(cvm.res[n.min], n.min=n.min))
#   })
# names(which.max(ridge_cvm_lambda_min['cvm',]))



ridge_cvmod_best <- ridge_cvmod[[optimal_run]]   
ridge_cvmod_best_pred <- data.table(ds.data[,c('Barcode',ds.incident),with=F], 
                                    fit_preval = ridge_cvmod_best$fit.preval[,n.min], 
                                    foldid = ridge_cvmod_best$foldid)


dim(dt <- merge(ridge_cvmod_best_pred, pheno[,-ds.incident,with=F], by='Barcode', sort=F))
dt[,(ds.incident):=lapply(.SD, as.numeric), .SDcol=ds.incident]


fml.microbiome <- list(
  
  as.formula(paste0('Surv(', ds.agediff,',',ds.incident,')~ scale(fit_preval)')) ,
  
  as.formula(paste0(paste0('Surv(', ds.agediff,',',ds.incident,')~'),
                    paste0(c( 'scale(fit_preval)',strat_sex), collapse = "+")  ))   ,
  
  as.formula(paste0(paste0('Surv(', ds.agediff,',',ds.incident,')~'),
                    paste0(c( 'scale(fit_preval)','scale(BL_AGE)',strat_sex), collapse = "+")  ))  ,
  
  as.formula(paste0('Surv(', ds.agediff,',',ds.incident,')~', 
                    paste0(c( paste0('scale(',num,')', collapse = ' + '), 
                              paste0(setdiff(fac, 'MEN'), collapse = ' + '),
                              strat_sex), collapse = "+")  )), 
  
  as.formula(paste0(paste0('Surv(', ds.agediff,',',ds.incident,')~'), 
                    paste0(c( 'scale(fit_preval)',
                              paste0('scale(',num,')', collapse = ' + '), 
                              paste0(setdiff(fac, 'MEN'), collapse = ' + '),
                              strat_sex), collapse = "+")  ))
)

names(fml.microbiome) <- c('GMscore', 
                           'GMscore_strat_sex', 
                           'GMscore_strat_sex_adjust_age', 
                           'crf_strat_sex', 
                           'GMscore_crf_strat_sex')

md.microbiome <- lapply(fml.microbiome, function(this.fml){
  output.fit <- coxph(this.fml, data = dt)
  output.fit$call$formula <- this.fml
  return(output.fit)
})
res.microbiome <- lapply(md.microbiome, get_multivar_cox_stats)





dim(dt <- merge(ridge_cvmod_best_pred, pheno[,-ds.incident,with=F], by='Barcode', sort=F))
dim(dt <- merge(dt, dt.prs, by='FID', sort=F))
dt[,(ds.incident):=lapply(.SD, as.numeric), .SDcol=ds.incident]



fml.integrative <- list(
  
  as.formula(paste0(paste0('Surv(', ds.agediff,',',ds.incident,')~'),
                    paste0(c( 'scale(BL_AGE)',strat_sex), collapse = "+")  ) )   ,
  
  as.formula(paste0(paste0('Surv(', ds.agediff,',',ds.incident,')~'),
                    paste0(c(paste0('scale(',prs,')'), strat_sex), collapse = "+")  ))   ,
  
  # 3.	cox(ds, time on study) ~ scale(microbiome score) + strat(sex)
  as.formula(paste0(paste0('Surv(', ds.agediff,',',ds.incident,')~'),
                    paste0(c( 'scale(fit_preval)',strat_sex), collapse = "+")  ))   ,
  
  as.formula(paste0(paste0('Surv(', ds.agediff,',',ds.incident,')~'),
                    paste0(c(paste0('scale(',prs,')'), 'scale(BL_AGE)', strat_sex), collapse = "+")  ))   ,
  
  as.formula(paste0(paste0('Surv(', ds.agediff,',',ds.incident,')~'),
                    paste0(c( 'scale(fit_preval)','scale(BL_AGE)', strat_sex), collapse = "+")  ))   ,
  
  as.formula(paste0(paste0('Surv(', ds.agediff,',',ds.incident,')~'),
                    paste0(c(paste0('scale(',prs,')'), 'scale(fit_preval)','scale(BL_AGE)', strat_sex), collapse = "+")  ))   ,
  
  as.formula(paste0('Surv(', ds.agediff,',',ds.incident,')~', 
                    paste0(c( paste0('scale(',num,')', collapse = ' + '), 
                              paste0(setdiff(fac, 'MEN'), collapse = ' + '),
                              strat_sex), collapse = "+")  )) ,
  
  as.formula(paste0(paste0('Surv(', ds.agediff,',',ds.incident,')~'), 
                    paste0(c( paste0('scale(',prs,')'), 'scale(fit_preval)',
                              paste0('scale(',num,')', collapse = ' + '), 
                              paste0(setdiff(fac, 'MEN'), collapse = ' + '),
                              strat_sex), collapse = "+")  ))
)
names(fml.integrative) <- c('age_stratSex', 'prs_stratSex', 'GMscore_stratSex', 
                            'age_stratSex_prs','age_stratSex_GMscore', 'age_stratSex_prs_GMscore',
                            'crf_stratSex','crf_stratSex_prs_GMscore')


md.integrative <- lapply(fml.integrative, function(this.fml){
  output.fit <- coxph(this.fml, data = dt)
  output.fit$call$formula <- this.fml
  return(output.fit)
})

res.integrative <- lapply(md.integrative, get_multivar_cox_stats)
res2excel <- rbindlist(res.integrative, idcol = "model")
res2excel[, Disease:=ds]
res2excel

dt.newmodelnm <- data.table(model=c('age_stratSex','prs_stratSex','GMscore_stratSex',
                                    'age_stratSex_prs','age_stratSex_GMscore','age_stratSex_prs_GMscore',
                                    'crf_stratSex','crf_stratSex_prs_GMscore'),
                            modelname=c('Age','PRS','Gut microbiome score',
                                        'Age+PRS','Age+Gut microbiome score','Age+PRS+Gut microbiome score',
                                        'Conventional risk factors','Conventional risk factors+PRS+Gut microbiome score'))



best.rocs <- 
  lapply(ridge_cvmod, function(cvmd){
    rocs <- roc.glmnet(cvmd$fit.preval, newy = y)
    best <- which.max(cvmd$cvm)
    return(rocs[[best]])
  })
plot(best.rocs[[1]], type = "l", main='', lwd = 1, col="grey")
for (i in 2:10){
  lines(best.rocs[[i]], lwd = 1,col = "grey")
}
lines(x = c(0,1), y = c(0,1), col='black', lty='dashed', lwd=2)
lines(best.rocs[[optimal_run]], lwd = 2,col = "red")











