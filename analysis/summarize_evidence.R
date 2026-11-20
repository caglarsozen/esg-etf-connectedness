for(f in sort(list.files('R',pattern='[.]R$',full.names=TRUE)))source(f)
out<-'outputs/evidence';dir.create(out,showWarnings=FALSE)
all<-inc<-stress_increment<-list()
for(bl in c(5,10,20)){
 p<-paste0('outputs/matched_market_ablation',if(bl==10)'' else paste0('_b',bl))
 all[[as.character(bl)]]<-read.csv(file.path(p,'MATCHED_TCI_DIFFERENCES.csv'))
 inc[[as.character(bl)]]<-read.csv(file.path(p,'EXPLORATORY_HORIZON_INCREMENT.csv'))
 boot<-readRDS(file.path(p,'BOOTSTRAP_DRAWS.rds'))
  stopifnot(identical(boot$rng_method,'single_per_block_stream'),nrow(boot$boot)==9999L)
  # Added 2026-11-20: report two-sided evidence for these exploratory level contrasts.
  # Preserve the original one-sided columns, draws and all null results.
  aa<-all[[as.character(bl)]]
  aa$p_two_sided<-vapply(aa$contrast,function(nm)
    unname(esg_basic_inference(boot$observed[nm],boot$boot[,nm],'two.sided')['p_value']),numeric(1))
  aa$p_holm_six_two_sided<-p.adjust(aa$p_two_sided,'holm')
  all[[as.character(bl)]]<-aa
 stress_increment[[as.character(bl)]]<-data.frame(block=bl,t(esg_basic_inference(boot$observed['H10_Stress']-boot$observed['H10_Calm'],boot$boot[,'H10_Stress']-boot$boot[,'H10_Calm'],'two.sided')),scope='post hoc stress-minus-calm contrast of exposure-control reduction',check.names=FALSE)
}
a<-do.call(rbind,all);b<-do.call(rbind,inc);a$p_holm_18<-p.adjust(a$p_value,'holm');b$p_holm_nine<-p.adjust(b$p_value,'holm')
joint<-p.adjust(c(a$p_value,b$p_value),'holm');a$p_holm_all27<-head(joint,nrow(a));b$p_holm_all27<-tail(joint,nrow(b))
cst<-do.call(rbind,stress_increment);cst$p_holm_three<-p.adjust(cst$p_value,'holm')
j30<-p.adjust(c(a$p_value,b$p_value,cst$p_value),'holm');a$p_holm_all30<-j30[seq_len(nrow(a))];b$p_holm_all30<-j30[nrow(a)+seq_len(nrow(b))];cst$p_holm_all30<-tail(j30,nrow(cst))
a$p_holm_18_two_sided<-p.adjust(a$p_two_sided,'holm')
j30two<-p.adjust(c(a$p_two_sided,b$p_value,cst$p_value),'holm')
a$p_holm_all30_two_sided<-j30two[seq_len(nrow(a))]
b$p_holm_all30_two_sided<-j30two[nrow(a)+seq_len(nrow(b))]
cst$p_holm_all30_two_sided<-tail(j30two,nrow(cst))
esg_csv(cst,file.path(out,'POSTHOC_STRESS_AMPLIFICATION.csv'))
esg_csv(a,file.path(out,'MATCHED_MARKET_BLOCK_SENSITIVITY.csv'));esg_csv(b,file.path(out,'EXPLORATORY_HORIZON_INCREMENT.csv'))
gates<-networks<-ppcs<-pilots<-list()
for(cut in c('20230228','20230929')){
 p<-file.path('outputs/fsv',paste0('anchored_dynamic_',cut,'_160000_40000_40'))
 gp<-file.path(p,'MCMC_GATE.csv');if(!file.exists(gp))next
 g<-read.csv(gp);gates[[cut]]<-g
 files<-c(network='network_uncertainty/POSTERIOR_PROPAGATED_SUMMARY.csv',ppc='POSTERIOR_PREDICTIVE_SUMMARY.csv',pilot='predictive_pilot/SUMMARY.csv')
 for(nm in names(files))if(file.exists(file.path(p,files[[nm]]))){
  z<-read.csv(file.path(p,files[[nm]]));z$cutoff<-g$cutoff
  if(nm=='network')networks[[cut]]<-z else if(nm=='ppc')ppcs[[cut]]<-z else pilots[[cut]]<-z
 }
}
for(nm in c('gates','networks','ppcs','pilots')){
 x<-get(nm);fname<-switch(nm,gates='LONG_ANCHORED_GATES.csv',networks='POSTERIOR_NETWORK_TWO_DATES.csv',ppcs='CONDITIONAL_PPC_TWO_DATES.csv',pilots='PREDICTIVE_PILOTS_TWO_DATES.csv')
 if(length(x))esg_csv(do.call(rbind,x),file.path(out,fname))
}
fsv_status<-if(length(gates)<2L)'INCOMPLETE_ARCHIVE' else if(all(vapply(
 gates,function(g)isTRUE(g$passed),logical(1))))'TWO_DIAGNOSTIC_GATES_PASSED' else 'MCMC_GATE_REJECTION'
ledger<-data.frame(experiment=c('Matched SPY exposure','Long FSV','Posterior network','Predictive pilot','Holding horizon','Earlier network sensitivities'),
 status=c('EXPLORATORY_CONDITIONAL',fsv_status,
 if(length(networks)==2L)'FIRST_STAGE_PATH_PROPAGATION' else 'INCOMPLETE_ARCHIVE',
 if(length(pilots)==2L)'SHORT_PREDICTIVE_PILOT' else 'INCOMPLETE_ARCHIVE',
 'NO_INCREMENTAL_ES_GAIN','ARCHIVED_LEGACY_BOOTSTRAP'),
 restriction=c('fitted-path conditional inference','two overlapping cutoffs and selected covariance dates',
 'conditional on plug-in HMM and VAR fits','21 days per cutoff and particle ancestor depletion',
 'four fixed arms compared with C','old seeding and independently refitted regimes'))
esg_csv(ledger,file.path(out,'EXPERIMENT_LEDGER.csv'))
print(a[a$block==10,c('contrast','estimate','lower','upper','p_holm_18_two_sided','p_holm_all30_two_sided')])
