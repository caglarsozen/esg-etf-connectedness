for(f in sort(list.files('R',pattern='[.]R$',full.names=TRUE)))source(f)
pr<-readRDS('outputs/PRESPECIFIED_PROTOCOL.rds');pr$root<-getwd()
pan<-revision_read('input/archived_returns.csv',pr$development_end)
model<-readLines('outputs/FROZEN_COVARIANCE_MODEL.txt',warn=FALSE)[1]
b<-revision_bundle(pan$R,pan$dates,nrow(pan$R),pr$models[match(model,pr$models$id),],pr)
cfg<-esg_config('pilot',getwd());cfg$bootstrap_B<-pr$bootstrap_B;cfg$var_min_obs<-40L
path<-'outputs/network_development';dir.create(path,showWarnings=FALSE)
tables<-esg_conn_tables(b$var,b$hmm$labels,colnames(pan$R),cfg$blocks,c(1,5,10,20))
for(nm in names(tables))esg_csv(tables[[nm]],file.path(path,paste0(nm,'.csv')))
esg_csv(data.frame(source=b$hmm$labels,b$hmm$P),file.path(path,'transition_matrix.csv'))
esg_csv(data.frame(date=b$dates[b$cv$valid],b$hmm$filtered),file.path(path,'regime_probabilities.csv'))
for(bl in pr$block_lengths){
  cat('Network bootstrap',bl,'\n')
  z<-esg_conn_bootstrap(b$packet,b$var,b$hmm$labels,colnames(pan$R),cfg$blocks,cfg,
    'development_observation_driven',block=bl)
  z$table$scope<-'conditional on fitted observation-driven variance path and state labels'
  esg_csv(z$table,file.path(path,paste0('bootstrap_',bl,'.csv')))
  esg_csv(z$role_switch,file.path(path,paste0('CE_role_switch_',bl,'.csv')))
}
# Diagnostic counterfactuals. Neither identifies a structural causal network.
counter<-list()
for(s in seq_along(b$var))for(h in c(1L,10L))for(kind in c('original','diagonal_innovations','diagonal_dynamics')){
  f<-b$var[[s]]
  if(kind=='diagonal_innovations')f$Sigma<-diag(diag(f$Sigma))
  if(kind=='diagonal_dynamics')f$A<-lapply(f$A,function(A)diag(diag(A)))
  counter[[length(counter)+1]]<-data.frame(regime=b$hmm$labels[s],horizon=h,counterfactual=kind,
    TCI=esg_conn(esg_fevd(f,h),colnames(pan$R))$TCI)
}
esg_csv(do.call(rbind,counter),file.path(path,'dependence_counterfactuals.csv'))
writeLines(c('Observation-driven log conditional variances, NOT a validated latent Factor-SV network.',
 'Development sample only. Descriptive fitted regimes; parameter estimation uses the whole development prefix.',
 'Bootstrap conditional on extracted variance path and state classification.',
 'H10 minus H1 is a horizon change, not an additive or identified causal dynamic share.',
 'Within-regime residual Ljung-Box values in diagnostics.csv are compressed-row diagnostics, not daily-time tests.'),
 file.path(path,'SCOPE.txt'))
