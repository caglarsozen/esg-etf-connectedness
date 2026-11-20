for(f in sort(list.files('R',pattern='[.]R$',full.names=TRUE)))source(f)
folders<-list.dirs('outputs/fsv',recursive=FALSE,full.names=TRUE);rows<-list()
for(folder in folders){
  gatepath<-file.path(folder,'MCMC_GATE.csv');if(!file.exists(gatepath))next
  gate<-read.csv(gatepath,check.names=FALSE);p<-file.path(folder,if(gate$passed)'QUALIFIED_PREFLIGHT.rds' else 'REJECTED_PREFLIGHT.rds')
  z<-readRDS(p);d<-revision_read('input/archived_returns.csv',as.Date(gate$cutoff))
  E<-matrix(NA_real_,nrow(d$R),ncol(d$R),dimnames=dimnames(d$R))
  for(t in z$valid)E[t,]<-as.numeric(d$R[t,,drop=FALSE]%*%esg_root(z$H[t,,],TRUE))
  rr<-revision_residual_checks(list(eps=E),d$dates,min(d$dates),max(d$dates))
  rr$scope<-'smoothed conditional posterior-mean whitening; not calibrated predictive PIT'
  esg_csv(rr,file.path(folder,'CONDITIONAL_RESIDUAL_CHECKS.csv'))
  gate$squared_dependence_rejections_Holm<-sum(rr$squared_Holm_p<.05)
  gate$scientific_status<-if(!gate$passed)'REJECTED_MCMC' else
    if(any(rr$squared_Holm_p<.05))'MCMC_PASS_MODEL_ADEQUACY_UNRESOLVED' else 'PREFLIGHT_ONLY_PREDICTIVE_VALIDATION_REQUIRED'
  rows[[length(rows)+1]]<-gate
}
if(length(rows))esg_csv(do.call(rbind,rows),'outputs/FSV_EVIDENCE_SUMMARY.csv')
