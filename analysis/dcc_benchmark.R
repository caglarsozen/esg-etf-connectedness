for(f in sort(list.files('R',pattern='[.]R$',full.names=TRUE)))source(f)
source('benchmarks/gaussian_dcc.R')
pr<-readRDS('outputs/PRESPECIFIED_PROTOCOL.rds');pr$root<-getwd()
pan<-revision_read('input/archived_returns.csv');R<-pan$R;dates<-pan$dates
out<-'outputs/dcc_benchmark';dir.create(out,showWarnings=FALSE)
ends<-esg_month_ends(dates);origins<-ends[dates[ends]>=pr$development_end&ends<nrow(R)]
W<-array(NA_real_,c(length(origins),ncol(R),1),dimnames=list(as.character(dates[origins]),colnames(R),'DCC_Gaussian_MinVar'))
rows<-list()
for(j in seq_along(origins)){
  t<-origins[j];cat('Gaussian DCC',as.character(dates[t]),'\n')
  loc<-file.path(out,paste0('fit_',dates[t],'.rds'))
  identity<-esg_hash(list(R=R[1:t,,drop=FALSE],source=unname(tools::md5sum('benchmarks/gaussian_dcc.R'))))
  if(file.exists(loc)){fit<-readRDS(loc);esg_assert(identical(fit$identity,identity),'DCC cache mismatch.')}else{
    fit<-revision_dcc(R[1:t,,drop=FALSE]);fit$identity<-identity;esg_save(fit,loc)}
  W[j,,1]<-esg_minvar(fit$forecast,pr$cap)
  rows[[j]]<-data.frame(origin=dates[t],a=fit$a,b=fit$b,convergence=fit$convergence,
    max_GARCH_persistence=max(fit$univariate$alpha+fit$univariate$beta))
}
esg_csv(do.call(rbind,rows),file.path(out,'convergence.csv'))
for(cost in c(0,5,10,25)){
  bt<-esg_backtest_targets(R,dates,origins,W,cost)
  esg_csv(data.frame(date=bt$dates,bt$returns,check.names=FALSE),file.path(out,paste0('returns_cost',cost,'.csv')))
  esg_csv(esg_performance_table(bt,cost),file.path(out,paste0('performance_cost',cost,'.csv')))
}
esg_save(list(targets=W,origins=origins),file.path(out,'TARGETS.rds'))
