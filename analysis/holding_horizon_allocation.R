for(f in sort(list.files('R',pattern='[.]R$',full.names=TRUE)))source(f)
source('experiments/holding_horizon.R')
pr<-revision_protocol(getwd());pan<-revision_read('input/archived_returns.csv',pr$development_end)
old<-readRDS('outputs/DEVELOPMENT_targets.rds');origins<-old$origins;spec<-old$model
out<-'outputs/holding_horizon';dir.create(out,showWarnings=FALSE)
id<-esg_hash(list(pan=pan,origins=origins,spec=spec,source=tools::md5sum(c('experiments/holding_horizon.R','analysis/holding_horizon_allocation.R'))))
cache<-file.path(out,'CHECKPOINT.rds');Ws<-list();scores<-list()
if(file.exists(cache)){saved<-readRDS(cache);stopifnot(identical(saved$id,id));Ws<-saved$Ws;scores<-saved$scores}
for(i in seq_along(origins))if(i>length(Ws)){
 t<-origins[i];cat('Holding-horizon origin',as.character(pan$dates[t]),'\n')
 b<-revision_bundle(pan$R,pan$dates,t,spec,pr);hs<-horizon_bundle_forecasts(b)
 W<-matrix(NA_real_,ncol(pan$R),5,dimnames=list(colnames(pan$R),c('C',names(hs))))
 W[,'C']<-old$targets[i,,'C']
 for(nm in names(hs)){
  sc<-esg_match_scenarios(b$eps,rep(1/nrow(b$eps),nrow(b$eps)),hs[[nm]]$H)
  W[,nm]<-esg_cvar_lp(sc$scenarios,sc$prob,pr$alpha,pr$cap)$weights
 }
 Ws[[i]]<-W
 HH<-c(list(C=b$forecast),lapply(hs,function(z)z$H));tt<-(t+1):min(t+21,nrow(pan$R))
 scores[[i]]<-do.call(rbind,lapply(names(HH),function(nm){S<-HH[[nm]]
  data.frame(origin=pan$dates[t],date=pan$dates[tt],model=nm,Gaussian_score=as.numeric(determinant(S,logarithm=TRUE)$modulus)+rowSums((pan$R[tt,,drop=FALSE]%*%solve(S))*pan$R[tt,,drop=FALSE]))
 }))
 esg_save(list(id=id,Ws=Ws,scores=scores),cache)
}
W<-aperm(array(unlist(Ws),c(ncol(pan$R),5,length(Ws))),c(3,1,2));dimnames(W)<-list(as.character(pan$dates[origins]),colnames(pan$R),colnames(Ws[[1]]))
esg_save(list(targets=W,origins=origins,protocol=pr),file.path(out,'TARGETS.rds'))
esg_csv(do.call(rbind,scores),file.path(out,'HOLDING_COVARIANCE_SCORES.csv'))
for(cost in c(0,5,10,25)){
 bt<-esg_backtest_targets(pan$R,pan$dates,origins,W,cost);keep<-bt$dates>=pr$assessment_start;rr<-bt$returns[keep,,drop=FALSE]
 esg_csv(data.frame(date=bt$dates,bt$returns,check.names=FALSE),file.path(out,paste0('returns_cost',cost,'.csv')))
 esg_csv(bt$trades,file.path(out,paste0('trades_cost',cost,'.csv')))
 perf<-do.call(rbind,lapply(colnames(rr),function(s)data.frame(strategy=s,cost_bps=cost,t(esg_metrics(rr[,s])),check.names=FALSE)))
 esg_csv(perf,file.path(out,paste0('performance_2023_cost',cost,'.csv')))
 if(cost==10){
  pairs<-c(lapply(colnames(rr)[-1],function(s)c(s,'C')),list(c('H_VAR','H_OWN'),c('H_STATE_VAR','H_STATE_OWN')))
  inf<-do.call(rbind,lapply(pr$block_lengths,function(bl)revision_pairs(rr,pairs,pr$bootstrap_B,bl)))
  esg_csv(inf,file.path(out,'PORTFOLIO_INFERENCE.csv'))
 }
}
writeLines(c('Four fixed exploratory arms; no use of returns after 2023-10-31.',
 '21-session average conditional daily covariance with fixed initial correlation, NOT exact compounded 21-day CVaR.',
 'Reestimated own/full VAR, pooled/state-dependent, at the same monthly origins as baseline.',
 'Dynamics indexed by lagged state; propagate moments then transition, consistent with daily lag packets.',
 'Conditional exponential moments approximate Gaussian mixtures when states differ.'),file.path(out,'SCOPE.txt'))
print(read.csv(file.path(out,'performance_2023_cost10.csv'))[,c('strategy','ES_loss_95')])
