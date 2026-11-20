# Run from this project's root. Does not analyze returns after 2023-10-31.
for(f in sort(list.files('R',pattern='[.]R$',full.names=TRUE)))source(f)
pr<-revision_protocol(getwd());dir.create('outputs',showWarnings=FALSE)
saveRDS(pr,'outputs/PRESPECIFIED_PROTOCOL.rds',version=2)
dput(pr,file='outputs/PRESPECIFIED_PROTOCOL.R')
pan<-revision_read('input/archived_returns.csv',pr$development_end);R<-pan$R;dates<-pan$dates
esg_csv(data.frame(start=min(dates),end=max(dates),days=nrow(R),assets=ncol(R),
  scope='development only; old pilot examined later historical period'),'outputs/data_scope.csv')
scores<-revision_scores(R,dates,pr$models,pr$warmup)
esg_csv(scores,'outputs/covariance_daily_scores.csv')
cal<-scores[scores$date>=pr$calibration_start&scores$date<=pr$calibration_end,]
rank<-aggregate(Gaussian_score~model,cal,mean);rank<-rank[order(rank$Gaussian_score),]
chosen<-rank$model[1];spec<-pr$models[match(chosen,pr$models$id),]
esg_csv(rank,'outputs/covariance_calibration_ranking.csv')
writeLines(chosen,'outputs/FROZEN_COVARIANCE_MODEL.txt')
assess<-scores[scores$date>=pr$assessment_start&scores$date<=pr$development_end,]
wide<-reshape(assess[,c('date','model','Gaussian_score')],idvar='date',timevar='model',direction='wide')
base<-wide$Gaussian_score.EWMA094;inf<-list()
for(m in pr$models$id[pr$models$id!='EWMA094'])for(bl in pr$block_lengths){
  delta<-wide[[paste0('Gaussian_score.',m)]]-base
  boot<-esg_with_seed(esg_seed(c('covscore',m,bl)),replicate(pr$bootstrap_B,mean(delta[esg_mbb_indices(length(delta),bl)])))
  inf[[length(inf)+1]]<-data.frame(model=m,benchmark='EWMA094',block=bl,
    t(esg_basic_inference(mean(delta),boot,'less')),check.names=FALSE)
}
inf<-do.call(rbind,inf);inf$p_holm<-ave(inf$p_value,inf$block,FUN=function(x)p.adjust(x,'holm'))
esg_csv(inf,'outputs/covariance_development_inference.csv')
for(m in unique(c('EWMA094',chosen))){
  cv<-revision_covariance(R,pr$models[match(m,pr$models$id),],pr$warmup)
  esg_csv(revision_residual_checks(cv,dates,pr$assessment_start,pr$development_end),
    paste0('outputs/residual_checks_',m,'.csv'))
}
cat('Frozen covariance:',chosen,'\n')
ends<-esg_month_ends(dates);origins<-ends[dates[ends]>=as.Date('2021-12-01')&ends<nrow(R)]
Ws<-list();diagnostics<-list()
for(i in seq_along(origins)){
  t<-origins[i];cat('Development allocation',as.character(dates[t]),'\n')
  # All calibration and assessment portfolios use fixed candidate definitions.
  # Model choice uses 2022 scores; 2022 portfolio results are tuning diagnostics.
  b<-revision_bundle(R,dates,t,spec,pr);a<-revision_targets(b,pr)
  Ws[[i]]<-a$weights
  diagnostics[[i]]<-data.frame(origin=dates[t],a$diagnostics,
    covariance_error=a$max_covariance_error,ambiguity_gap=a$ambiguity_gap,
    guard_limit=a$guard_limit,baseline_CVaR=a$baseline_risk['C'])
  esg_save(list(Ws=Ws,origins=origins[seq_len(i)],model=spec,protocol=pr),
    'outputs/DEVELOPMENT_checkpoint.rds')
}
W<-array(unlist(Ws),c(ncol(R),ncol(Ws[[1]]),length(Ws)))
W<-aperm(W,c(3,1,2));dimnames(W)<-list(as.character(dates[origins]),colnames(R),colnames(Ws[[1]]))
esg_save(list(targets=W,origins=origins,protocol=pr,model=spec),'outputs/DEVELOPMENT_targets.rds')
esg_csv(do.call(rbind,diagnostics),'outputs/scenario_diagnostics.csv')
for(cost in c(0,5,10,25)){
  bt<-esg_backtest_targets(R,dates,origins,W,cost)
  esg_csv(data.frame(date=bt$dates,bt$returns,check.names=FALSE),paste0('outputs/development_returns_cost',cost,'.csv'))
  esg_csv(bt$trades,paste0('outputs/development_trades_cost',cost,'.csv'))
  keep<-bt$dates>=pr$assessment_start
  perf<-do.call(rbind,lapply(colnames(bt$returns),function(s)data.frame(strategy=s,cost_bps=cost,
    t(esg_metrics(bt$returns[keep,s])),check.names=FALSE)))
  esg_csv(perf,paste0('outputs/performance_2023_cost',cost,'.csv'))
  if(cost==pr$cost){
    rr<-bt$returns[keep,,drop=FALSE]
    pairs<-list(c('C_R','C'),c('C_N','C'),c('C_WCOMP','C'),c('C_WMIX','C'),
      c('C_WMIX_N','C'),c('C','MinVar'),c('C','EqualWeight'))
    all<-do.call(rbind,lapply(pr$block_lengths,function(bl)revision_pairs(rr,pairs,pr$bootstrap_B,bl)))
    esg_csv(all,'outputs/portfolio_development_inference.csv')
  }
}
capture.output(sessionInfo(),file='outputs/R_SESSION.txt')
writeLines(c('Development computation completed. No post-2023-10-31 return was analyzed.',
 'See README.md for the scope of this assessment.'),'outputs/RUN_STATUS.txt')
