for(f in sort(list.files('R',pattern='[.]R$',full.names=TRUE)))source(f)
freeze<-readRDS('outputs/FROZEN_HISTORICAL_DESIGN.rds');pr<-freeze$protocol;pr$root<-getwd()
now<-tools::md5sum(names(freeze$source_md5))
esg_assert(identical(now,freeze$source_md5),'Frozen sources changed. Do not treat changed methods as the frozen design.')
pan<-revision_read('input/archived_returns.csv');R<-pan$R;dates<-pan$dates
esg_assert(identical(esg_hash(revision_read('input/archived_returns.csv',pr$development_end)),freeze$development_identity),
  'Development data changed after freezing.')
ends<-esg_month_ends(dates);origins<-ends[dates[ends]>=pr$development_end&ends<nrow(R)]
out<-'outputs/historical_evaluation';dir.create(out,showWarnings=FALSE)
W<-array(NA_real_,c(length(origins),ncol(R),length(freeze$arms)),
  dimnames=list(as.character(dates[origins]),colnames(R),freeze$arms))
ds<-list();spec<-freeze$model
for(j in seq_along(origins)){
  t<-origins[j];cat('Frozen historical allocation',as.character(dates[t]),'\n')
  b<-revision_bundle(R,dates,t,spec,pr);a<-revision_targets(b,pr)
  W[j,,colnames(a$weights)]<-a$weights
  ref<-revision_covariance(R[1:t,,drop=FALSE],pr$models[1,],pr$warmup)
  ix<-tail(ref$valid,pr$pool);sc<-esg_match_scenarios(ref$eps[ix,,drop=FALSE],rep(1/length(ix),length(ix)),ref$forecast)
  W[j,,'C_EWMA094']<-esg_cvar_lp(sc$scenarios,sc$prob,pr$alpha,pr$cap)$weights
  W[j,,'MinVar_EWMA094']<-esg_minvar(ref$forecast,pr$cap)
  ds[[j]]<-data.frame(origin=dates[t],first_return_date=dates[t+1],last_training_date=max(b$dates),
    a$diagnostics,ambiguity_gap=unname(a$ambiguity_gap),covariance_error=a$max_covariance_error)
  esg_save(list(targets=W,completed=j,source_md5=freeze$source_md5),file.path(out,'CHECKPOINT.rds'))
}
esg_save(list(targets=W,origins=origins,source_md5=freeze$source_md5),file.path(out,'TARGETS.rds'))
esg_csv(do.call(rbind,ds),file.path(out,'information_and_scenario_audit.csv'))
for(cost in c(0,5,10,25)){
  bt<-esg_backtest_targets(R,dates,origins,W,cost)
  esg_csv(data.frame(date=bt$dates,bt$returns,check.names=FALSE),file.path(out,paste0('returns_cost',cost,'.csv')))
  esg_csv(bt$trades,file.path(out,paste0('trades_cost',cost,'.csv')))
  esg_csv(esg_performance_table(bt,cost),file.path(out,paste0('performance_cost',cost,'.csv')))
  if(cost==pr$cost){
    pairs<-list(c('C_R','C'),c('C_N','C'),c('C_WCOMP','C'),c('C_WMIX','C'),c('C_WMIX_N','C'),
      c('C','MinVar'),c('C','EqualWeight'),c('C','C_EWMA094'),c('MinVar','MinVar_EWMA094'))
    tests<-do.call(rbind,lapply(pr$block_lengths,function(bl)revision_pairs(bt$returns,pairs,pr$bootstrap_B,bl)))
    esg_csv(tests,file.path(out,'paired_inference.csv'))
  }
}
score<-revision_scores(R,dates,pr$models[pr$models$id%in%c('EWMA094',spec$id),],pr$warmup)
score<-score[score$date>=pr$historical_test_start,]
esg_csv(score,file.path(out,'covariance_scores.csv'))
wide<-reshape(score[,c('date','model','Gaussian_score')],idvar='date',timevar='model',direction='wide')
delta<-wide[[paste0('Gaussian_score.',spec$id)]]-wide$Gaussian_score.EWMA094
tabs<-lapply(pr$block_lengths,function(bl){
  boot<-esg_with_seed(esg_seed(c('historical_score',bl)),replicate(pr$bootstrap_B,mean(delta[esg_mbb_indices(length(delta),bl)])))
  data.frame(block=bl,t(esg_basic_inference(mean(delta),boot,'less')),check.names=FALSE)})
esg_csv(do.call(rbind,tabs),file.path(out,'covariance_inference.csv'))
writeLines(c('COMPLETED: frozen historical re-evaluation; not untouched confirmation.',
 'No choice, candidate strength, cost, level, or network definition was changed using this run.',
 'Only the fixed current universe is analyzed; no survival-free universe claim.',
 'Gaussian-DCC, first-stage latent posterior uncertainty, and future independent confirmation are separate modules.'),file.path(out,'STATUS.txt'))
