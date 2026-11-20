# Additional exploratory development sensitivities. No allocation is reselected.
for(f in sort(list.files('R',pattern='[.]R$',full.names=TRUE)))source(f)
source('experiments/market_residuals.R')
pr<-readRDS('outputs/PRESPECIFIED_PROTOCOL.rds');pr$root<-getwd()
d<-revision_read('input/archived_returns.csv',pr$development_end)
adjusted<-revision_market_residuals(d$R,d$dates)
cut<-nrow(d$R)-50L;d2<-d$R;d2[(cut+1L):nrow(d2),]<-d2[(cut+1L):nrow(d2),]+.2
a2<-revision_market_residuals(d2,d$dates)
stopifnot(identical(adjusted$R[adjusted$rows<=cut,],a2$R[a2$rows<=cut,]),
          all(adjusted$latest_estimation_date<adjusted$dates))
out<-'outputs/network_robustness';dir.create(out,showWarnings=FALSE)
esg_csv(data.frame(check=c('market_adjustment_future_invariance','market_coefficients_strictly_lagged'),
                  pass=TRUE),file.path(out,'information_checks.csv'))
esg_csv(data.frame(date=adjusted$dates,latest_estimation_date=adjusted$latest_estimation_date,
                  adjusted$beta),file.path(out,'market_betas.csv'))
raw<-d$R[adjusted$rows,,drop=FALSE]
cases<-list(matched_raw=raw,matched_without_ESGV=raw[,colnames(raw)!='ESGV',drop=FALSE],
            matched_market_adjusted=adjusted$R)
id<-readLines('outputs/FROZEN_COVARIANCE_MODEL.txt',warn=FALSE)[1]
spec<-pr$models[match(id,pr$models$id),]
cfg<-esg_config('pilot',getwd());cfg$bootstrap_B<-pr$bootstrap_B;cfg$var_min_obs<-40L
cache_identity<-esg_hash(list(input=list(R=d$R,dates=d$dates),window=252L,
  case_names=names(cases),spec=spec,
  B=cfg$bootstrap_B,block=10L,seed=cfg$seed,
  source=tools::md5sum(c(sort(list.files('R',pattern='[.]R$',full.names=TRUE)),
                       'experiments/market_residuals.R','analysis/network_robustness.R'))))
identity_path<-file.path(out,'CACHE_IDENTITY.rds')
if(file.exists(identity_path)){
  esg_assert(identical(readRDS(identity_path),cache_identity),'Network sensitivity cache mismatch; preserve old outputs and recompute in a new folder.')
}else{
  esg_assert(!length(list.files(out,pattern='_bootstrap[.]csv$')),
    'Unversioned network bootstrap cache. Preserve it and use a clean output folder.')
  saveRDS(cache_identity,identity_path,version=2)
}
tables<-list()
for(nm in names(cases)){
  cat('Development network sensitivity:',nm,'\n')
  X<-cases[[nm]]
  cached<-file.path(out,paste0(nm,'_bootstrap.csv'))
  if(file.exists(cached)){
    z<-list(table=read.csv(cached,check.names=FALSE))
  }else{
  b<-revision_bundle(X,adjusted$dates,nrow(X),spec,pr)
  tab<-esg_conn_tables(b$var,b$hmm$labels,colnames(X),cfg$blocks,c(1,5,10,20))
  esg_csv(tab$tci,file.path(out,paste0(nm,'_tci.csv')))
  esg_csv(tab$diagnostics,file.path(out,paste0(nm,'_diagnostics.csv')))
  z<-esg_conn_bootstrap(b$packet,b$var,b$hmm$labels,colnames(X),cfg$blocks,cfg,
       paste0('development_robustness_',nm),block=10L)
  z$table$scope<-'exploratory; conditional on fitted observation-driven variance path and refitted state labels'
  esg_csv(z$table,file.path(out,paste0(nm,'_bootstrap.csv')))
  }
  t<-z$table[z$table$statistic%in%c('D_TCI','D_CE_NET','D_TECH_MINUS_FOSSIL'),]
  t$case<-nm;t$assets<-ncol(X);t$start<-min(adjusted$dates);t$end<-max(adjusted$dates)
  tables[[nm]]<-t
}
all<-do.call(rbind,tables);all$p_holm_across_nine<-p.adjust(all$p_value,'holm')
esg_csv(all,file.path(out,'primary_sensitivities.csv'))
writeLines(c('Exploratory sensitivities added after observing the initial development network.',
  'All three cases share a shortened calendar: 252 observations removed for rolling market estimation.',
  'Each case refits the covariance filter, HMM, regime labels and VAR; regimes need not denote identical dates.',
  'Market-adjusted residuals subtract contemporaneous SPY with intercept/beta estimated on the previous 252 days.',
  'SPY is excluded from that residual panel. These residuals are not investable ETF returns.',
  'The no-ESGV case changes composition on this same calendar; it does not extend the original data start.',
  'B=1999, block=10; Holm correction across all nine directional contrasts.',
  'Inference remains conditional on extracted variances and estimated classifications; no structural causality claim.'),
  file.path(out,'SCOPE.txt'))
print(all)
