# Development-only mechanism references. These are NOT calibrated hypothesis tests.
# Fixed filters/regime/VAR specification; no cross-asset volatility transmission
# is present in either generator. Common shocks and filtering remain possible.
for(f in sort(list.files('R',pattern='[.]R$',full.names=TRUE)))source(f)
pr<-readRDS('outputs/PRESPECIFIED_PROTOCOL.rds');pr$root<-getwd()
pan<-revision_read('input/archived_returns.csv',pr$development_end)
id<-readLines('outputs/FROZEN_COVARIANCE_MODEL.txt',warn=FALSE)[1]
spec<-pr$models[match(id,pr$models$id),];n<-nrow(pan$R);p<-ncol(pan$R)
train<-pan$dates<pr$calibration_start
S<-crossprod(pan$R[train,,drop=FALSE])/sum(train);L<-esg_root(S)
cfg<-esg_config('pilot',getwd())
getstat<-function(b){
  st<-esg_conn_statistics(b$var,b$hmm$labels,colnames(pan$R),cfg$blocks,10)
  h1<-esg_conn_statistics(b$var,b$hmm$labels,colnames(pan$R),cfg$blocks,1)
  c(D_TCI_H10=unname(st['D_TCI']),D_TCI_H1=unname(h1['D_TCI']),
    horizon_contrast=unname(st['D_TCI']-h1['D_TCI']),
    CE_NET_Stress=unname(st['CE_NET|Stress']),tech_minus_fossil=unname(st['D_TECH_MINUS_FOSSIL']))
}
observed<-getstat(revision_bundle(pan$R,pan$dates,n,spec,pr))
B<-as.integer(Sys.getenv('ESG_REFERENCE_REPS','100'));workers<-as.integer(Sys.getenv('ESG_REFERENCE_WORKERS','2'))
work<-function(job){
  generator<-job$generator;b<-job$replicate
  esg_with_seed(esg_seed(c('mechanism_reference',generator,b)),{
    Z<-matrix(rnorm(n*p),n,p);Y<-Z%*%L
    if(generator=='common_SV'){
      phi<-.97;sig<-.15;vv<-sig^2/(1-phi^2);h<-numeric(n)
      h[1]<-rnorm(1,0,sqrt(vv));for(t in 2:n)h[t]<-phi*h[t-1]+rnorm(1,0,sig)
      Y<-Y*exp((h-vv/2)/2)
    }
    colnames(Y)<-colnames(pan$R)
    result<-tryCatch(getstat(revision_bundle(Y,pan$dates,n,spec,pr)),error=function(e)e)
    if(inherits(result,'error'))return(data.frame(generator=generator,replicate=b,valid=FALSE,
      D_TCI_H10=NA,D_TCI_H1=NA,horizon_contrast=NA,CE_NET_Stress=NA,tech_minus_fossil=NA,error=conditionMessage(result)))
    data.frame(generator=generator,replicate=b,valid=TRUE,t(result),error='',check.names=FALSE)
  })
}
grid<-expand.grid(generator=c('iid_Gaussian','common_SV'),replicate=seq_len(B),stringsAsFactors=FALSE)
jobs<-lapply(seq_len(nrow(grid)),function(i)grid[i,,drop=FALSE])
if(.Platform$OS.type=='unix'&&workers>1)z<-parallel::mclapply(jobs,work,mc.cores=workers,mc.set.seed=FALSE) else z<-lapply(jobs,work)
rows<-do.call(rbind,z);out<-'outputs/mechanism_references';dir.create(out,showWarnings=FALSE)
esg_csv(rows,file.path(out,'replicates.csv'))
summary<-list()
for(g in unique(rows$generator))for(nm in names(observed)){
  rr<-rows[rows$generator==g,];x<-rr[[nm]][rr$valid]
  summary[[length(summary)+1]]<-data.frame(generator=g,statistic=nm,observed=observed[nm],
    valid=length(x),requested=B,median=median(x),q025=quantile(x,.025,names=FALSE),
    q975=quantile(x,.975,names=FALSE),reference_fraction_ge_observed=mean(x>=observed[nm]),
    scope='model-implied simulation reference, not a calibrated p-value',row.names=NULL)
}
esg_csv(do.call(rbind,summary),file.path(out,'summary.csv'))
writeLines(c('Both generators have no cross-asset volatility feedback coefficients.',
 'iid_Gaussian: independent time draws with a covariance estimated before 2022.',
 'common_SV: one common log-SV driver, phi=.97, innovation SD=.15, stationary initialization and unit mean variance.',
 'The Gaussian standardized innovations share the fixed pre-2022 covariance in both generators.',
 'The observation-driven covariance filter, HMM, daily-lag regime VAR, and GFEVD are re-estimated per replicate.',
 'This isolates common shocks/filter-induced dependence; reference fractions are NOT p-values.',
 'Parameters and model selection uncertainty are not calibrated by these illustrative generators.'),file.path(out,'SCOPE.txt'))
print(do.call(rbind,summary))
