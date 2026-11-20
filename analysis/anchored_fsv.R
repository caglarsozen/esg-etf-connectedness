# Long anchored experiment: default 160000 draws after 40000 burn-in, thinning 40.
for(f in sort(list.files('R',pattern='[.]R$',full.names=TRUE)))source(f)
source('experiments/parallel_fsv_engine.R')
draws<-as.integer(Sys.getenv('ESG_FSV_DRAWS','160000'))
burnin<-as.integer(Sys.getenv('ESG_FSV_BURNIN','40000'))
thin<-as.integer(Sys.getenv('ESG_FSV_THIN','40'))
spec<-revision_fsv_spec('anchored_dynamic',draws,burnin,thin)
cutoff<-as.Date(Sys.getenv('ESG_FSV_CUTOFF','2023-02-28'))
esg_assert(cutoff<=as.Date('2023-10-31'),'Preflight may not inspect historical test.')
pan<-revision_read('input/archived_returns.csv',cutoff)
out<-file.path('outputs','fsv',paste('anchored_dynamic',format(cutoff,'%Y%m%d'),draws,burnin,thin,sep='_'))
esg_need(c('factorstochvol','posterior'))
config_path<-file.path(out,'CONFIGURATION.rds')
if(file.exists(config_path)){
  stored<-readRDS(config_path)
  expected<-revision_fsv_worker_identity(pan$R,pan$dates,spec)$identity
  esg_assert(identical(stored$identity,expected),paste(
    'Existing chain archive uses a different data/specification/software identity.',
    'Use the recorded software versions or move this result folder aside before refitting:',out))
}
workers<-max(1L,min(4L,as.integer(Sys.getenv('ESG_WORKERS','2'))))
root<-getwd()
work<-function(ch)revision_fsv_fit_worker(pan$R,pan$dates,spec,out,only_chain=ch)
if(workers==1){ans<-lapply(1:4,work)}else if(.Platform$OS.type=='unix'){
  ans<-parallel::mclapply(1:4,work,mc.cores=workers,mc.preschedule=FALSE,mc.set.seed=FALSE)
}else{
  cl<-parallel::makePSOCKcluster(workers)
  parallel::clusterExport(cl,c('root','pan','spec','out'),envir=environment())
  parallel::clusterEvalQ(cl,{
    setwd(root)
    for(f in sort(list.files('R',pattern='[.]R$',full.names=TRUE)))source(f)
    source('experiments/parallel_fsv_engine.R')
    NULL
  })
  ans<-tryCatch(parallel::parLapply(cl,1:4,function(ch)
    revision_fsv_fit_worker(pan$R,pan$dates,spec,out,only_chain=ch)),
    finally=parallel::stopCluster(cl))
}
esg_assert(!any(vapply(ans,inherits,logical(1),'try-error')),'An MCMC worker failed; inspect logs.')
fit<-revision_fsv_fit_worker(pan$R,pan$dates,spec,out)
