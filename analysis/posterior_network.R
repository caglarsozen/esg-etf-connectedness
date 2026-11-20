# Propagate each coherent first-stage posterior path through freshly fitted HMM and VAR.
for(f in sort(list.files('R',pattern='[.]R$',full.names=TRUE)))source(f)
folder<-Sys.getenv('ESG_FSV_FOLDER');stopifnot(nzchar(folder))
gate<-read.csv(file.path(folder,'MCMC_GATE.csv'));stopifnot(gate$passed)
out<-file.path(folder,'network_uncertainty');dir.create(out,showWarnings=FALSE)
cfg<-esg_config('pilot',getwd());rows<-list();failures<-list();root<-getwd()
path_work<-function(ch){
 z<-readRDS(file.path(folder,paste0('coherent_paths_chain_',ch,'.rds')))
 cache<-file.path(out,paste0('path_statistics_chain_',ch,'.rds'))
 id<-esg_hash(list(z$identity,z$draws,source=tools::md5sum(c(
   sort(list.files('R',pattern='[.]R$',full.names=TRUE)),
   'analysis/posterior_network.R'))))
 if(file.exists(cache)){
  zz<-readRDS(cache);stopifnot(identical(zz$id,id));return(zz)
 }else{
 rr<-list();err<-list()
 for(j in seq_along(z$draws)){
  ans<-tryCatch({
   ff<-esg_features(z$H[,,,j]);colnames(ff$X)<-z$assets
   h<-esg_hmm_fit(ff$Z,3L,starts=3L,seed=esg_seed(c('posterior_network',gate$cutoff,ch,z$draws[j])),tol=1e-6)
   esg_assert(h$converged,'HMM did not converge')
   pk<-esg_lag_packets(ff$X,h$state,1L,z$dates);vf<-esg_var_from_packets(pk,3L,40L,.1)
   v<-esg_conn_statistics(vf,h$labels,z$assets,cfg$blocks,10L)
   v['role_switch']<-as.numeric(v['CE_NET|Calm']<0&&v['CE_NET|Stress']>0)
   data.frame(chain=ch,draw=z$draws[j],statistic=names(v),value=as.numeric(v))
  },error=function(e){err[[length(err)+1]]<<-data.frame(chain=ch,draw=z$draws[j],error=conditionMessage(e));NULL})
  if(!is.null(ans))rr[[length(rr)+1]]<-ans
 }
 zz<-list(id=id,rows=do.call(rbind,rr),failures=do.call(rbind,err))
 esg_save(zz,cache)
 }
 cat('Posterior network chain',ch,'completed\n')
 zz
}
workers<-max(1L,min(4L,as.integer(Sys.getenv('ESG_WORKERS','2'))))
if(workers==1L)ans<-lapply(1:4,path_work) else if(.Platform$OS.type=='unix'){
 ans<-parallel::mclapply(1:4,path_work,mc.cores=workers,mc.preschedule=FALSE,mc.set.seed=FALSE)
}else{
 cl<-parallel::makePSOCKcluster(workers)
 parallel::clusterExport(cl,c('root','folder','gate','out','cfg','path_work'),envir=environment())
 parallel::clusterEvalQ(cl,{
  setwd(root);for(f in sort(list.files('R',pattern='[.]R$',full.names=TRUE)))source(f)
  NULL
 })
 ans<-tryCatch(parallel::parLapply(cl,1:4,path_work),finally=parallel::stopCluster(cl))
}
stopifnot(!any(vapply(ans,inherits,logical(1),'try-error')))
rows<-lapply(ans,function(z)z$rows);failures<-lapply(ans,function(z)z$failures)
d<-do.call(rbind,rows);err<-do.call(rbind,failures)
esg_csv(d,file.path(out,'PATH_STATISTICS.csv'));if(!is.null(err))esg_csv(err,file.path(out,'FAILURES.csv'))
tab<-do.call(rbind,lapply(unique(d$statistic),function(nm){x<-d$value[d$statistic==nm]
 data.frame(statistic=nm,paths=length(x),mean=mean(x),median=median(x),q025=unname(quantile(x,.025)),
 q975=unname(quantile(x,.975)),fraction_positive=mean(x>0),
 scope='first-stage posterior propagation; downstream plug-in estimates; not a joint-model credible interval')
}))
esg_csv(tab,file.path(out,'POSTERIOR_PROPAGATED_SUMMARY.csv'))
writeLines(c('All successful and failed path fits are retained. Cutoffs overlap; these are not independent validations.',
 'HMM and VAR refitted per coherent posterior covariance path. No allocation uses these retrospective smoothed states.',
 'State labels use model volatility/correlation, not an external stress definition.'),file.path(out,'SCOPE.txt'))
print(tab[tab$statistic%in%c('D_TCI','D_TECH_MINUS_FOSSIL','STRESS_TECH_MINUS_FOSSIL','role_switch'),])
