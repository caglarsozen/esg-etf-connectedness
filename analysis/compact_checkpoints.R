# Postprocessing cache: retain all diagnostics and all 40 coherent paths per chain.
for(f in sort(list.files('R',pattern='[.]R$',full.names=TRUE)))source(f)
folder<-Sys.getenv('ESG_FSV_FOLDER');stopifnot(nzchar(folder),file.exists(file.path(folder,'MCMC_GATE.csv')))
for(ch in 1:4){
 path<-file.path(folder,paste0('sampling_checkpoint_chain_',ch,'.rds'));z<-readRDS(path)
 if(isTRUE(z$compacted)){rm(z);gc(FALSE);next}
 n<-dim(z$fit$logvar)[1];nd<-dim(z$fit$logvar)[3]
 times<-unique(c(1L,ceiling(n/2),n));take<-unique(round(seq(1,nd,length.out=min(40L,nd))))
 kt<-z$fit$logvar[times,,,drop=FALSE];kp<-z$fit$logvar[,,take,drop=FALSE]
 z$fit$logvar[setdiff(seq_len(n),times),,setdiff(seq_len(nd),take)]<-NA_real_
 stopifnot(identical(kt,z$fit$logvar[times,,,drop=FALSE]),identical(kp,z$fit$logvar[,,take,drop=FALSE]))
 fields<-if(is.environment(z$fit))ls(z$fit,all.names=TRUE) else names(z$fit)
 for(nm in setdiff(fields,c('y','facload','para','logvar','runningstore')))z$fit[[nm]]<-NULL
 z$compacted<-TRUE;z$retained_diagnostic_times<-times;z$retained_full_path_draws<-take
 z$scope<-'Postprocessing checkpoint, not every latent cell; all monitored/PPC/network values retained.'
 esg_save(z,path);rm(z,kt,kp);gc(FALSE);cat('Compacted checkpoint',ch,'\n')
}
