# Same posterior engine as R/09_fsv_research.R; independent per-chain workers.
revision_fsv_worker_identity <- function(R, dates, spec) {
  restriction <- matrix(FALSE,ncol(R),spec$factors)
  if (spec$anchored) {
    esg_assert('SPY'%in%colnames(R),'SPY anchor missing.')
    restriction[match('SPY',colnames(R)),2] <- TRUE
  }
  list(restriction=restriction,identity=esg_hash(list(
    R=R,dates=dates,spec=spec,restriction=restriction,
    R_version=as.character(getRversion()),FSV=as.character(packageVersion('factorstochvol')),
    SV=as.character(packageVersion('stochvol')),posterior=as.character(packageVersion('posterior')))))
}

revision_fsv_fit_worker <- function(R,dates,spec,output,only_chain=NULL) {
  esg_need(c('factorstochvol','posterior'));dir.create(output,recursive=TRUE,showWarnings=FALSE)
  esg_check_panel(R,dates);p<-ncol(R);nfull<-nrow(R);k<-spec$factors
  sc<-revision_fsv_scale(R,spec$scale_lambda,spec$warmup)
  valid<-if(spec$scaled)sc$valid else seq_len(nfull)
  scale_hist<-if(spec$scaled)sc$D[valid,,drop=FALSE] else matrix(.01,length(valid),p)
  scale_next<-if(spec$scaled)sc$next_scale else rep(.01,p)
  Y<-R[valid,,drop=FALSE]/scale_hist;n<-nrow(Y)
  # Zero conditional mean is fixed. No full-prefix demeaned transformation.
  fingerprint<-revision_fsv_worker_identity(R,dates,spec)
  restriction<-fingerprint$restriction
  identity<-fingerprint$identity
  esg_save(list(spec=spec,identity=identity,cutoff=max(dates),restriction=restriction),file.path(output,'CONFIGURATION.rds'))
  dput(spec,file=file.path(output,'SPECIFICATION.R'))
  hh<-array(0,c(n,p,p));fc<-matrix(0,p,p);trace<-list();ppcs<-list()
  timegrid<-unique(c(1L,ceiling(n/2),n));pairs<-which(upper.tri(matrix(0,p,p)),arr.ind=TRUE)
  summarizeS<-function(S){cc<-cov2cor(S);c(log(diag(S)),as.numeric(determinant(S,logarithm=TRUE)$modulus),cc[upper.tri(cc)])}
  snames<-c(paste0('logvariance_',colnames(R)),'logdet',paste0('corr_',colnames(R)[pairs[,1]],'_',colnames(R)[pairs[,2]]))
  for(ch in if(is.null(only_chain))seq_len(spec$chains) else only_chain){
    cache<-file.path(output,paste0('chain_',ch,'_',identity,'.rds'))
    auxiliary<-file.path(output,c(paste0('coherent_paths_chain_',ch,'.rds'),paste0('terminal_particles_chain_',ch,'.rds')))
    if(file.exists(cache)&&all(file.exists(auxiliary))){z<-readRDS(cache);cat('Verified chain cache',ch,'\n')} else {
      cat('Sampling',spec$candidate,as.character(max(dates)),'chain',ch,'/',spec$chains,'\n')
      args<-c(list(y=Y,factors=k,draws=spec$draws,burnin=spec$burnin,thin=spec$thin,
        restrict=if(spec$anchored)restriction else 'none',zeromean=TRUE,keeptime='all',
        heteroskedastic=c(spec$dynamic_idio,TRUE),runningstore=6,runningstorethin=10,
        runningstoremoments=2,signident=TRUE,signswitch=FALSE,interweaving=4,quiet=TRUE),
        spec[c('priorfacloadtype','priorfacload','priorng','priormu','priorphiidi','priorphifac',
          'priorsigmaidi','priorsigmafac','priorh0idi','priorh0fac')])
      if(!spec$dynamic_idio)args$priorhomoskedastic<-matrix(rep(spec$constant_U_prior,p),p,2,byrow=TRUE)
      sampling_path<-file.path(output,paste0('sampling_checkpoint_chain_',ch,'.rds'))
      if(file.exists(sampling_path)){
        saved<-readRDS(sampling_path);esg_assert(identical(saved$identity,identity),'Sampler identity mismatch.')
        fit<-saved$fit;rm(saved)
      }else{
        fit<-esg_with_seed(esg_seed(c(identity,ch),spec$seed),do.call(factorstochvol::fsvsample,args))
        esg_save(list(identity=identity,fit=fit),sampling_path)
      }
      writeLines('Sampler complete; postprocessing',file.path(output,paste0('chain_',ch,'_progress.txt')))
      nd<-dim(fit$facload)[3]
      esg_assert(all(dim(fit$logvar)[1:2]==c(n,p+k))&&nd>=400,'Unexpected stored posterior dimensions.')
      H<-array(0,c(n,p,p))
      for(tt in seq_len(n))H[tt,,]<-factorstochvol::runningcovmat(fit,tt,statistic='mean',type='cov')*outer(scale_hist[tt,],scale_hist[tt,])
      parnames<-if(spec$dynamic_idio)unlist(lapply(seq_len(p),function(j)paste0('idio_',colnames(R)[j],'_',c('mu','phi','sigma')))) else paste0('logU_',colnames(R))
      factor_names<-c(paste0('sorted_phi_',1:k),paste0('sorted_sigma_',1:k),paste0('sorted_log_loading_norm_',1:k))
      cn<-c(parnames,factor_names,unlist(lapply(timegrid,function(tt)paste0('t',tt,'_',snames))),paste0('forecast_',snames))
      tr<-matrix(NA_real_,nd,length(cn),dimnames=list(NULL,cn));forecast<-matrix(0,p,p)
      for(d in seq_len(nd)){
        L<-matrix(fit$facload[,,d],p,k);norm2<-colSums(L^2);o<-order(norm2,decreasing=TRUE)
        par<-if(spec$dynamic_idio)as.numeric(fit$para[,seq_len(p),d]) else fit$logvar[1,seq_len(p),d]
        vals<-c(par,fit$para[2,p+o,d],fit$para[3,p+o,d],log(norm2[o]))
        for(tt in timegrid){h<-fit$logvar[tt,,d];B<-sweep(L,2,exp(h[p+seq_len(k)]/2),'*')
          S<-(tcrossprod(B)+diag(exp(h[seq_len(p)]),p))*outer(scale_hist[tt,],scale_hist[tt,]);vals<-c(vals,summarizeS(S))}
        h<-fit$logvar[n,,d];mu<-fit$para[1,,d];phi<-fit$para[2,,d];sig<-fit$para[3,,d]
        nv<-exp(mu+phi*(h-mu)+.5*sig^2)
        if(!spec$dynamic_idio)nv[seq_len(p)]<-exp(h[seq_len(p)])
        B<-sweep(L,2,sqrt(nv[p+seq_len(k)]),'*')
        S<-(tcrossprod(B)+diag(nv[seq_len(p)],p))*outer(scale_next,scale_next)
        forecast<-forecast+S/nd;tr[d,]<-c(vals,summarizeS(S))
      }
      # Conditional posterior predictive discrepancy, distinct from forecast validation.
      # Coherent historical covariance paths use the SAME draw at every time.
      ppc<-list();take<-unique(round(seq(1,nd,length.out=min(40L,nd))))
      discrepancy<-function(E){
        lb<-vapply(seq_len(p),function(j)as.numeric(Box.test(E[,j]^2,lag=10,type='Ljung-Box')$statistic),numeric(1))
        c(max_squared_LB=max(lb),mean_squared_LB=mean(lb),max_abs=max(abs(E)))
      }
      paths<-array(NA_real_,c(n,p,p,length(take)))
      for(dd in take){
        L<-matrix(fit$facload[,,dd],p,k);E<-matrix(0,n,p)
        for(tt in seq_len(n)){
          h<-fit$logvar[tt,,dd];B<-sweep(L,2,exp(h[p+seq_len(k)]/2),'*')
          V<-(tcrossprod(B)+diag(exp(h[seq_len(p)]),p))*outer(scale_hist[tt,],scale_hist[tt,])
          paths[tt,,,match(dd,take)]<-V
          E[tt,]<-as.numeric(R[valid[tt],,drop=FALSE]%*%esg_root(V,TRUE))
        }
        obs<-discrepancy(E)
        rep<-esg_with_seed(esg_seed(c(identity,'ppc',ch,dd)),discrepancy(matrix(rnorm(n*p),n,p)))
        ppc[[length(ppc)+1]]<-data.frame(chain=ch,draw=dd,statistic=names(obs),observed=unname(obs),replicated=unname(rep))
      }
      esg_save(list(identity=identity,draws=take,H=paths,dates=dates[valid],assets=colnames(R)),
        file.path(output,paste0('coherent_paths_chain_',ch,'.rds')))
      esg_save(list(identity=identity,facload=fit$facload,para=fit$para,
        h=matrix(fit$logvar[n,,],p+k,nd),scale_next=scale_next,assets=colnames(R),cutoff=max(dates)),
        file.path(output,paste0('terminal_particles_chain_',ch,'.rds')))
      rm(paths)
      z<-list(H=H,forecast=forecast,trace=tr,stored_draws=nd,identity=identity,ppc=do.call(rbind,ppc))
      esg_save(z,cache);rm(fit);gc(FALSE)
    }
    esg_assert(identical(z$identity,identity),'Chain cache identity mismatch.')
    hh<-hh+z$H/spec$chains;fc<-fc+z$forecast/spec$chains;trace[[ch]]<-z$trace;ppcs[[ch]]<-z$ppc
  }
  if(!is.null(only_chain))return(invisible(list(chain=only_chain,identity=identity)))
  nd<-nrow(trace[[1]]);arr<-array(NA_real_,c(nd,spec$chains,ncol(trace[[1]])),dimnames=list(NULL,NULL,colnames(trace[[1]])))
  for(ch in seq_len(spec$chains))arr[,ch,]<-trace[[ch]]
  diagnostics<-as.data.frame(posterior::summarise_draws(posterior::as_draws_array(arr),
    'mean','sd','rhat','ess_bulk','ess_tail','mcse_mean'))
  diagnostics$pass<-with(diagnostics,is.finite(rhat)&rhat<=spec$max_rhat&
    is.finite(ess_bulk)&ess_bulk>=spec$min_ess&is.finite(ess_tail)&ess_tail>=spec$min_ess)
  if(length(Filter(Negate(is.null),ppcs))==spec$chains){
    pd<-do.call(rbind,ppcs);esg_csv(pd,file.path(output,'POSTERIOR_PREDICTIVE_DRAWS.csv'))
    ps<-do.call(rbind,lapply(unique(pd$statistic),function(nm){
      a<-pd[pd$statistic==nm,];data.frame(statistic=nm,posterior_paths=nrow(a),
        mean_observed=mean(a$observed),mean_replicated=mean(a$replicated),
        posterior_predictive_fraction=mean(a$replicated>=a$observed),
        scope='conditional Bayesian discrepancy; not a calibrated frequentist p-value or predictive PIT')
    }));esg_csv(ps,file.path(output,'POSTERIOR_PREDICTIVE_SUMMARY.csv'))
  }
  passed<-all(diagnostics$pass);esg_csv(diagnostics,file.path(output,'MCMC_DIAGNOSTICS.csv'))
  gate<-data.frame(candidate=spec$candidate,cutoff=max(dates),chains=spec$chains,
    requested_draws=spec$draws,thin=spec$thin,stored_draws_per_chain=nd,
    quantities=nrow(diagnostics),failed=sum(!diagnostics$pass),max_Rhat=max(diagnostics$rhat),
    min_bulk_ESS=min(diagnostics$ess_bulk),min_tail_ESS=min(diagnostics$ess_tail),passed=passed)
  esg_csv(gate,file.path(output,'MCMC_GATE.csv'))
  H<-array(NA_real_,c(nfull,p,p));H[valid,,]<-hh
  result<-list(H=H,forecast=fc,valid=valid,diagnostics=diagnostics,diagnostic_pass=passed,
    spec=spec,cutoff=max(dates),identity=identity,scope='research posterior; no downstream approval')
  esg_save(result,file.path(output,if(passed)'QUALIFIED_PREFLIGHT.rds' else 'REJECTED_PREFLIGHT.rds'))
  writeLines(c('All covariance pairs at three historical dates and the one-step conditional covariance mean are monitored.',
    'Unscaled anchored candidate uses fixed 0.01 unit conversion; all idiosyncratic variances remain dynamic.',
    'The SPY loading restriction, where active, is a substantive model assumption.',
    'Passing this gate is not model adequacy, structural identification, or portfolio superiority.'),file.path(output,'SCOPE.txt'))
  capture.output(sessionInfo(),file=file.path(output,'R_SESSION.txt'));print(gate)
  result
}
