# Explicit research specifications. Failed candidates are saved and rejected.
revision_fsv_scale <- function(R,lambda=.98,warmup=126L) {
  n<-nrow(R);p<-ncol(R);v<-colMeans(R[1:warmup,,drop=FALSE]^2)
  D<-matrix(NA_real_,n,p,dimnames=dimnames(R))
  for(tt in (warmup+1):n){D[tt,]<-sqrt(v);v<-lambda*v+(1-lambda)*R[tt,]^2}
  list(D=D,next_scale=sqrt(v),valid=(warmup+1):n)
}
revision_fsv_spec <- function(candidate='scaled_constant',draws=12000L,burnin=6000L,thin=4L) {
  esg_assert(candidate%in%c('scaled_constant','anchored_dynamic','anchored_scaled'),'Unknown FSV candidate.')
  list(candidate=candidate,factors=2L,chains=4L,draws=as.integer(draws),burnin=as.integer(burnin),
    thin=as.integer(thin),scaled=candidate!='anchored_dynamic',
    dynamic_idio=candidate=='anchored_dynamic',anchored=candidate!='scaled_constant',
    scale_lambda=.98,warmup=126L,max_rhat=1.01,min_ess=400,
    priorfacloadtype='rowwiseng',priorfacload=.1,priorng=c(1,1),
    priormu=c(0,10),priorphiidi=c(10,3),priorphifac=c(10,3),
    priorsigmaidi=if(candidate=='anchored_dynamic').2 else 1,
    priorsigmafac=1,priorh0idi='stationary',priorh0fac='stationary',
    constant_U_prior=c(shape=1.1,rate=.055),seed=20261113L)
}
revision_fsv_fit <- function(R,dates,spec,output) {
  esg_need(c('factorstochvol','posterior'));dir.create(output,recursive=TRUE,showWarnings=FALSE)
  esg_check_panel(R,dates);p<-ncol(R);nfull<-nrow(R);k<-spec$factors
  sc<-revision_fsv_scale(R,spec$scale_lambda,spec$warmup)
  valid<-if(spec$scaled)sc$valid else seq_len(nfull)
  scale_hist<-if(spec$scaled)sc$D[valid,,drop=FALSE] else matrix(.01,length(valid),p)
  scale_next<-if(spec$scaled)sc$next_scale else rep(.01,p)
  Y<-R[valid,,drop=FALSE]/scale_hist;n<-nrow(Y)
  # Zero conditional mean is fixed. No full-prefix demeaned transformation.
  restriction<-matrix(FALSE,p,k)
  if(spec$anchored){esg_assert('SPY'%in%colnames(R),'SPY anchor missing.');restriction[match('SPY',colnames(R)),2]<-TRUE}
  identity<-esg_hash(list(R=R,dates=dates,spec=spec,restriction=restriction,
    R_version=as.character(getRversion()),FSV=as.character(packageVersion('factorstochvol')),
    SV=as.character(packageVersion('stochvol')),posterior=as.character(packageVersion('posterior'))))
  esg_save(list(spec=spec,identity=identity,cutoff=max(dates),restriction=restriction),file.path(output,'CONFIGURATION.rds'))
  dput(spec,file=file.path(output,'SPECIFICATION.R'))
  hh<-array(0,c(n,p,p));fc<-matrix(0,p,p);trace<-list()
  timegrid<-unique(c(1L,ceiling(n/2),n));pairs<-which(upper.tri(matrix(0,p,p)),arr.ind=TRUE)
  summarizeS<-function(S){cc<-cov2cor(S);c(log(diag(S)),as.numeric(determinant(S,logarithm=TRUE)$modulus),cc[upper.tri(cc)])}
  snames<-c(paste0('logvariance_',colnames(R)),'logdet',paste0('corr_',colnames(R)[pairs[,1]],'_',colnames(R)[pairs[,2]]))
  for(ch in seq_len(spec$chains)){
    cache<-file.path(output,paste0('chain_',ch,'_',identity,'.rds'))
    if(file.exists(cache)){z<-readRDS(cache);cat('Verified chain cache',ch,'\n')} else {
      cat('Sampling',spec$candidate,as.character(max(dates)),'chain',ch,'/',spec$chains,'\n')
      args<-c(list(y=Y,factors=k,draws=spec$draws,burnin=spec$burnin,thin=spec$thin,
        restrict=if(spec$anchored)restriction else 'none',zeromean=TRUE,keeptime='all',
        heteroskedastic=c(spec$dynamic_idio,TRUE),runningstore=6,runningstorethin=10,
        runningstoremoments=2,signident=TRUE,signswitch=FALSE,interweaving=4,quiet=TRUE),
        spec[c('priorfacloadtype','priorfacload','priorng','priormu','priorphiidi','priorphifac',
          'priorsigmaidi','priorsigmafac','priorh0idi','priorh0fac')])
      if(!spec$dynamic_idio)args$priorhomoskedastic<-matrix(rep(spec$constant_U_prior,p),p,2,byrow=TRUE)
      fit<-esg_with_seed(esg_seed(c(identity,ch),spec$seed),do.call(factorstochvol::fsvsample,args))
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
      z<-list(H=H,forecast=forecast,trace=tr,stored_draws=nd,identity=identity)
      esg_save(z,cache);rm(fit);gc(FALSE)
    }
    esg_assert(identical(z$identity,identity),'Chain cache identity mismatch.')
    hh<-hh+z$H/spec$chains;fc<-fc+z$forecast/spec$chains;trace[[ch]]<-z$trace
  }
  nd<-nrow(trace[[1]]);arr<-array(NA_real_,c(nd,spec$chains,ncol(trace[[1]])),dimnames=list(NULL,NULL,colnames(trace[[1]])))
  for(ch in seq_len(spec$chains))arr[,ch,]<-trace[[ch]]
  diagnostics<-as.data.frame(posterior::summarise_draws(posterior::as_draws_array(arr),
    'mean','sd','rhat','ess_bulk','ess_tail','mcse_mean'))
  diagnostics$pass<-with(diagnostics,is.finite(rhat)&rhat<=spec$max_rhat&
    is.finite(ess_bulk)&ess_bulk>=spec$min_ess&is.finite(ess_tail)&ess_tail>=spec$min_ess)
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
    'Positive diagonal scaling transports known past-only variance dynamics; scale-estimation uncertainty is not integrated.',
    'The SPY loading restriction, where active, is a substantive model assumption.',
    'Passing this gate is not model adequacy, structural identification, or portfolio superiority.'),file.path(output,'SCOPE.txt'))
  capture.output(sessionInfo(),file=file.path(output,'R_SESSION.txt'));print(gate)
  result
}
