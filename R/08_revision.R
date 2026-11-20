# Revision research engine. All returns are simple decimals; forecasts are causal.
# Observation-driven alternatives are not latent Factor-SV estimates.
revision_protocol <- function(root=getwd()) {
  list(root=root, seed=20261113L,
    calibration_start=as.Date('2022-01-01'), calibration_end=as.Date('2022-12-31'),
    assessment_start=as.Date('2023-01-01'), development_end=as.Date('2023-10-31'),
    historical_test_start=as.Date('2023-11-01'),
    historical_test_status='previously examined in an EWMA pilot; historical evaluation',
    warmup=126L, cap=.35, alpha=.95, pool=1260L, cost=10,
    guard=.05, network_lambda=.03, bootstrap_B=1999L,
    block_lengths=c(5L,10L,20L), states=3L,
    models=data.frame(id=c('EWMA094','EWMA098','SHRINK005','SHRINK020','SHRINK050','SCALE_CORR','PCA_DYNAMIC'),
      kind=c('ewma','ewma','shrink','shrink','shrink','scale_corr','pca'),
      lambda=c(.94,.98,.94,.94,.94,.94,.94),
      shrink=c(0,0,.005,.02,.05,.02,.10),stringsAsFactors=FALSE))
}
revision_read <- function(file, through=NULL) {
  if (!file.exists(file)) stop('Return panel required at ', file, '. See README.md.', call. = FALSE)
  d<-read.csv(file,check.names=FALSE);dates<-as.Date(d[[1]])
  R<-as.matrix(d[,-1,drop=FALSE]);storage.mode(R)<-'double';esg_check_panel(R,dates)
  if(!is.null(through)){j<-dates<=through;R<-R[j,,drop=FALSE];dates<-dates[j]}
  list(R=R,dates=dates)
}
revision_covariance <- function(R,spec,warmup=126L) {
  n<-nrow(R);p<-ncol(R);esg_assert(n>warmup,'Insufficient history.')
  S<-esg_spd(crossprod(R[1:warmup,,drop=FALSE])/warmup)
  S0<-S;v<-diag(S);Q<-cov2cor(S);lam<-spec$lambda;rho<-spec$shrink
  # PCA basis is fixed using warmup only. This is an observed-data rotation,
  # with all PC variances dynamic, not a latent stochastic volatility model.
  U<-eigen(S,symmetric=TRUE)$vectors;pcv<-diag(t(U)%*%S%*%U)
  Cpc<-matrix(0,p,p);diag(Cpc)<-1
  H<-array(NA_real_,c(n,p,p));eps<-matrix(NA_real_,n,p,dimnames=dimnames(R))
  makeS<-function(){
    if(spec$kind=='ewma')return(S)
    if(spec$kind=='shrink')return((1-rho)*S+rho*diag(diag(S),p))
    if(spec$kind=='scale_corr'){
      C<-cov2cor(Q);C<-(1-rho)*C+rho*diag(p)
      return(C*outer(sqrt(v),sqrt(v)))
    }
    U%*%diag(pcv,p)%*%t(U)*(1-rho)+rho*S
  }
  for(tt in (warmup+1):n){
    St<-esg_spd(makeS());H[tt,,]<-St
    eps[tt,]<-as.numeric(R[tt,,drop=FALSE]%*%esg_root(St,TRUE))
    z<-R[tt,]/sqrt(v)
    Q<-.98*Q+.02*tcrossprod(z)
    v<-lam*v+(1-lam)*R[tt,]^2
    pcv<-lam*pcv+(1-lam)*as.numeric(R[tt,,drop=FALSE]%*%U)^2
    S<-lam*S+(1-lam)*tcrossprod(R[tt,])
  }
  list(H=H,forecast=esg_spd(makeS()),eps=eps,valid=(warmup+1):n,
    model=spec$id,spec=spec,initial_covariance=S0)
}
revision_scores <- function(R,dates,models,warmup) {
  out<-list()
  for(j in seq_len(nrow(models))){
    cat('Covariance model',models$id[j],'\n')
    cv<-revision_covariance(R,models[j,],warmup);n<-nrow(R)
    score<-logdet<-condition<-rep(NA_real_,n)
    for(tt in cv$valid){S<-cv$H[tt,,];logdet[tt]<-as.numeric(determinant(S,logarithm=TRUE)$modulus)
      score[tt]<-logdet[tt]+sum(R[tt,]*solve(S,R[tt,]));condition[tt]<-kappa(S,exact=TRUE)}
    out[[j]]<-data.frame(date=dates,model=models$id[j],Gaussian_score=score,
      log_determinant=logdet,condition_number=condition)
  };do.call(rbind,out)
}
revision_residual_checks <- function(cv,dates,start,end) {
  ix<-which(dates>=start&dates<=end&complete.cases(cv$eps))
  out<-do.call(rbind,lapply(seq_len(ncol(cv$eps)),function(j){
    x<-cv$eps[ix,j];u<-(x-mean(x))/sd(x)
    data.frame(asset=colnames(cv$eps)[j],n=length(x),mean=mean(x),sd=sd(x),
      squared_LB_p=Box.test(x^2,lag=10,type='Ljung-Box')$p.value,
      return_LB_p=Box.test(x,lag=10,type='Ljung-Box')$p.value,
      skewness=mean(u^3),excess_kurtosis=mean(u^4)-3,
      scope='one-step plug-in symmetric whitening; dependence diagnostic, not a calibrated PIT')
  }));out$squared_Holm_p<-p.adjust(out$squared_LB_p,'holm');out
}
# Unified epigraph LP for worst component CVaR and convex-hull mixture CVaR.
# For mixtures, the VaR threshold MUST be shared across components.
# For max_s CVaR_s it is separate. They coincide only in special cases.
revision_robust_cvar <- function(scenarios,probabilities,alpha=.95,cap=.35,
    ambiguity=c('components','mixtures'),guard_scenarios=NULL,guard_prob=NULL,
    guard_limit=NULL,network_score=NULL,network_lambda=0,risk_scale=NULL) {
  ambiguity<-match.arg(ambiguity);esg_need('lpSolve')
  K<-length(scenarios);p<-ncol(scenarios[[1]])
  esg_assert(K>0&&p*cap>=1&&length(probabilities)==K,'Invalid robust components.')
  S<-lapply(scenarios,function(x)as.matrix(x)*100)
  P<-lapply(probabilities,function(x)x/sum(x))
  esg_assert(all(vapply(seq_len(K),function(k)ncol(S[[k]])==p&&
    length(P[[k]])==nrow(S[[k]])&&all(is.finite(S[[k]]))&&
    all(is.finite(P[[k]])&P[[k]]>=0),logical(1))),'Invalid robust scenarios.')
  if(is.null(network_score))network_score<-rep(0,p)
  if(is.null(risk_scale))risk_scale<-max(esg_cvar(-rowMeans(S[[1]]/100),P[[1]],alpha),1e-6)
  esg_assert(length(network_score)==p&&all(is.finite(network_score))&&
    network_lambda>=0&&is.finite(risk_scale)&&risk_scale>0,'Invalid network objective.')
  # variables w, t_plus, t_minus, z+/z-, tail slacks, optional guard.
  obj<-c(network_lambda*100*risk_scale*network_score,1,-1)
  trip<-rbind(cbind(1,seq_len(p),1),cbind(1+seq_len(p),seq_len(p),1))
  rhs<-c(1,rep(cap,p));dir<-c('=',rep('<=',p))
  alloc<-function(m){old<-length(obj);obj<<-c(obj,rep(0,m));old+seq_len(m)}
  add<-function(cols,vals,right=0){r<-length(rhs)+1L
    trip<<-rbind(trip,cbind(r,cols,vals));rhs<<-c(rhs,right);dir<<-c(dir,'<=');r}
  common_z<-if(ambiguity=='mixtures')alloc(2) else NULL
  for(k in seq_len(K)){
    z<-if(is.null(common_z))alloc(2) else common_z;u<-alloc(nrow(S[[k]]))
    rows<-length(rhs)+seq_len(nrow(S[[k]]))
    trip<-rbind(trip,cbind(rep(rows,p),rep(seq_len(p),each=length(rows)),-as.vector(S[[k]])),
      cbind(rows,z[1],-1),cbind(rows,z[2],1),cbind(rows,u,-1))
    rhs<-c(rhs,rep(0,length(rows)));dir<-c(dir,rep('<=',length(rows)))
    add(c(z,u,p+1,p+2),c(1,-1,P[[k]]/(1-alpha),-1,1))
  }
  if(!is.null(guard_scenarios)){
    G<-as.matrix(guard_scenarios)*100;m<-nrow(G)
    if(is.null(guard_prob))guard_prob<-rep(1/m,m);guard_prob<-guard_prob/sum(guard_prob)
    esg_assert(ncol(G)==p&&length(guard_prob)==m&&all(is.finite(G))&&
      all(is.finite(guard_prob)&guard_prob>=0)&&length(guard_limit)==1&&
      is.finite(guard_limit),'Invalid robust guard.')
    z<-alloc(2);u<-alloc(m);rows<-length(rhs)+seq_len(m)
    trip<-rbind(trip,cbind(rep(rows,p),rep(seq_len(p),each=m),-as.vector(G)),
      cbind(rows,z[1],-1),cbind(rows,z[2],1),cbind(rows,u,-1))
    rhs<-c(rhs,rep(0,m));dir<-c(dir,rep('<=',m))
    add(c(z,u),c(1,-1,guard_prob/(1-alpha)),guard_limit*100)
  }
  fit<-lpSolve::lp('min',obj,const.dir=dir,const.rhs=rhs,dense.const=trip,scale=196)
  esg_assert(fit$status==0,paste('Robust LP failed',fit$status))
  w<-fit$solution[seq_len(p)];err<-max(abs(sum(w)-1),pmax(-w,0),pmax(w-cap,0))
  esg_assert(err<1e-6,'Robust weights violate constraints.')
  if(!is.null(guard_scenarios))esg_assert(esg_cvar(-as.numeric(guard_scenarios%*%w),guard_prob,alpha)<=guard_limit+1e-7,'Robust guard failed.')
  losses<-lapply(S,function(x)-as.numeric(x%*%w)/100)
  component_risk<-vapply(seq_len(K),function(k)esg_cvar(losses[[k]],P[[k]],alpha),numeric(1))
  # Evaluate the shared-threshold representation independently at LP threshold.
  risk<-if(ambiguity=='components')max(component_risk) else {
    zval<-(fit$solution[common_z[1]]-fit$solution[common_z[2]])/100
    max(vapply(seq_len(K),function(k)zval+sum(P[[k]]*pmax(losses[[k]]-zval,0))/(1-alpha),numeric(1)))
  }
  esg_assert(abs(risk+network_lambda*risk_scale*sum(w*network_score)-fit$objval/100)<1e-6,'Robust objective mismatch.')
  list(weights=w,risk=risk,component_risks=component_risk,objective=fit$objval/100,
    constraint_error=err,status=fit$status,ambiguity=ambiguity)
}
revision_bundle <- function(R,dates,origin,spec,protocol) {
  R<-R[seq_len(origin),,drop=FALSE];dates<-dates[seq_len(origin)]
  cv<-revision_covariance(R,spec,protocol$warmup);f<-esg_features(cv$H)
  valid<-cv$valid;X<-f$X[valid,,drop=FALSE];colnames(X)<-colnames(R);Z<-f$Z[valid,,drop=FALSE]
  h<-esg_hmm_fit(Z,k=protocol$states,starts=3L,seed=esg_seed(c('regime',tail(dates,1))),tol=1e-6)
  esg_assert(h$converged,'HMM failed convergence; no fallback.')
  nextZ<-esg_features(array(cv$forecast,c(1,ncol(R),ncol(R))))$Z
  nextp<-as.numeric(tail(esg_hmm_filter(h,rbind(Z,nextZ)),1));nextp<-nextp/sum(nextp)
  pk<-esg_lag_packets(X,h$state,1L,dates[valid])
  vf<-esg_var_from_packets(pk,protocol$states,min_obs=40,lambda=.1)
  cfg<-esg_config('pilot',protocol$root)
  network<-esg_network_score(vf,nextp,colnames(R),cfg$blocks,10L)
  pool<-tail(valid,protocol$pool)
  list(R=R,dates=dates,forecast=cv$forecast,eps=cv$eps[pool,,drop=FALSE],
    history_prob=h$filtered[match(pool,valid),,drop=FALSE],next_prob=nextp,
    network=network,hmm=h,var=vf,packet=pk,features=f,cv=cv)
}
revision_targets <- function(b,pr) {
  E<-b$eps;m<-nrow(E);p<-ncol(E);base<-esg_match_scenarios(E,rep(1/m,m),b$forecast)
  fit0<-esg_cvar_lp(base$scenarios,base$prob,pr$alpha,pr$cap)
  lim<-fit0$CVaR+pr$guard*max(fit0$CVaR,1e-6)
  scale<-max(esg_cvar(-rowMeans(base$scenarios),base$prob,pr$alpha),1e-6)
  wp<-esg_regime_probabilities(b$history_prob,b$next_prob,1,min_ess=80)
  mixture<-esg_match_scenarios(E,wp$probability,b$forecast)
  mixfit<-esg_cvar_lp(mixture$scenarios,mixture$prob,pr$alpha,pr$cap,
    guard_scenarios=base$scenarios,guard_prob=base$prob,guard_limit=lim)
  scenario<-list();prob<-list();diag<-list()
  # Include the baseline as an explicit component of the ambiguity set.
  scenario[[1]]<-base$scenarios;prob[[1]]<-base$prob
  for(s in seq_len(ncol(b$history_prob))){
    q<-b$history_prob[,s];q<-q/sum(q)
    # A fixed 20% uniform floor prevents an empty state dominating the LP.
    q<-.8*q+.2/m;sc<-esg_match_scenarios(E,q,b$forecast)
    scenario[[s+1]]<-sc$scenarios;prob[[s+1]]<-q
    diag[[s]]<-data.frame(state=s,ESS=1/sum(q^2),tail_ESS=(1-pr$alpha)/sum(q^2),
      covariance_error=sc$covariance_error)
  }
  component<-revision_robust_cvar(scenario,prob,pr$alpha,pr$cap,'components',
    base$scenarios,base$prob,lim)
  robust<-revision_robust_cvar(scenario,prob,pr$alpha,pr$cap,'mixtures',
    base$scenarios,base$prob,lim)
  net<-esg_cvar_lp(base$scenarios,base$prob,pr$alpha,pr$cap,b$network$score,pr$network_lambda,
    scale,base$scenarios,base$prob,lim)
  both<-revision_robust_cvar(scenario,prob,pr$alpha,pr$cap,'mixtures',
    base$scenarios,base$prob,lim,b$network$score,pr$network_lambda,scale)
  W<-cbind(C=fit0$weights,C_R=mixfit$weights,C_N=net$weights,
    C_WCOMP=component$weights,C_WMIX=robust$weights,C_WMIX_N=both$weights,
    MinVar=esg_minvar(b$forecast,pr$cap),EqualWeight=rep(1/p,p))
  risks<-apply(W,2,function(w)esg_cvar(-as.numeric(base$scenarios%*%w),base$prob,pr$alpha))
  list(weights=W,diagnostics=do.call(rbind,diag),baseline_risk=risks,
    guard_limit=lim,ambiguity_gap=robust$risk-component$risk,
    max_covariance_error=max(vapply(seq_along(scenario),function(k)
      max(abs(crossprod(scenario[[k]],scenario[[k]]*prob[[k]])-b$forecast))/max(diag(b$forecast)),numeric(1))))
}
revision_pairs <- function(R,pairs,B=1999L,block=10L,seed=20261113L) {
  fun<-list(ES95=function(x)esg_cvar(-x,alpha=.95),
    downside_squared=function(x)mean(pmin(x,0)^2),negative_CER=function(x)-252*(mean(x)-2.5*var(x)))
  out<-list();ii<-esg_with_seed(seed,replicate(B,esg_mbb_indices(nrow(R),block)))
  for(pair in pairs)for(nm in names(fun)){
    f<-fun[[nm]];a<-R[,pair[1]];b<-R[,pair[2]];obs<-f(a)-f(b)
    boot<-apply(ii,2,function(ix)f(a[ix])-f(b[ix]))
    out[[length(out)+1]]<-data.frame(proposed=pair[1],benchmark=pair[2],metric=nm,
      t(esg_basic_inference(obs,boot,'less')),block=block,check.names=FALSE)
  };out<-do.call(rbind,out)
  out$p_holm<-ave(out$p_value,out$metric,FUN=function(x)p.adjust(x,'holm'));out
}
