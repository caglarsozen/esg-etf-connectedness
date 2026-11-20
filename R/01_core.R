# All model and portfolio inputs are simple daily decimal returns.
esg_assert<-function(ok,message){if(!isTRUE(ok))stop(message,call.=FALSE);invisible(TRUE)}
esg_seed<-function(key,base=20261104L){
  b<-as.double(utf8ToInt(paste(key,collapse="|")))
  as.integer((base+sum(b*(seq_along(b)%%997+1)))%%2147483000+1)
}
esg_with_seed<-function(seed,expr){
  had<-exists(".Random.seed",.GlobalEnv,inherits=FALSE)
  if(had)old<-get(".Random.seed",.GlobalEnv)
  on.exit(if(had)assign(".Random.seed",old,.GlobalEnv) else
    if(exists(".Random.seed",.GlobalEnv,inherits=FALSE))rm(".Random.seed",envir=.GlobalEnv))
  set.seed(seed);force(expr)
}
esg_hash<-function(x){p<-tempfile();on.exit(unlink(p));saveRDS(x,p,version=2);unname(tools::md5sum(p))}
esg_csv<-function(x,path){dir.create(dirname(path),recursive=TRUE,showWarnings=FALSE)
  write.csv(x,path,row.names=FALSE,na="",fileEncoding="UTF-8");invisible(path)}
esg_save<-function(x,path){dir.create(dirname(path),recursive=TRUE,showWarnings=FALSE)
  p<-tempfile(tmpdir=dirname(path));saveRDS(x,p,version=2)
  if(file.exists(path))unlink(path)
  esg_assert(file.rename(p,path),paste("Cannot save",path));invisible(path)}
esg_need<-function(pkgs){
  no<-pkgs[!vapply(pkgs,requireNamespace,logical(1),quietly=TRUE)]
  if(length(no))stop("Missing packages: ",paste(no,collapse=", "),
                    ". Run source('install_dependencies.R').",call.=FALSE)
}
esg_spd<-function(S){
  S<-as.matrix(S);esg_assert(nrow(S)==ncol(S)&&all(is.finite(S)),"Invalid covariance.")
  S<-(S+t(S))/2;e<-eigen(S,symmetric=TRUE);floor<-max(max(e$values),1e-14)*1e-10
  esg_assert(min(e$values)>=-max(1e-10,100*floor),"Materially indefinite covariance.")
  if(min(e$values)<floor)S<-e$vectors%*%(pmax(e$values,floor)*t(e$vectors))
  (S+t(S))/2
}
esg_root<-function(S,inverse=FALSE){
  e<-eigen(esg_spd(S),symmetric=TRUE)
  v<-if(inverse)1/sqrt(e$values) else sqrt(e$values)
  e$vectors%*%(v*t(e$vectors))
}
esg_config<-function(mode="audit",root=getwd()){
  esg_assert(mode%in%c("audit","pilot","full"),"Mode must be audit/pilot/full.")
  full<-mode=="full"
  list(mode=mode,root=normalizePath(root,mustWork=TRUE),
    output=file.path(root,paste0("outputs/",toupper(mode))),seed=20261104L,
    price_file=file.path(root,"input","adjusted_prices.csv"),returns_file="",
    from=as.Date("2015-01-01"),through=as.Date("2026-06-25"),
    tickers=c("ESGU","SUSA","ESGV","ICLN","TAN","QCLN","SPY","XLE","XLK","GLD","TLT"),
    blocks=c(ESGU="ESG",SUSA="ESG",ESGV="ESG",ICLN="CleanEnergy",TAN="CleanEnergy",
      QCLN="CleanEnergy",SPY="BroadMarket",XLE="FossilEnergy",XLK="Technology",
      GLD="Gold",TLT="Bonds"),
    covariance=if(full)"factor_sv" else "ewma",ewma_lambda=.94,ewma_warmup=60L,
    pool_window=1260L,factors=3L,factor_candidates=if(full)1:4 else 2:3,
    factor_selection_origins=if(full)3L else 2L,
    fsv_chains=if(full)4L else 2L,fsv_draws=if(full)8000L else 500L,
    fsv_burnin=if(full)3000L else 300L,fsv_oos_chains=2L,
    fsv_oos_draws=if(full)4000L else 500L,fsv_oos_burnin=if(full)1500L else 300L,
    posterior_paths=if(full)200L else 20L,max_rhat=1.01,min_ess=400,strict_mcmc=full,
    states=3L,state_candidates=if(full)2:4 else 2:3,
    hmm_starts=if(full)8L else 3L,hmm_maxit=400L,hmm_tol=1e-7,
    hmm_var_floor=.025,hmm_min_occupancy=15,
    var_lags=c(1L,2L),var_penalties=c(.001,.01,.1,1),
    var_min_obs=if(full)80L else 40L,var_radius=.995,
    horizons=c(1L,5L,10L,20L),horizon=10L,
    bootstrap_B=if(full)999L else 99L,block_lengths=if(full)c(5L,10L,20L) else 10L,
    train_fraction=.65,development_fraction=.50,min_train=if(full)504L else 252L,
    validation_origins=if(full)8L else 4L,
    alpha=.95,alpha_grid=if(full)c(.95,.975,.99) else .95,
    regime_grid=c(0,.25,.5,1),network_grid=c(0,.01,.03,.1),
    min_scenario_ess=80,weight_cap=.35,guard_fraction=.05,
    costs_bps=if(full)c(0,5,10,25) else c(0,10),primary_cost_bps=10,
    utility_gamma=5,utility_margin=.002,
    include_dcc=full,run_no_esgv=full,run_channel_controls=full,
    run_factor_sensitivity=full,run_state_sensitivity=full,
    holdings_file=file.path(root,"input","holdings.csv"),
    external_stress_file=file.path(root,"input","stress_index.csv"))
}
esg_check_panel<-function(R,dates){
  esg_assert(is.matrix(R)&&is.numeric(R)&&nrow(R)==length(dates)&&ncol(R)>=2,
             "Invalid panel dimensions.")
  esg_assert(all(is.finite(R))&&all(R> -1),"Simple decimal returns must be finite and > -1.")
  esg_assert(!is.null(colnames(R))&&!anyDuplicated(colnames(R)),"Unique asset names required.")
  esg_assert(!anyNA(dates)&&!anyDuplicated(dates)&&all(diff(as.Date(dates))>0),
             "Sorted unique dates required.");invisible(TRUE)
}
esg_read_panel<-function(cfg){
  if(nzchar(cfg$returns_file)){
    d<-read.csv(cfg$returns_file,check.names=FALSE);dates<-as.Date(d[[1]])
    R<-as.matrix(d[,-1,drop=FALSE]);storage.mode(R)<-"double"
    keep<-dates>=cfg$from&dates<=cfg$through;R<-R[keep,,drop=FALSE];dates<-dates[keep]
    esg_check_panel(R,dates);return(list(R=R,dates=dates,source=cfg$returns_file))
  }
  if(file.exists(cfg$price_file))d<-read.csv(cfg$price_file,check.names=FALSE) else{
    esg_need("quantmod")
    xs<-lapply(cfg$tickers,function(a){
      x<-quantmod::getSymbols(a,src="yahoo",from=cfg$from,to=cfg$through+1,auto.assign=FALSE)
      x<-quantmod::Ad(x);colnames(x)<-a;x})
    px<-Reduce(function(a,b)merge(a,b,all=TRUE),xs)
    d<-data.frame(date=as.Date(zoo::index(px)),as.matrix(px),check.names=FALSE)
    esg_csv(d,cfg$price_file)
  }
  esg_assert(all(c("date",cfg$tickers)%in%names(d)),"Price CSV misses required ETFs.")
  dates<-as.Date(d$date);P<-as.matrix(d[,cfg$tickers,drop=FALSE]);storage.mode(P)<-"double"
  j<-dates>=cfg$from&dates<=cfg$through;dates<-dates[j];P<-P[j,,drop=FALSE]
  good<-is.finite(P)&P>0
  first<-vapply(seq_len(ncol(P)),function(j){a<-which(good[,j]);if(length(a))a[1] else NA_integer_},integer(1))
  last<-vapply(seq_len(ncol(P)),function(j){a<-which(good[,j]);if(length(a))tail(a,1) else NA_integer_},integer(1))
  esg_assert(!anyNA(first)&&!anyNA(last),"An ETF has no valid prices.")
  j<-seq.int(max(first),min(last));P<-P[j,,drop=FALSE];dates<-dates[j]
  j<-rowSums(is.finite(P))>0;P<-P[j,,drop=FALSE];dates<-dates[j]
  esg_assert(all(is.finite(P)&P>0),"Interior price gap: repair source; do not compress days.")
  esg_assert(nrow(P)>cfg$min_train+100,"Sample too short.")
  R<-P[-1,,drop=FALSE]/P[-nrow(P),,drop=FALSE]-1;dates<-dates[-1]
  esg_check_panel(R,dates);list(R=R,dates=dates,source=cfg$price_file)
}
esg_month_ends<-function(dates){z<-format(as.Date(dates),"%Y-%m");which(c(z[-1]!=head(z,-1),TRUE))}
esg_splits<-function(dates,cfg){
  ep<-esg_month_ends(dates)
  cutoff<-ep[which(ep>=max(cfg$min_train,floor(length(dates)*cfg$train_fraction)))[1]]
  esg_assert(is.finite(cutoff)&&cutoff<length(dates)-10,"No OOS period.")
  val<-tail(ep[ep>=max(cfg$min_train,floor(length(dates)*cfg$development_fraction))&
                 ep<cutoff],cfg$validation_origins)
  oos<-ep[ep>=cutoff&ep<length(dates)]
  esg_assert(length(val)>=2&&length(oos)>=2,"Too few origins.")
  list(validation=val,oos=oos,cutoff=cutoff,validation_end=dates[cutoff],
       first_oos_return=dates[cutoff+1])
}
esg_ewma<-function(R,lambda=.94,warmup=60L){
  n<-nrow(R);p<-ncol(R);esg_assert(n>warmup,"EWMA warmup exceeds sample.")
  H<-array(NA_real_,c(n,p,p));eps<-matrix(NA_real_,n,p,dimnames=dimnames(R))
  S<-esg_spd(cov(R[seq_len(warmup),,drop=FALSE]))
  for(tt in (warmup+1):n){
    H[tt,,]<-S;eps[tt,]<-as.numeric(R[tt,,drop=FALSE]%*%esg_root(S,TRUE))
    S<-esg_spd(lambda*S+(1-lambda)*tcrossprod(R[tt,]))
  }
  list(H=H,forecast=S,eps=eps,valid=(warmup+1):n)
}
esg_features<-function(H){
  n<-dim(H)[1];p<-dim(H)[2];X<-matrix(NA_real_,n,p)
  Z<-matrix(NA_real_,n,2,dimnames=list(NULL,c("log_variance","abs_correlation")))
  for(tt in seq_len(n)){
    S<-H[tt,,];if(any(!is.finite(S)))next
    S<-esg_spd(S);X[tt,]<-log(diag(S))
    C<-S/outer(sqrt(diag(S)),sqrt(diag(S)))
    Z[tt,]<-c(mean(X[tt,]),mean(abs(C[upper.tri(C)])))
  };list(X=X,Z=Z)
}
esg_cvar<-function(loss,prob=NULL,alpha=.95){
  loss<-as.numeric(loss);if(is.null(prob))prob<-rep(1/length(loss),length(loss))
  esg_assert(length(loss)==length(prob)&&all(is.finite(loss))&&
    all(is.finite(prob)&prob>=0)&&sum(prob)>0&&alpha>0&&alpha<1,"Invalid CVaR input.")
  prob<-prob/sum(prob);j<-order(loss,decreasing=TRUE);prob<-prob[j];loss<-loss[j]
  used<-pmin(prob,pmax(0,1-alpha-c(0,head(cumsum(prob),-1))))
  sum(used*loss)/(1-alpha)
}
esg_metrics<-function(r,gamma=5,rf_daily=0){
  r<-as.numeric(r);esg_assert(length(r)>1&&all(is.finite(r))&&all(r> -1),"Invalid realized returns.")
  e<-r-rf_daily;v<-sd(e);ds<-sqrt(mean(pmin(e,0)^2))
  w<-c(1,cumprod(1+r))
  c(CAGR=expm1(252*mean(log1p(r))),annual_mean=252*mean(r),
    volatility=sqrt(252)*sd(r),Sharpe=if(v>0)sqrt(252)*mean(e)/v else NA_real_,
    Sortino=if(ds>0)sqrt(252)*mean(e)/ds else NA_real_,max_drawdown=min(w/cummax(w)-1),
    VaR_loss_95=unname(quantile(-r,.95,type=1)),ES_loss_95=esg_cvar(-r,alpha=.95),
    ES_loss_975=esg_cvar(-r,alpha=.975),ES_loss_99=esg_cvar(-r,alpha=.99),
    CER_quadratic=252*(mean(r)-gamma*var(r)/2),
    downside_squared=mean(pmin(r,0)^2),observations=length(r))
}
