esg_project_weights<-function(v,cap=.35){
  n<-length(v);esg_assert(all(is.finite(v))&&n*cap>=1-1e-12,"Infeasible weight cap.")
  lo<-min(v)-cap-1;hi<-max(v)+1
  for(it in 1:100){mid<-(lo+hi)/2
    if(sum(pmin(cap,pmax(0,v-mid)))>1)lo<-mid else hi<-mid}
  w<-pmin(cap,pmax(0,v-(lo+hi)/2))
  esg_assert(abs(sum(w)-1)<1e-8,"Weight projection failed.");as.numeric(w)
}
esg_minvar<-function(S,cap=.35){
  esg_need("quadprog");S<-esg_spd(S);n<-ncol(S)
  esg_assert(n*cap>=1-1e-12,"Infeasible QP cap.")
  q<-quadprog::solve.QP(2*S/max(diag(S)),rep(0,n),
     cbind(rep(1,n),diag(n),-diag(n)),c(1,rep(0,n),rep(-cap,n)),meq=1)
  w<-as.numeric(q$solution)
  esg_assert(abs(sum(w)-1)<1e-6&&min(w)> -1e-6&&max(w)<cap+1e-6,"QP constraint error.")
  esg_project_weights(w,cap)
}
esg_cvar_lp<-function(scenarios,prob=NULL,alpha=.95,cap=.35,
  network_score=NULL,network_lambda=0,risk_scale=NULL,
  guard_scenarios=NULL,guard_prob=NULL,guard_limit=NULL){
  esg_need("lpSolve");R<-as.matrix(scenarios);m<-nrow(R);n<-ncol(R)
  if(is.null(prob))prob<-rep(1/m,m)
  esg_assert(all(is.finite(R))&&length(prob)==m&&all(is.finite(prob)&prob>=0)&&
               sum(prob)>0&&n*cap>=1&&alpha>0&&alpha<1,"Invalid scenario LP.")
  prob<-prob/sum(prob)
  if(is.null(network_score))network_score<-rep(0,n)
  if(is.null(risk_scale))risk_scale<-max(esg_cvar(-rowMeans(R),prob,alpha),1e-6)
  esg_assert(length(network_score)==n&&all(is.finite(network_score))&&network_lambda>=0,
             "Invalid network penalty.")
  R<-R*100;risk_scale<-risk_scale*100
  # Every LP variable is nonnegative. VaR is z_plus-z_minus, not restricted to >=0.
  obj<-c(network_lambda*risk_scale*network_score,1,-1,prob/(1-alpha))
  tailrows<-n+1+seq_len(m)
  trip<-rbind(cbind(1,seq_len(n),1),cbind(1+seq_len(n),seq_len(n),1),
    cbind(rep(tailrows,times=n),rep(seq_len(n),each=m),-as.vector(R)),
    cbind(tailrows,n+1,-1),cbind(tailrows,n+2,1),
    cbind(tailrows,n+2+seq_len(m),-1))
  rhs<-c(1,rep(cap,n),rep(0,m))
  if(!is.null(guard_scenarios)){
    G<-as.matrix(guard_scenarios)*100;mg<-nrow(G)
    if(is.null(guard_prob))guard_prob<-rep(1/mg,mg)
    esg_assert(ncol(G)==n&&all(is.finite(G))&&length(guard_prob)==mg&&
      all(is.finite(guard_prob)&guard_prob>=0)&&sum(guard_prob)>0&&
      length(guard_limit)==1&&is.finite(guard_limit),"Invalid CVaR guard.")
    guard_prob<-guard_prob/sum(guard_prob);last<-length(obj)
    rows<-length(rhs)+seq_len(mg);budget<-length(rhs)+mg+1
    trip<-rbind(trip,
      cbind(rep(rows,times=n),rep(seq_len(n),each=mg),-as.vector(G)),
      cbind(rows,last+1,-1),cbind(rows,last+2,1),cbind(rows,last+2+seq_len(mg),-1),
      cbind(budget,last+1,1),cbind(budget,last+2,-1),
      cbind(budget,last+2+seq_len(mg),guard_prob/(1-alpha)))
    obj<-c(obj,rep(0,mg+2));rhs<-c(rhs,rep(0,mg),guard_limit*100)
  }
  fit<-lpSolve::lp("min",objective.in=obj,const.dir=c("=",rep("<=",length(rhs)-1)),
                    const.rhs=rhs,dense.const=trip,scale=196)
  esg_assert(fit$status==0,paste("CVaR solver failed; no fallback. Status:",fit$status))
  w<-fit$solution[seq_len(n)]
  error<-max(abs(sum(w)-1),pmax(-w,0),pmax(w-cap,0))
  esg_assert(error<1e-5,"CVaR constraint error.");w<-esg_project_weights(w,cap)
  value<-esg_cvar(-as.numeric((R/100)%*%w),prob,alpha)
  achieved<-value+network_lambda*(risk_scale/100)*sum(network_score*w)
  esg_assert(abs(achieved-fit$objval/100)<max(1e-6,abs(achieved)*1e-4),"LP objective mismatch.")
  guarded<-if(is.null(guard_scenarios))NA_real_ else
    esg_cvar(-as.numeric((G/100)%*%w),guard_prob,alpha)
  if(!is.null(guard_scenarios))esg_assert(guarded<=guard_limit+1e-6,"CVaR guard violated.")
  list(weights=w,CVaR=value,objective=achieved,guard_CVaR=guarded,
       max_constraint_error=error,lp_status=fit$status)
}
esg_regime_probabilities<-function(history,next_probability,strength=0,min_ess=80){
  G<-as.matrix(history);m<-nrow(G)
  esg_assert(ncol(G)==length(next_probability)&&all(is.finite(G))&&all(G>=0)&&
               strength>=0&&strength<=1,"Invalid regime weights.")
  q<-next_probability/sum(next_probability);mass<-colSums(G)
  esg_assert(all(mass>1e-8),"No history for a forecast state.")
  weighted<-as.numeric(G%*%(q/mass));weighted<-weighted/sum(weighted)
  uniform<-rep(1/m,m);effective<-strength;target<-min(min_ess,.9*m)
  while(effective>1e-8){
    prob<-(1-effective)*uniform+effective*weighted
    if(1/sum(prob^2)>=target)break
    effective<-effective*.9
  }
  prob<-(1-effective)*uniform+effective*weighted
  list(probability=prob/sum(prob),effective_strength=effective,ESS=1/sum(prob^2))
}
esg_match_scenarios<-function(standardized_residuals,prob,target_covariance){
  E<-as.matrix(standardized_residuals);prob<-prob/sum(prob)
  esg_assert(nrow(E)==length(prob)&&all(is.finite(E)),"Invalid FHS pool.")
  centered<-sweep(E,2,colSums(E*prob),"-")
  C<-crossprod(centered,centered*prob)
  R<-centered%*%esg_root(C,TRUE)%*%esg_root(target_covariance)
  error<-max(abs(crossprod(R,R*prob)-target_covariance))/max(diag(target_covariance))
  esg_assert(error<1e-6,"FHS covariance mismatch.")
  esg_assert(all(R> -1),"A simulated simple return is below -100%; no tail clipping.")
  list(scenarios=R,prob=prob,covariance_error=error)
}
esg_trade<-function(pre,target,cost_bps=10){
  esg_assert(length(pre)==length(target)&&all(pre>=-1e-10)&&sum(pre)<=1+1e-8&&
               abs(sum(target)-1)<1e-8&&all(target>=0),"Invalid trade.")
  rate<-cost_bps/10000;esg_assert(rate>=0&&rate<.1,"Invalid per-dollar cost.")
  fraction<-if(rate==0)1 else uniroot(function(a)
    a+rate*sum(abs(a*target-pre))-1,c(0,1),tol=1e-13)$root
  traded<-sum(abs(fraction*target-pre))
  esg_assert(abs(fraction+rate*traded-1)<1e-10,"Self-financing fee equation failed.")
  list(remaining=fraction,cost_fraction=rate*traded,traded_notional=traded,turnover_half=.5*traded)
}
esg_backtest_targets<-function(R,dates,origins,targets,cost_bps=10){
  esg_check_panel(R,dates)
  esg_assert(length(dim(targets))==3&&dim(targets)[1]==length(origins)&&
    dim(targets)[2]==ncol(R)&&all(diff(origins)>0)&&min(origins)>=1&&max(origins)<nrow(R),
    "Invalid target schedule.")
  strategies<-dimnames(targets)[[3]];ns<-length(strategies)
  esg_assert(ns==dim(targets)[3]&&!anyDuplicated(strategies),"Strategy labels required.")
  returns<-matrix(NA_real_,nrow(R),ns,dimnames=list(NULL,strategies))
  pre<-matrix(0,ns,ncol(R));trades<-list()
  for(k in seq_along(origins)){
    first<-origins[k]+1;last<-if(k<length(origins))origins[k+1] else nrow(R)
    for(j in seq_len(ns)){
      w<-as.numeric(targets[k,,j]);fee<-esg_trade(pre[j,],w,cost_bps)
      trades[[length(trades)+1]]<-data.frame(signal_date=as.character(dates[origins[k]]),
        strategy=strategies[j],turnover_half=fee$turnover_half,
        traded_notional=fee$traded_notional,cost_fraction=fee$cost_fraction,
        initial_purchase=k==1,max_target_weight=max(w))
      for(tt in first:last){
        gross<-sum(w*R[tt,]);esg_assert(gross> -1,"Portfolio bankruptcy.")
        returns[tt,j]<-(if(tt==first)fee$remaining else 1)*(1+gross)-1
        w<-w*(1+R[tt,])/(1+gross)
        esg_assert(abs(sum(w)-1)<1e-8,"Weight drift error.")
      };pre[j,]<-w
    }
  }
  keep<-which(complete.cases(returns));r<-returns[keep,,drop=FALSE]
  list(returns=r,dates=as.Date(dates)[keep],trades=do.call(rbind,trades),
       wealth=apply(r,2,function(z)cumprod(1+z)),targets=targets,origins=origins)
}
esg_performance_table<-function(bt,cost_bps=10,gamma=5){
  do.call(rbind,lapply(colnames(bt$returns),function(s){
    m<-esg_metrics(bt$returns[,s],gamma)
    tr<-bt$trades[bt$trades$strategy==s&!bt$trades$initial_purchase,,drop=FALSE]
    data.frame(strategy=s,cost_bps=cost_bps,t(m),
      mean_turnover=if(nrow(tr))mean(tr$turnover_half) else NA_real_,check.names=FALSE)
  }))
}
esg_portfolio_candidates<-function(bundle,cfg,alpha=cfg$alpha){
  candidates<-expand.grid(regime=cfg$regime_grid,network=cfg$network_grid)
  candidates<-candidates[order(candidates$regime,candidates$network),,drop=FALSE]
  candidates$id<-paste0("r",sprintf("%.8g",candidates$regime),"_n",sprintf("%.8g",candidates$network))
  W<-matrix(NA_real_,nrow(candidates),ncol(bundle$R),dimnames=list(candidates$id,colnames(bundle$R)))
  base<-esg_match_scenarios(bundle$eps,rep(1/nrow(bundle$eps),nrow(bundle$eps)),bundle$forecast)
  fit0<-esg_cvar_lp(base$scenarios,base$prob,alpha,cfg$weight_cap)
  scale<-max(esg_cvar(-rowMeans(base$scenarios),base$prob,alpha),1e-6)
  limit<-fit0$CVaR+cfg$guard_fraction*max(fit0$CVaR,1e-6);diag<-list()
  for(eta in unique(candidates$regime)){
    pr<-esg_regime_probabilities(bundle$history_prob,bundle$next_prob,eta,cfg$min_scenario_ess)
    sc<-esg_match_scenarios(bundle$eps,pr$probability,bundle$forecast)
    for(j in which(candidates$regime==eta)){
      z<-if(eta==0&&candidates$network[j]==0)fit0 else
        esg_cvar_lp(sc$scenarios,sc$prob,alpha,cfg$weight_cap,bundle$network$score,
          candidates$network[j],scale,base$scenarios,base$prob,limit)
      W[j,]<-z$weights
      diag[[length(diag)+1]]<-data.frame(candidate=candidates$id[j],alpha=alpha,
        regime=eta,effective_regime_strength=pr$effective_strength,network=candidates$network[j],
        historical_vectors=nrow(bundle$eps),scenario_ESS=pr$ESS,
        effective_tail_vectors=(1-alpha)*pr$ESS,CVaR=z$CVaR,
        baseline_CVaR=fit0$CVaR,guard_limit=limit,
        realized_guard_CVaR=esg_cvar(-as.numeric(base$scenarios%*%z$weights),base$prob,alpha),
        covariance_match_error=sc$covariance_error,constraint_error=z$max_constraint_error)
    }
  };list(candidates=candidates,weights=W,diagnostics=do.call(rbind,diag))
}
esg_select_ablation<-function(bt,candidates,block=10L){
  loss<-pmin(bt$returns,0)^2;score<-colMeans(loss)
  groups<-split(seq_len(nrow(loss)),ceiling(seq_len(nrow(loss))/block))
  se<-vapply(seq_len(ncol(loss)),function(j){
    b<-vapply(groups,function(ix)mean(loss[ix,j]),numeric(1));if(length(b)>1)sd(b)/sqrt(length(b)) else 0
  },numeric(1))
  arms<-list(C=which(candidates$regime==0&candidates$network==0),
    C_R=which(candidates$network==0),C_N=which(candidates$regime==0),C_R_N=seq_len(nrow(candidates)))
  out<-lapply(names(arms),function(a){
    pool<-arms[[a]];best<-pool[which.min(score[pool])]
    good<-pool[score[pool]<=score[best]+se[best]+1e-18]
    good<-good[order(candidates$regime[good]+candidates$network[good],
                      candidates$network[good],candidates$regime[good])]
    j<-good[1]
    data.frame(arm=a,candidate=candidates$id[j],regime=candidates$regime[j],network=candidates$network[j],
      validation_downside_squared=score[j],validation_SE=se[j])
  })
  list(selected=do.call(rbind,out),all_scores=cbind(candidates,validation_score=score,validation_SE=se))
}
esg_paired_tests<-function(returns,cfg,key="oos",B=cfg$bootstrap_B,block=10L){
  R<-as.matrix(returns)
  contrasts<-list(c("C_R","C"),c("C_N","C"),c("C_R_N","C"),c("C_R_N","MinVar"),
    c("C_R_N","EqualWeight"),c("C_R_N","LowConnectedness"),c("C_R_N","DCC_MinVar"))
  contrasts<-Filter(function(z)all(z%in%colnames(R)),contrasts)
  fs<-list(downside_squared=function(r)mean(pmin(r,0)^2),
    ES95=function(r)esg_cvar(-r,alpha=.95),ES975=function(r)esg_cvar(-r,alpha=.975),
    ES99=function(r)esg_cvar(-r,alpha=.99),
    negative_Sharpe=function(r)-unname(esg_metrics(r)["Sharpe"]),
    negative_CER=function(r)-unname(esg_metrics(r,cfg$utility_gamma)["CER_quadratic"]))
  rows<-list()
  for(z in contrasts)for(metric in names(fs)){
    f<-fs[[metric]];obs<-f(R[,z[1]])-f(R[,z[2]])
    boot<-esg_with_seed(esg_seed(c(key,z,metric,block),cfg$seed),replicate(B,{
      ix<-esg_mbb_indices(nrow(R),block);f(R[ix,z[1]])-f(R[ix,z[2]])}))
    rows[[length(rows)+1]]<-data.frame(proposed=z[1],benchmark=z[2],metric=metric,
      t(esg_basic_inference(obs,boot,"less")),block_length=block,
      interpretation="negative favors proposed",check.names=FALSE)
  }
  out<-do.call(rbind,rows);out$p_holm<-ave(out$p_value,out$metric,FUN=function(x)p.adjust(x,"holm"));out
}
