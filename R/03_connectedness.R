# Preserve original daily lag packets before selecting states or bootstrapping.
esg_lag_packets<-function(X,state,p=1L,dates=NULL){
  X<-as.matrix(X);n<-nrow(X);q<-ncol(X)
  esg_assert(length(state)==n&&n>p+10,"Invalid lag input.")
  center<-colMeans(X);scale<-apply(X,2,sd)
  esg_assert(all(is.finite(scale)&scale>1e-10),"Degenerate VAR variable.")
  Xs<-sweep(sweep(X,2,center,"-"),2,scale,"/");tt<-(p+1):n
  Z<-cbind(intercept=1,do.call(cbind,lapply(seq_len(p),function(j)Xs[tt-j,,drop=FALSE])))
  list(Y=Xs[tt,,drop=FALSE],Z=Z,state=state[tt-1],time=tt,
    lag_time=sapply(seq_len(p),function(j)tt-j),p=p,center=center,scale=scale,
    dates=if(is.null(dates))NULL else as.Date(dates)[tt],assets=colnames(X),q=q)
}
esg_companion<-function(A){
  p<-length(A);q<-nrow(A[[1]])
  if(p==1)return(A[[1]])
  rbind(do.call(cbind,A),cbind(diag(q*(p-1)),matrix(0,q*(p-1),q)))
}
esg_var_fit_rows<-function(Y,Z,p,lambda=.01,radius=.995){
  Y<-as.matrix(Y);Z<-as.matrix(Z);n<-nrow(Y);q<-ncol(Y);m<-ncol(Z)
  esg_assert(n>q+5&&m==1+p*q&&all(is.finite(Y))&&all(is.finite(Z)),"Invalid VAR support.")
  Q<-crossprod(Z)+diag(c(0,rep(n*lambda,m-1)),m)
  B<-solve(Q,crossprod(Z,Y))
  getA<-function(B)lapply(seq_len(p),function(j)t(B[1+(j-1)*q+seq_len(q),,drop=FALSE]))
  A<-getA(B);rho0<-max(Mod(eigen(esg_companion(A),only.values=TRUE)$values))
  rho<-rho0;shrink<-1
  while(rho>=radius){
    B[-1,]<-B[-1,,drop=FALSE]*.97;shrink<-shrink*.97
    B[1,]<-colMeans(Y)-as.numeric(colMeans(Z[,-1,drop=FALSE])%*%B[-1,,drop=FALSE])
    A<-getA(B);rho<-max(Mod(eigen(esg_companion(A),only.values=TRUE)$values))
    esg_assert(shrink>1e-10,"VAR stability adjustment failed.")
  }
  E<-Y-Z%*%B;E<-sweep(E,2,colMeans(E),"-");S<-esg_spd(crossprod(E)/n)
  df<-sum(diag(solve(Q,crossprod(Z))))
  bic<-n*as.numeric(determinant(S,logarithm=TRUE)$modulus)+log(n)*q*df
  lb<-vapply(seq_len(q),function(j)Box.test(E[,j],lag=min(10L,floor(n/5)),
                     type="Ljung-Box")$p.value,numeric(1))
  list(B=B,A=A,Sigma=S,residual=E,n=n,p=p,lambda=lambda,radius_before=rho0,
    radius_after=rho,stability_multiplier=shrink,effective_df=df,BIC=bic,
    residual_ljung_box_p=lb)
}
esg_var_from_packets<-function(packet,k,min_obs=40L,lambda=.01,radius=.995,
                               row_index=seq_len(nrow(packet$Y))){
  Y<-packet$Y[row_index,,drop=FALSE];Z<-packet$Z[row_index,,drop=FALSE]
  states<-packet$state[row_index]
  lapply(seq_len(k),function(s){
    j<-which(states==s)
    esg_assert(length(j)>=max(min_obs,packet$q+6),paste("Too few genuine lag pairs in state",s))
    esg_var_fit_rows(Y[j,,drop=FALSE],Z[j,,drop=FALSE],packet$p,lambda,radius)
  })
}
esg_var_select<-function(X,state,k,cfg){
  rows<-list();best<-NULL
  for(p in cfg$var_lags)for(lambda in cfg$var_penalties){
    z<-tryCatch({
      pk<-esg_lag_packets(X,state,p)
      f<-esg_var_from_packets(pk,k,cfg$var_min_obs,lambda,cfg$var_radius)
      list(packet=pk,fit=f,score=sum(vapply(f,function(a)a$BIC,numeric(1))))
    },error=function(e)NULL)
    rows[[length(rows)+1]]<-data.frame(p=p,lambda=lambda,BIC_proxy=if(is.null(z))Inf else z$score)
    if(!is.null(z)&&(is.null(best)||z$score<best$score)){best<-z;best$p<-p;best$lambda<-lambda}
  }
  esg_assert(!is.null(best),"No supported daily regime VAR.")
  best$table<-do.call(rbind,rows);best
}
esg_fevd<-function(fit,horizon=10L){
  S<-fit$Sigma;q<-nrow(S);p<-length(fit$A);Psi<-vector("list",horizon);Psi[[1]]<-diag(q)
  if(horizon>1)for(h in 1:(horizon-1)){
    v<-matrix(0,q,q)
    for(l in seq_len(min(p,h)))v<-v+fit$A[[l]]%*%Psi[[h-l+1]]
    Psi[[h+1]]<-v
  }
  numerator<-matrix(0,q,q);denominator<-numeric(q)
  for(Phi in Psi){M<-Phi%*%S;numerator<-numerator+sweep(M^2,2,diag(S),"/")
    denominator<-denominator+diag(M%*%t(Phi))}
  raw<-sweep(numerator,1,denominator,"/");theta<-raw/rowSums(raw)
  esg_assert(all(is.finite(theta))&&max(abs(rowSums(theta)-1))<1e-8,"Invalid GFEVD.")
  theta
}
esg_conn<-function(theta,assets){
  off<-theta;diag(off)<-0;FROM<-100*rowSums(off);TO<-100*colSums(off)
  list(theta=theta,TCI=100*sum(off)/length(assets),
       asset=data.frame(asset=assets,FROM=FROM,TO=TO,NET=TO-FROM))
}
esg_block<-function(theta,assets,blocks){
  g<-unname(blocks[assets]);esg_assert(!anyNA(g),"Unclassified asset.");out<-list()
  for(source in unique(g))for(target in unique(g)){
    ij<-expand.grid(i=which(g==target),j=which(g==source))
    ij<-ij[ij$i!=ij$j,,drop=FALSE];if(!nrow(ij))next
    z<-theta[cbind(ij$i,ij$j)]
    out[[length(out)+1]]<-data.frame(source=source,target=target,mean_pair_share=100*mean(z),
      sum_share=100*sum(z),target_average_total_share=100*sum(z)/sum(g==target),
      system_share=100*sum(z)/length(assets),pairs=length(z))
  };do.call(rbind,out)
}
esg_conn_tables<-function(fits,labels,assets,blocks,horizons=c(1L,10L)){
  tci<-asset<-pair<-block<-diagnostics<-list()
  for(s in seq_along(fits)){
    f<-fits[[s]]
    diagnostics[[s]]<-data.frame(regime=labels[s],genuine_lag_pairs=f$n,lag=f$p,lambda=f$lambda,
      radius_before=f$radius_before,radius_after=f$radius_after,stability_multiplier=f$stability_multiplier,
      min_Ljung_Box_p=min(f$residual_ljung_box_p),effective_df=f$effective_df)
    for(h in horizons){
      cc<-esg_conn(esg_fevd(f,h),assets)
      tci[[length(tci)+1]]<-data.frame(regime=labels[s],horizon=h,TCI=cc$TCI,genuine_lag_pairs=f$n)
      asset[[length(asset)+1]]<-cbind(regime=labels[s],horizon=h,cc$asset)
      ij<-expand.grid(affected=assets,shock_from=assets,stringsAsFactors=FALSE)
      pair[[length(pair)+1]]<-cbind(regime=labels[s],horizon=h,ij,share=100*as.vector(cc$theta))
      block[[length(block)+1]]<-cbind(regime=labels[s],horizon=h,esg_block(cc$theta,assets,blocks))
    }
  }
  list(tci=do.call(rbind,tci),asset=do.call(rbind,asset),pairwise=do.call(rbind,pair),
       block=do.call(rbind,block),diagnostics=do.call(rbind,diagnostics))
}
esg_network_score<-function(fits,regime_prob,assets,blocks,horizon=10L){
  theta<-Reduce("+",Map(function(f,p)esg_fevd(f,horizon)*p,fits,regime_prob))
  g<-unname(blocks[assets])
  external<-vapply(seq_along(assets),function(j)mean(theta[g!=g[j],j]),numeric(1))
  width<-diff(range(external))
  score<-if(width<1e-12)rep(0,length(assets)) else (external-min(external))/width
  list(score=setNames(score,assets),external_raw=external,theta=theta)
}
esg_conn_statistics<-function(fits,labels,assets,blocks,horizon=10L){
  v<-numeric()
  for(s in seq_along(fits)){
    cc<-esg_conn(esg_fevd(fits[[s]],horizon),assets)
    v[paste("TCI",labels[s],sep="|")]<-cc$TCI
    v[paste("NET",labels[s],assets,sep="|")]<-cc$asset$NET
    ce<-which(unname(blocks[assets])=="CleanEnergy")
    if(length(ce))v[paste("CE_NET",labels[s],sep="|")]<-mean(cc$asset$NET[ce])
    b<-esg_block(cc$theta,assets,blocks)
    v[paste("BLOCK",labels[s],b$source,b$target,sep="|")]<-b$mean_pair_share
  }
  if(all(c("Calm","Stress")%in%labels)){
    v["D_TCI"]<-v["TCI|Stress"]-v["TCI|Calm"]
    if("CE_NET|Stress"%in%names(v))v["D_CE_NET"]<-v["CE_NET|Stress"]-v["CE_NET|Calm"]
    for(a in assets)v[paste0("D_NET|",a)]<-v[paste("NET","Stress",a,sep="|")]-
      v[paste("NET","Calm",a,sep="|")]
    keys<-c("BLOCK|Stress|Technology|CleanEnergy","BLOCK|Calm|Technology|CleanEnergy",
            "BLOCK|Stress|FossilEnergy|CleanEnergy","BLOCK|Calm|FossilEnergy|CleanEnergy")
    if(all(keys%in%names(v))){
      v["D_TECH_CE"]<-v[keys[1]]-v[keys[2]];v["D_FOSSIL_CE"]<-v[keys[3]]-v[keys[4]]
      v["D_TECH_MINUS_FOSSIL"]<-v["D_TECH_CE"]-v["D_FOSSIL_CE"]
      v["STRESS_TECH_MINUS_FOSSIL"]<-v[keys[1]]-v[keys[3]]
    }
  };v
}
esg_mbb_indices<-function(n,block=10L){
  esg_assert(n>=block&&block>=1,"Insufficient bootstrap support.")
  starts<-sample.int(n-block+1,ceiling(n/block),replace=TRUE)
  head(unlist(lapply(starts,function(a)seq.int(a,length.out=block))),n)
}
esg_basic_inference<-function(observed,boot,alternative="two.sided"){
  observed<-as.numeric(observed);b<-as.numeric(boot);b<-b[is.finite(b)]
  esg_assert(length(b)>=20,"Too few valid bootstrap draws.")
  root<-b-observed;ci<-observed-as.numeric(quantile(root,c(.975,.025),names=FALSE))
  pg<-(sum(root>=observed)+1)/(length(b)+1);pl<-(sum(root<=observed)+1)/(length(b)+1)
  pv<-if(alternative=="greater")pg else if(alternative=="less")pl else min(1,2*min(pg,pl))
  c(estimate=observed,lower=ci[1],upper=ci[2],p_value=pv,draws=length(b))
}
esg_conn_bootstrap<-function(packet,fits,labels,assets,blocks,cfg,key="main",
                             B=cfg$bootstrap_B,block=10L){
  observed<-esg_conn_statistics(fits,labels,assets,blocks,cfg$horizon)
  boot<-matrix(NA_real_,B,length(observed),dimnames=list(NULL,names(observed)))
  errors<-character()
  for(b in seq_len(B)){
    z<-tryCatch(esg_with_seed(esg_seed(c(key,"conn_boot",block,b),cfg$seed),{
      ix<-esg_mbb_indices(nrow(packet$Y),block)
      ff<-esg_var_from_packets(packet,length(labels),cfg$var_min_obs,
                               fits[[1]]$lambda,cfg$var_radius,ix)
      esg_conn_statistics(ff,labels,assets,blocks,cfg$horizon)
    }),error=function(e){errors<<-c(errors,conditionMessage(e));NULL})
    if(!is.null(z))boot[b,]<-z[names(observed)]
  }
  esg_assert(sum(complete.cases(boot))>=max(20,.9*B),"Too many failed connection bootstraps.")
  tab<-do.call(rbind,lapply(seq_along(observed),function(j){
    nm<-names(observed)[j];alt<-if(nm%in%c("D_TCI","D_CE_NET","D_TECH_CE",
      "D_FOSSIL_CE","D_TECH_MINUS_FOSSIL","STRESS_TECH_MINUS_FOSSIL"))"greater" else "two.sided"
    data.frame(statistic=nm,t(esg_basic_inference(observed[j],boot[,j],alt)),
      alternative=alt,block_length=block,scope="conditional on latent path and state labels",
      check.names=FALSE)
  }))
  family<-ifelse(tab$statistic%in%c("D_TCI","D_CE_NET","D_TECH_MINUS_FOSSIL"),"primary",
                 ifelse(grepl("^D_NET[|]",tab$statistic),"asset_changes","additional"))
  tab$p_holm<-ave(tab$p_value,family,FUN=function(x)p.adjust(x,"holm"))
  role<-if(all(c("CE_NET|Calm","CE_NET|Stress")%in%colnames(boot)))
    data.frame(joint_sign_frequency=mean(boot[,"CE_NET|Calm"]<0&boot[,"CE_NET|Stress"]>0,na.rm=TRUE),
      scope="bootstrap frequency, not posterior probability") else NULL
  list(table=tab,draws=boot,role_switch=role,failures=errors)
}
