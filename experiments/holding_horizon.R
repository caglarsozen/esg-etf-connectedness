# Coefficients are indexed by the LAGGED state, as in esg_lag_packets().
# Therefore update the current-state dynamics BEFORE the state transition.
horizon_own_fit<-function(packet,rows=seq_len(nrow(packet$Y)),lambda=.1){
 q<-ncol(packet$Y);B<-matrix(0,1+q,q)
 for(j in seq_len(q)){
  x<-packet$Z[rows,c(1,j+1),drop=FALSE];y<-packet$Y[rows,j]
  bb<-solve(crossprod(x)+diag(c(0,length(rows)*lambda),2),crossprod(x,y))
  slope<-max(-.995,min(.995,bb[2]));B[1,j]<-mean(y)-slope*mean(x[,2]);B[j+1,j]<-slope
 }
 E<-packet$Y[rows,,drop=FALSE]-packet$Z[rows,,drop=FALSE]%*%B
 list(B=B,A=list(t(B[-1,,drop=FALSE])),Sigma=esg_spd(crossprod(E)/nrow(E)),n=length(rows))
}
horizon_forecast<-function(S,packet,fits,probability,transition,days=21L){
 q<-nrow(S);K<-length(fits);probability<-probability/sum(probability)
 stopifnot(days>=1,nrow(transition)==K,max(abs(rowSums(transition)-1))<1e-8)
 x<-(log(diag(S))-packet$center)/packet$scale
 m<-lapply(probability,function(a)a*x);M<-lapply(probability,function(a)a*tcrossprod(x))
 probs<-probability;v<-matrix(NA_real_,days,q);v[1,]<-diag(S)
 if(days>1L)for(h in 2:days){
  next_m<-lapply(seq_len(K),function(s)rep(0,q));next_M<-lapply(seq_len(K),function(s)matrix(0,q,q))
  next_p<-as.numeric(probs%*%transition)
  for(s in seq_len(K)){
   A<-fits[[s]]$A[[1]];cc<-as.numeric(fits[[s]]$B[1,]);am<-as.numeric(A%*%m[[s]])
   lm<-probs[s]*cc+am
   lM<-A%*%M[[s]]%*%t(A)+tcrossprod(cc,am)+tcrossprod(am,cc)+probs[s]*(tcrossprod(cc)+fits[[s]]$Sigma)
   for(j in seq_len(K)){next_m[[j]]<-next_m[[j]]+transition[s,j]*lm;next_M[[j]]<-next_M[[j]]+transition[s,j]*lM}
  }
  probs<-next_p;m<-next_m;M<-next_M;v[h,]<-0
  for(s in seq_len(K))if(probs[s]>1e-12){
   mu<-m[[s]]/probs[s];V<-M[[s]]/probs[s]-tcrossprod(mu)
   lv<-packet$center+packet$scale*mu+.5*packet$scale^2*pmax(diag(V),0)
   esg_assert(all(is.finite(lv))&&max(lv)<10,'Unstable log-variance forecast; no clipping or fallback.')
   v[h,]<-v[h,]+probs[s]*exp(lv)
  }
 }
 # Joint first/second moments are exact for fitted switching linear dynamics.
 # Exponential moment is exact only for K=1; K>1 uses a conditional Gaussian approximation.
 avg<-colMeans(v);H<-cov2cor(S)*outer(sqrt(avg),sqrt(avg))
 list(H=esg_spd(H),variance_path=v,final_state_probability=probs,final_joint_mean=m,final_joint_second=M)
}
horizon_bundle_forecasts<-function(bundle,days=21L){
 pk<-bundle$packet;k<-length(bundle$var)
 own<-horizon_own_fit(pk);full<-esg_var_fit_rows(pk$Y,pk$Z,1L,lambda=.1)
 own_states<-lapply(seq_len(k),function(s)horizon_own_fit(pk,which(pk$state==s)))
 list(H_OWN=horizon_forecast(bundle$forecast,pk,list(own),1,matrix(1,1,1),days),
 H_VAR=horizon_forecast(bundle$forecast,pk,list(full),1,matrix(1,1,1),days),
 H_STATE_OWN=horizon_forecast(bundle$forecast,pk,own_states,bundle$next_prob,bundle$hmm$P,days),
 H_STATE_VAR=horizon_forecast(bundle$forecast,pk,bundle$var,bundle$next_prob,bundle$hmm$P,days))
}
