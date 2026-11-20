# Two conditionally independent Gaussian emissions. Explicit scaled EM,
# filtering, smoothing and Viterbi; no quantile fallback.
esg_hmm_fb<-function(Z,mu,variance,P,initial){
  Z<-as.matrix(Z);n<-nrow(Z);k<-nrow(mu)
  logB<-vapply(seq_len(k),function(s)
    -.5*rowSums(sweep(sweep(Z,2,mu[s,],"-")^2,2,variance[s,],"/"))-
      .5*sum(log(2*pi*variance[s,])),numeric(n))
  logB<-matrix(logB,n,k);shift<-apply(logB,1,max);B<-exp(logB-shift)
  A<-matrix(0,n,k);scales<-numeric(n)
  a<-initial*B[1,];scales[1]<-sum(a);esg_assert(scales[1]>0,"HMM underflow.")
  A[1,]<-a/scales[1]
  if(n>1)for(tt in 2:n){
    a<-as.numeric(A[tt-1,,drop=FALSE]%*%P)*B[tt,];scales[tt]<-sum(a)
    esg_assert(is.finite(scales[tt])&&scales[tt]>0,"HMM underflow.")
    A[tt,]<-a/scales[tt]
  }
  backward<-matrix(1,n,k)
  if(n>1)for(tt in (n-1):1)
    backward[tt,]<-as.numeric(P%*%(B[tt+1,]*backward[tt+1,]))/scales[tt+1]
  G<-A*backward;G<-G/rowSums(G);transitions<-matrix(0,k,k)
  if(n>1)for(tt in seq_len(n-1)){
    xi<-P*outer(A[tt,],B[tt+1,]*backward[tt+1,])
    transitions<-transitions+xi/sum(xi)
  }
  list(filtered=A,smoothed=G,transitions=transitions,loglik=sum(log(scales)+shift),logB=logB)
}
esg_hmm_viterbi<-function(logB,P,initial){
  n<-nrow(logB);k<-ncol(logB);D<-matrix(-Inf,n,k);psi<-matrix(1L,n,k)
  D[1,]<-log(pmax(initial,1e-300))+logB[1,]
  if(n>1)for(tt in 2:n)for(s in seq_len(k)){
    z<-D[tt-1,]+log(pmax(P[,s],1e-300));psi[tt,s]<-which.max(z)
    D[tt,s]<-max(z)+logB[tt,s]
  }
  state<-integer(n);state[n]<-which.max(D[n,])
  if(n>1)for(tt in (n-1):1)state[tt]<-psi[tt+1,state[tt+1]]
  state
}
esg_hmm_fit<-function(features,k=3L,starts=3L,seed=123L,maxit=400L,tol=1e-7,
                      variance_floor=.025,min_occupancy=15){
  raw<-as.matrix(features)
  esg_assert(nrow(raw)>max(40,k*min_occupancy)&&all(is.finite(raw)),"HMM support insufficient.")
  center<-colMeans(raw);scale<-apply(raw,2,sd)
  esg_assert(all(scale>1e-10),"HMM feature has zero variance.")
  Z<-sweep(sweep(raw,2,center,"-"),2,scale,"/");best<-NULL;failures<-character()
  for(start in seq_len(starts)){
    fit<-tryCatch(esg_with_seed(esg_seed(c("hmm",seed,start)),{
      km<-kmeans(Z,centers=k,nstart=1,iter.max=100);mu<-km$centers
      variance<-matrix(1,k,ncol(Z))
      for(s in seq_len(k)){
        z<-Z[km$cluster==s,,drop=FALSE]
        if(nrow(z)>2)variance[s,]<-pmax(apply(z,2,var),variance_floor)
      }
      P<-matrix(.1/(k-1),k,k);diag(P)<-.9;initial<-rep(1/k,k)
      old<- -Inf;converged<-FALSE
      for(iter in seq_len(maxit)){
        fb<-esg_hmm_fb(Z,mu,variance,P,initial)
        if(is.finite(old)&&abs(fb$loglik-old)<=tol*(1+abs(old))){converged<-TRUE;break}
        old<-fb$loglik;G<-fb$smoothed;mass<-colSums(G)
        esg_assert(all(mass>2),"Collapsed HMM state.")
        mu<-crossprod(G,Z)/mass
        for(s in seq_len(k))
          variance[s,]<-pmax(colSums(sweep(Z,2,mu[s,],"-")^2*G[,s])/mass[s],variance_floor)
        P<-pmax(fb$transitions,1e-8);P<-P/rowSums(P)
        initial<-pmax(G[1,],1e-8);initial<-initial/sum(initial)
      }
      fb<-esg_hmm_fb(Z,mu,variance,P,initial)
      esg_assert(all(colSums(fb$smoothed)>=min_occupancy),"Small HMM state.")
      o<-order(mu[,1],seq_len(k));mu<-mu[o,,drop=FALSE]
      variance<-variance[o,,drop=FALSE];P<-P[o,o,drop=FALSE];initial<-initial[o]
      fb<-esg_hmm_fb(Z,mu,variance,P,initial)
      labels<-if(k==2)c("Calm","Stress") else if(k==3)c("Calm","Transition","Stress") else
        c("Calm",paste0("Middle",seq_len(k-2)),"Stress")
      dimnames(P)<-list(labels,labels)
      colnames(fb$filtered)<-colnames(fb$smoothed)<-labels
      df<-k*(k-1)+(k-1)+2*k*ncol(Z)
      list(mu=mu,variance=variance,P=P,initial=initial,center=center,scale=scale,
        filtered=fb$filtered,smoothed=fb$smoothed,
        viterbi=esg_hmm_viterbi(fb$logB,P,initial),
        state=max.col(fb$filtered,ties.method="first"),labels=labels,
        loglik=fb$loglik,BIC= -2*fb$loglik+df*log(nrow(Z)),occupancy=colSums(fb$smoothed),
        duration=1/(1-diag(P)),converged=converged,iterations=iter,k=k,start=start)
    }),error=function(e){failures<<-c(failures,conditionMessage(e));NULL})
    if(!is.null(fit)&&(is.null(best)||fit$loglik>best$loglik))best<-fit
  }
  esg_assert(!is.null(best),paste("HMM failed:",paste(unique(failures),collapse="; ")))
  best$failed_starts<-failures;best
}
esg_hmm_filter<-function(model,features){
  Z<-sweep(sweep(as.matrix(features),2,model$center,"-"),2,model$scale,"/")
  esg_hmm_fb(Z,model$mu,model$variance,model$P,model$initial)$filtered
}
