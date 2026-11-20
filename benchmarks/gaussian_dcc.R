# Gaussian two-step DCC(1,1) QMLE. This is not rmgarch's multivariate-t DCC.
# All parameters are re-estimated on the prefix ending at the signal date.
dcc_simplex <- function(x){v<-exp(c(x,0)-max(c(x,0)));.999*v[1:2]/sum(v)}
dcc_start <- function(a,b)log(c(a,b)/(.999-a-b))
dcc_garch_path <- function(y,theta) {
  ab<-dcc_simplex(theta[2:3]);omega<-exp(theta[1]);n<-length(y)
  h<-numeric(n);h[1]<-mean(head(y^2,min(60,n)))
  for(t in 2:n)h[t]<-omega+ab[1]*y[t-1]^2+ab[2]*h[t-1]
  list(h=h,next_h=omega+ab[1]*y[n]^2+ab[2]*h[n],omega=omega,alpha=ab[1],beta=ab[2])
}
dcc_correlation_path <- function(Z,theta,last_only=TRUE) {
  ab<-dcc_simplex(theta);n<-nrow(Z);p<-ncol(Z)
  Qbar<-esg_spd(crossprod(Z)/n);Q<-Qbar;loss<-0
  for(t in seq_len(n)){
    C<-cov2cor(Q);U<-chol(C);wh<-backsolve(U,Z[t,],transpose=TRUE)
    if(t>60)loss<-loss+2*sum(log(diag(U)))+sum(wh^2)
    Q<-(1-sum(ab))*Qbar+ab[1]*tcrossprod(Z[t,])+ab[2]*Q
  }
  list(loss=loss/(n-60),forecast=cov2cor(Q),a=ab[1],b=ab[2])
}
revision_dcc <- function(R,window=756L) {
  R<-tail(R,window)*100;n<-nrow(R);p<-ncol(R)
  esg_assert(n>100,'DCC sample too short.');Z<-matrix(0,n,p);nextv<-numeric(p);diag<-list()
  for(j in seq_len(p)){
    y<-R[,j];s2<-mean(y^2)
    fn<-function(x){v<-dcc_garch_path(y,x)$h
      if(any(!is.finite(v))||min(v)<=0)return(1e10)
      mean((log(v)+y^2/v)[61:n])}
    fits<-lapply(list(c(.05,.9),c(.1,.8)),function(ab)
      optim(c(log(s2*(1-sum(ab))),dcc_start(ab[1],ab[2])),fn,method='BFGS',control=list(maxit=300,reltol=1e-8)))
    good<-which(vapply(fits,function(x)x$convergence==0&&is.finite(x$value),logical(1)))
    esg_assert(length(good)>0,paste('GARCH convergence failed',colnames(R)[j]))
    f<-fits[[good[which.min(vapply(fits[good],function(x)x$value,numeric(1)))]]]
    v<-dcc_garch_path(y,f$par);Z[,j]<-y/sqrt(v$h);nextv[j]<-v$next_h
    diag[[j]]<-data.frame(asset=colnames(R)[j],omega=v$omega,alpha=v$alpha,beta=v$beta,
      objective=f$value,convergence=f$convergence)
  }
  fn<-function(x){tryCatch(dcc_correlation_path(Z,x)$loss,error=function(e)1e10)}
  fits<-lapply(list(c(.02,.95),c(.05,.8)),function(ab)
    optim(dcc_start(ab[1],ab[2]),fn,method='BFGS',control=list(maxit=200,reltol=1e-7)))
  good<-which(vapply(fits,function(x)x$convergence==0&&is.finite(x$value)&&x$value<1e9,logical(1)))
  esg_assert(length(good)>0,'DCC correlation optimization failed.')
  f<-fits[[good[which.min(vapply(fits[good],function(x)x$value,numeric(1)))]]]
  a<-dcc_correlation_path(Z,f$par);S<-a$forecast*outer(sqrt(nextv),sqrt(nextv))/10000
  list(forecast=esg_spd(S),univariate=do.call(rbind,diag),a=a$a,b=a$b,
    objective=f$value,convergence=f$convergence)
}
