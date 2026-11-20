source('R/01_core.R');source('experiments/particle_prediction.R')
esg_with_seed(81,{
 p<-4L;N<-100L;L<-array(rnorm(p*2*N,sd=.3),c(p,2,N));h<-matrix(rnorm((p+2)*N,sd=.2),p+2,N)
 sc<-c(.01,.02,.015,.007);r<-c(.01,-.004,.008,-.005);w<-rep(1/p,p)
 den<-particle_density(r,L,h,sc,w)
 for(j in 1:N){B<-L[,,j]%*%diag(exp(h[p+1:2,j]/2));S<-(tcrossprod(B)+diag(exp(h[1:p,j])))*outer(sc,sc)
  expected<--.5*(p*log(2*pi)+as.numeric(determinant(S,logarithm=TRUE)$modulus)+sum(r*solve(S,r)))
  stopifnot(abs(expected-den$log_density[j])<1e-10,abs(sum(w*(S%*%w))-den$portfolio_variance[j])<1e-12)
 }
 stopifnot(abs(particle_logsumexp(c(-1000,-1000))-(-1000+log(2)))<1e-12)
 ix<-particle_resample(rep(1/N,N));stopifnot(identical(ix,seq_len(N)))
 pa<-array(0,c(3,p+2,N));pa[2,,]<-.8;pa[3,,]<-.2
 term<-list(L=L,para=pa,h=h,scale=sc);rr<-matrix(rnorm(6*p,sd=.01),6,p);dates<-as.Date('2023-01-01')+1:6
 a<-particle_pilot(term,rr,dates,18);b<-particle_pilot(term,rr[1:3,,drop=FALSE],dates[1:3],18)
 stopifnot(identical(a[1:3,],b),all(a$filtering_ESS>=1-1e-8&a$filtering_ESS<=N+1e-8),all(a$portfolio_PIT>=0&a$portfolio_PIT<=1))
 rr[4:6,]<-rr[4:6,]+.1;cc<-particle_pilot(term,rr,dates,18);stopifnot(identical(a[1:3,],cc[1:3,]))
})
writeLines('PASS: dense-vs-Woodbury density, portfolio variance, stable normalizer, resampling, prefix and future invariance, ESS/PIT bounds.','outputs/PARTICLE_CHECKS.txt')
cat('Particle checks passed\n')
