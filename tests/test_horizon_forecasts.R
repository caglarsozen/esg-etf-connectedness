for(f in sort(list.files('R',pattern='[.]R$',full.names=TRUE)))source(f)
source('experiments/holding_horizon.R')
checks<-list();check<-function(n,ok){stopifnot(isTRUE(ok));checks[[length(checks)+1]]<<-data.frame(check=n,passed=TRUE)}
S<-matrix(c(.0001,.00006,.00006,.0004),2);pk<-list(center=log(diag(S)),scale=c(.7,1.2))
A<-diag(c(.7,.4));cc<-c(.05,-.03);Q<-diag(c(.1,.2));fit<-list(A=list(A),B=rbind(cc,t(A)),Sigma=Q)
z<-horizon_forecast(S,pk,list(fit),1,matrix(1,1,1),3)
check('horizon_one',max(abs(horizon_forecast(S,pk,list(fit),1,matrix(1,1,1),1)$H-S))<1e-12)
check('Gaussian_second_step',max(abs(z$variance_path[2,]-exp(pk$center+pk$scale*cc+.5*pk$scale^2*diag(Q))))<1e-12)
m3<-cc+as.numeric(A%*%cc);V3<-A%*%Q%*%t(A)+Q
check('Gaussian_third_step',max(abs(z$variance_path[3,]-exp(pk$center+pk$scale*m3+.5*pk$scale^2*diag(V3))))<1e-12)
check('correlation_preserved',max(abs(cov2cor(z$H)-cov2cor(S)))<1e-10)
check('positive_definite',min(eigen(z$H,symmetric=TRUE)$values)>0)
T<-matrix(c(.8,.2,.3,.7),2,byrow=TRUE);prob<-c(.4,.6)
dup<-horizon_forecast(S,pk,list(fit,fit),prob,T,3)
check('duplicate_states_collapse',max(abs(dup$H-z$H))<1e-12)
fit2<-list(A=list(A*.5),B=rbind(-cc,t(A*.5)),Sigma=Q*2)
pk$center<-pk$center-c(.2,-.1)*pk$scale;x<-c(.2,-.1)
zz<-horizon_forecast(S,pk,list(fit,fit2),prob,T,2)
check('transition_probabilities',max(abs(zz$final_state_probability-as.numeric(prob%*%T)))<1e-12)
means<-list(cc+as.numeric(A%*%x),-cc+as.numeric((A*.5)%*%x));seconds<-list(tcrossprod(means[[1]])+Q,tcrossprod(means[[2]])+2*Q)
for(j in 1:2){
 mm<-prob[1]*T[1,j]*means[[1]]+prob[2]*T[2,j]*means[[2]]
 ss<-prob[1]*T[1,j]*seconds[[1]]+prob[2]*T[2,j]*seconds[[2]]
 check(paste0('lag_state_joint_mean_',j),max(abs(zz$final_joint_mean[[j]]-mm))<1e-12)
 check(paste0('lag_state_joint_second_',j),max(abs(zz$final_joint_second[[j]]-ss))<1e-12)
}
check('probability_sum',abs(sum(zz$final_state_probability)-1)<1e-12)
xx<-matrix(seq_len(80),40,2);st<-rep(1:2,20)
check('fit_uses_lagged_state',identical(esg_lag_packets(xx,st)$state,st[1:39]))
esg_csv(do.call(rbind,checks),'outputs/HORIZON_FORECAST_CHECKS.csv');cat(length(checks),'horizon checks passed\n')
