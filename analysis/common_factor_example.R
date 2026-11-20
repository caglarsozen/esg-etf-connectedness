for(f in sort(list.files('R',pattern='[.]R$',full.names=TRUE)))source(f)
# X_t = a*f_t + e_t; only the common factor f_t = phi*f_(t-1)+eta_t has dynamics.
# Observed process is generally VARMA. Below is its population VAR(1) projection.
out<-'outputs/common_factor_example';dir.create(out,showWarnings=FALSE)
scenarios<-list(Calm=list(phi=.8,vf=1,noise=c(4,1,1)),Stress=list(phi=.8,vf=2,noise=c(.1,1,1)))
a<-c(1,1,1);rows<-list()
for(nm in names(scenarios)){
 s<-scenarios[[nm]];G0<-s$vf*tcrossprod(a)+diag(s$noise);G1<-s$phi*s$vf*tcrossprod(a)
 A<-G1%*%solve(G0);Sigma<-G0-A%*%G0%*%t(A)
 cc<-esg_conn(esg_fevd(list(A=list(A),Sigma=Sigma),10),c('A','B','C'))
 rows[[nm]]<-data.frame(regime=nm,TCI=cc$TCI,node=cc$asset$asset,NET=cc$asset$NET)
 esg_csv(A,file.path(out,paste0(nm,'_projection_A.csv')))
}
esg_csv(do.call(rbind,rows),file.path(out,'EXACT_EXAMPLE.csv'));dput(scenarios,file=file.path(out,'PARAMETERS.R'))
writeLines(c('Known common-factor mechanism, not a new theorem or ETF-calibrated model.',
 'No structural cross-asset feedback. GFEVD belongs to the population VAR(1) projection; residual whiteness is not asserted.',
 'Two parameterizations illustrate possibility, not a data-generating explanation of observed regimes.'),file.path(out,'SCOPE.txt'))
print(do.call(rbind,rows))
