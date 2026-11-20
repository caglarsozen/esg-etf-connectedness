for(f in sort(list.files('R',pattern='[.]R$',full.names=TRUE)))source(f)
source('experiments/particle_prediction.R')
folder<-Sys.getenv('ESG_FSV_FOLDER');stopifnot(nzchar(folder))
gate<-read.csv(file.path(folder,'MCMC_GATE.csv'));stopifnot(gate$passed)
cutoff<-as.Date(gate$cutoff);pan<-revision_read('input/archived_returns.csv',as.Date('2023-10-31'))
tt<-head(which(pan$dates>cutoff),21L);stopifnot(length(tt)==21)
out<-file.path(folder,'predictive_pilot');dir.create(out,showWarnings=FALSE)
zz<-lapply(1:4,function(ch)readRDS(file.path(folder,paste0('terminal_particles_chain_',ch,'.rds'))))
take<-lapply(zz,function(z)unique(round(seq(1,ncol(z$h),length.out=2000L))))
p<-ncol(pan$R);N<-sum(lengths(take));L<-array(NA_real_,c(p,2,N));pa<-array(NA_real_,c(3,p+2,N));hh<-matrix(NA_real_,p+2,N)
pos<-0L
for(ch in 1:4){ix<-take[[ch]];target<-pos+seq_along(ix);L[,,target]<-zz[[ch]]$facload[,,ix];pa[,,target]<-zz[[ch]]$para[,,ix];hh[,target]<-zz[[ch]]$h[,ix];pos<-pos+length(ix)}
term<-list(L=L,para=pa,h=hh,scale=zz[[1]]$scale_next);rm(zz)
rows<-list()
# Preserve the archived Monte Carlo seeds for each cutoff and replicate.
pilot_seed_base <- 20258676L
for(replica in 1:2){
 z<-particle_pilot(term,pan$R[tt,,drop=FALSE],pan$dates[tt],esg_seed(c('predictive_path',cutoff,replica),pilot_seed_base))
 z$model<-'FSV_particle';z$replicate<-replica;rows[[length(rows)+1]]<-z
}
pr<-revision_protocol(getwd())
for(nm in c('EWMA094','EWMA098')){
 cv<-revision_covariance(pan$R,pr$models[pr$models$id==nm,],pr$warmup)
 z<-do.call(rbind,lapply(tt,function(t){S<-cv$H[t,,];r<-pan$R[t,]
  data.frame(date=pan$dates[t],negative_log_density=.5*(p*log(2*pi)+as.numeric(determinant(S,logarithm=TRUE)$modulus)+sum(r*solve(S,r))),
  portfolio_PIT=pnorm(mean(r)/sqrt(sum(S)/p^2)),filtering_ESS=NA_real_,unique_ancestors=NA_integer_,ancestors_before=NA_integer_,resampled=FALSE)
 }));z$model<-nm;z$replicate<-0L;rows[[length(rows)+1]]<-z
}
all<-do.call(rbind,rows);esg_csv(all,file.path(out,'DAILY_PREDICTIONS.csv'))
tab<-do.call(rbind,lapply(rows,function(z)data.frame(model=z$model[1],replicate=z$replicate[1],days=nrow(z),
 mean_negative_log_density=mean(z$negative_log_density),min_filtering_ESS=if(all(is.na(z$filtering_ESS)))NA else min(z$filtering_ESS),
 last_unique_ancestors=tail(z$unique_ancestors,1),lower_5pct_PIT_count=sum(z$portfolio_PIT<.05))))
esg_csv(tab,file.path(out,'SUMMARY.csv'))
writeLines(c('8000 particles, two fixed Monte Carlo replicates, 21 next trading days, no post-2023-10-31 observations.',
 'State propagated before outcome; mixture density and PIT use prior weights; weights updated only after scoring.',
 'Static parameters resampled without rejuvenation: severe ancestral depletion invalidates general superiority claims.',
 'Decimal-return density with full Gaussian constant; marginal equal-weight portfolio PIT is diagnostic, not an allocation result.'),file.path(out,'SCOPE.txt'))
print(tab)
