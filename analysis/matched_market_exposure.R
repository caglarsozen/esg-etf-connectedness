# Fixed exploratory contrast: identical assets, calendar, and raw-panel state labels.
for(f in sort(list.files('R',pattern='[.]R$',full.names=TRUE)))source(f)
source('experiments/market_residuals.R')
source('experiments/matched_bootstrap.R')
pr<-revision_protocol(getwd());pan<-revision_read('input/archived_returns.csv',pr$development_end)
adj<-revision_market_residuals(pan$R,pan$dates);assets<-colnames(adj$R)
raw<-pan$R[adj$rows,assets,drop=FALSE]
spec<-pr$models[pr$models$id=='EWMA098',]
br<-revision_bundle(raw,adj$dates,nrow(raw),spec,pr)
ca<-revision_covariance(adj$R,spec,pr$warmup);fa<-esg_features(ca$H)
Xa<-fa$X[ca$valid,,drop=FALSE];colnames(Xa)<-assets
pa<-esg_lag_packets(Xa,br$hmm$state,1L,adj$dates[ca$valid]);praw<-br$packet
stopifnot(identical(pa$dates,praw$dates),identical(pa$state,praw$state),
          identical(colnames(raw),colnames(adj$R)),all(adj$latest_estimation_date<adj$dates))
block<-as.integer(Sys.getenv('ESG_MATCHED_BLOCK','10'));stopifnot(block%in%c(5,10,20))
out<-paste0('outputs/matched_market_ablation',if(block==10)'' else paste0('_b',block))
dir.create(out,recursive=TRUE,showWarnings=FALSE)
B<-9999L;labels<-br$hmm$labels
all_rows<-seq_len(nrow(pa$Y))
matched_statistic<-function(packet,ix=all_rows,fast=TRUE)
 matched_tci_statistic(packet,labels,assets,ix,fast)
indices<-matched_block_indices(length(all_rows),block,B,pr$seed)
checked<-c(list(all_rows),indices[unique(round(seq(1,B,length.out=12)))])
max_error<-max(vapply(checked,function(ix)max(
 abs(matched_statistic(praw,ix)-matched_statistic(praw,ix,FALSE)),
 abs(matched_statistic(pa,ix)-matched_statistic(pa,ix,FALSE))),numeric(1)))
stopifnot(max_error<1e-9)
r0<-matched_statistic(praw,fast=FALSE);a0<-matched_statistic(pa,fast=FALSE);observed<-r0-a0
boot<-matrix(NA_real_,B,length(observed),dimnames=list(NULL,names(observed)))
for(b in seq_len(B)){
 ix<-indices[[b]]
 boot[b,]<-matched_statistic(praw,ix)-matched_statistic(pa,ix)
}
tab<-do.call(rbind,lapply(seq_along(observed),function(j)data.frame(contrast=names(observed)[j],
 raw_TCI=unname(r0[j]),market_adjusted_TCI=unname(a0[j]),
 t(esg_basic_inference(observed[j],boot[,j],'greater')),block=block,check.names=FALSE)))
tab$p_holm_six<-p.adjust(tab$p_value,'holm')
inc<-do.call(rbind,lapply(labels,function(s){
 j1<-paste0('H1_',s);j10<-paste0('H10_',s)
 data.frame(state=s,t(esg_basic_inference(observed[j10]-observed[j1],boot[,j10]-boot[,j1],'two.sided')),block=block,check.names=FALSE)
}))
esg_csv(tab,file.path(out,'MATCHED_TCI_DIFFERENCES.csv'));esg_csv(inc,file.path(out,'EXPLORATORY_HORIZON_INCREMENT.csv'))
esg_save(list(observed=observed,boot=boot,raw=r0,adjusted=a0,
 rng_method='single_per_block_stream',seed=pr$seed+block),file.path(out,'BOOTSTRAP_DRAWS.rds'))
esg_csv(data.frame(block=block,draws=B,distinct_index_vectors=length(indices),
 checked_packet_sets=length(checked),maximum_TCI_error=max_error,
 seed=pr$seed+block,RNG=paste(RNGkind(),collapse='/')),
 file.path(out,'RNG_AND_EQUIVALENCE.csv'))
esg_csv(data.frame(date=adj$dates[br$cv$valid],state=br$hmm$state),file.path(out,'FIXED_RAW_REGIMES.csv'))
dput(list(block=block,B=B,seed=pr$seed+block,RNG=RNGkind(),key='single_per_block_stream',lambda=.98,window=252L,
 warmup=126L,p=1L,ridge=.1,through=pr$development_end,raw_states_fixed=TRUE),file=file.path(out,'PROTOCOL.R'))
writeLines(c('Exploratory, conditional on estimated variance paths and common raw-panel regime labels.',
 paste('Effective lag-packet calendar:',min(pa$dates),max(pa$dates),'n=',nrow(pa$Y)),
 'Both panels exclude SPY and include the same ten assets. Contemporaneous SPY is removed using preceding-252-day OLS.',
 'Residuals are not investable returns. Bootstrap resamples paired genuine lag packets; first-stage filters and HMM are held fixed.',
 'This does not identify a causal common-factor share. H10-H1 is an exploratory nonadditive diagnostic.',
 'The superseded per-replicate seeding caused duplicate indices. This run uses one fixed RNG stream per block length.',
 'All resampled index vectors are distinct. Fast VAR inputs match the core VAR on the observed and 12 resampled packet sets.'),file.path(out,'SCOPE.txt'))
print(tab)
