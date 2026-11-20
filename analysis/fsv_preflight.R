# Optional expensive research module. Parameters are explicit; no fallback.
for(f in sort(list.files('R',pattern='[.]R$',full.names=TRUE)))source(f)
args<-commandArgs(trailingOnly=TRUE)
candidate<-if(length(args)>=1)args[1] else 'scaled_constant'
draws<-if(length(args)>=2)as.integer(args[2]) else 12000L
burnin<-if(length(args)>=3)as.integer(args[3]) else 6000L
thin<-if(length(args)>=4)as.integer(args[4]) else 4L
cutoffs<-if(length(args)>=5)as.Date(args[5]) else as.Date(c('2023-02-28','2023-09-29'))
spec<-revision_fsv_spec(candidate,draws,burnin,thin)
for(cutoff in as.character(cutoffs)){
  esg_assert(as.Date(cutoff)<=as.Date('2023-10-31'),'Preflight may not inspect historical test.')
  pan<-revision_read('input/archived_returns.csv',as.Date(cutoff))
  folder<-file.path('outputs','fsv',paste(candidate,gsub('-','',cutoff),draws,burnin,thin,sep='_'))
  fit<-revision_fsv_fit(pan$R,pan$dates,spec,folder)
  # Every candidate is retained, including failure. No unqualified network fit.
}
