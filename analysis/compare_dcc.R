for(f in sort(list.files('R',pattern='[.]R$',full.names=TRUE)))source(f)
pr<-readRDS('outputs/PRESPECIFIED_PROTOCOL.rds')
a<-read.csv('outputs/historical_evaluation/returns_cost10.csv',check.names=FALSE)
b<-read.csv('outputs/dcc_benchmark/returns_cost10.csv',check.names=FALSE)
esg_assert(identical(a$date,b$date),'DCC/main calendars differ.')
R<-cbind(as.matrix(a[,-1]),as.matrix(b[,-1,drop=FALSE]))
pairs<-list(c('C','DCC_Gaussian_MinVar'),c('MinVar','DCC_Gaussian_MinVar'),
  c('C_WMIX','DCC_Gaussian_MinVar'),c('C_WMIX_N','DCC_Gaussian_MinVar'))
z<-do.call(rbind,lapply(pr$block_lengths,function(bl)revision_pairs(R,pairs,pr$bootstrap_B,bl)))
z$family<-'four DCC contrasts per metric and block length; additional benchmark family'
esg_csv(z,'outputs/dcc_benchmark/paired_inference.csv')
