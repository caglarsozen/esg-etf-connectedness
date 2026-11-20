# Freeze the fixed research arms before historical evaluation. Not a preregistration.
for(f in sort(list.files('R',pattern='[.]R$',full.names=TRUE)))source(f)
pr<-readRDS('outputs/PRESPECIFIED_PROTOCOL.rds');pr$root<-getwd()
id<-readLines('outputs/FROZEN_COVARIANCE_MODEL.txt',warn=FALSE)[1]
pan<-revision_read('input/archived_returns.csv',pr$development_end)
files<-sort(c(list.files('R',pattern='[.]R$',full.names=TRUE),'analysis/historical_allocation.R'))
freeze<-list(protocol=pr,model=pr$models[match(id,pr$models$id),],
  development_identity=esg_hash(pan),source_md5=tools::md5sum(files),
  arms=c('C','C_R','C_N','C_WCOMP','C_WMIX','C_WMIX_N','MinVar','EqualWeight','C_EWMA094','MinVar_EWMA094'),
  created_utc=format(Sys.time(),tz='UTC',usetz=TRUE),
  status='retrospective protocol freeze; later period previously examined by an earlier design')
freeze_path<-'outputs/FROZEN_HISTORICAL_DESIGN.rds'
if(file.exists(freeze_path)){
  old<-readRDS(freeze_path)
  esg_assert(identical(old$source_md5,freeze$source_md5)&&
    identical(old$development_identity,freeze$development_identity)&&
    identical(old$model,freeze$model)&&identical(old$arms,freeze$arms),
    'A different historical design already exists. Preserve its outputs; use a separately labeled research run for changes.')
  freeze<-old
}else esg_save(freeze,freeze_path)
esg_csv(data.frame(file=names(freeze$source_md5),md5=unname(freeze$source_md5)),
  'outputs/FROZEN_SOURCE_IDENTITY.csv')
writeLines(c('Historical evaluation may proceed with these exact files.',
  'The earlier EWMA pilot already examined this period. It is not a virgin confirmatory sample.',
  'Any method change after inspection requires a new development label and new future validation.'),
  'outputs/FREEZE_SCOPE.txt')
cat('Frozen',id,'and',length(freeze$arms),'fixed arms.\n')
