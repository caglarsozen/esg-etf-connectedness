# Run the matched exposure, holding-horizon, and anchored FSV analyses.
esg_run_revision<-function(){
 for(p in c('factorstochvol','posterior','lpSolve','quadprog'))if(!requireNamespace(p,quietly=TRUE))stop('Run install_dependencies.R first: ',p)
 envnames<-c('ESG_MATCHED_BLOCK','ESG_FSV_CUTOFF','ESG_FSV_FOLDER');old<-Sys.getenv(envnames,unset=NA_character_)
 on.exit(for(i in seq_along(envnames))if(is.na(old[i]))Sys.unsetenv(envnames[i]) else do.call(Sys.setenv,setNames(list(old[i]),envnames[i])),add=TRUE)
 run<-function(f){message('Running ',f);source(f,local=new.env(parent=.GlobalEnv),encoding='UTF-8')}
 if(!file.exists('outputs/DEVELOPMENT_targets.rds'))run('analysis/develop_covariance.R')
 run('tests/test_horizon_forecasts.R');run('tests/test_particle_prediction.R')
 run('analysis/holding_horizon_allocation.R');run('analysis/common_factor_example.R')
 matched_work<-function(block){
  Sys.setenv(ESG_MATCHED_BLOCK=as.character(block))
  run('analysis/matched_market_exposure.R')
  invisible(TRUE)
 }
 workers<-max(1L,min(3L,as.integer(Sys.getenv('ESG_WORKERS','2'))))
 if(.Platform$OS.type=='unix'&&workers>1L){
  matched_runs<-parallel::mclapply(c(5,10,20),matched_work,mc.cores=workers,
    mc.preschedule=FALSE,mc.set.seed=FALSE)
  stopifnot(all(vapply(matched_runs,isTRUE,logical(1))))
 }else invisible(lapply(c(5,10,20),matched_work))
 run('analysis/decompose_tci.R')
 for(cut in c('2023-02-28','2023-09-29')){
  Sys.setenv(ESG_FSV_CUTOFF=cut);run('analysis/anchored_fsv.R')
  folder<-file.path('outputs/fsv',paste('anchored_dynamic',gsub('-','',cut),Sys.getenv('ESG_FSV_DRAWS','160000'),Sys.getenv('ESG_FSV_BURNIN','40000'),Sys.getenv('ESG_FSV_THIN','40'),sep='_'))
  Sys.setenv(ESG_FSV_FOLDER=folder)
  if(isTRUE(read.csv(file.path(folder,'MCMC_GATE.csv'))$passed)){
   run('analysis/posterior_network.R');run('analysis/particle_forecast_pilot.R')
  }else message('MCMC rejected; latent downstream analysis skipped at ',cut)
 }
 run('analysis/summarize_fsv.R')
 vals<-Sys.getenv(c('ESG_FSV_DRAWS','ESG_FSV_BURNIN','ESG_FSV_THIN'));defaults<-c('160000','40000','40')
 if(all(ifelse(vals=='',defaults,vals)==defaults)){
  run('analysis/summarize_evidence.R');run('analysis/plot_evidence.R')
  run('analysis/export_tables.R')
 }else message('Custom MCMC: inspect custom output folders; canonical report not rebuilt.')
 writeLines('Research computations completed. See README.md for the study scope.','outputs/REVISION_RUN_STATUS.txt')
}
esg_run_revision()
