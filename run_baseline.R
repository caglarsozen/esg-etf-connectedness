# Run from the repository root.
# Default: the completed observation-driven research pipeline. No anchored run.
# Optional scaled Factor-SV preflights: Sys.setenv(ESG_RUN_SCALED_FSV='true').
esg_run_all<-function(){
  if (!file.exists('input/archived_returns.csv'))
    stop('Provide input/archived_returns.csv as described in README.md.', call. = FALSE)
  for(p in c('lpSolve','quadprog'))
    if(!requireNamespace(p,quietly=TRUE))stop('Run source("install_dependencies.R") first: ',p)
  for(script in c('tests/test_revision.R','analysis/develop_covariance.R','analysis/development_network.R',
                 'analysis/freeze_historical_design.R','analysis/historical_allocation.R',
                 'tests/test_information.R','analysis/dcc_benchmark.R','analysis/compare_dcc.R',
                 'analysis/network_reference_models.R','analysis/network_robustness.R',
                 'analysis/plot_baseline.R')){
    message('Running ',script);source(script,local=new.env(parent=.GlobalEnv),encoding='UTF-8')
  }
  for(f in sort(list.files('R',pattern='[.]R$',full.names=TRUE)))source(f,local=.GlobalEnv)
  if(tolower(Sys.getenv('ESG_RUN_SCALED_FSV','false'))=='true')
    for(cutoff in as.Date(c('2023-02-28','2023-09-29'))){
    cutoff<-as.Date(cutoff,origin='1970-01-01')
    pan<-revision_read('input/archived_returns.csv',cutoff)
    revision_fsv_fit(pan$R,pan$dates,revision_fsv_spec('scaled_constant'),
      file.path('outputs','fsv',paste('scaled_constant',format(cutoff,'%Y%m%d'),12000,6000,4,sep='_')))
  }
  source('analysis/summarize_fsv.R',local=new.env(parent=.GlobalEnv),encoding='UTF-8')
  writeLines('Computations completed. See README.md for the study scope.',
             'outputs/FULL_RUN_STATUS.txt')
}
esg_run_all()
