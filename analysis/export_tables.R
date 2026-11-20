# Export manuscript table fragments from completed outputs, without new estimation.
local({
  out <- 'outputs/manuscript_tables'
  dir.create(out, recursive=TRUE, showWarnings=FALSE)
  read <- function(x) read.csv(x, check.names=FALSE)
  num <- function(x,d=4L) formatC(x,format='f',digits=d,decimal.mark='.')
  table_tex <- function(d,caption,label,note,file,align) {
    lines <- c('\\begin{table}[!htbp]','\\centering',
      paste0('\\caption{',caption,'}'),paste0('\\label{',label,'}'),
      '\\small',paste0('\\begin{tabular}{',align,'}'),'\\toprule',
      paste0(paste(names(d),collapse=' & '),' \\\\'),'\\midrule',
      apply(d,1,function(x)paste0(paste(x,collapse=' & '),' \\\\')),
      '\\bottomrule','\\end{tabular}','\\par\\smallskip',
      '\\begin{minipage}{0.97\\textwidth}',
      paste0('\\footnotesize\\textit{Note:} ',note),
      '\\end{minipage}','\\end{table}')
    writeLines(lines,file.path(out,file),useBytes=TRUE)
  }
  a <- read('outputs/evidence/MATCHED_MARKET_BLOCK_SENSITIVITY.csv')
  a <- a[a$block==10,]
  stopifnot(nrow(a)==6,all(is.finite(a$p_holm_all30_two_sided)))
  stress <- a[a$contrast=='H10_Stress',]
  amplification <- read('outputs/evidence/POSTHOC_STRESS_AMPLIFICATION.csv')
  amplification <- amplification[amplification$block==10,]
  stopifnot(nrow(stress)==1L,nrow(amplification)==1L)
  values <- c(MatchedStressRaw=num(stress$raw_TCI),
    MatchedStressAdjusted=num(stress$market_adjusted_TCI),
    MatchedStressDifference=num(stress$estimate),
    MatchedStressLow=num(stress$lower),MatchedStressHigh=num(stress$upper),
    MatchedStressHolm=num(stress$p_holm_all30_two_sided),
    MatchedAmplification=num(amplification$estimate),
    MatchedAmplificationLow=num(amplification$lower),
    MatchedAmplificationHigh=num(amplification$upper),
    MatchedAmplificationHolm=num(amplification$p_holm_all30_two_sided),
    MatchedDraws=format(stress$draws,big.mark=',',trim=TRUE,scientific=FALSE))
  writeLines(paste0('\\providecommand{\\',names(values),'}{',values,'}'),
    file.path(out,'MatchedValues.tex'))
  write.csv(a,file.path(out,'Table1_matched_exposure_numeric.csv'),row.names=FALSE)
  d <- data.frame(State=sub('^H[0-9]+_','',a$contrast),
    Horizon=sub('_.*','',a$contrast),Raw=num(a$raw_TCI),
    Adjusted=num(a$market_adjusted_TCI),Difference=num(a$estimate),
    Interval=paste0('[',num(a$lower),', ',num(a$upper),']'),
    Holm=num(a$p_holm_all30_two_sided),check.names=FALSE)
  names(d)[6:7] <- c('95\\% interval','Holm $p$')
  table_tex(d,'Matched common-market exposure comparison','tab:matched_exposure',
    paste('Ten identical non-SPY ETFs, common dates and raw-panel state labels.',
      'Connectedness is in percent. Differences are raw minus adjusted, in percentage points.',
      'The intervals are pointwise basic 95\\% moving-block intervals at block length 10.',
      'Inference conditions on the fitted first-stage paths and labels.',
      'Two-sided Holm adjustment covers all 30 reported level, horizon and stress-amplification contrasts across block lengths 5, 10 and 20.',
      'The exploratory network comparison uses the adjusted series for descriptive analysis.'),
    'Table1_matched_exposure.tex','llrrrlr')
  gate_file <- 'outputs/evidence/LONG_ANCHORED_GATES.csv'
  if (file.exists(gate_file)) {
  g <- read(gate_file)
  stopifnot(nrow(g)>=1)
  write.csv(g,file.path(out,'Table2_long_fsv_numeric.csv'),row.names=FALSE)
  d <- data.frame(Date=g$cutoff,Quantities=g$quantities,Failed=g$failed,
    Rhat=num(g$max_Rhat),Bulk=num(g$min_bulk_ESS,1),Tail=num(g$min_tail_ESS,1),
    check.names=FALSE)
  names(d)[4:6] <- c('Max. $\\widehat R$','Min. bulk ESS','Min. tail ESS')
  table_tex(d,'Long anchored dynamic factor-SV diagnostic gates','tab:long_fsv_gate',
    paste('Each cutoff uses four chains, 40,000 burn-in and 160,000 subsequent iterations, with thinning 40 (4,000 retained draws per chain).',
      'All idiosyncratic variances are dynamic.',
      'Thresholds are $\\widehat R\\leq1.01$ and bulk/tail ESS of at least 400.',
      'The monitored historical covariance quantities cover three dates.',
      'Passing the gate is not model adequacy or complete latent-path convergence; the two fitted prefixes overlap.',
      'Numerical solver warnings are retained.'),
    'Table2_long_fsv.tex','lrrrrr')
  }
  p <- read('outputs/holding_horizon/performance_2023_cost10.csv')
  p$ES_percent <- 100*p$ES_loss_95
  p$ES_difference_bp <- 1e4*(p$ES_loss_95-p$ES_loss_95[p$strategy=='C'])
  write.csv(p,file.path(out,'Table3_holding_allocation_numeric.csv'),row.names=FALSE)
  labels <- c(C='Covariance--CVaR',H_OWN='Own dynamics',H_VAR='Cross-variable dynamics',
    H_STATE_OWN='State-dependent own dynamics',H_STATE_VAR='State-dependent cross-variable dynamics')
  d <- data.frame(Strategy=unname(labels[p$strategy]),ES=num(p$ES_percent),
    Difference=num(p$ES_difference_bp),Days=p$observations,check.names=FALSE)
  names(d)[2:3] <- c('ES$_{0.95}$ (\\%)','Difference (bp)')
  table_tex(d,'Development expected shortfall of the additional holding-horizon rules','tab:holding_allocation',
    paste('The outcome is daily net portfolio return in the 2023 development period, with 10 basis points of proportional trading cost.',
      'The difference subtracts covariance--CVaR; lower is better.',
      'The additional arms use 21-session average conditional daily variances with a fixed initial correlation matrix,',
      'rather than exact compounded 21-day CVaR.',
      'All four additional rules have higher observed ES; this does not establish universal inferiority or equivalence.'),
    'Table3_holding_allocation.tex','lrrr')
  components <- read('outputs/evidence/MATCHED_TCI_COMPONENTS.csv')
  components <- components[components$horizon==10,]
  write.csv(components,file.path(out,'Table4_matched_tci_components_numeric.csv'),row.names=FALSE)
  d <- data.frame(State=components$state,
    Difference=num(components$difference),
    Covariance=num(components$innovation_covariance_component),
    Dynamics=num(components$lag_matrix_component),check.names=FALSE)
  names(d)[3:4] <- c('Innovation covariance','VAR lag matrix')
  table_tex(d,'Components of the matched ten-session TCI difference','tab:matched_tci_components',
    paste('Units are percentage points.',
      'The two fitted VARs are expressed in common log-variance coordinates.',
      'Each signed component averages both orders of substituting its input and the two components sum to the observed difference.',
      'The table describes fitted-model sensitivity. It provides point estimates without causal identification or component tests.'),
    'Table4_matched_tci_components.tex','lrrr')
  writeLines(c('Generated from the completed R output tables.',
    'LaTeX fragments require booktabs. Inspect formatting within the target manuscript.',
    'The accompanying numeric CSVs contain full-precision output.'),
    file.path(out,'README.txt'))
  cat('Manuscript table fragments and numeric source tables exported.\n')
})
