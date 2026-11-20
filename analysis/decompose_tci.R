# Descriptive two-input decomposition of the raw minus SPY-adjusted TCI.
# Keep the two fitted VARs in common log-variance coordinates before crossing
# their transition matrices and innovation covariances.
local({
  for (file in sort(list.files('R', pattern='[.]R$', full.names=TRUE))) source(file)
  source('experiments/market_residuals.R')
  protocol <- revision_protocol(getwd())
  panel <- revision_read('input/archived_returns.csv', protocol$development_end)
  adjusted <- revision_market_residuals(panel$R, panel$dates)
  assets <- colnames(adjusted$R)
  raw <- panel$R[adjusted$rows, assets, drop=FALSE]
  spec <- protocol$models[protocol$models$id == 'EWMA098', ]
  raw_bundle <- revision_bundle(raw, adjusted$dates, nrow(raw), spec, protocol)
  adjusted_covariance <- revision_covariance(adjusted$R, spec, protocol$warmup)
  adjusted_features <- esg_features(adjusted_covariance$H)$X
  valid <- adjusted_covariance$valid
  adjusted_packet <- esg_lag_packets(
    adjusted_features[valid, , drop=FALSE], raw_bundle$hmm$state, 1L,
    adjusted$dates[valid])
  raw_packet <- raw_bundle$packet
  stopifnot(identical(raw_packet$dates, adjusted_packet$dates),
    identical(raw_packet$state, adjusted_packet$state))
  adjusted_fits <- esg_var_from_packets(adjusted_packet, protocol$states, 40L, .1)

  # x_native = center + D x_standardized. Intercepts do not affect FEVD.
  to_native <- function(fit, scale) {
    list(
      A=lapply(fit$A, function(a) a * outer(scale, 1/scale)),
      Sigma=fit$Sigma * outer(scale, scale))
  }
  connectedness <- function(lag_fit, innovation_fit, horizon) {
    hybrid <- list(A=lag_fit$A, Sigma=innovation_fit$Sigma)
    esg_conn(esg_fevd(hybrid, horizon), assets)$TCI
  }
  cross_inputs <- function(raw_fit, adjusted_fit, horizon) {
    rr <- connectedness(raw_fit, raw_fit, horizon)
    ra <- connectedness(raw_fit, adjusted_fit, horizon)
    ar <- connectedness(adjusted_fit, raw_fit, horizon)
    aa <- connectedness(adjusted_fit, adjusted_fit, horizon)
    from_innovations <- ((rr-ra)+(ar-aa))/2
    from_dynamics <- ((rr-ar)+(ra-aa))/2
    stopifnot(abs(rr-aa-from_innovations-from_dynamics) < 1e-10)
    c(raw=rr, raw_lag_adjusted_innovations=ra,
      adjusted_lag_raw_innovations=ar, adjusted=aa, difference=rr-aa,
      innovation_covariance_component=from_innovations,
      lag_matrix_component=from_dynamics)
  }
  reference <- read.csv('outputs/matched_market_ablation/MATCHED_TCI_DIFFERENCES.csv')
  rows <- list()
  for (horizon in c(1L,10L)) for (state in seq_len(protocol$states)) {
    label <- raw_bundle$hmm$labels[state]
    raw_fit <- to_native(raw_bundle$var[[state]], raw_packet$scale)
    adjusted_fit <- to_native(adjusted_fits[[state]], adjusted_packet$scale)
    result <- cross_inputs(raw_fit, adjusted_fit, horizon)
    check <- reference[reference$contrast == paste0('H',horizon,'_',label),]
    stopifnot(nrow(check) == 1L,
      max(abs(result[c('raw','adjusted','difference')] -
        c(check$raw_TCI,check$market_adjusted_TCI,check$estimate))) < 1e-8)
    # A joint change of units leaves every crossed TCI invariant.
    units <- exp(seq(-.7,.6,length.out=length(assets)))
    rescaled <- cross_inputs(to_native(raw_fit, units), to_native(adjusted_fit, units), horizon)
    stopifnot(max(abs(result-rescaled)) < 1e-9)
    if (horizon == 1L) stopifnot(abs(result['lag_matrix_component']) < 1e-10)
    rows[[length(rows)+1L]] <- data.frame(state=label,horizon=horizon,t(result),
      scope='descriptive two-input decomposition of fitted TCI',check.names=FALSE)
  }
  result <- do.call(rbind, rows)
  out <- 'outputs/evidence'
  esg_csv(result, file.path(out,'MATCHED_TCI_COMPONENTS.csv'))
  writeLines(c(
    'Descriptive, added after the original matched exposure result.',
    'Both fitted VARs were transformed back to common log-variance units before crossing their inputs.',
    'The covariance contribution averages both orders of replacing the innovation covariance.',
    'The lag-matrix contribution averages both orders of replacing the VAR transition matrix.',
    'The two signed components add to the observed raw minus adjusted TCI difference.',
    'Hybrid fits are arithmetic substitutions. They neither identify causal shocks nor estimate investable strategies.',
    'The lag-matrix component can reflect omitted common factors as well as lag dependence.',
    'No hypothesis tests or new multiplicity family were added.',
    'Original TCI, additivity, common-coordinate rescaling and horizon-one checks passed.'),
    file.path(out,'MATCHED_TCI_COMPONENTS_SCOPE.txt'))
  print(result[,c('state','horizon','difference','innovation_covariance_component',
    'lag_matrix_component')],row.names=FALSE)
})
