# Inspect completed tables using base R only. Does not estimate or select models.
local({
  files <- c(
    market = 'outputs/evidence/MATCHED_MARKET_BLOCK_SENSITIVITY.csv',
    components = 'outputs/evidence/MATCHED_TCI_COMPONENTS.csv',
    allocation = 'outputs/holding_horizon/performance_2023_cost10.csv',
    ledger = 'outputs/evidence/EXPERIMENT_LEDGER.csv')
  if (!all(file.exists(files))) stop(
    'Run the analysis from the repository root to generate these tables. ',
    'Missing: ', paste(files[!file.exists(files)], collapse = ', '))
  z <- lapply(files, read.csv, check.names = FALSE)
  cat('\nMatched exposure adjustment: exploratory, conditional, two-sided inference\n')
  print(z$market[z$market$block == 10,
    c('contrast','raw_TCI','market_adjusted_TCI','estimate','lower','upper',
      'p_holm_all30_two_sided')], row.names = FALSE)
  cat('\nDescriptive components of the H10 TCI difference\n')
  print(z$components[z$components$horizon == 10L,
    c('state','difference','innovation_covariance_component','lag_matrix_component')],
    row.names=FALSE)
  gate_file <- 'outputs/evidence/LONG_ANCHORED_GATES.csv'
  if (file.exists(gate_file)) {
    cat('\nLong anchored dynamic FSV diagnostic gates\n')
    gate <- read.csv(gate_file,check.names=FALSE)
    print(gate[,c('cutoff','quantities','failed','max_Rhat','min_bulk_ESS',
      'min_tail_ESS','passed')], row.names=FALSE)
  } else cat('\nLong anchored diagnostic tables are pending in this snapshot.\n')
  cat('\nAllocation: daily net ES95, 2023 development, 10 bp trading cost\n')
  print(z$allocation[,c('strategy','ES_loss_95','observations')], row.names = FALSE)
  cat('\nEvidence scope\n')
  print(z$ledger, row.names = FALSE)
  cat('\nSee README.md for the study scope and replication instructions.\n')
})
