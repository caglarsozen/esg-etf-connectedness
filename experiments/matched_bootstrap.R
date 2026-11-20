# Draw from one stream per block length. Per-replicate calls to esg_seed()
# collided for the old matched bootstrap and are deliberately avoided here.
matched_block_indices <- function(n, block, draws, base_seed) {
  indices <- esg_with_seed(base_seed + block, replicate(
    draws, esg_mbb_indices(n, block), simplify=FALSE))
  stopifnot(length(indices)==draws,anyDuplicated(indices)==0L)
  indices
}

# The matched TCI only uses the fitted slopes and innovation covariance.
# This is esg_var_fit_rows' ridge and stability calculation without its
# unused residual-test and BIC work. Equivalence to the core fit is checked
# on the original packets and fixed resampled packets in the caller.
matched_var_inputs <- function(y, z, min_obs=40L, ridge=.1, radius=.995) {
  observations <- nrow(y)
  assets <- ncol(y)
  stopifnot(observations>=max(min_obs,assets+6L),ncol(z)==assets+1L)
  penalty <- diag(c(0,rep(observations*ridge,assets)),assets+1L)
  coefficients <- solve(crossprod(z)+penalty,crossprod(z,y))
  lag_matrix <- t(coefficients[-1,,drop=FALSE])
  max_root <- function(a) max(Mod(eigen(a,only.values=TRUE)$values))
  while(max_root(lag_matrix)>=radius) {
    coefficients[-1,] <- coefficients[-1,,drop=FALSE]*.97
    coefficients[1,] <- colMeans(y) -
      as.numeric(colMeans(z[,-1,drop=FALSE]) %*% coefficients[-1,,drop=FALSE])
    lag_matrix <- t(coefficients[-1,,drop=FALSE])
  }
  residuals <- y-z%*%coefficients
  residuals <- sweep(residuals,2,colMeans(residuals),'-')
  list(A=list(lag_matrix),Sigma=esg_spd(crossprod(residuals)/observations))
}

matched_tci_statistic <- function(packet, labels, assets, indices, fast=TRUE) {
  if (fast) {
    selected_states <- packet$state[indices]
    fits <- lapply(seq_along(labels),function(state) {
      rows <- indices[selected_states==state]
      matched_var_inputs(packet$Y[rows,,drop=FALSE],packet$Z[rows,,drop=FALSE])
    })
  } else {
    fits <- esg_var_from_packets(packet,length(labels),40L,.1,row_index=indices)
  }
  result <- numeric()
  for (horizon in c(1L,10L)) for (state in seq_along(labels)) {
    result[paste0('H',horizon,'_',labels[state])] <-
      esg_conn(esg_fevd(fits[[state]],horizon),assets)$TCI
  }
  result
}
