# Usage: Rscript --vanilla reviewer_sensitivity.R PROJECT_ROOT OUTPUT_DIR VIXCLS_CSV
# The VIX file is the FRED VIXCLS download with date and index-value columns.
# This supplement leaves the evidence and checkpoints in PROJECT_ROOT unchanged.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) {
  stop("Provide the repository directory, a separate output directory, and a FRED VIXCLS CSV.")
}
root <- normalizePath(args[1], mustWork = TRUE)
output <- file.path(normalizePath(dirname(args[2]), mustWork = TRUE), basename(args[2]))
vix_file <- normalizePath(args[3], mustWork = TRUE)
if (output == root || startsWith(output, paste0(root, .Platform$file.sep))) {
  stop("Use an output directory outside the repository.")
}
dir.create(output, recursive = TRUE, showWarnings = FALSE)
for (file in sort(list.files(file.path(root, "R"), "[.]R$", full.names = TRUE))) {
  source(file)
}
source(file.path(root, "experiments", "market_residuals.R"))
save_table <- function(x, name) esg_csv(x, file.path(output, paste0(name, ".csv")))
read_project <- function(path) read.csv(file.path(root, path), check.names = FALSE)

pr <- revision_protocol(root)
sample <- revision_read(file.path(root, "input", "archived_returns.csv"), pr$development_end)
spec <- pr$models[pr$models$id == "EWMA098", ]
stopifnot(readLines(file.path(root, "outputs", "FROZEN_COVARIANCE_MODEL.txt"), n = 1L) == spec$id)
adjusted <- revision_market_residuals(sample$R, sample$dates)
raw <- sample$R[adjusted$rows, colnames(adjusted$R), drop = FALSE]
raw_fit <- revision_bundle(raw, adjusted$dates, nrow(raw), spec, pr)
valid <- raw_fit$cv$valid
days <- adjusted$dates[valid]
model <- raw_fit$hmm
labels <- model$labels
frozen <- read_project("outputs/matched_market_ablation/FIXED_RAW_REGIMES.csv")
stopifnot(identical(as.character(days), frozen$date), identical(model$state, frozen$state))

cat("Review supplement: HMM, VIX, and network diagnostics\n")
hmm_fits <- setNames(lapply(2:4, function(k) {
  tryCatch(
    if (k == 3L) model else esg_hmm_fit(
      raw_fit$features$Z[valid, , drop = FALSE], k = k, starts = 3L,
      seed = esg_seed(c("regime", tail(adjusted$dates, 1L))), tol = 1e-6
    ),
    error = function(e) e
  )
}), 2:4)
comparison <- lapply(2:4, function(k) {
  fit <- hmm_fits[[as.character(k)]]
  if (inherits(fit, "error")) {
    return(data.frame(states = k, loglik = NA_real_, BIC = NA_real_, converged = FALSE,
                      best_start = NA_integer_, iterations = NA_integer_,
                      min_occupancy = NA_real_, note = conditionMessage(fit)))
  }
  data.frame(states = k, loglik = fit$loglik, BIC = fit$BIC, converged = fit$converged,
             best_start = fit$start, iterations = fit$iterations,
             min_occupancy = min(fit$occupancy), note = "")
})
save_table(do.call(rbind, comparison), "HMM_STATE_COUNT")
probs <- function(x, prefix) {
  colnames(x) <- paste0(prefix, "_", labels)
  x
}
stopifnot(max(abs(rowSums(model$filtered) - 1)) < 1e-12,
          max(abs(rowSums(model$smoothed) - 1)) < 1e-12)
save_table(data.frame(date = days, filtered_label = labels[model$state],
                      smoothed_label = labels[max.col(model$smoothed, ties.method = "first")],
                      viterbi_label = labels[model$viterbi],
                      probs(model$filtered, "filtered"), probs(model$smoothed, "smoothed")),
           "MATCHED_HMM_PROBABILITIES")
save_table(data.frame(from_state = labels, model$P, mean_duration = model$duration,
                      filtered_days = tabulate(model$state, length(labels)),
                      smoothed_occupancy = model$occupancy), "MATCHED_HMM_TRANSITIONS")
save_table(data.frame(state = labels, initial = model$initial,
                      log_variance_mean = model$mu[, 1], log_variance_variance = model$variance[, 1],
                      abs_correlation_mean = model$mu[, 2], abs_correlation_variance = model$variance[, 2]),
           "MATCHED_HMM_EMISSIONS")
save_table(data.frame(feature = colnames(raw_fit$features$Z), center = model$center, scale = model$scale),
           "MATCHED_HMM_FEATURE_SCALING")

vix <- read.csv(vix_file, na.strings = c("", "."), stringsAsFactors = FALSE)
stopifnot(ncol(vix) == 2L, names(vix)[2] == "VIXCLS", !anyDuplicated(vix[[1]]))
vix_value <- as.numeric(vix$VIXCLS[match(days, as.Date(vix[[1]]))])
stopifnot(length(vix_value) == length(days), all(is.finite(vix_value)))
save_table(do.call(rbind, lapply(seq_along(labels), function(s) {
  x <- vix_value[model$state == s]
  data.frame(state = labels[s], dates = length(x), VIX_median = median(x),
             VIX_q25 = unname(quantile(x, .25)), VIX_q75 = unname(quantile(x, .75)),
             VIX_ge_30 = sum(x >= 30), fraction_VIX_ge_30 = mean(x >= 30))
})), "EXTERNAL_VIX_BY_STATE")
save_table(data.frame(dates = length(days), start = min(days), end = max(days),
                      spearman_state_vs_vix = cor(model$state, vix_value, method = "spearman"),
                      spearman_filtered_stress_vs_vix = cor(model$filtered[, "Stress"], vix_value,
                                                            method = "spearman"),
                      filtered_smoothed_mode_agreement =
                        mean(model$state == max.col(model$smoothed, ties.method = "first"))),
           "EXTERNAL_VIX_SUMMARY")

vol_paths <- function(R) {
  cv <- revision_covariance(R, spec, pr$warmup)
  X <- esg_features(cv$H)$X[cv$valid, , drop = FALSE]
  colnames(X) <- colnames(R)
  X
}
matched_tables <- function(raw_path, adjusted_path, names_of_nodes, panel, hmm = model) {
  blocks <- esg_config("pilot", root)$blocks
  stopifnot(identical(colnames(raw_path), names_of_nodes),
            identical(colnames(adjusted_path), names_of_nodes))
  paths <- list(raw = raw_path, adjusted = adjusted_path)
  lapply(paths, function(X) {
    lagged <- esg_lag_packets(X, hmm$state, 1L, days)
    stopifnot(all(diff(lagged$time) == 1L))
    fits <- esg_var_from_packets(lagged, hmm$k, min_obs = 40L, lambda = .1)
    lapply(esg_conn_tables(fits, hmm$labels, names_of_nodes, blocks, horizons = c(1L, 10L)),
           function(x) data.frame(panel = panel, x))
  })
}
paired_tci <- function(x) {
  merged <- merge(x$raw$tci, x$adjusted$tci,
                  by = c("panel", "regime", "horizon", "genuine_lag_pairs"),
                  suffixes = c("_raw", "_adjusted"))
  stopifnot(nrow(merged) == nrow(x$raw$tci), nrow(merged) == nrow(x$adjusted$tci))
  merged$raw_minus_adjusted <- merged$TCI_raw - merged$TCI_adjusted
  merged
}
raw_path <- vol_paths(raw)
adjusted_path <- vol_paths(adjusted$R)
main <- matched_tables(raw_path, adjusted_path, colnames(raw), "SPY")
for (kind in names(main$raw)) {
  save_table(rbind(data.frame(representation = "raw", main$raw[[kind]]),
                   data.frame(representation = "adjusted", main$adjusted[[kind]])),
             paste0("MATCHED_NETWORK_", toupper(kind)))
}
primary <- read_project("outputs/evidence/MATCHED_TCI_COMPONENTS.csv")
aligned <- merge(paired_tci(main), primary, by.x = c("regime", "horizon"),
                 by.y = c("state", "horizon"))
stopifnot(nrow(aligned) == 6L, all(abs(aligned$TCI_raw - aligned$raw) < 1e-9),
          all(abs(aligned$TCI_adjusted - aligned$adjusted) < 1e-9))
state_sensitivity <- lapply(2:4, function(k) {
  fit <- hmm_fits[[as.character(k)]]
  if (inherits(fit, "error") || !fit$converged) {
    stop("State-count check needs a converged HMM for ", k, " states.")
  }
  x <- if (k == 3L) main else matched_tables(raw_path, adjusted_path, colnames(raw), "SPY", fit)
  data.frame(states = k, paired_tci(x))
})
save_table(do.call(rbind, state_sensitivity), "MATCHED_HMM_STATE_SENSITIVITY")

# Both comparison panels below have nine nodes. Fit both controls on past
# returns and remove their same-day fitted contribution descriptively.
factors <- c("SPY", "XLK")
assets <- setdiff(colnames(sample$R), factors)
dual_raw <- sample$R[adjusted$rows, assets, drop = FALSE]
dual_residual <- dual_raw
for (j in seq_along(adjusted$rows)) {
  t <- adjusted$rows[j]
  past <- seq.int(t - 252L, t - 1L)
  F <- cbind(1, sample$R[past, factors, drop = FALSE])
  fit <- lm.fit(F, sample$R[past, assets, drop = FALSE])
  stopifnot(fit$rank == ncol(F), all(is.finite(fit$coefficients)), max(past) < t)
  dual_residual[j, ] <- dual_raw[j, ] -
    as.numeric(c(1, sample$R[t, factors]) %*% fit$coefficients)
}
dual <- matched_tables(vol_paths(dual_raw), vol_paths(dual_residual), assets, "SPY_XLK")
save_table(paired_tci(dual), "MARKET_AND_TECHNOLOGY_ADJUSTMENT")

grid <- expand.grid(lag = c(1L, 2L), ridge = c(.03, .10, .30))
sens <- lapply(seq_len(nrow(grid)), function(j) {
  lag <- grid$lag[j]
  penalty <- grid$ridge[j]
  packs <- lapply(list(raw_path, adjusted_path), esg_lag_packets,
                  state = model$state, p = lag, dates = days)
  rows <- which(packs[[1]]$time >= 3L)
  fits <- lapply(packs, esg_var_from_packets, k = 3L, min_obs = 40L,
                 lambda = penalty, row_index = rows)
  do.call(rbind, lapply(seq_along(labels), function(s) {
    f <- fits[[1]][[s]]
    g <- fits[[2]][[s]]
    tci <- function(fit) esg_conn(esg_fevd(fit, 10L), colnames(raw))$TCI
    data.frame(lag = lag, ridge = penalty, state = labels[s], n = f$n,
               raw_TCI = tci(f), adjusted_TCI = tci(g), drop = tci(f) - tci(g),
               raw_radius_before = f$radius_before, raw_radius_after = f$radius_after,
               adjusted_radius_before = g$radius_before, adjusted_radius_after = g$radius_after)
  }))
})
save_table(do.call(rbind, sens), "MATCHED_LAG_RIDGE_SENSITIVITY")
calendar_residuals <- lapply(seq_along(labels), function(s) {
  packet_index <- which(raw_fit$packet$state == s)
  pairs <- which(diff(raw_fit$packet$time[packet_index]) == 1L)
  E <- raw_fit$var[[s]]$residual
  correlations <- vapply(seq_len(ncol(E)), function(j) {
    cor(E[pairs, j], E[pairs + 1L, j])
  }, numeric(1))
  data.frame(state = labels[s], asset = colnames(raw), adjacent_pairs = length(pairs),
             lag1_residual_correlation = correlations)
})
save_table(do.call(rbind, calendar_residuals), "MATCHED_CALENDAR_RESIDUAL_CORRELATIONS")

cat("Review supplement: higher-confidence CVaR allocation at fixed origins\n")
panel <- revision_read(file.path(root, "input", "archived_returns.csv"))
R <- panel$R
dates <- panel$dates
cv <- revision_covariance(R, spec, pr$warmup)
ends <- esg_month_ends(dates)
origins <- ends[dates[ends] >= pr$development_end & ends < nrow(R)]
levels <- c(C95 = .95, C975 = .975, C99 = .99)
arms <- c(names(levels), "MinVar", "EqualWeight")
W <- array(NA_real_, c(length(origins), ncol(R), length(arms)),
           dimnames = list(as.character(dates[origins]), colnames(R), arms))
scenario_audit <- vector("list", length(origins))
for (j in seq_along(origins)) {
  t <- origins[j]
  ix <- tail(cv$valid[cv$valid <= t], pr$pool)
  target <- cv$H[t + 1L, , ]
  sc <- esg_match_scenarios(cv$eps[ix, , drop = FALSE], rep(1 / length(ix), length(ix)), target)
  for (name in names(levels)) {
    W[j, , name] <- esg_cvar_lp(sc$scenarios, sc$prob, levels[name], pr$cap)$weights
  }
  W[j, , "MinVar"] <- esg_minvar(target, pr$cap)
  W[j, , "EqualWeight"] <- rep(1 / ncol(R), ncol(R))
  scenario_audit[[j]] <- data.frame(origin = dates[t], first_payoff = dates[t + 1L],
                                   last_residual = dates[max(ix)], observations = length(ix),
                                   tail_vectors_975 = length(ix) * .025, tail_vectors_99 = length(ix) * .01)
}
stopifnot(all(is.finite(W)), all(W >= 0), all(W <= pr$cap + 1e-9))
save_table(do.call(rbind, scenario_audit), "HIGHER_CVAR_SCENARIO_AUDIT")
costs <- c(0, 5, 10, 25)
bts <- setNames(lapply(costs, function(cost) esg_backtest_targets(R, dates, origins, W, cost)), costs)
for (cost in as.character(costs)) {
  base <- read_project(paste0("outputs/historical_evaluation/returns_cost", cost, ".csv"))
  bt <- bts[[cost]]
  stopifnot(identical(as.character(bt$dates), base$date), nrow(base) == 663L,
            max(abs(bt$returns[, "C95"] - base$C)) < 1e-10,
            max(abs(bt$returns[, "MinVar"] - base$MinVar)) < 1e-10,
            max(abs(bt$returns[, "EqualWeight"] - base$EqualWeight)) < 1e-10)
  save_table(data.frame(date = bt$dates, bt$returns), paste0("HIGHER_CVAR_RETURNS_COST", cost))
}
performance <- do.call(rbind, lapply(costs, function(cost) {
  esg_performance_table(bts[[as.character(cost)]], cost)
}))
save_table(performance, "HIGHER_CVAR_PERFORMANCE")
returns <- bts[["10"]]$returns
dcc <- read_project("outputs/dcc_benchmark/returns_cost10.csv")
stopifnot(identical(dcc$date, as.character(bts[["10"]]$dates)))
returns <- cbind(returns, DCC_Gaussian_MinVar = dcc$DCC_Gaussian_MinVar)
contrast <- function(a, b, metric, fun, indices) {
  measure <- function(ix) fun(returns[ix, a]) - fun(returns[ix, b])
  data.frame(proposed = a, benchmark = b, metric = metric,
             t(esg_basic_inference(measure(seq_len(nrow(returns))), apply(indices, 2, measure))),
             check.names = FALSE)
}
tests <- do.call(rbind, lapply(pr$block_lengths, function(block) {
  idx <- esg_with_seed(esg_seed(c("review_tail", block)), replicate(
    pr$bootstrap_B, esg_mbb_indices(nrow(returns), block)
  ))
  es <- function(level) function(x) esg_cvar(-x, alpha = level)
  sharpe <- function(x) sqrt(252) * mean(x) / sd(x)
  z <- rbind(contrast("C975", "C95", "ES975", es(.975), idx),
             contrast("C99", "C95", "ES99", es(.99), idx),
             contrast("C95", "DCC_Gaussian_MinVar", "Sharpe", sharpe, idx),
             contrast("C95", "MinVar", "Sharpe", sharpe, idx))
  z$block <- block
  z$holm <- p.adjust(z$p_value, method = "holm")
  z
}))
save_table(tests, "HIGHER_CVAR_PAIRED_INFERENCE")

# Vector graphics for the manuscript, using the main analysis tables.
matched_boot <- read_project("outputs/evidence/MATCHED_MARKET_BLOCK_SENSITIVITY.csv")
ten <- primary[primary$horizon == 10L, ]
ten <- ten[match(labels, ten$state), ]
stopifnot(identical(ten$state, labels), all(ten$horizon == 10L))
intervals <- matched_boot[matched_boot$block == 10L, ]
intervals <- intervals[match(paste0("H10_", labels), intervals$contrast), ]
hold <- read_project("outputs/holding_horizon/performance_2023_cost10.csv")
hold_order <- c("H_OWN", "H_VAR", "H_STATE_OWN", "H_STATE_VAR")
net_es_bp <- 1e4 * (hold$ES_loss_95[match(hold_order, hold$strategy)] -
                     hold$ES_loss_95[hold$strategy == "C"])
navy <- "#23405A"
teal <- "#2A776F"
orange <- "#A96632"
pdf(file.path(output, "MATCHED_EVIDENCE_FIGURE.pdf"), width = 7.2, height = 5.6,
    family = "Helvetica", pointsize = 11)
par(mfrow = c(2, 2), mar = c(4.0, 4.3, 2.4, 0.5), mgp = c(2.3, 0.65, 0), tcl = -.24)
barplot(rbind(ten$raw, ten$adjusted), beside = TRUE, names.arg = labels,
        col = c(navy, teal), border = NA, ylim = c(0, 78),
        ylab = "Total connectedness (%)", main = "A. Matched SPY adjustment", cex.main = 1.02)
legend("topleft", c("Raw", "SPY adjusted"), fill = c(navy, teal), border = NA,
       bty = "n", cex = .85, horiz = TRUE)
plot(intervals$estimate, seq_along(labels), xlim = c(0, 42), ylim = c(.6, 3.4), yaxt = "n",
     xlab = "Raw minus adjusted (percentage points)", ylab = "", pch = 19, col = navy,
     main = "B. Conditional 95% intervals", cex.main = 1.02)
axis(2, at = seq_along(labels), labels = labels, las = 1, cex.axis = .84)
abline(v = 0, col = "gray45", lty = 3)
segments(intervals$lower, seq_along(labels), intervals$upper, seq_along(labels), col = teal, lwd = 2)
points(intervals$estimate, seq_along(labels), pch = 19, col = navy)
barplot(rbind(ten$innovation_covariance_component, ten$lag_matrix_component),
        names.arg = labels, col = c(teal, navy), border = NA, ylim = c(0, 32),
        ylab = "Difference (percentage points)", main = "C. Fitted components", cex.main = 1.02)
legend("topleft", c("Innovation covariance", "Lag matrix"), fill = c(teal, navy),
       border = NA, bty = "n", cex = .81)
barplot(net_es_bp, names.arg = c("Own", "VAR", "State\nown", "State\nVAR"),
        col = orange, border = NA, ylim = c(0, 10), cex.names = .88,
        ylab = "ES increase (basis points)", main = "D. Net development ES", cex.main = 1.02)
invisible(dev.off())

pdf(file.path(output, "MATCHED_REGIME_PROBABILITIES.pdf"), width = 7.2, height = 5.35,
    family = "Helvetica", pointsize = 11)
par(mfrow = c(3, 1), mar = c(2.7, 4, 1.0, 0.6), mgp = c(2.2, .62, 0), tcl = -.24)
for (s in seq_along(labels)) {
  plot(days, model$smoothed[, s], type = "l", col = "grey57", lwd = 1.3,
       ylim = c(0, 1), xlab = "", ylab = labels[s], axes = FALSE)
  lines(days, model$filtered[, s], col = teal, lwd = 1.1, lty = 2)
  axis.Date(1, at = seq(as.Date("2020-07-01"), as.Date("2023-07-01"), by = "6 months"),
            format = "%b %Y", cex.axis = .82)
  axis(2, at = c(0, .5, 1), las = 1)
  box(col = "grey68")
  if (s == 1L) {
    legend("topright", c("Filtered", "Smoothed"), col = c(teal, "grey57"),
           lty = c(2, 1), lwd = 1.3, bg = "white", cex = .84, horiz = TRUE)
  }
}
invisible(dev.off())
capture.output(sessionInfo(), file = file.path(output, "REVIEW_SUPPLEMENT_R_SESSION.txt"))
writeLines(c("Completed reviewer sensitivity analysis.",
             "State counts, VIX and adjusted networks are exploratory descriptions of the fitted prefix.",
             "Higher-CVaR rules retain the archived returns, calendar, model, cap, scenario design, and costs.",
             "Paired inference is a historical comparison without an untouched evaluation period."),
           file.path(output, "REVIEW_SUPPLEMENT_STATUS.txt"))
