# Usage: Rscript --vanilla holdings/compute_holdings_overlap.R [holdings_directory]
# The disclosed market values and report dates are frozen under data/.
args <- commandArgs(trailingOnly = TRUE)
root <- normalizePath(if (length(args)) args[1] else "holdings", mustWork = TRUE)
data_dir <- file.path(root, "data")
output <- file.path(root, "results")
dir.create(output, showWarnings = FALSE)
read_input <- function(name) {
  read.csv(file.path(data_dir, paste0(name, ".csv")), stringsAsFactors = FALSE, check.names = FALSE)
}
save_result <- function(x, name) {
  write.csv(x, file.path(output, paste0(name, ".csv")), row.names = FALSE, na = "")
}
positions <- read_input("equity_positions")
disclosed <- read_input("disclosed_industry_weights")
filings <- read_input("filing_register")
funds <- c("ICLN", "TAN", "QCLN", "XLK", "SPY")
rounds <- c("Spring", "Autumn")
sectors <- c("Information Technology", "Industrials", "Utilities", "Consumer Discretionary",
             "Materials", "Energy", "Real Estate", "Financials", "Health Care",
             "Communication Services", "Consumer Staples", "Unclassified")
stopifnot(nrow(filings) == length(funds) * length(rounds),
          !anyDuplicated(filings[c("round", "fund")]),
          all(positions$market_value_usd >= 0), all(nzchar(positions$issuer_key)),
          all(positions$report_date == filings$report_date[match(positions$source_file, filings$source_file)]))

reconcile <- aggregate(market_value_usd ~ round + fund + report_date, positions, sum)
check <- merge(filings, reconcile, by = c("round", "fund", "report_date"))
stopifnot(nrow(check) == nrow(filings), all(check$market_value_usd == check$equity_value_usd))

issuer <- aggregate(market_value_usd ~ round + fund + report_date + issuer_key, positions, sum)
issuer <- issuer[issuer$market_value_usd > 0, ]
issuer$weight <- issuer$market_value_usd /
  ave(issuer$market_value_usd, issuer$round, issuer$fund, FUN = sum)
save_result(issuer, "issuer_weights")

sector_weights <- do.call(rbind, lapply(seq_len(nrow(filings)), function(j) {
  period <- filings$round[j]
  fund <- filings$fund[j]
  reported <- disclosed[disclosed$round == period & disclosed$fund == fund, ]
  if (nrow(reported)) {
    reported <- reported[reported$sector != "Non-equity residual", ]
    raw <- tapply(reported$reported_percent, reported$sector, sum)
    source <- "Published industry percentages"
    stopifnot(abs(sum(raw) - 100) <= 1)
  } else {
    pos <- positions[positions$round == period & positions$fund == fund, ]
    stopifnot(all(nzchar(pos$sector)))
    raw <- tapply(pos$market_value_usd, pos$sector, sum)
    source <- if (fund == "SPY") "Same-date Select Sector SPDR membership" else "Filing equity schedule"
  }
  stopifnot(all(names(raw) %in% sectors), all(raw >= 0), sum(raw) > 0)
  x <- setNames(numeric(length(sectors)), sectors)
  x[names(raw)] <- raw / sum(raw)
  data.frame(round = period, fund = fund, report_date = filings$report_date[j],
             sector = names(x), weight = unname(x), classification_source = source)
}))
stopifnot(all(abs(aggregate(weight ~ round + fund, sector_weights, sum)$weight - 1) < 1e-12))
save_result(sector_weights, "sector_weights")

pair_rows <- list()
lines <- list()
for (period in rounds) {
  for (pair in combn(funds, 2, simplify = FALSE)) {
    align <- function(fund) issuer[issuer$round == period & issuer$fund == fund,
                                   c("issuer_key", "weight")]
    shared <- merge(align(pair[1]), align(pair[2]), by = "issuer_key", suffixes = c("_a", "_b"))
    shared$matched_weight <- pmin(shared$weight_a, shared$weight_b)
    rows <- sector_weights[sector_weights$round == period, ]
    a <- setNames(rows$weight[rows$fund == pair[1]], rows$sector[rows$fund == pair[1]])
    b <- setNames(rows$weight[rows$fund == pair[2]], rows$sector[rows$fund == pair[2]])
    known <- setdiff(sectors, "Unclassified")
    bound <- sum(pmin(a[known], b[known]))
    uncertain <- a[["Unclassified"]] + b[["Unclassified"]]
    date_for <- function(fund) filings$report_date[filings$round == period & filings$fund == fund]
    pair_rows[[length(pair_rows) + 1L]] <- data.frame(
      round = period, fund_a = pair[1], fund_b = pair[2],
      date_a = date_for(pair[1]), date_b = date_for(pair[2]),
      date_gap_days = abs(as.numeric(as.Date(date_for(pair[1])) - as.Date(date_for(pair[2])))),
      shared_issuers = nrow(shared), issuer_overlap = sum(shared$matched_weight),
      weight_a_in_b = sum(shared$weight_a), weight_b_in_a = sum(shared$weight_b),
      sector_overlap_lower = bound, sector_overlap_upper = min(1, bound + uncertain)
    )
    if (nrow(shared)) {
      lines[[length(lines) + 1L]] <- data.frame(round = period, fund_a = pair[1],
                                                fund_b = pair[2], shared)
    }
  }
}
overlaps <- do.call(rbind, pair_rows)
stopifnot(nrow(overlaps) == 20, all(overlaps$date_gap_days <= 61),
          all(overlaps$issuer_overlap <= pmin(overlaps$weight_a_in_b, overlaps$weight_b_in_a) + 1e-12),
          all(overlaps$sector_overlap_lower <= overlaps$sector_overlap_upper + 1e-12))
save_result(overlaps, "pairwise_overlap")
save_result(do.call(rbind, lines), "shared_issuers")
save_result(check, "equity_total_reconciliation")

draw_overlap <- function(period, measure, heading) {
  data <- overlaps[overlaps$round == period, ]
  values <- matrix(NA_real_, length(funds), length(funds), dimnames = list(funds, funds))
  for (k in seq_len(nrow(data))) {
    i <- match(data$fund_a[k], funds)
    j <- match(data$fund_b[k], funds)
    values[i, j] <- 100 * data[[measure]][k]
  }
  shades <- colorRampPalette(c("#F1F5F4", "#88BBB1", "#216860", "#173842"))(101)
  plot.new()
  plot.window(c(.4, 5.7), c(.4, 5.7), xaxs = "i", yaxs = "i", asp = 1)
  title(main = heading, line = .35, cex.main = .95)
  for (i in seq_along(funds)) {
    for (j in seq_along(funds)) {
      value <- values[i, j]
      if (!is.na(value)) {
        rect(j - .49, 6 - i - .49, j + .49, 6 - i + .49,
             border = NA, col = shades[1 + round(value)])
        text(j, 6 - i, formatC(value, digits = 1, format = "f"),
             cex = 1, col = if (value >= 52) "white" else "#162928")
      }
    }
  }
  axis(1, at = 1:5, labels = funds, tick = FALSE, line = -.7, cex.axis = .8)
  axis(2, at = 5:1, labels = funds, tick = FALSE, line = -.5, las = 1, cex.axis = .8)
}
pdf(file.path(output, "holdings_sector_overlap.pdf"), width = 7.2, height = 6, pointsize = 11, useDingbats = FALSE)
par(mfrow = c(2, 2), mar = c(1.8, 2.5, 1.55, .2), oma = c(.4, .1, .2, .1))
draw_overlap("Spring", "issuer_overlap", "A. Issuer overlap, Apr-Jun")
draw_overlap("Autumn", "issuer_overlap", "B. Issuer overlap, Sep-Nov")
draw_overlap("Spring", "sector_overlap_lower", "C. Sector overlap, Apr-Jun")
draw_overlap("Autumn", "sector_overlap_lower", "D. Sector overlap, Sep-Nov")
invisible(dev.off())
cat(nrow(positions), "equity lines across", nrow(filings), "public filings reconciled exactly.\n")
cat(nrow(overlaps), "pairs across two filing rounds.\n")
