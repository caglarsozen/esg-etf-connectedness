# ESG ETF Connectedness

R replication code for *Common Market Exposure, Volatility Connectedness, and Tail-Risk Allocation across ESG, Clean-Energy, and Conventional ETFs*. The analysis examines how shared equity exposure affects measured volatility connectedness and whether regime and network information improves expected shortfall after trading costs.

The main analysis estimates SPY exposure from the preceding 252 trading sessions and compares volatility networks derived from raw and adjusted returns for the same ten ETFs, trading dates, and regime labels. The repository also covers factor stochastic-volatility diagnostics, conditional value-at-risk (CVaR) allocation, Gaussian dynamic conditional correlation (DCC) benchmarks, and historical fund overlap using two sets of 2023 filings.

## Inputs

Put the archived daily return panel from the replication materials at `input/archived_returns.csv`. It contains 1,949 observations from 21 September 2018 through 25 June 2026. Use a wide CSV with unique trading dates in ascending order and simple decimal returns in the following column order.

```text
date,ESGU,SUSA,ESGV,ICLN,TAN,QCLN,SPY,XLE,XLK,GLD,TLT
```

Dates use `YYYY-MM-DD`. The archived panel was derived from Yahoo Finance adjusted closes and has SHA-256 `84053ad560ecfe5b85c3652be6ee86060665584b5b51d09decfb044a8c16bb19`.

For a run with separately obtained prices, put adjusted closes in `input/adjusted_prices.csv` using the same header. From the repository root, the following R commands convert them to daily returns on a common trading calendar. The archived panel is the input for reproducing the study estimates.

```r
source("R/01_core.R")
config <- esg_config()
stopifnot(file.exists(config$price_file), !file.exists("input/archived_returns.csv"))
panel <- esg_read_panel(config)
esg_csv(data.frame(date = panel$dates, panel$R, check.names = FALSE),
        "input/archived_returns.csv")
```

## Run the analysis

Run these commands from the repository root. The study used R 4.5.0. The dependency script installs `lpSolve`, `quadprog`, `factorstochvol`, and `posterior` from CRAN.

```sh
Rscript --vanilla install_dependencies.R
Rscript --vanilla run_all.R
Rscript --vanilla read_results.R
```

`run_baseline.R` runs the observation-driven analysis and historical benchmarks. `run_extended.R` runs the matched network, holding-horizon, and anchored Factor-SV analyses. The full run executes both and writes tables, figures, and diagnostic records under `outputs/`.

The long Factor-SV fits use four chains at each of two cutoffs. Each chain uses 40,000 burn-in iterations and retains 4,000 of the subsequent 160,000 draws. Set `ESG_WORKERS` to choose the number of workers. The fits reuse checkpoints with a matching input panel, specification, and software fingerprint. For the scaled, constant-idiosyncratic-variance comparison, set `ESG_RUN_SCALED_FSV=true` before running the baseline.

The archived fits used `factorstochvol` 1.1.2, `stochvol` 3.2.9, and `posterior` 1.6.1. The allocation results used `lpSolve` 5.6.23 and `quadprog` 1.5-8.

```sh
Rscript --vanilla tests/test_revision.R
Rscript --vanilla tests/test_horizon_forecasts.R
Rscript --vanilla tests/test_particle_prediction.R
```

The matched analysis uses data through 31 October 2023. Its bootstrap conditions on the fitted variance paths and common regime labels. The later portfolio evaluation is historical and uses a period examined during development.

## Reviewer checks and fund overlap

After the full run, pass the repository root, an output directory outside the repository, and the FRED VIXCLS CSV to generate the regime, network, and tail-level sensitivity tables. [FRED VIXCLS](https://fred.stlouisfed.org/series/VIXCLS)

```sh
Rscript --vanilla replication/reviewer_sensitivity.R . ../review_results input/VIXCLS.csv
```

For the holdings analysis, add `equity_positions.csv` and `disclosed_industry_weights.csv` from the replication materials to `holdings/data/`. The calculator reconciles disclosed equity totals and writes issuer and sector overlap tables and a figure to `holdings/results/`. Fund-specific report dates and source links are in the [filing register](holdings/data/filing_register.csv).

```sh
Rscript --vanilla holdings/compute_holdings_overlap.R
```

## License

The source code is distributed under the [MIT License](LICENSE). Market data and filing materials are subject to their respective source terms.
