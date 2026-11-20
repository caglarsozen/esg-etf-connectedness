# Each day uses OLS coefficients estimated from the preceding 252 observations.
# Same-day SPY is removed for a descriptive network sensitivity, not a tradable return.
revision_market_residuals <- function(R, dates, window=252L, market='SPY') {
  stopifnot(market %in% colnames(R), nrow(R)>window, length(dates)==nrow(R))
  assets<-setdiff(colnames(R),market); tt<-(window+1L):nrow(R)
  E<-matrix(NA_real_,length(tt),length(assets),dimnames=list(NULL,assets))
  beta<-intercept<-E
  for(j in seq_along(tt)){
    t<-tt[j];past<-(t-window):(t-1L);x<-R[past,market];Y<-R[past,assets,drop=FALSE]
    xc<-x-mean(x);yc<-sweep(Y,2,colMeans(Y),'-')
    slope<-as.numeric(crossprod(xc,yc)/sum(xc^2))
    alpha<-colMeans(Y)-slope*mean(x)
    E[j,]<-R[t,assets]-alpha-slope*R[t,market]
    beta[j,]<-slope;intercept[j,]<-alpha
  }
  list(R=E,dates=dates[tt],rows=tt,beta=beta,intercept=intercept,
       latest_estimation_date=dates[tt-1L])
}
