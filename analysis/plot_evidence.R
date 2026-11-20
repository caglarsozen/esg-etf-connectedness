local({
  read <- function(name) read.csv(file.path('outputs/evidence',name))
  matched <- read('MATCHED_MARKET_BLOCK_SENSITIVITY.csv')
  h10 <- matched[matched$block==10 & grepl('^H10_',matched$contrast),]
  h10$state <- sub('H10_','',h10$contrast)
  parts <- read('MATCHED_TCI_COMPONENTS.csv')
  parts <- parts[parts$horizon==10,]
  allocation <- read.csv('outputs/holding_horizon/performance_2023_cost10.csv')
  stopifnot(identical(h10$state,parts$state),nrow(h10)==3,nrow(allocation)==5)
  dir.create('outputs/figures',showWarnings=FALSE)
  navy <- '#243E55'
  teal <- '#388E82'
  bronze <- '#BB7A46'
  render <- function() {
    par(mfrow=c(2,2),mar=c(4.2,5,3.2,1.5),las=1,mgp=c(2.8,.7,0))
    barplot(rbind(h10$raw_TCI,h10$market_adjusted_TCI),beside=TRUE,
      col=c(navy,teal),border=NA,ylim=c(0,77),names.arg=h10$state,
      ylab='Total connectedness (%)',main='A. Matched SPY adjustment')
    legend('topright',c('Raw','SPY adjusted'),fill=c(navy,teal),border=NA,
      bty='n',cex=.8,inset=c(0,.01))
    plot(h10$estimate,1:3,xlim=range(c(h10$lower,h10$upper,0)),ylim=c(.55,3.45),
      yaxt='n',pch=19,col=navy,xlab='Raw minus adjusted (percentage points)',
      ylab='',main='B. Pointwise conditional 95% intervals')
    axis(2,at=1:3,labels=h10$state,las=1,cex.axis=.82)
    abline(v=0,lty=2,col='grey50')
    segments(h10$lower,1:3,h10$upper,1:3,lwd=3,col=teal)
    points(h10$estimate,1:3,pch=19,col=navy)
    composition <- rbind(parts$innovation_covariance_component,parts$lag_matrix_component)
    barplot(composition,names.arg=parts$state,col=c(teal,navy),border=NA,
      ylim=c(0,30),ylab='Difference (percentage points)',
      main='C. Components of the fitted difference')
    legend('topright',c('Innovation covariance','VAR lag matrix'),fill=c(teal,navy),
      border=NA,bty='n',cex=.8)
    difference <- 1e4*(allocation$ES_loss_95-allocation$ES_loss_95[allocation$strategy=='C'])
    barplot(difference,col=c(navy,rep(bronze,4)),border=NA,
      names.arg=c('C','Own','VAR','State\nOwn','State\nVAR'),cex.names=.8,
      ylab='ES95 difference from C (bp)',main='D. Net development allocation',
      ylim=c(0,max(difference)*1.16))
    mtext('Daily loss at 10 bp trading cost',side=1,line=3.25,cex=.72)
  }
  png('outputs/figures/connectedness_evidence.png',width=3000,height=2150,res=245)
  render()
  dev.off()
  pdf('outputs/figures/connectedness_evidence.pdf',width=12,height=8.6,useDingbats=FALSE)
  render()
  dev.off()
})
