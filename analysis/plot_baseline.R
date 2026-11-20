# Scientific evidence figure from the completed R tables; no model selection.
dir.create('outputs/figures',showWarnings=FALSE)
a<-read.csv('outputs/network_development/bootstrap_10.csv')
a<-a[a$statistic=='D_TCI',c('estimate','lower','upper')]
b<-read.csv('outputs/network_robustness/primary_sensitivities.csv')
b<-b[b$statistic=='D_TCI',c('estimate','lower','upper')]
n<-rbind(a,b)
p<-read.csv('outputs/historical_evaluation/paired_inference.csv')
p<-p[p$block==10&p$metric=='ES95'&p$benchmark=='C',]
stopifnot(nrow(n)==4,nrow(p)==5)
png('outputs/figures/evidence_boundaries.png',width=3300,height=1650,res=300)
par(mfrow=c(1,2),mar=c(5,7.6,3.4,1.4),oma=c(2,0,1,0),family='sans',las=1,
    cex=.82,col.axis='#344353',col.lab='#263746',fg='#344353')
yl<-4:1
plot(NA,xlim=c(-4,34),ylim=c(.5,4.5),xlab='Stress - Calm TCI (percentage points)',
     ylab='',yaxt='n',bty='n',main='Network contrast depends on the sample',cex.main=1.05)
abline(v=0,lty=2,col='#9ca8b2');abline(h=yl,col='#e8edf1')
cc<-c('#14557e',rep('#6c7d89',3))
segments(n$lower,yl,n$upper,yl,col=cc,lwd=2)
points(n$estimate,yl,pch=19,col=cc,cex=1.2)
axis(2,at=yl,labels=c('Full development\n2018-2023','Matched raw\n2019-2023',
                    'Without ESGV\n2019-2023','Market-adjusted\n2019-2023'),tick=FALSE,cex.axis=.84)
par(mar=c(5,5.4,3.4,1.4))
yp<-5:1
plot(NA,xlim=c(-5.5,10.5),ylim=c(.5,5.5),xlab='ES95 difference from C (basis points)',
     ylab='',yaxt='n',bty='n',main='Additional allocation arms versus baseline',cex.main=1.05)
abline(v=0,lty=2,col='#9ca8b2');abline(h=yp,col='#e8edf1')
segments(10000*p$lower,yp,10000*p$upper,yp,col='#477b84',lwd=2)
points(10000*p$estimate,yp,pch=19,col='#477b84',cex=1.2)
axis(2,at=yp,labels=c('Regime mixture','Network penalty','Worst component',
                    'Worst mixture','Worst mixture + network'),tick=FALSE,cex.axis=.84)
mtext('Historical re-evaluation: 663 days; 10 bps transaction charges.',side=3,line=.3,cex=.69)
mtext('Pointwise 95% conditional block intervals (block = 10). Negative ES differences favor the added arm.',
      side=1,outer=TRUE,line=.4,cex=.75)
dev.off()
