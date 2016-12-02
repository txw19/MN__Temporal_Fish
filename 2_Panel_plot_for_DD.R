########## SCATTER PLOT 

# ###############################################################
# ##### Plot Popl'n average effect [** not plotted **]
# 
# predX <- seq(min(dat$length), max(dat$length), length=100) # fake data to predict
# bayesB <- out1$BUGSoutput$summary[c("mu.a","mu.b"),1] # Extract pop'n average coefficents
# zz <- seq(min(dat$wgt,na.rm=T),max(dat$wgt,na.rm=T),length=100) # fake y-axis
# GroupCoef <- matrix(out1$BUGSoutput$summary[1:(2*(length(unique(dat$Pool)) )),1], c(length(unique(dat$Pool)),2), byrow=F)
# 

### Population average effect
##############################
predX <- seq(min(dat2$Zdd), max(dat2$Zdd), length=25) # fake data to predict

# Obtain 90% CIs for fitted line
est.lineA <- matrix(NA, ncol=length(predX), nrow=out1$BUGSoutput$n.sims) #container for predicted values

for(i in 1:out1$BUGSoutput$n.sims ){
  for(t in 1:length(predX) ){
    est.lineA[i,t] <- out1$BUGSoutput$sims.list$mu.a[i] + out1$BUGSoutput$sims.list$mu.b[i] * predX[t]

  }
}

# 95% CIs and fitted values
y.pred <- apply(est.lineA, 2, mean )
upper.CIA <- apply(est.lineA, 2, quantile, probs=c(0.975) )
lower.CIA <- apply(est.lineA, 2, quantile, probs=c(0.025) )


#################################
### Group-specific
#################################

# Container for predicted values
est.lineB <- array(NA, c(out1$BUGSoutput$n.sims,length(predX),J) )

# Put each groups MCMC draws for all 2 parameters in its own list
group.params <- list()
for(m in 1:J){
  group.params[[m]] <- out1$BUGSoutput$sims.list$BB[,m,]
}


for(k in 1:J){ # loop over groups (J)
  for(i in 1:out1$BUGSoutput$n.sims ){  
    for(t in 1:length(predX)){
      est.lineB[i,t,k] <-  group.params[[k]][i,1] + group.params[[k]][i,2] * predX[t]
    }	  
  }
}

groupMean <- array(NA, c(1,length(predX),J) )
upper.CIB <- array(NA, c(1,length(predX),J) )
lower.CIB <- array(NA, c(1,length(predX),J) )

for(i in 1:J){
  
  # Means
  groupMean[,,i] <- apply(est.lineB[,,i], 2, mean )
  # 95% CIs for fitted values
  upper.CIB[,,i] <- apply(est.lineB[,,i], 2, quantile, probs=c(0.975) )
  lower.CIB[,,i] <- apply(est.lineB[,,i], 2, quantile, probs=c(0.025) )
}



################# PLOT ############################

res <- 6
name_figure <- "DD_nested.png"
png(filename = name_figure, height = 500*res, width = 800*res, res=72*res)
def.par <- par(no.readonly = TRUE)

size.labels = 1
size.text = 1
axissize <- 1
x.label = "Standardized DD"
y.label = expression(paste(log[e],'(CPE)+0.1' ))

z <- seq(min(dat2$logYEP,na.rm=T),max(dat2$logYEP,na.rm=T),length=25) 

nf <- layout(matrix(c(1:100),nrow=10,ncol=10,byrow=TRUE),  TRUE) 
layout.show(nf)
par(mar=c(0.0,0.1,0.1,0.1),oma=c(3,3,0,1),mai=c(0.0,0.05,0.05,0) )


# Group-specific plots

for(i in 1:J){
  
  plot(predX,z, ylim=c(min(dat2$logYEP,na.rm=T),max(dat2$logYEP,na.rm=T)),
       xlim=c(min(dat2$Zdd),max(dat2$Zdd)), axes=F, ylab='', xlab='', type='n')
  

  if( i <=88){
    axis(side=1,cex.axis=axissize , mgp=c(1,0,0),tck= -0.01, labels=F ) 
  } else {
    axis(side=1,cex.axis=axissize , mgp=c(1,0,0),tck= -0.01)
  }	
  
  if( i ==1 | i==11 |i==21 |i==31|i==41|i==51|i==61|i==71|i==81|i==91){
    axis(side=2,cex.axis=axissize , mgp=c(0,0.3,0),tck= -0.01, las=1)
  } else {
    axis(side=2,cex.axis=axissize , mgp=c(1,0,0),tck= -0.01, labels=F)
  }	
  
  # Add credible region
  i.for <- order(predX)
  i.back <- order(predX, decreasing = TRUE )
  x.polygon <- c( predX[i.for] , predX[i.back] )
  y.polygon <- c( lower.CIB[,,i][i.for] , upper.CIB[,,i][i.back] )
  polygon( x.polygon , y.polygon , col = "lightgray" , border = NA )
  
  points(dat2$Zdd[dat2$G==i], dat2$logYEP[dat2$G==i], cex=0.8, pch=16,col="black" )
  
  
  # Add posterior means
  lines(predX, groupMean[,,i],lwd=1, col='black',lty=1)
  
  # text(0.1,0.9,i,cex=0.8)
  #text(0.3,0.9,paste('(',sum2$bkt.n[i],')'),cex=0.8)
  box()
  
}

mtext(y.label, line = 1.5, side = 2, cex = size.text,outer=T)
mtext(x.label, line = 1.5, side = 1, cex = size.text, outer=T)


plot(predX,z, ylim=c(min(dat2$logYEP,na.rm=T),max(dat2$logYEP,na.rm=T)),
     xlim=c(min(dat2$Zdd),max(dat2$Zdd)), axes=F, ylab='', xlab='', type='n')
axis(side=1,cex.axis=axissize , mgp=c(1,0,0),tck= -0.01)
# Add credible region
i.for <- order(predX)
i.back <- order(predX, decreasing = TRUE )
x.polygon <- c( predX[i.for] , predX[i.back] )
y.polygon <- c( lower.CIA[i.for] , upper.CIA[i.back] )
polygon( x.polygon , y.polygon , col = "lightgray" , border = NA )

points(dat2$Zdd, dat2$logYEP, cex=0.8, pch=16,col="black" )

lines(predX, y.pred,lwd=2, lty=1, col='blue')
box()

par(def.par)
dev.off()
