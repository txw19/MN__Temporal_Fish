require(jagsUI)
require(doBy)

########################################################
# Format data
dat <- read.csv("MN_fish_timeseries.csv") 
head(dat)

## Clean data
# Sites with no DD data (some sites may have some missing DD)
noTemp <- aggregate(dat$gdd_wtr_5c,list(dat$DOW), mean, na.rm=T)
withTemp <- noTemp[!is.na(noTemp$x),]

# Only use sites that have modeled temp data
dat <- dat[dat$DOW %in% withTemp$Group.1,]
summary(dat)

# Remove sites with zero lake area
noArea <- aggregate(dat$LAKE_AREA_DOW_ACRES,list(dat$DOW), mean, na.rm=T)
withArea <- noArea[noArea$x!=0,]

dat <- dat[dat$DOW %in% withArea$Group.1,]
summary(dat)

# Select out lakes with n >= X
tbl1 <- table(dat$DOW)
dat <- dat[dat$DOW %in% names(tbl1)[tbl1 >=10],]

####### DATA PREP
Sites <- unique(dat$DOW)
nSites <- length(Sites)
nYears <- length(unique(dat$Year))
Years <- as.data.frame(seq(min(dat$Year),max(dat$Year),1))
Years$idx <- NA
colnames(Years) <- c("Year","idx")

# dat[dat$DOW==1000100,]

# Create matrix [i,t]: Response variable
ys <- matrix(nrow=nSites, ncol=nYears)
yx <- dat[c("DOW", "Year", "YEP")]
for(i in 1:nSites){
  siteidx <- Sites[i]
  yx1 <- subset(yx, DOW==siteidx)
  yx2 <- merge(Years, yx1, by="Year", all=T)
  ys[i,] <- t(yx2$YEP)
}
ys <- log(ys+0.1) # log-transform
rm(yx, yx1, yx2)


# Create DD predictor matrix [i,t]
dd <- matrix(nrow=nSites, ncol=nYears)
ddx <- dat[c("DOW", "Year", "gdd_wtr_5c")]
for(i in 1:nSites){
  siteidx <- Sites[i]
  ddx1 <- subset(ddx, DOW==siteidx)
  ddx2 <- merge(Years, ddx1, by="Year", all=T)
  dd[i,] <- t(ddx2$gdd_wtr_5c)
}
dd <- scale(dd) # z-scores by year
dd[is.na(dd)] <- 0 # assuming mean values for unsampled years
rm(ddx, ddx1, ddx2)

# Create Walleye predictor matrix [i,t]
wae <- matrix(nrow=nSites, ncol=nYears)
waex <- dat[c("DOW", "Year", "WAE")]
for(i in 1:nSites){
  siteidx <- Sites[i]
  waex1 <- subset(waex, DOW==siteidx)
  waex2 <- merge(Years, waex1, by="Year", all=T)
  wae[i,] <- t(waex2$WAE)
}
wae <- log(wae + 0.1) # log-transform
wae <- scale(wae) # z-scores by year
wae[is.na(wae)] <- 0 # assuming mean values for unsampled years
rm(waex, waex1, waex2)


# Create lake-level covariate vectors [i]: elevation, basin size, gradient
area <- as.numeric(by(dat$LAKE_AREA_DOW_ACRES, dat$DOW, mean))
area <- log(area)
area <- scale(area)
depth <- as.numeric(by(dat$MAX_DEPTH_FEET, dat$DOW, mean))
depth <- log(depth)
depth <- scale(depth)

########################################################
# Compile data
data <- list(ys=ys, # YEP CPE [i, t]
             dd=dd, # DD [i, t]
             wae=wae, # WAE CPE [i, t]
             area=area, # Lake area [i]
             depth=depth, # Lake depth [i]
             nSites=nSites, # count of sites
             nYears=nYears) # count of years

########################################################
# Write model
writeLines("
           model{
           for(i in 1:nSites){
              for(t in 1:nYears){
                 ys[i,t] ~ dnorm(mu[i,t], tau)

                    mu[i,t] <-  b1[i]*dd[i,t] + b2[i]*wae[i,t] + site.ran[i]
              }
           }


  # Priors: random slopes 
    for(i in 1:nSites){
      b1[i] ~ dnorm(mu.b, tau.b)
      b2[i] ~ dnorm(mu.b2, tau.b2)
    }
    # hyperpriors
    mu.b ~ dnorm(0, 0.001)
    tau.b <- pow(sigma.b, -2)
    sigma.b ~ dunif(0, 10)
    mu.b2 ~ dnorm(0, 0.001)
    tau.b2 <- pow(sigma.b, -2)
    sigma.b2 ~ dunif(0, 10)
  

          for(i in 1:nSites){
            site.ran[i] ~ dnorm(site.hat[i], tau.site)
              site.hat[i] <- alpha + g1 * area[i] + g2*depth[i] 

          }

           # Priors
            tau.site <- pow(sigma.site, -2)
            sigma.site ~ dunif(0,5)

           alpha ~ dunif(-5,5)
           g1 ~ dnorm(0, 0.0001)
           g2 ~ dnorm(0, 0.0001)
           sigma ~ dunif(0, 10)
           tau <- pow(sigma, -2)
           }
           ", con="model.txt")
modfile <- 'model.txt'

########################################################
# Initial values
inits <- function() {
  list(alpha=rnorm(1,0,1), 
       g1=rnorm(1,0,1), 
       g2=rnorm(1,0,1), 
       sigma=runif(1,0,3),
       mu.b=rnorm(1),
       mu.b2=rnorm(1),
       sigma.b = runif(1),
       sigma.b2 = runif(1)
       )
}

########################################################
# Set parameters to monitor
params <- c("b1","b2","g1","g2","alpha","sigma", "sigma.site","mu.b","mu.b2","sigma.b","sigma.b2") # ,"site.ran"
params.names <- c("Degree Days","Walleye CPE","Lake Area","Lake Depth","Alpha","Sigma","Sigma Site")

########################################################
# Run analysis
out.tw <- jags(data=data,
            inits=inits,
            parameters.to.save=params,
            model.file=modfile,
            n.chains=3,
            n.adapt=100,
            n.iter=12000,
            n.burnin=10000,
            n.thin=2)
save(out.tw, file="out_tw.R")


out2 <- out.tw$summary
# write.csv(out2, 'tw.out.csv', row.names = T)

which(out.tw$summary[, c("Rhat")] > 1.1)

summary(out.tw)

########################################################
# Posterior distribution plots
# Number of mcmc samples
n.keep <- out.tw$mcmc.info$n.samples
# Grab all three chains
output <- out.tw$samples
# Combine all three chains
output2 <- rbind(output[[1]], output[[2]], output[[3]])



par(mfrow=c(5,3), mar=c(3,4,1,0), oma=c(1,1,1,1), las=1)
for(i in 1:7){
  hist(output2[,i], main=params.names[i], xlab="", col="black")
  abline(v=0, col="red", lty=2)
}


####### MORE DETAILED PLOTS OF POSTEROR DISTRIBUTIONS

# 150 parameters 
# Number of pags
n.plot <- 18
# Number of figures per page
n.panel <- 20

# Create names for random site effects for plots
siteRan <- numeric()
for(i in 1:length(elev)){
  siteRan[i] <- paste("site",i)
}
# Combine random effect names with params.names
names <- c(params.names,siteRan,"Deviance")

# 333 + 17

## Create plots of posterior densities
for(hh in 1:n.plot){
  pdf(paste('yoy_size',hh,'.pdf'),width=8,height=12)
  

  if(hh<n.plot){
    n.fig=n.panel #number of pannels per plot
    layout(matrix(1:n.fig,nrow=5,ncol=4,byrow=T))
    par(mar=c(4,4,1,1),oma=c(2,2,1,1))
  }

  if(hh==n.plot){
    n.fig=ncol(output[[1]])-(n.plot-1)*n.panel
    layout(matrix(c(1:n.fig,rep(0,(n.panel-n.fig))),nrow=5,ncol=4,byrow=T))
    par(mar=c(4,4,1,1),oma=c(2,2,1,1))
  }
  
  for(pp in 1:n.fig){
    den.1=density(output[[1]][,(hh-1)*n.panel+pp])
    den.2=density(output[[2]][,(hh-1)*n.panel+pp])
    den.3=density(output[[3]][,(hh-1)*n.panel+pp])
    
    x.min=min(output2[,(hh-1)*n.panel+pp])
    x.max=max(output2[,(hh-1)*n.panel+pp])
    y.max=max(max(den.1$y),max(den.2$y),max(den.3$y))
    
    plot(den.1,type='l',lwd=1,col='red',xlim=c(x.min,x.max),ylim=c(0,y.max),ylab='',xlab=names[(hh-1)*n.panel+pp],main='')
    lines(den.2,type='l',lwd=1,col='green')
    lines(den.3,type='l',lwd=1,col='blue') #3 chains in total
    abline(v=0, lwd=1.5, col='blue', lty=2)
  }
  mtext('Density',2,outer=T)
  
  dev.off()
}


