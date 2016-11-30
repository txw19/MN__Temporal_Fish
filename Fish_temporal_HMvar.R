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
ys <- log(ys+1) # log-transform
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
wae <- log(wae + 1) # log-transform
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
             area=as.numeric(area), # Lake area [i]
             depth=as.numeric(depth), # Lake depth [i]
             nSites=nSites, # count of sites
             nYears=nYears) # count of years

########################################################
# Write model
writeLines("
           model{
           for(i in 1:nSites){
              for(t in 1:nYears){
                 ys[i,t] ~ dnorm(mu[i,t], tau[i])

                    mu[i,t] <- alpha + b1*dd[i,t] + b2*wae[i,t] + site.ran[i]
              }
           }


          for(i in 1:nSites){
          # Site random effect
            site.ran[i] ~ dnorm(site.hat[i], tau.site)
            site.hat[i] <- g1 * area[i] + g2*depth[i] 

          # Hierarchical variance model
          log.sigma[i] ~ dnorm(mu.sigma[i], inv.omega.sigma.squared)
          log(sigma[i]) <- log.sigma[i]
          tau[i] <- 1/pow(sigma[i], 2)
           
          mu.sigma[i] <- Mu.Sigma # can add predictors of the sd here
          med.sigma[i] <- exp(mu.sigma[i])

          }

           # Priors
            tau.site <- pow(sigma.site, -2)
            sigma.site ~ dunif(0,1)

           alpha ~ dnorm(0, 0.0001)
           b1 ~ dnorm(0, 0.0001)
           b2 ~ dnorm(0, 0.0001)
           g1 ~ dnorm(0, 0.0001)
           g2 ~ dnorm(0, 0.0001)


           Mu.Sigma ~ dnorm(0,0.001)
           Med.Sigma <- exp(Mu.Sigma)
           omega.sigma ~ dunif(0, 1)
           inv.omega.sigma.squared <- 1 / pow(omega.sigma, 2)
           }
           ", con="model.txt")
modfile <- 'model.txt'

########################################################
# Initial values
inits <- function() {
  list(alpha=rnorm(1), 
       b1=rnorm(1,0,1), 
       b2=rnorm(1,0,1), 
       g1=rnorm(1,0,1), 
       g2=rnorm(1,0,1), 
       Mu.Sigma=rnorm(1), 
       omega.sigma=runif(1),
       log.sigma=runif(nSites) )
}

########################################################
# Set parameters to monitor
params <- c("b1","b2","g1","g2","alpha", "sigma.site","Mu.Sigma") # "sigma","site.ran"
params.names <- c("Degree Days","Walleye CPE","Lake Area","Lake Depth","Alpha","Sigma Site", "Mu.Sigma")

########################################################
# Run analysis
out_HMvar<- jags(data=data,
            inits=inits,
            parameters.to.save=params,
            model.file=modfile,
            n.chains=3,
            n.adapt=100,
            n.iter=7000,
            n.burnin=5000,
            n.thin=3)
save(out_HMvar, file="out_HMvar.R")


# out2 <- out_HMvar$summary
# write.csv(out2, 'tw.out.csv', row.names = T)

which(out_HMvar$summary[, c("Rhat")] > 1.1)

summary(out_HMvar)
summary(out.tw)

########################################################
# Posterior distribution plots
# xlims=c(-0.8,0.8)
# par(mfrow=c(5,3), mar=c(3,4,1,0), oma=c(1,1,1,1), las=1)
# 
# for(i in 1:(length(params)-2)){
#   hist(c(out$samples[[1]][,i], out$samples[[2]][,i], out$samples[[3]][,i]), main=params.names[i], xlab="", col="black", xlim=xlims)
#   abline(v=0, col="red", lty=2)
# }


# Number of mcmc samples
n.keep <- out.tw$mcmc.info$n.samples
# Grab all three chains
output <- out.tw$samples
# Combine all three chains
output2 <- rbind(output[[1]], output[[2]], output[[3]])

# 350 parameters 
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


