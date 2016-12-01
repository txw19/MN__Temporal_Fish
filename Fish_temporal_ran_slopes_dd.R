require(jagsUI)
require(doBy)
# rm(list=ls())
########################################################
# Format data
dat <- read.csv("MN_GN_CPUE_all_lakes.csv") 
head(dat)
length(unique(dat$DOW))

temp <- read.csv("temperature_data_for_ty.csv") 
head(temp)

wqdat <- read.csv("wq_data_for_ty.csv")
head(wqdat)

## Clean data
# Select out lakes with n >= X
tbl1 <- table(dat$DOW)
dat <- dat[dat$DOW %in% names(tbl1)[tbl1 >=8],]
length(unique(dat$DOW))
head(dat)
range(dat$Year)

# Select temp data for lakes in dat
temp <- temp[temp$DOW %in% dat$DOW,]
length(unique(temp$DOW))

# Convert ice-off date to decimal day of year
temp$ice_off_date <- as.POSIXct(temp$ice_off_date, format="%m/%d/%Y")
temp$ice_off_julian <- as.numeric(strftime(temp$ice_off_date, format = "%j"))
head(temp)

# Not all lakes in dat have temp, so remove lakes in dat that are not in temp
dat <- dat[dat$DOW %in% temp$DOW,]
length(unique(dat$DOW))

# Grab water quality data for lakes in dat
wqdat <- wqdat[wqdat$DOW %in% dat$DOW,]
length(unique(dat$DOW))

# Look at years contained in all datasets - match them up
range(dat$Year)
range(temp$Year)
range(wqdat$SAMPLE_YEAR)

# Remove sample year 2016 from dat - other data sets only go to 2015
# and start all data sets at 1987 (latest date for cpe data)
dat <- dat[dat$Year < 2016, ]
range(dat$Year)

temp <- temp[temp$Year > 1986, ]
range(temp$Year)

wqdat <- wqdat[wqdat$SAMPLE_YEAR > 1986, ]
range(wqdat$SAMPLE_YEAR)

# Sort data frames by DOW and Year for simplicity
dat <- dat[order(dat$DOW, dat$Year), ]
temp <- temp[order(temp$DOW, temp$Year), ]
wqdat <- wqdat[order(wqdat$DOW, wqdat$SAMPLE_YEAR), ]

# max(table(dat$DOW, dat$Year))

# Fill in missing years for all lakes with NA in dat
datY <- expand.grid(DOW = unique(dat$DOW), Year = min(dat$Year):max(dat$Year))
head(datY)

# Grab YEP data
yep <- dat[c("DOW", "Year", "YEP")]

# Merge yep with datY that has missing data for years with no data for each lake
yep <- merge(yep,datY,all=TRUE)
dim(yep)
head(yep)
length(unique(yep$Year))

# Grab WAE data
wae <- dat[c("DOW", "Year", "WAE")]
wae <- merge(wae,datY,all=TRUE)
dim(wae)
head(wae)
length(unique(wae$Year))

# Convert to matrices [i, t]
nSites <- length(unique(yep$DOW))
nYears <- length(unique(yep$Year))
# Matrices of CPE data
# Response variable
yepM <- matrix(yep$YEP, nrow=nSites, ncol=nYears, byrow=T)
# Log-tranform CPE
yepM <- log(yepM + 0.1)
# Predictor
waeM <- matrix(wae$WAE, nrow=nSites, ncol=nYears, byrow=T)
# Log-tranform CPE
waeM <- log(waeM + 0.1)
waeM[is.na(waeM)] <- 0 # assuming mean values for unsampled years

# Create matrices of temperature and water quality data
head(temp)
dd <- matrix(temp$gdd_wtr_5c, nrow=nSites, ncol=nYears, byrow=T)
dd <- scale(dd)
sum(is.na(dd)) # Check for any missing values (NAs)

iceoff <- matrix(temp$ice_off_julian, nrow=nSites, ncol=nYears, byrow=T)
iceoff <- scale(iceoff)
sum(is.na(iceoff))

iceduration <- matrix(temp$ice_duration_days, nrow=nSites, ncol=nYears, byrow=T)
iceduration <- scale(iceduration)
sum(is.na(iceduration))

coefvar <- matrix(temp$coef_var_30.60, nrow=nSites, ncol=nYears, byrow=T)
coefvar <- scale(coefvar)
sum(is.na(coefvar))

# Create lake-level covariate vectors [i]: elevation, basin size, gradient
area <- as.numeric(by(dat$LAKE_AREA_GIS_ACRES, dat$DOW, mean))
area <- log(area)
area <- as.numeric(scale(area))
depth <- as.numeric(by(dat$MAX_DEPTH_FEET, dat$DOW, mean))
depth <- log(depth)
depth <- as.numeric(scale(depth))

########################################################
# Compile data
data <- list(ys=yepM, # YEP CPE [i, t]
             dd=dd, # DD [i, t]
             nSites=nSites, # count of sites
             nYears=nYears) # count of years

########################################################
# Write model
writeLines("
           model{
           for(i in 1:nSites){
              for(t in 1:nYears){
                 ys[i,t] ~ dnorm(mu[i,t], tau)

                    mu[i,t] <-  b1[i]*dd[i,t] + site.ran[i]
              }
           }


  # Priors: random slopes 
    for(i in 1:nSites){
      b1[i] ~ dnorm(mu.b, tau.b)
    }
    # hyperpriors
    mu.b ~ dnorm(0, 0.001)
    tau.b <- pow(sigma.b, -2)
    sigma.b ~ dunif(0, 20)

  

          for(i in 1:nSites){
            site.ran[i] ~ dnorm(site.hat[i], tau.site)
              site.hat[i] <- alpha  

          }

           # Priors
            tau.site <- pow(sigma.site, -2)
            sigma.site ~ dunif(0,5)

           alpha ~ dunif(-5,5)
           g1 ~ dnorm(0, 0.0001)
           g2 ~ dnorm(0, 0.0001)
           sigma ~ dunif(0, 20)
           tau <- pow(sigma, -2)
           }
           ", con="model.txt")
modfile <- 'model.txt'

########################################################
# Initial values
inits <- function() {
  list(alpha=rnorm(1,0,1), 
       sigma=runif(1,0,3),
       mu.b=rnorm(1),
       sigma.b = runif(1)
       )
}

########################################################
# Set parameters to monitor
params <- c("b1","alpha","sigma", "sigma.site","mu.b","sigma.b","site.ran") # 
# params.names <- c("Degree Days","Walleye CPE","Lake Area","Lake Depth","Alpha","Sigma","Sigma Site")

########################################################
# Run analysis
out.tw <- jags(data=data,
            inits=inits,
            parameters.to.save=params,
            model.file=modfile,
            n.chains=3,
            n.adapt=100,
            n.iter=34000,
            n.burnin=30000,
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
head(output2)
dim(output2)
# Exclude deviance column
output2 <- output2[,-200]

# Posterior means
postmeans <- apply(output2, 2, mean)
cis <- apply(output2, 2, quantile, c(0.025, 0.975))

# for plotting
slopes <- postmeans[1:97]
int <- postmeans["alpha"]
ran.ints <- postmeans[103:199]
popslope <- postmeans["mu.b"]

res <- 6
name_figure <- "DD_effects.png"
png(filename = name_figure, height = 500*res, width = 800*res, res=72*res)
def.par <- par(no.readonly = TRUE)

size.labels = 1
size.text = 1
axissize <- 1
x.label = "Standardized DD"
y.label = expression(paste(log[e],'(CPE)' ))

nf <- layout(matrix(c(1:100),nrow=10,ncol=10,byrow=TRUE),  TRUE) 
layout.show(nf)
par(mar=c(0.0,0.1,0.1,0.1),oma=c(3,3,0,1),mai=c(0.0,0.05,0.05,0) )

# Group-specific plots

for(i in 1:100){
  
  plot(dd[i,],yepM[i,], ylim=c(min(yepM,na.rm=T),max(yepM,na.rm=T)),
       xlim=c(min(dd),max(dd)), axes=F, ylab='', xlab='', type='n')
  
  points(dd[i,],yepM[i,], cex=0.8, pch=16,col='black' )
  
  
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
  
  
  # Add posterior means
  abline((ran.ints[i]),slopes[i], lwd=2 )
  
  # text(0.1,0.9,i,cex=0.8)
  #text(0.3,0.9,paste('(',sum2$bkt.n[i],')'),cex=0.8)
  box()
  
}

mtext(y.label, line = 1.5, side = 2, cex = size.text,outer=T)
mtext(x.label, line = 1.5, side = 1, cex = size.text, outer=T)

plot(dd,yepM, ylim=c(min(yepM,na.rm=T),max(yepM,na.rm=T)),
     xlim=c(min(dd),max(dd)), axes=F, ylab='', xlab='', type='n')

points(dd,yepM, cex=0.8, pch=16,col='black' )
abline(int,popslope, lwd=2, col='blue' )

axis(side=1,cex.axis=axissize , mgp=c(1,0,0),tck= -0.01)

box()
par(def.par)
dev.off()

# 
# plot(1:29, yepM[1,])
# plot(1:29, dd[1,])

