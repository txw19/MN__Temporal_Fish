require(R2jags)
require(doBy)
library(MCMCpack) # rwish function
library(plyr)
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

dim(dat)
dim(temp)
dim(wqdat)

# Merge temp data with fish data
dat2 <- join(dat, temp, by=c("DOW", "Year"), type='left', match='all')
dim(dat2)

# log-transform YEP adn WAE data and standardize WAE
dat2$logYEP <- log(dat2$YEP+0.1)
dat2$logWAE <- as.numeric(scale(dat2$gdd_wtr_5c))

dat2$Zdd <- as.numeric(scale(dat2$gdd_wtr_5c))

# Create lake-level covariate vectors [i]: elevation, basin size, gradient
area <- as.numeric(by(dat$LAKE_AREA_GIS_ACRES, dat$DOW, mean))
area <- log(area)
area <- as.numeric(scale(area))
depth <- as.numeric(by(dat$MAX_DEPTH_FEET, dat$DOW, mean))
depth <- log(depth)
depth <- as.numeric(scale(depth))

#################################################################
########## BUGS CODE ############################################
#################################################################
sink("model_test.txt")
cat("
    model {
    for (i in 1:n){
    y[i] ~ dnorm (y.hat[i], tau.y)
    
    y.hat[i] <- alpha[group[i]] + beta[group[i]] * x[i]  
    
    }
    
    tau.y <- pow(sigma.y, -2)
    sigma.y ~ dunif (0, 10)
    
    # Level-2 of the model
    for(j in 1:J){
    alpha[j] <- BB[j,1]
    beta[j] <- BB[j,2]
    
    BB[j,1:K] ~ dmnorm (BB.hat[j,], Tau.B[,])
    BB.hat[j,1] <- mu.a 
    BB.hat[j,2] <- mu.b 
    
    }
    
    
    mu.a ~ dnorm(0,0.0001)
    mu.b ~ dnorm(0,0.0001)
    
    
    # Model variance-covariance
    Tau.B[1:K,1:K] ~ dwish(W[,], df)
    df <- K+1
    Sigma.B[1:K,1:K] <- inverse(Tau.B[,])
    for (k in 1:K){
    for (k.prime in 1:K){
    rho.B[k,k.prime] <- Sigma.B[k,k.prime]/
    sqrt(Sigma.B[k,k]*Sigma.B[k.prime,k.prime])
    }
    sigma.B[k] <- sqrt(Sigma.B[k,k])
    }
    
    }
    ",fill=TRUE)
sink()

# Number of sites
J <-length(unique(dat$DOW))


# Create identity matrix for Wishart dist'n
#!!!!!!!Number of parameters to estimate (K)
K <- 2

W <- diag(K)

as.numeric(unique(dat$DOW))

# Make sure Pool ID is numeric and goes from 1-10
dat2$G <- as.numeric(as.factor(as.numeric(dat2$DOW)))
# unique(as.numeric(as.factor(as.numeric(dat$Pool))))

# load data
data <- list(y = dat2$logYEP , x = dat2$Zdd, group = dat2$G, n = dim(dat)[1],
             J = J, W = W, K = K )


# Initial values
inits <- function (){
  list (BB=array(c(rep(rnorm(1,0,1),J),rep(rnorm(1,0,1),J)), c(J,K)), 
        mu.a=rnorm(1,0,1),mu.b=rnorm(1,0,1),
        sigma.y=runif(1,0,10), 
        Tau.B=rwish(K+1,diag(K))	 )
}


# Parameters monitored
params1 <- c("BB","mu.a","mu.b", "sigma.y","sigma.B","rho.B")


# MCMC settings
ni <- 9000
nt <- 3
nb <- 5000
nc <- 3


####### RUN JAGS ON SINGLE CORE #############

start.time = Sys.time()         # Set timer (2.3 mins)
# Call JAGS from R 

out1 <- jags(data, inits, params1, "model_test.txt", n.chains = nc, 
             n.thin = nt, n.iter = ni, n.burnin = nb)

end.time = Sys.time()
elapsed.time = round(difftime(end.time, start.time, units='mins'), dig = 2)
cat('Posterior computed in ', elapsed.time, ' minutes\n\n', sep='') 
# Calculate computation time


# Find which parameters, if any, have Rhat > 1.1
which(out1$BUGSoutput$summary[, c("Rhat")] > 1.1)

# See what max Rhat value is
# max(out1$BUGSoutput$summary[, c("Rhat")])


# Summarize posteriors
print(out1, dig = 3)

# outNested <- out1$BUGSoutput$summary
# write.csv(outNested, 'outNested.csv', row.names = T)
