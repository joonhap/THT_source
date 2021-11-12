## 2d sensor network model from Ihler et al. (2005) and Ahn et al. (2013)c Distributed and adaptive darting monte carlo through regenerations
setwd('~/Research/multimodalhmc/experiments/sensor/')
if (Sys.info()['sysname'] == 'darwin') { options(device='quartz') }
if (Sys.info()['sysname'] == 'Linux') { options(device='X11') }

## parameters
N <- 8 # number of sensors
Nknown <- 3  # number of sensors with known location
L <- 1 # side length of the square in which sensors are placed
rounds <- 1 # how many rounds of meausrements were made?
R2_t <- (.3)^2 # the true parameter value (R^2) determining the probability of inter-sensor distance measurement
signu2_t <- (.02)^2 # the true additive noise variance (sigma_nu^2) for distance measurements

incseq <- function(a, b) {
    if (b >= a) return(a:b)
    if (b < a) return(NULL)
}

## data generation
set.seed(4838+100*N+Nknown+28237*(rounds-1))
sensor.loc <- matrix(runif((N+Nknown)*2, 0, L), ncol=2)
measurements <- array(NA, dim=c(N, N+Nknown, rounds))
true_distances <- matrix(0, N, N+Nknown)
for (i in 1:N) {
    for (j in incseq(i+1,N+Nknown)) {
        distance <- sqrt(sum((sensor.loc[i,] - sensor.loc[j,])^2))
        true_distances[i,j] <- distance; if (j <= N) { true_distances[j,i] <- distance }
        observed <- ifelse(runif(rounds) < exp(-.5*distance^2/R2_t), 1, 0)
        for (r in 1:rounds) {
            if (observed[r]) {
                measurements[i,j,r] <- distance + rnorm(1,0,sqrt(signu2_t))
                if (j <= N) { measurements[j,i,r] <- measurements[i,j,r] }
            }
        }
    }
}

## pairwise distances for a given sensor configuration
EstDistances <- function(xy) {
    ## EstSensorLoc: estimated sensor locations with respect to which the pairwise distances are to be computed
    EstSensorLoc <- matrix(xy, ncol=2)
    returnval <- matrix(0, N, N+Nknown)
    for (i in 1:N) {
        for (j in incseq(i+1,N)) {
            returnval[i,j] <- returnval[j,i] <- sqrt(sum((EstSensorLoc[i,]-EstSensorLoc[j,])^2))
        }
        for (j in incseq(N+1,N+Nknown)) {
            returnval[i,j] <- sqrt(sum((EstSensorLoc[i,]-sensor.loc[j,])^2))
        }
    }
    return(returnval)
}


## likelihood and prior models
llik <- function(xy, lR2=log(R2_t), lsignu2=log(signu2_t)) { # measurement loglikelihood
    ## x is the distance of N sensors of which the locations are unknown
    ## x is a vector of length 2N, where the first N numbers represent the x-coordinates
    ## lR2: log(R^2), lsignu2: log(sigma_nu^2)
    loglik <- 0
    xy.mat <- matrix(xy, ncol=2)
    for (i in 1:N) {
        for (j in incseq(i+1,N+Nknown)) {
            if (j <= N) { jsensor.loc <- xy.mat[j,] }
            if (j > N) { jsensor.loc <- sensor.loc[j,] }
            distance <- sqrt(sum((xy.mat[i,] - jsensor.loc)^2))
            loglikinc <- ifelse(is.na(measurements[i,j,]), log(1 - exp(-.5*distance^2*exp(-lR2))),
                -.5*distance^2*exp(-lR2) + dnorm(measurements[i,j,], distance, exp(.5*lsignu2), log=TRUE))
            loglik <- loglik + sum(loglikinc)
        }
    }
    return(loglik)
} 

lprior <- function(xy) { # log-density of uniform prior inside [0,L]X[0,L]
    if (all(xy>0 & xy<L)) {
        return(0)
    } else {
        return(-Inf)
    }
}

lpost <- function(xy, lR2, lsignu2) { llik(xy, lR2, lsignu2) + lprior(xy) }

gllik <- function(xy, lR2=log(R2_t), lsignu2=log(signu2_t)) { # gradient of measurement loglikelihood
    if (length(xy) != 2*N) { stop("The length of the input to gllik is not equal to 2N") }
    if (length(lR2)>1 || length(lsignu2)>1) { stop("lR2 and lsignu2 should not be a vector.") }
    xy.mat <- matrix(xy, ncol=2)
    gloglik.mat <- matrix(0, N, 2)
    for (i in 1:N) {
        for (j in incseq(i+1,N+Nknown)) {
            if (j <= N) { jsensor.loc <- xy.mat[j,] }
            if (j > N) { jsensor.loc <- sensor.loc[j,] }
            distance <- sqrt(sum((xy.mat[i,] - jsensor.loc)^2))
            nNAmeas <- sum(is.na(measurements[i,j,]))
            nmeas <- rounds - nNAmeas
            NAmeas_gl <- nNAmeas * (xy.mat[i,]-jsensor.loc)*exp(-lR2) * exp(-.5*distance^2*exp(-lR2))/(1 - exp(-.5*distance^2*exp(-lR2))) # gradient of log likelihood for unobserved data
            meas_gl <- nmeas * (-xy.mat[i,]+jsensor.loc)*exp(-lR2) - sum(distance-measurements[i,j,], na.rm=TRUE)/distance*exp(-lsignu2)*(xy.mat[i,]-jsensor.loc) # gradient of log likelihood for observed data
            gloglik.mat[i,] <- gloglik.mat[i,] + NAmeas_gl + meas_gl
            if (j <= N) { gloglik.mat[j,] <- gloglik.mat[j,] - NAmeas_gl - meas_gl }
        }
    }
    return(c(gloglik.mat))
}

glprior <- function(xy) { # gradient of lprior
    buffer_width <- 0.001*L # the virtual "thickness" of the wall
    ## for each component, the gradient is proportional to the amount the boundary is violated,
    ## but capped at 10.
    #### NOTE: this prior has been disabled by multiplying by zero. Instead, the velocity is reflected if the path hits the boundary.
    return(0*(ifelse(xy < 0, pmin(abs(xy)^2/buffer_width^2, Inf), 0) + ifelse(xy > L, pmax(-(xy-L)^2/buffer_width^2, -Inf), 0)))
}

glpost <- function(xy) { gllik(xy) + glprior(xy) }



#############################################################################
## distribution of R2 and signu2 conditional on xy (for the Gibbs sampler) ##
#############################################################################
## the prior for R is given by the exponential distribution with scale 0.5
Rprior_rate <- 1/.5
## the prior for signu is given by the exponential distribution with scale 0.05
signuprior_rate <- 1/.05

llR2 <- function(lR2, xy) {
    ## log density for log(R^2) given the sensor locations
    ## lR2: log(R^2), xy: sensor locations
    returnval <- -Rprior_rate*exp(lR2/2) + lR2/2
    xy.mat <- matrix(xy, ncol=2)
    for (i in 1:N) {
        for (j in incseq(i+1,N+Nknown)) {
            if (j <= N) { jsensor.loc <- xy.mat[j,] }
            if (j > N) { jsensor.loc <- sensor.loc[j,] }
            distance <- sqrt(sum((xy.mat[i,] - jsensor.loc)^2))
            nNAmeas <- sum(is.na(measurements[i,j,]))
            returnval <- returnval + nNAmeas * log(1 - exp(-.5*distance^2*exp(-lR2)))
            returnval <- returnval + (rounds-nNAmeas) * -.5*distance^2*exp(-lR2)
        }
    }
    return(returnval)
}

llsignu2 <- function(lsignu2, xy) {
    ## log density for log(signu^2) given the sensor locations
    ## lsignu2: log(signu^2), xy: sensor locations
    if (length(lsignu2)>1) { stop("lsignu2 should not be a vector.") }
    returnval <- -signuprior_rate*exp(lsignu2/2) + lsignu2/2
    xy.mat <- matrix(xy, ncol=2)
    for (i in 1:N) {
        for (j in incseq(i+1,N+Nknown)) {
            if (j <= N) { jsensor.loc <- xy.mat[j,] }
            if (j > N) { jsensor.loc <- sensor.loc[j,] }
            distance <- sqrt(sum((xy.mat[i,] - jsensor.loc)^2))
            returnval <- returnval + sum(dnorm(measurements[i,j,], mean=distance, sd=exp(.5*lsignu2), log=TRUE), na.rm=TRUE)
        }
    }
    return(returnval)
}

gllR2 <- function(lR2, xy) {
    ## gradient of log density for log(R^2) given the sensor locations
    ## lR2: log(R^2), xy: sensor locations
    if (length(xy) != 2*N) { stop("The length of the input to gllik is not equal to 2N") }
    xy.mat <- matrix(xy, ncol=2)
    returnval <- -.5*Rprior_rate*exp(lR2/2) + .5
    for (i in 1:N) {
        for (j in incseq(i+1,N+Nknown)) {
            if (j <= N) { jsensor.loc <- xy.mat[j,] }
            if (j > N) { jsensor.loc <- sensor.loc[j,] }
            distance <- sqrt(sum((xy.mat[i,] - jsensor.loc)^2))
            nNAmeas <- sum(is.na(measurements[i,j,]))
            tempval <- .5*distance^2*exp(-lR2)
            returnval <- returnval + nNAmeas * tempval*(-exp(-tempval))/(1-exp(-tempval)) + (rounds-nNAmeas) * .5*distance^2*exp(-lR2)
        }
    }
    return(returnval)
}

gllsignu2 <- function(lsignu2, xy) {
    ## gradient of log density for log(signu^2) given the sensor locations
    ## lsignu2: log(signu^2), xy: sensor locations
    if (length(xy) != 2*N) { stop("The length of the input to gllik is not equal to 2N") }
    if (length(lsignu2) > 1) { stop("lsignu2 should not be a vector.") }
    xy.mat <- matrix(xy, ncol=2)
    returnval <- -.5*signuprior_rate*exp(lsignu2/2) + .5
    for (i in 1:N) {
        for (j in incseq(i+1,N+Nknown)) {
            if (j <= N) { jsensor.loc <- xy.mat[j,] }
            if (j > N) { jsensor.loc <- sensor.loc[j,] }
            distance <- sqrt(sum((xy.mat[i,] - jsensor.loc)^2))
            returnval <- returnval - .5*(rounds-sum(is.na(measurements[i,j,]))) + .5*sum((distance-measurements[i,j,])^2, na.rm=TRUE)*exp(-lsignu2)
        }
    }
    return(returnval)
}
