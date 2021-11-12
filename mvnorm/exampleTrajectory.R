## simulates example trajectories with enhanced mass for a mixture of two high dimensional Gaussian distributions
setwd('~/Research/multimodalhmc/experiments/mvnorm/')
source('~/Research/multimodalhmc/experiments/leapfrog.R')
lplot <- function(...) plot(..., type='l')
library(tidyverse)

## bimodal target
logtarget <- function(x, give_log=TRUE) { ## target denstiy
    if (x.d != length(x)) {
        stop("the length of x does not equal the specified dimension (x.d)") }
    if (length(comp.sd)==1) { comp.sd <- rep(comp.sd, x.d) }
    lpdfs <- apply(modes, 1, function(mode) -sum((mode-x)^2/(2*comp.sd^2))) -
        x.d/2*log(2*pi) - sum(log(comp.sd))
    maxlpdf <- max(lpdfs)
    lmixpdf <- maxlpdf + log(sum(exp(lpdfs-maxlpdf))) - log(2)
    return(ifelse(give_log, lmixpdf, exp(lmixpdf)))
}

gradlt <- function(x) { ## gradient of the log of the target density
    if (x.d != length(x)) stop("the length of x does not equal the specificed dimension (x.d)")
    if (length(comp.sd)==1) { comp.sd <- rep(comp.sd, x.d) }
    lpdfs <- apply(modes, 1, function(mode) -sum((mode-x)^2/(2*comp.sd^2))) -
        x.d/2*log(2*pi) - sum(log(comp.sd))
    return(apply(rbind(sapply(1:2, function(nm) -(x-modes[nm,])/comp.sd^2/
            sum(exp(lpdfs-lpdfs[nm])))), 1, sum))
}

## unimodal target
unimodal <- TRUE
if (unimodal) {
    logtarget <- function(x, give_log=TRUE) { ## target density
        if (x.d != length(x)) {
            stop("the length of x does not equal the specified dimension (x.d)") }
        if (length(comp.sd)==1) { comp.sd <- rep(comp.sd, x.d) }
        lpdf <- -sum(x^2/(2*comp.sd^2)) - x.d/2*log(2*pi) - sum(log(comp.sd))
        return(ifelse(give_log, lpdf, exp(lpdf)))
    }
    gradlt <- function(x) { ## gradient of the log of the target density
        if (x.d != length(x)) stop("the length of x does not equal the specificed dimension (x.d)")
        if (length(comp.sd)==1) { comp.sd <- rep(comp.sd, x.d) }
        return(-x/comp.sd^2)
    }
}

## model
x.d <- 10
modeDist <- 20 ## the distance between two modes
set.seed(302385)
comp.sd <- runif(x.d, min=.5, max=4) # the sd of each mixture component of the target density
#comp.sd <- rep(1, x.d)
##set.seed(12403)
set.seed(Sys.time())
direction <- c(1, rep(0,x.d-1))
modes <- outer(c(-1,1), modeDist*direction/2)

## simulated path
massInv <- 1
massScaling <- 2000
simlen <- 5000
jsize <- 0.1
njumps <- 1
x.init <- modes[1,] + comp.sd*rnorm(x.d)
x.init <- comp.sd*rnorm(x.d)

dimension <- length(x.init)
m <- 0; x <- x.init
xvmat <- matrix(NA, simlen, 2*dimension)
v <- rnorm(dimension, 0, sqrt(massInv)); xvmat[1,] <- c(x,v)
n <- 1; n_acc <- 0; xold <- x; vold <- v
repeat {
    n <- n+1; if (n > simlen) { break }
    xvnew <- lf(xold, vold, gd=gradlt, jsize=jsize*sqrt(massScaling),
        njumps=njumps, massInv=massInv/massScaling)
    xvmat[n,] <- xvnew
    xnew <- xvnew[1:dimension]; vnew <- xvnew[dimension+1:dimension]
    Unew <- -logtarget(xnew)
    xold <- xnew; vold <- vnew
}
U <- -apply(xvmat[,1:x.d], 1, logtarget)
K <- apply(xvmat[,x.d+1:x.d], 1, function(v) .5*sum(v*v/massInv))
H <- U+K

##
save_pdf <- TRUE
if(save_pdf) {
    filename <- paste0('figures/Hamiltonian_exampleTrajectory_normal_','dim',x.d,'massScaling',massScaling,'_misspecified_M.pdf')
    pdf(filename, width=7, height=3)
    lplot(H, xlab='leapfrog step', ylab='Hamiltonian'); abline(h=H[1], col='orange')
    dev.off()
}

lplot(H, xlab='leapfrog steps', ylab='Hamiltonian'); abline(h=H[1], col='orange')

lplot(xvmat[,1]); points(which(H<=H[1]+log(2)), xvmat[H<=H[1]+log(2),1], pch='+', col='red', cex=2)

lplot(H, ylim=min(H)+c(0,200)); abline(h=H[1], lty=2, col='blue')

lplot(U+massScaling*K); abline(v=which(H<H[1]), col='red')
