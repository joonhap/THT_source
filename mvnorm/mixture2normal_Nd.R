## a probability distribution with a mixture of 2 normal distributions on an arbitrary dimension
#if (!exists(x.d)) {x.d <- 20} ## the dimension of the target density 
modeDist <- 400 ## the distance between two modes
comp.sd <- 1 # the sd of each mixture component of the target density
set.seed(12403)
direction <- {raw.direction <- rnorm(x.d); raw.direction/sqrt(sum(raw.direction^2))}
modes <- outer(c(-1,1), modeDist*direction/2)

target <- function(x, give_log=TRUE) { ## target denstiy
    if (x.d != length(x)) {
        stop("the length of x does not equal the specified dimension (x.d)") }
    lpdfs <- apply(modes, 1, function(mode) -sum((mode-x)^2)/(2*comp.sd^2)) -
        x.d/2*log(2*pi) - x.d*log(comp.sd)
    maxlpdf <- max(lpdfs)
    lmixpdf <- maxlpdf + log(sum(exp(lpdfs-maxlpdf))) - log(2)
    return(ifelse(give_log, lmixpdf, exp(lmixpdf)))
}

gd.target <- function(x) { ## gradient of the log of the target density
    if (x.d != length(x)) stop("the length of x does not equal the specificed dimension (x.d)")
    lpdfs <- apply(modes, 1, function(mode) -sum((mode-x)^2)/(2*comp.sd^2)) -
        x.d/2*log(2*pi) - x.d*log(comp.sd)
    return(apply(rbind(sapply(1:2, function(nm) -(x-modes[nm,])/comp.sd^2/
            sum(exp(lpdfs-lpdfs[nm])))), 1, sum))
}

## name the mode that a point is closest to
closest_mode <- function(x) {
    dists <- apply(modes, 1, function(v) sum((x-v)^2))
    return(which.min(dists))
}
