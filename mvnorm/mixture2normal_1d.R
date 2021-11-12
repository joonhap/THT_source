## a probability distribution with a mixture of two normal distributions on 1d
x.d <- 1 ## the dimension of the target density
modes <- 200*c(-1,1) # centers of the mixture components of the target density
o.comp.sd <- c(1,1) # the std.dev's of each mixture component of the target density (numeric if the two std.dev's are the same, a vector of length two if different)
if (length(o.comp.sd)==1) { o.comp.sd <- rep(o.comp.sd, 2) }
o.target <- function(x, give_log=TRUE) { ## target denstiy
    log_c1 <- dnorm(x, modes[1], o.comp.sd[1], log=TRUE)
    log_c2 <- dnorm(x, modes[2], o.comp.sd[2], log=TRUE)
    maxlog <- max(log_c1, log_c2)
    out <- maxlog + log(exp(log_c1 - maxlog) + exp(log_c2 - maxlog)) - log(2)
    return(ifelse(give_log, out, exp(out)))
}

gd.o.target <- function(x) { ## gradient of the log of the target density
    log_c1 <- dnorm(x, modes[1], o.comp.sd[1], log=TRUE)
    log_c2 <- dnorm(x, modes[2], o.comp.sd[2], log=TRUE)
    gd <- (-(x-modes[1])/o.comp.sd[1]^2)/(1+exp(log_c2-log_c1)) + (-(x-modes[2])/o.comp.sd[2]^2)/(1+exp(log_c1-log_c2))
    gd
}

