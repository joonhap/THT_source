source("../leapfrog.R")

## tempered Hamiltonian transitions
THT <- function(x.init, logtarget, gradlt, niter, jsize, njumps=1,
    massInv=1, massScalingLvls, rk, dk, k_support_maxMSF, acc_find, spmax, gamma=2) {
    ## x.init: initial condition
    ## logtarget: log target density function (R function)
    ## gradlt: gradient of log target density (R function)
    ## niter: number of iterations
    ## jsize: the size of leapfrog jumps
    ## njumps: the number of leapfrog jumps at each mass scaling factor level
    ## massInv: the inverse of the mass
    ## massScalingLvls: a vector of the levels of the mass scaling factor. The odd indexed (k=integer) values are the levels at which the velocity variable is drawn, and the even indexed (k=half-integer) values are the levels at which the Hamiltonian dynamics is simulated.
    ## rk: a function that randomly draws from the target distribution for k
    ## dk: a function that evaluates the probability mass function for k
    ## k_support_maxMSF: instead of specifying rk and dk, the distribution for k can be defined as uniform distribution on the integer k for which massScalingLvls[2*k+1] is below k_support_maxMSF. Ignored if rk and dk were supplied.
    ## acc_find: the number of acceptable proposals to find in each iteration
    ## spmax: max number of sequential proposals to try
    ## gamma: a time scale parameter. The leapfrog step size for each step is given by epsilon*exp(4*eta/(gamma+2)) where eta=log(massScalingFactor)/2.
    dimension <- length(x.init)
    nMSFlvls <- length(massScalingLvls)/2
    if (!(exists("rk") && exists("dk"))) {
        if (!exists("k_support_maxMSF")) {
            stop("Distribution on k is not specified.")
        } else {
            ksupport <- which(massScalingLvls[2*0:(nMSFlvls-1)+1]<k_support_maxMSF)-1
            len_ksupport <- length(ksupport)
            if (len_ksupport==0) { stop("k_support_maxMSF is too small. There is no k such that massScalingLvls[2*k+1]<k_supoprt_maxMSF") }
            rk <- function() { sample(ksupport, 1) }
            dk <- function(k) { (k %in% ksupport) }
        }
    }
    i <- 0; x <- x.init; U <- -logtarget(x)
    xvmat <- matrix(NA, niter+1, 2*dimension)
    xvmat[1,1:dimension] <- x
    na.vec <- numeric(niter) # number of acceptable candidates proposed
    k0.vec <- numeric(niter) # k at the start of iteration
    k.vec <- numeric(niter) # k at the end of iteration
    kacc.vec <- numeric(niter) # k for the next state (new k if accepted or starting k if rejected)
    repeat {
        Lambda <- runif(1); i <- i+1; if (i > niter) { break }
        ##cat(paste0("| (",i,") "))
        k <- rk() # initial mass scaling factor level (the index)
        k0.vec[i] <- k
        v <- rnorm(dimension, 0, sqrt(massInv/massScalingLvls[2*k+1]))
        H <- U + .5*sum(v*v/massInv)*massScalingLvls[2*k+1] -.5*dimension*log(massScalingLvls[2*k+1]) - log(dk(k))
        n <- 0; n_acc <- 0; xold <- x; vold <- v
        repeat {
            n <- n+1; if (n > spmax) { break }
            xvnew <- lf(xold, vold, gd=gradlt, jsize=jsize*(massScalingLvls[2*k+2])^(2/(gamma+2)),
                njumps=njumps, massInv=massInv/massScalingLvls[2*k+2])
            xnew <- xvnew[1:dimension]; vnew <- xvnew[dimension+1:dimension]
            xold <- xnew; vold <- vnew
            k <- k+1
            if (k>=nMSFlvls) { k <- k-nMSFlvls }
            if (is.finite(log(dk(k)))) { # if k is in the support of its target distribution
                Unew <- -logtarget(xnew)
                Hnew <- Unew+.5*sum(vnew*vnew/massInv)*massScalingLvls[2*k+1]-.5*dimension*log(massScalingLvls[2*k+1])
                if (is.nan(H) || is.nan(Unew)) { # numerical instability occurs
                    cat("i", i, "n", n, "\nx", x, "\nv", v, "\nxold", xold, "\nvold", vold,
                        "\nxnew", xnew, "\nvnew", vnew, "\nUnew", Unew, "H", H)  }
                if (Unew+.5*sum(vnew*vnew/massInv)*massScalingLvls[2*k+1]-.5*dimension*log(massScalingLvls[2*k+1])-log(dk(k)) - H < -log(Lambda)) {
                    n_acc <- n_acc + 1
                    cat(paste0("<-[",n_acc,"]\n"))
                }
                if (n_acc >= acc_find) {
                    x <- xnew; v <- vnew; U <- Unew; break
                }
            }
        }
        xvmat[i+1,] <- c(x,v)
        na.vec[i] <- n_acc
        k.vec[i] <- k
        kacc.vec[i] <- ifelse(n_acc >= acc_find, k, k0.vec[i])
    }
    return(list(xv=xvmat, na=na.vec, k=k.vec, k0=k0.vec, kacc=kacc.vec))
}

