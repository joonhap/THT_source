## Tempered Hamiltonian transitions for the sensor localization model
## velocity reflection happens when the path reaches the boundary.
source("../leapfrog.R")

THT_sensor <- function(x.init, logtarget, gradlt, niter, jsize, njumps=1,
    massInv=1, massScalingLvls, rk, dk, k_support_maxMSF, acc_find, spmax, gamma=2, randomize.eps=FALSE) {
    ## mass-scaled HMC
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
    xvmat <- matrix(NA, niter, 2*dimension)
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
            xvnew <- lf(xold, vold, gd=gradlt, jsize=jsize*(massScalingLvls[2*k+2])^(2/(gamma+2))*ifelse(randomize.eps, runif(1,.8,1.2), 1),
                njumps=njumps, massInv=massInv/massScalingLvls[2*k+2])
            xnew <- xvnew[1:dimension]; vnew <- xvnew[dimension+1:dimension]
            for (comp in 1:(2*N)) {
                if (xnew[comp] < 0) { vnew[comp] <- abs(vnew[comp]) }
                if (xnew[comp] > L) { vnew[comp] <- -abs(vnew[comp]) }
            } # reflection at the boundaries
            Unew <- -logtarget(xnew)
            if (is.nan(H) || is.nan(Unew)) {
                cat("i", i, "n", n, "\nx", x, "\nv", v, "\nxold", xold, "\nvold", vold,
                    "\nxnew", xnew, "\nvnew", vnew, "\nUnew", Unew, "H", H)  }
            xold <- xnew; vold <- vnew
            k <- k+1
            if (k>=nMSFlvls) { k <- k-nMSFlvls }
            if (Unew+.5*sum(vnew*vnew/massInv)*massScalingLvls[2*k+1]-.5*dimension*log(massScalingLvls[2*k+1])-log(dk(k)) - H < -log(Lambda)) {
                n_acc <- n_acc + 1
            }
            if (n_acc >= acc_find) {
                x <- xnew; v <- vnew; U <- Unew; break
            }
        }
        xvmat[i,] <- c(x,v)
        na.vec[i] <- n_acc
        k.vec[i] <- k
        kacc.vec[i] <- ifelse(n_acc >= acc_find, k, k0.vec[i])
    }
    return(list(xv=xvmat, na=na.vec, k=k.vec, k0=k0.vec, kacc=kacc.vec))
}



################################################
### Gibbs sampler for the hierarchical model ###
################################################
THT_sensor_hierarchical <- function(x.init, lR2.init, lsignu2.init, niter, jsize, njumps=1,
    massInv=1, massScalingLvls, rk, dk, k_support_maxMSF, acc_find, spmax, gamma=2, lR2jsize, lsignu2jsize) {
    ## mass-scaled HMC
    ## x.init: initial condition
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
    i <- 0; x <- x.init; lR2 <- lR2.init; lsignu2 <- lsignu2.init
    U <- -lpost(x,lR2,lsignu2)
    xvmat <- matrix(NA, niter, 2*dimension)
    na.vec <- numeric(niter) # number of acceptable candidates proposed
    k0.vec <- numeric(niter) # k at the start of iteration
    k.vec <- numeric(niter) # k at the end of iteration
    kacc.vec <- numeric(niter) # k for the next state (new k if accepted or starting k if rejected)
    lR2.vec <- numeric(niter) # log(R^2)
    lsignu2.vec <- numeric(niter) # log(sigma_nu^2)
    repeat {
        i <- i+1; if (i > niter) { break }
        #####################################################
        ## draw the sensor locations given R^2 and signu^2 ##
        #####################################################
        Lambda <- runif(1)
        k <- rk() # initial mass scaling factor level (the index)
        k0.vec[i] <- k
        v <- rnorm(dimension, 0, sqrt(massInv/massScalingLvls[2*k+1]))
        H <- U + .5*sum(v*v/massInv)*massScalingLvls[2*k+1] -.5*dimension*log(massScalingLvls[2*k+1]) - log(dk(k))
        n <- 0; n_acc <- 0; xold <- x; vold <- v
        repeat {
            n <- n+1; if (n > spmax) { break }
            gradlt <- function(val) gllik(val, lR2, lsignu2)
            xvnew <- lf(xold, vold, gd=gradlt, jsize=jsize*(massScalingLvls[2*k+2])^(2/(gamma+2)),
                njumps=njumps, massInv=massInv/massScalingLvls[2*k+2])
            xnew <- xvnew[1:dimension]; vnew <- xvnew[dimension+1:dimension]
            for (comp in 1:(2*N)) {
                if (xnew[comp] < 0) { vnew[comp] <- abs(vnew[comp]) }
                if (xnew[comp] > L) { vnew[comp] <- -abs(vnew[comp]) }
            } # reflection at the boundaries
            Unew <- -lpost(xnew,lR2,lsignu2)
            if (is.nan(H) || is.nan(Unew)) {
                cat("i", i, "n", n, "\nx", x, "\nv", v, "\nxold", xold, "\nvold", vold,
                    "\nxnew", xnew, "\nvnew", vnew, "\nUnew", Unew, "H", H)  }
            xold <- xnew; vold <- vnew
            k <- k+1
            if (k>=nMSFlvls) { k <- k-nMSFlvls }
            if (Unew+.5*sum(vnew*vnew/massInv)*massScalingLvls[2*k+1]-.5*dimension*log(massScalingLvls[2*k+1])-log(dk(k)) - H < -log(Lambda)) {
                n_acc <- n_acc + 1
            }
            if (n_acc >= acc_find) {
                x <- xnew; v <- vnew; U <- Unew; break
            }
        }
        xvmat[i,] <- c(x,v)
        na.vec[i] <- n_acc
        k.vec[i] <- k
        kacc.vec[i] <- ifelse(n_acc >= acc_find, k, k0.vec[i])
        ########################################################
        ## draw R^2 and sigma_nu^2 given the sensor locations ##
        ########################################################
        ## update log(R^2) using HMC
        v_lR2 <- rnorm(1)
        gradlt_lR2 <- function(val) gllR2(val, x)
        lR2_new <- lf(lR2, v_lR2, gd=gradlt_lR2, jsize=lR2jsize, njumps=30)
        if (-llR2(lR2_new[1], x) + lR2_new[2]^2/2 < -llR2(lR2, x) + v_lR2^2/2 - log(runif(1))) {
            lR2 <- lR2_new[1]
        }
        lR2.vec[i] <- lR2
        ## update log(sigma_nu^2) using HMC
        v_lsignu2 <- rnorm(1)
        gradlt_lsignu2 <- function(val) gllsignu2(val, x)
        lsignu2_new <- lf(lsignu2, v_lsignu2, gd=gradlt_lsignu2, jsize=lsignu2jsize, njumps=30)
        if (-llsignu2(lsignu2_new[1], x) + lsignu2_new[2]^2/2 < -llsignu2(lsignu2, x) + v_lsignu2^2/2 - log(runif(1))) {
            lsignu2 <- lsignu2_new[1]
        }
        lsignu2.vec[i] <- lsignu2
    }
    return(list(xv=xvmat, na=na.vec, k=k.vec, k0=k0.vec, kacc=kacc.vec, lR2=lR2.vec, lsignu2=lsignu2.vec))
}

