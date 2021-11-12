## mass-enhanced HMC
source("../leapfrog.R")

meHMC <- function(x.init, logtarget, gradlt, niter, jsize, njumps=1,
    massInv=1, massScaling=1, acc_find=1, spmax=1, randomize.eps=FALSE, prob) {
    ## x.init: initial condition
    ## logtarget: log target density function (R function)
    ## gradlt: gradient of log target density (R function)
    ## niter: number of iterations
    ## jsize: the size of leapfrog jumps
    ## njumps: the number of leapfrog jumps to obtain each successive proposal
    ## massInv: the inverse of the mass which which the velocity is drawn (vector of diagonal entries)
    ## massScaling: mass scaling factor for simulation of the Hamiltonian dynamics (can be supplied
    ## as a number or a function so that massScaling() outputs a draw of mass scaling factor.
    ## acc_find: the number of acceptable proposals to find in each iteration
    ## spmax: max number of sequential proposals to try
    if (hasArg(prob)) { nsettings <- length(prob) } else { nsettings <- 1 }
    isvec.jsize <- FALSE; isvec.massScaling <- FALSE
    isvec.acc_find <- FALSE; isvec.spmax <- FALSE # are there more than one options for jsize, massScaling, acc_find, or spmax?
    if (length(jsize)!=nsettings && length(jsize)>1) { stop("The length of jsize does not match that of prob.") }
    if (length(jsize)>1) { isvec.jsize <- TRUE }
    if (length(massScaling)!=nsettings && length(massScaling)>1) { stop("The length of massScaling does not match that of prob.") }
    if (length(massScaling)>1) { isvec.massScaling <- TRUE }
    if (length(acc_find)!=nsettings && length(acc_find)>1) { stop("The length of acc_find does not match that of prob.") }
    if (length(acc_find)>1) { isvec.acc_find <- TRUE }
    if (length(spmax)!=nsettings && length(spmax)>1) { stop("The length of spmax does not match that of prob.") }
    if (length(spmax)>1) { isvec.spmax <- TRUE }
    dimension <- length(x.init)
    m <- 0; x <- x.init; U <- -logtarget(x)
    xvmat <- matrix(NA, niter, 2*dimension)
    repeat {
        Lambda <- runif(1); m <- m+1; if (m > niter) { break }
        ## draw kernel settings from supplied options
        if (hasArg(prob)) { this.setting <- sample(1:nsettings, 1, prob=prob) }
        if (isvec.jsize) { this.jsize <- jsize[this.setting]
        } else { this.jsize <- jsize }
        if (isvec.massScaling) { this.massScaling <- massScaling[this.setting]
        } else { this.massScaling <- massScaling }
        if (isvec.acc_find) { this.acc_find <- acc_find[this.setting]
        } else { this.acc_find <- acc_find }
        if (isvec.spmax) { this.spmax <- spmax[this.setting]
        } else { this.spmax <- spmax }
        v <- rnorm(dimension, 0, sqrt(massInv)); H <- U + .5*sum(v*v/massInv)
        n <- 0; n_acc <- 0; xold <- x; vold <- v
        repeat {
            n <- n+1; if (n > this.spmax) { break }
            xvnew <- lf(xold, vold, gd=gradlt, jsize=this.jsize*ifelse(randomize.eps, runif(1,.8,1.2), 1),
                njumps=njumps, massInv=massInv/this.massScaling)
            xnew <- xvnew[1:dimension]; vnew <- xvnew[dimension+1:dimension]
            Unew <- -logtarget(xnew)
            if (is.nan(H) || is.nan(Unew)) {
                cat("m", m, "n", n, "\nx", x, "\nv", v, "\nxold", xold, "\nvold", vold,
                    "\nxnew", xnew, "\nvnew", vnew, "\nUnew", Unew, "H", H)
            }
            xold <- xnew; vold <- vnew
            Hinc <- Unew+.5*sum(vnew*vnew/massInv) - H
            if (Hinc < -log(Lambda)) {
                n_acc <- n_acc + 1
            }
            if (n_acc >= this.acc_find) {
                x <- xnew; v <- vnew; U <- Unew; break
            }
        }
        xvmat[m,] <- c(x,v)
    }
    return(xvmat)
}
