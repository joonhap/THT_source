## this file simulates example trajectories like the ones constructed by tempered Hamiltonian transitions.
## this file was written to create some plots in the manuscript
lplot <- function(...) plot(..., type='l')
source("../leapfrog.R")
set.seed(28357)

run <- function(sgn, jsize, x.init, sin_period=2000) {
    njumps <- 1
    massInv <- 1
    maxLogScale <- 24
    massScalingLvls <- exp(sgn*maxLogScale/2-maxLogScale/2*cos((0:(sin_period-1))*(2*pi)/sin_period))
    acc_find <- 1
    nMSFlvls <- length(massScalingLvls)/2
    simlen <- nMSFlvls+1

    xv <- matrix(NA, simlen, 2*x.d)
    Utrace <- numeric(simlen)
    Ktrace <- numeric(simlen)
    MSFtrace <- numeric(simlen)
    Lambda <- runif(1)
    x <- x.init
    Utrace[1] <- U <- -logtarget(x)
    k <- 0#round(nMSFlvls*.4) # initial mass scaling factor
    v <- rnorm(x.d, 0, sqrt(massInv/massScalingLvls[2*k+1])) # this is actually v-tilde in the manuscript
    Ktrace[1] <- .5*sum(v*v/massInv)
    xv[1,] <- c(x,v)
    MSFtrace[1] <- massScalingLvls[2*k+1]
    H <- U + .5*sum(v*v/massInv)*massScalingLvls[2*k+1] - x.d/2*log(massScalingLvls[2*k+1])
    n <- 1; xold <- x; vold <- v
    repeat{
        n <- n+1
        xvnew <- lf(xold, vold, gd=gradlt, jsize=jsize*(massScalingLvls[2*k+2])^(2/(power+2)),
            njumps=1, massInv=massInv/massScalingLvls[2*k+2])
        xv[n,] <- xvnew
        xnew <- xvnew[1:x.d]; vnew <- xvnew[x.d+1:x.d]
        Utrace[n] <- Unew <- -logtarget(xnew)
        Ktrace[n] <- .5*sum(vnew*vnew/massInv)
        k <- k+1
        if (k >= nMSFlvls) { k <- k-nMSFlvls }
        MSFtrace[n] <- massScalingLvls[2*k+1]
        xold <- xnew; vold <- vnew
        if (n >= simlen) { break }
    }
    Ktildetrace <- Ktrace * MSFtrace
    Kenergytrace <- Ktildetrace - x.d/2*log(MSFtrace)-log(massInv) # negative log of extended velocity density (up to additive const.)
    Htrace <- Utrace + Ktrace
    Htildetrace <- Utrace + Ktildetrace
    Henergytrace <- Utrace + Kenergytrace
    return(Henergytrace[simlen]-Henergytrace[1])
}

##U(x) propto |x|^power
x.d <- 1
logtarget <- function(x) { -sum(x^2)^(power.truth/2) }
gradlt <- function(x) { -power.truth*sum(x^2)^(power.truth/2-1)*x }

# eta and eta+c (c=const.) is numerically the same if we change eps to eps*exp((-2a+1)c) = eps*exp((gamma-2)/(gamma+2)c)
power <- power.truth <- 2
negeta <- replicate(500, { x.init <- rnorm(1); run(-1, .1, x.init, 2000) });1
poseta <- replicate(500, { x.init <- rnorm(1); run(1, .1, x.init, 2000) });1
boxplot(negeta, poseta)
poly2_result <- list(neg=negeta, pos=poseta)
save(poly2_result, file="data/temperedTrajectory_Fig_DelHdist_Jun1_2021_poly2.RData")

power <- power.truth <- 1.2
negeta <- replicate(500, { x.init <- rnorm(1); run(-1, exp(3)*.1, x.init, 2000) });1 # init v * init jsize (propto) exp((2a-1)eta_0) = exp((2-power)/(2+power)eta_0)
poseta <- replicate(500, { x.init <- rnorm(1); run(1, .1, x.init, 2000) });1
boxplot(negeta, poseta)
poly1.2_result <- list(neg=negeta, pos=poseta)
save(poly1.2_result, file="data/temperedTrajectory_Fig_DelHdist_Jun1_2021_poly1.2.RData")

power <- power.truth <- 3
negeta <- replicate(500, { x.init <- rnorm(1); run(-1, exp(-2.4)*.1, x.init) });1
poseta <- replicate(500, { x.init <- rnorm(1); run(1, .1, x.init) });1
boxplot(negeta, poseta)
poly3_result <- list(neg=negeta, pos=poseta)
save(poly3_result, file="data/temperedTrajectory_Fig_DelHdist_Jun1_2021_poly3.RData")

power <- 3; power.truth <- 2
negeta <- replicate(500, { x.init <- rnorm(1); run(-1, exp(-2.4)*1/exp(4.8), x.init, sin_period=round(60000)) });1
poseta <- replicate(500, { x.init <- rnorm(1); run(1, 1/exp(4.8), x.init, sin_period=round(60000)) });1
boxplot(negeta, poseta)
poly_true2_est3_result <- list(neg=negeta, pos=poseta)
save(poly_true2_est3_result, file="data/temperedTrajectory_Fig_DelHdist_Jun1_2021_poly_trueDeg2_est3.RData")

power <- 2; power.truth <- 3
negeta <- replicate(500, { x.init <- rnorm(1); run(-1, 1/exp(4.8), x.init, 60000)});1
poseta <- replicate(500, { x.init <- rnorm(1); run(1, 1/exp(4.8), x.init, 60000) });1
boxplot(negeta, poseta)
poly_true3_est2_result <- list(neg=negeta, pos=poseta)
save(poly_true3_est2_result, file="data/temperedTrajectory_Fig_DelHdist_Jun1_2021_poly_trueDeg3_est2.RData")

power <- 1.2; power.truth <- 2
negeta <- replicate(100, { x.init <- rnorm(1); run(-1, exp(3)*exp(-3), x.init, sin_period=2000) });1
poseta <- replicate(100, { x.init <- rnorm(1); run(1, exp(-3), x.init, sin_period=2000) });1
boxplot(negeta, poseta)
poly_true2_est1.2_result <- list(neg=negeta, pos=poseta)
save(poly_true2_est1.2_result, file="data/temperedTrajectory_Fig_DelHdist_Jun1_2021_poly_trueDeg2_est1.2.RData")

power <- 2; power.truth <- 1.2
negeta <- replicate(100, { x.init <- rnorm(1); run(-1, exp(-3), x.init, sin_period=20000) });1
poseta <- replicate(100, { x.init <- rnorm(1); run(1, exp(-3), x.init, sin_period=20000) });1
boxplot(negeta, poseta)
poly_true1.2_est2_result <- list(neg=negeta, pos=poseta)
save(poly_true1.2_est2_result, file="data/temperedTrajectory_Fig_DelHdist_Jun1_2021_poly_trueDeg1.2_est2.RData")


##bimodal 1-dim normal
source("../mvnorm/mixture2normal_1d.R"); logtarget <- o.target; gradlt <- gd.o.target
power <- 2; x.d <- 1
negeta <- replicate(500, { x.init <- modes[1] + o.comp.sd[1]*rnorm(1); run(-1, .1, x.init) })
poseta <- replicate(500, { x.init <- modes[1] + o.comp.sd[1]*rnorm(1); run(1, .1, x.init) })
bimodal_result <- list(neg=negeta, pos=poseta)
save(bimodal_result, file="data/temperedTrajectory_Fig_DelHdist_Jun1_2021_bimodal.RData")

