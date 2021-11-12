## run mass-enhanced HMC for the sensor self-localization model
source('sensor_model.R')
lplot <- function(...) plot(..., type='l')
x.d <- 2*N
logtarget <- lpost
gradlt <- glpost

######################################################################################
## Use tempered Hamiltonian transitions to sample from the posterior for this model ##
######################################################################################
source('THT_sensor.R')
njumps <- 1
massInv <- 1
jsize <- 0.001
sin_period <- 4000 # number of mass levels times 2 (integers (in this code, even numbers) and half integers (resp. odd numbers))
maxLogScale <- 8
massScalingLvls <- exp(maxLogScale/2-maxLogScale/2*cos((0:(sin_period-1))*(2*pi)/sin_period))
#ksupport <- c(0:(round(sin_period*.1)-1), (sin_period/2-round(sin_period*.1)):(sin_period/2-1))
ksupport <- c(0:30, sin_period/2-1:30)
len_ksupport <- length(ksupport)
rk <- function() { sample(ksupport, 1) }
dk <- function(k) { (k %in% ksupport) }
acc_find <- round(len_ksupport/3)
spmax <- round(sin_period/2*1.1)
gamma <- 2
niter <- 1000
set.seed(124920)

## run parallel chains
library(parallel)
nparallel <- 12
numCores <- nparallel
run_chain_THT <- function(iii) {
    set.seed(92938+82745*iii)
    x.init <- runif(x.d, 0, L)
    out <- THT_sensor(x.init, logtarget, gradlt, niter, jsize, njumps,
        massInv, massScalingLvls, rk, dk, acc_find=acc_find, spmax=spmax, gamma=gamma)
    ## return value : list(xv=xvmat, na=na.vec, k=k.vec, k0=k0.vec, kacc=kacc.vec)
    return(out)
}
starttime <- Sys.time()
chains <- mclapply(1:nparallel, run_chain_THT, mc.cores=numCores)
endtime <- Sys.time()
duration <- endtime-starttime
beepr::beep(4)


## diagnose the run
library(tidyverse)
chain_data <- do.call(rbind, lapply(chains, function(arg) arg$xv))[,1:(2*N)] %>%
    data.frame %>% `colnames<-`(c(paste0('X',1:N),paste0('Y',1:N))) %>%
    mutate(chain=rep(1:nparallel, each=niter)) %>%
    mutate(iter=rep(1:niter, times=nparallel)) %>%
    pivot_longer(cols=!c(iter,chain), names_to=c('coord','sensor'), names_pattern='([XY])(\\d+)', values_to='value') %>%
    pivot_wider(names_from="coord", values_from="value")

p <- ggplot(data=chain_data) +
    facet_wrap(.~chain, nrow=4) +
    geom_point(aes(x=X, y=Y, col=sensor), size=0.01) +
    geom_text(data=data.frame(sensor.loc) %>% `colnames<-`(c('X','Y')) %>% mutate(sensor=factor(1:(N+Nknown))),
        mapping=aes(x=X, y=Y, label=as.character(sensor)), size=3) +
    guides(color = guide_legend(override.aes = list(size=3)))

p

plot(unlist(chain_data%>%filter(sensor==6)%>%select("Y")), pch='.'); abline(v=niter*1:nparallel, lty=2)

pairs(chain_data%>%filter(chain==1)%>%pivot_wider(id_cols=c(chain,iter),
    names_from=sensor, names_prefix='Y', values_from=Y)%>%select(Y1:Y8), pch='.',
    xlim=c(0,1), ylim=c(0,1))

par(mfrow=c(2,2))
chainid <- 9
U <- -apply(chains[[chainid]]$xv[,1:x.d], 1, logtarget)
H <- U + .5*sum(chains[[chainid]]$xv[,x.d+1:x.d]^2/massInv)
lplot(U, ylim=min(U)+c(0,30))
lplot(H, ylim=min(H)+c(0,30))
lplot(chains[[chainid]]$xv[,3+N])
lplot(chains[[chainid]]$xv[,6+N]); par(mfrow=c(1,1))

## save the run
#filename <- paste0('chains_','sensor',N,'known',Nknown,'niter',niter,'jsize',jsize,'maxLogScale',maxLogScale,'sinPeriod',sin_period,'accfind',acc_find,'spmax',spmax,'gamma',gamma,'nparallel',nparallel)
#save(N, Nknown, sensor.loc, logtarget, gradlt, niter, jsize, njumps, massInv, maxLogScale, sin_period, massScalingLvls, acc_find, spmax, ksupport, gamma, nparallel, chains, file=paste0('data/', filename, '.RData'))
## load data file
#load(paste0('data/',filename,'.RData'))

## save the plots ##
## marginal distributions
save_marginal_dists <- FALSE
if (save_marginal_dists) {
    filenamecore <- "chains_sensor8known3niter1000jsize0.001maxLogScale8sinPeriod4000accfind20spmax2200gamma2nparallel12"
    load(paste0("data/",filenamecore,".RData"))
    library(tidyverse)
    chain_data <- do.call(rbind, lapply(chains, function(arg) arg$xv))[,1:(2*N)] %>%
        data.frame %>% `colnames<-`(c(paste0('X',1:N),paste0('Y',1:N))) %>%
        mutate(chain=rep(1:nparallel, each=niter)) %>%
        mutate(iter=rep(1:niter, times=nparallel)) %>%
        pivot_longer(cols=!c(iter,chain), names_to=c('coord','sensor'), names_pattern='([XY])(\\d+)', values_to='value') %>%
        pivot_wider(names_from="coord", values_from="value")
    symmetry_line_slope <- (sensor.loc[11,2]-(sensor.loc[10,2]+sensor.loc[9,2])/2)/(sensor.loc[11,1]-(sensor.loc[10,1]+sensor.loc[9,1])/2) # the line connecting sensor 11 and the midpoint of sensors 9 and 10 (slope)
    symmetry_line_intercept <- sensor.loc[11,2] - sensor.loc[11,1]*symmetry_line_slope # intercept
    p <- ggplot(data=chain_data) +
        facet_wrap(.~chain, nrow=4) +
        geom_point(aes(x=X, y=Y, col=sensor), size=0.008) +
        geom_abline(slope=symmetry_line_slope, intercept=symmetry_line_intercept, color='red', linetype='dashed') +
        geom_text(data=data.frame(sensor.loc) %>% `colnames<-`(c('X','Y')) %>% mutate(sensor=factor(1:(N+Nknown))),
            mapping=aes(x=X, y=Y, label=as.character(sensor)), size=3) +
        theme_bw() + guides(color = guide_legend(override.aes = list(size=3)))
    ggsave(filename=paste0('figure/marginal_',filenamecore,'.png'), plot=p, width=7, height=9, units="in", dpi=200)
}

## comparison of tempered Hamiltonian transitions (=mass-modulated HMC) and HMC in terms of posterior marginal samples
save_marginal_comparison <- FALSE
if (save_marginal_comparison) {
    ## load tempered Hamiltonian transitions
    THTfilenamecore <- "chains_sensor8known3niter1000jsize0.001maxLogScale8sinPeriod4000accfind20spmax2200gamma2nparallel12"
    load(paste0("data/",THTfilenamecore,".RData"))
    library(tidyverse)
    THTchain_data <- do.call(rbind, lapply(chains, function(arg) arg$xv))[,1:(2*N)] %>%
        data.frame %>% `colnames<-`(c(paste0('X',1:N),paste0('Y',1:N))) %>%
        mutate(chain=rep(1:nparallel, each=niter)) %>%
        mutate(iter=rep(1:niter, times=nparallel)) %>%
        pivot_longer(cols=!c(iter,chain), names_to=c('coord','sensor'), names_pattern='([XY])(\\d+)', values_to='value') %>%
        pivot_wider(names_from="coord", values_from="value") %>%
        mutate(algo="THT")
    ## load HMC
    HMCfilenamecore <- "chains_sensor8known3niter1000jsize0.001maxLogScale0sinPeriod4000accfind20spmax2200gamma2nparallel12"
    load(paste0("data/",HMCfilenamecore,".RData"))
    library(tidyverse)
    HMCchain_data <- do.call(rbind, lapply(chains, function(arg) arg$xv))[,1:(2*N)] %>%
        data.frame %>% `colnames<-`(c(paste0('X',1:N),paste0('Y',1:N))) %>%
        mutate(chain=rep(1:nparallel, each=niter)) %>%
        mutate(iter=rep(1:niter, times=nparallel)) %>%
        pivot_longer(cols=!c(iter,chain), names_to=c('coord','sensor'), names_pattern='([XY])(\\d+)', values_to='value') %>%
        pivot_wider(names_from="coord", values_from="value") %>%
        mutate(algo="HMC")
    ## combine two data frames
    combined_chain_data <- rbind(THTchain_data, HMCchain_data)
    ## marginal posterior samples
    symmetry_line_slope <- (sensor.loc[11,2]-(sensor.loc[10,2]+sensor.loc[9,2])/2)/(sensor.loc[11,1]-(sensor.loc[10,1]+sensor.loc[9,1])/2) # the line connecting sensor 11 and the midpoint of sensors 9 and 10 (slope)
    symmetry_line_intercept <- sensor.loc[11,2] - sensor.loc[11,1]*symmetry_line_slope # intercept
    p <- ggplot(data=combined_chain_data%>%filter(chain<=6)%>%mutate(algo=factor(algo,levels=c("THT","HMC")))) +
        facet_wrap(~algo+chain, nrow=4) +
        geom_point(aes(x=X, y=Y, col=sensor), size=0.008) +
        geom_abline(slope=symmetry_line_slope, intercept=symmetry_line_intercept, color='red', linetype='dashed') +
        geom_text(data=data.frame(sensor.loc) %>% `colnames<-`(c('X','Y')) %>% mutate(sensor=factor(1:(N+Nknown))),
            mapping=aes(x=X, y=Y, label=as.character(sensor)), size=3) + theme_bw() +
        guides(color = guide_legend(override.aes = list(size=3)))
    ggsave(filename='figure/marginal_comparison_chains_sensor8known3niter1000jsize0.001maxLogScale0_or_8sinPeriod4000accfind20spmax2200gamma2nparallel12.png', plot=p, width=7, height=9, units="in", dpi=200)
}

## paired distributions
save_paired_dists <- FALSE
if (save_paired_dists) {
    load("data/chains_sensor8known3niter1000jsize0.001maxLogScale8sinPeriod4000accfind20spmax2200gamma2nparallel12.RData")
    png(filename='figure/paired_chain1_sensor8known3niter1000jsize0.001maxLogScale8sinPeriod4000accfind20spmax2200gamma2.png', width=8, height=8, units="in", res=200)
    chain_data <- do.call(rbind, lapply(chains, function(arg) arg$xv))[,1:(2*N)] %>%
        data.frame %>% `colnames<-`(c(paste0('X',1:N),paste0('Y',1:N))) %>%
        mutate(chain=rep(1:nparallel, each=niter)) %>%
        mutate(iter=rep(1:niter, times=nparallel)) %>%
        pivot_longer(cols=!c(iter,chain), names_to=c('coord','sensor'), names_pattern='([XY])(\\d+)', values_to='value') %>%
        pivot_wider(names_from="coord", values_from="value")
    pairs(chain_data%>%filter(chain==1)%>%pivot_wider(id_cols=c(chain,iter),
        names_from=sensor, names_prefix='Y', values_from=Y)%>%select(Y1:Y8), pch='.',
        xlim=c(0,1), ylim=c(0,1))
    dev.off()
}

## traceplots
save_traceplots <- FALSE
if (save_traceplots) {
    load("data/chains_sensor8known3niter1000jsize0.001maxLogScale8sinPeriod4000accfind20spmax2200gamma2nparallel12.RData")
    source('sensor_model.R')
    lplot <- function(...) plot(..., type='l')
    x.d <- 2*N
    png(filename='figure/traceplot_chain1_Y3Y6_sensor8known3niter1000jsize0.001maxLogScale8sinPeriod4000accfind20spmax2200gamma2.png', width=6, height=4, units="in", res=200)
    par(mfrow=c(2,2), mar=c(4,4,0,0), oma=c(0,0,2,0.5))
    chainid <- 1
    linewidth <- .5; axislabel <- .8
    U <- -apply(chains[[chainid]]$xv[,1:x.d], 1, logtarget)
    lplot(U, ylim=min(U)+c(0,30), bty='L', xlab='', ylab='', lwd=linewidth); mtext('U', side=2, line=2, cex=axislabel); mtext(LETTERS[1], adj=0, line=0.5)
    lplot(chains[[chainid]]$xv[,3+N], bty='L', xlab='', ylab='', lwd=linewidth); mtext(expression(Y[3]), side=2, line=2, cex=axislabel); abline(h=sensor.loc[3,2], lty=2); mtext(LETTERS[2], adj=0, line=0.5)
    lplot(chains[[chainid]]$xv[,6+N], bty='L', xlab='iteration', ylab='', lwd=linewidth); mtext(expression(Y[6]), side=2, line=2, cex=axislabel); abline(h=sensor.loc[6,2], lty=2); mtext(LETTERS[3], adj=0, line=0.5)
    lplot(chains[[chainid]]$xv[,8+N], bty='L', xlab='iteration', ylab='', lwd=linewidth); mtext(expression(Y[8]), side=2, line=2, cex=axislabel); abline(h=sensor.loc[8,2], lty=2); mtext(LETTERS[4], adj=0, line=0.5)
    dev.off()
}



## pilot runs: construct Hamiltonian trajectories
rm(list=ls()[ls()!="start.point"])
lplot <- function(...) plot(..., type='l')
source("../leapfrog.R")
power <- 2
## see an example trajectory (Tempered Hamiltonian transitions)
jsize <- 0.001
massInv <- 1
sin_period <- 6000
maxLogScale <- 8
massScalingLvls <- exp(maxLogScale/2-maxLogScale/2*cos((0:(sin_period-1))*(2*pi)/sin_period))
nMSFlvls <- length(massScalingLvls)/2

xold <- start.point# <- chains[[5]]$xv[niter, 1:(2*N)]
#k <- sample(0:(nMSFlvls-1), 1) # initial mass scaling factor level (the index)
k <- 0
vold <- rnorm(x.d, 0, sqrt(massInv/massScalingLvls[2*k+1]))
U <- -lpost(xold)
H <- U + .5*sum(vold*vold/massInv*massScalingLvls[2*k+1]) -.5*x.d*log(massScalingLvls[2*k+1])
i <- 1
simlen <- nMSFlvls; Utraj <- numeric(simlen+1); Utraj[1] <- U
Ktildetraj <- numeric(simlen+1); Ktildetraj[1] <- .5*sum(vold*vold/massInv*massScalingLvls[2*k+1])
Kenergytraj <- numeric(simlen+1); Kenergytraj[1] <- Ktildetraj[1] - -.5*x.d*log(massScalingLvls[2*k+1])
MSFtraj <- numeric(simlen+1)
xv <- matrix(NA, nrow=simlen+1, ncol=4*N); xv[1,] <- c(xold, vold)
n <- 1
repeat {
    n <- n+1; if (n > simlen+1) { break }
    xvnew <- lf(xold, vold, gd=glpost, jsize=jsize*(massScalingLvls[2*k+2])^(2/(power+2)),
        njumps=1, massInv=massInv/massScalingLvls[2*k+2])
    xnew <- xvnew[1:x.d]; vnew <- xvnew[x.d+1:x.d]
    for (comp in 1:(2*N)) {
        if (xnew[comp] < 0) { vnew[comp] <- abs(vnew[comp]) }
        if (xnew[comp] > L) { vnew[comp] <- -abs(vnew[comp]) }
    } # reflection at the boundaries
    xv[n,] <- c(xnew, vnew)
    Unew <- -lpost(xnew)
    xold <- xnew; vold <- vnew
    k <- k+1
    if (k>=nMSFlvls) { k <- k-nMSFlvls }
    Utraj[n] <- Unew; Ktildetraj[n] <- .5*sum(vnew*vnew/massInv*massScalingLvls[2*k+1])
    Kenergytraj[n] <- Ktildetraj[n] - .5*x.d*log(massScalingLvls[2*k+1])
    MSFtraj[n] <- massScalingLvls[2*k+1]
}
Henergytraj <- Kenergytraj + Utraj

sids <- 6:10
lplot(NA, xlim=c(0,1), ylim=c(0,1)); text(sensor.loc[sids,], col=sids, labels=sids, cex=2)
for (sid in sids) {
    lines(xv[,sid], xv[,N+sid], col=sid)
    text(xv[1,sid], xv[1,N+sid], col=sid, labels="S", cex=2)
    text(xv[simlen+1,sid], xv[simlen+1,N+sid], col=sid, labels="E", cex=2)
}

lplot(Utraj+Kenergytraj); abline(h=Henergytraj[1], lty=2); lines(Utraj, col='red', lty=2); lines(Ktildetraj, col='cyan', lty=2); lines(MSFtraj/2, col='green', lty=3)
(Henergytraj[simlen+1] - Henergytraj[1])

power <- 2
lplot(xv[,x.d+1]*MSFtraj^(1/(power+2))); points(xv[,x.d+1]*MSFtraj^(1/(power+2)), pch='+', col='green')
lplot(xv[,1]*MSFtraj^(-1/(power+2)))



## summary of results and implications for tuning:
## the trajectory for velocity typically looks oscillatory. A good scale to look at v*alpha^(1/(power+2)), where power is the degree of the polynomial that the potential function U(x) grows at. Note that the leapfrog step size is typically scaled as base_step_size * alpha^(2/(power+2)).  (You know alpha is the mass scaling factor).
## jsize and the length of massScalingLvls should be tuned as follows:
## jsize should be tuned so that there at least some number of (say >~ 10) leapfrog points per velocity oscillation cycle.
## As the length of massScalingLvls increases, the total length of the trajectory increases.  The length of massScalingLvls should be tuned so that there are a good number of velocity oscillation cycles (say around 20 or more) during one period of massScalingLvls. If there are too few velocity oscillation cycles, (say 1~2), the Hamiltonian during the second half of the trajectory does not decrease as much as the increase in Hamiltonian during the first half, leading to exponentially small acceptance probability.

