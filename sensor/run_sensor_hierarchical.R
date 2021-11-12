## run mass-enhanced HMC for the sensor self-localization model
source('sensor_model.R')
lplot <- function(...) plot(..., type='l')
x.d <- 2*N
logtarget <- lpost
gradlt <- glpost

######################################
## Tempered Hamiltonian transitions ##
######################################
source('THT_sensor.R')
njumps <- 1
massInv <- 1
jsize <- 0.001
sin_period <- 4000 # number of mass levels times 2 (integers (in this code, even numbers) and half integers (resp. odd numbers))
maxLogScale <- 8
massScalingLvls <- exp(maxLogScale/2-maxLogScale/2*cos((0:(sin_period-1))*(2*pi)/sin_period))
##ksupport <- c(0:(round(sin_period*.1)-1), (sin_period/2-round(sin_period*.1)):(sin_period/2-1))
ksupport <- c(0:30, sin_period/2-1:30)
len_ksupport <- length(ksupport)
rk <- function() { sample(ksupport, 1) }
dk <- function(k) { (k %in% ksupport) }
acc_find <- round(len_ksupport/3)
spmax <- round(sin_period/2*1.1)
gamma <- 2
lR2jsize <- .02
lsignu2jsize <- .02
lR2.init <- log(.5^2)
lsignu2.init <- log(.05^2)
niter <- 1000
set.seed(124920)

## run parallel chains
library(parallel)
nparallel <- 12
numCores <- nparallel
run_chain_THT <- function(iii) {
    set.seed(92938+82745*iii)
    x.init <- runif(x.d, 0, L)
    out <- THT_sensor_hierarchical(x.init, lR2.init, lsignu2.init, niter, jsize, njumps,
        massInv, massScalingLvls, rk, dk, acc_find=acc_find, spmax=spmax, gamma=gamma, lR2jsize=lR2jsize, lsignu2jsize=lsignu2jsize)
    ## return value : list(xv=xvmat, na=na.vec, k=k.vec, k0=k0.vec, kacc=kacc.vec, lR2.vec, lsignu2.vec)
    return(out)
}
starttime <- Sys.time()
timeD <- system.time(chains <- mclapply(1:nparallel, run_chain_THT, mc.cores=numCores))
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

par(mfrow=c(3,2))
chainid <- 1
U <- -sapply(1:niter, function(iii) { logtarget(chains[[chainid]]$xv[iii,1:x.d], chains[[chainid]]$lR2[iii], chains[[chainid]]$lsignu2[iii]) })
logden <- -U + sapply(1:niter, function(iii) {
    llR2(chains[[chainid]]$lR2[iii], chains[[chainid]]$xv[iii,1:x.d]) +
    llsignu2(chains[[chainid]]$lsignu2[iii], chains[[chainid]]$xv[iii,1:x.d]) } )
lplot(U)
lplot(logden); abline(v=FirstDomModeIter[chainid], col='red')
lplot(chains[[chainid]]$xv[,3+N]); abline(v=FirstDomModeIter[chainid], col='red')
lplot(chains[[chainid]]$xv[,6+N])
lplot(exp(chains[[chainid]]$lR2/2))
lplot(exp(chains[[chainid]]$lsignu2/2), log='y'); par(mfrow=c(1,1))

lplot(exp(chains[[chainid]]$lsignu2/2)[400:1000])

## save the run
filename <- paste0('chains_hierarchical_','sensor',N,'known',Nknown,'rounds',rounds,'niter',niter,'jsize',jsize,'maxLogScale',maxLogScale,'sinPeriod',sin_period,'accfind',acc_find,'spmax',spmax,'gamma',gamma,'lR2jsize',lR2jsize,'lsignu2jsize',lsignu2jsize,'nparallel',nparallel)
save(N, Nknown, rounds, sensor.loc, logtarget, gradlt, niter, jsize, njumps, massInv, maxLogScale, sin_period, massScalingLvls, acc_find, spmax, ksupport, gamma, lR2jsize, lsignu2jsize, nparallel, chains, duration, timeD, file=paste0('data/', filename, '.RData'))
## load data file
#load(paste0('data/',filename,'.RData'))

## save the plots ##
## marginal distributions
save_marginal_dists <- FALSE
if (save_marginal_dists) {
    filenamecore <- "chains_hierarchical_sensor8known3rounds1niter1000jsize0.001maxLogScale8sinPeriod4000accfind20spmax2200gamma2lR2jsize0.02lsignu2jsize0.02nparallel12"
    load(paste0("data/",filenamecore,".RData"))
    library(tidyverse)
    chain_data <- do.call(rbind, lapply(chains, function(arg) arg$xv))[,1:(2*N)] %>%
        data.frame %>% `colnames<-`(c(paste0('X',1:N),paste0('Y',1:N))) %>%
        mutate(chain=rep(1:nparallel, each=niter)) %>%
        mutate(iter=rep(1:niter, times=nparallel)) %>%
        pivot_longer(cols=!c(iter,chain), names_to=c('coord','sensor'), names_pattern='([XY])(\\d+)', values_to='value') %>%
        pivot_wider(names_from="coord", values_from="value")
    p <- ggplot(data=chain_data) +
        facet_wrap(.~chain, nrow=4) +
        geom_point(aes(x=X, y=Y, col=sensor), size=0.008) +
        geom_text(data=data.frame(sensor.loc) %>% `colnames<-`(c('X','Y')) %>% mutate(sensor=factor(1:(N+Nknown))),
            mapping=aes(x=X, y=Y, label=as.character(sensor)), size=3) + theme_bw() +
        guides(color = guide_legend(override.aes = list(size=3)))
    ggsave(filename=paste0('figure/marginal_',filenamecore,'.png'), plot=p, width=8, height=8, units="in", dpi=200)
}

## paired distributions
save_paired_dists <- FALSE
if (save_paired_dists) {
    filenamecore <- "sensor8known3rounds1niter1000jsize0.001maxLogScale8sinPeriod4000accfind20spmax2200gamma2lR2jsize0.02lsignu2jsize0.005nparallel12"
    load(paste0("data/chains_hierarchical_",filenamecore,".RData"))
    png(filename=paste0('figure/paired_hierarchical_chain1_',filenamecore,'.png'), width=8, height=8, units="in", res=200)
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
    filenamecore <- "sensor8known3rounds1niter1000jsize0.001maxLogScale8sinPeriod4000accfind20spmax2200gamma2lR2jsize0.02lsignu2jsize0.02nparallel12"
    load(paste0("data/chains_hierarchical_",filenamecore,".RData"))
    source('~/Research/multimodalhmc/experiments/sensor/sensor_model.R')
    lplot <- function(...) plot(..., type='l')
    x.d <- 2*N
    png(filename=paste0('figure/traceplot_hierarchical_chain1_Y3Y6_',filenamecore,'.png'), width=6, height=6, units="in", res=200)
    par(mfrow=c(3,2), mar=c(4,4,0,0), oma=c(0,0,2,.5))
    chainid <- 1
    U <- -sapply(1:niter, function(iii) { logtarget(chains[[chainid]]$xv[iii,1:x.d], chains[[chainid]]$lR2[iii], chains[[chainid]]$lsignu2[iii]) })
    logden <- -U + sapply(1:niter, function(iii) {
        llR2(chains[[chainid]]$lR2[iii], chains[[chainid]]$xv[iii,1:x.d]) +
            llsignu2(chains[[chainid]]$lsignu2[iii], chains[[chainid]]$xv[iii,1:x.d]) } )
    linewidth=.5; axislabel=.8
    lplot(U, bty='L', lwd=linewidth, xlab='', ylab=''); mtext(side=2, line=2, text='U', cex=axislabel); mtext(text=LETTERS[1], adj=0, line=0.5)
    lplot(logden, bty='L', lwd=linewidth, xlab='', ylab=''); mtext(side=2, line=2, text='log posterior density', cex=axislabel); mtext(text=LETTERS[2], adj=0, line=0.5)
    lplot(exp(chains[[chainid]]$lR2/2), bty='L', lwd=linewidth, xlab='', ylab=''); mtext(side=2, line=2, text=expression(R), cex=axislabel); abline(h=sqrt(R2_t), lty=2); mtext(text=LETTERS[3], adj=0, line=0.5)
    lplot(exp(chains[[chainid]]$lsignu2/2), log='y', bty='L', lwd=linewidth, xlab='', ylab=''); mtext(side=2, line=2, text=expression(sigma[e]), cex=axislabel); abline(h=sqrt(signu2_t), lty=2); mtext(text=LETTERS[4], adj=0, line=0.5)
    lplot(chains[[chainid]]$xv[,3+N], bty='L', lwd=linewidth, xlab='iteration', ylab=''); mtext(side=2, line=2, text=expression(Y[3]), cex=axislabel); abline(h=sensor.loc[3,2], lty=2); mtext(text=LETTERS[5], adj=0, line=0.5)
    lplot(chains[[chainid]]$xv[,6+N], bty='L', lwd=linewidth, xlab='iteration', ylab=''); mtext(side=2, line=2, expression(Y[6]), cex=axislabel); abline(h=sensor.loc[6,2], lty=2); mtext(text=LETTERS[6], adj=0, line=0.5)
    dev.off()
}



## Rhat statistic ###
filenamecore <- "sensor8known3rounds1niter1000jsize0.001maxLogScale8sinPeriod4000accfind20spmax2200gamma2lR2jsize0.02lsignu2jsize0.02nparallel12"
load(paste0("data/chains_hierarchical_",filenamecore,".RData"))
lplot <- function(...) plot(..., type='l')
x.d <- 2*N
library(rstan)
outvar <- array(
    c(sapply(1:x.d, function(jj) { sapply(1:nparallel, function(val) { chains[[val]]$xv[,jj] }) }), # sensor locations (2*N)
    sapply(1:nparallel, function(val) { exp(chains[[val]]$lR2/2) }), # parameter R
    sapply(1:nparallel, function(val) { exp(chains[[val]]$lsignu2/2) })), # parameter sigma_nu
    dim=c(niter, nparallel, x.d+2))
burnin <- 200
Rhatvals <- apply(outvar[(burnin+1):niter,,], 3, Rhat)

## time to a dominant mode
Us <- -sapply(1:nparallel, function(chainid) {
    sapply(1:niter, function(iii) { logtarget(chains[[chainid]]$xv[iii,1:x.d], chains[[chainid]]$lR2[iii], chains[[chainid]]$lsignu2[iii]) })
    })
logdens <- -Us + sapply(1:nparallel, function(chainid) {
        sapply(1:niter, function(iii) { llR2(chains[[chainid]]$lR2[iii], chains[[chainid]]$xv[iii,1:x.d]) +
                llsignu2(chains[[chainid]]$lsignu2[iii], chains[[chainid]]$xv[iii,1:x.d]) }) })
DomModeCriterion <- mean(logdens[501:1000,]) - 2*sd(logdens[501:1000,])/sqrt(20) # a chain will be considered to have reached one of the two dominant modes if the average log density of the next twenty iterations have exceeded the DomModeCriterion
wFirstDomModeIter <- apply(logdens, 2,
    function(vec) { csum <- cumsum(vec);
        which.max((csum[21:niter]-csum[1:(niter-20)])/20 > DomModeCriterion) })
c(mean(FirstDomModeIter), sd(FirstDomModeIter))/niter * duration
modehoppings <- sapply(1:nparallel,
    function(chainid) {
        Y3 <- chains[[chainid]]$xv[,3+N]
        sum((Y3[(FirstDomModeIter[chainid]+1):niter]-.4)*(Y3[FirstDomModeIter[chainid]:(niter-1)]-.4)<0) })






#####################################
## run standard HMC for this model ##
#####################################
x.d <- 2*N
logtarget <- lpost
gradlt <- glpost
njumps <- 30
massInv <- 1
jsize <- 0.001
sin_period <- 2 # number of mass levels times 2 (integers (in this code, even numbers) and half integers (resp. odd numbers))
maxLogScale <- 0
massScalingLvls <- exp(maxLogScale/2-maxLogScale/2*cos((0:(sin_period-1))*(2*pi)/sin_period))
#ksupport <- c(0:(round(sin_period*.1)-1), (sin_period/2-round(sin_period*.1)):(sin_period/2-1))
ksupport <- 0
len_ksupport <- length(ksupport)
rk <- function() { sample(ksupport, 1) }
dk <- function(k) { (k %in% ksupport) }
acc_find <- 1
spmax <- 1
gamma <- 2
lR2jsize <- .02
lsignu2jsize <- .02
lR2.init <- log(.5^2)
lsignu2.init <- log(.05^2)
niter <- 35000
set.seed(124920)

## run parallel chains
library(parallel)
nparallel <- 12
numCores <- nparallel
run_chain_THT <- function(iii) {
    set.seed(92938+82745*iii)
    x.init <- runif(x.d, 0, L)
    out <- THT_sensor_hierarchical(x.init, lR2.init, lsignu2.init, niter, jsize, njumps,
        massInv, massScalingLvls, rk, dk, acc_find=acc_find, spmax=spmax, gamma=gamma, lR2jsize=lR2jsize, lsignu2jsize=lsignu2jsize)
    ## return value : list(xv=xvmat, na=na.vec, k=k.vec, k0=k0.vec, kacc=kacc.vec, lR2.vec, lsignu2.vec)
    return(out)
}
starttime <- Sys.time()
timeD <- system.time(chains <- mclapply(1:nparallel, run_chain_THT, mc.cores=numCores))
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

par(mfrow=c(3,2))
chainid <- 12
lplot(Us[,chainid], ylim=min(Us[,chainid])+c(0,100))
lplot(logdens[,chainid], ylim=max(logdens[,chainid])+c(-100,0))
lplot(chains[[chainid]]$xv[,3+N])
lplot(chains[[chainid]]$xv[,6+N])
lplot(exp(chains[[chainid]]$lR2/2))
lplot(exp(chains[[chainid]]$lsignu2/2), log='y'); par(mfrow=c(1,1))


lplot(exp(chains[[chainid]]$lsignu2/2)[400:1000])

## save the run
filename <- paste0('chains_hierarchical_','sensor',N,'known',Nknown,'rounds',rounds,'niter',niter,'jsize',jsize,'maxLogScale',maxLogScale,'sinPeriod',sin_period,'accfind',acc_find,'spmax',spmax,'gamma',gamma,'lR2jsize',lR2jsize,'lsignu2jsize',lsignu2jsize,'nparallel',nparallel)
save(N, Nknown, rounds, sensor.loc, logtarget, gradlt, niter, jsize, njumps, massInv, maxLogScale, sin_period, massScalingLvls, acc_find, spmax, ksupport, gamma, lR2jsize, lsignu2jsize, nparallel, chains, duration, timeD, file=paste0('data/', filename, '.RData'))
## load data file
#load(paste0('data/',filename,'.RData'))

## save the plots ##
## marginal distributions
save_marginal_dists <- 1
if (save_marginal_dists) {
    filenamecore <- "sensor8known3rounds1niter35000jsize0.001maxLogScale0sinPeriod2accfind1spmax1gamma2lR2jsize0.02lsignu2jsize0.02nparallel12"
    load(paste0("data/chains_hierarchical_",filenamecore,".RData"))
    library(tidyverse)
    chain_data <- do.call(rbind, lapply(chains, function(arg) arg$xv))[,1:(2*N)] %>%
        data.frame %>% `colnames<-`(c(paste0('X',1:N),paste0('Y',1:N))) %>%
        mutate(chain=rep(1:nparallel, each=niter)) %>%
        mutate(iter=rep(1:niter, times=nparallel)) %>%
        pivot_longer(cols=!c(iter,chain), names_to=c('coord','sensor'), names_pattern='([XY])(\\d+)', values_to='value') %>%
        pivot_wider(names_from="coord", values_from="value")
    p <- ggplot(data=chain_data) +
        facet_wrap(.~chain, nrow=4) +
        geom_point(aes(x=X, y=Y, col=sensor), size=0.008) +
        geom_text(data=data.frame(sensor.loc) %>% `colnames<-`(c('X','Y')) %>% mutate(sensor=factor(1:(N+Nknown))),
            mapping=aes(x=X, y=Y, label=as.character(sensor)), size=3) + theme_bw() +
        guides(color = guide_legend(override.aes = list(size=3)))
    ggsave(filename=paste0('figure/marginal_',filenamecore,'.png'), plot=p, width=8, height=8, units="in", dpi=200)
}

## paired distributions
save_paired_dists <- 0
if (save_paired_dists) {
    filenamecore <- "sensor8known3rounds1niter35000jsize0.001maxLogScale0sinPeriod2accfind1spmax1gamma2lR2jsize0.02lsignu2jsize0.02nparallel12"
    load(paste0("data/chains_hierarchical_",filenamecore,".RData"))
    png(filename=paste0('figure/paired_hierarchical_chain1_',filenamecore,'.png'), width=8, height=8, units="in", res=200)
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
save_traceplots <- 1
if (save_traceplots) {
    filenamecore <- "sensor8known3rounds1niter35000jsize0.001maxLogScale0sinPeriod2accfind1spmax1gamma2lR2jsize0.02lsignu2jsize0.02nparallel12"
    load(paste0("data/chains_hierarchical_",filenamecore,".RData"))
    source('~/Research/multimodalhmc/experiments/sensor/sensor_model.R')
    lplot <- function(...) plot(..., type='l')
    x.d <- 2*N
    ##m1 <- numeric(12); nhop <- numeric(12)
    #burnin <- 30
    #for (cid in 1:12) {
    #    m1[cid] <- sum(chains[[cid]]$xv[(burnin+1):niter,3+8]>.4)
    #    nhop[cid] <- sum((chains[[cid]]$xv[(burnin):(niter-1),3+8]-0.4)*(chains[[cid]]$xv[(burnin+1):niter,3+8]-0.4)<0)
    #}
    #print(mean(m1)); print(sd(m1)); print(mean(m1)+qt(.975,11)*sd(m1)*c(-1,1)/sqrt(12))
    #boxplot(nhop)
    png(filename=paste0('figure/traceplot_hierarchical_chain1_Y3Y6_',filenamecore,'.png'), width=6, height=6, units="in", res=200)
    par(mfrow=c(3,2), mar=c(4,4,0,0), oma=c(0,0,2,.5))
    chainid <- 1
    U <- -sapply(1:niter, function(iii) { logtarget(chains[[chainid]]$xv[iii,1:x.d], chains[[chainid]]$lR2[iii], chains[[chainid]]$lsignu2[iii]) })
    logden <- -U + sapply(1:niter, function(iii) {
        llR2(chains[[chainid]]$lR2[iii], chains[[chainid]]$xv[iii,1:x.d]) +
            llsignu2(chains[[chainid]]$lsignu2[iii], chains[[chainid]]$xv[iii,1:x.d]) } )
    linewidth=.5; axislabel=.8
    lplot(U, bty='L', lwd=linewidth, xlab='', ylab=''); mtext(side=2, line=2, text='U', cex=axislabel); mtext(text=LETTERS[1], adj=0, line=0.5)
    lplot(logden, bty='L', lwd=linewidth, xlab='', ylab=''); mtext(side=2, line=2, text='log posterior density', cex=axislabel); mtext(text=LETTERS[2], adj=0, line=0.5)
    lplot(exp(chains[[chainid]]$lR2/2), bty='L', lwd=linewidth, xlab='', ylab=''); mtext(side=2, line=2, text=expression(R), cex=axislabel); abline(h=sqrt(R2_t), lty=2); mtext(text=LETTERS[3], adj=0, line=0.5)
    lplot(exp(chains[[chainid]]$lsignu2/2), log='y', bty='L', lwd=linewidth, xlab='', ylab=''); mtext(side=2, line=2, text=expression(sigma[nu]), cex=axislabel); abline(h=sqrt(signu2_t), lty=2); mtext(text=LETTERS[4], adj=0, line=0.5)
    lplot(chains[[chainid]]$xv[,3+N], bty='L', lwd=linewidth, xlab='iteration', ylab=''); mtext(side=2, line=2, text=expression(Y[3]), cex=axislabel); mtext(text=LETTERS[5], adj=0, line=0.5)
    lplot(chains[[chainid]]$xv[,6+N], bty='L', lwd=linewidth, xlab='iteration', ylab=''); mtext(side=2, line=2, expression(Y[6]), cex=axislabel); mtext(text=LETTERS[6], adj=0, line=0.5)
    dev.off()
}


## Rhat statistic ###
filenamecore <- "sensor8known3rounds1niter35000jsize0.001maxLogScale0sinPeriod2accfind1spmax1gamma2lR2jsize0.02lsignu2jsize0.02nparallel12"
load(paste0("data/chains_hierarchical_",filenamecore,".RData"))
source('~/Research/multimodalhmc/experiments/sensor/sensor_model.R')
lplot <- function(...) plot(..., type='l')
x.d <- 2*N
library(rstan)
outvar <- array(
    c(sapply(1:x.d, function(jj) { sapply(1:nparallel, function(val) { chains[[val]]$xv[,jj] }) }), # sensor locations (2*N)
    sapply(1:nparallel, function(val) { exp(chains[[val]]$lR2/2) }), # parameter R
    sapply(1:nparallel, function(val) { exp(chains[[val]]$lsignu2/2) })), # parameter sigma_nu
    dim=c(niter, nparallel, x.d+2))
burnin <- 7000
Rhatvals <- apply(outvar[(burnin+1):niter,,], 3, Rhat)


## time to a dominant mode
Us <- -sapply(1:nparallel, function(chainid) {
    sapply(1:niter, function(iii) { logtarget(chains[[chainid]]$xv[iii,1:x.d], chains[[chainid]]$lR2[iii], chains[[chainid]]$lsignu2[iii]) })
    })
logdens <- -Us + sapply(1:nparallel, function(chainid) {
        sapply(1:niter, function(iii) { llR2(chains[[chainid]]$lR2[iii], chains[[chainid]]$xv[iii,1:x.d]) +
                llsignu2(chains[[chainid]]$lsignu2[iii], chains[[chainid]]$xv[iii,1:x.d]) }) })
DomModeCriterion # this is computed using the Markov chains constructed by tempered Hamitonian transitions
FirstDomModeIter <- apply(logdens, 2,
    function(vec) { csum <- cumsum(vec);
        which.max((csum[21:niter]-csum[1:(niter-20)])/20 > DomModeCriterion) }) # first time the chain reached one of the two dominant modes
c(mean(FirstDomModeIter), sd(FirstDomModeIter))/niter*duration


