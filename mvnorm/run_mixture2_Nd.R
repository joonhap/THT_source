library(tidyverse)

### compare mass-enhanced HMC to plain HMC
source("../meHMC/meHMC.R") ### load mass-enhanced HMC

x.d_list <- c(2, 20, 200)
niter <- 2000
mass <- 1
acc_find <- 1; spmax <- 1000
massScaling <- 100

THTResult <- lapply(x.d_list, function(x.d) {
    ### load target distribution
    source("mixture2normal_Nd.R", local=TRUE)

    set.seed(134653)
    set.seed(Sys.time())
    x.init <- modes[1,]+rnorm(x.d)*comp.sd

    THT.time <- as.numeric(system.time(
        THT.out <- msHMC(x.init=x.init, logtarget=target, gradlt=gd.target, niter=niter,
            jsize=.1, njumps=1, massInv=1/mass, massScaling=massScaling,
            acc_find=acc_find, spmax=spmax, randomize.eps=TRUE)
    )[3])
    pilot_HMC.time <- as.numeric(system.time(
        pilot_HMC.out <- msHMC(x.init=x.init, logtarget=target, gradlt=gd.target, niter=500,
            jsize=.1, njumps=5, massInv=1/mass, massScaling=1,
            acc_find=acc_find, spmax=1, randomize.eps=TRUE)
    )[3]) # test run to see how long an HMC iteration take on average
    niter_HMC <- round(THT.time / pilot_HMC.time * 500)
    HMC.time <- as.numeric(system.time(
        HMC.out <- msHMC(x.init=x.init, logtarget=target, gradlt=gd.target, niter=niter_HMC,
            jsize=.1, njumps=5, massInv=1/mass, massScaling=1,
            acc_find=acc_find, spmax=1, randomize.eps=TRUE)
    )[3]) # run HMC for roughly the same amount of time as for THT

    # name of the closest mode
    THT.cm <- apply(THT.out[,1:x.d], 1, closest_mode)
    HMC.cm <- apply(HMC.out[,1:x.d], 1, closest_mode)

    # parallel and a perpendicular projection
    perp <- { ran.vec <- rnorm(x.d);
    perp.raw <- ran.vec - sum(ran.vec * direction) * direction
    perp.raw / sqrt(sum(perp.raw^2)) } # a direction perpendicular to the inter-mode direction
    THT.parallel <- apply(THT.out[,1:x.d], 1, function(v) sum(v*direction))
    THT.perp <- apply(THT.out[,1:x.d], 1, function(v) sum(v*perp))
    HMC.parallel <- apply(HMC.out[,1:x.d], 1, function(v) sum(v*direction))
    HMC.perp <- apply(HMC.out[,1:x.d], 1, function(v) sum(v*perp))

    return(list(THT.out=THT.out, HMC.out=HMC.out, THT.time=THT.time,
        HMC.time=HMC.time, HMC.niter=niter_HMC, THT.cm=THT.cm, HMC.cm=HMC.cm,
        THT.parallel=THT.parallel, THT.perp=THT.perp,
        HMC.parallel=HMC.parallel, HMC.perp=HMC.perp))
})
description <- paste0('d=2,20,200, comparison of HMC and THT, where for the latter acc_find=', acc_find,
    'spmax=', spmax, 'massScaling=', massScaling, '. The plain HMC was run with njumps=5.')
save(THTResult, x.d_list, description, file=paste0('data/Feb25_mixture2_d2_20_200_massScaling_',massScaling,
    '_accfind',acc_find,'_spmax',spmax,'_results.RData'))
load('data/Feb25_mixture2_d2_20_200_massScaling_100_accfind1_spmax1000_results.RData')

THT.df <- rbind(
    data.frame(iter=1:2000, time=THTResult[[1]]$THT.time*(1:2000)/2000,
        x.d=2, acc_find=1, spmax=100, massScaling=30, method="mass-enhanced HMC",
        cm=THTResult[[1]]$THT.cm,
        parallel=THTResult[[1]]$THT.parallel,
        perp=THTResult[[1]]$THT.perp),
    data.frame(iter=1:THTResult[[1]][["HMC.niter"]],
        time=THTResult[[1]]$HMC.time*(1:THTResult[[1]]$HMC.niter)/THTResult[[1]]$HMC.niter,
        x.d=2, acc_find=1, spmax=1, massScaling=1, method="HMC",
        cm=THTResult[[1]]$HMC.cm,
        parallel=THTResult[[1]]$HMC.parallel,
        perp=THTResult[[1]]$HMC.perp),
    data.frame(iter=1:2000, time=THTResult[[2]]$THT.time*(1:2000)/2000,
        x.d=20, acc_find=1, spmax=100, massScaling=30, method="mass-enhanced HMC",
        cm=THTResult[[2]]$THT.cm,
        parallel=THTResult[[2]]$THT.parallel,
        perp=THTResult[[2]]$THT.perp),
    data.frame(iter=1:THTResult[[2]][["HMC.niter"]],
        time=THTResult[[2]]$HMC.time*(1:THTResult[[2]]$HMC.niter)/THTResult[[2]]$HMC.niter,
        x.d=20, acc_find=1, spmax=1, massScaling=1, method="HMC",
        cm=THTResult[[2]]$HMC.cm,
        parallel=THTResult[[2]]$HMC.parallel,
        perp=THTResult[[2]]$HMC.perp),
    data.frame(iter=1:2000, time=THTResult[[3]]$THT.time*(1:2000)/2000,
        x.d=200, acc_find=1, spmax=100, massScaling=30, method="mass-enhanced HMC",
        cm=THTResult[[3]]$THT.cm,
        parallel=THTResult[[3]]$THT.parallel,
        perp=THTResult[[3]]$THT.perp),
    data.frame(iter=1:THTResult[[3]][["HMC.niter"]],
        time=THTResult[[3]]$HMC.time*(1:THTResult[[3]]$HMC.niter)/THTResult[[3]]$HMC.niter,
        x.d=200, acc_find=1, spmax=1, massScaling=1, method="HMC",
        cm=THTResult[[3]]$HMC.cm,
        parallel=THTResult[[3]]$HMC.parallel,
        perp=THTResult[[3]]$HMC.perp)
)

## plots
## closest_mode plot
cm_gplot <- ggplot(THT.df, aes(time, cm)) +
    facet_grid(method~x.d, labeller=labeller(x.d=c('d=2','20','200')%>%`names<-`(c(2,20,200))), scale="free_x") +
    geom_line() + theme(text=element_text(size=16)) +
    xlab('runtime (s)') + ylab('closest density component') +
    scale_y_continuous(breaks=c(1,2))
cm_gplot
save_cmplot <- TRUE
if (save_cmplot) { ggsave(paste0('figures/closest_mode_x.d2,20,200_massScale',massScaling,'_accFind',acc_find,
    '_spmax', spmax, '_niter', niter,'_HMC_THT.png'), plot=cm_gplot, width=9, height=4.5, dpi=150) }
## sample points (parallel/perpendicular)
scatter_plot <- ggplot(THT.df, aes(parallel, perp)) +
    facet_grid(x.d~method, labeller=labeller(x.d=c('d=2','20','200')%>%`names<-`(c(2,20,200)))) +
    geom_point(size=0.08) + xlim(-1.4,1.4) + ylim(-0.4,0.4) + coord_fixed() +
    geom_point(aes(x=-1, y=0), shape=3, color='red') +
    geom_point(aes(x=1, y=0), shape=3, color='red') +
    theme(text=element_text(size=16)) + xlab('parallel component') + ylab('perpendicular component')
scatter_plot
save_scatterplot <- TRUE
if (save_scatterplot) { ggsave(paste0('figures/samples_x.d2,20,200_massScale',massScaling,'_accFind', acc_find,
    '_spmax', spmax, '_niter', niter, '_HMC_THT.png'), plot=scatter_plot, width=9, height=4.5, dpi=150) }




jumpiters <- which(closest_modes[1:(niter-1)]!=closest_modes[2:niter])+1; lapply(props[jumpiters],
    function(l) l[,"jumpmethod"])

par(ask=1); it <- jumpiters[3]; for(portion in 1:ceiling(dim(props[[it]])[1]/10)) { plotiters <- portion*
    10+-9:min(1,dim(props[[it]])[1]-portion*10); todraw<-props[[it]][plotiters,1:2];
    plot(todraw, xlim=c(-5,5), ylim=c(-5,5), type='p', cex=3, pch='.', col=props[[it]][plotiters,"jumpmethod"]);
    points(modes[,1], modes[,2], pch='+', col='red') }; par(ask=0)

print(dim(props[[it]])[1])
plot(props[[it]][,"jumpmethod"], col=props[[it]][,"acc"])

plot(Result[,1], Result[,2], xlim=c(-2.5,2.5), ylim=c(-2.5,2.5),
    pch='+', type='p'); points(modes[,1], modes[,2], pch='+', col='red')

closest_modes <- apply(Result[,1:x.d], 1, closest_mode); plot(closest_modes, type='l')

plcmp <- c(1,10); plot(Result[,plcmp[1]], Result[,plcmp[2]], pch='+',
    xlim=c(-2.5,2.5),ylim=c(-2.5,2.5)); points(modes[,plcmp[1]], modes[,plcmp[2]], pch='+', col='red')

dist_to_closest_mode <- sapply(1:niter, function(iter) {
    sum((Result[iter,1:x.d] - modes[closest_modes[iter],])^2)
    }); plot(dist_to_closest_mode, type='l')

plot(Result[,plcmp[1]], type='l')

iter <- niter; plot(modes[closest_modes[iter],], ylim=c(-2.5,2.5)); points(
    Result[iter,1:x.d], col='red')



