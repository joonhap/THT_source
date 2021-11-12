## this file creates example trajectories like those constructed by tempered Hamiltonian transitions
lplot <- function(...) plot(..., type='l')
source("../leapfrog.R")
#set.seed(Sys.time())
#set.seed(109241)

#########################
## LOAD VARIOUS MODELS ##
#########################
## mixture of d-dim normals
x.d <- 10000
source("../mvnorm/mixture2normal_Nd.R"); logtarget <- target; gradlt <- gd.target
x.init <- modes[1,] + comp.sd*rnorm(x.d); power <- 2; power.truth <- 2;  set.seed(912211); set.seed(Sys.time())
##bimodal 2-dim normal
#source("../mvnorm/mixture2normal_2d.R"); logtarget <- o.target; gradlt <- gd.o.target
#x.init <- modes[1,] + o.comp.sd*rnorm(2); power <- 2
##bimodal 1-dim normal
#source("../mvnorm/mixture2normal_1d.R"); logtarget <- o.target; gradlt <- gd.o.target
#x.init <- modes[1] + o.comp.sd[1]*rnorm(1); power <- 2
##U(x) propto |x|^power
#x.d <- 11; power <- 2;  x.init <- rnorm(x.d); power.truth <- 2
#logtarget <- function(x) { -sum(x^2)^(power.truth/2) }
#gradlt <- function(x) { -power.truth*sum(x^2)^(power.truth/2-1)*x }
#set.seed(323509)


jsize <- .1
njumps <- 1
massInv <- 1
sin_period <- 3000
maxLogScale <- 24 # this is equal to 4 times eta_* in the manuscript
massScalingLvls <- exp(maxLogScale/2-maxLogScale/2*cos((0:(sin_period-1))*(2*pi)/sin_period))
nMSFlvls <- length(massScalingLvls)/2
simlen <- nMSFlvls+1


xv <- matrix(NA, simlen, 2*x.d)
acc <- logical(simlen) # logical, is a proposal accepted?
last_accepted <- matrix(NA, simlen, x.d) # the last accepted proposal
Utrace <- numeric(simlen)
Ktrace <- numeric(simlen)
MSFtrace <- numeric(simlen)
acc[1] <- TRUE
last_accepted[1,] <- x.init
Lambda <- .5
x <- x.init
Utrace[1] <- U <- -logtarget(x)
k <- 0#sample(c(0:4, nMSFlvls-4:1), 1)#round(nMSFlvls*.4) # initial mass scaling factor
v <- rnorm(x.d, 0, sqrt(massInv/massScalingLvls[2*k+1])) #^(-1/(power+2))) # this is actually v-tilde in the manuscript
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
    if (Unew+.5*sum(vnew*vnew/massInv)*massScalingLvls[2*k+1]-x.d/2*log(massScalingLvls[2*k+1]) - H < -log(Lambda)) {
        acc[n] <- TRUE
        last_accepted[n,] <- xnew
    } else {
        last_accepted[n,] <- last_accepted[n-1,]
    }
    if (n >= simlen) { break }
}

Ktildetrace <- Ktrace * MSFtrace
Kenergytrace <- Ktildetrace - x.d/2*log(MSFtrace) # negative log of extended velocity density (up to additive const.)
Htrace <- Utrace + Ktrace
Htildetrace <- Utrace + Ktildetrace
Henergytrace <- Utrace + Kenergytrace
#closest_mode <- apply(xv[,1:x.d], 1, function(x) ifelse(sum((x-modes[1,])^2)<sum((x-modes[2,])^2), 1, 2))
#dist_to_clsmd <- sapply(1:simlen, function(n) sqrt(sum((xv[n,1:x.d]-modes[closest_mode[n],])^2)))

x_parallel <- xv[,1:x.d]%*%direction; lplot(x_parallel)
lplot(xv[,1]); points(xv[,1], pch='+', col='green', cex=0.8); abline(h=modes[,1], col='blue', lty=2)
lplot(Henergytrace, log='y'); abline(h=Henergytrace[1], col='blue', lty=2); lines(MSFtrace^(power/(power+2)), col='cyan')
lplot(log(Htildetrace)); abline(h=log(Htildetrace[1]), col='blue', lty=2)
lplot(log(MSFtrace)[],log(Henergytrace)[]); points(log(MSFtrace[c(1,simlen)]), log(Henergytrace[c(1,simlen)]), col=c('red', 'cyan'))

print(Henergytrace[1]); print(Henergytrace[simlen])

par(mfrow=c(2,1)); lplot(acc); lplot(log(MSFtrace)); par(mfrow=c(1,1))

vbar <- xv[,x.d+1]*MSFtrace^(1/(power+2)); lplot(vbar); points(vbar, col='green', pch='+')#^(1/(power.truth+2)))
xbar <- xv[,1]*MSFtrace^(-1/(power+2)); lplot(xbar); points(xbar, col='green', pch='+')#^(-1/(power.truth+2)))
corrected_timescale <- cumsum(jsize*massScalingLvls[c(2*1:nMSFlvls,2)]^(2/(power+2)-2/(power.truth+2)))
lplot(corrected_timescale, vbar); points(corrected_timescale, vbar, pch='+', col='green')
lplot(corrected_timescale, xbar); points(corrected_timescale, xbar, pch='+', col='green')

Kbar <- Ktrace*MSFtrace^(2/(power+2)); Ubar <- Utrace*MSFtrace^(-power.truth/(power+2)); Hbar <- Kbar + Ubar
lplot(Hbar, ylim=range(c(Hbar,Ubar,Kbar))); lines(Ubar, col='red'); lines(Kbar, col='cyan'); abline(h=Hbar[1], lty=2, col='blue')
print(Hbar[1]); print(Hbar[simlen])
#lplot(log(Ktildetrace+Utrace*MSFtrace^((power-2)/(power+2))))
#lplot(Ktildetrace); lines(Utrace, col='red')


## mixture of two d=10000 dimensional Gaussian.
trajectoryplot_2d <- function() {
    ##readme <- "used seed 912211 after x.init generation. jsize=.1, sin_period=3000, initial k=0, modeDist=400, comp.sd=1, x.d=10000, power=2"
    ##save(xv, simlen, readme, file='data/trajectory_bimodal_jsize0.1_sinperiod3000_d10000.RData')
    load('data/trajectory_bimodal_jsize0.1_sinperiod3000_d10000.RData')
    ## parallel and a perpendicular projection
    perp <- local({ ran.vec <- sin(1:x.d);
        perp.raw <- ran.vec - sum(ran.vec * direction) * direction
        perp.raw / sqrt(sum(perp.raw^2)) }) # a direction perpendicular to the inter-mode direction
    x_parallel <- c(xv[,1:x.d]%*%direction)
    x_perp <- c(xv[,1:x.d]%*%perp)
    cexsize <- 1
    parallel_range <- range(x_parallel); perp_range <- range(x_perp)
    plot(x_parallel, x_perp, type='l', pch='.', col='black', axes=FALSE, xlim=range(x_parallel), ylim=range(x_perp), xlab='', ylab='')
    axis(side=1, cex.axis=cexsize); axis(side=2, cex.axis=cexsize)
    points(x_parallel[simlen], x_perp[simlen], col='red', pch="+", cex=1.8); text(x_parallel[simlen], x_perp[simlen], labels=simlen-1, col='red', pos=2, cex=cexsize)
    points(x_parallel[1], x_perp[1], col='red', pch="+", cex=1.8); text(x_parallel[1], x_perp[1], labels=0, col='red', pos=2, cex=cexsize)
    labels <- c(400, 600, 800, 1000, 1200)
    points(x_parallel[labels+1], x_perp[labels+1], pch='+', col='blue')
    text(x_parallel[labels+1], x_perp[labels+1], labels=labels, col='blue', pos=4, cex=cexsize)
    mtext(side=1, line=2.2, text='parallel component', cex=cexsize)
    mtext(side=2, line=2.3, text='perpendicular component', cex=cexsize)
}
## run tempered Hamiltonian transitions for mixture of two Gaussians (d=10000)
source('THT.R')
##
set.seed(246209)
jsize <- .2
njumps <- 1
massInv <- 1
sin_period <- 2000
maxLogScale <- 24
massScalingLvls <- exp(maxLogScale/2-maxLogScale/2*cos((0:(sin_period-1))*(2*pi)/sin_period))
ksupport <- c(0:4, (sin_period/2-4):(sin_period/2-1))
#ksupport <- 0
len_ksupport <- length(ksupport)
rk <- function() { sample(ksupport, 1) }
dk <- function(k) { (k %in% ksupport) }
acc_find <- len_ksupport
spmax <- sin_period/2+1
gamma <- 2
niter <- 100
##
##out <- THT(x.init, logtarget, gradlt, niter=niter, jsize=jsize, njumps=njumps,
##    massInv=massInv, massScalingLvls=massScalingLvls, rk=rk, dk=dk, acc_find=acc_find, spmax=spmax, gamma=gamma)
##readme <- "jsize .2, sin_period 2000, maxLogScale 24, len_ksupport 9, gamma2, niter 100, seed 246209, x.d 10000"
##save(out, readme, file='data/modulatedHMC_bimodal_d10000_jsize0.2_sinperiod2000_niter100.RData')
load('data/modulatedHMC_bimodal_d10000_jsize0.1_sinperiod3000_niter100_2.RData')
## create a plot for an example trajectory and a traceplot for the parallel component
gen_trajectory_traceplot_bimodal_d10000 <- 0
if (gen_trajectory_traceplot_bimodal_d10000) {
    png('figures/bimodal_trajectory_traceplot_jsize0.1_sinperiod3000_d10000.png', width=9, height=4, units='in', res=150)
    par(mfrow=c(1,2), mar=c(3.5,3.5,1,1))
    trajectoryplot_2d()
    x_parallel <- c(out$xv[,1:x.d]%*%direction)
    lplot(0:niter, x_parallel, axes=FALSE, xlab='', ylab='')
    axis(side=1); axis(side=2)
    mtext(side=1, line=2, text='iteration')
    mtext(side=2, line=2.2, text='parallel component')
    dev.off()
}


## xbar, vbar, Hbar for gamma=1 and 3 and for the bimodal example.
gen_xvHbar_gamma13_bimodal <- 1
if(gen_xvHbar_gamma13_bimodal) {
    ##readme <- "xbarN, vbarN, KbarN, UbarN, HbarN are the xbar, vbar, etc. variables where power=power.truth=N. jsize=0.1, sin_period=1600, maxLogScale=24, seed=109241 for U(x)=|x|^gamma and seed=323509 for bimodal U"
    ##save(xbar1,vbar1,Kbar1,Ubar1,Hbar1,xbar3,vbar3,Kbar3,Ubar3,Hbar3,xbar_bm,vbar_bm,Hbar_bm,Kbar_bm,Ubar_bm, readme, file="data/xbar_vbar_Hbar_gamma_1_3_bimodal_Jun8_2021.RData")
    load('data/xbar_vbar_Hbar_gamma_1_3_bimodal_Jun8_2021.RData')
    simlen=1600/2+1
    library(tidyverse); library(cowplot); library(scales)
    linewidth <- .1
    xvdata <- rbind(
        data.frame(vbar=vbar1, xbar=xbar1, target='gamma1', step=1:simlen),
        data.frame(vbar=vbar3, xbar=xbar3, target='gamma3', step=1:simlen),
        data.frame(vbar=vbar_bm, xbar=xbar_bm, target='bimodal', step=1:simlen)
    ) %>%
    pivot_longer(cols=c(vbar, xbar), names_to='variable', values_to='value') %>%
    mutate(target= recode_factor(target, `gamma1`="gamma==1", `gamma3`="gamma==3", `bimodal`="bimodal")) %>%
    mutate(variable= recode_factor(variable, `xbar`="bar(x)", `vbar`="bar(v)"))
    xplot <- ggplot(xvdata%>%filter(variable=='bar(x)')) +
        facet_wrap(facets=vars(target), nrow=3, ncol=1, scales="free_y", label="label_parsed") +
        geom_line(aes(x=step, y=value, linetype=variable), size=linewidth) + labs(x='leapfrog step', linetype="") + theme_bw() +
        theme(axis.title.y=element_blank(), plot.margin=margin(t=0,b=0, unit="pt"),
            legend.margin=margin(t=0,b=0, unit="pt"), legend.position='top', #legend.key=element_blank(),
            strip.background=element_blank(), strip.text=element_blank()) +
        scale_linetype_discrete(labels= parse_format())
    vplot <- ggplot(xvdata%>%filter(variable=='bar(v)')) +
        facet_wrap(facets=vars(target), nrow=3, ncol=1, scales="free_y", label="label_parsed") +
        geom_line(aes(x=step, y=value, linetype=variable), size=linewidth) + labs(x='leapfrog step', linetype="") + theme_bw() +
        theme(axis.title.y=element_blank(), plot.margin=margin(t=0,b=0, unit="pt"),
            legend.margin=margin(t=0,b=0,unit="pt"), legend.position='top', #legend.key=element_blank(),
            strip.background=element_blank(), strip.text=element_blank()) +
        scale_linetype_discrete(labels= parse_format())
    KUHdata <- rbind(
        data.frame(Hbar=Hbar1, Kbar=Kbar1, Ubar=Ubar1, target='gamma1', step=1:simlen, facetLab="energy"),
        data.frame(Hbar=Hbar3, Kbar=Kbar3, Ubar=Ubar3, target='gamma3', step=1:simlen, facetLab="energy"),
        data.frame(Hbar=Hbar_bm, Kbar=Kbar_bm, Ubar=Ubar_bm, target='bimodal', step=1:simlen, facetLab="energy")
    ) %>%
    pivot_longer(cols=c(Hbar, Kbar, Ubar), names_to='energy', values_to='value') %>%
    mutate(target= recode_factor(target, `gamma1`="gamma==1", `gamma3`="gamma==3", `bimodal`="bimodal")) %>%
    mutate(energy= recode_factor(energy, `Hbar`="U(bar(x))+K(bar(v))", `Ubar`="U(bar(x))", `Kbar`="K(bar(v))"))
    KUHplot <- ggplot(KUHdata) +
        facet_wrap(facets=vars(target), nrow=3, ncol=1, scales="free_y", label="label_parsed", strip.position="right") +
        geom_line(aes(x=step, y=value, color=energy), size=linewidth) + xlab("leapfrog step") +
        ggthemes::scale_color_colorblind(labels=parse_format()) + theme_bw() +
        theme(axis.title.y=element_blank(), plot.margin=margin(b=0,t=0, unit="pt"),
            legend.position='top', legend.margin=margin(t=0,b=0,unit="pt"), legend.title=element_blank())
    combined_plot <- plot_grid(xplot, vplot, KUHplot, ncol=3, rel_widths=c(1,1,1.1))
    ggsave(filename="figures/xbar_vbar_gamma_1_3_bimodal_jsize0.1_sinPeriod1600_maxLogScale24.png", plot=combined_plot, width=8, height=6, units="in", dpi=250)
}


## vbar for misspecified gamma
gen_vbar_misspecified_gamma <- 1
if(gen_vbar_misspecified_gamma) {
    ##readme <- "vbar_trueA_estB, where A is the true value of gamma ('power.truth'), B is the estimated gamma ('power') used to scale leapfrog step sizes. For A=1,B=2, jsize=.4, sin_period=6000; for A=2,B=1, jsize=.01, sin_period=2000; for A=2,B=3, jsize=.1, sin_period=4000; for A=3,B=2, jsize=.01, sin_period=4000.  seed=323509 for all cases."
    ##save(vbar_true1_est2, vbar_true2_est1, vbar_true2_est3, vbar_true3_est2, simlen_true1_est2, simlen_true2_est1, simlen_true2_est3, simlen_true3_est2, readme, file="data/vbar_gamma_misspecified_Jun8_2021.RData")
    load('data/vbar_gamma_misspecified_Jun8_2021.RData')
    library(tidyverse); library(cowplot); library(scales)
    linewidth <- .1
    vdata <- rbind(
        data.frame(vbar=vbar_true1_est2, gamma=1, gammahat=2, step=1:simlen_true1_est2, label1=1, label2=1),
        data.frame(vbar=vbar_true2_est3, gamma=2, gammahat=3, step=1:simlen_true2_est3, label1=2, label2=2),
        data.frame(vbar=vbar_true2_est1, gamma=2, gammahat=1, step=1:simlen_true2_est1, label1=3, label2=3),
        data.frame(vbar=vbar_true3_est2, gamma=3, gammahat=2, step=1:simlen_true3_est2, label1=4, label2=4)
    ) %>%
        mutate(gamma= recode_factor(gamma, `1`="gamma==1", `2`="gamma==2", `3`="gamma==3")) %>%
        mutate(gammahat= recode_factor(gammahat, `1`="hat(gamma)==1", `2`="hat(gamma)==2", `3`="hat(gamma)==3")) %>%
        mutate(label1= recode_factor(label1, `1`="gamma==1*',  '*hat(gamma)*'='*2", `2`="gamma==2*',  '*hat(gamma)*'='*3", `3`="gamma==2*',  '*hat(gamma)*'='*1", `4`="gamma==3*',  '*hat(gamma)*'='*2")) %>%
        mutate(label2= recode_factor(label2, `1`="epsilon==0.4*', '*K*'='*3000", `2`="epsilon==0.1*', '*K*'='*2000", `3`="epsilon==0.01*', '*K*'='*1000", `4`="epsilon==0.01*', '*K*'='*2000"))
    vplot <- ggplot(vdata) +
        facet_wrap(facets=vars(label1,label2), nrow=1, ncol=4, scales="free", label="label_parsed") +
        geom_line(aes(x=step, y=vbar), size=linewidth) + xlab('leapfrog step') + ylab(expression(bar(v))) +
        geom_point(aes(x=step, y=vbar), size=linewidth, color='green') + xlab('leapfrog step') + ylab(expression(bar(v))) + theme_bw()
    ggsave(filename="figures/vbar_gamma_misspecified_maxLogScale24.png", plot=vplot, width=8, height=2.5, units="in", dpi=250)
}


## an example figure
gen_draft_fig <- 0 ; ## USED SEED 109241
if (gen_draft_fig) {
    ##save(jsize, maxLogScale, x.init, sin_period, massInv, njumps, massScalingLvls, simlen, xv, file='data/temperedTrajectory_bimodal_1d_norm_maxLogScale24jsize0.1.RData')
    load('data/temperedTrajectory_bimodal_1d_norm_maxLogScale24jsize0.1.RData')
    ##pdf(paste0('figures/temperedTrajectory_bimodal_1d_norm_maxLogScale',maxLogScale,'jsize',jsize,'.pdf'), width=11, height=4)
    png(paste0('figures/temperedTrajectory_bimodal_1d_norm_maxLogScale',maxLogScale,'jsize',jsize,'.png'), width=11, height=4, units="in", res=150)
    mag <- 1.3
    par(mar=c(3.6,3.5,1.5,1), mfrow=c(1,2))
    ## panel 1 (trajectory)
    lplot(xv[,1], xlab='', ylab='', axes=FALSE); points(xv[,1], pch='+', col='green', cex=0.8); abline(h=modes, col='blue', lty=2)
    mtext(side=1, line=2.5, text='leapfrog steps', cex=mag); axis(side=1, cex.axis=mag)
    mtext(side=2, line=2.5, text=expression(x), cex=mag); axis(side=2, cex.axis=mag)
    ## panel 2 (mass scaling factor)
    #lplot(MSFtrace, xlab='', ylab='', log='y', axes=FALSE)
    #mtext(side=1, line=2.5, text='leapfrog steps'); axis(side=1)
    #mtext(side=2, line=2.5, text=expression('mass scale factor, '*alpha)); axis(side=2)
    ## panel 3 (energy)
    lplot(Henergytrace, log='y', ylim=c(min(Htildetrace), max(Htildetrace)*8), xlab='', ylab='', axes=FALSE); lines(MSFtrace^.5, lty=2); abline(h=Henergytrace[1], col='blue', lty=2)
    mtext(side=1, line=2.5, text='leapfrog steps', cex=mag); axis(side=1, cex.axis=mag)
    axis(side=2, cex.axis=mag); #mtext(side=2, line=2.5, text=expression('Hamiltonian H(x,'*tilde(v)*')= U(x) + '*frac(1,2)*tilde(v)^T*'('*alpha*M*')'*tilde(v))
    legend('topleft', bty='n', lty=c(1,2), legend=c(expression('Hamiltonian H(x,k,'*tilde(v)*')'), expression(e^eta)), cex=mag, ncol=2)
    dev.off()
}

## vbar, xbar plot for various power(=power.truth)
gen_vxbar_plot <- 0
if(gen_vxbar_plot) {
    ##readme <- "xbarN, vbarN are the xbar, vbar variables where power=power.truth=N. jsize=0.2, sin_period=1000, maxLogScale=24"
    ##save(xbar2,vbar2,xbar1.2,vbar1.2,xbar3,vbar3, readme, file="data/xbar_vbar_gamma_1.2_2_3_Jun8_2021.RData")
    load('data/xbar_vbar_gamma_1.2_2_3_Jun8_2021.RData')
    library(tidyverse)
    data <- rbind(
        cbind(vbar=vbar2, xbar=xbar2, gamma=2, step=1:simlen),
        cbind(vbar=vbar1.2, xbar=xbar1.2, gamma=1.2, step=1:simlen),
        cbind(vbar=vbar3, xbar=xbar3, gamma=3, step=1:simlen)
    ) %>% as.data.frame %>% pivot_longer(cols=c(vbar, xbar), names_to='variable', values_to='value')
    p <- ggplot(data) +
        facet_wrap(facets=vars(variable,gamma), nrow=2, ncol=3, scales="free_y") +
        geom_line(aes(x=step, y=value))
    ggsave(filename="figures/xbar_vbar_gamma_1.2_2_3_jsize0.2_sinPeriod1000_maxLogScale24.png", plot=p, width=7, height=4.5, units="in", dpi=150)
}

