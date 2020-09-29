require(tidyverse)
require(ggplot2)
require(deSolve)
require(ggpubr)
require(patchwork)
require(triangle)
require(ggrepel)
require(ggpubr)
require(ggthemes)
require(extraDistr)
require(mvtnorm)

source("code/simulation_functions.R")
source("code/model_funcs.R")

setwd("~/Documents/GitHub/ct_dynamics_preprint/")

cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

#### Functions for Viral Load Distributions: ####
## Describes average viral load as a function of time since infection
viral_load_func <- function(pars, a, convert_ct=TRUE){
  viral_load <- viral_load_func_asymp(pars, obs_t=a, convert_ct)
  viral_load
}

## Function to give probability of observing x given age a and the viral kinetics curve
p_a_use <- function(x, a, pars, obs_model="normal", log_prob=FALSE, ct=TRUE, additional_detect) {
  if (!ct) {
    x <- pars["intercept"]-log2(10)*(x-pars["LOD"])
  }
  p <- p_a(x, a, pars, viral_load_func(pars, a, convert_ct=TRUE), obs_model, log_prob, additional_detect)
  p
}

## Proportion of individuals with detectable viral loads for a given age
prop_det_use <- function(a, pars, obs_model="normal", additional_detect) {
  phi <- prop_detectable(a, pars, viral_load_func(pars, a, convert_ct=TRUE), obs_model, additional_detect)
  phi
}

## Function to generate random viral loads given a single age a and the viral kinetics curve
rx_a <- function(n, a, pars, obs_model="normal", additional_detect=FALSE) {
  ct_mean <- viral_load_func(pars, a, convert_ct=TRUE)
  if(obs_model=="normal") {
    xs <- rnorm(n, mean=ct_mean, sd=pars["obs_sd"])
  } else {
    xs <- rgumbel(n, mu=ct_mean, sigma=pars["obs_sd"])
  }
  xs <- xs[xs < pars["intercept"]]
  if(additional_detect){
    t_switch <-  pars["t_switch"] + pars["desired_mode"] + pars["tshift"]
    if (a > t_switch) {
      ts <- a - t_switch
      additional_prob <- (1-pars["prob_detect"]*pars["t_unit"])^ts
      xs <- xs[rbinom(n, size=1, prob=additional_prob)==1]
    }
  }
  xs
}

## Function to generate random viral load trajectories:
rx <- function(pars, obs_model="normal", additional_detect=FALSE, a.range=1:35, corr=.8) {
  ct_mean <- viral_load_func(pars, a.range, convert_ct=TRUE)
  if (obs_model=="normal") {
    Sigma <- matrix(pars["obs_sd"]^2, nrow=length(a.range), ncol=length(a.range))
    for (i in 1:length(a.range)) {
      for (j in 1:length(a.range)) {
        Sigma[i,j] <- Sigma[i,j]*corr^(abs(i-j))
      }
    }
    xs <- rmvnorm(1, mean=ct_mean, sigma=Sigma)[1,]
    xs[xs >= pars["intercept"]] <- NA
    xs[xs < 1] <- 1
  } else {
    ind <- extraDistr::rgumbel(1, mu=0, sigma=pars["obs_sd"])
    xs <- ct_mean + ind
    xs[xs >= pars["intercept"]] <- NA
    xs[xs < 1] <- 1
  }
  if (additional_detect) {
    t_switch <-  pars["t_switch"] + pars["desired_mode"] + pars["tshift"]
    t_switch_int <- ceiling(t_switch)
    t_switch_inc <- t_switch_int-t_switch
    ts <- a.range[a.range >= t_switch_int]
    addl_det <- c(rbinom(1, 1, (1-pars["prob_detect"]*pars["t_unit"])^t_switch_inc),
                  rbinom(length(ts)-1, 1, (1-pars["prob_detect"]*pars["t_unit"])))
    last.day <- min(c(max(a.range)+1,ts[addl_det==0]))-1
    xs[a.range > last.day] <- NA
  }
  xs
}

#### Viral Load Kinetics Parameters: ####
pars <- c(beta = 0.1, tshift = 2, desired_mode = 2, viral_peak = 9, 
          obs_sd = 6, true_0 = -3, intercept = 40, LOD = 3, limit_quantification = 2, 
          incu = 5, t_switch = 16, level_switch = 3.7, 
          wane_rate2 = 30, prob_detect = 0.1, t_unit = 1)

#### Simulating the epidemic: ####

set.seed(0)
## Entire population size
population_n <- 1000000
## Sample size
sample_n <- 500
## Duration of epidemic in days
times <- 0:365

## Simulate the epidemic process
#epidemic_process <- simulate_epidemic_process(population_n,0.1,times)
#epidemic_process$plot
best_pars_fitting <- read.csv("pars/seir_fits_best_pars.csv")
best_pars <- as.numeric(best_pars_fitting[best_pars_fitting$institution == "NH 1",])
names(best_pars) <- colnames(best_pars_fitting)

best_pars["R0"] <- 2.5 ## To hard-code R0
best_pars["I0"] <- 1000 ## To hard-code number of infections at time 0

seir_pars <- c("R0"=best_pars["R0"],"gamma"=best_pars["gamma"],"sigma"=best_pars["sigma"],"I0"=best_pars["I0"],"recovered0"=0)
names(seir_pars) <- c("R0","gamma","sigma","I0","recovered0")
epidemic_process <- simulate_seir_process(population_n,seir_pars,times)

x <- seq(0,40,by=0.01)
lastday <- 35
incidence <- epidemic_process$incidence[1:200]

#### Calculating Growth Rates: ####
GR_daily <- log(incidence[2:200]/incidence[1:199])
GR_daily <- ifelse(is.finite(GR_daily), GR_daily, NA)
GR_full <- NULL
GR_25 <- NULL
for (i in (lastday+1):length(GR_daily)) {
  GR_full <- c(GR_full,mean(GR_daily[(i-lastday):(i-1)], na.rm=TRUE))
  GR_25 <- c(GR_25,mean(GR_daily[(i-25):(i-1)], na.rm=TRUE))
}
GR <- data.frame(Days=1:length(GR_daily), GR=c(rep(NA,lastday),GR_full), GR25=c(rep(NA,lastday),GR_25), DailyGR=GR_daily)

#### Calculating viral load and age distributions: ####
Ct_res <- matrix(0,nrow=length(incidence), ncol=length(x))
age_res <- matrix(0,nrow=length(incidence),ncol=lastday)

# Detectable probabilities over [lastday] days prior to test
detectable_props <- prop_det_use(1:lastday, pars, obs_model="gumbel", additional_detect=TRUE)
p_a_sums <- sapply(1:lastday, function(a) sum(p_a_use(x, a, pars, obs_model="gumbel", log_prob=FALSE,
                                                      ct=TRUE, additional_detect=TRUE)))

# Create matrices of relative frequencies of observing age of infection and viral load by day of testing:
for (i in 2:(dim(Ct_res)[1])) {
  past_inc <- incidence[(i-1):(max(i-lastday,1))]
  days <- 1:length(past_inc)
  age_res[i,days] <- past_inc*detectable_props[days]
  Ct_res[i,] <- sapply(X=x, FUN=function(xval) sum(age_res[i,days]*sapply(X=days, 
                                                                          FUN=function(a) p_a_use(xval, a, pars, obs_model="gumbel", log_prob=FALSE,
                                                                                                  ct=TRUE, additional_detect=TRUE)/p_a_sums[a])))
}


# Summaries of age and VL distributions:
age_res_std <- age_res/apply(age_res, 1, sum, na.rm=TRUE)
age_mean <- apply(age_res_std, 1, function(res) sum(res*(1:lastday), na.rm=TRUE))
age_res_std_csum <- t(apply(age_res_std, 1, FUN=function(res) cumsum(res)))
age_median <- apply(age_res_std_csum, 1, FUN=function(res) min((1:lastday)[res >= 0.5]))
plot(x=0:199, y=age_mean, type="l", lwd=2, col="blue")
lines(x=0:199, y=age_median, lwd=2, col="green", lty=2)

Ct_res_std <- Ct_res/apply(Ct_res, 1, sum, na.rm=TRUE)
Ct_mean <- apply(Ct_res_std, 1, function(res) sum(res*x, na.rm=TRUE))
Ct_mean <- ifelse(Ct_mean>=pars["intercept"], NA, Ct_mean)
Ct_res_std_csum <- t(apply(Ct_res_std, 1, FUN=function(res) cumsum(res)))
Ct_median <- apply(Ct_res_std_csum, 1, FUN=function(res) min(x[res >= 0.5]))
plot(x=0:199, y=Ct_mean, type="l", lwd=2, col="blue")
lines(x=0:199, y=Ct_median, lwd=2, col="green", lty=2)


#### Setting Testing Days to Highlight in Figures: ####
vertvals <- c(50,75,100,125,150)
GRs <- format(round(GR$GR[GR$Days %in% vertvals],3), nsmall=3)


#### Panel A: ####
DF.a <- data.frame(Days=c(0:199,GR$Days,GR$Days), Value=c(incidence,GR$GR/10+.01,GR$DailyGR/10+.01), 
                   Type=c(rep("Incidence",length(incidence)),rep("35-Day Growth Rate",length(GR$GR)),rep("1-Day Growth Rate", length(GR$DailyGR))))
DF.a$Typef <- factor(DF.a$Type, levels=c("Incidence", "1-Day Growth Rate", "35-Day Growth Rate"))

p.a <- ggplot(data=DF.a, 
              aes(x=Days,y=Value, group=Typef, color=Typef)) + geom_line(size=1.5) + 
  theme_light() +
  scale_color_manual(values=cbbPalette[c(1,7,2)]) + 
  scale_y_continuous(name="Per Capita Incidence", sec.axis=sec_axis(trans=~.*10-.1, name="Growth Rate")) + 
  labs(title="A. Per Capita Incidence and Growth Rate", x="Days Since Outbreak") + 
  theme(legend.title=element_blank(), legend.position="top",
        axis.title.x=element_blank(), axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        panel.grid.minor.x=element_blank()) + 
  geom_vline(xintercept=vertvals, col=cbbPalette[8], lty="dashed") + 
  scale_x_continuous(limits=c(25,175), breaks=seq(25,175,by=25))

#### Panel B: ####
set.seed(3712)
times <- sort(sample(x=0:199, size=sample_n, replace=TRUE, prob=incidence), decreasing=TRUE)
DF.b <- data.frame()
for (n in 1:sample_n) {
  Person <- n
  Days <- times[n] + 1:lastday
  Cts <- rx(pars, obs_model="gumbel", additional_detect=TRUE, a.range=1:lastday, corr=0.9)
  DF.b <- rbind(DF.b, data.frame(Person, Days, Cts))
}

p.b <- ggplot(data=DF.b) + 
  geom_tile(aes(x=Days,y=as.factor(Person),fill=Cts)) + 
  scale_fill_gradientn(name="Ct Values",
                       colours = c("white","#00468bff","#ffdd55ff","#ff751a","#ed0000ff","#ed0000ff","#ed0000ff"),
                       breaks=seq(0,40,by=8), labels=seq(0,40,by=8), na.value="white",
                       trans="reverse", limits=c(42,-2)) + 
  geom_vline(xintercept=vertvals, col=cbbPalette[8], lty="dashed") + 
  scale_x_continuous(limits=c(25,175), breaks=seq(25,175,by=25)) +
  theme_light() +
  theme(axis.text.y=element_blank(),
        axis.line.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title.x=element_blank(), 
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        panel.grid=element_blank(),
        legend.position=c(0.9,0.5),
        plot.background=element_blank()) +
  labs(title="B. Individual Cycle Threshold Trajectories", x="Days Since Outbreak", y="Individuals")

#### Panel C: ####
AgeDF <- age_res[1+vertvals,]
rownames(AgeDF) <- as.character(vertvals)
AgeDF2 <- NULL
for (i in 1:(dim(AgeDF)[1])) { ## Assigns [sample_n] points to test days & times of infection by coarsening cumulative distributions
  probs <- AgeDF[i,]*sample_n
  cum.probs <- cumsum(probs)
  cum.probs.r <- floor(cum.probs+.5)
  NumPerAge <- diff(c(0,cum.probs.r))
  ages <- rep(1:35, times=NumPerAge)
  cum.ages <- cumsum(ages)
  AgeDF2 <- rbind(AgeDF2, data.frame(TestDay=rep(vertvals[i], length(ages)), TSI=ages))
}
ADFmeans <- apply(age_res, 1, function(x) sum(x*(1:lastday), na.rm=TRUE)/sum(x, na.rm=TRUE))
ADFmeans <- data.frame(TestDay=0:199, MeanAge=ADFmeans)
ADFmedians <- data.frame(TestDay=0:199, MedianAge=age_median)

NumShifts <- length(vertvals)
age_distns <- age_res_std[vertvals+1,]
AgeDF.hists <- NULL
for (i in 1:NumShifts) {
  AgeDF.hists <- rbind(AgeDF.hists, data.frame(Pos=NumShifts-i, TestDay=vertvals[i], TSI=1:35, 
                                               Distn=age_distns[i,]))
}

p.c <- ggplot(data=AgeDF2, aes(x=TestDay,y=TSI)) +
  geom_jitter(height=.2, width=4, color=cbbPalette[1]) + 
  theme_light() + 
  labs(title="C. Time Since Infection", y="Days Since Infection", x="Days Since Outbreak") + 
  theme(legend.position="bottom",
        axis.title.x=element_blank(), axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        panel.grid.minor.x=element_blank(), legend.box="horizontal") + 
  geom_vline(xintercept=vertvals, col=cbbPalette[8], lty="dashed") + 
  scale_x_continuous(limits=c(25,175), breaks=seq(25,175,by=25)) + 
  geom_line(data=ADFmedians, aes(x=TestDay, y=MedianAge), color=cbbPalette[6], 
            linetype=1, inherit.aes=FALSE, size=1.5) + 
  guides(shape = FALSE, 
         color = guide_legend(override.aes=list(values=c("Median"=cbbPalette[6],"P25"=cbbPalette[6],"P75"=cbbPalette[6]),
                                                linetype=c(1,NA,NA), 
                                                shape=c(17,18,18)))) +
  coord_cartesian(ylim=c(-6,lastday)) + scale_y_continuous(expand=expansion(add=c(0,1)), breaks=seq(0,lastday,by=5))
for (i in 0:(NumShifts-1)) {
  ph <- ggplot() + geom_histogram(data=AgeDF.hists[AgeDF.hists$Pos==i,], aes(x=TSI, weight=Distn), breaks=seq(0.5,35.5,by=5)) +
    theme_minimal() + scale_x_continuous(breaks=seq(3,33,by=5), minor_breaks=NULL, expand=expansion(add=c(0,0))) +
    scale_y_continuous(expand=expansion(add=c(0,0))) + 
    theme(axis.title=element_blank(), 
          axis.text=element_blank(),
          # axis.ticks.length=unit(0, "pt"),
          axis.ticks=element_blank(),
          panel.grid.major=element_blank(), panel.grid.minor=element_blank())
  xval <- min(AgeDF.hists$TestDay[AgeDF.hists$Pos==i])
  p.c <- p.c + annotation_custom(grob=ggplotGrob(ph), xmin=xval-12, xmax=xval+12, 
                                 ymin=-6, ymax=2)
}

p.c.violin <- ggplot(data=AgeDF2, aes(x=TestDay,y=TSI)) +
  geom_jitter(height=.2, width=4, color=cbbPalette[1]) + 
  geom_violin(aes(x=TestDay, y=TSI, group=TestDay), scale="width",
              fill="grey70", alpha=0.75, color=NA) +
  theme_light() + 
  labs(title="C. Time Since Infection", y="Days Since Infection", x="Days Since Outbreak") + 
  theme(legend.position="bottom",
        axis.title.x=element_blank(), axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        panel.grid.minor.x=element_blank(), legend.box="horizontal") + 
  geom_vline(xintercept=vertvals, col=cbbPalette[8], lty="dashed") + 
  scale_x_continuous(limits=c(25,175), breaks=seq(25,175,by=25)) + 
  geom_line(data=ADFmedians, aes(x=TestDay, y=MedianAge), color=cbbPalette[6], 
            linetype=1, inherit.aes=FALSE, size=1.5) + 
  guides(shape = FALSE, 
         color = guide_legend(override.aes=list(values=c("Median"=cbbPalette[6],"P25"=cbbPalette[6],"P75"=cbbPalette[6]),
                                                linetype=c(1,NA,NA), 
                                                shape=c(17,18,18)))) +
  coord_cartesian(ylim=c(0,lastday)) + scale_y_continuous(expand=expansion(add=c(0,1)), breaks=seq(0,lastday,by=5))

#### Panel D: ####
CtDF <- Ct_res[1+vertvals,]
CtDF2 <- NULL
for (i in 1:(dim(CtDF)[1])) { ## Assigns [sample_n] points to VLs and test days by coarsening cumulative distribution
  probs <- CtDF[i,]*sample_n
  cum.probs <- cumsum(probs)
  cum.probs.r <- floor(cum.probs+.5)
  NumPerCt <- diff(c(0,cum.probs.r))
  Cts <- rep(x, times=NumPerCt)
  CtDF2 <- rbind(CtDF2, data.frame(TestDay=rep(vertvals[i], length(Cts)), Ct=Cts))
}

CtMeds <- data.frame(Days=0:199, Median=Ct_median)
CtMeds$Value <- -1*DF.h.2$Median/100+.3

Ct_distns <- Ct_res_std[vertvals+1,]
CtDF.hist <- NULL
for (i in 1:NumShifts) {
  CtDF.hist <- rbind(CtDF.hist, data.frame(Pos=NumShifts-i, TestDay=vertvals[i], Ct=x, 
                                           Distn=Ct_distns[i,]))
}
CtDF.hist$Ctrd <- ceiling(CtDF.hist$Ct) ## Coarsens viral load distribution by rounding each Ct up (approximates Ct discreteness?)
CtDF.hist2 <- as.data.frame(CtDF.hist %>% group_by(Pos, TestDay, Ctrd) %>% summarize(Distn=sum(Distn)))

p.d <- 
  ggplot() +
  geom_jitter(data=CtDF2, aes(x=TestDay, y=Ct), inherit.aes=FALSE, height=0, width=4, color=cbbPalette[1]) +
  geom_line(data=CtMeds, aes(x=Days, y=Median), inherit.aes=FALSE, size=1.5, color=cbbPalette[4]) +
  theme_light() +
  labs(tilte="D. Cycle Threshold Distribution", y="Cycle Threshold (Ct) Value", x="Days Since Outbreak") +
  scale_y_continuous(name="Cycle Threshold (Ct) Value", breaks=seq(10,40,by=5),
                     trans="reverse", expand=expansion(add=c(0,0))) +
  coord_cartesian(ylim=c(47,10)) +
  theme(legend.position="top",
        panel.grid.minor.x=element_blank()) + 
  geom_vline(xintercept=vertvals, col=cbbPalette[8], lty="dashed") + 
  scale_x_continuous(limits=c(25,175), breaks=seq(25,175,by=25)) 
for (i in 0:(NumShifts-1)) { ## Inset histograms
  ph <- ggplot() + geom_histogram(data=CtDF.hist2[CtDF.hist2$Pos==i,], aes(x=Ctrd, weight=Distn), breaks=c(-0.5,seq(4.5,40.5,by=4)), 
                                  inherit.aes=FALSE) + 
    theme_minimal() + scale_x_continuous(breaks=seq(2,38,by=4), minor_breaks=NULL, expand = c(0,0), trans="reverse")  +
    scale_y_continuous(expand=c(0,0)) +
    theme(axis.title=element_blank(), 
          axis.text=element_blank(),
          axis.ticks.length=unit(0, "pt"),
          panel.grid.major=element_blank(), panel.grid.minor=element_blank())
  xval <- vertvals[NumShifts-i]
  p.d <- p.d + annotation_custom(grob=ggplotGrob(ph), xmin=xval-15, xmax=xval+15, 
                                 ymin=-47, ymax=-40)
}

p.d.violin <- 
  ggplot() +
  geom_jitter(data=CtDF2, aes(x=TestDay, y=Ct), inherit.aes=FALSE, height=0, width=4, color=cbbPalette[1]) +
  geom_line(data=CtMeds, aes(x=Days, y=Median), inherit.aes=FALSE, size=1.5, color=cbbPalette[4]) +
  geom_violin(data=CtDF2, aes(x=TestDay, y=Ct, group=TestDay), scale="width",
              fill="grey70", alpha=0.75, color=NA) +
  theme_light() +
  labs(tilte="D. Cycle Threshold Distribution", y="Cycle Threshold (Ct) Value", x="Days Since Outbreak") +
  scale_y_continuous(name="Cycle Threshold (Ct) Value", breaks=seq(10,40,by=5),
                     trans="reverse", expand=expansion(add=c(0,0))) +
  coord_cartesian(ylim=c(42,10)) +
  theme(legend.position="top",
        panel.grid.minor.x=element_blank()) + 
  geom_vline(xintercept=vertvals, col=cbbPalette[8], lty="dashed") + 
  scale_x_continuous(limits=c(25,175), breaks=seq(25,175,by=25)) 

ggsave(filename="figs/fig1.png",
       plot=p.a+labs(title=NULL, tag="A")+theme(plot.tag.position="topleft", legend.position=c(.85,.75), legend.direction="vertical") +
         p.b+labs(title=NULL, tag="B")+theme(plot.tag.position="topleft") +
         p.c+labs(title=NULL, tag="C")+theme(plot.tag.position="topleft") +
         p.d+labs(title=NULL, tag="D")+theme(plot.tag.position="topleft") +
         plot_layout(guides="keep", ncol=1, nrow=4),
       width=8, height=12, units="in", dpi=400)

ggsave(filename="figs/fig1_violins.png",
       plot=p.a+labs(title=NULL, tag="A")+theme(plot.tag.position="topleft", legend.position=c(.85,.75), legend.direction="vertical") +
         p.b+labs(title=NULL, tag="B")+theme(plot.tag.position="topleft") +
         p.c.violin+labs(title=NULL, tag="C")+theme(plot.tag.position="topleft") +
         p.d.violin+labs(title=NULL, tag="D")+theme(plot.tag.position="topleft") +
         plot_layout(guides="keep", ncol=1, nrow=4),
       width=8, height=12, units="in", dpi=400)
