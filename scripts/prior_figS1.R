require(tidyverse)
require(ggplot2)
require(ggthemes)
require(ggsci)
require(deSolve)
require(ggpubr)
require(patchwork)
require(triangle)
require(ggrepel)
require(ggpubr)
require(extraDistr)
require(mvtnorm)
library(moments)

cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

setwd("~/Documents/GitHub/ct_dynamics_preprint/")

source("code/priors.R")
source("~/Documents/GitHub/virosolver_paper/code/plot_funcs.R")
devtools::load_all("~/Documents/GitHub/virosolver")

n_samp <- 10000

#### Viral Load Kinetics Parameters: ####
pars1 <- c(beta = 0.1, tshift = 2, desired_mode = 2, viral_peak = 8.87848688565733, 
           obs_sd = 5.94042141408726, true_0 = -3, intercept = 40, LOD = 3, 
           limit_quantification = 2, incu = 5, t_switch = 15.788514698239, 
           level_switch = 3.68421266367764, wane_rate2 = 29.6102449407044, 
           prob_detect = 0.116087743649125, t_unit = 1)
pars <- c(beta = 0.1, tshift = 2, desired_mode = 2, viral_peak = 9, 
          obs_sd = 6, true_0 = -3, intercept = 40, LOD = 3, limit_quantification = 2, 
          incu = 5, t_switch = 16, level_switch = 3.7, 
          wane_rate2 = 30, prob_detect = 0.1, t_unit = 1)

pars <- c(beta = 0.1, tshift = 0, desired_mode = 5, viral_peak = 10.6, 
          obs_sd = 1.42, true_0 =40, intercept = 40, LOD = 3, limit_quantification = 2, 
          incu = 5, t_switch = 11.5, level_switch = 37, 
          wane_rate2 = 10000, prob_detect = 0.041, t_unit = 1)

sds <- c(beta=0.5,tshift=0,desired_mode=0,viral_peak=3,
         obs_sd=2,true_0=0,intercept=0,LOD=0,limit_quantification=0,
         incu=0,t_switch=3,level_switch=0.5,wane_rate2=5,prob_detect=0.1,
         t_unit=0)

parTab <- read.csv("~/Documents/GitHub/virosolver_paper/pars/partab_gp_model_start.csv")
parTab <- read.csv("~/Documents/GitHub/virosolver_paper/pars/partab_fitted.csv")

pars <- parTab$values
names(pars) <- parTab$names
pars["level_switch"] <- 35
pars["wane_rate2"] <- 50
pars["prob_detect"] <- 0.13
pars["t_switch"] <- 10
pars["mod"] <- 0.8

parTab$fixed <- 1
free_pars <- parTab[parTab$fixed == 0 & parTab$names != "prob", "names"]

sds <- c("beta"=0.25,"obs_sd"=1,
         "viral_peak"=1,"wane_rate2"=3,
         "t_switch"=3,"level_switch"=1,
         "prob_detect"=0.03)
sds <- c("beta"=0.25,"obs_sd"=0,
         "viral_peak"=0,"wane_rate2"=0,
         "t_switch"=0,"level_switch"=0,
         "prob_detect"=0)

ages <- 1:100

trajs <- matrix(0, n_samp, length(ages))
vl_trajs <- matrix(pars["true_0"], n_samp,length(ages))
sampd_pars <- matrix(0, n_samp, length(pars))


for(j in 1:n_samp){
  tmp_pars <- pars
  for(i in 1:length(pars)) {
    par_name <- names(pars)[i]
    if(par_name %in% free_pars) {
      if(par_name == "prob_detect"){
        beta1_mean <- pars["prob_detect"]
        beta1_sd <- sds["prob_detect"]
        beta_alpha <- ((1-beta1_mean)/beta1_sd^2 - 1/beta1_mean)*beta1_mean^2
        beta_beta <- beta_alpha*(1/beta1_mean - 1)
        tmp_pars[i] <- -1
        while(tmp_pars[i] < 0){
          tmp_pars[i] <- rbeta(1,beta_alpha,beta_beta)
        }
        
      } else {
        tmp_pars[par_name] <- -1
        while(tmp_pars[par_name] < 0){
          tmp_pars[par_name] <- rnorm(1,pars[par_name],sds[par_name])
        }
      }
    }
  }
  names(tmp_pars) <- names(pars)
  sampd_pars[j,] <- tmp_pars
  vl <- viral_load_func(tmp_pars, ages, FALSE)
  vl1 <- viral_load_func(tmp_pars, ages, FALSE)
  
  t_switch <-  tmp_pars["t_switch"] + tmp_pars["desired_mode"] + tmp_pars["tshift"]
  
  pre_tswitch <- which(ages < t_switch)
  post_tswitch <- which(ages >= t_switch)
  obs <- vl
  obs <- extraDistr::rgumbel(length(pre_tswitch),vl[pre_tswitch], tmp_pars["obs_sd"])
  obs[pre_tswitch] <- extraDistr::rgumbel(length(pre_tswitch),vl[pre_tswitch], tmp_pars["obs_sd"])
  obs[post_tswitch] <- extraDistr::rgumbel(length(post_tswitch),vl[post_tswitch], tmp_pars["obs_sd"]*tmp_pars["mod"])
  vl_trajs[j,] <- obs
  trajs[j,] <- prop_detectable(ages, tmp_pars,vl)
}


trajs1 <- t(apply(trajs, 2, function(x) quantile(x, c(0.025,0.25,0.5,0.75,0.975))))
trajs1 <- as.data.frame(trajs1)
colnames(trajs1) <- c("lower","lower_mid","median","upper_mid","upper")
trajs1$t <- ages

vl_trajs1 <- t(apply(vl_trajs, 2, function(x) quantile(x, c(0.025,0.25,0.5,0.75,0.975))))
vl_trajs1 <- as.data.frame(vl_trajs1)
colnames(vl_trajs1) <- c("lower","lower_mid","median","upper_mid","upper")
vl_trajs1$t <- ages

assumed_incu_period <- 5
p_detect <- ggplot(trajs1) + 
  geom_ribbon(aes(x=t,ymin=lower,ymax=upper),alpha=0.25) +
  geom_ribbon(aes(x=t,ymin=lower_mid,ymax=upper_mid),alpha=0.5) +
  geom_line(aes(x=t,y=median)) +
  geom_point(data=data.frame(ages=c(-5, 0, 5, 10, 15, 20, 25, 30, 35,40,45,50) + assumed_incu_period,
                             desired_probs=c(0,1,0.9,0.72,0.49,0.29,0.14,0.04,0,0,0,0)),
             aes(x=ages,y=desired_probs)) +
  xlab("Days since infection") +
  ylab("Proportion detectable") +
  scale_x_continuous(limits=c(0,100),breaks=seq(0,100,by=5)) +
  geom_vline(xintercept=5,linetype="dashed") +
  export_theme


p_vl <- ggplot(vl_trajs1) + 
  geom_ribbon(aes(x=t,ymin=lower,ymax=upper),alpha=0.25) +
  geom_ribbon(aes(x=t,ymin=lower_mid,ymax=upper_mid),alpha=0.5) +
  geom_line(aes(x=t,y=median)) +
  xlab("Days since infection") +
  scale_y_continuous(trans="reverse") +
  scale_x_continuous(limits=c(0,100),breaks=seq(0,100,by=5)) +
  coord_cartesian(ylim=c(40,0)) +
  ylab("Viral load (log10 RNA copies/ml)") +
  geom_vline(xintercept=5,linetype="dashed") +
  export_theme

p_vl/p_detect

trajs_melt <- reshape2::melt(trajs)
colnames(trajs_melt) <- c("samp","t","prop_detectable")

#write_csv(trajs1, "~/Documents/figS1_for_joel.csv")
#write_csv(trajs_melt, "~/Documents/figS1_draws_for_joel.csv")
