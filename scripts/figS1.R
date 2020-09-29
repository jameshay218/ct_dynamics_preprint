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

source("code/model_funcs.R")
source("code/simulation_functions.R")
source("code/odin_funcs.R")
source("code/priors.R")
source("code/sim_data_funcs.R")
source("code/fig1_code.R")

#### Viral Load Kinetics Parameters: ####
pars <- c(beta = 0.1, tshift = 2, desired_mode = 2, viral_peak = 9, 
          obs_sd = 6, true_0 = -3, intercept = 40, LOD = 3, limit_quantification = 2, 
          incu = 5, t_switch = 16, level_switch = 3.7, 
          wane_rate2 = 30, prob_detect = 0.1, t_unit = 1)


#### Figure S1: ####

S1.df <- data.frame(Day=rep(1:lastday))
S1.df$VL <- viral_load_func_asymp(pars, S1.df$Day, FALSE)
S1.df$Ct <- viral_load_func_asymp(pars, S1.df$Day, TRUE)
S1.df$Ct.LI <- qgumbel(.975, mu=S1.df$Ct, sigma=pars["obs_sd"], lower.tail=TRUE, log.p=FALSE)
S1.df$Ct.UI <- qgumbel(.025, mu=S1.df$Ct, sigma=pars["obs_sd"], lower.tail=TRUE, log.p=FALSE)
S1.df$VL.LI <- (pars["intercept"]-S1.df$Ct.LI)/log2(10)+pars["LOD"]
S1.df$VL.UI <- (pars["intercept"]-S1.df$Ct.UI)/log2(10)+pars["LOD"]
S1.df$PropDet <- prop_det_use(S1.df$Day, pars, obs_model="gumbel", additional_detect=TRUE)

set.seed(2013)
violdays <- seq(5,lastday,by=5)
S1.rx.df <- NULL
for (i in violdays) {
  Cts <- rx_a(10000, i, pars, obs_model="gumbel", additional_detect=TRUE)
  Cts <- Cts[!is.na(Cts)]
  S1.rx.df <- rbind(S1.rx.df, data.frame(Day=rep(i,length(Cts)), Cts=Cts))
}
S1.rx.df <- merge(S1.rx.df, S1.df[c("Day","PropDet")], by="Day", all.x=TRUE, all.y=FALSE)

s1.a <- ggplot(data=S1.df) +
  geom_violin(data=S1.rx.df, aes(x=Day, y=Cts, group=Day, fill=PropDet)) +
  scale_fill_gradient(name="Proportion of Infections Detectable:", low="blue", high="red",
                      breaks=seq(0,1,by=.25), labels=seq(0,1,by=.25), limits=c(-.05,1.05)) +
  geom_line(aes(x=Day, y=Ct), col=cbbPalette[4]) + geom_point(aes(x=Day, y=Ct), col=cbbPalette[4]) +
  theme_light() + 
  scale_x_continuous(breaks=seq(0,lastday,by=5)) +
  scale_y_continuous(breaks=seq(10,60,by=10), trans="reverse", expand=expansion(add=c(0,0)),
                     sec.axis=sec_axis(trans=~(.*(-1)+40)/log2(10)+3, 
                                       name=expression("Viral Load ("*log["10"]*" RNA copies per mL)"),
                                       breaks=seq(-3,12,by=3))) +
  theme(legend.position="bottom") +
  coord_cartesian(ylim=c(65,5)) +
  labs(x="Time Since Infection (Days)", y="Cycle Threshold (Ct) Value") + 
  # geom_errorbar(aes(x=Day, ymin=Ct.UI, ymax=Ct.LI), col=cbbPalette[4]) +
  geom_hline(yintercept=pars["intercept"], lty="dashed")
s1.a

s1.b <- ggplot(data=S1.df) +
  geom_line(aes(x=Day, y=PropDet, col=PropDet)) + geom_point(aes(x=Day, y=PropDet, col=PropDet)) +
  scale_color_gradient(name="Proportion of Infections Detectable:", low="blue", high="red",
                       breaks=seq(0,1,by=.25), labels=seq(0,1,by=.25), limits=c(-.05,1.05),
                       guide=NULL) +
  theme_light() +
  scale_x_continuous(breaks=seq(0,lastday,by=5)) +
  labs(x="Time Since Infection (Days)", y="Proportion of Infections Detectable")
s1.b

ggsave(filename="schematics/FigS1.png",
       plot=s1.a+theme(plot.tag.position="topleft")+labs(tag="A  ") + 
         s1.b+theme(plot.tag.position="topleft")+labs(tag="B  ") + 
         plot_layout(nrow=2, ncol=1),
       width=6.51, height=6.51, units="in", dpi=600)
