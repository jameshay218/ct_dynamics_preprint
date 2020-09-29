
setwd("~/Documents/GitHub/ct_dynamics_preprint/")

#### Pulls in functions and data from the Fig 1 generation (note: takes some time to run): ####
source("scripts/fig1.R")

#### Parameters for the alternative observation model: ####
pars_alt <- c(beta = 0.1, tshift = 2, desired_mode = 2, viral_peak = 6,
              obs_sd = 4, true_0 = -3, intercept = 40, LOD = 3, limit_quantification = 2,
              incu = 5, t_switch = 16, level_switch = 3,
              wane_rate2 = 30, prob_detect = 0.05, t_unit = 1)


GR <- GR[,c("Days","GR")]

days <- 1:lastday

TSI_proportions <- age_res_std
TSI_medians <- age_median
TSI_skews <- apply(TSI_proportions, 1, 
                   function(x) sum(x*(days-sum(x*days))^3)/(sum(x*((days-sum(x*days))^2))^(3/2)))

Ct_vals <- seq(0,40,by=0.01)

Ct_1_proportions <- Ct_res_std
Ct_1_medians <- Ct_median
Ct_1_skews <- apply(Ct_1_proportions, 1, 
                   function(x) sum(x*(Ct_vals-sum(x*Ct_vals))^3)/(sum(x*((Ct_vals-sum(x*Ct_vals))^2))^(3/2)))

#### Calculating Cts for alternate observation model: ####
Ct_2_res <- matrix(0,nrow=length(incidence), ncol=length(Ct_vals))
detectable_props_2 <- prop_det_use(days, pars_alt, obs_model="gumbel", additional_detect=TRUE)
p_a_sums_2 <- sapply(days, function(a) sum(p_a_use(Ct_vals, a, pars_alt, obs_model="gumbel", log_prob=FALSE,
                                                      ct=TRUE, additional_detect=TRUE)))
for (i in 2:(dim(Ct_2_res)[1])) {
  past_inc <- incidence[(i-1):(max(i-lastday,1))]
  days_int <- 1:length(past_inc)
  age_res[i,days_int] <- past_inc*detectable_props_2[days_int]
  Ct_2_res[i,] <- sapply(X=Ct_vals, FUN=function(xval) sum(age_res[i,days_int]*sapply(X=days_int, 
                                                                          FUN=function(a) p_a_use(xval, a, pars_alt, obs_model="gumbel", log_prob=FALSE,
                                                                                                  ct=TRUE, additional_detect=TRUE)/p_a_sums_2[a])))
}

Ct_2_proportions <- Ct_2_res/apply(Ct_2_res, 1, sum, na.rm=TRUE)
Ct_2_res_std_csum <- t(apply(Ct_2_proportions, 1, FUN=function(res) cumsum(res)))
Ct_2_medians <- apply(Ct_2_res_std_csum, 1, FUN=function(res) min(x[res >= 0.5]))
Ct_2_skews <- apply(Ct_2_proportions, 1, 
                    function(x) sum(x*(Ct_vals-sum(x*Ct_vals))^3)/(sum(x*((Ct_vals-sum(x*Ct_vals))^2))^(3/2)))


DF.S2 <- data.frame(Days=rep(GR[,"Days"],6),GR=rep(GR[,"GR"],6),
                Property=rep(c("Median","Skew"),each=3*length(GR[,"Days"])),
                Variable=rep(rep(c("Time since infection", "Ct value", "Ct value"),each=length(GR[,"Days"])),2),
                Model=rep(rep(c("1","1","2"),each=length(GR[,"Days"])),2),
                Group=rep(rep(c(1,2,3),each=length(GR[,"Days"])),2),
                Value=c(TSI_medians[-1],Ct_1_medians[-1],Ct_2_medians[-1],
                        TSI_skews[-1],Ct_1_skews[-1],Ct_2_skews[-1]))

DF.S2 <- DF.S2[DF.S2$Days >= 25,]

p.S2 <- ggplot(data=DF.S2, 
               aes(x=GR, y=Value, group=Group, color=Variable, alpha=Model)) + 
  geom_path(lwd=1.5) + 
  scale_x_reverse(name="Growth Rate", breaks=seq(-.1,.1,by=.025), limits=c(.1,-.1)) +
  scale_y_continuous(name="Distribution Measure") +
  scale_alpha_manual(labels=NULL, values=c(1,0.5)) +
  scale_color_manual(labels=c("Ct value", NULL,"Time since infection"), values=cbbPalette[c(4,6,6)]) +
  guides(alpha=FALSE) +
  geom_point(data=DF.S2[DF.S2$Days %in% vertvals,], 
             aes(x=GR, y=Value, group=Group, shape=as.factor(Days)), size=2.5, color="black",
             bg=cbbPalette[8]) +
  theme_light() +
  scale_shape_manual(name="Test Day", values=21:25) + 
  facet_wrap(~Property, ncol=1, scales="free")
p.S2

ggsave(filename="figs/figS2.png",
       plot=p.S2,
       width=6, height=6, units="in", dpi=400)

