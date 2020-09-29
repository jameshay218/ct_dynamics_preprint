library(tidyverse)
library(patchwork)
library(EpiEstim)
library(zoo)
library(ggsci)
library(EpiNow2)

main_wd <- "~/Documents/GitHub/ct_inference_preprint/"
setwd(main_wd)

rerun <- FALSE

## Create an epiweek calendar
dates <- seq(as.Date("2020-01-01"),as.Date("2020-12-31"),by="1 day")
epiweeks <- lubridate::epiweek(dates)
epi_calendar <- tibble(date=dates,week=epiweeks)
epi_calendar <- epi_calendar %>% group_by(week) %>% mutate(first_day=min(date))

## NYT data
nyt_dat <- read_csv("data/us-states.csv")

nyt_dat <- nyt_dat %>% 
  filter(state=="Massachusetts") %>%
  group_by(state) %>%
  mutate(new_cases=cases-lag(cases, 1))

## Plot case counts
p1 <- nyt_dat %>% 
  ## Highlight anomalous cases
  mutate(new_cases1=ifelse(new_cases > 3500, lag(new_cases,1), new_cases)) %>%
  ## Get rolling mean
  mutate(roll_mean=rollmean(new_cases1, 7,fill=NA,align="right")) %>%
  ggplot()+ 
  geom_bar(aes(x=date,y=new_cases),stat="identity",alpha=0.25) +
  #geom_ribbon(aes(x=date,ymax=roll_mean,ymin=0),fill="grey70",alpha=0.5,col="grey20") +
  #scale_fill_npg() +
  scale_fill_manual(values=c("grey40","darkred")) +
  scale_y_continuous(expand=c(0,0),limits=c(0,5000)) +
  scale_x_date(limits=as.Date(c("2020-03-01", "2020-09-01"), "%Y-%m-%d"), breaks="7 days") +
  theme_classic()+
  theme(legend.position="none",
        panel.grid.minor=element_blank(),
        #axis.text.x=element_text(angle=45,hjust=1),
        axis.text.x=element_blank(),
        axis.title.x=element_blank(),
        axis.ticks.x=element_blank()) +
  xlab("Date") +
  ylab("New infections") +
  labs(tag="A")
p1

## Rt estimation on these case counts
rt_dat <- nyt_dat %>% 
  ungroup() %>%
  select(date, new_cases) %>%
  drop_na() %>%
  rename(I=new_cases,
         dates=date)

reporting_delay <- EpiNow2::bootstrapped_dist_fit(rlnorm(100, log(3), 0.5))
## Set max allowed delay to 30 days to truncate computation
reporting_delay$max <- 30

generation_time <- list(mean = EpiNow2::covid_generation_times[1, ]$mean,
                        mean_sd = EpiNow2::covid_generation_times[1, ]$mean_sd,
                        sd = EpiNow2::covid_generation_times[1, ]$sd,
                        sd_sd = EpiNow2::covid_generation_times[1, ]$sd_sd,
                        max = 30)

incubation_period <- list(mean = EpiNow2::covid_incubation_period[1, ]$mean,
                          mean_sd = EpiNow2::covid_incubation_period[1, ]$mean_sd,
                          sd = EpiNow2::covid_incubation_period[1, ]$sd,
                          sd_sd = EpiNow2::covid_incubation_period[1, ]$sd_sd,
                          max = 30)
reported_cases <- EpiNow2::example_confirmed[1:50]
reported_cases <- rt_dat %>% rename(confirm=I, date=dates)
if(rerun){
  estimates <- EpiNow2::epinow(reported_cases = reported_cases, 
                             generation_time = generation_time,
                             delays = list(incubation_period, reporting_delay),
                             horizon = 7, samples = 4000, warmup = 1000, 
                             cores = 4, chains = 4, verbose = TRUE, 
                             adapt_delta = 0.95)
  saveRDS(estimates,"results/ma_rt_fit.RData")

} else {
  estimates <- readRDS("results/ma_rt_fit.RData")
}

dat_infections <- estimates$estimates$summarised %>% filter(variable == "infections" & 
                                                              type=="estimate")
rt_dat <- estimates$estimates$summarised %>% filter(variable %in% c("R") & 
                                                      type=="estimate")
coeff <- 1000
p1 <- p1 + 
  geom_ribbon(data=dat_infections,aes(x=date,ymin=bottom,ymax=top),
              fill="red",alpha=0.5) +
  geom_line(data=dat_infections,aes(x=date,y=mean),col="red") +
  geom_ribbon(data=rt_dat,aes(x=date,ymin=bottom*coeff,ymax=top*coeff),
              fill="darkgreen",alpha=0.5) +
  geom_line(data=rt_dat,aes(x=date,y=mean*coeff)) +
  scale_y_continuous(expand=c(0,0),limits=c(0,5000),
                     sec.axis=sec_axis(~.*1/coeff, name="Rt"))


## Plot Ct data
bwh_data <- read_csv("data/BWH_COVID_Cts_deid_20200403-20200831.csv")

bwh_data_use <-  bwh_data %>% 
  filter(platform=="Panther" &
           first_pos %in% c(1,0)) %>%
  rename(date=coll_date) %>%
  left_join(epi_calendar)

bwh_data_week <- bwh_data_use %>% 
  group_by(date) %>% 
  summarize(median_ct=median(panther_Ct),
            skew_ct=moments::skewness(panther_Ct),
            n=n())

combined_dat <- bwh_data_week %>% 
  #rename(date=first_day) %>%
  full_join(estimates$estimates$summarised %>% as_tibble() %>% 
              filter(variable == "R") %>%
              mutate(date = date)
              ) %>%
  arrange(date) %>% 
  rename(mean_rt=mean) %>%
  select(date, median_ct,skew_ct, n, mean_rt) %>%
  filter(n >= 10)

combined_dat1 <- combined_dat %>% 
  pivot_longer(-c(date,n)) %>%
  left_join(epi_calendar) %>%
  group_by(first_day, name) %>% 
  summarize(mean=mean(value)) %>%
  pivot_wider(values_from=mean,names_from=name)
combined_dat1 <- combined_dat
  
median_cts <- smooth(combined_dat$median_ct)
skew_cts <- combined_dat$skew_ct
Rt <- combined_dat$mean_rt
ccf(median_cts, Rt,lag.max = 25)
ccf(skew_cts, Rt,lag.max = 25)


ggplot(combined_dat) +
  geom_line(aes(x=date,y=median_ct))

pA <- combined_dat1 %>%
  ggplot() +
  geom_smooth(aes(x=mean_rt,y=skew_ct),se=TRUE)+
  geom_point(aes(x=mean_rt,y=skew_ct),size=2) +
  geom_vline(xintercept=1,linetype="dashed") +
  theme_classic() +
  ylab("Skewness of Ct distribution") +
  xlab("Rt (posterior mean)") +
  labs(tag="C")
pA

pB <- combined_dat1 %>%
  ggplot() +
  geom_smooth(aes(x=mean_rt,y=median_ct),se=TRUE)+
  geom_point(aes(x=mean_rt,y=median_ct),size=2)+
  geom_vline(xintercept=1,linetype="dashed") +
  theme_classic() +
  scale_y_continuous(trans="reverse") +
  ylab("Median of Ct distribution") +
  xlab("Rt (posterior mean)") +
  labs(tag="D")
  
pC <- ggplot(combined_dat1 %>% rename(Rt=mean_rt)) +
  geom_point(aes(x=skew_ct,y=median_ct,col=Rt),alpha=0.9,size=2) +
  scale_color_gradient2(low="green",mid="blue",high="red",midpoint=1)+
  scale_y_continuous(trans="reverse") +
  theme_classic() +
  xlab("Skewness of Ct distribution") +
  ylab("Median of Ct distribution") +
  theme(legend.position=c(0.2,0.7)) +
  labs(tag="E")


p2 <- bwh_data_use %>% 
  ggplot() + 
  geom_violin(aes(x=first_day, y=panther_Ct,group=week),scale="width",
              fill="grey70",alpha=0.75,color=NA) +
  #geom_jitter(aes(x=first_day, y=panther_Ct,group=week),width=2,height=0,
  #            fill="grey70",size=0.25,alpha=0.25) +
  geom_dotplot(aes(x=first_day, y=panther_Ct,group=week),binaxis="y",
               binwidth=1,stackdir="center",binpositions="all",dotsize=0.25) +
  geom_smooth(data=bwh_data_use %>% group_by(first_day) %>% summarize(median_ct=median(panther_Ct)),
            aes(x=first_day,y=median_ct),col="blue",se=FALSE) +
  scale_y_continuous(trans="reverse",limits=c(45, 5),expand=c(0,0)) +
  geom_hline(yintercept=40,linetype="dashed") +
  theme_bw() +
  scale_x_date(limits=as.Date(c("2020-03-01", "2020-09-01"), "%Y-%m-%d"), breaks=unique(epi_calendar$first_day)) +
  theme_classic()+
  xlab("Date") +
  ylab("Ct value") +
  theme(legend.position="none",
        plot.title=element_blank(),
        panel.grid.minor =element_blank(),
        axis.text.x=element_text(angle=45,hjust=1)) +
  labs(tag="B")



main_p <- p1/p2/(pA|pB|pC)
main_p

ggsave("figs/brigham_data_daily_firstpos.png",main_p,
       width=8,height=8,units="in")
