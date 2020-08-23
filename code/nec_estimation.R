########## 
##########
# This code contains the analysis for the fitting of the NEC parameter as a
# part of the larger analysis of the effect of microplastics on daphnia
##########
##########
# AUTHOR: Cole B. Brookson
# DATE OF CREATION: 2020-07-04
##########
##########

library(here)
library(tidyverse)
library(rstan)
library(gdata)
library(bayesplot)
library(brms)
library(beepr)

### read in data
dir = 'C:/Users/brookson/Documents/Github/Schur-etal-Data' #private data repo

offspring = 
  read_csv(
    'C:/Users/brookson/Documents/Github/Schur-etal-Data/offspring_data.csv',
    guess_max = 20000)
offspring = offspring %>% 
  select(-`X1`)

### create data groupings
str(offspring)
offspring$generation = as.factor(offspring$generation)
levels(offspring$generation)[levels(offspring$generation)== 'F3'] = 3

for(i in unique(offspring$generation)) {
  
  temp = offspring[which(offspring$generation == i),]
  assign(paste('generation_f', i, sep = ''), temp)
  rm(temp)
} 

### pull the data and get into proper format
#note the response is the proportion of surviving/original
nec_fit_data_f0 = generation_f0 %>% 
  filter(treatment %in% c('HFC', 'PS_400', 'PS_2000', 'PS_10000'),
         day == 21) 

nec_fit_data_f0$treatment = as.factor(nec_fit_data_f0$treatment)
levels(nec_fit_data_f0$treatment)
levels(nec_fit_data_f0$treatment) = c('0', '10000', '2000', '400')
nec_fit_data_f0$treatment = 
  as.numeric(levels(nec_fit_data_f0$treatment))[nec_fit_data_f0$treatment]

nec_fit_data_f0 = nec_fit_data_f0 %>% 
  select(treatment, survival) %>% 
  rename(x = treatment, y = survival)

#put into a list
N = nrow(nec_fit_data_f0)
x = nec_fit_data_f0$x
y = (nec_fit_data_f0$y/100)*20
n = rep(20,100)
nec_fit_data_f0_sample = list(N=N,
                              n=n,                                              
                              x=x,
                              y=y)
### read in model
nec_model = stan_model(file = here('./code/stan_files/nec_estimation.stan'))

### fit model
nec_fit = sampling(nec_model, 
                   data = nec_fit_data_f0_sample, 
                   seed = 12); beep(3)

### diagnose model fit
nec_fit_summ = print(nec_fit, 
                     pars=c('a', 'b', 'g'),
                     probs=c(0.1, 0.5, 0.9), digits = 3)
nec_fit_params = c('a', 'b', 'g')

### model checks
nec_fit_output = rstan::extract(nec_fit_summ,
                                permuted=TRUE,
                                include=TRUE)
#look at the parameters
nec_fit_trace = stan_trace(nec_fit, nec_fit_params)
nec_fit_dens = mcmc_dens(nec_fit, nec_fit_params)
nec_fit_overlay = mcmc_dens_overlay(nec_fit, nec_fit_params)
nec_fit_violin = mcmc_violin(nec_fit, nec_fit_params, probs = c(0.1, 0.5, 0.9))
nec_fit_pairs = mcmc_pairs(nec_fit, nec_fit_params)

nec_fit_trace
nec_fit_dens
nec_fit_overlay
nec_fit_violin
nec_fit_pairs
# Data; density
y_rep = Output3$y_rep
y_rep = data.frame(y_rep)
length(y_rep)
y_rep_P = y_rep[-c(21:40)]
y_rep_P = as.matrix(y_rep_P)
y_rep_H = y_rep[-c(1:20)]
y_rep_H = as.matrix(y_rep_H)
y_df = as.data.frame(y)
y_P = y_df[-c(2)]
y_P = as.vector(y_P$V1,mode = "numeric")
y_H = y_df[-c(1)]
y_H = as.vector(y_H$V2,mode = "numeric")
color_scheme_set("green")
P_plot = ppc_dens_overlay(y_P, y_rep_P[1:200,])+
  theme_minimal()+
  xlab("Abundance")+
  ylab("Density")+
  labs(title="Empirical vs. Estimated Distribution of Parasite Abundance Data",
       subtitle="Showing 200 Draws from the Posterior")
P_plot # How estimates line up with data for P
color_scheme_set("brightblue")
H_plot = ppc_dens_overlay(y_H, y_rep_H[1:200,])+
  theme_minimal()+
  xlab("Abundance")+
  ylab("Density")+
  labs(title="Empirical vs. Estimated Distribution of Host Immune Cell Abundance Data",
       subtitle="Showing 200 Draws from the Posterior")
H_plot # How estimates line up with data for H
# Data; time series
Z_df = data.frame(y_rep)
Z_means = Z_df %>% summarise_each(funs(mean))
# Parse df by P and H
Z_means_P = Z_means[-c(21:40)]
ncol(Z_means_P) # 2922; good
Z_means_H = Z_means[-c(1:20)]  
ncol(Z_means_H) # 2922
# Invert dfs
Z_means_P = t(Z_means_P)
Z_means_H = t(Z_means_H)
# Re-name columns 
colnames(Z_means_P) = c("post_means_P")
colnames(Z_means_H) = c("post_means_H")
head(y) # Recall that y is the data we gave the model (our stochastic data)
y_df = data.frame(y)
colnames(y_df) = c("P","H")
# Merge dfs
Mega_df = cbind(Z_means_P,Z_means_H,y_df)
# Add time steps
ts = seq(1:20)
ts = data.frame(ts)
# Final df
Mega_df = cbind(Mega_df,ts)
Post_By_Data_Plot = ggplot(Mega_df,aes(x=ts))+
  geom_line(aes(y=P,colour="springgreen4"))+ # Data P
  geom_line(aes(y=H,colour="cornflowerblue"))+ # Data H
  geom_line(aes(y=post_means_P,colour="lightgreen"),size=1.5)+ # Post. P
  geom_line(aes(y=post_means_H,colour="lightskyblue1"),size=1.5)+ # Post. H
  scale_x_continuous(breaks=seq(0,3000,by=500))+ # Cleaning up the axes
  scale_y_continuous(breaks=seq(0,400,by=50))+
  scale_color_manual(name="Legend",
                     labels=c("Host immune cell abundance from data",
                              "Estimated parasite abundance (av. over 2000 posterior draws)",
                              "Estimated host immune cell abundance (av. over 2000 posterior draws)",
                              "Parasite abundance from data"),
                     values=c("cornflowerblue","lightgreen","lightskyblue1","springgreen4"))+
  xlab("Time Step")+
  ylab("Population Abundance")+
  ggtitle("Posterior Estimates Plotted Against Data")+
  theme_minimal()
Post_By_Data_Plot

