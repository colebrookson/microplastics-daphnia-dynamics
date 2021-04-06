########## 
##########
# This code contains some data preparation and model fitting for an analysis
# investigating the effects of microplastic pollution on Daphnia magna
# presented in Brookson et al. ()
##########
##########
# AUTHOR: Cole B. Brookson
# DATE OF CREATION: 2021-04-06
##########
##########

# set-up =======================================================================

library(here)
library(tidyverse)
library(rstan)
library(gdata)
library(bayesplot)
library(brms)
library(beepr)
library(parallel)
library(R2WinBUGS)
library(deSolve)
detectCores()


dir = 'C:/Users/brookson/Documents/Github/Schur-etal-Data/' #private data repo

growth_data = read_csv(paste0(dir, 'growth_data_gen0.csv'))
reproduction_data = read_csv(paste0(dir, 'offspring_data.csv'), 
                             guess_max = 15000)
unique(reproduction_data$treatment)

# get the data required ========================================================

# keep only the LFC and generation 0
reproduction_data = reproduction_data %>% 
  dplyr::select(-`X1`) %>% 
  filter(treatment == 'PS_2000') 

# make the length into mm
growth_data = growth_data %>% 
  select(-`x1`) %>% 
  filter(treatment == 'MP2000') %>% #keep only ps2000 treatment
  filter(test == 1) %>% 
  filter(day < 23) %>%
  mutate(length_mm = length*0.001)  

# prep data ====================================================================

# create data matrix
l_y_obs_data = as.matrix(growth_data %>% 
                           filter(day < 22) %>% 
                           select(day, length_mm, replicate) %>% 
                           pivot_wider(id_cols = c(replicate, day, length_mm),
                                       names_from = replicate,
                                       values_from = length_mm) %>% 
                           select_if(~ !any(is.na(.))) %>% 
                           select(-day))

# identify replicates that didn't survive
reps_to_remove = reproduction_data %>% 
  select(day, replicate, offspring) %>% 
  group_by(day, replicate) %>% 
  filter(offspring == 'X')
unique(reps_to_remove$replicate)

# remove reps that didn't surviv
reproduction_data_x_removed = reproduction_data %>% 
  filter(!replicate %in% reps_to_remove$replicate)

# make numeric
reproduction_data_x_removed$offspring = 
  as.numeric(reproduction_data_x_removed$offspring)

# turn NA's into zeros
reproduction_data_x_removed$offspring[which(
  is.na(reproduction_data_x_removed$offspring))] = 0

# create data matrix
r_y_data = as.matrix(reproduction_data_x_removed %>% 
                       filter(day < 22) %>%  
                       pivot_wider(id_cols = c(replicate, day, offspring),
                                   names_from = replicate,
                                   values_from = offspring, 
                                   values_fn = max))
r_y_data = r_y_data[,2:ncol(r_y_data)]
r_y_data_cumul = matrix(nrow = nrow(r_y_data), ncol = ncol(r_y_data))
for(i in 1:nrow(r_y_data_cumul)) {
  for(j in 1:ncol(r_y_data_cumul)){
    
    r_y_data_cumul[i,j] = sum(r_y_data[1:i,j])
    
  }
}

#declare all variables
funct = function(t, y, p){
  cq = y[1]
  
  d_cq = ke*(2000-cq)
  
  return(list(d_cq))
}
ke = 0.5
parms = c(ke)

N0 = 0
TT = seq(1,22,1) 
results = lsoda(N0,TT,funct,parms)

N_obs = 10
N_mis = 12
ii_obs = c(1, 3, 6, 8, 10, 13, 15, 17, 20, 22)
ii_mis = c(2, 4, 5, 7, 9, 11, 12, 14, 16, 18, 19, 21) 
ll_init = 0.2
cq = results[,2]
l_y_obs = l_y_obs_data[,1]#this is taking first replicate (column)
ts = 1:nrow(r_y_data)
r_y = r_y_data_cumul[,1] #this is taking first replciate (column )

gen_0_ps2000_data = list(
  N_obs = N_obs, 
  N_mis = N_mis,
  ll_init = array(ll_init),
  #cq_init = array(cq_init),
  ii_obs = ii_obs,
  ii_mis = ii_mis, 
  l_y_obs = l_y_obs,
  ts = ts,
  r_y = r_y,
  cq = cq)

# fit model ====================================================================
gen_0_ps2000_onerep_fit = stan(file = 
                                here('./code/stan_files/organism_costs_model_gen0_ps2000.stan'),
                              data = gen_0_ps2000_data,
                              chains = 4,
                              cores = 8,
                              warmup = 5000,
                              iter = 10000,
                              seed = 12,
                              verbose = TRUE,
                              #open_progress = TRUE,
                              control = list(adapt_delta = 0.999,
                                             max_treedepth = 12)); beep(3)
saveRDS(gen_0_ps2000_onerep_fit, 
        here('/output/intermediate-objects/gen_0_ps2000_onerep_fit.RDS'))
gen_0_ps2000_onerep_fit = readRDS(here('/output/intermediate-objects/gen_0_ps2000_onerep_fit.RDS'))
print(gen_0_ps2000_onerep_fit, 
                         pars=c("theta_ll[1]", "cstar", "NEC",
                                "Lp", "Rm", "Lm", "tau_l", "tau_r"),
                         probs=c(0.1, 0.5, 0.9), digits = 3)

parms = c("theta_ll[1]", 
          #"theta_cq[1]", 
          "cstar", "NEC",
          "Lp", "Rm", "Lm", "tau_l", "tau_r")
gen_0_ps2000_output = rstan::extract(gen_0_ps2000_onerep_fit,
                                    permuted=TRUE,include=TRUE)

gen_0_ps2000_onerep_fit_Trace = stan_trace(gen_0_ps2000_onerep_fit,parms)
gen_0_ps2000_onerep_fit_Dens = mcmc_dens(gen_0_ps2000_onerep_fit,parms)
gen_0_ps2000_onerep_fit_Overlay = mcmc_dens_overlay(gen_0_ps2000_onerep_fit,parms)
gen_0_ps2000_onerep_fit_Violin = mcmc_violin(gen_0_ps2000_onerep_fit,parms,
                                              probs = c(0.1, 0.5, 0.9))
gen_0_ps2000_onerep_fit_Pairs = mcmc_pairs(gen_0_ps2000_onerep_fit,parms)
gen_0_ps2000_onerep_fit_Trace
gen_0_ps2000_onerep_fit_Dens
gen_0_ps2000_onerep_fit_Overlay
gen_0_ps2000_onerep_fit_Violin
gen_0_ps2000_onerep_fit_Pairs

# Data; density
y_rep_ps2000 = gen_0_ps2000_output$r_y_rep
y_rep_ps2000 = data.frame(y_rep_ps2000)
y_rep_ps2000 = as.matrix(y_rep_ps2000)
y_df_ps2000 = data.frame(r_y)
y_df_ps2000 = as.vector(y_df_ps2000$r_y, mode = 'numeric')
color_scheme_set("green")
P_plot_ps2000 = ppc_dens_overlay(as.vector(y_df_ps2000),
                                  y_rep_ps2000[1:200,])+
  theme_bw()+
  theme(
    panel.grid = element_blank(),
    legend.position = 'none'
  )+
  labs(title="Empirical vs. Estimated Distribution \nof Cumulative Reproduction Data (ps2000 Treatment)",
       subtitle="Showing 200 Draws from the Posterior")
ggsave(here('output/model-fit-figs/ps2000_dens_overlay.png'), P_plot_ps2000)

# Data; time series
Z_df_ps2000 = data.frame(gen_0_ps2000_output$r_y_rep)
Z_df_ps2000_means = Z_df_ps2000 %>% summarize(across(`X1`:`X22`, mean))
Z_df_ps2000_means = t(Z_df_ps2000_means)
# Parse df by P and H
# Invert dfs
# Re-name columns 
# Merge dfs
Mega_df_ps2000 = cbind(Z_df_ps2000_means,
                        data.frame(r_y), ts = seq(1,22,1))
# Add time steps
# Final df
#Mega_df = cbind(Mega_df,ts)
Post_By_Data_Plot_ps2000 = ggplot(Mega_df_ps2000,aes(x=ts))+
  geom_line(aes(y=Z_df_ps2000_means), colour="black", size = 1.5)+ # Data P
  geom_point(aes(y=r_y), shape = 21, fill="cornflowerblue", size = 3.2)+ 
  xlab("Time Step")+
  ylab("Cumulative Reproduction")+
  ggtitle("Posterior Estimates Plotted Against Data \n(ps2000)")+
  theme_bw()+
  theme(
    panel.grid = element_blank(),
    legend.position = 'none',
    axis.title = element_text(size = 17),
    axis.text = element_text(size = 12),
    plot.title = element_text(size = 20)
  )
ggsave(here('output/model-fit-figs/ps2000_model_data.png'), 
       Post_By_Data_Plot_ps2000)


